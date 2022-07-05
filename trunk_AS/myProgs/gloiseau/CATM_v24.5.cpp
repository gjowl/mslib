#include <iostream>
#include <fstream>
#include <unistd.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"
#include "HelixDimer.h"
#include "HelixDimerCluster.h"
#include "hbondInfo.h"
#include "catmFunctions.h"
#include "catmOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

string programName = "CATM";
string programDescription = "This program repacks two helices from a set starting position: electrostatics are on in this version";
string programAuthor = "Benjamin K. Mueller, Sabareesh Subramaniam; functions, options, and classes organized into separate files by Gilbert Loiseau";
string programVersion = "0.0.24.5";
string programDate = "21 April 2022";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

static SysEnv ENV;

struct compareHelixDimers {
	bool operator () (HelixDimer* lhs, HelixDimer* rhs) {return *lhs < *rhs;}
};

int main(int argc, char *argv[]) {

	time(&startTime);

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

  	Options opt = CATMParseOptions(argc, argv, defaults);
  	if (opt.errorFlag) {
  		cerr << endl;
  		cerr << "The program terminated with errors:" << endl;
  		cerr << endl;
  		cerr << opt.errorMessages << endl;
  		cerr << endl;
  		cerr << opt.OPerrors << endl;

  		CATMUsage();
  		exit(1);
  	}

	string logFile;
	setupOutputDirectory(opt, logFile);

	ofstream fout;
	fout.open(logFile.c_str());

	if(!fout.is_open()) {
		cerr << "Unable to open " << logFile << endl;
		exit(0);
	}

	printOptions(opt, fout);

	/******************************************************************************
	 *                     === READ IN GEOMETRY FILE ===
	 ******************************************************************************/
	vector<string> fileVec;
	readGeometryFile(opt.helixGeoFile,fileVec);

	// Import sequence, determine length
	string originalTMSeq = opt.fullSequence.substr(opt.startResNum-opt.fullSequenceStart,opt.endResNum - opt.startResNum + 1);
	unsigned int sequenceLength = originalTMSeq.length();
	// cant do sequences of size less than 4
	if(sequenceLength < 4) {
		cerr << "Sequence " << originalTMSeq << "is too small (should be >= 4 AA long)" << endl;
		exit(0);
	}

	string modelledTMSeq ="" ; // replace P with A
	string prolineMask = "";
	for(int i = 0; i < originalTMSeq.length(); i++) {
		char AA = originalTMSeq[i];
		if(AA == 'P') {
			modelledTMSeq += "A";
			prolineMask += "1";
		} else {
			modelledTMSeq += AA;
			prolineMask += "0";
		}
	}
	/**********************************************************************************
	*
	*    printProteinOutFile
	*
	**********************************************************************************/

	ofstream pout;
	string poutName  = opt.pdbOutputDir + "/" + opt.uniprotAccession + ".out";
	pout.open(poutName.c_str());
	if(!pout.is_open()) {
		cerr << "Unable to open " << poutName << endl;
		exit(0);
	}

	pout << "uniprotName " << opt.uniprotName << endl;
	pout << "uniprotAccession " << opt.uniprotAccession << endl;
	pout << "sequence " << opt.fullSequence << endl;
	pout << "sequenceStart " << opt.fullSequenceStart << endl;
	pout << "originalSeq " << originalTMSeq << endl;
	pout << "modelledSeq " << modelledTMSeq << endl;
	pout << "prolineMask " << prolineMask << endl;
	pout << "resStart " << opt.startResNum << endl;
	pout << "resEnd " << opt.endResNum << endl;
	pout << "tmStart " << opt.tmStart << endl;
	pout << "tmEnd " << opt.tmEnd << endl;
	pout << "CATMversion " << programVersion << endl;

	//cout << modelledTMSeq << endl;

	string sequence = convertToPolymerSequence(modelledTMSeq,1); // so that the 4th residue will be the middle one (35th) on the GLY 69 backbone
	PolymerSequence PS(sequence);

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", true);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << sequence << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	vector<Position*>& positions = sys.getPositions();

	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd);
	if(!cRead.read()) {
		fout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();
	AtomPointerVector& glyAPV = cRead.getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	// Objects used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(opt.seed);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	hb.buildInteractions(30);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", true);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);

	/******************************************************************************
	*             Delete hydrogen bonds in the 3 (+3) terminal residues
	* 	      Residue numbers 1,2,3 and last,last-1,last-2
	*******************************************************************************/
	if(opt.deleteTerminalBonds) {
		fout << "\nDeleting the following bonds from the N- and C-terminal extension residues:\n" << endl;
		for (uint i=0; i< opt.deleteTerminalInteractions.size(); i++){
			fout << opt.deleteTerminalInteractions[i] << endl;
		}
		deleteTerminalBondInteractions(sys,opt);
	}

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/

	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";

	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	helicalAxis.saveCoor("originState");

	// Declare SelfPairManager and Set Seed
	SelfPairManager spm;
	spm.seed(RNG1.getSeed());

	// Set weights
	Eset->setWeight("CHARMM_ELEC", opt.weight_elec);
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);
	fout << "ELEC weight: " << Eset->getWeight("CHARMM_ELEC") << " VDW weight: " << Eset->getWeight("CHARMM_VDW") << " HB weight: " << Eset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;
	fout << endl;

	// TEMPORARY SET VARIABLE
		vector<string> variablePositionIdList;
		for (uint k=0; k < sys.positionSize(); k++) {
			Position &pos = sys.getPosition(k);

			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				variablePositionIdList.push_back(pos.getPositionId());
			}
		}
		sys.setVariablePositions(variablePositionIdList);

	spm.setSystem(&sys);
	spm.updateWeights();
	spm.setVerbose(opt.verbose);

	/******************************************************************************
	 *              === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	map<string,double> monomerEnergyByTerm;
	double monomerEnergy;

	if (opt.inputMonomerE) {
		monomerEnergy = opt.monoE_vdw + opt.monoE_hbond + opt.monoE_solv + opt.monoE_solvRef;
		monomerEnergyByTerm["CHARMM_VDW"] = opt.monoE_vdw;
		monomerEnergyByTerm["CHARMM_IMM1"] = opt.monoE_solv;
		monomerEnergyByTerm["CHARMM_IMM1REF"] = opt.monoE_solvRef;
		monomerEnergyByTerm["SCWRL4_HBOND"] = opt.monoE_hbond;
	} else {
		monomerEnergy = computeMonomerEnergy(sys, trans, opt, helicalAxis, RNG1, monomerEnergyByTerm, fout, opt.greedyCycles, opt.MCCycles, opt.MCMaxRejects);
	}

	fout << setiosflags(ios::fixed) << setprecision(4) << "Monomer Energy: " << monomerEnergy << endl;
	fout << endl;

	/******************************************************************************
	 *                  === READ IN RULES ===
	 ******************************************************************************/
	map<int, string> rulesFileMap;
	readRulesFile(opt.rulesFile, rulesFileMap);

	/******************************************************************************
	 *                  === NECESSARY VARIABLES ===
	 ******************************************************************************/
	// To store the structures in case we are clustering
	vector<HelixDimer*> structures;

	// We will need a format converter to convert to PDB names
	FormatConverter fc;

	map<string,double> energyByTerm;

	// Only calculate energies if necessary
	bool needToCalcEnergiesForThread = true;
	vector<vector<vector<vector<bool> > > > savedEnergyFlagTable;

	// Store info of top MC Repack structures
	vector<vector<double> > topMCRepackVector;

	// Global lowest energy found (if above monomer we won't save anyways)
	double globalLowestE = monomerEnergy;

	/******************************************************************************
	 *              === LOOP OVER ALL POSSIBLE INTERFACE THREADS ===
	 ******************************************************************************/
	for(int j = opt.threadStart; j <= opt.threadEnd; j++) {

		chainA.renumberChain(j);
		chainB.renumberChain(j);

		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();

		/******************************************************************************
		 *                     === LOAD ROTAMERS ===
		 ******************************************************************************/
		for (uint k=0; k < sys.positionSize(); k++) {
			Position &pos = sys.getPosition(k);

			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!sysRot.loadRotamers(&pos, pos.getResidueName(),"SL95.00")) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}

		sys.saveAltCoor("currentTMFragment");

		needToCalcEnergiesForThread = true;

		/******************************************************************************
		 *              === LOOP OVER ALL LINES IN HELIX GEOMETRY FILE ===
		 ******************************************************************************/
		for (uint i=0; i<fileVec.size()-1; i++) {

			chainA.renumberChain(j);
			chainB.renumberChain(j);

			helicalAxis.applySavedCoor("originState");
			sys.applySavedCoor("currentTMFragment");

			// Parse parameter file line
			fileVec[i] = MslTools::trim(fileVec[i], "\t\n\r");
			vector<string> parsedGeoInformation = MslTools::tokenize(fileVec[i], " ");
			int centroidIDNum = MslTools::toInt(parsedGeoInformation[0]);

			fout << centroidIDNum << " A:" << j << " B:" << j << endl;

			if (parsedGeoInformation.size() % 5 != 0) {
				cerr << "File line: " << i << " is of incompatible length" << endl;
				exit(1);
			}

			// index axialRot crossingAngle zShift xShift ChainDonor DonorResNum AcceptorResNum DonorAtom DistBondBreaks
			// 00001 75 -35.0 1.0 6.8 A 7 8 HA2 9.6
			// vector length = 5, 0 hbonds
			// vector length = 10, 1 hbonds
			// vector length = 15, 2 hbonds

			/******************************************************************************
			 *                     === HYDROGEN BOND COUNT CHECK ===
			 ******************************************************************************/
			double xShiftStart = 0;
			if(!hydrogenBondCheck(sys, opt, parsedGeoInformation, xShiftStart)) {
				fout << "less than " << opt.hbondCheckNumber << " hbonds" << endl;

				// renumber residues - rethreading
				//reThreadResidues(positions,-1);
				fout << endl;
				continue;
			}

			/******************************************************************************
			 *                     === CHECK SEQUENCE RULES ===
			 ******************************************************************************/
			if(!rulesCheck(sys, parsedGeoInformation[0], rulesFileMap)) {
				fout << "sequence does not conform with given rules" << endl;

				// renumber residues - rethreading
				//reThreadResidues(positions,-1);
				fout << endl;
				continue;
			}

			/******************************************************************************
			 *                  === CALCULATE ENERGIES TO SAVE ===
			 ******************************************************************************/
			//if(needToCalcEnergiesForThread) {
			//	spm.setOnTheFly(0);
			//	spm.calculateEnergies();
			//	savedEnergyFlagTable = createSavedEnergyFlagTable(sys);
			//	needToCalcEnergiesForThread = false;
			//}

			/******************************************************************************
			 *                     === INITIAL STARTING POSITION ===
			 ******************************************************************************/
			double xShift = xShiftStart;
			double crossingAngle = MslTools::toDouble(parsedGeoInformation[2]);
			double axialRotation = MslTools::toDouble(parsedGeoInformation[1]);
			double zShift = MslTools::toDouble(parsedGeoInformation[3]);

			fout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
			fout << " xShift: " << xShift << endl;

			// Transform helices to initial starting position
			transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
			moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

			// Optimizatize Initial Starting Position
			repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);

			sys.setActiveRotamers(spm.getMinStates()[0]);
			double currentEnergy = spm.getMinBound()[0];
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");

			/******************************************************************************
			 *              === CHECK AGAINST MONOMER ENERGY ===
			 ******************************************************************************/
			// if the outermost energy is still > opt.energyCutOff (default 100) it is likely we will not get <0 energy with this model

			fout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
			//sys.writePdb("first_" + MslTools::doubleToString(xShift) + ".pdb");
			if(currentEnergy-monomerEnergy > opt.energyCutOff) {
				fout << "discard model since deltaE Exceeds " << opt.energyCutOff << " at dOut" << endl;
				fout << endl;
				continue;
			}

			/******************************************************************************
			 *                     === X SHIFT REPACKS ===
			 ******************************************************************************/
			double bestEnergy = currentEnergy;
			double savedXShift = xShift;
			double previousEnergy = monomerEnergy;
			double deltaXShift = -0.1;
			double xShiftEnd = MslTools::toDouble(parsedGeoInformation[4]);

			while (xShift >= xShiftEnd) {

				xShift += deltaXShift;

				// Move the helix
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );

				// Run Optimization
				repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);

				vector<unsigned int> MCOFinal;
				MCOFinal = spm.getMinStates()[0];
				sys.setActiveRotamers(MCOFinal);

				currentEnergy = spm.getMinBound()[0];

				if (currentEnergy < bestEnergy) {
					bestEnergy = currentEnergy;
					savedXShift = xShift;
					sys.saveAltCoor("savedBestState");
					helicalAxis.saveAltCoor("BestAxis");
				}

				fout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

				// If energy increase twice in a row, and it is above the global lowest energy, quit
				if (currentEnergy < globalLowestE) {
					globalLowestE = currentEnergy;
				}
				if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
					fout << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
					break;
				}
				else {
					previousEnergy = currentEnergy;
				}

			}
			fout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl;

			vector<double> currentStuctureInfo;
			currentStuctureInfo.push_back(bestEnergy); // [0]
			currentStuctureInfo.push_back(j); // thread [1]
			currentStuctureInfo.push_back(savedXShift); // [2]
			currentStuctureInfo.push_back(crossingAngle); // [3]
			currentStuctureInfo.push_back(axialRotation); // [4]
			currentStuctureInfo.push_back(zShift); // [5]
			currentStuctureInfo.push_back(centroidIDNum); // [6]

			addStructureToTopMCRepackList(topMCRepackVector, currentStuctureInfo, fout, opt.numberOfStructuresToMCRepack);
			fout << endl;

		} // Geo Loop End
	} // Threading Loop End


	/******************************************************************************
	 *      === LOCAL BACKBONE MONTE CARLO REPACKS FOR TOP STRUCTURES ===
	 ******************************************************************************/
	fout << "==================================================================" << endl;
	fout << "Performing Monte Carlo Repacks on Top " << topMCRepackVector.size() << " Structures." << endl;
	fout << "==================================================================" << endl;
	fout << endl;
	for (uint j=0; j<topMCRepackVector.size(); j++) {

		chainA.renumberChain(topMCRepackVector[j][1]);
		chainB.renumberChain(topMCRepackVector[j][1]);

		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();

		/******************************************************************************
		 *                     === LOAD ROTAMERS ===
		 ******************************************************************************/
		for (uint k=0; k < sys.positionSize(); k++) {
			Position &pos = sys.getPosition(k);

			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!sysRot.loadRotamers(&pos, pos.getResidueName(),"SL95.00")) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}

		/******************************************************************************
		 *                     === INITIAL STARTING POSITION ===
		 ******************************************************************************/
		int thread = (int) topMCRepackVector[j][1];
		double xShift = topMCRepackVector[j][2];
		double crossingAngle = topMCRepackVector[j][3];
		double axialRotation = topMCRepackVector[j][4];
		double zShift = topMCRepackVector[j][5];
		int centroidIDNum = (int) topMCRepackVector[j][6];

		fout << "Monte Carlo Repack " << j+1 << " of " << topMCRepackVector.size() << endl;
		fout << setiosflags(ios::fixed) << setprecision(0) << centroidIDNum << " A: " << thread << " B: " << thread << endl;
		fout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " xShift: " << xShift << endl;


		// Transform helices to initial starting position
		transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
		moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

		// Optimizatize Initial Starting Position
		repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);

		sys.setActiveRotamers(spm.getMinStates()[0]);
		double currentEnergy = spm.getMinBound()[0];
		sys.saveAltCoor("savedBestState");
		helicalAxis.saveAltCoor("BestAxis");

		/******************************************************************************
		 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
		 ******************************************************************************/
		double bestEnergy = currentEnergy;
		time_t startTimeMC, endTimeMC;
		double diffTimeMC;
		time(&startTimeMC);

		if (opt.MCCycles > 0) {
			//MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);
			MonteCarloManager MCMngr(opt.MCStartTemp, opt.MCEndTemp, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects);

			MCMngr.setEner(bestEnergy);

			while(!MCMngr.getComplete()) {

				sys.applySavedCoor("savedBestState");
				helicalAxis.applySavedCoor("BestAxis");

				int moveToPreform = RNG1.getRandomInt(3);

				double deltaXShift = 0.0;
				double deltaZShift = 0.0;
				double deltaCrossingAngle = 0.0;
				double deltaAxialRotation = 0.0;

				//======================================
				//====== Z Shift (Crossing Point) ======
				//======================================
				if (moveToPreform == 0) {
					//deltaZShift = getStandardNormal(RNG1) * 0.1;
					deltaZShift = getStandardNormal(RNG1) * opt.deltaZ;
					backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
				} else if (moveToPreform == 1) {
				//===========================
				//===== Axial Rotation ======
				//===========================
					//deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
					deltaAxialRotation = getStandardNormal(RNG1) * opt.deltaAx;
					backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
				} else if (moveToPreform == 2) {
				//==================================
				//====== Local Crossing Angle ======
				//==================================
					//deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
					deltaCrossingAngle = getStandardNormal(RNG1) * opt.deltaCross;
					backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
				} else if (moveToPreform == 3) {
				//==============================================
				//====== X shift (Interhelical Distance) =======
				//==============================================
					//deltaXShift = getStandardNormal(RNG1) * 0.1;
					deltaXShift = getStandardNormal(RNG1) * opt.deltaX;
					backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
				}

				// Run Optimization
				repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);

				vector<unsigned int> MCOFinal = spm.getMinStates()[0];
				sys.setActiveRotamers(MCOFinal);
				currentEnergy = spm.getMinBound()[0];

				if (!MCMngr.accept(currentEnergy)) {
					//fout << "state rejected   energy: " << currentEnergy << endl;
				}
				else {
					// check number of hydrogen bonds and if it is less than hydrogen bond option dont accept
					vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);
					if(interHelicalHbonds.size() < opt.hbondCheckNumber) {
						// dont accept if we lost hydrogen bonds
						fout << "state rejected, lost Hbonds" << endl;
					} else {
						bestEnergy = currentEnergy;
						sys.saveAltCoor("savedBestState");
						helicalAxis.saveAltCoor("BestAxis");

						xShift = xShift + deltaXShift;
						crossingAngle = crossingAngle + deltaCrossingAngle;
						axialRotation = axialRotation + deltaAxialRotation;
						zShift = zShift +  deltaZShift;

						fout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
					}
				}
			}
		}
		time(&endTimeMC);
		diffTimeMC = difftime (endTimeMC, startTimeMC);
		fout << "Monte Carlo repack complete. Time: " << diffTimeMC << endl << endl;

		sys.applySavedCoor("savedBestState");

		double finalEnergy = sys.calcEnergy()-monomerEnergy;
		if(opt.printTermEnergies) {
			energyByTerm = getEnergyByTerm(sys.getEnergySet());

			for(map<string,double>::iterator it = monomerEnergyByTerm.begin(); it != monomerEnergyByTerm.end(); it++) {
				if(energyByTerm.find(it->first) != energyByTerm.end()) {
					energyByTerm[it->first] -= monomerEnergyByTerm[it->first];
				} else {
					// impossible
					fout << "ERROR 1242 Missing Energy term " << it->first << " Exiting" << endl;
					exit(0);
				}
			}
		}
		if (finalEnergy > 0) { // TODO Check number of Hbonds???
			// renumber residues - rethreading
			//fout << "Energy " << finalEnergy << " above zero - not written" << endl;
			//fout << endl;
			////reThreadResidues(positions,-1);
//
			if (opt.onlySaveNegativeStructures){
				fout << "Energy " << finalEnergy << " above zero - not written" << endl;
				fout << endl;
				continue;
			} else {
				fout << "Energy " << finalEnergy << " above zero - written" << endl;
				fout << endl;
			}
		}
//TODO: I think if all of the below was a function that took a system, it would be possible to just take the best positive structure

		// Renumber structures (all should go 1 - N, or from starting point determined by user)
		//int currStartNum = positions[0]->getResidueNumber();
		chainA.renumberChain(opt.startResNum);
		chainB.renumberChain(opt.startResNum);
		//renumberResidues(sys,opt.startResNum);

		// get vector of Hbonds, check size
		vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);

		// Store everything
		// create class structure to store data
		HelixDimer * st = new HelixDimer("",finalEnergy,thread);

		st->addAtoms(sys.getAtomPointers());
		// save the structure interface residues, energy, hbond list, CA of closest approach mark it with 2 in the mask
		map<string, unsigned int> interfaceMap = interfaceResidueCheck(apvChainA, apvChainB);

		//unsigned int closestCA = CAOfClosestApproach(chainA,chainB);
		string interfaceString = "";
		for (uint r=0; r < chainA.positionSize(); r++) {
			sequence = sequence + MslTools::getOneLetterCode(chainA.getPosition(r).getResidueName());

			if (interfaceMap.find(chainA.getIdentity(r).getIdentityId()) != interfaceMap.end()) {
				interfaceString = interfaceString + "1";
			}
			else {
				interfaceString = interfaceString + "0";
			}
		}

		//mark residue with CA of closest approach with 2
		//interfaceString.replace(closestCA - opt.startResNum,1,"2");

		st->setHelixDimerDetails(xShift,zShift,axialRotation,crossingAngle,interfaceString,prolineMask,interHelicalHbonds,centroidIDNum);
		if(opt.printTermEnergies) {
			st->setDeltaEnergyByTerm(energyByTerm);
		}

		if (opt.clusterSolutions) {
			// Save structure to vector of structures
			fout << "clustering structure with Energy: " << finalEnergy << endl;
			fout << endl;
			structures.push_back(st);

			AtomPointerVector& axisAtoms = helicalAxis.getAtomPointers();
			for(AtomPointerVector::iterator it =  axisAtoms.begin(); it != axisAtoms.end(); it++) {
				Atom* a = new Atom(* *it);
				st->addAxisAtom(a);
			}

		}
		else {
			char tmp[1000];


			sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.crd", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
			CRDWriter crd;
			string crdName = opt.pdbOutputDir+"/"+string(tmp);
			crd.open(crdName);
			crd.write(sys.getAtomPointers());
			fout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
			crd.close();

			// TODO convert to PDB names before printing
			sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.pdb", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
			fout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
			sys.writePdb(opt.pdbOutputDir+"/"+string(tmp));

			sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.txt", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
			fout << "wrote text file: " << string(tmp) << endl;

			ofstream txt;
			txt.open(tmp);
			if(txt.is_open()) {
				txt << "resStart " << opt.startResNum << endl;
				txt << "resEnd " << opt.endResNum << endl;
				txt << "TMStart " << opt.startResNum << endl;
				txt << "TMEnd " << opt.endResNum << endl;
				txt << "originalSeq " << originalTMSeq << endl;
				txt << "modelledSeq " << modelledTMSeq << endl;
				st->printHelixDimerDetails(txt);
				delete st;

			} else {
				cerr << "Unable to open " << tmp << endl;
			}


			fout << sys.getEnergySummary();
			fout << endl;
		}

	}


	/******************************************************************************
	 *                   === CLUSTER CATM MODEL SOLUTIONS ===
	 ******************************************************************************/
	if (opt.clusterSolutions) {
		sort(structures.begin(),structures.end(),compareHelixDimers());
		// convert all structures to PDB names
		/*
		for(int i = 0; i < structures.size(); i++) {
			How do we handle changed hydrogen bond names?
			AtomPointerVector & a = structures[i]->getAtomPointers();
			for(int j = 0; j < a.size(); j++) {
				int resNum = a[j]->getResidueNumber();
				fc.setPdbFromCharmm(*a[j],"22",resNum == opt.startResNum,resNum == opt.endResNum);
			}
		}
		*/
		vector<HelixDimerCluster*> clusters;
		clusterSolutions(clusters,structures,opt.rmsdCutoff,originalTMSeq,modelledTMSeq,opt);
		//cout << "Output " << clusters.size() << " HelixDimerClusters at rmsd " << _rmsdCutoff <<  endl;
		for(int i = 1; i <= clusters.size(); i++) {
			clusters[i-1]->printHelixDimerClusterCrds(i,opt.printAllCrds,opt.printAxes);
			clusters[i-1]->printDetails(i);
			if(opt.printTermEnergies) {
				clusters[i-1]->printTermEnergies(i);
			}
		}


		for(int i = 1; i <= clusters.size(); i++) {
			//clusters[i-1]->convertToPdbNames();
			clusters[i-1]->printHelixDimerClusterPdbs(i);
			clusters[i-1]->makePse(i);
		}

		// TODO delete HelixDimerCluster* form the vector clusters

	}

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	fout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
}


// Helper functions
void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

