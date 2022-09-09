#include <sstream>
#include <iterator>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
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
#include "ResidueSelection.h"
#include "SasaCalculator.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "localBackboneMove";
string programDescription = "This is an updated version of geomRepack: geomRepack does local backbone moves, then mutates sequences, but doesn't\n\
 	make moves for the mutated sequences. This will take a given backbone and do the local repacks for it. I've copied the code from two programs:\n\
 	readPDBAndCalcEnergy.cpp and geomRepack.cpp and just readjusted it for this program to work properly.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "27 May 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

//TODO: I think I have an idea for how to run this. I can queue everything in a file that's trimmed from the original csv file with clashing information,
// and queue it up as variables separated by commas. I just need to find the variables I need, get rid of all of the ones I don't, and then fix any constant
// variables and hard code it here. So the only thing I need to do is copy stuff over, make an option.cpp and options.h, functions.cpp and functions.h, and 
// test until it works. Should only take a couple of days.

// All of the below is copied from readPDBAndCalcEnergy
int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);

	string date(buffer);
	
    /******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
	//Add in some default options that can easily be changed here
	Options opt = parseOptions(argc, argv, defaults);

	if (opt.errorFlag) {
		outputErrorMessage(opt);
		exit(1);
	}

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	// summary file output
	ofstream sout;
	// monomer file output
	ofstream mout;
	// error file output
	ofstream err;
	// rerun config output
	ofstream rerun;

	setupOutputDirectory(opt);

	string soutfile = opt.outputDir + "/energy.csv";
	string moutfile = opt.outputDir + "/monomer.out";
	string errfile  = opt.outputDir + "/errors.out";
	string rerunfile = opt.outputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	rerun << opt.rerunConf << endl;
	rerun.close();

	err << date << endl;

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
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	System pdb;
	pdb.readPdb(opt.infile);

	// Set up chain A and chain B atom pointer vectors
	Chain & chain1A = pdb.getChain("A");
	Chain & chain1B = pdb.getChain("B");
	AtomPointerVector & apvChain1A = chain1A.getAtomPointers();
	AtomPointerVector & apvChain1B = chain1B.getAtomPointers();

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Transform to chosen geometry
	// Change this function to give other geometric parameters for the independent helices for the heterodimers
	// This is directly from CATM
	transformation(apvChain1A, apvChain1B, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
    System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	CSB.setBuildNonBondedInteractions(false);
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
    /******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// assign the coordinates of our system to the given geometry that was assigned without energies using System pdb
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	// initialize the object for loading rotamers into our system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	int firstPos = 0;
    int lastPos = sys.positionSize();
    deleteTerminalHydrogenBondInteractions(sys,firstPos,lastPos);

	// Up to here is from readDPBAndCalcEnergy
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// Assign number of rotamers by residue burial
	loadRotamersBySASABurial(sys, sysRot, opt);
	CSB.updateNonBonded(10,12,50);

	//slightly convoluted, but something with the helical axis just wouldn't get correct if I didn't tranform it, but if I did it messed with my structure...so I worked around above. I don't think it should affect my results in anyway other than putting my helices in the membrane, and it's consistent with what I did for the results in my design code

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 

	// Optimize Initial Starting Position 
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	
	//Switch to current sequence
	//TODO: get sequence from the pdb here
	string sequence = opt.sequence;
	double bestEnergy = sys.calcEnergy();
	
	//Repack dimer
	repackSideChains(spm, 10);
	vector<uint> startStateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(startStateVec);
	double currentEnergy = spm.getStateEnergy(startStateVec);
	sys.setActiveRotamers(startStateVec);

	
	/******************************************************************************
	 *                    === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
    map<string,double> monomerEnergyByTerm;
    double monomer = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);
    
	double startDimer = sys.calcEnergy();
    cout << "Monomer Energy: " << monomer << endl;
    cout << "Dimer-Monomer: " << startDimer-monomer << endl;
    
	double monomerEnergy = monomer;

	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	cout << "Starting Geometry" << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << opt.xShift << " crossingAngle: " << opt.crossingAngle << " axialRotation: " << opt.axialRotation << " zShift: " << opt.zShift << endl << endl;
	cout << "Current Best Energy: " << bestEnergy-monomerEnergy << endl;
	cout << "Interaction Energies: " << endl;
	cout << spm.getSummary(startStateVec) << endl;
	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
		
	double bestRepackEnergy;
	vector<uint> bestRepackState;

	double bestXShift = opt.xShift;
	double bestAxialRotation = opt.axialRotation;
	double bestZShift = opt.zShift;
	double bestCrossingAngle = opt.crossingAngle;
	
	double prevBestEnergy = startDimer;
	sys.applySavedCoor("startingState");
	helicalAxis.applySavedCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
		
	double xShift = opt.xShift;
	double crossingAngle = opt.crossingAngle;
	double axialRotation = opt.axialRotation;
	double zShift = opt.zShift;
	
	if (opt.verbose){
		cout << "======================================" << endl;
		cout << "Performing Local Monte Carlo Repack   " << endl;
		cout << "======================================" << endl;
	}
	vector<unsigned int> MCOBest = startStateVec;
		
	//MonteCarloManager MCMngr(opt.MCStartTemp, opt.MCEndTemp, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects);
	MonteCarloManager MCMngr(100, 0.5, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects);
	// MonteCarloManager MCMngr(1000, 0.5, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects, 10, 0.01);
	//MonteCarloManager MCMngr(1000, 0.5, 1000, opt.MCCurve, opt.MCMaxRejects);
	//MonteCarloManager MCMngr(1000, 0.5, 10, 2, 2);//same amount as in monomer
		
	MCMngr.setEner(prevBestEnergy);
		
	unsigned int counter = 0;
	while(!MCMngr.getComplete()) {
		
		sys.applySavedCoor("savedBestState");
		helicalAxis.applySavedCoor("BestAxis");
		
		int moveToPreform = RNG.getRandomInt(3);
		
		double deltaXShift = 0.0;
		double deltaZShift = 0.0;
		double deltaCrossingAngle = 0.0;
		double deltaAxialRotation = 0.0; 
		
		//======================================
		//====== Z Shift (Crossing Point) ======
		//======================================
		if (moveToPreform == 0) {
			//deltaZShift = getStandardNormal(RNG1) * 0.1;
			deltaZShift = getStandardNormal(RNG) * opt.deltaZ;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(RNG) * opt.deltaAx;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(RNG) * opt.deltaCross;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(RNG1) * 0.1;
			deltaXShift = getStandardNormal(RNG) * opt.deltaX;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
		}
		
		// Run Optimization
		repackSideChains(spm, 10);

		vector<unsigned int> MCOFinal = spm.getMinStates()[0];
		currentEnergy = spm.getMinBound()[0];
		
		if (!MCMngr.accept(currentEnergy)) {
			if (opt.verbose){
				cout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
		} else {
			prevBestEnergy = currentEnergy;
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");
		
			xShift = xShift + deltaXShift;
			crossingAngle = crossingAngle + deltaCrossingAngle;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift + deltaZShift;
			MCOBest = MCOFinal;
		
			if (opt.verbose){
				cout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
			counter++;
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
		
	sys.applySavedCoor("savedBestState");
	double dimerEnergy = spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-monomerEnergy;
	double vdw = spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	double hbond = spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	double imm1 = spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	double dimerDiff = dimerEnergy-startDimer;
	double xShiftDiff = opt.xShift-xShift;
	double angleDiff = opt.crossingAngle-crossingAngle;
	double axialRotDiff = opt.axialRotation-axialRotation;
	double zShiftDiff = opt.zShift-zShift;

	// Print out info to the summary csv file
	//vector<string> cols{'Sequence','Total','Dimer','Monomer','DimerDiff','VDWDimer','VDWDiff','HBONDDiff','IMM1Dimer','IMM1Diff','startXShift','xShift','xShiftDiff','startCrossingAngle','crossingAngleDiff','startAxialRotation','startZShift','zShift','zShiftDiff'};
	//for (uint i = 0; i < cols.size(); i++){
	//	if (i < cols.size()-1){
	//		sout << cols[i] << ',';
	//	} else {
	//		sout << cols[i] << endl;
	//	}
	//}
	//sout << 'Sequence,Total,Dimer,Monomer,DimerDiff,VDWDimer,VDWDiff,HBONDDiff,IMM1Dimer,IMM1Diff,startXShift,xShift,xShiftDiff,startCrossingAngle,crossingAngleDiff,startAxialRotation,startZShift,zShift,zShiftDiff' << endl;
	sout << sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << monomer << ',' << dimerDiff << ',' << vdw << ',' << hbond << ',' << imm1 << ',' << opt.xShift << ',' << xShift << ',' << xShiftDiff << ',' << opt.crossingAngle << ',' << crossingAngle << ',' << angleDiff << ',' << opt.axialRotation << ',' << axialRotation << ',' << axialRotDiff << ',' << opt.zShift << ',' << zShift << ',' << zShiftDiff << endl;
	cout << sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << monomer << ',' << dimerDiff << ',' << vdw << ',' << hbond << ',' << imm1 << ',' << opt.xShift << ',' << xShift << ',' << xShiftDiff << ',' << opt.crossingAngle << ',' << crossingAngle << ',' << angleDiff << ',' << opt.axialRotation << ',' << axialRotation << ',' << axialRotDiff << ',' << opt.zShift << ',' << zShift << ',' << zShiftDiff << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

    //outputs a pdb file for the structure 
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(opt.outputDir + "/" + opt.sequence + "_backboneOptimized.pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

    // close all of the output file writers
    sout.close();
    mout.close();
    err.close();
}