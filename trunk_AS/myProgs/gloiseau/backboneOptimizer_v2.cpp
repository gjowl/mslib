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
#include "backboneOptimizerFunctions_v2.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "backboneOptimizer";
string programDescription = "This is an updated version of backboneOptimizer: backboneOptimizer does local backbone moves on given backbones\n\
 	from an input pdb file and an input configuration file that gives the positions that are interfacial. This version takes code from seqDesign\n\
 	to identify interfacial amino acids and place the backbone into a given geometry rather than needing inputs (helpful for running on chtc.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "7 July 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

// Options Functions
void prepareSystem(Options &_opt, System &_sys, string _polySeq, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, CartesianPoint &_xAxis,
 CartesianPoint &_zAxis, Transforms &_trans);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
 CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels,
 string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, 
 vector<uint> &_allInterfacePositions, vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq);
vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);

int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	// Start the timer
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
	Options opt = backboneOptimizerParseOptions(argc, argv, defaults);

	if (opt.errorFlag) {
		backboneOptimizerOutputErrorMessage(opt);
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

	// function that defines the output directory
	setupOutputDirectoryChtc(opt);

	// setup the output files
	string soutfile = opt.outputDir + "/energy.csv";
	string moutfile = opt.outputDir + "/monomer.out";
	string errfile  = opt.outputDir + "/errors.out";
	string rerunfile = opt.outputDir + "/rerun.config";

	// open the output files
	sout.open(soutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	// write the rerun config file
	rerun << opt.rerunConf << endl;
	// close the rerun config file
	rerun.close();

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
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	/******************************************************************************
	 *                       === GENERATE POLYMER SEQUENCE ===
	 ******************************************************************************/
	// polymer sequences have: chain, starting position of chain residue, three letter AA code
	string polySeq = generatePolymerSequence("L", opt.backboneLength, opt.thread);
	PolymerSequence PS(polySeq);
	
	// set up the system
	System sys;
	prepareSystem(opt, sys, polySeq, helicalAxis, axisA, axisB, ori, xAxis, zAxis, trans);

	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	/******************************************************************************
	 *                          === PRINT GEOMETRY ===
	 ******************************************************************************/
	// Output the starting geometry
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << opt.xShift << endl;
	cout << "crossingAngle: " << opt.crossingAngle << endl;
	cout << "axialRotation: " << opt.axialRotation << endl;
	cout << "zShift:        " << opt.zShift << endl;

	//String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.ids);
	cout << "Amino acids for design: LEU " << alternateIds << endl;
	
	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// Variables to output from defineInterfaceAndRotamerSampling function
	string rotamerLevels;
	string variablePositionString;
	string rotamerSamplingString;
	// vector of the positions that will be linked
	vector<int> linkedPositions;
	// vector of positions at the interface excluding termini positions
	vector<uint> interfacePositions;
	// vector of positions at the interface including the terminal positions
	vector<uint> allInterfacePositions;
	// vector of rotamer level for each position
	vector<int> rotamerSamplingPerPosition;

	// Defines the interfacial positions and the number of rotamers to give each position
	// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
	// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
	// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
	defineInterfaceAndRotamerSampling(opt, PS, rotamerLevels, polySeq, variablePositionString, rotamerSamplingString, linkedPositions, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, sout, axis);
	
	
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	// setup random number generator object
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
	
	// Switch to current sequence
	string sequence = opt.sequence;
	
	// Repack dimer
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
	
	// calculate the energy of the system
	double startDimer = sys.calcEnergy();
    cout << "Monomer Energy: " << monomer << endl;
    cout << "Dimer-Monomer: " << startDimer-monomer << endl;
    
	double monomerEnergy = monomer;

	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	/******************************************************************************
	 *                     === X SHIFT REPACKS ===
	 ******************************************************************************/
	double xTranslate = 7;
	double bestEnergy = currentEnergy-monomerEnergy;
	double xShift = opt.xShift+xTranslate/2;
	double savedXShift = xShift;
	double previousEnergy = monomerEnergy;
	double deltaXShift = -0.1;
	double xShiftEnd = 6.5;

	// Global lowest energy found (if above monomer we won't save anyways)
	double globalLowestE = monomerEnergy;

	// while loop for x-shifts	
	while (xShift >= xShiftEnd) {

		// add the xShift change to the current xShift
		xShift += deltaXShift;

		// Move the helix
		backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );

		// Run Optimization
		repackSideChains(spm, opt.greedyCycles);

		vector<unsigned int> MCOFinal;
		MCOFinal = spm.getMinStates()[0];
		sys.setActiveRotamers(MCOFinal);

		currentEnergy = spm.getMinBound()[0];

		// Check if this is the lowest energy found so far
		if (currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			savedXShift = xShift;
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");
		}
		cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

		// If energy increase twice in a row, and it is above the global lowest energy, quit
		if (currentEnergy < globalLowestE) {
			globalLowestE = currentEnergy;
		}
		if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
			cout << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
			break;
		} else {
			previousEnergy = currentEnergy;
		}
	}
	cout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl;
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	PDBWriter writer1;
	writer1.open(opt.outputDir + "/" + opt.sequence + "_xShifted.pdb");
	writer1.write(sys.getAtomPointers(), true, false, true);
	writer1.close();

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	cout << "Starting Geometry" << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << xShift << " crossingAngle: " << opt.crossingAngle << " axialRotation: " << opt.axialRotation << " zShift: " << opt.zShift << endl << endl;
	cout << "Current Best Energy: " << bestEnergy-monomerEnergy << endl;
	cout << "Interaction Energies: " << endl;
	cout << spm.getSummary(startStateVec) << endl;

	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
		
	double bestRepackEnergy;
	vector<uint> bestRepackState;

	double bestXShift = xShift;
	double bestAxialRotation = opt.axialRotation;
	double bestZShift = opt.zShift;
	double bestCrossingAngle = opt.crossingAngle;
	
	double prevBestEnergy = startDimer;
	sys.applySavedCoor("savedBestState");
	helicalAxis.applySavedCoor("BestAxis");
		
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
		
	MCMngr.setEner(prevBestEnergy);
		
	writer1.open(opt.outputDir + "/" + opt.sequence + "_beforeRepack.pdb");
	writer1.write(sys.getAtomPointers(), true, false, true);
	writer1.close();
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
			writer1.open(opt.outputDir + "/" + opt.sequence + "_duringRepack.pdb");
			writer1.write(sys.getAtomPointers(), true, false, true);
			writer1.close();
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

	// calculate the solvent accessible surface area
	SasaCalculator sasa(sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	// Print out info to the summary csv file
	//sout << 'Sequence,Total,Dimer,Monomer,DimerDiff,VDWDimer,VDWDiff,HBONDDiff,IMM1Dimer,IMM1Diff,startXShift,xShift,xShiftDiff,startCrossingAngle,crossingAngleDiff,startAxialRotation,startZShift,zShift,zShiftDiff' << endl;
	sout << sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << monomer << ',' << dimerDiff << ',' << dimerSasa << ',' << vdw << ',' << hbond << ',' << imm1 << ',' << opt.xShift << ',' << xShift << ',' << xShiftDiff << ',' << opt.crossingAngle << ',' << crossingAngle << ',' << angleDiff << ',' << opt.axialRotation << ',' << axialRotation << ',' << axialRotDiff << ',' << opt.zShift << ',' << zShift << ',' << zShiftDiff << endl;
	cout << sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << monomer << ',' << dimerDiff << ',' << dimerSasa << ',' << vdw << ',' << hbond << ',' << imm1 << ',' << opt.xShift << ',' << xShift << ',' << xShiftDiff << ',' << opt.crossingAngle << ',' << crossingAngle << ',' << angleDiff << ',' << opt.axialRotation << ',' << axialRotation << ',' << axialRotDiff << ',' << opt.zShift << ',' << zShift << ',' << zShiftDiff << endl;
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

//Functions
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, CartesianPoint &_xAxis,
 CartesianPoint &_zAxis, Transforms &_trans) {
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	_sys.readPdb(_opt.backboneFile);

	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, _axisA, _axisB, _ori, _xAxis, _zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

void prepareSystem(Options &_opt, System &_sys, string _polySeq, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, CartesianPoint &_xAxis,
 CartesianPoint &_zAxis, Transforms &_trans){
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	System startGly69;
	setGly69ToStartingGeometry(_opt,startGly69,_helicalAxis,_axisA,_axisB,_ori,_xAxis,_zAxis,_trans);
	
	// declare system
	CharmmSystemBuilder CSB(_sys,_opt.topFile,_opt.parFile,_opt.solvFile);
	if (_opt.useElec == false){
		CSB.setBuildTerm("CHARMM_ELEC", false);
	} else {
		CSB.setBuildTerm("CHARMM_ELEC", true);
	}
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	//
	CSB.setIMM1Params(15, 10);
	//
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	PolymerSequence PL(_polySeq);
	if(!CSB.buildSystem(PL)) {
		cerr << "Unable to build system from " << _polySeq << endl;
		exit(0);
	}

	// get chain A and B from the _system
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// assign the coordinates of our system to the given geometry 
	_sys.assignCoordinates(startGly69.getAtomPointers(),false);
	_sys.buildAllAtoms();

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(_sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(_sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = _sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	//Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	int firstPos = 0;
    int lastPos = _sys.positionSize();
    deleteTerminalHydrogenBondInteractions(_sys,firstPos,lastPos);

	// Up to here is from readDPBAndCalcEnergy
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// Assign number of rotamers by residue burial
	loadRotamersBySASABurial(_sys, sysRot, _opt);
	CSB.updateNonBonded(10,12,50);
	
	// as of 2022-7-5: not sure if the above works or needs to be reworked
	PDBWriter writer1;
	writer1.open(_opt.outputDir + "/" + _opt.sequence + "_inputGeometry.pdb");
	writer1.write(_sys.getAtomPointers(), true, false, true);
	writer1.close();

}

void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq,
 string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis){
	// Declare system
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setBuildNonBondedInteractions(false);
	//CSB.setBuildNoTerms();

	if(!CSB.buildSystem(_PS)) {
		cout << "Unable to build system from " << _PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);

	CSB.updateNonBonded(10,12,50);

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(_opt.backboneFile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate

	Chain & chainA = pdb.getChain("A");
	Chain & chainB = pdb.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	PDBReader readAxis;
	if(!readAxis.read(_axis)) {
		cout << "Unable to read axis" << endl;
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

	// Objects used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Transformation to zShift, axialRotation, crossingAngle, and short xShift to identify potential interacting positions
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, trans);//one of the shortest distances given from my pdb searc
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	//TODO: add in a comparison for the monomer at each position here instead
	string backboneSeq = generateString("L", _opt.backboneLength);
	vector<pair <int, double> > resiBurial = calculateResidueBurial(sys, _opt, backboneSeq);
	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});

	//cout << "Determining interfacial residues by residue burial..." << endl;
	int levelCounter = 0;
	vector<int> interfacePositions;

	// Output variable Set up
	backboneSeq = generateBackboneSequence("L", _opt.backboneLength, _opt.useAlaAtCTerminus);
	string variablePositionString = generateString("0", _opt.backboneLength);
	string rotamerLevels = generateString("0", _opt.backboneLength);

	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	int highestRotamerLevel = numberOfRotamerLevels-1;
	//make into a function
	int numAAs = 0;
	double lvlavg = 0;
	vector<double> avgs;
	//cout << "lvl " << levelCounter << endl;
	for (uint i = 0; i < resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(resiBurial.size());
		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
			levelCounter++;
			lvlavg = lvlavg/numAAs;
			avgs.push_back(lvlavg);
			lvlavg=0;
			numAAs=0;
			//cout << "lvl " << levelCounter << endl;
		}
		//cout << resiBurial[i].first << ": " << resiBurial[i].second << endl;
		lvlavg = lvlavg+resiBurial[i].second;
		numAAs++;
		int backbonePosition = resiBurial[i].first;
		Position &pos = sys.getPosition(backbonePosition);
		string posRot = _opt.sasaRepackLevel[levelCounter];
		int resiNum = pos.getResidueNumber();
		int posNum = resiNum-_opt.thread;
		string add;

		if (levelCounter < _opt.interfaceLevel){
			add = "Add all Ids at this pos";
			interfacePositions.push_back(resiNum);
			if (backbonePosition > 2 && backbonePosition < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 3 and 4 here instead of 4 and 5 to prevent changes at the interface like others
				variablePositionString.replace(variablePositionString.begin()+posNum, variablePositionString.begin()+posNum+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		} else {
			add = "Only 1 ID";
		}
		rotamerLevels.replace(rotamerLevels.begin()+posNum, rotamerLevels.begin()+posNum+1, MslTools::intToString(levelCounter));
	}
	lvlavg = lvlavg/numAAs;
	avgs.push_back(lvlavg);
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.ids, interfacePositions);

	vector<int> rotamerSamplingPerPosition = getRotamerSampling(rotamerLevels);
	vector<int> linkedPositions = getLinkedPositions(rotamerSamplingPerPosition, _opt.interfaceLevel, highestRotamerLevel);

	//String for the positions of the sequences that are considered interface for positions amd high rotamers
	string rotamerSamplingString = getInterfaceString(rotamerSamplingPerPosition, _opt.backboneLength);

	// Vector for linked positions in "A,25 B,25" format


	// Define referenced output variables
	_rotamerLevels = rotamerLevels;
	_polySeq = polySeq;
	_rotamerSamplingPerPosition = rotamerSamplingPerPosition;
	_variablePositionString = variablePositionString;
	_rotamerSamplingString = rotamerSamplingString;
	_linkedPositions = linkedPositions;
	_out << endl;
	_out << "PolyLeu Backbone:   " << backboneSeq << endl;
	_out << "Variable Positions: " << variablePositionString << endl;
	_out << "Rotamers Levels:    " << rotamerSamplingString << endl;
	_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition);
	_allInterfacePositions = getAllInterfacePositions(_opt, rotamerSamplingPerPosition);

	//cout << endl << "Averages:\t";
	//for (uint a=0; a<avgs.size(); a++){
	//	cout << avgs[a] << "\t";
	//}
	//cout << endl;

	//for (uint i=0; i<_interfacePositions.size(); i++){
	//	cout << _interfacePositions[i] << ",";
	//}
	//cout << endl;
	cout << endl;
	cout << "PolyLeu Backbone:   " << backboneSeq << endl;
	cout << "Variable Positions: " << variablePositionString << endl;
	cout << "Rotamers Levels:    " << rotamerSamplingString << endl;

	int numPosAtInterface = 0;
	for(string::iterator it = variablePositionString.begin(); it != variablePositionString.end();it++ ) {
		stringstream ss;
		ss << *it;
		string num = ss.str();
		if (MslTools::toInt(num) == 1){
			numPosAtInterface++;
		}
	}
}

// calculate the SASA burial for each residue
std::vector<pair <int, double> > calculateResidueBurial (System &_sys) {
	/*
	  SASA reference:
	  Protein Engineering vol.15 no.8 pp.659–667, 2002
	  Quantifying the accessible surface area of protein residues in their local environment
	  Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
	*/
	map<string,double> refSasa;
	refSasa["G"] = 83.91;
	refSasa["A"] = 116.40;
	refSasa["S"] = 125.68;
	refSasa["C"] = 141.48;
	refSasa["P"] = 144.80;
	refSasa["T"] = 148.06;
	refSasa["D"] = 155.37;
	refSasa["V"] = 162.24;
	refSasa["N"] = 168.87;
	refSasa["E"] = 187.16;
	refSasa["Q"] = 189.17;
	refSasa["I"] = 189.95;
	refSasa["L"] = 197.99;
	refSasa["H"] = 198.51;
	refSasa["K"] = 207.49;
	refSasa["M"] = 210.55;
	refSasa["F"] = 223.29;
	refSasa["Y"] = 238.30;
	refSasa["R"] = 249.26;
	refSasa["W"] = 265.42;

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator b(_sys.getAtomPointers());
	//b.setProbeRadius(3.0);
	b.calcSasa();
	for (uint i = 0; i < _sys.positionSize()/2; i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdA = _sys.getPosition(i).getPositionId();
		//string posIdB = _sys.getPosition(i+_sys.positionSize()/2).getPositionId();
		string oneLetter = MslTools::getOneLetterCode(_sys.getPosition(i).getResidueName());
		double resiSasa = b.getResidueSasa(posIdA);
		double burial = resiSasa / refSasa[oneLetter];
		residueBurial.push_back(pair<int,double>(i, burial));
		//residueBurial.push_back(pair<string,double>(posIdB, burial));
		//cout << posId << "\t" << resiSasa << "\t" << burial << endl;
	}

	return residueBurial;
}

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq) {
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);
	PolymerSequence PS(polySeq);

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
	}

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	System gly69;
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	gly69.readPdb(_opt.backboneFile);

	AtomPointerVector& glyAPV = gly69.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator dimerSasa(_sys.getAtomPointers());
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	dimerSasa.calcSasa();
	monoSasa.calcSasa();
	dimerSasa.setTempFactorWithSasa(true);

	for (uint i = 0; i < monoSys.positionSize(); i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdMonomer = monoSys.getPosition(i).getPositionId();
		string posIdDimer = _sys.getPosition(i).getPositionId();
		double resiSasaMonomer = monoSasa.getResidueSasa(posIdMonomer);
		double resiSasaDimer = dimerSasa.getResidueSasa(posIdDimer);
		double burial = resiSasaDimer/resiSasaMonomer;
		residueBurial.push_back(pair<int,double>(i, burial));

		//set sasa for each residue in the b-factor
		AtomSelection selA(_sys.getPosition(i).getAtomPointers());
		AtomSelection selB(_sys.getPosition(i+_opt.backboneLength).getAtomPointers());
		AtomPointerVector atomsA = selA.select("all");
		AtomPointerVector atomsB = selB.select("all");
		for (AtomPointerVector::iterator k=atomsA.begin(); k!=atomsA.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
		for (AtomPointerVector::iterator k=atomsB.begin(); k!=atomsB.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
	}

	PDBWriter writer;
	writer.open(_opt.outputDir+"/interfaceSASA.pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	return residueBurial;
}

vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=0; k<_opt.backboneLength; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=3; k<_opt.backboneLength-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}
