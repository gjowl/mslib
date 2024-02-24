#include <sstream>
#include <iterator>
#include <unistd.h>

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
#include "backboneOptimizerOptions.h"
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

time_t startTime, endTime;
double diffTime;

// BBOptions Functions
void prepareSystem(BBOptions &_opt, System &_sys, System &_startGeom, string &_polySeq);
void setGly69ToStartingGeometry(BBOptions &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
 CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
vector<pair <int, double>> getResiBurial(System &_sys, BBOptions &_opt, string _sequence);
void defineRotamerLevels(System &_sys, BBOptions &_opt, vector<pair <int, double>> _resiBurial, 
 vector<uint> &_interfacePositions, string &_variablePositionString, string &_rotamerLevels);
void defineInterfaceAndRotamerSampling(BBOptions &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
 string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, 
 vector<uint> &_allInterfacePositions, vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, BBOptions &_opt, string _seq);
vector<uint> getAllInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition);
void localXShiftDocking(System &_sys, BBOptions &_opt, double &_bestEnergy, double _monomerEnergy, 
 SelfPairManager &_spm, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB,
 AtomPointerVector &_apvChainA, AtomPointerVector &_apvChainB, Transforms &_trans, double &_savedXShift);
void localBackboneRepack(BBOptions &_opt, System &_sys, double _savedXShift, SelfPairManager &_spm,
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA, AtomPointerVector &_apvChainB,
 Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm, double _monomerEnergy, ofstream &_out);
void monteCarloRepack(BBOptions &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 map<string,double> _monomerEnergyByTerm, double _monomerEnergy, ofstream &_out);
void checkOptionErrors(BBOptions &_opt);	

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k < _sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
}


/***********************************
 *help functions
 ***********************************/
void usage();
void help(BBOptions defaults);
void outputErrorMessage(BBOptions &_opt);
void outputWarningMessage(BBOptions &_opt);

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
	BBOptions defaults;
	//Add in some default options that can easily be changed here
	BBOptions opt = BBParseOptions(argc, argv, defaults);
	checkOptionErrors(opt);

	/******************************************************************************
	 *                       === SETUP OUTPUTS ===
	 ******************************************************************************/
	ofstream sout;  // summary file output
	ofstream mout;  // monomer file output
	ofstream err;   // error file output
	ofstream rerun; // rerun config output

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

	// write and close the rerun config file
	rerun << opt.rerunConf << endl;
	rerun.close();
	
	// get the starting geometries; convert to parallelogram axialRot and Z (Mueller 2014; Fig. S1)
	convertToAxAndZForTranformation(opt);
	
	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	helicalAxis.saveCoor("originState");
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	/******************************************************************************
	 *                       === SETUP STARTING GEOMETRY ===
	 ******************************************************************************/
	// get the polymer sequence for determining the interface
	string polySeq = generatePolymerSequence("V", opt.sequence.length(), opt.thread);
	PolymerSequence PS(polySeq);

	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,axisA,axisB,ori,xAxis,zAxis,trans);

	//helicalAxis.applySavedCoor("originState");
	// set up the system for the input sequence
	System sys;
	string polySeq1 = convertToPolymerSequenceNeutralPatch(opt.sequence, opt.thread);
	prepareSystem(opt, sys, startGeom, polySeq1);
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?
	// make sure to add the above to my code for design; increases stability with imm1 by quite a bit, ~10 from these runs of old CATM structures
	// maybe I should try to redesign gpa???
	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Output the starting geometry
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << opt.xShift << endl;
	cout << "crossingAngle: " << opt.crossingAngle << endl;
	cout << "axialRotation: " << opt.axialRotation << endl;
	cout << "zShift:        " << opt.zShift << endl;

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Assign number of rotamers by residue burial
	loadRotamers(sys, sysRot, "SL97.00");
	//loadRotamersBySASABurial(sys, sysRot, opt, rotamerSamplingPerPosition);
	//cout << sys.getAtomPointers() << endl;
	
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 

	/******************************************************************************
	 *                    === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
    map<string,double> monomerEnergyByTerm;
    double monomerEnergy = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);
	
	// Optimize Initial Starting Position 
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	//spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);// changed to make more similar to CATM
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	
	// Repack dimer
	repackSideChains(spm, opt.greedyCycles);
	vector<uint> startStateVec = spm.getMinStates()[0];
	double currentEnergy = spm.getMinBound()[0];
	sys.setActiveRotamers(startStateVec);

	// Output the starting energies	
	cout << spm.getSummary(startStateVec) << endl;

	// calculate the energy of the system
	double startDimer = sys.calcEnergy();
    cout << "Monomer Energy: " << monomerEnergy << endl;
    cout << "Dimer-Monomer: " << startDimer-monomerEnergy << endl;
	exit(0);
    
	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	/******************************************************************************
	 *                     === X SHIFT REPACKS ===
	 ******************************************************************************/
	double bestEnergy = currentEnergy-monomerEnergy;
	double savedXShift = opt.xShift;

	// xShift repack to get close to best xShift
	localXShiftDocking(sys, opt, bestEnergy, monomerEnergy, spm, helicalAxis, axisA, axisB, apvChainA, apvChainB, trans, savedXShift);
	//helicalAxis.applySavedCoor("originState");
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?
	
	repackSideChains(spm, opt.greedyCycles);
	vector<uint> xShiftStateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(xShiftStateVec);
	currentEnergy = spm.getMinBound()[0]-monomerEnergy;

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	cout << "Starting Geometry" << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << savedXShift << " crossingAngle: " << opt.crossingAngle << " axialRotation: " << opt.axialRotation << " zShift: " << opt.zShift << endl << endl;
	cout << "Current Best Energy: " << bestEnergy-monomerEnergy << endl;

	double bestRepackEnergy = sys.calcEnergy();
	vector<uint> bestRepackState;

	// Local backbone monte carlo repacks
	localBackboneRepack(opt, sys, savedXShift, spm, helicalAxis, axisA, axisB, apvChainA, apvChainB, trans, RNG, monomerEnergyByTerm, monomerEnergy, sout);

	// output the total time
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
    // close all of the output file writers
    sout.close();
    mout.close();
    err.close();
}

//Functions
//TODO: to add this into the design code, I think I just need to add in a mask here to the optimizer
// set a number of cycles
void localBackboneRepack(BBOptions &_opt, System &_sys, double _savedXShift, SelfPairManager &_spm,
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA, AtomPointerVector &_apvChainB,
 Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm, double _monomerEnergy, ofstream &_out){
	//for (int i=0; i < _opt.numRepacks; i++){
		if (_opt.verbose){
			cout << "==============================================" << endl;
			cout << " Performing Local Monte Carlo Backbone Repack " << endl;
			cout << "==============================================" << endl;
		}
		// set the system to the original xShift state
		_sys.applySavedCoor("savedBestState");
		_helicalAxis.applySavedCoor("BestAxis");
		// save the current state as the repack state
		_sys.saveAltCoor("savedRepackState");
		_helicalAxis.saveAltCoor("BestRepack");
		// get the best repack energy
		double prevBestEnergy = _sys.calcEnergy();

		// do backbone geometry repacks
		monteCarloRepack(_opt, _sys, _savedXShift, _spm, _helicalAxis, _axisA, _axisB, _apvChainA, _apvChainB, _trans, _RNG, prevBestEnergy, 
		 _monomerEnergyByTerm, _monomerEnergy, _out);
    
 	    //outputs a pdb file for the structure 
		// Initialize PDBWriter
		PDBWriter writer;
		writer.open(_opt.outputDir + "/backboneOptimized.pdb");
		writer.write(_sys.getAtomPointers(), true, false, true);
		writer.close();
	//}
}

void monteCarloRepack(BBOptions &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 map<string,double> _monomerEnergyByTerm, double _monomerEnergy, ofstream &_out){
	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);

	// starting geometry
	double xShift = _savedXShift;	
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	// Monte Carlo Repack Manager Setup
	MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MCMngr.setEner(_prevBestEnergy);

	vector<uint> startStateVec = _spm.getMinStates()[0];
	vector<unsigned int> MCOBest = startStateVec;
		
	unsigned int counter = 0;
	double currentEnergy = _prevBestEnergy;
	double startDimer = _prevBestEnergy;
	
	PDBWriter writer;
	// loop through the MC cycles for backbone repacks
	while(!MCMngr.getComplete()) {
		double startTemp = MCMngr.getCurrentT();

		_sys.applySavedCoor("savedRepackState");
		_helicalAxis.applySavedCoor("BestRepack");
		
		int moveToPreform = _RNG.getRandomInt(3);
		
		double deltaXShift = 0.0;
		double deltaZShift = 0.0;
		double deltaCrossingAngle = 0.0;
		double deltaAxialRotation = 0.0; 
		
		//======================================
		//====== Z Shift (Crossing Point) ======
		//======================================
		if (moveToPreform == 0) {
			//deltaZShift = getStandardNormal(RNG1) * 0.1;
			deltaZShift = getStandardNormal(_RNG) * _opt.deltaZ;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(_RNG) * _opt.deltaAx;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * _opt.deltaCross;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * _opt.deltaX;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		}
		
		// Run _optimization
		repackSideChains(_spm, _opt.greedyCycles);

		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0];
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE
		
		if (!MCMngr.accept(currentEnergy)) {
			if (_opt.verbose){
				cout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy-_monomerEnergy << endl;
			}
		} else {
			_prevBestEnergy = currentEnergy;
			_sys.saveAltCoor("savedRepackState");
			_helicalAxis.saveAltCoor("BestRepack");
		
			xShift = xShift + deltaXShift;
			crossingAngle = crossingAngle + deltaCrossingAngle;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift + deltaZShift;
			MCOBest = MCOFinal;

			// if accept, decrease the value of the moves by the sigmoid function
			double endTemp = MCMngr.getCurrentT();
			//double beta = 1.0 / (1.0 + exp(6.0)); // R(endT)
			//double gamma = 1.0 / (1.0 + exp(-6.0)) - beta; // R(startT) - R(endT)
			//double decreaseMultiplier = ((startTemp-endTemp)/(gamma) * (R * MCMngr.getCurrentCycle() - beta) + endTemp)/startTemp;
			//cout << "Start Temp: " << startTemp << " End Temp: " << endTemp << " R: " << R << " beta: " << beta << " gamma: " << gamma << endl;
			double decreaseMultiplier = endTemp/startTemp;
			_opt.deltaAx *= decreaseMultiplier;
			_opt.deltaX *= decreaseMultiplier;
			_opt.deltaZ *= decreaseMultiplier;
			_opt.deltaCross *= decreaseMultiplier;
			//cout << decreaseMultiplier << endl; // although this rounds when I print, the value still holds when I multiply, so should be good here
			//cout << _opt.deltaAx << " " << _opt.deltaX << " " << _opt.deltaZ << " " << _opt.deltaCross << endl;
			
			if (_opt.verbose){
				cout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-_monomerEnergy << endl;
			}
			counter++;
			writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	writer.close();
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
		
	_sys.applySavedCoor("savedRepackState");
	currentEnergy = _spm.getMinBound()[0];
	cout << "Best: " << _spm.getSummary(MCOBest) << endl;
	// TODO: make this into a map that saves all of these to be output
	double dimerEnergy = _sys.calcEnergy();
	double finalEnergy = dimerEnergy-_monomerEnergy;
	double vdw = _spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	double hbond = _spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	double imm1 = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1");
	double imm1Ref = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	double dimerDiff = dimerEnergy-startDimer;

	// get monomer energies from monomerEnergyTerm map
	double monomerVdw = _monomerEnergyByTerm.find("CHARMM_VDW")->second;
	double monomerHbond = _monomerEnergyByTerm.find("SCWRL4_HBOND")->second;
	double monomerImm1 = _monomerEnergyByTerm.find("CHARMM_IMM1")->second;
	double monomerImm1Ref =_monomerEnergyByTerm.find("CHARMM_IMM1REF")->second;

	// calculate the solvent accessible surface area
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	// Print out info to the summary csv file
	if (_opt.useElec == false){
		_out << "Sequence, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		_out << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		_out << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		cout << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		//cout << ',' << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		cout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		cout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		cout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	} else {
		double elec = _spm.getStateEnergy(MCOBest, "CHARMM_ELEC");
		double monomerElec =_monomerEnergyByTerm.find("CHARMM_ELEC")->second;
		_out << "Sequence, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, elec, monoElec, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		_out << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		_out << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		cout << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		//cout << ',' << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		cout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		cout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		cout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	}
}

void setGly69ToStartingGeometry(BBOptions &_opt, System &_sys, System &_helicalAxis,
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

void prepareSystem(BBOptions &_opt, System &_sys, System &_startGeom, string &_polySeq){	
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
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
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
	_sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	_sys.buildAllAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(_sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	CSB.updateNonBonded(10,12,50);
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
	if (_opt.useElec == true){
		Eset->setTermActive("CHARMM_ELEC", true);
		Eset->setWeight("CHARMM_ELEC", _opt.weight_elec);
	}

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// updated on 2023-8-26 to be more similar to sequence design
    deleteTerminalBondInteractions(_sys,_opt.deleteTerminalInteractions);

	// Up to here is from readDPBAndCalcEnergy
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	//CSB.updateNonBonded(10,12,50);

	// TODO: maybe calculate the energy here, then see if there's clashing. If so, then move helices away until no clashing, keeping the same
	// other coordinates? Just for simplicity for now. And if I want to implement this in design, adding this in will likely be a good idea.	
	// as of 2022-7-5: not sure if the above works or needs to be reworked
	PDBWriter writer1;
	writer1.open(_opt.outputDir + "/" + _opt.sequence + "_inputGeometry.pdb");
	writer1.write(_sys.getAtomPointers(), true, false, true);
	writer1.close();
}

//vector<pair <int, double>> getResiBurial(System &_sys, BBOptions &_opt, string _sequence){
//	vector<pair <int, double> > resiBurial = calculateResidueBurial(_sys, _opt, _sequence);
//	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
//		return left.second < right.second;
//	});
//	return resiBurial;
//}
//
//void defineRotamerLevels(System &_sys, BBOptions &_opt, vector<pair <int, double>> _resiBurial, 
// vector<uint> &_interfacePositions, string &_variablePositionString, string &_rotamerLevels){
//	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
//	int levelCounter = 0;
//	int numAAs = 0;
//	double lvlavg = 0;
//	vector<double> avgs;
//	//cout << "lvl " << levelCounter << endl;
//	for (uint i = 0; i < _resiBurial.size(); i++) {
//		double sasaPercentile = double(i) / double(_resiBurial.size());
//		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
//			levelCounter++;
//			lvlavg = lvlavg/numAAs;
//			avgs.push_back(lvlavg);
//			lvlavg=0;
//			numAAs=0;
//			//cout << "lvl " << levelCounter << endl;
//		}
//		//cout << _resiBurial[i].first << ": " << resiBurial[i].second << endl;
//		lvlavg = lvlavg+_resiBurial[i].second;
//		numAAs++;
//		int backbonePosition = _resiBurial[i].first;
//		Position &pos = _sys.getPosition(backbonePosition);
//		string posRot = _opt.sasaRepackLevel[levelCounter];
//		int resiNum = pos.getResidueNumber();
//		int posNum = resiNum-_opt.thread;
//		string add;
//
//		if (levelCounter < _opt.interfaceLevel){
//			add = "Add all Ids at this pos";
//			_interfacePositions.push_back(resiNum);
//			if (backbonePosition > 2 && backbonePosition < _opt.sequence.length()-4){//backbone position goes from 0-20, so numbers need to be 3 and 4 here instead of 4 and 5 to prevent changes at the interface like others
//				_variablePositionString.replace(_variablePositionString.begin()+posNum, _variablePositionString.begin()+posNum+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
//			}
//		} else {
//			add = "Only 1 ID";
//		}
//		_rotamerLevels.replace(_rotamerLevels.begin()+posNum, _rotamerLevels.begin()+posNum+1, MslTools::intToString(levelCounter));
//	}
//	lvlavg = lvlavg/numAAs;
//	avgs.push_back(lvlavg);
//}
//
//// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
//// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
//// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
//void defineInterfaceAndRotamerSampling(BBOptions &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq,
// string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
// vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out){
//	// Declare system
//	System sys;
//	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
//	CSB.setBuildTerm("CHARMM_ELEC", false);
//	CSB.setBuildTerm("CHARMM_ANGL", false);
//	CSB.setBuildTerm("CHARMM_BOND", false);
//	CSB.setBuildTerm("CHARMM_DIHE", false);
//	CSB.setBuildTerm("CHARMM_IMPR", false);
//	CSB.setBuildTerm("CHARMM_U-BR", false);
//
//	CSB.setBuildNonBondedInteractions(false);
//	//CSB.setBuildNoTerms();
//
//	if(!CSB.buildSystem(_PS)) {
//		cout << "Unable to build system from " << _PS << endl;
//		exit(0);
//	} else {
//		//fout << "CharmmSystem built for sequence" << endl;
//	}
//
//	/******************************************************************************
//	 *                     === INITIAL VARIABLE SET UP ===
//	 ******************************************************************************/
//	EnergySet* Eset = sys.getEnergySet();
//	// Set all terms active, besides Charmm-Elec
//	Eset->setAllTermsActive();
//	Eset->setTermActive("CHARMM_ELEC", false);
//	Eset->setTermActive("CHARMM_ANGL", false);
//	Eset->setTermActive("CHARMM_BOND", false);
//	Eset->setTermActive("CHARMM_DIHE", false);
//	Eset->setTermActive("CHARMM_IMPR", false);
//	Eset->setTermActive("CHARMM_U-BR", false);
//	Eset->setTermActive("CHARMM_VDW", true);
//	Eset->setTermActive("SCWRL4_HBOND", true);
//
//	// Set weights
//	Eset->setWeight("CHARMM_VDW", 1);
//	Eset->setWeight("SCWRL4_HBOND", 1);
//
//	CSB.updateNonBonded(10,12,50);
//
//	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
//	sysRot.defineRotamerSamplingLevels();
//
//	// Add hydrogen bond term
//	HydrogenBondBuilder hb(sys, _opt.hbondFile);
//	hb.buildInteractions(50);//when this is here, the HB weight is correct
//
//	/******************************************************************************
//	 *                     === COPY BACKBONE COORDINATES ===
//	 ******************************************************************************/
//	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
//	sys.buildAllAtoms();
//
//	// hardcoded Valine (two C beta carbons to better approximate interface interactions)
//	string backboneSeq = generateString("V", _opt.sequence.length());
//	vector<pair <int, double>> resiBurial = getResiBurial(sys, _opt, backboneSeq);
//
//	// Output variable Set up
//	string variablePositionString = generateString("0", _opt.sequence.length());
//	string rotamerLevels = generateString("0", _opt.sequence.length());
//
//	defineRotamerLevels(sys, _opt, resiBurial, _interfacePositions, variablePositionString, rotamerLevels);
//	backboneSeq = generateBackboneSequence("V", _opt.sequence.length(), _opt.useAlaAtCTerminus);
//
//	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
//	int highestRotamerLevel = numberOfRotamerLevels-1;
//
//	vector<int> rotamerSamplingPerPosition = getRotamerSampling(rotamerLevels);
//	vector<int> linkedPositions = getLinkedPositions(rotamerSamplingPerPosition, _opt.interfaceLevel, highestRotamerLevel);
//
//	//String for the positions of the sequences that are considered interface for positions amd high rotamers
//	string rotamerSamplingString = getInterfaceString(rotamerSamplingPerPosition, _opt.sequence.length());
//
//	// Define referenced output variables
//	_rotamerLevels = rotamerLevels;
//	//_polySeq = polySeq;
//	_rotamerSamplingPerPosition = rotamerSamplingPerPosition;
//	_variablePositionString = variablePositionString;
//	_rotamerSamplingString = rotamerSamplingString;
//	_linkedPositions = linkedPositions;
//
//	// output 
//	_out << "Sequence:           " << _opt.sequence << endl;
//	_out << "Rotamers Levels:    " << rotamerSamplingString << endl;
//	_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition);
//	_allInterfacePositions = getAllInterfacePositions(_opt, rotamerSamplingPerPosition);
//
//	cout << "Sequence:           " << _opt.sequence << endl;
//	cout << "Rotamers Levels:    " << rotamerSamplingString << endl;
//
//	int numPosAtInterface = 0;
//	for(string::iterator it = variablePositionString.begin(); it != variablePositionString.end();it++ ) {
//		stringstream ss;
//		ss << *it;
//		string num = ss.str();
//		if (MslTools::toInt(num) == 1){
//			numPosAtInterface++;
//		}
//	}
//	//cout << "Interface Positions: " << _interfacePositions.size() << endl;
//	//for (uint i=0; i<_interfacePositions.size(); i++){
//	//	cout << _interfacePositions[i] << ",";
//	//}
//	//cout << endl;
//}

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, BBOptions &_opt, string _seq) {
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
		AtomSelection selB(_sys.getPosition(i+_opt.sequence.length()).getAtomPointers());
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

vector<uint> getAllInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=0; k<_opt.sequence.length(); k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=3; k<_opt.sequence.length()-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

void localXShiftDocking(System &_sys, BBOptions &_opt, double &_bestEnergy, double _monomerEnergy, 
 SelfPairManager &_spm, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB,
 AtomPointerVector &_apvChainA, AtomPointerVector &_apvChainB, Transforms &_trans, double &_savedXShift) {
	double deltaXShift = -0.1; // xShift changes
	double xShiftEnd = 6; // xShift to stop at
	double xShift = _opt.xShift; // current xShift
	double previousEnergy = _monomerEnergy; // previous energy to compare to
	double globalLowestE = _monomerEnergy; // Global lowest energy found (if above monomer we won't save anyways)
	double currentEnergy = _bestEnergy; // current energy
	// while loop for x-shifts	
	while (xShift >= xShiftEnd) {
		// add the xShift change to the current xShift
		xShift += deltaXShift;
		// Move the helix
		backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaXShift, 3 );

		// Run Optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal;
		MCOFinal = _spm.getMinStates()[0];
		_sys.setActiveRotamers(MCOFinal);

		// get current energy
		currentEnergy = _spm.getMinBound()[0];

		// Check if this is the lowest energy found so far
		if (currentEnergy < _bestEnergy) {
			_bestEnergy = currentEnergy;
			_savedXShift = xShift;
			_sys.saveAltCoor("savedBestState");
			_helicalAxis.saveAltCoor("BestAxis");
		}
		cout << "xShift: " << xShift << " energy: " << currentEnergy-_monomerEnergy << endl;

		// If energy increase twice in a row, and it is above the global lowest energy, quit
		if (currentEnergy < globalLowestE) {
			globalLowestE = currentEnergy;
		}
		if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
			cout << "Energy increasing above global lowest energy... (currently " << globalLowestE-_monomerEnergy << ")" << endl;
			break;
		} else {
			previousEnergy = currentEnergy;
		}
	}
	cout << "Best Energy at x shift: " << _bestEnergy-_monomerEnergy << " at " << _savedXShift << endl;
	PDBWriter writer1;
	writer1.open(_opt.outputDir + "/" + _opt.sequence + "_xShifted.pdb");
	writer1.write(_sys.getAtomPointers(), true, false, true);
	writer1.close();
}

/***********************************
 *help functions
 ***********************************/
void usage() {
	cout << endl;
	cout << "Run as :" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputWarningMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program has the following warning:" << endl;
		cerr << endl;
		cerr << _opt.warningMessages << endl;
		cerr << endl;
}
void outputErrorMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << _opt.errorMessages << endl;
		cerr << endl;
		cerr << _opt.OPerrors << endl;
		usage();
}

//TODO: finish writing up this help
void help(BBOptions defaults) {
	cout << "This program runs as:" << endl;
	cout << " % seqDesign " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --numberOfStructuresToMCRepack <int> --energyCutOff <double> --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --MCStartTemp <double> --MCEndTemp <double> --MCCurve <CONSTANT-0, LINEAR-1, EXPONENTIAL-2, SIGMOIDAL-3, SOFT-4>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --thread <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   --weight_hbond <double> --weight_vdw <double> --weight_solv <double> --weight_seqEntropy <double>" << endl;
	cout << "   --sasaRepackLevel <rotLevel> (in format SL95.00; 4 levels used by default) --interfaceLevel <int> " << endl << endl;
	cout << "Template Configuration file (copy and paste the below into a file.config and run code as bin/seqDesign --config file.config" << endl;
	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;

	cout << endl << "#Energy term weights" << endl;
	cout << setw(20) << "weight_vdw " << defaults.weight_vdw << endl;
	cout << setw(20) << "weight_hbond " << defaults.weight_hbond << endl;
	cout << setw(20) << "weight_solv " << defaults.weight_solv << endl;
	cout << endl;
}

// check through error options and exit if too many
void checkOptionErrors(BBOptions &_opt){	
	if (_opt.errorFlag) {
		outputErrorMessage(_opt);
		exit(1);
	} else if (!_opt.errorFlag && !_opt.warningFlag && _opt.errorMessages != ""){
		outputWarningMessage(_opt);
		usage();
		exit(0);
	}
}