#include <iostream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>
#include <chrono>

// MSL Functions
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "CRDReader.h"
#include "SasaCalculator.h"

// My functions
#include "versatileFunctions.h"
#include "designFunctions.h"
#include "designOptions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "forwardFolding";
string programDescription = "Runs an input sequence through a set of given backbones and determines the best geometry";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "25 January 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime, spmTime;
auto clockTime = chrono::system_clock::now();

// initialize time variables
time_t rawtime;
struct tm * timeinfo;
char buffer[80];


// Functions
/*
	I have left the most important functions within this code in case they need to be looked at or changed.
	Many of the auxiliary functions are found within the designFunctions.cpp file.
*/

//void switchSequence(System &_sys, Options &_opt, string _sequence = "AAALLLLLLLLLLLLLLLAAA");
void switchSequence(System &_sys, Options &_opt, string _sequence);

// geometry setup functions
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS, bool _useBaseline);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, Transforms &_trans);
void getStartingGeometry(Options &_opt, ofstream &_sout);
void checkForClashing(System &_startGeom, Options &_opt, vector<uint> _interfacePositions, ofstream &_sout);

// define interface and rotamer level functions
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<uint> &_rotamerSamplingPerPosition);
void defineRotamerLevels(Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString);
void useInputInterface(Options &_opt, string &_variablePositionString, string &_rotamerLevels, vector<int> &_interfacePositions);

// backbone optimization functions
void backboneOptimizer(Options &_opt, System &_startGeom, string _sequence, vector<uint> _bestState, map<string, map<string,double>> &_sequenceEnergyMap,
 uint _rep, System &_helicalAxis, AtomPointerVector &_axisA,  AtomPointerVector &_axisB, vector<uint> _rotamerSampling, Transforms &_trans,
 RandomNumberGenerator &_RNG, ofstream &_sout);
double backboneOptimizeMonteCarlo(Options &_opt, System &_sys, SelfPairManager &_spm, map<string, map<string,double>> &_sequenceEnergyMap,
 string _sequence, vector<uint> &_bestState, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 uint _rep, ofstream &_sout);
void getCurrentMoveSizes(Options &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);

// output functions
void outputFiles(Options &_opt, double _seed, vector<uint> _rotamerSamplingPerPosition,
 map<string,map<string,double>> _sequenceEnergyMap, ofstream &_sout);
string getRunParameters(Options &_opt, vector<double> _densities);
void outputTime(auto _clockTime, string _descriptor, ofstream &_sout);

// axial rotation and zShift conversion from given input ax' and z' to ax and z (Mueller 2014; Fig. S1)
void convertToRelativeAxAndZ(double _axialRotation, double _zShift, double &_relativeAx, double &_relativeZ);
void convertToAxAndZForTranformation(Options &_opt);
map<string, double> getGeometryMap(Options &opt, string _descriptor);
void addGeometryToEnergyMap(map<string, double> _geometryMap, map<string, map<string, double>> &_energyMap, string _bestSequence);

void computeMonomerEnergy(System &_sys, System &_helicalAxis, Options &_opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq,
 RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
double calculateInterfaceSequenceEntropy(string _sequence, map<string, double> _sequenceEntropyMap, vector<uint> _interfacePositions,
 double &_sequenceProbability, double _seqEntropyWeight);
void outputStartDesignInfo(Options &_opt);

// help functions
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

/******************************************
 *  =======  BEGIN MAIN =======
 ******************************************/
int main(int argc, char *argv[]){
	// start the timer for the program
	time(&startTime);

	// setup time with the current time
	time_t startRealTime = chrono::system_clock::to_time_t(clockTime); 
	cout << "Program Start time: " << ctime(&startRealTime) << endl;

	// get date
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);
	
	// parse command line options
	Options opt = parseOptions(argc, argv);
	if (opt.errorFlag) {
		outputErrorMessage(opt);
		exit(1);
	} else if (!opt.errorFlag && !opt.warningFlag && opt.errorMessages != ""){
		outputErrorMessage(opt);
		usage();
		exit(0);
	}
	
	setupDesignDirectory(opt); // makes a directory in the directory that you run from, with design_<runNumber> as the name
	
	// setup the rerun config file
	ofstream rerun; // rerun config output
	string rerunfile = opt.outputDir + "/rerun.config"; // rerun config file name
	rerun.open(rerunfile.c_str()); // open rerun config file
	rerun << opt.rerunConf << endl; // write rerun config file with config options
	rerun.close(); // close rerun config file

	// setup the summary output file
	ofstream sout; // summary file output
	string soutfile = opt.outputDir + "/summary.out"; // summary file name
	sout.open(soutfile.c_str()); // open summary file
	sout << date << endl; // output the date and time to the summary file

	// setup the error file
	ofstream err; // error file output
	string errfile  = opt.outputDir + "/errors.out"; // error file name
	err.open(errfile.c_str()); // open error file

	// get the starting geometries; convert to parallelogram axialRot and Z (Mueller 2014; Fig. S1)
	convertToAxAndZForTranformation(opt);

	// output starting design information
	outputStartDesignInfo(opt);

	// setup for saving the best sets of geometries 
	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	// AtomPointerVector for the helical axis for each chain
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// compute monomer energy
	computeMonomerEnergy(sys, helicalAxis, opt, trans, sequenceEnergyMapBest, bestSequence, RNG, sout, err);

	// start multithreading through different backbones here
	/******************************************************************************
	 *      === COPY BACKBONE COORDINATES AND TRANSFORM TO INPUT GEOMETRY ===
	 ******************************************************************************/
	// get the starting geometry using polyglycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,axisA,axisB,trans);

	// I think I can rid of the interfacial stuff?	
	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	vector<uint> interfacePositions; // vector of positions at the interface excluding termini positions
	vector<uint> allInterfacePositions; // vector of positions at the interface including the terminal positions
	vector<uint> rotamerSamplingPerPosition; // vector of rotamer level for each position

	// Defines the interfacial positions and the number of rotamers to give each position
	PolymerSequence interfacePolySeq = getInterfacialPolymerSequence(opt, startGeom, allInterfacePositions, interfacePositions,
	 rotamerSamplingPerPosition);
	
	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	prepareSystem(opt, sys, startGeom, interfacePolySeq, opt.useBaseline);
	checkForClashing(sys, opt, interfacePositions, sout); // checks an interface for clashes; if too clashing (>10 VDW) with alanine at those positions, exit
	
	// initialize the object for loading rotamers into system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	switchSequence(sys, opt, opt.sequence);

	/******************************************************************************
	 *           === VARIABLES FOR SAVING ENERGIES AND SEQUENCES ===
	 ******************************************************************************/
	map<string, map<string,double>> sequenceEnergyMapBest; // energyMap to hold all energies for output into a summary file
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile); // Get sequence entropy map

	/******************************************************************************
	 *                    === GET THE BEST STARTING STATE ===
	 ******************************************************************************/
	// Initialize RNG with seed (time or given seed number)
	RandomNumberGenerator RNG;
	if (opt.useTimeBasedSeed){
		RNG.setTimeBasedSeed();
	} else {
		RNG.setSeed(opt.seed);
	}
	
	// if runSCMF is true, get the best starting state using SelfConsistentMeanField
	// add in getting energy using self pair manager here
	if (opt.runSCMF == true){
		runSCMFToGetStartingSequence(sys, opt, RNG, rotamerSamplingPerPosition, 
		 interfacePositions, sequenceEnergyMapBest, sequenceEntropyMap, sout);
	}

	// get the starting sequence	
	Chain &chainA = sys.getChain("A");
	string seq = convertPolymerSeqToOneLetterSeq(chainA);
	cout << "Sequence before stateMC: " << seq << endl;

	/******************************************************************************
	 *   === MONTE CARLO TO SEARCH FOR BEST SEQUENCES AND BACKBONE OPTIMIZE ===
	 ******************************************************************************/
	map<string, map<string,double>> sequenceEnergyMapFinalSeqs; // energyMap to hold all energies for output into a summary file
	string bestSequence = seq;
	string prevSequence;
	map<string, double> startGeometries = getGeometryMap(opt, "start");
	addGeometryToEnergyMap(startGeometries, sequenceEnergyMapBest, bestSequence);

	for (uint i=0; i<opt.backboneSearchCycles; i++){
		vector<uint> bestState;
		// run state Monte Carlo to get a random sequence from the best state
		//searchForBestSequence(startGeom, opt, interfacePolySeq, sequenceEnergyMapBest, sequenceEntropyMap, bestState, bestSequence,
		// allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, i, sout, err);
		// add in the starting geometries to the map	
		addGeometryToEnergyMap(startGeometries, sequenceEnergyMapBest, bestSequence);
		// TODO: going to set up this today; save x number of sequences from the search and then optimize them with threading
		// if the sequence is different, optimize the backbone
		
		// optimize the backbone for the sequence
		backboneOptimizer(opt, startGeom, bestSequence, bestState, sequenceEnergyMapBest, i, helicalAxis, axisA, axisB,
		 rotamerSamplingPerPosition, trans, RNG, sout);
		// get the best sequence from the energy map
		sequenceEnergyMapFinalSeqs[bestSequence] = sequenceEnergyMapBest[bestSequence];	
		// reset the energy map
		sequenceEnergyMapBest.clear(); // energyMap to hold all energies for output into a summary file
	}

	/******************************************************************************
	 *                   === WRITE OUT ENERGY AND DESIGN FILES ===
	 ******************************************************************************/
	double seed = RNG.getSeed();
	outputFiles(opt, seed, rotamerSamplingPerPosition, sequenceEnergyMapFinalSeqs, sout);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

	err.close();
	sout.close();
}

//Functions
double calculateInterfaceSequenceEntropy(string _sequence, map<string, double> _sequenceEntropyMap, vector<uint> _interfacePositions,
 double &_sequenceProbability, double _seqEntropyWeight){
	//Get residue name for each interfacial identity
	vector<string> seqVector;
	int numInterfacials = _interfacePositions.size();
	for (uint i=0; i<numInterfacials; i++){
		stringstream tmp;
		tmp << _sequence[_interfacePositions[i]];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
		//cout << _interfacialPositionsList[i] << " : " << resName << endl;
	}
	map<string,int> seqCountMap = getAACountMap(seqVector); // get the count of each amino acid in the sequence
	double numberOfPermutations = calcNumberOfPermutations(seqCountMap, numInterfacials); // calculate the number of permutations for the sequence
	_sequenceProbability = calculateSequenceProbability(seqCountMap, _sequenceEntropyMap, numberOfPermutations); // calculate the sequence probability
	double seqEntropy = -log(_sequenceProbability)*0.592*_seqEntropyWeight; // multiply by the weight and -1
	//double seqEntropy = _seqEntropyWeight/_sequenceProbability; // divide the weight by the probability to get the entropy
	return seqEntropy;
}

map<string, double> getGeometryMap(Options &_opt, string _descriptor){
	map<string, double> geometryMap;
	geometryMap[_descriptor+"XShift"] = _opt.xShift;
	geometryMap[_descriptor+"CrossingAngle"] = _opt.crossingAngle;
	geometryMap[_descriptor+"AxialRotation"] = _opt.axialRotation;
	geometryMap[_descriptor+"ZShift"] = _opt.zShift;
	double relativeAx;
	double relativeZ;
	convertToRelativeAxAndZ(_opt.axialRotation, _opt.zShift, relativeAx, relativeZ);
	geometryMap[_descriptor+"AxialRotationPrime"] = relativeAx;
	geometryMap[_descriptor+"ZShiftPrime"] = relativeZ;
	return geometryMap;
}

void outputGeometry(Options &_opt, double _xShift, double _crossingAngle, double _axialRotation, double _zShift, ofstream &_sout){
	double relativeAxBefore;
	double relativeZBefore;
	convertToRelativeAxAndZ(_opt.axialRotation, _opt.zShift, relativeAxBefore, relativeZBefore);
	double relativeAxAfter;
	double relativeZAfter;
	convertToRelativeAxAndZ(_axialRotation, _zShift, relativeAxAfter, relativeZAfter);
	// convert axialRotation and zShift for output (from parallelogram from Ben's paper for interface to square)
	_sout << "xShift;        Before: " << _opt.xShift << "; After: " << _xShift << endl;
	_sout << "crossingAngle; Before: " << _opt.crossingAngle << "; After: " << _crossingAngle << endl;
	_sout << "axialRotation; Before: " << _opt.axialRotation << "; After: " << _axialRotation << endl;
	_sout << "zShift;        Before: " << _opt.zShift << "; After: " << _zShift << endl << endl;
	_sout << "axialRotationPrime; Before: " << relativeAxBefore << "; After: " << relativeAxAfter << endl;
	_sout << "zShiftPrime;        Before: " << relativeZBefore << "; After: " << relativeZAfter << endl << endl;
}

void addGeometryToEnergyMap(map<string, double> _geometryMap, map<string, map<string, double>> &_energyMap, string _bestSequence){
	for (auto &it : _geometryMap){
		_energyMap[_bestSequence][it.first] = it.second;
	}
}

void convertToRelativeAxAndZ(double _axialRotation, double _zShift, double &_relativeAx, double &_relativeZ){
	// use the positive axial rotation for conversion since our axial rotations are negative
	//double axialRotation = abs(_axialRot);
	double axialRotation = _axialRotation;
	// use equations from Mueller 2014; Fig. S1
	_relativeAx = (10*axialRotation/9)+(200*_zShift/27);
	_relativeAx = 100+_relativeAx;
	_relativeZ = (10*_zShift/9)+(0.15*axialRotation/9);
}

void convertToAxAndZForTranformation(Options &_opt){
	double axialRotation = abs(_opt.axialRotation);
	double zShift = _opt.zShift;
	// convert input axialRotation and zShift to transformation for interfacial parallelogram (Mueller 2014; Fig S1) 
	_opt.axialRotation = axialRotation+(20*zShift/3);
	_opt.axialRotation = -_opt.axialRotation;
	_opt.zShift = zShift+(0.015*axialRotation);
}

void getStartingGeometry(Options &_opt, ofstream &_sout){
	// use the given xShift, crossing angle, axial rotation, and zShift
	_sout << "***STARTING GEOMETRY:***" << endl;
	_sout << "xShift:        " << _opt.xShift << endl;
	_sout << "crossingAngle: " << _opt.crossingAngle << endl;
	convertToAxAndZForTranformation(_opt);
	_sout << "axialRotation: " << _opt.axialRotation << endl;
	_sout << "zShift:        " << _opt.zShift << endl << endl;
	_sout << "axRotAndZShift:" << _opt.density << endl << endl;
}

// switches to the starting sequence (if given, otherwise set to polyleu)
void switchSequence(System &_sys, Options &_opt, string _sequence){
	if (_sequence == ""){
		_sequence = generateBackboneSequence(_opt.backboneAA, _opt.backboneLength, _opt.useAlaAtTermini);
	}
	setActiveSequence(_sys, _sequence);
}

// The below checks for clashing at the interface by looking at the energy for a polyAla interface sequence
void checkForClashing(System &_startGeom, Options &_opt, vector<uint> _interfacePositions, ofstream &_sout){
	// change the ends of sequence to ala
	string polyLeu;
	for (uint i = 0; i < _opt.backboneLength; i++){
		if (i < 4 || i > _opt.backboneLength-5){
			polyLeu += "A";
		} else {
			polyLeu += "L";
		}
	}
	if (_opt.xShift <= 7.5){
		for (uint i=0; i< _interfacePositions.size(); i++){
			if (i%2==0){
				polyLeu[_interfacePositions[i]] = 'G';
			} else {
				polyLeu[_interfacePositions[i]] = 'A';
			}
		}
	} else {
		for (uint i=0; i< _interfacePositions.size(); i++){
			polyLeu[_interfacePositions[i]] = 'A';
		}
	}

	// convert sequence to polymer sequence
	string backboneSeq = convertToPolymerSequence(polyLeu, _opt.thread);
	PolymerSequence PS(backboneSeq);
	
	// declare system
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << PS << endl;
		exit(0);
	}
	
	// assign the coordinates of our system to the given geometry 
	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	sys.buildAllAtoms();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = sys.getEnergySet();

	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	CSB.updateNonBonded(10,12,50);

	AtomSelection sel(sys.getAtomPointers());
	sel.select("chainA, chain A");
	sel.select("chainB, chain B");
	double energy = sys.calcEnergy("chainA", "chainB");
	if (energy > 10){
		_sout << "Clashing at the interface; energy = " << energy << endl;
		cout << "Clashing at the interface; energy = " << energy << endl;
		exit(0);
	} 
	
	// if no clashing, then output the pdb
	// convert axial rotation to positive using absolute value, for outputting
	double relativeAx;
	double relativeZ;
	convertToRelativeAxAndZ(_opt.axialRotation, _opt.zShift, relativeAx, relativeZ);
	string geometry = "x"+MslTools::doubleToString(_opt.xShift)+"_cross"+MslTools::doubleToString(_opt.crossingAngle)
	 +"_ax"+MslTools::doubleToString(relativeAx)+"_z"+MslTools::doubleToString(relativeZ)+"_vdW"+to_string(energy)+".pdb";
	writePdb(sys, _opt.outputDir, geometry);
}

// function to prepare the system for design:
// - sets up CharmmSystemBuilder
// - assigns to given coordinates
// - Adds in the correct energies for the EnergySet
// - deletes terminal hydrogen bonds
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS, bool _useBaseline){	
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
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(_PS)) {
		cerr << "Unable to build system from " << _PS << endl;
		exit(0);
	}
	
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
	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
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

	_sys.calcEnergy();
	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
    deleteTerminalBondInteractions(_sys,_opt.deleteTerminalInteractions);

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	// add in an estimate of monomer energy for each sequence (calculated using CHARMM_VDW, SCWRL4_HBOND, CHARMM_IMM1, and CHARMM_IMM1REF)
	if (_useBaseline){
		buildBaselines(_sys, _opt.selfEnergyFile, _opt.pairEnergyFile);
	}
}

// sets the gly69 backbone to starting geometry
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, Transforms &_trans) {
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
	
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, _axisA, _axisB, ori, xAxis, zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

/****************************************************************
 *           === BACKBONE OPTIMIZER FUNCTIONS ===
 ****************************************************************/
// backbone optimization function, setting up the system with a single sequence to run the optimization
void backboneOptimizer(Options &_opt, System &_startGeom, string _sequence, vector<uint> _bestState, map<string, map<string,double>> &_sequenceEnergyMap,
 uint _rep, System &_helicalAxis, AtomPointerVector &_axisA,  AtomPointerVector &_axisB, vector<uint> _rotamerSampling, Transforms &_trans,
 RandomNumberGenerator &_RNG, ofstream &_sout){
	// output the start time
	outputTime(clockTime, "Backbone optimize replicate " + to_string(_rep) + " Start", _sout);

	// get the polymer sequence for the input sequence
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PS(polySeq);

	// set up the system for the input sequence
	System sys;
	prepareSystem(_opt, sys, _startGeom, PS, false);
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// load the rotamers
	loadRotamers(sys, sysRot, _opt.SL);

	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// Setup the SelfPairManager object
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	spm.runGreedyOptimizer(_opt.greedyCycles);

	// get monomer energy
	double monomerEnergy = _sequenceEnergyMap[_sequence]["Monomer"];
	// running into a weird issue with polyala, so just getting the best state which should work
	if (_opt.useAlaAtTermini == true){
		_bestState = spm.getMinStates()[0];
	}
	double currentEnergy = spm.getStateEnergy(_bestState)-monomerEnergy;
	double dimer = spm.getStateEnergy(_bestState); // best state is too long when rotamer level 60, but works otherwise
	double calcDimer = sys.calcEnergy();
	sys.setActiveRotamers(_bestState);
	cout << "Energy Before Backbone Optimization: " << currentEnergy << endl;

	// save the pre backbone optimize energy and the replicate number in the 
	_sequenceEnergyMap[_sequence]["preOptimizeEnergy"] = currentEnergy;
	_sequenceEnergyMap[_sequence]["replicateNumber"] = _rep;
	
	// add the geometry to the map pre backbone optimization
	map<string, double> preOptimizeGeometry = getGeometryMap(_opt, "preOptimize");
	addGeometryToEnergyMap(preOptimizeGeometry, _sequenceEnergyMap, _sequence);
	// TODO: thread the below, run multiple, and get the one with the best energy
	double finalEnergy = backboneOptimizeMonteCarlo(_opt, sys, spm, _sequenceEnergyMap, _sequence, _bestState, _helicalAxis, _axisA, _axisB, apvChainA, apvChainB, _trans, _RNG, monomerEnergy, _rep, _sout);
	// add the geometry to the map post backbone optimization
	map<string, double> endGeometry = getGeometryMap(_opt, "end");
	addGeometryToEnergyMap(endGeometry, _sequenceEnergyMap, _sequence);

	// write the pdb for the final state of the backbone optimization
	string designNumber = _opt.runNumber+"_"+to_string(_rep);
	writePdb(sys, _opt.outputDir, designNumber);

	// assign the coordinates of our system to the given geometry 
	_startGeom.assignCoordinates(sys.getAtomPointers(),false);
	_startGeom.buildAllAtoms();

	// output the energy of the system post backbone repack
	double repackEnergy = _sequenceEnergyMap[_sequence]["Total"];
	cout << "Energy After BBoptimize: " << repackEnergy << endl;

	// output the end time
	outputTime(clockTime, "Backbone optimize replicate " + to_string(_rep) + " End", _sout);
}

// monte carlo backbone optimization function
double backboneOptimizeMonteCarlo(Options &_opt, System &_sys, SelfPairManager &_spm, map<string, map<string,double>> &_sequenceEnergyMap,
 string _sequence, vector<uint> &_bestState, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 uint _rep, ofstream &_sout){
	// Setup backbone repack file
	ofstream bbout;
	string bboutfile  = _opt.outputDir + "/bbRepack_" + to_string(_rep) + ".out";
	bbout.open(bboutfile.c_str());

	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	// starting geometry
	double xShift = _opt.xShift;	
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	// final geometry
	double finalXShift = _opt.xShift;
	double finalCrossingAngle = _opt.crossingAngle;
	double finalAxialRotation = _opt.axialRotation;
	double finalZShift = _opt.zShift;

	bbout << "***STARTING GEOMETRY***" << endl;
	outputGeometry(_opt, xShift, crossingAngle, axialRotation, zShift, bbout);
	// calculate starting sasa
	SasaCalculator startSasa(_sys.getAtomPointers());
	startSasa.calcSasa();
	double sasa = startSasa.getTotalSasa();
	double monomerSasa = _sequenceEnergyMap[_sequence]["MonomerSasa"];
	_sequenceEnergyMap[_sequence]["PreBBOptimizeSasa"] = sasa;
	cout << "Pre BBOptimize SASA: " << sasa << endl;
	double bestSasa = sasa-monomerSasa;

	// Monte Carlo Repack Manager Setup
	MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects, _opt.backboneConvergedSteps, _opt.backboneConvergedE);

	vector<unsigned int> MCOBest = _bestState;
	
	unsigned int counter = 0;
	_sys.setActiveRotamers(_bestState);
	double currentEnergy = _spm.getStateEnergy(_bestState)-_monomerEnergy;
	double dimer = _spm.getStateEnergy(_bestState);
	double calcDimer = _sys.calcEnergy();
	cout << "Starting Energy: " << currentEnergy << endl;
	cout << "Starting Dimer Energy: " << dimer << endl;
	cout << "Monomer Energy: " << _monomerEnergy << endl;
	cout << "Calculated Dimer Energy: " << calcDimer << endl;
	_sequenceEnergyMap[_sequence]["TotalPreOptimize"] = currentEnergy;
	_sequenceEnergyMap[_sequence]["VDWDimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_VDW");
	_sequenceEnergyMap[_sequence]["IMM1DimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1")+_spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	_sequenceEnergyMap[_sequence]["HBONDDimerPreOptimize"] = _spm.getStateEnergy(_bestState, "SCWRL4_HBOND");
	
	double bestEnergy = currentEnergy;
	double prevBestEnergy = currentEnergy;
	if (_opt.compareSasa){
		MCMngr.setEner(bestSasa);
	} else {
		MCMngr.setEner(currentEnergy);
	}
	//double startDimer = _prevBestEnergy;

	// setup variables for shifts: ensures that they start from the proper values for every repack and not just the final value from the initial repack
	bool decreaseMoveSize = _opt.decreaseMoveSize;
	double deltaX = _opt.deltaX;
	double deltaCross = _opt.deltaCross;
	double deltaAx = _opt.deltaAx;
	double deltaZ = _opt.deltaZ;

	// uncomment the below and add the end of the accept to output the geometry trajectory pdb
	//PDBWriter writer;
	//writer.open(_opt.outputDir + "/bbRepack_"+to_string(_rep)+".pdb");
	// loop through the MC cycles for backbone repacks
	bbout << "Starting Repack Cycles" << endl; 
	while(!MCMngr.getComplete()) {
		// get the current temperature of the MC
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
			deltaZShift = getStandardNormal(_RNG) * deltaZ;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(_RNG) * deltaAx;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaAxialRotation, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * deltaCross;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaCrossingAngle, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * deltaX;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaXShift, moveToPreform);
		}

		// Run optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0]-_monomerEnergy;
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE
		if (_opt.compareSasa){
			SasaCalculator startSasa(_sys.getAtomPointers());
			startSasa.calcSasa();
			double sasa = startSasa.getTotalSasa()-monomerSasa;
			if (!MCMngr.accept(sasa)) {
				bbout << "MCReject   xShift: " << finalXShift+deltaXShift << " crossingAngle: " << finalCrossingAngle+deltaCrossingAngle << " axialRotation: " << finalAxialRotation+deltaAxialRotation << " zShift: " << finalZShift+deltaZShift << " sasa: " << sasa << " energy: " << currentEnergy << endl;
			} else {
				bestEnergy = currentEnergy;
				bestSasa = sasa;
				_sys.saveAltCoor("savedRepackState");
				_helicalAxis.saveAltCoor("BestRepack");
		
				finalXShift = finalXShift + deltaXShift;
				finalCrossingAngle = finalCrossingAngle + deltaCrossingAngle;
				finalAxialRotation = finalAxialRotation + deltaAxialRotation;
				finalZShift = finalZShift + deltaZShift;
				MCOBest = MCOFinal;

				// if accept, decrease the value of the moves by the sigmoid function
				if (_opt.decreaseMoveSize == true){
					double endTemp = MCMngr.getCurrentT();
					getCurrentMoveSizes(_opt, startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, decreaseMoveSize);
				}
				bbout << "MCAccept " << counter <<  " xShift: " << finalXShift << " crossingAngle: " << finalCrossingAngle << " axialRotation: " << finalAxialRotation << " zShift: " << finalZShift << " sasa: " << sasa << " energy: " << currentEnergy << endl;
				counter++;
			}
		} else {
			if (counter == 0){
				_sequenceEnergyMap[_sequence]["firstRepackEnergy"] = currentEnergy;
			}
			if (!MCMngr.accept(currentEnergy)) {
				bbout << "MCReject   xShift: " << finalXShift+deltaXShift << " crossingAngle: " << finalCrossingAngle+deltaCrossingAngle << " axialRotation: " << finalAxialRotation+deltaAxialRotation << " zShift: " << finalZShift+deltaZShift << " energy: " << currentEnergy << endl;
			} else {
				bestEnergy = currentEnergy;
				_sys.saveAltCoor("savedRepackState");
				_helicalAxis.saveAltCoor("BestRepack");
		
				finalXShift = finalXShift + deltaXShift;
				finalCrossingAngle = finalCrossingAngle + deltaCrossingAngle;
				finalAxialRotation = finalAxialRotation + deltaAxialRotation;
				finalZShift = finalZShift + deltaZShift;
				MCOBest = MCOFinal;

				// if accept, decrease the value of the moves by the sigmoid function
				if (_opt.decreaseMoveSize == true){
					double endTemp = MCMngr.getCurrentT();
					getCurrentMoveSizes(_opt, startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, decreaseMoveSize);
				}
				bbout << "MCAccept " << counter <<  " xShift: " << finalXShift << " crossingAngle: " << finalCrossingAngle << " axialRotation: " << finalAxialRotation << " zShift: " << finalZShift << " energy: " << currentEnergy << endl;
				counter++;
				//writer.write(_sys.getAtomPointers(), true, false, true);
			}
		}
	}
	//writer.close();
	bbout << "End Repack Cycles" << endl << endl; 
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_bestState = MCOBest;	
	_sys.applySavedCoor("savedRepackState");
	double dimerEnergy = _spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-_monomerEnergy;
	
	// Output change in geometry
	bbout << "***REPACK GEOMETRY***" << endl;
	outputGeometry(_opt, finalXShift, finalCrossingAngle, finalAxialRotation, finalZShift, bbout);
	bbout << "Energy;        Before: " << prevBestEnergy << "; After: " << bestEnergy << endl << endl;

	SasaCalculator endDimerSasa(_sys.getAtomPointers());
	endDimerSasa.calcSasa();
	double endSasa = endDimerSasa.getTotalSasa();
	_sequenceEnergyMap[_sequence]["OptimizeSasa"] = endSasa;
	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;
	_sequenceEnergyMap[_sequence]["VDWDimerOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_VDW");
	_sequenceEnergyMap[_sequence]["IMM1DimerOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1")+_spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	_sequenceEnergyMap[_sequence]["HBONDDimerOptimize"] = _spm.getStateEnergy(_bestState, "SCWRL4_HBOND");

	// sets the updated backbone parameters
	_opt.xShift = finalXShift;
	_opt.crossingAngle = finalCrossingAngle;
	_opt.axialRotation = finalAxialRotation;
	_opt.zShift = finalZShift;
	bbout << MCMngr.getReasonCompleted() << endl;	
	bbout << "Monte Carlo repack complete. Time: " << diffTimeMC/60 << "min" << endl << endl;
	if (finalEnergy > 100){
        bbout << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		_sout << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		exit(0);
	}
	return finalEnergy;
}

void outputFiles(Options &_opt, double _seed, vector<uint> _rotamerSamplingPositionVector,
 map<string,map<string,double>> _sequenceEnergyMap, ofstream &_sout){
	// Setup vector to hold energy file lines
	vector<string> energyLines;
	// get the run parameters
	string t = "\t";
	stringstream enerTerms;
	// For loop to setup the energy file
	uint i = 0;
	string rotamerValues = convertVectorUintToString(_rotamerSamplingPositionVector); // string of rotamer sampling number (if 4 rotamer levels, 0-3 for each position)
	for (auto &seq : _sequenceEnergyMap){
		stringstream seqLine;
		string sequence = seq.first;
		// get the interface sequence
		string interfaceSequence = getInterfaceSequence(_opt.interfaceLevel, rotamerValues, sequence);
		seqLine << sequence << t << rotamerValues << t << _opt.interface << t << interfaceSequence << t << _seed << t;
		map<string,double> energyMap = _sequenceEnergyMap[sequence];
		// For adding in strings to a line for the energy file; looping through the terms instead of my input terms this time; sort later
		for (auto &ener: energyMap){
			string energyTerm = ener.first;
			double energy = energyMap[energyTerm];
			string term = MslTools::doubleToString(energy)+t;
			seqLine << term;
			if (i == 0){
				enerTerms << energyTerm << t;
			}
		}
		string line = seqLine.str();
		energyLines.push_back(line);
		i++;
	}
	ofstream eout;
	string eoutfile = _opt.outputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	eout << "Sequence" << t << "RotamerValues" << t << "Interface" << t << "InterfaceSequence" << t << "Seed" << t;
	eout << enerTerms.str() << endl;
	_sout << "Sequence" << t << "RotamerValues" << t << "Interface" << t << "InterfaceSequence" << t << "Seed" << t;
	_sout << enerTerms.str() << endl;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
		_sout << energyLines[i] << endl;
	}
	eout.close();
}

void outputTime(auto _start, string _descriptor, ofstream &_sout){
	auto end = chrono::system_clock::now();
	time_t endTimeFormatted = chrono::system_clock::to_time_t(end); 
	chrono::duration<double> elapsedTime = end-_start;
	_sout.precision(3);// rounds output to 3 decimal places
	cout.precision(3);// rounds output to 3 decimal places
	_sout << _descriptor << " Time: " << ctime(&endTimeFormatted);
	_sout << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
	cout << _descriptor << " Time: " << ctime(&endTimeFormatted);
	cout << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
}

// help functions
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputErrorMessage(Options &_opt){
	cout << endl;
	cout << "The program terminated with errors:" << endl;
	cout << endl;
	cout << _opt.errorMessages << endl;
	cout << endl;
	cout << _opt.OPerrors << endl;
	usage();
}

void help(Options defaults) {
	cout << "This program runs as:" << endl;
	cout << " % seqDesign " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --numberOfStructuresToMCRepack <int> --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --MCStartTemp <double> --MCEndTemp <double> --MCCurve <CONSTANT-0, LINEAR-1, EXPONENTIAL-2, SIGMOIDAL-3, SOFT-4>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --thread <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   --weight_hbond <double> --weight_vdw <double> --weight_solv <double> --weight_seqEntropy <double>" << endl;
	cout << "   --sasaRepackLevel <rotLevel> (in format SL95.00; 4 levels used by default) --interfaceLevel <int> " << endl << endl;
	cout << "Template Configuration file (copy and paste the below into a file.config and run code as bin/seqDesign --config file.config" << endl;
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "outputDir " << defaults.outputDir << endl;

	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "geometryDensityFile " << defaults.geometryDensityFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;
	cout << setw(20) << "backboneFile " << defaults.backboneFile << endl;
	cout << setw(20) << "selfEnergyFile " << defaults.selfEnergyFile << endl;
	cout << setw(20) << "pairEnergyFile " << defaults.pairEnergyFile << endl;

	cout << "#Geometry and Transformation parameters" << endl;
	cout << setw(20) << "xShift" << defaults.xShift << endl;
	cout << setw(20) << "crossingAngle" << defaults.crossingAngle << endl;
	cout << setw(20) << "axialRotation" << defaults.axialRotation << endl;
	cout << setw(20) << "zShift" << defaults.zShift << endl;
	cout << setw(20) << "thread" << defaults.thread << endl;
	cout << setw(20) << "backboneLength " << defaults.backboneLength << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;
	cout << setw(20) << "deleteTerminalBonds" << defaults.deleteTerminalBonds << endl;
	cout << setw(20) << "useSasaBurial" << defaults.useSasaBurial << endl;
	cout << setw(20) << "getGeoFromPDBData" << false << endl;//Since we already have the geometry output here, default to false in the rerun config
	cout << setw(20) << "runDEESingles" << defaults.runDEESingles << endl;
	cout << setw(20) << "runDEEPairs" << defaults.runDEEPairs << endl;
	cout << setw(20) << "runSCMF" << defaults.runSCMF << endl;
	//TODO: set this up so that instead of running through that, I just have the output sequence and state below
	if (defaults.runDEESingles == false && defaults.runDEEPairs == false && defaults.runSCMF == false){

	}

	if (defaults.useSasaBurial == true){
		//cout << "#Load Rotamers based on SASA scores" << endl;
		for (uint i=0; i<defaults.sasaRepackLevel.size()-1; i++){
			cout << setw(20) << "sasaRepackLevel" << defaults.sasaRepackLevel[i] << endl;
		}
		cout << setw(20) << "interfaceLevel" << defaults.interfaceLevel << endl;
	} else {
		//cout << "#Load Rotamers by interface and non interfacial positions" << endl;
		cout << setw(20) << "SL" << defaults.SL << endl;
	}

	cout << "#MonteCarlo Paramenters" << endl;
	cout << setw(20) << "MCCycles" << defaults.MCCycles << endl;
	cout << setw(20) << "MCMaxRejects" << defaults.MCMaxRejects << endl;
	cout << setw(20) << "MCStartTemp" << defaults.MCStartTemp << endl;
	cout << setw(20) << "MCEndTemp" << defaults.MCEndTemp << endl;
	cout << setw(20) << "MCCurve" << defaults.MCCurve << endl;
	cout << setw(20) << "greedyCycles" << defaults.greedyCycles << endl;

	cout << "#Alternate IDs" << endl;
	for (uint i=0; i<defaults.Ids.size()-1; i++){
		cout << setw(20) << defaults.Ids[i] << endl;
	}

	cout << endl << "#Rerun Seed" << endl;
	cout << setw(20) << "seed" << defaults.seed << endl;

	cout << endl << "#Energy term weights" << endl;
	cout << setw(20) << "weight_vdw " << defaults.weight_vdw << endl;
	cout << setw(20) << "weight_hbond " << defaults.weight_hbond << endl;
	cout << setw(20) << "weight_solv " << defaults.weight_solv << endl;
	cout << setw(20) << "weight_seqEntropy " << defaults.weight_seqEntropy << endl;
	cout << endl;
}

void computeMonomerEnergy(System &_sys, System &_helicalAxis, Options &_opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq,
 RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {
	// output the start time
	outputTime(clockTime, "Compute Monomer Energy Start", _sout);

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);
	CSBMono.setBuildTerm("CHARMM_IMM1REF", true);
	CSBMono.setBuildTerm("CHARMM_IMM1", true);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(30);
	
	CSBMono.updateNonBonded(10,12,50);
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsInactive();
	monoEset->setTermActive("CHARMM_IMM1REF", true);
	monoEset->setTermActive("CHARMM_IMM1", true);
	monoEset->setTermActive("CHARMM_VDW", true);
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all bonding near the termini of our helices for a list of interactions
    deleteTerminalBondInteractions(_sys,_opt.deleteTerminalInteractions);

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.updateWeights();
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	_trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), _trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	_trans.translate(chainA, interDistVect);

	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	_trans.translate(chainA, move5Down);
	double bestZ = -5.0;

	monoSys.calcEnergy();

	// Repack side chains
	monoSpm.setOnTheFly(1);
	monoSpm.calculateEnergies();
        monoSpm.runGreedyOptimizer(_opt.greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	_helicalAxis.saveAltCoor("BestMonomerAxis");
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		//double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			//_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	_helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		_helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		_trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		_trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		for(int j=0; j<=3; j++) { // test at 4 rotations 0, 90, 180 and 270 degrees
			//==================================
			//====== Axial Rot ======
			//==================================
			monoAxialRotation = j * 90.0;

			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			//_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			//monoSys.writePdb("mono_" + MslTools::doubleToString(monoTilt) + "_" + MslTools::doubleToString(monoAxialRotation) + ".pdb");

			if(currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				bestTilt = monoTilt;
				bestRotation = monoAxialRotation;
				monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
				monoSys.saveAltCoor("savedBestMonomer");
				_helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			_trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

		}
	}

	//MonteCarloManager MCMngr(1000.0, 0.5, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		_helicalAxis.applySavedCoor("BestMonomerAxis");

		int moveToPreform = _RNG.getRandomInt(2);

		double deltaZShift = 0.0;
		double deltaTilt = 0.0;
		double deltaAxialRotation = 0.0;

		//======================================
		//====== Z Shift ======
		//======================================
		if (moveToPreform == 0) {
			deltaZShift = getStandardNormal(_RNG) * 1.0;
			CartesianPoint translateA = axisA(1).getCoor() - axisA(0).getCoor(); // vector minus helical center
			translateA = translateA.getUnit() * deltaZShift; // unit vector of helical _axis times the amount to shift by
			_trans.translate(chainA, translateA);
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			_trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			_trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			_trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_fout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_fout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_fout << "state rejected   energy: " << currentEnergy << endl;
		} else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;
		}
		counter++;
	}

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	//Calculate Monomer energy for output
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getMinBound()[0]*2;

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	monoSasa.calcSasa();
	double monomerSasa = monoSasa.getTotalSasa();
	double totalMonomerSasa = monomerSasa*2;

	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	_sequenceEnergyMap[_seq]["MonomerSasa"] = totalMonomerSasa;
	double dimerEnergy = _sequenceEnergyMap[_seq]["Dimer"];
	cout << "Dimer Energy: " << _seq << ": " << dimerEnergy << endl;
	double totalEnergy = dimerEnergy-monomerEnergy;
	_sout << "-Dimer - Monomer = " << dimerEnergy << " - " << monomerEnergy << " = " << totalEnergy << endl;
	cout << "-Dimer - Monomer = " << dimerEnergy << " - " << monomerEnergy << " = " << totalEnergy << endl;
	_sequenceEnergyMap[_seq]["preRepackTotal"] = totalEnergy;

	// cal;ulate the energy of the monomer for positions 4-18
	if (_opt.useAlaAtTermini){
		vector<double> selfVec = calcBaselineEnergies(monoSys, _opt.thread, _opt.backboneLength);
		vector<double> pairVec = calcPairBaselineEnergies(monoSys, _opt.thread, _opt.backboneLength);
		double self = sumDoubleVector(selfVec);
		double pair = sumDoubleVector(pairVec);
		_sequenceEnergyMap[_seq]["MonomerWithoutAlaEnds"] = 2*(self+pair);
	}
	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	_helicalAxis.clearSavedCoor("BestMonomerAxis");
	_helicalAxis.clearSavedCoor("bestZ");

	// output end time
	outputTime(clockTime, "Compute Monomer Energy End", _sout);
}

void runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _rotamerSamplingPositionVector,
 vector<uint> _interfacePositions, map<string, map<string,double>> &_sequenceEnergyMap, map<string, double> _sequenceEntropyMap,
 ofstream &_out){
	// output the starting time
	outputTime(clockTime, "Self Consistent Mean Field Start", _out);

	// link the interfacial positions (for quicker calculation of initial sequence for homodimers)
	if (_opt.linkInterfacialPositions){
		vector<vector<string>> linkedPos = convertToLinkedFormat(_sys, _interfacePositions, _opt.backboneLength);
		_sys.setLinkedPositions(linkedPos);
	}
	
	// SelfPairManager setup
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&_sys);
	spm.setVerbose(false);
	spm.setOnTheFly(true);
	spm.setMCOptions(1000, 0.5, 5000, 3, 10, 1000, 0.01);//changed to sigmoid and added up to 5000
	spm.saveEnergiesByTerm(true); //added back in on 09_21_2021 to get the vdw and hbond energies
	spm.calculateEnergies();

	//Setup running SCMF
	cout << "Running Self Consistent Mean Field" << endl;
	_out << "Running Self Consistent Mean Field" << endl;
	spm.setRunSCMF(true);
	spm.setRunSCMFBiasedMC(true);
	spm.setRunUnbiasedMC(false);

	// run and find a sequence using the chosen parameters (MCOptions, SCMF, DEE, etc.)
	spm.runOptimizer();

	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
	_sys.setActiveRotamers(bestState);
	_sys.saveAltCoor("SCMFState");
	string startSequence = convertPolymerSeqToOneLetterSeq(_sys.getChain("A")); //used for outputting starting sequence
	string rotamerSamplingString = convertVectorUintToString(_rotamerSamplingPositionVector); // string of rotamer sampling number (if 4 rotamer levels, 0-3 for each position)
	string interfaceSeq = getInterfaceSequence(_opt.interfaceLevel, rotamerSamplingString, startSequence);

	// output spm run optimizer information
	spmRunOptimizerOutput(spm, _sys, interfaceSeq, _out);
	
	//Add energies for initial sequences into the sequenceEnergyMap
	pair<string,vector<uint>> startSequenceStatePair = make_pair(startSequence, bestState);
	getEnergiesForStartingSequence(_opt, spm, startSequence, bestState, _interfacePositions, _sequenceEnergyMap, _sequenceEntropyMap);
	getSasaForStartingSequence(_sys, startSequence, bestState, _sequenceEnergyMap);
	
	// reset the energy set
	resetEnergySet(_sys, _opt.energyTermList);
	vector<uint> unlinkedState = unlinkBestState(bestState, _interfacePositions, _opt.backboneLength);
	
	// output the ending time
	outputTime(clockTime, "Self Consistent Mean Field End", _out);
}

void outputStartDesignInfo(Options &_opt){
	// output the starting geometry	
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << _opt.xShift << endl;
	cout << "crossingAngle: " << _opt.crossingAngle << endl;
	cout << "axialRotation: " << _opt.axialRotation << endl;
	cout << "zShift:        " << _opt.zShift << endl << endl;

	// output the alternate ids to be used for design	
	string alternateIds = getAlternateIdString(_opt.Ids);
	cout << "Variable amino acids at the interface: " << alternateIds << endl << endl;
}