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

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "bbOptimizer";
string programDescription = "Reads a PDB, does a backbone repack, and calculates the energy. Can also take an input sequence using input pdb as a backbone";
// - code that can either take in a pdb or set a homodimer sequence to a given geometry and does a 
// local backbone optimization to get the most stable energy for the sequence
// - Can input a sequence that is different from the given pdb and it should replace the sequence
// - I use this code after sequence design to locally optimize the geometry (also done during the
// design process, but this does a longer locally optimized repack)
// - To ensure a local backbone optimization, it decreases the potential geometric shifts 
// (ie. if xShift starts allowing 0.5 changes, after the first time it makes that change, it decreases
//  the potential next xShift towards a minimum using a monte carlo step decrease)
string programAuthor = "Gilbert Loiseau";
string programVersion = "3";
string programDate = "23 Feb 2024";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;
auto clockTime = chrono::system_clock::now();

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/
struct Options{
	// input files
	string backboneCrd; //initial coordinates for helix backbones: crd file
	string backboneFile; //initial coordinates for helix backbones: pdb file
	string outputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string pdbFile; //input pdb file
	string helicalAxis; //file with helical axis information

	string sequence; //sequence of the protein
	int thread; //thread number on the gly69 helix: depending on the sequence length, proper thread number can be used to center the helix in the membrane
	string interface; //interface to design on; if empty, use SASA

	// Monte Carlo parameters
    int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	// booleans
	bool verbose; //TRUE: write energy outputs and calculations throughout the run to the terminal OR FALSE: only write outputs to output files
	bool deleteTerminalBonds; //TRUE: delete hydrogen bonds at the termini of sequences to not be considered in hydrogen bonding score OR FALSE: keep hydrogen bonds at termini
	bool useSasaBurial; //TRUE: use solvent accessible surface area to designated the number of rotamers at each position on the dimer OR FALSE: input set number of rotamers for both interface and non-interface

	// use different energy parameters
	bool useIMM1;
	bool useElec;

	// repack parameters
	int greedyCycles;

	// load rotamers useSasa = false
	string SL; //number of rotamers
	std::vector<string> sasaRepackLevel; //vector of levels
	int interfaceLevel; // level for the interface

    // energy weights
	double weight_vdw; //weight of vdw energy contribution to total energy: default = 1
	double weight_IMM1; //weight of solvation energy contribution to total energy: default = 1
	double weight_IMM1REF; //weight of solvation energy contribution to total energy: default = 1
	double weight_elec; //weight of electrostatic energy contribution to total energy: default = 1
	double weight_hbond; //weight of hbond energy contribution to total energy: default = 1

	// Starting Geometry
	double xShift; //distance between helices
	double crossingAngle; //crossing angle between helices
	double axialRotation; //rotation of helices
	double zShift; //position of the crossing point between helices
	double density; //density of axialRotation and zShift from clashing tests at different geometries
	bool negAngle; //designates if the crossing angle is negative or positive
	bool negRot; //designates if the axial rotation is negative or positive
	
	// Backbone Monte Carlo parameters
	int backboneMCCycles;
	int backboneMCMaxRejects;
	double backboneMCStartTemp;
	double backboneMCEndTemp;
	int backboneMCCurve;
	int backboneConvergedSteps;
	double backboneConvergedE;

	// backbone repack variables
	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;
	double deltaXLimit; 
	double deltaCrossLimit;
	double deltaAxLimit;
	double deltaZLimit;
	bool decreaseMoveSize;
	int backboneRepackCycles;

	// alternate identities
	vector<string> alternateIds; //alternate AA identities for interfacial positions
    
    // terminal interactions
    vector<string> deleteTerminalInteractions; //terminal interactions to be considered in the energy calculation
	vector<string> energyTermList; // list of energy terms to evaluate and save (not yet fully implemented for usage, only for saving)

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPErrors; //the errors from the option parser

	// run parameters
	string rerunConf; // data for a configuration file that would rerun the job as the current run
	string configfile;
	int seed;
};

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === BACKBONE OPTIMIZATION FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void backboneOptimizer(Options &_opt, RandomNumberGenerator &_RNG, string _sequence, PolymerSequence _PS, AtomPointerVector &_startPdb, map<string,double> _startGeometry,
 map<string, map<string, double>> &_sequenceEnergyMapFinal, ofstream &_eout);
void backboneOptimizeMonteCarlo(Options &_opt, System &_sys, SelfPairManager &_spm, map<string,double> &_sequenceEnergyMap,
 string _sequence, vector<uint> _bestState, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 string _outputDir, int _rep, ofstream &_sout);

// decrease move size functions
void getCurrentMoveSizes(double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 double _deltaXLimit, double _deltaCrossLimit, double _deltaAxLimit, double _deltaZLimit, bool &_decreaseMoveSize, int _moveToPerform);
double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);

// axial rotation and zShift conversion from given input ax' and z' to ax and z (Mueller 2014; Fig. S1)
void convertToRelativeAxAndZ(double _axialRotation, double _zShift, double &_relativeAx, double &_relativeZ);

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === FUNCTIONS FROM CATM ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
// geometry functions
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
double getStandardNormal(RandomNumberGenerator& RNG);

// compute the monomer energy  
// REDACTED VERSION: void computeMonomerEnergy(System &_sys, System &_helicalAxis, Options &_opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq,
// RandomNumberGenerator &_RNG, ofstream &_sout);
// NEW VERSION: transforms object is made in the function
void computeMonomerEnergy(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, ofstream &_sout);

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === HELPER FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
// function that adds the energy terms to the system, sets up the energy map, and transforms the system to the given geometry
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);

// delete terminal interactions for a given list of interactions (CHARMM_IMM1, CHARMM_IMM1REF, CHARMM_VDW, SCWRL4_HBOND, etc.)
void deleteTerminalBondInteractions(System &_sys, vector<string> &_deleteTerminalInteractionList);

// set the active identity for each position to the identity in the given sequence (only for homodimers)
string extractSequence(System &_sys);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB,
 map<string,double> _geometry, Transforms &_trans);

// side chain functions
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);

// usage function
void usage();

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === OUTPUT FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
map<string, double> getGeometryMap(map<string,double> _geometry, string _descriptor);
void addEnergiesToMap(Options &_opt, SelfPairManager &_spm, map<string, double> &_energyMap, vector<uint> _bestState, string _descriptor);
void addGeometryToEnergyMap(map<string, double> _geometryMap, map<string, double> &_energyMap);
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1);
void outputFiles(Options &_opt, double _seed, map<string,map<string,double>> _sequenceEnergyMap, ofstream &_sout);
void outputTime(auto _start, string _descriptor, ofstream &_sout);
string setupOutputDirectory(string _dirName);
void outputErrorMessage(string _message, string _optionParserErrors);

// function for writing a pdb
void writePdb(System &_sys, string _outputDir, string _pdbName);

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === OPTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
Options parseOptions(int _argc, char * _argv[]);

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
	//Add in some default options that can easily be changed here
	Options opt = parseOptions(argc, argv);

	if (opt.errorFlag) {
		outputErrorMessage(opt.errorMessages, opt.OPErrors);
		exit(1);
	}
	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	// summary file output
	ofstream sout;
	// error file output
	ofstream err;
	// rerun config output
	ofstream rerun;

	string outputDir = setupOutputDirectory(opt.outputDir);

	string soutfile = outputDir + "/summary.out";
	string errfile  = outputDir + "/errors.out";
	string rerunfile = outputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	rerun << opt.rerunConf << endl;
	rerun.close();
	cout << 1 << endl;	
	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(opt.pdbFile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	AtomPointerVector & apv = pdb.getAtomPointers();
	cout << 2 << endl;	

	Chain & chainA = pdb.getChain("A");
	Chain & chainB = pdb.getChain("B");
	cout << 3 << endl;	

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// write the pdb for the final state of the backbone optimization
	writePdb(pdb, outputDir, "initialPdb"); // Issue with the PDBs from the design program; fixed now, but initialPdb is still useful for debugging if the problem occurs again
	cout << 4 << endl;	

	string originalSeq = extractSequence(pdb);
	// replace the first and last positions with the patch residues
	// mutate the sequence here; should I add in a check to make sure that the sequence is calculate normally too?
	// TODO: add in option useAlaAtTermini to switch AAs to Ala at termini
	string seq;
	if (opt.sequence == ""){
		seq = originalSeq;
		seq.replace(seq.begin(), seq.begin()+3, "AAA");
		seq.replace(seq.end()-3, seq.end(), "AAA");
	} else {
		seq = opt.sequence;
	}

	/******************************************************************************
	 *                   === INITIALIZE POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polySeq = convertToPolymerSequenceNeutralPatch(seq, opt.thread);
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	
	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	// AtomPointerVector for the helical axis for each chain
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Setup random number generator
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 

	// initialize the sequenceEnergyMap
	map<string, map<string, double>> sequenceEnergyMap;

	// add in the backbone geometry repack for the sequence (adapted from my sequence design code)
	map<string,double> geometry;
	// add the backbone geometry to the geometry map
	geometry["xShift"] = opt.xShift;
	geometry["crossingAngle"] = opt.crossingAngle;
	geometry["axialRotation"] = opt.axialRotation;
	geometry["zShift"] = opt.zShift;
	
	backboneOptimizer(opt, RNG, seq, PS, apv, geometry, sequenceEnergyMap, sout);

	// write out the summary file
	double seed = RNG.getSeed();
	outputFiles(opt, seed, sequenceEnergyMap, sout);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

    // close all of the output file writers
    sout.close();
    err.close();
	exit(0);
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === BACKBONE OPTIMIZATION FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void backboneOptimizer(Options &_opt, RandomNumberGenerator &_RNG, string _sequence, PolymerSequence _PS, AtomPointerVector &_startPdb, map<string,double> _startGeometry,
 map<string, map<string, double>> &_sequenceEnergyMapFinal, ofstream &_eout){
	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	// AtomPointerVector for the helical axis for each chain
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// get the starting geometry using polyglycine
	System startGeom;
	setGly69ToStartingGeometry(_opt,startGeom,helicalAxis,axisA,axisB,_startGeometry,trans);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	prepareSystem(_opt, sys, startGeom, _PS);

	sys.assignCoordinates(_startPdb,false);
	sys.buildAllAtoms();
	
	double dimer = sys.calcEnergy();
    cout << "Dimer Energy: " << dimer << endl;
	sys.printEnergySummary();
	// initialize the object for loading rotamers into system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	loadRotamers(sys, sysRot, _opt.SL);
	sys.buildAllAtoms();
	
	// Optimize Initial Starting Position
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	//repackSideChains(spm, _opt.greedyCycles);
	repackSideChains(spm, 50);
	vector<uint> stateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(stateVec);
	
	dimer = sys.calcEnergy();
    cout << "Dimer Energy: " << dimer << endl;
    sys.printEnergySummary();
	writePdb(sys, _opt.outputDir, "startPdb");
	/******************************************************************************
	 *   === MONTE CARLO TO SEARCH FOR BEST SEQUENCES AND BACKBONE OPTIMIZE ===
	 ******************************************************************************/
	// define the best state
	vector<uint> bestState = stateVec;
	
	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// set the starting geometry
	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingState");
	double currentEnergy = sys.calcEnergy();

	string repackDir = _opt.outputDir;
	//computeMonomerEnergy(sys, helicalAxis, _opt, trans, _sequenceEnergyMapFinal, _sequence, _RNG, _eout);
	computeMonomerEnergy(sys, _opt, _RNG,  _sequenceEnergyMapFinal, _sequence, _eout);

	// get the monomer energy from the sequence energy map
	double monomerEnergy = _sequenceEnergyMapFinal[_sequence]["Monomer"];
	map<string,double> monomerEnergyTerms;
	// loop through the terms in the sequence energy map
	for (auto & term : _sequenceEnergyMapFinal[_sequence]){
		// check if the term name contains "Monomer"
		if (term.first.find("Monomer") != string::npos){
			// if it does, add it to the monomer energy terms map
			monomerEnergyTerms[term.first] = term.second;
		}
	}
	_sequenceEnergyMapFinal.clear();
	// loop through the number of backbone repack cycles
	for (uint i=0; i<_opt.backboneRepackCycles; i++){
		// output the start time
		outputTime(clockTime, "Backbone optimize replicate " + to_string(i) + " Start", _eout);

		// initialize the sequence energy map
		map<string, double> energyMap;

		// add the starting geometry to the energy map
		map<string, double> startGeometry = getGeometryMap(_startGeometry, "preOptimize");
		addGeometryToEnergyMap(startGeometry, energyMap);

		// save the pre backbone optimize energy and the replicate number in the 
		energyMap["preOptimizeEnergy"] = currentEnergy;
		energyMap["replicateNumber"] = i;
		
		// backbone optimizer function
		backboneOptimizeMonteCarlo(_opt, sys, spm, energyMap, _sequence, bestState, helicalAxis, axisA, axisB, apvChainA, apvChainB,
		 trans, _RNG, monomerEnergy, repackDir, i, _eout);
		
		// append the sequence energy map to the final sequence energy map by looping through the sequence energy map
		for (auto &energy : energyMap){
			_sequenceEnergyMapFinal[_sequence+"_"+to_string(i)][energy.first] = energy.second;
		}
		// add the monomer energy terms to the sequence energy map
		for (auto &term : monomerEnergyTerms){
			_sequenceEnergyMapFinal[_sequence+"_"+to_string(i)][term.first] = term.second;
		}

		// reset the energy map
		energyMap.clear(); // reset the energyMap for the next cycle

		// reset the system to the starting state
		sys.applySavedCoor("startingState");
		helicalAxis.applySavedCoor("startingState");
		
		// output the end time
		outputTime(clockTime, "Backbone optimize replicate " + to_string(i) + " End", _eout);
	}
}

// monte carlo backbone optimization function
void backboneOptimizeMonteCarlo(Options &_opt, System &_sys, SelfPairManager &_spm, map<string,double> &_energyMap,
 string _sequence, vector<uint> _bestState, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 string _outputDir, int _rep, ofstream &_sout){
	// Setup backbone repack file
	ofstream bbout;
	string bboutfile = _outputDir + "/" + to_string(_rep) + "_repack.out";
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
	// calculate starting sasa
	SasaCalculator startSasa(_sys.getAtomPointers());
	startSasa.calcSasa();
	double sasa = startSasa.getTotalSasa();
	_energyMap["PreBBOptimizeSasa"] = sasa;

	// Monte Carlo Repack Manager Setup
	MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects, _opt.backboneConvergedSteps, _opt.backboneConvergedE);

	vector<unsigned int> MCOBest = _bestState;
	unsigned int counter = 0;
	_sys.setActiveRotamers(_bestState);
	double currentEnergy = _spm.getStateEnergy(_bestState)-_monomerEnergy;
	double dimer = _spm.getStateEnergy(_bestState);
	double calcDimer = _sys.calcEnergy();
	_energyMap["TotalPreOptimize"] = currentEnergy;
	addEnergiesToMap(_opt, _spm, _energyMap, _bestState, "PreOptimize");
	// 2024-2-24: change to the below to split the IMM1 and IMM1REF energies
	//_energyMap["VDWDimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_VDW");
	//_energyMap["IMM1DimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1")+_spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	//_energyMap["CHARMM_IMM1_DimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1");
	//_energyMap["CHARMM_IMM1REF_DimerPreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	//_energyMap["HBONDDimerPreOptimize"] = _spm.getStateEnergy(_bestState, "SCWRL4_HBOND");
	//if (_opt.useElec == true){
	//	_energyMap["CHARMM_ELEC_PreOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_ELEC");
	//}
	
	double bestEnergy = currentEnergy;
	double prevBestEnergy = currentEnergy;
	MCMngr.setEner(currentEnergy);

	// setup variables for shifts: ensures that they start from the proper values for every repack and not just the final value from the initial repack
	bool decreaseMoveSize = _opt.decreaseMoveSize;
	double deltaX = _opt.deltaX;
	double deltaCross = _opt.deltaCross;
	double deltaAx = _opt.deltaAx;
	double deltaZ = _opt.deltaZ;

	_sys.saveAltCoor("savedRepackState");
	_helicalAxis.saveAltCoor("BestRepack");
	// uncomment the below and lines ... ; it makes the outputs quite large, but you'll have outputs for the entire optimization
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
		if (counter == 0){
			_energyMap["firstRepackEnergy"] = currentEnergy;
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
				getCurrentMoveSizes(startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, _opt.deltaXLimit,
				 _opt.deltaCrossLimit, _opt.deltaAxLimit, _opt.deltaZLimit, decreaseMoveSize, moveToPreform);
			}
			bbout << "MCAccept " << counter <<  " xShift: " << finalXShift << " crossingAngle: " << finalCrossingAngle << " axialRotation: " << finalAxialRotation << " zShift: " << finalZShift << " energy: " << currentEnergy << endl;
			counter++;
			//writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	//writer.close();
	// TODO: add in here to run a final optimization where the parameters decrease to 0
	bbout << "End Repack Cycles" << endl << endl; 
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_sys.applySavedCoor("savedRepackState");
	_helicalAxis.applySavedCoor("BestRepack");
	// write the pdb for the final state of the backbone optimization
	writePdb(_sys, _outputDir, to_string(_rep));
	double dimerEnergy = _spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-_monomerEnergy;
	
	// Output change in geometry
	bbout << "***REPACK GEOMETRY***" << endl;
	bbout << "Energy;        Before: " << prevBestEnergy << "; After: " << bestEnergy << endl << endl;

	// sets the updated backbone parameters
	map<string,double> geometry;
	geometry["xShift"] = finalXShift;
	geometry["crossingAngle"] = finalCrossingAngle;
	geometry["axialRotation"] = finalAxialRotation;
	geometry["zShift"] = finalZShift;

	// add the geometry to the map post backbone optimization
	map<string, double> endGeometry = getGeometryMap(geometry, "end");
	addGeometryToEnergyMap(endGeometry, _energyMap);

	SasaCalculator endDimerSasa(_sys.getAtomPointers());
	endDimerSasa.calcSasa();
	double endSasa = endDimerSasa.getTotalSasa();
	_energyMap["OptimizeSasa"] = endSasa;
	_energyMap["Total"] = finalEnergy;

	addEnergiesToMap(_opt, _spm, _energyMap, _bestState, "Optimize");
	// 2024-2-24: change to the below to split the IMM1 and IMM1REF energies
	//_energyMap["VDWDimerOptimize"] = _spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	//_energyMap["IMM1DimerOptimize"] = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+_spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	//_energyMap["CHARMM_IMM1_DimerOptimize"] = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1");
	//_energyMap["CHARMM_IMM1REF_DimerOptimize"] = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	//_energyMap["HBONDDimerOptimize"] = _spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	
	bbout << MCMngr.getReasonCompleted() << endl;	
	bbout << "Monte Carlo repack complete. Time: " << diffTimeMC/60 << "min" << endl << endl;

	// clear the saved repack state
	_sys.clearSavedCoor("savedRepackState");
	_helicalAxis.clearSavedCoor("BestRepack");
}

// adjust the move sizes based on the monte carlo function (decrease by multiplying by the change in temperature)
void getCurrentMoveSizes(double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 double _deltaXLimit, double _deltaCrossLimit, double _deltaAxLimit, double _deltaZLimit, bool &_decreaseMoveSize, int _moveToPerform) {
	// define the temperature change
	double decreaseMultiplier = _endTemp/_currTemp;
	// setup the decrease booleans; if the move size has reached its minimum, don't decrease it further and set to false
	bool decreaseZ = false;
	bool decreaseAx = false;
	bool decreaseCross = false;
	bool decreaseX = false;

	// move to perform is the move that was just performed from the backboneOptimizer function
	if (_moveToPerform == 0){
		decreaseZ = true;
		_deltaZ = decreaseMoveSize(_deltaZ, _deltaZLimit, decreaseMultiplier, decreaseZ);
	} else if (_moveToPerform == 1) {
		decreaseAx = true;
		_deltaAx = decreaseMoveSize(_deltaAx, _deltaAxLimit, decreaseMultiplier, decreaseAx);
	} else if (_moveToPerform == 2) {
		decreaseCross = true;
		_deltaCross = decreaseMoveSize(_deltaCross, _deltaCrossLimit, decreaseMultiplier, decreaseCross);
	} else if (_moveToPerform == 3) {
		decreaseX = true;
		_deltaX = decreaseMoveSize(_deltaX, _deltaXLimit, decreaseMultiplier, decreaseX);
	}

	// if all of the move sizes have reached their minimum, set the decrease move size boolean to false	
	if (decreaseX == false && decreaseCross == false && decreaseAx == false && decreaseZ == false){
		_decreaseMoveSize = false;
	}
}

// function that decreases the acceptable move size based on the monte carlo function temperature change (decreaseMultiplier)
// only decreases the move size if the move size is greater than the move limit
// as the temperature decreases, the acceptable move size decreases

double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease) {
	// get the change between the current acceptable move size and the move limit (ex: xShift 0.5-0.1)
	double diffMoveSize = _moveSize-_moveLimit;
	// get the current move size decrease based on the change in temperature (ex: xShift 0.4-(0.4*995/1000) = 0.4-0.398 = 0.002)
	double moveDecrease = diffMoveSize-(diffMoveSize*_decreaseMultiplier);// edited on 2022-10-6: it now works properly; before it decreased to the move limit
	// get the adjusted move size (the change in acceptable move size minus the move decrease; ex: xShift 0.5-0.002 = 0.498 acceptable)
	double newMoveSize = _moveSize - moveDecrease;
	// check to see if the new move size is greater than the accepted move limit (ex: xShift 0.498 > 0.1; 0.498 is the new maximum xShift change)
	if (newMoveSize > _moveLimit){
		return newMoveSize;
	} else {
		// if the move size is less, set the move size to the move limit and set the decrease boolean to false 
		// to stop accepting any decreases in that move
		_decrease = false;
		return _moveSize;
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



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === OUTPUT FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
map<string, double> getGeometryMap(map<string,double> _geometry, string _descriptor){
	// get geometries from the initial geometry map
	double xShift = _geometry["xShift"];
	double crossingAngle = _geometry["crossingAngle"];
	double axialRotation = _geometry["axialRotation"];
	double zShift = _geometry["zShift"];

	// add the descriptor to the output geometry map
	map<string, double> geometryMap;
	geometryMap[_descriptor+"_xShift"] = xShift;
	geometryMap[_descriptor+"_crossingAngle"] = crossingAngle;
	geometryMap[_descriptor+"_axialRotation"] = axialRotation;
	geometryMap[_descriptor+"_zShift"] = zShift;

	// convert axialRotation and zShift for output (from parallelogram from Ben)
	double relativeAx;
	double relativeZ;
	convertToRelativeAxAndZ(axialRotation, zShift, relativeAx, relativeZ);
	geometryMap[_descriptor+"_axialRotationPrime"] = relativeAx;
	geometryMap[_descriptor+"_zShiftPrime"] = relativeZ;
	return geometryMap;
}

// 2024-2-24: change energy outputs to the below to split the IMM1 and IMM1REF energies, and make it able to accomodate any additional energy terms
void addEnergiesToMap(Options &_opt, SelfPairManager &_spm, map<string, double> &_energyMap, vector<uint> _bestState, string _descriptor){
	// loop through the energy terms in the options and get the energy for each term
	for (auto &it : _opt.energyTermList){
		string energyLabel = it.substr(7,it.length()) + '_' + _descriptor;//Removes the CHARMM_ and SCWRL4_ before energyTerm names and adds descriptor to the end
		_energyMap[energyLabel] = _spm.getStateEnergy(_bestState, it);
	}
}

void addGeometryToEnergyMap(map<string, double> _geometryMap, map<string, double> &_energyMap){
	for (auto &it : _geometryMap){
		_energyMap[it.first] = it.second;
	}
}

// output energies by term into a referenced energyMap
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1){
	if (_includeIMM1 == false){//No IMM1 Energy (for the Monte Carlos, both dimer and monomer)
		for (uint i=0; i<_energyTermList.size(); i++){
			string energyTerm = _energyTermList[i]; //CHARMM_ and SCWRL4_ terms
			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;//Removes the CHARMM_ and SCWRL4_ before energyTerm names
			if (energyTerm.find("IMM1") != string::npos){
				continue;
			} else {
				if (_energyDescriptor.find("Monomer") != string::npos){
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
				}
			}
		}
		if (_energyDescriptor.find("Monomer") != string::npos){
			//skip if monomer; could add calc baseline here at some point
		} else {
			_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
			_energyMap["DimerSelfBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE");
			_energyMap["DimerPairBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
		}
	} else if (_includeIMM1 == true){//IMM1 Energies
		for (uint i=0; i<_energyTermList.size(); i++){
			string energyTerm = _energyTermList[i];
			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;
			// 2024-2-24: change to the below to split the IMM1 and IMM1REF energies
			//if (_energyDescriptor.find("Monomer") != string::npos){
			//	if (energyTerm.find("IMM1") != string::npos){
			//		_energyMap["IMM1Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF"))*2;
			//	} else {
			//		_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
			//	}
			//} else {
			//	if (energyTerm.find("IMM1") != string::npos){
			//		_energyMap["IMM1Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF");
			//	} else {
			//		_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
			//}
			if (_energyDescriptor.find("Monomer") != string::npos){
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1_Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1"))*2;
					_energyMap["IMM1REF_Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF"))*2;
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
				}
			} else {
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1_Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1");
					_energyMap["IMM1REF_Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF");
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
				}
			}
		}
		//_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
	}
}

void outputFiles(Options &_opt, double _seed, map<string,map<string,double>> _sequenceEnergyMap, ofstream &_sout){
	// Setup vector to hold energy file lines
	vector<string> energyLines;
	// get the run parameters
	string t = ",";
	stringstream enerTerms;
	// For loop to setup the energy file
	uint i = 0;
	for (auto &seq : _sequenceEnergyMap){
		stringstream seqLine;
		string geometry = seq.first;
		// get the interface sequence
		seqLine << geometry << t << _opt.interface << t << _seed << t;
		map<string,double> energyMap = _sequenceEnergyMap[geometry];
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
	eout << "Geometry" << t << "Interface" << t << "Seed" << t;
	eout << enerTerms.str() << endl;
	_sout << "Geometry" << t << "Interface" << t << "Seed" << t;
	_sout << enerTerms.str() << endl;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
		_sout << energyLines[i] << endl;
	}
	eout.close();
}

// setups up the output directory in the current working directory, good for running on chtc and biochem servers
string setupOutputDirectory(string _dirName){
	_dirName = string(get_current_dir_name()) + "/" + _dirName;
	string cmd = "mkdir -p " + _dirName;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	return _dirName;
}

void outputErrorMessage(string _message, string _optionParserErrors){
	cout << endl;
	cout << "The program terminated with errors:" << endl;
	cout << endl;
	cout << _message << endl;
	cout << endl;
	cout << _optionParserErrors << endl;
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

void writePdb(System &_sys, string _outputDir, string _pdbName){
	PDBWriter writer;
	writer.open(_outputDir + "/" + _pdbName + ".pdb");
	writer.write(_sys.getAtomPointers(), true, false, false);
	writer.close();
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === HELPER FUNCTIONS ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS){	
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
	Eset->setWeight("CHARMM_IMM1", _opt.weight_IMM1);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_IMM1REF);

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
}

void deleteTerminalBondInteractions(System &_sys, vector<string> &_deleteTerminalInteractionList){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		int firstResidueNumber = 0; // this is kinda done stupidly but it works; removes hbonds first and last 4 positions
		int lastResidueNumber = positions[positions.size()-1]->getResidueNumber();
		for(int i = firstResidueNumber; i < firstResidueNumber+4; i++) { // changed on 2024-2-29 to address hbonding AAs at the interface that don't have a way to hbond in the monomer
			// rid of hbonds from first 4 positions
			if(firstResidueNumber <= i) {
				atoms += positions[i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			// rid of hbonds from last 4 positions
			if(lastResidueNumber > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	for (uint i=0; i<_deleteTerminalInteractionList.size(); i++){
		pESet->deleteInteractionsWithAtoms(atoms,_deleteTerminalInteractionList[i]);
	}
}

// set the active identity for each position to the identity in the given sequence (only for homodimers)
string extractSequence(System &_sys){
	// initialize the sequence string
	string sequence = "";
	// get the first chain from the system
	Chain &chain = _sys.getChain(0);
	// loop through the chain
	for (uint i=0; i<chain.positionSize(); i++){
		// get the ith position in the system
		Position &pos = chain.getPosition(i);
		// get the residue name of the ith position
		string res = pos.getResidueName();
		// convert the residue name to one letter code
		string aa = MslTools::getOneLetterCode(res);
		// add the one letter code to the sequence string
		sequence += aa;
	}
	return sequence;
}

string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		if (it == _seq.begin() || it == _seq.end()-1){
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (it == _seq.begin()){
				if(resName == "HIS") {
					ps = ps + " HSE-ACE";
				} else {
					ps = ps + " " + resName + "-ACE";
				}
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

// sets the gly69 backbone to starting geometry
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB,
 map<string,double> _geometry, Transforms &_trans) {
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

	double xShift = _geometry["xShift"];
	double zShift = _geometry["zShift"];
	double axialRotation = _geometry["axialRotation"];
	double crossingAngle = _geometry["crossingAngle"];

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, _axisA, _axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}


void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
}

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k < _sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
        // check if there is more than one identity for this position
        if (pos.identitySize() > 1) {
            for (uint i=0; i<pos.identitySize(); i++){
                // set the active identity
			    pos.setActiveIdentity(i);
		        if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
		        	if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
		        		cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
		        	}
		        }
            }
        } else if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    === FUNCTIONS FROM CATM ===
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {
	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, zShift);
	_trans.translate(_apV, interDistVect);
	_trans.translate(_axis, interDistVect);
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
	Matrix m(3,3,0.0);
	m[0][0] = -1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = -1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;

	// Rotate chain B around Z axis
	Transforms trans;
	trans.rotate(_apvB, m);
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {

	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {

	 if (moveType == 0) {
		// Z Shift
		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 1) {
		// Axial Rotation
		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else 	if (moveType == 2) {
		// Crossing Angle
		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 3) {
		// XShift
		// Helix A interhelical distance
		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center
		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);
		_trans.translate(_axisA, translateA);

		// Helix B interhelical distance
		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else {
		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
	}
}

// Just add 10 U(0,1) uniform random variables, offset by 0.5 to make mean = 0 and divide by variance = (10 * var(U(0,1)))
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

void computeMonomerEnergy(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, ofstream &_sout) {
	// output the start time
	outputTime(clockTime, "Compute Monomer Energy Start", _sout);
	
	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	
	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

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
	CSBMono.buildSystemFromPDB(inputChain.getAllAtomPointers());

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
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_IMM1);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_IMM1REF);
	
	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all bonding near the termini of our helices for a list of interactions
    deleteTerminalBondInteractions(monoSys,_opt.deleteTerminalInteractions);
	
	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);
	monoSys.buildAllAtoms();

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
	AtomPointerVector &chainA = monoSys.getAllAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	trans.translate(chainA, interDistVect);

	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	trans.translate(chainA, move5Down);
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
	helicalAxis.saveAltCoor("BestMonomerAxis");
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		trans.translate(chainA, zUnitVector);

		//double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			//helicalAxis.saveAltCoor("BestMonomerAxis");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
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
				helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

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
		helicalAxis.applySavedCoor("BestMonomerAxis");

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
			trans.translate(chainA, translateA);
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
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
			helicalAxis.saveAltCoor("BestMonomerAxis");
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
	helicalAxis.applySavedCoor("BestMonomerAxis");
	monoSpm.runGreedyOptimizer(_opt.greedyCycles);
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
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "_Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	_sequenceEnergyMap[_seq]["MonomerSasa"] = totalMonomerSasa;

	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	helicalAxis.clearSavedCoor("BestMonomerAxis");
	helicalAxis.clearSavedCoor("bestZ");
	monoSys.printEnergySummary();
	writePdb(monoSys, _opt.outputDir, "monomer");

	// output end time
	outputTime(clockTime, "Compute Monomer Energy End", _sout);
}

void computeMonomerEnergy(System &_sys, System &_helicalAxis, Options &_opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq,
 RandomNumberGenerator &_RNG, ofstream &_sout) {
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
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_IMM1);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_IMM1REF);
	
	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all bonding near the termini of our helices for a list of interactions
    deleteTerminalBondInteractions(monoSys,_opt.deleteTerminalInteractions);
	
	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);
	monoSys.buildAllAtoms();

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
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "_Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	_sequenceEnergyMap[_seq]["MonomerSasa"] = totalMonomerSasa;

	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	_helicalAxis.clearSavedCoor("BestMonomerAxis");
	_helicalAxis.clearSavedCoor("bestZ");
	monoSys.printEnergySummary();
	writePdb(monoSys, _opt.outputDir, "monomer");

	// output end time
	outputTime(clockTime, "Compute Monomer Energy End", _sout);
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    === OPTIONS ===

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
Options parseOptions(int _argc, char * _argv[]){

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file
	 ******************************************/

	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuration file:
	 *
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/seqDesign.config
	 *
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	// optional
	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_IMM1");
	opt.allowed.push_back("weight_IMM1REF");

	// booleans
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("energyTermList");
	opt.allowed.push_back("deleteTerminalBonds");
	opt.allowed.push_back("deleteTerminalInteractions");
	
	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("density");
	opt.allowed.push_back("negAngle");
	opt.allowed.push_back("negRot");
	opt.allowed.push_back("thread");

	// load rotamers
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("useSasaBurial");
	opt.allowed.push_back("interfaceLevel");
	opt.allowed.push_back("SL");

	//Input Files
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("backboneFile");
	opt.allowed.push_back("pdbFile");
	opt.allowed.push_back("helicalAxis");
	opt.allowed.push_back("configfile");

	//
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("thread");

	//
	opt.allowed.push_back("seed");
	opt.allowed.push_back("greedyCycles");
	
	//monomer options
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");

	// Backbone Monte Carlo variables
	opt.allowed.push_back("backboneMCCycles");
	opt.allowed.push_back("backboneMCMaxRejects");
	opt.allowed.push_back("backboneMCStartTemp");
	opt.allowed.push_back("backboneMCEndTemp");
	opt.allowed.push_back("backboneMCCurve");
	opt.allowed.push_back("backboneConvergedSteps");
	opt.allowed.push_back("backboneConvergedE");
	opt.allowed.push_back("backboneRepackCycles");

	// Shift Size
	opt.allowed.push_back("deltaX");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaXLimit");
	opt.allowed.push_back("deltaCrossLimit");
	opt.allowed.push_back("deltaAxLimit");
	opt.allowed.push_back("deltaZLimit");
	opt.allowed.push_back("decreaseMoveSize");

	//Begin Parsing through the options
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		opt.errorMessages += "No options given!\n";
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPErrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	//Energy Terms to Output
	opt.energyTermList = OP.getStringVector("energyTermList");
	if (OP.fail()) {
		//TODO: write in an error that will tell you if the above is the case
		opt.energyTermList.push_back("CHARMM_VDW");
		opt.energyTermList.push_back("SCWRL4_HBOND");
		opt.energyTermList.push_back("CHARMM_IMM1");
		opt.energyTermList.push_back("CHARMM_IMM1REF");
	}
	opt.deleteTerminalBonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalBonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalInteractions = OP.getMultiString("deleteTerminalInteractions");
	if (OP.fail()) {
		opt.deleteTerminalInteractions.push_back("SCWRL4_HBOND");
		opt.warningMessages += "deleteTerminalInteractions not specified, defaulting to delete SCWRL4_HBOND\n";
		opt.warningFlag = true;
	}
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	
	// starting geometry
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified\n";
		opt.warningFlag = true;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified\n";
		opt.warningFlag = true;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified\n";
		opt.warningFlag = true;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified\n";
		opt.warningFlag = true;
	}
	opt.density = OP.getDouble("density");
	if (OP.fail()) {
		opt.warningMessages += "density not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.density = 0;
	}
	opt.negAngle = OP.getBool("negAngle");
	if (OP.fail()) {
		opt.warningMessages += "negAngle not specified using false\n";
		opt.warningFlag = true;
		opt.negAngle = false;
	}
	if (opt.negAngle == true){
		opt.crossingAngle = -opt.crossingAngle;
	}
	opt.negRot = OP.getBool("negRot");
	if (OP.fail()) {
		opt.warningMessages += "negRot not specified using false\n";
		opt.warningFlag = true;
		opt.negRot = false;
	}
	if (opt.negRot == true){
		opt.axialRotation = -opt.axialRotation;
		//opt.axialRotation = -opt.axialRotation;//for CATM geometries?
	}
	
	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL95.00\n";
		opt.SL = "SL95.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}
	opt.useSasaBurial = OP.getBool("useSasaBurial");
	if (OP.fail()) {
		opt.warningMessages += "useSasaBurial not specified, default true\n";
		opt.warningFlag = true;
		opt.useSasaBurial = true;
	}
	opt.sasaRepackLevel = OP.getMultiString("sasaRepackLevel");
	if (OP.fail()) {
		opt.warningMessages += "sasaRepacklevel not specified! Default to one level at " + opt.SL;
		opt.sasaRepackLevel.push_back(opt.SL);
	}
	opt.interfaceLevel = OP.getInt("interfaceLevel");
	if (OP.fail()) {
		opt.warningMessages += "interfaceLevel not specified using 1\n";
		opt.warningFlag = true;
		opt.interfaceLevel = 1;
	}

    //Weights
	opt.weight_vdw = OP.getDouble("weight_vdw");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_vdw not specified, default 1.0\n";
		opt.weight_vdw = 1.0;
	}
	opt.weight_hbond = OP.getDouble("weight_hbond");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_hbond not specified, default 1.0\n";
		opt.weight_hbond = 1.0;
	}
	opt.weight_IMM1 = OP.getDouble("weight_IMM1");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_IMM1 not specified, default 1.0\n";
		opt.weight_IMM1 = 1.0;
	}
	opt.weight_IMM1REF = OP.getDouble("weight_IMM1REF");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_IMM1REF not specified, default 1.0\n";
		opt.weight_IMM1REF = 1.0;
	}

    //Parameter files
	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(envVar)) {
			opt.topFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.parFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_SOLV";
		if(SYSENV.isDefined(envVar)) {
			opt.solvFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(SYSENV.isDefined(envVar)) {
			opt.rotLibFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + ", defaulting to " + SYSENV.getEnv(envVar) + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}
	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine backboneCrd";
		opt.errorFlag = true;
	}
	opt.backboneFile = OP.getString("backboneFile");
	if (OP.fail()) {
		opt.warningMessages += "backboneFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.backboneFile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	opt.hbondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hbondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hbondFile not specified using " + opt.hbondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hbondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages += "pdbFile not specified";
		opt.errorFlag = true;
	}

	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.errorMessages += "helicalAxis file not specified";
		opt.errorFlag = true;
	}
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine outputDir, default to current directory\n";
		opt.warningFlag = true;
	}
	// backbone repack variables
	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaX = 0.5;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 3.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 3.0;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 3.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 3.0;
	}
	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.5;
	}
	opt.deltaXLimit = OP.getDouble("deltaXLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaXLimit not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaXLimit = 0.1;
	}
	opt.deltaCrossLimit = OP.getDouble("deltaCrossLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaCrossLimit not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaCrossLimit = 1.0;
	}
	opt.deltaAxLimit = OP.getDouble("deltaAxLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaAxLimit not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaAxLimit = 1.0;
	}
	opt.deltaZLimit = OP.getDouble("deltaZLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaZLimit not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaZLimit = 0.1;
	}
	opt.decreaseMoveSize = OP.getBool("decreaseMoveSize");
	if (OP.fail()) {
		opt.warningMessages += "decreaseMoveSize not specified using true\n";
		opt.warningFlag = true;
		opt.decreaseMoveSize = true;
	}

	// sequence parameters
	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.warningMessages += "sequence not specified, defaulting to polyleu\n";
		opt.warningFlag = true;
		opt.sequence = "";
	}
	opt.interface = OP.getString("interface");
	if (OP.fail()) {
		opt.warningMessages += "interface not specified\n";
		opt.warningFlag = true;
		opt.interface = "";
	}
	// if sequence is specified, define interface for sequence
	if (opt.sequence != "" && opt.interface == "") {
		opt.interface = "";
		for (uint i=0; i<opt.sequence.length(); i++) {
			// assumes polyLeu backbone
			if (i > 3 || i < opt.sequence.length()-5) {
				// if not Leu, then interface
				if (opt.sequence[i] != 'L') {
					opt.interface += "1";
				} else {
					opt.interface += "0";
				}
			} else {
				opt.interface += "0";
			}
		}
	}
	if (opt.interface != "") {
		if (opt.sequence != ""){
			if (opt.interface.length() != opt.sequence.length()) {
				opt.errorMessages += "interface string and backbone length must be the same length\n";
				opt.errorFlag = true;
			}
		}
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 23\n";
		opt.warningFlag = true;
		opt.thread = 23;
	}

	// alternate identities
	opt.alternateIds = OP.getStringVector("alternateIds");
	if (OP.fail()) {
		opt.warningMessages += "No alternate AA identities given. If error, make sure they are space separated\n";
		opt.warningFlag = true;
	}
	// Monomer Options
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "MCCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.MCCycles = 10;
	}
	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "MCMaxRejects not specified using 2\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 2;
	}

	// Backbone Monte Carlo parameters
	opt.backboneMCCycles = OP.getInt("backboneMCCycles");
	if (OP.fail()) {
		opt.backboneMCCycles = 100;
		opt.warningMessages += "Number of backboneMC cycles not specified, default to 100\n";
		opt.warningFlag = true;
	}
	opt.backboneMCMaxRejects = OP.getInt("backboneMCMaxRejects");
	if (OP.fail()) {
		opt.backboneMCMaxRejects = 5;
		opt.warningMessages += "Number of backboneMC max rejects not specified, default to using 5\n";
		opt.warningFlag = true;
	}
	opt.backboneMCStartTemp = OP.getDouble("backboneMCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCStartTemp not specified using 100.0\n";
		opt.warningFlag = true;
		opt.backboneMCStartTemp = 100.0;
	}
	opt.backboneMCEndTemp = OP.getDouble("backboneMCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.backboneMCEndTemp = 0.5;
	}
	opt.backboneMCCurve = OP.getInt("backboneMCCurve");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.backboneMCCurve = 2;
	}
	opt.backboneConvergedSteps = OP.getInt("backboneConvergedSteps");
	if (OP.fail()) {
		opt.backboneConvergedSteps = opt.backboneMCCycles/2;
		opt.warningMessages += "backboneConvergedSteps not specified using half of given cycles (" + to_string(opt.backboneConvergedSteps) + "\n";
		opt.warningFlag = true;
	}
	opt.backboneConvergedE = OP.getDouble("backboneConvergedE");
	if (OP.fail()) {
		opt.warningMessages += "backboneConvergedE not specified using 0.001\n";
		opt.warningFlag = true;
		opt.backboneConvergedE = 0.001;
	}
	opt.backboneRepackCycles = OP.getInt("backboneRepackCycles");
	if (OP.fail()) {
		opt.backboneRepackCycles = 5;
		opt.warningMessages += "Number of backbone search cycles not specified, default to 5\n";
		opt.warningFlag = true;
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
