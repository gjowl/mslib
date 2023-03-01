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
string programName = "seqDesign";//TODO: better name
string programDescription = "Designs sequences for backbone geometries extracted from the PDB";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "18 August 2022";
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

// sequence search functions
void searchForBestSequence(System &_startGeom, Options &_opt, PolymerSequence &_PS,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
 vector<unsigned int> &_bestState, string &_bestSequence, vector<uint> &_allInterfacialPositionsList,
 vector<uint> &_interfacialPositionsList, vector<uint> &_rotamerSampling,
 RandomNumberGenerator &_RNG, int _rep, ofstream &_sout, ofstream &_err);
void sequenceSearchMonteCarlo(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<uint> &_rotamerSampling, int _rep, ofstream &_sout, ofstream &_err);
map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap, map<string,double> _sequenceEntropyMap,
 vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList, vector<uint> _rotamerSampling);
string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap, string _energyTerm);
string getSequenceUsingMetropolisCriteria(map<string,map<string,double>> &_sequenceEnergyMap, RandomNumberGenerator &_RNG, double _currTemp, string _energyTerm);
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, double _prevEnergy, string _currSeq, vector<uint> _currVec,
 vector<uint> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap);
//void switchSequence(System &_sys, Options &_opt, string _sequence = "AAALLLLLLLLLLLLLLLAAA");
void switchSequence(System &_sys, Options &_opt, string _sequence);
// function that runs a SelfConsistentMeanField algorithm to find the best starting sequence, only taking into account monomer energies
void runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _rotamerSamplingPositionVector,
 vector<uint> _interfacePositions, map<string, map<string,double>> &_sequenceEnergyMap, map<string, double> _sequenceEntropyMap,
 ofstream &_out);

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

// string functions
string getAlternateIdString(vector<string> _alternateIds);

// help functions
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

// directory setup functions
void setupDesignDirectory(Options &_opt){
	_opt.outputDir = string(get_current_dir_name()) + "/design_" + _opt.runNumber;
	//_opt.outputDir = "/exports/home/gloiseau/mslib/trunk_AS/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (Options &_opt, System &_startGeom, string _seq) {
	// polymer sequences have: chain, starting position of chain residue, three letter AA code
	string polySeq = convertToPolymerSequence(_seq, _opt.thread);
	PolymerSequence PS(polySeq);

	// Declare system for dimer
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

	if(!CSB.buildSystem(PS)) {
		cout << "Unable to build system from " << PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);

	CSB.updateNonBonded(10,12,50);

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	sys.buildAllAtoms();

	string monoPolySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);
	PolymerSequence MPS(monoPolySeq);

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
	if (!CSBMono.buildSystem(MPS)){
		cerr << "Unable to build system from " << monoPolySeq << endl;
	}

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(_opt.backboneCrd);
	if(!cRead.read()) {
		cerr << "Unable to read " << _opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator dimerSasa(sys.getAtomPointers());
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	dimerSasa.calcSasa();
	monoSasa.calcSasa();
	dimerSasa.setTempFactorWithSasa(true);

	for (uint i = 0; i < monoSys.positionSize(); i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdMonomer = monoSys.getPosition(i).getPositionId();
		string posIdDimer = sys.getPosition(i).getPositionId();
		double resiSasaMonomer = monoSasa.getResidueSasa(posIdMonomer);
		double resiSasaDimer = dimerSasa.getResidueSasa(posIdDimer);
		double burial = resiSasaDimer/resiSasaMonomer;
		residueBurial.push_back(pair<int,double>(i, burial));

		//set sasa for each residue in the b-factor
		AtomSelection selA(sys.getPosition(i).getAtomPointers());
		AtomSelection selB(sys.getPosition(i+_opt.backboneLength).getAtomPointers());
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
	// sort in descending order of burial (most buried first)
	sort(residueBurial.begin(), residueBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});
	return residueBurial;
}

// parse config file for given options
Options parseOptions(int _argc, char * _argv[]);

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

	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	// AtomPointerVector for the helical axis for each chain
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	/******************************************************************************
	 *      === COPY BACKBONE COORDINATES AND TRANSFORM TO INPUT GEOMETRY ===
	 ******************************************************************************/
	// get the starting geometry using polyglycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,axisA,axisB,trans);
	
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
	loadRotamers(sys, sysRot, opt.useSasaBurial, opt.sasaRepackLevel, opt.SL, rotamerSamplingPerPosition);
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

	// loop through sequence search and backbone optimization
	// realized that my sequence search works well for a backbone, and my optimization works well for a sequence, so why am I doing them concurrently?
	// instead, I can do more backbone cycles for multiple sequence searches? May be easier to explain than iterative sequence search and optimization
	// is it possible to pass an energy set to multiple systems?
	// TODO: I don't think I can do this; I would need to be able to copy the energy set to a new system, which I don't think I can do
	for (uint i=0; i<opt.backboneSearchCycles; i++){
		vector<uint> bestState;
		// run state Monte Carlo to get a random sequence from the best state
		searchForBestSequence(startGeom, opt, interfacePolySeq, sequenceEnergyMapBest, sequenceEntropyMap, bestState, bestSequence,
		 allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, i, sout, err);
		// add in the starting geometries to the map	
		addGeometryToEnergyMap(startGeometries, sequenceEnergyMapBest, bestSequence);
		// TODO: going to set up this today; save x number of sequences from the search and then optimize them with threading
		// if the sequence is different, optimize the backbone
		if (prevSequence != bestSequence){ 
			// switch the sequence to the best sequence from the search monte carlo
			switchSequence(sys, opt, bestSequence);
			// compute monomer energy
			computeMonomerEnergy(sys, helicalAxis, opt, trans, sequenceEnergyMapBest, bestSequence, RNG, sout, err);
			// optimize the backbone for the sequence
			backboneOptimizer(opt, startGeom, bestSequence, bestState, sequenceEnergyMapBest, i, helicalAxis, axisA, axisB,
			 rotamerSamplingPerPosition, trans, RNG, sout);
			// get the best sequence from the energy map
			sequenceEnergyMapFinalSeqs[bestSequence] = sequenceEnergyMapBest[bestSequence];	
			// reset the energy map
			sequenceEnergyMapBest.clear(); // energyMap to hold all energies for output into a summary file
		} else {
			// end loop if same sequence from seq monte carlo
			i = opt.backboneSearchCycles;	
		}
		prevSequence = bestSequence;// save the current best sequence as the previous sequence for comparison before backbone optimization
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

/****************************************************************
 *             === SEQUENCE SEARCH FUNCTIONS ===
 ****************************************************************/
void searchForBestSequence(System &_startGeom, Options &_opt, PolymerSequence &_PS,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
 vector<unsigned int> &_bestState, string &_bestSequence, vector<uint> &_allInterfacialPositionsList,
 vector<uint> &_interfacialPositionsList, vector<uint> &_rotamerSampling,
 RandomNumberGenerator &_RNG, int _rep, ofstream &_sout, ofstream &_err){
	// output the start time
	outputTime(clockTime, "Monte Carlo Sequence Search Replicate " + to_string(_rep) + " Start", _sout);

	// initialize the system	
	System sys;
	prepareSystem(_opt, sys, _startGeom, _PS, _opt.useBaseline);

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	loadRotamers(sys, sysRot, _opt.useSasaBurial, _opt.sasaRepackLevel, _opt.SL, _rotamerSampling);

	// Setup self pair manager to calculate the self and pair energies for each sequence
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);

	// calculate the self and pair energies for each amino acid and rotamer combo in the sequence	
	outputTime(clockTime, "  - SelfPairManager energy calculation " + to_string(_rep) + " Start", _sout);
	spm.calculateEnergies();
	outputTime(clockTime, "  - SelfPairManager energy calculation " + to_string(_rep) + " End", _sout);
	
	// monte carlo for finding the best sequence
	sequenceSearchMonteCarlo(sys, _opt, spm, _RNG, _bestState, _bestSequence, _sequenceEnergyMap,  
	_sequenceEntropyMap, _allInterfacialPositionsList, _interfacialPositionsList, _rotamerSampling, _rep, _sout, _err);

	outputTime(clockTime, "Monte Carlo Sequence Search Replicate " + to_string(_rep) + " End", _sout);
	// output the best sequence energy
	cout << "Best Sequence: " << _bestSequence << endl;
	cout << sys.getEnergySummary() << endl;	
}

void sequenceSearchMonteCarlo(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<uint> &_rotamerSampling, int _rep, ofstream &_sout, ofstream &_err){
	// Setup time variables
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);

	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	switchSequence(_sys, _opt, _bestSequence); // switch sequence to the best sequence (from SCMF or input sequence)
	vector<vector<bool>> mask = getActiveMask(_sys); // get the active mask for the current sequence
	_spm.runGreedyOptimizer(_opt.greedyCycles, mask); // run the greedy optimizer to get the best state for the current sequence
	_bestState = _spm.getMinStates()[0]; // set the best state to the best state from the greedy optimizer
	_sys.setActiveRotamers(_bestState); // set the active rotamers to the best state
	double bestEnergy = _spm.getStateEnergy(_bestState); // save the energy of the best state
	
	// calculate sequence entropy of the starting sequence and add to the energy map
	double sequenceProbability = 1;
	double sequenceEntropy = calculateInterfaceSequenceEntropy(_bestSequence, _sequenceEntropyMap, _interfacialPositionsList, sequenceProbability, _opt.weight_seqEntropy); // calculate the sequence entropy for the interface

	// add sequence entropy to map
	double baseline = _spm.getStateEnergy(_bestState, "BASELINE")+_spm.getStateEnergy(_bestState, "BASELINE_PAIR");
	double energyAndEntropyTotal = bestEnergy+sequenceEntropy;
	_sequenceEnergyMap[_bestSequence]["Dimerw/Baseline"] = bestEnergy;
	_sequenceEnergyMap[_bestSequence]["Dimer"] = bestEnergy-baseline;
	_sequenceEnergyMap[_bestSequence]["Baseline"] = baseline;
	_sequenceEnergyMap[_bestSequence]["SequenceProbability"] = sequenceProbability;
	_sequenceEnergyMap[_bestSequence]["SequenceEntropy"] = sequenceEntropy;
	_sequenceEnergyMap[_bestSequence]["Totalw/Entropy"] = energyAndEntropyTotal;
	bestEnergy = energyAndEntropyTotal;

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	MC.setEner(bestEnergy);

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	int acceptCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);
	cout << "Sequence at start of stateMC: " << prevStateSeq << endl;

	// Sequence Search Energy Landscape file
	ofstream lout;
	string loutfile  = _opt.outputDir + "/sequenceSearchEnergyLandscape_" + to_string(_rep) + ".out";
	lout.open(loutfile.c_str());
	lout << "Accept/Reject\tCycle\tPrevSequence\tCurrSequence\tPrevEnergy\tCurrEnergy\tPrevEntropy\tCurrEntropy\tCurrentTemp\tAccepts" << endl;

	// initialize energy variables for the MonteCarlo
	string bestSeq = prevStateSeq;
	double prevStateEntropy = sequenceEntropy;
	map<string,double> bestSequenceEnergyMap;
	bool randomMCRun = false;
	bool decreaseMCRun = false;

	// Monte Carlo while loop for finding the best sequences
	string energyTermToEvaluate = "Totalw/Entropy";
	double bestProb = sequenceProbability;
	while (!MC.getComplete()){
		// get the sequence entropy probability for the current best sequence
		map<string,vector<uint>> sequenceVectorMap;
		map<string,map<string,double>> sequenceEnergyMap = mutateRandomPosition(_sys, _opt, _spm, _RNG, bestSeq, bestEnergy, 
		 sequenceVectorMap, _sequenceEntropyMap, _allInterfacialPositionsList, _interfacialPositionsList, _rotamerSampling);
		
		// get the best sequence and energy for the current mutation position (picks sequence with energy including vdw, hbond, imm1, baseline, sequence entropy)
		//string currSeq = getBestSequenceInMap(sequenceEnergyMap);
		double currTemp = MC.getCurrentT();
		string currSeq = getSequenceUsingMetropolisCriteria(sequenceEnergyMap, _RNG, currTemp, energyTermToEvaluate);
		// extract energies from the sequence energy map for the chosen sequence
		double currEnergy = sequenceEnergyMap[currSeq][energyTermToEvaluate]; // energy total for current sequence (dimer+baseline+sequence entropy)
		double currStateEntropy = sequenceEnergyMap[currSeq]["SequenceEntropy"];
		// compare the entropy between the current sequence and previous sequence 
		vector<uint> currStateVec = sequenceVectorMap[currSeq]; // current state vector for current sequence (best rotamers and sequence identities)

		//MC.setEner(bestEnergyTotal);
		// MC accept and reject conditions
		double acceptTemp = MC.getCurrentT();
		if (!MC.accept(currEnergy)){
			setActiveSequence(_sys, bestSeq);
			_sys.setActiveRotamers(prevStateVec); // set rotamers to the previous state
			lout << "Reject\t" << cycleCounter << "\t" << prevStateSeq << "\t" << currSeq << "\t" << bestEnergy << "\t" << currEnergy << "\t";
			lout << prevStateEntropy << "\t" << currStateEntropy << "\t" << acceptTemp << "\t" << acceptCounter << endl;
		} else {
			acceptCounter++; // increase the accept counter by 1
			//bestProb = currProb;
			setActiveSequence(_sys, currSeq);
			_sys.setActiveRotamers(currStateVec); // set rotamers to the current state
			map<string,double> energyMap = sequenceEnergyMap[currSeq];
			lout << "Accept\t" << cycleCounter << "\t" << prevStateSeq << "\t" << currSeq << "\t" << bestEnergy << "\t" << currEnergy << "\t";
			lout << prevStateEntropy << "\t" << currStateEntropy << "\t" << acceptTemp << "\t" << acceptCounter << endl;
			bestEnergy = currEnergy;
			_bestState = currStateVec; // set the 
			prevStateVec = currStateVec; // set the previous state vector to be the current state vector
			prevStateSeq = currSeq; // set the previous sequence to be the current sequence
			bestSeq = currSeq; // set the best sequence to the newly accepted current sequence
			bestEnergy = sequenceEnergyMap[bestSeq][energyTermToEvaluate]; // set the best energy to the current energy (vdw, hbond, imm1, baseline)
			sequenceEnergyMap[bestSeq]["acceptCycleNumber"] = acceptCounter; // gets the accept cycle number for the current sequence
			sequenceEnergyMap[bestSeq]["acceptTemp"] = acceptTemp; // gets the accept cycle number for the current sequence
			bestSequenceEnergyMap = sequenceEnergyMap[bestSeq]; // saves the current sequence energy map 
			prevStateEntropy = currStateEntropy;
		}
		//Reset the MC to run x more cycles
		if (MC.getComplete() == true && randomMCRun == false){
			//MC.reset(3649, 3649, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
			//MC.reset(547, 547, 100, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);//Approximately 25% likely to accept within 5kcal
			MC.reset(_opt.MCResetTemp, _opt.MCResetTemp, _opt.MCResetCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);//Approximately 25% likely to accept within 5kcal
			randomMCRun = true;
		} 
		// run more cycles of MC to decrease
		if (MC.getComplete() == true && randomMCRun == true && decreaseMCRun == false){
			//MC.reset(547, _opt.MCEndTemp, 100, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
			MC.reset(_opt.MCResetTemp, _opt.MCEndTemp, _opt.MCResetCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);//Approximately 25% likely to accept within 5kcal
			decreaseMCRun = true;
		}
		cycleCounter++;
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);
	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();

	_bestSequence = bestSeq;
	_sequenceEnergyMap[bestSeq] = bestSequenceEnergyMap; // saves the best sequence energy map into the global energy map
	bestEnergy = _sequenceEnergyMap[bestSeq]["Dimerw/Baseline"]; // set the best energy to the current energy (vdw, hbond, imm1, baseline)
	
	// output the best sequence pdb
	string pdbOutput = "bestSequence_"+to_string(_rep);
	_sys.setActiveRotamers(_bestState);
	writePdb(_sys, _opt.outputDir, "seqSearch");
	cout << "End monte carlo sequence search #" << _rep << ": " << diffTimeSMC/60 << "min" << endl;
	cout << "Monte Carlo ended at Temp: " << MC.getCurrentT() << endl << endl;
	cout << "- Best Sequence: " << bestSeq << "; Energy w/ baseline: " << bestEnergy << endl;
	_sout << "- End monte carlo sequence search #" << _rep << ": " << diffTimeSMC/60 << "min; ";
	_sout << "Monte Carlo ended at Temp: " << MC.getCurrentT() << endl << endl;
	_sout << "- Best Sequence: " << bestSeq << "; Energy w/ baseline: " << bestEnergy << endl;

	vector<vector<unsigned int>> statePositionIdRotamerIndeces = _spm.getStatePositionIdentityRotamerIndeces(_bestState);
	vector<uint> rotamerState;
	// offset the rotamer indeces whenever the residue is GLY, ALA, or PRO since they only have one rotamer
	// the below is definitely a hack, but I couldn't figure out how to do it and was too frustrated after hours of thinking
	if (_opt.useAlaAtTermini){
		uint offset = 3;
		for (uint i=3; i<_sys.getPositions().size()-3; i++){
			Position &pos = _sys.getPosition(i);
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				rotamerState.push_back(statePositionIdRotamerIndeces[i-offset][2]);
				//cout << statePositionIdRotamerIndeces[i-offset][2] << ",";
			} else {
				if (i > 17 && i < 24 || pos.getResidueName() == "ALA"){
					offset++;
				}
			}
		}
	} else {
		for (uint i=0; i<statePositionIdRotamerIndeces.size(); i++){
			Position &pos = _sys.getPosition(i);
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				rotamerState.push_back(statePositionIdRotamerIndeces[i][2]);
				cout << statePositionIdRotamerIndeces[i][2] << ",";
			} 
		}
	}
	cout << endl;	
	_bestState = rotamerState;
}

map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap, map<string,double> _sequenceEntropyMap,
 vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList, vector<uint> _rotamerSampling){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _interfacialPositionsList.size()-1);
	int interfacePosA = _interfacialPositionsList[rand];
	int interfacePosB = interfacePosA+_bestSeq.length();

	// Get the random position from the system
	Position &randPosA = _sys.getPosition(interfacePosA);
	Position &randPosB = _sys.getPosition(interfacePosB);
	string posIdA = randPosA.getPositionId();
	string posIdB = randPosB.getPositionId();

	// variable setup for current state
	map<string,map<string,double>> sequenceEnergyMap;
	vector<thread> threads;
	for (uint i=0; i<_opt.Ids.size(); i++){
		// pick an identity for each thread 
		int idNum = i;
		// generate polymer sequence for each identity at the corresponding chosen position
		string id = _opt.Ids[idNum];
		// input into the thread function for calculating energies
		string currAA = MslTools::getThreeLetterCode(_bestSeq.substr(interfacePosA, 1));
		if (currAA != id){
			// replace the id at the position in bestSeq with the current id to get current sequence
			string currSeq = _bestSeq;
			string oneLetterId = MslTools::getOneLetterCode(id);
			currSeq.replace(interfacePosA, 1, oneLetterId);
			// switch to the mutated sequence
			setActiveSequence(_sys, currSeq);
			// Set a mask and run a greedy to get the best state for the current sequence
			vector<vector<bool>> mask = getActiveMask(_sys);
			_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
			vector<uint> currVec = _spm.getMinStates()[0];
			_sys.setActiveRotamers(currVec);
			_sequenceStateMap[currSeq] = currVec;
			// start threading and calculating energies for each identity; changed the below interface from allInterface to just the given interface since these are the ones that change and ala at ends as of 2022-10-18
			threads.push_back(thread{energyFunction, ref(_opt), ref(_spm), _bestSeq, _bestEnergy, currSeq, currVec, ref(_rotamerSampling), ref(_interfacialPositionsList), 
			 ref(sequenceEnergyMap), ref(_sequenceEntropyMap)});
		}
		setActiveSequence(_sys, _bestSeq);
	} 
	// join all the threads (wait for them all to finish before continuing)
	for (auto& th : threads){
		th.join();
	}
	return sequenceEnergyMap;
}

// threaded function that gets energies for a sequence and appends to a map
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, double _prevEnergy, string _currSeq, vector<uint> _currVec,
 vector<uint> &_rotamerSampling, vector<uint> &_interfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap){
	// variable setup
	map<string,double> energyMap;

	// Compute dimer energy
	outputEnergiesByTerm(_spm, _currVec, energyMap, _opt.energyTermList, "Dimer", true);
	double currEnergy = _spm.getStateEnergy(_currVec);

	// initialize variables for this thread
	double sequenceProbability = 1;
	double sequenceEntropy = calculateInterfaceSequenceEntropy(_currSeq, _sequenceEntropyMap, _interfacePositions, sequenceProbability, _opt.weight_seqEntropy); // calculate the sequence entropy for the interface

	// output info
	double baseline = _spm.getStateEnergy(_currVec, "BASELINE")+_spm.getStateEnergy(_currVec, "BASELINE_PAIR");
	double energyAndEntropyTotal = currEnergy+sequenceEntropy;
	energyMap["Dimerw/Baseline"] = currEnergy;
	energyMap["Dimer"] = currEnergy-baseline;
	energyMap["Baseline"] = baseline;
	energyMap["SequenceProbability"] = sequenceProbability;
	energyMap["SequenceEntropy"] = sequenceEntropy;
	energyMap["Totalw/Entropy"] = energyAndEntropyTotal;
	_seqEnergyMap[_currSeq] = energyMap;
}

string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap, string _energyTerm){
	string currSeq;
	double currEnergyComparison;
	uint i=0;
	for (auto& seq : _sequenceEnergyMap){
		if (i==0){
			currSeq = seq.first;
			currEnergyComparison = seq.second[_energyTerm];
			i++;
		} else {
			if (seq.second[_energyTerm] < currEnergyComparison){
				string prevSeq = currSeq;
				currSeq = seq.first;
				currEnergyComparison = seq.second[_energyTerm];
			}
		}
	}
	return currSeq;
}

string getSequenceUsingMetropolisCriteria(map<string,map<string,double>> &_sequenceEnergyMap, RandomNumberGenerator &_RNG, double _currTemp, string _energyTerm){
	// get a random sequence from the map
	int rand = _RNG.getRandomInt(0, _sequenceEnergyMap.size()-1);	
	map<string,map<string,double>>::iterator it = _sequenceEnergyMap.begin();
	advance(it, rand);
	string currSeq = it->first;
	double lastEnergy = it->second[_energyTerm];
	// create a vector to hold sequences that are already chosen; allows me to randomly traverse the map
	vector<string> chosenSeqs;
	chosenSeqs.push_back(it->first);
	// loop through map size
	for (uint i=0; i<_sequenceEnergyMap.size()-1; i++){
		// get a random sequence from the map
		int rand = _RNG.getRandomInt(0, _sequenceEnergyMap.size()-1);	
		map<string,map<string,double>>::iterator it = _sequenceEnergyMap.begin();
		// check if random sequence has already been chosen
		if (find(chosenSeqs.begin(), chosenSeqs.end(), it->first) != chosenSeqs.end()){
			// if it has been chosen, get a new random sequence
			while (find(chosenSeqs.begin(), chosenSeqs.end(), it->first) != chosenSeqs.end()){
				rand = _RNG.getRandomInt(0, _sequenceEnergyMap.size()-1);	
				it = _sequenceEnergyMap.begin();
				advance(it, rand);
			}
		}
		// add the chosen sequence to the list of chosen sequences
		chosenSeqs.push_back(it->first);
		// advance to the random sequence
		double energy = it->second[_energyTerm];
		/***********************************
		 *  Metropolis criterion
		 ***********************************/
		double exp = pow(M_E,(lastEnergy - energy) / (MslTools::R * _currTemp));
		//cout << "currSeq: " << currSeq << ";energy: " << energy << ";lastEnergy: " << lastEnergy << ";exp: " << exp << endl;
		if (exp > 1) {
			lastEnergy = energy;
			currSeq = it->first;
			//cout << "currSeq: " << currSeq << ";energy: " << energy << ";lastEnergy: " << lastEnergy << endl;
		} else {
			double randomN = _RNG.getRandomDouble();
			//cout << exp << " < " << randomN << endl;
			if (exp > randomN) {
				lastEnergy = energy;
				currSeq = it->first;
			}
		}
	}
	return currSeq;
}


/****************************************************************
 *       === IDENTIFY INTERFACE USING RESIDUE BURIAL ===
 ****************************************************************/
// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
//TODO: fix this to start...make this entire thing less complicated, don't need as many vectors, just two or three (rotamer levels, interface, and allInterface)
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<uint> &_rotamerSamplingPerPosition){
	
	// generate a backboneSequence to determine the interface positions using residue burial (defaults to poly-Valine sequence)
	string polyVal = generateString("V", _opt.backboneLength);

	// save into vector of backbone positions and residue burial pairs
	vector<pair <int, double> > resiBurial = calculateResidueBurial(_opt, _startGeom, polyVal);
	vector<int> interfacePositions;
	
	// if sequence is not empty, use polyLeu
	// I feel like heterodimerization would be best as an option in this case, but I would have to set everything up with an option, which isn't ideal
	// for now, I think my initial idea of starting the code as an alternate hetero version is best, and then I can set the homodimerization to be the 
	// same if need be?
	string backboneSeq;
	if (_opt.sequence == ""){
		backboneSeq = generateBackboneSequence(_opt.backboneAA, _opt.backboneLength, _opt.useAlaAtTermini);
	} else {
		backboneSeq = _opt.sequence;
	}

	// setup 0 string to represent variable positions and rotamer levels
	string variablePositionString = generateString("0", _opt.backboneLength);
	string rotamerLevels = generateString("0", _opt.backboneLength);

	// define rotamer levels for each position based on residue burial
	defineRotamerLevels(_opt, resiBurial, interfacePositions, rotamerLevels, variablePositionString);

	// checks if interface is defined; if so, check the position and set those positions to the highest rotamer level
	if (_opt.interface != ""){
		interfacePositions.clear(); // reset the interface from the defineRotamerLevels function
		useInputInterface(_opt, variablePositionString, rotamerLevels, interfacePositions);
	}	

	// save the rotamer levels for all positions  
	vector<uint> rotamerSamplingPerPosition = convertStringToVectorUint(rotamerLevels); // converts the rotamer sampling for each position as a vector
	
	// define the interface positions for the sequence search	
	if (_opt.interface != ""){
		_interfacePositions.clear();
		_allInterfacePositions.clear();
		for (uint k=0; k<interfacePositions.size(); k++){
			int pos = interfacePositions[k];
			_interfacePositions.push_back(pos-_opt.thread);
		}
		for (uint k=0; k<backboneSeq.length(); k++){
			if (k < 3 || k > backboneSeq.length()-5){
				if (rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
					_allInterfacePositions.push_back(k);
				}
			} else {
		 		if (_opt.interface[k] == '1'){
					_allInterfacePositions.push_back(k);
				}
			}
		}
	} else {
		_interfacePositions = getInterfacePositions(_opt.interfaceLevel, rotamerSamplingPerPosition, backboneSeq.length());
		_allInterfacePositions = getAllInterfacePositions(_opt.interfaceLevel, rotamerSamplingPerPosition, backboneSeq.length());
	}

	//String for the positions of the sequences that are considered interface for positions amd high rotamers
	string rotamerSamplingString = convertVectorUintToString(rotamerSamplingPerPosition);
	cout << "Backbone:           " << backboneSeq << endl;
	cout << "Variable Positions: " << variablePositionString << endl;
	cout << "Rotamers Levels:    " << rotamerSamplingString << endl;

	// Define referenced output variables
	_rotamerSamplingPerPosition = rotamerSamplingPerPosition;

	// makes the polymer sequence to return basedon the interface positions
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.Ids, interfacePositions);
	PolymerSequence PS(polySeq);
	return PS;
}

// define the rotamer levels for each position based on residue burial
void defineRotamerLevels(Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString){
	// initialize variables
	int levelCounter = 0;
	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	// loop through the residue burial values calculated above for each position (these are ordered by most to least buried)
	for (uint i = 0; i < _resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(_resiBurial.size()); // calculate the SASA percentile for this position
		// if percentile is greater, move on to the next rotamer level
		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
			levelCounter++;
		}
		int positionNumber = _resiBurial[i].first; // get the backbone position for this position
		int resiNum = positionNumber+_opt.thread;
		// check if the current interface level is below the accepted option for interface levels (SASA repack level)
		if (levelCounter < _opt.interfaceLevel) {
			_interfacePositions.push_back(resiNum);
			// check to see if the position is found within the core of protein (i.e. not the first 3 residues or the last 4 residues)
			if (positionNumber > 2 && positionNumber < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 2 and 4 here instead of 3 and 5 to prevent changes at the interface like others
				// replace 0 with 1 for variable positions that are found at the interface
				_variablePositionString.replace(_variablePositionString.begin()+positionNumber, _variablePositionString.begin()+positionNumber+1, "1");
			}
		}
		_rotamerLevels.replace(_rotamerLevels.begin()+positionNumber, _rotamerLevels.begin()+positionNumber+1, MslTools::intToString(levelCounter));
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

	// clear the saved repack state
	sys.clearSavedCoor("savedRepackState");
	_helicalAxis.clearSavedCoor("BestRepack");

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

	_sys.saveAltCoor("savedRepackState");
	_helicalAxis.saveAltCoor("BestRepack");
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
				if (decreaseMoveSize == true){
					double endTemp = MCMngr.getCurrentT();
					getCurrentMoveSizes(startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, _opt.deltaXLimit,
					 _opt.deltaCrossLimit, _opt.deltaAxLimit, _opt.deltaZLimit, decreaseMoveSize);
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
					getCurrentMoveSizes(startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, _opt.deltaXLimit,
					_opt.deltaCrossLimit, _opt.deltaAxLimit, _opt.deltaZLimit, decreaseMoveSize);
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
	_helicalAxis.applySavedCoor("BestRepack");
	repackSideChains(_spm, _opt.greedyCycles);
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

// define interface functions
void useInputInterface(Options &_opt, string &_variablePositionString, string &_rotamerLevels, vector<int> &_interfacePositions){
	for (uint i=3; i<_opt.interface.length()-4; i++){
		if (_opt.interface[i] == '0'){
			if (_variablePositionString[i] != '0'){
				_variablePositionString.replace(_variablePositionString.begin()+i, _variablePositionString.begin()+i+1, "0");
			}
		} else if (_opt.interface[i] == '1'){
			if (_rotamerLevels[i] != '0'){
				_rotamerLevels.replace(_rotamerLevels.begin()+i, _rotamerLevels.begin()+i+1, "0");
			}
			if (_variablePositionString[i] == '0'){
				_variablePositionString.replace(_variablePositionString.begin()+i, _variablePositionString.begin()+i+1, "1");
			}
		}
	}
	for (uint i=0; i<_variablePositionString.length(); i++){
		if (_variablePositionString[i] == '1'){
			_interfacePositions.push_back(i+_opt.thread);
		}
	}
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
    deleteTerminalBondInteractions(monoSys,_opt.deleteTerminalInteractions);

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

string getAlternateIdString(vector<string> _alternateIds){
	string alternateIdsString = "";
	for (uint i=0; i<_alternateIds.size(); i++){
		if (i == _alternateIds.size()-1){
			alternateIdsString += _alternateIds[i];
		} else {
			alternateIdsString += _alternateIds[i] += " ";
		}
	}
	return alternateIdsString;
}

/****************************************
 *
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/
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

	//opt.required.push_back("");
	//opt.allowed.push_back("");
	
	// Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("geometryDensityFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("backboneFile");
	opt.allowed.push_back("selfEnergyFile");
	opt.allowed.push_back("pairEnergyFile");
	opt.allowed.push_back("sequenceEntropyFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("helicalAxis");

	// sequence parameters
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	opt.allowed.push_back("interface");
	opt.allowed.push_back("numberOfSequencesToSave");

	// booleans
	opt.allowed.push_back("getGeoFromPDBData");
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("useSasaBurial");
	opt.allowed.push_back("useTimeBasedSeed");
	opt.allowed.push_back("deleteTerminalBonds");
	opt.allowed.push_back("linkInterfacialPositions");
	opt.allowed.push_back("useAlaAtTermini");
	opt.allowed.push_back("useBaseline");
	opt.allowed.push_back("getRandomAxRotAndZShift");

	// repack parameters
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");


	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("endResNum");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("density");
	opt.allowed.push_back("negAngle");
	opt.allowed.push_back("negRot");
	opt.allowed.push_back("thread");

	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	opt.allowed.push_back("MCConvergedSteps");
	opt.allowed.push_back("MCConvergedE");
	opt.allowed.push_back("MCResetTemp");
	opt.allowed.push_back("MCResetCycles");

	// Backbone Monte Carlo variables
	opt.allowed.push_back("backboneMCCycles");
	opt.allowed.push_back("backboneMCMaxRejects");
	opt.allowed.push_back("backboneMCStartTemp");
	opt.allowed.push_back("backboneMCEndTemp");
	opt.allowed.push_back("backboneMCCurve");
	opt.allowed.push_back("backboneConvergedSteps");
	opt.allowed.push_back("backboneConvergedE");
	opt.allowed.push_back("backboneSearchCycles");

	// use different energy parameters
	opt.allowed.push_back("useIMM1");
	opt.allowed.push_back("useElec");
	opt.allowed.push_back("compareSasa");
	
	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_solv");
	opt.allowed.push_back("weight_elec");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_seqEntropy");

	// Alternate Identities
	opt.allowed.push_back("Ids");

	// MonteCarlo Arguments
	opt.allowed.push_back("numStatesToSave");

	// SelfPairManager Arguments
	opt.allowed.push_back("runDEESingles");
	opt.allowed.push_back("runDEEPairs");
	opt.allowed.push_back("runSCMF");

	// Energy Terms to Output
	opt.allowed.push_back("energyTermList");
	opt.allowed.push_back("deleteTerminalInteractions");

	//Rotamers
	// Load Rotamers from SASA values (from sgfc)
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("interfaceLevel");
	opt.allowed.push_back("SL");
	
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

	opt.allowed.push_back("runNumber");

	// Begin Parsing through the options
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";
	
	if (OP.countOptions() == 0){
		opt.errorMessages += "No options given!\n";
		return opt;
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	
	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
		}
	}

	// Input Files
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
	opt.geometryDensityFile = OP.getString("geometryDensityFile");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine geometryDensityFile, defaulting to original density file\n";
		opt.warningFlag = true;
		opt.geometryDensityFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_09_28_geometryDensityFile.txt";
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
	opt.backboneFile = OP.getString("backboneFile");
	if (OP.fail()) {
		opt.warningMessages += "backboneFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.backboneFile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	opt.selfEnergyFile = OP.getString("selfEnergyFile");
	if (OP.fail()) {
		opt.warningMessages += "selfEnergyFile not specified, default \n";
		opt.warningFlag = true;
		opt.selfEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanSelf_par.txt";
	}
	opt.pairEnergyFile = OP.getString("pairEnergyFile");
	if (OP.fail()) {
		opt.warningMessages += "pairEnergyFile not specified, default \n";
		opt.warningFlag = true;
		opt.pairEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanPair_par.txt";
	}
	opt.sequenceEntropyFile = OP.getString("seqEntropyFile");
	if (OP.fail()) {
		opt.warningMessages += "seqEntropyFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.sequenceEntropyFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/sequenceEntropies.txt";
	}
	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.warningMessages += "helicalAxis not specified\n";
		opt.warningFlag = true;
		opt.helicalAxis = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/parameterFiles/helicalAxis.pdb";
	}
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine outputDir, default to current directory\n";
		opt.warningFlag = true;
	}
	
	// sequence parameters
	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.warningMessages += "sequence not specified, defaulting to polyleu\n";
		opt.warningFlag = true;
		opt.sequence = "";
	}
	opt.backboneAA = OP.getString("backboneAA");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneAA not specified, default to V\n";
		opt.backboneAA = "V";
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 21\n";
		opt.backboneLength = 21;
	}
	opt.numberOfSequencesToSave = OP.getInt("numberOfSequencesToSave");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "numberOfSequencesToSave not specified, default to 5\n";
		opt.numberOfSequencesToSave = 5;
	}
	// override backboneLength if a sequence is specified
	if (opt.sequence != "") {
		opt.backboneLength = opt.sequence.length();
	}
	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(1) + "\n";
		opt.warningFlag = true;
		opt.startResNum = 1;
	}
	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.startResNum+opt.backboneLength) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.startResNum+opt.backboneLength;
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
		if (opt.interface.length() != opt.backboneLength) {
			opt.errorMessages += "interface string and backbone length must be the same length\n";
			opt.errorFlag = true;
		}
	}

	// booleans
	opt.getGeoFromPDBData = OP.getBool("getGeoFromPDBData");
	if (OP.fail()) {
		opt.warningMessages += "getGeoFromPDBData not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = false;
	}
	opt.deleteTerminalBonds = OP.getBool("deleteTerminalBonds");
	if (OP.fail()) {
		opt.deleteTerminalBonds = true;
		opt.warningMessages += "deleteTerminalBonds not specified using true\n";
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
	opt.linkInterfacialPositions = OP.getBool("linkInterfacialPositions");
	if (OP.fail()) {
		opt.warningMessages += "linkInterfacialPositions not specified using true\n";
		opt.warningFlag = true;
		opt.linkInterfacialPositions = true;
	}
	opt.useTimeBasedSeed = OP.getBool("useTimeBasedSeed");
	if (OP.fail()) {
		opt.warningMessages += "useTimeBasedSeed not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useTimeBasedSeed = false;
	}
	opt.useAlaAtTermini = OP.getBool("useAlaAtTermini");
	if (OP.fail()) {
		opt.warningMessages += "useAlaAtTermini not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useAlaAtTermini = false;
	}
	opt.useBaseline = OP.getBool("useBaseline");
	if (OP.fail()) {
		opt.warningMessages += "useBaseline not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useBaseline = false;
	}
	opt.getRandomAxRotAndZShift = OP.getBool("getRandomAxRotAndZShift");
	if (OP.fail()) {
		opt.warningMessages += "getRandomAxRotAndZShift not specified, defaulting to false";
		opt.warningFlag = true;
		opt.getRandomAxRotAndZShift = false;
	}
	//use different energy parameters
	opt.useIMM1 = OP.getBool("useIMM1");
	if (OP.fail()) {
		opt.warningMessages += "useIMM1 not specified, defaulting to true\n";
		opt.warningFlag = true;
		opt.useIMM1 = true;
	}
	opt.useElec = OP.getBool("useElec");
	if (OP.fail()) {
		opt.warningMessages += "useElec not specified using false\n";
		opt.warningFlag = true;
		opt.useElec = false;
	}
	opt.compareSasa = OP.getBool("compareSasa");
	if (OP.fail()) {
		opt.warningMessages += "compareSasa not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.compareSasa = false;
	}
	
	// starting geometry
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
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
		opt.axialRotation = opt.axialRotation-100;
		//opt.axialRotation = -opt.axialRotation;//for CATM geometries?
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.thread = 0;
	}

	//Load Rotamers using SASA values (from sgfc)
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
	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL95.00\n";
		opt.SL = "SL95.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}

	//Monte Carlo parameters
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MCResetCycles not specified, using 100\n";
		opt.warningFlag = true;
		opt.MCCycles = 100;
	}
	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.MCMaxRejects = 10;
		opt.warningMessages += "Number of MC max rejects not specified, default to using 10\n";
		opt.warningFlag = true;
	}
	opt.MCStartTemp = OP.getDouble("MCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCStartTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.MCStartTemp = 1000.0;
	}
	opt.MCEndTemp = OP.getDouble("MCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.MCEndTemp = 0.5;
	}
	opt.MCCurve = OP.getInt("MCCurve");
	if (OP.fail()) {
		opt.warningMessages += "MCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.MCCurve = 2;
	}
	opt.MCConvergedSteps = OP.getInt("MCConvergedSteps");
	if (OP.fail()) {
		opt.warningMessages += "MCConvergedSteps not specified using 10\n";
		opt.warningFlag = true;
		opt.MCConvergedSteps = 10;
	}
	//opt.MCConvergedE = OP.getDouble("MCConvergedE");
	//if (OP.fail()) {
	//	opt.warningMessages += "MCConvergedE not specified using 0.01\n";
	//	opt.warningFlag = true;
	//	opt.MCConvergedE = 0.01;
	//}
	opt.MCResetTemp = OP.getDouble("MCResetTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCResetTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.MCResetTemp = 3649; // 50% likelihood of accepting a solution within 5kcals
	}
	opt.MCResetCycles = OP.getInt("MCResetCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MCResetCycles not specified, using 100\n";
		opt.warningFlag = true;
		opt.MCResetCycles = 100;
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
	opt.backboneSearchCycles = OP.getInt("backboneSearchCycles");
	if (OP.fail()) {
		opt.backboneSearchCycles = 5;
		opt.warningMessages += "Number of backbone search cycles not specified, default to 5\n";
		opt.warningFlag = true;
	}

	// repack parameters	
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
	}

	// energy weights
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
	opt.weight_solv = OP.getDouble("weight_solv");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_solv not specified, default 1.0\n";
		opt.weight_solv = 1.0;
	}
	opt.weight_seqEntropy = OP.getDouble("weight_seqEntropy");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_seqEntropy not specified, default 1.0\n";
		opt.weight_seqEntropy = 1.0;
	}
	opt.weight_elec = OP.getDouble("weight_elec");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_elec not specified, default 1.0\n";
		opt.weight_elec = 1.0;
	}

	// alternate identities
	opt.Ids = OP.getStringVector("Ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	// String for the alternateIds at the interface
	if (opt.xShift <= 7.5){
		opt.Ids.push_back("GLY");
	}

	//SelfPairManager Optimization Options
	opt.runDEESingles = OP.getBool("runDEESingles");
	if (OP.fail()) {
		opt.warningMessages += "runDEESingles not specified, defaulting to false";
		opt.warningFlag = true;
		opt.runDEESingles = false;
	}
	opt.runDEEPairs = OP.getBool("runDEEPairs");
	if (OP.fail()) {
		opt.warningMessages += "runDEEPairs not specified, defaulting to false";
		opt.warningFlag = true;
		opt.runDEEPairs = false;
	}
	opt.runSCMF = OP.getBool("runSCMF");
	if (OP.fail()) {
		opt.warningMessages += "runSCMF not specified, defaulting to true";
		opt.warningFlag = true;
		opt.runSCMF = true;
	}

	//Energy Terms to Output
	opt.energyTermList = OP.getStringVector("energyTermList");
	if (OP.fail()) {
		//This works, but I think if you ever want to output more terms in the future, need to add them to the terms above
		//TODO: write in an error that will tell you if the above is the case
		opt.energyTermList.push_back("CHARMM_VDW");
		opt.energyTermList.push_back("SCWRL4_HBOND");
		opt.energyTermList.push_back("CHARMM_IMM1");
		opt.energyTermList.push_back("CHARMM_IMM1REF");
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
	
	// run parameters	
	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()) {
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}
	opt.rerunConf = OP.getConfFile();

	return opt;
}