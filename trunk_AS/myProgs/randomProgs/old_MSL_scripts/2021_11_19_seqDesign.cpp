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
#include "BaselineEnergyBuilder.h"
#include "BaselineInteraction.h"
#include "BaselinePairInteraction.h"
#include "BaselineOuterPairInteraction.h"
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"
#include "SasaCalculator.h"
//#include "seqDesign.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "seqDesign";
string programDescription = "Designs sequences for backbone geometries extracted from the PDB, optimizing specifically for vdW energies";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "22 September 2021";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************
 *  
 *  =======  OPTIONS =======
 *
 ******************************************/
struct Options{
	// input files
	string backboneCrd;
	string pdbOutputDir;
	string topFile;
	string parFile;
	string geometryDensityFile;
	string solvFile;
	string hbondFile;
	string rotLibFile;
	string infile;
	string selfEnergyFile;
	string pairEnergyFile;
	string seqEntropyFile;
	string AACompositionPenaltyFile;

	// sequence parameters 
	string sequence;
	int sequenceLength;
	string backboneAA;
	int backboneLength;

	// booleans
	bool getGeoFromPDBData; //TRUE: randomly choose a dimeric geometry from the membrane protein pdb landscape OR FALSE: use a given dimer geometry
	bool verbose; //TRUE: write more outputs throughout the run
	bool deleteTerminalHbonds; //TRUE: delete hydrogen bonds at the termini OR FALSE: keep hydrogen bonds at termini
	bool linkInterfacialPositions; //TRUE: keeps interfacial positions linked when searching for the best states in stateMC (less memory) OR FALSE: unlinks positions (memory intensive)
	bool useSasa; //TRUE: solvent accessible surface area used to choose the number rotamers at each position OR FALSE: give a set number of rotamers to interface vs non-interface
	bool useTimeBasedSeed; //TRUE: use time based seed for RandomNumberGenerator OR FALSE: use given seed
	bool energyLandscape; //TRUE: collect all sequences and their respective monomer and dimer energies
	bool useBaselineComposition; //TRUE: use Rosetta BaselineAAComposition to navigate sequence space or FALSE: only use vdw, hbond, and sequence entropy

	// repack parameters
	int greedyCycles;
	int seed;

	// load rotamers useSasa = false
	string SL; //number of rotamers
	string SLInterface; //number of rotamers for interfacial AAs
	// load rotamers useSasa = true
	std::vector<string> sasaRepackLevel;
	int interfaceLevel;
	
	// tm
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int sequenceStart;

	// Starting Geometry
	double xShift;
	double zShift;
	double crossingAngle;
	double axialRotation;
	bool transform;

	// transformation variables
	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;

	// crossing point
	int thread;
	int bbThread;
	
	// Monte Carlo parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	double repackEnergyCutoff;
	double vdwEnergyCutoff;
	
	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	double weight_seqEntropy;
	
	// alternate identities
	vector<string> Ids;

	// state Monte Carlo Options
	int numStatesToSave;

	// energy terms to output
	vector<string> monomerEnergyTerms;
	vector<string> monomerIMM1EnergyTerms;
	vector<string> dimerEnergyTerms;
	vector<string> energyLandscapeTerms;
	vector<string> energyTermsToOutput;
	vector<string> energyTermList;

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

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
	string runNumber;
	bool useIMM1;

	//SelfPairManager Options
	bool runDEESingles;
	bool runDEEPairs;
	bool runSCMF;
};

/**************************************************
 *  
 *  =======  INTERNAL FUNCTIONS  =======
 *
 **************************************************/
// parse config file for given options
Options parseOptions(int _argc, char * _argv[], Options defaults);

/***********************************
 *help functions
 ***********************************/
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

/***********************************
 *geometry
 ***********************************/
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
void xShiftTransformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, double _xShift, Transforms & _trans);
void readGeometryFile(string _filename, vector<string>& _fileVec);
void getGeometry(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_densities, ofstream &_out);

/***********************************
 *string output
 ***********************************/
// TODO: if possible, make some of these more multipurpose
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum);
string convertPolymerSeqToOneLetterSeq(Chain &_chain);
string generateString(string _backbone, int _length);
string generateBackboneSequence(string _backbone, int _length);
string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum);
string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum);
string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string getInterfaceString(vector<int> _interface, int _seqLength);
string getAlternateIdString(vector<string> _alternateIds); 
string getInterfaceSequence(Options &_opt, string _interface, string _sequence);

/***********************************
 *repack functions
 ***********************************/
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
std::vector < std::vector < bool > > getActiveMask (System &_sys);

/***********************************
 *define interface and rotamer sampling
 ***********************************/
vector<int> getRotamerSampling(string _rotamerLevels);
vector<int> getLinkedPositions(vector<int> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel);
vector<uint> getVariablePositions(vector<int> &_interfacialPositions);
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<int> &_interfacialPositions, int _backboneLength);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq);
void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_variablePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis);

/***********************************
 *output file functions
 ***********************************/
void setupDesignDirectory(Options &_opt, string _date);
void printConfigFile(Options & _opt, ofstream & _out);
void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns);
void makeRepackConfig(Options &_opt, string _sequence, string _designDir, string _designNumber, string _pdbPath, string _crdPath, map<string,double> _energyMap);
void makeDockingConfig(Options &_opt, string _sequence, vector<uint> _state, string _designNumber, string _pdbPath, map<string,double> _energyMap, vector<int> _rotamerSampling);
void outputRepackFile(Options &_opt, vector<string> _dockingDesigns);
void outputDesignFiles(Options &_opt, string _interface, vector<int> _rotamerSampling, vector<pair<string,vector<uint>>> _sequenceStatePair, map<string,map<string,double>> _sequenceEnergyMap, vector<double> _densities);

/***********************************
 *load rotamer functions
 ***********************************/
void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot);
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling);
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<int> _interface);
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling);//Uses rotamer sampling defined by SASA values to load rotamers by position

/***********************************
 *baseline energy helper functions
 ***********************************/
vector<double> calcBaselineEnergies(System &_sys, int _seqLength);
vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength);
double sumEnergyVector(vector<double> _energies);

/***********************************
 *calculate energies
 ***********************************/
void computeDimerEnergy(System &_sys, Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_sequence, vector<uint> &_stateVec, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, int _seqNumber, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err);
void computeDimerEnergiesLinked(System &_sys, Options &_opt, map<string,map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err);
void computeDimerEnergies(System &_sys, Options &_opt, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> _rotamerSamplingPerPosition, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
void computeMonomerEnergyNoIMM1(Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);

/***********************************
 *other helper functions
 ***********************************/
void saveEnergyDifference(Options _opt, map<string,map<string,double>> &_sequenceEnergyMap, string _sequence);
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap, vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1);
void outputEnergiesByTermLinked(EnergySet *_Eset, map<string,double> &_energyMap, vector<string> _energyTermList, string _energyDescriptor);
void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt);
map<string, double> readSingleParameters(string _baselineFile);
map<string,map<string,map<uint, double>>> readPairParameters(string _baselineFile);

/***********************************
 *energy builders
 ***********************************/
void buildSelfInteractions(System &_sys, map<string, double> &_selfMap);
void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap);
void buildSequenceEntropy(System &_sys, map<string, double> &_sequenceEntropyMap, double _weight);

/***********************************
 *stateMC helper functions
 ***********************************/
//gets a random position and chooses a random rotamer
void randomRotamerChange(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _variablePositions, vector<unsigned int> &_stateVec);
void randomRotamerChangeNonLinked(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, map<int,map<int,pair<uint,uint>>> &_posRotLimitMap, vector<uint> _variablePositions, vector<unsigned int> &_stateVec);
void sameSequenceChecker(string &_newSeq, vector<string> &_seqs);
bool sameSequenceChecker(string &_newSeq, double &_newEnergy, vector<uint> &_state, vector<pair<double,string>> &_enerSeqPair, vector<pair<double,vector<uint>>> &_energyStateVec);
void saveSequence(Options &_opt, vector<pair<double,string>> &_energyVector, vector<pair<double,vector<uint>>> &_energyStateVec, string _sequence, vector<uint> _state, double _energy);
map<int,map<int,pair<uint, uint>>> setupRotamerPositionMap(System &_sys, vector<uint> _interfacialPositionsList);
void unlinkBestState(Options &_opt, vector<uint> &_bestState, vector<int> _linkedPositions, int _backboneLength);
bool convertStateMapToSequenceMap(System &_sys, vector<pair<double,vector<uint>>> &_energyStateVec, map<vector<uint>, map<string,double>> &_stateEnergyMap, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair);

/***********************************
 *sequence entropy functions
 ***********************************/
map<string,int> getAACountMap(vector<string> _seq);
double calcNumberOfPermutations(map<string,int> _seqAACounts, int _seqLength);
void internalAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, int _seqLength);
void sequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, int _seqLength);
double calculateSequenceProbability(map<string,int> &_seqCountMap, map<string,double> &_entropyMap, double _numberOfPermutations);
void calculateSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq, map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_internalSEProb, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal);

// other functions
double getStandardNormal(RandomNumberGenerator& RNG);
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);
void addSequencesToVector(vector<pair<double,string>> &_energyVector, vector<string> &_allSeqs);

// Linked version of the state monte carlo
void stateMCLinked(System &_sys, SelfPairManager &_spm, Options &_opt, PolymerSequence &_PS, map<string, map<string,double>> &_sequenceEnergyMap, vector<unsigned int> &_bestState, vector<string> &_seqs, vector<string> &_allSeqs, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	//initialize baseline map for sequence entropy calculation: uses entropy to give an energy score for the likelihood that a sequence is found within the membrane
	map<string, double> seqEntMap = readSingleParameters(_opt.seqEntropyFile);

	//initialize baselineAAComposition energy map (based on Rosetta baselineAAComposition to prevent unlikely sequences by adding energy (ex. if more than 2SER in sequence, add 1000 energy score for each additional PHE)
	//if (_opt.useBaselineComposition){
	//	BaselineAAComposition bac(&_sys);
	//	bac.readPenaltyFile(_opt.AACompositionPenaltyFile);
	//	Eset->addInteraction(new BaselineAAComposition(bac));
	//}

	/******************************************************************************
	 *                   === SETUP FOR STATE MONTE CARLO ===
	 ******************************************************************************/
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);
	
	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, 50);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _spm.getStateEnergy(_bestState);

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	vector<unsigned int> currStateVec = _bestState;
	MC.setEner(bestEnergy);
	
	// initialize map for accepting energies
	map<string,double> stateMCEnergies;
	map<vector<uint>, map<string,double>> stateEnergyMap;

	// initialize energy variables for the MonteCarlo
	double prevStateSEProb = 0;
	double totEnergy = 0;

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);

	// initialize energy vectors
	vector<pair<double,string>> energyVector;
	vector<pair<double,vector<uint>>> energyStateVec;

	// Setup EnergyLandscape file
	ofstream lout;
	string loutfile  = _opt.pdbOutputDir + "/stateMCEnergyLandscape.out";
	lout.open(loutfile.c_str());
	lout << "***STARTING GEOMETRY:***" << endl;
	lout << "xShift: " << _opt.xShift << endl;
	lout << "crossingAngle: " << _opt.crossingAngle << endl;
	lout << "axialRotation: " << _opt.axialRotation << endl;
	lout << "zShift: " << _opt.zShift << endl << endl;
	lout << "Number of MCCycles: " << _opt.MCCycles << endl;
	lout << "Design Number\tTotal Cycles\tSequence\tTotal\tDimer\tBaseline\tVDWDimer\tHBONDDimer\tEnergyw/seqEntropy\tCycle" << endl;

	/******************************************************************************
	 *                      === BEGIN STATE MONTE CARLO ===
	 ******************************************************************************/
	cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	_sout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		
		// Get energy term and probability for the first sequence (these terms can then get replaced by future sequences that are accepted by MC)
		if (cycleCounter == 0){
			outputEnergiesByTerm(_spm, prevStateVec, stateMCEnergies, _opt.energyTermList, "DimerNoIMM1", false);
			stateMCEnergies["DimerNoIMM1"] = bestEnergy;

			stateEnergyMap[prevStateVec] = stateMCEnergies;
			totEnergy = stateMCEnergies["EnergyBeforeLocalMC"];
			sequences[totEnergy] = prevStateSeq;
			double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
			saveSequence(_opt, energyVector, energyStateVec, prevStateSeq, prevStateVec, prevVDW);
		}
		// Reset the energy map to save energies from new state after changing the rotamer
		stateMCEnergies.clear();
		randomRotamerChange(_sys, _opt, _RNG, _interfacialPositionsList, currStateVec);

		// Set a mask and run a greedy to get the best state for the sequence
		_sys.setActiveRotamers(currStateVec);
		vector<vector<bool>> mask = getActiveMask(_sys);
		_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
		currStateVec = _spm.getMinStates()[0];

		// Get the sequence for the random state
		_sys.setActiveRotamers(currStateVec);
		string currStateSeq = convertPolymerSeqToOneLetterSeq(chain);

		// Compute dimer energy
		outputEnergiesByTerm(_spm, currStateVec, stateMCEnergies, _opt.energyTermList, "DimerNoIMM1", false);
		double currStateEnergy = _spm.getStateEnergy(currStateVec);
		stateMCEnergies["DimerNoIMM1"] = currStateEnergy;
		
		// Convert the energy term (which actually saves the probability of the sequence in the whole _system)
		// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)
		double bestEnergyTotal;
		double currEnergyTotal;
		double currStateSEProb;
		double internalSEProb;
		calculateSequenceEntropy(_opt, prevStateSeq, currStateSeq, seqEntMap, prevStateSEProb, currStateSEProb, internalSEProb, bestEnergy, currStateEnergy, bestEnergyTotal, currEnergyTotal);
		MC.setEner(bestEnergyTotal);

		// MC accept and reject conditions
		double currVDW = _spm.getStateEnergy(currStateVec, "CHARMM_VDW");
		double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
		if (!MC.accept(currEnergyTotal) && currVDW > prevVDW+5){
			_sys.setActiveRotamers(prevStateVec);
			currStateVec = prevStateVec;

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			double currStateEnergyNoComposition = currStateEnergy;
			if (_opt.useBaselineComposition){
				double compositionEnergy = _spm.getStateEnergy(currStateVec, "BASELINE_COMPOSITION");
				currStateEnergyNoComposition = currStateEnergy-compositionEnergy;
			}
			saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currVDW);

			//TODO: make these a separate function or put in comments  for them
			bestEnergy = currStateEnergy;
			MC.setEner(currEnergyTotal);
			prevStateSEProb = currStateSEProb;
			prevStateSeq = currStateSeq;
			prevStateVec = currStateVec;
			_sys.setActiveRotamers(currStateVec);
			stateMCEnergies["EnergyBeforeLocalMC"] = currStateEnergyNoComposition+stateMCEnergies["BASELINE"]-stateMCEnergies["BASELINE_PAIR"];
			stateMCEnergies["EnergyBeforeLocalMCw/seqEntropy"] = bestEnergyTotal-currEnergyTotal;
			stateMCEnergies["SequenceProbability"] = currStateSEProb;
			stateMCEnergies["InternalSequenceProbability"] = internalSEProb;
			stateEnergyMap[currStateVec] = stateMCEnergies;

			if (_opt.energyLandscape){
				map<string,double> energyMap = stateEnergyMap.at(currStateVec);
				lout << _opt.runNumber << "\t" << _opt.MCCycles << "\t" << currStateSeq << "\t";
				for (uint j=0; j<_opt.energyLandscapeTerms.size(); j++){
					lout << energyMap.at(_opt.energyLandscapeTerms[j]) << "\t";
				}
				lout << cycleCounter << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currStateSeq << "; PrevE=  " << bestEnergy << " : CurrE= " << currStateEnergy << "; PrevVDW: " << prevVDW << " : CurrVDW: " << currVDW << endl;
			}
		cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to 
		if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
			MC.reset(3649, 3649, 100, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		}
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();

	addSequencesToVector(energyVector, _allSeqs);
	convertStateMapToSequenceMap(_sys, energyStateVec, stateEnergyMap, _sequenceEnergyMap, _sequenceStatePair);
	//bool sequencesClashing = checkForClashing(_allSeqs,  _sequenceEnergyMap){

	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_sout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_seqs = _allSeqs;
}

void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS, map<string, map<string,double>> &_sequenceEnergyMap, vector<unsigned int> &_bestState, vector<string> &_seqs, vector<string> &_allSeqs, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setBuildNonBondedInteractions(false);
	if(!CSB.buildSystem(_PS)) {
		cerr << "Unable to build system from " << _PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	Chain & chain = sys.getChain("A");

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(_sys.getAtomPointers(),false);
	sys.buildAllAtoms();
	CSB.updateNonBonded(10,12,50);

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct
	
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
	Eset->setTermActive("SCWRL4_HBOND", true);
	
	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
	map<string, double> selfMap = readSingleParameters(_opt.selfEnergyFile);
	map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_opt.pairEnergyFile);
	buildSelfInteractions(sys, selfMap);
	buildPairInteractions(sys, pairMap);

	//initialize baseline map for sequence entropy calculation: uses entropy to give an energy score for the likelihood that a sequence is found within the membrane
	map<string, double> seqEntMap = readSingleParameters(_opt.seqEntropyFile);

	//initialize baselineAAComposition energy map (based on Rosetta baselineAAComposition to prevent unlikely sequences by adding energy (ex. if more than 2SER in sequence, add 1000 energy score for each additional PHE)
	if (_opt.useBaselineComposition){
		BaselineAAComposition bac(&sys);
		bac.readPenaltyFile(_opt.AACompositionPenaltyFile);
		Eset->addInteraction(new BaselineAAComposition(bac));
	}

	/******************************************************************************
	 *              === LOAD ROTAMERS AND CHOOSE TO LINK INTERFACE ===
	 ******************************************************************************/
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);
	CSB.updateNonBonded(10,12,50);

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	sys.buildAllAtoms();
	
	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	//spm.updateWeights();
	spm.setOnTheFly(false);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();//TODO: length 25 seems to break here or a little further. Likely takes up too much memory on server?  I know I shouldn't be recalculating...but it's the only way it's worked so far. Find a way to make it work without recalculationg, likely would have to unlink positions (searched but couldn't find) and recalculate; or update the original spm with a smaller system so it takes up less memory? That's the only thing I can think of for here unless it's in the spm code somewhere; WHENEVER I WROTE THIS ABOUT MEMORY I WAS RIGHT...I just forgot about using top to evaluate

	//Setup the map for limiting the number of rotamers per position (allows me to keep AAs the same without needing the rotamers to be symmetric)
	map<int,map<int,pair<uint, uint>>> posRotLimitMap;
	posRotLimitMap = setupRotamerPositionMap(sys, _interfacialPositionsList);

	/******************************************************************************
	 *                   === SETUP FOR STATE MONTE CARLO ===
	 ******************************************************************************/
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);
	
	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, 50);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	sys.setActiveRotamers(_bestState);
	double bestEnergy = spm.getStateEnergy(_bestState);

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	vector<unsigned int> currStateVec = _bestState;
	MC.setEner(bestEnergy);
	
	// initialize map for accepting energies
	map<string,double> stateMCEnergies;
	map<vector<uint>, map<string,double>> stateEnergyMap;

	// initialize energy variables for the MonteCarlo
	double prevStateSEProb = 0;
	double totEnergy = 0;

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	//Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);

	// initialize energy vectors
	vector<pair<double,string>> energyVector;
	vector<pair<double,vector<uint>>> energyStateVec;

	// Setup EnergyLandscape file
	ofstream lout;
	string loutfile  = _opt.pdbOutputDir + "/stateMCEnergyLandscape.out";
	lout.open(loutfile.c_str());
	lout << "***STARTING GEOMETRY:***" << endl;
	lout << "xShift: " << _opt.xShift << endl;
	lout << "crossingAngle: " << _opt.crossingAngle << endl;
	lout << "axialRotation: " << _opt.axialRotation << endl;
	lout << "zShift: " << _opt.zShift << endl << endl;
	lout << "Number of MCCycles: " << _opt.MCCycles << endl;
	lout << "Design Number\tTotal Cycles\tSequence\tTotal\tDimer\tBaseline\tVDWDimer\tHBONDDimer\tEnergyw/seqEntropy\tCycle" << endl;

	/******************************************************************************
	 *                      === BEGIN STATE MONTE CARLO ===
	 ******************************************************************************/
	cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	_sout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		
		// Get energy term and probability for the first sequence (these terms can then get replaced by future sequences that are accepted by MC)
		if (cycleCounter == 0){
			outputEnergiesByTerm(spm, prevStateVec, stateMCEnergies, _opt.energyTermList, "DimerNoIMM1", false);
			stateMCEnergies["DimerNoIMM1"] = bestEnergy;

			stateEnergyMap[prevStateVec] = stateMCEnergies;
			totEnergy = stateMCEnergies["EnergyBeforeLocalMC"];
			sequences[totEnergy] = prevStateSeq;
			double prevVDW = spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
			saveSequence(_opt, energyVector, energyStateVec, prevStateSeq, prevStateVec, prevVDW);
		}
		// Reset the energy map to save energies from new state after changing the rotamer
		stateMCEnergies.clear();
		randomRotamerChangeNonLinked(sys, _opt, _RNG, posRotLimitMap, _interfacialPositionsList, currStateVec);//this was previously _sys...could that have been screwing up the rotamer changes...? looking at the function it's likely a no, but still glad I caught it

		// Set a mask and run a greedy to get the best state for the sequence
		sys.setActiveRotamers(currStateVec);
		vector<vector<bool>> mask = getActiveMask(sys);
		spm.runGreedyOptimizer(_opt.greedyCycles, mask);
		currStateVec = spm.getMinStates()[0];

		// Get the sequence for the random state
		sys.setActiveRotamers(currStateVec);
		string currStateSeq = convertPolymerSeqToOneLetterSeq(chain);

		// Compute dimer energy
		outputEnergiesByTerm(spm, currStateVec, stateMCEnergies, _opt.energyTermList, "DimerNoIMM1", false);
		double currStateEnergy = spm.getStateEnergy(currStateVec);
		stateMCEnergies["DimerNoIMM1"] = currStateEnergy;
		
		// Convert the energy term (which actually saves the probability of the sequence in the whole system)
		// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)
		double bestEnergyTotal;
		double currEnergyTotal;
		double currStateSEProb;
		double internalSEProb;
		calculateSequenceEntropy(_opt, prevStateSeq, currStateSeq, seqEntMap, prevStateSEProb, currStateSEProb, internalSEProb, bestEnergy, currStateEnergy, bestEnergyTotal, currEnergyTotal);
		MC.setEner(bestEnergyTotal);

		// MC accept and reject conditions
		double currVDW = spm.getStateEnergy(currStateVec, "CHARMM_VDW");
		double prevVDW = spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
		if (!MC.accept(currEnergyTotal) && currVDW > prevVDW+5){
			sys.setActiveRotamers(prevStateVec);
			currStateVec = prevStateVec;

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			double currStateEnergyNoComposition = currStateEnergy;
			if (_opt.useBaselineComposition){
				double compositionEnergy = spm.getStateEnergy(currStateVec, "BASELINE_COMPOSITION");
				currStateEnergyNoComposition = currStateEnergy-compositionEnergy;
			}
			saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currVDW);

			//TODO: make these a separate function or put in comments  for them
			bestEnergy = currStateEnergy;
			MC.setEner(currEnergyTotal);
			prevStateSEProb = currStateSEProb;
			prevStateSeq = currStateSeq;
			prevStateVec = currStateVec;
			sys.setActiveRotamers(currStateVec);
			stateMCEnergies["EnergyBeforeLocalMC"] = currStateEnergyNoComposition+stateMCEnergies["BASELINE"]-stateMCEnergies["BASELINE_PAIR"];
			stateMCEnergies["EnergyBeforeLocalMCw/seqEntropy"] = bestEnergyTotal-currEnergyTotal;
			stateMCEnergies["SequenceProbability"] = currStateSEProb;
			stateMCEnergies["InternalSequenceProbability"] = internalSEProb;
			stateEnergyMap[currStateVec] = stateMCEnergies;

			if (_opt.energyLandscape){
				map<string,double> energyMap = stateEnergyMap.at(currStateVec);
				lout << _opt.runNumber << "\t" << _opt.MCCycles << "\t" << currStateSeq << "\t";
				for (uint j=0; j<_opt.energyLandscapeTerms.size(); j++){
					lout << energyMap.at(_opt.energyLandscapeTerms[j]) << "\t";
				}
				lout << cycleCounter << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currStateSeq << "; PrevE=  " << bestEnergy << " : CurrE= " << currStateEnergy << "; PrevVDW: " << prevVDW << " : CurrVDW: " << currVDW << endl;
			}
		cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to 
		if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
			MC.reset(3649, 3649, 100, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
			//cout << "MonteCarlo Reset!" << endl;
		}

	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();

	addSequencesToVector(energyVector, _allSeqs);
	convertStateMapToSequenceMap(sys, energyStateVec, stateEnergyMap, _sequenceEnergyMap, _sequenceStatePair);
	//bool sequencesClashing = checkForClashing(_allSeqs,  _sequenceEnergyMap){

	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_sout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_seqs = _allSeqs;
}

vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=4; k<_opt.backboneLength-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

/******************************************
 *  
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
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
	ofstream sout;
	ofstream err;
	ofstream rerun;

	setupDesignDirectory(opt, date);

	string soutfile = opt.pdbOutputDir + "/summary.out";
	string errfile  = opt.pdbOutputDir + "/errors.out";
	string rerunfile = opt.pdbOutputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	sout << date << endl;
	err << date << endl;
	rerun << "#" << date << endl;

	//TODO: if I decide to use the same sequence geometry and no run SCMF, put this at the end of the code instead
	printConfigFile(opt, rerun);
	
	/******************************************************************************
	 *               === LOAD RANDOM GEOMETRY FROM GEOMETRY FILE ===
	 ******************************************************************************/
	RandomNumberGenerator RNG;
	if (opt.useTimeBasedSeed){
		RNG.setTimeBasedSeed();
	} else {
		RNG.setSeed(opt.seed); 
	}

	vector<double> densities;
	if (opt.getGeoFromPDBData){
		getGeometry(opt, RNG, densities, sout);
	}

	// Output Geometry
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << opt.xShift << "\tDensity: " << densities[0] << endl;
	cout << "crossingAngle: " << opt.crossingAngle << "\tDensity: " << densities[0] << endl;
	cout << "axialRotation: " << opt.axialRotation << "\tDensity: " << densities[1] << endl;
	cout << "zShift:        " << opt.zShift << "\tDensity: " << densities[2] << endl << endl;

	//String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.Ids);
	cout << "Amino acids for design: LEU " << alternateIds << endl;

	/******************************************************************************
	 *                         === GENERATE POLYGLY ===
	 ******************************************************************************/
	string polySeq = generatePolymerSequence(opt.backboneAA, opt.backboneLength, opt.thread);
	PolymerSequence PS(polySeq);

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
		err << "Unable to read axis" << endl;
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
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// Variables to output from defineInterfaceAndRotamerSampling function
	string rotamerLevels;
	string variablePositionString;
	string rotamerSamplingString;
	vector<int> linkedPositions;
	vector<uint> interfacePositions;
	vector<int> rotamerSamplingPerPosition;

	// Defines the interfacial positions and the number of rotamers to give each position
	defineInterfaceAndRotamerSampling(opt, PS, rotamerLevels, polySeq, variablePositionString, rotamerSamplingString, linkedPositions, interfacePositions, rotamerSamplingPerPosition, sout, axis);

	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(opt.infile);
	
	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = pdb.getChain("A");
	Chain & chainB = pdb.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// Set up object used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Transform to chosen geometry 
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile, opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	CSB.setBuildNonBondedInteractions(false);
	
	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	PolymerSequence PL(polySeq);
	if(!CSB.buildSystem(PL)) {
		err << "Unable to build system from " << polySeq << endl;
		exit(0);
	}

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();
	
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
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
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1REF", 0);
	Eset->setWeight("CHARMM_IMM1", 0);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,opt);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	checkIfAtomsAreBuilt(sys, err);

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
	map<string, double> selfMap = readSingleParameters(opt.selfEnergyFile);
	map<string,map<string,map<uint,double>>> pairMap = readPairParameters(opt.pairEnergyFile);
	buildSelfInteractions(sys, selfMap);
	buildPairInteractions(sys, pairMap);

	sys.calcEnergy();
	vector<vector<string>> linkedPos = convertToLinkedFormat(sys, linkedPositions, opt.backboneLength);
	sys.setLinkedPositions(linkedPos);

	//make into a function
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs
	Eset->setTermActive("CHARMM_IMM1REF", false);
	Eset->setTermActive("CHARMM_IMM1", false);
	sys.calcEnergy();

	/******************************************************************************
	 *                        === SETUP SPM AND RUN SCMF ===
	 ******************************************************************************/
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.setRunDEE(opt.runDEESingles, opt.runDEEPairs);
	spm.setOnTheFly(false);
	spm.setMCOptions(5000, 0.5, 50000, 3, 5000, 5000, 0.01);//changed to sigmoid and added up to 5000
	spm.saveEnergiesByTerm(true); //added back in on 09_21_2021 to get the vdw and hbond energies 
	spm.calculateEnergies();

	//Setup running SCMF or UnbiasedMC
	if (opt.runSCMF == true){
		cout << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	} else {
		cout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	time(&spmStart);
	spm.runOptimizer();
	time(&spmEnd);
	
	spmTime = difftime (spmEnd, spmStart);
	vector<unsigned int> initialState = spm.getSCMFstate();
	vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
	double bestEnergy = spm.getStateEnergy(bestState);
	double hbondEnergy = spm.getStateEnergy(bestState, "SCWRL4_HBOND");
	double vdwEnergy = spm.getStateEnergy(bestState, "CHARMM_VDW");

	sout << endl << "End SelfPairManager Optimization: " << spmTime << "s" << endl;
	cout << endl << "End SelfPairManager Optimization: " << spmTime << "s" << endl;
	
	sys.setActiveRotamers(initialState);
	string initialSeq = convertPolymerSeqToOneLetterSeq(sys.getChain("A"));
	sys.setActiveRotamers(bestState);
	string SCMFBestSeq = convertPolymerSeqToOneLetterSeq(sys.getChain("A"));
	sys.calcEnergy();
	opt.sequence = SCMFBestSeq;//used for outputting starting sequence

	// TODO: make into a function?
	string seqInterface = getInterfaceSequence(opt, rotamerSamplingString, opt.sequence);
	int numInterfacialPositions = linkedPos.size();
	sout << "                   " << "Sequence\tTotalEnergy\tVDW\tHBOND\tInterface\tInterfacialPositions\tSASARepackLevels\tInterfaceLevels" << endl;
	sout << "Starting Sequence: " << opt.sequence << "\t" << bestEnergy << "\t" << vdwEnergy << "\t" << hbondEnergy << "\t" << seqInterface << "\t" << numInterfacialPositions << "\t" << opt.sasaRepackLevel.size() << "\t" << opt.interfaceLevel << "\t" << initialSeq << "\t" << SCMFBestSeq << endl;
	sout << "Interface:         " << variablePositionString << endl;
	sout << "Rotamer Levels:    " << rotamerLevels << endl;
	
	/******************************************************************************
	 *           === METHODS FOR DETERMINING ALTERNATE SEQUENCES ===
	 ******************************************************************************/
	vector<string> seqs;
	vector<string> allSeqs;

	//Initialize energyMap to hold all energies for output into a summary file
	map<string, map<string,double>> sequenceEnergyMap;

	// Initialize sequence and state pair vector: each sequence will be tied to it's state with the proper rotamers
	vector<pair<string,vector<uint>>> sequenceStatePair;

	/******************************************************************************
	 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
	 ******************************************************************************/
	// Unlink the best state from SCMF if not using linked positions during the state Monte Carlo
	if (!opt.linkInterfacialPositions){
		unlinkBestState(opt, bestState, rotamerSamplingPerPosition, opt.backboneLength);
		stateMCUnlinked(sys, opt, PL, sequenceEnergyMap, bestState, seqs, allSeqs, sequenceStatePair, interfacePositions, rotamerSamplingPerPosition, linkedPos, RNG, sout, err);
	} else {
		stateMCLinked(sys, spm, opt, PL, sequenceEnergyMap, bestState, seqs, allSeqs, sequenceStatePair, interfacePositions, rotamerSamplingPerPosition, linkedPos, RNG, sout, err);
	}

	/******************************************************************************
	 *            === CALCULATE MONOMER ENERGIES OF EACH SEQUENCE ===
	 ******************************************************************************/
	computeMonomerEnergies(opt, trans, sequenceEnergyMap, seqs, RNG, sout, err);

	/******************************************************************************
	 *               === LOCAL REPACKS ON EACH UNIQUE SEQUENCE ===
	 ******************************************************************************/
	computeDimerEnergies(sys, opt, sequenceEnergyMap, sequenceStatePair, rotamerSamplingPerPosition, linkedPos, RNG, sout, err);

	/******************************************************************************
	 *                   === WRITE OUT ENERGY AND DESIGN FILES ===
	 ******************************************************************************/
	outputDesignFiles(opt, rotamerSamplingString, rotamerSamplingPerPosition, sequenceStatePair, sequenceEnergyMap, densities);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

	err.close();
	sout.close();
}

/*********************************
 *  
 *  ======= FUNCTIONS =======
 *
 *********************************/

/***********************************
 *help functions
 ***********************************/
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputErrorMessage(Options &_opt){
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << _opt.errorMessages << endl;
		cerr << endl;
		cerr << _opt.OPerrors << endl;
		usage();
}

void help(Options defaults) {
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
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "pdbOutputDir " << defaults.pdbOutputDir << endl;

	cout << setw(20) << "tmStart " << defaults.tmStart << endl;
	cout << setw(20) << "tmEnd " << defaults.tmEnd << endl;

	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "geometryDensityFile " << defaults.geometryDensityFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;
	cout << setw(20) << "infile " << defaults.infile << endl;
	cout << setw(20) << "selfEnergyFile " << defaults.selfEnergyFile << endl;
	cout << setw(20) << "pairEnergyFile " << defaults.pairEnergyFile << endl;
	cout << setw(20) << "AACompositionPenaltyFile " << defaults.AACompositionPenaltyFile << endl << endl;

	cout << "#Geometry and Transformation parameters" << endl;
	cout << setw(20) << "xShift" << defaults.xShift << endl;
	cout << setw(20) << "crossingAngle" << defaults.crossingAngle << endl;
	cout << setw(20) << "axialRotation" << defaults.axialRotation << endl;
	cout << setw(20) << "zShift" << defaults.zShift << endl;
	cout << setw(20) << "thread" << defaults.thread << endl;
	cout << setw(20) << "backboneLength " << defaults.backboneLength << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;
	cout << setw(20) << "deleteTerminalHbonds" << defaults.deleteTerminalHbonds << endl;
	cout << setw(20) << "useSasa" << defaults.useSasa << endl;
	cout << setw(20) << "getGeoFromPDBData" << false << endl;//Since we already have the geometry output here, default to false in the rerun config
	cout << setw(20) << "runDEESingles" << defaults.runDEESingles << endl;
	cout << setw(20) << "runDEEPairs" << defaults.runDEEPairs << endl;
	cout << setw(20) << "runSCMF" << defaults.runSCMF << endl;
	//TODO: set this up so that instead of running through that, I just have the output sequence and state below
	if (defaults.runDEESingles == false && defaults.runDEEPairs == false && defaults.runSCMF == false){

	}

	if (defaults.useSasa == true){
		//cout << "#Load Rotamers based on SASA scores" << endl;
		for (uint i=0; i<defaults.sasaRepackLevel.size()-1; i++){
			cout << setw(20) << "sasaRepackLevel" << defaults.sasaRepackLevel[i] << endl;
		}
		cout << setw(20) << "interfaceLevel" << defaults.interfaceLevel << endl;
	} else {
		//cout << "#Load Rotamers by interface and non interfacial positions" << endl;
		cout << setw(20) << "SL" << defaults.SL << endl;
		cout << setw(20) << "SLInterface" << defaults.SLInterface << endl;
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

/***********************************
 *geometry
 ***********************************/
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
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by

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

void xShiftTransformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, double _xShift, Transforms & _trans) {
	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void readGeometryFile(string _filename, vector<string>& _fileVec) {
	ifstream file;
	file.open(_filename.c_str()); 
	if(!file.is_open()) {
		cerr << "Unable to open " << _filename << endl;
		exit(0);
	}

	string parameterList;

	while(file) {
		getline(file, parameterList);
		_fileVec.push_back(parameterList);
	}
	file.close();
}

void getGeometry(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_densities, ofstream &_out){
	// Setup file reader
	Reader reader(_opt.geometryDensityFile);
	reader.open();
	if(!(reader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.geometryDensityFile << endl;
		exit(0);
	}
	vector<string> lines = reader.getAllLines();

	// Extract the geometries from a random line of the geometry file
	int geometryLine = _RNG.getRandomInt(1,lines.size()-1);
	vector<string> tokens = MslTools::tokenize(lines[geometryLine], "\t");//xShift, crossingAngle, axialRotation, zShift, angleDistDensity, axialRotationDensity, zShiftDensity
	_opt.xShift = MslTools::toDouble(tokens[0]);
	_opt.crossingAngle = MslTools::toDouble(tokens[1]);
	_opt.axialRotation = MslTools::toDouble(tokens[2]);
	_opt.zShift = MslTools::toDouble(tokens[3]);
	double angleDistDensity = MslTools::toDouble(tokens[4]);
	double axialRotationDensity = MslTools::toDouble(tokens[5]);
	double zShiftDensity = MslTools::toDouble(tokens[6]);
	_densities.push_back(angleDistDensity);
	_densities.push_back(axialRotationDensity);
	_densities.push_back(zShiftDensity);

	// Output to summary file
	_out << "***STARTING GEOMETRY:***" << endl;
	_out << "xShift:        " << _opt.xShift << "\tDensity: " << angleDistDensity << endl;
	_out << "crossingAngle: " << _opt.crossingAngle << "\tDensity: " << angleDistDensity << endl;
	_out << "axialRotation: " << _opt.axialRotation << "\tDensity: " << axialRotationDensity << endl;
	_out << "zShift:        " << _opt.zShift << "\tDensity: " << zShiftDensity << endl << endl;

	_out << "Geometry: " << _opt.xShift << "\t" << _opt.crossingAngle << "\tDensity: " << angleDistDensity << endl;
}

/***********************************
 *string output functions
 ***********************************/
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

string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
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
	return "A" + ps;
}

string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getCurrentIdentity().getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}

string generateString(string _backbone, int _length) {
	string str = "";
	for (uint i=0; i<_length; i++){
		str = str + _backbone;
	}
	return str;
}

string generateBackboneSequence(string _backbone, int _length) {
	string str = "";
	//2021-09-21: add in an alanine cap to allow for more variable positions at the leucine region
	for (uint i=0; i<_length; i++){
		//if (i<4 || i>_length-5){
		if (i<4){
			str = str + "L";
		}
		else if (i==_length-3 || i==_length-1){
		//else if (i==_length-7 || i==_length-5){
			str = str + "I";
		} else {
		      str = str + _backbone;
		}
	}
	return str;
}

string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
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
			counter++;
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (_interfacialPositions[counter] == 1){
				ps = ps + " [";
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
				for (uint i=0; i<_alternateIds.size(); i++){
					if(_alternateIds[i] == "HIS") {
						ps = ps + " HSE";
					} else {
						ps = ps + " " + _alternateIds[i];
					}
				}
				ps = ps + "] ";
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
			}
			counter++;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}

string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum) {
	string ps = "";
	string resName = MslTools::getThreeLetterCode(_backboneAA);
	if(resName == "HIS") {
		resName = "HSE";
	}
	for (uint i=0; i<_backboneLength; i++){
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum) {
	string ps = "";
	for (uint i=0; i<_sequence.length(); i++){
		stringstream tmp;
		tmp << _sequence[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		if(resName == "HIS") {
			resName = "HSE";
		}
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}

string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	int startPos = _startResNum;
	int endPos = _startResNum+_seq.length();
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		int pos = it-_seq.begin()+_startResNum;
		if (it == _seq.begin() || it == _seq.end()-1){
		//if (it == _seq.begin()){
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
			counter++;
		} else if (pos < startPos+4 || pos > endPos-5){
		//} else if (pos < startPos+4 || pos > endPos-9){//2021-09-21: accounts for adding in Ala at the ends of the sequence; doesn't allow LILI to change
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			//cout << pos << endl;
			if (find(_interfacialPositions.begin(), _interfacialPositions.end(), pos) != _interfacialPositions.end()){	
				ps = ps + " [";
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
				for (uint i=0; i<_alternateIds.size(); i++){
					if(_alternateIds[i] == "HIS") {
						ps = ps + " HSE";
					} else {
						ps = ps + " " + _alternateIds[i];
					}
				}
				ps = ps + "] ";
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
			}
			counter++;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

string getInterfaceString(vector<int> _interface, int _seqLength){
	string interfaceString = "";
	for (uint i=0; i<_interface.size(); i++){
		if (i == _seqLength){
			i = _interface.size();
		} else {
			interfaceString += MslTools::intToString(_interface[i]);
		}
	}
	return interfaceString;
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

string getInterfaceSequence(Options &_opt, string _interface, string _sequence){
	string interfaceSequence = "";
	for(string::iterator it = _interface.begin(); it != _interface.end(); it++) {
		stringstream ss;
		ss << *it;
		int pos = it-_interface.begin();
		int tmp = MslTools::toInt(ss.str());
		if (tmp > _opt.interfaceLevel-1){//interfaceLevel counts from 1 but rotamerLevel coutns from 0
			interfaceSequence = interfaceSequence + "-";
		} else {
			interfaceSequence = interfaceSequence + _sequence[pos];
		}
	}
	return interfaceSequence;
}

/***********************************
 *repack functions
 ***********************************/
void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
}

//Code Samson made a while back that should get each active ID and set a mask for anything that isn't active
std::vector < std::vector < bool > > getActiveMask (System &_sys) {
	_sys.updateVariablePositions();
	std::vector <unsigned int> residueState;
	std::vector < std::vector<unsigned int> > resRots(_sys.getMasterPositions().size());
	std::vector < std::vector<bool> > resMask(_sys.getMasterPositions().size());
	//Initialize residue state at the current active identity for each position
	for (unsigned int i = 0; i < _sys.getMasterPositions().size(); i++) {
		Position &pos = _sys.getPosition(_sys.getMasterPositions()[i]);
		unsigned int activeRes = pos.getActiveIdentity();
		residueState.push_back(activeRes);
		
		resRots[i] = std::vector<unsigned int> (pos.identitySize());
		for (unsigned int j = 0; j < pos.identitySize(); j++) {
			resRots[i][j] = pos.getTotalNumberOfRotamers(j);
		}
	}

	for (unsigned int i = 0; i < residueState.size(); i++) {
		unsigned int activeResidue = residueState[i];
		if (activeResidue >= resRots[i].size()) {
			cerr << "ERROR: the current residue number exceeds the number of residues for position " << i << endl;
			exit(100);
		}
		for (unsigned int j = 0; j < resRots[i].size(); j++) {
			if (j==activeResidue) {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(true);
				}
			} else {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(false);
				}
			}
		}
	
		//Sanity check for presence of true rotamers

		bool trueRots = false;
		for (unsigned int j = 0; j < resMask[i].size(); j++) {
			if (resMask[i][j]) {
				trueRots = true;
			}
		}
		if (!trueRots) {
			cerr << "ERROR AT POSITION: " << i << endl;
			cerr << "Current Residue: " << activeResidue << endl;
			cerr << "resRots at this position: " << endl;
			for (uint k = 0; k < resRots[i].size(); k++) {
				cerr << resRots[i][k] << " ";
			}
			cerr << endl;
			cerr << "resMask at this position: " << endl;
			for (uint k = 0; k < resMask[i].size(); k++) {
				cerr << resMask[i][k] << " ";
			}
			cerr << endl;	
			exit(9123);
		}
	}
	return resMask;
}

/***********************************
 *define interface and rotamer sampling
 ***********************************/
vector<int> getRotamerSampling(string _rotamerLevels){
	vector<int> rotamerSampling;
	for (uint n=0; n<2; n++){
		for (uint i=0; i<_rotamerLevels.size(); i++){
			stringstream ss;
			ss << _rotamerLevels[i];
			rotamerSampling.push_back(MslTools::toInt(ss.str()));
		}
	}
	return rotamerSampling;
}

vector<int> getLinkedPositions(vector<int> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel){
	vector<int> positionsToLink;
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
			positionsToLink.push_back(1);
		} else {
			positionsToLink.push_back(0);
		}
	}
	return positionsToLink;
}

//Identify which positions are found at the identified interface (example: Sequence LLLLIGLLIGLLIGLLLL would be 000011001100110000 where positions at interface are 1, others are 0)
vector<int> getLinked(vector<int> _rotamerSampling, int _backboneLength, int _interfaceLevel, int _highestRotamerLevel){
	vector<int> linkedPositions;
	for (uint i=0; i<_backboneLength; i++){
		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
			linkedPositions.push_back(i);
		}
	}
	return linkedPositions;
}

// Convert positions to string for setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions) which uses "A,19" "B,19" format!
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<int> &_interfacialPositions, int _backboneLength){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_backboneLength; k++){//TODO: if you want to make it possible to do sequences that are not linked, need to change this function!
		if (_interfacialPositions[k] == 1){
			vector<string> tempPos;
	
			Position &posA = _sys.getPosition(k);
			Position &posB = _sys.getPosition(k+_backboneLength);
	
			string A = posA.toString();
			string B = posB.toString();
	
			string delimiter = " ";
			
			size_t p = 0;
			p = A.find(delimiter);
	
			tempPos.push_back(A.substr(0, p));
			tempPos.push_back(B.substr(0, p));
	
			stringPositions.push_back(tempPos);
		}
	}
	return stringPositions;
}

/***********************************
 *  calculate residue burial
 ***********************************/
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
	writer.open(_opt.pdbOutputDir+"/interfaceSASA.pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	return residueBurial;
}

void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis){
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
	pdb.readPdb(_opt.infile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
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
	string backboneSeq = generateString(_opt.backboneAA, _opt.backboneLength);
	vector<pair <int, double> > resiBurial = calculateResidueBurial(sys, _opt, backboneSeq);
	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});

	//cout << "Determining interfacial residues by residue burial..." << endl;
	int levelCounter = 0;
	vector<int> interfacePositions;

	// Output variable Set up
	backboneSeq = generateBackboneSequence("L", _opt.backboneLength);
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
		if (levelCounter < _opt.interfaceLevel){//TODO: change this to an option so I can choose how many levels I want to have multiple IDs
			add = "Add all Ids at this pos";
			interfacePositions.push_back(resiNum);
			if (backbonePosition > 3 && backbonePosition < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 3 and 4 here instead of 4 and 5 to prevent changes at the interface like others
				variablePositionString.replace(variablePositionString.begin()+posNum, variablePositionString.begin()+posNum+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		} else {
			add = "Only 1 ID";
		}
		rotamerLevels.replace(rotamerLevels.begin()+posNum, rotamerLevels.begin()+posNum+1, MslTools::intToString(levelCounter));
	}
	lvlavg = lvlavg/numAAs;
	avgs.push_back(lvlavg);
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.Ids, interfacePositions);

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
	_out << "Starting Sequence:  " << backboneSeq << endl;
	_out << "Variable Positions: " << variablePositionString << endl;
	_out << "Rotamers Levels:    " << rotamerSamplingString << endl;
	_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition);

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
	cout << "Starting Sequence:  " << backboneSeq << endl;
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

/***********************************
 *output file functions
 ***********************************/
void setupDesignDirectory(Options &_opt, string _date){
	//_opt.pdbOutputDir = _opt.pdbOutputDir + "/" + _date;
	_opt.pdbOutputDir = get_current_dir_name();
	_opt.pdbOutputDir = _opt.pdbOutputDir + "/10_30_2021";
	string cmd = "mkdir -p " + _opt.pdbOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	_opt.pdbOutputDir = _opt.pdbOutputDir + "/design_" + _opt.runNumber;
	cmd = "mkdir -p " + _opt.pdbOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

void printConfigFile(Options & _opt, ofstream & _out) {
	_out << "#Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
	_out << "#All other parameters were run with the defaults as specified in the program options" << endl;

	_out << setw(20) << "backboneCrd " << _opt.backboneCrd << endl;
	_out << setw(20) << "pdbOutputDir " << _opt.pdbOutputDir << endl;

	_out << setw(20) << "tmStart " << _opt.tmStart << endl;
	_out << setw(20) << "tmEnd " << _opt.tmEnd << endl;

	_out << "#Input Files" << endl;
	_out << setw(20) << "topFile " << _opt.topFile << endl;
	_out << setw(20) << "parFile " << _opt.parFile << endl;
	_out << setw(20) << "geometryDensityFile " << _opt.geometryDensityFile << endl;
	_out << setw(20) << "rotLibFile " << _opt.rotLibFile << endl;
	_out << setw(20) << "solvFile " << _opt.solvFile << endl;
	_out << setw(20) << "backboneCrd " << _opt.backboneCrd << endl;
	_out << setw(20) << "hbondFile " << _opt.hbondFile << endl;
	_out << setw(20) << "infile " << _opt.infile << endl;
	_out << setw(20) << "selfEnergyFile " << _opt.selfEnergyFile << endl;
	_out << setw(20) << "pairEnergyFile " << _opt.pairEnergyFile << endl;
	_out << setw(20) << "AACompositionPenaltyFile " << _opt.AACompositionPenaltyFile << endl << endl;

	_out << "#Geometry and Transformation parameters" << endl;
	_out << setw(20) << "xShift" << _opt.xShift << endl;
	_out << setw(20) << "crossingAngle" << _opt.crossingAngle << endl;
	_out << setw(20) << "axialRotation" << _opt.axialRotation << endl;
	_out << setw(20) << "zShift" << _opt.zShift << endl;
	_out << setw(20) << "thread" << _opt.thread << endl;
	_out << setw(20) << "backboneLength " << _opt.backboneLength << endl;

	_out << "#Booleans" << endl;
	_out << setw(20) << "verbose " << _opt.verbose << endl;
	_out << setw(20) << "deleteTerminalHbonds" << _opt.deleteTerminalHbonds << endl;
	_out << setw(20) << "useSasa" << _opt.useSasa << endl;
	_out << setw(20) << "getGeoFromPDBData" << false << endl;//Since we already have the geometry output here, default to false in the rerun config
	_out << setw(20) << "runDEESingles" << _opt.runDEESingles << endl;
	_out << setw(20) << "runDEEPairs" << _opt.runDEEPairs << endl;
	_out << setw(20) << "runSCMF" << _opt.runSCMF << endl;

	if (_opt.useSasa == true){
		_out << "#Load Rotamers based on SASA scores" << endl;
		for (uint i=0; i<_opt.sasaRepackLevel.size()-1; i++){
			_out << setw(20) << "sasaRepackLevel" << _opt.sasaRepackLevel[i] << endl;
		}
		_out << setw(20) << "interfaceLevel" << _opt.interfaceLevel << endl;
	} else {
		_out << "#Load Rotamers by interface and non interfacial positions" << endl;
		_out << setw(20) << "SL" << _opt.SL << endl;
		_out << setw(20) << "SLInterface" << _opt.SLInterface << endl;
	}

	_out << "#MonteCarlo Paramenters" << endl;
	_out << setw(20) << "MCCycles" << _opt.MCCycles << endl;
	_out << setw(20) << "MCMaxRejects" << _opt.MCMaxRejects << endl;
	_out << setw(20) << "MCStartTemp" << _opt.MCStartTemp << endl;
	_out << setw(20) << "MCEndTemp" << _opt.MCEndTemp << endl;
	_out << setw(20) << "MCCurve" << _opt.MCCurve << endl;
	_out << setw(20) << "greedyCycles" << _opt.greedyCycles << endl;

	_out << "#Alternate IDs" << endl;
	for (uint i=0; i<_opt.Ids.size()-1; i++){
		_out << setw(20) << _opt.Ids[i] << endl;
	}

	_out << endl << "#Rerun Seed" << endl;
	_out << setw(20) << "seed" << _opt.seed << endl;

	_out << setw(20) << "weight_vdw " << _opt.weight_vdw << endl;
	_out << setw(20) << "weight_hbond " << _opt.weight_hbond << endl;
	_out << setw(20) << "weight_solv " << _opt.weight_solv << endl;
	_out << setw(20) << "weight_seqEntropy " << _opt.weight_seqEntropy << endl;
}

void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns){
	ofstream eout;
	string eoutfile = _opt.pdbOutputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	string tab = "\t";

	// Outputs
	eout << tab << endl;
	for (uint i=0; i<_opt.energyTermsToOutput.size(); i++){
		eout << _opt.energyTermsToOutput[i] << tab;
		cout << _opt.energyTermsToOutput[i] << tab;
	}
	eout << "Sequence" << tab;
	eout << "InterfaceSeq" << tab;
	eout << "xShift" << tab;
	eout << "crossingAngle" << tab;
	eout << "axialRotation" << tab;
	eout << "zShift" << tab;
	eout << "angleDistDensity" << tab;
	eout << "axialRotationDensity" << tab;
	eout << "zShiftDensity" << tab;
	eout << "repackLevels" << tab;
	eout << "interfaceLevels" << tab;
	eout << "backboneLength" << tab;
	eout <<	"PDBPath" << endl; 
	
	cout << "Sequence" << tab;
	cout << "InterfaceSeq" << tab;
	cout << "xShift" << tab;
	cout << "crossingAngle" << tab;
	cout << "axialRotation" << tab;
	cout << "zShift" << tab;
	cout << "angleDistDensity" << tab;
	cout << "axialRotationDensity" << tab;
	cout << "zShiftDensity" << tab;
	cout << "repackLevels" << tab;
	cout << "interfaceLevels" << tab;
	cout << "backboneLength" << tab;
	cout <<	"PDBPath" << endl; 

	for (uint i=0; i<_allDesigns.size() ; i++){
		eout << _allDesigns[i] << endl;
		cout << _allDesigns[i] << endl;
	}
	eout.close();
}

//void makeRepackConfig(Options &_opt, string _sequence, vector<uint> _state, string _designNumber, string _pdbPath, ofstream &_globalRepackOut){
void makeRepackConfig(Options &_opt, string _sequence, string _designDir, vector<uint> _state, string _designNumber, string _pdbPath, string _crdPath, map<string,double> _energyMap, vector<int> _rotamerSampling){
	ofstream dout;
	string doutfile = _opt.pdbOutputDir + "/" + _designNumber +"/repack.config";
	dout.open(doutfile.c_str());
	string designFileDir = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/";//Running on chtc server and needed these to access the proper directories after transferring them back to our server
	
	dout << "#Input Files" << endl;
	dout << "topFile            " << designFileDir + _opt.topFile << endl;
	dout << "parFile            " << designFileDir + _opt.parFile << endl;
	dout << "solvFile           " << designFileDir + _opt.solvFile << endl;
	dout << "hbondFile          " << designFileDir + _opt.hbondFile << endl;
	dout << "rotLibFile         " << designFileDir + _opt.rotLibFile << endl;
	dout << "designPdb          " << _pdbPath << endl;
	dout << "designCrd          " << _crdPath << endl;
	dout << "outputDir          " << _designDir << endl << endl;
	
	dout << "#Geometry" << endl;
	dout << "xShift             " << _opt.xShift << endl;
	dout << "crossingAngle      " << _opt.crossingAngle << endl;
	dout << "axialRotation      " << _opt.axialRotation << endl;
	dout << "zShift             " << _opt.zShift << endl;
	dout << "thread             " << _opt.thread << endl;
	dout << "tmStart            " << _opt.tmStart << endl;
	dout << "tmEnd              " << _opt.tmEnd << endl;
	
	dout << "#Design Parameters" << endl;
	dout << "sequence              " << _sequence << endl;
	dout << "rotamerSamplingString ";
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (i<_rotamerSampling.size()-1){
			dout << _rotamerSampling[i];
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << "seed               " << _opt.seed << endl;
	dout << "state              ";
	for (uint i=0; i<_state.size(); i++){
		if (i<_state.size()-1){
			dout << _state[i] << " ";
		} else {
			dout << _state[i] << endl;
		}
	}
	dout << endl;
	dout << "#Rotamer Sampling Vector" << endl;
	dout << "rotamerSampling           ";
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (i<_rotamerSampling.size()-1){
			dout << _rotamerSampling[i] << " ";
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}

	dout << "designNumber       " << _designNumber << endl;
	dout << endl;

	// SASA Rotamer Values
	dout << "#SASA Rotamers" << endl;
	for (uint i=0; i<_opt.sasaRepackLevel.size(); i++){
		dout << "sasaRepackLevel  " << _opt.sasaRepackLevel[i] << endl;
	}
	dout << endl;

	// Weights
	dout << "#Weights" << endl;
	dout << "weight_vdw         " << _opt.weight_vdw << endl;
	dout << "weight_hbond       " << _opt.weight_hbond << endl;
	dout << "weight_solv        " << _opt.weight_solv << endl << endl;

	// Monomer Energies
	dout << "#Monomer Energies" << endl;
	dout << "monomer            " << _energyMap.at("Monomer") << endl;
	dout << "monoVdw            " << _energyMap.at("VDWMonomer") << endl;
	dout << "monoHbond          " << _energyMap.at("HBONDMonomer") << endl;
	dout << "monoIMM1           " << _energyMap.at("IMM1Monomer") << endl << endl;

	// Energies before repack
	dout << "#Dimer Energies" << endl;
	dout << "dimer            " << _energyMap.at("Dimer") << endl;
	dout << "dimerVdw         " << _energyMap.at("VDWDimer") << endl;
	dout << "dimerHbond       " << _energyMap.at("HBONDDimer") << endl;
	dout << "dimerIMM1        " << _energyMap.at("IMM1Dimer") << endl << endl;
	
	dout << "#Total Energy" << endl;
	dout << "Total            " << _energyMap.at("Total") << endl;
	dout.close();
}

void makeDockingConfig(Options &_opt, string _sequence, string _designDir, string _designNumber, string _pdbPath, map<string,double> _energyMap, vector<int> _rotamerSampling){
	ofstream dout;
	string outputDir = _opt.pdbOutputDir + "/" + _designNumber;
	string doutfile = outputDir + "/docking.config";
	dout.open(doutfile.c_str());
	dout << "#Design Parameters" << endl;
	dout << "chainSeq           " << _sequence << " " << _sequence << endl;
	dout << "rotamerSamplingString ";
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (i<_rotamerSampling.size()-1){
			dout << _rotamerSampling[i];
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << "chainStartNum      " << _opt.thread << endl << endl;
	dout << "outputDir          " << _designDir << endl;
	dout << "seed               " << _opt.seed << endl;

	string designFileDir = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/";
	dout << "#Input Files" << endl;
	dout << "topFile            " << designFileDir + _opt.topFile << endl;
	dout << "parFile            " << designFileDir + _opt.parFile << endl;
	dout << "solvFile           " << designFileDir + _opt.solvFile << endl;
	dout << "hbondFile          " << designFileDir + _opt.hbondFile << endl;
	dout << "rotLibFile         " << designFileDir + _opt.rotLibFile << endl;
	dout << "bbqFile            " << designFileDir << "PiscesBBQTable.txt" << endl;
	dout << "designPdb          " << _pdbPath << endl;//In case I ever wnat to do an RMSD calculation against the orignal

	dout << "#Geometry" << endl;
	dout << "xShift             " << _opt.xShift << endl;
	dout << "crossingAngle      " << _opt.crossingAngle << endl;
	dout << "axialRotation      " << _opt.axialRotation << endl;
	dout << "zShift             " << _opt.zShift << endl;
	dout << endl;

	dout << "#Rotamer Sampling Vector" << endl;
	dout << "rotamerSampling    ";
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (i<_rotamerSampling.size()-1){
			dout << _rotamerSampling[i] << " ";
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << "#Weights" << endl;
	dout << "weight_vdw         " << _opt.weight_vdw << endl;
	dout << "weight_hbond       " << _opt.weight_hbond << endl;
	dout << "weight_solv        " << _opt.weight_solv << endl << endl;

	// Monomer Energies
	dout << "#Monomer Energies" << endl;
	dout << "monomer            " << _energyMap.at("Monomer") << endl;
	dout << "monoVdw            " << _energyMap.at("VDWMonomer") << endl;
	dout << "monoHbond          " << _energyMap.at("HBONDMonomer") << endl;
	dout << "monoIMM1           " << _energyMap.at("IMM1Monomer") << endl << endl;

	// Energies before repack
	dout << "#Dimer Energies" << endl;
	dout << "dimer            " << _energyMap.at("Dimer") << endl;
	dout << "dimerVdw         " << _energyMap.at("VDWDimer") << endl;
	dout << "dimerHbond       " << _energyMap.at("HBONDDimer") << endl;
	dout << "dimerIMM1        " << _energyMap.at("IMM1Dimer") << endl << endl;
	
	dout << "#Total Energy" << endl;
	dout << "Total            " << _energyMap.at("Total") << endl;
	dout.close();
}

void outputRepackFile(Options &_opt, vector<string> _dockingDesigns){
	ofstream dout;
	string doutfile = _opt.pdbOutputDir + "/designsToRepack.out";
	dout.open(doutfile.c_str());
	for (uint i=0; i<_dockingDesigns.size() ; i++){
		dout << _dockingDesigns[i] << endl;
	}
	dout.close();
}

void outputDesignFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, vector<pair<string,vector<uint>>> _sequenceStatePair, map<string,map<string,double>> _sequenceEnergyMap, vector<double> _densities){
	// Setup vector to hold energy file lines
	vector<string> energyLines;

	// Setup the parameters for this specific run
	string tab = "\t";
	string xShift = MslTools::doubleToString(_opt.xShift);
	string crossingAngle = MslTools::doubleToString(_opt.crossingAngle);
	string axialRotation = MslTools::doubleToString(_opt.axialRotation);
	string zShift = MslTools::doubleToString(_opt.zShift);
	string angleDistDensity = MslTools::doubleToString(_densities[0]);
	string axialRotationDensity = MslTools::doubleToString(_densities[1]);
	string zShiftDensity = MslTools::doubleToString(_densities[2]);
	string thread = MslTools::intToString(_opt.thread);
	string repackLevels = MslTools::intToString(_opt.sasaRepackLevel.size());
	string interfaceLevel = MslTools::intToString(_opt.interfaceLevel);
	string bbLength = MslTools::intToString(_opt.backboneLength);
	//string runParameters = xShift+tab+crossingAngle+tab+axialRotation+tab+zShift+tab+thread+tab+repackLevels+tab+interfaceLevel;
	string runParameters = xShift+tab+crossingAngle+tab+axialRotation+tab+zShift+tab+angleDistDensity+tab+axialRotationDensity+tab+zShiftDensity+tab+repackLevels+tab+interfaceLevel+tab+bbLength;
	// For loop to setup the energy file
	for (uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		vector<uint> state = _sequenceStatePair[i].second;
		map<string,double> energyMap = _sequenceEnergyMap.at(sequence);
		// For adding in strings to a line for the energy file
		string tmp = "Sequence Info: ";
		for (uint j=0; j<_opt.energyTermsToOutput.size(); j++){
			string energyTerm = _opt.energyTermsToOutput[j];
			double energy = energyMap.at(energyTerm);
			tmp.append(MslTools::doubleToString(energy));
			tmp.append(tab);
		}

		//Add in path to design PDB and make repack and docking configuration files
		string designNumber = MslTools::doubleToString(energyMap.at("DesignNumber"));
		//string designDir = _opt.pdbOutputDir + "/" + designNumber;
		string designDir = "/data02/gloiseau/Sequence_Design_Project/vdwSequenceDesign/sequenceDesign/10_30_2021/design_" + _opt.runNumber + "/"+ designNumber;
		string pdbPath = designDir + "/" + designNumber + ".pdb";
		string crdPath = designDir + "/" + designNumber + ".crd";
		makeRepackConfig(_opt, sequence, designDir, state, designNumber, pdbPath, crdPath, energyMap, _rotamerSamplingPerPosition);
		makeDockingConfig(_opt, sequence, designDir, designNumber, pdbPath, energyMap, _rotamerSamplingPerPosition);
	
		// Append other important features to the end energy files lines
		string interfaceSequence = getInterfaceSequence(_opt,_interface, sequence);
		tmp.append(sequence);
		tmp.append(tab);
		tmp.append(interfaceSequence);
		tmp.append(tab);
		tmp.append(runParameters);
		tmp.append(tab);
		tmp.append(pdbPath);
		energyLines.push_back(tmp);
		tmp.clear();
	}
	outputEnergyFile(_opt, _interface, energyLines);
}

/***********************************
 *load rotamer functions
 ***********************************/
void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), "SL90.00")) {//lower rotamer level because I did baselines at this level
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}	
}

void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling){
	//Repack side chains based on sasa scores
	for (uint i = 0; i < _rotamerSampling.size()/2; i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+_opt.backboneLength);
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _opt.sasaRepackLevel[_rotamerSampling[i]];
				if (_opt.verbose){
					cout << posA.getPositionId() << ", " << posA.getResidueName() << ":" << posRot << endl;
					cout << posB.getPositionId() << ", " << posB.getResidueName() << ":" << posRot << endl;
				}
				if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot)) { 
						cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
					}
					if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot)) { 
						cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
					}
				}
			}
		} else {
			string posRot = _opt.sasaRepackLevel[_rotamerSampling[i]];
			if (_opt.verbose){
				cout << posA.getPositionId() << ", " << posA.getResidueName() << ": " << posRot << endl;
				cout << posB.getPositionId() << ", " << posB.getResidueName() << ": " << posRot << endl;
			}
			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot)) { 
					cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
				}
				if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot)) { 
					cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
				}
			}
		}
	}
}

//below function only loads rotamers onto the interfacial positions by interfacialPositions (01 where 0 = non-interfacial and 1 = interfacial)
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<int> _interface){
	for (uint k=0; k<_interface.size(); k++) {
		if (_interface[k] < _numRotamerLevels){
			Position &pos = _sys.getPosition(k);
			if (pos.identitySize() > 1){
				for (uint j=0; j < pos.getNumberOfIdentities(); j++){
					pos.setActiveIdentity(j);
					if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
						if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
							cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
						}
					}
					pos.setActiveIdentity(0);
				}
			} else {
				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
					}
				}
			}
		}
	}
}

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.identitySize() > 1){
			for (uint j=0; j < pos.getNumberOfIdentities(); j++){
				pos.setActiveIdentity(j);
				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
					}
				}
				pos.setActiveIdentity(0);
			}
		} else {
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}
	}	
}

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling){
	if (_opt.useSasa){
		if (_opt.verbose){
			cout << "Load rotamers by difference in residue burial..." << endl;
		}
		loadRotamersBySASABurial(_sys, _sysRot, _opt, _rotamerSampling);
	} else {
		if (_opt.verbose){
			cout << "Load rotamers..." << endl;
			cout << "Non-interface: " << _opt.SL << endl;
			cout << "Interface:     " << _opt.SLInterface << endl;
		}
		loadRotamers(_sys, _sysRot, _opt.SL);
		loadInterfacialRotamers(_sys, _sysRot, _opt.SLInterface, _opt.sasaRepackLevel.size(), _rotamerSampling);
	}
}

/***********************************
 *baseline energy helper functions
 ***********************************/
//Function to calculate the self energies of a chain
vector<double> calcBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	//cout << "Self Energies" << endl;
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
		//cout << number << ": " << resi << endl;
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to calculate the pair energies of a chain
vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+1);
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+1);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			} else {
				j = _seqLength;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to get the sum of a vector of doubles, typically energies
double sumEnergyVector(vector<double> _energies){
	double ener = 0;
	for (uint i=0; i<_energies.size(); i++){
		ener = ener + _energies[i];
	}
	return ener;
}

/***********************************
 *calculate energies
 ***********************************/
void computeDimerEnergiesLinked(System &_sys, Options &_opt, map<string,map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err) {
	
	EnergySet* Eset = _sys.getEnergySet();
	_sys.calcEnergy();

	Eset->eraseTerm("BASELINE");
	Eset->eraseTerm("BASELINE_PAIR");
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	for(uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		vector<uint> stateVec = _sequenceStatePair[i].second;
		_sys.setActiveRotamers(stateVec);
		
		//Calculate Dimer energy for output
		double dimerEnergy = _sys.calcEnergy();
		double finalEnergy = dimerEnergy-_sequenceEnergyMap[sequence]["Monomer"];//TODO: make this more elegant by pulling from the map 
		_sout << "Sequence: " << sequence << endl;
		_sout << "Dimer Energy: " << dimerEnergy << endl;
		_sout << "Final Energy = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl << endl;
		
		map<string,double> &energyMap = _sequenceEnergyMap[sequence];
		outputEnergiesByTermLinked(Eset, energyMap, _opt.energyTermList, "Dimer");
		
		//Setup directory for individual design
		string repackDir = _opt.pdbOutputDir + "/" + MslTools::intToString(i);
		string cmd = "mkdir -p " + repackDir;
		if (system(cmd.c_str())){
			_err << "Unable to make directory" << endl;
			exit(0);
		}
		_sequenceEnergyMap[sequence]["Dimer"] = dimerEnergy;
		_sequenceEnergyMap[sequence]["Total"] = finalEnergy;
		_sequenceEnergyMap[sequence]["DesignNumber"] = i;
		
		_sys.setActiveRotamers(stateVec);
		PDBWriter designWriter;
		designWriter.open(repackDir + "/" + MslTools::intToString(i) + ".pdb");
		designWriter.write(_sys.getAtomPointers(), true, false, true);
		designWriter.close();
		_writer.write(_sys.getAtomPointers(), true, false, true);
		saveEnergyDifference(_opt, _sequenceEnergyMap, sequence);
		
		//Write CRD if below energy cutoff
		if (finalEnergy < _opt.repackEnergyCutoff){
			CRDWriter writer;
			writer.open(repackDir + "/" + MslTools::intToString(i) + ".crd");
			writer.write(_sys.getAtomPointers(), false);
			writer.close();
		}
	}
}

void computeDimerEnergies(System &_sys, Options &_opt, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> _rotamerSamplingPerPosition, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	// Setup timer for LocalMC Cycles
	time_t startTimeLMC, endTimeLMC;
	double diffTimeLMC;
	time(&startTimeLMC);
	
	// Initialize PDBWriter for design
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/allDesigns.pdb");

	_sout << "Calculating dimer energies..." << endl;
	cout << "Calculating dimer energies..." << endl;
	if (_opt.linkInterfacialPositions){	
		computeDimerEnergiesLinked(_sys, _opt, _sequenceEnergyMap, _sequenceStatePair, _rotamerSamplingPerPosition, _linkedPos, _RNG, writer, _sout, _err);
	} else {
		for (uint i=0; i<_sequenceStatePair.size(); i++){
			string sequence = _sequenceStatePair[i].first;
			vector<uint> state = _sequenceStatePair[i].second;
			cout << "Sequence " << i+1 << ": " << sequence << endl;
			_sout << "Sequence " << i+1 << ": " << sequence << endl;

			computeDimerEnergy(_sys, _opt, _sequenceEnergyMap, sequence, state, _rotamerSamplingPerPosition, _linkedPos, i, _RNG, writer, _sout, _err);
			_sequenceStatePair[i].second = state;
		}
	}
	
	time(&endTimeLMC);
	diffTimeLMC = difftime (endTimeLMC, startTimeLMC);
	writer.close();
	cout << "End dimer Calculations: " << diffTimeLMC << "s" << endl << endl;
	_sout << "End dimer Calculations: " << diffTimeLMC << "s" << endl << endl;
}

void computeDimerEnergy(System &_sys, Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_sequence, vector<uint> &_stateVec, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, int _seqNumber, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err) {
	
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PL(polySeq);

	// Declare new system
	System sys;
	CharmmSystemBuilder CSB(sys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	CSB.setBuildNonBondedInteractions(false);
	if(!CSB.buildSystem(PL)) {
		_err << "Unable to build system from " << PL << endl;
		exit(0);
	}

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(_sys.getAtomPointers(),false);
	sys.buildAllAtoms();

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	Eset->setAllTermsActive();
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

	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	Eset->setWeight("CHARMM_IMM1", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);
	
	/******************************************************************************
	 *              === LOAD ROTAMERS AND CHOOSE TO LINK INTERFACE ===
	 ******************************************************************************/
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);
	CSB.updateNonBonded(10,12,50);

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	// Setup spm and calculate energies
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();

	//Calculate Dimer energy for output
	repackSideChains(spm, 10);
	_stateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(_stateVec);
	
	double dimerEnergy = spm.getStateEnergy(_stateVec);
	double finalEnergy = dimerEnergy-_sequenceEnergyMap[_sequence]["Monomer"];//TODO: make this more elegant by pulling from the map 
	_sout << "Dimer Energy: " << dimerEnergy << endl;
	_sout << "Final Energy w/ IMM1 = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;
	cout << "Dimer Energy: " << dimerEnergy << endl;
	cout << "Final Energy = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;
	
	map<string,double> &energyMap = _sequenceEnergyMap[_sequence];
	outputEnergiesByTerm(spm, _stateVec, energyMap, _opt.energyTermList, "Dimer", true);
	
	//Setup directory for individual design
	string repackDir = _opt.pdbOutputDir + "/" + MslTools::intToString(_seqNumber);
	string cmd = "mkdir -p " + repackDir;
	if (system(cmd.c_str())){
		_err << "Unable to make directory" << endl;
		exit(0);
	}
	_sequenceEnergyMap[_sequence]["Dimer"] = dimerEnergy;
	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;
	_sequenceEnergyMap[_sequence]["DesignNumber"] = _seqNumber;
	
	sys.setActiveRotamers(_stateVec);
	PDBWriter designWriter;
	designWriter.open(repackDir + "/" + MslTools::intToString(_seqNumber) + ".pdb");
	designWriter.write(sys.getAtomPointers(), true, false, true);
	designWriter.close();
	_writer.write(sys.getAtomPointers(), true, false, true);
	saveEnergyDifference(_opt, _sequenceEnergyMap, _sequence);
	
	//Write CRD if below energy cutoff
	if (finalEnergy < _opt.repackEnergyCutoff){
		CRDWriter writer;
		writer.open(repackDir + "/" + MslTools::intToString(_seqNumber) + ".crd");
		writer.write(sys.getAtomPointers(), false);
		writer.close();
	}
}

void computeMonomerEnergyNoIMM1(Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {
	
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

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);
	
	/*****************************************************************************
	 *                 === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, "SL95.00");
	CSBMono.updateNonBonded(10,12,50);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	/*****************************************************************************
	 *            === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ===
	 ******************************************************************************/
	repackSideChains(monoSpm, 10);
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getStateEnergy(stateVec)*2;

	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "MonomerNoIMM1", false);
	_sequenceEnergyMap[_seq]["MonomerNoIMM1"] = monomerEnergy;
	vector<double> selfVec = calcBaselineEnergies(monoSys, _opt.backboneLength);
	vector<double> pairVec = calcPairBaselineEnergies(monoSys, _opt.backboneLength);
	//double self = sumEnergyVector(selfVec);
	//double pair = sumEnergyVector(pairVec);
	_sout << "Monomer Energy No IMM1: " << monomerEnergy << endl;
	//cout << "Self Energy:  " << self << endl;
	//cout << "Pair Energy:  " << pair << endl;
	//cout << "Total Energy: " << self+pair << endl;
}

void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {
	
	//string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, _opt.thread);//fixed monomer calculation issue on 05_12_2021
	string polySeq = generateMonomerPolymerSequenceFromSequence(_seq, _opt.thread);
	PolymerSequence PS(polySeq);

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

	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
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

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);
	
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);
	monoEset->setTermActive("CHARMM_IMM1REF", true);
	monoEset->setTermActive("CHARMM_IMM1", true);
	monoEset->setTermActive("CHARMM_VDW", true);
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);
	monoEset->setWeight("CHARMM_IMM1REF", 1);
	monoEset->setWeight("CHARMM_IMM1", 1);
	
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
		cout << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, "SL95.00");
	CSBMono.updateNonBonded(10,12,50);
	
	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	repackSideChains(monoSpm, _opt.greedyCycles);
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);

	monoSys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");

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
	monoSys.calcEnergy();

	// move center of mass to origin
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
	helicalAxis.saveAltCoor("BestMonomerAxis");
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
				helicalAxis.saveAltCoor("BestMonomerAxis");
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
			helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			//_fout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
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
	double monomerEnergy = monoSpm.getStateEnergy(stateVec)*2;
	_sout << "Monomer Energy w/ IMM1: " << monomerEnergy << endl;
	cout << "Monomer Energy w/ IMM1: " << monomerEnergy << endl;
	
	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	
	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	helicalAxis.clearSavedCoor("BestMonomerAxis");
	helicalAxis.clearSavedCoor("bestZ");
}

void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	cout << "Calculating monomer energies..." << endl;
	_sout << "Calculating monomer energies..." << endl;
	for (uint m=0; m<_seqs.size(); m++){
		cout << "Sequence " << m+1 << ": " << _seqs[m] << endl;
		_sout << "Sequence " << m+1 << ": " << _seqs[m] << endl;
		computeMonomerEnergyNoIMM1(_opt, _sequenceEnergyMap, _seqs[m], _RNG, _sout, _err);
		computeMonomerEnergyIMM1(_opt, _trans, _sequenceEnergyMap, _seqs[m], _RNG, _sout, _err);
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	cout << "End monomer calculations: " << diffTimeMono << "s" << endl << endl;
	_sout << "End monomer calculations: " << diffTimeMono << "s" << endl << endl;
}


/***********************************
 *other helper functions
 ***********************************/
void saveEnergyDifference(Options _opt, map<string,map<string,double>> &_sequenceEnergyMap, string _sequence){
	map<string,double> &energyMap = _sequenceEnergyMap[_sequence];
	energyMap["Baseline-Monomer"] = energyMap["Baseline"] + energyMap["MonomerNoIMM1"];
	energyMap["HBONDDiff"] = energyMap["HBONDDimer"] - energyMap["HBONDMonomer"];
	energyMap["VDWDiff"] = energyMap["VDWDimer"] - energyMap["VDWMonomer"];
	energyMap["IMM1Diff"] = energyMap["IMM1Dimer"] - energyMap["IMM1Monomer"];
}

void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap, vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1){
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
			if (_energyDescriptor.find("Monomer") != string::npos){
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF"))*2;
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
				}
			} else {
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF");
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
				}
			}
		}
	}
}

void outputEnergiesByTermLinked(EnergySet *_Eset, map<string,double> &_energyMap, vector<string> _energyTermList, string _energyDescriptor){
	for (uint i=0; i<_energyTermList.size(); i++){
		string energyTerm = _energyTermList[i];
		string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;
		if (energyTerm.find("IMM1") != string::npos){
			_energyMap["IMM1Dimer"] = _Eset->getTermEnergy("CHARMM_IMM1")+_Eset->getTermEnergy("CHARMM_IMM1REF");
		} else {
			_energyMap[energyLabel] = _Eset->getTermEnergy(energyTerm);
		}
	}
}

void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	int frontExt = _opt.tmStart;
	int endExt = _opt.endResNum;
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(frontExt > i) {
				atoms += positions[i]->getAtomPointers();
			}
			if(endExt > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
			}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");
}

//TODO: add the below functions to baseline files?
map<string, double> readSingleParameters(string _baselineFile){
	Reader selfReader(_baselineFile);
	selfReader.open();
	map<string, double> selfEnergies;

	if(!(selfReader.is_open())){
		cerr << "WARNING: Unable to open " << _baselineFile << endl;
		exit(0);
	}

	vector<string> lines = selfReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 2){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) Energy(double)";
			continue;
		}
		selfEnergies[MslTools::toUpper(tokens[0])] = MslTools::toDouble(tokens[1]);
		//cout << tokens[0] << " " << tokens[1] << " " << tokens[2] << " = " << tokens[3] << endl;
	}
	
	selfReader.close();
	return selfEnergies;
}

map<string,map<string,map<uint, double>>> readPairParameters(string _baselineFile){
	Reader pairReader(_baselineFile);
	pairReader.open();
	map<string,map<string,map<uint, double>>> pairEnergies;

	if(!(pairReader.is_open())){
		cerr << "WARNING: Unable to open " << _baselineFile << endl;
		exit(0);
	}

	vector<string> lines = pairReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], " ");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) ResName(string) Distance(uint) Energy(double)";
			continue;
		}
		if (tokens[0].compare(tokens[1]) == 0){//Added in on 03-18-2021: apparently the code that I was using to flip the AAs in buildPairInteractions isn't good, but this adds the flips to the map which works better and is cleaner
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		} else {
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
			pairEnergies[MslTools::toUpper(tokens[1])][MslTools::toUpper(tokens[0])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		}
	}
	
	pairReader.close();
	return pairEnergies;
}

/***********************************
 *stateMC helper functions
 ***********************************/
void randomRotamerChange(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _variablePositions, vector<unsigned int> &_stateVec){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _variablePositions.size()-1);
	int pos = _variablePositions[rand];

	// Get the random position from the system
	Position &randPos = _sys.getPosition(pos);
	string posId = randPos.getPositionId();
	int numRots = randPos.getTotalNumberOfRotamers();
	
	// Get a random integer to pick the rotamer
	int randRot = _RNG.getRandomInt(0, numRots-1);
	_stateVec[pos] = randRot;
	randPos.setActiveRotamer(randRot);
}

//This version is for unlinked positions, but keeps the AAs the same, allowing me to change rotamers for both spots but keep the same sequence (posRotLimitMap has limits for each AA)
void randomRotamerChangeNonLinked(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, map<int,map<int,pair<uint,uint>>> &_posRotLimitMap, vector<uint> _variablePositions, vector<unsigned int> &_stateVec){
	// Get a random position from the variable positions
	int rand = _RNG.getRandomInt(0, _variablePositions.size()-1);
	int posA = _variablePositions[rand];
	int posB = posA+_opt.backboneLength;

	// Get the random position from the system
	Position &randPosA = _sys.getPosition(posA);
	Position &randPosB = _sys.getPosition(posB);
	string posIdA = randPosA.getPositionId();
	string posIdB = randPosB.getPositionId();

	// Get rotamer upper and lower limits for the position from the map (keeps AA the same on both chains)
	int idNum = _RNG.getRandomInt(0, randPosA.getNumberOfIdentities()-1);
	map<int,pair<uint,uint>> posRotLimits = _posRotLimitMap.at(posA);
	pair<uint,uint> rotLimits = posRotLimits.at(idNum);

	// Get a random integer to pick the rotamer
	int randRotA = _RNG.getRandomInt(rotLimits.first, rotLimits.second-1);
	int randRotB = _RNG.getRandomInt(rotLimits.first, rotLimits.second-1);
	_stateVec[posA] = randRotA;
	_stateVec[posB] = randRotB;
	randPosA.setActiveRotamer(randRotA);
	randPosB.setActiveRotamer(randRotB);
}

//Checks through a vector of sequences to see if the new sequence has already been found
void sameSequenceChecker(string &_newSeq, vector<string> &_seqs){
	bool sameSeq = false;
	if (_seqs.size() == 0){
		_seqs.push_back(_newSeq);
	} else {
		for (int j=0; j<_seqs.size(); j++){
			if (_seqs[j] == _newSeq){
				sameSeq = true;
				j = _seqs.size()-1;
			} else if (j==_seqs.size()-1){
				if (sameSeq == false){
					_seqs.push_back(_newSeq);
					j = _seqs.size();
				}
			}
		}
	}
}

//new functions for saving sequences in order and not having repeats (makes sure I will save x best sequences rather than x-duplicates)
bool sameSequenceChecker(string &_newSeq, double &_newEnergy, vector<uint> &_state, vector<pair<double,string>> &_enerSeqPair, vector<pair<double,vector<uint>>> &_energyStateVec){
	bool sameSeq = false;
	if (_enerSeqPair.size() == 0){
		sameSeq = false;
	} else {
		for (int j=0; j<_enerSeqPair.size(); j++){
			if (_enerSeqPair[j].second == _newSeq){
				sameSeq = true;
				if (_newEnergy < _enerSeqPair[j].first){
					_enerSeqPair[j].first = _newEnergy;
					_energyStateVec[j].first = _newEnergy;
					_energyStateVec[j].second = _state;
				}
				j = _enerSeqPair.size()-1;
			}
		}
	}
	return sameSeq;
}

void saveSequence(Options &_opt, vector<pair<double,string>> &_energyVector, vector<pair<double,vector<uint>>> &_energyStateVec, string _sequence, vector<uint> _state, double _energy){
	bool sameSeq = sameSequenceChecker(_sequence, _energy, _state, _energyVector, _energyStateVec);
	if (sameSeq == false){
		if (_energyVector.size() < _opt.numStatesToSave){
			_energyVector.push_back(make_pair(_energy, _sequence));
			_energyStateVec.push_back(make_pair(_energy, _state));
		} else {
			sort(_energyVector.begin(), _energyVector.end());
			sort(_energyStateVec.begin(), _energyStateVec.end());
			int lim = _energyVector.size()-1;
			if (_energyVector[lim].first > _energy){
				_energyVector[lim].first = _energy;
				_energyVector[lim].second = _sequence;
				_energyStateVec[lim].first = _energy;
				_energyStateVec[lim].second = _state;
			}
		}
	}
}

//Setup the map for limititing the number of rotamers per position (allows me to keep AAs the same without needing the rotamers to be symmetric)
map<int,map<int,pair<uint, uint>>> setupRotamerPositionMap(System &_sys, vector<uint> _interfacialPositionsList){
	map<int,map<int,pair<uint,uint>>> posRotLimitMap;
	for (uint i=0; i<_interfacialPositionsList.size(); i++){
		Position &pos = _sys.getPosition(_interfacialPositionsList[i]);
		uint startNum = 0;
		uint endNum = 0;
		map<int,pair<uint,uint>> currIdMap;
		for(uint j=0; j<pos.identitySize(); j++){
			Residue &resi = pos.getIdentity(j);
			string id = resi.getResidueName();
			int numRots = pos.getTotalNumberOfRotamers(j);
			//cout << id << ": " << numRots << endl;
			endNum = endNum + numRots;
			currIdMap.insert(make_pair(j, make_pair(startNum, endNum)));
			//cout << _interfacialPositionsList[i] << ":" << j << "," << id << ": " << startNum << "; " << endNum << endl;
			startNum = startNum + numRots;
		}
		posRotLimitMap.insert(make_pair(_interfacialPositionsList[i], currIdMap));
	}
	return posRotLimitMap;
}

//makes the best state applicable for unlinked positions by duplicating the rotamer at each interfacial position on the opposite chain
void unlinkBestState(Options &_opt, vector<uint> &_bestState, vector<int> _rotamerSampling, int _backboneLength){
	vector<int> linkedPositions = getLinked(_rotamerSampling, _backboneLength, _opt.interfaceLevel, _opt.sasaRepackLevel.size()-1);
	// Inserts the current rotamer at interfacial positions at positions on second chain
	//for (uint i=0; i<linkedPositions.size(); i++){
	//	cout << linkedPositions[i] << ",";
	//}
	//cout << endl;

	for (uint i=_backboneLength; i<_backboneLength*2; i++){
		vector<int>::iterator itr;
		itr = find(linkedPositions.begin(), linkedPositions.end(), i-_backboneLength);
		if (itr != linkedPositions.end()){
			vector<uint>::iterator itPos;
			itPos = _bestState.begin()+i;
			_bestState.insert(itPos, _bestState[i-_backboneLength]);
		}
	}
	//for (uint i=0; i<_bestState.size(); i++){
	//	cout << _bestState[i] << ",";
	//}
}

bool convertStateMapToSequenceMap(System &_sys, vector<pair<double,vector<uint>>> &_energyStateVec, map<vector<uint>, map<string,double>> &_stateEnergyMap, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair){
	//convert stateEnergyMap to sequenceEnergyMap
	map<vector<uint>,map<string,double>>::iterator itr;
	bool nonClashing = false;
	Chain & chain = _sys.getChain("A");
	// Setup sequence and state pair vector
	cout << "Sequence               Energy No IMM1" << endl;
	for (uint i=0; i<_energyStateVec.size(); i++){
		vector<uint> state = _energyStateVec[i].second;
		_sys.setActiveRotamers(state);
		string currSeq = convertPolymerSeqToOneLetterSeq(chain);
		_sequenceEnergyMap[currSeq] = _stateEnergyMap.at(state);
		_sequenceStatePair.push_back(make_pair(currSeq, state));
		cout << currSeq << ": " << _sequenceEnergyMap[currSeq]["DimerNoIMM1"] << endl;
		if (_sequenceEnergyMap[currSeq]["VDWDimerNoIMM1"] < 0){
			nonClashing = true;
		}
	}
	return nonClashing;
}

/***********************************
 *sequence entropy functions
 ***********************************/
//while loop to reset the sequence and make another mutation
map<string,int> getAACountMap(vector<string> _seq){
	map<string,int> AAcounts;
	for (uint i=0; i<_seq.size(); i++){
		try{
			if (AAcounts.count(_seq[i]) > 0){
				AAcounts.at(_seq[i])++;
			} else {
				AAcounts[_seq[i]] = 1;
			}
		}
		catch(const out_of_range& e){
			continue;
		}
	}
	return AAcounts;
}

double calcNumberOfPermutations(map<string,int> _seqAACounts, int _seqLength){
	//This function calculates number of permutations using following equation: n!/(n-r)! where n is number of positions and r is number of AAs	
	double numPermutation = 1;
	double permutationDenominator = 1;
	for(uint i=_seqLength; i>1; i--){
		numPermutation = numPermutation*i;
	}
	map<string,int>::iterator itr;
	for(itr = _seqAACounts.begin(); itr != _seqAACounts.end(); itr++){
		for (uint j=itr->second; j>1; j--){
			permutationDenominator = permutationDenominator*j;
		}
	}
	numPermutation = numPermutation/permutationDenominator;

	return numPermutation;
}

void internalAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, int _seqLength){
	vector<string> seqVector;
	for (uint i=4; i<_seqLength-4; i++){
		stringstream tmp;
		tmp << _seq[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, _seqLength-8);
}

void sequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, int _seqLength){
	vector<string> seqVector;
	for (uint i=0; i<_seqLength; i++){
		stringstream tmp;
		tmp << _seq[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, _seqLength);
}

double calculateSequenceProbability(map<string,int> &_seqCountMap, map<string,double> &_entropyMap, double _numberOfPermutations){
	//Find AA in count map, get counts of AA, then calculate the sequence probability by multiplying each AAs membrane probability contribution
	double seqProb = 1;
	map<string,int>::iterator it;
	for (it=_seqCountMap.begin(); it != _seqCountMap.end(); it++){
		double memProb = _entropyMap.at(it->first);
		int count = it->second;
		seqProb = seqProb*(pow(memProb, count));
	}
	seqProb = seqProb*_numberOfPermutations;
	return seqProb;
}

//TODO: maybe just read the sequence entropy map here? 
void calculateSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq, map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_internalSEProb, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal){
	map<string,int> prevSeqCountMap;
	map<string,int> currSeqCountMap;
	map<string,int> internalSeqCountMap;
	double prevNumberOfPermutations;
	double currNumberOfPermutations;
	double internalNumberOfPermutations;
	
	//double prevInternalSEProb;
	//Calculate the sequence entropy for both sequences
	sequenceEntropySetup(_prevSeq, prevSeqCountMap, prevNumberOfPermutations, _opt.backboneLength);
	sequenceEntropySetup(_currSeq, currSeqCountMap, currNumberOfPermutations, _opt.backboneLength);
	internalAASequenceEntropySetup(_currSeq, internalSeqCountMap, internalNumberOfPermutations, _opt.backboneLength);
	_prevSEProb = calculateSequenceProbability(prevSeqCountMap, _entropyMap, prevNumberOfPermutations);
	_currSEProb = calculateSequenceProbability(currSeqCountMap, _entropyMap, currNumberOfPermutations);
	_internalSEProb = calculateSequenceProbability(internalSeqCountMap, _entropyMap, internalNumberOfPermutations);

	//Calculate the probability of each sequence compared to the other
	//double totSEProb = _prevSEProb+_currSEProb;
	//On 10_26_21: after changing to the 2021_10_25_seqEntropies.txt file with more accurate membrane composition from Liu et al. 2002, I started seeing a lot of GxxxG sequences. I also checked my runs from before this date, and there were qute a bit of GxxxG sequences as well. Uncomment the below if you would like to try running using the internal probabilities, which should give a bit more flexibiilty to choosing AAs with lower likelilhoods of being in the membrane
	map<string,int> prevInternalSeqCountMap;
	double prevInternalNumberOfPermutations;
	internalAASequenceEntropySetup(_prevSeq, prevInternalSeqCountMap, prevInternalNumberOfPermutations, _opt.backboneLength);
	double prevInternalSEProb = calculateSequenceProbability(prevInternalSeqCountMap, _entropyMap, prevInternalNumberOfPermutations);
	double totSEProb = prevInternalSEProb+_internalSEProb;
	//double prevSeqProp = _prevSEProb/totSEProb;
	//double currSeqProp = _currSEProb/totSEProb;
	double prevSeqProp = prevInternalSEProb/totSEProb;
	double currSeqProp = _internalSEProb/totSEProb;
	

	//Convert the probability of each sequence to an entropy term that can be added to the original energy to directly compare two sequences
	double prevEner = -log(prevSeqProp)*0.592*_opt.weight_seqEntropy;//Multiplies the energy (log proportion times RT in kcal/mol) by the sequence entropy weight (weight of 1 gives entropy same weight as other terms)
	double currEner = -log(currSeqProp)*0.592*_opt.weight_seqEntropy;
	
	//Calculate energy total for best sequence vs current sequence
	//The below includes the baseline energy, which is an estimate of monomer energy
	_bestEnergyTotal = _bestEnergy+prevEner;
	_currEnergyTotal = _currEnergy+currEner;//TODO: I would like to include these in analysis somehow, but can't thikn of a way to do so

	//Output the terms if verbose
	if (_opt.verbose){
		cout << "Prev Prob:    " << _prevSEProb << endl;
		cout << "New Prob:     " << _currSEProb << endl;
		cout << "Prev Seq Proportion: " << prevSeqProp << endl;
		cout << "New Seq Proportion:  " << currSeqProp << endl;
		cout << "PrevEner =    " << prevEner << endl;
		cout << "NewEner =     " << currEner << endl;
		cout << "Diff =        " << (prevEner-currEner) << endl;
		cout << "Best Energy: " << _bestEnergyTotal << endl;
		cout << "New Energy: " << _currEnergyTotal << endl;
	}
}

// other random functions
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

void checkIfAtomsAreBuilt(System &_sys, ofstream &_err){
	for (uint i=0; i<_sys.atomSize(); i++){
		Atom atom = _sys.getAtom(i);
		if (!atom.hasCoor()){
			_err << "Atom " << i << " was not assigned coordinates; program termination";
			cout << "Atom " << i << " was not assigned coordinates; program termination";
			break;
		} else {
			continue;
		}
	}
}

void addSequencesToVector(vector<pair<double,string>> &_energyVector, vector<string> &_allSeqs){
	for (uint i=0; i<_energyVector.size(); i++){
		//cout << i << ": " << _energyVector[i].first << " " << _energyVector[i].second << endl;
		sameSequenceChecker(_energyVector[i].second, _allSeqs);
	}
}


/***********************************
 *energy builders
 ***********************************/
void buildSelfInteractions(System &_sys, map<string, double> &_selfMap){
	EnergySet* ESet = _sys.getEnergySet();

	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);
				string baseId = res.getResidueName();
				if (p-positions.begin() < 4 || p-positions.begin() > positions.size()-5){//On 03_18_2021 I found this error; position.size() is weird, so need to use 5 instead of 4
					baseId = baseId.append("-OUT");
				}
				try{
					double ener = _selfMap.at(baseId);
					Atom *a = &res.getAtom("CA");
					ESet->addInteraction(new BaselineInteraction(*a,ener));
				}
				catch (const out_of_range& e){
					continue;		
				}
			}
		}
	}
}

void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	EnergySet* ESet = _sys.getEnergySet();
	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j < (*p)->identitySize(); j++){
				Residue &res1 = (*p)->getIdentity(j);
				string baseId1 = res1.getResidueName();
				if (p-positions.begin() < 4){
					baseId1 = baseId1.append("-ACE");
				}
				if (p-positions.begin() > positions.size()-5){//
					baseId1 = baseId1.append("-CT2");
				}
				//cout << "Identity " << j << ": " << baseId1 << endl;
				for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					//cout << "Position 2: " << p2-positions.begin() << endl;
					if (d <= 10){
						//cout << "Distance: " << d << endl;
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							if (p2-positions.begin() < 4){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-5){
								baseId2 = baseId2.append("-CT2");
							}
							try{
								map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
								map<uint,double> AA2 = AA1.at(baseId2);
								double ener = AA2.at(d);
								Atom *a = &res1.getAtom("CA");
								Atom *b = &res2.getAtom("CA");
								ESet->addInteraction(new BaselinePairInteraction(*a,*b,-1*ener));//I forgot that this needs to be the opposite sign to actually counteract the energies of vdW and hydrogen bonding; switched it here but should eventually just switch in my baseline parameter file
							}
							catch (const out_of_range& e){
								continue;		
							}
						}
					}
				}
			}
		}
	}
}

void buildSequenceEntropy(System &_sys, map<string, double> &_sequenceEntropyMap, double _weight){
	EnergySet* ESet = _sys.getEnergySet();

	for(uint i = 0; i < 1; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		//cout << "Chain: " << thisChain.toString() << endl;
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			//cout << "Position: " << p-positions.begin() << endl;
			for (uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);//TODO: figure out if this is an error from 12_4_2020; if so, none of the runs were working for a reason...? changed it from i to j which I think it should be
				string baseId = res.getResidueName();
				//cout << "Identity " << j << ": " << baseId << endl;
				try{
					//cout << baseId << " worked" << endl;
					double ener = _sequenceEntropyMap.at(baseId);
					ener = -ener*(log2(ener))*_weight;
					//cout << "Energy: " << ener << endl;
					Atom *a = &res.getAtom("CA");
					ESet->addInteraction(new BaselineSequenceEntropy(*a,ener));
				}
				catch (const out_of_range& e){
					continue;
				}
			}
		}
	}
}

/****************************************
 *  
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults){
	
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

	//opt.allowed.push_back("");
	// optional
	opt.allowed.push_back("getGeoFromPDBData");
	
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("transform");

	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	
	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaX");

	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	opt.allowed.push_back("weight_seqEntropy");

	//Rotamers
	opt.allowed.push_back("SL");
	opt.allowed.push_back("SLInterface");
	
	//
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	// Cutoffs
	opt.allowed.push_back("dockingEnergyCutoff");
	opt.allowed.push_back("vdWEnergyCutoff");
	
	opt.allowed.push_back("printAllCrds");
	opt.allowed.push_back("printAxes");
	opt.allowed.push_back("printTermEnergies");
	opt.allowed.push_back("deleteTerminalHbonds");
	opt.allowed.push_back("linkInterfacialPositions");

	//Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("geometryDensityFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("selfEnergyFile");
	opt.allowed.push_back("pairEnergyFile");
	opt.allowed.push_back("seqEntropyFile");
	opt.allowed.push_back("AACompositionPenaltyFile");
	opt.allowed.push_back("configfile");
	
	opt.allowed.push_back("thread");
	opt.allowed.push_back("bbThread");

	//Alternate 
	opt.allowed.push_back("Ids");

	//Command Line Arguments
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("useIMM1");

	//MonteCarlo Arguments
	opt.allowed.push_back("numStatesToSave");

	//RNG Arguments
	opt.allowed.push_back("useTimeBasedSeed");
	
	opt.allowed.push_back("energyLandscape");
	opt.allowed.push_back("useBaselineComposition");
	
	//SelfPairManager Arguments
	opt.allowed.push_back("runDEESingles");
	opt.allowed.push_back("runDEEPairs");
	opt.allowed.push_back("runSCMF");
	
	//Energy Terms to Output
	opt.allowed.push_back("monomerEnergyTerms");
	opt.allowed.push_back("monomerIMM1EnergyTerms");
	opt.allowed.push_back("dimerEnergyTerms");
	opt.allowed.push_back("energyLandscapeTerms");
	opt.allowed.push_back("energyTermsToOutput");

	opt.allowed.push_back("energyTermList");

	//Load Rotamers from SASA values (from sgfc)
	opt.allowed.push_back("useSasa");
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("interfaceLevel");

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
		opt.OPerrors = OP.getErrors();
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
		}
	}
	
	opt.getGeoFromPDBData = OP.getBool("getGeoFromPDBData");
	if (OP.fail()) {
		opt.warningMessages += "getGeoFromPDBData not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = false;
	}
	
	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()) {
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}

	opt.useIMM1 = OP.getBool("useIMM1");
	if (OP.fail()) {
		opt.warningMessages += "useIMM1 not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.useIMM1 = false;
	}
	
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}

	opt.sequence = OP.getString("sequence");
	if(OP.fail()) {
		opt.warningMessages += "sequence not specified using L\n";
		opt.warningFlag = true;
		opt.sequence = "L";
	}

	opt.sequenceLength = OP.getInt("sequenceLength");
	if(OP.fail()) {
		opt.warningMessages += "sequenceLength not specified using 21\n";
		opt.warningFlag = true;
		opt.sequenceLength = 21;
	}

	
	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}

	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.sequenceLength;
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.backboneLength) + "\n";
		opt.warningFlag = true;
	}

	opt.sequenceStart = OP.getInt("sequenceStart");
	if (OP.fail()) {
		opt.warningMessages += "sequenceStart not specified using 1\n";
		opt.warningFlag = true;
		opt.sequenceStart = 1;
	}

	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(opt.tmStart) + "\n";
		opt.warningFlag = true;
		opt.startResNum = opt.tmStart;
	}

	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.tmEnd;
	}

	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting to 9.2\n";
		opt.warningFlag = true;
		opt.xShift = 9.2;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to 2\n";
		opt.warningFlag = true;
		opt.zShift = 2;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 40\n";
		opt.warningFlag = true;
		opt.axialRotation = 40;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to 25\n";
		opt.warningFlag = true;
		opt.crossingAngle = 25;
	}

	opt.transform = OP.getBool("transform");
	if (OP.fail()) {
		opt.warningMessages += "transform not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.transform = false;
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.thread = 0;
	}
	if (opt.thread == 0) {
		if (opt.sequenceLength > opt.backboneLength){
			opt.bbThread = opt.sequenceLength - opt.backboneLength + 1;
		} else {
			opt.bbThread = opt.backboneLength - opt.sequenceLength + 1;
		}
	}

	//Load Rotamers using SASA values (from sgfc)
	opt.useSasa = OP.getBool("useSasa");
	if (OP.fail()) {
		opt.warningMessages += "useSasa not specified, default true\n";
		opt.warningFlag = true;
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

	//Monte Carlo variables
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.errorMessages += "Number of MC cycles not specified!\n";
		opt.errorFlag = true;
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

	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.1;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 1.0;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 1.0;
	}
	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaX = 0.1;
	}

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.linkInterfacialPositions = OP.getBool("linkInterfacialPositions");
	if (OP.fail()) {
		opt.warningMessages += "linkInterfacialPositions not specified using true for less memory intensive version of code\n";
		opt.warningFlag = true;
		opt.linkInterfacialPositions = true;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
	}
	
	opt.repackEnergyCutoff = OP.getDouble("repackEnergyCutoff");
	if (OP.fail()) {
		opt.warningMessages += "repackEnergyCutoff not specified using 5.0\n";
		opt.warningFlag = true;
		opt.repackEnergyCutoff = 5.0;
	}
	opt.vdwEnergyCutoff = OP.getDouble("vdwEnergyCutoff");
	if (OP.fail()) {
		opt.warningMessages += "vdwEnergyCutoff not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.vdwEnergyCutoff = 0;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}

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
	opt.weight_seqEntropy = opt.weight_seqEntropy*100;//Default allows 1 to be weighted equally to other energy terms (I should convert to actual weighting conversion used with other energy terms)

	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}
	opt.SLInterface = OP.getString("SLInterface");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	} else {
		opt.SLInterface = "SL"+opt.SLInterface;
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
		opt.warningMessages += "backboneLength not specified, default to 35\n";
		opt.backboneLength = 35;
	}

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

	opt.infile = OP.getString("infile");
	if (OP.fail()) { 
		opt.warningMessages += "infile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.infile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
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
	opt.seqEntropyFile = OP.getString("seqEntropyFile");
	if (OP.fail()) { 
		opt.warningMessages += "seqEntropyFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.seqEntropyFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/seqEntropies.txt";
	}
	opt.AACompositionPenaltyFile = OP.getString("AACompositionPenaltyFile");
	if (OP.fail()) { 
		opt.warningMessages += "AACompositionPenaltyFile not specified, default \n";
		opt.warningFlag = true;
		opt.AACompositionPenaltyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/AACompositionPenalties.out";
	}

	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine pdbOutputDir";
		opt.errorFlag = true;
	}

	opt.Ids = OP.getStringVector("Ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}

	//MonteCarlo Options
	opt.numStatesToSave = OP.getInt("numStatesToSave");
	if (OP.fail()){
		opt.errorMessages += "numStatesToSave not specified, defaulting to 5";
		opt.warningFlag = true;
		opt.numStatesToSave = 5;
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

	opt.useTimeBasedSeed = OP.getBool("useTimeBasedSeed");
	if (OP.fail()) {
		opt.warningMessages += "useTimeBasedSeed not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useTimeBasedSeed = false;
	}

	opt.energyLandscape = OP.getBool("energyLandscape");
	if (OP.fail()) {
		opt.warningMessages += "energyLandscape not specified, defaulting to false";
		opt.warningFlag = true;
		opt.energyLandscape = false;
	}

	opt.useBaselineComposition = OP.getBool("useBaselineComposition");
	if (OP.fail()) {
		opt.warningMessages += "useBaselineComposition not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useBaselineComposition = false;
	}

	//Energy Terms to Output
	opt.monomerEnergyTerms = OP.getStringVector("monomerEnergyTerms");
	if (OP.fail()) {
		opt.monomerEnergyTerms.push_back("Monomer");
		opt.monomerEnergyTerms.push_back("VDWMonomer");
		opt.monomerEnergyTerms.push_back("HbondMonomer");
		opt.monomerEnergyTerms.push_back("MonomerSelfBaseline");
		opt.monomerEnergyTerms.push_back("MonomerPairBaseline");
	}
	opt.monomerIMM1EnergyTerms = OP.getStringVector("monomerIMM1EnergyTerms");
	if (OP.fail()) {
		opt.monomerIMM1EnergyTerms.push_back("Monomerw/IMM1");
		opt.monomerIMM1EnergyTerms.push_back("VDWMonomerw/IMM1");
		opt.monomerIMM1EnergyTerms.push_back("HbondMonomerw/IMM1");
		opt.monomerIMM1EnergyTerms.push_back("IMM1Monomer");
	}
	opt.dimerEnergyTerms = OP.getStringVector("dimerEnergyTerms");
	if (OP.fail()) {
		opt.dimerEnergyTerms.push_back("Dimer");
		opt.dimerEnergyTerms.push_back("HbondDimer");
		opt.dimerEnergyTerms.push_back("VDWDimer");
		opt.dimerEnergyTerms.push_back("IMM1Dimer");
	}
	opt.energyLandscapeTerms = OP.getStringVector("energyLandscapeTerms");
	if (OP.fail()) {
		opt.energyLandscapeTerms.push_back("EnergyBeforeLocalMC");
		opt.energyLandscapeTerms.push_back("DimerNoIMM1");
		opt.energyLandscapeTerms.push_back("Baseline");
		opt.energyLandscapeTerms.push_back("VDWDimerNoIMM1");
		opt.energyLandscapeTerms.push_back("HBONDDimerNoIMM1");
		opt.energyLandscapeTerms.push_back("EnergyBeforeLocalMCw/seqEntropy");
	}
	opt.energyTermsToOutput = OP.getStringVector("energyTermsToOutput");
	if (OP.fail()) {
		opt.energyTermsToOutput.push_back("Total");
		opt.energyTermsToOutput.push_back("Dimer");
		opt.energyTermsToOutput.push_back("Monomer");
		opt.energyTermsToOutput.push_back("VDWDimer");
		opt.energyTermsToOutput.push_back("VDWMonomer");
		opt.energyTermsToOutput.push_back("VDWDiff");
		opt.energyTermsToOutput.push_back("HBONDDimer");
		opt.energyTermsToOutput.push_back("HBONDMonomer");
		opt.energyTermsToOutput.push_back("HBONDDiff");
		opt.energyTermsToOutput.push_back("IMM1Dimer");
		opt.energyTermsToOutput.push_back("IMM1Monomer");
		opt.energyTermsToOutput.push_back("IMM1Diff");
		opt.energyTermsToOutput.push_back("MonomerNoIMM1");
		opt.energyTermsToOutput.push_back("DimerNoIMM1");
		opt.energyTermsToOutput.push_back("Baseline");
		opt.energyTermsToOutput.push_back("Baseline-Monomer");
		opt.energyTermsToOutput.push_back("VDWDimerNoIMM1");
		opt.energyTermsToOutput.push_back("VDWMonomerNoIMM1");
		opt.energyTermsToOutput.push_back("HBONDDimerNOIMM1");
		opt.energyTermsToOutput.push_back("HBONDMonomerNOIMM1");
		opt.energyTermsToOutput.push_back("DimerSelfBaseline");
		opt.energyTermsToOutput.push_back("DimerPairBaseline");
		opt.energyTermsToOutput.push_back("DesignNumber");
	}
	opt.energyTermList = OP.getStringVector("energyTermList");
	if (OP.fail()) {
		//This works, but I think if you ever want to output more terms in the future, need to add them to the terms above
		//TODO: write in an error that will tell you if the above is the case
		opt.energyTermList.push_back("CHARMM_VDW");
		opt.energyTermList.push_back("SCWRL4_HBOND");
		opt.energyTermList.push_back("CHARMM_IMM1");
		opt.energyTermList.push_back("CHARMM_IMM1REF");
	}
	
	opt.rerunConf = OP.getConfFile();

	return opt;
}
