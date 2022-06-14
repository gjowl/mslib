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
#include "BaselineEnergyBuilder.h"
#include "BaselineInteraction.h"
#include "BaselineIMM1Interaction.h"
#include "BaselinePairInteraction.h"
#include "BaselineOuterPairInteraction.h"
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"
#include "SasaCalculator.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "geomRepack";
string programDescription = "This pairs with seqDesign: takes designed sequences at their given geometry and does a local geometric repack on them";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "19 August 2021";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************
 *        ======= OPTIONS =======
 ******************************************/
struct Options{//TODO: divide these into necessary and non-necessary
	// input files
	string topFile;
	string parFile;
	string solvFile;
	string hbondFile;
	string rotLibFile;
	string designPdb;
	string designCrd;
	string backboneCrd;
	string outputDir;
	string sequenceEntropyFile;
	string selfEnergyFile;
	string pairEnergyFile;

	// Geometry
	double xShift;
	double crossingAngle;
	double axialRotation;
	double zShift;
	int thread;

	// tm
	int tmStart;
	int tmEnd;

	// Design parameters
	string sequence;
	string rotamerSamplingString;
	vector<uint> state;
	vector<int> rotamerSamplingVector;

	// load rotamers from SASA values (from sgfc)
	bool keepOriginalRotamer;
	std::vector<string> sasaRepackLevel;

	// Monte Carlo Parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;
	int numRepacks;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	double monomer;
	double monoVdw;
	double monoHbond;
	double monoIMM1;
	
	// input previous dimerEnergy
	double dimer;
	double dimerVdw;
	double dimerHbond;
	double dimerIMM1;
	double total;
	
	// Shift Size
	double deltaX; 
	double deltaCross;
	double deltaAx;
	double deltaZ;

	// Other options
	bool verbose;
	int greedyCycles;
	int seed;
	int interfaceLevel;
	string runNumber;
	bool estimateMonomerWithBaseline;

	vector<string> energyTermsToOutput;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help //figure out how to add stuff to help

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

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
	bool useIMM1;
};

//TODO: edit this
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}
/******************************************
 *        ======= FUNCTIONS =======
 ******************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt);
vector<pair<int,double>> calculateResidueBurial (System &_sys);
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
double getStandardNormal(RandomNumberGenerator& RNG);

void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {

	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());
	//fout << x << " " << y << " " << pt << " " << caApV.size() << endl;

	// old code
	//for(int i = 0; i < _apV.size(); i++) {
	//	CartesianPoint& pt = _apV[i]->getCoor();
	//	pt.setZ(pt.getZ() +  zShift);
	//}

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, zShift);
	_trans.translate(_apV, interDistVect);
	_trans.translate(_axis, interDistVect);
	

}

string generatePolymerSequenceFromSequence(string _sequence, int _startResNum) {
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
	return "A" + ps + "\nB" + ps;
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

string mutateSequence(string _sequence, string _AA, uint _posToSwitch){
	stringstream ss;
	for (uint i=0; i<_sequence.length(); i++){
		if (i != _posToSwitch){
			ss << _sequence[i];
		} else {
			ss << _AA;
		}
	}
	string mutantSeq = ss.str();
	return mutantSeq;
}

vector<string> getMutantInterfacialSequences(Options &_opt, vector<string> &_AAList){
	vector<string> mutantSeqList;
	for (uint i=3; i<_opt.rotamerSamplingString.length()-4; i++){
		char tmp1 = _opt.rotamerSamplingString[i];
		stringstream ss1;
		ss1 << tmp1;
		string AA1 = ss1.str();
		if (AA1 == "0"){
			uint position = i;
			char tmp = _opt.sequence[i];
			stringstream ss;
			ss << tmp;
			string AA = ss.str();
			for (uint j=0; j<_AAList.size(); j++){
				if (AA != _AAList[j]){
					string currSeq = mutateSequence(_opt.sequence, _AAList[j], position);
					mutantSeqList.push_back(currSeq);
					cout << currSeq << endl;
				}
			}
		}
	}
	return mutantSeqList;
}

//TODO: make this into a multithread; add in a way to take the repack sequence

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

void interfaceSequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<int> _interfacialPositionsList){
	//Get residue name for each interfacial identity
	vector<string> seqVector;
	int numInterfacials = _interfacialPositionsList.size()-1;
	for (uint i=0; i<numInterfacials; i++){
		//Position &pos = _sys.getPosition(_interfacialPositionsList[i]);
		//Residue &resi = pos.getCurrentIdentity();
		//string id = resi.getResidueName();
		stringstream tmp;
		tmp << _seq[_interfacialPositionsList[i]];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
		//cout << _interfacialPositionsList[i] << " : " << resName << endl;
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, numInterfacials);
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
	//cout << "pSeq = " << seqProb << endl;
	seqProb = seqProb*_numberOfPermutations;
	//cout << "comb = " << _numberOfPermutations << endl;
	//cout << "pRes = " << seqProb << endl;
	return seqProb;
}


double getInterfaceSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap, vector<int> _interfacialPositionsList){
	map<string,int> AACountMap;
	double numberOfPermutations;
	interfaceSequenceEntropySetup(_sequence, AACountMap, numberOfPermutations, _interfacialPositionsList);
	double seqProb = calculateSequenceProbability(AACountMap, _entropyMap, numberOfPermutations);
	return seqProb;
}


void insertMonomerEnergiesIntoEnergyMap(Options &_opt, SelfPairManager &_spm, string _sequence, vector<unsigned int> &_stateVector, map<string, double> &_energyMap){
	_energyMap["Monomer"]      = _spm.getStateEnergy(_stateVector);
	_energyMap["VDWMonomer"]   = _spm.getStateEnergy(_stateVector, "CHARMM_VDW");
	_energyMap["HBONDMonomer"] = _spm.getStateEnergy(_stateVector, "SCWRL4_HBOND");
	_energyMap["IMM1Monomer"]  = _spm.getStateEnergy(_stateVector, "CHARMM_IMM1REF")+_spm.getStateEnergy(_stateVector, "CHARMM_IMM1");
}

vector<int> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<int> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=0; k<_opt.sequence.length(); k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

void insertGeometriesIntoMap(Options &_opt, map<string, double> &_energyMap, double _xShift, double _crossingAngle, double _axialRotation, double _zShift){
	_energyMap["startXShift"]         = _opt.xShift;
	_energyMap["startCrossingAngle"]  = _opt.crossingAngle;
	_energyMap["startAxialRotation"]  = _opt.axialRotation;
	_energyMap["startZShift"]         = _opt.zShift;
	_energyMap["xShift"]              = _xShift;
	_energyMap["crossingAngle"]       = _crossingAngle;
	_energyMap["axialRotation"]       = _axialRotation;
	_energyMap["zShift"]              = _zShift;
	_energyMap["xShiftDiff"]          = _opt.xShift-_xShift;
	_energyMap["crossingAngleDiff"]   = _opt.crossingAngle-_crossingAngle;
	_energyMap["axialRotationDiff"]   = _opt.axialRotation-_axialRotation;
	_energyMap["zShiftDiff"]          = _opt.zShift-_zShift; 
}

void insertDimerEnergiesIntoMap(Options &_opt, SelfPairManager &_spm, string _sequence, vector<unsigned int> &_stateVector, map<string, double> &_energyMap, map<string, double> &_entropyMap){
	_energyMap["Dimer"]      = _spm.getStateEnergy(_stateVector);
	_energyMap["VDWDimer"]   = _spm.getStateEnergy(_stateVector, "CHARMM_VDW");
	_energyMap["HBONDDimer"] = _spm.getStateEnergy(_stateVector, "SCWRL4_HBOND");
	_energyMap["IMM1Dimer"]  = _spm.getStateEnergy(_stateVector, "CHARMM_IMM1REF")+_spm.getStateEnergy(_stateVector, "CHARMM_IMM1");
	
	vector<int> interfacePositions = getAllInterfacePositions(_opt, _opt.rotamerSamplingVector);
	double interfaceSEProb = getInterfaceSequenceEntropyProbability(_opt, _sequence, _entropyMap, interfacePositions);
	_energyMap["SequenceProbability"] = interfaceSEProb;
}

void getStartingStateEnergies(Options &_opt, SelfPairManager &_spm, vector<unsigned int> &_stateVector, map<string, double> &_energyMap, map<string, double> &_entropyMap){
	insertDimerEnergiesIntoMap(_opt, _spm, _opt.sequence, _stateVector, _energyMap, _entropyMap);
	_energyMap["Monomer"]      = _opt.monomer;
	_energyMap["Total"]        = _energyMap["Dimer"]-_opt.monomer;
	_energyMap["DimerDiff"]    = _energyMap["Dimer"]-_opt.dimer;
	_energyMap["VDWDiff"]      = _energyMap["VDWDimer"]-_opt.dimerVdw;
	_energyMap["HBONDDiff"]    = _energyMap["HBONDDimer"]-_opt.dimerHbond;
	_energyMap["IMM1Diff"]     = _energyMap["IMM1Dimer"]-_opt.dimerIMM1;
}

void repackMutantSequence(Options &_opt, vector<string> &_acceptedSeqs, string _sequence, AtomPointerVector &_designApv, map<string,map<string,double>> &_energyMap, double _bestRepackEnergy, double _monomerEnergy, double _xShift, double _crossingAngle, double _axialRotation, double _zShift){
	/******************************************************************************
	 *                   === INITIALIZE POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	//string polySeq = generatePolymerSequenceFromSequence(_opt.sequence, _opt.thread);
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
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	//CRDReader cRead;
	//cRead.open(_crdFile); 
	bool noCrdFile = false;
	AtomPointerVector designApv;
	//cRead.close();
	
	System pdb;
	pdb.readPdb(_opt.designPdb);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
	Chain & chainA1 = pdb.getChain("A");
	Chain & chainB1 = pdb.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA1 = chainA1.getAtomPointers();
	AtomPointerVector & apvChainB1 = chainB1.getAtomPointers();
	
	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA1, apvChainB1, axisA, axisB, ori, xAxis, zAxis, _zShift, _axialRotation, _crossingAngle, _xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	if (noCrdFile == true){
		designApv = pdb.getAtomPointers();
	}
	cout << designApv << endl;

	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile,_opt.solvFile);
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
	//CSB.buildSystemFromPDB(designApv);
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from pdb" << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(_designApv,false);
	sys.buildAllAtoms();
	
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                          === SETUP ENERGY SET ===
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
	Eset->setTermActive("CHARMM_IMM1", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	
	/******************************************************************************
	 *             === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// Assign number of rotamers by residue burial
	loadRotamersBySASABurial(sys, sysRot, _opt);
	CSB.updateNonBonded(10,12,50);

	//slightly convoluted, but something with the helical axis just wouldn't get correct if I didn't tranform it, but if I did it messed with my structure...so I worked around above. I don't think it should affect my results in anyway other than putting my helices in the membrane, and it's consistent with what I did for the results in my design code

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 

	// Optimize Initial Starting Position (using Baseline to get back to original result)
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
	string sequence = _sequence;
	//TODO: calculate energy here and suctract from monomer
	
	//Repack dimer
	repackSideChains(spm, 10);
	vector<uint> startStateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(startStateVec);
	double currentEnergy = spm.getStateEnergy(startStateVec);
	sys.setActiveRotamers(startStateVec);

	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	double bestEnergy = currentEnergy;
	cout << "Monomer     = " << _monomerEnergy << endl;
	cout << "Dimer       = " << currentEnergy << endl;
	cout << "CurrentBest = " << bestEnergy-_monomerEnergy << endl;

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/

	cout << "Starting Geometry" << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << _opt.xShift << " crossingAngle: " << _opt.crossingAngle << " axialRotation: " << _opt.axialRotation << " zShift: " << _opt.zShift << endl << endl;
	cout << "Current Best Energy: " << bestEnergy-_monomerEnergy << endl;
	cout << "Interaction Energies: " << endl;
	cout << spm.getSummary(startStateVec) << endl;
	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	for (uint i=0; i<_opt.numRepacks; i++){
		double prevBestEnergy = bestEnergy;
		sys.applySavedCoor("startingState");
		helicalAxis.applySavedCoor("startingAxis");
		sys.saveAltCoor("savedBestState");
		helicalAxis.saveAltCoor("BestAxis");
		
		double xShift = _xShift;
		double crossingAngle = _crossingAngle;
		double axialRotation = _axialRotation;
		double zShift = _zShift;
	
		if (_opt.verbose){
			cout << "======================================" << endl;
			cout << "Performing Local Monte Carlo Repack " << i << endl;
			cout << "======================================" << endl;
		}

		vector<unsigned int> MCOBest = startStateVec;
		
		//MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
		MonteCarloManager MCMngr(0.5, 0.5, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects, 10, 0.01);
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
				deltaZShift = getStandardNormal(RNG) * _opt.deltaZ;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
			} else if (moveToPreform == 1) {
			//===========================
			//===== Axial Rotation ======
			//===========================
				//deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
				deltaAxialRotation = getStandardNormal(RNG) * _opt.deltaAx;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
			} else if (moveToPreform == 2) {
			//==================================
			//====== Local Crossing Angle ======
			//==================================
				//deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
				deltaCrossingAngle = getStandardNormal(RNG) * _opt.deltaCross;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
			} else if (moveToPreform == 3) {
			//==============================================
			//====== X shift (Interhelical Distance) =======
			//==============================================
				//deltaXShift = getStandardNormal(RNG1) * 0.1;
				deltaXShift = getStandardNormal(RNG) * _opt.deltaX;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
			}
		
			// Run Optimization
			repackSideChains(spm, 10);

			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			currentEnergy = spm.getMinBound()[0];
		
			if (!MCMngr.accept(currentEnergy)) {
				if (_opt.verbose){
					cout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy-_monomerEnergy << endl;
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
		
				if (_opt.verbose){
					cout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-_monomerEnergy << endl;
				}
				counter++;
			}
		}
		
		sys.applySavedCoor("savedBestState");

		double dimerEnergy = spm.getStateEnergy(MCOBest);
		double dimerEnergyDiff = dimerEnergy-_opt.dimer;
		double finalEnergy = dimerEnergy-_monomerEnergy;
		double vdw = spm.getStateEnergy(MCOBest, "CHARMM_VDW");
		double hbond = spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
		double imm1 = spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");

		time(&endTimeMC);
		diffTimeMC = difftime (endTimeMC, startTimeMC);
		
		string tab = "\t";
		// Print out info to the summary file
		cout << "Geometry " << i << endl;
		cout << "xShift:               " << xShift << endl;
		cout << "crossingAngle:        " << crossingAngle << endl;
		cout << "axialRotation:        " << axialRotation << endl;
		cout << "zShift:               " << zShift << endl << endl;
		cout << "Monomer Energy:       " << _monomerEnergy << endl;
		cout << "Energy Before Repack: " << _opt.dimer << endl << endl;
		size_t spacerSize = 20;
		string spacer(spacerSize, ' ');
		//TODO: add in a way to save these so that I can output them at the end (or just only save these and not the above geometries?)
		cout << spacer + " Final Energy" << tab << "VDW" << tab << "HBOND" << tab << "IMM1" << tab << "startXShift" << tab << "startCrossingAngle" << tab << "startAxialRotation" << tab << "startZShift" << tab << "currxShift" << tab << "currCrossingAngle" << tab << "currAxialRotation" << tab << "currZShift" << tab << "thread" << tab << "Time" << endl;
		cout << "Final Repack Energy " << i << ": " << finalEnergy << tab << vdw << tab << hbond << tab << imm1 << tab << _xShift << tab << _crossingAngle << tab << _axialRotation << tab << _zShift << tab << xShift << tab << crossingAngle << tab << axialRotation << tab << zShift << _opt.thread  << tab << diffTimeMC << endl << endl;

		// Add energies of the sequence to the energy map
		//map<string,double> &energyMap = _seqEnergyMap[sequence];
		//outputEnergiesByTerm(spm, MCOBest, energyMap, _opt.energyTermList, "", 1);
		//saveEnergyDifference(_opt, _seqEnergyMap, sequence);
		//seqEnergyMap[sequence][""] = dimerEnergy;
		//seqEnergyMap[sequence]["Total"] = finalEnergy;
		if (finalEnergy < _bestRepackEnergy-5){
			cout << _sequence << " more stable: " << finalEnergy << endl;
			_acceptedSeqs.push_back(_sequence);
			PDBWriter designWriter;
			designWriter.open(_opt.outputDir + "/mutant_"+_sequence+".pdb");//names the repack by design (seed) and sequence number (designNumber)
			designWriter.write(sys.getAtomPointers(), true, false, true);
			designWriter.close();
		} else if (finalEnergy > _bestRepackEnergy+5){
			cout << _sequence << " less stable: " << finalEnergy << endl;
			_acceptedSeqs.push_back(_sequence);
			PDBWriter designWriter;
			designWriter.open(_opt.outputDir + "/mutant_"+_sequence+".pdb");//names the repack by design (seed) and sequence number (designNumber)
			designWriter.write(sys.getAtomPointers(), true, false, true);
			designWriter.close();
		}
	}

	//TODO: for below: get the apv for a crd file, get the file, and enter that info above to get that
}

double computeDimerEnergy(System &_sys, Options& _opt, map<string,double> &_energyMap, string &_sequence, RandomNumberGenerator &_RNG, PDBWriter &_writer, map<string, double> &_entropyMap, string _mutantOutputDir, ofstream &_sout, ofstream &_err) {

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
	loadRotamersBySASABurial(sys, sysRot, _opt);
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
	vector<uint> stateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(stateVec);
	
	double dimerEnergy = spm.getStateEnergy(stateVec);
	
	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator sasa(sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	PDBWriter designWriter;
	designWriter.open(_mutantOutputDir + "/"+_sequence+".pdb");//names the repack by design (seed) and sequence number (designNumber)
	designWriter.write(sys.getAtomPointers(), true, false, true);
	_writer.write(sys.getAtomPointers(), true, false, true);
	designWriter.close();
	
	insertDimerEnergiesIntoMap(_opt, spm, _sequence, stateVec, _energyMap, _entropyMap);
	_energyMap["DimerSasa"] = dimerSasa;
	return dimerEnergy;
}

void insertStateEnergiesIntoMap(string _sequence, map<string, map<string, double>> &_energyMap, SelfPairManager &_spm, vector<uint> &_state, vector<string> _energies){
	_energyMap[_sequence]["Total"] = _spm.getStateEnergy(_state);
	for (uint i=0; i<_energies.size(); i++){
		string energyTerm = _energies[i];
		double energy = _spm.getStateEnergy(_state, energyTerm);
		_energyMap[_sequence][energyTerm] = energy;
	}
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

void buildSelfInteractions(System &_sys, map<string, double> &_selfMap){
	EnergySet* ESet = _sys.getEnergySet();

	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);
				string baseId = res.getResidueName();
				if (p-positions.begin() < 3 || p-positions.begin() > positions.size()-5){//On 03_18_2021 I found this error; position.size() is weird, so need to use 5 instead of 4; on 11_20_2021 saw that a lot of clashing occurs at hte 4th position, so changed this to only use 3
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
				if (p-positions.begin() < 1){
					baseId1 = baseId1.append("-ACE");
				}
				//Changed this on 11_24_2021 for the designFiles/2021_11_22_IMM1Self and Pair baselines 
				//if (p-positions.begin() > positions.size()-5){//
				//	baseId1 = baseId1.append("-CT2");
				//}
				//cout << "Identity " << j << ": " << baseId1 << endl;
				for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					//cout << "Position 2: " << p2-positions.begin() << endl;
					if (d <= 10){
						//cout << "Distance: " << d << endl;
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							//if (p2-positions.begin() < 3){
							//	baseId2 = baseId2.append("-ACE");
							//}
							if (p2-positions.begin() > positions.size()-2){
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


double computeMonomerEnergy(Options& _opt, Transforms &_trans, map<string,double> &_energyMap, string _sequence, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {
	//string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, _opt.thread);//fixed monomer calculation issue on 05_12_2021
	string polySeq = generateMonomerPolymerSequenceFromSequence(_sequence, _opt.thread);
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

	double monomerEnergy = 0;
	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();
	
	if (_opt.estimateMonomerWithBaseline){
		map<string, double> selfMap = readSingleParameters(_opt.selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_opt.pairEnergyFile);
		buildSelfInteractions(monoSys, selfMap);
		buildPairInteractions(monoSys, pairMap);

		monomerEnergy = -(monoSys.calcEnergy()*2);
		_energyMap["Monomer"] = monomerEnergy;
		cout << "Monomer: " << monomerEnergy << endl;
	} else {
		/******************************************************************************
		 *                 === LOAD ROTAMERS AND HYDROGEN BONDING ===
		 ******************************************************************************/
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

		// Clear saved coordinates
		monoSys.clearSavedCoor("savedBestMonomer");
		monoSys.clearSavedCoor("bestZ");
		helicalAxis.clearSavedCoor("BestMonomerAxis");
		helicalAxis.clearSavedCoor("bestZ");
	
		insertMonomerEnergiesIntoEnergyMap(_opt, monoSpm, _sequence, stateVec, _energyMap);
	}
	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	monoSasa.calcSasa();
	double monomerSasa = monoSasa.getTotalSasa();
	double totalMonomerSasa = monomerSasa*2;

	_energyMap["MonomerSasa"] = totalMonomerSasa;
	return monomerEnergy;
}

vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=3; k<_opt.sequence.length()-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
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

void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns){
	ofstream eout;
	string eoutfile = _opt.outputDir + "/optimizationEnergyFile.csv";
	eout.open(eoutfile.c_str());
	string tab = "\t";

	// Outputs
	eout << tab;
	for (uint i=0; i<_opt.energyTermsToOutput.size(); i++){
		eout << _opt.energyTermsToOutput[i] << tab;
		cout << _opt.energyTermsToOutput[i] << tab;
	}
	eout << "Starting Sequence" << tab;
	eout << "Sequence" << tab;
	eout << "InterfaceSeq" << tab;
	eout << "backboneLength" << tab;
	eout << "runNumber" << tab;
	eout <<	"PDBPath" << endl; 
	
	cout << "Starting Sequence" << tab;
	cout << "Sequence" << tab;
	cout << "InterfaceSeq" << tab;
	cout << "backboneLength" << tab;
	cout << "runNumber" << tab;
	cout <<	"PDBPath" << endl; 

	for (uint i=0; i<_allDesigns.size() ; i++){
		eout << _allDesigns[i] << endl;
		cout << _allDesigns[i] << endl;
	}
	eout.close();
}


void outputDesignFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, vector<string> &_sequences, map<string,map<string,double>> _sequenceEnergyMap){
	// Setup vector to hold energy file lines
	vector<string> energyLines;

	// Setup the parameters for this specific run
	string tab = "\t";
	string startXShift = MslTools::doubleToString(_opt.xShift);
	string startCrossingAngle = MslTools::doubleToString(_opt.crossingAngle);
	string startAxialRotation = MslTools::doubleToString(_opt.axialRotation);
	string startZShift = MslTools::doubleToString(_opt.zShift);
	string thread = MslTools::intToString(_opt.thread);
	string repackLevels = MslTools::intToString(_opt.sasaRepackLevel.size());
	string interfaceLevel = MslTools::intToString(_opt.interfaceLevel);
	string bbLength = MslTools::intToString(_opt.sequence.length());
	//string runParameters = xShift+tab+crossingAngle+tab+axialRotation+tab+zShift+tab+thread+tab+repackLevels+tab+interfaceLevel;
	string runParameters = bbLength+tab+_opt.runNumber;
	// For loop to setup the energy file
	string startSequence = _sequences[0];
	for (uint i=0; i<_sequences.size(); i++){
		string sequence = _sequences[i];
		map<string,double> energyMap = _sequenceEnergyMap.at(sequence);
		// For adding in strings to a line for the energy file
		string tmp = "Backbone Optimization: ";
		for (uint j=0; j<_opt.energyTermsToOutput.size(); j++){
			string energyTerm = _opt.energyTermsToOutput[j];
			//cout << energyTerm << endl;
			double energy = energyMap.at(energyTerm);
			//cout << sequence << ": " << energyTerm << " = " << energy << endl;
			tmp.append(MslTools::doubleToString(energy));
			tmp.append(tab);
		}

		//Add in path to design PDB and make repack and docking configuration files
		string designDir = _opt.outputDir + "/"+ sequence + ".pdb";
	
		// Append other important features to the end energy files lines
		string interfaceSequence = getInterfaceSequence(_opt,_interface, sequence);
		tmp.append(startSequence);
		tmp.append(tab);
		tmp.append(sequence);
		tmp.append(tab);
		tmp.append(interfaceSequence);
		tmp.append(tab);
		tmp.append(runParameters);
		tmp.append(tab);
		tmp.append(designDir);
		energyLines.push_back(tmp);
		tmp.clear();
	}
	outputEnergyFile(_opt, _interface, energyLines);
}

void addSequenceToMap(Options &_opt, map<string,map<string,double>> &_sequenceEnergyMap, map<string,double> &_energyMap, double _finalEnergy, string _sequence, ofstream &_out){
	if (_finalEnergy > 0){
			cout << _sequence << " clashing: " << _finalEnergy << endl;
			_out << _sequence << " clashing: " << _finalEnergy << endl;
			_sequenceEnergyMap[_sequence] = _energyMap;
			_sequenceEnergyMap[_sequence]["Total"] = _finalEnergy;
	} else {
		cout << _sequence << " not clashing: " << _finalEnergy << endl;
		_out << _sequence << " not clashing: " << _finalEnergy << endl;
		_sequenceEnergyMap[_sequence] = _energyMap;
		_sequenceEnergyMap[_sequence]["Total"] = _finalEnergy;
	}
}




void setupOutputDirectory(Options &_opt){
	_opt.outputDir = _opt.outputDir + "/backboneOptimization";//Had problems with dynamic date because of runs starting one day and ending another
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}



void saveEnergyDifference(Options _opt, map<string,map<string,double>> &_sequenceEnergyMap, string _sequence){
	string startSequence = _opt.sequence;
	_sequenceEnergyMap[_sequence]["DimerDiff"] = _sequenceEnergyMap[_sequence]["Dimer"] - _sequenceEnergyMap[startSequence]["Dimer"];
	_sequenceEnergyMap[_sequence]["HBONDDiff"] = _sequenceEnergyMap[_sequence]["HBONDDimer"] - _sequenceEnergyMap[startSequence]["HBONDDimer"];
	_sequenceEnergyMap[_sequence]["VDWDiff"]   = _sequenceEnergyMap[_sequence]["VDWDimer"]   - _sequenceEnergyMap[startSequence]["VDWDimer"];
	_sequenceEnergyMap[_sequence]["IMM1Diff"]  = _sequenceEnergyMap[_sequence]["IMM1Dimer"]  - _sequenceEnergyMap[startSequence]["IMM1Dimer"];

}

/******************************************
 *        ======= BEGIN MAIN =======
 ******************************************/
int main(int argc, char *argv[]){
//TODO: 
//	1. input new config file for a sequence from seqDesign
//	2. transform to the appropriate geometry and state
//	3. check energy to see if it's the same
//	4.  save to the same directory as already originally saved in (make sure I think of a way to keep this continuous. It would be great if I could give my code an option so that it will atuo do this or not; and that this had an option to be connected directly to seqDesgin so that it can be used by lots of people (If I want to be able to setup a server for this, which I think I do:
//		-Make it so that people choose a density type for the kde files, choose the types of AAs they want to use; choose lengths of helices, choose using baseline or not, etc
//		-The tricky part is going to be making this also editable for lots of people if they ever want to edit the code if it breaks on something I'm not thinkin about...add comments in within the next month (it would actually be pretty incredible if I actually am able to say hat I've got a webserver up and running for this by December and also that I have data on the other part too...that would be really clutch and is definitely doable
//	5. after that all works, make a way to take all the config files and generate a submit file so that it's a pretty automatic process to run the seqDesign, then run the repack after (and save in the same folder)
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%Y_%m_%d",timeinfo);
	string date(buffer);
	
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);

	if (opt.errorFlag) {
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		usage();
		exit(1);
	}

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	setupOutputDirectory(opt);
	
	ofstream rerun;
	string rerunfile = opt.outputDir + "/optimizationRerun.config";
	rerun.open(rerunfile.c_str());
	rerun << opt.rerunConf << endl;
	rerun.close();
	
	ofstream sout;
	ofstream err;
	
	string soutfile = opt.outputDir + "/summary.out";
	string errfile  = opt.outputDir + "/error.err";

	sout.open(soutfile.c_str());
	err.open(errfile.c_str());

	sout << date << endl;
	err << date << endl;


	/******************************************************************************
	 *                   === INITIALIZE POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polySeq = convertToPolymerSequenceNeutralPatch(opt.sequence, opt.thread);
	//string polySeq = generatePolymerSequenceFromSequence(opt.sequence, opt.thread);
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
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.designCrd); 
	bool noCrdFile = false;
	AtomPointerVector designApv;
	if(!cRead.read()) {
		//cerr << "Unable to read " << opt.designCrd << endl;
		noCrdFile = true;
		//exit(0);
	} else {
		designApv = cRead.getAtomPointers();//*/
	}
	cRead.close();
	
	System pdb;
	pdb.readPdb(opt.designPdb);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
	Chain & chainA1 = pdb.getChain("A");
	Chain & chainB1 = pdb.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA1 = chainA1.getAtomPointers();
	AtomPointerVector & apvChainB1 = chainB1.getAtomPointers();
	
	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA1, apvChainB1, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	if (noCrdFile == true){
		designApv = pdb.getAtomPointers();
	}

	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
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
	//CSB.buildSystemFromPDB(designApv);
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from pdb" << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(designApv,false);
	sys.buildAllAtoms();
	
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                          === SETUP ENERGY SET ===
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
	Eset->setTermActive("CHARMM_IMM1", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	
	/******************************************************************************
	 *             === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,opt);

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

	// Optimize Initial Starting Position (using Baseline to get back to original result)
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
	string sequence = opt.sequence;
	double bestEnergy = opt.dimer;
	
	//Repack dimer
	repackSideChains(spm, 10);
	vector<uint> startStateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(startStateVec);
	double currentEnergy = spm.getStateEnergy(startStateVec);
	sys.setActiveRotamers(startStateVec);

	sys.saveAltCoor("startingState");
	helicalAxis.saveAltCoor("startingAxis");
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	double monomerEnergy = opt.monomer;

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
		

	PDBWriter designWriter;
	designWriter.open(opt.outputDir + "/optimizedBackbones.pdb");//names the repack by design (seed) and sequence number (designNumber)

	double bestRepackEnergy;
	vector<uint> bestRepackState;

	double bestXShift = opt.xShift;
	double bestAxialRotation = opt.axialRotation;
	double bestZShift = opt.zShift;
	double bestCrossingAngle = opt.crossingAngle;
	for (uint i=0; i<opt.numRepacks; i++){
		double prevBestEnergy = opt.dimer;
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
			cout << "Performing Local Monte Carlo Repack " << i << endl;
			cout << "======================================" << endl;
		}

		vector<unsigned int> MCOBest = startStateVec;
		
		//MonteCarloManager MCMngr(opt.MCStartTemp, opt.MCEndTemp, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects);
		MonteCarloManager MCMngr(0.5, 0.5, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects, 10, 0.01);
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
		
		sys.applySavedCoor("savedBestState");

		double dimerEnergy = spm.getStateEnergy(MCOBest);
		double dimerEnergyDiff = dimerEnergy-opt.dimer;
		double finalEnergy = dimerEnergy-monomerEnergy;
		double vdw = spm.getStateEnergy(MCOBest, "CHARMM_VDW");
		double hbond = spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
		double imm1 = spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");

		time(&endTimeMC);
		diffTimeMC = difftime (endTimeMC, startTimeMC);
		
		string tab = "\t";
		// Print out info to the summary file
		sout << "Optimization " << i << endl;
		sout << "Start Geometry" << endl;
		sout << "startXShift:         " << opt.xShift << endl;
		sout << "startCrossingAngle:  " << opt.crossingAngle << endl;
		sout << "startAxialRotation:  " << opt.axialRotation << endl;
		sout << "startZShift:         " << opt.zShift << endl << endl;
		sout << "Final Geometry" << endl;
		sout << "xShift:              " << xShift << endl;
		sout << "crossingAngle:       " << crossingAngle << endl;
		sout << "axialRotation:       " << axialRotation << endl;
		sout << "zShift:              " << zShift << endl << endl;
		sout << "Monomer Energy:      " << monomerEnergy << endl;
		sout << "Before Optimization" << endl;
		sout << "Dimer Energy:        " << opt.dimer << endl;
		sout << "Total Energy:        " << opt.total << endl << endl;
		sout << "Before Optimization" << endl;
		sout << "Dimer Energy:        " << dimerEnergy << endl;
		sout << "Total Energy:        " << finalEnergy << endl << endl;

		// Write an individual pdb for each sequence
		if (i ==0){
			bestRepackEnergy = finalEnergy;
			bestRepackState = MCOBest;
			sys.saveAltCoor("bestRepack");
			bestXShift = xShift;
			bestAxialRotation = axialRotation;
			bestZShift = zShift;
			bestCrossingAngle = crossingAngle;
		} else if (finalEnergy < bestRepackEnergy) {
			bestRepackEnergy = finalEnergy;
			bestRepackState = MCOBest;
			sys.saveAltCoor("bestRepack");
			bestXShift = xShift;
			bestAxialRotation = axialRotation;
			bestZShift = zShift;
			bestCrossingAngle = crossingAngle;
		}
		designWriter.write(sys.getAtomPointers(), true, false, true);

	}
	designWriter.close();
	sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

	sys.applySavedCoor("bestRepack");
	PDBWriter pdbWriter;
	pdbWriter.open(opt.outputDir + "/bestOptimizedBackbone.pdb");
	pdbWriter.write(sys.getAtomPointers(), true, false, true);
	pdbWriter.close();

	string repackFile = opt.outputDir + "/bestRepack.pdb";

	vector<string> energyList;
	energyList.push_back("CHARMM_VDW");
	energyList.push_back("SCWRL4_HBOND");
	energyList.push_back("CHARMM_IMM1REF");
	energyList.push_back("CHARMM_IMM1");

	//TODO: should have been an option in the design code, but for now I'll add these in as defaults and allow it to change in the future
	vector<string> AAList;
	AAList.push_back("A");
	AAList.push_back("L");
	AAList.push_back("I");
	AAList.push_back("F");
	AAList.push_back("S");
	AAList.push_back("T");
	AAList.push_back("W");
	AAList.push_back("Y");

	map<string, map<string, double>> sequenceEnergyMap;
	vector<string> mutants = getMutantInterfacialSequences(opt, AAList);

	//shuffle the mutant sequence list
	random_shuffle(mutants.begin(), mutants.end());
	
	AtomPointerVector &bestApv = sys.getAtomPointers();
	
	// Get sequence entropy map
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile);

	map<string, double> startEnergyMap;
	insertGeometriesIntoMap(opt, startEnergyMap, bestXShift, bestCrossingAngle, bestAxialRotation, bestZShift);
	getStartingStateEnergies(opt, spm, bestRepackState, startEnergyMap, sequenceEntropyMap);
	SasaCalculator sasa(sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();
	startEnergyMap["DimerSasa"] = dimerSasa;
	startEnergyMap["MonomerSasa"] = 0;
	sequenceEnergyMap[opt.sequence] = startEnergyMap;
	PDBWriter mutantWriter;
	mutantWriter.open(opt.outputDir + "/allMutants.pdb");

	//Make directory for individual mutant pdbs
	string mutantOutputDir = opt.outputDir + "/mutants";
	string cmd = "mkdir -p " + mutantOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	for (uint i=0; i<mutants.size(); i++){
		string mutantSeq = mutants[i];
		map<string, double> energyMap;
		double monomerEnergy = computeMonomerEnergy(opt, trans, energyMap, mutantSeq, RNG, sout, err);
		double dimerEnergy = computeDimerEnergy(sys, opt, energyMap, mutantSeq, RNG, mutantWriter, sequenceEntropyMap, mutantOutputDir, sout, err);
		double finalEnergy = dimerEnergy-monomerEnergy;
		cout << "Dimer-Monomer = " << dimerEnergy << " - " << monomerEnergy << " = " << finalEnergy << endl;
		insertGeometriesIntoMap(opt, energyMap, bestXShift, bestCrossingAngle, bestAxialRotation, bestZShift);
		addSequenceToMap(opt, sequenceEnergyMap, energyMap, finalEnergy, mutantSeq, sout);
		cout << "Total = " << sequenceEnergyMap[mutantSeq]["Total"] << endl;
		saveEnergyDifference(opt, sequenceEnergyMap, mutantSeq);
	}
	mutantWriter.close();

	mutants.push_back(opt.sequence);
	rotate(mutants.rbegin(), mutants.rbegin()+1, mutants.rend());
	//TODO: save best geometry to a file and then use that file to load that geometry in repacks for mutants
	outputDesignFiles(opt, opt.rotamerSamplingString, opt.rotamerSamplingVector, mutants, sequenceEnergyMap);

	//TODO: add in a way to mutate through all of the potential combinations at positions that are potentially higher in clashing
	//1. I think the easy way to do this is to make the above repack code into a function, mutate all individually and add them to a list, then randomly pick and choose others to mutate until I get a desired number
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	sout.close();
	err.close();
}

/******************************************
 *  ======= FUNCTIONS DEFINITIONS =======
 ******************************************/
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
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
			}
			else{
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
		}
		else{
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

void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	int frontExt = _opt.tmStart;
	int endExt = _opt.tmEnd;//TODO: I changed this from endResNum...it should still work but if energies are weird this is likely why
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

vector<pair<int,double>> calculateResidueBurial(System &_sys) {
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

//TODO: add in the new load rotamers term from seqDesign and output properly depending on if I want to load by burial or not
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt){
	//Repack side chains based on sasa scores
	for (uint i = 0; i < _opt.rotamerSamplingVector.size()/2; i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+_opt.sequence.length());
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _opt.sasaRepackLevel[_opt.rotamerSamplingVector[i]];
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
			string posRot = _opt.sasaRepackLevel[_opt.rotamerSamplingVector[i]];
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

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	
	/* Faster code 
	for (uint i=0; i < _apvA.size(); i++) {
			_apvB[i]->copyAllCoor(*_apvA[i]);
			vector<CartesianPoint*>& bCoors = _apvB[i]->getAllCoor();

			for(uint j = 0; j < bCoors.size(); j++) {
				bCoors[j]->setX(0 - bCoors[j]->getX());
				bCoors[j]->setY(0 - bCoors[j]->getY());
			}
					
		}
	*/

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

void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); 
	_spm.runGreedyOptimizer(_greedyCycles);
}

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

void getSasaDifference(vector<pair<string,vector<uint>>> &_sequenceStatePair, map<string, map<string,double>> &_sequenceEnergyMap){
	for (uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		
		double monomerSasa = _sequenceEnergyMap[sequence]["MonomerSasa"];
		double dimerSasa = _sequenceEnergyMap[sequence]["DimerSasa"];

		// Calculate the amount of surface area lost from monomer to dimer
		double interfaceSasa = monomerSasa-dimerSasa;
		_sequenceEnergyMap[sequence]["InterfaceSasa"] = interfaceSasa;
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
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//Required Options
	opt.allowed.push_back("designCrd");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("designPdb");
	opt.allowed.push_back("configfile");

	//Design Parameters
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("state");
	opt.allowed.push_back("rotamerSamplingString");
	opt.allowed.push_back("rotamerSamplingVector");

	//Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("sequenceEntropyFile");
	opt.allowed.push_back("selfEnergyFile");
	opt.allowed.push_back("pairEnergyFile");
	
	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("thread");
	
	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	opt.allowed.push_back("numRepacks");
	
	//input monomerEnergy
	opt.allowed.push_back("monomer");
	opt.allowed.push_back("monoVdw");
	opt.allowed.push_back("monoHbond");
	opt.allowed.push_back("monoIMM1");
	
	//input previous dimerEnergy
	opt.allowed.push_back("dimer");
	opt.allowed.push_back("dimerVdw");
	opt.allowed.push_back("dimerHbond");
	opt.allowed.push_back("dimerIMM1");
	opt.allowed.push_back("total");

	//Shift Size
	opt.allowed.push_back("deltaX");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaZ");

	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	//Load Rotamers from SASA values (from sgfc)
	opt.allowed.push_back("sasaRepack");
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("keepOriginalRotamer");
	
	//Other Options
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	opt.allowed.push_back("estimateMonomerWithBaseline");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");
	
	opt.allowed.push_back("energyTermsToOutput");
	opt.allowed.push_back("runNumber");
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

	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()) {
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}
	// Design Parameters
	opt.sequence = OP.getString("sequence");
	if(OP.fail()) {
		opt.errorMessages += "sequence not specified using L\n";
		opt.errorFlag = true;
	}
	opt.rotamerSamplingString = OP.getString("rotamerSamplingString");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify state, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.rotamerSamplingVector = OP.getIntVector("rotamerSampling");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify state, make sure they are space separated\n";
		opt.errorFlag = true;
	}

	// Geometry
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting to 9.2\n";
		opt.warningFlag = true;
		opt.xShift = 9.2;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to 25\n";
		opt.warningFlag = true;
		opt.crossingAngle = 25;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 40\n";
		opt.warningFlag = true;
		opt.axialRotation = 40;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to 2\n";
		opt.warningFlag = true;
		opt.zShift = 2;
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 25\n";
		opt.warningFlag = true;
		opt.thread = 25;
	}
	//TODO: maybe not in this code, but I feel like changing the thread of a sequence to see how it interacts at different threads could be helpful?

	//Load Rotamers using SASA values (from sgfc)
	opt.sasaRepackLevel = OP.getMultiString("sasaRepackLevel");
	if (OP.fail()) {
		opt.warningMessages += "sasaRepacklevel not specified! Default to one level at SL90.00";
		opt.sasaRepackLevel.push_back("SL90.00");
	}
	opt.keepOriginalRotamer = OP.getBool("keepOriginalRotamer");
	if (OP.fail()) {
		opt.warningMessages += "keepOriginalRotamer not specified, default true\n";
		opt.warningFlag = true;
	}

	//Monte Carlo variables
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC cycles not specified, default to 100\n";
		opt.warningFlag = true;
		opt.MCCycles = 100;
	}

	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC max rejects not specified, default to using 5\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 5;
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
		opt.warningMessages += "MCCurve not specified using SIGMOID(3)\n";
		opt.warningFlag = true;
		opt.MCCurve = 3;
	}
	
	opt.numRepacks = OP.getInt("numRepacks");
	if (OP.fail()) {
		opt.warningMessages += "Number of repack cycles not specified, default to 1\n";
		opt.warningFlag = true;
		opt.numRepacks = 1;
	}

	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaX = 0.5;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 5.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 5.0;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 4.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 4.0;
	}
	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.5;
	}

	// monomer energies
	opt.monomer = OP.getDouble("monomer");
	if (OP.fail()) {
		opt.errorMessages += "monomer not specified\n";
		opt.errorFlag = true;
	}
	opt.monoVdw = OP.getDouble("monoVdw");
	if (OP.fail()) {
		opt.errorMessages += "monoVdw not specified\n";
		opt.errorFlag = true;
	}
	opt.monoHbond = OP.getDouble("monoHbond");
	if (OP.fail()) {
		opt.errorMessages += "monoHbond not specified\n";
		opt.errorFlag = true;
	}
	opt.monoIMM1 = OP.getDouble("monoIMM1");
	if (OP.fail()) {
		opt.errorMessages += "monoIMM1 not specified\n";
		opt.errorFlag = true;
	}

	// dimer energies
	opt.dimer = OP.getDouble("dimer");
	if (OP.fail()) {
		opt.errorMessages += "dimer not specified\n";
		opt.errorFlag = true;
	}
	opt.dimerVdw = OP.getDouble("dimerVdw");
	if (OP.fail()) {
		opt.errorMessages += "dimerVdw not specified\n";
		opt.errorFlag = true;
	}
	opt.dimerHbond = OP.getDouble("dimerHbond");
	if (OP.fail()) {
		opt.errorMessages += "dimerHbond not specified\n";
		opt.errorFlag = true;
	}
	opt.dimerIMM1 = OP.getDouble("dimerIMM1");
	if (OP.fail()) {
		opt.errorMessages += "dimerIMM1 not specified\n";
		opt.errorFlag = true;
	}
	opt.total = OP.getDouble("Total");
	if (OP.fail()) {
		opt.errorMessages += "total not specified\n";
		opt.errorFlag = true;
	}

	// Weights
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
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine outputDir";
		opt.errorFlag = true;
	}
	opt.selfEnergyFile = OP.getString("selfEnergyFile");
	if (OP.fail()) { 
		opt.warningMessages += "selfEnergyFile not specified, default /exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_22_IMM1SelfBaseline.txt\n";
		opt.warningFlag = true;
		opt.selfEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_22_IMM1SelfBaseline.txt";
	}
	opt.pairEnergyFile = OP.getString("pairEnergyFile");
	if (OP.fail()) { 
		opt.warningMessages += "pairEnergyFile not specified, default /exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_22_IMM1PairBaseline.txt\n";
		opt.warningFlag = true;
		opt.pairEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_22_IMM1PairBaseline.txt";
	}
	opt.sequenceEntropyFile = OP.getString("sequenceEntropyFile");
	if (OP.fail()) { 
		opt.warningMessages += "sequenceEntropyFile not specified, default to /exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_21_seqEntropies.txt\n";
		opt.warningFlag = true;
		opt.sequenceEntropyFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_11_21_seqEntropies.txt";
	}
	opt.designPdb = OP.getString("designPdb");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine designPdb";
		opt.errorFlag = true;
	}
	opt.designCrd = OP.getString("designCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine designCrd";
		opt.errorFlag = true;
	}
	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine backboneCrd";
		opt.warningFlag = true;
		opt.backboneCrd = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/69-gly-residue-helix.crd";
	}

	//Other Options
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using true\n";
		opt.warningFlag = true;
		opt.verbose = true;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}
	opt.interfaceLevel = OP.getInt("interfaceLevel");
	if (OP.fail()) {
		opt.warningMessages += "interfaceLevel not specified using 2\n";
		opt.warningFlag = true;
		opt.interfaceLevel = 2;
	}
	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}

	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.sequence.length();
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
	}

	opt.estimateMonomerWithBaseline = OP.getBool("estimateMonomerWithBaseline");
	if (OP.fail()) {
		opt.warningMessages += "estimateMonomerWithBaseline not specified using true\n";
		opt.warningFlag = true;
		opt.estimateMonomerWithBaseline = true;
	}

	opt.energyTermsToOutput = OP.getStringVector("energyTermsToOutput");
	if (OP.fail()) {
		if (opt.estimateMonomerWithBaseline == true){
			opt.energyTermsToOutput.push_back("Total");
			opt.energyTermsToOutput.push_back("Dimer");
			opt.energyTermsToOutput.push_back("Monomer");
			opt.energyTermsToOutput.push_back("DimerDiff");
			opt.energyTermsToOutput.push_back("VDWDimer");
			opt.energyTermsToOutput.push_back("VDWDiff");
			opt.energyTermsToOutput.push_back("HBONDDimer");
			opt.energyTermsToOutput.push_back("HBONDDiff");
			opt.energyTermsToOutput.push_back("IMM1Dimer");
			opt.energyTermsToOutput.push_back("IMM1Diff");
			opt.energyTermsToOutput.push_back("SequenceProbability");
			opt.energyTermsToOutput.push_back("DimerSasa");
			opt.energyTermsToOutput.push_back("MonomerSasa");
			opt.energyTermsToOutput.push_back("startXShift");
			opt.energyTermsToOutput.push_back("xShift");
			opt.energyTermsToOutput.push_back("xShiftDiff");
			opt.energyTermsToOutput.push_back("startCrossingAngle");
			opt.energyTermsToOutput.push_back("crossingAngle");
			opt.energyTermsToOutput.push_back("crossingAngleDiff");
			opt.energyTermsToOutput.push_back("startAxialRotation");
			opt.energyTermsToOutput.push_back("axialRotation");
			opt.energyTermsToOutput.push_back("axialRotationDiff");
			opt.energyTermsToOutput.push_back("startZShift");
			opt.energyTermsToOutput.push_back("zShift");
			opt.energyTermsToOutput.push_back("zShiftDiff");
		} else {
			opt.energyTermsToOutput.push_back("Total");
			opt.energyTermsToOutput.push_back("Dimer");
			opt.energyTermsToOutput.push_back("Monomer");
			opt.energyTermsToOutput.push_back("DimerDiff");
			opt.energyTermsToOutput.push_back("VDWDimer");
			opt.energyTermsToOutput.push_back("VDWMonomer");
			opt.energyTermsToOutput.push_back("VDWDiff");
			opt.energyTermsToOutput.push_back("HBONDDimer");
			opt.energyTermsToOutput.push_back("HBONDMonomer");
			opt.energyTermsToOutput.push_back("HBONDDiff");
			opt.energyTermsToOutput.push_back("IMM1Dimer");
			opt.energyTermsToOutput.push_back("IMM1Monomer");
			opt.energyTermsToOutput.push_back("IMM1Diff");
			opt.energyTermsToOutput.push_back("SequenceProbability");
			opt.energyTermsToOutput.push_back("DimerSasa");
			opt.energyTermsToOutput.push_back("MonomerSasa");
		}
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
