#include <iostream>
#include <fstream>
#include <sstream>

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
#include "ResidueSelection.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "generateBaselineIMM1_noREF";
string programDescription = "This program is based on generateBaselinevdW and designs randomized sequences to calculate baselines for self and pair energies for future sequence design of vdW dependent sequences; this is a better organized version of the original generateSelfPairBaselinevdW, getting rid of any extraneous code; the CHARMM_IMM1REF energy is not used in this code (gives a large standard deviation to my self energy data)";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "22 September 2020";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

/******************************************************************************************************************************************************************************/

struct Options{
	string sequence;
	string backboneAA;
	int backboneLength;
	int seqNumber;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	// TODO: think I have to change this to a vector for specifically the interfacial residues 
	int startResNum;
	int endResNum;
	int sequenceStart;

	bool deleteTerminalHbonds;
	
	string SL; //number of rotamers

	// transformation
	double XShift;
	double ZShift;
	double crossAng;
	double axialRot;
	bool transform;
	int thread;
	int bbThread;
	
	// input files
	string helixGeoFile;
	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string monoRotLibFile;
	string infile;
	string rulesFile;

	// side-chain repack variable
	int mcCycles;
	int mcMaxRejects;
	double mcStartTemp;
	double mcEndTemp;
	int mcCurve;

	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;

	bool verbose;
	int greedyCycles;
	int seed;

	int numberOfStructuresToMCRepack;
	double energyCutOff;
	
	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	bool inputMonomerE;
	int monoE_vdw;
	int monoE_hbond;
	int monoE_solv;
	int monoE_solvRef;
	bool printTermEnergies;

	int start;
	int end;

	// alternate identities and weights
	vector<string> ids;
	vector<int> weights;

	//Pair parameters
	int pairDist;

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
	string runNumber;

	string configfile;
	string datafile;
};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

map<string, unsigned int> interfaceResidueCheck(AtomPointerVector & _chainA, AtomPointerVector & _chainB) {
	map<string, unsigned int> atomsWithIn4AOfPosition;
	for (uint i=0; i < _chainA.size(); i++) {
		if (_chainA[i]->getName() != "CA" && _chainA[i]->getName() != "C" && _chainA[i]->getName() != "N" && _chainA[i]->getName() != "HN" && _chainA[i]->getName() != "O") {
			for (uint j=0; j < _chainB.size(); j++) {
				if (_chainB[j]->getName() != "CA" && _chainB[j]->getName() != "C" && _chainB[j]->getName() != "N" && _chainB[j]->getName() != "HN" && _chainB[j]->getName() != "O") {
					if (_chainA[i]->distance(*_chainB[j]) < 4.0) {
						if (atomsWithIn4AOfPosition.find(_chainA[i]->getIdentityId()) != atomsWithIn4AOfPosition.end()) {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()]++;
						}
						else {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()] = 1;
						}
					}
				}
			}
		}
	}
	
	//for(map<string, unsigned int>::iterator it=atomsWithIn4AOfPosition.begin(); it != atomsWithIn4AOfPosition.end(); it++) {
	//	_fout << "pos: " << it->first << " count: " << it->second << endl;
	//}

	return atomsWithIn4AOfPosition;
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
	return "A" + ps;
}//changed this for only one chain


string convertToPolymerSequence(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}//changed this for only one chain

string generatePolyLeu(string _backboneAA, int _sequenceLength) {
	string polyLeu = "";
	for (uint i=0; i<_sequenceLength; i++){
		polyLeu = polyLeu + _backboneAA;
	}
	return polyLeu;
}

string generateMultiIDPolymerSequence(string _seq, vector<string> _alternateIds) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
		if (it-_seq.begin() == 0 || it-_seq.begin() == 35){
			ps = ps + " [";
			for (uint i=0; i<_alternateIds.size(); i++){
				if(_alternateIds[i] == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + _alternateIds[i];
				}
			}
			ps = ps + "] ";
			counter++;
		}
	}
	ps = ":{" + MslTools::intToString(0) + "} " + ps;
	//return "A" + ps + "\nB" + ps;
	return "A" + ps;
}

string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _varPos) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if (_varPos[counter] == 1){
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
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

string generatePolymerSequence(string _backboneAA, int _sequenceLength, int _startResNum) {
	string ps = "";
	string resName = MslTools::getThreeLetterCode(_backboneAA);
	if(resName == "HIS") {
		resName = "HSE";
	}
	for (uint i=0; i<_sequenceLength; i++){
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	//return "A" + ps + "\nB" + ps;
	return "A" + ps;
}

//void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable) {
void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {

	_spm.setOnTheFly(1);
	//_spm.recalculateNonSavedEnergies(_savedEnergyFlagTable);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
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

struct HbondInfo{
	string chain;
	string resA;
	string resB;
	string atomName;
	double distBreak;

	HbondInfo(string _c, string _resA, string _resB, string _aName, string _dist) {
		chain = _c;
		resA = _resA;
		resB = _resB;
		atomName = _aName;
		distBreak = MslTools::toDouble(_dist);
	}

	bool isValid(System & _sys) {
		string posA = chain + "," + resA;
		string posB = "";
		if(chain == "A") {
			posB = "B," + resB;
		} else {
			posB = "A," + resB;
		}
		if(!_sys.positionExists(posA) || !_sys.positionExists(posB) ) {
			return false;
		}

		if(atomName == "HA1") {
			if(!_sys.identityExists(posA + ",GLY" )) {
				return false;
			}
		}
		return true;
	}

	void print(ofstream & pout) {
		pout << chain << " " << resA << " " << resB << " " << atomName << " " << distBreak << endl;
	}
};

bool hydrogenBondCheck(System & _sys, vector<string> _parsedGeoInformation, double & _xShiftStart) {
	
	vector<HbondInfo*> hbondList;
	for(unsigned k = 5; k < _parsedGeoInformation.size(); k+=5 ) {
		HbondInfo * t = new HbondInfo(_parsedGeoInformation[k],_parsedGeoInformation[k+1],_parsedGeoInformation[k+2],_parsedGeoInformation[k+3],_parsedGeoInformation[k+4]);
		if(t->isValid(_sys)) {
			hbondList.push_back(t);
		} else {
			delete t;
		}
	}


	/*
	for(int i = 0; i < hbondList.size(); i++) {
		hbondList[i]->print(fout);
	}
	*/
	if(hbondList.size() < 4) {
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return false;
	} else {
		_xShiftStart = hbondList[hbondList.size() - 4]->distBreak - 0.1;
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return true;
	}

}//*/

void deleteTerminalHydrogenBondInteractions(System &_sys, Options& _opt) {
	// look at all hbond interactions and remove those
	// remove any interaction that has a donor or acceptor from residues 1 2 3 and n-2 n-1 n on each chain that are not part of the TM

	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize(); // number of positions in chain A = chain B
	// compute the extensions at the beginning and the end
	//int frontExt = _opt.tmStart - _opt.startResNum;
	//int endExt = _opt.endResNum - _opt.tmEnd;
	int frontExt = _opt.tmStart;
	int endExt = _opt.endResNum;
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain& thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(frontExt > i) {
				atoms += positions[i]->getAtomPointers(); 
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			if(endExt > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers(); 
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");
}

map<string,double> getEnergyByTerm(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

map<string,double> getEnergyByTermDoubled(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  2.0* _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

void readRulesFile(string _fileName, map<int, string> & _rulesFileMap) {
	ifstream file;
	file.open(_fileName.c_str());
	if(!file.is_open()) {
		cerr << "Unable to open " << _fileName << endl;
		exit(0);
	}

	string tmpRulesLine;

	while(file) {
		getline(file, tmpRulesLine);
		if (tmpRulesLine.length() > 1) {
			vector<string> token = MslTools::tokenizeAndTrim(tmpRulesLine,",");
			_rulesFileMap[MslTools::toInt(token[0])] = tmpRulesLine;
		}
	}
	file.close();
}

vector<string> getInterHelicalHbonds(EnergySet* & _ESet) {
	unsigned int numHbonds = 0;
	// Why are we doing this?
	//_ESet->setAllTermsInactive();
	//_ESet->setTermActive("SCWRL4_HBOND", true);
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];

	vector<string> hbonds;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			//_fout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
			hbonds.push_back(atoms[0]->getAtomOfIdentityId() + ":" + atoms[2]->getAtomOfIdentityId() + "=" + MslTools::doubleToString(atoms[0]->distance(*atoms[2])) + "," + MslTools::doubleToString(e));
			numHbonds++;
		}
	}
	// Why are we doing this?
	//_ESet->setAllTermsActive();
	//_ESet->setTermActive("CHARMM_ELEC", false);
	return hbonds;
}

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

vector<vector<vector<vector<bool> > > > createSavedEnergyFlagTable(System & _sys) {
	vector<vector<vector<vector<bool> > > > savedEnergyFlagTable; // true if the energy is computed already and available
	vector<unsigned int> variablePos = _sys.getVariablePositions();

	for (int i=0; i<variablePos.size(); i++) {
		savedEnergyFlagTable.push_back(vector<vector<vector<bool> > >());
		//cout << i;

		for (int j=0; j < _sys.getPosition(variablePos[i]).getTotalNumberOfRotamers(); j++) {
			savedEnergyFlagTable[i].push_back(vector<vector<bool> >());
			//cout << " " << i << "/" << j << "-";

			for (int k=0; k < i; k++) {
				if (i==k) {
					continue;
				}
				savedEnergyFlagTable[i][j].push_back(vector<bool>());

				for (int l=0; l < _sys.getPosition(variablePos[k]).getTotalNumberOfRotamers(); l++) {
			//		cout << " " << k << "/" << l;
					if (_sys.getPosition(variablePos[i]).getChainId() == _sys.getPosition(variablePos[k]).getChainId()) {
						savedEnergyFlagTable[i][j][k].push_back(true);
			//			cout << "=T";
					}
					else {
						savedEnergyFlagTable[i][j][k].push_back(false);
			//			cout << "=F";
					}
				}
			}
		}
		//cout << endl;
	}

	return savedEnergyFlagTable;
}

double computeMonomerEnergy(System & _sys, Transforms & _trans, Options& _opt, System & _helicalAxis, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);	
	
	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	//CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	SystemRotamerLoader monoRot(monoSys, _opt.monoRotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hBondFile);
	monohb.buildInteractions(30);

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

	deleteTerminalHydrogenBondInteractions(monoSys, _opt);
	
	/******************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	SelfPairManager monoSpm;
	monoSpm.seed(_opt.seed); 
	monoSpm.setVerbose(_opt.verbose);

	for (uint k=0; k < monoSys.positionSize(); k++) {
		Position &pos = monoSys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!monoRot.loadRotamers(&pos, pos.getResidueName(), "SL97.00")) { 
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}


	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	_fout << "Monomer - VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;

	monoSpm.setSystem(&monoSys);
	monoSpm.updateWeights();

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
	//moveZCenterOfCAMassToOrigin(chainA, _helicalAxis.getAtomPointers(), _trans);
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
        monoSpm.runGreedyOptimizer(_greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	_helicalAxis.saveAltCoor("BestMonomerAxis");
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
	_fout << "current Z: -5 Energy: " << currentEnergy << endl; 

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0); 
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
		_fout << "current Z: " << currentZ << " Energy: " << currentEnergy << endl;

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
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
			monoSpm.runGreedyOptimizer(_greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			//_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy << endl;
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

	MonteCarloManager MCMngr(1000.0, 0.5, _MCCycles, MonteCarloManager::EXPONENTIAL, _MCMaxRejects);
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
			monoSpm.runGreedyOptimizer(_greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_fout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_fout << "state rejected   energy: " << currentEnergy << endl;
		}
		else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			_fout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}

		counter++;
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	_fout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;


	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");
	_fout << monoEset->getSummary();
	_fout << endl;

	// print the monomer
	// TODO:change seq monomer maybe to seq# for each sequence
	string monoOutCrdFile  = _opt.pdbOutputDir + "/seq_monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}
	
	string monoOutPdbFile  = _opt.pdbOutputDir + "/seq_monomer.pdb";
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	// Store monomer energy by term
	if(_opt.printTermEnergies) {
		monoSys.calcEnergy();
		//_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix
		_monomerEnergyByTerm = getEnergyByTerm(monoSys.getEnergySet()); 

		ofstream meOut;
		string meOutName  = _opt.pdbOutputDir + "/seq_monomer.energy";
		meOut.open(meOutName.c_str());
		if(!meOut.is_open()) {
			cerr << "Unable to open " << meOutName << endl;
			exit(0);
		}
		for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
			meOut << it->first << " " << it->second << endl;
		}
		meOut.close();

	}

	//double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	double finalEnergy = monoSpm.getMinBound()[0]; 
	return finalEnergy;

}


Options parseOptions(int _argc, char * _argv[], Options defaults);

//skeleton code for adding given sequence into different portions of a polyLeu chain
string generateSequence(int _startResNum, int _sequenceLength, int _backboneLength, string _backboneAA, string &_sequence, bool _backbone=true){
	string seq = "";
	for(uint i=0; i<_sequenceLength; i++){
		//somehow multiply the amount of _backboneAA, then add in the _sequence, then add in leftover _backboneAA
		//might be able to switch this to a while loop like below?
		seq += _backboneAA;	
	}
	//seq += _sequence; use at some point if I want to do something like threading through the sequence
	if (_backbone == true){
		while (seq.length() < _backboneLength){
			seq += _backboneAA;
			//cout << seq << endl;
		}
	}
	return seq;
}

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Datafile: " << _op.datafile << endl << endl;

	_fout << "Options from Datafile" << endl;
	_fout << setiosflags(ios::fixed) << setprecision(3) << "XShift: " << _op.XShift << " ZShift: " << _op.ZShift << " Axial Rotation: " << _op.axialRot << " Crossing Angle: " << _op.crossAng << endl;
	
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "Warning Messages: " << _op.warningMessages << endl << endl;

	_fout << "Other Parameters" << endl;
	_fout << "backboneCrd " << _op.backboneCrd << endl;
	//_fout << "logFile " << _op.logFile << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	//_fout << "fullSequence " << _op.fullSequence << endl;
	_fout << "tmStart " << _op.tmStart << endl;
	_fout << "tmEnd " << _op.tmEnd << endl;

	_fout << "helixGeoFile " << _op.helixGeoFile << endl;
	_fout << "rulesFile " << _op.rulesFile << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "solvFile " << _op.solvFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;
	_fout << "monoRotLibFile " << _op.monoRotLibFile << endl;

	_fout << "MCCycles " << _op.mcCycles << endl;
	_fout << "MCMaxRejects " << _op.mcMaxRejects << endl;
	_fout << "MCStartTemp " << _op.mcStartTemp << endl;
	_fout << "MCEndTemp " << _op.mcEndTemp << endl;
	_fout << "MCCurve " << _op.mcCurve << endl;

	//_fout << "deltaZ " << _op.deltaZ << endl;
	//_fout << "deltaAx " << _op.deltaAx << endl;
	//_fout << "deltaCross " << _op.deltaCross << endl;
	//_fout << "deltaX " << _op.deltaX << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "numberOfStructuresToMCRepack " << _op.numberOfStructuresToMCRepack << endl;
	_fout << "energyCutOff " << _op.energyCutOff << endl;

	//_fout << "uniprotName " << _op.uniprotName << endl;
	//_fout << "uniprotAccession " << _op.uniprotAccession << endl;

	_fout << "monoE_vdw " << _op.monoE_vdw << endl;
	_fout << "monoE_solv " << _op.monoE_solv << endl;
	_fout << "monoE_solvRef" << _op.monoE_solvRef << endl;
	_fout << "monoE_hbond" << _op.monoE_hbond << endl;

	_fout << "deleteTerminalHbonds " << _op.deleteTerminalHbonds << endl;

	_fout << "fullSequenceStart " << _op.sequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;
	_fout << "weight_solv " << _op.weight_solv << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

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
		}
		else{
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}
	}	
}

void identifyInterface(System &_sys, vector<int> &_pos, vector<double> &_dists){
	int count = 0;
	for (uint k=0; k<_sys.positionSize()/2; k++) {
		//cout << _sys.positionSize() << endl;
		Position &pos1 = _sys.getPosition(k);
		Position &pos2 = _sys.getPosition(k+_sys.positionSize()/2);
		
		Atom &c1 = pos1.getAtom("CA");
		Atom &c2 = pos2.getAtom("CA");
		
		double dist;
		double most;
		double x = 0;
		dist = c1.distance(c2);
		//cout << "Dist " << k << ": " << dist << endl;
		if (_pos.size() < 12){
			_dists.push_back(dist);
			_pos.push_back(k);
		}
		else {
			count++;
			if (count <= 5){
				most = _dists[0];
				x = 0;
				//cout << "Current: " << endl;
				for (uint j=0; j<8; j++){
					//cout << _pos[j] << ": " << _dists[j] << endl;
					if (most < _dists[j]){
						most = _dists[j];
						x = j;
					}
				}
				if (dist <= most){
					//cout << "Before: " << x << ": " << _dists[x];
					_dists[x] = dist;
					_pos[x] = k;
					//cout << " ; After: " << _pos[x] << ": " << _dists[x] << endl;
					count = 0;
				}
			}
			else{
				cout << "Interface Identified!" << endl;
				k = _sys.positionSize()/2;
			}
		}
	}
}

vector<int> interface01(System &_sys, vector<int> &_pos){
	vector<int> varPos;
	for (uint k=0; k<_sys.positionSize(); k++){
		varPos.push_back(0);
	}
	for (uint j=0; j<_pos.size(); j++){
		varPos[_pos[j]] = 1;
		varPos[_pos[j]+_sys.positionSize()/2] = 1;
	}
	return varPos;
}

void addIdentities(CharmmSystemBuilder &_CSB, System &_sys, vector<string> &_ids, vector<int> &_varPos, bool multipleIds=true){
	if (_ids.size() == 1 || multipleIds != true){
		for (uint k=0; k<_varPos.size(); k++){
			if (_varPos[k] == 1){
				Position &pos = _sys.getPosition(k);
				_CSB.addIdentity(pos, _ids[0]);
			}
			else{
				continue;
			}
		}
		cout << "Id added!" << endl;
	}
	else{
		for (uint k=0; k<_varPos.size(); k++){
			if (_varPos[k] == 1){
				Position &pos = _sys.getPosition(k);
				_CSB.addIdentity(pos, _ids);
			}
			else{
				continue;
			}
		}
		cout << "Ids added!" << endl;
	}
}

vector<string> weightAAs(vector<string> &_AAs, vector<int> &_weights, int _seqLength, RandomNumberGenerator &_RNG){
	vector<string> weightedAAs;
	for (uint j=0; j<_AAs.size(); j++){
		for (uint k=0; k<_weights[j]; k++){
			weightedAAs.push_back(_AAs[j]);
		}
	}
	random_shuffle(weightedAAs.begin(), weightedAAs.end(), _RNG);
	return weightedAAs;
}

string randomAASeqLeuEnds(vector<string> _weightedAAs, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	for (uint i=0; i<_seqLength; i++){
		if (i < 4 || i > _seqLength-5){
			seq += "L";
		}
		else{
			int a=0;
			a = _RNG.getRandomInt(_weightedAAs.size()-1);
			//cout << i << ": " << a << endl;
			if (a == -1){
				while (a == -1){
					a = _RNG.getRandomInt(_weightedAAs.size()-1);
				}
			}//I don't think I need this as long as I change a to an int
			//cout << "AA: " << weightedAAs[a] << endl;
			//randAA.push_back(a);
			seq += _weightedAAs[a];
		}
	}
	return seq;
}

string randomAASequence(vector<string> _weightedAAs, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	for (uint i=0; i<_seqLength; i++){
		int a=0;
		a = _RNG.getRandomInt(_weightedAAs.size()-1);
		//cout << i << ": " << a << endl;
		if (a == -1){
			while (a == -1){
				a = _RNG.getRandomInt(_weightedAAs.size()-1);
			}
		}//I don't think I need this as long as I change a to an int
		//cout << "AA: " << weightedAAs[a] << endl;
		//randAA.push_back(a);
		seq += _weightedAAs[a];
	}
	return seq;
}

vector<double> calcBaselineIMM1Energies(System &_sys, int _seqLength, double &_totalProt){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	vector<Position*>& positions = _sys.getPositions();

	for (vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
		cout << p-positions.begin() << endl;
		if (p-positions.begin() == 0 || p-positions.begin() == 35){
			for (uint i=0; i < (*p)->identitySize(); i++){
				(*p)->setActiveIdentity(i);
				string residue = "resi, chain A and resi ";
				string number = to_string(p-positions.begin()+1);
				sel.select(residue += number);
				double resi = _sys.calcEnergy("resi");
				ener.push_back(resi);
				_totalProt += resi;
				sel.clearStoredSelections();
			}
		}
	}
	return ener;
}

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, double &_totalProt, double &_interactions, double &_refInteractions){
	vector <double> ener;
	EnergySet * eset = _sys.getEnergySet();
	double imm1refener = 0;
	double imm1ener = 0;
	double interref = 0;
	double inter = 0;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=0; i<_seqLength; i++){
		sel.select("allProt, all");
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		cout << residue << endl;
		//double prot = _sys.calcEnergy("resi", "allProt");
		//cout << eset->getSummary() << endl;
		//imm1refener = imm1refener + eset->getTermEnergy("CHARMM_IMM1REF");
		//imm1ener = imm1ener + eset->getTermEnergy("CHARMM_IMM1")/2;
		//interref = interref + eset->getTermNumberOfInteractionsCalculated("CHARMM_IMM1REF");
		//inter = inter + eset->getTermNumberOfInteractionsCalculated("CHARMM_IMM1")/2;
		double resi = _sys.calcEnergy("resi");
		cout << eset->getSummary() << endl;
		cout << eset->getTermEnergy("CHARMM_IMM1") << endl;
		imm1ener = imm1ener + eset->getTermEnergy("CHARMM_IMM1");
		inter = inter + eset->getTermNumberOfInteractionsCalculated("CHARMM_IMM1");
		cout << resi << endl;
		//ener.push_back(eset->getTermEnergy("CHARMM_IMM1REF") + eset->getTermEnergy("CHARMM_IMM1")/2);
		ener.push_back(resi);
		_totalProt += resi;
		sel.clearStoredSelections();
	}
	cout << "REF: " << interref << "; " << imm1refener << endl;
	cout << "IMM1: " << inter << "; " << imm1ener << endl;
	_interactions = inter;
	_refInteractions = interref;

	return ener;
}

vector<double> calcPairBaselineEnergies(System &_sys, Options _opt, int _seqLength, double &_totalProt, double &_interactions, double &_refInteractions){
	vector <double> ener;
	EnergySet * eset = _sys.getEnergySet();
	AtomSelection sel(_sys.getAtomPointers());
	double imm1refener = 0;
	double imm1ener = 0;
	double interref = 0;
	double inter = 0;
	//string r1CA   = "1rCA, resi 3 and name CA";
	//string r1CB   = "1rCB, resi 3 and name CB";
	//string r1CG   = "1rCG, resi 3 and name CG";
	//string r1CD1  = "1rCD1, resi 3 and name CD1";
	//string r1CD2  = "1rCD2, resi 3 and name CD2";
	//string r1N    = "1rN, resi 3 and name N";
	//string r1HN   = "1rHN, resi 3 and name HN";
	//string r1C    = "1rC, resi 3 and name C";
	//string r1O    = "1rO, resi 3 and name O";
	//string r1HA   = "1rHA, resi 3 and name HA";
	//string r1HG   = "1rHG, resi 3 and name HG";
	//string r1HB1  = "1rHB1, resi 3 and name HB1";
	//string r1HB2  = "1rHB2, resi 3 and name HB2";
	//string r1HD11 = "1rHD11, resi 3 and name HD11";
	//string r1HD12 = "1rHD12, resi 3 and name HD12";
	//string r1HD13 = "1rHD13, resi 3 and name HD13";
	//string r1HD21 = "1rHD21, resi 3 and name HD21";
	//string r1HD22 = "1rHD22, resi 3 and name HD22";
	//string r1HD23 = "1rHD23, resi 3 and name HD23";
	//
	//string r2CA   = "2rCA, resi 2 and name CA";
	//string r2CB   = "2rCB, resi 2 and name CB";
	//string r2CG   = "2rCG, resi 2 and name CG";
	//string r2CD1  = "2rCD1, resi 2 and name CD1";
	//string r2CD2  = "2rCD2, resi 2 and name CD2";
	//string r2N    = "2rN, resi 2 and name N";
	//string r2HN   = "2rHN, resi 2 and name HN";
	//string r2C    = "2rC, resi 2 and name C";
	//string r2O    = "2rO, resi 2 and name O";
	//string r2HA   = "2rHA, resi 2 and name HA";
	//string r2HG   = "2rHG, resi 2 and name HG";
	//string r2HB1  = "2rHB1, resi 2 and name HB1";
	//string r2HB2  = "2rHB2, resi 2 and name HB2";
	//string r2HD11 = "2rHD11, resi 2 and name HD11";
	//string r2HD12 = "2rHD12, resi 2 and name HD12";
	//string r2HD13 = "2rHD13, resi 2 and name HD13";
	//string r2HD21 = "2rHD21, resi 2 and name HD21";
	//string r2HD22 = "2rHD22, resi 2 and name HD22";
	//string r2HD23 = "2rHD23, resi 2 and name HD23";

	//sel.select(r1CA);
	//sel.select(r1CB);
	//sel.select(r1CG);
	//sel.select(r1CD1);
	//sel.select(r1CD2);
	//sel.select(r1N);
	//sel.select(r1HN);
	//sel.select(r1C);
	//sel.select(r1O);
	//sel.select(r1HA);
	//sel.select(r1HG);
	//sel.select(r1HB1);
	//sel.select(r1HB2);
	//sel.select(r1HD11);
	//sel.select(r1HD12);
	//sel.select(r1HD13);
	//sel.select(r1HD21);
	//sel.select(r1HD22);
	//sel.select(r1HD23);
	//
	//sel.select(r2CA);
	//sel.select(r2CB);
	//sel.select(r2CG);
	//sel.select(r2CD1);
	//sel.select(r2CD2);
	//sel.select(r2N);
	//sel.select(r2HN);
	//sel.select(r2C);
	//sel.select(r2O);
	//sel.select(r2HA);
	//sel.select(r2HG);
	//sel.select(r2HB1);
	//sel.select(r2HB2);
	//sel.select(r2HD11);
	//sel.select(r2HD12);
	//sel.select(r2HD13);
	//sel.select(r2HD21);
	//sel.select(r2HD22);
	//sel.select(r2HD23);
	//vector<string> atomNames1;
	//vector<string> atomNames2;
	//atomNames1.push_back("1rCA");
	//atomNames1.push_back("1rCB");
	//atomNames1.push_back("1rCG");
	//atomNames1.push_back("1rCD1");
	//atomNames1.push_back("1rCD2");
	//atomNames1.push_back("1rN");
	//atomNames1.push_back("1rHN");
	//atomNames1.push_back("1rC");
	//atomNames1.push_back("1rO");
	//atomNames1.push_back("1rHA");
	//atomNames1.push_back("1rHG");
	//atomNames1.push_back("1rHB1");
	//atomNames1.push_back("1rHB2");
	//atomNames1.push_back("1rHD11");
	//atomNames1.push_back("1rHD12");
	//atomNames1.push_back("1rHD13");
	//atomNames1.push_back("1rHD21");
	//atomNames1.push_back("1rHD22");
	//atomNames1.push_back("1rHD23");
	//
	//atomNames2.push_back("2rCA");
	//atomNames2.push_back("2rCB");
	//atomNames2.push_back("2rCG");
	//atomNames2.push_back("2rCD1");
	//atomNames2.push_back("2rCD2");
	//atomNames2.push_back("2rN");
	//atomNames2.push_back("2rHN");
	//atomNames2.push_back("2rC");
	//atomNames2.push_back("2rO");
	//atomNames2.push_back("2rHA");
	//atomNames2.push_back("2rHG");
	//atomNames2.push_back("2rHB1");
	//atomNames2.push_back("2rHB2");
	//atomNames2.push_back("2rHD11");
	//atomNames2.push_back("2rHD12");
	//atomNames2.push_back("2rHD13");
	//atomNames2.push_back("2rHD21");
	//atomNames2.push_back("2rHD22");
	//atomNames2.push_back("2rHD23");
	//double bbener = 0;
	//vector<double> e;
	//vector<string> p;
	//for (uint x=0; x<atomNames1.size(); x++){
	//	for (uint y=0; y<atomNames2.size(); y++){
	//		string a1 = atomNames1[x];
	//		string a2 = atomNames2[y];
	//		bbener = bbener + _sys.calcEnergy(a1, a2);
	//		string s = a1 + ", " + a2;
	//		e.push_back(_sys.calcEnergy(a1, a2));
	//		p.push_back(s);
	//		cout << a1 << ", " << a2 << ": " << bbener << endl;
	//	}
	//}
	//for (uint x=0; x<e.size(); x++){
	//	if (e[x] > 0){
	//		cout << p[x] << ": " << e[x] << endl;
	//	}
	//}TODO: this is gross but was for troubleshooting; probably could have just looped through a selection of a residue for it's atoks? but came to me too late :(
	_sys.calcEnergy();
	cout << eset->getSummary() << endl;

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+1);
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength; j++){
			int dist = j-i;
			if (dist <= _opt.pairDist){
			//if (dist >= 0){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+1);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				imm1ener = imm1ener + eset->getTermEnergy("CHARMM_IMM1");
				inter = inter + eset->getTermNumberOfInteractionsCalculated("CHARMM_IMM1");
				//cout << eset->getSummary() << endl;
				cout << residue << "; " << resi1 << endl;
				//ener.push_back(eset->getTermEnergy("CHARMM_IMM1REF") + eset->getTermEnergy("CHARMM_IMM1")/2);
				ener.push_back(pair);
				_totalProt += pair;
				//cout << residue << "; " << resi1 << " : " << pair << endl;
			}
			else{
				j = _seqLength;
			}
		}
	}
	cout << "REF: " << interref << "; " << imm1refener << endl;
	cout << "IMM1: " << inter << "; " << imm1ener << endl;
	_interactions = _interactions + inter;
	_refInteractions = _refInteractions + interref;
	sel.clearStoredSelections();
	return ener;
}

double sumEnergyVector(vector<double> _energies){
	double ener = 0;
	for (uint i=0; i<_energies.size(); i++){
		ener = ener + _energies[i];
	}
	return ener;
}

void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, double _totalEnergy, double _hbond, double _vdw, int _seqNumber, bool _selfEner, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	if (_selfEner == true){
		_out << "AA:      Position:      Energy" << endl;//could be interesting to add rotamer number to this
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		_out << "Pair:    Distance:      Energy" << endl;
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				stringstream ss;
				ss <<_seq[i] << _seq[j];
				string p = ss.str();
				int d = j-i;
				pairs.push_back(p);
				dists.push_back(d);
			}
		}
		string spacer = ":         ";
		for (uint i=0; i<pairs.size(); i++){
			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << endl;
		}
	}
	_out << "Total Energy: " << _totalEnergy << endl;
	_out << "H-bond Energy: " << _hbond << endl;
	_out << "VDW Energy: " << _vdw << endl << endl;
}

void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, double _totalEnergy, int _seqNumber, bool _selfEner, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	if (_selfEner == true){
		_out << "AA:      Position:      Energy:" << endl;//could be interesting to add rotamer number to this
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		_out << "Pair:    Distance:      Energy:" << endl;
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				stringstream ss;
				ss <<_seq[i] << _seq[j];
				string p = ss.str();
				int d = j-i;
				pairs.push_back(p);
				dists.push_back(d);
			}
		}
		string spacer = ":         ";
		for (uint i=0; i<pairs.size(); i++){
			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << spacer << endl;
		}
	}
	_out << "Total Energy: " << _totalEnergy << endl;
}

//void printEnerFile(string _seq, vector<double> &_ener, int _seqNumber, bool _selfEner, ofstream &_out){
//	if (_selfEner == true){
//		if (_seqNumber == 0){
//			_out << "AA:      Position:      Energy" << endl;//could be interesting to add rotamer number to this
//		}
//		for (uint i=0; i<_seq.length(); i++){
//			_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
//		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
//	}
//	else{
//		if (_seqNumber == 0){
//			_out << "Pair:      Distance:      Energy" << endl;//could be interesting to add rotamer number to this
//		}
//		vector<string> pairs;
//		vector<double> dists;
//		for (uint i=0; i<_seq.length(); i++){
//			for (uint j=i+1; j<_seq.length(); j++){
//				stringstream ss;
//				ss <<_seq[i] << _seq[j];
//				string p = ss.str();
//				int d = j-i;
//				pairs.push_back(p);
//				dists.push_back(d);
//			}
//		}
//		string spacer = ":         ";
//		for (uint i=0; i<pairs.size(); i++){
//			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << endl;
//		}
//	}
//}

void printEnerFile(string _seq, Options _opt, vector<double> &_ener, int _seqNumber, bool _selfEner, ofstream &_out){
	if (_selfEner == true){
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "AA:      Position:      Energy:" << endl;//could be interesting to add rotamer number to this
		}
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "Pair: Distance: Position1: Position2: Energy:" << endl;//could be interesting to add rotamer number to this
		}
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				int d = j-i;
				if (d <= _opt.pairDist){
					stringstream ss;
					ss <<_seq[i] << _seq[j];
					string p = ss.str();
					pairs.push_back(p);
					dists.push_back(d);
				}
			}
		}
		string sp = ": ";
		int pos1 = 0;

		for (uint i=0; i<pairs.size(); i++){
			if (dists[i] == 1){
				pos1++;
			}
			int pos2 = pos1 + dists[i];
			_out << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << endl;
			//cout << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << endl;
		}
	}
}
//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

vector<vector<string>> positionToString(System &_sys, vector<int> &_varPos){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_varPos.size()/2; k++){
		//cout << "string" << endl;
		if (_varPos[k] == 1){
			vector<string> tempPos;

			Position &posA = _sys.getPosition(k);
			Position &posB = _sys.getPosition(k+_varPos.size()/2);

			string A = posA.toString();
			string B = posB.toString();

			string delimiter = " ";
			
			size_t p = 0;
			p = A.find(delimiter);

			tempPos.push_back(A.substr(0, p));
			tempPos.push_back(B.substr(0, p));

			stringPositions.push_back(tempPos);
		}
		else{
			continue;
		}
	}
	return stringPositions;
}

vector<double> getSubtractBetweenVectors(vector<double> &_v1, vector<double> &_v2){
	vector<double> v3;
	if (_v1.size() != _v2.size()){
		cout << "Vectors of different sizes; terminate program." << endl;
	}
	else{
		for (uint i=0; i<_v1.size(); i++){
			v3.push_back(_v1[i]-_v2[i]);
		}
	}
	return v3;
}

vector<double> getAdditionBetweenVectors(vector<double> &_v1, vector<double> &_v2){
	vector<double> v3;
	if (_v1.size() != _v2.size()){
		cout << "Vectors of different sizes; terminate program." << endl;
	}
	else{
		for (uint i=0; i<_v1.size(); i++){
			v3.push_back(_v1[i]+_v2[i]);
		}
	}
	return v3;
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
	
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%Y_%m_%d",timeinfo);
	string date(buffer);
	
	cout << date << endl;

	time (&startTime);
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		usage();
		exit(1);
	}

	/**********************************************************************************
	*
	*    printOutFiles
	*
	**********************************************************************************/
	ofstream sout;
	ofstream pout;
	ofstream ssout;
	ofstream spout;
	
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string soutName  = dir + "/IMM1_Self_Energies_" + opt.runNumber + ".out";
	string poutName  = dir + "/IMM1_Pair_Energies_" + opt.runNumber + ".out";
	string ssoutName = dir + "/Sequences_Self_" + opt.runNumber + ".out";
	string spoutName = dir + "/Sequences_Pair_" + opt.runNumber + ".out";
	
	sout.open(soutName.c_str());
	pout.open(poutName.c_str());
	ssout.open(ssoutName.c_str());
	spout.open(spoutName.c_str());
	
	/******************************************************************************
	 *                    === INITIALIZE POLYGLY AS BACKBONE ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd); 
	if(!cRead.read()) {
		cout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	PDBWriter writer;
	writer.open(dir + "/Run" + opt.runNumber + ".pdb");
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	int seqDiscard = 0;
	int seqAccept = 0;
	uint a = 0;
	int seqNumber = opt.seqNumber;
	//vector<string> completeSequences;//TODO: if I want to have a check to be sure that no sequences are the same (currently have a check outside of this code)
	
	//Random Number Generator
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);
	
	/******************************************************************************
	 *                     === GIVE SEQUENCES WEIGHTS ===
	 ******************************************************************************/
	vector<string> weightedAAs = weightAAs(opt.ids, opt.weights, opt.backboneLength, RNG);
	
	if (opt.backboneLength > weightedAAs.size()){
		cout << "ERROR: Length of protein is too long for total number of AAs given." << endl;
		exit(1);
	}

	while (seqAccept < seqNumber){
		/******************************************************************************
		 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
		 ******************************************************************************/
		cout << "Number Accepted: " << seqAccept << endl;
		cout << "Number Discarded: " << seqDiscard << endl;
		a = seqAccept;
		
		/******************************************************************************
		 *                         === RANDOMIZE SEQUENCES ===
		 ******************************************************************************/
		string seq = randomAASeqLeuEnds(weightedAAs, opt.backboneLength, RNG);
		string polySeq = randomAASeqLeuEnds(weightedAAs, opt.backboneLength, RNG);
		//string seq = generatePolyLeu("L", opt.backboneLength);
		//cout << seq << endl;
		//string polySeq = generateMultiIDPolymerSequence(seq, opt.ids);//To have neutral patches at the ends of my sequences to prevent any hydrogen bonding that could affect energy score
		polySeq = convertToPolymerSequenceNeutralPatch(polySeq, 1);
		//string polySeq = generatePolymerSequence("L", opt.backboneLength, opt.thread);
		PolymerSequence PS(polySeq);
		cout << PS << endl;
	
		/******************************************************************************
		 *                           === DECLARE SYSTEM ===
		 ******************************************************************************/
		System sys;
		CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
		CSB.setBuildTerm("CHARMM_ELEC", false);
		CSB.setBuildTerm("CHARMM_ANGL", false);
		CSB.setBuildTerm("CHARMM_BOND", false);
		CSB.setBuildTerm("CHARMM_DIHE", false);
		CSB.setBuildTerm("CHARMM_IMPR", false);
		CSB.setBuildTerm("CHARMM_U-BR", false);
		CSB.setBuildTerm("CHARMM_IMM1", true);
		
		CSB.setSolvent("MEMBRANE");
		CSB.setIMM1Params(15, 10);
	
		CSB.setBuildNonBondedInteractions(false);
		if(!CSB.buildSystem(PS)) {
			cerr << "Unable to build system from " << polySeq << endl;
			exit(0);
		} else {
			//fout << "CharmmSystem built for sequence" << endl;
		}
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();
		
		CSB.updateNonBonded(10,12,50); //eventually change these to options for cuton, cutoff, and cutnb
		sys.buildAllAtoms();
	
		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
		sysRot.defineRotamerSamplingLevels();
		
		// Add hydrogen bond term
		//HydrogenBondBuilder hb(sys, opt.hBondFile);
		//hb.buildInteractions(50);//set this to same as cutnb
		
		/******************************************************************************
		 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
		 ******************************************************************************/
		for (uint i=0; i<sys.atomSize(); i++){
			Atom at = sys.getAtom(i);
			if (!at.hasCoor()){
				cout << "Atom " << i << " was not assigned coordinates; program termination";
				break;
			}
			else{
				continue;
			}
		}
		cout << "All atoms have coordinates" << endl;
		writer.write(sys.getAtomPointers(), true, false, true);
	
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
		Eset->setTermActive("CHARMM_VDW", false);
		Eset->setTermActive("SCWRL4_HBOND", false);
		Eset->setTermActive("CHARMM_IMM1REF", false);
		
		// Set weights
		//Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
		//Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
		Eset->setWeight("CHARMM_IMM1", 1);
	
		/******************************************************************************
		 *            === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
		 ******************************************************************************/
		deleteTerminalHydrogenBondInteractions(sys,opt);
	
		/******************************************************************************
		 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
		 ******************************************************************************/
		loadRotamers(sys, sysRot, opt.SL);
		
		SelfPairManager spm;
		spm.seed(RNG.getSeed());
		spm.setSystem(&sys);
		spm.setVerbose(true);
		//spm.getMinStates()[0];
		spm.updateWeights();
		spm.setOnTheFly(true);
		//spm.calculateEnergies();
		
		//repackSideChains(spm, 10);
		
		/******************************************************************************
		 *                  === SET SYSTEM TO BEST SPM ROTAMERS ===
		 ******************************************************************************/
		//sys.setActiveRotamers(spm.getMinStates()[0]);
		sys.calcEnergy();
		double totalEnergy = sys.calcEnergy();
		cout << Eset->getSummary() << endl;
		double IMM1 = Eset->getTermEnergy("CHARMM_IMM1") ;
		
		/******************************************************************************
		 *               === CALCULATE ENERGIES FOR EACH POSITION ===
		 ******************************************************************************/
		double totalProt = 0;
		double totalSelfPair = 0;
		double totalPair = 0;
		double interactions = 0;
		double refInteractions = 0;
		cout <<"Sequence Length: " << seq.length() << endl;
		vector<double> AAener = calcBaselineEnergies(sys, seq.length(), totalSelfPair, interactions, refInteractions);
		vector<double> pairAAener = calcPairBaselineEnergies(sys, opt, seq.length(), totalSelfPair, interactions, refInteractions);
		
		/******************************************************************************
		 *                 === OUTPUT ENERGIES FOR EACH TERM ===
		 ******************************************************************************/
		double baselineEner     = sumEnergyVector(AAener);
		double baselinePairEner = sumEnergyVector(pairAAener);
		double totalIMM1Energy  = baselineEner + baselinePairEner;
	
		double numInteractions = interactions + refInteractions;
		cout << "Total Energy: " << totalEnergy << endl;
		cout << "BaselineIMM1 Energy: " << baselineEner << endl;
		cout << "BaselinePairIMM1 Energy: " << baselinePairEner << endl;
		cout << "Baseline Addition: " << totalIMM1Energy << endl;
		cout << "BaselineIMM1 Energy: " << totalSelfPair << endl;
		cout << "BaselinePairIMM1 Energy: " << totalPair << endl;
		cout << "Ref Interactions: " << numInteractions << endl;
		cout << "IMM1 Interactions: " << numInteractions << endl;
		cout << "Total Interactions: " << numInteractions << endl;
		sys.calcEnergy();
		cout << Eset->getSummary() << endl;
		
		/******************************************************************************
		 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
		 ******************************************************************************/
		printSeqFile(PS, seq, AAener, totalEnergy, seqAccept, true, ssout);
		printEnerFile(seq, opt, AAener, seqAccept, true, sout);
		//printSeqFile(PS, seq, pairAAener, totalEnergy, seqAccept, false, spout);//Commented out because repeating runs from 12_03_2020 which I already have data for pairs (only need the IMM1 self energies)
		//printEnerFile(seq, opt, pairAAener, seqAccept, false, pout);
		seqAccept++;
		sys.reset();
	}

	time (&endTime);
	diffTime = difftime(endTime, startTime);
	ssout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	spout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	writer.close();
	sout.close();
	pout.close();
	ssout.close();
	spout.close();
}

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
	 *  /exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/helixGenerator.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//opt.required.push_back("");
	//opt.allowed.push_back("");

	//opt.allowed.push_back("");
	
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	opt.allowed.push_back("seqNumber");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("threadBool");
	
	//transformation
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("transform");

	//side-chain repack variable
	opt.allowed.push_back("mcCycles");
	opt.allowed.push_back("mcMaxRejects");
	opt.allowed.push_back("mcStartTemp");
	opt.allowed.push_back("mcEndTemp");
	opt.allowed.push_back("mcCurve");

	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaX");

	opt.allowed.push_back("SL");
	
	// energy weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("start");
	opt.allowed.push_back("end");
	
	opt.allowed.push_back("ener");
	
	opt.allowed.push_back("ivalues");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	opt.allowed.push_back("numberOfStructuresToMCRepack");
	opt.allowed.push_back("energyCutOff");
	
	opt.allowed.push_back("uniprotName");
	opt.allowed.push_back("uniprotAccession");
	
	// input monomerEnergy
	opt.allowed.push_back("inputMonomerE");
	opt.allowed.push_back("monoE_vdw");
	opt.allowed.push_back("monoE_hbond");
	opt.allowed.push_back("monoE_solv");
	opt.allowed.push_back("monoE_solvRef");
	opt.allowed.push_back("printTermEnergies");

	// input files
	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("solvFile");
	
	// alternate ids and weights
	opt.allowed.push_back("ids");
	opt.allowed.push_back("weights");

	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("pairDist");

	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		cerr << "No options given!" << endl;
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()){
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}

	opt.pairDist = OP.getInt("pairDist");
	if (OP.fail()){
		opt.warningMessages += "pairDist not specified, using 15\n";
		opt.warningFlag = true;
		opt.pairDist = 15;
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
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}

	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages += "sequence (1 letter aa) not specified\n";
		opt.errorFlag = true;
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

	opt.rulesFile = OP.getString("rulesFile");
	if (OP.fail()) {
		opt.rulesFile = "/data01/sabs/tmRepacks/GLY_69_Homo_2/tmRules/rules_10kcals_vdw_only/tmRules.out";
		opt.warningMessages += "rulesFile not specified using " + opt.rulesFile + "\n";
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

	opt.XShift = OP.getDouble("XShift");
	if (OP.fail()) {
		opt.warningMessages += "XShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.XShift = 6.7;
	}
	opt.ZShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.ZShift = 0;
	}
	opt.axialRot = OP.getDouble("axialRot");
	if (OP.fail()) {
		opt.warningMessages += "axialRot not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.axialRot = 0;
	}
	opt.crossAng = OP.getDouble("crossAng");
	if (OP.fail()) {
		opt.warningMessages += "crossAng not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.crossAng = -40;
	}
	opt.transform = OP.getBool("transform");
	if (OP.fail()) {
		opt.warningMessages += "transform not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.transform = false;
	}
		
	opt.mcCycles = OP.getInt("mcCycles");
	if (OP.fail()) {
		opt.errorMessages += "Number of MC cycles not specified!\n";
		opt.errorFlag = true;
	}

	opt.mcMaxRejects = OP.getInt("mcMaxRejects");
	if (OP.fail()) {
		opt.mcMaxRejects = 10;
		opt.warningMessages += "Number of MC max rejects not specified, default to using 10\n";
		opt.warningFlag = true;
	}

	opt.mcStartTemp = OP.getDouble("MCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCStartTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.mcStartTemp = 1000.0;
	}
	opt.mcEndTemp = OP.getDouble("MCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.mcEndTemp = 0.5;
	}
	opt.mcCurve = OP.getInt("MCCurve");
	if (OP.fail()) {
		opt.warningMessages += "MCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.mcCurve = 2;
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
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.numberOfStructuresToMCRepack = OP.getInt("numberOfStructuresToMCRepack");
	if (OP.fail()) {
		opt.warningMessages += "numberOfStructuresToMCRepack not specified using 20\n";
		opt.warningFlag = true;
		opt.numberOfStructuresToMCRepack = 20;
	}
	opt.energyCutOff = OP.getDouble("energyCutOff");
	if (OP.fail()) {
		opt.warningMessages += "energyCutOff not specified using 100.0\n";
		opt.warningFlag = true;
		opt.energyCutOff = 100.0;
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

	opt.SL = OP.getString("rotLevel");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	}

	opt.backboneAA = OP.getString("backboneAA");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneAA not specified, default to glycine\n";
		opt.backboneAA = "G";
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 35\n";
		opt.backboneLength = 35;
	}
	opt.seqNumber = OP.getInt("seqNumber");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "seqNumber not specified, default to 10\n";
		opt.seqNumber = 10;
	}

	opt.start = OP.getInt("start");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "Start point not specified, default to 0\n";
		opt.start = 0;
	}
	opt.end = OP.getInt("end");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "End point not specified, default to 60\n";
		opt.end = 60;
	}
	
	opt.helixGeoFile = OP.getString("helixGeoFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "helixGeoFile not specified, default to /data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
		opt.helixGeoFile = "/data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
	}

	opt.inputMonomerE = OP.getBool("inputMonomerE");
	if (OP.fail()) {
		opt.warningMessages += "monomer energy will be calculated\n";
		opt.warningFlag = true;
		opt.inputMonomerE = true;
	}
	opt.monoE_vdw = OP.getDouble("monoE_vdw");
	if (OP.fail()) {
		opt.monoE_vdw = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_hbond = OP.getDouble("monoE_hbond");
	if (OP.fail()) {
		opt.monoE_hbond = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_solv = OP.getDouble("monoE_solv");
	if (OP.fail()) {
		opt.monoE_solv = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_solvRef = OP.getDouble("monoE_solvRef");
	if (OP.fail()) {
		opt.monoE_solvRef= 1000000; //Default large, easy to spot error.
	}
	opt.printTermEnergies = OP.getBool("printTermEnergies");
	if (OP.fail()) {
		opt.printTermEnergies = true;
		opt.warningMessages += "printTermEnergies not specified using true\n";
		opt.warningFlag = true;
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
	opt.monoRotLibFile = OP.getString("monoRotLibFile");
	if (OP.fail()) {
		opt.warningMessages += "monoRotLibFile not specified using " + opt.rotLibFile + "\n";
		opt.warningFlag = true;
		opt.monoRotLibFile = opt.rotLibFile;
	}

	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine backboneCrd";
		opt.errorFlag = true;
	}
	
	opt.hBondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hBondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hbondFile not specified using " + opt.hBondFile + "\n";
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
	
	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine pdbOutputDir";
		opt.errorFlag = true;
	}

	opt.ids = OP.getStringVector("ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.weights = OP.getIntVector("weights");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA weights, defaulting each weight to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(10);
		}
	}
	if (opt.weights.size() != opt.ids.size()){
		opt.errorMessages += "Unable to identify alternate AA weights, make sure to correspond a weight to each AA\n";
		opt.errorFlag = true;
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
