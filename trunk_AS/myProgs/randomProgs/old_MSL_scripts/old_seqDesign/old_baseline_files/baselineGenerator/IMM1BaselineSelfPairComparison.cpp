#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

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
string programName = "IMM1BaselineSelfPairComparison";
string programDescription = "This program is a version of baselineSelfPairComparison_02 that calculates IMM1 baseline comparisons.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "4 December 2020";
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
	int pairDist;

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
	vector<double> baselines;
	vector<double> baselinesInner;
	vector<double> baselinesvdw;
	vector<double> baselinesvdwInner;

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
	return "A" + ps + "\nB" + ps;
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

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
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

vector<double> calcBaselineEnergies(System &_sys, int _seqLength){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		//cout << _sys.getPosition(i) << ": " << resi << endl;
		ener.push_back(resi);
	}
	sel.clearStoredSelections();
	return ener;
}
				
//vector<double> calcPairBaselineEnergies(System &_sys, Options _opt, int _seqLength, double &_totalProt){
vector<double> calcPairBaselineEnergies(System &_sys, Options _opt, int _seqLength){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
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

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+1);
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength; j++){
			int dist = j-i;
			if (dist <= _opt.pairDist){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+1);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
				//_totalProt += pair;
				//cout << residue << "; " << resi1 << " : " << pair << endl;
			//cout << _sys.getPosition(i) << ": " << _sys.getPosition(j) << ": " << pair << endl;
			}
			else{
				j = _seqLength;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

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

double seqSelfEnergy(string &_seq, map<string, double> _selfEner){
	double ener = 0;
	for (uint i=0; i<_seq.length(); i++){
		string aa(1, _seq[i]);
		for (map<string, double>::const_iterator it = _selfEner.begin(); it != _selfEner.end(); ++it){
			if (it->first == aa){
				ener = ener + it->second;
			}
		}
	}
	return ener;
}

double seqSelfPositionEnergy(string &_seq, map<string, vector<double>> _selfEner){
	double ener = 0;
	for (uint i=0; i<_seq.length(); i++){
		string aa(1, _seq[i]);
		for (map<string, vector<double>>::const_iterator it = _selfEner.begin(); it != _selfEner.end(); ++it){
			if (it->first == aa){
				ener = ener + it->second[i];		
			}
		}
	}
	return ener;
}

double seqPairEnergy(string &_seq, Options _opt, map<string, vector<double>> _pairEner){
	double ener = 0;
	vector<string> pairs1;
	vector<string> pairs2;
	vector<int> dists;
	for (uint i=0; i<_seq.length(); i++){
		for (uint j=i+1; j<_seq.length(); j++){
			stringstream ss1;
			stringstream ss2;
			string aa1(1, _seq[i]);
			string aa2(1, _seq[j]);
			ss1 << aa1 << aa2;
			ss2 << aa2 << aa1;
			string key1 = ss1.str();
			string key2 = ss2.str();
			pairs1.push_back(key1);
			pairs2.push_back(key2);
			dists.push_back(j-i);
		}
	}
	for (uint i=0; i<pairs1.size(); i++){
		for (map<string, vector<double>>::const_iterator it = _pairEner.begin(); it != _pairEner.end(); ++it){
			if (it->first == pairs1[i] || it->first == pairs2[i]){
				ener = ener + it->second[dists[i]-1];
				it == _pairEner.end();
			}
		}
	}
	return ener;
}

vector<double> printSelfBaselineEnergy(System &_sys, map<string, map<int,double>> _selfEner){
	vector<double> ener;

	for (uint i=0; i<_sys.chainSize(); i++){
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for (vector<Position*>::iterator p = positions.begin(); p!= positions.end(); p++){
			Residue &res = (*p)->getCurrentIdentity();
			string baseId = MslTools::getOneLetterCode(res.getResidueName());
			uint pos = p-positions.begin()+1;
			try{
				double tmp = _selfEner.at(baseId).at(pos);
				//cout << baseId << ": " << pos << ": " << tmp << endl;
				ener.push_back(tmp);
			}
			catch (const out_of_range& e){
				continue;
			}
		}
	}
	return ener;
}

vector<double> printPairBaselineEnergy(System &_sys, Options _opt, map<string,map<int,map<int,double>>> &_pairMap){
	vector<double> ener;
	vector<string> pairs1;
	vector<string> pairs2;
	vector<int> dists;
	for (uint i=0; i<_sys.chainSize(); i++){
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for (vector<Position*>::iterator p = positions.begin(); p!= positions.end(); p++){
			Residue &res1 = (*p)->getCurrentIdentity();
			string baseId1 = MslTools::getOneLetterCode(res1.getResidueName());
			int pos = p-positions.begin()+1;
			for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
				uint d = p2-p;
				if (d <= 10){
					Residue &res2 = (*p2)->getCurrentIdentity();
					string baseId2 = MslTools::getOneLetterCode(res2.getResidueName());
					stringstream ss;
					ss << baseId1 << baseId2;
					string key = ss.str();
					double energy = _pairMap.at(key).at(pos).at(d);
					//cout << key << ": " << pos << ": " << d << ": " << energy << endl;
					ener.push_back(energy);
				//	uint order = baseId1.compare(baseId2);
				//	if (order == 0){
				//		try{
				//			stringstream ss;
				//			ss << baseId1 << baseId2;
				//			string key = ss.str();
				//			cout << key << ": " << pos << ": " << d << endl;
				//			double energy = _pairMap.at(key).at(pos).at(d);
				//	cout << key << ": " << pos << ": " << d << ": " << energy << endl;
				//			ener.push_back(energy);
				//		}
				//		catch (const out_of_range& e){
				//			continue;
				//		}
				//	}
				//	else{
				//		set<string> sortAAs;
				//		sortAAs.insert(baseId1);
				//		sortAAs.insert(baseId2);
				//		vector<string> sort;
				//		for (set<string>::iterator s = sortAAs.begin(); s != sortAAs.end(); ++s){
				//			string a = *s;
				//			sort.push_back(a);
				//		}
				//		try{
				//			stringstream ss;
				//			ss << sort[0] << sort[1];
				//			string key = ss.str();
				//			cout << key << ": " << pos << ": " << d << endl;
				//			double energy = _pairMap.at(key).at(pos).at(d);
				//	cout << key << ": " << pos << ": " << d << ": " << energy << endl;
				//			ener.push_back(energy);
				//		}
				//		catch (const out_of_range& e){
				//			continue;
				//		}
				//	}
				}
			}
		}
	}
	return ener;
}

map<string, map<int,double>> makeBaselineMap(string _file){
	map<string, map<int,double>> m;
	map<int,double> m2;
	vector<string> line;
	string AA = "";
	vector<int> pos;
	vector<double> ener;
	ifstream fs;
	fs.open(_file.c_str());
	if (!fs.is_open()){
		cerr << "Could not open baseline file" << endl;
	}
	else{
		string s;
		while(getline(fs, s)){
			istringstream iss(s);
			copy((istream_iterator<string>(iss)), istream_iterator<string>(), back_inserter(line));
		}
		fs.close();
	}
	uint count = 0;
	string prevAA = "";
	//separate the input of each line into two separate vectors
	for (uint i=0; i<line.size(); i++){
		if (i == 0){
			AA = line[i];
			prevAA = AA;
		}
		if (prevAA == AA){
			if (count == 0){
				AA = line[i];
				count++;
			}
			else{
				if (count == 1){
					pos.push_back(MslTools::toInt(line[i]));
					count++;
				}
				else{
					ener.push_back(MslTools::toDouble(line[i]));
					count = 0;
				}
			}

		}
		else{
			for (uint j=0; j<pos.size(); j++){
				m2.insert(make_pair(pos[j], ener[j]));
				//cout << prevAA << ": " << pos[j] << ": " << ener[j] << endl;
			}
			m.insert(make_pair(prevAA, m2));
			prevAA = AA;
			pos.clear();
			ener.clear();
			m2.clear();
			i = i-1;//so that it doesn't skip a value on the line which it would as currently written
		}
		if (i == line.size()-1){
			for (uint j=0; j<pos.size(); j++){
				m2.insert(make_pair(pos[j], ener[j]));
				//cout << prevAA << ": " << pos[j] << ": " << ener[j] << endl;
			}
			m.insert(make_pair(prevAA, m2));
		}
	}

	return m;
}

map<string, map<int,map<int,double>>> makeBaselineMapPair(string _file){
	map<string, map<int,map<int,double>>> m;
	map<int,map<int,double>> m1;
	map<int,double> m2;
	vector<string> line;
	string AA = "";
	int pos = 0;
	vector<double> dist;
	vector<double> ener;
	ifstream fs;
	fs.open(_file.c_str());
	if (!fs.is_open()){
		cerr << "Could not open baseline file" << endl;
	}
	else{
		string s;
		while(getline(fs, s)){
			istringstream iss(s);
			copy((istream_iterator<string>(iss)), istream_iterator<string>(), back_inserter(line));
		}
		fs.close();
	}
	int count = 0;
	string prevAA = "";
	int prevPos = 0;

	//separate the input of each line into separate and values vectors
	for (uint i=0; i<line.size(); i++){
		if (i == 0){
			AA = line[i];
			prevAA = AA;
			pos = MslTools::toInt(line[i+1]);
			prevPos = pos;
		}
		if (prevAA == AA){
			if (count == 0){
				AA = line[i];
				count++;
			}
			else{
				if (prevPos == pos){
					if (count == 1){
						pos = MslTools::toInt(line[i]);
						count++;
					}
					else{
						if (count == 2){
							dist.push_back(MslTools::toInt(line[i]));
							count++;
						}
						else{
							ener.push_back(MslTools::toDouble(line[i]));
							count = 0;
						}
					}
				}
				else{
					for (uint j=0; j<dist.size(); j++){
						m2.insert(make_pair(dist[j], ener[j]));
						//cout << prevAA << ": "<< prevPos << ": " << dist[j] << ": " << ener[j] << endl;
					}
					m1.insert(make_pair(prevPos, m2));
					prevPos = pos;
					dist.clear();
					ener.clear();
					m2.clear();
					i = i-1;//so that it doesn't skip a value on the line which it would as currently written
				}
			}
		}
		else{
			for (uint k=0; k<dist.size(); k++){
				m2.insert(make_pair(dist[k], ener[k]));
				//cout << prevAA << ": "<< prevPos << ": " << dist[j] << ": " << ener[j] << endl;
			}
			m1.insert(make_pair(prevPos, m2));
			m.insert(make_pair(prevAA, m1));
			prevAA = AA;
			pos = MslTools::toInt(line[i]);
			prevPos = pos;
			dist.clear();
			ener.clear();
			m1.clear();
			m2.clear();
			count++;
		}
		if (i == line.size()-1){
			for (uint j=0; j<dist.size(); j++){
				m2.insert(make_pair(dist[j], ener[j]));
				//cout << prevAA << ": "<< prevPos << ": " << dist[j] << ": " << ener[j] << endl;
			}
			m1.insert(make_pair(prevPos, m2));
			m.insert(make_pair(prevAA, m1));
		}
	}
	return m;
}

double sumEnergyVector(vector<double> _energies){
	double ener = 0;
	for (uint i=0; i<_energies.size(); i++){
		ener = ener + _energies[i];
	}
	return ener;
}

void printEnerFile(string _seq, Options _opt, vector<double> &_ener, vector<double> &_bEner, int _seqNumber, bool _selfEner, ofstream &_out){
	if (_selfEner ==true){
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "AA:     Position:     Energy:     Baseline" << endl;
		}
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << ":    " << i+1 << ":        " << _ener[i] << ":      " << _bEner[i] << endl;
		}
	}
	else{
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "Pair: Distance: Position1: Position2: Energy: Baseline" << endl;//could be interesting to add rotamer number to this
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
			_out << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << _bEner[i] << endl;
		}
	}
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
	ofstream bsout;
	ofstream sout;
	ofstream pout;
	
	opt.pdbOutputDir = opt.pdbOutputDir + "/" + date + "/";
	string dir = opt.pdbOutputDir;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string bsoutName  = dir + "/baselineEnergySummary_" + opt.runNumber + ".out";
	string soutName   = dir + "/selfEnergyComparison_" + opt.runNumber + ".out";
	string poutName   = dir + "/pairEnergyComparison_" + opt.runNumber + ".out";
	
	bsout.open(bsoutName.c_str());
	sout.open(soutName.c_str());
	pout.open(poutName.c_str());
	
	/******************************************************************************
	 *                     === READ IN BASELINE FILES AS MAPS ===
	 ******************************************************************************/
	map<string, map<int,double>> selfNoPos       = makeBaselineMap("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/IMM1_meanSelf_noRef.txt");
	double selfOuterLeu = 1.2329;
	map<string, map<int, map<int,double>>> pairMap   = makeBaselineMapPair("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/IMM1_meanPair.txt");
	
	//for (map<string, map<int, map<int,double>>>::const_iterator i = pairMap.begin(); i != pairMap.end(); i++){
	//	//cout << i->first << ": " << i->second.first << ": " << i->second->second->first << ": " << i->second->second->second << endl;
	//	map<int, map<int, double>> pm2 = i->second;
	//	for (map<int, map<int,double>>::const_iterator j = pm2.begin(); j != pm2.end(); j++){
	//		map<int,double> pm3 = j->second;
	//		for (map<int, double>::const_iterator k = pm3.begin(); k != pm3.end(); k++){
	//			cout << i->first << ": " << j->first << ": " << k->first << ": " << k->second << endl;
	//		}
	//	}
	//}

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
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
	writer.open(dir + "/polyLeu.pdb");
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	//vector<int> highEnergySequences;
	int seqDiscard = 0;
	int seqAccept = 0;
	uint a = 0;
	int seqNumber = opt.seqNumber;
	
	//Random Number Generator
	RandomNumberGenerator RNG;
	//RNG.setTimeBasedSeed();
	RNG.setSeed(opt.seed);
	
	/******************************************************************************
	 *                     === GIVE SEQUENCES WEIGHTS ===
	 ******************************************************************************/
	vector<string> weightedAAs = weightAAs(opt.ids, opt.weights, opt.backboneLength, RNG);
	
	if (opt.backboneLength > weightedAAs.size()){
		cout << "ERROR: Length of protein is too long for total number of AAs given." << endl;
		exit(1);
	}
	
	/******************************************************************************
	 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
	 ******************************************************************************/
	while (seqAccept < seqNumber){
		cout << "Number Accepted: " << seqAccept << endl;
		cout << "Number Discarded: " << seqDiscard << endl;
		/******************************************************************************
		 *                     === RANDOMIZE SEQUENCES ===
		 ******************************************************************************/
		string seq = randomAASeqLeuEnds(weightedAAs, opt.backboneLength, RNG);
		cout << seq << endl;
		string polySeq = convertToPolymerSequenceNeutralPatch(seq, 1);
		//string polySeq = convertToPolymerSequence(seq, 1);
		PolymerSequence PS(polySeq);
		cout << PS << endl;

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
		CSB.setIMM1Params(15,10);
		
		CSB.setBuildNonBondedInteractions(false);	
		if(!CSB.buildSystem(PS)) {
			cerr << "Unable to build system from " << polySeq << endl;
			exit(0);
		} else {
			//fout << "CharmmSystem built for sequence" << endl;
		}
		
		/******************************************************************************
		 *                     === CALCULATE PAIR ENERGIES ===
		 ******************************************************************************/
		vector<double> selfBaselines = printSelfBaselineEnergy(sys, selfNoPos);
		double seqSelfEner = sumEnergyVector(selfBaselines);
		vector<double> pairBaselines = printPairBaselineEnergy(sys, opt, pairMap);
		double seqPairEner = sumEnergyVector(pairBaselines);
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();
		
		CSB.updateNonBonded(10,12,50);	
		
		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
		sysRot.defineRotamerSamplingLevels();
		
		// Add hydrogen bond term
		HydrogenBondBuilder hb(sys, opt.hBondFile);
		hb.buildInteractions(30);
	
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
		 *           === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
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
		spm.getMinStates()[0];
		spm.updateWeights();
		spm.setOnTheFly(true);
		spm.calculateEnergies();
		
		repackSideChains(spm, 10);
		
		/******************************************************************************
		 *                  === SET SYSTEM TO BEST SPM ROTAMERS ===
		 ******************************************************************************/
		sys.setActiveRotamers(spm.getMinStates()[0]);
		sys.calcEnergy();
		double totalEnergy = sys.calcEnergy();
		//double vdw = Eset->getTermEnergy("CHARMM_VDW");
		//double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
	
		/******************************************************************************
		 *             === EXTRACT SELF PAIR ENERGIES FOR PRINTING ===
		 ******************************************************************************/
		vector<double> selfEnerVec = calcBaselineEnergies(sys, seq.length()); //reapplying this function from generateSelfPairBaselinevdW_01.cpp to extract self and pair energies for e`ach sequence	
		vector<double> pairEnerVec = calcPairBaselineEnergies(sys, opt, seq.length()); //reapplying this function from generateSelfPairBaselinevdW_01.cpp to extract self and pair energies for e`ach sequence	
		/******************************************************************************
		 *               === CALCULATE ENERGIES FOR EACH POSITION ===
		 ******************************************************************************/
		//set up way to select for the 22 inner residues (selection seems a little weird to me, so I think I'm going to try the alternate way subtracting from the outermeans first)
		totalEnergy = sys.calcEnergy();
		cout << "Total Energy: " << totalEnergy << endl;
		cout << sys.getEnergySummary() << endl;

		double selfEner = sumEnergyVector(selfEnerVec);
		double pairEner = sumEnergyVector(pairEnerVec);
		
		double seqBaselineEnergy = seqSelfEner + seqPairEner;
		cout << "Total Energy: " << totalEnergy << endl;
		cout << "Self Energy: " << selfEner << endl;
		cout << "Pair Energy: " << pairEner << endl;
		cout << "Baseline Self Energy: " << seqSelfEner << endl;
		cout << "Baseline Pair Energy: " << seqPairEner << endl;
		cout << "SelfPair Energy: " << seqBaselineEnergy << endl;
		cout << "Total - SelfPair Energy: " << totalEnergy-seqBaselineEnergy << endl;
		
		double selfPairDif = totalEnergy - seqBaselineEnergy;

		/******************************************************************************
		 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
		 ******************************************************************************/
		if (seqAccept == 0 && opt.seed == 0){
			bsout << "Monomer: SelfPair: Difference:" << endl;
		}
		bsout << totalEnergy << ":    " << seqBaselineEnergy << ":   " << selfPairDif << endl;
		if (pairEnerVec.size() > pairBaselines.size()){
			while (pairEnerVec.size() != pairBaselines.size()){
				pairBaselines.push_back(0);//for some reason my code isn't getting consistent numbers of baselines; this is a cheat until I fix the functions (doesn't work; check YA on line 1570 to be sure if my function fixes work or not
			}
		}
		else{
			if (pairEnerVec.size() < pairBaselines.size()){
				cout << "Pair size is not equal to baseline size: Pair " << pairEnerVec.size() << "vs Baseline " << pairBaselines.size() << "." << endl;
				break;
			}
		}
		//printEnerFile(seq, opt, selfEnerVec, selfBaselines, seqAccept, true, sout);
		printEnerFile(seq, opt, pairEnerVec, pairBaselines, seqAccept, false, pout);
		seqAccept++;
		sys.reset();
	}
	time (&endTime);
	diffTime = difftime(endTime, startTime);
	cout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	writer.close();
	bsout.close();
	sout.close();
	pout.close();
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
	opt.allowed.push_back("baselines");
	opt.allowed.push_back("baselinesInner");
	opt.allowed.push_back("baselinesvdw");
	opt.allowed.push_back("baselinesvdwInner");
	
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
	opt.baselines = OP.getDoubleVector("baselines");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA baselines, defaulting each baseline to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(0);
		}
	}
	opt.baselinesInner = OP.getDoubleVector("baselinesInner");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA baselines, defaulting each baseline to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(0);
		}
	}
	opt.baselinesvdw = OP.getDoubleVector("baselinesvdw");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA baselinesvdw, defaulting each baseline to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(0);
		}
	}
	opt.baselinesvdwInner = OP.getDoubleVector("baselinesvdwInner");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA baselinesvdw, defaulting each baseline to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(0);
		}
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
