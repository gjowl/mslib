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
#include "BaselineEnergyBuilder.h"
#include "BaselineInteraction.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "vdwSequenceDesign";
string programDescription = "This program designs sequences reliant on van der Waals forces onto the interfacial residues of a given GASright geometry";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "11 October 2019";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************************************************************************************************************************************************/

struct Options{
	string sequence;
	int sequenceLength;
	string backboneAA;
	int backboneLength;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int sequenceStart;

	int threadStart;
	int threadEnd;
	bool threadBool;
	
	bool deleteTerminalHbonds;
	
	string SL; //number of rotamers

	// transformation
	double xShift;
	double zShift;
	double crossingAngle;
	double axialRotation;
	bool transform;
	int thread;
	int bbThread;

	// input files
	string helixGeoFile;
	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string baselineFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string monoRotLibFile;
	string infile;
	string rulesFile;

	// side-chain repack variable
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;

	bool verbose;
	int greedyCycles;
	int seed;

	int numberOfStructuresToMCRepack;
	double energyCutOff;
	
	// weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	bool inputMonomerE;
	int monoE_vdw;
	int monoE_hbond;
	int monoE_solv;
	int monoE_solvRef;

	// clustering options (optional)
	double rmsdCutoff;
	bool printAllCrds;
	bool printAxes;
	bool printTermEnergies;

	// alternate identities
	vector<string> Ids;
	int numPositions;

	// protein information (optional)
	string uniprotName;
	string uniprotAccession;

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

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
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
	return "A" + ps + "\nB" + ps;
}

string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}

string generatePolyLeu(string _backboneAA, int _sequenceLength) {
	string polyLeu = "";
	for (uint i=0; i<_sequenceLength; i++){
		polyLeu = polyLeu + _backboneAA;
	}
	return polyLeu;
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

bool rulesCheck(System & _sys, string _geoIndex, map<int, string> _rulesMap) {

	int idx = MslTools::toInt(_geoIndex);
	if (_rulesMap.find(idx) == _rulesMap.end()) {
		return true;
	}

	vector<string> parsedRules = MslTools::tokenize(_rulesMap[idx], ",");
	// read over all the rules for a given model
	for (uint m = 1; m < parsedRules.size(); m++) {
		int delimiter  = parsedRules[m].find_first_of(":");
		string position = parsedRules[m].substr(0,1)+","+parsedRules[m].substr(1,delimiter-1); // makes string A,35
		string rule = parsedRules[m].substr(delimiter+1,(parsedRules[m].length()-delimiter-1)); 
		if(!_sys.positionExists(position)) {
			continue;
		}
		//fout << "position " << position << " must be " << rule << " , " << " is actually " << _sys.getPosition(position).getResidueName() << endl;

		// Check if there is a ! "not"
		size_t findNot = rule.find("!");

		// if there is a "not"
		if (findNot == 0) { // ! will only appear in position 0
			rule = rule.substr(2);
			rule = MslTools::trim(rule, "]"); // residues that cannot exist in given position
			// look at position given in rule	
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					return false;
				}
			}
			
		}
		else { // no "not" found
			rule = rule.substr(1);
			rule = MslTools::trim(rule, "]"); // residues that must exist in given position
			// look at position given in rule	
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			bool found = false;
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					found = true;
					break;
				}
			}
			if(!found) {
				return false;
			}
		}
	}
	
	return true;
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

Options parseOptions(int _argc, char * _argv[], Options defaults);

void printOptions(Options & _op, ofstream & _fout) {
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
	_fout << "baselineFile " << _op.baselineFile << endl;
	_fout << "solvFile " << _op.solvFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;
	_fout << "monoRotLibFile " << _op.monoRotLibFile << endl;

	_fout << "MCCycles " << _op.MCCycles << endl;
	_fout << "MCMaxRejects " << _op.MCMaxRejects << endl;
	_fout << "MCStartTemp " << _op.MCStartTemp << endl;
	_fout << "MCEndTemp " << _op.MCEndTemp << endl;
	_fout << "MCCurve " << _op.MCCurve << endl;

	_fout << "deltaZ " << _op.deltaZ << endl;
	_fout << "deltaAx " << _op.deltaAx << endl;
	_fout << "deltaCross " << _op.deltaCross << endl;
	_fout << "deltaX " << _op.deltaX << endl;

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

	_fout << "rmsdCutoff " << _op.rmsdCutoff << endl;
	_fout << "printAllCrds " << _op.printAllCrds << endl;
	_fout << "printAxes " << _op.printAxes << endl;
	_fout << "printTermEnergies " << _op.printTermEnergies << endl;
	_fout << "deleteTerminalHbonds " << _op.deleteTerminalHbonds << endl;

	_fout << "fullSequenceStart " << _op.sequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "threadStart " << _op.threadStart << endl;
	_fout << "threadEnd " << _op.threadEnd << endl;

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

//below function only loads rotamers onto the interfacial positions by varPos (01 where 0 = non-interfacial and 1 = interfacial)
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, vector<int> _interface){
	for (uint k=0; k<_interface.size(); k++) {
		if (_interface[k] == 1){
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
}

void loadInterfacialPeripheralRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, vector<int> _interface){
	for (uint k=0; k<_interface.size(); k++) {
		if (_interface[k] == 1){
			Position &posA = _sys.getPosition(k-1);
			Position &posB = _sys.getPosition(k+1);
			if (posA.identitySize() > 1){
				for (uint j=0; j < posA.getNumberOfIdentities(); j++){
					posA.setActiveIdentity(j);
					if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
						if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), _SL)) {
							cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
						}
					}
					posA.setActiveIdentity(0);
				}
			}
			else{
				if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
					}
				}
			}
			if (posB.identitySize() > 1){
				for (uint j=0; j < posB.getNumberOfIdentities(); j++){
					posB.setActiveIdentity(j);
					if (posB.getResidueName() != "GLY" && posB.getResidueName() != "ALA" && posB.getResidueName() != "PRO") {
						if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), _SL)) {
							cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
						}
					}
					posB.setActiveIdentity(0);
				}
			}
			else{
				if (posB.getResidueName() != "GLY" && posB.getResidueName() != "ALA" && posB.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
					}
				}
			}
		}	
	}
}

void identifyInterface(System &_sys, vector<int> &_pos, vector<double> &_dists, int _numPositions){
	int count = 0;
	for (uint k=0; k<_sys.getChain("A").positionSize(); k++) {
		//lout << _sys.positionSize() << endl;
		Position &pos1 = _sys.getChain("A").getPosition(k);
		Position &pos2 = _sys.getChain("B").getPosition(k);
	//cout << pos1 << endl;	
	//cout << pos2 << endl;	
		Atom &c1 = pos1.getAtom("CA");
		Atom &c2 = pos2.getAtom("CA");
		
		double dist;
		double most;
		double x = 0;
		dist = c1.distance(c2);
		//lout << "Dist " << k << ": " << dist << endl;
		if (_pos.size() < _numPositions){
			_dists.push_back(dist);
			_pos.push_back(k);
		}
		else {
			count++;
			if (count <= _numPositions/2){
				most = _dists[0];
				x = 0;
				//lout << "Current: " << endl;
				for (uint j=0; j<8; j++){
					//lout << _pos[j] << ": " << _dists[j] << endl;
					if (most < _dists[j]){
						most = _dists[j];
						x = j;
					}
				}
				if (dist <= most){
					//lout << "Before: " << x << ": " << _dists[x];
					_dists[x] = dist;
					_pos[x] = k;
					//lout << " ; After: " << _pos[x] << ": " << _dists[x] << endl;
					count = 0;
				}
			}
			else{
				cout << "Interface Identified!" << endl;
				k = _sys.getChain("A").positionSize();
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
		varPos[_pos[j]+_sys.getChain("A").positionSize()] = 1;
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

//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

vector<vector<string>> positionToString(System &_sys, vector<int> &_varPos){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_varPos.size()/2; k++){
		//lout << "string" << endl;
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

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}


double computeMonomerEnergy(System & _sys, Transforms & _trans, Options& _opt, System & _helicalAxis, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects, bool _useIMM1) {

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);	
	
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
	monoEset->setTermActive("SCWRL4_HBOND", false);

	if (_useIMM1 == false){
		monoEset->setTermActive("CHARMM_IMM1REF", false);
		monoEset->setTermActive("CHARMM_IMM1", false);
		monoEset->setWeight("CHARMM_IMM1REF", 0);
		monoEset->setWeight("CHARMM_IMM1", 0);
	}

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
	monoEset->setWeight("SCWRL4_HBOND", 0);
	monoEset->setWeight("CHARMM_IMM1REF", 1);
	monoEset->setWeight("CHARMM_IMM1", 1);
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
	_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0); 
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

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
			_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
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
	string monoOutCrdFile  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}
	
	string monoOutPdbFile  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.pdb";
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
		_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

		ofstream meOut;
		string meOutName  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.energy";
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

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;
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

System localMC(System &_pdb, System &_helicalAxis, Options &_opt, string _sequence, double &_monomer, ofstream &_out, PDBWriter &_writer){

	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequence(_sequence, _opt.thread);

	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile,_opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", false);
	CSB.setBuildTerm("CHARMM_IMM1", false);
	
	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);
	
	if(!CSB.buildSystem(polymerSeq)) {
		cerr << "Unable to build system from " << polymerSeq << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");
		
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

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
	Eset->setTermActive("SCWRL4_HBOND", false);
	
	/******************************************************************************
	 *             === SETUP IMM1 ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	if (_opt.useIMM1 == true){
		_out << "Using IMM1 Energies for Monomer Comparison" << endl << endl;
		Eset->setTermActive("CHARMM_IMM1REF", true);
		Eset->setTermActive("CHARMM_IMM1", true);
		Eset->setWeight("CHARMM_IMM1REF", 1);
		Eset->setWeight("CHARMM_IMM1", 1);
	}
	else{
		_out << "Not using IMM1 Energies for Monomer Comparison" << endl << endl;
		Eset->setTermActive("CHARMM_IMM1REF", false);
		Eset->setTermActive("CHARMM_IMM1", false);
		Eset->setWeight("CHARMM_IMM1REF", 0);
		Eset->setWeight("CHARMM_IMM1", 0);
	}
	
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 0);

	_out << "Energy Weights: " << endl;
	_out << "VDW weight: " << Eset->getWeight("CHARMM_VDW") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
	_out << endl;
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	
	beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
	beb.readParameters(_opt.baselineFile);
	//beb.printParameters();// Ensures they are being read properly
	
	beb.buildInteractions();	

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.wipeAllCoordinates();
	//sys.assignCoordinates(glyAPV,false);
	sys.assignCoordinates(_pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);//don't think I need this?
	
	//writer.write(sys.getAtomPointers(), true, false, true);
	//writer2.write(pdb.getAtomPointers(), true, false, true);
	
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 
	//RNG.setTimeBasedSeed();

	CSB.updateNonBonded();
	sys.buildAllAtoms();
	
	loadRotamers(sys, sysRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm1;
	spm1.seed(RNG.getSeed());
	spm1.setSystem(&sys);
	spm1.setVerbose(true);
	spm1.getMinStates()[0];
	spm1.updateWeights();
	spm1.setOnTheFly(true);
	spm1.calculateEnergies();

	repackSideChains(spm1, _opt.greedyCycles);

	sys.setActiveRotamers(spm1.getMinStates()[0]);
	sys.saveAltCoor("savedBestState");
	_helicalAxis.saveAltCoor("BestAxis");
	cout << sys.calcEnergy() << endl;
	
	Eset->eraseTerm("BASELINE");//gets rid of the baseline energy term
	cout << sys.calcEnergy() << endl;
	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	//Was previously generating a new axis, but realized I still have one from my initial transformation to geometry, so I'm just gonna use that one (energies were strange so hopefully this helps?)
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();
	
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	//xShiftTransformation(apvChainA, apvChainB, axisA, axisB, _opt.xShift, trans);//this fixes the axis problem I was having, allowing MC to complete
	
	/******************************************************************************
	 *              === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	map<string,double> monomerEnergyByTerm;
	double monomerEnergy;
	
	_out << "Monomer calculations..." << endl;
	cout << "Monomer calculations..." << endl;
	monomerEnergy = computeMonomerEnergy(sys, trans, _opt, _helicalAxis, RNG, monomerEnergyByTerm, _out, _opt.greedyCycles, _opt.MCCycles, _opt.MCMaxRejects, _opt.useIMM1);
	_monomer = monomerEnergy;
	_out << "End monomer calculations! Monomer Energy: " << monomerEnergy << endl;
	cout << "End monomer calculations! Monomer Energy: " << monomerEnergy << endl;

	//sys.setActiveRotamers(spm.getMinStates()[0]);

	//double currentEnergy = spm.getMinBound()[0];

	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	double xShift = _opt.xShift;
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	_out << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	_out << " xShift: " << xShift << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	cout << " xShift: " << xShift << endl;

	// Transform helices to initial starting position
	//transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), trans);

	cout << sys.calcEnergy() << endl;
	// Optimizatize Initial Starting Position
	//repackSideChains(spm, _opt.greedyCycles);
	sys.setActiveRotamers(spm1.getMinStates()[0]);
	cout << sys.calcEnergy() << endl;
	//double currentEnergy = spm1.getMinBound()[0];
	double currentEnergy = sys.calcEnergy();
	cout << currentEnergy << endl;

	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.calculateEnergies();
	
	/******************************************************************************
	 *                     === X SHIFT REPACKS ===
	 ******************************************************************************/
	double bestEnergy = currentEnergy;
	double savedXShift = xShift;
	double previousEnergy = monomerEnergy;
	double deltaXShift = -0.1;
	double globalLowestE = monomerEnergy;
	//double xShiftEnd = MslTools::toDouble(parsedGeoInformation[4]);
	double xShiftEnd = 9.5;//TODO: probably should change this to an option, or xShift-0.5 or something of the like

	while (xShift >= xShiftEnd) {
	
		xShift += deltaXShift;
		// Move the helix
		backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );
		// Run Optimization
		repackSideChains(spm, _opt.greedyCycles);
	
		vector<unsigned int> MCOFinal;
		MCOFinal = spm.getMinStates()[0];
		sys.setActiveRotamers(MCOFinal);
	
		cout << Eset->getSummary() << endl;
		currentEnergy = spm.getMinBound()[0];
	
		if (currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			savedXShift = xShift;
			sys.saveAltCoor("savedBestState");
			_helicalAxis.saveAltCoor("BestAxis");
		}
	
		_out << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
		cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

		// If energy increase twice in a row, and it is above the global lowest energy, quit
		if (currentEnergy < globalLowestE) {
			globalLowestE = currentEnergy;
		}
		if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
			_out << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
			break;	
		}
		else {
			previousEnergy = currentEnergy;
		}
	
	}
	_out << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl; 
	cout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl; 
	xShift = savedXShift;

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	_out << "====================================" << endl;
	_out << "Performing Local Monte Carlo Repacks" << endl;
	_out << "====================================" << endl;
	_out << endl;
	cout << "====================================" << endl;
	cout << "Performing Local Monte Carlo Repacks" << endl;
	cout << "====================================" << endl;
	cout << endl;

//	double bestEnergy = currentEnergy;
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	//bestEnergy = currentEnergy;

	vector<unsigned int> MCOBest = spm.getMinStates()[0];
	
	if (_opt.MCCycles > 0) {
		//MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);
		MonteCarloManager MCMngr(0, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);

		MCMngr.setEner(bestEnergy);
		
		while(!MCMngr.getComplete()) {

			sys.applySavedCoor("savedBestState");
			_helicalAxis.applySavedCoor("BestAxis");

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
			repackSideChains(spm, _opt.greedyCycles);
	
			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
			currentEnergy = sys.calcEnergy();
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				_out << "state rejected   energy: " << currentEnergy << endl;
			}
			else {
				bestEnergy = currentEnergy;
				sys.saveAltCoor("savedBestState");
				_helicalAxis.saveAltCoor("BestAxis");

				xShift = xShift + deltaXShift;
				crossingAngle = crossingAngle + deltaCrossingAngle;
				axialRotation = axialRotation + deltaAxialRotation;
				zShift = zShift + deltaZShift;
				MCOBest = MCOFinal;

				_out << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_out << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-monomerEnergy << endl;
	cout << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-monomerEnergy << endl;
	_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

	sys.applySavedCoor("savedBestState");
	_writer.write(sys.getAtomPointers(),true,false,true);
	
	double finalEnergy = sys.calcEnergy()-monomerEnergy;
	return sys;
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
	
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);
	
	time(&startTime);
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);
	
	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	ofstream fout;
	ofstream eout;
	ofstream mout;
	ofstream pout;
	
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		eout << "Unable to make directory" << endl;
		//exit(0);
	}

	string outfile = dir + "/vdw_" + opt.runNumber + ".out";
	string errfile = dir + "/vdw_" + opt.runNumber + ".err";
	string moutfile = dir + "/monomer_" + opt.runNumber + ".out";
	string poutfile = dir + "/options.out";

	fout.open(outfile.c_str());
	eout.open(errfile.c_str());
	mout.open(moutfile.c_str());
	pout.open(poutfile.c_str());

	fout << date << endl;
	eout << date << endl;
	mout << date << endl;
	pout << date << endl;

	printOptions(opt, pout);
	
	if (opt.errorFlag) {
		eout << endl;
		eout << "The program terminated with errors:" << endl;
		eout << endl;
		eout << opt.errorMessages << endl;
		eout << endl;
		eout << opt.OPerrors << endl;

		usage();
		exit(1);
	}
	
	/******************************************************************************
	 *                     === READ IN GEOMETRY FILE ===
	 ******************************************************************************/
	vector<string> fileVec;
	readGeometryFile(opt.helixGeoFile, fileVec);

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd); 
	if(!cRead.read()) {
		eout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === GENERATE POLYGLY ===
	 ******************************************************************************/
	string polyGly = generatePolymerSequence("G", opt.sequenceLength, opt.thread);
	
	PolymerSequence PS(polyGly);

	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	PDBWriter writer;
	PDBWriter writer2;
	PDBWriter writer3;
	PDBWriter writer4;
	
	writer.open(dir + "/polyGly_" + opt.runNumber + ".pdb");
	writer2.open(dir + "/helixAlignToCoord_" + opt.runNumber + ".pdb");
	writer3.open(dir + "/polyLeu_" + opt.runNumber + ".pdb");
	if (opt.runSCMF){
		writer4.open(dir + "/SCMFstate_" + opt.runNumber + ".pdb");
	}
	else{
		writer4.open(dir + "/BestUnbiasedMCState_" + opt.runNumber + ".pdb");
	}

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	CSB.setBuildNoTerms();
	
	if(!CSB.buildSystem(PS)) {
		eout << "Unable to build system from " << polyGly << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	vector<Position*>& positions = sys.getPositions();

	//writer.write(sys.getAtomPointers(), true, false, true);//for some reason this only builds part of chainB? even though when printed out the atomPointer are correct?
	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(opt.infile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
	sys.wipeAllCoordinates();
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	writer.write(sys.getAtomPointers(), true, false, true);//for some reason this only builds part of chainB? even though when printed out the atomPointer are correct?
	//writer.write(apvChainA, true, false, true);//for some reason this only builds part of chainB? even though when printed out the atomPointer are correct?
	//writer.write(apvChainB, true, false, true);//for some reason this only builds part of chainB? even though when printed out the atomPointer are correct?
	string seqA = convertPolymerSeqToOneLetterSeq(chainA);
	string seqB = convertPolymerSeqToOneLetterSeq(chainB);

	//writer.close();
	
	//cout << sys.getAtomPointers() << endl;

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
		eout << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	//writer2.write(sys.getAtomPointers(), true, false, true);
	
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);//TODO: changed this line to reflect the config options, eventually should probably change back to datafile (or extracted geometry from pdb)
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	sys.buildAllAtoms();
	writer2.write(sys.getAtomPointers(), true, false, true);

	/******************************************************************************
	 *                   === IDENTIFY INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	vector<int> pos;
	vector<double> dist;
	vector<vector<vector<double>>> posDistVector;
	identifyInterface(sys, pos, dist, opt.numPositions);
	vector<int> organizedPos;
	vector<double>organizedDist;
	//Organizes the position and distance combo from lowest position to largest
	for (uint i=0; i<pos.size(); i++){
		organizedPos.push_back(pos[i]);
		organizedDist.push_back(dist[i]);
		if (organizedPos.size() > 1){
			int count = organizedPos.size()-1;
			double tempDist = 0;
			int tempPos = 0;
			while (count > 0){
				if (organizedPos[count] < organizedPos[count-1]){
					tempPos = organizedPos[count];
					tempDist = organizedDist[count];
					organizedPos[count] = organizedPos[count-1];
					organizedDist[count] = organizedDist[count-1];
					organizedPos[count-1] = tempPos;
					organizedDist[count-1] = tempDist;
					count -= 1;
				}
				else{
					count -= 1;
				}
			}
		}
	}
	// Prints out organized positions and distance combos
	//for (uint i=0; i < organizedPos.size(); i++){
	//	cout << organizedPos[i] << ": " << organizedDist[i] << endl;
	//}

	vector<int> varPos = interface01(sys, organizedPos);
	//
	//for (uint i=0; i<varPos.size(); i++){
	//	if (i == 21){
	//		cout << endl;
	//	}
	//	cout << varPos[i];
	//}
	//cout << endl;
	/******************************************************************************
	 *     === INITIALIZE POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	string polyLeu = generatePolyLeu("L", opt.sequenceLength);

	string polyLeuPS = generateMultiIDPolymerSequence(polyLeu, opt.thread, opt.Ids, varPos);
	PolymerSequence PL(polyLeuPS);

	/******************************************************************************
	 *                   === DECLARE SYSTEM FOR POLYLEU ===
	 ******************************************************************************/
	System sysL;
	CharmmSystemBuilder CSBL(sysL,opt.topFile,opt.parFile,opt.solvFile);
	CSBL.setBuildTerm("CHARMM_ELEC", false);
	CSBL.setBuildTerm("CHARMM_ANGL", false);
	CSBL.setBuildTerm("CHARMM_BOND", false);
	CSBL.setBuildTerm("CHARMM_DIHE", false);
	CSBL.setBuildTerm("CHARMM_IMPR", false);
	CSBL.setBuildTerm("CHARMM_U-BR", false);
	CSBL.setBuildTerm("CHARMM_IMM1REF", false);
	CSBL.setBuildTerm("CHARMM_IMM1", false);//having these two lines here has the energies back to normal, but can I call the energies later if they're not built?
	
	CSBL.setSolvent("MEMBRANE");
	CSBL.setIMM1Params(15, 10);

	if(!CSBL.buildSystem(PL)) {
		eout << "Unable to build system from " << polyLeu << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainAL = sysL.getChain("A");
	Chain & chainBL = sysL.getChain("B");
	
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainAL = chainAL.getAllAtomPointers();
	AtomPointerVector & apvChainBL = chainBL.getAllAtomPointers();
	vector<Position*>& positionsL = sysL.getPositions();

	SystemRotamerLoader sysRot(sysL, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hbL(sysL, opt.hBondFile);
	hbL.buildInteractions(30);//when this is here, the HB weight is correct
	
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* EsetL = sysL.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	EsetL->setAllTermsActive();
	EsetL->setTermActive("CHARMM_ELEC", false);
	EsetL->setTermActive("CHARMM_ANGL", false);
	EsetL->setTermActive("CHARMM_BOND", false);
	EsetL->setTermActive("CHARMM_DIHE", false);
	EsetL->setTermActive("CHARMM_IMPR", false);
	EsetL->setTermActive("CHARMM_U-BR", false);
	EsetL->setTermActive("CHARMM_IMM1REF", false);
	EsetL->setTermActive("CHARMM_IMM1", false);
	EsetL->setTermActive("SCWRL4_HBOND", false);
	
	// Set weights
	EsetL->setWeight("CHARMM_VDW", opt.weight_vdw);
	//EsetL->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	//EsetL->setWeight("SCWRL4_HBOND", 0);
	EsetL->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	EsetL->setWeight("CHARMM_IMM1", opt.weight_solv);
	fout << "Energy Weights: " << endl;
	fout << "VDW weight: " << EsetL->getWeight("CHARMM_VDW") << " IMM1REF weight: " << EsetL->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << EsetL->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
	fout << endl;

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	
	beb.setSystem(sysL);//had to have this before readParameters to get it to work! So it works now
	beb.readParameters(opt.baselineFile);
	//beb.printParameters();// Ensures they are being read properly
	
	beb.buildInteractions();	

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	//chainAL.renumberChain(opt.thread);
	//chainBL.renumberChain(opt.thread);
	sysL.wipeAllCoordinates();
	//sysL.assignCoordinates(glyAPV,false);//this seems to be fucking things up; maybe just assign to previously transformed coordinates in sys?
	sysL.assignCoordinates(sys.getAtomPointers(),false);
	sysL.buildAllAtoms();
	
	//writer.write(sysL.getAtomPointers(), true, false, true);//for some reason this only builds part of chainB? even though when printed out the atomPointer are correct?

	EsetL->saveEnergySubset("CHARMM_VDW");
	fout << "Total interactions calculated: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	//cout << EsetL->getSummary() << endl;//for troubleshooting; make sure to calculate energies before
	
	//loadInterfacialRotamers(sysL, sysRot, opt.SL, varPos);
	//loadInterfacialPeripheralRotamers(sysL, sysRot, opt.SL, varPos);
	loadRotamers(sysL, sysRot, opt.SL);
	//fout << sysL.getAllAtomPointers() << endl;

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms transL; 
	transL.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	transL.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	//transformation(apvChainAL, apvChainBL, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, transL);
	//moveZCenterOfCAMassToOrigin(sysL.getAllAtomPointers(), helicalAxis.getAllAtomPointers(), transL);
	
	/******************************************************************************
	 *                   === SET INTERFACIAL LINKED POSITIONS ===
	 ******************************************************************************/
	vector<vector<string>> pos2String;
	pos2String = positionToString(sys, varPos);
	//for (uint i=0; i<pos2String.size(); i++){
	//	for (uint j=0; j<pos2String[i].size(); j++){
	//		cout << i << "; " << j << ": " << pos2String[i][j] << endl;
	//	}
	//}
	
	sysL.setLinkedPositions(pos2String);//having this uncommented should make it go faster (homodimer)

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	//Random Number Generator
	RandomNumberGenerator RNG;
	//RNG.setSeed(opt.seed);//this isn't doing anything? 
	RNG.setTimeBasedSeed();//for my future runs, this is should be uncommented
	
	CSBL.updateNonBonded();
	sysL.buildAllAtoms();
	//writer3.write(sysL.getAtomPointers(), true, false, true);//these were here for troubleshooting
	fout << "Setting up SelfPairManager for optimization of sequence" << endl;
	SelfPairManager spm;
	spm.seed(RNG.getSeed());//I think in the future for my future runs, this should just be 0 (makes spm use a time based seed)
	spm.setSystem(&sysL);
	//spm.setVerbose(opt.verbose);
	spm.setVerbose(true);
	//spm.setRunDEE(opt.runDEESingles, opt.runDEEPairs);//TODO: figure out why this wasn't working
	spm.setRunDEE(true, true);
	spm.setOnTheFly(true);
	
	//Setup running SCMF or UnbiasedMC
	if (opt.runSCMF){
		fout << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	}
	else{
		fout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	fout << "Starting SelfPairManager Optimization..." << endl << endl;
	fout << "SPM Seed: " << RNG.getSeed() << endl;
	sysL.calcEnergy();
	fout << "Total interactions calced: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	fout << EsetL->getSummary() << endl;
	time(&spmStart);
	spm.runOptimizer();
	//spm.calculateEnergies();//TODO: run this whenever using just greedy
	//spm.runGreedyOptimizer(1);//for testing
	time(&spmEnd);
	
	fout << endl << "End SelfPairManager Optimization" << endl;
	spmTime = difftime (spmEnd, spmStart);
	fout << "SelfPairManager runOptimizer time: " << spmTime << " seconds" << endl;
	
	/******************************************************************************
	 *             === SETUP IMM1 ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************///decide if I need this or not (probably not?)
	if (opt.useIMM1 == true){
		fout << "Using IMM1 Energies for Monomer Comparison" << endl << endl;
		EsetL->setTermActive("CHARMM_IMM1REF", true);
		EsetL->setTermActive("CHARMM_IMM1", true);
	}
	else{
		fout << "Not using IMM1 Energies for Monomer Comparison" << endl << endl;
	}
	EsetL->eraseTerm("BASELINE");//gets rid of the baseline energy term
	
	sysL.setActiveRotamers(spm.getBestSCMFBiasedMCState());
	string seq1 = convertPolymerSeqToOneLetterSeq(chainAL);//make sure these are in the proper order: I was struggling to figure out what was wrong with my localMC code until I switched this
	fout << "System loading best SCMF Biased MC: " << seq1 << endl;
	fout << "SCMF Energy before local Monte Carlo: " << sysL.calcEnergy() << endl;

	double monomerEnergy;

	//The below compares the monomeric state and does a local repack, getting it to the best energetic state comparing dimer to monomer
	System finalPdb = localMC(sysL, helicalAxis, opt, seq1, monomerEnergy, fout, writer4);
	fout << "SCMF Energy after local Monte Carlo: " << finalPdb.calcEnergy() << endl;

	if (opt.runSCMF){
		fout << "Final SCMF Energy: " << finalPdb.calcEnergy()-monomerEnergy << endl;
		fout << finalPdb.getEnergySet()->getSummary() << endl;

	}
	else{
		sysL.setActiveRotamers(spm.getBestUnbiasedMCState());
		fout << "BestUnbiasedMC Energy: " << sysL.calcEnergy() << endl;
		fout << EsetL->getSummary() << endl;
	}

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	fout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

	writer.close();
	writer2.close();
	writer3.close();
	writer4.close();
	fout.close();
	mout.close();
	eout.close();
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
	Options data;

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
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("threadBool");

	opt.allowed.push_back("xShift");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("transform");
	//opt.allowed.push_back("thread");
	//opt.allowed.push_back("bbThread");

	//localMC repack variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	
	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaX");

	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("SL");
	
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	opt.allowed.push_back("numberOfStructuresToMCRepack");
	opt.allowed.push_back("energyCutOff");
	
	opt.allowed.push_back("uniprotName");
	opt.allowed.push_back("uniprotAccession");

	opt.allowed.push_back("inputMonomerE");
	opt.allowed.push_back("monoE_vdw");
	opt.allowed.push_back("monoE_hbond");
	opt.allowed.push_back("monoE_solv");
	opt.allowed.push_back("monoE_solvRef");

	opt.allowed.push_back("rmsdCutoff");
	opt.allowed.push_back("printAllCrds");
	opt.allowed.push_back("printAxes");
	opt.allowed.push_back("printTermEnergies");
	opt.allowed.push_back("deleteTerminalHbonds");

	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("baselineFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	
	opt.allowed.push_back("thread");
	opt.allowed.push_back("bbThread");
	
	opt.allowed.push_back("Ids");
	opt.allowed.push_back("numPositions");

	//Command Line Arguments
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("useIMM1");

	//SelfPairManager Arguments
	opt.allowed.push_back("runDEESingles");
	opt.allowed.push_back("runDEEPairs");
	opt.allowed.push_back("runSCMF");
	
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

	opt.useIMM1 = OP.getBool("useIMM1");
	if (OP.fail()) {
		opt.warningMessages += "useIMM1 not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.useIMM1 = false;
	}
	
	opt.printAllCrds = OP.getBool("printAllCrds");
	if (OP.fail()) {
		opt.warningMessages += "printAllCrds not specified using false\n";
		opt.warningFlag = true;
		opt.printAllCrds = false;
	}
	opt.printAxes = OP.getBool("printAxes");
	if (OP.fail()) {
		opt.warningMessages += "printAxes not specified using false\n";
		opt.warningFlag = true;
		opt.printAxes = false;
	}
	opt.printTermEnergies = OP.getBool("printTermEnergies");
	if (OP.fail()) {
		opt.printTermEnergies = true;
		opt.warningMessages += "printTermEnergies not specified using true\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}
	opt.rmsdCutoff = OP.getDouble("rmsdCutoff");
	if (OP.fail()) {
		opt.rmsdCutoff = 2.0;
		opt.warningMessages += "rmsdCutoff not specified using 2.0\n";
		opt.warningFlag = true;
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
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
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

	opt.threadStart = OP.getInt("threadStart");
	if (OP.fail()) {
		// excluded 2 residues in the beginning
		// opt.threadStart = 35 - ( (tmEnd-tmStart + 1) - 3 );  
		opt.threadStart = opt.startResNum + 37 - opt.endResNum;  
		opt.warningMessages += "threadStart not specified using " + MslTools::intToString(opt.threadStart) + "\n";
		opt.warningFlag = true;
	}
	opt.threadEnd = OP.getInt("threadEnd");
	if (OP.fail()) {
		opt.threadEnd = 33;
		opt.warningMessages += "threadEnd not specified using " + MslTools::intToString(opt.threadEnd) + "\n";
		opt.warningFlag = true;
	}
	opt.threadBool = OP.getBool("threadBool");
	if (OP.fail()) {
		opt.warningMessages += "threadBool not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.threadBool = false;
	}

	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.xShift = 6.7;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.zShift = 0;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.axialRotation = 0;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.crossingAngle = -40;
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
		}
		else{
			opt.bbThread = opt.backboneLength - opt.sequenceLength + 1;
		}
	}
	//opt.thread = opt.bbThread;
		
	//localMC repack variables
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
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
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

	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	}
	else{
		opt.SL = "SL"+opt.SL;
	}

	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 35\n";
		opt.backboneLength = 35;
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

	opt.baselineFile = OP.getString("baselineFile");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine baselineFile, no baseline forces used\n";
		//opt.errorFlag = true;
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
	
	opt.uniprotName = OP.getString("uniprotName");
	if (OP.fail()) {
		opt.uniprotName = "PROTEIN_UNK";
		opt.warningMessages += "uniprotName not specified using " + opt.uniprotName + "\n";
		opt.warningFlag = true;
	}

	opt.uniprotAccession = OP.getString("uniprotAccession");
	if (OP.fail()) {
		opt.uniprotAccession = "P00000";
		opt.warningMessages += "uniprotAccession not specified using " + opt.uniprotAccession + "\n";
		opt.warningFlag = true;
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

	opt.numPositions = OP.getInt("numPositions");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify number of positions\n";
		opt.errorFlag = true;
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
	
	opt.rerunConf = OP.getConfFile();

	return opt;
	return data;
}
