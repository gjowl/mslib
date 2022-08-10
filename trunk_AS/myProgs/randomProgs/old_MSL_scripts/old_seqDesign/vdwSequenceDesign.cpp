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

time_t startTime, endTime;
double diffTime;

/******************************************************************************************************************************************************************************/

struct Options{
	string sequence;
	string backboneAA;
	int backboneLength;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	// TODO: think I have to change this to a vector for specifically the interfacial residues 
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

	int start;
	int end;
	double ener;
	vector<int> ivalues;

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
};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
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

void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable) {

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
		if (_pos.size() < _numPositions){
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
	
	cout << date << endl;

	time(&startTime);
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
	*    printProteinOutFile
	*
	**********************************************************************************/
	ofstream pout;
	opt.pdbOutputDir = opt.pdbOutputDir + "/" + date + "/tester";
	string dir = opt.pdbOutputDir;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		//exit(0);
	}
	string poutName = dir + "/tester.out";
	
	pout.open(poutName.c_str());
	
	printOptions(opt, pout);

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
		cout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === GENERATE POLYGLY ===
	 ******************************************************************************/
	string polyGly = generatePolymerSequence("G", 30, opt.thread);
	
	PolymerSequence PS(polyGly);

	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	PDBWriter writer;
	PDBWriter writer2;
	PDBWriter writer3;
	PDBWriter writer4;
	PDBWriter writer5;
	PDBWriter writer6;
	
	writer.open(dir + "/polyGly.pdb");
	writer2.open(dir + "/helixAlignToCoord.pdb");
	writer3.open(dir + "/polyLeu.pdb");
	writer4.open(dir + "/MCOFinalStates.pdb");
	writer5.open(dir + "/SCMFstate.pdb");
	writer6.open(dir + "/BestUnbiasedMCState.pdb");
	
	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	CSB.setBuildNoTerms();

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << polyGly << endl;
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

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();

	writer.write(sys.getAtomPointers(), true, false, true);

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
	
	writer2.write(sys.getAtomPointers(), true, false, true);
	
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);//TODO: changed this line to reflect the config options, eventually should probably change back to datafile (or extracted geometry from pdb)
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
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

	/******************************************************************************
	 *     === INITIALIZE POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	string polyLeu = generatePolyLeu("L", 30);

	string polyLeuPS = generateMultiIDPolymerSequence(polyLeu, opt.thread, opt.Ids, varPos);

	PolymerSequence PL(polyLeuPS);
	
	/******************************************************************************
	 *                   === DECLARE SYSTEM FOR POLYLEU ===
	 ******************************************************************************/
	System sysL;
	CharmmSystemBuilder CSBL(sysL,opt.topFile,opt.parFile);
	CSBL.setBuildTerm("CHARMM_ELEC", false);
	CSBL.setBuildTerm("CHARMM_ANGL", false);
	CSBL.setBuildTerm("CHARMM_BOND", false);
	CSBL.setBuildTerm("CHARMM_DIHE", false);
	CSBL.setBuildTerm("CHARMM_IMPR", false);
	CSBL.setBuildTerm("CHARMM_U-BR", false);

	if(!CSBL.buildSystem(PL)) {
		cerr << "Unable to build system from " << polyLeu << endl;
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
	//EsetL->setWeight("CHARMM_IMM1REF", 0);
	//EsetL->setWeight("CHARMM_IMM1", 0);
	cout << "VDW weight: " << EsetL->getWeight("CHARMM_VDW") << " IMM1REF weight: " << EsetL->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << EsetL->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
	pout << endl;

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
	sysL.assignCoordinates(glyAPV,false);
	sysL.buildAllAtoms();
	
	EsetL->saveEnergySubset("CHARMM_VDW");
	cout << "System Energy: " << sysL.calcEnergy() << endl;//this calculates the current energy of all of the interactions
	cout << "Total interactions calculated: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	cout << EsetL->getSummary() << endl;
	
	loadInterfacialRotamers(sysL, sysRot, opt.SL, varPos);
	loadInterfacialPeripheralRotamers(sysL, sysRot, opt.SL, varPos);
	//loadRotamers(sysL, sysRot, opt.SL);
	cout << sysL.getAllAtomPointers() << endl;


	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms transL; 
	transL.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	transL.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	transformation(apvChainAL, apvChainBL, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, transL);
	moveZCenterOfCAMassToOrigin(sysL.getAllAtomPointers(), helicalAxis.getAllAtomPointers(), transL);
	
	/******************************************************************************
	 *                   === SET INTERFACIAL LINKED POSITIONS ===
	 ******************************************************************************/
	vector<vector<string>> pos2String;
	pos2String = positionToString(sys, varPos);
	
	sysL.setLinkedPositions(pos2String);//having this uncommented should make it go faster (homodimer)

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	//Random Number Generator
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);
	
	CSB.updateNonBonded();
	sysL.buildAllAtoms();
	writer3.write(sysL.getAtomPointers(), true, false, true);//not sure why I added these but don't think they're necessary, probably just troubleshooting
	
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sysL);
	spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.setRunDEE(true, true);
	spm.updateWeights();
	cout << "System Energy: " << sysL.calcEnergy() << endl;
	
	cout << "Total interactions calced: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	cout << "VDW Energy: " << EsetL->getTermEnergy("CHARMM_VDW") << endl;
	cout << "HBOND Energy: " << EsetL->getTermEnergy("SCWRL4_HBOND") << endl;
	cout << "IMM1REF Energy: " << EsetL->getTermEnergy("CHARMM_IMM1REF") << endl;
	cout << "IMM1 Energy: " << EsetL->getTermEnergy("CHARMM_IMM1") << endl;
	cout << EsetL->getSummary() << endl;
	
	writer3.write(sysL.getAtomPointers(), true, false, true);
	//pout << "PR005: " << endl;
	//pout << sysL.getAllAtomPointers() << endl;

	writer5.write(sysL.getAtomPointers(), true, false, true);
	spm.setOnTheFly(true);
	//spm.calculateEnergies();//already done in repackSideChains and runOptimizer
	//spm.runGreedyOptimizer(1);
	
	spm.setRunSCMF(false);
	spm.setRunUnbiasedMC(false);
	spm.runOptimizer();
	spm.runGreedyOptimizer(10);

	/******************************************************************************
	 *               === PRINT OUT SPM RUN INFO AND WRITE PDBS ===
	 ******************************************************************************/
	vector<unsigned int> MCOFinal;
	cout << "MCOFinal Energies: " << endl;
	for (uint i=0; i<spm.getMinStates().size(); i++){
		MCOFinal = spm.getMinStates()[i];
		sysL.setActiveRotamers(MCOFinal);
		pout << "PR00" << i+6 << ": " << endl;
		pout << sysL.getAllAtomPointers() << endl;
		//sysL.buildAllAtoms();
		cout << sysL.calcEnergy() << endl;
		cout << EsetL->getSummary() << endl;
		writer4.write(sysL.getAtomPointers(), true, false, true);
	}
	cout << "Min Energy: " << spm.getMinBound()[0] << endl;
	sysL.setActiveRotamers(spm.getSCMFstate());
	cout << EsetL->getSummary() << endl;
	cout << "System Energy: " << sysL.calcEnergy() << endl;
	cout << "VDW Energy: " << EsetL->getTermEnergy("CHARMM_VDW") << endl;
	cout << "HBOND Energy: " << EsetL->getTermEnergy("SCWRL4_HBOND") << endl;
	cout << "IMM1REF Energy: " << EsetL->getTermEnergy("CHARMM_IMM1REF") << endl;
	cout << "IMM1 Energy: " << EsetL->getTermEnergy("CHARMM_IMM1") << endl;
	cout << EsetL->getSummary() << endl;
	pout << "System Energy: " << sysL.calcEnergy() << endl;
	pout << "VDW Energy: " << EsetL->getTermEnergy("CHARMM_VDW") << endl;
	pout << "HBOND Energy: " << EsetL->getTermEnergy("SCWRL4_HBOND") << endl;
	pout << "IMM1REF Energy: " << EsetL->getTermEnergy("CHARMM_IMM1REF") << endl;
	pout << "IMM1 Energy: " << EsetL->getTermEnergy("CHARMM_IMM1") << endl;
	pout << EsetL->getSummary() << endl;
	//cout << spm.getSummary(spm.getSCMFstate()) << endl;
	writer5.write(sysL.getAtomPointers(), true, false, true);
	//cout << "SPM Energy: " << spm.getStateEnergy(spm.getSCMFstate()) << endl;//this portion doesn't work for vdw for some reason
	cout << endl;
	sysL.setActiveRotamers(spm.getBestUnbiasedMCState());
	cout << EsetL->getSummary() << endl;
	cout << "System Energy: " << sysL.calcEnergy() << endl;
	cout << "VDW Energy: " << EsetL->getTermEnergy("CHARMM_VDW") << endl;
	cout << "HBOND Energy: " << EsetL->getTermEnergy("SCWRL4_HBOND") << endl;
	cout << "IMM1REF Energy: " << EsetL->getTermEnergy("CHARMM_IMM1REF") << endl;
	cout << "IMM1 Energy: " << EsetL->getTermEnergy("CHARMM_IMM1") << endl;
	cout << EsetL->getSummary() << endl;
	//cout << spm.getSummary(spm.getBestUnbiasedMCState()) << endl;
	writer6.write(sysL.getAtomPointers(), true, false, true);
	cout << "System Energy: " << sysL.calcEnergy() << endl;

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

	writer.close();
	writer2.close();
	writer3.close();
	writer4.close();
	writer5.close();
	writer6.close();
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

	opt.allowed.push_back("mcCycles");
	opt.allowed.push_back("mcMaxRejects");
	opt.allowed.push_back("mcStartTemp");
	opt.allowed.push_back("mcEndTemp");
	opt.allowed.push_back("mcCurve");

	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaX");

	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("SL");
	
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
	opt.bbThread = opt.thread + opt.sequence.length() - opt.backboneLength + 1;
	opt.thread = opt.bbThread;
		
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
	
	opt.ener = OP.getInt("ener");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "Energy not specified, default to -15\n";
		opt.end = -15;
	}

	opt.ivalues = OP.getMultiInt("ivalues");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "i+x values not specified, default to 1, 3, and 4\n";
		opt.ivalues.push_back(1);
		opt.ivalues.push_back(3);
		opt.ivalues.push_back(4);
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
	
	opt.rerunConf = OP.getConfFile();

	return opt;
	return data;
}
