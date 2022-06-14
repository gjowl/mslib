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
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "seqDesign";
string programDescription = "This program aims to use changes I made to SPM to get it to save more than just one sequence state, and uses new baselines that previous design programs did not use (self + pair)";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "14 October 2020";
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
	int pairDist;

	int numberOfStructuresToMCRepack;
	double energyCutOff;
	
	// weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	double weight_seqEntropy;
	
	// input monomerEnergy
	bool inputMonomerE;
	int monoE_vdw;
	int monoE_hbond;
	int monoE_solv;
	int monoE_solvRef;

	// clustering options (optional)
	bool printAllCrds;
	bool printAxes;
	bool printTermEnergies;

	// alternate identities
	vector<string> Ids;
	vector<int> varPos;
	int numPositions;

	//
	bool runMCAfterSPM;
	vector<string> extraAAs;

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
	return "A" + ps + "\nB" + ps;
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

//TODO: make a separate alternateIds for this loop that only has AAs that I want to see in the extra monteCarlo
void addAAIdentities(CharmmSystemBuilder &_CSB, System &_sys, vector<string> _alternateIds, vector<int> _varPos, int _startResNum, int _seqLength) {
	int counter = 0;
	for(uint i = _startResNum; i<_startResNum+_seqLength; i++) {
		Position &posA = _sys.getChain("A").getPosition(i);
		Position &posB = _sys.getChain("B").getPosition(i);
		if (_varPos[counter] == 1){
			for (uint j=0; j<_alternateIds.size(); j++){
				_CSB.addIdentity(posA, _alternateIds[j]);
				_CSB.addIdentity(posB, _alternateIds[j]);
			}
		}
		else{
			continue;
		}
		counter++;
	}
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
			counter++;
		}
		else{
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (_varPos[counter] == 1){
				cout << _startResNum + counter << endl;
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

	_fout << "monoE_vdw " << _op.monoE_vdw << endl;
	_fout << "monoE_solv " << _op.monoE_solv << endl;
	_fout << "monoE_solvRef" << _op.monoE_solvRef << endl;
	_fout << "monoE_hbond" << _op.monoE_hbond << endl;

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
	_fout << "weight_seqEntropy " << _op.weight_seqEntropy << endl;

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

void loadRotamersCurrId(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
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
	for (uint k=4; k<_sys.getChain("A").positionSize()-4; k++) {
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
	}//TODO: fix this code; doesn't give me the interfacial positions?
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

vector<int> getVariablePos(vector<int> &_varPos){
	vector<int> pos;
	for (int k=0; k<_varPos.size()/2; k++){
		if (_varPos[k] == 1){
			pos.push_back(k);
		}
		else{
			continue;
		}
	}
	return pos;
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
	//CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	//CSBMono.setSolvent("MEMBRANE");
	//CSBMono.setIMM1Params(15, 10);
	
	CSBMono.setBuildNonBondedInteractions(false);
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
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", 1);
//	if (_useIMM1 == false){
//		monoEset->setTermActive("CHARMM_IMM1REF", false);
//		monoEset->setTermActive("CHARMM_IMM1", false);
//		monoEset->setWeight("CHARMM_IMM1REF", 0);
//		monoEset->setWeight("CHARMM_IMM1", 0);
//	}
//	else{
//		monoEset->setTermActive("CHARMM_IMM1REF", true);
//		monoEset->setTermActive("CHARMM_IMM1", true);
//		monoEset->setWeight("CHARMM_IMM1REF", 1);
//		monoEset->setWeight("CHARMM_IMM1", 1);
//	}

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSBMono.updateNonBonded(10,12,50);
	monoSys.buildAllAtoms();

	loadRotamers(monoSys, monoRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(true);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.calculateEnergies();

	repackSideChains(monoSpm, _opt.greedyCycles);

	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestState");
	_helicalAxis.saveAltCoor("BestAxis");
	cout << monoEset->getSummary() << endl;

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
	cout << monoSys.calcEnergy() << endl;

	//Select first and last resi of sequence
	//AtomSelection sel(monoSys.getAtomPointers());
	AtomSelection sel(chainA);
	string resi1 = "resi1, chain A and resi 1";
	string resi2 = "resi2, chain A and resi 21";//TODO: if this works, make it so that it's the final residue of chain using an option
	sel.select(resi1);
	sel.select(resi2);

	monoEset->deleteInteractionsWithinSelection("resi1");
	monoEset->deleteInteractionsWithinSelection("resi2");
	monoSys.calcEnergy();
	cout << monoSys.calcEnergy() << endl;
	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, _helicalAxis.getAtomPointers(), _trans);
	//AtomSelection sel(chainA);
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
	cout << monoSys.calcEnergy() << endl;
	
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

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0); 
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_greedyCycles);
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

			//_fout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
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
	string monoOutCrdFile  = _opt.pdbOutputDir + "/monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}
	
	string monoOutPdbFile  = _opt.pdbOutputDir + "/monomer.pdb";
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
		string meOutName  = _opt.pdbOutputDir + "/monomer.energy";
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

System localMC(System &_pdb, System &_helicalAxis, Options &_opt, string _sequence, double &_monomer, ofstream &_out, ofstream &_dout, ofstream &_sout, vector<double> &_xShifts, vector<double> &_crossingAngles, vector<double> &_axialRotations, vector<double> &_zShifts, int _seqNumber, string _dir, PDBWriter &_writer){
	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	//PDBWriter writer;
	//writer.open(_dir + "/sequence_" + to_string(_seqNumber) + ".pdb");

	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	cout << polymerSeq << endl;	

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
	if(!CSB.buildSystem(polymerSeq)) {
		cerr << "Unable to build system from " << polymerSeq << endl;
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
	sys.assignCoordinates(_pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	HydrogenBondBuilder hb(sys, _opt.hBondFile);
	hb.buildInteractions(30);//when this is here, the HB weight is correct
	
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
	
	/******************************************************************************
	 *             === SETUP ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	
	beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
	//beb.readParameters(_opt.baselineFile);
	
	//beb.buildInteractions();	

	// Objects used for transformations
	/******************************************************************************
	 *                     === SET UP TRANSFORMS OBJECT ===
	 ******************************************************************************/
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSB.updateNonBonded(10,12,50);
	sys.buildAllAtoms();

	loadRotamers(sys, sysRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.calculateEnergies();

	repackSideChains(spm, _opt.greedyCycles);

	sys.setActiveRotamers(spm.getMinStates()[0]);
	sys.saveAltCoor("savedBestState");
	_helicalAxis.saveAltCoor("BestAxis");
	
	//Eset->eraseTerm("BASELINE_PAIR");//gets rid of the baseline energy term
	cout << "Energy without Baseline: " << sys.calcEnergy() << endl;
	cout << Eset->getSummary() << endl;

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
	_dout << "Monomer calculations..." << endl;
	cout << "Monomer calculations..." << endl;
	monomerEnergy = computeMonomerEnergy(sys, trans, _opt, _helicalAxis, RNG, monomerEnergyByTerm, _out, _opt.greedyCycles, _opt.MCCycles, _opt.MCMaxRejects, false);
	_monomer = monomerEnergy;
	_out << "End monomer calculations! Monomer Energy: " << monomerEnergy << endl;
	_dout << "End monomer calculations! Monomer Energy: " << monomerEnergy << endl;
	cout << "End monomer calculations! Monomer Energy: " << _monomer << endl;

	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	double xShift = _opt.xShift;
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	_dout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	_dout << " xShift: " << xShift << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	cout << " xShift: " << xShift << endl;

	// Transform helices to initial starting position //Already set by setting coordinates earlier
	//transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), trans);

	cout << sys.calcEnergy() << endl;
	// Optimizatize Initial Starting Position
	sys.setActiveRotamers(spm.getMinStates()[0]);
	sys.calcEnergy();

	//double currentEnergy = spm1.getMinBound()[0];
	
	double currentEnergy = sys.calcEnergy();
	cout << currentEnergy << endl;
	//SelfPairManager spm;
	//spm.seed(RNG.getSeed());
	//spm.setSystem(&sys);
	//spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.calculateEnergies();
	repackSideChains(spm, _opt.greedyCycles);
	
	/******************************************************************************
	 *                     === X SHIFT REPACKS ===
	 ******************************************************************************/
	double bestEnergy = currentEnergy;
	double savedXShift = xShift;
	double previousEnergy = monomerEnergy;
	double deltaXShift = -0.1;
	double globalLowestE = monomerEnergy;
	//double xShiftEnd = MslTools::toDouble(parsedGeoInformation[4]);
	double xShiftEnd = xShift-0.5;//TODO: probably should change this to an option, or xShift-0.5 or something of the like

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
	
		//_out << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
		cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

		// If energy increase twice in a row, and it is above the global lowest energy, quit
		if (currentEnergy < globalLowestE) {
			globalLowestE = currentEnergy;
		}
		if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
			//_out << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
			break;	
		}
		else {
			previousEnergy = currentEnergy;
		}
	
	}
	_dout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl; 
	cout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl; 
	xShift = savedXShift;

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	_dout << "====================================" << endl;
	_dout << "Performing Local Monte Carlo Repacks" << endl;
	_dout << "====================================" << endl;
	_dout << endl;
	cout << "====================================" << endl;
	cout << "Performing Local Monte Carlo Repacks" << endl;
	cout << "====================================" << endl;
	cout << endl;

//	double bestEnergy = currentEnergy;
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	//bestEnergy = currentEnergy;
	_dout << "Best Energy before: " << bestEnergy << endl;
	_dout << Eset->getSummary() << endl;

	vector<unsigned int> MCOBest = spm.getMinStates()[0];
	
	if (_opt.MCCycles > 0) {
		//MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);
		MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);

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
			sys.calcEnergy();
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				//_out << "state rejected   energy: " << currentEnergy-monomerEnergy << endl;
				//_writer.write(sys.getAtomPointers(),true,false,true);
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

				cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_dout << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-monomerEnergy << endl;
	cout << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-monomerEnergy << endl;
	_dout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

	sys.applySavedCoor("savedBestState");
	_writer.write(sys.getAtomPointers(),true,false,true);
	double finalEnergy = sys.calcEnergy()-monomerEnergy;
	cout << "Final Energy: " << finalEnergy << endl;
	_sout << "Current Geometry " << endl;
	_sout << "xShift: " << xShift << endl;
	_sout << "crossingAngle: " << crossingAngle << endl;
	_sout << "axialRotation: " << axialRotation << endl;
	_sout << "zShift: " << zShift << endl << endl;
	_xShifts.push_back(xShift);
	_crossingAngles.push_back(crossingAngle);
	_axialRotations.push_back(axialRotation);
	_zShifts.push_back(zShift);
	_sout << "Energy Summary" << endl;
	_sout << "Monomer Energy: " << monomerEnergy << endl;
	_sout << "Dimer Energy: " << sys.calcEnergy() << endl;
	_sout << "Final Energy: " << finalEnergy << endl << endl;
	_dout << Eset->getSummary() << endl;
	sys.clearSavedCoor("savedBestState");
	_helicalAxis.clearSavedCoor("BestAxis");
	return sys;
}

void sameSequenceChecker(string &_sequence, vector<unsigned int> _stateVec, vector<string> &_seqs, vector<vector<unsigned int>> &_uniqueSeqs){
	if (_seqs.size() == 0){
		_seqs.push_back(_sequence);
		_uniqueSeqs.push_back(_stateVec);
		cout << "Unique Seq : " << _sequence << endl;
	}
	else {
		bool sameseq = false;
		for (int j=0; j<_seqs.size(); j++){
			int compare = _seqs[j].compare(_sequence);
			cout << "Seq " << j << ": " << _seqs[j] << "; " << _sequence << endl;
			cout << "Compare: " << compare << endl;
			if (compare == 0){
				sameseq = true;
				j == _seqs.size()-1;
			}
			cout << "Same Seq: " << sameseq << endl;
			if (j==_seqs.size()-1){
				if (sameseq == false){
					_seqs.push_back(_sequence);
					_uniqueSeqs.push_back(_stateVec);
					cout << "Unique Seq " << ": " << _sequence << endl;
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
		pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		//cout << tokens[0] << " " << tokens[1] << " " << tokens[2] << " = " << tokens[3] << endl;
	}
	
	pairReader.close();
	return pairEnergies;
}

void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	EnergySet* ESet = _sys.getEnergySet();
	
	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		//cout << "Chain: " << thisChain.toString() << endl;
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			//cout << "Position 1: " << p-positions.begin() << endl;
			for (uint j=0; j < (*p)->identitySize(); j++){
				Residue &res1 = (*p)->getIdentity(j);
				string baseId1 = res1.getResidueName();
				if (p-positions.begin() < 4){
					baseId1 = baseId1.append("-ACE");
				}
				if (p-positions.begin() > positions.size()-4){
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
							string baseId2 = res2.getResidueName();//TODO: check Charmm file top and par to see 
							if (p2-positions.begin() < 4){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-4){
								baseId2 = baseId2.append("-CT2");
							}
							//cout << "Identity " << k << ": " << baseId2 << endl;
							uint order = baseId1.compare(baseId2);
							if (order == 0){
								//cout << "Same AAs!" << endl;
								try{
									map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
									//cout << baseId1 << " worked" << endl;
									map<uint,double> AA2 = AA1.at(baseId1);
									//cout << baseId1 << " worked" << endl;
									double ener = AA2.at(d);
									//cout << ener << " worked" << endl;
									Atom *a = &res1.getAtom("CA");
									Atom *b = &res2.getAtom("CA");
									ESet->addInteraction(new BaselinePairInteraction(*a,*b,-1*ener));//I forgot that this needs to be the opposite sign to actually counteract the energies of vdW and hydrogen bonding; switched it here but should eventually just switch in my baseline parameter file
									//cout << baseId1 << ", " << baseId2 << ": " << ener << endl;//TODO: figure out how getResidueName() works; it's outputting the first and last LEU as just LEU, not LEU-ACE and LEU-CT2 respectively; I could hard code it here but then I'd also have to do so in the BaselinePairInteraction object
								}//TODO: do I need a distance?
								catch (const out_of_range& e){
									continue;		
								}
							}
							else{
								//cout << "Different AAs!" << endl;
								set<string> sortAAs;
								sortAAs.insert(baseId1);
								sortAAs.insert(baseId2);
								vector<string> sort;
								for (set<string>::iterator s = sortAAs.begin(); s != sortAAs.end(); ++s){
									string a = *s;
									sort.push_back(a);
								}
								try{
									map<string,map<uint,double>> AA1 = _pairMap.at(sort[0]);
									//cout << sort[0] << " worked" << endl;
									map<uint,double> AA2 = AA1.at(sort[1]);
									//cout << sort[1] << " worked" << endl;
									double ener = AA2.at(d);
									//cout << ener << "worked" << endl;
									Atom *a = &res1.getAtom("CA");
									Atom *b = &res2.getAtom("CA");
									ESet->addInteraction(new BaselinePairInteraction(*a,*b,-1*ener));//I forgot that this needs to be the opposite sign to actually counteract the energies of vdW and hydrogen bonding; switched it here but should eventually just switch in my baseline parameter file
									//cout << sort[0] << ", " << sort[1] << ": " << ener << endl;//TODO: figure out how getResidueName() works; it's outputting the first and last LEU as just LEU, not LEU-ACE and LEU-CT2 respectively; I could hard code it here but then I'd also have to do so in the BaselinePairInteraction object
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
}//TODO: maybe try a getEnergy for a combo to see what the interaction is and if it's actually being saved?

void buildSelfInteractions(System &_sys, map<string, double> &_selfMap){
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
				if (p-positions.begin() < 4 || p-positions.begin() > positions.size()-4){
					baseId = baseId.append("-OUT");
				}
				//cout << "Identity " << j << ": " << baseId << endl;
				try{
					//cout << baseId << " worked" << endl;
					double ener = _selfMap.at(baseId);
					//cout << ener << " worked" << endl;
					Atom *a = &res.getAtom("CA");
					ESet->addInteraction(new BaselineInteraction(*a,ener));//TODO:: make this more consistent in the code: right now the parameter file flips the sign but should probably just put it here
				}
				catch (const out_of_range& e){
					continue;		
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
					cout << baseId << " worked" << endl;
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

double baselineSelfEnergyOutput(System &_sys, map<string, double> &_selfMap){
	double totalEner = 0;

	for(uint i = 0; i < 1; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();

		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			Residue &res = (*p)->getCurrentIdentity();
			string baseId = res.getResidueName();
			if (p-positions.begin() < 1){
				baseId = baseId.append("-ACE");
			}
			if (p-positions.begin() > positions.size()-2){
				baseId = baseId.append("-CT2");
			}
			try{
				double ener = _selfMap.at(baseId);
				totalEner = totalEner + ener;
			}
			catch (const out_of_range& e){
				continue;		
			}
		}
	}
	return totalEner;
}

double baselinePairEnergyOutput(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	double totalEner = 0;

	for(uint i = 0; i < 1; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();

		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			Residue &res1 = (*p)->getCurrentIdentity();
			string baseId1 = res1.getResidueName();
			if (p-positions.begin() < 4){
				baseId1 = baseId1.append("-ACE");
			}
			if (p-positions.begin() > positions.size()-4){
				baseId1 = baseId1.append("-CT2");
			}
			for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
				uint d = p2-p;
				if (d <= 10){
					Residue &res2 = (*p2)->getCurrentIdentity();
					string baseId2 = res2.getResidueName();
					if (p2-positions.begin() < 4){
						baseId2 = baseId2.append("-ACE");
					}
					if (p2-positions.begin() > positions.size()-4){
						baseId2 = baseId2.append("-CT2");
					}
					uint order = baseId1.compare(baseId2);
					if (order == 0){
						//cout << "Same AAs!" << endl;
						try{
							map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
							//cout << baseId1 << "worked" << endl;
							map<uint,double> AA2 = AA1.at(baseId1);
							//cout << baseId1 << "worked" << endl;
							double ener = -1*AA2.at(d);
							//cout << ener << "worked" << endl;
							totalEner = totalEner + ener;
							//cout << baseId1 << ", " << baseId2 << ": " << ener << endl;//TODO: figure out how getResidueName() works; it's outputting the first and last LEU as just LEU, not LEU-ACE and LEU-CT2 respectively; I could hard code it here but then I'd also have to do so in the BaselinePairInteraction object
						}
						catch (const out_of_range& e){
							continue;		
						}
					}
					else{
						//cout << "Different AAs!" << endl;
						set<string> sortAAs;
						sortAAs.insert(baseId1);
						sortAAs.insert(baseId2);
						vector<string> sort;
						for (set<string>::iterator s = sortAAs.begin(); s != sortAAs.end(); ++s){
							string a = *s;
							sort.push_back(a);
						}
						try{
							map<string,map<uint,double>> AA1 = _pairMap.at(sort[0]);
							//cout << sort[0] << "worked" << endl;
							map<uint,double> AA2 = AA1.at(sort[1]);
							//cout << sort[1] << "worked" << endl;
							double ener = -1*AA2.at(d);
							//cout << ener << "worked" << endl;
							totalEner = totalEner + ener;
							//cout << sort[0] << ", " << sort[1] << ": " << ener << endl;//TODO: figure out how getResidueName() works; it's outputting the first and last LEU as just LEU, not LEU-ACE and LEU-CT2 respectively; I could hard code it here but then I'd also have to do so in the BaselinePairInteraction object
						}
						catch (const out_of_range& e){
							continue;		
						}
					}
				//	if (order == 0 || order < 0){
				//		try{
				//			map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
				//			map<uint,double> AA2 = AA1.at(baseId2);
				//			double ener = -1*AA2.at(d);
				//			totalEner = totalEner + ener;
				//		}
				//		catch (const out_of_range& e){
				//			continue;		
				//		}
				//	}
				//	else{
				//		try{
				//			map<string,map<uint,double>> AA1 = _pairMap.at(baseId2);
				//			map<uint,double> AA2 = AA1.at(baseId1);
				//			double ener = -1*AA2.at(d);
				//			totalEner = totalEner + ener;
				//		}	
				//		catch (const out_of_range& e){
				//			continue;		
				//		}
				//	}
				}
			}
		}
	}
	return totalEner;
}

vector<double> printSelfBaselineEnergy(string &_seq, map<string, double> _selfEner, double _outerSelfEner){
	vector<double> ener;
	for (uint i=0; i<_seq.length(); i++){
		string aa(1, _seq[i]);
		if (i < 4 || i > _seq.length()-4){//TODO: change this to an option for how many outer AAs there are
			ener.push_back(_outerSelfEner);
		}
		else{
			for (map<string, double>::const_iterator it = _selfEner.begin(); it != _selfEner.end(); ++it){
				if (it->first == aa){
					ener.push_back(it->second);
				}
			}
		}
	}
	return ener;
}

vector<double> printPairBaselineEnergy(string &_seq, Options _opt, map<string, vector<double>> _innerPairs, map<string, vector<double>> _outerPairs, ofstream &_out){
	vector<double> ener;
	vector<string> pairs1;
	vector<string> pairs2;
	vector<int> dists;
	for (uint i=0; i<_seq.length(); i++){
		for (uint j=i+1; j<_seq.length(); j++){
			int d = j-i;
			if (d <= _opt.pairDist){
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
	}
	int pos = 0;//uses the distance between pairs to increment position and choose either inner or outer pair energy maps to choose energy from
	for (uint i=0; i<pairs1.size(); i++){
		if (pairs1[i] != pairs2[i]){
			if (dists[i] == 1){
				pos++;
			}
			if (pos < 2 || pos > _seq.length()-2){//TODO: change these to an option for number of outer AAs
				for (map<string, vector<double>>::const_iterator it = _outerPairs.begin(); it != _outerPairs.end(); ++it){
					if (it->first == pairs1[i] || it->first == pairs2[i]){
						ener.push_back(it->second[dists[i]-1]);
						it == _outerPairs.end();
					}
				}
			}
			else{
				for (map<string, vector<double>>::const_iterator it = _innerPairs.begin(); it != _innerPairs.end(); ++it){
					if (it->first == pairs1[i] || it->first == pairs2[i]){
						ener.push_back(it->second[dists[i]-1]);
						it == _innerPairs.end();
					}
				}
			}
		}
		else{
			if (dists[i] == 1){
				pos++;
			}
			if (pos < 2 || pos > _seq.length()-2){//TODO: change these to an option for number of outer AAs
				for (map<string, vector<double>>::const_iterator it = _outerPairs.begin(); it != _outerPairs.end(); ++it){
					if (it->first == pairs1[i]){
						ener.push_back(it->second[dists[i]-1]);
						it == _outerPairs.end();
					}
				}
			}
			else{
				for (map<string, vector<double>>::const_iterator it = _innerPairs.begin(); it != _innerPairs.end(); ++it){
					if (it->first == pairs1[i]){
						ener.push_back(it->second[dists[i]-1]);
						it == _innerPairs.end();
					}
				}
			}
		}
	}

	double pos1 = 24;
	for (uint i=0; i<pairs1.size(); i++){
		double pos2 = pos1 + dists[i];
		if (dists[i] == 1){
			pos1++;
		}
		_out << pairs1[i] << pos1 << ", " << pos2 << ": " << ener[i] << endl;
	}
	_out << endl;
	pairs1.clear();
	pairs2.clear();
	dists.clear();
	return ener;
}

double calcEnergyOfSequence(System &_pdb, Options &_opt, string _sequence, int _seqNumber, double &_newSeqEnt, map<string, double> _selfMap, map<string, double> _seqEntMap, map<string,map<string,map<uint, double>>> _pairMap){
	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	//PDBWriter writer;
	//writer.open(_dir + "/sequence_" + to_string(_seqNumber) + ".pdb");

	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = generatePolymerSequenceFromSequence(_sequence, _opt.thread);
	cout << polymerSeq << endl;	

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
	if(!CSB.buildSystem(polymerSeq)) {
		cerr << "Unable to build system from " << polymerSeq << endl;
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
	sys.assignCoordinates(_pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	HydrogenBondBuilder hb(sys, _opt.hBondFile);
	hb.buildInteractions(30);//when this is here, the HB weight is correct
	
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
	Eset->setTermActive("BASELINE_ENTROPY_NORM", true);
	
	/******************************************************************************
	 *             === SETUP ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
	
	buildSelfInteractions(sys, _selfMap);
	buildPairInteractions(sys, _pairMap);
	
	BaselineSequenceEntropyNormalized bs(&sys);
	bs.setMap(_seqEntMap);
	Eset->addInteraction(new BaselineSequenceEntropyNormalized(bs));

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSB.updateNonBonded(10,12,50);
	sys.buildAllAtoms();
	
	sys.calcEnergy();

	loadRotamers(sys, sysRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.calculateEnergies();

	repackSideChains(spm, _opt.greedyCycles);

	sys.setActiveRotamers(spm.getMinStates()[0]);
	sys.saveAltCoor("savedBestState");
	
	sys.calcEnergy();
	cout << Eset->getSummary() << endl;
	_newSeqEnt = Eset->getTermEnergy("BASELINE_ENTROPY_NORM");
	return sys.calcEnergy();

}

map<double, string> sequenceMC(System &_sys, Options &_opt, vector<unsigned int> &_bestState, vector<vector<unsigned int>> &_deeAliveRotamers, vector<string> &_seqs, vector<vector<unsigned int>> &_uniqueSeqs, vector<int> &_varPosList, RandomNumberGenerator &_RNG, map<string, double> _selfMap, map<string, double> _seqEntMap, map<string,map<string,map<uint, double>>> _pairMap){
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MC.setRandomNumberGenerator(&_RNG);

	map<double, string> sequences;

	EnergySet *ESet = _sys.getEnergySet();
	// Start from most probable state
	vector<unsigned int> currentState;
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _sys.calcEnergy();
	cout << ESet->getSummary() << endl;
	_sys.saveAltCoor("bestState");

	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	vector<unsigned int> prevStateVec = _bestState;
	vector<unsigned int> stateVec = _bestState;
	
	currentState = _bestState;
	MC.setEner(bestEnergy);

	int cycleCounter = 0;
	Chain & chain = _sys.getChain("A");

	double enerToSave = 0;
	double prevSeqEnt = 0;
	double prevSEProb = 0;
	string prevSeq = convertPolymerSeqToOneLetterSeq(chain);

	while (!MC.getComplete()){
		//TODO: add in some form of changing to a random state here by states left from AAs...?
		try{
			double convertE = 0;
			
			if (_opt.verbose){
				cout << "Cycle #" << cycleCounter << "" << endl;
				cout << "Starting Seq: " << prevSeq << endl;
			}
			
			if (cycleCounter == 0){
				prevSeqEnt = ESet->getTermEnergy("BASELINE_ENTROPY_NORM");
				prevSEProb = exp(-prevSeqEnt/0.592);
			}

			cout << "Prev Seq Ent: " << prevSeqEnt << endl;
			cout << "Prev Prob: " << prevSEProb << endl;

			int rand = _RNG.getRandomInt(0, _varPosList.size()-1);
			int pos = _varPosList[rand];
			int rand1 = _RNG.getRandomInt(ids.size());
			string posId = _sys.getPosition(pos).getPositionId();
			string randId;
			randId = ids[rand1];
			
			string posA = "A," + MslTools::intToString(pos);
			string posB = "B," + MslTools::intToString(pos);
			//string prevId = _sys.getIdentity(pos);
			stringstream tmp;
			tmp << _sys.getIdentity(pos);
			string prevId = tmp.str();
			cout << _sys.getIdentity(pos) << endl;
			_sys.setActiveIdentity(posId, randId);
			cout << _sys.getIdentity(pos) << endl;

			string currSeq = convertPolymerSeqToOneLetterSeq(chain);
			cout << "Previous Seq: " << prevSeq << endl;
			cout << "Current Seq: " << currSeq << endl;

			if (currSeq.compare(prevSeq) == 0){
				while (currSeq.compare(prevSeq) == 0){
					rand1 = _RNG.getRandomInt(ids.size());
					randId = ids[rand1];
					_sys.setActiveIdentity(posId, randId);
					currSeq = convertPolymerSeqToOneLetterSeq(chain);
				}
			}
			double newSeqEnt = 0;
			double oligomerEnergy = calcEnergyOfSequence(_sys, _opt, currSeq, cycleCounter, newSeqEnt, _selfMap, _seqEntMap, _pairMap);

			//		Check to make sure the prevEnt thing is working properly

			cout << ESet->isTermActive("BASELINE_ENTROPY_NORM") << endl;
			//Resets the energy after every move to the calculated energy
			
			currSeq = convertPolymerSeqToOneLetterSeq(chain);
			if (_opt.verbose){
				cout << "Current Seq: " << currSeq << endl;
				//cout << "Current state: ";
				//for (int i=0; i<stateVec.size(); i++){
				//	cout << stateVec[i] << ",";
				//}
				//cout << endl;
			}
		
			//Get the proportionality by energy from old to new sequence

		
			//Calculate the Sequence Entropy Energy by proportionality of the old sequence to new sequence
			double newSEProb = exp(-newSeqEnt/0.592);
			double totSEProb = prevSEProb+newSEProb;
			double prevSeqProp = prevSEProb/totSEProb;
			double newSeqProp = newSEProb/totSEProb;
			double prevEner = -log(prevSeqProp)*0.592*100;
			double newEner = -log(newSeqProp)*0.592*100;
		
			if (_opt.verbose){
			cout << "Previous Seq: " << prevSeq << endl;
			cout << "Current Seq: " << currSeq << endl;
			cout << "Prev Seq Ent: " << prevSeqEnt << endl;
			cout << "Prev Prob: " << prevSEProb << endl;
			cout << "New Seq Ent: " << newSeqEnt << endl;
			cout << "New Prob: " << newSEProb << endl;

				cout << "Prev Seq Proportion: " << prevSeqProp << endl;
				cout << "New Seq Proportion: " << newSeqProp << endl;
				cout << "PrevEner = " << prevEner << endl;
				cout << "NewEner = " << newEner << endl;
				cout << "Diff = " << (prevEner-newEner) << endl;
			}
			
			double bestEnergyTotal = bestEnergy-prevSeqEnt+prevEner;//this needs to be run everytime...I think I need a new variable
			MC.setEner(bestEnergyTotal);
			double oligomerEnergyTotal = oligomerEnergy-newSeqEnt+newEner;
			
			cout << "Best Energy: " << bestEnergyTotal << endl;
			cout << "New Energy: " << oligomerEnergyTotal << endl;

			if (cycleCounter == 0){
				sequences[bestEnergy-prevSeqEnt] = prevSeq;
			}

			if (!MC.accept(oligomerEnergyTotal)){
				if (cycleCounter == 0){
					sameSequenceChecker(prevSeq, stateVec, _seqs, _uniqueSeqs);
				}
				//_sys.setActiveRotamers(prevStateVec);
				_sys.applySavedCoor("bestState");
				//stateVec = prevStateVec;
				_sys.setActiveRotamer(prevId, 1);

				if (_opt.verbose){
					cout << "State not accepted, E= " << oligomerEnergyTotal << "\n";
				}
			} else {
				sameSequenceChecker(currSeq, stateVec, _seqs, _uniqueSeqs);
				//prevStateVec = stateVec;
				_sys.saveAltCoor("bestState");
				sequences[oligomerEnergy-newSeqEnt] = currSeq;
				bestEnergy = oligomerEnergy;
				MC.setEner(oligomerEnergyTotal);
				prevSeqEnt = newSeqEnt;
				prevSEProb = newSEProb;
				prevSeq = currSeq;

				if (_opt.verbose) {
					cout << "State accepted, E= " << oligomerEnergyTotal << "\n";
				}
			}
			cycleCounter++;
		}
		catch (const out_of_range& e){
			_sys.setActiveRotamers(prevStateVec);
			cout << "Out of bounds exception caught!" << endl;
			continue;	
		}
	}
	cout << "Best State: ";
	for (int i=0; i<_bestState.size(); i++){
		cout << _bestState[i] << ",";
	}
	cout << endl;
	cout << "Energy: " << bestEnergy << endl;
	return sequences;
}

//hardcoded pair Baseline code that needs to be put into BaselineEnergyBuilder 
//The below are redacted and switched to readSingleParameter and readPairParameter
map<string, double> makeBaselineMap(string _file){
	map<string, double> m;
	vector<string> line;
	vector<string> AA;
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
	//separate the input of each line into two separate vectors
	for (uint i=0; i<line.size(); i++){
		if (i % 2 == 0){
			AA.push_back(line[i]);
		}
		else{
			ener.push_back(MslTools::toDouble(line[i]));
		}
	}
	for (uint i=0; i<AA.size(); i++){
		m.insert(make_pair(AA[i], ener[i]));
	}

	return m;
}


map<string, vector<double>> makeBaselineMapPair(string _file){
	map<string, vector<double>> m;
	vector<string> line;
	vector<string> AA;
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
	//separate the input of each line into two separate vectors
	for (uint i=0; i<line.size(); i++){
		if (i % 2 == 0){
			AA.push_back(line[i]);
		}
		else{
			ener.push_back(MslTools::toDouble(line[i]));
		}
	}
	
	string tmp;
	vector<double> es;
	for (uint i=0; i<AA.size(); i++){
		if  (i == 0){
			tmp = AA[i];
			es.push_back(ener[i]);
		}
		else{
			if (AA[i] == tmp){
				es.push_back(ener[i]);
			}
			else{
				m.insert({tmp, es});
				//reset the string and vector and push in the first value
				es.clear();
				tmp = AA[i];
				es.push_back(ener[i]);
			}
			if (i == AA.size()-1){
				m.insert({tmp, es});
				es.clear();
			}
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

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	ofstream fout;
	ofstream eout;
	ofstream sout;
	ofstream bmout;
	ofstream sdout;
	ofstream out3;
	
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string ddir = dir + "/design_" + opt.runNumber;
	cmd = "mkdir -p " + ddir;

	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string foutfile  = ddir + "/design.out";
	string eoutfile  = ddir + "/errors.out";
	string soutfile  = ddir + "/seq_summary.out";
	string bmoutfile = ddir + "/baselineMonomerComparison.out";
	string sdoutfile = ddir + "/seq_detailed.out";
	string outfile3 = ddir + "/pair_interactions.out";

	fout.open(foutfile.c_str());
	eout.open(eoutfile.c_str());
	sout.open(soutfile.c_str());
	bmout.open(bmoutfile.c_str());
	sdout.open(sdoutfile.c_str());
	out3.open(outfile3.c_str());

	eout << date << endl;
	sout << date << endl;
	bmout << date << endl;
	sdout << date << endl;

	//printOptions(opt, pout);
	
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
	 *                        === SETUP SUMMARY FILE ===
	 ******************************************************************************/
	sout << "Starting Geometry: " << endl;
	sout << "xShift: " << opt.xShift << endl;
	sout << "crossingAngle: " << opt.crossingAngle << endl;
	sout << "axialRotation: " << opt.axialRotation << endl;
	sout << "zShift: " << opt.zShift << endl << endl;
	bmout << "Starting Geometry: " << endl;
	bmout << "xShift: " << opt.xShift << endl;
	bmout << "crossingAngle: " << opt.crossingAngle << endl;
	bmout << "axialRotation: " << opt.axialRotation << endl;
	bmout << "zShift: " << opt.zShift << endl << endl;
	sdout << "Starting Geometry: " << endl;
	sdout << "xShift: " << opt.xShift << endl;
	sdout << "crossingAngle: " << opt.crossingAngle << endl;
	sdout << "axialRotation: " << opt.axialRotation << endl;
	sdout << "zShift: " << opt.zShift << endl << endl;
	fout << "Starting Geometry: " << endl;
	fout << "xShift: " << opt.xShift << endl;
	fout << "crossingAngle: " << opt.crossingAngle << endl;
	fout << "axialRotation: " << opt.axialRotation << endl;
	fout << "zShift: " << opt.zShift << endl << endl;
	
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
	 ******************************************************************************///for troubleshooting
	PDBWriter writer;
	writer.open(ddir + "/design_" + opt.runNumber + ".pdb");

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

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(opt.infile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
	sys.wipeAllCoordinates();
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	//writer.write(sys.getAtomPointers(), true, false, true);//for troubleshooting
	string seqA = convertPolymerSeqToOneLetterSeq(chainA);
	string seqB = convertPolymerSeqToOneLetterSeq(chainB);

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

	// Transformation to command line arguments for zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	sys.buildAllAtoms();

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

	vector<int> varPos = interface01(sys, organizedPos);

	/******************************************************************************
	 *     === INITIALIZE POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	string polyLeu = generatePolyLeu("L", opt.backboneLength);
	string polyLeuPS = generateMultiIDPolymerSequence(polyLeu, opt.thread, opt.Ids, varPos);
	cout << polyLeuPS << endl;
	PolymerSequence PL(polyLeuPS);

	/******************************************************************************
	 *                   === BEGINNING OUTPUT FOR OUT FILE ===
	 ******************************************************************************/
	fout << "Sequence: " << polyLeu << endl;
	fout << "DiffPos:  ";
	for (uint i=0; i<varPos.size(); i++){
		if (i == 21){
			i = varPos.size();
			fout << endl;
		}
		fout << varPos[i];
	}
	fout << endl;

	fout << "RotLevel: " << opt.SL << endl;
	fout << "Alternate Ids: ";
	for (uint i=0; i<opt.Ids.size(); i++){
		if (i == opt.Ids.size()-1){
			fout << endl;
		}
		fout << opt.Ids[i] << " ";
	}
	fout << endl;
	
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

	CSBL.setBuildNonBondedInteractions(false);
	if(!CSBL.buildSystem(PL)) {
		eout << "Unable to build system from " << polyLeuPS << endl;
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

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sysL.assignCoordinates(sys.getAtomPointers(),false);
	sysL.buildAllAtoms();

	CSBL.updateNonBonded(10,12,50);
	
	SystemRotamerLoader sysRot(sysL, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hbL(sysL, opt.hBondFile);
	hbL.buildInteractions(50);//when this is here, the HB weight is correct
	
	/******************************************************************************
	 *                 === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	for (uint i=0; i<sysL.atomSize(); i++){
		Atom at = sysL.getAtom(i);
		if (!at.hasCoor()){
			eout << "Atom " << i << " was not assigned coordinates; program termination";
			break;
		}
		else{
			continue;
		}
	}
	cout << "All atoms have coordinates" << endl;

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
	EsetL->setTermActive("CHARMM_VDW", true);
	EsetL->setTermActive("SCWRL4_HBOND", true);
	
	// Set weights
	EsetL->setWeight("CHARMM_VDW", 1);
	EsetL->setWeight("SCWRL4_HBOND", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sysL,opt);

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	beb.setSystem(sysL);//had to have this before readParameters to get it to work! So it works now

	//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
	map<string, double> selfMap = readSingleParameters("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/2020_10_07_meanSelf_par.txt");
	map<string, double> seqEntMap = readSingleParameters("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/seqEntropies.txt");
	map<string,map<string,map<uint,double>>> pairMap = readPairParameters("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/2020_10_07_meanPair_par.txt");
	map<string,map<string,map<uint,double>>> pairMapOut = readPairParameters("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/2020_10_07_meanPair_par.txt");
	buildSelfInteractions(sysL, selfMap);
	buildPairInteractions(sysL, pairMap);
	//buildSequenceEntropy(sysL, seqEntMap, opt.weight_seqEntropy);
	//TODO: I got the entropy map, and built it...do I have to do anything else before I run it to see if it works?

	//initialize baselineAAComposition energy map (based on Rosetta baselineAAComposition to prevent unlikely sequences by adding energy (ex. if more than 2PHE in sequence, add 100 energy score for each additional PHE)

	for (map<string,double>::const_iterator it = seqEntMap.begin(); it != seqEntMap.end(); it++){
		cout << it->first << ": " << it->second << endl;
	}
	//BaselinePermutation bp(&sysL);
	//bp.setMap(seqEntMap);
	//for (map<string,double>::const_iterator it = seqEntMap.begin(); it != seqEntMap.end(); it++){
	//	cout << it->first << ": " << it->second << endl;
	//}
//	EsetL->addInteraction(new BaselineSequenceEntropyNormalized(bs));
	//BaselineAAComposition bac(&sysL);
	//bac.readPenaltyFile("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/AACompositionPenalties.out");
	//map<string,double> penaltyMap = bac.getMap();
	//for (map<string,double>::const_iterator it = penaltyMap.begin(); it != penaltyMap.end(); it++){
	//	cout << it->first << ": " << it->second << endl;
	//}
	//EsetL->addInteraction(new BaselineAAComposition(bac));

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms transL; 
	transL.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	transL.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	vector<vector<string>> linkedPos = positionToString(sysL, varPos);
	sysL.setLinkedPositions(linkedPos);
	
	loadRotamers(sysL, sysRot, opt.SL);

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	//Random Number Generator
	RandomNumberGenerator RNG;
	RNG.setTimeBasedSeed();
	
	CSBL.updateNonBonded();
	sysL.buildAllAtoms();
	fout << "Setting up SelfPairManager for optimization of sequence" << endl;
	cout << EsetL->getSummary() << endl;

	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sysL);
	spm.setVerbose(true);
	spm.setRunDEE(opt.runDEESingles, opt.runDEEPairs);
	//spm.setRunDEE(false, false);
	spm.setOnTheFly(true);
	
	//Setup running SCMF or UnbiasedMC
	if (opt.runSCMF == true){
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
	cout << "Total interactions calced: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	cout << EsetL->getSummary() << endl;
	time(&spmStart);
	spm.runOptimizer();
	time(&spmEnd);
	
	fout << endl << "End SelfPairManager Optimization" << endl;
	cout << endl << "End SelfPairManager Optimization" << endl;
	spmTime = difftime (spmEnd, spmStart);
	//TODO: would be nice to output the run times for DEE and SCMF
	fout << "SelfPairManager runOptimizer time: " << spmTime << " seconds" << endl;
	//fout << "States after DEE: " << spm.getDEEAliveRotamers().size() << endl;
	fout << "States after SCMF: " << spm.getBestSCMFBiasedMCStates().size() << endl;
	cout << "SelfPairManager runOptimizer time: " << spmTime << " seconds" << endl;
	//cout << "States after DEE: " << spm.getDEEAliveRotamers().size() << endl;
	cout << "States after SCMF: " << spm.getBestSCMFBiasedMCStates().size() << endl;

	//TODO: On 1-22-21: run was started for DEE with all AAs around 2 and is still currently running at 10
	//TODO: Depending on how runs look from 12-18-20, may need to input a monte carlo here to go through sequences after getting the best state and just checking energy as a way to save 100 sequences to try at this geometry
	//As of 1-8-20 I'm still waiting on runs from the SL90 rotamer level to see if any better results; otherwise, could be nice to add this in and then try these instead
	//TODO: add in a way to run and save the DEE for each of these runs in case I want to run them again?

	/******************************************************************************
	 *           === METHODS FOR DETERMINING ALTERNATE SEQUENCES ===
	 ******************************************************************************/
	vector<string> seqs;
	vector<vector<unsigned int>> uniqueSeqs;
	map<double, string> unorderedSeqMap;
	map<string, double> orderedSeqMap;
	
	if (opt.runMCAfterSPM){
		cout << "Run MC on sequence after SPM: " << opt.runMCAfterSPM << endl;
		/******************************************************************************
		 *             === ADD IDENTITIES FOR AAS WITH MANY ROTAMERS ===
		 ******************************************************************************/
		//addAAIdentities(CSBL, sysL, opt.extraAAs, varPos, opt.thread, opt.sequenceLength);	
		//TODO: I thought of a potentially quicker way to run these: for the AAs with high rotamers, I can add those AAs in AFTER the SelfPairManager Optimization. These can then be potential solutions that will be found during the monteCarlo above
		//TODO: make this a function	
		/******************************************************************************
		 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
		 ******************************************************************************/
		//DEE alive rotamers to be randomly oriented in
		vector<vector<unsigned int>> deeAliveRotamers = spm.getDEEAliveRotamers();
		vector<int> varPosList = getVariablePos(varPos);
		vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
		
	//buildSelfInteractions(sysL, selfMap);
	//buildPairInteractions(sysL, pairMap);
		BaselineSequenceEntropyNormalized bs(&sysL);
		bs.setMap(seqEntMap);
		EsetL->addInteraction(new BaselineSequenceEntropyNormalized(bs));
		sysL.calcEnergy();
		unorderedSeqMap = sequenceMC(sysL, opt, bestState, deeAliveRotamers, seqs, uniqueSeqs, varPosList, RNG, selfMap, seqEntMap, pairMap);
		
	} else{
		/******************************************************************************
		 *              === RETRIEVE ALL UNIQUE SEQUENCES FROM DEE SCMF  ===
		 ******************************************************************************/
		if (spm.getBestSCMFBiasedMCStates().size() > 0){
			cout << "Number of good SCMF states: " << spm.getBestSCMFBiasedMCStates().size() << endl;
			for (int i=0; i<spm.getBestSCMFBiasedMCStates().size(); i++){
				if (i < 100){
					vector <unsigned int> stateVec = spm.getBestSCMFBiasedMCStates()[i];
					sysL.setActiveRotamers(stateVec);//add in sequence checker here
					sysL.buildAtoms();
					string seq = convertPolymerSeqToOneLetterSeq(chainAL);
					cout << "Sequence " << i << ": " << seq << endl;
					sysL.calcEnergy();
					//cout << EsetL->getSummary() << endl;
					sameSequenceChecker(seq, stateVec, seqs, uniqueSeqs);
				}
				else{
					i = spm.getBestSCMFBiasedMCStates().size();
				}
			}
		}
		else{
			cout << "SCMF didn't run" << endl;
		}
	}
	cout << "Total Sequences: " << seqs.size() << endl;

	//sort(unorderedSeqMap.begin(), unorderedSeqMap.end());
	for (map<double,string>::const_iterator it = unorderedSeqMap.begin(); it != unorderedSeqMap.end(); it++){
		cout << it->first << ": " << it->second << endl;
	}

	/******************************************************************************
	 *               === LOCAL REPACKS ON EACH UNIQUE SEQUENCE ===
	 ******************************************************************************/
	vector<double> monomerEnergies;
	vector<double> dimerEnergies;
	vector<double> baselineEnergies;
	vector<double> baselineEnergiesAdd;
	vector<double> baselineEnergiesAdd1;
	vector<double> baselineEnergiesAdd2;
	vector<double> xShifts;
	vector<double> crossingAngles;
	vector<double> axialRotations;
	vector<double> zShifts;
//TODO: unedited
	for (int i=0; i<unorderedSeqMap.size(); i++){
		if (i < 20){ // Currently limiting number of sequences to 20
			if (seqs[i].compare("LLLLLLLLLLLLLLLLLLLLL") == 0){
				i++;
			}
			sysL.setActiveRotamers(uniqueSeqs[i]);
			sysL.buildAtoms();
			sysL.calcEnergy();
			cout << EsetL->getSummary() << endl;
			cout << chainAL.toString() << endl;
			cout << "Seq #" << i << ": " << seqs[i] << endl;
			
			//for troubleshooting
			writer.write(sysL.getAtomPointers(), true, false, true);
			vector<double> selfBaselines = printSelfBaselineEnergy(seqs[i], self1, selfOuterLeu);
			vector<double> selfBaselines1 = printSelfBaselineEnergy(seqs[i], self1, 0.476915);
			out3 << seqs[i] << endl;
			vector<double> pairBaselines = printPairBaselineEnergy(seqs[i], opt, pair1, pair2, out3);
			out3 << seqs[i] << "1" << endl;
			vector<double> pairBaselines1 = printPairBaselineEnergy(seqs[i], opt, pair1, pair1, out3);
			map<string, vector<Interaction*>> *es =  EsetL->getEnergyTerms();
			vector<Interaction*> intSelf = es->at("BASELINE");
			vector<Interaction*> ints = es->at("BASELINE_PAIR");
			
			/******************************************************************************
			 *                      === SETUP OUTPUT FILES ===
			 ******************************************************************************/
			fout << "System loading SCMF Sequence #" << i << ": " << seqs[i] << endl;
			fout << EsetL->getSummary() << endl;
			cout << EsetL->getSummary() << endl;
			sdout << "SCMF Energy before local Monte Carlo: " << sysL.calcEnergy() << endl;
			sdout << EsetL->getSummary() << endl;
			
			//The below compares the monomeric state and does a local repack, getting it to the best energetic state by comparing dimer to monomer energies
			double monomerEnergy;
			System finalPdb = localMC(sysL, helicalAxis, opt, seqs[i], monomerEnergy, fout, sdout, sout, xShifts, crossingAngles, axialRotations, zShifts,  i, ddir, writer);
			sout << "END SEQUENCE " << i << endl;
			monomerEnergies.push_back(monomerEnergy);
			fout << "Energy after local Monte Carlo Sequence #" << i << ": " << finalPdb.calcEnergy() << endl;
		
			//Calculate the Self and Pair energies for troubleshooting
			double self = baselineSelfEnergyOutput(sysL, selfMap);
			double pair = baselinePairEnergyOutput(sysL, pairMap);//TODO: this is only partially complete: make sure to add in the outer baselines to the pairBaselineMap file
			baselineEnergies.push_back(2*(self + pair));

			//Calculate Dimer energy for output
			//TODO: maybe change this to make it easier for me to extract into a new file (if need be; may already be good enough)
			double dimerEnergy = finalPdb.calcEnergy()-monomerEnergy;
			dimerEnergies.push_back(dimerEnergy);
			fout << "Dimer Energy: " << dimerEnergy << " for Seq #" << i << ": " << seqs[i] << endl << endl;
			fout << finalPdb.getEnergySet()->getSummary() << endl;
		}
		else{
			i = uniqueSeqs.size();
		}
	}


	/******************************************************************************
	 *         === WRITE OUT BASELINE TO MONOMER COMPARISON BY SEQUENCE ===
	 ******************************************************************************/
	int seqNumber = MslTools::toInt(opt.runNumber);
	if (seqNumber == 0){
 		bmout << "Dimer: Monomer: Baseline: startXShift: startCrossingAngle: startAxialRotation: startZShift: xShift: crossingAngle: axialRotation: zShift: Sequence:" << endl;//TODO: change this output to some sort of sstring with the %s like in BaselinePairInteractions
		seqNumber++;//TODO: make a more permanent solution for this; temporary just to get runs through on 2020_10_16
	}
	else{
		for (uint i=0; i<monomerEnergies.size(); i++){
			bmout << "Final Sequence Info: " << dimerEnergies[i] << ": " << monomerEnergies[i] << ": " << baselineEnergies[i] << ": " << opt.xShift << ": " << opt.crossingAngle << ": " << opt.axialRotation << ": " << opt.zShift << ": " << xShifts[i] << ": " << crossingAngles[i] << ": " << axialRotations[i] << ": " << zShifts[i] << ": " << baselineEnergies[i] << ": " << seqs[i] << endl;
 		}
	}
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	fout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

	writer.close();
	fout.close();
	eout.close();
	//pout.close();
	sout.close();
	bmout.close();
	sdout.close();
	//out1.close();
	//out2.close();
	out3.close();
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

	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	opt.allowed.push_back("weight_seqEntropy");
	
	opt.allowed.push_back("SL");
	
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	opt.allowed.push_back("numberOfStructuresToMCRepack");
	opt.allowed.push_back("energyCutOff");
	
	opt.allowed.push_back("inputMonomerE");
	opt.allowed.push_back("monoE_vdw");
	opt.allowed.push_back("monoE_hbond");
	opt.allowed.push_back("monoE_solv");
	opt.allowed.push_back("monoE_solvRef");

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
	opt.allowed.push_back("varPos");
	opt.allowed.push_back("numPositions");

	opt.allowed.push_back("runMCAfterSPM");
	opt.allowed.push_back("extraAAs");

	//Command Line Arguments
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("useIMM1");

	//SelfPairManager Arguments
	opt.allowed.push_back("runDEESingles");
	opt.allowed.push_back("runDEEPairs");
	opt.allowed.push_back("runSCMF");
	opt.allowed.push_back("pairDist");
	
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

	opt.pairDist = OP.getInt("pairDist");
	if (OP.fail()){
		opt.warningMessages += "pairDist not specified, using 15\n";
		opt.warningFlag = true;
		opt.pairDist = 15;
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
		opt.warningMessages += "greedyCycles not specified using 10\n";
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
	opt.weight_seqEntropy = OP.getDouble("weight_seqEntropy");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_seqEntropy not specified, default 0.0\n";
		opt.weight_solv = 0.0;
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
	opt.varPos = OP.getIntVector("varPos");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA positions, make sure they are space separated\n";
		opt.errorFlag = true;
	}

	opt.numPositions = OP.getInt("numPositions");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify number of positions\n";
		opt.errorFlag = true;
	}
	
	opt.extraAAs = OP.getStringVector("extraAAs");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify additional AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.runMCAfterSPM = OP.getBool("runMCAfterSPM");
	if (OP.fail()){
		opt.errorMessages += "runMCAfterSPM not specified, defaulting to false";
		opt.warningFlag = true;
		opt.runMCAfterSPM = false;
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
}
