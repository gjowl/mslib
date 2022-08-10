//#include <iostream>
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
#include "BaselineOuterPairInteraction.h"
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"
#include "SasaCalculator.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "seqDesign";
string programDescription = "This is the most updated version of seqDesign: aims design sequences from geometric data from the PDB optimizing specifically for vdW energies";
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
	bool useGeoFromPDBData;

	// tm
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int sequenceStart;

	bool deleteTerminalHbonds;
	
	string SL; //number of rotamers
	string SLInterface; //number of rotamers

	// transformation
	double xShift;
	double zShift;
	bool leftHanded;
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

	// energy terms to output
	vector<string> monomerEnergyTerms;
	vector<string> monomerIMM1EnergyTerms;
	vector<string> dimerEnergyTerms;
	vector<string> calcEnergyOfSequenceTerms;
	vector<string> sequenceMCEnergyTerms;
	vector<string> enerLandscapeTerms;
	vector<string> energyTermsToOutput;

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

	_fout << "thread " << _op.thread << endl;

	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;
	_fout << "weight_solv " << _op.weight_solv << endl;
	_fout << "weight_seqEntropy " << _op.weight_seqEntropy << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), "SL80.00")) {//lower rotamer level because I did baselines at this level
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
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

//Interface identified using occluded surface area
void identifyInterface(System &_sys, Options &_opt, PolymerSequence _PS, vector<int> &_pos, string _outputDir, string _axis){
	// Declare system
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildNoTerms();
	
	if(!CSB.buildSystem(_PS)) {
		cout << "Unable to build system from " << _PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	sys.assignCoordinates(_sys.getAtomPointers(), false);

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

	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, 6.4, trans);//one of the shortest distances given from my pdb searc
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	sys.buildAllAtoms();

	// Setup output files
	ofstream out;
	string outfile  = _outputDir + "/Interface_SASA.out";
	out.open(outfile.c_str());
	out << "Using SASA values, interfacial residues are chosen for any SASA with less than 85 percent of highest SASA" << endl << endl;
	
	PDBWriter writer;
	writer.open(_outputDir + "/polyAla_for_SASA.pdb");
	writer.write(sys.getAtomPointers(),true,false,true);
	writer.close();

	// Setup SASA calculator
	SasaCalculator sa;
	sa.addAtoms(sys.getAtomPointers());

	sa.calcSasa();
	double totSasa = 0;
	cout << "Total SASA: " << sa.getTotalSasa() << endl;
	
	double worstSasa = 0;
	vector<double> sasas;
	vector<int> pos;
	vector<string> posIds;

	// Get all SASA values for each residue on one chain
	out << "SASA values per residue:" << endl;
	for (uint k=1; k<sys.getChain("A").positionSize()-1; k++) {
		Position &pos1 = sys.getChain("A").getPosition(k);
		
		Atom &c1 = pos1.getAtom("CA");
		string posId = c1.getPositionId();

		double sasa = sa.getResidueSasa(posId);
		out << posId << ": " << sasa << endl;
		cout << posId << ": " << sasa << endl;
		totSasa = totSasa+sasa;
		//cout << posId << ": " << sasa << endl;
		// TODO: add a cutoff for sasa (right now I'm thinking 85% of the highest sasa value because it works well for the polyAla of the first sequence I looked at) 
		sasas.push_back(sasa);
		pos.push_back(k);
		posIds.push_back(posId);
		if (k == 4){
			worstSasa = sasa;
		}
		else if (k < sys.getChain("A").positionSize()-4){
			if (sasa > worstSasa){
				worstSasa = sasa;
			}
		}
	}
	cout << "Added sasa: " << totSasa << endl;
	cout << "Avg SASA: " << totSasa/19 << endl;

	// Get the interface by values less than 0.85*worst SASA (worked well for the first trial I did to not get back facing residues)
	double cutoff = totSasa/19;// now using the average of internal positions
	out << "Residues chosen with less than " << cutoff << " SASA:" << endl;
	cout << "Residues chosen with less than " << cutoff << " SASA:" << endl;
	for (uint j=0;j<sasas.size(); j++){
		if (sasas[j] < cutoff){
			_pos.push_back(pos[j]);
			out << posIds[j] << endl;
			cout << posIds[j] << endl;
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
		Atom &c1 = pos1.getAtom("CB");
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

//boolean added to make multiple interface01 (I need to keep Leu in the first 4AAs because of baseline, but some are considered interfacial and should have more rotamers)
vector<int> interface01(System &_sys, vector<int> &_pos, bool _onlyInterface=true){
	vector<int> varPos;
	for (uint k=0; k<_sys.positionSize(); k++){
		varPos.push_back(0);
	}
	for (uint j=0; j<_pos.size(); j++){
		varPos[_pos[j]] = 1;
		varPos[_pos[j]+_sys.getChain("A").positionSize()] = 1;
	}
	if (_onlyInterface == false){
		for (uint k=0; k<4; k++){
			varPos[k] = 0;
			varPos[k+_sys.getChain("A").positionSize()] = 0;
		}
		for (uint k=_sys.getChain("A").positionSize()-4; k<_sys.getChain("A").positionSize(); k++){
			varPos[k] = 0;
			varPos[k+_sys.getChain("A").positionSize()] = 0;
		}
	}
	return varPos;
}

//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

vector<vector<string>> positionToString(System &_sys, vector<int> &_varPos){
	vector<vector<string>> stringPositions;
	if (_varPos.size()%2){
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
	}
	else{
		for (uint k=0; k<_varPos.size()/2+1; k++){
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

	}
	return stringPositions;
}

vector<int> getVariablePos(vector<int> &_varPos){
	vector<int> pos;
	if (_varPos.size()%2){
		for (int k=0; k<_varPos.size()/2; k++){
			if (_varPos[k] == 1){
				pos.push_back(k);
			}
			else{
				continue;
			}
		}
	}
	else{
		for (int k=0; k<_varPos.size()/2+1; k++){
			if (_varPos[k] == 1){
				pos.push_back(k);
			}
			else{
				continue;
			}
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

vector<double> calcBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
	}
	sel.clearStoredSelections();
	return ener;
}

vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+1);//TODO: fix this so it uses the thread
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+1);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			}
			else{
				j = _seqLength;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, int _seqStart){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+_seqStart);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
	}
	sel.clearStoredSelections();
	return ener;
}

vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength, int _seqStart){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+_seqStart);//TODO: fix this so it uses the thread
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+_seqStart);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			}
			else{
				j = _seqLength;
			}
		}
	}
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

double computeMonomerEnergyNoMoves(System & _sys, Options& _opt, SelfPairManager &_monoSpm, map<string,double> &_seqMCEnergies, string &_seq, RandomNumberGenerator &_RNG, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {
	
	map<string, double> monomerEnergies;

	for (uint i=0; i<_opt.monomerEnergyTerms.size(); i++){
		monomerEnergies[_opt.monomerEnergyTerms[i]] = 0;
	}

	
	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);	
	
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);//03_26_21: found my monomer issue (maybe?); In the baselineSelfPairComparison code, this had to be 1 to line up with the baselines
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

	//CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

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

	CSBMono.updateNonBonded(10,12,50);

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hBondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	PDBWriter writer;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	writer.open(dir + "/monomer.pdb");
	for (uint i=0; i<monoSys.atomSize(); i++){
		Atom at = monoSys.getAtom(i);
		if (!at.hasCoor()){
			cout << "Atom " << i << " was not assigned coordinates; program termination";
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
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", 1);

	//TODO: I think I would like to use this just in case the sequence actually needs to be compared, but this would mess up the baseline comparison. Figure out what can be done instead	
	//BaselineAAComposition bac(&monoSys);
	//bac.readPenaltyFile("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/AACompositionPenalties.out");
	//monoEset->addInteraction(new BaselineAAComposition(bac));
	
	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);
	
	/*****************************************************************************
	 *                 === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//loadRotamers(monoSys, monoRot, _opt.SL);
	loadMonomerRotamers(monoSys, monoRot);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	_monoSpm.setSystem(&monoSys);

	repackSideChains(_monoSpm, 10);

	/*****************************************************************************
	 *            === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ===
	 ******************************************************************************/
	//_monomerMinState.push_back(_monoSpm.getMinStates()[0]);
	monoSys.setActiveRotamers(_monoSpm.getMinStates()[0]);
	monoSys.calcEnergy();
	cout << monoEset->getSummary() << endl;
	
	string monoOutPdbFile  = _opt.pdbOutputDir + "/monomer.pdb";
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);

	cout << "Monomer time: " << diffTimeMono << endl;
	double finalEnergy = 2.0 * _monoSpm.getMinBound()[0]; // double the energy for 2 helices
	double vdw = monoEset->getTermEnergy("CHARMM_VDW")*2.0;
	double hbond = monoEset->getTermEnergy("SCWRL4_HBOND")*2.0;
	vector<double> self = calcBaselineEnergies(monoSys, 21);
	vector<double> pair = calcPairBaselineEnergies(monoSys, 21);
	double monomerSelf = sumEnergyVector(self)*2.0;
	double monomerPair = sumEnergyVector(pair)*2.0;

	monomerEnergies["VDWMonomer"] = vdw;
	monomerEnergies["HbondMonomer"] = hbond;
	monomerEnergies["MonomerSelfBaseline"] = monomerSelf;
	monomerEnergies["MonomerPairBaseline"] = monomerPair;
	monomerEnergies["Monomer"] = finalEnergy;

	map<string, double>::iterator itr;
	for (itr = monomerEnergies.begin(); itr != monomerEnergies.end(); ++itr){
		cout << itr->first << ": " << itr->second << endl;
		_seqMCEnergies.at(itr->first); //checks if the string exists in the map
		_seqMCEnergies[itr->first] = itr->second;
	}
	return finalEnergy;
}

double computeMonomerEnergy(System & _sys, Options& _opt, Transforms & _trans, map<string,double> &_monomerIMM1Energies, string _seq, RandomNumberGenerator &_RNG, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {

	_monomerIMM1Energies.clear();

	for (uint i=0; i<_opt.monomerIMM1EnergyTerms.size(); i++){
		_monomerIMM1Energies[_opt.monomerIMM1EnergyTerms[i]] = 0;
	}

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);	
	
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

	CSBMono.updateNonBonded(10,12,50);

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hBondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	PDBWriter writer;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	writer.open(dir + "/monomer.pdb");
	for (uint i=0; i<monoSys.atomSize(); i++){
		Atom at = monoSys.getAtom(i);
		if (!at.hasCoor()){
			cout << "Atom " << i << " was not assigned coordinates; program termination";
			break;
		}
		else{
			continue;
		}
	}
	cout << "All atoms have coordinates" << endl;

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
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSBMono.updateNonBonded(10,12,50);
	monoSys.buildAllAtoms();

	//loadRotamers(monoSys, monoRot, _opt.SL);
	loadMonomerRotamers(monoSys, monoRot);
	
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
	helicalAxis.saveAltCoor("BestAxis");
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

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), _trans);
	AtomSelection sel(chainA);
	//string resi1 = "resi1, chain A and resi 1";
	//string resi2 = "resi2, chain A and resi 21";//TODO: if this works, make it so that it's the final residue of chain using an option
	//sel.select(resi1);
	//sel.select(resi2);

	//monoEset->deleteInteractionsWithinSelection("resi1");
	//monoEset->deleteInteractionsWithinSelection("resi2");
	monoSys.calcEnergy();
	cout <<" Monomer before membrane: " << monoSys.calcEnergy() << endl;
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
	cout << "Monomer into membrane: " << monoSys.calcEnergy() << endl;
	
	// Repack side chains
	monoSpm.setOnTheFly(1);
	monoSpm.calculateEnergies();
        monoSpm.runGreedyOptimizer(_greedyCycles);

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
				helicalAxis.saveAltCoor("BestMonomerAxis");
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
			helicalAxis.saveAltCoor("BestMonomerAxis");
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
	cout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");
	monoSys.calcEnergy();
	cout << monoEset->getSummary();
	cout << endl;

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

	map<string,double> monomerEnergyByTerm;
	// Store monomer energy by term
	if(_opt.printTermEnergies) {
		monoSys.calcEnergy();
		monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

		ofstream meOut;
		string meOutName  = _opt.pdbOutputDir + "/monomer_IMM1.energy";
		meOut.open(meOutName.c_str());
		if(!meOut.is_open()) {
			cerr << "Unable to open " << meOutName << endl;
			exit(0);
		}
		for(map<string,double>::iterator it = monomerEnergyByTerm.begin(); it != monomerEnergyByTerm.end(); it++) {
			meOut << it->first << " " << it->second << endl;
		}
		meOut.close();
	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	double imm1 = (monoEset->getTermEnergy("CHARMM_IMM1")+monoEset->getTermEnergy("CHARMM_IMM1REF"))*2.0;
	double vdw = monoEset->getTermEnergy("CHARMM_VDW")*2.0;
	double hbond = monoEset->getTermEnergy("SCWRL4_HBOND")*2.0;
	//vector<double> self = calcBaselineEnergies(monoSys, 21, _opt.thread);
	//vector<double> pair = calcPairBaselineEnergies(monoSys, 21, _opt.thread);
	//double monomerSelf = sumEnergyVector(self)*2.0;
	//double monomerPair = sumEnergyVector(pair)*2.0;

	_monomerIMM1Energies["VDWMonomerw/IMM1"] = vdw;
	_monomerIMM1Energies["HbondMonomerw/IMM1"] = hbond;
	_monomerIMM1Energies["IMM1Monomer"] = imm1;
	//_monomerIMM1Energies["MonomerSelfBaselinew/IMM1"] = monomerSelf;
	//_monomerIMM1Energies["MonomerPairBaselinew/IMM1"] = monomerPair;
	_monomerIMM1Energies["Monomerw/IMM1"] = finalEnergy;

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

double computeDimerIMM1Energy(System & _sys, Transforms & _trans, Options& _opt, int seqNumber, System &_helicalAxis, RandomNumberGenerator &_RNG, ofstream &_fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

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
	CSB.buildSystemFromPDB(inputChain.getAtomPointers());

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(_sys.getAtomPointers(),false);
	sys.buildAllAtoms();

	SystemRotamerLoader rot(sys, _opt.rotLibFile);
	rot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hBondFile);
	hb.buildInteractions(30);

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
	Eset->setTermActive("CHARMM_VDW", false);
	Eset->setTermActive("SCWRL4_HBOND", false);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	Eset->setWeight("CHARMM_VDW", 0);
	Eset->setWeight("SCWRL4_HBOND", 0);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	Eset->setWeight("CHARMM_IMM1", 1);

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSB.updateNonBonded(10,12,50);

	// Optimize Initial Starting Position (using Baseline to get back to original result)

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chains = sys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	_trans.translate(axisB, moveAxisBOneAngstrom);
	
	sys.calcEnergy();
	cout << sys.calcEnergy() << endl;

	// move center of mass to origin
	moveZCenterOfCAMassToOrigin(chains, _helicalAxis.getAtomPointers(), _trans);
	AtomSelection sel(chains);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	_trans.translate(chains, interDistVect);

	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	_trans.translate(chains, move5Down);
	double bestZ = -5.0;

	sys.calcEnergy();
	cout << sys.calcEnergy() << endl;
	
	// Repack side chains
	double currentEnergy = sys.calcEnergy();
	double bestEnergy = currentEnergy;
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chains, zUnitVector);

		//double currentZ = -5.0 + ((i+1)*1.0); 
		currentEnergy = sys.calcEnergy();
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			sys.saveAltCoor("savedBestDimer");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}
	cout << "Best Z shift = " << bestZ << endl;

	sys.applySavedCoor("SavedBestDimer");

	double finalEnergy = sys.calcEnergy(); // double the energy for 2 helices
	return finalEnergy;
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

System localMC(System &_pdb, System &_helicalAxis, Options &_opt, int _seqNumber, map<string,double> &_localMCEnergies, string _sequence, vector<int> _varPosRotamers, double _monomerEnergy, vector<double> &_xShifts, vector<double> &_crossingAngles, vector<double> &_axialRotations, vector<double> &_zShifts, ofstream &_dout, ofstream &_sout, string _dir, PDBWriter &_writer){
	_localMCEnergies.clear();

	for (uint i=0; i<_opt.dimerEnergyTerms.size();i++){
		_localMCEnergies[_opt.dimerEnergyTerms[i]] = 0;
	}
	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************///for troubleshooting
	PDBWriter writer;
	writer.open(_dir + "/design_" + MslTools::intToString(_seqNumber) + ".pdb");

	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	//string polymerSeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	string polymerSeq = generatePolymerSequenceFromSequence(_sequence, _opt.thread);
	cout << polymerSeq << endl;	

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
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);
	
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
	Eset->setTermActive("CHARMM_IMM1", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	
	/******************************************************************************
	 *             === SETUP ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);
	
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
	loadInterfacialRotamers(sys, sysRot, _opt.SLInterface, _varPosRotamers);

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
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	double xShift = _opt.xShift;
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	_dout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	_dout << " xShift: " << xShift << endl << endl;
	_sout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	_sout << " xShift: " << xShift << endl << endl;
	cout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	cout << " xShift: " << xShift << endl;

	cout << sys.calcEnergy() << endl;
	// Optimizatize Initial Starting Position
	sys.setActiveRotamers(spm.getMinStates()[0]);
	sys.calcEnergy();

	//double currentEnergy = spm1.getMinBound()[0];
	
	double currentEnergy = sys.calcEnergy();
	cout << currentEnergy << endl;
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
	double previousEnergy = _monomerEnergy;
	double deltaXShift = -0.1;
	double globalLowestE = _monomerEnergy;
	//double xShiftEnd = MslTools::toDouble(parsedGeoInformation[4]);
	//double xShiftEnd = xShift-0.5;//TODO: probably should change this to an option, or xShift-0.5 or something of the like
	double xShiftEnd = 7.5;//TODO: probably should change this to an option, or xShift-0.5 or something of the like
	
	if (crossingAngle < 0){
		xShiftEnd = 6.4;
	}

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
		cout << "xShift: " << xShift << " energy: " << currentEnergy-_monomerEnergy << endl;

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
	_dout << "Best Energy at x shift: " << bestEnergy-_monomerEnergy << " at " << savedXShift << endl; 
	cout << "Best Energy at x shift: " << bestEnergy-_monomerEnergy << " at " << savedXShift << endl; 
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

	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
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
				_dout << "state rejected   energy: " << currentEnergy-_monomerEnergy << endl;
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

				cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-_monomerEnergy << endl;
				_dout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-_monomerEnergy << endl;
			}
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_dout << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-_monomerEnergy << endl;
	cout << "Best MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << bestEnergy-_monomerEnergy << endl;
	_dout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	_sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

	sys.applySavedCoor("savedBestState");
	_writer.write(sys.getAtomPointers(),true,false,true);
	writer.write(sys.getAtomPointers(),true,false,true);
	double finalEnergy = sys.calcEnergy()-_monomerEnergy;
	cout << "Final Energy: " << finalEnergy << endl;
	_dout << "Current Geometry " << endl;
	_dout << "xShift: " << xShift << endl;
	_dout << "crossingAngle: " << crossingAngle << endl;
	_dout << "axialRotation: " << axialRotation << endl;
	_dout << "zShift: " << zShift << endl << endl;
	_sout << "Current Geometry " << endl;
	_sout << "xShift: " << xShift << endl;
	_sout << "crossingAngle: " << crossingAngle << endl;
	_sout << "axialRotation: " << axialRotation << endl;
	_sout << "zShift: " << zShift << endl << endl;
	double vdw = Eset->getTermEnergy("CHARMM_VDW");
	double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
	double imm1 = Eset->getTermEnergy("CHARMM_IMM1REF")+Eset->getTermEnergy("CHARMM_IMM1");
	//_vdw.push_back(Eset->getTermEnergy("CHARMM_VDW"));
	//_hbond.push_back(Eset->getTermEnergy("SCWRL4_HBOND"));
	_localMCEnergies["VDWDimer"] = vdw;
	_localMCEnergies["HbondDimer"] = hbond;
	_localMCEnergies["IMM1Dimer"] = imm1;
	_localMCEnergies["Dimer"] = sys.calcEnergy();
	_xShifts.push_back(xShift);
	_crossingAngles.push_back(crossingAngle);
	_axialRotations.push_back(axialRotation);
	_zShifts.push_back(zShift);
	_sout << "Energy Summary Below" << endl;
	_sout << "Monomer Energy: " << _monomerEnergy << endl;
	_sout << "Dimer Energy: " << sys.calcEnergy() << endl;
	_sout << "Final Energy: " << finalEnergy << endl << endl;
	_sout << Eset->getSummary() << endl << endl;
	_dout << "Energy Summary Below" << endl;
	_dout << "Monomer Energy: " << _monomerEnergy << endl;
	_dout << "Dimer Energy: " << sys.calcEnergy() << endl;
	_dout << "Final Energy: " << finalEnergy << endl << endl;
	_dout << Eset->getSummary() << endl;
	sys.clearSavedCoor("savedBestState");
	_helicalAxis.clearSavedCoor("BestAxis");
	writer.close();
	return sys;
}

bool sameSequenceChecker(string &_sequence, vector<string> &_seqs){
	bool sameseq = false;
	if (_seqs.size() == 0){
		_seqs.push_back(_sequence);
	}
	else {
		for (int j=0; j<_seqs.size(); j++){
			cout << _seqs[j] << " : " << _sequence << ":" << sameseq << endl;
			if (_seqs[j] == _sequence){
				sameseq = true;
				j = _seqs.size()-1;
			}
			else if (j==_seqs.size()-1){
			cout << _seqs[j] << " : " << _sequence << ":" << sameseq << endl;
				if (sameseq == false){
					_seqs.push_back(_sequence);
					cout << "Unique Seq " << ": " << _sequence << endl;
					j = _seqs.size();
				}
			}
		}
	}
	return sameseq;
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
		}
		else{
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
			pairEnergies[MslTools::toUpper(tokens[1])][MslTools::toUpper(tokens[0])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		}
	}
	
	pairReader.close();
	return pairEnergies;
}

//Based on the above function: Reads a Gaussian Kde file to be used to search for a density estimate to choose a spot in the geometric space for the run geometry
// This one is specifically for rotation and z shift
void getGaussianKdeValues(string _file, Options &_opt, double _kdeValue, double &_xShift, double &_crossingAngle, double &_axialRotation, double &_zShift, ofstream &_out){
	// Setup kde file reader
	Reader kdeReader(_file);
	kdeReader.open();
	map<double,vector<double>> kdeEnergies;

	if(!(kdeReader.is_open())){
		cerr << "WARNING: Unable to open " << _file << endl;
		exit(0);
	}

	vector<string> lines = kdeReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 7){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: KDE(double) xShift(double) Angle(double) Rot1(double) Rot2(double) Z1(double) Z2(double)";
			continue;
		}
		vector<double> geometries;
		geometries.push_back(MslTools::toDouble(tokens[1]));
		geometries.push_back(MslTools::toDouble(tokens[2]));
		geometries.push_back(MslTools::toDouble(tokens[3]));
		geometries.push_back(MslTools::toDouble(tokens[4]));
		geometries.push_back(MslTools::toDouble(tokens[5]));
		geometries.push_back(MslTools::toDouble(tokens[6]));
		kdeEnergies[MslTools::toDouble(tokens[0])] = geometries;
	}

	map<double,vector<double>>::iterator low;
	low = kdeEnergies.lower_bound(_kdeValue);
	_out << "Found kde value: " << low->first << endl;

	RandomNumberGenerator ht;
	ht.setSeed(_opt.seed);

	int r = ht.getRandomInt(2,3);
	int z = ht.getRandomInt(4,5);

	_out << "xShift: " << kdeEnergies.at(low->first)[0] << endl;
	_out << "crossingAngle: " << kdeEnergies.at(low->first)[1] << endl;
	_out << "axialRotation: " << r-1 << "	" << kdeEnergies.at(low->first)[r] << endl; // subtract to get which column it's a part of
	_out << "zShift: " << z-3 << "	" << kdeEnergies.at(low->first)[z] << endl;
	_xShift = kdeEnergies.at(low->first)[0];
	_crossingAngle = kdeEnergies.at(low->first)[1];
	_axialRotation = kdeEnergies.at(low->first)[r];
	_zShift = kdeEnergies.at(low->first)[z];

	kdeReader.close();
}

// This one is for angle and distance
void getGaussianKdeValues(string _file, Options &_opt, double _kdeValue, double &_xShift, double &_crossingAngle, ofstream &_out){
	Reader kdeReader(_file);
	kdeReader.open();
	map<double,vector<double>> kdeEnergies;

	if(!(kdeReader.is_open())){
		cerr << "WARNING: Unable to open " << _file << endl;
		exit(0);
	}

	vector<string> lines = kdeReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 3){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: KDE(double) xShift(double) crossingAngle(double)";
			continue;
		}
		vector<double> geometries;
		geometries.push_back(MslTools::toDouble(tokens[1]));
		geometries.push_back(MslTools::toDouble(tokens[2]));
		kdeEnergies[MslTools::toDouble(tokens[0])] = geometries;
	}

	map<double,vector<double>>::iterator low;
	low = kdeEnergies.lower_bound(_kdeValue);

	_out << low->first << "	";
	_xShift = kdeEnergies.at(low->first)[0];
	_crossingAngle = kdeEnergies.at(low->first)[1];


	kdeReader.close();
}

void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	EnergySet* ESet = _sys.getEnergySet();
	ofstream pout;
	string dir = "exports/home/gloiseau/mslib/trunk_AS";
	
	for(uint i = 0; i < 1; i++) {
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
							string baseId2 = res2.getResidueName();//TODO: check Charmm file top and par to see 
							if (p2-positions.begin() < 4){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-5){//On 03_18_2021 I found this error; position.size() is weird, so need to use 5 instead of 4
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
}//TODO: maybe try a getEnergy for a combo to see what the interaction is and if it's actually being saved?

void buildSelfInteractions(System &_sys, map<string, double> &_selfMap){
	EnergySet* ESet = _sys.getEnergySet();

	for(uint i = 0; i < 1; i++) {
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

vector<double> printSelfBaselineEnergy(string &_seq, map<string, double> _selfEner, double _outerSelfEner){
	vector<double> ener;
	for (uint i=0; i<_seq.length(); i++){
		stringstream ss;
		string aa(1, _seq[i]);
		ss << aa;
		string key = ss.str();
		if (i < 4 || i > _seq.length()-5){//TODO: change this to an option for how many outer AAs there are
			ener.push_back(_outerSelfEner);
			cout << key << ": " << _outerSelfEner << endl;
		}
		else{
			for (map<string, double>::const_iterator it = _selfEner.begin(); it != _selfEner.end(); ++it){
				//cout << it->first << endl;
				if (it->first == key){
					ener.push_back(it->second);
					cout << key << ": " << it->second << endl;
				}
			}
		}
	}
	return ener;
}

vector<double> printPairBaselineEnergy(string &_seq, Options _opt, map<string, vector<double>> _innerPairs, map<string, vector<double>> _outerPairs){
	vector<double> ener;
	vector<string> pairs1;
	vector<string> pairs2;
	vector<int> dists;
	//ofstream pout;
	//string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	//string poutName  = dir + "/pairBaseline_1.out";	
	//pout.open(poutName.c_str());
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
	cout << "Sequence Length: " <<  _seq.length() << endl;
	for (uint i=0; i<pairs1.size(); i++){
		//cout << pos << endl;
		if (pairs1[i] != pairs2[i]){
			if (dists[i] == 1){
				pos++;
			}
			if (pos < 4 || pos > _seq.length()-5){
				for (map<string, vector<double>>::const_iterator it = _outerPairs.begin(); it != _outerPairs.end(); ++it){
					if (it->first == pairs1[i] || it->first == pairs2[i]){
						cout << pairs1[i] << " : " << pairs2[i] << "; " << dists[i] << " = " << it->second[dists[i]-1] << endl;
						ener.push_back(it->second[dists[i]-1]);
						it == _outerPairs.end();
					}
				}
			}
			else{
				for (map<string, vector<double>>::const_iterator it = _innerPairs.begin(); it != _innerPairs.end(); ++it){
					if (it->first == pairs1[i] || it->first == pairs2[i]){
						cout << pairs1[i] << " : " << pairs2[i] << "; " << dists[i] << " = " << it->second[dists[i]-1] << endl;
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
			if (pos < 4 || pos > _seq.length()-5){
				for (map<string, vector<double>>::const_iterator it = _outerPairs.begin(); it != _outerPairs.end(); ++it){
					if (it->first == pairs1[i]){
						cout << pairs1[i] << " : " << pairs2[i] << "; " << dists[i] << " = " << it->second[dists[i]-1] << endl;
						ener.push_back(it->second[dists[i]-1]);
						it == _outerPairs.end();
					}
				}
			}
			else{
				for (map<string, vector<double>>::const_iterator it = _innerPairs.begin(); it != _innerPairs.end(); ++it){
					if (it->first == pairs1[i]){
						cout << pairs1[i] << " : " << pairs2[i] << "; " << dists[i] << " = " << it->second[dists[i]-1] << endl;
						ener.push_back(it->second[dists[i]-1]);
						it == _innerPairs.end();
					}
				}
			}
		}
	}
	pairs1.clear();
	pairs2.clear();
	dists.clear();
	//pout.close();
	return ener;
}

SelfPairManager getMonomerSpm(System &_sys, Options& _opt, PolymerSequence &_PS, RandomNumberGenerator &_RNG){
	// For calculating time for the monomerSpm
	time_t startTimeMonoSpm, endTimeMonoSpm;
	double diffTimeMonoSpm;
	time(&startTimeMonoSpm);
	
	Chain & inputChain = _sys.getChain(0);

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
	if (!CSBMono.buildSystem(_PS)){
		cerr << "Unable to build system from " << _PS << endl;
	}

	//CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

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

	CSBMono.updateNonBonded(10,12,50);

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hBondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	PDBWriter writer;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	writer.open(dir + "/monomer.pdb");
	for (uint i=0; i<monoSys.atomSize(); i++){
		Atom at = monoSys.getAtom(i);
		if (!at.hasCoor()){
			cout << "Atom " << i << " was not assigned coordinates; program termination";
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
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", 1);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);
	
	/*****************************************************************************
	 *                 === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(true);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.calculateEnergies();
	
	time(&endTimeMonoSpm);
	diffTimeMonoSpm = difftime (endTimeMonoSpm, startTimeMonoSpm);

	return monoSpm;
}


double calcEnergyOfSequence(System &_pdb, Options &_opt, SelfPairManager &_spm, map<string,double> &_seqMCEnergies, vector<int> _varPosRotamers, string _sequence, int _seqNumber, double &_newSeqEnt, double &_baseline, double &_self, double &_pair, map<string, double> _selfMap, map<string, double> _seqEntMap, map<string,map<string,map<uint, double>>> _pairMap){
	
	map<string, double> seqEnergies;
	for (uint i=0; i<_opt.calcEnergyOfSequenceTerms.size(); i++){
		seqEnergies[_opt.calcEnergyOfSequenceTerms[i]] = 0;
	}
	
	time_t startTimeCalc, endTimeCalc;
	double diffTimeCalc;
	time(&startTimeCalc);
	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	//PDBWriter writer;
	//writer.open(_dir + "/sequence_" + to_string(_seqNumber) + ".pdb");

	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	//cout << polymerSeq << endl;	

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

	BaselineAAComposition bac(&sys);
	bac.readPenaltyFile("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/AACompositionPenalties.out");
	Eset->addInteraction(new BaselineAAComposition(bac));
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	//RNG.setSeed(_opt.seed); 
	RNG.setTimeBasedSeed();

	CSB.updateNonBonded(10,12,50);
	sys.buildAllAtoms();
	
	//sys.calcEnergy();

	loadRotamers(sys, sysRot, _opt.SL);
	loadInterfacialRotamers(sys, sysRot, _opt.SLInterface, _varPosRotamers);
	// Optimize Initial Starting Position (using Baseline to get back to original result)
	//SelfPairManager spm;
	//spm.seed(RNG.getSeed());
	//spm.setSystem(&sys);
	//spm.setVerbose(true);
	//spm.getMinStates()[0];
	//spm.updateWeights();
	//spm.setOnTheFly(true);
	//spm.calculateEnergies();

	//repackSideChains(spm, _opt.greedyCycles);
	//sys.setActiveRotamers(spm.getMinStates()[0]);
	//sys.saveAltCoor("savedBestState");
	
	_spm.setSystem(&sys);
	repackSideChains(_spm, _opt.greedyCycles);
	sys.setActiveRotamers(_spm.getMinStates()[0]);
	sys.calcEnergy();
	cout << Eset->getSummary() << endl;
	// TODO: troubleshoot this: repackSideChains gets me the right minimum because it recalculates energies...but the energies are already calculated previous. Can I somehow take those energies? it's slightly faster now, but probably could be even faster with this edit
	//_spm.runGreedyOptimizer(_opt.greedyCycles);

	//cout << "State Pre Calc: ";
	//for (uint i=0; i< _spm.getMinStates()[0].size(); i++){
	//	cout << _spm.getMinStates()[0][i] << "," ;
	//}
	//cout << endl;
	//cout << Eset->getSummary() << endl;
	//sys.setActiveRotamers(_spm.getMinStates()[0]);

	//sys.calcEnergy();
	//cout << "State: ";
	//for (uint i=0; i< _spm.getMinStates()[0].size(); i++){
	//	cout << _spm.getMinStates()[0][i] << "," ;
	//}
	//cout << Eset->getSummary() << endl;
	//
	//repackSideChains(_spm, _opt.greedyCycles);
	//sys.setActiveRotamers(_spm.getMinStates()[0]);

	//sys.calcEnergy();
	//cout << "State with Repack Calcs: ";
	//for (uint i=0; i< _spm.getMinStates()[0].size(); i++){
	//	cout << _spm.getMinStates()[0][i] << "," ;
	//}
	//cout << Eset->getSummary() << endl;
	//cout << Eset->getSummary() << endl;
	_newSeqEnt = Eset->getTermEnergy("BASELINE_ENTROPY_NORM");
	double self = Eset->getTermEnergy("BASELINE");
	double pair = Eset->getTermEnergy("BASELINE_PAIR");
	_self = self*2;
	_pair = pair*2;
	//cout << "SelfCalc = " << selfCalc << "; PairCalc = " << pairCalc << endl;
	//cout << "Self = " << self << "; Pair: " << pair << endl;
	
	_baseline = _self+_pair;
	seqEnergies["Baseline"] = (self+pair)*2;
	seqEnergies["DimerSelfBaseline"] = self*2;
	seqEnergies["DimerPairBaseline"] = pair*2;
	double nonBaselineEnergy = Eset->getTermEnergy("CHARMM_VDW")+Eset->getTermEnergy("SCWRL4_HBOND");
	seqEnergies["EnergyBeforeLocalMC"] = nonBaselineEnergy; //TODO: this fixes comparison to dimer energy without adding other terms, but not with IMM1. Is it possible for me to get that?

	map<string, double>::iterator itr;
	for (itr = seqEnergies.begin(); itr != seqEnergies.end(); ++itr){
		cout << itr->first << ": " << itr->second << endl;
		_seqMCEnergies.at(itr->first); //checks if the string exists in the map
		_seqMCEnergies[itr->first] = itr->second;
	}

	time(&endTimeCalc);
	diffTimeCalc = difftime (endTimeCalc, startTimeCalc);
	cout << "Time: " << diffTimeCalc << "s" << endl;
	return sys.calcEnergy();
}

double calcIMM1EnergyOfSequence(System &_pdb, Options &_opt, SelfPairManager &_spm, string _sequence, int _seqNumber, double &_newSeqEnt, double &_baseline, double &_self, double &_pair, map<string, double> _selfMap, map<string, double> _seqEntMap, map<string,map<string,map<uint, double>>> _pairMap){
	time_t startTimeCalc, endTimeCalc;
	double diffTimeCalc;
	time(&startTimeCalc);	
	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);

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
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);
	
	CSB.setBuildNonBondedInteractions(false);
	if(!CSB.buildSystem(polymerSeq)) {
		cerr << "Unable to build system from " << polymerSeq << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
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
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);
	Eset->setTermActive("BASELINE_ENTROPY_NORM", true);
	
	/******************************************************************************
	 *             === SETUP ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	Eset->setWeight("CHARMM_IMM1", 1);

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
	
	//sys.calcEnergy();

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

	time(&endTimeCalc);
	diffTimeCalc = difftime (endTimeCalc, startTimeCalc);
	cout << "Time: " << diffTimeCalc << "s" << endl;
	return sys.calcEnergy();
}

map<string,vector<double>> sequenceMC(System &_sys, Options &_opt, Transforms &_trans, SelfPairManager &_spm, SelfPairManager &_monoSpm, vector<unsigned int> _bestState, vector<string> &_seqs, vector<vector<unsigned int>> &_uniqueSeqs, vector<int> &_varPosList, vector<int> &_varPosRotamers, RandomNumberGenerator &_RNG, map<string, double> _selfMap, map<string, double> _seqEntMap, map<string,map<string,map<uint, double>>> _pairMap, ofstream &_ddout, ofstream &_dsout, string _dir){

	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);
	
	// Setup output energy landscape file
	ofstream elout;
	string eloutfile = _dir + "/energyLandscapeFile.out";
	elout.open(eloutfile.c_str());
	
	// Setup MonteCarloManager
	//MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, 20, MonteCarloManager::EXPONENTIAL, 20);//20 cycles and 20 rejects gets me 30K sequences; good first trial run
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, 20, MonteCarloManager::EXPONENTIAL, 10);
	//MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, 10, MonteCarloManager::EXPONENTIAL, 5);// as a trial for my map code
	MC.setRandomNumberGenerator(&_RNG);

	EnergySet *ESet = _sys.getEnergySet();
	Chain & chain = _sys.getChain("A");

	// Start from most probable state
	vector<unsigned int> currentState;
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _sys.calcEnergy();
	cout << ESet->getSummary() << endl;
	_sys.saveAltCoor("bestState");

	//initialize map for accepting energies
	map<string,vector<double>> seqMCOutput;
	vector<double> empty;
	
	//initialize map for accepting energies
	map<string,double> seqMCEnergies;

	for (uint i=0; i<_opt.sequenceMCEnergyTerms.size(); i++){
		seqMCEnergies[_opt.sequenceMCEnergyTerms[i]] = 0;
		seqMCOutput[_opt.sequenceMCEnergyTerms[i]] = empty;
		cout << _opt.sequenceMCEnergyTerms[i] << endl;
	}
	
	map<string, map<string, map<string,double>>> enerLandscapeData;
	//TODO: make maps like the above for all functions that output energies with the correct energy terms
	//then make sure to output them like the below part of this code (read energy terms, save them in a map, then output them to this map)

	double baseline = 0;
	double baselineSelf = 0;
	double baselinePair = 0;
	double monomerEnergy = 0;
	vector<vector<unsigned int>> monomerMinState;

	
	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	vector<unsigned int> stateVec = _bestState;
	currentState = _bestState;
	MC.setEner(bestEnergy);

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	//double enerToSave = 0;
	double prevMonomerEnergy = 0;
	double prevSeqEnt = 0;
	double currSeqEnt = 0;
	double prevSEProb = 0;
	string prevSeq = convertPolymerSeqToOneLetterSeq(chain);
	double prevSeqNum = 0;
	double currSeqNum = 0;

	//_seqs.push_back(prevSeq);

	vector<string> seqs;
	vector<string> allSeqs;
	
	while (!MC.getComplete()){
		try{
			if (_opt.verbose){
				cout << "Cycle #" << cycleCounter << "" << endl;
				cout << "Starting Seq: " << prevSeq << endl;
				_ddout << "Cycle #" << cycleCounter << endl;
				_ddout << "Starting Seq: " << prevSeq << endl;
			}
		
			// Get energy term and probability for the first sequence (these terms can then get replaced by future sequences that are accepted by MC)
			if (cycleCounter == 0){
				prevSeqEnt = ESet->getTermEnergy("BASELINE_ENTROPY_NORM");
				prevSEProb = exp(-prevSeqEnt/0.592);
				seqs.push_back(prevSeq);
				prevMonomerEnergy = computeMonomerEnergyNoMoves(_sys, _opt, _monoSpm, seqMCEnergies, prevSeq, _RNG, _opt.greedyCycles, _opt.MCCycles, _opt.MCMaxRejects);
				
				baselineSelf = ESet->getTermEnergy("BASELINE");
				baselinePair = ESet->getTermEnergy("BASELINE_PAIR");
				baselineSelf = baselineSelf*2;
				baselinePair = baselinePair*2;
				baseline = (baselineSelf+baselinePair);

				seqMCEnergies["Baseline"] = baseline;
				seqMCEnergies["DimerSelfBaseline"] = baselineSelf;
				seqMCEnergies["DimerPairBaseline"] = baselinePair;
				double nonBaselineEnergy = ESet->getTermEnergy("CHARMM_VDW")+ESet->getTermEnergy("SCWRL4_HBOND");
				seqMCEnergies["EnergyBeforeLocalMC"] = nonBaselineEnergy; //TODO: this fixes comparison to dimer energy without adding other terms, but not with IMM1. Is it possible for me to get that?
			
				map<string, double>::iterator itr;
				for (itr = seqMCEnergies.begin(); itr != seqMCEnergies.end(); ++itr){
					cout << itr->first << ": " << itr->second << endl;
					vector<double> energies = seqMCOutput.at(itr->first); //checks if the string exists in the map
					energies.push_back(itr->second);
					seqMCOutput[itr->first] = energies;
				}
				_seqs.push_back(prevSeq);
				allSeqs.push_back(prevSeq);
				//enerLandscapeData[prevSeq+":"+MslTools:intToString(_seqs.size())] = seqMCEnergies;
			}

			//TODO: make this into a function to get random sequence
			// Get a random integer to pick through the AA identities for a new random sequence to compare
			int rand = _RNG.getRandomInt(0, _varPosList.size()-1);
			int pos = _varPosList[rand];
			cout << "Position: " << pos << endl;
			int rand1 = _RNG.getRandomInt(0, ids.size()-1);
			string posId = _sys.getPosition(pos).getPositionId();
			string res = _sys.getPosition(pos).getResidueName();
			string randId;
			randId = ids[rand1];

			string resOneL = MslTools::getOneLetterCode(randId.substr(0,3));//For some reason some of my AAs have spaces between them, likely because that's how the Msl getStringVector() reads it	
			cout << "New Id: " << randId << " : " << resOneL << endl;

			string prevResID = MslTools::getOneLetterCode(res);
			cout << "Old Id: " << res << " : " << prevResID << endl;
			cout << prevResID << endl;


			string currSeq = "";
			for (uint i=0; i<prevSeq.length(); i++){
				if (i == pos){
					currSeq += resOneL;
				}
				else{
					currSeq += prevSeq[i];
				}
			}
			cout << prevResID << " at " << pos+_opt.thread << " switched to " << randId << endl;
			cout << "Prev seq:     " << prevSeq << endl;
			cout << "Switched seq: " << currSeq << endl;

			// Check the new random sequence to make sure it's not the same as the previous sequence
			bool same = false;
			for (uint s=0; s<allSeqs.size(); s++){
				if (currSeq == allSeqs[s]){
					same = true;
					s = allSeqs.size();
					currSeq.clear();
				}
			}
			int c = 0;
			if (same == true){
				while (same == true){
					cout << "while loop: " << c << endl;
					rand = _RNG.getRandomInt(0, _varPosList.size()-1);
					pos = _varPosList[rand];
					posId = _sys.getPosition(pos).getPositionId();
					res = _sys.getPosition(pos).getResidueName();
					rand1 = _RNG.getRandomInt(0, ids.size()-1);
					randId = ids[rand1];
					resOneL = MslTools::getOneLetterCode(randId.substr(0,3));//For some reason some of my AAs have spaces between them, likely because that's how the Msl getStringVector() reads it

					cout << "New Id: " << randId << " : " << resOneL << endl;
					prevResID = MslTools::getOneLetterCode(res);
					cout << "Old Id: " << res << " : " << prevResID << endl;
					for (uint i=0; i<prevSeq.length(); i++){
						if (i == pos){
							currSeq += resOneL;
						}
						else{
							currSeq += prevSeq[i];
						}
					}
					cout << "Previous seq: " << prevSeq << endl;
					cout << "Switched seq: " << currSeq << endl;
					cout << "Seq #" << c << ": " << currSeq << endl;

					same = false;
					for (uint s=0; s<allSeqs.size(); s++){
						if (currSeq == allSeqs[s]){
							same = true;
							s = allSeqs.size();
							currSeq.clear();
							c++;
						}
					}
					if (same == false){
						allSeqs.push_back(currSeq);
					}
					else if (c > 1000){//It has gotten stuck on a sequence that switching to any AA
						same = false;
						//MC.getComplete() = true;
					}
				}
			}//TODO: I wanted to check through all possible solutions, but coding it seems quite difficult. Maybe do that in the future
			
			// Calculate the energy of the random sequence (needed a separate function for it; spm defaulted to best sequence rather than just optimizing rotamers)
			baseline = 0;
			baselineSelf = 0;
			baselinePair = 0;
			currSeqEnt = 0;
			
			// TODO: assuming this SPM works well, should I also just make one monomerSPM and then referenee it in the same way?
			double currEnergy = calcEnergyOfSequence(_sys, _opt, _spm, seqMCEnergies, _varPosRotamers, currSeq, cycleCounter, currSeqEnt, baseline, baselineSelf, baselinePair, _selfMap, _seqEntMap, _pairMap);
	
			// Compute monomer energy
			monomerEnergy = computeMonomerEnergyNoMoves(_sys, _opt, _monoSpm, seqMCEnergies, currSeq, _RNG, _opt.greedyCycles, _opt.MCCycles, _opt.MCMaxRejects);
			
			// Converts the energy term (which actually saves the probability of the sequence in the whole system)
			// to the proper comparison of proportion to energy between individual sequences (done outside of the actual eneryg term)
			// TODO: understand how this works again (it's been a bit and I haven't had to explain it yet)
			double currSEProb = exp(-currSeqEnt/0.592);
			double totSEProb = prevSEProb+currSEProb;
			double prevSeqProp = prevSEProb/totSEProb;
			double currSeqProp = currSEProb/totSEProb;
			double prevEner = -log(prevSeqProp)*0.592*100;
			double currEner = -log(currSeqProp)*0.592*100;
		
			if (_opt.verbose){
				cout << "Previous Sequence: " << prevSeq << endl;
				cout << "Current Sequence: " << currSeq << endl;
				//cout << "Prev Seq Ent: " << prevSeqEnt << endl;
				//cout << "Prev Prob: " << prevSEProb << endl;
				//cout << "New Seq Ent: " << currSeqEnt << endl;
				//cout << "New Prob: " << currSEProb << endl;
				//cout << "Prev Seq Proportion: " << prevSeqProp << endl;
				//cout << "New Seq Proportion: " << currSeqProp << endl;
				//cout << "PrevEner = " << prevEner << endl;
				//cout << "NewEner = " << currEner << endl;
				//cout << "Diff = " << (prevEner-currEner) << endl;
				_ddout << "Previous Seq: " << prevSeq << endl;
				_ddout << "Current Seq:  " << currSeq << endl;
				_ddout << "Baseline:     " << baseline << endl;
				_ddout << "Monomer:      " << monomerEnergy << endl;
				_ddout << "Prev Seq Ent: " << prevSeqEnt << endl;
				_ddout << "Prev Prob:    " << prevSEProb << endl;
				_ddout << "New Seq Ent:  " << currSeqEnt << endl;
				_ddout << "New Prob:     " << currSEProb << endl;
				_ddout << "Prev Seq Proportion: " << prevSeqProp << endl;
				_ddout << "New Seq Proportion:  " << currSeqProp << endl;
				_ddout << "PrevEner =    " << prevEner << endl;
				_ddout << "NewEner =     " << currEner << endl;
				_ddout << "Diff =        " << (prevEner-currEner) << endl;
			}
		
			// Sets the energy to compare acceptance to as the total energy with the seqEntropy term converted from the energy term
			// TODO: determine if I should use baselines here or just the energy? Or the actual mnomer too since that is an option here now
			// The below includes the baseline energy, which is basically monomer energy
			double bestEnergyTotal = bestEnergy-prevSeqEnt+prevEner;
			double currEnergyTotal = currEnergy-currSeqEnt+currEner;//TODO: I would like to include these in analysis somehow, but can't thikn of a way to do so
			MC.setEner(bestEnergyTotal);
			
			cout << "Best Energy: " << bestEnergyTotal << endl;
			cout << "New Energy: " << currEnergyTotal << endl;
			_ddout << "Best Energy: " << bestEnergyTotal << endl;
			_ddout << "New Energy: " << currEnergyTotal << endl;

			//TODO: make a sequence energy landscape: save all the energies of your random sequences as well as the sequences themselves in a energyLandscape file
			//Also while at it, make sure to make a sequence checker for a list of sequences so it doesn't do that same sequence over and over again
			// Adds the first sequence and energy to the sequences list
			//
			vector<double> enerLandscapeVector;
			enerLandscapeVector.push_back(bestEnergyTotal);
			enerLandscapeVector.push_back(currEnergyTotal);
			enerLandscapeVector.push_back(prevSeqProp);
			enerLandscapeVector.push_back(currSeqProp);
			enerLandscapeVector.push_back(prevEner);
			enerLandscapeVector.push_back(currEner);
			enerLandscapeVector.push_back(prevSeqNum);
			enerLandscapeVector.push_back(currSeqNum);
			for (uint i=0; i<_opt.enerLandscapeTerms.size(); i++){
				enerLandscapeData[prevSeq][currSeq][_opt.enerLandscapeTerms[i]] = enerLandscapeVector[i];//TODO: add in something that says not to work if these are not the same size
			}
			if (cycleCounter == 0){
				sequences[bestEnergy-prevSeqEnt] = prevSeq;
			}

			// MC accept and reject conditions
			if (!MC.accept(currEnergyTotal)){
				_sys.applySavedCoor("bestState");
				currSeqNum++;
				//_sys.setActiveRotamer(prevId, 1);

				if (_opt.verbose){
					cout << "State not accepted, E= " << currEnergyTotal << endl << endl;
					_ddout << "State not accepted, E= " << currEnergyTotal << endl << endl;
				}
			} else {

				bool sameseq = sameSequenceChecker(currSeq, _seqs);
				cout << sameseq << endl;
				if (sameseq == false){
					prevSeqNum++;
					currSeqNum = 1;
					_sys.saveAltCoor("bestState");
					sequences[currEnergy-currSeqEnt] = currSeq;
					bestEnergy = currEnergy;
					MC.setEner(currEnergyTotal);
					prevSeqEnt = currSeqEnt;
					prevSEProb = currSEProb;
					prevSeq = currSeq;
					prevMonomerEnergy = monomerEnergy;
					//TODO: this currently pushes back all sequences (even if they are the same because my code doesn't account too well for that yet; this will make the code overall take longer because it will do localMC checks on the same sequence (if it were a map it wouldn't matter, but I didn't implement that properly yet); fix one of these two things
					seqs.push_back(currSeq);
					map<string, double>::iterator itr;
					for (itr = seqMCEnergies.begin(); itr != seqMCEnergies.end(); ++itr){
						vector<double> energies = seqMCOutput.at(itr->first); //checks if the string exists in the map
						energies.push_back(itr->second);
						seqMCOutput[itr->first] = energies;
						cout << itr->first << ": " << energies.size() << ": " << itr->second << endl;
					}
				}

				if (_opt.verbose) {
					cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << " E= " << currEnergy << " : MCEner= " << currEnergyTotal << endl << endl;//includes baselineComposition
					_ddout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << " E= " << currEnergy << " : MCEner= " << currEnergyTotal << endl << endl;//includes baselineComposition
					_dsout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << " E= " << currEnergy << " : MCEner= " << currEnergyTotal << endl << endl;
				}
			}
			cycleCounter++;
		}
		catch (const out_of_range& e){
			//_sys.setActiveRotamers(prevStateVec);
			cout << "Out of bounds exception caught! Needs to be fixed!" << endl;
			_ddout << "Out of bounds exception caught! Needs to be fixed!" << endl;
			break;	
		}
	}

	// Print Energy Landscape

	elout << "\tSequence\t";
	for (uint i=0; i<_opt.enerLandscapeTerms.size(); i++){
		string enerTerm = _opt.enerLandscapeTerms[i];
		elout << enerTerm << "\t";
	}
	elout << endl;

	map<string,map<string,map<string,double>>>::iterator seqItr1;
	map<string,map<string,double>>::iterator seqItr2;
	for(seqItr1 = enerLandscapeData.begin(); seqItr1 != enerLandscapeData.end(); ++seqItr1){
		cout << seqItr1->first << "\t";
		for(seqItr2 = seqItr1->second.begin(); seqItr2 != seqItr1->second.end(); ++seqItr2){
			elout << "EnerLandscape: " << seqItr1->first << "\t" << seqItr2->first << "\t";
			for (uint i=0; i<_opt.enerLandscapeTerms.size(); i++){
				string enerTerm = _opt.enerLandscapeTerms[i];
				double ener = seqItr2->second.at(enerTerm);
				if (i != _opt.enerLandscapeTerms.size()-1){
					cout << ener << "\t";
					elout << ener << "\t";
				}
				else{
					cout << ener << endl;
					elout << ener << endl;
				}
			}
		}
	}

	cout << "Best Energy: " << bestEnergy << endl;
	cout << "#Seqs: " << seqs.size() << endl;
	_ddout << "Best Energy: " << bestEnergy << endl;
	_ddout << "#Seqs: " << seqs.size() << endl;
	_dsout << "Best Energy: " << bestEnergy << endl;
	_dsout << "#Seqs: " << seqs.size() << endl;
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);
	cout << "SMC Time: " << diffTimeSMC << "s" << endl;
	_ddout << "SMC Time: " << diffTimeSMC << "s" << endl << endl;
	_dsout << "SMC Time: " << diffTimeSMC << "s" << endl << endl;
	return seqMCOutput;
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
	ofstream sout;
	ofstream dout;
	ofstream eout;

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

	string foutfile = ddir + "/finalDesignInfo.out";
	string soutfile = ddir + "/design_summary.out";
	string doutfile = ddir + "/design_details.out";
	string eoutfile  = ddir + "/errors.out";

	fout.open(foutfile.c_str());
	eout.open(eoutfile.c_str());
	sout.open(soutfile.c_str());
	dout.open(doutfile.c_str());

	eout << date << endl;
	fout << date << endl;
	dout << date << endl;
	sout << date << endl;

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
	 *   === CHOOSE XSHIFT AND CROSSING ANGLE FROM DISTRIBUTION OF PDB DATA ===
	 ******************************************************************************/
	if (opt.useGeoFromPDBData){
		RandomNumberGenerator geoRNG;
		geoRNG.setSeed(MslTools::toInt(opt.runNumber));
		//TODO: automate this so that if I decide to use different kdes I can just plug the file in
		double ad_normSD = 0.00107102*3;
		double ad_normMean = 0.001041667;
		
		double ad_randNorm = geoRNG.getRandomDouble(0,1);//This is good in principle, but it's usually giving large values...how do I get it to be more likely to give small values?
		cout << ad_randNorm << endl;
		ad_randNorm = (ad_randNorm*ad_normSD)+ad_normMean;
		cout << "Rand Norm " << ": " << ad_randNorm << endl;
		
		string ad_kdeFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/kde.txt";//I think this is the file, but make sure to check
	
		// Setup output files
		ofstream gkout;
		string gkoutfile = ddir + "/gaussian_kde_values.out";
		gkout.open(gkoutfile.c_str());
		gkout << "Randomly chosen kde values: xShift: Angle: Rot: Z" << endl;
		gkout << "Random Kde Values:	" << ad_randNorm << endl;
		getGaussianKdeValues(ad_kdeFile, opt, ad_randNorm, opt.xShift, opt.crossingAngle, gkout);
		gkout << "Geometry:	" << opt.xShift << "	" << opt.crossingAngle << "	" << opt.axialRotation << "	" << opt.zShift << endl;	
		gkout.close();
	
		/******************************************************************************
		 * === CHOOSE AXIAL ROTATION AND ZSHIFT FROM DENSITY HISTOGRAMS OF PDB DATA ===
		 ******************************************************************************/
		string otherGeometriesFile; 
		if (opt.crossingAngle < 0){
			otherGeometriesFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/2021_05_13_geoBinsRight.txt";
		}
		else{
			otherGeometriesFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/2021_05_13_geoBinsLeft.txt";
		}
		Reader axZreader(otherGeometriesFile);
		axZreader.open();
		map<vector<double>,vector<double>> geoBinEnergies;
	
		if(!(axZreader.is_open())){
			cerr << "WARNING: Unable to open " << otherGeometriesFile << endl;
			exit(0);
		}
	
		vector<string> lines = axZreader.getAllLines();
	
		for (int i=1; i<lines.size(); i++){
			vector<string> tokens = MslTools::tokenize(lines[i], "\t");
			if(tokens.size() < 1){
				continue;
			}
			if(tokens.size() < 6){
				cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: xbin1(double) xbin2(double) axRotbin1(double) axRotbin2(double) zShiftbin1(double) zShiftbin2(double)";
				continue;
			}
			vector<double> xBins;
			vector<double> axZ;
			for (uint j=0; j<tokens.size(); j++){
				if (j<2){
					xBins.push_back(MslTools::toDouble(tokens[j]));
				}
				else{
					axZ.push_back(MslTools::toDouble(tokens[j]));
				}
			}
			geoBinEnergies[xBins] = axZ;
		}
	
		map<vector<double>, vector<double>>::iterator geoItr;
		for (geoItr = geoBinEnergies.begin(); geoItr != geoBinEnergies.end(); ++geoItr){
			if (opt.xShift >= geoItr->first[0] && opt.xShift <= geoItr->first[1]){
				if (geoItr->second.size()){
					opt.axialRotation = geoRNG.getRandomDouble(geoItr->second[0],geoItr->second[1]);
					opt.zShift = geoRNG.getRandomDouble(geoItr->second[2],geoItr->second[3]);
				}
				else{
					opt.zShift = geoRNG.getRandomDouble(geoItr->second[2],geoItr->second[3]);
					int rand = geoRNG.getRandomInt(1,2);
					if (rand == 1){
						opt.axialRotation = geoRNG.getRandomDouble(geoItr->second[0],geoItr->second[1]);
					}
					else{
						opt.axialRotation = geoRNG.getRandomDouble(geoItr->second[4],geoItr->second[5]);
					}
				}
			}
			else{
				continue;
			}
		}
	}
	//TODO: read through file, see what the xShift is, get values for both of the others between two values

	/******************************************************************************
	 *                        === SETUP SUMMARY FILE ===
	 ******************************************************************************/
	cout  << "***STARTING GEOMETRY:***" << endl;
	cout  << "xShift: " << opt.xShift << endl;
	cout  << "crossingAngle: " << opt.crossingAngle << endl;
	cout  << "axialRotation: " << opt.axialRotation << endl;
	cout  << "zShift: " << opt.zShift << endl << endl;
	fout  << "***STARTING GEOMETRY:***" << endl;
	fout  << "xShift: " << opt.xShift << endl;
	fout  << "crossingAngle: " << opt.crossingAngle << endl;
	fout  << "axialRotation: " << opt.axialRotation << endl;
	fout  << "zShift: " << opt.zShift << endl << endl;
	sout << "***STARTING GEOMETRY:***" << endl;
	sout << "xShift: " << opt.xShift << endl;
	sout << "crossingAngle: " << opt.crossingAngle << endl;
	sout << "axialRotation: " << opt.axialRotation << endl;
	sout << "zShift: " << opt.zShift << endl << endl;
	dout << "***STARTING GEOMETRY:***" << endl;
	dout << "xShift: " << opt.xShift << endl;
	dout << "crossingAngle: " << opt.crossingAngle << endl;
	dout << "axialRotation: " << opt.axialRotation << endl;
	dout << "zShift: " << opt.zShift << endl << endl;
	
	/******************************************************************************
	 *                     === READ IN GEOMETRY FILE ===
	 ******************************************************************************/
	vector<string> fileVec;
	readGeometryFile(opt.helixGeoFile, fileVec);

	/******************************************************************************
	 *                         === GENERATE POLYGLY ===
	 ******************************************************************************/
	string polyGly = generatePolymerSequence("A", opt.sequenceLength, opt.thread);
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

	//vector<Position*>& positions = sys.getPositions();

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
	 *                   === IDENTIFY INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	vector<int> pos;
	vector<double> dist;
	vector<vector<vector<double>>> posDistVector;
	identifyInterface(sys, opt, PS, pos, ddir, axis);

	vector<int> varPosRotamers = interface01(sys, pos);
	vector<int> varPos = interface01(sys, pos, false);

	//String for the positions of the sequences that are considered interface for positions
	string interfaceDiffAA = "";
	for (uint i=0; i<varPos.size(); i++){
		if (i == 21){
			i = varPos.size();
			dout << endl;
		}
		else{
			dout << varPos[i];
			interfaceDiffAA += MslTools::intToString(varPos[i]);
		}
	}
	cout << "VarPos: " << varPos.size() << ": " << interfaceDiffAA << endl;
	
	//String for the positions of the sequence that are considered interface with higher rotamers
	string interfacePos = "";
	for (uint i=0; i<varPosRotamers.size(); i++){
		if (i == 21){
			i = varPosRotamers.size();
			dout << endl;
		}
		else{
			dout << varPosRotamers[i];
			interfacePos += MslTools::intToString(varPosRotamers[i]);
		}
	}
	cout << "VarPosRots: " << interfacePos << endl;

	//String for the alternateIds at the interface
	string alternateIds = "";
	for (uint i=0; i<opt.Ids.size(); i++){
		if (i == opt.Ids.size()-1){
			alternateIds += opt.Ids[i];
			dout << endl;
		}
		else{
			alternateIds += opt.Ids[i] += " ";
		}
	}
	cout << alternateIds << endl;

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	sys.buildAllAtoms();

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
	dout << "Sequence: " << polyLeu << endl;
	dout << "DiffPos:  " << interfaceDiffAA << endl;
	dout << "DiffRots: " << interfacePos << endl;
	dout << "Alternate Ids: " << alternateIds << endl;
	dout << "RotLevel: " << opt.SL << endl;
	dout << "RotLevelInterface: " << opt.SLInterface << endl << endl;
	cout << "Sequence: " << polyLeu << endl;
	cout << "DiffPos:  " << interfaceDiffAA << endl;
	cout << "DiffRots: " << interfacePos << endl;
	cout << "Alternate Ids: " << alternateIds << endl;
	cout << "RotLevel: " << opt.SL << endl;
	cout << "RotLevelInterface: " << opt.SLInterface << endl << endl;

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
	//AtomPointerVector & apvChainAL = chainAL.getAllAtomPointers();
	//AtomPointerVector & apvChainBL = chainBL.getAllAtomPointers();
	//vector<Position*>& positionsL = sysL.getPositions();

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
	buildSelfInteractions(sysL, selfMap);
	buildPairInteractions(sysL, pairMap);

	//initialize baselineAAComposition energy map (based on Rosetta baselineAAComposition to prevent unlikely sequences by adding energy (ex. if more than 2PHE in sequence, add 100 energy score for each additional PHE)
	BaselineAAComposition bac(&sysL);
	bac.readPenaltyFile("/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/AACompositionPenalties.out");
	EsetL->addInteraction(new BaselineAAComposition(bac));

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
	loadInterfacialRotamers(sysL, sysRot, opt.SLInterface, varPosRotamers);

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	//Random Number Generator
	RandomNumberGenerator RNG;
	RNG.setTimeBasedSeed();
	
	CSBL.updateNonBonded();
	sysL.buildAllAtoms();
	dout << "***SEQUENCE OPTIMIZATION***" << endl;
	sout << "***SEQUENCE OPTIMIZATION***" << endl;
	cout << "***SEQUENCE OPTIMIZATION***" << endl;
	cout << EsetL->getSummary() << endl;

	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sysL);
	spm.setVerbose(true);
	spm.setRunDEE(opt.runDEESingles, opt.runDEEPairs);
	spm.setOnTheFly(true);
	
	//Setup running SCMF or UnbiasedMC
	if (opt.runSCMF == true){
		dout << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	}
	else{
		dout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	dout << "Starting SelfPairManager Optimization..." << endl << endl;
	dout << "SPM Seed: " << RNG.getSeed() << endl;
	sysL.calcEnergy();
	dout << "Total interactions calced: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	cout << "Total interactions calced: " << EsetL->getTotalNumberOfInteractionsCalculated() << endl;
	cout << EsetL->getSummary() << endl;
	time(&spmStart);
	spm.runOptimizer();
	time(&spmEnd);
	
	dout << endl << "End SelfPairManager Optimization" << endl;
	cout << endl << "End SelfPairManager Optimization" << endl;
	spmTime = difftime (spmEnd, spmStart);
	
	//TODO: would be nice to output the run times for DEE and SCMF
	dout << "SelfPairManager runOptimizer time: " << spmTime << " seconds" << endl;
	//dout << "States after DEE: " << spm.getDEEAliveRotamers().size() << endl;
	dout << "States after SCMF: " << spm.getBestSCMFBiasedMCStates().size() << endl << endl;
	fout << "SelfPairManager runOptimizer time: " << spmTime << " seconds" << endl;
	//cout << "States after DEE: " << spm.getDEEAliveRotamers().size() << endl;
	fout << "States after SCMF: " << spm.getBestSCMFBiasedMCStates().size() << endl << endl;

	// Setup Monomer SelfPairManager: should increase the speed of computing the monomer energies
	SelfPairManager monoSpm = getMonomerSpm(sys, opt, PL, RNG);

	/******************************************************************************
	 *           === METHODS FOR DETERMINING ALTERNATE SEQUENCES ===
	 ******************************************************************************/
	vector<string> seqs;

	//Initialize energyMap to hold all energies for output into a summary file
	//TODO: I need to fix these to make it save energies for each sequence rather than this weird thing; it definitely works for now tested for seq and monomerEnergies with mutateSequenceDesign (05_11_2021), but would like to make it better
	//map<string, map<string,double>> seqEnergyMap;
	map<string,vector<double>> energyMap;//TODO: these vector energy maps make me nervous, but I can't see a reason why they'd be wrong
	// Would be better if I just paired each string map with a sequence string map...maybe do that at some point, but not today (05_11_2021 errors in the monomer code really messed me up)
	vector<double> empty; //TODO: is there a more elegant way to do this?
	for (uint i=0; i<opt.energyTermsToOutput.size(); i++){
		energyMap[opt.energyTermsToOutput[i]] = empty;	
	}
	vector<vector<unsigned int>> uniqueSeqs;

	//Initialize maps for outputs for sequenceMC and localMC
	map<string,vector<double>> seqMCOutput;
	map<string,vector<double>> monomerIMM1Output;
	for (uint i=0; i<opt.monomerIMM1EnergyTerms.size(); i++){
		monomerIMM1Output[opt.monomerIMM1EnergyTerms[i]] = empty;
	}
	map<string,vector<double>> localMCOutput;
	for (uint i=0; i<opt.dimerEnergyTerms.size();i++){
		localMCOutput[opt.dimerEnergyTerms[i]] = empty;
	}

	map<string, double> localMCEnergies;
	vector<double> monomerEnergies;
	vector<vector<vector<unsigned int>>> monomerMinState;
	
	if (opt.runMCAfterSPM){
		/******************************************************************************
		 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
		 ******************************************************************************/
		cout << "Run MC on sequence after SPM: true" << endl;
		sout << "Run MC on sequence after SPM: true" << endl;
		
		//DEE alive rotamers to be randomly oriented in
		//vector<vector<unsigned int>> deeAliveRotamers = spm.getDEEAliveRotamers();
		vector<int> varPosList = getVariablePos(varPos);
		vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
		
		BaselineSequenceEntropyNormalized bs(&sysL);
		bs.setMap(seqEntMap);
		EsetL->addInteraction(new BaselineSequenceEntropyNormalized(bs));
		sysL.calcEnergy();
		cout << "BSE WORKED!" << endl;

		//TODO: add an option for a multiplier for emphasizing larger/smaller differences between sequence entropies (currently defaults to 100)
		seqMCOutput = sequenceMC(sysL, opt, trans, spm, monoSpm, bestState, seqs, uniqueSeqs, varPosList, varPosRotamers, RNG, selfMap, seqEntMap, pairMap, dout, sout, ddir);

		cout << "Monomer energies size: " << monomerEnergies.size() << endl;
		cout << "map size: " << seqMCOutput.size() << endl;

		map<string, double> monIMM1Energies;
		for (uint m=0; m<seqs.size(); m++){

			// TODO: I've already changed the below function, but I think it's finally time for me to move things over to a map...do that and then check how well it workd
			//double monIMM1 = computeMonomerIMM1Energy(sysL, trans, opt, monomerMinState, m, helicalAxis, RNG, sout, opt.greedyCycles, opt.MCCycles, opt.MCMaxRejects);
			double monIMM1 = computeMonomerEnergy(sysL, opt, trans, monIMM1Energies, seqs[m], RNG, opt.greedyCycles, opt.MCCycles, opt.MCMaxRejects);
			//monomerIMM1Energies.push_back(monIMM1);
			cout << seqs[m] << " Monomer Energy: " << monIMM1 << endl;
			monomerEnergies.push_back(monIMM1);
			cout << monIMM1Energies.size() << endl;
			map<string, double>::iterator itr;
			for (itr = monIMM1Energies.begin(); itr != monIMM1Energies.end(); ++itr){
				cout << itr->first << ": " << itr->second << endl;
				vector<double> energies = monomerIMM1Output.at(itr->first); //checks if the string exists in the map
				energies.push_back(itr->second);
				monomerIMM1Output[itr->first] = energies;
			}
		}
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

	//TODO: add in a way to rid of sequences that are outside of the baseline-monomer range that I'm interested in (filtering step...needs to get rid of all parts, could saving things in a map be better somehow...?

	/******************************************************************************
	 *               === LOCAL REPACKS ON EACH UNIQUE SEQUENCE ===
	 ******************************************************************************/
	vector<double> xShifts;
	vector<double> crossingAngles;
	vector<double> axialRotations;
	vector<double> zShifts;

	dout << "***LOCAL REPACKS ON UNIQUE SEQUENCES***" << endl << endl;
	sout << "***LOCAL REPACKS ON UNIQUE SEQUENCES***" << endl << endl;
	
	// Setup timer for LocalMC Cycles
	time_t startTimeLMC, endTimeLMC;
	double diffTimeLMC;
	time(&startTimeLMC);
	//TODO: change this to go through the energies of dimers from the sequenceMC and do them in order
	for(uint i=0; i<seqs.size(); i++){
		/******************************************************************************
		 *                      === SETUP OUTPUT FILES ===
		 ******************************************************************************/
		cout << "Sequence #" << i << ": " << seqs[i] << endl;
		fout << "Sequence #" << i << ": " << seqs[i] << endl;
		sout << "Sequence #" << i << ": " << seqs[i] << endl;
		dout << "Sequence #" << i << ": " << seqs[i] << endl;
		
		//The below compares the monomeric state and does a local repack, getting it to the best energetic state by comparing dimer to monomer energies
		System finalPdb = localMC(sysL, helicalAxis, opt, i, localMCEnergies, seqs[i], varPosRotamers, monomerEnergies[i], xShifts, crossingAngles, axialRotations, zShifts, dout, sout, ddir, writer);
		dout << "Energy after local Monte Carlo: " << finalPdb.calcEnergy() << endl;

		map<string, double>::iterator itr;
		for (itr = localMCEnergies.begin(); itr != localMCEnergies.end(); ++itr){
			cout << itr->first << ": " << itr->second << endl;
			vector<double> energies = localMCOutput.at(itr->first); //checks if the string exists in the map
			energies.push_back(itr->second);
			localMCOutput[itr->first] = energies;
		}

		//Calculate Dimer energy for output
		double dimerEnergy = finalPdb.calcEnergy();
		double totalEnergy = dimerEnergy-monomerEnergies[i];
		dout << "Dimer Energy: " << dimerEnergy << endl << endl;
		sout << "Dimer Energy: " << dimerEnergy << endl << endl;
	}
	time(&endTimeLMC);
	diffTimeLMC = difftime (endTimeLMC, startTimeLMC);
	cout << "End Local Monte Carlo Repacks: " << diffTimeLMC << "s" << endl;
	sout << "End Local Monte Carlo Repacks: " << diffTimeLMC << "s" << endl;
	dout << "End Local Monte Carlo Repacks: " << diffTimeLMC << "s" << endl;

	/******************************************************************************
	 *           === INPUT ENERGIES INTO ENERGY MAP FOR OUTPUT FILE ===
	 ******************************************************************************/
	map<string, vector<double>>::iterator itr;
	for (itr = seqMCOutput.begin(); itr != seqMCOutput.end(); ++itr){
		vector<double> energies = seqMCOutput.at(itr->first); //checks if the string exists in the map
		cout << itr->first << " size: " << energies.size() << endl;
		energyMap.at(itr->first);
		energyMap[itr->first] = energies;
	}
	for (itr = monomerIMM1Output.begin(); itr != monomerIMM1Output.end(); ++itr){
		vector<double> energies = monomerIMM1Output.at(itr->first); //checks if the string exists in the map
		cout << itr->first << " size: " << energies.size() << endl;
		energyMap.at(itr->first);
		energyMap[itr->first] = energies;
	}
	for (itr = localMCOutput.begin(); itr != localMCOutput.end(); ++itr){
		vector<double> energies = localMCOutput.at(itr->first); //checks if the string exists in the map
		cout << itr->first << " size: " << energies.size() << endl;
		energyMap.at(itr->first);
		energyMap[itr->first] = energies;
	}
	cout << "# seqs: " << seqs.size() << endl;

	/******************************************************************************
	 *         === WRITE OUT BASELINE TO MONOMER COMPARISON BY SEQUENCE ===
	 ******************************************************************************/
	int runNumber = MslTools::toInt(opt.runNumber);
	fout << "\t";	
	cout << "\t";	
	for (uint i=0; i<opt.energyTermsToOutput.size(); i++){
		fout << opt.energyTermsToOutput[i] << "\t";
		cout << opt.energyTermsToOutput[i] << "\t";
	}
	fout << "startXShift\tstartCrossingAngle\tstartAxialRotation\tstartZShift\txShift\tcrossingAngle\taxialRotation\tzShift\tSequence\tPDBPath\tThread\tInterfacePositions" << endl;
	string tab = "\t";
	cout << "# Sequences: " << seqs.size() << endl;
	map<string,vector<double>>::iterator energyItr;
	for (uint i=0; i<seqs.size(); i++){
		for (uint j=0; j<opt.energyTermsToOutput.size(); j++){
			if (opt.energyTermsToOutput[j].compare("Total") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << "Sequence Info: " << energyMap.at("Dimer")[i]-energyMap.at("Monomerw/IMM1")[i] << "\t";
				fout << "Sequence Info: " << energyMap.at("Dimer")[i]-energyMap.at("Monomerw/IMM1")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("Baseline-Monomer") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("Baseline")[i]+energyMap.at("Monomer")[i] << "\t";
				fout << energyMap.at("Baseline")[i]+energyMap.at("Monomer")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("HbondDifference") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("HbondDimer")[i]-energyMap.at("HbondMonomerw/IMM1")[i] << "\t";
				fout << energyMap.at("HbondDimer")[i]-energyMap.at("HbondMonomerw/IMM1")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("VDWDifference") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("VDWDimer")[i]-energyMap.at("VDWMonomerw/IMM1")[i] << "\t";
				fout << energyMap.at("VDWDimer")[i]-energyMap.at("VDWMonomerw/IMM1")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("IMM1Difference") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("IMM1Dimer")[i]-energyMap.at("IMM1Monomer")[i] << "\t";
				fout << energyMap.at("IMM1Dimer")[i]-energyMap.at("IMM1Monomer")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("SelfDifference") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("MonomerSelfBaseline")[i]+energyMap.at("DimerSelfBaseline")[i] << "\t";
				fout << energyMap.at("MonomerSelfBaseline")[i]+energyMap.at("DimerSelfBaseline")[i] << "\t";
			}
			else if (opt.energyTermsToOutput[j].compare("PairDifference") == 0){
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at("MonomerPairBaseline")[i]+energyMap.at("DimerPairBaseline")[i] << "\t";
				fout << energyMap.at("MonomerPairBaseline")[i]+energyMap.at("DimerPairBaseline")[i] << "\t";
			}
			else{
				cout << opt.energyTermsToOutput[j] << ": "; 
				cout << energyMap.at(opt.energyTermsToOutput[j])[i] << "\t";
				fout << energyMap.at(opt.energyTermsToOutput[j])[i] << "\t";
			}
		}
		cout << opt.xShift << tab << opt.crossingAngle << tab << opt.axialRotation << tab << opt.zShift << tab << xShifts[i] << tab << crossingAngles[i] << tab << axialRotations[i] << tab << zShifts[i] << tab << seqs[i] << tab << ddir << "/design_" << i << ".pdb" << tab << opt.thread << tab << interfacePos << endl << endl;
		fout << opt.xShift << tab << opt.crossingAngle << tab << opt.axialRotation << tab << opt.zShift << tab << xShifts[i] << tab << crossingAngles[i] << tab << axialRotations[i] << tab << zShifts[i] << tab << seqs[i] << tab << ddir << "/design_" << i << ".pdb" << tab << opt.thread << tab << interfacePos << endl << endl;
	}

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	dout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	fout  << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	cout  << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

	writer.close();
	fout.close();
	eout.close();
	dout.close();
	sout.close();
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
	// optional
	opt.allowed.push_back("useGeoFromPDBData");
	
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("xShift");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("leftHanded");
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
	opt.allowed.push_back("SLInterface");
	
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
	
	//Energy Terms to Output
	opt.allowed.push_back("monomerEnergyTerms");
	opt.allowed.push_back("monomerIMM1EnergyTerms");
	opt.allowed.push_back("dimerEnergyTerms");
	opt.allowed.push_back("calcEnergyOfSequenceTerms");
	opt.allowed.push_back("sequenceMCEnergyTerms");
	opt.allowed.push_back("enerLandscapeTerms");
	opt.allowed.push_back("energyTermsToOutput");

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
	
	opt.useGeoFromPDBData = OP.getBool("useGeoFromPDBData");
	if (OP.fail()) {
		opt.warningMessages += "useGeoFromPDBData not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.useGeoFromPDBData = false;
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
	opt.leftHanded = OP.getBool("leftHanded");
	if (OP.fail()){
		opt.warningMessages += "leftHanded not specified, defaulting to true\n";
		opt.warningFlag = true;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (opt.leftHanded == false){
		opt.crossingAngle = -opt.crossingAngle;
	}
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
	opt.SLInterface = OP.getString("SLInterface");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	}
	else{
		opt.SLInterface = "SL"+opt.SLInterface;
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

	//Energy Terms to Output
	opt.monomerEnergyTerms = OP.getStringVector("monomerEnergyTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify monomer energy terms, defaulting to Monomer, VDWMonomer, HbondMonomer, IMM1Monomer, MonomerSelfBaseline, MonomerPairBaseline\n";
		opt.warningFlag = true;
	}
	opt.monomerIMM1EnergyTerms = OP.getStringVector("monomerIMM1EnergyTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify monomerIMM1 energy terms, defaulting to Monomerw/IMM1, VDWMonomerw/IMM1, HbondMonomerw/IMM1, Monomerw/IMM1SelfBaseline, Monomerw/IMM1PairBaseline\n";
		opt.warningFlag = true;
	}
	opt.dimerEnergyTerms = OP.getStringVector("dimerEnergyTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify dimer energy terms, defaulting to Dimer, VDWDimer, HbondDimer, IMM1Dimer, DimerSelfBaseline, DimerPairBaseline\n";
		opt.warningFlag = true;
	}
	opt.calcEnergyOfSequenceTerms = OP.getStringVector("calcEnergyOfSequenceTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify dimer energy terms, defaulting to Dimer, VDWDimer, HbondDimer, IMM1Dimer, DimerSelfBaseline, DimerPairBaseline\n";
		opt.warningFlag = true;
	}
	opt.sequenceMCEnergyTerms = OP.getStringVector("sequenceMCEnergyTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify sequence energy terms, defaulting to Monomer, VDWMonomer, HbondMonomer, IMM1Monomer, MonomerSelfBaseline, MonomerPairBaseline, DimerSelfBaseline, DimerPairBaseline, EnergyBeforeLocalMC\n";
		opt.warningFlag = true;
	}
	opt.enerLandscapeTerms = OP.getStringVector("enerLandscapeTerms");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify sequence energy terms, defaulting to PrevEnergyTotal, CurrEnergyTotal, PrevSeqProp, CurrSeqProp, PrevSeqEntropy, CurrSeqEntropy\n";
		opt.warningFlag = true;
	}
	opt.energyTermsToOutput = OP.getStringVector("energyTermsToOutput");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify sequence energy terms, defaulting to Monomer, VDWMonomer, HbondMonomer, IMM1Monomer, MonomerSelfBaseline, MonomerPairBaseline, DimerSelfBaseline, DimerPairBaseline, EnergyBeforeLocalMC\n";
		opt.warningFlag = true;
	}
	
	opt.rerunConf = OP.getConfFile();

	return opt;
}
