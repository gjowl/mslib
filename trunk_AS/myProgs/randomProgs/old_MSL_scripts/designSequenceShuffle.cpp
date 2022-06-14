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
string programName = "designSequenceShuffle";
string programDescription = "This program is derived from generateRandomSequence; it is meant to go through all the permutations AA sequences from the designed sequence by shuffling the AAs and/or makding single mutations at each of the positions";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "2 November 2019";
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

	// recursively iterate to get the best sequence
	bool iterate;

	// localMC repack variable
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

	// alternate identities
	//vector<string> ids;
	string ids;
	int numPositions;

	// command line arguments
	string runNumber;
	bool useIMM1;
	string insertPdb;
	
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

string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}

void standardizeResidueList(Chain &_chain, vector<string> &_residues){
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getResidueName();
		bool oldResi = false;
		for (uint j=0; j<_residues.size(); j++){
			if (resName == _residues[j]){
				oldResi = true;
				j = _residues.size();
			}
		}
		if (!oldResi){
			_residues.push_back(resName);
			//cout << "New Resi: " << resName << endl;
		}
	}
}

string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		ps = ps + " [";
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
		for (uint i=0; i<_alternateIds.size(); i++){
			if(_alternateIds[i] != resName){
				if(_alternateIds[i] == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + _alternateIds[i];
				}
			}
			else{
				continue;
			}
		}
		ps = ps + "] ";
		counter++;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
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
				if(_alternateIds[i] != resName){
					if(_alternateIds[i] == "HIS") {
						ps = ps + " HSE";
					} else {
						ps = ps + " " + _alternateIds[i];
					}
				}
				else{
					continue;
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

vector<string> singleSeqMutants(string _sequence, string _mutations){
	string origSeq = _sequence;
	vector<string> mutantSequences;
	mutantSequences.push_back(_sequence);//so the original sequence to be compared is at the top
	for (uint i=0; i<_sequence.length(); i++){
		for (uint j=0; j<_mutations.length(); j++){
			char id = _mutations[j];
			_sequence[i] = id;
			if (_sequence != origSeq){
				cout << _sequence << endl;
				mutantSequences.push_back(_sequence);
			}
			else{
				continue;
			}
		}
		_sequence = origSeq;
	}
	return mutantSequences;
}

//alternate version of above to only do mutations on interfacial positions
vector<string> singleSeqMutants(string _sequence, string _mutations, vector<int> _varPos){
	string origSeq = _sequence;
	vector<string> mutantSequences;
	mutantSequences.push_back(_sequence);//so the original sequence to be compared is at the top
	for (uint i=0; i<_sequence.length(); i++){
		if (_varPos[i] == 1){
			for (uint j=0; j<_mutations.length(); j++){
				char id = _mutations[j];
				_sequence[i] = id;
				if (_sequence != origSeq){
					cout << _sequence << endl;
					mutantSequences.push_back(_sequence);
				}
				else{
					continue;
				}
			}
		}
		_sequence = origSeq;
	}
	return mutantSequences;
}

void recursiveSequenceSwitches(vector<string> _mutantSequences, string _sequence, int _id1, int _id2){
	if (_id1 == _id2-1){
		cout << _sequence << endl;
		_mutantSequences.push_back(_sequence);
		return;
	}

	for (uint j=_id1; j<_id2; j++){
		swap(_sequence[_id1], _sequence[j]);

		recursiveSequenceSwitches(_mutantSequences, _sequence, _id1+1, _id2);

		swap(_sequence[_id1], _sequence[j]);
	}
}

string seqSwitcher(string _sequence, int _res1, int _res2){
	string seq = _sequence;
	swap(seq[_res1], seq[_res2]);
	//if (seq != _sequence){
		cout << seq << endl;
	//}
	return seq;
}

void sequenceShuffler(vector<string> _mutantSequences, string _sequence, vector<int> _varPos){
	vector<int> interfaceNumbers;
	for (uint i=0; i<_sequence.length(); i++){
		if (_varPos[i] == 1){
			interfaceNumbers.push_back(i);
		}
	}

	int numAltPos = interfaceNumbers.size()-1;
	for (uint j=0; j<numAltPos; j++){
		_mutantSequences.push_back(seqSwitcher(_sequence, interfaceNumbers[j], interfaceNumbers[numAltPos]));
		if (j == numAltPos-1){
			j=-1;//get's added back up to 0 to restart the recursion at the end
			numAltPos = numAltPos-1;
		}
		if (numAltPos == 0){
			return;
		}
	}//end's up getting 28 sequences with 8 interfacial residues (which by my count is correct)
}

void resetSwitches(System &_sys, vector<string> _resNames, vector<int> _resNums){
	for (uint i=0; i<_resNums.size(); i++){
		Position &pos = _sys.getChain("A").getPosition(_resNums[i]);
		pos.setActiveIdentity(_resNames[_resNums[i]]);
	}
}

void idSwitch(System &_sys, int _id1, int _id2){
	Position &pos1 = _sys.getChain("A").getPosition(_id1);	
	Position &pos2 = _sys.getChain("A").getPosition(_id2);	
	
	string resName1 = pos1.getResidueName();
	string resName2 = pos2.getResidueName();

	if (resName1 != resName2){
		pos1.setActiveIdentity(resName2);
		pos2.setActiveIdentity(resName1);
	}
}

//adapted from some permutation code (found online, but basically the same thing I did in Java class before with tree traversal?)
void recursiveIdSwitch(System &_sys, vector<string> &_sequences, vector<double> &_energies, int _id1, int _id2){
	if (_id1 == _id2-1){
		string seq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));
		//cout << "Dead end: " << seq << endl;
		_sequences.push_back(seq);
		_energies.push_back(_sys.calcEnergy());
		return;
	}

	for (uint j=_id1; j<_id2; j++){
		idSwitch(_sys, _id1, j);

		recursiveIdSwitch(_sys, _sequences, _energies, _id1+1, _id2);

		idSwitch(_sys, _id1, j);
	}
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

Options parseOptions(int _argc, char * _argv[], Options defaults);

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Datafile: " << _op.datafile << endl << endl;

	_fout << "Options from Datafile" << endl;
	_fout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << _op.xShift << " zShift: " << _op.zShift << " Axial Rotation: " << _op.axialRotation << " Crossing Angle: " << _op.crossingAngle << endl;
	
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

	_fout << "MCCycles " << _op.MCCycles << endl;
	_fout << "MCMaxRejects " << _op.MCMaxRejects << endl;
	_fout << "MCStartTemp " << _op.MCStartTemp << endl;
	_fout << "MCEndTemp " << _op.MCEndTemp << endl;
	_fout << "MCCurve " << _op.MCCurve << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "numberOfStructuresToMCRepack " << _op.numberOfStructuresToMCRepack << endl;
	_fout << "energyCutOff " << _op.energyCutOff << endl;

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

void switchAA(CharmmSystemBuilder &_CSB, System &_sys, vector<int> &_varPos){
	for (uint k=0; k<_varPos.size(); k++){
		if (_varPos[k] == 1){
			Position &pos = _sys.getPosition(k);
			pos.setActiveIdentity(1);
			cout << "Id switched at position " << k << "!" << endl;
			//for (uint j=0; j<_ids.size(); j++){
			//	CSB.addIdentity(pos, _ids[j]);
			//}
		}
		else{
			continue;
		}
	}
}

void switchAA(CharmmSystemBuilder &_CSB, System &_sys, vector<int> &_rand, vector<string> &_AAs, bool _removeIdentities){
	if (_removeIdentities == false){
		for (uint j=0; j<30; j++){
			Position &pos = _sys.getPosition(j);
			if (_AAs[j] != "LEU"){
				pos.setActiveIdentity(_AAs[_rand[j]]);
				//cout << "Id switched at position " << k << "!" << endl;
				//for (uint j=0; j<_ids.size(); j++){
				//	CSB.addIdentity(pos, _ids[j]);
				//}
			}
			else{
				continue;
			}
		}
	}
	else{//TODO: change this tomorrow to make it work
		for (uint j=0; j<30; j++){
			Position &pos = _sys.getPosition(j);
			for (uint k=0; k<4; k++){
				pos.setActiveIdentity(k);
				if (pos.getCurrentIdentity().toString() != _AAs[_rand[j]]){
					pos.removeIdentity(pos.getCurrentIdentity().toString());
				}
				else{
					uint current = k;
					if (k != 3){
						k++;
						pos.setActiveIdentity(k);
						pos.removeIdentity(pos.getCurrentIdentity().toString());
					}
					else{
						pos.setActiveIdentity(current);
					}
				}
			}
		}
	}
}//this doesn't work as of 2019-7-24 but could be nice to fix at some point for something with alternate IDs

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

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, double &_totalProt){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	sel.select("allProt, all");
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		double prot = _sys.calcEnergy("resi", "allProt")/2;
		double resi = _sys.calcEnergy("resi")/2; 
		double protEnergy = prot+resi;
		ener.push_back(protEnergy);
		_totalProt += protEnergy;
	}
	return ener;
}
				
void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, double _totalEnergy, double _hbond, double _vdw, int _seqNumber, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	_out << "AA      Position      Energy" << endl;//could be interesting to add rotamer number to this
	for (uint i=0; i<_seq.length(); i++){
		_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
	}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	_out << "Total Energy: " << _totalEnergy << endl;
	_out << "H-bond Energy: " << _hbond << endl;
	_out << "VDW Energy: " << _vdw << endl << endl;
}

void printEnerFile(string _seq, vector<double> &_ener, int _seqNumber, ofstream &_out){
	if (_seqNumber == 0){
		_out << "AA      Position      Energy" << endl;//could be interesting to add rotamer number to this
	}
	for (uint i=0; i<_seq.length(); i++){
		_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
	}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
}

//Goal of below is to recursively mutate the sequence with single mutants until I find the very best possible combination
//I think to start, might be best to just try this with the top sequence?
//TODO: add in MC repack
void iterateMutantSequences(Options &_opt, System &_pdb, string _sequence, double _originalEnergy, double _newEnergy, ofstream &_out){
	if (_newEnergy > _originalEnergy){
		return;
	}
	vector<string> singleMuts = singleSeqMutants(_sequence, _opt.ids);
	//TODO: decide whether I should also add in shuffles
	/******************************************************************************
	 *             === ITERATE THROUGH ALL SINGLE MUTANT SEQUENCES ===
	 ******************************************************************************/
	for (uint i=0; i<singleMuts.size(); i++){

		/******************************************************************************
		 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
		 ******************************************************************************/
		string polymerSeq = convertToPolymerSequence(singleMuts[i], 28);//TODO: make this opt.thread

		/******************************************************************************
		 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
		 ******************************************************************************/
		System sys;
		//CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
		CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
		CSB.setBuildTerm("CHARMM_ELEC", false);
		CSB.setBuildTerm("CHARMM_ANGL", false);
		CSB.setBuildTerm("CHARMM_BOND", false);
		CSB.setBuildTerm("CHARMM_DIHE", false);
		CSB.setBuildTerm("CHARMM_IMPR", false);
		CSB.setBuildTerm("CHARMM_U-BR", false);
		
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
		Eset->setTermActive("CHARMM_IMM1REF", false);
		Eset->setTermActive("CHARMM_IMM1", false);
		Eset->setTermActive("SCWRL4_HBOND", false);
		
		// Set weights
		Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
		//EsetL->setWeight("SCWRL4_HBOND", opt.weight_hbond);
		//EsetL->setWeight("SCWRL4_HBOND", 0);
		Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
		Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);
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
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		//sys.assignCoordinates(glyAPV,false);
		sys.assignCoordinates(_pdb.getAtomPointers(),false);
		sys.buildAllAtoms();
		
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
		 *          === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ENERGIES ===
		 ******************************************************************************/
		sys.setActiveRotamers(spm.getMinStates()[0]);
		double currentEnergy = sys.calcEnergy();
		cout << singleMuts[i] << endl;
		cout << currentEnergy << endl;//this worked! It outputs the proper energy now; not sure what was wrong with it before with the multiple identities, but it works properly now!
		cout << Eset->getSummary() << endl;//for troubleshooting; make sure to calculate energies before

		double previousEnergy = _originalEnergy;
		if (currentEnergy < previousEnergy){
			previousEnergy = currentEnergy;
			iterateMutantSequences(_opt, _pdb, singleMuts[i], previousEnergy, currentEnergy, _out);
		}
		else{
			iterateMutantSequences(_opt, _pdb, singleMuts[i], previousEnergy, currentEnergy, _out);
		}
		cout << "It worked!" << endl;//seems like this works as I expect, so might try running it on server tomorrow?
		//The above should continue on until the previousEnergy is no longer less
		//TODO: add in repackSidechains here to make sure best orientation
		//TODO: add in if statement for using IMM1 or not (can probably just take from vdwSequenceDesign_bcc
		//TODO: figure out why some of these results are more stable? Should I be comparing to monomer here? Should I first compare to baseline then monomer as well?
		//TODO: check through the sequence for any ids that aren't considered and add them to the ids
		//fout << sysL.getAllAtomPointers() << endl;
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

	string monoOutPdbFile  = _opt.pdbOutputDir + "/" + _opt.runNumber + "_monomer.pdb";//TODO: could change this to name?
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
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

void localMC(System &_pdb, Options &_opt, string _sequence, ofstream &_out){
	double xShift = _opt.xShift;
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	
	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequence(_sequence, 28);

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
	Eset->setTermActive("CHARMM_IMM1REF", false);
	Eset->setTermActive("CHARMM_IMM1", false);
	Eset->setTermActive("SCWRL4_HBOND", false);
	
	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	//EsetL->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	//EsetL->setWeight("SCWRL4_HBOND", 0);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	_out << "Energy Weights: " << endl;
	_out << "VDW weight: " << Eset->getWeight("CHARMM_VDW") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
	_out << endl;
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	//BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	//
	//beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
	//beb.readParameters(opt.baselineFile);
	////beb.printParameters();// Ensures they are being read properly
	//
	//beb.buildInteractions();	

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

	xShiftTransformation(apvChainA, apvChainB, axisA, axisB, _opt.xShift, trans);//this fixes the axis problem I was having, allowing MC to complete
	
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
	/******************************************************************************
	 *      === LOCAL BACKBONE MONTE CARLO REPACKS FOR TOP STRUCTURES ===
	 ******************************************************************************/
	_out << "====================================" << endl;
	_out << "Performing Local Monte Carlo Repacks" << endl;
	_out << "====================================" << endl;
	_out << endl;

	// Optimize Initial Starting Position
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(true);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.calculateEnergies();

	repackSideChains(spm, _opt.greedyCycles);

	/******************************************************************************
	 *              === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	map<string,double> monomerEnergyByTerm;
	double monomerEnergy;
	
	_out << "Monomer calculations..." << endl;
	monomerEnergy = computeMonomerEnergy(sys, trans, _opt, helicalAxis, RNG, monomerEnergyByTerm, _out, _opt.greedyCycles, _opt.MCCycles, _opt.MCMaxRejects, _opt.useIMM1);

	sys.setActiveRotamers(spm.getMinStates()[0]);

	double currentEnergy = spm.getMinBound()[0];
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	
	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	double bestEnergy = currentEnergy;
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);	

	if (_opt.MCCycles > 0) {
		//MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);
		MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
		MCMngr.setEner(bestEnergy);
		
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
			repackSideChains(spm, _opt.greedyCycles);
	
			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				_out << "state rejected   energy: " << currentEnergy << endl;
			}
			else {
				bestEnergy = currentEnergy;
				sys.saveAltCoor("savedBestState");
				helicalAxis.saveAltCoor("BestAxis");

				xShift = xShift + deltaXShift;
				crossingAngle = crossingAngle + deltaCrossingAngle;
				axialRotation = axialRotation + deltaAxialRotation;
				zShift = zShift + deltaZShift;

				_out << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;

	sys.applySavedCoor("savedBestState");

	double finalEnergy = sys.calcEnergy()-monomerEnergy;
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

	time(&startTime);
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	ofstream pout;
	ofstream rout;
	ofstream nout;
	ofstream fout;
	ofstream eout;
	
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		eout << "Unable to make directory" << endl;
		//exit(0);
	}
	string poutName = dir + "/" + date + "_options.out";
	string routName = dir + "/" + date + "_ResidueEnergies.out";
	string noutName = dir + "/" + date + "_Sequences.out";
	string foutName = dir + "/" + date + "_generalOut.out";
	string errfile = dir + "/" + date + "_test.err";
	
	pout.open(poutName.c_str());
	rout.open(routName.c_str());
	nout.open(noutName.c_str());
	fout.open(noutName.c_str());
	eout.open(errfile.c_str());

	cout << dir << endl;
	//printOptions(opt, pout);

	pout << date << endl;
	rout << date << endl;
	nout << date << endl;
	fout << date << endl;
	eout << date << endl;

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
	//readGeometryFile(opt.helixGeoFile, fileVec);

	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	PDBWriter writer;
	PDBWriter writer2;
	//PDBWriter writer3;
	//PDBWriter writer4;
	//PDBWriter writer5;
	////PDBWriter writer6;
	
	writer.open(dir + "/" + date + "_polySeq.pdb");
	writer2.open(dir + "/" + date + "_reference.pdb");
	//writer3.open(dir + "/" + date + "_polySeq-greedy.pdb");
	//writer4.open(dir + "/" + date + "_MCOFinalStates.pdb");
	//writer5.open(dir + "/" + date + "_SCMFstate.pdb");
	//writer6.open(dir + "/BestUnbiasedMCState.pdb");
	//TODO: starting here, start a for loop to generate and write energies of multiple sequences
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	//int seqNumber = opt.seqNumber;
	//vector<string> completeSequences;
	vector<string> str;
	vector<double> doub;

	/******************************************************************************
	 *             === READ PDB AND PREPARE NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polyLeu = generatePolyLeu("L", 4);

	string poly = convertToPolymerSequence(polyLeu, 28);//opt.thread must be the same
	//string poly = generateMultiIDPolymerSequence(polyLeu, opt.thread, opt.ids);
	//PolymerSequence PS(poly);
	
	System pdb;
	pdb.readPdb("/exports/home/gloiseau/Generated_PDBs/sequenceDesign/11_02_2019/SCMFstate_1.pdb", true);//TODO: make this an option
	//pdb.readPdb(opt.insertPdb, true);

	pdb.buildAtoms();
	Chain &chain1 = pdb.getChain(0);

	string seq1 = convertPolymerSeqToOneLetterSeq(chain1);
	cout << "Sequence: " << seq1 << endl;

	/******************************************************************************
	 *                   === IDENTIFY INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	vector<int> pos;
	vector<double> dist;
	vector<vector<vector<double>>> posDistVector;
	identifyInterface(pdb, pos, dist, opt.numPositions);
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

	vector<int> varPos = interface01(pdb, organizedPos);
	
	for (uint i=0; i<varPos.size(); i++){
		if (i == 21){
			cout << endl;
		}
		cout << varPos[i];
	}
	cout << endl;
	
	/******************************************************************************
	 *               === SINGLE MUTANTS SEQUENCES OF INTERFACE ===
	 ******************************************************************************/
	vector<string> singleMutants = singleSeqMutants(seq1, opt.ids, varPos);
	vector<string> mutSeq;

	cout << endl << endl;
	//recursively get all of the sequence permutations
	//recursiveSequenceSwitches(mutSeq, seq1, 0, seq1.length());
	sequenceShuffler(mutSeq, seq1, varPos);

	/******************************************************************************
	 *                       === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	double originalEnergy = 0;
	double bestEnergy = 0;
	double currentEnergy = 0;

	vector<string> betterSequences;
	vector<double> betterEnergies;
	//add in anything else that I feel might be helpful for the below code
	//
	//TODO: add in compute monomer energy for each of these single mutants, then if the overall energy is less, use that monomer energy in the localMC
	//nevermind already using baselines and anything better than that is what I'll use to evaluate if I should do monomer later
	
	/******************************************************************************
	 *             === ITERATE THROUGH ALL SINGLE MUTANT SEQUENCES ===
	 ******************************************************************************/
//	for (uint i=0; i<singleMutants.size(); i++){
//
//		/******************************************************************************
//		 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
//		 ******************************************************************************/
//		string polymerSeq = convertToPolymerSequence(singleMutants[i], 28);
//
//		/******************************************************************************
//		 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
//		 ******************************************************************************/
//		System sys;
//		//CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
//		CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
//		CSB.setBuildTerm("CHARMM_ELEC", false);
//		CSB.setBuildTerm("CHARMM_ANGL", false);
//		CSB.setBuildTerm("CHARMM_BOND", false);
//		CSB.setBuildTerm("CHARMM_DIHE", false);
//		CSB.setBuildTerm("CHARMM_IMPR", false);
//		CSB.setBuildTerm("CHARMM_U-BR", false);
//		
//		CSB.setSolvent("MEMBRANE");
//		CSB.setIMM1Params(15, 10);
//		
//		if(!CSB.buildSystem(polymerSeq)) {
//			cerr << "Unable to build system from " << polymerSeq << endl;
//			exit(0);
//		} else {
//			//fout << "CharmmSystem built for sequence" << endl;
//		}
//	
//		Chain & chainA = sys.getChain("A");
//		Chain & chainB = sys.getChain("B");
//			
//		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
//		sysRot.defineRotamerSamplingLevels();
//	
//		/******************************************************************************
//		 *                     === INITIAL VARIABLE SET UP ===
//		 ******************************************************************************/
//		EnergySet* Eset = sys.getEnergySet();
//		// Set all terms active, besides Charmm-Elec
//		Eset->setAllTermsActive();
//		Eset->setTermActive("CHARMM_ELEC", false);
//		Eset->setTermActive("CHARMM_ANGL", false);
//		Eset->setTermActive("CHARMM_BOND", false);
//		Eset->setTermActive("CHARMM_DIHE", false);
//		Eset->setTermActive("CHARMM_IMPR", false);
//		Eset->setTermActive("CHARMM_U-BR", false);
//		Eset->setTermActive("CHARMM_IMM1REF", false);
//		Eset->setTermActive("CHARMM_IMM1", false);
//		Eset->setTermActive("SCWRL4_HBOND", false);
//		
//		// Set weights
//		Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
//		//EsetL->setWeight("SCWRL4_HBOND", opt.weight_hbond);
//		//EsetL->setWeight("SCWRL4_HBOND", 0);
//		Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
//		Eset->setWeight("CHARMM_IMM1", opt.weight_solv);
//		fout << "Energy Weights: " << endl;
//		fout << "VDW weight: " << Eset->getWeight("CHARMM_VDW") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
//		fout << endl;
//	
//		/******************************************************************************
//		 *                     === ADD IN BASELINE ENERGIES ===
//		 ******************************************************************************/
//		BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
//		
//		beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
//		beb.readParameters(opt.baselineFile);
//		//beb.printParameters();// Ensures they are being read properly
//		
//		beb.buildInteractions();	
//
//		/******************************************************************************
//		 *                     === COPY BACKBONE COORDINATES ===
//		 ******************************************************************************/
//		sys.wipeAllCoordinates();
//		//sys.assignCoordinates(glyAPV,false);
//		sys.assignCoordinates(pdb.getAtomPointers(),false);
//		sys.buildAllAtoms();
//		
//		//writer.write(sys.getAtomPointers(), true, false, true);
//		//writer2.write(pdb.getAtomPointers(), true, false, true);
//	
//		/******************************************************************************
//		 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
//		 ******************************************************************************/
//		//Random Number Generator/
//		RandomNumberGenerator RNG;
//		RNG.setSeed(opt.seed); 
//		//RNG.setTimeBasedSeed();
//
//		CSB.updateNonBonded();
//		sys.buildAllAtoms();
//		
//		loadRotamers(sys, sysRot, opt.SL);
//		
//		SelfPairManager spm;
//		spm.seed(RNG.getSeed());
//		spm.setSystem(&sys);
//		spm.setVerbose(true);
//		spm.getMinStates()[0];
//		spm.updateWeights();
//		spm.setOnTheFly(true);
//		spm.calculateEnergies();
//		
//		repackSideChains(spm, 10);
//		
//		/******************************************************************************
//		 *          === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ENERGIES ===
//		 ******************************************************************************/
//		sys.setActiveRotamers(spm.getMinStates()[0]);
//		currentEnergy = sys.calcEnergy();
//		cout << singleMutants[i] << endl;
//		cout << currentEnergy << endl;//this worked! It outputs the proper energy now; not sure what was wrong with it before with the multiple identities, but it works properly now!
//		cout << Eset->getSummary() << endl;//for troubleshooting; make sure to calculate energies before
//
//		if (i == 0){
//			originalEnergy = currentEnergy;
//		}
//		else if (currentEnergy < originalEnergy){
//			betterSequences.push_back(singleMutants[i]);
//			betterEnergies.push_back(currentEnergy);
//		}
//		//double totalEnergy = sys.calcEnergy();
//		//double vdw = Eset->getTermEnergy("CHARMM_VDW");
//		//double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
//
//		//TODO: add in repackSidechains here to make sure best orientation
//		//TODO: add in if statement for using IMM1 or not (can probably just take from vdwSequenceDesign_bcc
//		//TODO: figure out why some of these results are more stable? Should I be comparing to monomer here? Should I first compare to baseline then monomer as well?
//		//TODO: check through the sequence for any ids that aren't considered and add them to the ids
//		//fout << sysL.getAllAtomPointers() << endl;
//	}
//	for (uint i=0; i<betterSequences.size(); i++){
//		cout << betterSequences[i] << ": " << betterEnergies[i] << endl;
//		if (opt.iterate){
//			iterateMutantSequences(opt, pdb, betterSequences[i], betterEnergies[i], betterEnergies[i], fout);
//		}
//		//I think this should work now? recursively?
//		//TODO: save the sequences that are best from these better sequences and then run shuffle on those as well
//		//I also should eventually ask Alessandro about this; good from a best sequence standpoint, but bad from a standardization standpoint...?
//	}
//	for (uint i=0; i<betterSequences.size(); i++){
//		localMC(pdb, opt, betterSequences[i], fout);
//	}

	//BEGIN TESTING
	/******************************************************************************
	 *            === CONVERT MUTANT SEQUENCE TO POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polymerSeq = convertToPolymerSequence(singleMutants[0], 28);//maybe change this code to do up to 10 of the sequences?

	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	//CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
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
		
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

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
	Eset->setTermActive("CHARMM_IMM1REF", false);
	Eset->setTermActive("CHARMM_IMM1", false);
	Eset->setTermActive("SCWRL4_HBOND", false);
	
	// Set weights
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	//EsetL->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	//EsetL->setWeight("SCWRL4_HBOND", 0);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);
	fout << "Energy Weights: " << endl;
	fout << "VDW weight: " << Eset->getWeight("CHARMM_VDW") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;//fixed the problem of getting improper weights, but this doesn't help me get the energies to work...
	fout << endl;

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	BaselineEnergyBuilder beb; //(sysL, "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/par_baseline.txt");
	
	beb.setSystem(sys);//had to have this before readParameters to get it to work! So it works now
	beb.readParameters(opt.baselineFile);
	//beb.printParameters();// Ensures they are being read properly
	
	beb.buildInteractions();	

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.wipeAllCoordinates();
	//sys.assignCoordinates(glyAPV,false);
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();
	
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//Random Number Generator/
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 
	//RNG.setTimeBasedSeed();

	CSB.updateNonBonded();
	sys.buildAllAtoms();
	
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
	 *          === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ENERGIES ===
	 ******************************************************************************/
	sys.setActiveRotamers(spm.getMinStates()[0]);
	currentEnergy = sys.calcEnergy();
	cout << singleMutants[0] << endl;
	cout << currentEnergy << endl;//this worked! It outputs the proper energy now; not sure what was wrong with it before with the multiple identities, but it works properly now!
	cout << Eset->getSummary() << endl;//for troubleshooting; make sure to calculate energies before

	originalEnergy = currentEnergy;
	
	cout << singleMutants[0] << ": " << currentEnergy << endl;

	localMC(pdb, opt, singleMutants[0], fout);
	//END TESTING

	//double totalEnergy = sys.calcEnergy();
	//double vdw = Eset->getTermEnergy("CHARMM_VDW");
	//double hbond = Eset->getTermEnergy("SCWRL4_HBOND");

	//using this to check how many better sequences there are
	//The below does a recursive search for the best energies

	/******************************************************************************
	 *                     === MAKING SINGLE MUTANTS ===
	 ******************************************************************************/
	cout << "Making Single Mutants..." << endl;

	//TODO: Need to make a way to go through all combinations of mutants; I think this is pretty much done (except for making sure the calculate energy changes)
	
	//for (uint i=0; i<sys.getChain("A").positionSize(); i++){
	//	Position &posA = sys.getChain("A").getPosition(i);
	//	string resName = posA.getResidueName();
	//	int origID = 0;
	//	for (uint j=0; j<opt.ids.size(); j++){
	//		//Is there a way to make the order of identities the same order as what is given? I think so, but I'll have to check; if so, this should make it a little easier
	//		//need IDchecker to get the original ID
	//		posA.setActiveIdentity(j);
	//		//cout << sys.getAtomPointers() << endl;
	//		seq = convertPolymerSeqToOneLetterSeq(chainA);
	//		if (seq != origSeq){
	//			cout << seq << endl;
	//			cout << sys.calcEnergy() << endl;//and save these in a vector to be printed later next to the sequence
	//			//extract what the sequence would be
	//		}
	//		else{
	//			continue;
	//		}
	//	}
	//	posA.setActiveIdentity(resName);//another way would be to save the original resName, then just iterate through the vector until I get that (which wouldn't need proper order?)
	//}

	//need to reset identity back to what it originally was
	//as far as I can tell from early testing, the above code is doing what I need it to do; gonna have to figure out how to build each properly again (I think just buildAllAtoms)
	
	/******************************************************************************
	 *                     === SHUFFLE SEQUENCE ===
	 ******************************************************************************/
	cout << "Shuffling..." << endl;
	//seq = convertPolymerSeqToOneLetterSeq(chainA);
	//cout << seq << endl;

	//Position &pos1 = sys.getChain("A").getPosition(0);
	//Position &pos2 = sys.getChain("A").getPosition(1);
	//Position &pos3 = sys.getChain("A").getPosition(2);
	//Position &pos4 = sys.getChain("A").getPosition(3);
	//pos2.setActiveIdentity("GLY");
	//pos3.setActiveIdentity("MET");
	//pos4.setActiveIdentity("SER");
	//vector<string> resNames;
	//vector<int> resNums;

	//for (uint i=0; i<sys.getChain("A").positionSize(); i++){
	//	Position &pos = sys.getChain("A").getPosition(i);
	//	resNames.push_back(pos.getResidueName());
	//}

	//vector<string> sequences;
	//vector<double> energies;

	//seq = convertPolymerSeqToOneLetterSeq(chainA);
	//sequences.push_back(seq);
	//energies.push_back(sys.calcEnergy());

	//idSwitch(sys, sequences, energies, 1, 2);
	//resetSwitches(sys, 0, resName1, resName2, resName3);
	//idSwitch(sys, sequences, energies, 0, 1);
	//idSwitch(sys, sequences, energies, 1, 2);
	//resetSwitches(sys, 0, resName1, resName2, resName3);
	//idSwitch(sys, sequences, energies, 0, 2);
	//idSwitch(sys, sequences, energies, 1, 2);
	//the above sequence of events works for 3AAs; how do I get it to work for 4?

	//TODO: Recursion! Need to make a way to recursively go through the entire sequence and all of the permutations of the sequence
	//figure out if it matters whether or not I reset the sequence or not, or if I can get through it all without doing that

	//recursiveIdSwitch(sys, sequences,energies, 0, sys.getChain("A").positionSize()); 

	//for (uint i=0; i<sequences.size(); i++){
	//	cout << sequences[i] << ": " << energies[i] << endl;
	//}
	//
	/******************************************************************************
	 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
	 ******************************************************************************/
	//while (seqAccept < seqNumber){
	//	cout << "Number Accepted: " << seqAccept << endl;
	//	cout << "Number Discarded: " << seqDiscard << endl;
	//	a = seqAccept;
	//	
	//	/******************************************************************************
	//	 *                     === RANDOMIZE SEQUENCES ===
	//	 ******************************************************************************/
	//	string seq = randomAASequence(weightedAAs, opt.backboneLength, RNG);
	//	cout << seq << endl;
	//	string polySeq = convertToPolymerSequence(seq, 1);
	//	PolymerSequence PS(polySeq);
	//	cout << PS << endl;

	//	//TODO: Have to add in a way to load the pdb (done! From generalCBContactMapper)
	//	/******************************************************************************
	//	 *                     === READ IN PDB FILE ===
	//	 ******************************************************************************/
	//	PDBReader reader("/exports/home/gloiseau/5x9x.pdb");
	//	
	//	PDBWriter writer;
	//	writer.open("/exports/home/gloiseau/Generated_PDBs/5x9x_copy.pdb");
	//	//string poutName = dir + "/" + opt.datafile + ".out";
	//	
	//	//pout.open("/exports/home/gloiseau/Generated_PDBs/5x9x.out");

	//	/******************************************************************************
	//	 *                     === DECLARE SYSTEM ===
	//	 ******************************************************************************/
	//	System pdb;
	//	pdb.readPdb("/exports/home/gloiseau/5x9x.pdb", true);
	//	
	//	System sys;
	//	//CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	//	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	//	CSB.setBuildTerm("CHARMM_ELEC", false);
	//	CSB.setBuildTerm("CHARMM_ANGL", false);
	//	CSB.setBuildTerm("CHARMM_BOND", false);
	//	CSB.setBuildTerm("CHARMM_DIHE", false);
	//	CSB.setBuildTerm("CHARMM_IMPR", false);
	//	CSB.setBuildTerm("CHARMM_U-BR", false);
	//	
	//	CSB.setSolvent("MEMBRANE");
	//	CSB.setIMM1Params(15, 10);
	//	
	//	if(!CSB.buildSystemFromPDB(pdb.getAtomPointers())) {
	//		cerr << "Unable to build system from " << polySeq << endl;
	//		exit(0);
	//	} else {
	//		//fout << "CharmmSystem built for sequence" << endl;
	//	}
	//	
	//	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	//	sysRot.defineRotamerSamplingLevels();
	//	
	//	// Add hydrogen bond term
	//	HydrogenBondBuilder hb(sys, opt.hBondFile);
	//	hb.buildInteractions(30);
	//	
	//	/******************************************************************************
	//	 *                     === INITIAL VARIABLE SET UP ===
	//	 ******************************************************************************/
	//	EnergySet* Eset = sys.getEnergySet();
	//	// Set all terms active, besides Charmm-Elec
	//	Eset->setAllTermsActive();
	//	Eset->setTermActive("CHARMM_ELEC", false);
	//	Eset->setTermActive("CHARMM_ANGL", false);
	//	Eset->setTermActive("CHARMM_BOND", false);
	//	Eset->setTermActive("CHARMM_DIHE", false);
	//	Eset->setTermActive("CHARMM_IMPR", false);
	//	Eset->setTermActive("CHARMM_U-BR", false);
	//	
	//	// Set weights
	//	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	//	//Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	//	//Eset->setWeight("CHARMM_IMM1REF", 1);
	//	//Eset->setWeight("CHARMM_IMM1", 1);
	//	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	//	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);
	//	
	//	/******************************************************************************
	//	 *                     === COPY BACKBONE COORDINATES ===
	//	 ******************************************************************************/
	//	sys.wipeAllCoordinates();
	//	sys.assignCoordinates(glyAPV,false);
	//	sys.buildAllAtoms();
	//	
	//	writer.write(sys.getAtomPointers(), true, false, true);
	//	
	//	//TODO: Here is where I have to add in my sequence stuff

	//	/******************************************************************************
	//	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	//	 ******************************************************************************/
	//	CSB.updateNonBonded();
	//	sys.buildAllAtoms();
	//	
	//	loadRotamers(sys, sysRot, opt.SL);
	//	
	//	SelfPairManager spm;
	//	spm.seed(RNG.getSeed());
	//	spm.setSystem(&sys);
	//	spm.setVerbose(true);
	//	spm.getMinStates()[0];
	//	spm.updateWeights();
	//	spm.setOnTheFly(true);
	//	spm.calculateEnergies();
	//	
	//	repackSideChains(spm, 10);
	//	
	//	/******************************************************************************
	//	 *                  === SET SYSTEM TO BEST SPM ROTAMERS ===
	//	 ******************************************************************************/
	//	sys.setActiveRotamers(spm.getMinStates()[0]);
	//	sys.calcEnergy();
	//	double totalEnergy = sys.calcEnergy();
	//	double vdw = Eset->getTermEnergy("CHARMM_VDW");
	//	double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
	//	
	//	if (totalEnergy < 0 && vdw < 0){
	//		/******************************************************************************
	//		 *               === CALCULATE ENERGIES FOR EACH POSITION ===
	//		 ******************************************************************************/
	//		double totalProt = 0;
	//		vector<double> AAener = calcBaselineEnergies(sys, seq.length(), totalProt);
	//		writer4.write(sys.getAtomPointers(), true, false, true);
	//	
	//		totalEnergy = sys.calcEnergy();
	//		cout << "Total Energy: " << totalEnergy << endl;
	//		cout << sys.getEnergySummary() << endl;
	//		
	//		int finalEnergy = totalEnergy;
	//		int finalTotal = totalProt;
	//	
	//		cout << "Total - ener: " << finalEnergy-finalTotal << endl;//to prove that the energies are the same since it was giving me some very small number that wasn't 0
	//		
	//		/******************************************************************************
	//		 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
	//		 ******************************************************************************/
	//		printSeqFile(PS, seq, AAener, totalEnergy, hbond, vdw, a, nout);
	//		printEnerFile(seq, AAener, a, rout);
	//		seqAccept++;
	//		sys.reset();
	//		AAener = doub;
	//	}
	//	else{
	//		if (totalEnergy > 0){
	//			cout << "Total Energy " << totalEnergy << " is too high; go to next random sequence." << endl;
	//		}
	//		if (vdw > 0){
	//			cout << "VDW Energy " << vdw << " is too high; go to next random sequence." << endl;
	//		}
	//		seqDiscard++;
	//		sys.reset();
	//	}
	//}
	//if (seqDiscard > 0){
	//	cout << "Total sequences not accepted: " << seqDiscard << endl;
	//	pout << "Total sequences not accepted: " << seqDiscard << endl;
	//}
	writer.close();
	writer2.close();
	//writer3.close();
	//writer4.close();
	//writer5.close();
	pout.close();
	nout.close();
	rout.close();
	fout.close();
	eout.close();
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
	opt.allowed.push_back("baselineFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("solvFile");

	// iterate
	opt.allowed.push_back("iterate");

	// alternate ids and weights
	opt.allowed.push_back("ids");
	opt.allowed.push_back("numPositions");

	//Command Line Arguments
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("useIMM1");
	opt.allowed.push_back("insertPdb");
	
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
			exit(1);
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

	opt.insertPdb = OP.getString("insertPdb");
	if (OP.fail()) {
		opt.warningMessages += "insertPdb not specified, using first one I used: /exports/home/gloiseau/Generated_PDBs/sequenceDesign/11_02_2019/SCMFstate_2.pdb\n";
		opt.warningFlag = true;
		opt.insertPdb = "/exports/home/gloiseau/Generated_PDBs/sequenceDesign/11_02_2019/SCMFstate_2.pdb";
	}

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

	opt.transform = OP.getBool("transform");
	if (OP.fail()) {
		opt.warningMessages += "transform not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.transform = false;
	}
	//original structure parameters
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified using 0.0\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.0;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified using 0.0\n";
		opt.warningFlag = true;
		opt.axialRotation = 0.0;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified using 25.0\n";
		opt.warningFlag = true;
		opt.crossingAngle = 25.0;
	}
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified using 10.0\n";
		opt.warningFlag = true;
		opt.xShift = 10.0;
	}
	
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

	opt.ids = OP.getString("ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	
	opt.numPositions = OP.getInt("numPositions");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify number of positions, defaulting to 8\n";//This number must be the same as from vdwSequenceDesign_bcc
		opt.warningFlag = true;
		opt.numPositions = 8;
	}
	
	opt.iterate = OP.getBool("iterate");
	if (OP.fail()) {
		opt.warningMessages += "iterate not specified defaulting to false\n";
		opt.warningFlag = true;
		opt.iterate = false;
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
