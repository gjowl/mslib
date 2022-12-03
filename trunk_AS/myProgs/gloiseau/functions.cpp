#include <sstream>
#include <iterator>
#include <unistd.h>
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv ENV;

// General functions that many of my programs use
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

// converts a polymer sequence to a string for one chain (for homodimer sequences)
string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getCurrentIdentity().getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}


/***********************************
 *functions from designFunctions
 ***********************************/
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

string convertVectorUintToString(vector<uint> _inputVector){
	string outputString = "";
	for (uint i=0; i<_inputVector.size(); i++){
		outputString += MslTools::intToString(_inputVector[i]);
	}
	return outputString;
}

// get the positions that will be linked on the interface (will have same AA identity and rotamer for self consistent mean field)
vector<uint> getLinkedPositions(vector<uint> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel){
	vector<uint> positionsToLink;
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
			positionsToLink.push_back(1);
		} else {
			positionsToLink.push_back(0);
		}
	}
	return positionsToLink;
}

// define the rotamer level for each position in the backbone
vector<uint> convertStringToVectorUint(string _inputString){
	vector<uint> outputVec;
	for (uint i=0; i<_inputString.size(); i++){
		stringstream ss;
		ss << _inputString[i];
		uint stringToInt = MslTools::toUnsignedInt(ss.str());
		outputVec.push_back(stringToInt);
	}
	return outputVec;
}

// get a backbone sequence with an alanine cap at the beginning and end as an option
string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) {
	// initial start of sequence
	string str = "";
	//2021-09-21: add in an alanine cap to allow for more variable positions at the leucine region
	for (uint i=0; i<_length-3; i++){
		if (i<3){
			if (_useAlaCap == true){
				str = str + "A";
			} else {
				str = str + _backboneAA;
			}
		} else {
			str = str + _backboneAA;
		}
	}
	// Adds in the LILI at the end of the sequence which is necessary for our TOXCAT plasmids
	if (_useAlaCap == true){
		str = str + "AAA";
	} else {
		str = str + "ILI";
	}
	return str;
}

// generate string for backbone sequence
string generateString(string _backbone, int _length) {
	string str = "";
	for (uint i=0; i<_length; i++){
		str = str + _backbone;
	}
	return str;
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
		} else if (pos < startPos+3 || pos > endPos-5){
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