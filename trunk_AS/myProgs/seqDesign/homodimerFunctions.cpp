/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/22
 * @Email:  gjowl04@gmail.com
 * @Filename: homodimerFunctions.cpp
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/22
 */


#include <sstream>
#include <iterator>
#include <unistd.h>
#include "homodimerFunctions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

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

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA,
	AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis,
	double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
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

	//_out << "Geometry: " << _opt.xShift << "\t" << _opt.crossingAngle << "\tDensity: " << angleDistDensity << endl;
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
 *define interface and rotamer sampling
 ***********************************/
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

//Identify which positions are found at the identified interface
//Example: Sequence:  LLLLIGLLIGLLIGLLLL
//         Interface: 000011001100110000
// Positions at interface are 1 and non-interfacial are 0
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

//makes the best state applicable for unlinked positions by duplicating the rotamer at each interfacial position on the opposite chain
void unlinkBestState(Options &_opt, vector<uint> &_bestState, vector<int> _rotamerSampling, int _backboneLength){
	vector<int> linkedPositions = getLinked(_rotamerSampling, _backboneLength, _opt.interfaceLevel, _opt.sasaRepackLevel.size()-1);

	for (uint i=_backboneLength; i<_backboneLength*2; i++){
		vector<int>::iterator itr;
		itr = find(linkedPositions.begin(), linkedPositions.end(), i-_backboneLength);
		if (itr != linkedPositions.end()){
			vector<uint>::iterator itPos;
			itPos = _bestState.begin()+i;
			_bestState.insert(itPos, _bestState[i-_backboneLength]);
		}
	}
}

/***********************************
 *stateMC helper functions
 ***********************************/
void randomPointMutation(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _variablePositions, vector<string> &_ids){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _variablePositions.size()-1);
	int pos = _variablePositions[rand];

	// Get a random integer to pick through the AA identities
	int randIdNum = _RNG.getRandomInt(0, _ids.size()-1);
	string posId = _sys.getPosition(pos).getPositionId();
	string randId;
	randId = _ids[randIdNum];

	string res = _sys.getPosition(pos).getResidueName();

	_sys.setActiveIdentity(posId, randId);
}

//
void randomPointMutationUnlinked(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, vector<uint> _variablePositions, vector<string> &_ids){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _variablePositions.size()-1);
	int posA = _variablePositions[rand];
	int posB = posA+_opt.backboneLength;

	// Get the random position from the system
	Position &randPosA = _sys.getPosition(posA);
	Position &randPosB = _sys.getPosition(posB);
	string posIdA = randPosA.getPositionId();
	string posIdB = randPosB.getPositionId();

	// Get a random integer to pick through the AA identities
	int randIdNum = _RNG.getRandomInt(0, _ids.size()-1);
	string randId;
	randId = _ids[randIdNum];

	string res = _sys.getPosition(posA).getResidueName();

	//if (_opt.designHomodimer){
	_sys.setActiveIdentity(posIdA, randId);
	_sys.setActiveIdentity(posIdB, randId);
	//TODO: add to hetero code
	//} else {
	//	_sys.setActiveIdentity(posIdA, randId);
	//	// get a second AA identity to choose for the other helix
	//	randIdNum = _RNG.getRandomInt(0, _ids.size()-1);
	//	randId = _ids[randIdNum];
	//	// set the position on helix B as another random AA identity
	//	_sys.setActiveIdentity(posIdB, randId);
	//}
}

//Checks through a vector of sequences to see if the new sequence has already been found
void checkSequenceVector(string &_newSeq, vector<string> &_seqs){
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
bool sameSequenceChecker(string &_newSeq, double &_newEnergy, vector<uint> &_state,
vector<pair<double,string>> &_enerSeqPair, vector<pair<double,vector<uint>>> &_energyStateVec){
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

void saveSequence(Options &_opt, vector<pair<double,string>> &_energyVector,
vector<pair<double,vector<uint>>> &_energyStateVec, string _sequence, vector<uint> _state, double _energy){
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

double getMapValueFromKey(map<string,double> &_map, string _key){
	double value = 0;
	map<string,double>::const_iterator itr = _map.find(_key);
	if (itr == _map.end()){
		cerr << _key << " not in map" << endl;
		exit(0);
	} else {
		value = itr->second;
	}
	return value;
}

void saveSequence(Options &_opt, RandomNumberGenerator &_RNG, map<vector<uint>,
map<string,double>> &_stateEnergyMap, vector<pair<double,string>> &_energyVector,
vector<pair<double,vector<uint>>> &_energyStateVec, string _sequence, vector<uint> _state, double _energy, ofstream &_out){
	bool sameSeq = sameSequenceChecker(_sequence, _energy, _state, _energyVector, _energyStateVec);
	if (sameSeq == false){
		if (_energyVector.size() < _opt.numStatesToSave){
			_energyVector.push_back(make_pair(_energy, _sequence));
			_energyStateVec.push_back(make_pair(_energy, _state));
		} else {
			//choose random saved state and replace if better
			int seqNumber = _RNG.getRandomInt(0, _energyStateVec.size()-1);
			int vecPos = seqNumber;
			vector<uint> lowStateProb = _energyStateVec[seqNumber].second;
			string lowSeqProb = _energyVector[seqNumber].second;
			map<string,double> energyMap = _stateEnergyMap.at(lowStateProb);
			double prevSEProb = getMapValueFromKey(energyMap, "SequenceProbability");
			string prevSeq = _energyVector[vecPos].second;

			////Get lowest prob from stateMap
			//string lowSeqProb;
			//vector<uint> lowStateProb;
			//double prevSEProb = 1;
			//uint vecPos = 0;
			//for (uint i=0; i<_energyStateVec.size(); i++){
			//	vector<uint> state = _energyStateVec[i].second;
			//	string seq = _energyVector[i].second;
			//	map<string,double> energyMap = _stateEnergyMap.at(state);
			//	double prob = getMapValueFromKey(energyMap, "SequenceProbability");
			//	if (prob < prevSEProb){
			//		prevSEProb = prob;
			//		lowSeqProb = seq;
			//		lowStateProb = state;
			//		vecPos = i;
			//	}
			//}

			//Compare lowest energy with seqProb to the new sequence
			map<string,double> prevSeqEnergyMap = _stateEnergyMap.at(lowStateProb);
			map<string,double> currSeqEnergyMap = _stateEnergyMap.at(_state);
			double prevEnergy = getMapValueFromKey(prevSeqEnergyMap, "Dimer");
			double prevVDW = getMapValueFromKey(prevSeqEnergyMap, "VDWDimer");
			double currEnergy = getMapValueFromKey(currSeqEnergyMap, "Dimer");
			double currVDW = getMapValueFromKey(currSeqEnergyMap, "VDWDimer");
			double currSEProb = getMapValueFromKey(currSeqEnergyMap, "SequenceProbability");

			double totSEProb = prevSEProb+currSEProb;
			double prevSeqProp = prevSEProb/totSEProb;
			double currSeqProp = currSEProb/totSEProb;

			//Convert the probability of each sequence to an entropy term that can be added to the original energy to directly compare two sequences
			double prevEntropy = -log(prevSeqProp)*0.592*_opt.weight_seqEntropy;
			double currEntropy = -log(currSeqProp)*0.592*_opt.weight_seqEntropy;

			//Calculate energy total for best sequence vs current sequence
			//The below includes the baseline energy, which is an estimate of monomer energy
			double prevEnergyTotal = prevEnergy+prevEntropy;
			double currEnergyTotal = currEnergy+currEntropy;

			if (currSeqProp > prevSeqProp+0.9 && prevVDW > currVDW){//Approximately the strength of 5 strong vdW interactions greater (looked to me like Trp had some strong interactions, but since always losing out, loosened the accept parameters)
				cout.precision(6);
				cout << "Switching Sequence: " << prevSeq << " to " << _sequence << ". VDW: " << prevVDW << " > " << currVDW << " & Sequence Probability comparison > 0.9: " << currSeqProp << " > " << prevSeqProp << endl;
				_out << "Switching Sequence: " << prevSeq << " to " << _sequence << ". VDW: " << prevVDW << " > " << currVDW << " & Sequence Probability comparison > 0.9: " << currSeqProp << " > " << prevSeqProp << endl;
				_energyVector[vecPos].first = _energy;
				_energyVector[vecPos].second = _sequence;
				_energyStateVec[vecPos].first = _energy;
				_energyStateVec[vecPos].second = _state;
			}
		}
	}
}

bool convertStateMapToSequenceMap(System &_sys, vector<pair<double,vector<uint>>> &_energyStateVec,
map<vector<uint>, map<string,double>> &_stateEnergyMap, map<string, map<string,double>> &_sequenceEnergyMap,
vector<pair<string,vector<uint>>> &_sequenceStatePair, ofstream &_out){
	//convert stateEnergyMap to sequenceEnergyMap
	map<vector<uint>,map<string,double>>::iterator itr;
	bool nonClashing = false;
	Chain & chain = _sys.getChain("A");
	// Setup sequence and state pair vector
	cout << endl << "Accepted Design Sequences" << endl;
	_out << endl << "Accepted Design Sequences" << endl;
	cout << "Sequence               Energy " << endl;
	_out << "Sequence               Energy " << endl;
	for (uint i=0; i<_energyStateVec.size(); i++){
		vector<uint> state = _energyStateVec[i].second;
		_sys.setActiveRotamers(state);
		string currSeq = convertPolymerSeqToOneLetterSeq(chain);
		_sequenceEnergyMap[currSeq] = _stateEnergyMap.at(state);
		_sequenceStatePair.push_back(make_pair(currSeq, state));
		cout << currSeq << ": " << _sequenceEnergyMap[currSeq]["Dimer"] << endl;
		_out << currSeq << ": " << _sequenceEnergyMap[currSeq]["Dimer"] << endl;
		if (_sequenceEnergyMap[currSeq]["VDWDimer"] < 0){
			nonClashing = true;
		}
	}
	return nonClashing;
}

void addSequencesToVector(vector<pair<double,string>> &_energyVector, vector<string> &_allSeqs){
	for (uint i=0; i<_energyVector.size(); i++){
		//cout << i << ": " << _energyVector[i].first << " " << _energyVector[i].second << endl;
		checkSequenceVector(_energyVector[i].second, _allSeqs);
	}
}

void getEnergiesForStartingSequence(Options &_opt, SelfPairManager &_spm, string _startSequence,
vector<unsigned int> &_stateVector, map<string, map<string, double>> &_sequenceEnergyMap, map<string, double> &_entropyMap){
	if (_opt.useBaseline){
		double baseline = _spm.getStateEnergy(_stateVector, "BASELINE")+_spm.getStateEnergy(_stateVector, "BASELINE_PAIR");
		_sequenceEnergyMap[_startSequence]["Dimer"] = _spm.getStateEnergy(_stateVector)-baseline;
		_sequenceEnergyMap[_startSequence]["Baseline"] = baseline;
	} else {
		_sequenceEnergyMap[_startSequence]["Dimer"] = _spm.getStateEnergy(_stateVector);
	}
	_sequenceEnergyMap[_startSequence]["VDWDimer"] = _spm.getStateEnergy(_stateVector, "CHARMM_VDW");
	_sequenceEnergyMap[_startSequence]["HBONDDimer"] = _spm.getStateEnergy(_stateVector, "SCWRL4_HBOND");
	_sequenceEnergyMap[_startSequence]["IMM1Dimer"] = _spm.getStateEnergy(_stateVector, "CHARMM_IMM1REF")+_spm.getStateEnergy(_stateVector, "CHARMM_IMM1");

	map<string,int> internalSeqCountMap;
	map<string,int> seqCountMap;
	double internalNumberOfPermutations;
	double numberOfPermutations;
	sequenceEntropySetup(_startSequence, seqCountMap, numberOfPermutations, _opt.backboneLength);
	internalAASequenceEntropySetup(_startSequence, internalSeqCountMap, internalNumberOfPermutations, _opt.backboneLength);
	double internalSEProb = calculateSequenceProbability(internalSeqCountMap, _entropyMap, internalNumberOfPermutations);
	double SEProb = calculateSequenceProbability(seqCountMap, _entropyMap, numberOfPermutations);

	_sequenceEnergyMap[_startSequence]["SequenceProbability"] = internalSEProb;
}

void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1){
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
		//_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
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

//TODO: add below functions to other header files
/***********************************
 *energy builders
 ***********************************/
void buildBaselineIMM1Interactions(System &_sys, map<string, double> &_selfMap){
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
					ESet->addInteraction(new BaselineIMM1Interaction(*a,ener));
				}
				catch (const out_of_range& e){
					continue;
				}
			}
		}
	}
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

/***********************************
 * SequenceEntropyFunctions
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

void interfaceAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<uint> _interfacialPositionsList){
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



double getInterfaceSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap, vector<uint> _interfacialPositionsList){
	map<string,int> AACountMap;
	double numberOfPermutations;
	interfaceAASequenceEntropySetup(_sequence, AACountMap, numberOfPermutations, _interfacialPositionsList);
	double seqProb = calculateSequenceProbability(AACountMap, _entropyMap, numberOfPermutations);
	return seqProb;
}

double getSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap){
	map<string,int> AACountMap;
	double numberOfPermutations;
	internalAASequenceEntropySetup(_sequence, AACountMap, numberOfPermutations, _opt.backboneLength);
	double seqProb = calculateSequenceProbability(AACountMap, _entropyMap, numberOfPermutations);
	return seqProb;
}

void calculateInterfaceSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq,
map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_prevEntropy,
double &_currEntropy, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal, vector<uint> _interfacePositionsList){

	_prevSEProb = getInterfaceSequenceEntropyProbability(_opt, _prevSeq, _entropyMap, _interfacePositionsList);
	_currSEProb = getInterfaceSequenceEntropyProbability(_opt, _currSeq, _entropyMap, _interfacePositionsList);

	//Calculate the probability of each sequence compared to the other
	//On 10_26_21: after changing to the 2021_10_25_seqEntropies.txt file with more accurate membrane composition from Liu et al. 2002, I started seeing a lot of GxxxG sequences. I also checked my runs from before this date, and there were qute a bit of GxxxG sequences as well. Uncomment the below if you would like to try running using the internal probabilities, which should give a bit more flexibiilty to choosing AAs with lower likelilhoods of being in the membrane
	double totSEProb = _prevSEProb+_currSEProb;
	double prevSeqProp = _prevSEProb/totSEProb;
	double currSeqProp = _currSEProb/totSEProb;

	//Convert the probability of each sequence to an entropy term that can be added to the original energy to directly compare two sequences
	_prevEntropy = -log(prevSeqProp)*0.592*_opt.weight_seqEntropy;//Multiplies the energy (log proportion times RT in kcal/mol) by the sequence entropy weight (weight of 1 gives entropy same weight as other terms)
	_currEntropy = -log(currSeqProp)*0.592*_opt.weight_seqEntropy;

	//Calculate energy total for best sequence vs current sequence
	//The below includes the baseline energy, which is an estimate of monomer energy
	_bestEnergyTotal = _bestEnergy+_prevEntropy;
	_currEnergyTotal = _currEnergy+_currEntropy;

	//Output the terms if verbose
	if (_opt.verbose){
		cout << "Prev Prob:    " << _prevSEProb << endl;
		cout << "New Prob:     " << _currSEProb << endl;
		cout << "Prev Seq Proportion: " << prevSeqProp << endl;
		cout << "New Seq Proportion:  " << currSeqProp << endl;
		cout << "PrevEner =    " << _prevEntropy << endl;
		cout << "NewEner =     " << _currEntropy << endl;
		cout << "Diff =        " << (_prevEntropy-_currEntropy) << endl;
		cout << "Best Energy: " << _bestEnergyTotal << endl;
		cout << "New Energy: " << _currEnergyTotal << endl;
	}
}

void calculateSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq,
map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_prevEntropy,
double &_currEntropy, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal){

	_prevSEProb = getSequenceEntropyProbability(_opt, _prevSeq, _entropyMap);
	_currSEProb = getSequenceEntropyProbability(_opt, _currSeq, _entropyMap);

	//Calculate the probability of each sequence compared to the other
	//double totSEProb = _prevSEProb+_currSEProb;
	//On 10_26_21: after changing to the 2021_10_25_seqEntropies.txt file with more accurate membrane composition from Liu et al. 2002, I started seeing a lot of GxxxG sequences. I also checked my runs from before this date, and there were qute a bit of GxxxG sequences as well. Uncomment the below if you would like to try running using the internal probabilities, which should give a bit more flexibiilty to choosing AAs with lower likelilhoods of being in the membrane
	double totSEProb = _prevSEProb+_currSEProb;
	double prevSeqProp = _prevSEProb/totSEProb;
	double currSeqProp = _currSEProb/totSEProb;

	//Convert the probability of each sequence to an entropy term that can be added to the original energy to directly compare two sequences
	_prevEntropy = -log(prevSeqProp)*0.592*_opt.weight_seqEntropy;//Multiplies the energy (log proportion times RT in kcal/mol) by the sequence entropy weight (weight of 1 gives entropy same weight as other terms)
	_currEntropy = -log(currSeqProp)*0.592*_opt.weight_seqEntropy;

	//Calculate energy total for best sequence vs current sequence
	//The below includes the baseline energy, which is an estimate of monomer energy
	_bestEnergyTotal = _bestEnergy+_prevEntropy;
	_currEnergyTotal = _currEnergy+_currEntropy;

	//Output the terms if verbose
	if (_opt.verbose){
		cout << "Prev Prob:    " << _prevSEProb << endl;
		cout << "New Prob:     " << _currSEProb << endl;
		cout << "Prev Seq Proportion: " << prevSeqProp << endl;
		cout << "New Seq Proportion:  " << currSeqProp << endl;
		cout << "PrevEner =    " << _prevEntropy << endl;
		cout << "NewEner =     " << _currEntropy << endl;
		cout << "Diff =        " << (_prevEntropy-_currEntropy) << endl;
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

void internalAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, int _seqLength){
	vector<string> seqVector;
	//TODO: make these numbers variable fr ohow long I want each of the termini cutoffs
	for (uint i=3; i<_seqLength-4; i++){
		stringstream tmp;
		tmp << _seq[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, _seqLength-7);
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
	//cout << "pSeq = " << seqProb << endl;
	seqProb = seqProb*_numberOfPermutations;
	//cout << "comb = " << _numberOfPermutations << endl;
	//cout << "pRes = " << seqProb << endl;
	return seqProb;
}
/////

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
 * MonteCarlo Functions
 ***********************************/
void stateMCLinked(System &_sys, SelfPairManager &_spm, Options &_opt, PolymerSequence &_PS,
map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
vector<unsigned int> &_bestState, vector<string> &_seqs, vector<string> &_allSeqs,
vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<uint> &_interfacialPositionsList,
vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	//initialize baselineAAComposition energy map (based on Rosetta baselineAAComposition to prevent unlikely sequences by adding energy (ex. if more than 2SER in sequence, add 1000 energy score for each additional PHE)

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
	double bestEnergyTotal = 0;
	double currEnergyTotal = 0;
	double currStateSEProb = 0;
	double prevStateSEProb = 0;
	double prevStateEntropy = 0;
	double currStateEntropy = 0;
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
	lout << "PrevSequence\tCurrSequence\tTotal\tDimer\tBaseline\tVDW\tHBOND\tEnergyDifferencew/prevSeq\tPrevEnergy\tCurrEnergy\tPrevEntropy\tCurrEntropy\tCurrentTemp" << endl;

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
			outputEnergiesByTerm(_spm, prevStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
			stateMCEnergies["Dimer"] = bestEnergy;

			stateEnergyMap[prevStateVec] = stateMCEnergies;
			totEnergy = stateMCEnergies["EnergyBeforeLocalMC"];
			sequences[totEnergy] = prevStateSeq;
			double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
			saveSequence(_opt, energyVector, energyStateVec, prevStateSeq, prevStateVec, prevVDW);
		}
		// Reset the energy map to save energies from new state after changing the rotamer
		stateMCEnergies.clear();
		randomPointMutation(_sys, _opt, _RNG, _interfacialPositionsList, ids);

		// Set a mask and run a greedy to get the best state for the sequence
		//_sys.setActiveRotamers(currStateVec);
		vector<vector<bool>> mask = getActiveMask(_sys);
		_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
		currStateVec = _spm.getMinStates()[0];

		// Get the sequence for the random state
		_sys.setActiveRotamers(currStateVec);
		string currStateSeq = convertPolymerSeqToOneLetterSeq(chain);

		// Compute dimer energy
		outputEnergiesByTerm(_spm, currStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
		double currStateEnergy = _spm.getStateEnergy(currStateVec);
		stateMCEnergies["Dimer"] = currStateEnergy;

		// Convert the energy term (which actually saves the probability of the sequence in the whole _system)
		// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)
		calculateSequenceEntropy(_opt, prevStateSeq, currStateSeq, _sequenceEntropyMap,
			prevStateSEProb, currStateSEProb, prevStateEntropy, currStateEntropy, bestEnergy,
			currStateEnergy, bestEnergyTotal, currEnergyTotal);
		MC.setEner(bestEnergyTotal);

		// MC accept and reject conditions
		double currVDW = _spm.getStateEnergy(currStateVec, "CHARMM_VDW");
		double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
		double currHBOND = _spm.getStateEnergy(currStateVec, "SCWRL4_HBOND");
		double prevHBOND = _spm.getStateEnergy(prevStateVec, "SCWRL4_HBOND");
		if (!MC.accept(currEnergyTotal)){
			_sys.setActiveRotamers(prevStateVec);
			currStateVec = prevStateVec;

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			//TODO: make these a separate function or put in comments  for them
			bestEnergy = currStateEnergy;
			MC.setEner(currEnergyTotal);
			prevStateSEProb = currStateSEProb;
			string prevStateSeq1 = prevStateSeq;
			prevStateSeq = currStateSeq;
			prevStateVec = currStateVec;
			_sys.setActiveRotamers(currStateVec);

			//addEnergiesToStateMap(stateMCEnergies, currStateVec, ){

			//}
			double EnergyBeforeLocalMC = currStateEnergy-(_spm.getStateEnergy(currStateVec, "BASELINE")+_spm.getStateEnergy(currStateVec, "BASELINE_PAIR"));
			stateMCEnergies["EnergyBeforeLocalMC"] = EnergyBeforeLocalMC;
			stateMCEnergies["Dimer"] = EnergyBeforeLocalMC;
			stateMCEnergies["Baseline"] = _spm.getStateEnergy(currStateVec, "BASELINE")+_spm.getStateEnergy(currStateVec, "BASELINE_PAIR");
			stateMCEnergies["EnergyBeforeLocalMCw/seqEntropy"] = bestEnergyTotal-currEnergyTotal;
			stateMCEnergies["SequenceProbability"] = currStateSEProb;
			stateEnergyMap[currStateVec] = stateMCEnergies;

			saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currStateEnergy);
			double prevEnergy = bestEnergyTotal;

			if (_opt.energyLandscape){
				map<string,double> energyMap = stateEnergyMap.at(currStateVec);
				lout << prevStateSeq1 << "\t" << currStateSeq << "\t";
				for (uint j=0; j<_opt.energyLandscapeTerms.size(); j++){
					lout << energyMap.at(_opt.energyLandscapeTerms[j]) << "\t";
				}
				lout << bestEnergyTotal << "\t" << currEnergyTotal << "\t" << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currStateSeq << "; PrevE=  " << bestEnergy << " : CurrE= " << currStateEnergy << "; PrevVDW: " << prevVDW << " : CurrVDW: " << currVDW << endl;
			}
		cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to
		if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
			MC.reset(3649, 3649, 500, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		}
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();

	addSequencesToVector(energyVector, _allSeqs);
	convertStateMapToSequenceMap(_sys, energyStateVec, stateEnergyMap, _sequenceEnergyMap, _sequenceStatePair, _sout);
	getDimerSasaScores(_sys, _sequenceStatePair, _sequenceEnergyMap);

	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_sout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_seqs = _allSeqs;
}



//Make it so that this will get all the info I need instead of having to run more code later
void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS,
map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
vector<unsigned int> &_bestState, vector<string> &_seqs, vector<string> &_allSeqs,
vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<uint> &_allInterfacialPositionsList,
vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
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
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
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
	if (_opt.useBaseline){
		//addBaselineToSelfPairManager();//TODO: was going to make this a function but I feel like it already exists in spm, so going to wait when I have time to look that up
		//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
		map<string, double> selfMap = readSingleParameters(_opt.selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_opt.pairEnergyFile);
		buildSelfInteractions(sys, selfMap);
		buildPairInteractions(sys, pairMap);
		sys.calcEnergy();
	}

	/******************************************************************************
	 *              === LOAD ROTAMERS AND CHOOSE TO LINK INTERFACE ===
	 ******************************************************************************/
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);

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
	if (_opt.verbose){
		cout << "Calculating self and pair energy terms..." << endl;
	}
	spm.calculateEnergies();

	/******************************************************************************
	 *                   === SETUP FOR STATE MONTE CARLO ===
	 ******************************************************************************/
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);

	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, 10);
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
	double bestEnergyTotal = 0;
	double currEnergyTotal = 0;
	double currStateSEProb = 0;
	double prevStateSEProb = 0;
	double prevStateEntropy = 0;
	double currStateEntropy = 0;
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
	lout << "PrevSequence\tCurrSequence\tTotal\tDimer\tBaseline\tVDW\tHBOND\tEnergyDifferencew/prevSeq\tPrevEnergy\tCurrEnergy\tPrevEntropy\tCurrEntropy\tCurrentTemp" << endl;
	if (_opt.verbose){
		cout << "***STARTING GEOMETRY:***" << endl;
		cout << "xShift: " << _opt.xShift << endl;
		cout << "crossingAngle: " << _opt.crossingAngle << endl;
		cout << "axialRotation: " << _opt.axialRotation << endl;
		cout << "zShift: " << _opt.zShift << endl << endl;
		cout << "Number of MCCycles: " << _opt.MCCycles << endl;
	}

	/******************************************************************************
	 *                      === BEGIN STATE MONTE CARLO ===
	 ******************************************************************************/
	//TODO: I think if I ever want to run this on multiple cores, I should set it up to run those here, then just save x number from each core
	cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	_sout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		// Get energy term and probability for the first sequence (these terms can then get replaced by future sequences that are accepted by MC)
		if (cycleCounter == 0){
			outputEnergiesByTerm(spm, prevStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
			stateMCEnergies["Dimer"] = bestEnergy;
			double prevSeqProb = getSequenceEntropyProbability(_opt, prevStateSeq, _sequenceEntropyMap);
			stateMCEnergies["SequenceProbability"] = prevSeqProb;

			stateEnergyMap[prevStateVec] = stateMCEnergies;
			totEnergy = stateMCEnergies["EnergyBeforeLocalMC"];

			sequences[totEnergy] = prevStateSeq;
			double prevVDW = spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
			//saveSequence(_opt, energyVector, energyStateVec, prevStateSeq, prevStateVec, bestEnergy);
		}
		// Reset the energy map to save energies from new state after changing the rotamer
		stateMCEnergies.clear();
		cout << 1 << endl;
		randomPointMutationUnlinked(sys, _opt, _RNG, _interfacialPositionsList, ids);
		cout << 2 << endl;

		// Set a mask and run a greedy to get the best state for the sequence
		//sys.setActiveRotamers(currStateVec);
		vector<vector<bool>> mask = getActiveMask(sys);
		spm.runGreedyOptimizer(_opt.greedyCycles, mask);
		currStateVec = spm.getMinStates()[0];

		// Get the sequence for the random state
		sys.setActiveRotamers(currStateVec);
		string currStateSeq = convertPolymerSeqToOneLetterSeq(chain);

		// Compute dimer energy
		outputEnergiesByTerm(spm, currStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
		double currStateEnergy = spm.getStateEnergy(currStateVec);

		// Convert the energy term (which actually saves the probability of the sequence in the whole system)
		// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)

		//TODO: I just realized that for heterodimers, I may need to completely remake some of these functions as with the below only taking the sequence of one helix; I may make a hetero and homo functions list?
		calculateInterfaceSequenceEntropy(_opt, prevStateSeq, currStateSeq, _sequenceEntropyMap, prevStateSEProb, currStateSEProb, prevStateEntropy, currStateEntropy, bestEnergy, currStateEnergy, bestEnergyTotal, currEnergyTotal, _allInterfacialPositionsList);
		MC.setEner(bestEnergyTotal);

		// MC accept and reject conditions
		double currVDW = spm.getStateEnergy(currStateVec, "CHARMM_VDW");
		double prevVDW = spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
		double currHBOND = spm.getStateEnergy(currStateVec, "SCWRL4_HBOND");
		double prevHBOND = spm.getStateEnergy(prevStateVec, "SCWRL4_HBOND");
		if (!MC.accept(currEnergyTotal)){
			sys.setActiveRotamers(prevStateVec);
			currStateVec = prevStateVec;

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			//TODO: make these a separate function or put in comments  for them
			bestEnergy = currStateEnergy;
			MC.setEner(currEnergyTotal);
			prevStateSEProb = currStateSEProb;
			string prevStateSeq1 = prevStateSeq;
			prevStateSeq = currStateSeq;
			prevStateVec = currStateVec;
			sys.setActiveRotamers(currStateVec);

			outputEnergiesByTerm(spm, currStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
			double EnergyBeforeLocalMC = currStateEnergy-(spm.getStateEnergy(currStateVec, "BASELINE")+spm.getStateEnergy(currStateVec, "BASELINE_PAIR"));
			stateMCEnergies["EnergyBeforeLocalMC"] = EnergyBeforeLocalMC;
			stateMCEnergies["Dimer"] = EnergyBeforeLocalMC;
			stateMCEnergies["Baseline"] = spm.getStateEnergy(currStateVec, "BASELINE")+spm.getStateEnergy(currStateVec, "BASELINE_PAIR");
			stateMCEnergies["EnergyBeforeLocalMCw/seqEntropy"] = bestEnergyTotal-currEnergyTotal;
			stateMCEnergies["SequenceProbability"] = currStateSEProb;
			stateEnergyMap[currStateVec] = stateMCEnergies;

			//TODO: change this so I just save energies in the same place and easily can get vdw, hbond, etc. for each saved sequence
			//saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currStateEnergy);
			if (_opt.weight_seqEntropy == 0){
				saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currVDW);
			} else {
				saveSequence(_opt, _RNG, stateEnergyMap, energyVector, energyStateVec, currStateSeq, currStateVec, currStateEnergy, _sout);
			}
			double prevEnergy = bestEnergyTotal;

			if (_opt.energyLandscape){
				map<string,double> energyMap = stateEnergyMap.at(currStateVec);
				lout << prevStateSeq1 << "\t" << currStateSeq << "\t";
				for (uint j=0; j<_opt.energyLandscapeTerms.size(); j++){
					lout << energyMap.at(_opt.energyLandscapeTerms[j]) << "\t";
				}
				lout << bestEnergyTotal << "\t" << currEnergyTotal << "\t" << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currStateSeq << "; PrevE=  " << prevEnergy << " : CurrE= " << currEnergyTotal << "; PrevVDW: " << prevVDW << " : CurrVDW: " << currVDW << "EnergyDifference" << bestEnergyTotal-currEnergyTotal << "; CurrTemp: " << MC.getCurrentT() << endl;
			}
		cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to
		if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
			MC.reset(3649, 3649, 500, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		}
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();

	addSequencesToVector(energyVector, _allSeqs);
	convertStateMapToSequenceMap(sys, energyStateVec, stateEnergyMap, _sequenceEnergyMap, _sequenceStatePair, _sout);
	getDimerSasaScores(sys, _sequenceStatePair, _sequenceEnergyMap);

	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_sout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_seqs = _allSeqs;
}


/***********************************
 *Other functions
 ***********************************/
void getTotalEnergyAndWritePdbs(System &_sys, Options &_opt, map<string, map<string,double>> &_sequenceEnergyMap,
string _sequence, vector<uint> _stateVec, vector<int> _rotamerSampling, RandomNumberGenerator &_RNG, int _seqNumber,
PDBWriter &_writer, ofstream &_sout, ofstream &_err){
	double dimerEnergy = _sequenceEnergyMap[_sequence]["Dimer"];
	double finalEnergy = dimerEnergy-_sequenceEnergyMap[_sequence]["Monomer"];
	_sout << "Dimer Energy: " << dimerEnergy << endl;
	_sout << "Final Energy w/ IMM1 = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;
	cout << "Dimer Energy: " << dimerEnergy << endl;
	cout << "Final Energy = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;
	_sequenceEnergyMap[_sequence]["Dimer"] = dimerEnergy;
	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;
	_sequenceEnergyMap[_sequence]["SequenceNumber"] = _seqNumber;

	//Setup directory for individual design
	string repackDir = _opt.pdbOutputDir + "/" + _sequence;
	string cmd = "mkdir -p " + repackDir;
	if (system(cmd.c_str())){
		_err << "Unable to make directory" << endl;
		exit(0);
	}

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

	loadRotamers(sys, sysRot, _opt, _rotamerSampling);
	CSB.updateNonBonded(10,12,50);

	sys.setActiveRotamers(_stateVec);
	PDBWriter designWriter;
	designWriter.open(repackDir + "/" + _sequence + ".pdb");
	designWriter.write(sys.getAtomPointers(), true, false, true);
	designWriter.close();
	_writer.write(sys.getAtomPointers(), true, false, true);
	saveEnergyDifference(_opt, _sequenceEnergyMap, _sequence);

	CRDWriter writer;
	writer.open(repackDir + "/" + _sequence + ".crd");
	writer.write(_sys.getAtomPointers(), false);
	writer.close();
}

void getDimerSasaScores(System &_sys, vector<pair<string,vector<uint>>> &_sequenceStatePair, map<string, map<string,double>> &_sequenceEnergyMap){
	for (uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		vector<uint> state = _sequenceStatePair[i].second;

		_sys.setActiveRotamers(state);

		//Setup SasaCalculator to calculate the monomer SASA
		SasaCalculator sasa(_sys.getAtomPointers());
		sasa.calcSasa();
		double dimerSasa = sasa.getTotalSasa();

		_sequenceEnergyMap[sequence]["DimerSasa"] = dimerSasa;
	}
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

//TODO: Use this as a template for making this code more clean: when I don't have the name in the config file, it still works...how can I make it like that for everything? Should I just have one function that will get all of the values that need to be returned from a much larger map that may have more than that? And then if the value isn't there, call out illegal value that needs to be added in?
void getSasaForStartingSequence(System &_sys, string _sequence, vector<uint> _state, map<string, map<string,double>> &_sequenceEnergyMap){
	_sys.setActiveRotamers(_state);

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	_sequenceEnergyMap[_sequence]["DimerSasa"] = dimerSasa;
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

//Other functions
void saveEnergyDifference(Options _opt, map<string,map<string,double>> &_sequenceEnergyMap, string _sequence){
	map<string,double> &energyMap = _sequenceEnergyMap[_sequence];
	energyMap["Baseline-Monomer"] = energyMap["Baseline"] + energyMap["Monomer"];
	energyMap["HBONDDiff"] = energyMap["HBONDDimer"] - energyMap["HBONDMonomer"];
	energyMap["VDWDiff"] = energyMap["VDWDimer"] - energyMap["VDWMonomer"];
	energyMap["IMM1Diff"] = energyMap["IMM1Dimer"] - energyMap["IMM1Monomer"];
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
