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
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

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

//bool convertStateMapToSequenceMap(System &_sys, vector<pair<double,vector<uint>>> &_energyStateVec,
//map<vector<uint>, map<string,double>> &_stateEnergyMap, map<string, map<string,double>> &_sequenceEnergyMap,
//vector<pair<string,vector<uint>>> &_sequenceStatePair, ofstream &_out){
//	//convert stateEnergyMap to sequenceEnergyMap
//	map<vector<uint>,map<string,double>>::iterator itr;
//	bool nonClashing = false;
//	Chain & chain = _sys.getChain("A");
//	// Setup sequence and state pair vector
//	cout << endl << "Accepted Design Sequences" << endl;
//	_out << endl << "Accepted Design Sequences" << endl;
//	cout << "Sequence               Energy " << endl;
//	_out << "Sequence               Energy " << endl;
//	for (uint i=0; i<_energyStateVec.size(); i++){
//		vector<uint> state = _energyStateVec[i].second;
//		_sys.setActiveRotamers(state);
//		string currSeq = convertPolymerSeqToOneLetterSeq(chain);
//		_sequenceEnergyMap[currSeq] = _stateEnergyMap.at(state);
//		_sequenceStatePair.push_back(make_pair(currSeq, state));
//		cout << currSeq << ": " << _sequenceEnergyMap[currSeq]["Dimer"] << endl;
//		_out << currSeq << ": " << _sequenceEnergyMap[currSeq]["Dimer"] << endl;
//		if (_sequenceEnergyMap[currSeq]["VDWDimer"] < 0){
//			nonClashing = true;
//		}
//	}
//	return nonClashing;
//}

void addSequencesToVector(vector<pair<double,string>> &_energyVector, vector<string> &_allSeqs){
	for (uint i=0; i<_energyVector.size(); i++){
		//cout << i << ": " << _energyVector[i].first << " " << _energyVector[i].second << endl;
		checkSequenceVector(_energyVector[i].second, _allSeqs);
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

//Other functions
void saveEnergyDifference(Options _opt, map<string,map<string,double>> &_sequenceEnergyMap, string _sequence){
	map<string,double> &energyMap = _sequenceEnergyMap[_sequence];
	energyMap["Baseline-Monomer"] = energyMap["Baseline"] + energyMap["Monomer"];
	energyMap["HBONDDiff"] = energyMap["HBONDDimer"] - energyMap["HBONDMonomer"];
	energyMap["VDWDiff"] = energyMap["VDWDimer"] - energyMap["VDWMonomer"];
	energyMap["IMM1Diff"] = energyMap["IMM1Dimer"] - energyMap["IMM1Monomer"];
}

//TODO: add the below functions to baseline files?

