/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/12
 * @Email:  gjowl04@gmail.com
 * @Filename: design.cpp
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/22
 */

#include <sstream>
#include <iterator>
#include <unistd.h>
#include "designFunctions.h"
#include "versatileFunctions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

/***********************************
 * load rotamers 
 ***********************************/
void loadRotamersByBurial(System &_sys, SystemRotamerLoader &_sysRot, vector<string> _repackLevels, vector<uint> _rotamerSampling){
	//Repack side chains based on sasa scores
	// get backbone length
	int backboneLength = _sys.getChain("A").positionSize();
	for (uint i = 0; i < _rotamerSampling.size(); i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+backboneLength);
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _repackLevels[_rotamerSampling[i]];
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
			string posRot = _repackLevels[_rotamerSampling[i]];
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

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, bool _useInterfaceBurial, vector<string> _repackLevels, string _SL,
 vector<uint> _rotamerSampling){
	// if using the SASA to identify the interface, then load the rotamers by the SASA burial
	if (_useInterfaceBurial){
		loadRotamersByBurial(_sys, _sysRot, _repackLevels, _rotamerSampling);
	} else {
		loadRotamers(_sys, _sysRot, _SL);
	}
}

/***************************************
 *define interface and rotamer sampling
 ***************************************/
//Identify which positions are found at the identified interface
//Example: Sequence:  LLLLIGLLIGLLIGLLLL
//         Interface: 000011001100110000
// Positions at interface are 1 and non-interfacial are 0
// Convert positions to string for setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions) which uses "A,19" "B,19" format!
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<uint> &_interfacePositions, int _backboneLength){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_interfacePositions.size(); k++){
		int interfacePosition = _interfacePositions[k];
		// search through the interface positions vector to get the interface
		vector<string> tempPos;

		// get the positions at the backbone position k
		Position &posA = _sys.getPosition(interfacePosition);
		Position &posB = _sys.getPosition(interfacePosition+_backboneLength);

		// convert position to string
		string A = posA.toString();
		string B = posB.toString();

		// find the space delimiter	
		string delimiter = " ";
		size_t p = 0;
		p = A.find(delimiter);

		// add the string without anything after the delimiter to the tempPos vector
		tempPos.push_back(A.substr(0, p));
		tempPos.push_back(B.substr(0, p));
			
		// add the tempPos vector to the stringPositions vector
		stringPositions.push_back(tempPos);
	}
	return stringPositions;
}

//makes the best state applicable for unlinked positions by duplicating the rotamer at each interfacial position on the opposite chain
vector<uint> unlinkBestState(vector<uint> _bestState, vector<uint> _interfacePositions, int _backboneLength){
	int originalSize = _bestState.size();
	for (uint i=0; i<originalSize; i++){
		_bestState.push_back(_bestState[i]);
	}
	return _bestState;
}

/***********************************
 *string output functions
 ***********************************/
// get the string of an interface sequence in 00010001 format where 1 is an interface residue and 0 is a non-interface residue
string getInterfaceSequence(int _interfaceLevelLimit, string _interface, string _sequence){
	string interfaceSequence = "";
	for(string::iterator it = _interface.begin(); it != _interface.end(); it++) {
		stringstream ss;
		ss << *it;
		int pos = it-_interface.begin();
		int tmp = MslTools::toInt(ss.str());
		if (tmp > _interfaceLevelLimit-1){//interfaceLevel counts from 1 but rotamerLevel coutns from 0
			interfaceSequence = interfaceSequence + "-";
		} else {
			interfaceSequence = interfaceSequence + _sequence[pos];
		}
	}
	return interfaceSequence;
}

/***********************************
 *  calculate residue burial
 ***********************************/
// function from Samson that calculates the residue burial of a system by comparing it to the typical burial of an amino acid
std::vector<pair <int, double> > calculateResidueBurial (System &_sys) {
	/*
	  SASA reference:
	  Protein Engineering vol.15 no.8 pp.659–667, 2002
	  Quantifying the accessible surface area of protein residues in their local environment
	  Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
	*/
	map<string,double> refSasa;
	refSasa["G"] = 83.91;
	refSasa["A"] = 116.40;
	refSasa["S"] = 125.68;
	refSasa["C"] = 141.48;
	refSasa["P"] = 144.80;
	refSasa["T"] = 148.06;
	refSasa["D"] = 155.37;
	refSasa["V"] = 162.24;
	refSasa["N"] = 168.87;
	refSasa["E"] = 187.16;
	refSasa["Q"] = 189.17;
	refSasa["I"] = 189.95;
	refSasa["L"] = 197.99;
	refSasa["H"] = 198.51;
	refSasa["K"] = 207.49;
	refSasa["M"] = 210.55;
	refSasa["F"] = 223.29;
	refSasa["Y"] = 238.30;
	refSasa["R"] = 249.26;
	refSasa["W"] = 265.42;

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator b(_sys.getAtomPointers());
	//b.setProbeRadius(3.0);
	b.calcSasa();
	for (uint i = 0; i < _sys.positionSize()/2; i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdA = _sys.getPosition(i).getPositionId();
		//string posIdB = _sys.getPosition(i+_sys.positionSize()/2).getPositionId();
		string oneLetter = MslTools::getOneLetterCode(_sys.getPosition(i).getResidueName());
		double resiSasa = b.getResidueSasa(posIdA);
		double burial = resiSasa / refSasa[oneLetter];
		residueBurial.push_back(pair<int,double>(i, burial));
		//residueBurial.push_back(pair<string,double>(posIdB, burial));
		//cout << posId << "\t" << resiSasa << "\t" << burial << endl;
	}

	return residueBurial;
}

// get a vector of all interfacial positions, including the ends
vector<uint> getAllInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength){
	vector<uint> variableInterfacePositions;
	for (uint k=0; k<_backboneLength; k++){//TODO: make this not hardcoded to skip RAS
		if (_rotamerSamplingPerPosition[k] < _interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

// get a vector of interface positions that doesn't include the ends of the sequence
vector<uint> getInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength){
	vector<uint> variableInterfacePositions;
	for (uint k=3; k<_backboneLength-5; k++){//TODO: make this not hardcoded to skip RAS
		if (_rotamerSamplingPerPosition[k] < _interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

/***********************************
 *baseline energy helper functions
 ***********************************/
//Function to calculate the non-termini self energies of a chain
vector<double> calcBaselineEnergies(System &_sys, int _thread, int _backboneLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	//cout << "Self Energies" << endl;
	for (uint i=_thread+3; i<_thread+_backboneLength-5; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
		//cout << number << ": " << resi << endl;
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to calculate the non-termini pair energies of a chain
vector<double> calcPairBaselineEnergies(System &_sys, int _thread, int _backboneLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=_thread+3; i<_thread+_backboneLength-5; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i);
		sel.select(residue += num1);
		for (uint j=i+1; j<_thread+_backboneLength-5;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			} else {
				j = _thread+_backboneLength-5;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

// build in the baseline energies to the system EnergySet
void buildBaselines(System &_sys, string _selfEnergyFile, string _pairEnergyFile){
		map<string, double> selfMap = readSingleParameters(_selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_pairEnergyFile);
		buildSelfInteractions(_sys, selfMap);
		buildPairInteractions(_sys, pairMap);
}
// reads in the self energy baseline file
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

// reads in the pair energy baseline file
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

// builds the self energy interactions into energy set
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
				//cout << "Adding self energy for " << baseId << " Position: " << p-positions.begin() << endl;
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

// builds the pair energy interactions into energy set
void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	EnergySet* ESet = _sys.getEnergySet();
	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j < (*p)->identitySize(); j++){
				Residue &res1 = (*p)->getIdentity(j);
				string baseId1 = res1.getResidueName();
				if (p-positions.begin() < 3){
					baseId1 = baseId1.append("-ACE");
				}
				//Changed this on 11_24_2021 for the designFiles/2021_11_22_IMM1Self and Pair baselines
				if (p-positions.begin() > positions.size()-5){//
					baseId1 = baseId1.append("-CT2");
				}
				for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					if (d <= 10){
						//cout << "Distance: " << d << endl;
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							if (p2-positions.begin() < 3){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-5){
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

/***********************************
 *calculate energies
 ***********************************/
void getSasaForStartingSequence(System &_sys, string _sequence, vector<uint> _state, map<string, map<string,double>> &_sequenceEnergyMap){
	_sys.setActiveRotamers(_state);

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	_sequenceEnergyMap[_sequence]["DimerSasa"] = dimerSasa;
}

void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, ofstream &_out){
	// output this information about the SelfPairManager run below
	// vector for the initial SCMF state
	vector<unsigned int> initialState = _spm.getSCMFstate();
	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = _spm.getBestSCMFBiasedMCState();
	
	// checks to see if the sequence changed between the initialState and the bestState
	_sys.setActiveRotamers(initialState);
	string initialSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));
	_sys.setActiveRotamers(bestState);
	string SCMFBestSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));

	// energy terms from the startingState
	double bestEnergy = _spm.getStateEnergy(bestState);
	double hbondEnergy = _spm.getStateEnergy(bestState, "SCWRL4_HBOND");
	double vdwEnergy = _spm.getStateEnergy(bestState, "CHARMM_VDW");

	// outputs	
	_out << "Initial Sequence:   " << initialSeq << endl;
	_out << "Best Sequence:      " << SCMFBestSeq << endl;
	_out << "Interface Sequence: " << _interfaceSeq << endl;
	_out << "Total Energy:       " << bestEnergy << endl;
	_out << "VDW:                " << vdwEnergy << endl;
	_out << "HBOND:              " << hbondEnergy << endl;
}

/***********************************
 *sequence entropy functions
 ***********************************/
// gets the sequence entropy value for the current sequence and adds it to the energy
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
	//cout << "Number Positions = " << _seqLength << endl;
	for(uint i=_seqLength; i>1; i--){
		numPermutation = numPermutation*i;
	}
	// was running the below before 2022-10-19; realized it doesn't give me the n!/(n-r)!, so not sure why I thought it was
	// as of 2022-11-3: running this again per Alessandro; it should be n!/(numAA1!*numAA2!*...)
	map<string,int>::iterator itr;
	for(itr = _seqAACounts.begin(); itr != _seqAACounts.end(); itr++){
		//cout << itr->first << ": " << itr->second << endl;
		for (uint j=itr->second; j>1; j--){
			permutationDenominator = permutationDenominator*j;
		}
	}
	//double denomFactorial = _seqLength-_seqAACounts.size();
	//cout << "Denominator Factorial = " << denomFactorial << "!" << endl;
	//cout << "Num unique AAs = " << _seqAACounts.size();
	//for (uint i=denomFactorial; i>1; i--){
	//	permutationDenominator = permutationDenominator*i;
	//}
	//cout << "Number Permutations: " << numPermutation << endl;
	//cout << "Permutation Denominator: " << permutationDenominator << endl;
	numPermutation = numPermutation/permutationDenominator;
	//cout << "Number Permutations/Denominator: " << numPermutation << endl;

	return numPermutation;
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


// output energies by term into a referenced energyMap
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

