/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/


#include "BaselineSequenceEntropyNormalized.h"
#include "System.h"
#include "Residue.h"
#include "Reader.h"

using namespace MSL;
using namespace std;
//TODO: I think I should make an option to normalize my sequenceEntropies if not using all AAs
const string BaselineSequenceEntropyNormalized::typeName = "BASELINE_ENTROPY_NORM";

BaselineSequenceEntropyNormalized::BaselineSequenceEntropyNormalized() {
	//setup(NULL,NULL,NULL,NULL,NULL,NULL,0.0);
}

//TODO: add an option here to change the energy to something that changes with linear or exponential (or add to my design code)

BaselineSequenceEntropyNormalized::BaselineSequenceEntropyNormalized(double _energy){
	setup(sys, type);
}

BaselineSequenceEntropyNormalized::BaselineSequenceEntropyNormalized(System *_sys, string _type){
	setup(_sys, _type);
}

BaselineSequenceEntropyNormalized::BaselineSequenceEntropyNormalized(BaselineSequenceEntropyNormalized & _interaction) {
	copy(_interaction);
}

BaselineSequenceEntropyNormalized::~BaselineSequenceEntropyNormalized() {
}


void BaselineSequenceEntropyNormalized::setup(System *_sys, string _type){
	sys = _sys;
	type = _type;
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	pAtoms.push_back(&sys->getAtom(0));//TODO: figure out how to make this an atom pointer; this is the problem it needs an atom to latch onto as active
}

void BaselineSequenceEntropyNormalized::copy(BaselineSequenceEntropyNormalized & _interaction) {
	sys = _interaction.getSysPointers();
	type = _interaction.getType();
	pMap = _interaction.getMap();
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	pAtoms.push_back(&sys->getAtom(0));
}

void BaselineSequenceEntropyNormalized::printParameters() {
	if(pAtoms.size() > 0 && pAtoms[0] && params.size() > 0) {
		cout << " ResName " << pAtoms[0]->getResidueName() << endl;
		cout << " atomName " << pAtoms[0]->getName() << endl;
		cout << " energy " << params[0] << endl;
	} else {
		cout << "No atoms defined in interaction" << endl;
	}
}

double BaselineSequenceEntropyNormalized::getEnergy(std::vector<double> *_paramDerivatives) {
	//return getEnergy(params[0],_paramDerivatives);
	return getEnergy();
}
double BaselineSequenceEntropyNormalized::getEnergy(double _dummy, std::vector<double> *_paramDerivatives) {
	// The gradient is zero
	if(_paramDerivatives) {
		_paramDerivatives->resize(pAtoms.size() * 3,0.0);
	}
	return _dummy;

}
vector<double> BaselineSequenceEntropyNormalized::getEnergyGrad(){
	std::vector<double> foo(pAtoms.size() * 3,0.0);
	return foo;
}


void BaselineSequenceEntropyNormalized::setCurrentSeq(){
	Chain &chain = sys->getChain("A");
	vector<Position*>& positions = chain.getPositions();
	for (vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
		Residue &res = (*p)->getCurrentIdentity();
		string resName = res.getResidueName();
		pSeq.push_back(resName);
	}
}

map<string,int> BaselineSequenceEntropyNormalized::getAACountMap(){
	map<string,int> AAcounts;
	for (uint i=0; i<pSeq.size(); i++){
		try{
			if (AAcounts.count(pSeq[i]) > 0){
				AAcounts.at(pSeq[i])++;
			}
			else{
				AAcounts[pSeq[i]] = 1;
			}
		}
		catch(const out_of_range& e){
			continue;
		}
	}
	return AAcounts;
}

//Redacted after chatting with Alessandro on 2-9-2021; new version better for accounting for total combinations to compare between sequences
/*void BaselineSequenceEntropyNormalized::calcNumeratorPermutations(double &_numPermutation, int _uniqueAAs){
	for(uint i=8; i>1; i--){
		_numPermutation = _numPermutation*i;
	}
}
void BaselineSequenceEntropyNormalized::calcDenominatorPermutations(double &_denPermutation, int _AAcount){
	for(uint i=_AAcount; i>1; i--){
		_denPermutation = _denPermutation*i;
	}

}*/

void BaselineSequenceEntropyNormalized::calcNumberPermutations(double &_numPermutation, vector<int> _counts){
	//This function calculates number of permutations using following equation: n!/(n-r)! where n is number of positions and r is number of AAs	
	_numPermutation = 1;
	double permutationDenominator = 1;
	for(uint i=21; i>1; i--){
		_numPermutation = _numPermutation*i;
	}
	for (uint i=0; i<_counts.size(); i++){
		for (uint j=_counts[i]; j>1; j--){
			permutationDenominator = permutationDenominator*j;
		}
	}
	_numPermutation = _numPermutation/permutationDenominator;
}

void BaselineSequenceEntropyNormalized::calcTotalCombinations(int _uniqueAAs){
	totCombinations = 1;
	totCombinations = pow(8, _uniqueAAs);//TODO: fix this so that it takes into account number of positions and uniqueAAs...
}

void BaselineSequenceEntropyNormalized::setPolyLeuEner(double _energy){
	polyLeuEner = _energy;
}

void BaselineSequenceEntropyNormalized::setPolyLeuProb(double _prob){
	polyLeuProb = _prob;
}

void BaselineSequenceEntropyNormalized::setTotalProb(double _prob){
	totalProb = 1/_prob;//normalizes (?) to energy of polyLeu being 0	
}

void BaselineSequenceEntropyNormalized::buildInteraction(){
	pSeq.clear();
	energy = 0;
	seqProb = 1;
	setCurrentSeq();
	map<string,int> countMap = getAACountMap();
	vector<int> counts;
	int numUniqueAAs = 0;
	numPermutations = 1;

	//Find AA in count map, get counts of AA, then calculate the sequence probability by multiplying each AAs membrane probability contribution
	map<string,int>::iterator it;
	for (it=countMap.begin(); it != countMap.end(); it++){
		int count = it->second;
		double memProb = pMap.at(it->first);
		//cout << "AA: " << it->first << endl;
		//cout << "Count: " << it->second << endl;
		//cout << "Prob: " << pMap.at(it->first) << endl;
		seqProb = seqProb*(pow(memProb, count));
		numUniqueAAs++;
		counts.push_back(count);
		//calcDenominatorPermutations(denPermutation, count);
	}
	//A fix for the seqProb: only count the interfacial positions
	//seqProb = seqProb/(pow(0.2, 13));//does give me slightly larger energy differences, but still not substantial; gonna need to think of another way
	//TODO: what if I still use these values, but to recalculate another value: like a straight up energy comparison between sequences to see there differences, rather than against the entire range of sequence energies: I can take the value from a good sequence, compare it, and only if the overall energy is better will I be able to change sequence (or maybe even better, only if the population change would be significant enough between the two sequences); so something like a 2.5 for example could be a significant enough change to warrant giving a good energy, where something worse get no additional energy: leads to a push for more diverse sequences (can then do this for multiple sequences that are saved until the energies are maxed out

	//TODO: add a way to normalize values (I am currently precomputing them externally, but should just compute them when given, depending on the list of available amino acids for design 
	
	//Calculate the number of permutations (need to know number of different AAs (if only 1, then 1 permutation, if 2, then 2, if 3, then 6, etc.
	calcNumberPermutations(numPermutations, counts);
	calcTotalCombinations(numUniqueAAs);//TODO: I don't think I need this anymore because of the below
	
	//Below Redacted: used to do this way, but after chatting with Alessandro on 2-9-2021, might be better to just use totPermutations (makes easier to compare between just two sequences)
	//seqProb = seqProb*(totPermutations/totCombinations);
	
	seqProb = seqProb*numPermutations;
	//seqProb = seqProb*(log2(seqProb));//Converts the sequence probability to an entropy
	//cout << "Seq Prob: " << seqProb << endl;
	//energy = seqProb*(log2(seqProb));//Converts the sequence probability to an entropy

	//if (numPermutations == 1){
	//	setPolyLeuEner(0);
	//	setPolyLeuProb(seqProb);
	//	setTotalProb(polyLeuProb);//TODO: make sure this works/assuming if polyLeu energy is 0
	////	energy = polyLeuEner/energy*100;//TODO: need to find a way to convert energy difference to something a bit more substantial
	//}
	//else{
	//	seqProb = seqProb/totalProb;//I think this works if I assume that the energy of polyLeu is 0 (e^0 = 1)
	//	//TODO: this technically works, but the energy change is really tiny (which makes sense, because the sequence probability change is small
	//	//1. Figure out how to make changes in sequence probability matter more (likely by decreasing the totalProb somehow?
	//	//energy = polyLeuEner/energy*100;//gives the percentage difference in states for new sequence versus original sequence (the more different, the better)
	//	energy = 0.6*log(seqProb);//reverse boltzmann calculation for energy from probability of sequence
	//}
	energy = -0.592*log(seqProb);//reverse boltzmann calculation for energy from probability of sequence
	//cout << "Prob: " << seqProb << endl;
	//cout << "Energy: " << energy << endl;
	//seqProb = seqProb/totalProb;//I think this works if I assume that the energy of polyLeu is 0 (e^0 = 1)
	//energy = 0.6*log(seqProb);//reverse boltzmann calculation for energy from probability of sequence
	//TODO: need to figure out how to get either p for an energy of 0/1 (or any value?), and then calculate the sum of all energies, to use as a totalProb 
	//energy = seqProb*(totPermutations/totCombinations);
	//cout << "Poly Leu Energy: " << polyLeuEner << endl;
	//cout << "Energy: " << energy << endl;
	if (params.size() == 0){
		params.push_back(energy);
	}
	else{
		params[0] = energy;
	}
}

//The below function is what is called by the EnergySet to get the energies for each interaction; I'm just going to get one energy from the whole term when it calls it, based around sequence
double BaselineSequenceEntropyNormalized::getEnergy() {
	//sys->getEnergySet()->eraseTerm("BASELINE_COMPOSITION");
	//sys->getEnergySet()->addInteraction(new BaselineSequenceEntropyNormalized(energy));
	buildInteraction();
	return(params[0]);
}//TODO: this function will not work for what I need with composition: I need to know how many AAs there are and calculate accordingly
// I think NOW that I know how both MSL and Rosetta work, I actually need to rethink this and implement the code here. If I implement it in the design code, I can't access things that need to be changed here. This getEnergy() should basically just calculate the total energy of the interaction rather than add it up, using functions from there. I think that's the best way to do it
// So instead I should take in the system and AAs and penalties in the design code to save them as private variables here, then call them with this getEnergy() into another function that actually calculates the energy by checking the sequence, adding them up, and return that sum instead...question is, if I'm referencing the system, is the system ACTUALLY undergoing changes within the SPM, or something else? This could be something I run into, but I think it's likely just referencing the system


