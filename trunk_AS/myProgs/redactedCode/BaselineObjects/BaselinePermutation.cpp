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


#include "BaselinePermutation.h"
#include "System.h"
#include "Residue.h"
#include "Reader.h"

using namespace MSL;
using namespace std;

const string BaselinePermutation::typeName = "BASELINE_PERMUTATION";

BaselinePermutation::BaselinePermutation() {
	//setup(NULL,NULL,NULL,NULL,NULL,NULL,0.0);
}

//TODO: add an option here to change the energy to something that changes with linear or exponential (or add to my design code)

BaselinePermutation::BaselinePermutation(double _energy){
	setup(sys, type, factor);
}

BaselinePermutation::BaselinePermutation(System *_sys, string _type, double _factor){
	setup(_sys, _type, _factor);
}

BaselinePermutation::BaselinePermutation(BaselinePermutation & _interaction) {
	copy(_interaction);
}

BaselinePermutation::~BaselinePermutation() {
}


void BaselinePermutation::setup(System *_sys, string _type, double _factor){
	sys = _sys;
	type = _type;
	factor = _factor;
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	pAtoms.push_back(&sys->getAtom(0));//TODO: figure out how to make this an atom pointer; this is the problem it needs an atom to latch onto as active
}

void BaselinePermutation::copy(BaselinePermutation & _interaction) {
	sys = _interaction.getSysPointers();
	type = _interaction.getType();
	factor = _interaction.getFactor();
	pMap = _interaction.getMap();
	penFracs = _interaction.getPenaltyFractions();
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	pAtoms.push_back(&sys->getAtom(0));
}

void BaselinePermutation::printParameters() {
	if(pAtoms.size() > 0 && pAtoms[0] && params.size() > 0) {
		cout << " ResName " << pAtoms[0]->getResidueName() << endl;
		cout << " atomName " << pAtoms[0]->getName() << endl;
		cout << " energy " << params[0] << endl;
	} else {
		cout << "No atoms defined in interaction" << endl;
	}
}

double BaselinePermutation::getEnergy(std::vector<double> *_paramDerivatives) {
	//return getEnergy(params[0],_paramDerivatives);
	return getEnergy();
}
double BaselinePermutation::getEnergy(double _dummy, std::vector<double> *_paramDerivatives) {
	// The gradient is zero
	if(_paramDerivatives) {
		_paramDerivatives->resize(pAtoms.size() * 3,0.0);
	}
	return _dummy;

}
vector<double> BaselinePermutation::getEnergyGrad(){
	std::vector<double> foo(pAtoms.size() * 3,0.0);
	return foo;
}


void BaselinePermutation::setCurrentSeq(){
	Chain &chain = sys->getChain("A");
	vector<Position*>& positions = chain.getPositions();
	for (vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
		Residue &res = (*p)->getCurrentIdentity();
		string resName = res.getResidueName();
		pSeq.push_back(resName);
	}
}

double BaselinePermutation::calcMultiplicationFactor(double _count, double _factor, string type){
	//TODO: add in the default for factor and all other types
	double f = 0;
	if (type == ""){
		f = _factor*_count;
		_factor = f;
	}
	//if (type == "LINEAR"){
	//	f = factor*count;
	//}
	//if (type == "QUADRATIC"){
	//	f = factor*count;
	//}
	//if (type == "EXPONENTIAL"){
	//	f = factor*count;
	//}
	return _factor;
}

map<string,int> BaselinePermutation::getAACountMap(){
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

void BaselinePermutation::buildInteraction(){
	pSeq.clear();
	energy = 0;
	setCurrentSeq();
	double numerator = 1;
	double denominator = 1;
	map<string,int> countMap = getAACountMap();

	//TODO: find AA in count map, get counts of AA, then figure out what to do with count (exponent?) and existence in membrane value
	map<string,int>::iterator it;
	for (it=countMap.begin(); it != countMap.end(); it++){
		int count = it->second;
		double memProb = pMap.at(it->first);
	
		numerator = numerator*memProb;
		denominator = denominator*(pow(memProb, count));
	}
	energy = numerator/denominator;
	cout << energy << endl;
	if (params.size() == 0){
		params.push_back(energy);
	}
	else{
		params[0] = energy;
	}
}

//The below function is what is called by the EnergySet to get the energies for each interaction; I'm just going to get one energy from the whole term when it calls it, based around sequence
double BaselinePermutation::getEnergy() {
	//sys->getEnergySet()->eraseTerm("BASELINE_COMPOSITION");
	//sys->getEnergySet()->addInteraction(new BaselinePermutation(energy));
	buildInteraction();
	return(params[0]);
}//TODO: this function will not work for what I need with composition: I need to know how many AAs there are and calculate accordingly
// I think NOW that I know how both MSL and Rosetta work, I actually need to rethink this and implement the code here. If I implement it in the design code, I can't access things that need to be changed here. This getEnergy() should basically just calculate the total energy of the interaction rather than add it up, using functions from there. I think that's the best way to do it
// So instead I should take in the system and AAs and penalties in the design code to save them as private variables here, then call them with this getEnergy() into another function that actually calculates the energy by checking the sequence, adding them up, and return that sum instead...question is, if I'm referencing the system, is the system ACTUALLY undergoing changes within the SPM, or something else? This could be something I run into, but I think it's likely just referencing the system


