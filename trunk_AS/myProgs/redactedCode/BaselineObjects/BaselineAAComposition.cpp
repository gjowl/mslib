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


#include "BaselineAAComposition.h"
#include "System.h"
#include "Residue.h"
#include "Reader.h"

using namespace MSL;
using namespace std;

const string BaselineAAComposition::typeName = "BASELINE_COMPOSITION";

BaselineAAComposition::BaselineAAComposition() {
	//setup(NULL,NULL,NULL,NULL,NULL,NULL,0.0);
}

//TODO: add an option here to change the energy to something that changes with linear or exponential (or add to my design code)

BaselineAAComposition::BaselineAAComposition(double _energy){
	setup(sys, type, factor);
}

BaselineAAComposition::BaselineAAComposition(System *_sys, string _type, double _factor){
	setup(_sys, _type, _factor);
}

BaselineAAComposition::BaselineAAComposition(BaselineAAComposition & _interaction) {
	copy(_interaction);
}

BaselineAAComposition::~BaselineAAComposition() {
}


void BaselineAAComposition::setup(System *_sys, string _type, double _factor){
	sys = _sys;
	type = _type;
	factor = _factor;
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	defaultPropertyMap();
	buildInteraction();
	pAtoms.push_back(&sys->getAtom(0));//TODO: figure out how to make this an atom pointer; this is the problem it needs an atom to latch onto as active
}

void BaselineAAComposition::copy(BaselineAAComposition & _interaction) {
	sys = _interaction.getSysPointers();
	type = _interaction.getType();
	factor = _interaction.getFactor();
	penAAs = _interaction.getPenaltyAAs();
	pMap = _interaction.getMap();
	penFracs = _interaction.getPenaltyFractions();
	propertyVecMap = _interaction.getPropertyMap();
	hasProperty = _interaction.getPropertyBooleans();
	function.clear();
	pSeq.clear();
	params.clear();
	pAtoms.clear();
	buildInteraction();
	pAtoms.push_back(&sys->getAtom(0));
}

void BaselineAAComposition::printParameters() {
	if(pAtoms.size() > 0 && pAtoms[0] && params.size() > 0) {
		cout << " ResName " << pAtoms[0]->getResidueName() << endl;
		cout << " atomName " << pAtoms[0]->getName() << endl;
		cout << " energy " << params[0] << endl;
	} else {
		cout << "No atoms defined in interaction" << endl;
	}
}

double BaselineAAComposition::getEnergy(std::vector<double> *_paramDerivatives) {
	//return getEnergy(params[0],_paramDerivatives);
	return getEnergy();
}
double BaselineAAComposition::getEnergy(double _dummy, std::vector<double> *_paramDerivatives) {
	// The gradient is zero
	if(_paramDerivatives) {
		_paramDerivatives->resize(pAtoms.size() * 3,0.0);
	}
	return _dummy;

}
std::vector<double> BaselineAAComposition::getEnergyGrad(){
	std::vector<double> foo(pAtoms.size() * 3,0.0);
	return foo;
}

void BaselineAAComposition::setCurrentSeq(){
	Chain &chain = sys->getChain("A");
	vector<Position*>& positions = chain.getPositions();
	for (vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
		Residue &res = (*p)->getCurrentIdentity();
		string resName = res.getResidueName();
		pSeq.push_back(resName);
	}
}

void BaselineAAComposition::readPenaltyFile(string _penaltyFile){
	Reader r(_penaltyFile);
	r.open();
	map<string, double> penalties;
	if(!(r.is_open())){
		cerr << "WARNING: Unable to open" << _penaltyFile << endl;
		exit(0);
	}

	vector<string> lines = r.getAllLines();

	vector<string> AA;
	vector<double> ener;
	int penCounter = 0;
	hasProperty.clear();

	//find a way to skip the first 6 lines
	for (uint i=6; i<lines.size()-1; i++){
		vector<string> tokens = MslTools::tokenize(lines[i], " ");
		if (tokens.size() == 0){
			continue;
		}
		else{
			if (tokens[0] == "PENALTY_DEFINITION"){
				penAAs.push_back(AA);
				penFracs.push_back(ener);
				hasProperty.push_back(false); //defaults to false and will change if property is found
				continue;
			}
			else if (tokens[0] == "END_PENALTY_DEFINITION"){
				for (uint k=0; k<AA.size(); k++){
					pMap[AA[k]] = ener[k];
				}
				AA.clear();
				ener.clear();
				penCounter++;
			}
			else{
				for (uint j=1; j<tokens.size(); j++){
					if (tokens[0] == "TYPE"){
						AA.push_back(MslTools::toUpper(tokens[j]));
						penAAs[penCounter].push_back(MslTools::toUpper(tokens[j]));
					}
					else if (tokens[0] == "PROPERTY" || tokens[0] == "NON-PROPERTY"){
						AA = propertyVecMap[MslTools::toUpper(tokens[j])];
						penAAs[penCounter] = propertyVecMap[MslTools::toUpper(tokens[j])];
						hasProperty[penCounter] = true;
					}
					if (tokens[0] == "PENALTIES"){
						ener.push_back(MslTools::toDouble(tokens[j]));
					}
					if (tokens[0] == "FRACTION"){
						penFracs[penCounter].push_back(MslTools::toDouble(tokens[j]));
					}
					if (tokens[0] == "FUNCTION"){
						function = MslTools::toUpper(tokens[j]);
					}
				}
			}
		}
	}
}

double BaselineAAComposition::calcMultiplicationFactor(double _count, double &_factor, string _function){
	//TODO: add in the default for factor and all other functions
	double f = 0;
	if (_function == ""){
		f = _factor*_count;
		_factor = f;
	}
	//if (_function == "LINEAR"){
	//	f = factor*count;
	//}
	//if (_function == "QUADRATIC"){
	//	f = factor*count;
	//}
	//if (_function == "EXPONENTIAL"){
	//	f = factor*count;
	//}
	return _factor;
}

void BaselineAAComposition::buildInteraction(){
	pSeq.clear();
	energy = 0;
	setCurrentSeq();
	
	for (uint p=0; p<penAAs.size(); p++){
		if (hasProperty.size() != 0){
			if (hasProperty[p] == false){
				for (uint i=0; i<penAAs[p].size(); i++){
					double count = 0;
					double f = factor;
					for (uint j=0; j<pSeq.size(); j++){
						if (penAAs[p][i] == pSeq[j]){
							count++;
							//cout << penAAs[p][i] << ": " << pSeq[j] << endl;
							//cout << count << endl;
							//cout << i << ": " << j << endl;
						}
					}
					if (count > penFracs[p][i]){
						f = calcMultiplicationFactor(count-penFracs[p][i], factor, function);
						energy += pMap.at(penAAs[p][i])*f;
						//cout << f << "*" << energy << endl;
					}
				}
			}
			else{
				double count = 0;
				double f = factor;
				for (uint i=0; i<penAAs[p].size(); i++){
					double f = factor;
					for (uint j=0; j<pSeq.size(); j++){
						if (penAAs[p][i] == pSeq[j]){
							count++;
							//cout << penAAs[p][i] << ": " << pSeq[j] << endl;
							//cout << count << endl;
							//cout << i << ": " << j << endl;
						}
					}
				}
				if (count > penFracs[p][0]){
					f = calcMultiplicationFactor(count-penFracs[p][0], f, function);
					energy += pMap.at(penAAs[p][0])*f;
				}
			}
			
		}
	}
	if (params.size() == 0){
		params.push_back(energy);
	}
	else{
		params[0] = energy;
	}
}

//The below function is what is called by the EnergySet to get the energies for each interaction; I'm just going to get one energy from the whole term when it calls it, based around sequence
double BaselineAAComposition::getEnergy() {
	//sys->getEnergySet()->eraseTerm("BASELINE_COMPOSITION");
	//sys->getEnergySet()->addInteraction(new BaselineAAComposition(energy));
	buildInteraction();
	return(params[0]);
}//TODO: this function will not work for what I need with composition: I need to know how many AAs there are and calculate accordingly
// I think NOW that I know how both MSL and Rosetta work, I actually need to rethink this and implement the code here. If I implement it in the design code, I can't access things that need to be changed here. This getEnergy() should basically just calculate the total energy of the interaction rather than add it up, using functions from there. I think that's the best way to do it
// So instead I should take in the system and AAs and penalties in the design code to save them as private variables here, then call them with this getEnergy() into another function that actually calculates the energy by checking the sequence, adding them up, and return that sum instead...question is, if I'm referencing the system, is the system ACTUALLY undergoing changes within the SPM, or something else? This could be something I run into, but I think it's likely just referencing the system


