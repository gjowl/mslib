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


#include "BaselineSequenceEntropy.h"

using namespace MSL;
using namespace std;

const string BaselineSequenceEntropy::typeName = "BASELINE_SEQ_ENTROPY";

BaselineSequenceEntropy::BaselineSequenceEntropy() {
	setup(NULL,0.0);
}

BaselineSequenceEntropy::BaselineSequenceEntropy(Atom & _d1, double _energy) {
	setup (&_d1,_energy);
}

BaselineSequenceEntropy::BaselineSequenceEntropy(const BaselineSequenceEntropy & _interaction) {
	setup(NULL, 0.0);
	copy(_interaction);
}

BaselineSequenceEntropy::~BaselineSequenceEntropy() {
}


void BaselineSequenceEntropy::setup(Atom * _pA1,double _energy) {
	pAtoms.push_back(_pA1);
	params.push_back(_energy);
}

void BaselineSequenceEntropy::copy(const BaselineSequenceEntropy & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;	
}

void BaselineSequenceEntropy::printParameters() {
	if(pAtoms.size() > 0 && pAtoms[0] && params.size() > 0) {
		cout << " ResName " << pAtoms[0]->getResidueName() << endl;
		cout << " atomName " << pAtoms[0]->getName() << endl;
		cout << " energy " << params[0] << endl;
	} else {
		cout << "No atoms defined in interaction" << endl;
	}
}

double BaselineSequenceEntropy::getEnergy(std::vector<double> *_paramDerivatives) {
	return getEnergy(params[0],_paramDerivatives);
}
double BaselineSequenceEntropy::getEnergy(double _dummy, std::vector<double> *_paramDerivatives) {
	// The gradient is zero
	if(_paramDerivatives) {
		_paramDerivatives->resize(pAtoms.size() * 3,0.0);
	}
	return _dummy;

}
std::vector<double> BaselineSequenceEntropy::getEnergyGrad(){
	std::vector<double> foo(pAtoms.size() * 3,0.0);
	return foo;
}

double BaselineSequenceEntropy::getEnergy() {//TODO: likely can add some sort of baseline cutoff here: If it is above the cutoff, add a large amount of energy
	return(params[0]);
}

