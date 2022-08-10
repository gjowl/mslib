/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2015 The MSL Developer Group (see README.TXT)
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

#include "SpringInteraction.h"

using namespace MSL;
using namespace std;


const string SpringInteraction::typeName = "SPRING_INTERACTION";

SpringInteraction::SpringInteraction() {
	setup(NULL, NULL, 0.0, 0.0, true, true);
}

SpringInteraction::SpringInteraction(Atom & _a1, Atom & _a2, double _Kb, double _b0, bool _extensionOn, bool _compressionOn) {
	setup (&_a1, &_a2, _Kb, _b0, _extensionOn, _compressionOn);
}

SpringInteraction::SpringInteraction(const SpringInteraction & _interaction) {
	setup(NULL, NULL, 0.0, 0.0, true, true);
	copy(_interaction);
}

SpringInteraction::~SpringInteraction() {
}




void SpringInteraction::setup(Atom * _pA1, Atom * _pA2, double _Kb, double _b0, bool _extensionOn, bool _compressionOn) {
	pAtoms = vector<Atom*> (2, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2);	
	params = vector<double>(2, 0.0);
	setParams(_Kb, _b0);
	extension = _extensionOn;
	compression = _compressionOn;
}

void SpringInteraction::copy(const SpringInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	extension = _interaction.extension;
	compression = _interaction.compression;
}

std::vector<double> SpringInteraction::getEnergyGrad(Atom& a1, Atom& a2, double Kb, double b0) {
	std::vector<double> dd;
	double distance = CartesianGeometry::distanceDerivative(a1.getCoor(), a2.getCoor(),&dd);
	if ((distance < params[1] && compression) || (distance > params[1] && extension)) {
		CharmmEnergy::instance()->springGrad(dd, distance, Kb, b0);
	}
	return dd;
}
