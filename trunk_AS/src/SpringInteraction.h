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

#ifndef SPRINGINTERACTION_H
#define SPRINGINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "CharmmBondInteraction.h"

/***********************************************************************
 *  This is essentially a CharmmBondInteraction (inherited) but it adds
 *  the option to only return an energy only if the spring is compressed
 *  or extended.
 *  Both compression and extension are on by default. Extension and
 *  compression can be turned off using setExtension(false) and 
 * setCompression(false).
 ***********************************************************************/

namespace MSL { 

class SpringInteraction: public CharmmBondInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		SpringInteraction();
		SpringInteraction(Atom & _a1, Atom & _a2, double _Kb, double _b0, bool _extensionOn=true, bool _compressionOn=true);

		// should implement an operator= as well 
		SpringInteraction(const SpringInteraction & _interaction);
		~SpringInteraction();

		// the energy for compression and extension can be turned of
		void setCompression(bool _compressionOn); // turn on or off compression
		void setExtension(bool _extensionOn); // turn on or off extension
/*
		/ * setting and getting the parameters * /
		void setParams(std::vector<double> _params);
		void setParams(double _Kb, double _b0);
		double getMinD() const;
		double getConstant() const;
		
		double getEnergy();
		double getEnergy(std::vector<double> *_dd);
*/
		double getEnergy(double _distance, std::vector<double> *_dd=NULL);
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, double Kb, double b0);
		std::string getName() const;
/*
		std::vector<double> getEnergyGrad();

		friend std::ostream & operator<<(std::ostream &_os, SpringInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::pair<double,std::vector<double> > partialDerivative();
*/		
	private:
		void setup(Atom * _a1, Atom * _a2, double _Kb, double _b0, bool _extensionOn=true, bool _compressionOn=true);
		void copy(const SpringInteraction & _interaction);
		bool extension;
		bool compression;

		//static const unsigned int type = 2;
		static const std::string typeName;

};

inline double SpringInteraction::getEnergy(double _distance,std::vector<double> *_dd) { 
	// check if it is compression or extension and check if these were activated before
	// returning the energy
	if ((_distance < params[1] && compression) || (_distance > params[1] && extension)) {
		return CharmmEnergy::instance()->spring(_distance, params[0], params[1],_dd); 
	} else {
		return 0.0;
	}
}
inline std::string SpringInteraction::getName() const {return typeName;}
/*
inline double SpringInteraction::getEnergy() { return getEnergy(pAtoms[0]->distance(*pAtoms[1])); }
inline double SpringInteraction::getEnergy(std::vector<double> *_dd) { 
	if(_dd) {
		
		ouble distance =  CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),_dd);
		return getEnergy(distance,_dd); 
	} 
	return getEnergy();
}
inline void SpringInteraction::setParams(std::vector<double> _params) { if (_params.size() != 2) {std::cerr << "ERROR 49123: invalid number of parameters in inline void SpringInteraction::setParams(std::vector<double> _params)" << std::endl; exit(49123);} params = _params;}
inline void SpringInteraction::setParams(double _Kb, double _b0) {params[0] = _Kb; params[1] = _b0;}
inline double SpringInteraction::getMinD() const {return params[1];};
inline double SpringInteraction::getConstant() const {return params[0];};
inline std::string SpringInteraction::toString() { 
	char c [1000]; 
	sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %20.6f",typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], pAtoms[0]->distance(*pAtoms[1]), getEnergy()); 
	return (std::string)c; 
}
//inline unsigned int SpringInteraction::getType() const {return type;}
inline std::vector<double> SpringInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],*pAtoms[1],params[0],params[1]);
}
inline std::pair<double,std::vector<double> > SpringInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
	return partials;
}
*/

}
#endif

