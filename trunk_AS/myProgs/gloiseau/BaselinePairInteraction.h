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


#ifndef BASELINEPAIRINTERACTION_H
#define BASELINEPAIRINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "Atom.h"


namespace MSL { 
	class BaselinePairInteraction: public TwoBodyInteraction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			BaselinePairInteraction();
			BaselinePairInteraction(Atom & _d1, Atom & _d2, double _energy);

			// should implement an operator= as well 
			BaselinePairInteraction(const BaselinePairInteraction & _interaction);
			~BaselinePairInteraction();

			double getEnergy(double _param, std::vector<double> *_paramDerivatives=NULL);
			double getEnergy(std::vector<double> *_paramDerivatives);
			std::vector<double> getEnergyGrad();//THESE ARE VIRTUAL FUNCTIONS IN TWO BODY INTERACTION AND MUST BE HERE FOR IT TO INITIALIZE

			double getEnergy();
			std::string toString() ;
			void printParameters();

			//unsigned int getType() const;
			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, BaselinePairInteraction & _term) {_os << _term.toString(); return _os;};
			std::pair<double,std::vector<double> > partialDerivative();

					
		private:
			void setup(Atom * _d1, Atom * _d2, double _energy);
			void copy(const BaselinePairInteraction & _interaction);
			static const std::string typeName;
	};

	inline std::string BaselinePairInteraction::toString() { 
		if(pAtoms.size() >= 2 && pAtoms[0] != NULL && pAtoms[1] != NULL) {
			char c [1000]; 
			sprintf(c, "%s %s %s %s %s %9.4f", typeName.c_str(), pAtoms[0]->toString().c_str(),pAtoms[1]->toString().c_str(),pAtoms[0]->getResidueName().c_str(),pAtoms[1]->getResidueName().c_str(),params[0]); 
			return (std::string)c;//TODO: check if this printer workds properly 
		}
		return "";
	};
	inline std::string BaselinePairInteraction::getName() const {return typeName;}
	inline std::pair<double,std::vector<double> > BaselinePairInteraction::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = 0.0;
		partials.second.resize(pAtoms.size() * 3, 0.0);
		return partials;
	}
}

#endif

