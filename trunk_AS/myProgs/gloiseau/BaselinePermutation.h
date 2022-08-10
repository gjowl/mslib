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


#ifndef BASELINEPERMUTATION_H
#define BASELINEPERMUTATION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "Atom.h"
#include "System.h"


namespace MSL { 
	class BaselinePermutation: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from OneBodyPermutation (a prototype object
		 *   for the interaction of one atom)
		 *******************************************************/

		public:
			BaselinePermutation();
			BaselinePermutation(double _energy);
			//BaselinePermutation(Atom & _d1, double _energy);
			BaselinePermutation(System *sys, string type="", double factor=1);

			// should implement an operator= as well 
			BaselinePermutation(BaselinePermutation & _interaction);
			~BaselinePermutation();

			double getEnergy(double _param, std::vector<double> *_paramDerivatives=NULL);
			double getEnergy(std::vector<double> *_paramDerivatives);
			std::vector<double> getEnergyGrad();

			void setCurrentSeq();
			double calcMultiplicationFactor(double _count, double _factor, std::string _type);
			void buildInteraction();

			double getEnergy();
			string getType();
			System* getSysPointers();
			double getFactor();
			void printCurrentSeq();
			void setFactor(int _factor);
			void setMap(std::map<std::string,double>);
			std::map<std::string,int> getAACountMap();
			std::vector<double> getPenaltyFractions();
			std::map<std::string,double> getMap();
			std::string toString();
			void printParameters();

			//unsigned int getType() const;
			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, BaselinePermutation & _term) {_os << _term.toString(); return _os;};
			std::pair<double,std::vector<double> > partialDerivative();
					
		private:
			void setup(System *sys, string type, double factor);
			void copy(BaselinePermutation & _interaction);
			static const std::string typeName;
			std::map<std::string,double> pMap;
			std::string type;
			std::vector<double> penFracs;
			System *sys;
			double energy;
			double factor;
			std::vector<std::string> pSeq;


	};

	inline std::string BaselinePermutation::toString() { 
		if(pAtoms.size() && pAtoms[0]) {
			char c [1000]; 
			sprintf(c, "%s %s %s %9.4f", typeName.c_str(), pAtoms[0]->toString().c_str(),pAtoms[0]->getResidueName().c_str(),params[0]); 
			return (std::string)c; 
		}
		return "";
	};
	inline std::string BaselinePermutation::getName() const {return typeName;}
	inline std::pair<double,std::vector<double> > BaselinePermutation::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = 0.0;
		partials.second.resize(pAtoms.size() * 3, 0.0);
		return partials;
	}

	//inline double BaselinePermutation::getEnergy(){ return energy;}
	inline double BaselinePermutation::getFactor(){ return factor;}
	inline void BaselinePermutation::setFactor(int _factor){ factor = _factor;}
	inline std::vector<double> BaselinePermutation::getPenaltyFractions(){ return penFracs;}
	inline std::map<std::string,double> BaselinePermutation::getMap(){ return pMap;}
	inline std::string BaselinePermutation::getType(){ return type;}
	inline System* BaselinePermutation::getSysPointers(){ return sys;}
	inline void BaselinePermutation::setMap(std::map<std::string,double> _pMap){ pMap = _pMap;}

}

#endif

