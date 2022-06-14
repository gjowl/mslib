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


#ifndef BASELINESEQUENCEENTROPYNORMALIZED_H
#define BASELINESEQUENCEENTROPYNORMALIZED_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "Atom.h"
#include "System.h"


namespace MSL { 
	class BaselineSequenceEntropyNormalized: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from OneBodySequenceEntropyNormalized (a prototype object
		 *   for the interaction of one atom)
		 *******************************************************/

		public:
			BaselineSequenceEntropyNormalized();
			BaselineSequenceEntropyNormalized(double _energy);
			//BaselineSequenceEntropyNormalized(Atom & _d1, double _energy);
			BaselineSequenceEntropyNormalized(System *sys, string type="");

			// should implement an operator= as well 
			BaselineSequenceEntropyNormalized(BaselineSequenceEntropyNormalized & _interaction);
			~BaselineSequenceEntropyNormalized();

			double getEnergy(double _param, std::vector<double> *_paramDerivatives=NULL);
			double getEnergy(std::vector<double> *_paramDerivatives);
			std::vector<double> getEnergyGrad();

			void setCurrentSeq();
			void calcNumberPermutations(double &_numPermutation, vector<int> _counts);
			//void calcNumeratorPermutations(double &_numPermutation, int _uniqueAAs);
			//void calcDenominatorPermutations(double &_denomPermutation, int _AAcount);
			void calcTotalCombinations(int _uniqueAAs);
			void buildInteraction();

			double getEnergy();
			string getType();
			System* getSysPointers();
			double getNumPermutations();
			double getTotCombinations();
			void printCurrentSeq();
			void setFactor(int _factor);
			void setMap(std::map<std::string,double>);
			void setPolyLeuEner(double _energy);
			void setPolyLeuProb(double _prob);
			void setTotalProb(double _prob);
			std::map<std::string,int> getAACountMap();
			std::vector<double> getPenaltyFractions();
			std::map<std::string,double> getMap();
			std::string toString();
			void printParameters();
			void setEnergy(double _energy);

			//unsigned int getType() const;
			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, BaselineSequenceEntropyNormalized & _term) {_os << _term.toString(); return _os;};
			std::pair<double,std::vector<double> > partialDerivative();
					
		private:
			void setup(System *sys, string type);
			void copy(BaselineSequenceEntropyNormalized & _interaction);
			static const std::string typeName;
			std::map<std::string,double> pMap;
			std::string type;
			std::vector<double> penFracs;
			System *sys;
			double energy;
			double seqProb;
			double numPermutations = 1;
			double totCombinations = 1;
			double polyLeuEner = 0;
			double polyLeuProb = 0;
			double totalProb = 0;
			std::vector<std::string> pSeq;


	};

	inline std::string BaselineSequenceEntropyNormalized::toString() { 
		if(pAtoms.size() && pAtoms[0]) {
			char c [1000]; 
			sprintf(c, "%s %s %s %9.4f", typeName.c_str(), pAtoms[0]->toString().c_str(),pAtoms[0]->getResidueName().c_str(),params[0]); 
			return (std::string)c; 
		}
		return "";
	};
	inline std::string BaselineSequenceEntropyNormalized::getName() const {return typeName;}
	inline std::pair<double,std::vector<double> > BaselineSequenceEntropyNormalized::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = 0.0;
		partials.second.resize(pAtoms.size() * 3, 0.0);
		return partials;
	}

	//inline double BaselineSequenceEntropyNormalized::getEnergy(){ return energy;}
	inline double BaselineSequenceEntropyNormalized::getNumPermutations(){ return numPermutations;}
	inline double BaselineSequenceEntropyNormalized::getTotCombinations(){ return totCombinations;}
	inline std::vector<double> BaselineSequenceEntropyNormalized::getPenaltyFractions(){ return penFracs;}
	inline std::map<std::string,double> BaselineSequenceEntropyNormalized::getMap(){ return pMap;}
	inline std::string BaselineSequenceEntropyNormalized::getType(){ return type;}
	inline System* BaselineSequenceEntropyNormalized::getSysPointers(){ return sys;}
	inline void BaselineSequenceEntropyNormalized::setMap(std::map<std::string,double> _pMap){ pMap = _pMap;}
	//inline void setEnergy(double _energy){ energy = _energy};

}

#endif

