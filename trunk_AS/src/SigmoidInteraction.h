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

#ifndef SIGMOIDINTERACTION_H
#define SIGMOIDINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"


namespace MSL {



class SigmoidInteraction: public TwoBodyInteraction {
	/***********************************************************************
	 * This is a sigmoidal energy function which will be used to model
	 * restraints between positions. See Ref from Baker Lab:
	 * 
	 * "Robust and accurate prediction of residue-residue interactions across 
	 *  protein interfaces using evolutionary information"
	 *
	 *  Ovchinnikov, Kamisetty, and Baker 2014. eLife 2014;3:e02030
	 *  DOI: http://dx.doi.org/10.7554/eLife.02030
	 *
	 *
	 ***********************************************************************/
	
	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/
	public:
		SigmoidInteraction();
		SigmoidInteraction(Atom & _a1, Atom & _a2, double _weight, double _slope, double  _cutoff, double _intercept);

		SigmoidInteraction (const SigmoidInteraction &_interaction);
		~SigmoidInteraction();

		//set and get the parameters
		void setParams(std::vector<double> _params);
		void setParams(double _weight, double _slope, double _cutoff, double _intercept);

		double getWeight() const;
		double getSlope() const;
		double getCutoff() const;
		double getIntercept() const;

		double getEnergy();
		double getEnergy(double _distance, std::vector<double> *_dd=NULL);
		double getEnergy(std::vector<double> *_dd);
		double getEnergy( double _dist, double _weight, double _slope, double _cutoff, double _intercept, std::vector<double> *_dd = NULL);

		//Mask the interaction
		bool getMask();
		void setMask(bool _mask);
		
		
		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom &_a1, Atom&_a2, double _weight, double _slope, double _cutoff, double _intercept);

		friend std::ostream & operator << (std::ostream &_os, SigmoidInteraction &_term) {_os << _term.toString(); return _os;};
		std::string toString();

		std::string getName() const;
		std::pair<double, std::vector<double> > partialDerivative();

	private:
		void setup (Atom *_a1, Atom *_a2, double _weight, double _slope, double _cutoff, double _intercept);
		void copy (const SigmoidInteraction &_interaction);
		static const std::string typeName;
		bool mask;
};


//Inline functions
inline void SigmoidInteraction::setParams(std::vector<double> _params) { if (_params.size() != 4) {std::cerr << "ERROR 623846: invalid number of parameters in inline void SigmoidInteraction::setParams(std::vector<double> _params)" << std::endl; exit(623846);} params = _params;}

inline void SigmoidInteraction::setParams(double _weight, double _slope, double _cutoff, double _intercept) {params[0] = _weight; params[1] = _slope; params[2] = _cutoff; params[3] = _intercept;}

inline double SigmoidInteraction::getWeight() const {return params[0];}
inline double SigmoidInteraction::getSlope() const {return params[1];}
inline double SigmoidInteraction::getCutoff() const {return params[2];}
inline double SigmoidInteraction::getIntercept() const {return params[3];}


inline double SigmoidInteraction::getEnergy() {
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]), params[0], params[1], params[2], params[3]); 
}

inline double SigmoidInteraction::getEnergy(std::vector<double> *_dd) {
	if(_dd) {
		double distance = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(), pAtoms[1]->getCoor(), _dd);
		return getEnergy(distance, _dd);
	}
	return getEnergy();

}
inline double SigmoidInteraction::getEnergy(double _dist, std::vector<double> *_dd) {
	return getEnergy(_dist, params[0], params[1], params[2], params[3], _dd);
}

inline double SigmoidInteraction::getEnergy( double _dist, double _weight, double _slope, double _cutoff, double _intercept, std::vector<double> *_dd) {
	
	if (mask == true) {
		return 0;
	}
	double diff = _dist - _cutoff;
	if (_dd != NULL) {
		double e_x = exp(_slope*diff);
		double p = _slope * e_x / pow((1 + e_x),2) + _weight;

		for (int i = 0; i < (*_dd).size(); i++) {
			(*_dd)[i] *= p;
		}
	}

	//Esig = (weight)/(1+exp(-slope(d-cutoff))) + intercept
	return _weight/(1+exp( -_slope*diff )) + _intercept;

}






inline std::string SigmoidInteraction::toString() { 
	char c[1000];
	sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], params[2], params[3], pAtoms[0]->distance(*pAtoms[1]), getEnergy());
	return (std::string)c;
}
inline std::string SigmoidInteraction::getName() const{return typeName;}
inline std::vector<double> SigmoidInteraction::getEnergyGrad() {
	return getEnergyGrad(*pAtoms[0], *pAtoms[1], params[0], params[1], params[2], params[3]);
}

inline std::vector<double> SigmoidInteraction::getEnergyGrad(Atom &_a1, Atom &_a2, double _weight, double _slope, double _cutoff, double _intercept) {
	std::vector<double> dd;
	double distance = CartesianGeometry::distanceDerivative(_a1.getCoor(), _a2.getCoor(), &dd);
	double e_x = exp(_slope*distance);
	double p = _slope * e_x / pow((1 + e_x),2) + _weight;

	for (int i = 0; i < dd.size(); i++) {
		if (mask == true) {
			dd[i] = 0;
		} else {
			dd[i] *= p;
		}
	}

	return dd;


}

inline std::pair<double, std::vector<double> > SigmoidInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
	return partials;
}


//Mask functions
inline bool SigmoidInteraction::getMask() {return mask;}
inline void SigmoidInteraction::setMask(bool _mask) {mask = _mask;}





	
	
}

#endif
