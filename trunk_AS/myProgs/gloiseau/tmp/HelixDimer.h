#ifndef HELIXDIMER_H
#define HELIXDIMER_H

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"

class HelixDimer : public AtomContainer {
	public:
	HelixDimer(string _name, double _energy, int _thread);

	void print();
	void setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds, int _centroidIDNum);

	void setDeltaEnergyByTerm(map<string,double> _deltaEbyTerm) {deltaEByTerm = _deltaEbyTerm;}

	map<string,double>& getDeltaEnergyByTerm(){return deltaEByTerm;}
	void addAxisAtom(Atom* _atom){axisAtoms.push_back(_atom);}
	AtomPointerVector& getAxes(){return axisAtoms;}

	bool operator< (HelixDimer& _b) const {
		if(energy < _b.energy) {
			return true;
		} else {
			return false;
		}
	}

	void printHelixDimerDetails(ofstream & _fout);
	int getNumHbonds() {return hbonds.size();}
	int getThread() {return thread;}
	int getCentId() {return centroidIDNum;}

	string getName() {return name;}
	double getEnergy() {return energy;}
	double getXShift() {return xshift;}
	double getZShift() {return zshift;}
	double getAxialRotation() {return axialrot;}
	double getRelZShift() {return relZshift;}
	double getRelAxialRotation() {return relAxialrot;}
	double getCrossingAngle() {return crossang;}
	string getInterface() {return interface;}
	string getProlineMask() {return prolineMask;}

	// Format: A,11,CA:B,12,O=2.6;
	vector<string>& getHbonds() {return hbonds;}


	private:

	string name;
	double energy;
	double xshift;
	double zshift;
	double relZshift;
	double axialrot;
	double relAxialrot;
	double crossang;
	string interface;
	int thread;
	string prolineMask;
	int centroidIDNum;
	vector<string> hbonds;
	map<string,double> deltaEByTerm;
	AtomPointerVector axisAtoms;
};

# endif
