#ifndef MULTICODE_H
#define MULTICODE_H

#include <sstream>
#include <thread>

#include "System.h"
#include "AtomSelection.h"
#include "SelfPairManager.h"
#include "Transforms.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"

using namespace std;
using namespace MSL;

/***********************************
 * System Functions
 ***********************************/
//load rotamers for non-interfacial positions
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);
// delete terminal hydrogen bonds
void deleteTerminalBondInteractions(System &_sys, vector<string> &_deleteTerminalInteractionList);
// sequence search functions
void setActiveSequence(System &_sys, string _sequence);
// makes sure that atoms of the system are built; errors out if not
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);
// checks to ensure that the atoms are built
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);
// gets a mask for use in calculating the energy of only the active sequence
std::vector < std::vector < bool > > getActiveMask (System &_sys);
// writes a pdb file
void writePdb(System &_sys, string _outputDir, string _pdbName);

/***********************************
 * EnergySet Functions
 ***********************************/
// resets the energy set of a system for a list of energy terms
void resetEnergySet(System &_sys, vector<string> _energyTermList);
// gets the energy of a system for a list of energy terms
map<string,double> getEnergyByTerm(EnergySet* _eSet);
// gets the energy of a system for a list of energy terms and doubles them (for monomers)
map<string,double> getEnergyByTermDoubled(EnergySet* _eSet);


// gets the sum of a vector of doubles
double sumDoubleVector(vector<double> _vector);
// runs a greedy to quickly repack sidechains
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);

/***********************************
* geometry move functions
 ***********************************/
double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);
double getStandardNormal(RandomNumberGenerator& RNG);
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);

#endif