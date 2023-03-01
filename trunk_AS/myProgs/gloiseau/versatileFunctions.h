#ifndef VERSATILEFUNC_H
#define VERSATILEFUNC_H

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

// Set of versatile functions that can be used in many of my programs
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
void getCurrentMoveSizes(double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 double _deltaXLimit, double _deltaCrossLimit, double _deltaAxLimit, double _deltaZLimit, bool &_decreaseMoveSize, int _moveToPerform);
double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);
double getStandardNormal(RandomNumberGenerator& RNG);
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);

/***********************************
* general functions
 ***********************************/
// string functions
string convertToPolymerSequence(string _seq, int _startResNum);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum);
string convertVectorUintToString(vector<uint> _inputVector);
string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) ;
string generateString(string _backbone, int _length) ;
string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string convertPolymerSeqToOneLetterSeq(Chain &_chain);

// vector functions
vector<uint> convertStringToVectorUint(string _rotamerLevels);

// void functions
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
  vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1);

#endif