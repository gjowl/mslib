#ifndef FUNCTIONS_H
#define FUNCTIONS_H

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

// Set of general functions that is used in most of my programs
double getStandardNormal(RandomNumberGenerator& RNG);
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
map<string,double> getEnergyByTerm(EnergySet* _eSet);
map<string,double> getEnergyByTermDoubled(EnergySet* _eSet);
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);
// decrease the move size for an accepted move
double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);

// runs a greedy to quickly repack sidechains
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
//load rotamers for a monomer
void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot);
//load rotamers for non-interfacial positions
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);
//below function only loads rotamers onto the interfacial positions by interfacialPositions (01 where 0 = non-interfacial and 1 = interfacial)
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<int> _interface);

//functions from seqDesign
string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum);
string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum);
string getAlternateIdString(vector<string> _alternateIds);
string convertVectorUintToString(vector<uint> _inputVector);
vector<uint> convertStringToVectorUint(string _rotamerLevels);
string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) ;
string generateString(string _backbone, int _length) ;
string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
std::vector < std::vector < bool > > getActiveMask (System &_sys);
string convertToPolymerSequence(string _seq, int _startResNum);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum);
string convertPolymerSeqToOneLetterSeq(Chain &_chain);
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
  vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1);
double sumEnergyVector(vector<double> _energies);
void resetEnergySet(System &_sys, vector<string> _energyTermList);
void writePdb(System &_sys, string _outputDir, string _pdbName);

# endif