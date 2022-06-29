#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <sstream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
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
#include "ResidueSelection.h"
#include "options.h"

/***********************************
 *output file functions
 ***********************************/
// function to setup the output directory
void setupOutputDirectory(Options &_opt);

//
void deleteTerminalHydrogenBondInteractions(System &_sys, int _firstResiNum, int _lastResiNum);
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL);
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);

/***********************************
 *functions from geomRepack
 ***********************************/
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt);
string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum);

/***********************************
 *energy output functions
 ***********************************/
map<string,double> getEnergyByTerm(EnergySet* _eSet);//get the energy of each term and input into a map
map<string,double> getEnergyByTermDoubled(EnergySet* _eSet);//double the energy of each term and input into a map: for monomer

/***********************************
 *repack functions
 ***********************************/
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);

/***********************************
 *monomer functions
 ***********************************/
double computeMonomerEnergy(System & _sys, Options& _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout);
double getStandardNormal(RandomNumberGenerator& RNG);

/***********************************
 *help functions
 ***********************************/
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

// parse config file for given options
Options parseOptions(int _argc, char * _argv[], Options defaults);


# endif