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