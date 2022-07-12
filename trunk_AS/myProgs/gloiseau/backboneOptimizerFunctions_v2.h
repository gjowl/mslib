#ifndef BACKBONEOPTIMIZERFUNCTIONSV2_H
#define BACKBONEOPTIMIZERFUNCTIONSV2_H

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
#include "backboneOptimizerOptions.h"
#include "SasaCalculator.h"

using namespace std;
using namespace MSL;

/***********************************
 *version 2 functions
 ***********************************/
string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum);

/***********************************
 *output file functions
 ***********************************/
// function to setup the output directory
void setupOutputDirectory(Options &_opt);
void setupOutputDirectoryChtc(Options &_opt);

//
void deleteTerminalHydrogenBondInteractions(System &_sys, int _firstResiNum, int _lastResiNum);

/***********************************
 *functions from geomRepack
 ***********************************/
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt);

/***********************************
 *repack functions
 ***********************************/
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);

/***********************************
 *monomer functions
 ***********************************/
double computeMonomerEnergy(System & _sys, Options& _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout);

/***********************************
 *help functions
 ***********************************/
void backboneOptimizerUsage();
void backboneOptimizerHelp(Options defaults);
void backboneOptimizerOutputErrorMessage(Options &_opt);

// parse config file for given options
Options optimizerParseOptions(int _argc, char * _argv[], Options defaults);

# endif