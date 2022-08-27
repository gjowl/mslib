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
 *output file functions
 ***********************************/
// function to setup the output directory
void setupOutputDirectory(BBOptions &_opt);
void setupOutputDirectoryChtc(BBOptions &_opt);

//
void deleteTerminalInteractions(System &_sys, BBOptions &_opt, int _firstResiNum, int _lastResiNum);

/***********************************
 *functions from geomRepack
 ***********************************/
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, BBOptions &_opt);
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, BBOptions &_opt, vector<int> &_rotamerSampling);

/***********************************
 *repack functions
 ***********************************/
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);

/***********************************
 *monomer functions
 ***********************************/
double computeMonomerEnergy(System & _sys, BBOptions &_opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout);

/***********************************
 *option parser
 ***********************************/
// parse config file for given options
BBOptions BBParseOptions(int _argc, char * _argv[], BBOptions defaults);

# endif