#ifndef BACKBONEOPTIMIZERFUNCTIONS_H
#define BACKBONEOPTIMIZERFUNCTIONS_H

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

// delete interactions at the termini
void deleteTerminalInteractions(System &_sys, BBOptions &_opt, int _firstResiNum, int _lastResiNum);


/***********************************
 *geometry setup functions
 ***********************************/
// builds the system and sets up the energy terms
void prepareSystem(BBOptions &_opt, System &_sys, System &_startGeom, string &_polySeq);
// set starting geometry
void setGly69ToStartingGeometry(BBOptions &_opt, System &_sys, System &_helicalAxis, Transforms &_trans, double _crossingAngle, double _xShift);
// set starting geometry in threadDockAndRepack
void setGly69ToStartingGeometry(BBOptions &_opt, map<string,double> _geometries, System &_sys, System &_helicalAxis, Transforms &_trans);

/***********************************
 *repack functions
 ***********************************/
// runs a greedy to quickly repack sidechains
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
// thread function that drives the docking and repacking for each thread and geometry combination 
void threadDockAndRepack(BBOptions &_opt, System &_helicalAxis, int _thread, int _repackNumber, double _crossingAngle, double _monomerEnergy, 
 map<string,double> _monomerEnergyByTerm, ofstream &_eout);
// get a random axial rotation and zShift 
void getAxialRotAndZShift(BBOptions &_opt, map<string,double> &_geometries, vector<double> &_densities);
// docking function
void localXShiftDocking(System &_sys, BBOptions &_opt, double &_bestEnergy, double _monomerEnergy, SelfPairManager &_spm, 
 System &_helicalAxis, Transforms &_trans, int _thread, double &_savedXShift, ofstream &_out);
// setup function for monteCarloRepack
void localBackboneRepack(BBOptions &_opt, System &_sys, map<string,double> _geometries, vector<double> _densities, SelfPairManager &_spm,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm,
 double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout);
// backbone repack using monte carlo
void monteCarloRepack(BBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string,double> _geometries, vector<double> _densities,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 map<string,double> _monomerEnergyByTerm, double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout);
// change the move size for the 
void getCurrentMoveSizes(BBOptions &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);
// increase the move size for an accepted move
double increaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);


/***********************************
 *monomer functions
 ***********************************/
double computeMonomerEnergy(System & _sys, BBOptions &_opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout);

/***********************************
 *option parser
 ***********************************/
// parse config file for given options
BBOptions BBParseOptions(int _argc, char * _argv[]);

# endif