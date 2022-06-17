/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/12
 * @Email:  gjowl04@gmail.com
 * @Filename: design.h
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/15
 */

#ifndef DESIGN_H
#define DESIGN_H

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
#include "BaselineEnergyBuilder.h"
#include "BaselineInteraction.h"
#include "designOptions.h"

// Design Functions: TODO figure out which of these are necessary; a lot of this code I just added internally I believe
#include "BaselineIMM1Interaction.h"
#include "BaselinePairInteraction.h"
#include "BaselineOuterPairInteraction.h"
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"
#include "SasaCalculator.h"

using namespace std;
using namespace MSL;

/***********************************
 *geometry
 ***********************************/
// moves the helices to the origin in 3D space (x=0, y=0, z=0)
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);

/***********************************
 *string output
 ***********************************/
// TODO: if possible, make some of these more multipurpose
string generateBackboneSequence(string _backbone, int _length);
string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum);
string getInterfaceString(vector<int> _interface, int _seqLength);
string getAlternateIdString(vector<string> _alternateIds);

/***********************************
 *repack functions
 ***********************************/
// runs a greedy to quickly repack sidechains
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
// Code Samson made a while back that should get each active ID and set a mask for anything that isn't active
// allows for repacking side chains of only active ID...


/***********************************
 *define interface and rotamer sampling
 ***********************************/
// gets a vector of the rotamer level for each position on the protein
vector<int> getRotamerSampling(string _rotamerLevels);
//
vector<uint> getVariablePositions(vector<int> &_interfacialPositions);
//
std::vector<pair <int, double> > calculateResidueBurial (System &_sys);
// Calculate Residue Burial for use in identifying the interfacial positions and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq);
//TODO: change all of the original interfacePositions to variablePositions and the allInterfacePositions to interfacePositions
void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_interfacePositions, vector<uint> &_variablePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis);

/***********************************
 *output file functions
 ***********************************/
// function to setup the output directory
void setupDesignDirectory(Options &_opt, string _date);
// function for outputting and writing and energy file
void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns);
// function that writes a config file for local geometric repacks by geomRepack.cpp
void makeRepackConfig(Options &_opt, string _sequence, string _designDir, string _designNumber, string _pdbPath, string _crdPath, map<string,double> _energyMap);
// function that writes a config file for running docking on our
void makeDockingConfig(Options &_opt, string _sequence, vector<uint> _state, string _pdbPath, map<string,double> _energyMap, vector<int> _rotamerSampling);
// function that writes and outputs the summarr, configuration, and energy files from this design run
void outputDesignFiles(Options &_opt, string _interface, vector<int> _rotamerSampling, vector<pair<string,vector<uint>>> _sequenceStatePair, map<string,map<string,double>> _sequenceEnergyMap, vector<double> _densities);





/***********************************
 *baseline energy helper functions
 ***********************************/
// calculates the self energies of a structure
vector<double> calcBaselineEnergies(System &_sys, int _seqLength);
// calculates the pair energies of a structure
vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength);
// adds up all of the values in a vector<double> to get the total energy
double sumEnergyVector(vector<double> _energies);

/***********************************
 *calculate energies
 ***********************************/
// function to setup calculating the dimer energies of a vector<pair<string(sequence),vector<uint>(rotamer state)>>
void computeDimerEnergies(System &_sys, Options &_opt, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> _rotamerSamplingPerPosition, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// helper function for computeDimerEnergies: calculates energy for unlinked positions (different rotamers allowed at each position on opposite helices)
void computeDimerEnergy(System &_sys, Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_sequence, vector<uint> &_stateVec, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, int _seqNumber, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err);
// helper function for computeDimerEnergies: calculates energy for linked positions (same rotamer at each position on opposite helices)
void computeDimerEnergiesLinked(System &_sys, Options &_opt, map<string,map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err);

// function to setup calculating the monomer energies of a vector<pair<string(sequence),vector<uint>(rotamer state)>>
void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// helper function for computeMonomerEnergies: calculates energy for monomer without solvation energy
void computeMonomerEnergyNoIMM1(Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// helper function for computeMonomerEnergies: calculates energy for monomer with solvation energy
void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);


/***********************************
 *other helper functions
 ***********************************/


// other functions

// helper function to ensure that all atoms are assigned coordinates (was having problems with this in the past and so I kept it in for troubleshooting)
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);


vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);

/***********************************
 *help functions
 ***********************************/
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

// parse config file for given options
Options parseOptions(int _argc, char * _argv[], Options defaults);

#endif
