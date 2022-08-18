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
#include <thread>

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
string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string getInterfaceString(vector<int> _interface, int _seqLength);
string getAlternateIdString(vector<string> _alternateIds);

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

vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
/***********************************
 *output file functions
 ***********************************/
// function to setup the output directory
void setupDesignDirectory(Options &_opt, string _date);
// function for outputting and writing and energy file
void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns);
// function that writes a config file for local geometric repacks by geomRepack.cpp
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
// baseline builder function; must have a selfEnergyFile and a pairEnergyFile in the options
void buildBaselines(System &_sys, Options &_opt);

/***********************************
 *sequence entropy functions
 ***********************************/
void calculateInterfaceSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq,
map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_prevEntropy,
double &_currEntropy, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal, vector<uint> _interfacePositionsList);
// counts the number of each AA in a sequence and outputs a map<string(AA),int(numAA)> ...
map<string,int> getAACountMap(vector<string> _seq);
// calculates the number of permutations possible for a sequence
double calcNumberOfPermutations(map<string,int> _seqAACounts, int _seqLength);
// sets up calculating sequence entropy for an interface (gets the AA counts for the interface
// and calculates the number of permutations for those AAs and the number of interfacial positions)
void interfaceAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<uint> _interfacialPositionsList);
double getInterfaceSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap, vector<uint> _interfacialPositionsList);

/***********************************
 *calculate energies
 ***********************************/
// function to setup calculating the monomer energies of a vector<pair<string(sequence),vector<uint>(rotamer state)>>
void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// helper function for computeMonomerEnergies: calculates energy for monomer without solvation energy
void computeMonomerEnergyNoIMM1(Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// helper function for computeMonomerEnergies: calculates energy for monomer with solvation energy
void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err);
// delete terminal hydrogen bonds
void deleteTerminalHydrogenBondInteractions(System &_sys, int _firstResiNum, int _lastResiNum);
// makes sure that atoms of the system are built; errors out if not
void checkIfAtomsAreBuilt(System &_sys, ofstream &_err);

// gets the interfacial positions from a vector
vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);

/***********************************
* sequence search functions
 ***********************************/
// function that runs a SelfConsistentMeanField algorithm to find the best starting sequence, only taking into account monomer energies
vector<uint> runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, string _rotamerSamplingString,
 string _variablePositionString, vector<string> _seqs, vector<uint> _interfacialPositions, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string, vector<uint>> &_sequenceVectorMap, map<string, double> _sequenceEntropyMap, ofstream &_out);
// outputs for runSCMFToGetStartingSequence
void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, string _variablePosString, double _spmTime, ofstream &_out);
// redacted old version without multithreading
void searchForBestSequences(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs, vector<uint> &_bestState,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> _sequenceEntropyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, 
 vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, ofstream &_out, ofstream &_err);
// gets the sasa score for sequences found during monte carlo search
void getDimerSasa(System &_sys, map<string, vector<uint>> &_sequenceVectorMap, map<string, map<string,double>> &_sequenceEnergyMap);
void getEnergiesForStartingSequence(Options &_opt, SelfPairManager &_spm, string _startSequence,
vector<uint> &_stateVector, vector<uint> _interfacialPositions, map<string, map<string, double>> &_sequenceEnergyMap, map<string, double> &_entropyMap);

// parse config file for given options
Options parseOptions(int _argc, char * _argv[], Options defaults);

#endif
