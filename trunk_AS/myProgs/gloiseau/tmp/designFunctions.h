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

// Design Functions: TODO figure out which of these are necessary; a lot of this code I just added internally I believe
#include "BaselinePairInteraction.h"
#include "SasaCalculator.h"

using namespace std;
using namespace MSL;

/***********************************
 *load rotamer functions
 ***********************************/
//Uses rotamer sampling defined by given interface values to load rotamers by position
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, bool _useInterfaceBurial, vector<string> _repackLevels, string _SL,
 vector<uint> _rotamerSampling);
//load rotamers for dimer with different rotamer levels at each position
void loadRotamersByBurial(System &_sys, SystemRotamerLoader &_sysRot, vector<string> _repackLevels, vector<uint> _rotamerSampling);
//load rotamers for interfacial positions
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<uint> _interface);

/***********************************
 *define interface and rotamer sampling
 ***********************************/
//
vector<uint> getLinkedPositions(vector<uint> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel);
// Convert positions to string for setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions) which uses "A,19" "B,19" format!
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<uint> &_interfacialPositions, int _backboneLength);
//unlinks a state vector, replicating the state vector of the linked structure (i.e. for a sequence with 5 AAs and 5 rotamers each: linked: 0,1,2,3,4; unlinked: 0,1,2,3,4,0,1,2,3,4)
vector<uint> unlinkBestState(vector<uint> _bestState, vector<uint> _interfacePositions, int _backboneLength);

/***********************************
 *geometry
 ***********************************/
// moves the helices to the origin in 3D space (x=0, y=0, z=0)
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans);

/***********************************
 *string output functions
 ***********************************/
// TODO: if possible, make some of these more multipurpose
string generateBackboneSequence(string _backbone, int _length);
string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
string getInterfaceString(vector<int> _interface, int _seqLength);
string getAlternateIdString(vector<string> _alternateIds);
string getInterfaceSequence(int _interfaceLevelLimit, string _interface, string _sequence);

/***********************************
 *define interface and rotamer sampling
 ***********************************/
//
vector<uint> getVariablePositions(vector<uint> &_interfacialPositions);
//
std::vector<pair <int, double> > calculateResidueBurial (System &_sys);
vector<uint> getAllInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength);
vector<uint> getInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength);

/***********************************
 *baseline energy helper functions
 ***********************************/
// calculates the self energies of a structure
vector<double> calcBaselineEnergies(System &_sys, int _thread, int _backboneLength);
// calculates the pair energies of a structure
vector<double> calcPairBaselineEnergies(System &_sys, int _thread, int _backboneLength);
// adds up all of the values in a vector<double> to get the total energy
double sumEnergyVector(vector<double> _energies);
// baseline builder function; must have a selfEnergyFile and a pairEnergyFile in the options
void buildBaselines(System &_sys, string _selfEnergyFile, string _pairEnergyFile);
map<string, double> readSingleParameters(string _baselineFile);
map<string,map<string,map<uint, double>>> readPairParameters(string _baselineFile);
void buildSelfInteractions(System &_sys, map<string, double> &_selfMap);
void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap);

/***********************************
 *sequence entropy functions
 ***********************************/
// counts the number of each AA in a sequence and outputs a map<string(AA),int(numAA)> ...
map<string,int> getAACountMap(vector<string> _seq);
// calculates the number of permutations possible for a sequence
double calcNumberOfPermutations(map<string,int> _seqAACounts, int _seqLength);
// sets up calculating sequence entropy for an interface (gets the AA counts for the interface
// and calculates the number of permutations for those AAs and the number of interfacial positions)
double calculateSequenceProbability(map<string,int> &_seqCountMap, map<string,double> &_entropyMap, double _numberOfPermutations);

//
void getSasaForStartingSequence(System &_sys, string _sequence, vector<uint> _state, map<string, map<string,double>> &_sequenceEnergyMap);

// gets the interfacial positions from a vector
vector<uint> getAllInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(int _interfaceLevel, vector<uint> &_rotamerSamplingPerPosition);

/***********************************
* sequence search functions
 ***********************************/
// outputs for runSCMFToGetStartingSequence
void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, ofstream &_out);
// gets the sasa score for sequences found during monte carlo search
void getDimerSasa(System &_sys, map<string, vector<uint>> &_sequenceVectorMap, map<string, map<string,double>> &_sequenceEnergyMap);
void getEnergiesForStartingSequence(SelfPairManager &_spm, string _startSequence, vector<uint> &_stateVector,
 vector<uint> _interfacialPositions, map<string, map<string, double>> &_sequenceEnergyMap, map<string, double> &_entropyMap,
 bool _useBaseline);

// for uncommenting and adding to the cpp file from functions.cpp
//string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum);
//string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum);
//string getAlternateIdString(vector<string> _alternateIds);
//string convertVectorUintToString(vector<uint> _inputVector);
//vector<uint> convertStringToVectorUint(string _rotamerLevels);
//string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) ;
//string generateString(string _backbone, int _length) ;
//string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions);
//std::vector < std::vector < bool > > getActiveMask (System &_sys);
//string convertToPolymerSequence(string _seq, int _startResNum);
//string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
//string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum);
//string convertPolymerSeqToOneLetterSeq(Chain &_chain);
//void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
//  vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1);
//double sumEnergyVector(vector<double> _energies);
//void resetEnergySet(System &_sys, vector<string> _energyTermList);
//void writePdb(System &_sys, string _outputDir, string _pdbName);

#endif
