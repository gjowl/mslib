#ifndef CATMFUNCTIONS_H
#define CATMFUNCTIONS_H

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
#include "catmOptions.h"
#include "functions.h"

//In no particular order at this point
void printOptions(catmOptions & _op, ofstream & _fout);
unsigned int CAOfClosestApproach(Chain & _chainA, Chain & _chainB);
map<string, unsigned int> interfaceResidueCheck(AtomPointerVector & _chainA, AtomPointerVector & _chainB);
string convertToPolymerSequence(string _seq, int _startResNum);
void reThreadResidues(vector<Position*> & positions, int offset);
bool hydrogenBondCheck(System & _sys, catmOptions &_opt, vector<string> _parsedGeoInformation, double & _xShiftStart);
void renumberResidues(System& _sys, int _startResNum);
bool rulesCheck(System & _sys, string _geoIndex, map<int, string> _rulesMap);
void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable);
void addStructureToTopMCRepackList(vector<vector<double> > & _topMCRepackVector, vector<double> & _currentStuctureInfo, ofstream & _fout, unsigned int _numberOfStructuresToMCRepack);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
void deleteTerminalBondInteractions(System &_sys, catmOptions & _opt);
map<string,double> getEnergyByTerm(EnergySet* _eSet);
map<string,double> getEnergyByTermDoubled(EnergySet* _eSet);
void readRulesFile(string _fileName, map<int, string> & _rulesFileMap);
vector<string> getInterHelicalHbonds(EnergySet* & _ESet);
void readGeometryFile(string _filename, vector<string>& _fileVec);
void clusterSolutions(vector<HelixDimerCluster*>& clusters,vector<HelixDimer*>& _structures, double _rmsdCutoff, string _origSeq, string _builtSeq, catmOptions & _opt);
vector<vector<vector<vector<bool> > > > createSavedEnergyFlagTable(System & _sys);
double computeMonomerEnergy(System & _sys, Transforms & _trans, catmOptions & _opt, System & _helicalAxis, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects);
void setupOutputDirectory(catmOptions &_opt, string &_logFile);

# endif