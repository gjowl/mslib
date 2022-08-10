#include <iostream>
//#include <sstream>
//#include <iterator>
//#include <unistd.h>
#include <thread>
//#include <chrono>
//#include <functional>

// MSL Functions
//#include "System.h"
//#include "CharmmSystemBuilder.h"
//#include "SystemRotamerLoader.h"
//#include "OptionParser.h"
//#include "SelfPairManager.h"
//#include "MslTools.h"
//#include "DeadEndElimination.h"
//#include "SelfConsistentMeanField.h"
//#include "HydrogenBondBuilder.h"
//#include "Transforms.h"
//#include "RandomNumberGenerator.h"
//#include "MonteCarloManager.h"
//#include "AtomSelection.h"
//#include "AtomContainer.h"
//#include "FormatConverter.h"
//#include "CRDReader.h"
//#include "CRDWriter.h"
//#include "SysEnv.h"
//#include "ResidueSelection.h"
//#include "BaselineEnergyBuilder.h"
//#include "BaselineInteraction.h"
//
//// Design Functions: TODO figure out which of these are necessary; a lot of this code I just added internally I believe
//// *There are apparently many redundant functions in designFunctions.h that may be present in many of the other headers
//#include "BaselineIMM1Interaction.h"
//#include "BaselinePairInteraction.h"
//#include "BaselineOuterPairInteraction.h"
//#include "BaselineAAComposition.h"
//#include "BaselineSequenceEntropy.h"
//#include "BaselineSequenceEntropyNormalized.h"
//#include "BaselinePermutation.h"
//#include "SasaCalculator.h"
//#include "designFunctions.h"
//#include "homodimerFunctions.h"
//#include "designOptions.h"
//#include "functions.h"

//using namespace MSL;
using namespace std;

//static SysEnv SYSENV;
//string programName = "seqDesign";//TODO: better name
//string programDescription = "Designs sequences for backbone geometries extracted from the PDB, optimizing specifically for vdW energies";
//string programAuthor = "Gilbert Loiseau";
//string programVersion = "2";
//string programDate = "11 February 2022";
//string mslVersion = MSLVERSION;
//string mslDate = MSLDATE;

//time_t startTime, endTime, spmStart, spmEnd;
//double diffTime, spmTime;

// Functions
//void buildBaselines(System &_sys, Options &_opt);
//void getBestSequenceAndGeometry(System &_startGeom, Options &_opt, RandomNumberGenerator &_RNG,
// vector<uint> _bestState, vector<uint> _interfacePositions, vector<int> _rotamerSampling,
// map<string,double> _sequenceEntropyMap, map<string,map<string,double>> &_sequenceEnergyMap, 
// double _initialEnergy, ofstream &_out);
//void energyFunction(System &_startGeom, Options &_opt, SelfPairManager &_spm,  PolymerSequence _PS, string _id, string _posIdA, string _posIdB,
// string _prevSeq, vector<uint> _prevVec, vector<int> _rotamerSampling, vector<uint> _interfacialPositions, map<string,map<string,double>> _seqEnergyMap,
// map<string,double> _sequenceEntropyMap);
//void searchForBestSequences(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs, vector<uint> &_bestState,
// map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> _sequenceEntropyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, 
// vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, ofstream &_out, ofstream &_err);
//vector<uint> runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, string _rotamerSamplingString, string _variablePositionString, 
// vector<string> _seqs, map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> _sequenceEntropyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, ofstream &_out);
//void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, string _variablePosString, double _spmTime, ofstream &_out);
//vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
//vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
//std::vector<pair <int, double> > calculateResidueBurial (System &_sys) ;
//std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq) ;
////TODO: change all of the original interfacePositions to variablePositions and the allInterfacePositions to interfacePositions
//void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);
//PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
// string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
// vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
//void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
// AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
// CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
//void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS,
//map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
//vector<unsigned int> &_bestState, vector<string> &_seqs, vector<string> &_allSeqs,
//vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<uint> &_allInterfacialPositionsList,
//vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, ofstream &_out, ofstream &_err);
//string generateMultiIDAtSinglePosPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, int _interfacialPosition);
//
///***********************************
// *help functions
// ***********************************/
//void usage();
//void help(Options defaults);
//void outputErrorMessage(Options &_opt);

void foo(int i){
	cout << "foo " << i << endl;
	this_thread::sleep_for(std::chrono::seconds(1));
}
static bool finished = false;

void DoWork(int x)
{
    while (!finished)
    {
        cout << x << endl;
        this_thread::sleep_for(chrono::milliseconds(1000));
    }
}

/******************************************
 *
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){

	thread t(DoWork, 5);
	if (t.joinable()){
		cout << "thread joinable" << endl;
	} else {
		cout << "thread not joinable" << endl;
	}
	//t = thread(energyFunction, ref(_startGeom), ref(_opt), ref(spm), PS, "ALA", posIdA, posIdB, startSeq, stateVec,
	// _rotamerSampling, _interfacePositions, _sequenceEnergyMap, _sequenceEntropyMap); 
	if (t.joinable()){
		cout << "thread joinable" << endl;
	} else {
		cout << "thread not joinable" << endl;
	}
    cin.get();
	finished = true;
	t.join();
	if (t.joinable()){
		cout << "thread joinable" << endl;
	} else {
		cout << "thread not joinable" << endl;
	}
	cout << "thread finished" << endl;
}
