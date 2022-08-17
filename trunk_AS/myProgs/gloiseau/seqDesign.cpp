#include <iostream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>
#include <chrono>
#include <functional>

// MSL Functions
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
// *There are apparently many redundant functions in designFunctions.h that may be present in many of the other headers
#include "BaselineIMM1Interaction.h"
#include "BaselinePairInteraction.h"
#include "BaselineOuterPairInteraction.h"
#include "BaselineAAComposition.h"
#include "BaselineSequenceEntropy.h"
#include "BaselineSequenceEntropyNormalized.h"
#include "BaselinePermutation.h"
#include "SasaCalculator.h"
#include "designFunctions.h"
#include "homodimerFunctions.h"
#include "designOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "seqDesign";//TODO: better name
string programDescription = "Designs sequences for backbone geometries extracted from the PDB, optimizing specifically for vdW energies";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "11 February 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime, spmTime;

// Functions
map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, vector<uint> _bestState, double _bestEnergy, double _bestSeqProb, map<string,vector<uint>> &_sequenceStateMap,
 map<string,double> _sequenceEntropyMap, vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList,
 vector<int> _rotamerSampling);
string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap);
void buildBaselines(System &_sys, Options &_opt);
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, vector<uint> _prevVec, double _prevEnergy, double _prevSEProb, 
 string _currSeq, vector<uint> _currVec, vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap);
void searchForBestSequences(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs, vector<uint> &_bestState,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> _sequenceEntropyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, 
 vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, ofstream &_out, ofstream &_err);
void searchForBestSequencesUsingThreads(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, map<string, vector<uint>> &_sequenceVectorMap,
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<int> &_rotamerSampling, int _rep, ofstream &_out, ofstream &_err);
vector<uint> runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, string _rotamerSamplingString,
 string _variablePositionString, vector<string> _seqs, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string, vector<uint>> &_sequenceVectorMap, map<string, double> _sequenceEntropyMap, ofstream &_out);
void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, string _variablePosString, double _spmTime, ofstream &_out);
vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq);
//TODO: change all of the original interfacePositions to variablePositions and the allInterfacePositions to interfacePositions
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
 string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
 CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
 vector<unsigned int> &_bestState, string &_bestSequence, vector<string> &_seqs, vector<string> &_allSeqs,
 map<string,vector<uint>> &_sequenceVectorMap, vector<uint> &_allInterfacialPositionsList,
 vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, int _rep, ofstream &_out, ofstream &_err);
string generateMultiIDAtSinglePosPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, int _interfacialPosition);
void localBackboneRepack(Options &_opt, System &_startGeom, string _sequence, uint _rep, double _savedXShift, System &_helicalAxis, AtomPointerVector &_axisA,
 AtomPointerVector &_axisB, vector<int> _rotamerSampling, Transforms &_trans, RandomNumberGenerator &_RNG, ofstream &_out);
void monteCarloRepack(Options &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 uint _rep, ofstream &_out);
double computeMonomerEnergy(System & _sys, Options & _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout);
void getDimerSasa(System &_sys, map<string, vector<uint>> &_sequenceVectorMap, map<string, map<string,double>> &_sequenceEnergyMap);
void outputFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, map<string,map<string,double>> _sequenceEnergyMap,
 vector<double> _densities);
string getRunParameters(Options &_opt, vector<double> _densities);

/***********************************
 *help functions
 ***********************************/
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

/******************************************
 *
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){
	// initialize time variables
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	// start the timer for the program
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	// setup time and date
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);

	/******************************************************************************
	 *                 === PARSE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
	//Add in some default options that can easily be changed here
	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		outputErrorMessage(opt);
		exit(1);
	} else if (!opt.errorFlag && !opt.warningFlag && opt.errorMessages != ""){
		outputErrorMessage(opt);
		usage();
		exit(0);
	}

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	// summary file output
	ofstream sout;
	// error file output
	ofstream err;
	// rerun config output
	ofstream rerun;

	// makes a directory in the directory that you run from, with design_<runNumber> as the name
	setupDesignDirectory(opt, date);

	// setup output files
	string soutfile = opt.pdbOutputDir + "/summary.out";
	string errfile  = opt.pdbOutputDir + "/errors.out";
	string rerunfile = opt.pdbOutputDir + "/rerun.config";
	
	// open output files
	sout.open(soutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	// setup the rerun config file
	rerun << opt.rerunConf << endl;
	rerun.close();

	// output the date and time to the summary file
	sout << date << endl;

	/******************************************************************************
	 *               === LOAD RANDOM GEOMETRY FROM GEOMETRY FILE ===
	 ******************************************************************************/
	// Initialize RNG with seed (time or given seed number)
	// *Look up easy RNG C++
	RandomNumberGenerator RNG;
	if (opt.useTimeBasedSeed){
		RNG.setTimeBasedSeed();
	} else {
		RNG.setSeed(opt.seed);
	}
	
	// *Change this function for selecting geometries to design
	vector<double> densities;
	if (opt.getGeoFromPDBData){
		getGeometry(opt, RNG, densities, sout);
		cout << "xShift:        " << opt.xShift << "\tDensity: " << densities[0] << endl;
		cout << "crossingAngle: " << opt.crossingAngle << "\tDensity: " << densities[0] << endl;
		cout << "axialRotation: " << opt.axialRotation << "\tDensity: " << densities[1] << endl;
		cout << "zShift:        " << opt.zShift << "\tDensity: " << densities[2] << endl << endl;
	} else {
		cout << "***STARTING GEOMETRY:***" << endl;
		cout << "xShift:        " << opt.xShift << endl;
		cout << "crossingAngle: " << opt.crossingAngle << endl;
		cout << "axialRotation: " << opt.axialRotation << endl;
		cout << "zShift:        " << opt.zShift << endl << endl;
		// TODO: temp fix; instead should search for closest density to the one given
		densities.push_back(0);
		densities.push_back(0);
		densities.push_back(0);
		densities.push_back(0);
	}

	//String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.Ids);
	cout << "Amino acids for design: LEU " << alternateIds << endl;

	/******************************************************************************
	 *                      === GENERATE POLYMER SEQUENCE ===
	 ******************************************************************************/
	// polymer sequences have: chain, starting position of chain residue, three letter AA code
	string polySeq = generatePolymerSequence(opt.backboneAA, opt.backboneLength, opt.thread);
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                         === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,axisA,axisB,ori,xAxis,zAxis,trans);
	
	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// Variables to output from defineInterfaceAndRotamerSampling function
	string rotamerLevels;
	string variablePositionString;
	string rotamerSamplingString;
	// vector of the positions that will be linked
	vector<int> linkedPositions;
	// vector of positions at the interface excluding termini positions
	vector<uint> interfacePositions;
	// vector of positions at the interface including the terminal positions
	vector<uint> allInterfacePositions;
	// vector of rotamer level for each position
	vector<int> rotamerSamplingPerPosition;

	// Defines the interfacial positions and the number of rotamers to give each position
	// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
	// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
	// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
	PolymerSequence interfacePolySeq = getInterfacialPolymerSequence(opt, startGeom, PS, rotamerLevels, variablePositionString, rotamerSamplingString,
	 linkedPositions, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, sout);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	prepareSystem(opt, sys, startGeom, interfacePolySeq);

	Chain chainA = sys.getChain("A");
	Chain chainB = sys.getChain("B");

	AtomPointerVector &apvA = chainA.getAtomPointers();
	AtomPointerVector &apvB = chainB.getAtomPointers();

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// check to verify that all atoms have coordinates
	checkIfAtomsAreBuilt(sys, err);

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	// This should already be ready to go for two independent hetero monomers
	// Would need to rerun this if wanting to expand library of AAs
	if (opt.useBaseline){
		//addBaselineToSelfPairManager();//TODO: was going to make this a function but I feel like it already exists in spm, so going to wait when I have time to look that up
		//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
		buildBaselines(sys, opt);
	}

	// link the interfacial positions (don't do this for hetero)
	if (opt.linkInterfacialPositions){
		vector<vector<string>> linkedPos = convertToLinkedFormat(sys, linkedPositions, opt.backboneLength);
		sys.setLinkedPositions(linkedPos);
	}
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	//decided to try loading low number of rotamers for this instead of high
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	//CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs
	
	/******************************************************************************
	 *           === VARIABLES FOR SAVING ENERGIES AND SEQUENCES ===
	 ******************************************************************************/
	// Initialize vector to hold designed sequences
	vector<string> seqs;
	// Initialize vector to hold all designed sequences from the stateMC
	vector<string> allSeqs;
	//Initialize energyMap to hold all energies for output into a summary file
	map<string, map<string,double>> sequenceEnergyMapBest;
	// Initialize sequence and state pair vector: each sequence will be tied to it's state with the proper rotamers
	map<string, vector<uint>> sequenceVectorMap;
	// Get sequence entropy map
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile);

	/******************************************************************************
	 *                        === SETUP SPM AND RUN SCMF ===
	 ******************************************************************************/
	//TODO: make it so that this part is optional: if I submit a sequence to start, don't even do this. Just find the interface, greedy for best sequence, and continue
	vector<uint> bestState = runSCMFToGetStartingSequence(sys, opt, RNG, rotamerSamplingString, variablePositionString,
	 seqs, sequenceEnergyMapBest, sequenceVectorMap, sequenceEntropyMap, sout);

	// set system to the best sequence or input sequence 
	sys.setActiveRotamers(bestState);
	double bestEnergy = sys.calcEnergy();
	
	/******************************************************************************
	 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
	 ******************************************************************************/
	// Unlink the best state from SCMF if not using linked positions during the state Monte Carlo
	if(opt.linkInterfacialPositions){
		unlinkBestState(opt, bestState, rotamerSamplingPerPosition, opt.backboneLength);
	}
	string bestSequence;
	for (uint i=0; i<3; i++){
		stateMCUnlinked(sys, opt, interfacePolySeq, sequenceEnergyMapBest, sequenceEntropyMap, bestState, bestSequence, seqs, allSeqs,
			sequenceVectorMap, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, i, sout, err);
		localBackboneRepack(opt, sys, bestSequence, i, opt.xShift, helicalAxis, axisA, axisB, rotamerSamplingPerPosition, trans, RNG, sout);
	}
	///******************************************************************************
	// *            === CALCULATE MONOMER ENERGIES OF EACH SEQUENCE ===
	// ******************************************************************************/
	computeMonomerEnergies(opt, trans, sequenceEnergyMapBest, seqs, RNG, sout, err);
	//getSasaDifference(sequenceStatePair, sequenceEnergyMap);

	///******************************************************************************
	// *              === CALCULATE TOTAL ENERGIES AND WRITE PDBS ===
	// ******************************************************************************/
	//// Initialize PDBWriter for designs
	//PDBWriter writer;
	//writer.open(opt.pdbOutputDir + "/allDesigns.pdb");

	//// Output these energy calculations to the summary file
	//sout << "Calculating Final Energies..." << endl;
	//uint i = 0;
	//for (auto &seq: sequenceVectorMap){
	//	string sequence = seq.first;
	//	vector<uint> state = seq.second;
	//	sout << "Sequence " << i+1 << ": " << sequence << endl;
	//	int seqNumber = i;
	//	// calculates the total energy difference between monomer and dimer and outputs individual pdbs for each sequence
	//	getTotalEnergyAndWritePdbs(sys, opt, sequenceEnergyMapBest, sequence, state, rotamerSamplingPerPosition, RNG, seqNumber, writer, sout, err);
	//	i++;
	//}
	//writer.close();

	///******************************************************************************
	// *                   === WRITE OUT ENERGY AND DESIGN FILES ===
	// ******************************************************************************/
	outputFiles(opt, rotamerSamplingString, rotamerSamplingPerPosition, sequenceEnergyMapBest, densities);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

	err.close();
	sout.close();
}

//Functions
// mutate a random position in the sequence to all other given amino acids and save the energetics
string getRunParameters(Options &_opt, vector<double> _densities){
	stringstream ss;
	string t = "\t";
	ss << _opt.xShift << t << _opt.crossingAngle << t << _opt.axialRotation << t << _opt.zShift << t << _densities[0] << t << _densities[1];
	ss << t << _densities[2] << t << _opt.thread << t << _opt.sasaRepackLevel.size() << t << _opt.interfaceLevel << t << _opt.backboneLength;
	string runParameters = ss.str(); 
	return runParameters;
}

void outputFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, map<string,map<string,double>> _sequenceEnergyMap,
 vector<double> _densities){
	// Setup vector to hold energy file lines
	vector<string> energyLines;
	// get the run parameters
	string runParameters = getRunParameters(_opt, _densities);
	string t = "\t";
	stringstream enerTerms;
	// For loop to setup the energy file
	uint i = 0;
	for (auto &seq : _sequenceEnergyMap){
		stringstream seqLine;
		string sequence = seq.first;
		// get the interface sequence
		string interfaceSequence = getInterfaceSequence(_opt,_interface, sequence);
		seqLine << sequence << t << interfaceSequence << t;
		map<string,double> energyMap = _sequenceEnergyMap[sequence];
		// For adding in strings to a line for the energy file
		for (uint j=0; j<_opt.energyTermsToOutput.size(); j++){
			string energyTerm = _opt.energyTermsToOutput[j];
			double energy = energyMap[energyTerm];
			//cout << sequence << ": " << energyTerm << " = " << energy << endl;
			string term = MslTools::doubleToString(energy)+t;
			seqLine << term;
			if (i == 0){
				enerTerms << energyTerm << t;
			}
		}
		//Add in path to design PDB and make repack and docking configuration files
		string seqNumber = MslTools::doubleToString(energyMap.at("SequenceNumber"));
		seqLine << seqNumber << t << runParameters;
		string line = seqLine.str();
		energyLines.push_back(line);
		i++;
	}
	ofstream eout;
	string eoutfile = _opt.pdbOutputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	eout << "Sequence" << t << "InterfaceSequence" << t;
	eout << enerTerms.str();
	eout << "Baseline" << t << "Sequence" << t << "InterfaceSeq" << t << "xShift" << t << "crossingAngle" << t; 
	eout << "axialRotation" << t << "zShift" << t << "angleDistDensity" << t << "axialRotationDensity" << t;
	eout << "zShiftDensity" << t << "repackLevels" << t << "interfaceLevels" << t << "backboneLength" << t;
	cout << "Sequence" << t << "InterfaceSequence" << t;
	cout << enerTerms.str();
	cout << "Baseline" << t << "Sequence" << t << "InterfaceSeq" << t << "xShift" << t << "crossingAngle" << t; 
	cout << "axialRotation" << t << "zShift" << t << "angleDistDensity" << t << "axialRotationDensity" << t;
	cout << "zShiftDensity" << t << "repackLevels" << t << "interfaceLevels" << t << "backboneLength" << t;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
		cout << energyLines[i] << endl;
	}
	eout.close();
}


map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, vector<uint> _bestState, double _bestEnergy, double _bestSeqProb, map<string,vector<uint>> &_sequenceStateMap,
 map<string,double> _sequenceEntropyMap, vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList,
 vector<int> _rotamerSampling){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _interfacialPositionsList.size()-1);
	int interfacePosA = _interfacialPositionsList[rand];
	int interfacePosB = interfacePosA+_opt.backboneLength;

	// Get the random position from the system
	Position &randPosA = _sys.getPosition(interfacePosA);
	Position &randPosB = _sys.getPosition(interfacePosB);
	string posIdA = randPosA.getPositionId();
	string posIdB = randPosB.getPositionId();

	// variable setup for current state
	map<string,map<string,double>> sequenceEnergyMap;
	vector<thread> threads;
	for (uint i=0; i<_opt.Ids.size()-1; i++){
		// pick an identity for each thread 
		int idNum = i;
		// generate polymer sequence for each identity at the corresponding chosen position
		string id = _opt.Ids[idNum];
		// input into the thread function for calculating energies
		string currAA = MslTools::getThreeLetterCode(_bestSeq.substr(interfacePosA, 1));
		if (currAA != id){
			// replace the id at the position in bestSeq with the current id to get current sequence
			string currSeq = _bestSeq;
			string oneLetterId = MslTools::getOneLetterCode(id);
			currSeq.replace(interfacePosA, 1, oneLetterId);
			// switch the position to the given id
			_sys.setActiveIdentity(posIdA, id);
			_sys.setActiveIdentity(posIdB, id);
			// Set a mask and run a greedy to get the best state for the current sequence
			vector<vector<bool>> mask = getActiveMask(_sys);
			_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
			vector<uint> currVec = _spm.getMinStates()[0];
			_sequenceStateMap[currSeq] = currVec;
			// start threading and calculating energies for each identity
			threads.push_back(thread{energyFunction, ref(_opt), ref(_spm), _bestSeq, _bestState, _bestEnergy, _bestSeqProb, currSeq, currVec, ref(_rotamerSampling), ref(_allInterfacialPositionsList), 
			 ref(sequenceEnergyMap), ref(_sequenceEntropyMap)});
		}
	} 
	// join all the threads (wait for them all to finish before continuing)
	for (auto& th : threads){
		th.join();
	}
	return sequenceEnergyMap;
}

string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap){
	string currSeq;
	double currEnergyComparison;
	uint i=0;
	for (auto& seq : _sequenceEnergyMap){
		if (i==0){
			currSeq = seq.first;
			currEnergyComparison = seq.second["energyComparison"];
			i++;
		} else {
			if (seq.second["energyComparison"] > currEnergyComparison){
				string prevSeq = currSeq;
				currSeq = seq.first;
				currEnergyComparison = seq.second["energyComparison"];
			}
		}
	}
	return currSeq;
}

void buildBaselines(System &_sys, Options &_opt){
		map<string, double> selfMap = readSingleParameters(_opt.selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_opt.pairEnergyFile);
		buildSelfInteractions(_sys, selfMap);
		buildPairInteractions(_sys, pairMap);
}

void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, vector<uint> _prevVec, double _prevEnergy, double _prevSEProb, 
 string _currSeq, vector<uint> _currVec, vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap){
	// variable setup
	map<string,double> energyMap;

	// Compute dimer energy
	outputEnergiesByTerm(_spm, _currVec, energyMap, _opt.energyTermList, "Dimer", true);
	double currEnergy = _spm.getStateEnergy(_currVec);

	// Convert the energy term (which actually saves the probability of the sequence in the whole system)
	// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)
	// initialize variables for this thread
	double currEnergyTotal = 0;
	double currSEProb = 0;
	double currEntropy = 0;
	double bestEnergyTotal = 0;
	double prevEntropy = 0;

	//TODO: I just realized that for heterodimers, I may need to completely remake some of these functions as with the below only taking the sequence of one helix; I may make a hetero and homo functions list?
	calculateInterfaceSequenceEntropy(_opt, _prevSeq, _currSeq, _sequenceEntropyMap, _prevSEProb,
	 currSEProb, prevEntropy, currEntropy, _prevEnergy, currEnergy, bestEnergyTotal,
	 currEnergyTotal, _allInterfacePositions);

	// output info
	outputEnergiesByTerm(_spm, _currVec, energyMap, _opt.energyTermList, "Dimer", true);
	double baseline = _spm.getStateEnergy(_currVec, "BASELINE")+_spm.getStateEnergy(_currVec, "BASELINE_PAIR");
	energyMap["EnergyBeforeLocalMC"] = currEnergy-baseline;
	energyMap["Dimer"] = currEnergy-baseline;
	energyMap["Baseline"] = baseline;
	double enerAndSeqEntropy = bestEnergyTotal-currEnergyTotal;
	energyMap["energyComparison"] = enerAndSeqEntropy;
	energyMap["SequenceProbability"] = currSEProb;
	energyMap["bestEnergyTotal"] = bestEnergyTotal;
	energyMap["currEnergyTotal"] = currEnergyTotal;
	energyMap["entropyDiff"] = currEntropy-prevEntropy;
	energyMap["currEntropy"] = currEntropy;
	energyMap["prevEntropy"] = prevEntropy;
	_seqEnergyMap[_currSeq] = energyMap;
}

string generateMultiIDAtSinglePosPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, int _interfacialPosition) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	int startPos = _startResNum;
	int endPos = _startResNum+_seq.length();
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		int pos = it-_seq.begin()+_startResNum;
		if (it == _seq.begin() || it == _seq.end()-1){
		//if (it == _seq.begin()){
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (it == _seq.begin()){
				if(resName == "HIS") {
					ps = ps + " HSE-ACE";
				} else {
					ps = ps + " " + resName + "-ACE";
				}
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
			counter++;
		} else if (pos < startPos+3 || pos > endPos-5){
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			// may need to change the below for polySeq at single pos
			if (pos == _interfacialPosition){
				ps = ps + " [";
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
				for (uint i=0; i<_alternateIds.size(); i++){
					if(_alternateIds[i] != resName){
						if(_alternateIds[i] == "HIS") {
							ps = ps + " HSE";
						} else {
							ps = ps + " " + _alternateIds[i];
						}
					}
				}
				ps = ps + "] ";
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + resName;
				}
			}
			counter++;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

vector<uint> runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, string _rotamerSamplingString,
 string _variablePositionString, vector<string> _seqs, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string, vector<uint>> &_sequenceVectorMap, map<string, double> _sequenceEntropyMap, ofstream &_out){
	// Setup time variables
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);

	// SelfPairManager setup
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&_sys);
	spm.setVerbose(false);
	spm.setRunDEE(_opt.runDEESingles, _opt.runDEEPairs);
	spm.setOnTheFly(false);
	spm.setMCOptions(1000, 0.5, 5000, 3, 10, 1000, 0.01);//changed to sigmoid and added up to 5000
	spm.saveEnergiesByTerm(true); //added back in on 09_21_2021 to get the vdw and hbond energies
	spm.calculateEnergies();

	//Setup running SCMF or UnbiasedMC
	if (_opt.runSCMF == true){
		cout << "Running Self Consistent Mean Field" << endl;
		_out << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	} else {
		cout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		_out << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	// run and find a sequence using the chosen parameters (MCOptions, SCMF, DEE, etc.)
	spm.runOptimizer();
	time(&endTime);
	diffTime = difftime (endTime, startTime);

	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
	_sys.setActiveRotamers(bestState);
	string startSequence = convertPolymerSeqToOneLetterSeq(_sys.getChain("A")); //used for outputting starting sequence
	string interfaceSeq = getInterfaceSequence(_opt, _rotamerSamplingString, startSequence);
	_sequenceVectorMap[startSequence] = bestState;	
	// output spm run optimizer information
	spmRunOptimizerOutput(spm, _sys, interfaceSeq, _variablePositionString, diffTime, _out);
	
	//Add energies for initial sequences into the sequenceEnergyMap
	_seqs.insert(_seqs.begin(), startSequence);
	pair<string,vector<uint>> startSequenceStatePair = make_pair(startSequence, bestState);
	getEnergiesForStartingSequence(_opt, spm, startSequence, bestState, _sequenceEnergyMap, _sequenceEntropyMap);
	getSasaForStartingSequence(_sys, startSequence, bestState, _sequenceEnergyMap);

	return bestState;
}

void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, string _variablePosString, double _spmTime, ofstream &_out){
	// output this information about the SelfPairManager run below
	// vector for the initial SCMF state
	vector<unsigned int> initialState = _spm.getSCMFstate();
	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = _spm.getBestSCMFBiasedMCState();
	
	// checks to see if the sequence changed between the initialState and the bestState
	_sys.setActiveRotamers(initialState);
	string initialSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));
	_sys.setActiveRotamers(bestState);
	string SCMFBestSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));

	// energy terms from the startingState
	double bestEnergy = _spm.getStateEnergy(bestState);
	double hbondEnergy = _spm.getStateEnergy(bestState, "SCWRL4_HBOND");
	double vdwEnergy = _spm.getStateEnergy(bestState, "CHARMM_VDW");

	// outputs	
	_out << "Initial Sequence:   " << initialSeq << endl;
	_out << "Best Sequence:      " << SCMFBestSeq << endl;
	_out << "Interface Sequence: " << _interfaceSeq << endl;
	_out << "Interface:          " << _variablePosString << endl;
	_out << "Total Energy:       " << bestEnergy << endl;
	_out << "VDW:                " << vdwEnergy << endl;
	_out << "HBOND:              " << hbondEnergy << endl;
	_out << endl << "End SelfPairManager Optimization: " << _spmTime << "s" << endl;
}

vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=0; k<_opt.backboneLength; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++)void defineInterfaceAndRotamerSampling(Options &_opt, PolymerSequence _PS, string &_rotamerLevels, string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions, vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out, string _axis){{
	for (uint k=3; k<_opt.backboneLength-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

std::vector<pair <int, double> > calculateResidueBurial (System &_sys) {
	/*
	  SASA reference:
	  Protein Engineering vol.15 no.8 pp.659–667, 2002
	  Quantifying the accessible surface area of protein residues in their local environment
	  Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
	*/
	map<string,double> refSasa;
	refSasa["G"] = 83.91;
	refSasa["A"] = 116.40;
	refSasa["S"] = 125.68;
	refSasa["C"] = 141.48;
	refSasa["P"] = 144.80;
	refSasa["T"] = 148.06;
	refSasa["D"] = 155.37;
	refSasa["V"] = 162.24;
	refSasa["N"] = 168.87;
	refSasa["E"] = 187.16;
	refSasa["Q"] = 189.17;
	refSasa["I"] = 189.95;
	refSasa["L"] = 197.99;
	refSasa["H"] = 198.51;
	refSasa["K"] = 207.49;
	refSasa["M"] = 210.55;
	refSasa["F"] = 223.29;
	refSasa["Y"] = 238.30;
	refSasa["R"] = 249.26;
	refSasa["W"] = 265.42;

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator b(_sys.getAtomPointers());
	//b.setProbeRadius(3.0);
	b.calcSasa();
	for (uint i = 0; i < _sys.positionSize()/2; i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdA = _sys.getPosition(i).getPositionId();
		//string posIdB = _sys.getPosition(i+_sys.positionSize()/2).getPositionId();
		string oneLetter = MslTools::getOneLetterCode(_sys.getPosition(i).getResidueName());
		double resiSasa = b.getResidueSasa(posIdA);
		double burial = resiSasa / refSasa[oneLetter];
		residueBurial.push_back(pair<int,double>(i, burial));
		//residueBurial.push_back(pair<string,double>(posIdB, burial));
		//cout << posId << "\t" << resiSasa << "\t" << burial << endl;
	}

	return residueBurial;
}

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq) {
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);
	PolymerSequence PS(polySeq);

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
	}

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(_opt.backboneCrd);
	if(!cRead.read()) {
		cerr << "Unable to read " << _opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator dimerSasa(_sys.getAtomPointers());
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	dimerSasa.calcSasa();
	monoSasa.calcSasa();
	dimerSasa.setTempFactorWithSasa(true);

	for (uint i = 0; i < monoSys.positionSize(); i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdMonomer = monoSys.getPosition(i).getPositionId();
		string posIdDimer = _sys.getPosition(i).getPositionId();
		double resiSasaMonomer = monoSasa.getResidueSasa(posIdMonomer);
		double resiSasaDimer = dimerSasa.getResidueSasa(posIdDimer);
		double burial = resiSasaDimer/resiSasaMonomer;
		residueBurial.push_back(pair<int,double>(i, burial));

		//set sasa for each residue in the b-factor
		AtomSelection selA(_sys.getPosition(i).getAtomPointers());
		AtomSelection selB(_sys.getPosition(i+_opt.backboneLength).getAtomPointers());
		AtomPointerVector atomsA = selA.select("all");
		AtomPointerVector atomsB = selB.select("all");
		for (AtomPointerVector::iterator k=atomsA.begin(); k!=atomsA.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
		for (AtomPointerVector::iterator k=atomsB.begin(); k!=atomsB.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
	}

	PDBWriter writer;
	writer.open(_opt.pdbOutputDir+"/interfaceSASA.pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	return residueBurial;
}

// function to prepare the system for design
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS){	
	// declare system
	CharmmSystemBuilder CSB(_sys,_opt.topFile,_opt.parFile,_opt.solvFile);
	if (_opt.useElec == false){
		CSB.setBuildTerm("CHARMM_ELEC", false);
	} else {
		CSB.setBuildTerm("CHARMM_ELEC", true);
	}
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(_PS)) {
		cerr << "Unable to build system from " << _PS << endl;
		exit(0);
	}

	// get chain A and B from the _system
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// assign the coordinates of our system to the given geometry 
	_sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	_sys.buildAllAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(_sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = _sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	//Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	int firstPos = 0;
    int lastPos = _sys.positionSize();
    deleteTerminalHydrogenBondInteractions(_sys,firstPos,lastPos);

	// Up to here is from readDPBAndCalcEnergy
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	CSB.updateNonBonded(10,12,50);

	// TODO: maybe calculate the energy here, then see if there's clashing. If so, then move helices away until no clashing, keeping the same
	// other coordinates? Just for simplicity for now. And if I want to implement this in design, adding this in will likely be a good idea.	
	//// as of 2022-7-5: not sure if the above works or needs to be reworked
	//PDBWriter writer1;
	//writer1.open(_osys.getAtomPointers(), true, false, true);
	//writer1.close();
}

//Make it so that this will get all the info I need instead of having to run more code later
void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS,
map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
vector<unsigned int> &_bestState, string &_bestSequence, vector<string> &_seqs, vector<string> &_allSeqs,
map<string,vector<uint>> &_sequenceVectorMap, vector<uint> &_allInterfacialPositionsList,
vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, int _rep, ofstream &_out, ofstream &_err){
	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	prepareSystem(_opt, sys, _sys, _PS);

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	if (_opt.useBaseline){
		//addBaselineToSelfPairManager();//TODO: was going to make this a function but I feel like it already exists in spm, so going to wait when I have time to look that up
		//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
		buildBaselines(sys, _opt);
	}

	/******************************************************************************
	 *              === LOAD ROTAMERS AND CHOOSE TO LINK INTERFACE ===
	 ******************************************************************************/
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);

	/******************************************************************************
	 *                        === SETUP SPM AND RUN DEE ===
	 ******************************************************************************/
	sys.buildAllAtoms(); //currently unsure if this is needed

	// Setup time variables
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	//spm.updateWeights();
	spm.setOnTheFly(false);
	spm.saveEnergiesByTerm(true);
	if (_opt.verbose){
		cout << "Calculating self and pair energy terms..." << endl;
	}
	spm.calculateEnergies();
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	_out << "Time to calculate energies: " << diffTime << endl;

	searchForBestSequencesUsingThreads(sys, _opt, spm, _RNG, _allSeqs, _bestState, _bestSequence, _sequenceEnergyMap, _sequenceVectorMap, 
	_sequenceEntropyMap, _allInterfacialPositionsList, _interfacialPositionsList, _rotamerSampling, _rep, _out, _err);
	
	_seqs.push_back(_bestSequence);
}

void searchForBestSequencesUsingThreads(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, map<string, vector<uint>> &_sequenceVectorMap,
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<int> &_rotamerSampling, int _rep, ofstream &_out, ofstream &_err){
	// Setup time variables
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);

	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects, _opt.MCConvergedSteps, 0.01);
	//MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _spm.getStateEnergy(_bestState);

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	MC.setEner(bestEnergy);

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);

	// Sequence Search Energy Landscape file
	ofstream lout;
	string loutfile  = _opt.pdbOutputDir + "/sequenceSearchEnergyLandscape.out";
	lout.open(loutfile.c_str());
	lout << "***STARTING GEOMETRY:***" << endl;
	lout << "xShift: " << _opt.xShift << endl;
	lout << "crossingAngle: " << _opt.crossingAngle << endl;
	lout << "axialRotation: " << _opt.axialRotation << endl;
	lout << "zShift: " << _opt.zShift << endl << endl;
	lout << "Number of MCCycles: " << _opt.MCCycles << endl;
	lout << "PrevSequence\tCurrSequence\tPrevEnergy\tCurrEnergy\tPrevEntropy\tCurrEntropy\tCurrentTemp" << endl;
	if (_opt.verbose){
		cout << "***STARTING GEOMETRY:***" << endl;
		cout << "xShift: " << _opt.xShift << endl;
		cout << "crossingAngle: " << _opt.crossingAngle << endl;
		cout << "axialRotation: " << _opt.axialRotation << endl;
		cout << "zShift: " << _opt.zShift << endl << endl;
		cout << "Number of MCCycles: " << _opt.MCCycles << endl;
	}

	/******************************************************************************
	 *                      === BEGIN STATE MONTE CARLO ===
	 ******************************************************************************/
	//TODO: I think if I ever want to run this on multiple cores, I should set it up to run those here, then just save x number from each core
	// - could I make local changes to geometry starting from here, then get it to run on multiple cores from here?
	// - what would multiple cores do here? Could I locally test a variety of sequences using multiple cores, and save x number per each core?
	//		 Say save 10 and run 10 replicates, that gets me 100 sequences right away?
	// initialize energy variables for the MonteCarlo
	double bestEnergyTotal = 0;
	double currEnergyTotal = 0;
	double currStateSEProb = 0;
	double prevStateSEProb = 0;
	double prevStateEntropy = 0;
	double currStateEntropy = 0;
	double totEnergy = 0;
	string bestSeq = prevStateSeq;
	map<string,map<string,double>> allSequenceEnergyMap;
	cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	_out << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		// get the sequence entropy probability for the current best sequence
		double bestSeqProb = getSequenceEntropyProbability(_opt, bestSeq, _sequenceEntropyMap);
		map<string,vector<uint>> sequenceVectorMap;
		map<string,map<string,double>> sequenceEnergyMap = mutateRandomPosition(_sys, _opt, _spm, _RNG, bestSeq, prevStateVec, bestEnergy, bestSeqProb, 
		 sequenceVectorMap, _sequenceEntropyMap, _allInterfacialPositionsList, _interfacialPositionsList, _rotamerSampling);

		// get the best sequence and energy for the current mutation position
		string currSeq = getBestSequenceInMap(sequenceEnergyMap);
		
		// TODO: add in check step to see if sequence is already found in the best sequence list and skip if so
		double currEnergyTotal = sequenceEnergyMap[currSeq]["currEnergyTotal"];
		double bestEnergyTotal = sequenceEnergyMap[currSeq]["bestEnergyTotal"];
		vector<uint> currStateVec = sequenceVectorMap[currSeq];
		//cout << "Prev Seq: " << bestSeq << ": " << bestEnergyTotal << endl;
		//cout << "Curr Seq: " << currSeq << ": " << currEnergyTotal << endl;
		MC.setEner(bestEnergyTotal);
		// MC accept and reject conditions
		if (!MC.accept(currEnergyTotal)){
			_sys.setActiveRotamers(prevStateVec);

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			_sys.setActiveRotamers(currStateVec);
			prevStateVec = currStateVec;
			bestSeq = currSeq;
			bestEnergy = sequenceEnergyMap[bestSeq]["Dimer"];
			sequenceEnergyMap[bestSeq]["acceptCycleNumber"] = cycleCounter;
			_sequenceVectorMap[bestSeq] = currStateVec;
			cout << "Best sequence: " << bestSeq << endl;
			cout << "Best sequence Info:" << endl;
			cout << "Baseline          " << sequenceEnergyMap[bestSeq]["Baseline"] << endl;
			cout << "Entropy           " << sequenceEnergyMap[bestSeq]["entropyDiff"] << endl;
			cout << "CurrEntropy       " << sequenceEnergyMap[bestSeq]["currEntropy"] << endl;
			cout << "PrevEntropy       " << sequenceEnergyMap[bestSeq]["prevEntropy"] << endl;
			allSequenceEnergyMap[bestSeq] = sequenceEnergyMap[bestSeq];
			double prevEnergy = bestEnergyTotal;

			if (_opt.energyLandscape){
				map<string,double> energyMap = sequenceEnergyMap[currSeq];
				lout << prevStateSeq << "\t" << currSeq << "\t" << bestEnergyTotal << "\t" << currEnergyTotal << "\t";
				lout << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << "; PrevE=  " << prevEnergy << " : CurrE= " << currEnergyTotal;
				cout << "; CurrTemp: " << MC.getCurrentT() << endl;
			}
			cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to
		//if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
		//	MC.reset(3649, 3649, 500, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		//}
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();
	//addSequencesToVector(energyVector, _allSeqs);
	getDimerSasa(_sys, _sequenceVectorMap, _sequenceEnergyMap);
	uint i=0;
	double ener = 0;
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/allDesigns_" + to_string(_rep) + ".pdb");
	for (auto &seq: allSequenceEnergyMap){
		_out << "Best Sequence #" << i << ": " << seq.first << "; Energy: " << seq.second["Dimer"] << endl;
		if (i == 0){
			_sys.setActiveRotamers(_sequenceVectorMap[seq.first]);
			writer.write(_sys.getAtomPointers(), true, false, true);
			_bestSequence = seq.first;
			ener = seq.second["SequenceProbability"];
		} else if (seq.second["SequenceProbability"] > ener){
			_sys.setActiveRotamers(_sequenceVectorMap[seq.first]);
			writer.write(_sys.getAtomPointers(), true, false, true);
			_bestSequence = seq.first;
			ener = seq.second["SequenceProbability"];
		}
		_sequenceEnergyMap[seq.first] = seq.second;
		i++;
	}
	writer.close();
	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_out << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_out << "Optimization end at Temp: " << MC.getCurrentT() << endl;
	//TODO: maybe set up something like the following:
	// - save sequences for x cycles
	// - do a backbone repack for the best sequence
	// - recalculate the energy for all sequences at this new geometry
	// - continue on to the next cycle with the current sequence
	// repeat for x times 10 cycles
}

void getDimerSasa(System &_sys, map<string, vector<uint>> &_sequenceVectorMap, map<string, map<string,double>> &_sequenceEnergyMap){
	for (auto &seq : _sequenceVectorMap){
		string sequence = seq.first;
		vector<uint> state = seq.second;

		_sys.setActiveRotamers(state);

		//Setup SasaCalculator to calculate the monomer SASA
		SasaCalculator sasa(_sys.getAtomPointers());
		sasa.calcSasa();
		double dimerSasa = sasa.getTotalSasa();

		_sequenceEnergyMap[sequence]["DimerSasa"] = dimerSasa;
	}
}

// old version of code without threading; used for my initial runs for CHIP1 in winter 2021
void searchForBestSequences(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs, vector<uint> &_bestState,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> _sequenceEntropyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, 
 vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, ofstream &_out, ofstream &_err){
	// Setup time variables
	time_t startTimeSMC, endTimeSMC;
	double diffTimeSMC;
	time(&startTimeSMC);

	// Setup MonteCarloManager
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, 10);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _spm.getStateEnergy(_bestState);

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	vector<unsigned int> currStateVec = _bestState;
	MC.setEner(bestEnergy);

	// initialize map for accepting energies
	map<string,double> stateMCEnergies;
	map<vector<uint>, map<string,double>> stateEnergyMap;
	
	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;
	ids.push_back("LEU");

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);

	// initialize energy vectors
	vector<pair<double,string>> energyVector;
	vector<pair<double,vector<uint>>> energyStateVec;

	// Sequence Search Energy Landscape file
	ofstream lout;
	string loutfile  = _opt.pdbOutputDir + "/sequenceSearchEnergyLandscape.out";
	lout.open(loutfile.c_str());
	lout << "***STARTING GEOMETRY:***" << endl;
	lout << "xShift: " << _opt.xShift << endl;
	lout << "crossingAngle: " << _opt.crossingAngle << endl;
	lout << "axialRotation: " << _opt.axialRotation << endl;
	lout << "zShift: " << _opt.zShift << endl << endl;
	lout << "Number of MCCycles: " << _opt.MCCycles << endl;
	lout << "PrevSequence\tCurrSequence\tTotal\tDimer\tBaseline\tVDW\tHBOND\tEnergyDifferencew/prevSeq\tPrevEnergy\tCurrEnergy\tPrevEntropy\tCurrEntropy\tCurrentTemp" << endl;
	if (_opt.verbose){
		cout << "***STARTING GEOMETRY:***" << endl;
		cout << "xShift: " << _opt.xShift << endl;
		cout << "crossingAngle: " << _opt.crossingAngle << endl;
		cout << "axialRotation: " << _opt.axialRotation << endl;
		cout << "zShift: " << _opt.zShift << endl << endl;
		cout << "Number of MCCycles: " << _opt.MCCycles << endl;
	}

	/******************************************************************************
	 *                      === BEGIN STATE MONTE CARLO ===
	 ******************************************************************************/
	//TODO: I think if I ever want to run this on multiple cores, I should set it up to run those here, then just save x number from each core
	// - could I make local changes to geometry starting from here, then get it to run on multiple cores from here?
	// - what would multiple cores do here? Could I locally test a variety of sequences using multiple cores, and save x number per each core?
	//		 Say save 10 and run 10 replicates, that gets me 100 sequences right away?
	//TODO: to get the backbone optimization working, I would need to implement it here and only with one set of sidechains, not the whole set
	// -how do I decide when to start that optimization?
	// - can I make it an option?
	// initialize energy variables for the MonteCarlo
	double bestEnergyTotal = 0;
	double currEnergyTotal = 0;
	double currStateSEProb = 0;
	double prevStateSEProb = 0;
	double prevStateEntropy = 0;
	double currStateEntropy = 0;
	double totEnergy = 0;

	cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	_out << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		// Get energy term and probability for the first sequence (these terms can then get replaced by future sequences that are accepted by MC)
		if (cycleCounter == 0){
			outputEnergiesByTerm(_spm, prevStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
			stateMCEnergies["Dimer"] = bestEnergy;
			double prevSeqProb = getSequenceEntropyProbability(_opt, prevStateSeq, _sequenceEntropyMap);
			stateMCEnergies["SequenceProbability"] = prevSeqProb;

			stateEnergyMap[prevStateVec] = stateMCEnergies;
			totEnergy = stateMCEnergies["EnergyBeforeLocalMC"];

			sequences[totEnergy] = prevStateSeq;
			double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
			//saveSequence(_opt, energyVector, energyStateVec, prevStateSeq, prevStateVec, bestEnergy);
		}
		// Reset the energy map to save energies from new state after changing the rotamer
		stateMCEnergies.clear();
		randomPointMutationUnlinked(_sys, _opt, _RNG, _interfacialPositionsList, ids);

		// Set a mask and run a greedy to get the best state for the sequence
		//_sys.setActiveRotamers(currStateVec);
		vector<vector<bool>> mask = getActiveMask(_sys);
		_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
		currStateVec = _spm.getMinStates()[0];

		// Get the sequence for the random state
		_sys.setActiveRotamers(currStateVec);
		string currStateSeq = convertPolymerSeqToOneLetterSeq(chain);

		// Compute dimer energy
		outputEnergiesByTerm(_spm, currStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
		double currStateEnergy = _spm.getStateEnergy(currStateVec);

		// Convert the energy term (which actually saves the probability of the sequence in the whole _system)
		// to the proper comparison of proportion to energy between individual sequences (done outside of the actual energy term)

		//TODO: I just realized that for heterodimers, I may need to completely remake some of these functions as with the below only taking the sequence of one helix; I may make a hetero and homo functions list?
		calculateInterfaceSequenceEntropy(_opt, prevStateSeq, currStateSeq, _sequenceEntropyMap, prevStateSEProb, 
		 currStateSEProb, prevStateEntropy, currStateEntropy, bestEnergy, currStateEnergy, bestEnergyTotal,
		 currEnergyTotal, _allInterfacialPositionsList);
		MC.setEner(bestEnergyTotal);

		// MC accept and reject conditions
		double currVDW = _spm.getStateEnergy(currStateVec, "CHARMM_VDW");
		double prevVDW = _spm.getStateEnergy(prevStateVec, "CHARMM_VDW");
		double currHBOND = _spm.getStateEnergy(currStateVec, "SCWRL4_HBOND");
		double prevHBOND = _spm.getStateEnergy(prevStateVec, "SCWRL4_HBOND");
		
		if (!MC.accept(currEnergyTotal)){
			_sys.setActiveRotamers(prevStateVec);
			currStateVec = prevStateVec;

			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			//TODO: make these a separate function or put in comments  for them
			bestEnergy = currStateEnergy;
			MC.setEner(currEnergyTotal);
			prevStateSEProb = currStateSEProb;
			string prevStateSeq1 = prevStateSeq;
			prevStateSeq = currStateSeq;
			prevStateVec = currStateVec;
			_sys.setActiveRotamers(currStateVec);
	
			outputEnergiesByTerm(_spm, currStateVec, stateMCEnergies, _opt.energyTermList, "Dimer", true);
			double EnergyBeforeLocalMC = currStateEnergy-(_spm.getStateEnergy(currStateVec, "BASELINE")+_spm.getStateEnergy(currStateVec, "BASELINE_PAIR"));
			stateMCEnergies["EnergyBeforeLocalMC"] = EnergyBeforeLocalMC;
			stateMCEnergies["Dimer"] = EnergyBeforeLocalMC;
			stateMCEnergies["Baseline"] = _spm.getStateEnergy(currStateVec, "BASELINE")+_spm.getStateEnergy(currStateVec, "BASELINE_PAIR");
			stateMCEnergies["EnergyBeforeLocalMCw/seqEntropy"] = bestEnergyTotal-currEnergyTotal;
			stateMCEnergies["SequenceProbability"] = currStateSEProb;
			stateEnergyMap[currStateVec] = stateMCEnergies;

			//TODO: change this so I just save energies in the same place and easily can get vdw, hbond, etc. for each saved sequence
			//saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currStateEnergy);
			if (_opt.weight_seqEntropy == 0){
				saveSequence(_opt, energyVector, energyStateVec, currStateSeq, currStateVec, currVDW);
			} else {
				saveSequence(_opt, _RNG, stateEnergyMap, energyVector, energyStateVec, currStateSeq, currStateVec, currStateEnergy, _out);
			}
			double prevEnergy = bestEnergyTotal;

			if (_opt.energyLandscape){
				map<string,double> energyMap = stateEnergyMap.at(currStateVec);
				lout << prevStateSeq1 << "\t" << currStateSeq << "\t";
				for (uint j=0; j<_opt.energyLandscapeTerms.size(); j++){
					lout << energyMap.at(_opt.energyLandscapeTerms[j]) << "\t";
				}
				lout << bestEnergyTotal << "\t" << currEnergyTotal << "\t" << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			}
			if (_opt.verbose){
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currStateSeq << "; PrevE=  " << prevEnergy << " : CurrE= " << currEnergyTotal << "; PrevVDW: " << prevVDW << " : CurrVDW: " << currVDW << "EnergyDifference" << bestEnergyTotal-currEnergyTotal << "; CurrTemp: " << MC.getCurrentT() << endl;
			}
		cycleCounter++;
		}
		//Reset the MC to run 100 more cycles to
		if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
			MC.reset(3649, 3649, 500, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		}
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();
	addSequencesToVector(energyVector, _allSeqs);
	convertStateMapToSequenceMap(_sys, energyStateVec, stateEnergyMap, _sequenceEnergyMap, _sequenceStatePair, _out);
	getDimerSasaScores(_sys, _sequenceStatePair, _sequenceEnergyMap);
	cout << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
	_out << "End sequence optimization by membrane composition: " << diffTimeSMC << "s" << endl << endl;
}

PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
 string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out){
	// Declare system
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setBuildNonBondedInteractions(false);
	//CSB.setBuildNoTerms();

	if(!CSB.buildSystem(_PS)) {
		cout << "Unable to build system from " << _PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);

	CSB.updateNonBonded(10,12,50);

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct
	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	sys.buildAllAtoms();

	//TODO: add in a comparison for the monomer at each position here instead
	string backboneSeq = generateString(_opt.backboneAA, _opt.backboneLength);
	vector<pair <int, double> > resiBurial = calculateResidueBurial(sys, _opt, backboneSeq);
	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});

	//cout << "Determining interfacial residues by residue burial..." << endl;
	int levelCounter = 0;
	vector<int> interfacePositions;

	// Output variable Set up
	backboneSeq = generateBackboneSequence("L", _opt.backboneLength, _opt.useAlaAtCTerminus);
	string variablePositionString = generateString("0", _opt.backboneLength);
	string rotamerLevels = generateString("0", _opt.backboneLength);

	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	int highestRotamerLevel = numberOfRotamerLevels-1;
	//make into a function
	int numAAs = 0;
	double lvlavg = 0;
	vector<double> avgs;
	//cout << "lvl " << levelCounter << endl;
	for (uint i = 0; i < resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(resiBurial.size());
		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
			levelCounter++;
			lvlavg = lvlavg/numAAs;
			avgs.push_back(lvlavg);
			lvlavg=0;
			numAAs=0;
			//cout << "lvl " << levelCounter << endl;
		}
		//cout << resiBurial[i].first << ": " << resiBurial[i].second << endl;
		lvlavg = lvlavg+resiBurial[i].second;
		numAAs++;
		int backbonePosition = resiBurial[i].first;
		Position &pos = sys.getPosition(backbonePosition);
		string posRot = _opt.sasaRepackLevel[levelCounter];
		int resiNum = pos.getResidueNumber();
		int posNum = resiNum-_opt.thread;
		string add;

		if (levelCounter < _opt.interfaceLevel){
			add = "Add all Ids at this pos";
			interfacePositions.push_back(resiNum);
			if (backbonePosition > 2 && backbonePosition < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 3 and 4 here instead of 4 and 5 to prevent changes at the interface like others
				variablePositionString.replace(variablePositionString.begin()+posNum, variablePositionString.begin()+posNum+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		} else {
			add = "Only 1 ID";
		}
		rotamerLevels.replace(rotamerLevels.begin()+posNum, rotamerLevels.begin()+posNum+1, MslTools::intToString(levelCounter));
	}
	lvlavg = lvlavg/numAAs;
	avgs.push_back(lvlavg);
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.Ids, interfacePositions);
	PolymerSequence PS(polySeq);

	vector<int> rotamerSamplingPerPosition = getRotamerSampling(rotamerLevels);
	vector<int> linkedPositions = getLinkedPositions(rotamerSamplingPerPosition, _opt.interfaceLevel, highestRotamerLevel);

	//String for the positions of the sequences that are considered interface for positions amd high rotamers
	string rotamerSamplingString = getInterfaceString(rotamerSamplingPerPosition, _opt.backboneLength);

	// Vector for linked positions in "A,25 B,25" format


	// Define referenced output variables
	_rotamerLevels = rotamerLevels;
	_rotamerSamplingPerPosition = rotamerSamplingPerPosition;
	_variablePositionString = variablePositionString;
	_rotamerSamplingString = rotamerSamplingString;
	_linkedPositions = linkedPositions;
	_out << endl;
	_out << "PolyLeu Backbone:   " << backboneSeq << endl;
	_out << "Variable Positions: " << variablePositionString << endl;
	_out << "Rotamers Levels:    " << rotamerSamplingString << endl;
	_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition);
	_allInterfacePositions = getAllInterfacePositions(_opt, rotamerSamplingPerPosition);

	//cout << endl << "Averages:\t";
	//for (uint a=0; a<avgs.size(); a++){
	//	cout << avgs[a] << "\t";
	//}
	//cout << endl;

	//for (uint i=0; i<_interfacePositions.size(); i++){
	//	cout << _interfacePositions[i] << ",";
	//}
	//cout << endl;
	cout << endl;
	cout << "PolyLeu Backbone:   " << backboneSeq << endl;
	cout << "Variable Positions: " << variablePositionString << endl;
	cout << "Rotamers Levels:    " << rotamerSamplingString << endl;

	int numPosAtInterface = 0;
	for(string::iterator it = variablePositionString.begin(); it != variablePositionString.end();it++ ) {
		stringstream ss;
		ss << *it;
		string num = ss.str();
		if (MslTools::toInt(num) == 1){
			numPosAtInterface++;
		}
	}
	return PS;
}

void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, CartesianPoint &_xAxis,
 CartesianPoint &_zAxis, Transforms &_trans) {
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	_sys.readPdb(_opt.backboneFile);

	//// read the initial helical axis coordinates	
	//PDBReader readAxis;
	//if(!readAxis.read(_axis)) {
	//	cout << "Unable to read axis" << endl;
	//	exit(0);
	//}

	//// setup the helical axis
	//System helicalAxis;
	//helicalAxis.addAtoms(readAxis.getAtomPointers());
	
	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, _axisA, _axisB, _ori, _xAxis, _zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

// sets the gly69 backbone to starting geometry

// help functions
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputErrorMessage(Options &_opt){
	cout << endl;
	cout << "The program terminated with errors:" << endl;
	cout << endl;
	cout << _opt.errorMessages << endl;
	cout << endl;
	cout << _opt.OPerrors << endl;
	usage();
}

void help(Options defaults) {
	cout << "This program runs as:" << endl;
	cout << " % seqDesign " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --numberOfStructuresToMCRepack <int> --energyCutOff <double> --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --MCStartTemp <double> --MCEndTemp <double> --MCCurve <CONSTANT-0, LINEAR-1, EXPONENTIAL-2, SIGMOIDAL-3, SOFT-4>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --thread <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   --weight_hbond <double> --weight_vdw <double> --weight_solv <double> --weight_seqEntropy <double>" << endl;
	cout << "   --sasaRepackLevel <rotLevel> (in format SL95.00; 4 levels used by default) --interfaceLevel <int> " << endl << endl;
	cout << "Template Configuration file (copy and paste the below into a file.config and run code as bin/seqDesign --config file.config" << endl;
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "pdbOutputDir " << defaults.pdbOutputDir << endl;
	cout << setw(20) << "tmStart " << defaults.tmStart << endl;
	cout << setw(20) << "tmEnd " << defaults.tmEnd << endl;

	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "geometryDensityFile " << defaults.geometryDensityFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "backboneCrd " << defaults.backboneCrd << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;
	cout << setw(20) << "backboneFile " << defaults.backboneFile << endl;
	cout << setw(20) << "selfEnergyFile " << defaults.selfEnergyFile << endl;
	cout << setw(20) << "pairEnergyFile " << defaults.pairEnergyFile << endl;
	cout << setw(20) << "AACompositionPenaltyFile " << defaults.AACompositionPenaltyFile << endl << endl;

	cout << "#Geometry and Transformation parameters" << endl;
	cout << setw(20) << "xShift" << defaults.xShift << endl;
	cout << setw(20) << "crossingAngle" << defaults.crossingAngle << endl;
	cout << setw(20) << "axialRotation" << defaults.axialRotation << endl;
	cout << setw(20) << "zShift" << defaults.zShift << endl;
	cout << setw(20) << "thread" << defaults.thread << endl;
	cout << setw(20) << "backboneLength " << defaults.backboneLength << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;
	cout << setw(20) << "deleteTerminalHbonds" << defaults.deleteTerminalHbonds << endl;
	cout << setw(20) << "useSasa" << defaults.useSasa << endl;
	cout << setw(20) << "getGeoFromPDBData" << false << endl;//Since we already have the geometry output here, default to false in the rerun config
	cout << setw(20) << "runDEESingles" << defaults.runDEESingles << endl;
	cout << setw(20) << "runDEEPairs" << defaults.runDEEPairs << endl;
	cout << setw(20) << "runSCMF" << defaults.runSCMF << endl;
	//TODO: set this up so that instead of running through that, I just have the output sequence and state below
	if (defaults.runDEESingles == false && defaults.runDEEPairs == false && defaults.runSCMF == false){

	}

	if (defaults.useSasa == true){
		//cout << "#Load Rotamers based on SASA scores" << endl;
		for (uint i=0; i<defaults.sasaRepackLevel.size()-1; i++){
			cout << setw(20) << "sasaRepackLevel" << defaults.sasaRepackLevel[i] << endl;
		}
		cout << setw(20) << "interfaceLevel" << defaults.interfaceLevel << endl;
	} else {
		//cout << "#Load Rotamers by interface and non interfacial positions" << endl;
		cout << setw(20) << "SL" << defaults.SL << endl;
		cout << setw(20) << "SLInterface" << defaults.SLInterface << endl;
	}

	cout << "#MonteCarlo Paramenters" << endl;
	cout << setw(20) << "MCCycles" << defaults.MCCycles << endl;
	cout << setw(20) << "MCMaxRejects" << defaults.MCMaxRejects << endl;
	cout << setw(20) << "MCStartTemp" << defaults.MCStartTemp << endl;
	cout << setw(20) << "MCEndTemp" << defaults.MCEndTemp << endl;
	cout << setw(20) << "MCCurve" << defaults.MCCurve << endl;
	cout << setw(20) << "greedyCycles" << defaults.greedyCycles << endl;

	cout << "#Alternate IDs" << endl;
	for (uint i=0; i<defaults.Ids.size()-1; i++){
		cout << setw(20) << defaults.Ids[i] << endl;
	}

	cout << endl << "#Rerun Seed" << endl;
	cout << setw(20) << "seed" << defaults.seed << endl;

	cout << endl << "#Energy term weights" << endl;
	cout << setw(20) << "weight_vdw " << defaults.weight_vdw << endl;
	cout << setw(20) << "weight_hbond " << defaults.weight_hbond << endl;
	cout << setw(20) << "weight_solv " << defaults.weight_solv << endl;
	cout << setw(20) << "weight_seqEntropy " << defaults.weight_seqEntropy << endl;
	cout << endl;
}

void localBackboneRepack(Options &_opt, System &_startGeom, string _sequence, uint _rep, double _savedXShift, System &_helicalAxis, AtomPointerVector &_axisA,
 AtomPointerVector &_axisB, vector<int> _rotamerSampling, Transforms &_trans, RandomNumberGenerator &_RNG, ofstream &_out){
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PS(polySeq);
	
	// set up the system for the input sequence
	System sys;
	prepareSystem(_opt, sys, _startGeom, PS);
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	//decided to try loading low number of rotamers for this instead of high
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);

	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// build in baseline as an estimate for monomer energies
	if (_opt.useBaseline){
		//addBaselineToSelfPairManager();//TODO: was going to make this a function but I feel like it already exists in spm, so going to wait when I have time to look that up
		//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
		buildBaselines(sys, _opt);
	}

	// set the system to the original xShift state
	//sys.applySavedCoor("savedBestState");
	//_helicalAxis.applySavedCoor("BestAxis");

	//// save the current state as the repack state
	//sys.saveAltCoor("savedRepackState");
	//_helicalAxis.saveAltCoor("BestRepack");
	
	// get the best repack energy
	double prevBestEnergy = sys.calcEnergy();

	// Setup time variables
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);

	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	//spm.updateWeights();
	spm.setOnTheFly(false);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	spm.runGreedyOptimizer(_opt.greedyCycles);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	_out << "Time to calculate energies for backbone repack: " << diffTime << endl;
	//TODO: seg faults below	
	// do backbone geometry repacks
	//for (int i=0; i < _opt.numRepacks; i++){
		if (_opt.verbose){
			cout << "===============================================" << endl;
			cout << "Performing Local Monte Carlo Backbone Repack   " << endl;
			cout << "===============================================" << endl;
		}
		//monteCarloRepack(_opt, sys, _savedXShift, spm, _helicalAxis, _axisA, _axisB, apvChainA, apvChainB, _trans, _RNG, prevBestEnergy, 
		//i, _out);
		monteCarloRepack(_opt, sys, _savedXShift, spm, _helicalAxis, _axisA, _axisB, apvChainA, apvChainB, _trans, _RNG, prevBestEnergy, 
		_rep, _out);
		// Initialize PDBWriter
		PDBWriter writer;
		writer.open(_opt.pdbOutputDir + "/backboneOptimized_" + to_string(_rep) + ".pdb");
		writer.write(sys.getAtomPointers(), true, false, true);
		writer.close();
	//}
	cout << "Backbone repack Worked!" << endl;
	// TODO: set the startGeom to the repacked geom at sys
	// assign the coordinates of our system to the given geometry 
	_startGeom.assignCoordinates(sys.getAtomPointers(),false);
	_startGeom.buildAllAtoms();

     //outputs a pdb file for the structure 
	
}

void monteCarloRepack(Options &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 uint _rep, ofstream &_out){
	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	// starting geometry
	double xShift = _savedXShift;	
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	// Monte Carlo Repack Manager Setup
	//MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects, _opt.backboneConvergedSteps, _opt.backboneConvergedE);
	MCMngr.setEner(_prevBestEnergy);

	vector<uint> startStateVec = _spm.getMinStates()[0];
	vector<unsigned int> MCOBest = startStateVec;
		
	unsigned int counter = 0;
	double currentEnergy = _prevBestEnergy;
	double startDimer = _prevBestEnergy;

	PDBWriter writer;
	// loop through the MC cycles for backbone repacks
	while(!MCMngr.getComplete()) {
		
		_sys.applySavedCoor("savedRepackState");
		_helicalAxis.applySavedCoor("BestRepack");
		
		int moveToPreform = _RNG.getRandomInt(3);
		
		double deltaXShift = 0.0;
		double deltaZShift = 0.0;
		double deltaCrossingAngle = 0.0;
		double deltaAxialRotation = 0.0; 
		
		//======================================
		//====== Z Shift (Crossing Point) ======
		//======================================
		if (moveToPreform == 0) {
			//deltaZShift = getStandardNormal(RNG1) * 0.1;
			deltaZShift = getStandardNormal(_RNG) * _opt.deltaZ;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(_RNG) * _opt.deltaAx;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * _opt.deltaCross;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * _opt.deltaX;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		}
		
		// Run _optimization
		repackSideChains(_spm, 10);

		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0];
		
		if (!MCMngr.accept(currentEnergy)) {
			if (_opt.verbose){
				cout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy << endl;
			}
		} else {
			_prevBestEnergy = currentEnergy;
			_sys.saveAltCoor("savedRepackState");
			_helicalAxis.saveAltCoor("BestRepack");
		
			xShift = xShift + deltaXShift;
			crossingAngle = crossingAngle + deltaCrossingAngle;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift + deltaZShift;
			MCOBest = MCOFinal;
		
			if (_opt.verbose){
				cout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy << endl;
			}
			counter++;
			writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	writer.close();
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
		
	_sys.applySavedCoor("savedRepackState");
	// TODO: make this into a map that saves all of these to be output
	double baseline = _spm.getStateEnergy(MCOBest, "BASELINE")+_spm.getStateEnergy(MCOBest, "BASELINE_PAIR");
	double dimerEnergy = _spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-baseline;
	double vdw = _spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	double hbond = _spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	double imm1 = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+_spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	double dimerDiff = dimerEnergy-startDimer;
	cout << "Energy #" << _rep << ": " << finalEnergy << endl;
	// calculate the solvent accessible surface area
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	// Output change in geometry
	_out << "***STARTING GEOMETRY***" << endl;
	_out << "xShift:        " << _opt.xShift << endl;
	_out << "crossingAngle: " << _opt.crossingAngle << endl;
	_out << "axialRotation: " << _opt.axialRotation << endl;
	_out << "zShift:        " << _opt.zShift << endl << endl;

	// sets the new backbone parameters
	_opt.xShift = xShift;
	_opt.crossingAngle = crossingAngle;
	_opt.axialRotation = axialRotation;
	_opt.zShift = zShift;
	
	_out << "***AFTER REPACK GEOMETRY***" << endl;
	_out << "xShift:        " << _opt.xShift << endl;
	_out << "crossingAngle: " << _opt.crossingAngle << endl;
	_out << "axialRotation: " << _opt.axialRotation << endl;
	_out << "zShift:        " << _opt.zShift << endl << endl;
	
	// Print out info to the summary csv file
	//_out << to_string(_rep) << ',' << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
	//_out << dimerDiff << ',' << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',';
	//_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
	//_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
	//cout << to_string(_rep) << ',' << _opt.sequence << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
	//cout << dimerDiff << ',' << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',';
	//cout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
	//cout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
}

double computeMonomerEnergy(System & _sys, Options & _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout){

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	// Objects used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";

	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	int firstPos = 0;
    int lastPos = monoSys.positionSize();
	deleteTerminalHydrogenBondInteractions(monoSys, firstPos, lastPos);

	/******************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	SelfPairManager monoSpm;
	monoSpm.seed(_opt.seed);
	monoSpm.setVerbose(_opt.verbose);

	for (uint k=0; k < monoSys.positionSize(); k++) {
		Position &pos = monoSys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!monoRot.loadRotamers(&pos, pos.getResidueName(), "SL97.00")) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	_mout << "Monomer - VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;

	monoSpm.setSystem(&monoSys);
	monoSpm.updateWeights();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	trans.translate(chainA, interDistVect);


	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	trans.translate(chainA, move5Down);
	double bestZ = -5.0;

	monoSys.calcEnergy();

	// Repack side chains
	monoSpm.setOnTheFly(1);
	monoSpm.calculateEnergies();
    monoSpm.runGreedyOptimizer(_opt.greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	helicalAxis.saveAltCoor("BestMonomerAxis");
	_mout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		_mout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		for(int j=0; j<=3; j++) { // test at 4 rotations 0, 90, 180 and 270 degrees
			//==================================
			//====== Axial Rot ======
			//==================================
			monoAxialRotation = j * 90.0;

			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			_mout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			//monoSys.writePdb("mono_" + MslTools::doubleToString(monoTilt) + "_" + MslTools::doubleToString(monoAxialRotation) + ".pdb");

			if(currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				bestTilt = monoTilt;
				bestRotation = monoAxialRotation;
				monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
				monoSys.saveAltCoor("savedBestMonomer");
				helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

		}
	}
	
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	// MonteCarloManager MCMngr(1000.0, 0.5, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		helicalAxis.applySavedCoor("BestMonomerAxis");

		int moveToPreform = _RNG.getRandomInt(2);

		double deltaZShift = 0.0;
		double deltaTilt = 0.0;
		double deltaAxialRotation = 0.0;

		//======================================
		//====== Z Shift ======
		//======================================
		if (moveToPreform == 0) {
			deltaZShift = getStandardNormal(_RNG) * 1.0;
			CartesianPoint translateA = axisA(1).getCoor() - axisA(0).getCoor(); // vector minus helical center
			translateA = translateA.getUnit() * deltaZShift; // unit vector of helical _axis times the amount to shift by
			trans.translate(chainA, translateA);
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_mout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_mout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_mout << "state rejected   energy: " << currentEnergy << endl;
		}
		else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			_mout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}

		counter++;
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	_mout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");
	_mout << monoEset->getSummary();
	_mout << endl;

	// print the monomer
	string monoOutCrdFile  = _opt.pdbOutputDir + "/monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}

	string monoOutPdbFile  = _opt.pdbOutputDir + "/monomer.pdb";
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	// Store monomer energy by term
	monoSys.calcEnergy();
	_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

	for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
		_mout << it->first << " " << it->second << endl;
	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;
}
