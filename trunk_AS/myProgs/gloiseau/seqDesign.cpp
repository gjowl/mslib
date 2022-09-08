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
#include "designOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "seqDesign";//TODO: better name
string programDescription = "Designs sequences for backbone geometries extracted from the PDB, optimizing specifically for vdW energies";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "18 August 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime, spmTime;
auto start = chrono::system_clock::now();

// Functions
/*
	I have left the most important functions within this code in case they need to be looked at or changed.
	Many of the auxiliary functions are found within the designFunctions.cpp file.
*/
// sequence search functions
void stateMCUnlinked(System &_sys, Options &_opt, PolymerSequence &_PS,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
 vector<unsigned int> &_bestState, string &_bestSequence, vector<string> &_seqs, vector<string> &_allSeqs,
 map<string,vector<uint>> &_sequenceVectorMap, vector<uint> &_allInterfacialPositionsList,
 vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, int _rep, ofstream &_out, ofstream &_err);
map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, vector<uint> _bestState, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap,
 map<string,double> _sequenceEntropyMap, vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList,
 vector<int> _rotamerSampling);
string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap);
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, vector<uint> _prevVec, double _prevEnergy, 
 string _currSeq, vector<uint> _currVec, vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap);
void searchForBestSequencesUsingThreads(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, map<string, vector<uint>> &_sequenceVectorMap,
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<int> &_rotamerSampling, int _rep, ofstream &_out, ofstream &_err);

// geometry setup functions
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
 string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
 CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
void defineRotamerLevelsByResidueBurial(System &_sys, Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString);

// backbone repack functions
void localBackboneRepack(Options &_opt, System &_startGeom, string _sequence, uint _rep, double _savedXShift, System &_helicalAxis, AtomPointerVector &_axisA,
 AtomPointerVector &_axisB, vector<int> _rotamerSampling, Transforms &_trans, RandomNumberGenerator &_RNG, ofstream &_out);
void monteCarloRepack(Options &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 uint _rep, ofstream &_out);
void getCurrentMoveSizes(Options &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);

// output functions
void outputFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, map<string,map<string,double>> _sequenceEnergyMap,
 vector<double> _densities);
string getRunParameters(Options &_opt, vector<double> _densities);

void outputTime(auto _start, string _descriptor, bool _beginFunction, ofstream &_out){
	auto end = chrono::system_clock::now();
	time_t endTimeFormatted = chrono::system_clock::to_time_t(end); 
	chrono::duration<double> elapsedTime = end-_start;
	if (_beginFunction == true){
		_out << _descriptor << " started. Time: " << ctime(&endTimeFormatted);
		_out << "Elapsed time of program: " << elapsedTime.count() << "s/" << elapsedTime.count()/60 << "min" << endl << endl;
		cout << _descriptor << " started. Time: " << ctime(&endTimeFormatted);
		cout << "Elapsed time of program: " << elapsedTime.count() << "s/" << elapsedTime.count()/60 << "min" << endl << endl;
	} else {
		_out << _descriptor << " finished. Time: " << ctime(&endTimeFormatted);
		_out << "Elapsed time of program: " << elapsedTime.count() << "s/" << elapsedTime.count()/60 << "min" << endl << endl;
		cout << _descriptor << " finished. Time: " << ctime(&endTimeFormatted);
		cout << "Elapsed time of program: " << elapsedTime.count() << "s/" << elapsedTime.count()/60 << "min" << endl << endl;
	}
}

void useInputInterface(Options &_opt, string &_variablePositionString, string &_rotamerLevels, vector<int> &_interfacePositions){
	for (uint i=3; i<_opt.interface.length()-4; i++){
		if (_opt.interface[i] == '0'){
			if (_variablePositionString[i] != '0'){
				_variablePositionString.replace(_variablePositionString.begin()+i, _variablePositionString.begin()+i+1, "0");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		} else if (_opt.interface[i] == '1'){
			if (_rotamerLevels[i] != '0'){
				_rotamerLevels.replace(_rotamerLevels.begin()+i, _rotamerLevels.begin()+i+1, "0");
			}
			if (_variablePositionString[i] == '0'){
				_variablePositionString.replace(_variablePositionString.begin()+i, _variablePositionString.begin()+i+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		}
	}
	for (uint i=0; i<_variablePositionString.length(); i++){
		if (_variablePositionString[i] == '1'){
			_interfacePositions.push_back(i+_opt.thread);
		}
	}
}

// help functions
void usage();
void help(Options defaults);
void outputErrorMessage(Options &_opt);

/******************************************
 *
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){
	time_t startRealTime = chrono::system_clock::to_time_t(start); 
	cout << "Program Start time: " << ctime(&startRealTime) << endl;
	
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
	Options opt = parseOptions(argc, argv);
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
	ofstream sout; // summary file output
	ofstream err; // error file output
	ofstream rerun; // rerun config output

	setupDesignDirectory(opt, date); // makes a directory in the directory that you run from, with design_<runNumber> as the name

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

	sout << date << endl; // output the date and time to the summary file

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
		if (opt.getRandomAxRotAndZShift){
			getAxialRotAndZShift(opt, RNG, densities, sout);
			cout << "***STARTING GEOMETRY:***" << endl;
			cout << "xShift:        " << opt.xShift << endl;
			cout << "crossingAngle: " << opt.crossingAngle << endl;
			cout << "axialRotation: " << opt.axialRotation << "\tDensity: " << densities[1] << endl;
			cout << "zShift:        " << opt.zShift << "\tDensity: " << densities[2] << endl << endl;
		} else {
			cout << "***STARTING GEOMETRY:***" << endl;
			cout << "xShift:        " << opt.xShift << endl;
			cout << "crossingAngle: " << opt.crossingAngle << endl;
			cout << "axialRotation: " << opt.axialRotation << endl;
			cout << "zShift:        " << opt.zShift << endl << endl;
			densities.push_back(0);
			densities.push_back(0);
			densities.push_back(0);
		}
	}

	//String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.Ids);
	cout << "Amino acids for design: " << alternateIds << endl;

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
	// get the starting geometry using polyglycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,axisA,axisB,ori,xAxis,zAxis,trans);
	//PDBWriter writer;
	//writer.open(opt.pdbOutputDir + "/startGeom.pdb");
	//writer.write(startGeom.getAtomPointers(), true, false, true);
	//writer.close();
	
	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	string rotamerLevels; // string of rotamer levels for each position
	string variablePositionString; // string of variable positions for each position
	string rotamerSamplingString; // string of rotamer sampling number (if 4 rotamer levels, 0-3 for each position)
	vector<int> linkedPositions; // vector of the positions that will be linked
	vector<uint> interfacePositions; // vector of positions at the interface excluding termini positions
	vector<uint> allInterfacePositions; // vector of positions at the interface including the terminal positions
	vector<int> rotamerSamplingPerPosition; // vector of rotamer level for each position
	
	// Defines the interfacial positions and the number of rotamers to give each position
	// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
	// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
	// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
	PolymerSequence interfacePolySeq = getInterfacialPolymerSequence(opt, startGeom, PS, rotamerLevels, variablePositionString, rotamerSamplingString,
	 linkedPositions, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, sout);

	outputTime(start, "Identify Interface", false, sout);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	prepareSystem(opt, sys, startGeom, interfacePolySeq);
	checkIfAtomsAreBuilt(sys, err); // check to verify that all atoms have coordinates
	// redacted on 2022-9-1: printed the pdbs of the helical axis and doing this moves
	// it a small bit more down that is unnecessary and slightly out of membrane 
	//moveZCenterOfCAMassToOrigin(sys.getAllAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	Chain &chainA = sys.getChain("A");
	int seqLength = chainA.positionSize();
	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	// add in an estimate of monomer energy for each sequence (calculate using only CHARMM_VDW, SCWRL4_HBOND, CHARMM_IMM1, and CHARMM_IMM1REF)
	if (opt.useBaseline){
		buildBaselines(sys, opt);
	}

	// link the interfacial positions (for quicker calculation of initial sequence for homodimers)
	//if (opt.linkInterfacialPositions && !opt.hetero){
	if (opt.linkInterfacialPositions){
		vector<vector<string>> linkedPos = convertToLinkedFormat(sys, linkedPositions, seqLength);
		sys.setLinkedPositions(linkedPos);
	}
	
	// initialize the object for loading rotamers into system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// load rotamers for each amino acid at each position into the system
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	//CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs
	/******************************************************************************
	 *           === VARIABLES FOR SAVING ENERGIES AND SEQUENCES ===
	 ******************************************************************************/
	vector<string> seqs; // vector to hold designed sequences
	vector<string> allSeqs; //vector to hold all designed sequences from the stateMC
	map<string, map<string,double>> sequenceEnergyMapBest; // energyMap to hold all energies for output into a summary file
	map<string, vector<uint>> sequenceVectorMap; //sequence and vector map: sequence will be tied vector with AAs and rotamers
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile); // Get sequence entropy map

	/******************************************************************************
	 *                        === SETUP SPM AND RUN SCMF ===
	 ******************************************************************************/
	//TODO: make it so that this part is optional: if I submit a sequence to start, don't even do this. Just find the interface, greedy for best sequence, and continue
	outputTime(start, "Self Consistent Mean Field", true, sout);
	vector<uint> bestState = runSCMFToGetStartingSequence(sys, opt, RNG, rotamerSamplingString, variablePositionString,
	 seqs, allInterfacePositions, sequenceEnergyMapBest, sequenceVectorMap, sequenceEntropyMap, sout);
	outputTime(start, "Self Consistent Mean Field", false, sout);
	
	// set system to the best sequence or input sequence 
	sys.setActiveRotamers(bestState);
	double bestEnergy = sys.calcEnergy();
	sys.getEnergySet()->eraseTerm("CHARMM_VDW");
	sys.getEnergySet()->eraseTerm("CHARMM_IMM1");
	sys.getEnergySet()->eraseTerm("CHARMM_IMM1REF");
	sys.getEnergySet()->eraseTerm("SCWRL4_HBOND");

	/******************************************************************************
	 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
	 ******************************************************************************/
	// Unlink the best state from SCMF if not using linked positions during the state Monte Carlo
	if(opt.linkInterfacialPositions){
		unlinkBestState(opt, bestState, rotamerSamplingPerPosition, seqLength);
	}
	string bestSequence;
	for (uint i=0; i<3; i++){
		outputTime(start, "Sequence search replicate " + to_string(i), true, sout);
		stateMCUnlinked(sys, opt, interfacePolySeq, sequenceEnergyMapBest, sequenceEntropyMap, bestState, bestSequence, seqs, allSeqs,
			sequenceVectorMap, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, i, sout, err);
		outputTime(start, "Sequence search replicate " + to_string(i), false, sout);
		outputTime(start, "Backbone repack replicate " + to_string(i), true, sout);
		localBackboneRepack(opt, sys, bestSequence, i, opt.xShift, helicalAxis, axisA, axisB, rotamerSamplingPerPosition, trans, RNG, sout);
		outputTime(start, "Backbone repack replicate " + to_string(i), false, sout);
	}
	// TODO: make a decision; should I try to calculate the energies for every sequence on the final backbone? Or should I just
	// say what the geometry is?
	// I think first I need to look at the geometries, the pds, etc. to make sure it looks reasonable
	/******************************************************************************
	 *            === CALCULATE MONOMER ENERGIES OF EACH SEQUENCE ===
	 ******************************************************************************/
	outputTime(start, "Monomer energy calculations ", true, sout);
	computeMonomerEnergies(opt, trans, sequenceEnergyMapBest, seqs, RNG, sout, err);
	outputTime(start, "Monomer energy calculations ", false, sout);
	//getSasaDifference(sequenceStatePair, sequenceEnergyMap);

	///******************************************************************************
	// *              === CALCULATE TOTAL ENERGIES AND WRITE PDBS ===
	// ******************************************************************************/
	//// Initialize PDBWriter for designs
	//PDBWriter writer;
	//writer.open(opt.pdbOutputDir + "/allDesigns.pdb");
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
	ss <<  _densities[0] << t << _densities[1] << t << _densities[2] << t << _opt.thread << t << _opt.sasaRepackLevel.size() << t << _opt.interfaceLevel << t << _opt.backboneLength;
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
		seqLine << sequence << t << _interface << t << interfaceSequence << t;
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
		string seqNumber = MslTools::doubleToString(energyMap.at("SequenceNumber"));
		seqLine << seqNumber << t << runParameters;
		string line = seqLine.str();
		energyLines.push_back(line);
		i++;
	}
	ofstream eout;
	string eoutfile = _opt.pdbOutputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	eout << "Sequence" << t << "Interface" << t << "InterfaceSequence" << t;
	eout << enerTerms.str();
	eout << "angleDistDensity" << t << "axialRotationDensity" << t << "zShiftDensity" << t << "repackLevels" << t << "interfaceLevels" << t << "backboneLength" << endl;
	cout << "Sequence" << t << "Interface" << t << "InterfaceSequence" << t;
	cout << enerTerms.str();
	eout << "angleDistDensity" << t << "axialRotationDensity" << t << "zShiftDensity" << t << "repackLevels" << t << "interfaceLevels" << t << "backboneLength" << endl;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
		cout << energyLines[i] << endl;
	}
	eout.close();
}

map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, vector<uint> _bestState, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap,
 map<string,double> _sequenceEntropyMap, vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList,
 vector<int> _rotamerSampling){
	// Get a random integer to pick through the variable positions
	int rand = _RNG.getRandomInt(0, _interfacialPositionsList.size()-1);
	int interfacePosA = _interfacialPositionsList[rand];
	int interfacePosB = interfacePosA+_bestSeq.length();

	// Get the random position from the system
	Position &randPosA = _sys.getPosition(interfacePosA);
	Position &randPosB = _sys.getPosition(interfacePosB);
	string posIdA = randPosA.getPositionId();
	string posIdB = randPosB.getPositionId();

	// variable setup for current state
	map<string,map<string,double>> sequenceEnergyMap;
	vector<thread> threads;
	for (uint i=0; i<_opt.Ids.size(); i++){
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
			threads.push_back(thread{energyFunction, ref(_opt), ref(_spm), _bestSeq, _bestState, _bestEnergy, currSeq, currVec, ref(_rotamerSampling), ref(_allInterfacialPositionsList), 
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
			currEnergyComparison = seq.second["currEnergyTotal"];
			i++;
		} else {
			if (seq.second["currEnergyTotal"] < currEnergyComparison){
				string prevSeq = currSeq;
				currSeq = seq.first;
				currEnergyComparison = seq.second["currEnergyTotal"];
			}
		}
	}
	return currSeq;
}

// threaded function that gets energies for a sequence and appends to a map
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, vector<uint> _prevVec, double _prevEnergy, 
 string _currSeq, vector<uint> _currVec, vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap){
	// variable setup
	map<string,double> energyMap;

	// Compute dimer energy
	outputEnergiesByTerm(_spm, _currVec, energyMap, _opt.energyTermList, "Dimer", true);
	double currEnergy = _spm.getStateEnergy(_currVec);

	// initialize variables for this thread
	double currEnergyTotal = 0;
	double bestEnergyTotal = 0;
	double currSEProb = 0;
	double prevSEProb = 0;
	double currEntropy = 0;
	double prevEntropy = 0;

	//TODO: I just realized that for heterodimers, I may need to completely remake some of these functions as with the below only taking the sequence of one helix; I may make a hetero and homo functions list?
	calculateInterfaceSequenceEntropy(_opt, _prevSeq, _currSeq, _sequenceEntropyMap, prevSEProb,
	 currSEProb, prevEntropy, currEntropy, _prevEnergy, currEnergy, bestEnergyTotal,
	 currEnergyTotal, _allInterfacePositions);

	// output info
	outputEnergiesByTerm(_spm, _currVec, energyMap, _opt.energyTermList, "Dimer", true);
	double baseline = _spm.getStateEnergy(_currVec, "BASELINE")+_spm.getStateEnergy(_currVec, "BASELINE_PAIR");
	double enerAndSeqEntropy = bestEnergyTotal-currEnergyTotal;
	energyMap["Dimerw/Baseline"] = currEnergy;
	energyMap["Dimer"] = currEnergy-baseline;
	energyMap["Baseline"] = baseline;
	energyMap["energyComparison"] = enerAndSeqEntropy;
	energyMap["SequenceProbability"] = currSEProb;
	energyMap["bestEnergyTotal"] = bestEnergyTotal;
	energyMap["currEnergyTotal"] = currEnergyTotal;
	energyMap["entropyDiff"] = prevEntropy-currEntropy;
	energyMap["currEntropy"] = currEntropy;
	energyMap["prevEntropy"] = prevEntropy;
	_seqEnergyMap[_currSeq] = energyMap;
}

// function to prepare the system for design:
// - sets up CharmmSystemBuilder
// - assigns to given coordinates
// - Adds in the correct energies for the EnergySet
// - deletes terminal hydrogen bonds
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
	//Chain & chainA = _sys.getChain("A");
	//Chain & chainB = _sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	//AtomPointerVector & apvChainA = chainA.getAtomPointers();
	//AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
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
	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);

	if (_opt.useElec == true){
		Eset->setTermActive("CHARMM_ELEC", true);
		Eset->setWeight("CHARMM_ELEC", _opt.weight_elec);
	}

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	int firstPos = 0;// should this be based on the thread and the
    int lastPos = _sys.positionSize();
    deleteTerminalBondInteractions(_sys,_opt,firstPos,lastPos);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	CSB.updateNonBonded(10,12,50);
	// TODO: maybe calculate the energy here, then see if there's clashing. If so, then move helices away until no clashing, keeping the same
	// other coordinates? Just for simplicity for now. And if I want to implement this in design, adding this in will likely be a good idea.	
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
	
	// add in baseline energies
	if (_opt.useBaseline){
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
	outputTime(start, "Self and Pair of all amino acids calculation " + to_string(_rep), true, _out);
	spm.calculateEnergies();
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	_out << "Time to calculate energies: " << diffTime  << "s" << endl;
	outputTime(start, "Self and Pair of all amino acids calculation " + to_string(_rep), false, _out);

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
	//MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects, _opt.MCConvergedSteps, _opt.backboneConvergedE);
	MonteCarloManager MC(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MC.setRandomNumberGenerator(&_RNG);

	// Start from most probable state
	_sys.setActiveRotamers(_bestState);
	double bestEnergy = _spm.getStateEnergy(_bestState);

	// State variable setup
	vector<unsigned int> prevStateVec = _bestState;
	MC.setEner(bestEnergy);

	// Alternate AA Ids for each of the interfacial positions
	vector<string> ids = _opt.Ids;

	// Variables setup for MC while loop
	map<double, string> sequences;
	int cycleCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);

	// Sequence Search Energy Landscape file
	ofstream lout;
	string loutfile  = _opt.pdbOutputDir + "/sequenceSearchEnergyLandscape_" + to_string(_rep) + ".out";
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
	
	string bestSeq = prevStateSeq;
	map<string,map<string,double>> allSequenceEnergyMap;
	//cout << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	//_out << "Finding " << _opt.numStatesToSave << " sequences using membrane composition (State MonteCarlo)..." << endl;
	// Monte Carlo while loop for finding the best sequences
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "" << endl;
			cout << "Starting Seq: " << prevStateSeq << endl;
		}
		// get the sequence entropy probability for the current best sequence
		map<string,vector<uint>> sequenceVectorMap;
		map<string,map<string,double>> sequenceEnergyMap = mutateRandomPosition(_sys, _opt, _spm, _RNG, bestSeq, prevStateVec, bestEnergy, 
		 sequenceVectorMap, _sequenceEntropyMap, _allInterfacialPositionsList, _interfacialPositionsList, _rotamerSampling);
		
		// get the best sequence and energy for the current mutation position (picks sequence with energy including vdw, hbond, imm1, baseline, sequence entropy)
		string currSeq = getBestSequenceInMap(sequenceEnergyMap);
		// TODO: add in check step to see if sequence is already found in the best sequence list and skip if so
		double currEnergyTotal = sequenceEnergyMap[currSeq]["currEnergyTotal"]; // energy total for current sequence
		double bestEnergyTotal = sequenceEnergyMap[currSeq]["bestEnergyTotal"]; // energy total for previous sequence (details in energyFunction)
		double prevStateEntropy = sequenceEnergyMap[currSeq]["prevEntropy"];
		double currStateEntropy = sequenceEnergyMap[currSeq]["currEntropy"];
		vector<uint> currStateVec = sequenceVectorMap[currSeq]; // current state vector for current sequence (best rotamers and sequence identities)
		MC.setEner(bestEnergyTotal);
		
		// MC accept and reject conditions
		if (!MC.accept(currEnergyTotal)){
			_sys.setActiveRotamers(prevStateVec); // set rotamers to the previous state
			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			_sys.setActiveRotamers(currStateVec); // set rotamers to the current state
			if (_opt.energyLandscape){
				map<string,double> energyMap = sequenceEnergyMap[currSeq];
				lout << prevStateSeq << "\t" << currSeq << "\t" << bestEnergyTotal << "\t" << currEnergyTotal << "\t";
				lout << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			}
			prevStateVec = currStateVec; // set the previous state vector to be the current state vector
			prevStateSeq = currSeq; // set the previous sequence to be the current sequence
			bestSeq = currSeq; // set the best sequence to the newly accepted current sequence
			bestEnergy = sequenceEnergyMap[bestSeq]["Dimerw/Baseline"]; // set the best energy to the current energy (vdw, hbond, imm1, baseline)
			sequenceEnergyMap[bestSeq]["acceptCycleNumber"] = cycleCounter; // gets the accept cycle number for the current sequence
			sequenceEnergyMap[bestSeq]["xShift"] = _opt.xShift; // gets the accept cycle number for the current sequence
			sequenceEnergyMap[bestSeq]["crossingAngle"] = _opt.crossingAngle; // gets the accept cycle number for the current sequence
			sequenceEnergyMap[bestSeq]["axialRotation"] = _opt.axialRotation; // gets the accept cycle number for the current sequence
			sequenceEnergyMap[bestSeq]["zShift"] = _opt.zShift; // gets the accept cycle number for the current sequence
			_sequenceVectorMap[bestSeq] = currStateVec; // saves the current state vector to the sequence vector map
			allSequenceEnergyMap[bestSeq] = sequenceEnergyMap[bestSeq]; // saves the current sequence energy map to the all sequence energy map
			
			if (_opt.verbose){
				double prevEnergy = bestEnergyTotal; 
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << "; PrevE=  " << prevEnergy << " : CurrE= " << currEnergyTotal;
				cout << "; CurrTemp: " << MC.getCurrentT() << endl;
				cout << "Best sequence: " << bestSeq << endl;
				cout << "Best sequence Info:" << endl;
				cout << "Baseline          " << sequenceEnergyMap[bestSeq]["Baseline"] << endl;
				cout << "Entropy           " << sequenceEnergyMap[bestSeq]["entropyDiff"] << endl;
				cout << "CurrEntropy       " << sequenceEnergyMap[bestSeq]["currEntropy"] << endl;
				cout << "PrevEntropy       " << sequenceEnergyMap[bestSeq]["prevEntropy"] << endl << endl;
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
	double entropy = 0;
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/allDesigns_" + to_string(_rep) + ".pdb");
	// TODO: should I switch this to the energy total?
	for (auto &seq: allSequenceEnergyMap){
		_out << "Best Sequence #" << i << ": " << seq.first << "; Energy: " << seq.second["Dimer"] << endl;
		if (i == 0){
			_sys.setActiveRotamers(_sequenceVectorMap[seq.first]);
			writer.write(_sys.getAtomPointers(), true, false, true);
			_bestSequence = seq.first;
			ener = seq.second["currEnergyTotal"];
			entropy = seq.second["entropyDiff"];
		} else if (seq.second["entropyDiff"] > entropy && seq.second["currEnergyTotal"] < ener){
			_sys.setActiveRotamers(_sequenceVectorMap[seq.first]);
			writer.write(_sys.getAtomPointers(), true, false, true);
			_bestSequence = seq.first;
			ener = seq.second["currEnergyTotal"];
			entropy = seq.second["entropyDiff"];
		}
		_sequenceEnergyMap[seq.first] = seq.second;
		i++;
	}
	//TODO: save the top x sequences during the last cycle after all backbone repacks finished; save the rest in a energy landscape file
	writer.close();
	cout << "End monte carlo sequence search #" << _rep << ": " << diffTimeSMC << "s" << endl;
	_out << "End monte carlo sequence search #" << _rep << ": " << diffTimeSMC << "s" << endl;
	_out << "Monte Carlo ended at Temp: " << MC.getCurrentT() << endl << endl;
	//TODO: maybe set up something like the following:
	// - save sequences for x cycles
	// - recalculate the energy for all sequences at this new geometry
	// - continue on to the next cycle with the current sequence
	// repeat for x times 10 cycles
}

void getDimerSasa(System &_sys, map<string, vector<uint>> &_sequenceVectorMap, map<string, map<string,double>> &_sequenceEnergyMap){
	for (auto &seq : _sequenceEnergyMap){
		string sequence = seq.first;
		vector<uint> state = _sequenceVectorMap[sequence];

		_sys.setActiveRotamers(state);

		//Setup SasaCalculator to calculate the monomer SASA
		SasaCalculator sasa(_sys.getAtomPointers());
		sasa.calcSasa();
		double dimerSasa = sasa.getTotalSasa();

		_sequenceEnergyMap[sequence]["DimerSasa"] = dimerSasa;
	}
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
	Eset->setAllTermsInactive();
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
	// save into vector of backbone positions and residue burial pairs
	vector<pair <int, double> > resiBurial = calculateResidueBurial(sys, _opt, backboneSeq);
	// sort in descending order of burial
	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});

	//cout << "Determining interfacial residues by residue burial..." << endl;
	vector<int> interfacePositions;

	// Output variable Set up
	backboneSeq = generateBackboneSequence("L", _opt.backboneLength, _opt.useAlaAtCTerminus);
	string variablePositionString = generateString("0", backboneSeq.length());
	string rotamerLevels = generateString("0", backboneSeq.length());
	defineRotamerLevelsByResidueBurial(sys, _opt, resiBurial, interfacePositions, rotamerLevels, variablePositionString);

	// checks if interface is defined; if so, check the position and set those positions to the highest rotamer level
	if (_opt.interface != ""){
		interfacePositions.clear(); // reset the interface from the defineRotamerLevelsByResidueBurial function
		useInputInterface(_opt, variablePositionString, rotamerLevels, interfacePositions);
	}	

	// makes the polymer sequence to return basedon the interface positions
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.Ids, interfacePositions);
	PolymerSequence PS(polySeq);
	cout << PS << endl;

	vector<int> rotamerSamplingPerPosition = getRotamerSampling(rotamerLevels); // converts the rotamer sampling for each position as a vector
	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	int highestRotamerLevel = numberOfRotamerLevels-1;
	vector<int> linkedPositions = getLinkedPositions(rotamerSamplingPerPosition, _opt.interfaceLevel, highestRotamerLevel);

	//String for the positions of the sequences that are considered interface for positions amd high rotamers
	string rotamerSamplingString = getInterfaceString(rotamerSamplingPerPosition, backboneSeq.length());

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
	_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition, backboneSeq.length());
	_allInterfacePositions = getAllInterfacePositions(_opt, rotamerSamplingPerPosition, backboneSeq.length());

	if (_opt.interface != ""){
		_interfacePositions.clear();
		_allInterfacePositions.clear();
		for (uint k=0; k<interfacePositions.size(); k++){//TODO: make this not hardcoded to skip RAS
			int pos = interfacePositions[k];
			_interfacePositions.push_back(pos-_opt.thread);
		}
		for (uint k=0; k<backboneSeq.length(); k++){//TODO: make this not hardcoded to skip RAS
			if (k < 3 || k > backboneSeq.length()-5){
				if (rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
					_allInterfacePositions.push_back(k);
				}
			} else {
		 		if (_opt.interface[k] == '1'){
					_allInterfacePositions.push_back(k);
				}
			}
		}
	}

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

void defineRotamerLevelsByResidueBurial(System &_sys, Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString){
	// initialize variables
	int levelCounter = 0;
	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	int highestRotamerLevel = numberOfRotamerLevels-1;
	// loop through the residue burial values calculated above for each position
	cout << "Interface: " << _opt.interface << endl;
	for (uint i = 0; i < _resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(_resiBurial.size()); // calculate the SASA percentile for this position
		// if percentile is greater, move on to the next rotamer level
		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
			levelCounter++;
		}
		int backbonePosition = _resiBurial[i].first; // get the backbone position for this position
		Position &position = _sys.getPosition(backbonePosition); // get the position object for this position
		string positionRotLevel = _opt.sasaRepackLevel[levelCounter]; // get the rotamer level for this position
		int resiNum = position.getResidueNumber(); // get the residue number for this position
		int positionNumber = resiNum-_opt.thread; // position number taking thread into account

		// After a couple of tests, the below should be working now. But if I run into interface problems in the future, come here first
		//TODO: would be nice to make this code ...better
		// check if the current interface level is below the accepted option for interface levels (SASA repack level)
		if (levelCounter < _opt.interfaceLevel) {
			_interfacePositions.push_back(resiNum);
			// check to see if the position is found within the core of protein (i.e. not the first 3 residues or the last 4 residues)
			if (backbonePosition > 5 && backbonePosition < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 3 and 4 here instead of 4 and 5 to prevent changes at the interface like others
				// replace 0 with 1 for variable positions that are found at the interface
				_variablePositionString.replace(_variablePositionString.begin()+positionNumber, _variablePositionString.begin()+positionNumber+1, "1");//TODO: I just added this if statement in. It may or may not work properly because of the numbers (I think it starts at 0 rather than 1 unlike many of the other parts where I hardcode these for baselines
			}
		}
		_rotamerLevels.replace(_rotamerLevels.begin()+positionNumber, _rotamerLevels.begin()+positionNumber+1, MslTools::intToString(levelCounter));
		
	}
}

// sets the gly69 backbone to starting geometry
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, CartesianPoint &_xAxis,
 CartesianPoint &_zAxis, Transforms &_trans) {
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	_sys.readPdb(_opt.backboneFile);

	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, _axisA, _axisB, _ori, _xAxis, _zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

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
	cout << setw(20) << "deleteTerminalBonds" << defaults.deleteTerminalBonds << endl;
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
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	//decided to try loading low number of rotamers for this instead of high
	loadRotamers(sys, sysRot, _opt.SL);

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
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	//spm.calculateEnergies();
	spm.runGreedyOptimizer(_opt.greedyCycles);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	_out << "Time to calculate energies for backbone repack: " << diffTime << endl;
	// do backbone geometry repacks
	if (_opt.verbose){
		cout << "==============================================" << endl;
		cout << " Performing Local Monte Carlo Backbone Repack " << endl;
		cout << "==============================================" << endl;
	}
	monteCarloRepack(_opt, sys, _savedXShift, spm, _helicalAxis, _axisA, _axisB, apvChainA, apvChainB, _trans, _RNG, prevBestEnergy, 
	_rep, _out);

	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/backboneOptimized_" + to_string(_rep) + ".pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?
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
	MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	//MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects, _opt.backboneConvergedSteps, _opt.backboneConvergedE);
	MCMngr.setEner(_prevBestEnergy);

	vector<uint> startStateVec = _spm.getMinStates()[0];
	vector<unsigned int> MCOBest = startStateVec;
		
	unsigned int counter = 0;
	double currentEnergy = _prevBestEnergy;
	double startDimer = _prevBestEnergy;

	// setup variables for shifts: ensures that they start from the proper values for every repack and not just the final value from the initial repack
	bool decreaseMoveSize = _opt.decreaseMoveSize;
	double deltaX = _opt.deltaX;
	double deltaCross = _opt.deltaCross;
	double deltaAx = _opt.deltaAx;
	double deltaZ = _opt.deltaZ;

	PDBWriter writer;
	// loop through the MC cycles for backbone repacks
	while(!MCMngr.getComplete()) {
		// get the current temperature of the MC
		double startTemp = MCMngr.getCurrentT();
		
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
			deltaZShift = getStandardNormal(_RNG) * deltaZ;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(_RNG) * deltaAx;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * deltaCross;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * deltaX;
			backboneMovement(_apvChainA, _apvChainB, _axisA, _axisB, _trans, deltaZShift, moveToPreform);
		}
		
		// Run _optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0];
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE
		
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
			
			// if accept, decrease the value of the moves by the sigmoid function
			if (_opt.decreaseMoveSize == true){
				double endTemp = MCMngr.getCurrentT();
				getCurrentMoveSizes(_opt, startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, decreaseMoveSize);
			}
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
	
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	//TODO: there may be a better way to resolve this, but as of 2022-9-8, I want to get the most data I can before a lab meeting, so putting this here
	if (finalEnergy > 100){
		_out << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		exit(0);
	}
}

void getCurrentMoveSizes(Options &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize) {
	double decreaseMultiplier = _endTemp/_currTemp;
	bool decreaseX = true;
	bool decreaseCross = true;
	bool decreaseAx = true;
	bool decreaseZ = true;
	_deltaX = decreaseMoveSize(_deltaX, _opt.deltaXLimit, decreaseMultiplier, decreaseX);
	_deltaCross = decreaseMoveSize(_deltaCross, _opt.deltaCrossLimit, decreaseMultiplier, decreaseCross);
	_deltaAx = decreaseMoveSize(_deltaAx, _opt.deltaAxLimit, decreaseMultiplier, decreaseAx);
	_deltaZ = decreaseMoveSize(_deltaZ, _opt.deltaZLimit, decreaseMultiplier, decreaseZ);	
	if (decreaseX == false && decreaseCross == false && decreaseAx == false && decreaseZ == false){
		_decreaseMoveSize = false;
	}
}