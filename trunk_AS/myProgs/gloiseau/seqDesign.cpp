#include <iostream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>
#include <chrono>

// MSL Functions
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "CRDReader.h"
#include "SasaCalculator.h"

// My functions
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
void stateMCUnlinked(System &_startGeom, Options &_opt, PolymerSequence &_PS,
 map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
 vector<unsigned int> &_bestState, string &_bestSequence, vector<string> &_seqs, vector<string> &_allSeqs,
 map<string,vector<uint>> &_sequenceVectorMap, vector<uint> &_allInterfacialPositionsList,
 vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, int _rep, ofstream &_out, ofstream &_err);
void searchForBestSequencesUsingThreads(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG, vector<string> &_allSeqs,
 vector<uint> &_bestState, string &_bestSequence, map<string, map<string,double>> &_sequenceEnergyMap, map<string, vector<uint>> &_sequenceVectorMap,
 map<string,double> _sequenceEntropyMap, vector<uint> &_allInterfacialPositionsList, vector<uint> &_interfacialPositionsList,
 vector<int> &_rotamerSampling, int _rep, ofstream &_out, ofstream &_err);
map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap, map<string,double> _sequenceEntropyMap,
 vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList, vector<int> _rotamerSampling);
string getBestSequenceInMap(map<string,map<string,double>> &_sequenceEnergyMap);
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, double _prevEnergy, string _currSeq, vector<uint> _currVec,
 vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
 map<string,double> &_sequenceEntropyMap);
void switchSequence(System &_sys, Options &_opt);

// geometry setup functions
void prepareSystem(Options &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);
void setGly69ToStartingGeometry(Options &_opt, System &_sys, System &_helicalAxis,
 AtomPointerVector &_axisA, AtomPointerVector &_axisB, CartesianPoint &_ori, 
 CartesianPoint &_xAxis, CartesianPoint &_zAxis, Transforms &_trans);
void getStartingGeometry(Options &_opt, vector<double> &_densities, RandomNumberGenerator &_RNG, ofstream &_out);
void checkForClashing(System &_startGeom, Options &_opt, vector<uint> _interfacePositions);

// define interface and rotamer level functions
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, string &_variablePositionString,
 string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
void defineRotamerLevels(Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString);
void useInputInterface(Options &_opt, string &_variablePositionString, string &_rotamerLevels, vector<int> &_interfacePositions);
void setActiveSequence(System &_sys, string _sequence);

// backbone repack functions
void localBackboneRepack(Options &_opt, System &_startGeom, string _sequence, map<string, map<string,double>> &_sequenceEnergyMap, uint _rep, double _savedXShift,
 System &_helicalAxis, AtomPointerVector &_axisA,  AtomPointerVector &_axisB, vector<int> _rotamerSampling, Transforms &_trans,
 RandomNumberGenerator &_RNG, ofstream &_out);
double monteCarloRepack(Options &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 uint _rep, ofstream &_out);
void getCurrentMoveSizes(Options &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);

// output functions
void outputFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, map<string,map<string,double>> _sequenceEnergyMap,
 vector<double> _densities, ofstream &_sout);
string getRunParameters(Options &_opt, vector<double> _densities);
void outputTime(auto _start, string _descriptor, bool _beginFunction, ofstream &_out);

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
	RandomNumberGenerator RNG;
	if (opt.useTimeBasedSeed){
		RNG.setTimeBasedSeed();
	} else {
		RNG.setSeed(opt.seed);
	}
	
	// get the starting geometries
	vector<double> densities;
	getStartingGeometry(opt, densities, RNG, sout); // loads the starting geometry from the geometry file

	// String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.Ids);
	cout << "Amino acids for design: " << alternateIds << endl;

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
	
	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	string variablePositionString; // string of variable positions for each position
	string rotamerSamplingString; // string of rotamer sampling number (if 4 rotamer levels, 0-3 for each position)
	vector<int> linkedPositions; // vector of the positions that will be linked
	vector<uint> interfacePositions; // vector of positions at the interface excluding termini positions
	vector<uint> allInterfacePositions; // vector of positions at the interface including the terminal positions
	vector<int> rotamerSamplingPerPosition; // vector of rotamer level for each position
	
	// Defines the interfacial positions and the number of rotamers to give each position
	PolymerSequence interfacePolySeq = getInterfacialPolymerSequence(opt, startGeom, variablePositionString, rotamerSamplingString,
	 linkedPositions, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, sout);

	outputTime(start, "Identify Interface", false, sout);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	prepareSystem(opt, sys, startGeom, interfacePolySeq);
	checkIfAtomsAreBuilt(sys, err); // check to verify that all atoms have coordinates

	Chain &chainA = sys.getChain("A");
	int seqLength = chainA.positionSize();

	checkForClashing(sys, opt, interfacePositions); // checks a given sequence for clashes; used to identify geometries for design

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	// add in an estimate of monomer energy for each sequence (calculated using only CHARMM_VDW, SCWRL4_HBOND, CHARMM_IMM1, and CHARMM_IMM1REF)
	if (opt.useBaseline){
		buildBaselines(sys, opt);
	}

	// link the interfacial positions (for quicker calculation of initial sequence for homodimers)
	if (opt.linkInterfacialPositions){
		vector<vector<string>> linkedPos = convertToLinkedFormat(sys, linkedPositions, seqLength);
		sys.setLinkedPositions(linkedPos);
	}
	string seq = convertPolymerSeqToOneLetterSeq(chainA);
	cout << "Sequence: " << seq << endl;
	
	// initialize the object for loading rotamers into system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// load rotamers for each amino acid at each position into the system
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	//CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs

	switchSequence(sys, opt);

	/******************************************************************************
	 *           === VARIABLES FOR SAVING ENERGIES AND SEQUENCES ===
	 ******************************************************************************/
	vector<string> seqs; // vector to hold designed sequences
	vector<string> allSeqs; //vector to hold all designed sequences from the stateMC
	map<string, map<string,double>> sequenceEnergyMapBest; // energyMap to hold all energies for output into a summary file
	map<string, vector<uint>> sequenceVectorMap; //sequence and vector map: sequence will be tied vector with AAs and rotamers
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile); // Get sequence entropy map
	
	/******************************************************************************
	 *                    === GET THE BEST STARTING STATE ===
	 ******************************************************************************/
	vector<uint> bestState;
	if (opt.sequence == ""){
		outputTime(start, "Self Consistent Mean Field", true, sout);
		bestState = runSCMFToGetStartingSequence(sys, opt, RNG, rotamerSamplingString, variablePositionString,
		 seqs, allInterfacePositions, sequenceEnergyMapBest, sequenceVectorMap, sequenceEntropyMap, sout);
		outputTime(start, "Self Consistent Mean Field", false, sout);
		// set system to the best sequence or input sequence 
		sys.setActiveRotamers(bestState);
		seq = convertPolymerSeqToOneLetterSeq(chainA);
		double bestEnergy = sys.calcEnergy();
		cout << "Sequence: " << seq << "; Energy: " << bestEnergy << endl;
		// reset the energy set
		sys.getEnergySet()->eraseTerm("CHARMM_VDW");
		sys.getEnergySet()->eraseTerm("CHARMM_IMM1");
		sys.getEnergySet()->eraseTerm("CHARMM_IMM1REF");
		sys.getEnergySet()->eraseTerm("SCWRL4_HBOND");
		// Unlink the best state from SCMF if not using linked positions during the state Monte Carlo
		if(opt.linkInterfacialPositions){
			unlinkBestState(opt, bestState, rotamerSamplingPerPosition, seqLength);
		}
	}

	seq = convertPolymerSeqToOneLetterSeq(chainA);
	cout << "Sequence before stateMC: " << seq << endl;
	/******************************************************************************
	 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
	 ******************************************************************************/
	map<string, map<string,double>> sequenceEnergyMapFinalSeq; // energyMap to hold all energies for output into a summary file
	string bestSequence;
	string prevSequence;
	map<string, double> geometries;
	PDBWriter helicalAxisWriter;
	for (uint i=0; i<3; i++){
		helicalAxisWriter.open(opt.pdbOutputDir+"/helicalAxis"+to_string(i)+".pdb");
		helicalAxisWriter.write(helicalAxis.getAtomPointers(), true, false, true);
		helicalAxisWriter.close();
		// run state Monte Carlo to get a random sequence from the best state
		outputTime(start, "Sequence search replicate " + to_string(i), true, sout);
		stateMCUnlinked(sys, opt, interfacePolySeq, sequenceEnergyMapBest, sequenceEntropyMap, bestState, bestSequence, seqs, allSeqs,
			sequenceVectorMap, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, i, sout, err);
		outputTime(start, "Sequence search replicate " + to_string(i), false, sout);
		
		// TODO: search through and see if the best sequence is found in the final energy map (won't have to do if improvements for more than one cycle are miniscule)
		if (prevSequence != bestSequence){
			// compute monomer energy
			outputTime(start, "Compute Monomer Energy " + to_string(i), true, sout);
			computeMonomerEnergyIMM1(opt, trans, sequenceEnergyMapBest, bestSequence, RNG, sout, err);
			outputTime(start, "Compute Monomer Energy " + to_string(i), false, sout);
			// get the current energy of the sequence by subtracting monomer from dimer
			double monomerEnergy = sequenceEnergyMapBest[bestSequence]["Monomer"];
			double dimerEnergy = sequenceEnergyMapBest[bestSequence]["Dimer"];
			cout << "Dimer Energy: " << bestSequence << ": " << dimerEnergy << endl;
			// set as the start energy
			sequenceEnergyMapBest[bestSequence]["startEnergy"] = dimerEnergy-monomerEnergy;
			sequenceEnergyMapBest[bestSequence]["geometryNumber"] = i;
			// add the monomer energy to the backbone repack to compare instead of using baseline energy (more accurate)
			outputTime(start, "Backbone repack replicate " + to_string(i), true, sout);
			localBackboneRepack(opt, sys, bestSequence, sequenceEnergyMapBest, i, opt.xShift, helicalAxis, axisA, axisB,
			 rotamerSamplingPerPosition, trans, RNG, sout);
			outputTime(start, "Backbone repack replicate " + to_string(i), false, sout);
			// set the end geometries after the repack for the sequence
			sequenceEnergyMapBest[bestSequence]["endXShift"] = opt.xShift;
			sequenceEnergyMapBest[bestSequence]["endCrossingAngle"] = opt.crossingAngle;
			sequenceEnergyMapBest[bestSequence]["endAxialRotation"] = opt.axialRotation;
			sequenceEnergyMapBest[bestSequence]["endZShift"] = opt.zShift;
			// get the best sequence from the energy map
			sequenceEnergyMapFinalSeq[bestSequence] = sequenceEnergyMapBest[bestSequence];	
			// reset the energy map
			sequenceEnergyMapBest.clear(); // energyMap to hold all energies for output into a summary file
		} else {
			// end loop if same sequence from seq monte carlo
			i = 3;	
		}
		prevSequence = bestSequence;
	}

	/******************************************************************************
	 *                   === WRITE OUT ENERGY AND DESIGN FILES ===
	 ******************************************************************************/
	outputFiles(opt, rotamerSamplingString, rotamerSamplingPerPosition, sequenceEnergyMapFinalSeq, densities, sout);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

	err.close();
	sout.close();
}

//Functions
void getStartingGeometry(Options &_opt, vector<double> &_densities, RandomNumberGenerator &_RNG, ofstream &_out){
	if (_opt.getGeoFromPDBData){
		// randomly choose a geometry based on densities from membrane protein pdb data
		getGeometry(_opt, _RNG, _densities, _out);
		cout << "xShift:        " << _opt.xShift << "\tDensity: " << _densities[0] << endl;
		cout << "crossingAngle: " << _opt.crossingAngle << "\tDensity: " << _densities[0] << endl;
		cout << "axialRotation: " << _opt.axialRotation << "\tDensity: " << _densities[1] << endl;
		cout << "zShift:        " << _opt.zShift << "\tDensity: " << _densities[2] << endl << endl;
	} else if (_opt.getRandomAxRotAndZShift){
			// use a given xShift and crossing angle, and randomly choose axial rotation and zShift
			getAxialRotAndZShift(_opt, _RNG, _densities, _out);
			cout << "***STARTING GEOMETRY:***" << endl;
			cout << "xShift:        " << _opt.xShift << endl;
			cout << "crossingAngle: " << _opt.crossingAngle << endl;
			cout << "axialRotation: " << _opt.axialRotation << "\tDensity: " << _densities[1] << endl;
			cout << "zShift:        " << _opt.zShift << "\tDensity: " << _densities[2] << endl << endl;
	} else {
		// use the given xShift, crossing angle, axial rotation, and zShift
		cout << "***STARTING GEOMETRY:***" << endl;
		cout << "xShift:        " << _opt.xShift << endl;
		cout << "crossingAngle: " << _opt.crossingAngle << endl;
		cout << "axialRotation: " << _opt.axialRotation << endl;
		cout << "zShift:        " << _opt.zShift << endl << endl;
		_densities.push_back(0);
		_densities.push_back(0);
		_densities.push_back(0);
	}
}
// switches to the starting sequence (if given, otherwise set to polyleu)
void switchSequence(System &_sys, Options &_opt){
	if (_opt.sequence != ""){
		setActiveSequence(_sys, _opt.sequence);
	} else {
		string polyLeu = "LLLLLLLLLLLLLLLLLLILI";
		setActiveSequence(_sys, polyLeu);
	}
}

// The below checks for clashing at the interface by looking at the energy for a polyAla interface sequence
void checkForClashing(System &_startGeom, Options &_opt, vector<uint> _interfacePositions){
	// declare system
	// change the ends of sequence to ala
	string polyLeu;
	for (uint i = 0; i < _opt.backboneLength; i++){
		if (i < 4 || i > _opt.backboneLength-5){
			polyLeu += "A";
		} else {
			polyLeu += "L";
		}
	}
	if (_opt.xShift <= 7.5){
		for (uint i=0; i< _interfacePositions.size(); i++){
			if (i%2==0){
				polyLeu[_interfacePositions[i]] = 'G';
			} else {
				polyLeu[_interfacePositions[i]] = 'A';
			}
		}
	} else {
		for (uint i=0; i< _interfacePositions.size(); i++){
			polyLeu[_interfacePositions[i]] = 'A';
		}
	}

cout << "polyLeu: " << polyLeu << endl;
	string backboneSeq = convertToPolymerSequence(polyLeu, _opt.thread);
	PolymerSequence PS(backboneSeq);
cout << "backboneSeq: " << PS << endl;
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << PS << endl;
		exit(0);
	}
	
	// assign the coordinates of our system to the given geometry 
	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	sys.buildAllAtoms();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	CSB.updateNonBonded(10,12,50);

	AtomSelection sel(sys.getAtomPointers());
	sel.select("chainA, chain A");
	sel.select("chainB, chain B");
	double energy = sys.calcEnergy("chainA", "chainB");
	if (energy > 10){
		cout << "Clashing at the interface; energy = " << energy << endl;
		exit(0);
	} 
	
	// if no clashing, then output the pdb
	// convert axial rotation to positive using absolute value, for outputting
	double absAxRot = abs(_opt.axialRotation);
	cout << absAxRot << endl;
	double adjustedAx = (10*absAxRot/9)-(200*_opt.zShift/27);
	double adjustedZ = (10*_opt.zShift/9)-(0.15*absAxRot/9);
	string geometry = "/x"+MslTools::doubleToString(_opt.xShift)+"_cross"+MslTools::doubleToString(_opt.crossingAngle)
	 +"_ax"+MslTools::doubleToString(adjustedAx)+"_z"+MslTools::doubleToString(adjustedZ)+"_vdW"+to_string(energy)+".pdb";
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + geometry);
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();
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

// sequence search functions
void stateMCUnlinked(System &_startGeom, Options &_opt, PolymerSequence &_PS,
map<string, map<string,double>> &_sequenceEnergyMap, map<string,double> &_sequenceEntropyMap,
vector<unsigned int> &_bestState, string &_bestSequence, vector<string> &_seqs, vector<string> &_allSeqs,
map<string,vector<uint>> &_sequenceVectorMap, vector<uint> &_allInterfacialPositionsList,
vector<uint> &_interfacialPositionsList, vector<int> &_rotamerSampling, RandomNumberGenerator &_RNG, int _rep, ofstream &_out, ofstream &_err){
	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH NEW POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	prepareSystem(_opt, sys, _startGeom, _PS);

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

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	//spm.updateWeights();
	spm.setOnTheFly(false);
	spm.saveEnergiesByTerm(true);
	
	outputTime(start, "Self and Pair of all amino acids calculation " + to_string(_rep), true, _out);
	spm.calculateEnergies();
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
	vector<vector<bool>> mask = getActiveMask(_sys);
	_spm.runGreedyOptimizer(_opt.greedyCycles, mask);
	_bestState = _spm.getMinStates()[0];
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
	int acceptCounter = 0;
	Chain & chain = _sys.getChain("A");
	string prevStateSeq = convertPolymerSeqToOneLetterSeq(chain);
	
	cout << "Sequence at start of stateMC: " << prevStateSeq << endl;

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
	// initialize energy variables for the MonteCarlo
	string bestSeq = prevStateSeq;
	map<string,double> bestSequenceEnergyMap;
	// Monte Carlo while loop for finding the best sequences
	while (!MC.getComplete()){
		if (_opt.verbose){
			cout << "Cycle #" << cycleCounter << "; Acceptances: " << acceptCounter << endl;
			//cout << "Starting Seq: " << prevStateSeq << endl;
		}
		// get the sequence entropy probability for the current best sequence
		map<string,vector<uint>> sequenceVectorMap;
		map<string,map<string,double>> sequenceEnergyMap = mutateRandomPosition(_sys, _opt, _spm, _RNG, bestSeq, bestEnergy, 
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
		
		cout << "Best Sequence: " << currSeq << "; Energy: " << currEnergyTotal << "; Entropy: " << currStateEntropy << endl;
		cout << "Prev Sequence: " << prevStateSeq << "; Energy: " << bestEnergyTotal << "; Entropy: " << prevStateEntropy << endl;

		// MC accept and reject conditions
		if (!MC.accept(currEnergyTotal)){
			_sys.setActiveRotamers(prevStateVec); // set rotamers to the previous state
			if (_opt.verbose){
				cout << "State not accepted, E= " << currEnergyTotal << "; PrevE= " << bestEnergyTotal << endl;
			}
		} else {
			_sys.setActiveRotamers(currStateVec); // set rotamers to the current state
			map<string,double> energyMap = sequenceEnergyMap[currSeq];
			lout << prevStateSeq << "\t" << currSeq << "\t" << bestEnergyTotal << "\t" << currEnergyTotal << "\t";
			lout << prevStateEntropy << "\t" << currStateEntropy << "\t" << MC.getCurrentT() << endl;
			prevStateVec = currStateVec; // set the previous state vector to be the current state vector
			prevStateSeq = currSeq; // set the previous sequence to be the current sequence
			bestSeq = currSeq; // set the best sequence to the newly accepted current sequence
			bestEnergy = sequenceEnergyMap[bestSeq]["Dimerw/Baseline"]; // set the best energy to the current energy (vdw, hbond, imm1, baseline)
			sequenceEnergyMap[bestSeq]["acceptCycleNumber"] = acceptCounter; // gets the accept cycle number for the current sequence
			// saves the geometry into the map; maybe make this a function or add it somewhere more appropriate
			sequenceEnergyMap[bestSeq]["xShift"] = _opt.xShift; 
			sequenceEnergyMap[bestSeq]["crossingAngle"] = _opt.crossingAngle;
			sequenceEnergyMap[bestSeq]["axialRotation"] = _opt.axialRotation;
			sequenceEnergyMap[bestSeq]["zShift"] = _opt.zShift; 
			_sequenceVectorMap[bestSeq] = currStateVec; // saves the current state vector to the sequence vector map
			bestSequenceEnergyMap = sequenceEnergyMap[bestSeq]; // saves the current sequence energy map to the all sequence energy map
			
			if (_opt.verbose){
				double prevEnergy = bestEnergyTotal; 
				cout << "Cycle#" << cycleCounter << " State accepted, Sequence: " << currSeq << "; PrevE=  " << prevEnergy << " : CurrE= " << currEnergyTotal;
				cout << "; CurrTemp: " << MC.getCurrentT() << endl;
				cout << "Best sequence: " << bestSeq << endl;
				cout << "Best sequence Info:" << endl;
				cout << "CurrVdw           " << sequenceEnergyMap[bestSeq]["VDWDimer"] << endl;
				cout << "CurrIMM1          " << sequenceEnergyMap[bestSeq]["IMM1Dimer"] << endl;
				cout << "CurrHBOND         " << sequenceEnergyMap[bestSeq]["HBONDDimer"] << endl;
				cout << "CurrEntropy       " << sequenceEnergyMap[bestSeq]["currEntropy"] << endl;
				cout << "PrevEntropy       " << sequenceEnergyMap[bestSeq]["prevEntropy"] << endl << endl;
			}
			acceptCounter++;
		}
		//Reset the MC to run 100 more cycles to
		//if (MC.getComplete() == true && MC.getCurrentT() < 546.4){
		//	MC.reset(3649, 3649, 500, MonteCarloManager::EXPONENTIAL, 10);//Approximately 50% likely to accept within 5kcal, and 25% likely to accept within 10kcal
		//}
		cycleCounter++;
	}
	time(&endTimeSMC);
	diffTimeSMC = difftime (endTimeSMC, startTimeSMC);
	_bestSequence = bestSeq;

	lout << "Time: " << diffTimeSMC << "s" << endl;
	lout.close();
	_allSeqs.clear();
	//addSequencesToVector(energyVector, _allSeqs);
	getDimerSasa(_sys, _sequenceVectorMap, _sequenceEnergyMap);
	_sequenceEnergyMap[bestSeq] = bestSequenceEnergyMap; // saves the best sequence energy map into the global energy map
	
	// output the best sequence pdb 
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/bestSequence_" + to_string(_rep) + ".pdb");
	_sys.setActiveRotamers(_sequenceVectorMap[_bestSequence]);
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	_sys.saveAltCoor("bestSequence_" + to_string(_rep));

	cout << "End monte carlo sequence search #" << _rep << ": " << diffTimeSMC/60 << "min" << endl;
	_out << "End monte carlo sequence search #" << _rep << ": " << diffTimeSMC/60 << "min" << endl;
	_out << "Monte Carlo ended at Temp: " << MC.getCurrentT() << endl << endl;
}

map<string,map<string,double>> mutateRandomPosition(System &_sys, Options &_opt, SelfPairManager &_spm, RandomNumberGenerator &_RNG,
 string _bestSeq, double _bestEnergy, map<string,vector<uint>> &_sequenceStateMap, map<string,double> _sequenceEntropyMap,
 vector<uint> _allInterfacialPositionsList, vector<uint> _interfacialPositionsList, vector<int> _rotamerSampling){
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
			threads.push_back(thread{energyFunction, ref(_opt), ref(_spm), _bestSeq, _bestEnergy, currSeq, currVec, ref(_rotamerSampling), ref(_allInterfacialPositionsList), 
			 ref(sequenceEnergyMap), ref(_sequenceEntropyMap)});
		}
	} 
	// join all the threads (wait for them all to finish before continuing)
	for (auto& th : threads){
		th.join();
	}
	return sequenceEnergyMap;
}

void setActiveSequence(System &_sys, string _sequence){
	// Set the active sequence for the system
	for (uint i=0; i<_sys.getPositions().size()/2; i++){
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+_sequence.length());
		string posIdA = posA.getPositionId();
		string posIdB = posB.getPositionId();
		string aa = MslTools::getThreeLetterCode(_sequence.substr(i, 1));
		_sys.setActiveIdentity(posIdA, aa);
		_sys.setActiveIdentity(posIdB, aa);
	}
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
void energyFunction(Options &_opt, SelfPairManager &_spm, string _prevSeq, double _prevEnergy, string _currSeq, vector<uint> _currVec,
 vector<int> &_rotamerSampling, vector<uint> &_allInterfacePositions, map<string,map<string,double>> &_seqEnergyMap,
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

// This takes poly-val helix to calculate the residue burial of every position and based on the burial and number
// of 'SASA interface level' decides rotamer level to assign to the position and also decides which of these positions are 'interfacial'
// PS is the actual polymerSeq object whereas polySeq is the string version of the polymerSeq
PolymerSequence getInterfacialPolymerSequence(Options &_opt, System &_startGeom, string &_variablePositionString,
 string &_rotamerSamplingString, vector<int> &_linkedPositions, vector<uint> &_allInterfacePositions,
 vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out){
	
	// generate a backboneSequence to determine the interface positions using residue burial (defaults to polyVal)
	string backboneSeq = generateString(_opt.backboneAA, _opt.backboneLength);

	// save into vector of backbone positions and residue burial pairs
	vector<pair <int, double> > resiBurial = calculateResidueBurial(_opt, _startGeom, backboneSeq);
	vector<int> interfacePositions;
	
	// if sequence is not empty, use polyLeu
	if (_opt.sequence == ""){
		backboneSeq = generateBackboneSequence("L", _opt.backboneLength, _opt.useAlaAtCTerminus);
	} else {
		backboneSeq = _opt.sequence;
	}

	// setup 0 string to represent variable positions and rotamer levels
	string variablePositionString = generateString("0", backboneSeq.length());
	string rotamerLevels = generateString("0", backboneSeq.length());

	// define rotamer levels for each position based on residue burial
	defineRotamerLevels(_opt, resiBurial, interfacePositions, rotamerLevels, variablePositionString);

	// checks if interface is defined; if so, check the position and set those positions to the highest rotamer level
	if (_opt.interface != ""){
		interfacePositions.clear(); // reset the interface from the defineRotamerLevels function
		useInputInterface(_opt, variablePositionString, rotamerLevels, interfacePositions);
	}	

	//	
	vector<int> rotamerSamplingPerPosition = getRotamerSampling(rotamerLevels); // converts the rotamer sampling for each position as a vector
	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	int highestRotamerLevel = numberOfRotamerLevels-1;
	vector<int> linkedPositions = getLinkedPositions(rotamerSamplingPerPosition, _opt.interfaceLevel, highestRotamerLevel);

	//String for the positions of the sequences that are considered interface for positions amd high rotamers
	string rotamerSamplingString = getInterfaceString(rotamerSamplingPerPosition, backboneSeq.length());

	// Define referenced output variables
	_rotamerSamplingPerPosition = rotamerSamplingPerPosition;
	_variablePositionString = variablePositionString;
	_rotamerSamplingString = rotamerSamplingString;
	_linkedPositions = linkedPositions;

	// define the interface positions for the sequence search	
	if (_opt.interface != ""){
		_interfacePositions.clear();
		_allInterfacePositions.clear();
		for (uint k=0; k<interfacePositions.size(); k++){
			int pos = interfacePositions[k];
			_interfacePositions.push_back(pos-_opt.thread);
		}
		for (uint k=0; k<backboneSeq.length(); k++){
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
	} else {
		_interfacePositions = getInterfacePositions(_opt, rotamerSamplingPerPosition, backboneSeq.length());
		_allInterfacePositions = getAllInterfacePositions(_opt, rotamerSamplingPerPosition, backboneSeq.length());
	}

	_out << endl;
	_out << "Backbone:           " << backboneSeq << endl;
	_out << "Variable Positions: " << variablePositionString << endl;
	_out << "Rotamers Levels:    " << rotamerSamplingString << endl;
	cout << "Backbone:           " << backboneSeq << endl;
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

	// makes the polymer sequence to return basedon the interface positions
	string polySeq = generateMultiIDPolymerSequence(backboneSeq, _opt.thread, _opt.Ids, interfacePositions);
	PolymerSequence PS(polySeq);
	cout << PS << endl;

	return PS;
}

void defineRotamerLevels(Options &_opt, vector<pair <int, double> > &_resiBurial, vector<int> &_interfacePositions,
 string &_rotamerLevels, string &_variablePositionString){
	// initialize variables
	int levelCounter = 0;
	int numberOfRotamerLevels = _opt.sasaRepackLevel.size();
	//int highestRotamerLevel = numberOfRotamerLevels-1;
	// loop through the residue burial values calculated above for each position
	cout << "Interface: " << _opt.interface << endl;
	for (uint i = 0; i < _resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(_resiBurial.size()); // calculate the SASA percentile for this position
		// if percentile is greater, move on to the next rotamer level
		if (sasaPercentile > (levelCounter+1)/double(numberOfRotamerLevels)) {
			levelCounter++;
		}
		int positionNumber = _resiBurial[i].first; // get the backbone position for this position
		int resiNum = positionNumber+_opt.thread;
		// check if the current interface level is below the accepted option for interface levels (SASA repack level)
		if (levelCounter < _opt.interfaceLevel) {
			_interfacePositions.push_back(resiNum);
			// check to see if the position is found within the core of protein (i.e. not the first 3 residues or the last 4 residues)
			if (positionNumber > 2 && positionNumber < _opt.backboneLength-4){//backbone position goes from 0-20, so numbers need to be 2 and 4 here instead of 3 and 5 to prevent changes at the interface like others
				// replace 0 with 1 for variable positions that are found at the interface
				_variablePositionString.replace(_variablePositionString.begin()+positionNumber, _variablePositionString.begin()+positionNumber+1, "1");
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

void localBackboneRepack(Options &_opt, System &_startGeom, string _sequence, map<string, map<string,double>> &_sequenceEnergyMap, uint _rep, double _savedXShift,
 System &_helicalAxis, AtomPointerVector &_axisA,  AtomPointerVector &_axisB, vector<int> _rotamerSampling, Transforms &_trans,
 RandomNumberGenerator &_RNG, ofstream &_out){
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PS(polySeq);

	// set up the system for the input sequence
	System sys;
	prepareSystem(_opt, sys, _startGeom, PS);
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// load the rotamers
	loadRotamers(sys, sysRot, _opt.SL);

	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
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
	// get monomer energy
	double monomerEnergy = _sequenceEnergyMap[_sequence]["Monomer"];

	double finalEnergy = monteCarloRepack(_opt, sys, _savedXShift, spm, _helicalAxis, _axisA, _axisB, apvChainA, apvChainB, _trans, _RNG, monomerEnergy, _rep, _out);

	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/backboneOptimized_" + to_string(_rep) + ".pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

	// assign the coordinates of our system to the given geometry 
	_startGeom.assignCoordinates(sys.getAtomPointers(),false);
	_startGeom.buildAllAtoms();
}

double monteCarloRepack(Options &_opt, System &_sys, double &_savedXShift, SelfPairManager &_spm, 
 System &_helicalAxis, AtomPointerVector &_axisA, AtomPointerVector &_axisB, AtomPointerVector &_apvChainA,
 AtomPointerVector &_apvChainB, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy,
 uint _rep, ofstream &_out){
	// Setup backbone repack file
	ofstream bbout;
	string bboutfile  = _opt.pdbOutputDir + "/bbRepack_" + to_string(_rep) + ".out";
	bbout.open(bboutfile.c_str());

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
	//MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects);

	vector<uint> startStateVec = _spm.getMinStates()[0];
	vector<unsigned int> MCOBest = startStateVec;
		
	unsigned int counter = 0;
	_sys.setActiveRotamers(startStateVec);
	double currentEnergy = _sys.calcEnergy()-_monomerEnergy;
	double bestEnergy = currentEnergy;
	double prevBestEnergy = currentEnergy;
	MCMngr.setEner(currentEnergy);
	//double startDimer = _prevBestEnergy;

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
		
		// Run optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0]-_monomerEnergy;
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE
		
		if (!MCMngr.accept(currentEnergy)) {
			bbout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy << endl;
		} else {
			bestEnergy = currentEnergy;
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
			bbout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy << endl;
			counter++;
			writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	writer.close();
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
		
	_sys.applySavedCoor("savedRepackState");
	double dimerEnergy = _spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-_monomerEnergy;
	
	bbout << "Energy #" << _rep << ": " << finalEnergy << endl;

	// Output change in geometry
	bbout << "***GEOMETRY SHIFT***" << endl;
	bbout << "xShift;        Before: " << _opt.xShift << "; After: " << xShift << endl;
	bbout << "crossingAngle; Before: " << _opt.crossingAngle << "; After: " << crossingAngle << endl;
	bbout << "axialRotation; Before: " << _opt.axialRotation << "; After: " << axialRotation << endl;
	bbout << "zShift;        Before: " << _opt.zShift << "; After: " << zShift << endl;
	bbout << "Energy;        Before: " << prevBestEnergy << "; After: " << bestEnergy << endl << endl;
	_out << "***GEOMETRY SHIFT***" << endl;
	_out << "xShift;        Before: " << _opt.xShift << "; After: " << xShift << endl;
	_out << "crossingAngle; Before: " << _opt.crossingAngle << "; After: " << crossingAngle << endl;
	_out << "axialRotation; Before: " << _opt.axialRotation << "; After: " << axialRotation << endl;
	_out << "zShift;        Before: " << _opt.zShift << "; After: " << zShift << endl << endl;
	_out << "Energy;        Before: " << prevBestEnergy << "; After: " << bestEnergy << endl << endl;

	// sets the updated backbone parameters
	_opt.xShift = xShift;
	_opt.crossingAngle = crossingAngle;
	_opt.axialRotation = axialRotation;
	_opt.zShift = zShift;
	bbout << MCMngr.getReasonCompleted() << endl;	
	_out << MCMngr.getReasonCompleted() << endl;	
	bbout << "Monte Carlo repack complete. Time: " << diffTimeMC/60 << "min" << endl << endl;
	//TODO: there may be a better way to resolve this, but as of 2022-9-8, I want to get the most data I can before a lab meeting, so putting this here
	if (finalEnergy > 100){
        bbout << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		_out << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		exit(0);
	}
	return finalEnergy;
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

// define interface functions
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

string getRunParameters(Options &_opt, vector<double> _densities){
	stringstream ss;
	string t = "\t";
	ss <<  _densities[0] << t << _densities[1] << t << _densities[2] << t << _opt.thread << t << _opt.sasaRepackLevel.size() << t << _opt.interfaceLevel << t << _opt.backboneLength;
	string runParameters = ss.str();
	return runParameters;
}

void outputFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, map<string,map<string,double>> _sequenceEnergyMap,
 vector<double> _densities, ofstream &_sout){
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
		// For adding in strings to a line for the energy file; looping through the terms instead of my input terms this time; sort later
		for (auto &ener: energyMap){
			string energyTerm = ener.first;
			double energy = energyMap[energyTerm];
			string term = MslTools::doubleToString(energy)+t;
			seqLine << term;
			if (i == 0){
				enerTerms << energyTerm << t;
			}
		}
		string line = seqLine.str();
		energyLines.push_back(line);
		i++;
	}
	ofstream eout;
	string eoutfile = _opt.pdbOutputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	eout << "Sequence" << t << "Interface" << t << "InterfaceSequence" << t;
	eout << enerTerms.str() << endl;
	_sout << "Sequence" << t << "Interface" << t << "InterfaceSequence" << t;
	_sout << enerTerms.str() << endl;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
		_sout << energyLines[i] << endl;
	}
	eout.close();
}

void outputTime(auto _start, string _descriptor, bool _beginFunction, ofstream &_out){
	auto end = chrono::system_clock::now();
	time_t endTimeFormatted = chrono::system_clock::to_time_t(end); 
	chrono::duration<double> elapsedTime = end-_start;
	_out.precision(3);// rounds output to 3 decimal places
	cout.precision(3);// rounds output to 3 decimal places
	if (_beginFunction == true){
		_out << _descriptor << " started. Time: " << ctime(&endTimeFormatted);
		_out << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
		cout << _descriptor << " started. Time: " << ctime(&endTimeFormatted);
		cout << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
	} else {
		_out << _descriptor << " finished. Time: " << ctime(&endTimeFormatted);
		_out << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
		cout << _descriptor << " finished. Time: " << ctime(&endTimeFormatted);
		cout << "Elapsed time of program: " << elapsedTime.count()/60 << "min" << endl << endl;
	}
}