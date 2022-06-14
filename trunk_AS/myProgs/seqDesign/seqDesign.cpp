/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/11
 * @Email:  gjowl04@gmail.com
 * @Filename: seqDesign.cpp
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/22
 */
#include <sstream>
#include <iterator>
#include <unistd.h>

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

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************
 *
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
	//Add in some default options that can easily be changed here
	Options opt = parseOptions(argc, argv, defaults);

	if (opt.errorFlag) {
		outputErrorMessage(opt);
		exit(1);
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

	setupDesignDirectory(opt, date);

	string soutfile = opt.pdbOutputDir + "/summary.out";
	string errfile  = opt.pdbOutputDir + "/errors.out";
	string rerunfile = opt.pdbOutputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	rerun << opt.rerunConf << endl;
	rerun.close();

	sout << date << endl;
	err << date << endl;

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


	// TODO:...
	// *Change this function for selecting geometries to design
	vector<double> densities;
	if (opt.getGeoFromPDBData){
		getGeometry(opt, RNG, densities, sout);
	}

	// Output the starting geometry
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << opt.xShift << "\tDensity: " << densities[0] << endl;
	cout << "crossingAngle: " << opt.crossingAngle << "\tDensity: " << densities[0] << endl;
	cout << "axialRotation: " << opt.axialRotation << "\tDensity: " << densities[1] << endl;
	cout << "zShift:        " << opt.zShift << "\tDensity: " << densities[2] << endl << endl;

	//String for the alternateIds at the interface
	string alternateIds = getAlternateIdString(opt.Ids);
	cout << "Amino acids for design: LEU " << alternateIds << endl;

	/******************************************************************************
	 *                       === GENERATE POLYMER SEQUENCE ===
	 ******************************************************************************/
	// polymer sequences have: chain, starting position of chain residue, three letter AA code
	string polySeq = generatePolymerSequence(opt.backboneAA, opt.backboneLength, opt.thread);
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	 // Sets the helices to the origin
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
		err << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

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
	defineInterfaceAndRotamerSampling(opt, PS, rotamerLevels, polySeq, variablePositionString, rotamerSamplingString, linkedPositions, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, sout, axis);

	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	System pdb;
	pdb.readPdb(opt.infile);

	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = pdb.getChain("A");
	Chain & chainB = pdb.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Transform to chosen geometry
	// Change this function to give other geometric parameters for the independent helices for the heterodimers
	// This is directly from CATM
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	/******************************************************************************
	 *  === DECLARE SYSTEM FOR POLYLEU WITH ALTERNATE IDENTITIES AT INTERFACE ===
	 ******************************************************************************/
	// This makes no fucking sense but you make another system to mirror the poly-gly geometry for the sequence that you actually want to run

	// Initialize system for dimer to design
	System sys;
	// Initialize CharmmSystemBuilder to build energy terms for design
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile, opt.solvFile);
	// Explicit definitions of which terms to use and which to not use
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	PolymerSequence PL(polySeq);
	if(!CSB.buildSystem(PL)) {
		err << "Unable to build system from " << polySeq << endl;
		exit(0);
	}

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// assign the coordinates of our system to the given geometry that was assigned without energies using System pdb
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	// initialize the object for loading rotamers into our system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	deleteTerminalHydrogenBondInteractions(sys,opt);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// check to verify that all atoms have coordinates
	// See if this is integratable into the system that is built previously for the helix of interest
	checkIfAtomsAreBuilt(sys, err);

	/******************************************************************************
	 *                     === ADD IN BASELINE ENERGIES ===
	 ******************************************************************************/
	// This should already be ready to go for two independent hetero monomers
	// Would need to rerun this if wanting to expand library of AAs

	if (opt.useBaseline){
		//addBaselineToSelfPairManager();//TODO: was going to make this a function but I feel like it already exists in spm, so going to wait when I have time to look that up
		//initialize baseline maps to calculate energy for each sequence for comparison in baselineMonomerComparison_x.out
		map<string, double> selfMap = readSingleParameters(opt.selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(opt.pairEnergyFile);
		buildSelfInteractions(sys, selfMap);
		buildPairInteractions(sys, pairMap);
	}

	// TODO: change this; I feel like this should be an option if linked is true
	vector<vector<string>> linkedPos = convertToLinkedFormat(sys, linkedPositions, opt.backboneLength);
	sys.setLinkedPositions(linkedPos);

	// Get sequence entropy map
	map<string, double> sequenceEntropyMap = readSingleParameters(opt.sequenceEntropyFile);

	//decided to try loading low number of rotamers for this instead of high
	loadRotamers(sys, sysRot, opt, rotamerSamplingPerPosition);
	CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs
	sys.calcEnergy();

	/******************************************************************************
	 *                        === SETUP SPM AND RUN SCMF ===
	 ******************************************************************************/
	//TODO: make it so that this part is optional: if I submit a sequence to start, don't even do this. Just find the interface, greedy for best sequence, and continue
	// DEE takes forever, apparently this is coded awfully (blame Alessandro)
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.setRunDEE(opt.runDEESingles, opt.runDEEPairs);
	spm.setOnTheFly(false);
	spm.setMCOptions(1000, 0.5, 5000, 3, 10, 1000, 0.01);//changed to sigmoid and added up to 5000
	spm.saveEnergiesByTerm(true); //added back in on 09_21_2021 to get the vdw and hbond energies
	spm.calculateEnergies();

	//Setup running SCMF or UnbiasedMC
	if (opt.runSCMF == true){
		cout << "Running Self Consistent Mean Field" << endl;
		sout << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	} else {
		cout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		sout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	// run and find a sequence using the chosen parameters (MCOptions, SCMF, DEE, etc.)
	time(&spmStart);
	spm.runOptimizer();
	time(&spmEnd);
	spmTime = difftime (spmEnd, spmStart);

	// vector for the initial SCMF state
	// This is fucking stupid
	vector<unsigned int> initialState = spm.getSCMFstate();
	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
	// define startingState vector as the bestState from the SelfPairManager optimizer
	vector<unsigned int> startingState = bestState;

	// energy terms from the startingState
	double bestEnergy = spm.getStateEnergy(bestState);
	double hbondEnergy = spm.getStateEnergy(bestState, "SCWRL4_HBOND");
	double vdwEnergy = spm.getStateEnergy(bestState, "CHARMM_VDW");

	// checks to see if the sequence changed between the initialState and the bestState
	sys.setActiveRotamers(initialState);
	string initialSeq = convertPolymerSeqToOneLetterSeq(sys.getChain("A"));
	sys.setActiveRotamers(bestState);
	string SCMFBestSeq = convertPolymerSeqToOneLetterSeq(sys.getChain("A"));
	sys.calcEnergy();
	string startSequence = SCMFBestSeq;//used for outputting starting sequence

	// TODO: make into a function?
	// output this information about the SelfPairManager run below
	string seqInterface = getInterfaceSequence(opt, rotamerSamplingString, startSequence);
	cout << "Initial Sequence:   " << initialSeq << endl;
	cout << "Best Sequence:      " << startSequence << endl;
	cout << "Interface Sequence: " << seqInterface << endl;
	cout << "Interface:          " << variablePositionString << endl;
	cout << "Total Energy:       " << bestEnergy << endl;
	cout << "VDW:                " << vdwEnergy << endl;
	cout << "HBOND:              " << hbondEnergy << endl;
	sout << "Initial Sequence:   " << initialSeq << endl;
	sout << "Best Sequence:      " << startSequence << endl;
	sout << "Interface Sequence: " << seqInterface << endl;
	sout << "Interface:          " << variablePositionString << endl;
	sout << "Total Energy:       " << bestEnergy << endl;
	sout << "VDW:                " << vdwEnergy << endl;
	sout << "HBOND:              " << hbondEnergy << endl;
	sout << endl << "End SelfPairManager Optimization: " << spmTime << "s" << endl;
	cout << endl << "End SelfPairManager Optimization: " << spmTime << "s" << endl;

	/******************************************************************************
	 *           === METHODS FOR DETERMINING ALTERNATE SEQUENCES ===
	 ******************************************************************************/
	// Initialize vector to hold designed sequences
	vector<string> seqs;

	// Initialize vector to hold all designed sequences from the stateMC
	vector<string> allSeqs;

	//Initialize energyMap to hold all energies for output into a summary file
	map<string, map<string,double>> sequenceEnergyMap;

	// Initialize sequence and state pair vector: each sequence will be tied to it's state with the proper rotamers
	vector<pair<string,vector<uint>>> sequenceStatePair;

	/******************************************************************************
	 *      === MONTE CARLO TO RANDOMIZE SEQUENCES FROM BEST SCMF STATE ===
	 ******************************************************************************/
	// Unlink the best state from SCMF if not using linked positions during the state Monte Carlo
	// Will likely need to change this function for heterodimers but try and use it first
	// TODO: should be working now, but going to test in lab tomorrow. need to attach an example config file to my github
	unlinkBestState(opt, bestState, rotamerSamplingPerPosition, opt.backboneLength);
	stateMCUnlinked(sys, opt, PL, sequenceEnergyMap, sequenceEntropyMap, bestState, seqs, allSeqs,
			sequenceStatePair, allInterfacePositions, interfacePositions, rotamerSamplingPerPosition, RNG, sout, err);

	// I moved this on 12_6_2021: I noticed that I am getting too many threonines in my starting seqeunces, so I decided to move the baseline down here (this likely won't work for stateMCLinked, so I'll have to deal with that

	//Add energies for initial sequences into the sequenceEnergyMap
	seqs.insert(seqs.begin(), startSequence);
	pair<string,vector<uint>> startSequenceStatePair = make_pair(startSequence, startingState);
	sequenceStatePair.insert(sequenceStatePair.begin(), startSequenceStatePair);
	getEnergiesForStartingSequence(opt, spm, startSequence, startingState, sequenceEnergyMap, sequenceEntropyMap);
	getSasaForStartingSequence(sys, startSequence, startingState, sequenceEnergyMap);

	/******************************************************************************
	 *            === CALCULATE MONOMER ENERGIES OF EACH SEQUENCE ===
	 ******************************************************************************/
	computeMonomerEnergies(opt, trans, sequenceEnergyMap, seqs, RNG, sout, err);
	getSasaDifference(sequenceStatePair, sequenceEnergyMap);

	/******************************************************************************
	 *              === CALCULATE TOTAL ENERGIES AND WRITE PDBS ===
	 ******************************************************************************/
	// Initialize PDBWriter for designs
	PDBWriter writer;
	writer.open(opt.pdbOutputDir + "/allDesigns.pdb");

	// Output these energy calculations to the summary file
	sout << "Calculating Final Energies..." << endl;
	for (uint i=0; i<sequenceStatePair.size(); i++){
		string sequence = sequenceStatePair[i].first;
		vector<uint> state = sequenceStatePair[i].second;
		sout << "Sequence " << i+1 << ": " << sequence << endl;
		int seqNumber = i;
		// calculates the total energy difference between monomer and dimer and outputs individual pdbs for each sequence
		getTotalEnergyAndWritePdbs(sys, opt, sequenceEnergyMap, sequence, state, rotamerSamplingPerPosition, RNG, seqNumber, writer, sout, err);
	}
	writer.close();

	/******************************************************************************
	 *                   === WRITE OUT ENERGY AND DESIGN FILES ===
	 ******************************************************************************/
	outputDesignFiles(opt, rotamerSamplingString, rotamerSamplingPerPosition, sequenceStatePair, sequenceEnergyMap, densities);

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;

	err.close();
	sout.close();
}
