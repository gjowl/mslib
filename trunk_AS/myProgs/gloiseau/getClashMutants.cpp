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
#include "SasaCalculator.h"
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
#include "versatileFunctions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = ""; 
string programDescription = "Reads a PDB into MSL, mutates the interfacial positions to a specified AA (typically Ile), and calculating the energy. Used for identifying clash mutations.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "13 Feb 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

double getSasaAtPosition(System &_pdb, vector<string> _chainIds, int _position, double &_totalSasa){
    double positionSasa = 0;
	SasaCalculator sasa(_pdb.getAtomPointers());
	sasa.calcSasa();
    for (uint i=0; i<_chainIds.size(); i++){
        // get the residue on the chain
        string chainResi = _chainIds[i] + ',' + MslTools::intToString(_position);
        // calculate the SASA of the mutated pdb
        double resiSasa = sasa.getResidueSasa(chainResi);
        positionSasa += resiSasa;
    }
    // calculate the total SASA of the mutated pdb
    _totalSasa = sasa.getTotalSasa();
    return positionSasa;
}

void setAminoAcidAtPosition(System &_pdb, vector<Chain*> _chains, int _position, int _chainPosition, string _aa){
    for (uint j=0; j<_chains.size(); j++){
        // get the position on the chain
        Position& pos = _chains[j]->getPosition(_position);
        string chain = _chains[j]->getChainId();
        string posId = chain+','+MslTools::intToString(_chainPosition);
        // switch the identity to alanine
        _pdb.setActiveIdentity(posId,_aa);
    }
}

void setupDirectory(string &_outputDir){
	_outputDir = string(get_current_dir_name()) + "/" + _outputDir;
	//_opt.outputDir = "/exports/home/gloiseau/mslib/trunk_AS/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

void computeMonomerEnergy(System &_sys, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, ofstream &_sout, string _topFile,
    string _parFile, string _solvFile, string _helicalAxisFile, string _backboneFile, string _rotLibFile, string _hbondFile, string _outputDir,
    int _greedyCycles, string _SL, vector<string> _deleteTerminalInteractions, vector<string> _energyTermList, int _thread);

int main(int argc, char *argv[]){
    // parse through the command line arguments
    OptionParser OP;
	OP.readArgv(argc, argv);
	OP.autoExtendOptions();
	// check if there is a config file
	if(OP.getString("config") != "NA"){
		OP.readFile(OP.getString("config"));
	}
    int thread = OP.getInt("thread"); // thread of the helix
	string pdbFile = OP.getString("pdbFile"); // input pdb file
    string topFile = OP.getString("topFile"); // topology file
    string parFile = OP.getString("parFile"); // parameter file
    string solvFile = OP.getString("solvFile"); // solvent file
    string hbondFile = OP.getString("hbondFile"); // hbond file
    string rotLibFile = OP.getString("rotLibFile"); // rotamer library file
    string backboneFile = OP.getString("backboneFile"); // backbone file
    string helicalAxis = OP.getString("helicalAxis"); // helical axis file
    string outputDir = OP.getString("outputDir"); // output directory
    string mutantAA = OP.getString("mutantAA"); // mutant amino acid
    vector<string> deleteTerminalInteractions = OP.getStringVector("deleteTerminalInteractions"); // delete terminal bond interactions
	vector<string> energyTermList = OP.getStringVector("energyTermList");
    vector<int> positionList = OP.getIntVector("positionList"); // list of positions to mutate
    int greedyCycles = OP.getInt("greedyCycles"); // number of greedy cycles
    string SL = OP.getString("SL"); // sampling level

	// summary file output
	ofstream sout;

	// output the command line arguments
    cout << "Command line arguments:" << endl;
	cout << "thread: " << thread << endl;
	cout << "pdbFile: " << pdbFile << endl;
	cout << "topFile: " << topFile << endl;
	cout << "parFile: " << parFile << endl;
	cout << "solvFile: " << solvFile << endl;
	cout << "backboneFile: " << backboneFile << endl;
	cout << "helicalAxisFile: " << helicalAxis << endl;
	
    // setup the output directory
    setupDirectory(outputDir);

    // read in the input pdb file
    System startGeom;
    startGeom.readPdb(pdbFile);

	// read in the input pdb file
    System pdb;
    CharmmSystemBuilder CSB(pdb,topFile,parFile,solvFile);
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
    CSB.buildSystemFromPDB(pdbFile);

	string startSequence = extractSequence(pdb);

    // save the starting state of the pdb (already repacked, don't need to load energy terms for another repack)
    pdb.saveAltCoor("start");

    // get the chains from the pdb
    vector<Chain*> chains = pdb.getChains();
    vector<Position*> positions = pdb.getPositions();

    // get the start position and convert to a number
    string startPosition = positions[0]->getPositionId(1);
    int startPos = MslTools::toInt(startPosition);

    // get the chains from the pdb
    for (uint i=0; i<positions.size(); i++){
        // get the position
        Position* pos = positions[i];
        string posId = pos->getPositionId();
        // add identity to the position
        Residue prevResi = positions[i]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        // if the first or last position of a chain, remove the identity (saw issues with unbuilt atoms in the center of the pdb at 0,0,0)
        if (i == 0 || i == chains[0]->positionSize()-1 || i == chains[0]->positionSize() || i == positions.size()-1){
            CSB.removeIdentity(posId,resi);
            CSB.addIdentity(posId,resi);
        } else {
            CSB.addIdentity(posId,mutantAA);
        }
    }
    pdb.assignCoordinates(startGeom.getAtomPointers(), true);
    pdb.buildAllAtoms();
	string newstartSequence = extractSequence(pdb);
    cout << "sequence after add identity: " << newstartSequence << endl;

	// Add hydrogen bond term
	HydrogenBondBuilder hb(pdb, hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct
	
    // initialize the object for loading rotamers into system
	SystemRotamerLoader sysRot(pdb, rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	loadRotamers(pdb, sysRot, SL);
	newstartSequence = extractSequence(pdb);
    cout << "sequence after load rotamers: " << newstartSequence << endl;

	setActiveSequence(pdb, startSequence);
	newstartSequence = extractSequence(pdb);
    cout << "sequence after load rotamers: " << newstartSequence << endl;
	CSB.updateNonBonded(10,12,50);
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = pdb.getEnergySet();
	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	Eset->setWeight("CHARMM_IMM1", 1);

	pdb.calcEnergy();
	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
    deleteTerminalBondInteractions(pdb,deleteTerminalInteractions);

    // write the pdb
	PDBWriter writer;
	writer.open(outputDir + "/" + startSequence + ".pdb");
	writer.write(pdb.getAtomPointers(), true, false, false);
    writer.close();

	// initialize the sequenceEnergyMap
	map<string, map<string, double>> sequenceEnergyMap;
    
    // initialize the random number generator
	RandomNumberGenerator RNG;
	RNG.setSeed(0); 

    // get the chain ids
    vector<string> chainIds;
    for (uint i=0; i<chains.size(); i++){
        chainIds.push_back(chains[i]->getChainId());
    }

	// first compute the energy of the initial pdb
	// Optimize Initial Starting Position
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&pdb);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	// Set a mask and run a greedy to get the best state for the current sequence
	vector<vector<bool>> mask = getActiveMask(pdb);
	spm.runGreedyOptimizer(100, mask);
	vector<uint> stateVec = spm.getMinStates()[0];
	pdb.setActiveRotamers(stateVec);

	// calculate the energy and add it to the sequenceEnergyMap
	double dimerEnergy = pdb.calcEnergy();
	sequenceEnergyMap[startSequence]["Dimer"] = dimerEnergy;

	// loop through the energy terms
	for (uint i=0; i<energyTermList.size(); i++){
		// get the ith energy term
		string energyTerm = energyTermList[i];
		// get the energy for the ith energy term
		double energy = Eset->getTermEnergy(energyTerm);
		// add the term to the energy map
		sequenceEnergyMap[startSequence][energyTerm] = energy;
	}

	// print the energy
	pdb.printEnergySummary();

	// compute the monomer energy of the initial pdb
	computeMonomerEnergy(pdb, sequenceEnergyMap, startSequence, sout, topFile, parFile, solvFile, helicalAxis, backboneFile, rotLibFile, hbondFile, outputDir, greedyCycles,
        SL, deleteTerminalInteractions, energyTermList, thread);

	// loop through the positions to mutate and calculate the energy
    for (uint i=0; i < positionList.size(); i++){
        cout << positionList.size() << endl;
        cout << "Position: " << positionList[i] << endl;
        int pos = positionList[i];
        int chainPos = pos+startPos;
		cout << "Chain Position: " << chainPos << endl;
        // get the previous identity of the position
        Residue prevResi = positions[pos]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        cout << "Position " << pos << " is " << resi << endl;
        // switch the identity to the mutant amino acid
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, mutantAA);
        pdb.buildAllAtoms();

        Residue currResi = positions[pos]->getCurrentIdentity();
        string resi1 = currResi.getResidueName();
        double posSasa = positions[pos]->getSasa();
        cout << "Position " << pos << " is " << resi1 << endl;
        // get the current sequence
	    string sequence = extractSequence(pdb);
        cout << "Sequence: " << sequence << endl;

		// Set a mask and run a greedy to get the best state for the current sequence
		vector<vector<bool>> mask = getActiveMask(pdb);
		spm.runGreedyOptimizer(100, mask);
	    vector<uint> stateVec = spm.getMinStates()[0];
	    pdb.setActiveRotamers(stateVec);

	    // write the pdb
	    PDBWriter writer;
	    writer.open(outputDir + "/" + sequence + ".pdb");
	    writer.write(pdb.getAtomPointers(), true, false, false);
        writer.close();

	    // calculate the energy and add it to the sequenceEnergyMap
	    double dimerEnergy = pdb.calcEnergy();
	    sequenceEnergyMap[sequence]["Dimer"] = dimerEnergy;

	    // loop through the energy terms
	    EnergySet *Eset = pdb.getEnergySet();
	    for (uint i=0; i<energyTermList.size(); i++){
	    	// get the ith energy term
	    	string energyTerm = energyTermList[i];
	    	// get the energy for the ith energy term
	    	double energy = Eset->getTermEnergy(energyTerm);
	    	// add the term to the energy map
	    	sequenceEnergyMap[sequence][energyTerm] = energy;
	    }

	    // print the energy
	    pdb.printEnergySummary();
        cout << "Pdb chain size:" << pdb.chainSize() << endl;

	    // compute the monomer energy of the initial pdb
	    computeMonomerEnergy(pdb, sequenceEnergyMap, sequence, sout, topFile, parFile, solvFile, helicalAxis, backboneFile, rotLibFile, hbondFile, outputDir, greedyCycles,
            SL, deleteTerminalInteractions, energyTermList, thread);

        // switch the identity back to the original amino acid
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, resi);
        pdb.buildAllAtoms();
    }
    
    // write the sasa map to a file
    ofstream sequenceEnergyFile;
    sequenceEnergyFile.open(outputDir + "/" + "energyFile.txt");
    // write the header by iterating through the first sequence
    map<string, double> firstSequence = sequenceEnergyMap.begin()->second;
    sequenceEnergyFile << "Sequence,";
    for (auto it2=firstSequence.begin(); it2!=firstSequence.end(); it2++){
        sequenceEnergyFile << it2->first << ",";
    }
    sequenceEnergyFile << endl;
    // write the output sasa values
    for (auto it=sequenceEnergyMap.begin(); it!=sequenceEnergyMap.end(); it++){
        sequenceEnergyFile << it->first << ",";
        for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
            sequenceEnergyFile << it2->second << ",";
        }
        sequenceEnergyFile << endl;
    }
    sequenceEnergyFile.close();
}

// compute monomer energy given a sequence and a system
void computeMonomerEnergy(System &_sys, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, ofstream &_sout, string _topFile,
    string _parFile, string _solvFile, string _helicalAxisFile, string _backboneFile, string _rotLibFile, string _hbondFile, string _outputDir,
    int _greedyCycles, string _SL, vector<string> _deleteTerminalInteractions, vector<string> _energyTermList, int _thread) {

	// initialize the polymer sequence
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, _thread);	
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	
	/*******************************************
	 *       === HELICAL AXIS SET UP ===
	 *******************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(_helicalAxisFile);

	// Setup random number generator
	RandomNumberGenerator RNG;
	RNG.setSeed(0); // use random seed 

	//Chain &inputChain = test.getChain(0);
	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _topFile, _parFile, _solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);
	CSBMono.setBuildTerm("CHARMM_IMM1REF", true);
	CSBMono.setBuildTerm("CHARMM_IMM1", true);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	//if(!CSBMono.buildSystem(PS)) {
	//	cerr << "Unable to build system from " << PS << endl;
	//	exit(0);
	//}
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());
	
	// assign the coordinates of our system to the given geometry 
	//monoSys.assignCoordinates(inputChain.getAtomPointers(),false);
	//monoSys.buildAtoms();
	SystemRotamerLoader monoRot(monoSys, _rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _hbondFile);
	monohb.buildInteractions(50);
	
	CSBMono.updateNonBonded(10,12,50);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsInactive();
	monoEset->setTermActive("CHARMM_IMM1REF", true);
	monoEset->setTermActive("CHARMM_IMM1", true);
	monoEset->setTermActive("CHARMM_VDW", true);
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);
	monoEset->setWeight("CHARMM_IMM1REF", 1);
	monoEset->setWeight("CHARMM_IMM1", 1);
	
	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all bonding near the termini of our helices for a list of interactions
    deleteTerminalBondInteractions(monoSys,_deleteTerminalInteractions);
	
	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.updateWeights();
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();
    //monoSpm.runGreedyOptimizer(_opt.greedyCycles);
	repackSideChains(monoSpm, 50);
	vector<uint> testVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(testVec);
	monoSys.calcEnergy();
	monoSys.printEnergySummary();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAllAtomPointers();

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
    monoSpm.runGreedyOptimizer(_greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	helicalAxis.saveAltCoor("BestMonomerAxis");
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		trans.translate(chainA, zUnitVector);

		//double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

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
			monoSpm.runGreedyOptimizer(_greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			//_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
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

	//MonteCarloManager MCMngr(1000.0, 0.5, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		helicalAxis.applySavedCoor("BestMonomerAxis");

		int moveToPreform = RNG.getRandomInt(2);

		double deltaZShift = 0.0;
		double deltaTilt = 0.0;
		double deltaAxialRotation = 0.0;

		//======================================
		//====== Z Shift ======
		//======================================
		if (moveToPreform == 0) {
			deltaZShift = getStandardNormal(RNG) * 1.0;
			CartesianPoint translateA = axisA(1).getCoor() - axisA(0).getCoor(); // vector minus helical center
			translateA = translateA.getUnit() * deltaZShift; // unit vector of helical _axis times the amount to shift by
			trans.translate(chainA, translateA);
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(RNG) * 20.0;
			trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(RNG) * 10;
			trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_fout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_fout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_fout << "state rejected   energy: " << currentEnergy << endl;
		} else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;
		}
		counter++;
	}

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	//Calculate Monomer energy for output
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");
	monoSys.calcEnergy();
	monoSys.printEnergySummary();
	monoSpm.runGreedyOptimizer(_greedyCycles);
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getMinBound()[0]*2;
	cout << "Monomer Energy: " << monomerEnergy << endl;
	monoSys.calcEnergy();
	monomerEnergy = monoSpm.getMinBound()[0]*2;
	cout << "Monomer Energy: " << monomerEnergy << endl;

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	monoSasa.calcSasa();
	double monomerSasa = monoSasa.getTotalSasa();
	double totalMonomerSasa = monomerSasa*2;
	monoSys.printEnergySummary();

	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _energyTermList, "Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	_sequenceEnergyMap[_seq]["MonomerSasa"] = totalMonomerSasa;
	monoSys.printEnergySummary();

	AtomSelection sele(monoSys.getAtomPointers());
	//monoEset->eraseTerm("CHARMM_IMM1");
	//monoEset->eraseTerm("CHARMM_IMM1REF");
	//monoEset->eraseTerm("SCWRL4_HBOND");
	monoSys.printEnergySummary();
	for (uint i=_thread; i<_seq.length()+_thread; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sele.select(residue += number);
		double resi = monoSys.calcEnergy("resi"); 
		cout << "Residue " << i << " vdw energy: " << resi << endl;
	}
	monoSys.calcEnergy();
	monoSys.printEnergySummary();
	sele.clearStoredSelections();
	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	helicalAxis.clearSavedCoor("BestMonomerAxis");
	helicalAxis.clearSavedCoor("bestZ");
	monoSys.printEnergySummary();
	writePdb(monoSys, _outputDir, "monomer");
}
