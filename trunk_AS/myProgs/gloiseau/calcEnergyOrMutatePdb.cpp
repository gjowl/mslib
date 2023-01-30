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
#include "versatileFunctions.h"
#include "calcPdbOptions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "calcEnergyOrMutatePdb"; //I think this should get the filename?
string programDescription = "Reads a PDB into MSL and then calculates the energy";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "19 April 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

string setupOutputDirectory(string _dirName);
void outputErrorMessage(string _message, string _optionParserErrors);
void deleteTerminalBondInteractions(System &_sys, vector<string> &_deleteTerminalInteractionList);
double computeMonomerEnergy(System & _sys, Options& _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout);
Options parseCalcEnergyOptions(int _argc, char * _argv[]);

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
	//Add in some default options that can easily be changed here
	Options opt = parseCalcEnergyOptions(argc, argv);

	if (opt.errorFlag) {
		outputErrorMessage(opt.errorMessages, opt.OPErrors);
		exit(1);
	}

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	// summary file output
	ofstream sout;
	// monomer file output
	ofstream mout;
	// error file output
	ofstream err;
	// rerun config output
	ofstream rerun;

	string outputDir = setupOutputDirectory(opt.outputDirName);

	string soutfile = outputDir + "/energy.csv";
	string moutfile = outputDir + "/monomer.out";
	string errfile  = outputDir + "/errors.out";
	string rerunfile = outputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	rerun << opt.rerunConf << endl;
	rerun.close();
	/******************************************************************************
	 *                   === INITIALIZE POLYMER SEQUENCE ===
	 ******************************************************************************/
	string seq = "AAAYYLLGTLLGYLLSTLAAA";
	string polySeq = convertToPolymerSequenceNeutralPatch(seq, 23);
	//string polySeq = generatePolymerSequenceFromSequence(opt.sequence, opt.thread);
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	System pdb;
	pdb.readPdb(opt.pdbFile);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	AtomPointerVector & apv = pdb.getAtomPointers();

	Chain & chainA1 = pdb.getChain("A");
	Chain & chainB1 = pdb.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA1 = chainA1.getAtomPointers();
	AtomPointerVector & apvChainB1 = chainB1.getAtomPointers();
	
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
	helicalAxis.readPdb(opt.helicalAxis);

	// AtomPointerVector for the helical axis for each chain
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	
	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA1, apvChainB1, axisA, axisB, ori, xAxis, zAxis, 0, 0, 0, 0, trans);
	//moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
    System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
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

	// read in the pdbfile

    //CSB.buildSystemFromPDB(opt.pdbFile);
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from pdb" << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(apv,false);
	// TODO: can I mutate the system here? or do I have to do a bunch of other stuff to get it to work (read sequence, geometry, etc)
	// CSB addIdentity. Add identity to the system at the desired positions (ends for now, but maybe accept a list?)
	// check if the terminal residues are the same as the identity
	//for (uint i = 0; i < sys.chainSize(); i++){
	//	// get the chain
	//	Chain &chain = sys.getChain(i);
	//	// get the positions
	//	vector<Position*>& positions = chain.getPositions();
	//	// loop through the positions
	//	for (uint j=0; j<positions.size(); j++){
	//		if (j < 3 || j > positions.size() - 4){
	//			// get the position
	//			Position &pos = *positions[i];
	//			// loop through the alternate identity list
	//			for (uint k=0; k< opt.alternateIds.size(); k++){
	//				string id = opt.alternateIds[k];
	//				CSB.addIdentity(pos, id);
	//			}
	//		}
	//	}
	//}

    /******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// assign the coordinates of our system to the given geometry that was assigned without energies using System pdb
	sys.buildAllAtoms();
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	
	PDBWriter writer1;
	writer1.open(outputDir + "/start.pdb");
	writer1.write(sys.getAtomPointers(), true, false, true);
	writer1.close();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	CSB.updateNonBonded(10,12,50);
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);

	sys.calcEnergy();
    cout << sys.getEnergySummary() << endl;
	
	// removes all bonding near the termini of our helices for a list of interactions
    deleteTerminalBondInteractions(sys,opt.deleteTerminalInteractions);
	
	// initialize the object for loading rotamers into our system
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
    // loading rotamers
	loadRotamers(sys, sysRot, opt.SL);

	CSB.updateNonBonded(10,12,50);
	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 

	// Optimize Initial Starting Position
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
    
    repackSideChains(spm, opt.greedyCycles);
	vector<uint> stateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(stateVec);
	
	double dimer = spm.getStateEnergy(stateVec);
    cout << "Dimer Energy: " << dimer << endl;
	sys.calcEnergy();
    cout << Eset->getSummary() << endl;

    if (opt.verbose){
        map<string,double> energyByTerm;
	    energyByTerm = getEnergyByTerm(sys.getEnergySet()); // must double the energy, as only computed energy for 1 helix
	    for(map<string,double>::iterator it = energyByTerm.begin(); it != energyByTerm.end(); it++) {
	    	cout << it->first << " " << it->second << endl;
	    }
    }

    map<string,double> monomerEnergyByTerm;
    double monomer = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);
	double finalEnergy = dimer-monomer;
	cout << "Final Energy: " << finalEnergy << endl;
	// added in way to get the differences between monomer and dimer for different energy terms
	double vdw = spm.getStateEnergy(spm.getMinStates()[0], "CHARMM_VDW");
	double hbond = spm.getStateEnergy(spm.getMinStates()[0], "SCWRL4_HBOND");
	double imm1 = spm.getStateEnergy(spm.getMinStates()[0], "CHARMM_IMM1")+spm.getStateEnergy(spm.getMinStates()[0], "CHARMM_IMM1REF");
	auto it = monomerEnergyByTerm.find("CHARMM_VDW");
	double vdwMonomer = it->second;
	it = monomerEnergyByTerm.find("SCWRL4_HBOND");
	double hbondMonomer = it->second;
	it = monomerEnergyByTerm.find("CHARMM_IMM1");
	double imm1Monomer = it->second;
	it = monomerEnergyByTerm.find("CHARMM_IMM1REF");
	imm1Monomer = imm1Monomer+it->second;
	// calculate differences: TODO turn this into a function?
	double vdwDiff = vdw-vdwMonomer;
	double hbondDiff = hbond-hbondMonomer;
	double imm1Diff = imm1-imm1Monomer;
	double dimerDiff = vdwDiff+hbondDiff+imm1Diff;

	sout << "energy,dimerEnergy,monomerEnergy,vdwDiff,hbondDiff,imm1Diff" << endl;
	sout << finalEnergy << ',' << dimer << ',' << monomer << ',' << vdwDiff << ',' << hbondDiff << ',' << imm1Diff << endl;
	cout << "energy,dimerEnergy,monomerEnergy,vdwDiff,hbondDiff,imm1Diff" << endl;
	cout << finalEnergy << ',' << dimer << ',' << monomer << ',' << vdwDiff << ',' << hbondDiff << ',' << imm1Diff << endl;

    //outputs a pdb file for the structure (already have the pdb, but the sidechains may be in different positions now)
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(outputDir + "/sideChainRepack.pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

    // close all of the output file writers
    sout.close();
    mout.close();
    err.close();
	exit(0);
}

// 
string setupOutputDirectory(string _dirName){
	_dirName = string(get_current_dir_name()) + "/" + _dirName;
	string cmd = "mkdir -p " + _dirName;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	return _dirName;
}

void outputErrorMessage(string _message, string _optionParserErrors){
	cout << endl;
	cout << "The program terminated with errors:" << endl;
	cout << endl;
	cout << _message << endl;
	cout << endl;
	cout << _optionParserErrors << endl;
}

double computeMonomerEnergy(System & _sys, Options& _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout){

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
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

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
	CSBMono.setBuildTerm("CHARMM_IMM1REF", true);
	CSBMono.setBuildTerm("CHARMM_IMM1", true);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

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

	deleteTerminalBondInteractions(monoSys, _opt.deleteTerminalInteractions);

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
    monoSpm.runGreedyOptimizer(_opt.greedyCycles);
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getStateEnergy(stateVec)*2;
	cout << "Monomer Energy w/ IMM1: " << monomerEnergy << endl;
	_mout << monoEset->getSummary();

	// print the monomer
	//string monoOutCrdFile  = _opt.outputDir + "/monomer.crd";
	//CRDWriter monoCrd;
	//monoCrd.open(monoOutCrdFile);
	//if(!monoCrd.write(monoSys.getAtomPointers())) {
	//	cerr << "Unable to write " << monoOutCrdFile << endl;
	//	exit(0);
	//}

	//string monoOutPdbFile  = _opt.outputDir + "/monomer.pdb";
	//PDBWriter monoPdb;
	//monoPdb.setConvertFormat("CHARMM22","PDB2");
	//monoPdb.open(monoOutPdbFile);
	//if(!monoPdb.write(monoSys.getAtomPointers())) {
	//	cerr << "Unable to write " << monoOutPdbFile << endl;
	//	exit(0);
	//}

	// Store monomer energy by term
	//if(_opt.verbose) {
	//	monoSys.calcEnergy();
	//	_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

	//	for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
	//		_mout << it->first << " " << it->second << endl;
	//	}
	//}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

Options parseCalcEnergyOptions(int _argc, char * _argv[]){

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file
	 ******************************************/

	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuration file:
	 *
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/seqDesign.config
	 *
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	// optional
	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("deleteTerminalBonds");
	opt.allowed.push_back("deleteTerminalInteractions");
	opt.allowed.push_back("SL");

	//Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDirName");
	opt.allowed.push_back("pdbFile");
	opt.allowed.push_back("helicalAxis");
	opt.allowed.push_back("configfile");

	//
	opt.allowed.push_back("seed");
	opt.allowed.push_back("greedyCycles");
	
	//monomer options
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");

	//Begin Parsing through the options
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		opt.errorMessages += "No options given!\n";
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPErrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	opt.deleteTerminalBonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalBonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalInteractions = OP.getMultiString("deleteTerminalInteractions");
	if (OP.fail()) {
		opt.deleteTerminalInteractions.push_back("SCWRL4_HBOND");
		opt.warningMessages += "deleteTerminalInteractions not specified, defaulting to delete SCWRL4_HBOND\n";
		opt.warningFlag = true;
	}
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	
	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL95.00\n";
		opt.SL = "SL95.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}

    //Weights
	opt.weight_vdw = OP.getDouble("weight_vdw");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_vdw not specified, default 1.0\n";
		opt.weight_vdw = 1.0;
	}
	opt.weight_hbond = OP.getDouble("weight_hbond");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_hbond not specified, default 1.0\n";
		opt.weight_hbond = 1.0;
	}
	opt.weight_solv = OP.getDouble("weight_solv");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_solv not specified, default 1.0\n";
		opt.weight_solv = 1.0;
	}

    //Parameter files
	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(envVar)) {
			opt.topFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.parFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_SOLV";
		if(SYSENV.isDefined(envVar)) {
			opt.solvFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(SYSENV.isDefined(envVar)) {
			opt.rotLibFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + ", defaulting to " + SYSENV.getEnv(envVar) + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}

	opt.hbondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hbondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hbondFile not specified using " + opt.hbondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hbondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages += "pdbFile not specified";
		opt.errorFlag = true;
	}

	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.errorMessages += "helicalAxis file not specified";
		opt.errorFlag = true;
	}
	opt.outputDirName = OP.getString("outputDirName");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine outputDirName";
		opt.errorFlag = true;
	}

	// alternate identities
	opt.alternateIds = OP.getStringVector("alternateIds");
	if (OP.fail()) {
		opt.warningMessages += "No alternate AA identities given. If error, make sure they are space separated\n";
		opt.warningFlag = true;
	}
	// Monomer Options
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "MCCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.MCCycles = 10;
	}
	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "MCMaxRejects not specified using 2\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 2;
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
