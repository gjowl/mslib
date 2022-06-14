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
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "readPDBAndCalcEnergy"; //I think this should get the filename?
string programDescription = "Reads a PDB into MSL and then calculates the energy";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "19 April 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

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
	// monomer file output
	ofstream mout;
	// error file output
	ofstream err;
	// rerun config output
	ofstream rerun;

	setupOutputDirectory(opt);

	string soutfile = opt.pdbOutputDir + "/energy.csv";
	string moutfile = opt.pdbOutputDir + "/monomer.out";
	string errfile  = opt.pdbOutputDir + "/errors.out";
	string rerunfile = opt.pdbOutputDir + "/rerun.config";

	sout.open(soutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	rerun << opt.rerunConf << endl;
	rerun.close();

	err << date << endl;
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

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

    CSB.buildSystemFromPDB(opt.pdbFile);
	CSB.setBuildNonBondedInteractions(false);

    /******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// assign the coordinates of our system to the given geometry that was assigned without energies using System pdb
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
	int firstPos = 0;
    int lastPos = sys.positionSize();
    deleteTerminalHydrogenBondInteractions(sys,firstPos,lastPos);
	cout << "Energy: " << sys.calcEnergy() << endl;
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// check to verify that all atoms have coordinates
	// See if this is integratable into the system that is built previously for the helix of interest
	checkIfAtomsAreBuilt(sys, err);

    //TODO: add in a check to see if parallel or antiparallel (distance of one termini to another > 20 angstroms or something)
    // if so, can still calculate, but make a note of it
	
    //loading rotamers
	loadRotamers(sys, sysRot, "SL95.00");
	CSB.updateNonBonded(10,12,50);//This for some reason updates the energy terms and makes the IMM1 terms active (still need to check where, but did a couple of calcEnergy and outputs

//TODO: make below calculation of energy into a function
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
	
	cout << 3 << endl;
    if (opt.verbose){
        map<string,double> energyByTerm;
	    energyByTerm = getEnergyByTerm(sys.getEnergySet()); // must double the energy, as only computed energy for 1 helix
	    for(map<string,double>::iterator it = energyByTerm.begin(); it != energyByTerm.end(); it++) {
	    	cout << it->first << " " << it->second << endl;
	    }
    }

	cout << 4 << endl;
    map<string,double> monomerEnergyByTerm;
    double monomer = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);
	double finalEnergy = dimer-monomer;
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

	sout << opt.pdbName << ',' << finalEnergy << ',' << dimer << ',' << monomer << ',' << dimerDiff << ',' << vdw << ',' << vdwMonomer << ',' << vdwDiff << ',' << hbond << ',' << hbondMonomer << ',' << hbondDiff << ',' << imm1 << ',' << imm1Monomer << ',' << imm1Diff << endl;
    //outputs a pdb file for the structure (already have the pdb, but the sidechains may be in different positions now)
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(opt.pdbOutputDir + "/" + opt.pdbName + "_sideChainRepack.pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

    // close all of the output file writers
    sout.close();
    mout.close();
    err.close();
}