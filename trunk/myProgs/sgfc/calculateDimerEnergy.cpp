#include <iostream>
#include <stdio.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "SysEnv.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "HelixGenerator.h"
#include "PDBReader.h"

using namespace MSL;
using namespace std;

static SysEnv ENV;
string programName = "calculateDimerEnergy";
string programDescription = "This program accepts a PDB file with two protein chains as input and calculates the energy.  It then moves the chains 500A apart and recalculates the energy to determine the energy of dimerization.";
string programAuthor = "Samson Condon";
string programVersion = "0.0.1";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

/*========================================== BEGIN MAIN ==========================================*/
int main(int argc, char *argv[]) {

	/******************************************************************************
	 *                     === SYSTEM SETUP ===
	 ******************************************************************************/
	System sys;
	Transforms trans;
	SelfPairManager spm;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	SystemRotamerLoader sysRot (sys, ENV.getEnv("MSL_ROTLIB"));
	// Set up Charmm System Builder
	string topFile = ENV.getEnv("MSL_CHARMM_TOP");
	string parFile = ENV.getEnv("MSL_CHARMM_PAR");
	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");

	
	// Read in PDB File

	if(!CSB.buildSystemFromPDB("/exports/home/scondon/Downloads/Model2.pdb")) {
		cerr << "Unable to build system" << endl;
		exit(0);
	}

	sys.buildAllAtoms();
	
	HydrogenBondBuilder hb(sys, ENV.getEnv("MSL_HBOND_CA_PAR"));
	hb.buildInteractions(30);
	CSB.updateNonBonded(9.0,10.0,11.0);
	cout << sys << endl;
	cout << "==================================" << endl;
	cout << "         Dimer Energy             " << endl;
	cout << "==================================" << endl;
	sys.calcEnergy();
	sys.printEnergySummary();
	double dimerE = sys.calcEnergy();


	//Begin computing monomer Energy
	//Load Rotamers and set up SPM
//	
//	/******************************************************************************
//	 *                  === LOAD ROTAMERS & SET-UP SPM ===
//	 ******************************************************************************/
//	for (uint k=0; k < sys.positionSize(); k++) {
//		Position &pos = sys.getPosition(k);
//
//		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
//			if (!sysRot.loadRotamers(&pos, pos.getResidueName(), "SL95.00")) { 
//				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
//			}
//		}
//	}
//	cout << "Rotamers Loaded!" << endl;
////	spm.saveEnergiesByTerm(true);
//	spm.setSystem(&sys);	
	
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector &chainA = sys.getChain("A").getAtomPointers();
	AtomPointerVector &chainB = sys.getChain("B").getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES 500 A APART ===
	 ******************************************************************************/
	double xShift = 500.0;
	cout << "Translating chains by 500 A" << endl;

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*xShift/2.0), 0.0, 0.0);
	trans.translate(chainA, interDistVect);
	interDistVect.setCoor((xShift*2.0), 0.0, 0.0);
	trans.translate(chainB, interDistVect);
	cout << "Chains separated!" << endl;
	cout << "Repacking side chains for monomers " << endl;
	// Repack side chains
	//spm.calculateEnergies();
	//spm.runGreedyOptimizer(1);
	sys.writePdb("monomer.pdb");
	double monomerE = sys.calcEnergy();
	//double monomerE = spm.getMinBound()[0];

	double dimerizationE = dimerE - monomerE;
	cout << "Dimer Energy : " << dimerE << endl;
	cout << "Monomer Energy: " << monomerE << endl;
	cout << "Energy of Dimerization: " << dimerizationE << endl;
	
	return 1;
}
