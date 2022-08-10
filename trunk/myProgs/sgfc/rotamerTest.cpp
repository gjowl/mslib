/* This program will take a residue and switch between rotamers */

#include <EnergySet.h>
#include <RotamerLibrary.h>
#include <System.h>
#include <HelixGenerator.h>
#include <CharmmSystemBuilder.h>
#include <SystemRotamerLoader.h>

using namespace MSL;
using namespace std;


int main () {
	//Files to use
		string topFile = "/exports/home/scondon/mslib/trunk/toppar/charmm22.top";
		string parFile = "/exports/home/scondon/mslib/trunk/toppar/charmm22.par";
		string solvFile = "/library/charmmTopPar/solvpar22.inp";
		string hBondFile = "/exports/home/scondon/mslib/trunk/toppar/scwrl4hb/par_hbond_CA_2.txt";
		string bbqTable = "/exports/home/scondon/mslib/trunk/tables/PeptidesBBQTable.txt";
		string rotLibFile = "/exports/home/scondon/mslib/trunk/library/EBL_11-2011_CHARMM22.txt";
	
	//Make an alpha helix from sequence
	System sys;
	PolymerSequence polSeq ("A: SER THR LEU LEU PHE TYR LYS");
	CharmmSystemBuilder charmm (sys,topFile, parFile);
	charmm.buildSystem(polSeq);
	
	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
		exit(0);
		}
	sys.buildAtoms();
	AtomPointerVector sysPointers = sys.getAtomPointers();
	cout << sys << endl;
/*
	System sys;

	if (sys.readPdb("/data04/scondon/TM_coevolution/experiments/2014-01-21/winnerMC-0001.pdb")) {
		cout << "PDB read successfully!" << endl;
	}

	CharmmSystemBuilder CSB (sys, topFile, parFile);
	if (CSB.buildSystemFromPDB("/data04/scondon/TM_coevolution/experiments/2014-01-21/winnerMC-0001.pdb")) {
		cout << "PDB built successfully!" << endl;
	}
*/
	SystemRotamerLoader SRL;

	RotamerLibrary rotLib;
	if (rotLib.readFile(rotLibFile)) {
		cout << "Rot Lib Read successfully!" << endl;
	}
	SRL.setSystem(sys);
	SRL.setRotamerLibrary(&rotLib);

	if (!SRL.loadRotamers ("A,3", "LEU", 10)) {
		cout << "Rotamers loading unsuccessful!" << endl;
	}

	EnergySet *energy = sys.getEnergySet();
	energy->calcEnergy();
	energy->printSummary();

	for (uint i = 0; i<sys.getTotalNumberOfRotamers("A,3"); i++) {
		sys.setActiveRotamer("A,3", i);
		energy->calcEnergy();
		energy->printSummary();
		if (sys.writePdb("rotTest_" + MslTools::intToString(i) + ".pdb") ) {
			cout << "rotTest_" + MslTools::intToString(i) + ".pdb" << endl;
			cout << "Pdb written successfully" << endl;
		}
	}



	/* Begin pseudo-greedy algorithm: loop through each position in Chain A and cycle through the rotamers.  If the rotamer energy is lower, set that rotamer as active.  When all the rotamers have been tested move to the next position*/



	return 0;
}
