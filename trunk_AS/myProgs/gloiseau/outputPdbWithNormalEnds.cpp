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

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "outputPdbWithNormalEnds";
string programDescription = "converts a pdb file from msl with alternate ends (ACE and ...) to a pdb file with normal ends for use in rosetta docking";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "3 March 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

int main(int argc, char *argv[]){
    // get the current working directory
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    string outputDir = string(cwd) + "/";

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
    string backboneFile = OP.getString("backboneFile"); // backbone file
    string helicalAxisFile = OP.getString("helicalAxisFile"); // helical axis file

	// output the command line arguments
    cout << "Command line arguments:" << endl;
	cout << "thread: " << thread << endl;
	cout << "pdbFile: " << pdbFile << endl;
	cout << "topFile: " << topFile << endl;
	cout << "parFile: " << parFile << endl;
	cout << "solvFile: " << solvFile << endl;
	cout << "backboneFile: " << backboneFile << endl;
	cout << "helicalAxisFile: " << helicalAxisFile << endl;

	// read in the input pdb file
	System pdb;
	pdb.readPdb(pdbFile,true);
	// get the sequence from the pdb
	string sequence = extractSequence(pdb);
    // make the polymer sequence
	string polySeq = convertToPolymerSequence(sequence, thread);
    PolymerSequence PS(polySeq);
	cout << PS << endl;

    // initialize system
    System sys;
	CharmmSystemBuilder CSB(sys,topFile,parFile,solvFile);
	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);
	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << PS << endl;
		exit(0);
	}

	// assign the coordinates of our system to the given geometry 
	sys.assignCoordinates(pdb.getAtomPointers(),false);
	sys.buildAllAtoms();

	// write the pdb
	writePdb(sys, outputDir, "dimer");
}