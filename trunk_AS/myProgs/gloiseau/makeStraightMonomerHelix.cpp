#include <iostream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>
#include <chrono>

// MSL Functions
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "AtomSelection.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "MonteCarloManager.h"
#include "CRDReader.h"
#include "SysEnv.h"
#include "versatileFunctions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "makeStraightMonomerHelix";
string programDescription = "Makes a straight monomer helix for rosetta";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "28 February 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime, spmTime;
auto clockTime = chrono::system_clock::now();

// initialize time variables
time_t rawtime;
struct tm * timeinfo;
char buffer[80];

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
    string sequence = OP.getString("sequence"); // sequence of the helix
    string topFile = OP.getString("topFile"); // topology file
    string parFile = OP.getString("parFile"); // parameter file
    string solvFile = OP.getString("solvFile"); // solvent file
    string backboneFile = OP.getString("backboneFile"); // backbone file
    string helicalAxisFile = OP.getString("helicalAxisFile"); // helical axis file

    // make the polymer sequence
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(sequence, thread);
    PolymerSequence PS(polySeq);

    // load the backbone coordinates
	System gly69;
	gly69.readPdb(backboneFile,true);

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
	sys.assignCoordinates(gly69.getAtomPointers(),false);
	sys.buildAllAtoms();

    // load the helical axis coordinates
    System helicalAxis;
    helicalAxis.readPdb(helicalAxisFile);

    // setup the Transforms object
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

    // move the system to the origin
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
	writePdb(sys, outputDir, "monomer");
}