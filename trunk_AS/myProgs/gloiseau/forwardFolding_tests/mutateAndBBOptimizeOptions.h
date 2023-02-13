#ifndef MUTATEANDBBOPTIMIZEOPTIONS_H 
#define MUTATEANDBBOPTIMIZEOPTIONS_H

#include <sstream>
#include <vector>

using namespace std;

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/
 struct MBBOptions{
	// input files
	string outputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string inputPdbFile; //initial coordinates for helix backbone: pdb file

	// Geometry
	double xShift;
	double crossingAngle;
	double axialRotation;
	double zShift;
	int thread;

	// Monte Carlo Parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;
    int convergedSteps;
    int convergedE;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	double weight_elec;
	
	// Shift Size: typically just use defaults
	double deltaX; 
	double deltaCross;
	double deltaAx;
	double deltaZ;
	double deltaXLimit; 
	double deltaCrossLimit;
	double deltaAxLimit;
	double deltaZLimit;
	bool decreaseMoveSize;

	// Other options
	bool verbose;
	int greedyCycles;
	int seed;
	int interfaceLevel;
	string runNumber;
	bool negAngle;
	bool negRot;
	
	// load rotamers
	string SL; //number of rotamers

	//version 2: I think it should work for both versions
	bool useElec;
	bool useIMM1;
	bool useAlaAtTermini;
	string helicalAxis;
	int backboneLength;
	vector<string> energyTermList;
	vector<string> deleteTerminalInteractions;

    vector<string> Ids;
    string interface;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help //figure out how to add stuff to help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
};

# endif