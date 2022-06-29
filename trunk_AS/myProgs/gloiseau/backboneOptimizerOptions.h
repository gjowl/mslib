#ifndef BACKBONEOPTIMIZEROPTIONS_H
#define BACKBONEOPTIMIZEROPTIONS_H

#include <sstream>
#include <vector>

using namespace std;

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/
 struct Options{
	// copied all from geomRepack; only keep the ones I need
	// input files
	string outputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string pdbFile; //initial coordinates for helix backbone: pdb file

	// Geometry
	double xShift;
	double crossingAngle;
	double axialRotation;
	double zShift;
	int thread;

	// tm
	int tmStart;
	int tmEnd;

	// Design parameters
	string sequence;
	string rotamerSamplingString;
	vector<uint> state;
	vector<int> rotamerSamplingVector;

	// load rotamers from SASA values (from sgfc)
	bool keepOriginalRotamer;
	std::vector<string> sasaRepackLevel;

	// Monte Carlo Parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;
	int numRepacks;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	double monomer;
	
	// Shift Size: typically just use defaults
	double deltaX; 
	double deltaCross;
	double deltaAx;
	double deltaZ;

	// Other options
	bool verbose;
	int greedyCycles;
	int seed;
	int interfaceLevel;
	string runNumber;
	bool negAngle;
	bool negRot;


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
	bool useIMM1;



};

# endif