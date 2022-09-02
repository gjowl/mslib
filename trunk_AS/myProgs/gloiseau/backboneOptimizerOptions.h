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
 struct BBOptions{
	// copied all from geomRepack; only keep the ones I need
	// input files
	string outputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string pdbFile; //initial coordinates for helix backbone: pdb file
	string geometryDensityFile; //geometries to choose for design with...density

	// Geometry
	double xShift;
	double crossingAngle;
	double axialRotation;
	double zShift;
	int thread;
	int threadStart;
	int threadEnd;

	// tm
	int tmStart;
	int tmEnd;

	// Design parameters
	string sequence;

	// load rotamers from SASA values (from sgfc)
	bool keepOriginalRotamer;

	// Monte Carlo Parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	double weight_elec;
	
	// input monomerEnergy
	double monomer;
	
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
	bool useAlaAtCTerminus;
	string backboneFile;
	string helicalAxis;
	int backboneLength;
	vector<string> deleteTerminalInteractions;
	string uniprotAccession;
	bool dockHelices;
	vector<double> crossAngle;
	bool getRandomAxAndZ;
	double energyCutoff;

	// 

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