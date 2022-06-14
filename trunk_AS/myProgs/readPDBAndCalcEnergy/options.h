#ifndef OPTIONS_H
#define OPTIONS_H

#include <sstream>
#include <vector>

using namespace std;

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/
 struct Options{
	// input files
	string backboneCrd; //initial coordinates for helix backbones: crd file
	string pdbOutputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string pdbFile; //initial coordinates for helix backbone: pdb file
	string pdbName; //name of pdb for naming directory and other output files

	//
	int seed;
	int greedyCycles;
	
	//monomer options
	int MCCycles;
	int MCMaxRejects;

    //booleans
	bool verbose; //TRUE: write energy outputs and calculations throughout the run to the terminal OR FALSE: only write outputs to output files
	bool deleteTerminalHbonds; //TRUE: delete hydrogen bonds at the termini of sequences to not be considered in hydrogen bonding score OR FALSE: keep hydrogen bonds at termini

	// energy weights
	double weight_vdw; //weight of vdw energy contribution to total energy: default = 1
	double weight_hbond;//weight of hbond energy contribution to total energy: default = 1
	double weight_solv;//weight of solvation energy contribution to total energy: default = 1

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

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
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
 };
 # endif