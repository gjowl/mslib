#ifndef CALCPDBOPTIONS_H
#define CALCPDBOPTIONS_H

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
	string outputDirName; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string pdbFile; //input pdb file
	string helicalAxis; //file with helical axis information

	string sequence; //sequence of the protein
	int thread; //thread number on the gly69 helix: depending on the sequence length, proper thread number can be used to center the helix in the membrane

	// Monte Carlo parameters
    int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	// booleans
	bool verbose; //TRUE: write energy outputs and calculations throughout the run to the terminal OR FALSE: only write outputs to output files
	bool deleteTerminalBonds; //TRUE: delete hydrogen bonds at the termini of sequences to not be considered in hydrogen bonding score OR FALSE: keep hydrogen bonds at termini

	// use different energy parameters
	bool useIMM1;
	bool useElec;
	bool compareSasa; //TRUE: use SASA to compare backboneOptimize states OR FALSE: use energy

	// repack parameters
	int greedyCycles;

	// load rotamers useSasa = false
	string SL; //number of rotamers

    // energy weights
	double weight_vdw; //weight of vdw energy contribution to total energy: default = 1
	double weight_solv; //weight of solvation energy contribution to total energy: default = 1
	double weight_elec; //weight of electrostatic energy contribution to total energy: default = 1
	double weight_hbond; //weight of hbond energy contribution to total energy: default = 1

	// alternate identities
	vector<string> alternateIds; //alternate AA identities for interfacial positions
    
    // terminal interactions
    vector<string> deleteTerminalInteractions; //terminal interactions to be considered in the energy calculation

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

	string OPErrors; //the errors from the option parser

	// run parameters
	string rerunConf; // data for a configuration file that would rerun the job as the current run
	string configfile;
	int seed;
};

#endif

