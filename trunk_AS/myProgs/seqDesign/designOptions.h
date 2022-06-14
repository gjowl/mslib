/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/13
 * @Email:  gjowl04@gmail.com
 * @Filename: design_options.h
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022-02-19
 */

#ifndef DESIGNOPTIONS_H
#define DESIGNOPTIONS_H

#include <sstream>
#include <vector>

using namespace std;

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/
 //TODO: simplify these options and have comments for each
struct Options{
	// input files
	string backboneCrd; //initial coordinates for helix backbones: crd file
	string pdbOutputDir; //output directory for all files
	string topFile; //topology file (default CHARMM22: defines distances between atoms)
	string parFile; //parameter file (defines Hbonding distances)
	string geometryDensityFile; //geometries to choose for design with...density
	string solvFile; //solvation file
	string hbondFile; //hydrogen bonding energy file
	string rotLibFile; //rotamer library file
	string infile; //initial coordinates for helix backbone: pdb file
	string selfEnergyFile; //self energy file estimates for building baseline energies
	string pairEnergyFile; //pair energy file estimates for building baseline energies
	string sequenceEntropyFile; //sequence entropy file defines average propensity of each AA from my pdb analysis
	string AACompositionPenaltyFile; //file with penalty definitions that penalizes energy of a structure based on it's AA composition

	// sequence parameters
	string backboneAA; //backbone amino acid (default to L)
	int backboneLength; //length of sequence for design (default to 21; code still needs to be reworked for other lengths; imm1 energy problem?)

	// booleans: changing these will ..TODO: add more here
	bool getGeoFromPDBData; //TRUE: randomly choose a dimeric geometry from the membrane protein pdb landscape OR FALSE: use a given dimer geometry
	bool verbose; //TRUE: write energy outputs and calculations throughout the run to the terminal OR FALSE: only write outputs to output files
	bool deleteTerminalHbonds; //TRUE: delete hydrogen bonds at the termini of sequences to not be considered in hydrogen bonding score OR FALSE: keep hydrogen bonds at termini
	bool linkInterfacialPositions; //TRUE: keep interfacial positions linked (same amino acid and rotamer) when searching for the best states in stateMC (less memory) OR FALSE: unlink positions (memory intensive)
	bool designHomodimer; //TRUE: design a homodimer sequence keeping ids same between positions on both helices OR FALSE: design a heterodimeric sequence
	bool useSasa; //TRUE: use solvent accessible surface area to designated the number of rotamers at each position on the dimer OR FALSE: input set number of rotamers for both interface and non-interface
	bool useTimeBasedSeed; //TRUE: use time based seed for all RandomNumberGenerator functions OR FALSE: use a given seed
	bool energyLandscape; //TRUE: collect all sequences and their respective monomer and dimer energies ..TODO: add more here
	bool useAlaAtCTerminus; //TRUE: use ALA at C terminus of sequence FALSE: use LEU at C terminus ..TODO: do I need this?
	bool useBaseline; //TRUE: calculate and use baseline values generated as estimates of the monomer sequence OR FALSE: don't use baselines to estimate the monomer

	// repack parameters
	int greedyCycles;
	int seed;

	// load rotamers useSasa = false
	string SL; //number of rotamers
	string SLInterface; //number of rotamers for interfacial AAs

	// load rotamers useSasa = true
	std::vector<string> sasaRepackLevel; //vector of levels
	int interfaceLevel; // level for the interface
	// Example:
	/*
	The number of given levels determines how many interfacial splits there are by normalized SASA value and sorted.

	sasaRepackLevel.push_back("SL95.00")
	sasaRepackLevel.push_back("SL95.00")
	sasaRepackLevel.push_back("SL85.00")
	sasaRepackLevel.push_back("SL60.00")
	sasaRepackLevel.size() = 4
	interfaceLevel = 2

	All positions below level 2 are consired interfacial.
	Since there are 4 levels, there are 4 splits, each being 25% of the total SASA value. Interfacial positions are the positions with the highest burial
	that occur in the first two splits. SASA value is added up for the first level until it passes 25%. Those positions are considered part of level 1. Values
	are continued to add up until reaching 50%, or level 2. And so on. In this example, all the positions that end up having a SASA value in level 1 and 2 are considered interface.

	*/

	// tm: do I need these?
	int tmStart; //specifies starting residue number
	int tmEnd; //specifies ending residue number

	// the actual AAs being modeled
	int startResNum; //starting residue number
	int endResNum; //end residue number

	// Starting Geometry
	double xShift; //distance between helices
	double zShift; //position of the crossing point between helices
	double crossingAngle; //crossing angle between helices
	double axialRotation; //rotation of helices

	// crossing point
	int thread; //crossing point...(more detail about this and gly69?)

	// Monte Carlo parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	bool useIMM1;

	// energy weights
	double weight_vdw; //weight of vdw energy contribution to total energy: default = 1
	double weight_hbond;//weight of hbond energy contribution to total energy: default = 1
	double weight_solv;//weight of solvation energy contribution to total energy: default = 1
	double weight_seqEntropy;//weight of sequence entropy contribution to total energy: default = ..

	// alternate identities
	vector<string> Ids; //alternate AA identities for interfacial positions

	// state Monte Carlo Options
	int numStatesToSave; //number of sequences to save for each design run

	//SelfPairManager Options
	bool runDEESingles;
	bool runDEEPairs;
	bool runSCMF;

	// energy terms to output: maybe rid of and just default them?
	vector<string> monomerEnergyTerms;
	vector<string> monomerIMM1EnergyTerms;
	vector<string> dimerEnergyTerms;
	vector<string> energyLandscapeTerms;
	vector<string> energyTermsToOutput;
	vector<string> energyTermList;

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
	string runNumber;
};

#endif
