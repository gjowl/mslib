#ifndef CATMOPTIONS_H
#define CATMOPTIONS_H

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
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
#include "HelixDimerCluster.h"
#include "hbondInfo.h"

using namespace std;

/******************************************
 *
 *  =======  OPTIONS =======
 *
 ******************************************/

struct catmOptions {
	// Required
	string fullSequence;

	// optional
	string backboneCrd;
	string pdbOutputDir;

	int tmStart;
	int tmEnd;

	string helixGeoFile;

	string rulesFile;

	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string monoRotLibFile;

	// side-chain repack variable (optional)
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;
	double deltaXLimit; 
	double deltaCrossLimit;
	double deltaAxLimit;
	double deltaZLimit;
	bool decreaseMoveSize;

	bool verbose;
	int greedyCycles;
	int seed;

	int numberOfStructuresToMCRepack;
	double energyCutOff;

	// protein information (optional)
	string uniprotName;
	string uniprotAccession;

	// input monomerEnergy
	bool inputMonomerE;
	double monoE_vdw;
	double monoE_solv;
	double monoE_solvRef;
	double monoE_hbond;

	// clustering options (optional)
	double rmsdCutoff;
	bool clusterSolutions;
	bool printAllCrds;
	bool printAxes;
	bool printTermEnergies;

	int fullSequenceStart;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;

	int threadStart;
	int threadEnd;

	bool deleteTerminalBonds;
	vector<string> deleteTerminalInteractions;

	// weights
	double weight_elec;
	double weight_vdw;
	double weight_hbond;
	double weight_solv;

	int hbondCheckNumber;
	bool onlySaveNegativeStructures;

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