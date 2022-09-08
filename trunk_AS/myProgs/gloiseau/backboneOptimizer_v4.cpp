#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>

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
#include "SasaCalculator.h"
#include "backboneOptimizerFunctions.h"
#include "backboneOptimizerOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "backboneOptimizer_v3";
string programDescription = "This is an updated version of backboneOptimizer: backboneOptimizer does local backbone moves on given backbones\n\
 	from a given input geometry. Functions in this are more improved, and it uses multithreading to try multiple threaded positions.\n\
	It also acts as a forward folding program, where it will try to minimize the energy of a given sequence and geometry.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "4";
string programDate = "30 August 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

/***********************************
 *help functions
 ***********************************/
void usage();
void help(BBOptions defaults);
void outputErrorMessage(BBOptions &_opt);
void outputWarningMessage(BBOptions &_opt);
void checkOptionErrors(BBOptions &_opt);	

int main(int argc, char *argv[]){
	
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	// Start the timer
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);

	string date(buffer);
    
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	//Add in some default options that can easily be changed here
	BBOptions opt = BBParseOptions(argc, argv);
	checkOptionErrors(opt);

	/******************************************************************************
	 *                       === SETUP OUTPUTS ===
	 ******************************************************************************/
	ofstream sout;  // summary file output
	ofstream eout;  // summary file output
	ofstream mout;  // monomer file output
	ofstream err;   // error file output
	ofstream rerun; // rerun config output

	// function that defines the output directory
	setupOutputDirectory(opt);

	// setup the output files
	string soutfile = opt.outputDir + "/summary.out";
	string eoutfile = opt.outputDir + "/allEnergies.csv";
	string moutfile = opt.outputDir + "/monomer.out";
	string errfile  = opt.outputDir + "/errors.out";
	string rerunfile = opt.outputDir + "/rerun.config";

	// open the output files
	sout.open(soutfile.c_str());
	eout.open(eoutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	// write and close the rerun config file
	rerun << opt.rerunConf << endl;
	rerun.close();

	// write the header for the energy file
	eout << "Energy,Thread,xShift,CrossingAngle,AxialRotation,zShift,axialRotDensity,zShiftDensity" << endl;

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	helicalAxis.saveCoor("originState");
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	/******************************************************************************
	 *                    === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,trans,opt.crossingAngle,opt.xShift);
	helicalAxis.applySavedCoor("originState");

	// compute monomer energy
	// takes this input system geometry as a starting point for monomer calculations
	System sys;
	string polySeq = convertToPolymerSequence(opt.sequence, opt.thread);
	prepareSystem(opt, sys, startGeom, polySeq);
	// setup random number generator for monomer energy calculations
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 
    map<string,double> monomerEnergyByTerm;
    double monomerEnergy = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);

	// loop over the number of repacks to perform for each thread
	for (uint n=0; n<opt.numRepacks; n++){
		// make directory for repack summaries
		string cmd = "mkdir -p " + opt.outputDir + "/RepackSummaries_" + to_string(n);
		if (system(cmd.c_str())){
			cout << "Unable to make directory" << endl;
			exit(0);
		}
		for (uint i=opt.threadStart; i<opt.threadEnd; i++){
			vector<thread> threads;
			//cout << "Docking Thread: " << i << endl;
			for (uint j=0; j<opt.crossAngle.size();j++){
				double crossingAngle = opt.crossAngle[j]; 
				//cout << "Angle: " << crossingAngle << endl;
				threads.push_back(thread(threadDockAndRepack, ref(opt), ref(helicalAxis), i, n, crossingAngle, monomerEnergy,
				 ref(monomerEnergyByTerm), ref(eout)));
			}
			for (auto &t : threads){
				t.join();
			}
		}
	}
	
	// output the total time
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
    // close all of the output file writers
    mout.close();
    err.close();
	sout.close();
	eout.close();
}

//Functions
vector<uint> getAllInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=0; k<_opt.sequence.length(); k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=3; k<_opt.sequence.length()-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

/***********************************
 *help functions
 ***********************************/
void usage() {
	cout << endl;
	cout << "Run as :" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputWarningMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program has the following warning:" << endl;
		cerr << endl;
		cerr << _opt.warningMessages << endl;
		cerr << endl;
}
void outputErrorMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << _opt.errorMessages << endl;
		cerr << endl;
		cerr << _opt.OPerrors << endl;
		usage();
}

//TODO: finish writing up this help
void help(BBOptions defaults) {
	cout << "This program runs as:" << endl;
	cout << " % seqDesign " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --numberOfStructuresToMCRepack <int> --energyCutOff <double> --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --MCStartTemp <double> --MCEndTemp <double> --MCCurve <CONSTANT-0, LINEAR-1, EXPONENTIAL-2, SIGMOIDAL-3, SOFT-4>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --thread <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   --weight_hbond <double> --weight_vdw <double> --weight_solv <double> --weight_seqEntropy <double>" << endl;
	cout << "   --sasaRepackLevel <rotLevel> (in format SL95.00; 4 levels used by default) --interfaceLevel <int> " << endl << endl;
	cout << "Template Configuration file (copy and paste the below into a file.config and run code as bin/seqDesign --config file.config" << endl;
	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;

	cout << endl << "#Energy term weights" << endl;
	cout << setw(20) << "weight_vdw " << defaults.weight_vdw << endl;
	cout << setw(20) << "weight_hbond " << defaults.weight_hbond << endl;
	cout << setw(20) << "weight_solv " << defaults.weight_solv << endl;
	cout << endl;
}

// check through error options and exit if too many
void checkOptionErrors(BBOptions &_opt){	
	if (_opt.errorFlag) {
		outputErrorMessage(_opt);
		exit(1);
	} else if (!_opt.errorFlag && !_opt.warningFlag && _opt.errorMessages != ""){
		outputWarningMessage(_opt);
		usage();
		exit(0);
	}
}