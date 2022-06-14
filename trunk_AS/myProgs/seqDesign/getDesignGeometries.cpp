#include <sstream>
#include <iterator>

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

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "getDesignGeometries";
string programDescription = "Outputs a weighted geometry file for seqDesign.cpp to choose a design geometry";
string programAuthor = "Gilbert Loiseau";
string programVersion = "0";
string programDate = "27 September 2021";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************
 *  
 *  ======= OPTIONS =======
 *
 ******************************************/
struct Options{
	string configfile;
	string angleDistKdeFile;
	string zShiftKdeFile;
	string axialRotationKdeFile;
	string outputDir;

	// Values for bin size of each geometry: helps to make the code explore the entire geometric space rather than fixed bin points
	double angleBinSize = 5.128205;
	double xShiftBinSize = 0.26087;
	double axialRotationBinSize = 2.564103;
	double zShiftBinSize = 0.26087;

	// Limits for density to choose geometries from
	double angleDistDensityLimit;
	double zShiftDensityLimit;
	double axialRotationDensityLimit;

	int numberOfGeometries;
	int seed;
	
	bool verbose;
	bool useTimeBasedSeed;
	
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
};

/******************************************
 *  
 *  ======= FUNCTIONS =======
 *
 ******************************************/

Options parseOptions(int _argc, char * _argv[], Options defaults);

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void getAngleDist(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_xShifts, vector<double> &_crossingAngles, vector<double> &_densities){
	//Read Kde File
	Reader kdeReader(_opt.angleDistKdeFile);
	kdeReader.open();
	map<double,vector<double>> geometryDensityMap;
	if(!(kdeReader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.angleDistKdeFile << endl;
		exit(0);
	}

	//Break up kde file into lines
	vector<string> lines = kdeReader.getAllLines();
	for (int i=1; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: KDE(double) xShift(double) crossingAngle(double) Percentile(double)";
			continue;
		}

		// Extract geometry from text file
		vector<double> geometries;
		geometries.push_back(MslTools::toDouble(tokens[1]));
		geometries.push_back(MslTools::toDouble(tokens[2]));
		
		// Extract density from text file and add density and geometry to map
		double density = MslTools::toDouble(tokens[3]);
		geometryDensityMap[density] = geometries;
	}

	// Randomly choose geometries from the map
	double xBinDiff = _opt.xShiftBinSize/2;
	double angleBinDiff = _opt.angleBinSize/2;
	for (uint i=0; i<_opt.numberOfGeometries; i++){
		// Get a geometry within the specified density limit
		map<double,vector<double>>::iterator low;
		double tmp = _RNG.getRandomDouble(0,_opt.angleDistDensityLimit);	
		low = geometryDensityMap.lower_bound(tmp);
		double density = low->first;

		double acceptCriteria = 1-density/2;
		double acceptRNG = _RNG.getRandomDouble(0,_opt.angleDistDensityLimit);
		if (acceptCriteria > acceptRNG){
		// Get the bin that chosen density occupies
		double xBin = geometryDensityMap.at(density)[0];
		double angleBin = geometryDensityMap.at(density)[1];
		
		// Define the bin with limits based on bin size
		double upperXShiftBinLimit = xBin+xBinDiff;
		double lowerXShiftBinLimit = xBin-xBinDiff;
		double upperAngleBinLimit = angleBin+angleBinDiff;
		double lowerAngleBinLimit = angleBin-angleBinDiff;
		
		// Choose a random point within the bin
		double xShift = _RNG.getRandomDouble(lowerXShiftBinLimit, upperXShiftBinLimit);
		double crossingAngle = _RNG.getRandomDouble(lowerAngleBinLimit, upperAngleBinLimit);

		//Set xShift and crossingAngle
		_xShifts.push_back(xShift);
		_crossingAngles.push_back(crossingAngle);
		_densities.push_back(density);
		} else {
			i--;
		}
	}
	kdeReader.close();
}

// This one is for zShift and axialRotation
void getAxialRotations(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_axialRotations, vector<double> &_densities){
	//Read Kde File
	Reader kdeReader(_opt.axialRotationKdeFile);
	kdeReader.open();
	map<double,vector<double>> geometryDensityMap;
	if(!(kdeReader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.axialRotationKdeFile << endl;
		exit(0);
	}

	//Break up kde file into lines
	vector<string> lines = kdeReader.getAllLines();
	for (int i=1; i<lines.size(); i++){//i=1 to skip first line
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: KDE(double) Rot1(double) Rot2(double) Percentile(double)";
			continue;
		}
		// Extract geometry from text file
		vector<double> rots;
		rots.push_back(MslTools::toDouble(tokens[1]));
		rots.push_back(MslTools::toDouble(tokens[2]));
		
		// Extract density from text file and add density and geometry to map
		double density = MslTools::toDouble(tokens[3]);
		geometryDensityMap[density] = rots;
	}
	
	// Randomly choose geometries from the map
	double axialRotationBinDiff = _opt.axialRotationBinSize/2;
	for (uint i=0; i<_opt.numberOfGeometries; i++){
		// Get a geometry within the specified density limit
		map<double,vector<double>>::iterator low;
		double tmp = _RNG.getRandomDouble(0,_opt.axialRotationDensityLimit);	
		low = geometryDensityMap.lower_bound(tmp);
		double density = low->first;

		double acceptCriteria = 1-density;
		double acceptRNG = _RNG.getRandomDouble(0,_opt.angleDistDensityLimit);
		if (acceptCriteria > acceptRNG){
		//Pick from one of the two geometries (zShift and axialRotation have two choices)
		int n = _RNG.getRandomInt(0,1);

		// Get the bin that chosen probability occupies
		double axialRotationBin = geometryDensityMap.at(density)[n];
		
		// Define the bin with limits based on bin size
		double upperAxialRotationBinLimit = axialRotationBin+axialRotationBinDiff;
		double lowerAxialRotationBinLimit = axialRotationBin-axialRotationBinDiff;
		
		// Choose a random point within the bin
		double axialRotation = _RNG.getRandomDouble(lowerAxialRotationBinLimit, upperAxialRotationBinLimit);

		//Set xShift and crossingAngle
		_axialRotations.push_back(axialRotation);
		_densities.push_back(density);
		} else {
			i--;
		}
	}	
	kdeReader.close();
}

void getZShifts(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_zShifts, vector<double> &_densities){
	//Read Kde File
	Reader kdeReader(_opt.zShiftKdeFile);
	kdeReader.open();
	map<double,vector<double>> geometryDensityMap;
	if(!(kdeReader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.zShiftKdeFile << endl;
		exit(0);
	}

	//Break up kde file into lines
	vector<string> lines = kdeReader.getAllLines();
	for (int i=1; i<lines.size(); i++){//i=1 to skip first line
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: KDE(double) geom1(double) geom2(double) Percentile(double)";
			continue;
		}
		// Extract geometry from text file
		vector<double> zShifts;
		zShifts.push_back(MslTools::toDouble(tokens[1]));
		zShifts.push_back(MslTools::toDouble(tokens[2]));

		// Extract density from text file and add density and geometry to map
		double density = MslTools::toDouble(tokens[3]);
		geometryDensityMap[density] = zShifts;
	}

	// Randomly choose geometries from the map
	double zShiftBinDiff = _opt.zShiftBinSize/2;
	for (uint i=0; i<_opt.numberOfGeometries; i++){
		// Get a geometry within the specified density limit
		map<double,vector<double>>::iterator low;
		double tmp = _RNG.getRandomDouble(0,_opt.zShiftDensityLimit);	
		low = geometryDensityMap.lower_bound(tmp);
		double density = low->first;

		double acceptCriteria = 1-density;
		double acceptRNG = _RNG.getRandomDouble(0,_opt.angleDistDensityLimit);
		if (acceptCriteria > acceptRNG){
		//Pick from one of the two geometries (zShift and axialRotation have two choices)
		int n = _RNG.getRandomInt(0,1);

		// Get the bin that chosen probability occupies
		double zShiftBin = geometryDensityMap.at(density)[n];
		
		// Define the bin with limits based on bin size
		double upperZShiftBinLimit = zShiftBin+zShiftBinDiff;
		double lowerZShiftBinLimit = zShiftBin-zShiftBinDiff;
		
		// Choose a random point within the bin
		double zShift = _RNG.getRandomDouble(lowerZShiftBinLimit, upperZShiftBinLimit);

		//Set xShift and crossingAngle
		_zShifts.push_back(zShift);
		_densities.push_back(density);
		} else {
			i--;
		}
	}
	kdeReader.close();
}

/******************************************
 *  
 *  ======= BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%Y_%m_%d",timeinfo);
	string date(buffer);
	
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
	Options opt = parseOptions(argc, argv, defaults);
	
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		exit(1);
	}

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	ofstream out;
	string cmd = "mkdir -p " + opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string outfile = opt.outputDir + "/" + date + "_geometryDensityFile.csv";
	out.open(outfile.c_str());

	/******************************************************************************
	 *               === GEOMETRY FROM DISTRIBUTION OF PDB DATA ===
	 ******************************************************************************/
	// Random Number Generator to randomly choose point in geometric space
	RandomNumberGenerator RNG;
	if (opt.useTimeBasedSeed){
		RNG.setTimeBasedSeed();
	} else {
		RNG.setSeed(opt.seed); 
	}

	// Values for bin size of each geometry: helps to make the code explore the entire geometric space rather than fixed bin points (TODO: calculate in the future)
	vector<double> xShifts;
	vector<double> crossingAngles;
	vector<double> axialRotations;
	vector<double> zShifts;
	vector<double> angleDistDensities;
	vector<double> axialRotationDensities;
	vector<double> zShiftDensities;
	getAngleDist(opt, RNG, xShifts, crossingAngles, angleDistDensities);
	getAxialRotations(opt, RNG, axialRotations, axialRotationDensities);
	getZShifts(opt, RNG, zShifts, zShiftDensities);

	out << "xShift\tcrossingAngle\taxialRotation\tzShift\tangleDistDensity\taxialRotationDensity\tzShiftDensity" << endl;
	for (uint i=0; i<opt.numberOfGeometries; i++){
		out << xShifts[i] << "\t" << crossingAngles[i] << "\t" << axialRotations[i] << "\t" << zShifts[i] << "\t" << angleDistDensities[i] << "\t" << axialRotationDensities[i] << "\t" << zShiftDensities[i] << endl;
	}

	out.close();
}

/****************************************
 *  
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults){
	
	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuration file:
	 *
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/seqDesign.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//opt.required.push_back("");
	//opt.allowed.push_back("");

	// optional
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("angleDistKdeFile");
	opt.allowed.push_back("zShiftKdeFile");
	opt.allowed.push_back("axialRotationKdeFile");

	opt.allowed.push_back("angleBinSize");
	opt.allowed.push_back("xShiftBinSize");
	opt.allowed.push_back("zShiftBinSize");
	opt.allowed.push_back("axialRotationBinSize");
	
	// Limits for density to choose geometries from
	opt.allowed.push_back("angleDistDensityLimit");
	opt.allowed.push_back("zShiftDensityLimit");
	opt.allowed.push_back("axialRotationDensityLimit");
	
	opt.allowed.push_back("numberOfGeometries");
	opt.allowed.push_back("seed");

	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("useTimeBasedSeed");
	//Begin Parsing through the options
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		opt.errorMessages += "No options given!\n";
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
		}
	}
	
	// input files
	opt.angleDistKdeFile = OP.getString("angleDistKdeFile");
	if (OP.fail()) { 
		opt.warningMessages += "angleDistKdeFile not specified, default \n";
		opt.warningFlag = true;
		opt.angleDistKdeFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_09_25_angleDistKde.txt";
	}
	opt.axialRotationKdeFile = OP.getString("axialRotationKdeFile");
	if (OP.fail()) { 
		opt.warningMessages += "rotkdeFile not specified, default \n";
		opt.warningFlag = true;
		opt.axialRotationKdeFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/designFiles/2021_09_25_axialRotationKde.txt";
	}
	opt.zShiftKdeFile = OP.getString("zShiftKdeFile");
	if (OP.fail()) { 
		opt.warningMessages += "zkdeFile not specified, default \n";
		opt.warningFlag = true;
		opt.zShiftKdeFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/designFiles/2021_09_25_zShiftKde.txt";
	}

	// output dir
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine outputDir";
		opt.errorFlag = true;
	}

	// Values for bin size of each geometry: helps to make the code explore the entire geometric space rather than fixed bin points
	opt.angleBinSize = OP.getDouble("angleBinSize");
	if (OP.fail()) {
		opt.warningMessages += "angleBinSize not specified, defaulting to 5.128205\n";
		opt.warningFlag = true;
		opt.angleBinSize = 5.128205;
	}
	opt.xShiftBinSize = OP.getDouble("xShiftBinSize");
	if (OP.fail()) {
		opt.warningMessages += "xShiftBinSize not specified, defaulting to 0.26087\n";
		opt.warningFlag = true;
		opt.xShiftBinSize = 0.26087;
	}
	opt.axialRotationBinSize = OP.getDouble("axialRotationBinSize");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 2.564103\n";
		opt.warningFlag = true;
		opt.axialRotationBinSize = 2.564103;
	}
	opt.zShiftBinSize = OP.getDouble("zShiftBinSize");
	if (OP.fail()) {
		opt.warningMessages += "zShiftBinSize not specified, defaulting to 0.26087\n";
		opt.warningFlag = true;
		opt.zShiftBinSize = 0.26087;
	}

	// Density limits to choose geometries from
	opt.angleDistDensityLimit = OP.getDouble("angleDistDensityLimit");
	if (OP.fail()) {
		opt.warningMessages += "angleDistDensityLimit not specified, defaulting to 1\n";
		opt.warningFlag = true;
		opt.angleDistDensityLimit = 1;
	}
	opt.axialRotationDensityLimit = OP.getDouble("axialRotationDensityLimit");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 1\n";
		opt.warningFlag = true;
		opt.axialRotationDensityLimit = 1;
	}
	opt.zShiftDensityLimit = OP.getDouble("zShiftDensityLimit");
	if (OP.fail()) {
		opt.warningMessages += "zShiftDensityLimit not specified, defaulting to 1\n";
		opt.warningFlag = true;
		opt.zShiftDensityLimit = 1;
	}
	//Number of Geometries
	opt.numberOfGeometries = OP.getInt("numberOfGeometries");
	if (OP.fail()) {
		opt.warningMessages += "numberOfGeometries not specified, defaulting to 100\n";
		opt.warningFlag = true;
		opt.numberOfGeometries = 100;
	}

	// booleans
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.useTimeBasedSeed = OP.getBool("useTimeBasedSeed");
	if (OP.fail()) {
		opt.warningMessages += "useTimeBasedSeed not specified using false\n";
		opt.warningFlag = true;
		opt.useTimeBasedSeed = false;
	}

	// seed
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}
	
	return opt;
}
