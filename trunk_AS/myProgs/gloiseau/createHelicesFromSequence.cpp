#include <iostream>
#include <time.h>

#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "OptionParser.h"		
#include "MslTools.h"			
#include "Transforms.h"			
#include "System.h"			
#include "SysEnv.h"			
#include "RandomNumberGenerator.h"	
#include "HelixGenerator.h"

using namespace MSL;
using namespace std;
static SysEnv SYSENV;
string programName = "createHelicesFromSequence";
string programDescription = "This program creates helices from a sequence";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "20 August 2022";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

string convertToFullPolymerSequence(std::vector<string> &_seqs, int _numChains, std::vector<int> &_chainStartNum);
void createHelixBackbone (std::vector<string> &_chainSeqs, int _numChains, std::vector<int> &_chainStartNum, AtomPointerVector &_apv, HelixGenerator &_hg);
string convertToThreeLetterSequence(string _seq);
string printParameters (std::vector <double> _params);
string transformBundleToStartingPosition (System &_sys, Transforms &_trans, int _numChains, std::vector<double> &_params);
void transformHelix (double _axialRotation, double _zShift, double _crossingAngle, double _xShift, AtomPointerVector &_helix, Transforms &_trans);

struct Options {
	//Parameters required for each set of chains
	std::vector<string> chainSeq;
	int numChains;
	std::vector<int> chainStartNum;
	std::vector<string> chainRotLevels;

	//Output Parameters
	string outputName;
	string outputDir;
	
	//Parameter Files
	string topFile;
	string parFile;
	string hbondFile;
	string rotLibFile;
	string bbqFile;

	//Orientation Parameter ranges 
	double axialRotStart;
	double axialRotEnd;
	double zShiftStart;
	double zShiftEnd;
	double crossAngleStart;
	double crossAngleEnd;
	double xShiftStart;
	double xShiftEnd;

    //Other Parameters
    int seed;

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

Options parseOptionsChfs(int _argc, char * _argv[]);
void usage();
void version();
void help();

/******************************************
 *  
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main (int argc, char *argv[]) {

	/******************************************************************************
	 *                  === SET UP OPTIONS AND BOOKKEEPING ===
	 ******************************************************************************/
	Options opt = parseOptionsChfs(argc, argv);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		usage();
		exit(1);
	}

	opt.outputDir = opt.outputDir + "/" + opt.outputName;
	string cmd = "mkdir -p " + opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}

	// Make log file if performing verbose
	string logFile = opt.outputDir + "/" + opt.outputName + ".log";
	ofstream lout;
	lout.open(logFile.c_str());
	if(!lout.is_open()) {
		cerr << "Unable to open " << logFile << endl;
		exit(0);
	}
	lout << opt.rerunConf << endl;
	
	/******************************************************************************
	 *                  === DECLARE SYSTEM ===
	 ******************************************************************************/	
	System sys;
	HelixGenerator hg(opt.bbqFile);
	//if(opt.useCatmBackbone) {
	//	hg.setHelixParameters(3.79, 90.77 / 180.0 * M_PI, 50.64 / 180.0 * M_PI);
	//}

	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);
	lout << "RNG Seed: " << RNG.getSeed() << endl;
	
	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");
	//CSB.setVdwRescalingFactor(opt.vdwRadius);
	HydrogenBondBuilder hb(sys, opt.hbondFile);

	/******************************************************************************
	 *                  === BUILD TM BUNDLE SYSTEM ===
	 ******************************************************************************/
	AtomPointerVector centered_zAlignedAPV;
	createHelixBackbone (opt.chainSeq, opt.chainSeq.size(), opt.chainStartNum, 
						centered_zAlignedAPV, hg);

	lout << "Backbone Generated!" << endl;
	
    string fullPS = convertToFullPolymerSequence(opt.chainSeq, opt.chainSeq.size(), 
			opt.chainStartNum);
	lout << fullPS << endl;
	PolymerSequence PS;
	PS.setSequence(fullPS);
	CSB.buildSystem(PS);
	lout << "Sequence built! Preparing to assign coordinates..." << endl;
	sys.wipeAllCoordinates();
	//sys.assignCoordinates(centered_zAlignedAPV, false);
	sys.assignCoordinates(centered_zAlignedAPV, true);
	sys.buildAllAtoms();
	lout << "Coordinates assigned to helix bundle!\n";

	/******************************************************************************
	 *               === CALCULATE STARTING POSITION ===
	 ******************************************************************************/
	double axialRot;
	double zShift;
	double crossingAngle;
	double xShift;

	axialRot = RNG.getRandomDouble(opt.axialRotStart, opt.axialRotEnd);
	zShift = RNG.getRandomDouble(opt.zShiftStart, opt.zShiftEnd);
	crossingAngle = RNG.getRandomDouble(opt.crossAngleStart, opt.crossAngleEnd);
	xShift = RNG.getRandomDouble(opt.xShiftStart, opt.xShiftEnd);

	vector<double> params;
	params.push_back(axialRot);
	params.push_back(zShift);
	params.push_back(crossingAngle);
	params.push_back(xShift);
	lout << "Starting Orientation Parameters: " << endl;
	lout << printParameters(params)<< endl;
	
	/******************************************************************************
	 *               === MOVE TO STARTING POSITION ===
	 ******************************************************************************/
	//Transform the bundle to initial coordinates
	transformBundleToStartingPosition (sys, trans, opt.chainSeq.size(), params);
	string startCoorFile = opt.outputDir + "/" + opt.outputName + ".pdb";
	sys.writePdb(startCoorFile);
    lout << "Wrote " << startCoorFile << endl;
    lout.close();
    exit(0);
}

string convertToFullPolymerSequence(std::vector<string> &_seqs, int _numChains, std::vector<int> &_chainStartNum) {
	string chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	string fullps = "";
	vector<string> threeLetterSeqs;
	for (unsigned int i = 0; i < _seqs.size(); i++) {
		threeLetterSeqs.push_back(convertToThreeLetterSequence(_seqs[i]));
	}
	for (unsigned int j = 0; j < _numChains; j++) {
		string chainID(1, chainIDs[j]);
		string startPos = "{" + MslTools::intToString(_chainStartNum[j]) + "}";
		fullps += chainID + ": " + startPos + " " + threeLetterSeqs[j] + "\n";
	}
	
	return fullps;
}

string convertToThreeLetterSequence(string _seq) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
	}
	return ps;
}

void createHelixBackbone (std::vector<string> &_chainSeqs, int _numChains, std::vector<int> &_chainStartNum, AtomPointerVector &_apv, HelixGenerator &_hg) {
	//For each set of chains
	for (unsigned int i = 0; i < _chainSeqs.size(); i++) {
		int chainLength = _chainSeqs[i].length();
		//For the number of chains within the set
		for (unsigned int c = 0; c < _numChains; c++) {

			//create a helix and assign its sequence
			AtomPointerVector helix;
			_hg.generateHelix(helix, chainLength, true, true);

			//Convert APV chain IDs and start numbers to match full polymer sequence
			for (unsigned int r = 0; r < helix.size(); r++) {
				string oldChain = helix[r]->getChainId();
				char newChain = int(oldChain[0]) + i;
				string chainID;
				chainID += newChain;
				helix[r]->setChainId(chainID);
				//Figure out how to do chain renumbering from this... 
				helix[r]->setResidueNumber(helix[r]->getResidueNumber() + _chainStartNum[i]);

				Atom *tmpAtom = new Atom(helix(r));
				_apv.push_back(tmpAtom);
				tmpAtom = NULL;
			}
		}
		//apply rotational symmetry to the set
	}
	
}

string printParameters (std::vector <double> _params) {
	std::vector<string> headers;
	headers.push_back("axialRot");
	headers.push_back("zShift");
	headers.push_back("crossingAngle");
	headers.push_back("xShift");

	std::stringstream ss;
	ss.precision(3);
	ss << fixed;
	ss << "ChainBundle\t";
	ss << endl;
	for (int i =0 ; i < headers.size(); i++) {
		ss << headers[i] << "\t";
		ss << _params[i] << "\t";
		ss << endl;
	}
	string niceParams = ss.str();
	return niceParams;
}

string transformBundleToStartingPosition (System &_sys, Transforms &_trans, int _numChains, std::vector<double> &_params) {
	//move helices to start position
	for (unsigned int i = 0; i < _numChains; i++) {
		//transform the first chain in the bundle
		AtomPointerVector bundleChain = _sys.getChain(i).getAtomPointers();
		transformHelix(_params[0], _params[1], _params[2], _params[3], bundleChain, _trans);
	}	
	
	double angle = 0.0;
	//apply rotational symmetry to the rest of the chains in the bundle
	for (unsigned int j = 1; j < _numChains; j++) {
		angle += 360.0/_numChains;			
		//Set coordinates of buddy chains to the first chain
		AtomPointerVector chainPrime = _sys.getChain(j).getAtomPointers();
		//Apply rotational symmetry to the buddy chain
		_trans.Zrotate(chainPrime, angle);
	}
	
	string paramString = printParameters(_params);
	return paramString;
}

void transformHelix (double _axialRotation, double _zShift, double _crossingAngle, double _xShift, AtomPointerVector &_helix, Transforms &_trans) {
	//set up coordinates for helix
	//assume that helix is parallel to z axis and centered on the z-shift
	CartesianPoint ori(0.0, 0.0, _zShift);
	CartesianPoint xax(1.0, 0.0, _zShift);
	CartesianPoint yax(0.0, 1.0, _zShift);
	CartesianPoint zax(0.0, 0.0, _zShift+1.0);

	//Transform the axial coordinates to match z shift
	
	_trans.rotate(_helix, _axialRotation, ori, zax); //axial rotation
	_trans.rotate(_helix, _crossingAngle, ori, xax); //crossing angle
	_trans.Xtranslate(_helix, -1*_xShift);	         //xshift
}


void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
//	cout << "This program runs as:" << endl;
//	cout << " % genHomoUniverse --helixPdbFile <pdbfile> --helicalAxisPdbFile <filename> --output <filename>" << endl;
//	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
//	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
//	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
//	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
//	cout << endl;
}

/******************************************
 *  
 *  =======  OPTIONS FUNCTIONS =======
 *
 ******************************************/
Options parseOptionsChfs(int _argc, char * _argv[]) {

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

	vector<string> required;
	vector<string> allowed;

	//Required
    // sequence parameters
	opt.allowed.push_back("chainSeq");
    opt.allowed.push_back("numChains");
    opt.allowed.push_back("chainStartNum");

    // output parameters
    opt.allowed.push_back("outputDir");
    opt.allowed.push_back("outputName");

    // paramter files
    opt.allowed.push_back("topFile");
    opt.allowed.push_back("parFile");
    opt.allowed.push_back("hbondFile");
    opt.allowed.push_back("rotLibFile");
    opt.allowed.push_back("bbqFile");

    // geometry parameters   
	opt.allowed.push_back("axialRotStart");
	opt.allowed.push_back("axialRotEnd");
	opt.allowed.push_back("zShiftStart");
	opt.allowed.push_back("zShiftEnd");
	opt.allowed.push_back("crossAngleStart");
	opt.allowed.push_back("crossAngleEnd");
	opt.allowed.push_back("xShiftStart");
	opt.allowed.push_back("xShiftEnd");

    // other parameters
    opt.allowed.push_back("seed");
	opt.allowed.push_back("configfile");

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
	if (opt.help) {
		help();
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

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/	
	//Required
	opt.chainSeq = OP.getStringVector("chainSeq");//for hetero, get multiple seqs from config
	if (OP.fail()) {
		opt.errorMessages += "chainSeq not specified!\n";
		opt.errorFlag = true;
	}
	opt.numChains = OP.getInt("numChains");
	if (OP.fail()) {
		opt.warningMessages += "numChains not specified, default to 1!\n";
		opt.warningFlag = true;
        opt.numChains = 1;
	}
	opt.chainStartNum = OP.getIntVector("chainStartNum");//for hetero, have multiple start numbers (may not need, could just always start them at the same number; I think this is here for her to keep track of her protein start numbers)
	if (OP.fail()) {	
		opt.warningMessages += "chainStartNum not specified, default to 1\n";
		opt.warningFlag = true;
        opt.chainStartNum.push_back(0);
	}
    if (opt.numChains > 1 && opt.chainSeq.size() != opt.numChains) {
        string chainSeq = opt.chainSeq[0];
        int chainNum = opt.chainStartNum[0];
        for (int i=1; i<opt.numChains; i++) {
            opt.chainSeq.push_back(chainSeq);
            opt.chainStartNum.push_back(chainNum);
        }
    }

    // output parameters
    opt.outputDir = OP.getString("outputDir");
    if (OP.fail()) {
        opt.errorMessages += "outputDir not specified!\n";
        opt.errorFlag = true;
    }
    opt.outputName = OP.getString("outputName");
    if (OP.fail()) {
        opt.errorMessages += "outputName not specified!\n";
        opt.errorFlag = true;
    }

    // paramter files
    opt.topFile = OP.getString("topFile");
    if (OP.fail()) {
        opt.errorMessages += "topFile not specified!\n";
        opt.errorFlag = true;
    }
    opt.parFile = OP.getString("parFile");
    if (OP.fail()) {
        opt.errorMessages += "parFile not specified!\n";
        opt.errorFlag = true;
    }
    opt.hbondFile = OP.getString("hbondFile");
    if (OP.fail()) {
        opt.errorMessages += "hbondFile not specified!\n";
        opt.errorFlag = true;
    }
    opt.rotLibFile = OP.getString("rotLibFile");
    if (OP.fail()) {
        opt.errorMessages += "rotLibFile not specified!\n";
        opt.errorFlag = true;
    }
    opt.bbqFile = OP.getString("bbqFile");
    if (OP.fail()) {
        opt.errorMessages += "bbqFile not specified!\n";
        opt.errorFlag = true;
    }
    
	//Geometric Parameters
	opt.axialRotStart = OP.getDouble("axialRotStart");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "axialRotStart not specified! Default set to -180.\n";
		opt.axialRotStart = -180;
	}
	opt.axialRotEnd = OP.getDouble("axialRotEnd");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "axialRotEnd not specified! Default set to 180.\n";
		opt.axialRotEnd = 180;
	}
	opt.zShiftStart = OP.getDouble("zShiftStart");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "zShiftStart not specified! Default set to -15.\n";
		opt.zShiftStart = -15;
	}
	opt.zShiftEnd = OP.getDouble("zShiftEnd");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "zShiftEnd not specified! Default set to 15.\n";
		opt.zShiftEnd = 15;
	}
	opt.crossAngleStart = OP.getDouble("crossAngleStart");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "crossAngleStart not specified! Default set to -50.\n";
		opt.crossAngleStart = -50;
	}
	opt.crossAngleEnd = OP.getDouble("crossAngleEnd");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "crossAngleEnd not specified! Default set to 50.\n";
		opt.crossAngleEnd = 50;
	}
	opt.xShiftStart = OP.getDouble("xShiftStart");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "xShiftStart not specified! Default set to 6.\n";
		opt.xShiftStart = 6;
	}
	opt.xShiftEnd = OP.getDouble("xShiftEnd");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "xShiftEnd not specified! Default set to 11.\n";
		opt.xShiftEnd = 11;
	}

    // other parameters
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified, default to 0\n";
		opt.warningFlag = true;
		opt.seed = 0;
	}
    opt.rerunConf = OP.getConfFile();
	return opt;
}