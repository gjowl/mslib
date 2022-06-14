// This program reads the output file from AL2CO and maps
// the conservation values onto each position.  If the 
// protein is a homo-oligomer, the conservation values
// are mapped onto each chain.
#include <iostream>
#include <stdio.h>

#include "System.h"
#include "OptionParser.h"
#include "MslTools.h"
#include "Atom.h"
#include "PDBReader.h"

using namespace MSL;
using namespace std;
string programName = "addConservationBFactor";
string programDescription = "This program reads the output from al2co and alters the B-factors in a matching PDB file. For homo-oligomers, conservation is mapped to each chain.";
string programAuthor = "Samson Condon";
string programVersion = "0.0.1";
string programDate = "2015 June 16";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

struct Options {
	string pdbFile;
	string al2coFile;
	string output;
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
void usage();
void version();
void help(Options defaults);
Options parseOptions(int _argc, char * _argv[], Options defaults);
void printOptions(Options& _opt);

//TODO check if position in pdb file matches al2co sequence


int main (int argc, char *argv[]) {
	Options defaults;
	Options opt = parseOptions (argc, argv, defaults);
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

	printOptions(opt);

	//Read in and parse conservation file
	ifstream al2coFile (opt.al2coFile.c_str());
	string line;
	
	std::vector<string> residues;
	std::vector<double> conservation;

	if (al2coFile.is_open()) {
		while (getline (al2coFile, line)) {
			if (!isdigit(line[0])) {
				break;
			}
			//Example AL2CO output:
			//1     -      -1.000    *
			//2     P       0.860
			//3     A      -1.352
			std::vector<string> lineInfo = MslTools::tokenize(line);
			if (lineInfo[1] == "-") {
				continue;
			} else {
				residues.push_back(lineInfo[1]);
				conservation.push_back(MslTools::toDouble(lineInfo[2]));
			}
				

		}

	}
	al2coFile.close();
	

	//Read PDB
	System sys;
	sys.readPdb(opt.pdbFile);
	cout << sys << endl;

	//Loop over chains
	for (uint i = 0; i < sys.chainSize(); i++) {
		//loop over positions
		Chain &chaini = sys.getChain(i);
		for (uint j = 0; j < chaini.positionSize(); j++) {
			Position &posij = chaini.getPosition(j);
			//loop over atoms
			for (uint k = 0; k < posij.atomSize(); k++) {
				//set B-factor to conservation value
				Atom &atomijk = posij.getAtom(k);
				atomijk.setTempFactor(conservation[j]);
			}
		}
	}
	sys.writePdb(opt.output);



}


Options parseOptions(int _argc, char * _argv[], Options defaults) {

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
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	vector<string> allowed;

	opt.required.push_back("pdbFile");
	opt.required.push_back("al2coFile");
	opt.required.push_back("output");


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
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
//	if (OP.fail()) {
//		opt.help = OP.getBool("h");
//	}

	if (opt.help) {
		help(defaults);
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

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages = "pdb file not specified";
		opt.errorFlag = true;
	}
	opt.al2coFile = OP.getString("al2coFile");
	if (OP.fail()) {
		opt.errorMessages = "conservation file not specified";
		opt.errorFlag = true;
	}
	opt.output = OP.getString("output");
	if (OP.fail()) {
		opt.errorMessages = "output pdb file not specified";
		opt.errorFlag = true;
	}

	return opt;
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


void help(Options defaults) {
	cout << "This program runs as:" << endl;
	cout << " %addConservationBFactor --pdbFile <pdbfile> --al2coFile <conservation file>" << endl;
	cout << "--output <output PDB file>" << endl;
	cout << endl;
}
void printOptions(Options& _opt) {
	cout << "Program " << programName << " v." << programVersion << "," << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl; 
	cout << "pdbFile            " <<  _opt.pdbFile << endl;
	cout << "al2coFile          " <<  _opt.al2coFile << endl;
	cout << "output             " <<  _opt.output << endl;
}
