
#include <iostream>
#include <stdio.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "SysEnv.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "HelixGenerator.h"


using namespace MSL;
using namespace std;

static SysEnv SYSENV;

string programName = "genHomoUniverse";
string programDescription = "This program creates helix dimers of specified geometries and measures interhelical hydrogen bond energy. Finds the Dmin and also the Douts for each hydrogen bond at Dmin.";
string programAuthor = "Sabareesh Subramaniam and Samson Condon";
string programVersion = "1.0.6";
string programDate = "16 April 2015";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {
	string pdbFile;
	string sequence;
	int startResNum;
	string helicalAxisPdbFile;
	string output;

	double xShiftStart;
	double zShiftStart;
	double axialRotStart;
	double crossingAngleStart;
	//vector<int> rotCount;	
	double xShiftEnd;
	double zShiftEnd;
	double axialRotEnd;
	double crossingAngleEnd;
	//vector<int> rotCount;	
	double xShiftSteps;
	double zShiftSteps;
	double axialRotSteps;
	double crossingAngleSteps;
	double axZShift; // to move the z along with the axialRotation along the trapezoid
	double zAxShift; // to move the axialRotation along with the zShift along the trapezoid
	//vector<int> rotCount;	
	bool unitCellCoord; // true if using unit cell coordinates w' and Z' instead of Cartesian coordinates w and Z
	string topFile;
	string parFile;
	string hBondFile;

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

void printOptions(Options& _opt) {
	cout << "Program " << programName << " v." << programVersion << "," << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl; 
	cout << "pdbFile            " <<  _opt.pdbFile << endl;
	cout << "sequence           " << _opt.sequence << endl;
	cout << "helicalAxisPdbFile " <<  _opt.helicalAxisPdbFile << endl;
	cout << "output             " <<  _opt.output << endl;
	cout << "topFile            " <<  _opt.topFile << endl;
	cout << "parFile            " <<  _opt.parFile << endl;
	cout << "hBondFile          " <<  _opt.hBondFile << endl;
	cout << "xShiftStart        " <<  _opt.xShiftStart << endl;
	cout << "zShiftStart        " <<  _opt.zShiftStart << endl;
	cout << "axialRotStart      " <<  _opt.axialRotStart << endl;
	cout << "crossingAngleStart " <<  _opt.crossingAngleStart << endl;
	cout << "xShiftEnd          " <<  _opt.xShiftEnd << endl;
	cout << "zShiftEnd          " <<  _opt.zShiftEnd << endl;
	cout << "axialRotEnd        " <<  _opt.axialRotEnd << endl;
	cout << "crossingAngleEnd   " <<  _opt.crossingAngleEnd << endl;
	cout << "xShiftSteps        " <<  _opt.xShiftSteps << endl;
	cout << "zShiftSteps        " <<  _opt.zShiftSteps << endl;
	cout << "axialRotSteps      " <<  _opt.axialRotSteps << endl;
	cout << "crossingAngleSteps " <<  _opt.crossingAngleSteps << endl;
	cout << "axZShift           " <<  _opt.axZShift << endl; 
	cout << "zAxShift           " <<  _opt.zAxShift << endl;
	cout << "unitCellCoord      " << _opt.unitCellCoord << endl;
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);

void transformCoiledCoil(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShift,double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, AtomPointerVector& _axisA, AtomPointerVector& _axisB);

string convertToPolymerSequence(string _seq, int _startResNum);
unsigned int getNumberOfInterHelicalHbonds(vector<Interaction*>& hbondInteractions);
map<string,double> getInterHelicalHbondInfo(vector<Interaction*>& hbondInteractions);
string printInfo(double axialRotate,double crossingAngle,double zShift,double relAxRot,double relZShift,vector<double> xShifts,vector<double> energies,vector<map<string,double> > hBonds);
bool createHelices (System &_sys, CharmmSystemBuilder &_csb, Transforms &_tr, HelixGenerator &_hg, PolymerSequence &_ps, string _sequence, int _startResNum);

/*========================================== BEGIN MAIN ==========================================*/
int main(int argc, char *argv[]) {

	time(&startTime);	
	
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);
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

	/******************************************************************************
	 *                     === SYSTEM SETUP ===
	 ******************************************************************************/
	
	// Declare System
	System sys;
	HelixGenerator hg;
	PolymerSequence PS;

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile);
	
	//createHelices (sys, CSB, trans, hg, PS, opt.sequence, opt.startResNum);
	//cout << sys << endl;
	//sys.writePdb("filler.pdb");
	
	// Read in PDB File
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");
	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
		cerr << "Unable to build system" << endl;
		exit(0);
	}


	// Build system
	sys.buildAllAtoms();
	cout << sys << endl;

	// Add hydrogen bonds
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	//HydrogenBondBuilder hb(sys, "/data00/bkmueller/dataFiles/hbondlist_nonCanon_adj.txt");
	hb.buildInteractions(30);
	cout << sys << endl;
	sys.calcEnergy();
	sys.printEnergySummary();



	// Redirect Output
	string filename = opt.output;
	if(!freopen (filename.c_str(), "w", stdout)) {
		cerr << "ERROR_freopen: Did not open output file " << opt.output << endl;
		exit(0);
	}
	// Set up APVs
	AtomPointerVector &chainA = sys.getChain("A").getAtomPointers();
	AtomPointerVector &chainB = sys.getChain("B").getAtomPointers();



	/******************************************************************************
	 *                     === MONTE CARLO SET UP ===
	 ******************************************************************************/


	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxisPdbFile);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/

	AtomSelection sel(sys.getAtomPointers());
	sel.select("chainA,chain A");
	sel.select("chainB,chain B");

	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	
	sys.saveCoor("initialState");
	helicalAxis.saveCoor("initialState");

	EnergySet* pESet = sys.getEnergySet();
	vector<Interaction*> hbondInteractions = (*(pESet->getEnergyTerms()))["SCWRL4_HBOND"];
	sys.saveEnergySubset("interHelical","chainA", "chainB");

	//cout << "Number of hbond interactions "<< hbondInteractions.size() << endl;
	// Reference points to set up Helical starting postions
	CartesianPoint origin(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);


	vector<double> xShifts;
	for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
		xShifts.push_back(xShift);
	}


	for(double axialRotate = opt.axialRotStart; axialRotate < opt.axialRotEnd; axialRotate += opt.axialRotSteps) {
		double zShiftStart = opt.zShiftStart + axialRotate * opt.axZShift; // opt.axZShift = -0.015, converts relative Z' to absolute Z

		//cout << "Zrel: " << opt.zShiftStart << endl;
		//cout << "Arel: " << axialRotate << endl;
		//cout << "-0.015: " << opt.axZShift << endl;
		//cout << "Arel * -0.015: " << axialRotate * opt.axZShift << endl;
		//cout <<  "Zabs: " << opt.zShiftStart + axialRotate * opt.axZShift << endl;

		double zShiftEnd = opt.zShiftEnd + axialRotate * opt.axZShift; // opt.axZShift = -0.015
		double relZ = opt.zShiftStart;
		for( double zShift = zShiftStart; zShift < zShiftEnd; zShift += opt.zShiftSteps ,relZ += opt.zShiftSteps) {
			double shiftedAxialRotate = axialRotate + (zShift - axialRotate * opt.axZShift)  * opt.zAxShift; //convert relative w' to absolute w
			//cout << "Zabs: " << zShift << endl;
			//cout << "(Zabs - Arel * -0.015) * -6.6666667: " << (zShift - axialRotate * opt.axZShift)  * opt.zAxShift << endl;
			//cout << "Aabs: " << shiftedAxialRotate <<  endl;
			//continue;
			
			//double zTranslation = zShift;
			//double axialRotation = axialRotate;
			//if(opt.unitCellCoord) {
			//	//convert unit cell coordinates to absolute coordinates for the transformation
			//	zTranslation = zShift - 0.015 * axialRotate;
			//	axialRotation = axialRotate - 20.0 * zShift / 3.0;
			//}

			for(double crossingAngle = opt.crossingAngleStart; crossingAngle < opt.crossingAngleEnd; crossingAngle += opt.crossingAngleSteps ) {

				vector<map<string,double> >  hBonds;
				vector<double> energies;
				sys.applySavedCoor("initialState");
				//helicalAxis.applySavedCoor("initialState");
				//transformCoiledCoil(chainA,chainB,zShift,crossingAngle,shiftedAxialRotate,opt.xShiftStart - opt.xShiftSteps ,trans,origin,zAxis,xAxis,axisA,axisB);
				transformCoiledCoil(chainA,chainB,zShift,crossingAngle,shiftedAxialRotate,opt.xShiftStart - opt.xShiftSteps ,trans,origin,zAxis,xAxis,axisA,axisB);
				CartesianPoint xShiftASize (opt.xShiftSteps/-2.0,0,0);
				CartesianPoint xShiftBSize (opt.xShiftSteps/2.0,0,0);
				
				for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
					// do the xShiftSize move on both chains
					trans.translate(chainA,xShiftASize);
					trans.translate(chainB,xShiftBSize);

					double thisEnergy = sys.calcEnergyOfSubset("interHelical");
					//unsigned int numHbonds = getNumberOfInterHelicalHbonds(hbondInteractions);
					hBonds.push_back(getInterHelicalHbondInfo(hbondInteractions));
					energies.push_back(thisEnergy);
					//cout << "axialRotation: " << axialRotate  << " crossingAngle: " << crossingAngle << " zShift: " << zShift <<  " xShift: " << xShift << " numCA_HBonds: " << numHbonds << endl;
				}
					string name = printInfo(axialRotate,crossingAngle,relZ,shiftedAxialRotate,zShift,xShifts,energies,hBonds);
					
				/*
				if(!sys.writePdb(name) ) {
					cerr << "Unable to write " << name << endl;
					exit(0);
				}
				*/
				//char name[100];
				//sprintf(name,"model_%02.0f_%02.0f_%+06.3f",shiftedAxialRotate,crossingAngle,zShift);
				//cout << name << endl;

			}
		}
	}



	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

}
/*========================================== END MAIN ==========================================*/


void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->setCoor(_apvA[i]->getCoor());	
	}

	// Rotation matrix for 180 degrees
	//      ____________
	//     |            |
	//     | -1   0   0 |
	//     |            |
	//     |  0  -1   0 |
	//     |            |
	//     |  0   0  -1 |
	//     |____________|
	
	Matrix m(3,3,0.0);
	m[0][0] = -1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = -1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;

	// Rotate chain B around Z axis
	Transforms trans; 
	trans.rotate(_apvB, m);
}

void transformCoiledCoil(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShift,double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, AtomPointerVector& _axisA, AtomPointerVector& _axisB) {
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotate, _origin, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _origin, _xAxis);
	//_trans.rotate(_axisA,  (_crossingAngle/2.0), _origin, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((_xShift/2.0) * -1.0, 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	//_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	//c2Symmetry(_axisA, _axisB);
}

unsigned int getNumberOfInterHelicalHbonds(vector<Interaction*>& hbondInteractions) {
	unsigned int numHbonds = 0;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			numHbonds++;
		}
	}
	return numHbonds;
}

map<string,double> getInterHelicalHbondInfo(vector<Interaction*>& hbondInteractions) {
	map<string,double> info;

	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			info[atoms[0]->getAtomId() + ":" + atoms[2]->getAtomId()] = e;
		}
	}
	return info;
}

string printInfo(double axialRotate,double crossingAngle,double zShift,double relAxRot,double relZShift,vector<double> xShifts,vector<double> energies,vector<map<string,double> > hBonds) {
	int minIdx = 0; // into xShifts, energies,hBonds
	for(int i = 1; i < energies.size(); i++) {
		if(energies[i] < energies[minIdx]) {
			minIdx = i;
		}
	}

	// get the hbonds at minIdx and see where they get broken
	map<string,double> hmin = hBonds[minIdx];
	map<string,double> douts;
	double energy = 0;

	for(map<string,double>::iterator it = hmin.begin(); it != hmin.end(); it++) {
		string name = it->first;
		douts[name] = xShifts[minIdx];
		energy += it->second;

		for(int i = minIdx+1; i < hBonds.size(); i++) {
			if(hBonds[i].find(name) == hBonds[i].end()) {
				douts[name] = xShifts[i];
				break;
			}
		}
	}
	char name[100];
	sprintf(name,"%03.3f\t" "%03.1f\t" "%03.3f\t" "%03.3f\t" "%03.3f\t" "%03.1f", axialRotate, crossingAngle, zShift, relAxRot, relZShift, xShifts[minIdx]);
	cout << name << "\t" << hBonds[minIdx].size() ;
	cout << "\t" << energies[minIdx] << "\t" << energy;

	for(int i = minIdx + 1; i < xShifts.size(); i++) {
		char str[100];
		sprintf(str,"%02.1f",xShifts[i]);
		for(map<string,double>::iterator it = hmin.begin(); it != hmin.end(); it++) {
			if(douts[it->first] == xShifts[i]) {
				char bondE[100];
				sprintf(bondE,"%+010.6f",it->second);
				cout << "\t" << it->first << "\t" << str << "\t" << bondE;
			}
		}
	}
	cout << endl;
	return string(name);

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
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("startResNum");
	opt.required.push_back("output");

	opt.required.push_back("xShiftStart");
	opt.required.push_back("xShiftEnd");
	opt.required.push_back("xShiftSteps");

	opt.required.push_back("zShiftStart");
	opt.required.push_back("zShiftEnd");
	opt.required.push_back("zShiftSteps");

	opt.required.push_back("axialRotStart");
	opt.required.push_back("axialRotEnd");
	opt.required.push_back("axialRotSteps");

	opt.required.push_back("crossingAngleStart");
	opt.required.push_back("crossingAngleEnd");
	opt.required.push_back("crossingAngleSteps");
	opt.required.push_back("axZShift");
	opt.required.push_back("zAxShift");

	//opt.required.push_back("rotCount");
	//
	opt.required.push_back("helicalAxisPdbFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("hBondFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("unitCellCoord");

	//opt.equivalent.push_back(vector<string>());
	//opt.equivalent.back().push_back("v");
	//opt.equivalent.back().push_back("version");
	//opt.equivalent.push_back(vector<string>());
	//opt.equivalent.back().push_back("h");
	//opt.equivalent.back().push_back("help");

	//opt.defaultArgs.push_back("configfile");


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

	/*****************************************
	 *  OUTPUT DIR AND FILES
	 *****************************************/


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages = "pdb file not specified";
		opt.errorFlag = true;
	}
	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages = "pdb file not specified";
		opt.errorFlag = true;
	}
	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.startResNum = 1;
		opt.warningMessages += "starting residue number not specified; defaulting to 1";
		opt.warningFlag = true;
	}
	opt.helicalAxisPdbFile = OP.getString("helicalAxisPdbFile");
	if (OP.fail()) {
		opt.errorMessages = "helicalAxisPdbFile file not specified";
		opt.errorFlag = true;
	}
	opt.output = OP.getString("output");
	if (OP.fail()) {
		opt.errorMessages = "output not specified";
		opt.errorFlag = true;
	}
	opt.xShiftStart = OP.getDouble("xShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "xShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.xShiftEnd = OP.getDouble("xShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "xShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.xShiftSteps = OP.getDouble("xShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "xShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.zShiftStart = OP.getDouble("zShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "zShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.zShiftEnd = OP.getDouble("zShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "zShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.zShiftSteps = OP.getDouble("zShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "zShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.axialRotStart = OP.getDouble("axialRotStart");
	if (OP.fail()) {
		opt.errorMessages = "axialRotStart not specified";
		opt.errorFlag = true;
	}
	opt.axialRotEnd = OP.getDouble("axialRotEnd");
	if (OP.fail()) {
		opt.errorMessages = "axialRotEnd not specified";
		opt.errorFlag = true;
	}
	opt.axialRotSteps = OP.getDouble("axialRotSteps");
	if (OP.fail()) {
		opt.errorMessages = "axialRotSteps not specified";
		opt.errorFlag = true;
	}

	opt.crossingAngleStart = OP.getDouble("crossingAngleStart");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleStart not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleEnd = OP.getDouble("crossingAngleEnd");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleEnd not specified";
		opt.errorFlag = true;
	}
	opt.axZShift = OP.getDouble("axZShift");
	if (OP.fail()) {
		opt.errorMessages = "axZShift not specified";
		opt.errorFlag = true;
	}

	opt.zAxShift = OP.getDouble("zAxShift");
	if (OP.fail()) {
		opt.errorMessages = "zAxShift not specified";
		opt.errorFlag = true;
	}

	opt.crossingAngleSteps = OP.getDouble("crossingAngleSteps");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleSteps not specified";
		opt.errorFlag = true;
	}

	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(envVar)) {
			opt.topFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}

	}

	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.parFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hBondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hBondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hBondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.unitCellCoord = OP.getBool("unitCellCoord");
	if (OP.fail()) {
		opt.unitCellCoord = false;
		opt.warningMessages += "Unit Cell coordinates not specified using --unitCellCoord. Assuming Cartesian coordinates for axRot and Z\n";
		opt.warningFlag = true;
	}

	// return the Options structure
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
	cout << " % genHomoUniverse --pdbFile <pdbfile> --helicalAxisPdbFile <filename> --output <filename>" << endl;
	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
	cout << endl;
}


string convertToPolymerSequence(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
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
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}

//function to create two idealized alpha helices of defined sequence that are in an orientation compatible for universe analysis
bool createHelices (System &_sys, CharmmSystemBuilder &_csb, Transforms &_tr, HelixGenerator &_hg, PolymerSequence &_ps, string _sequence, int _startResNum) {
	string polymerSequence = convertToPolymerSequence (_sequence, _startResNum);

	int chainMidpoint = polymerSequence.length() / 2 + 1;	//The alpha carbon at the chain midpoint will be set to lie on the X-axis
	string midpointCAID = "A " + MslTools::intToString(chainMidpoint) + " CA";

	_ps.setSequence(polymerSequence);
	_csb.buildSystem(_ps);

	AtomPointerVector *helix = _hg.generateHelix(_sequence.length());
	//TODO aign helix with Z-axis
	CartesianPoint O(0.0,0.0,0.0);
	CartesianPoint X(1.0,0.0,0.0);
	CartesianPoint Y(0.0,1.0,0.0);
	CartesianPoint Z(0.0,0.0,1.0);
	
	CartesianPoint geoCenter = helix->getGeometricCenter();

	_tr.orient(*helix, geoCenter, Z, O, X);

	//if (!align(&Helix, reference, Z) ) {
	//
	//}
	//TODO put the middle alpha carbon on the X-axis (2.271, 0.000, 0.000)

	_sys.assignCoordinates(*helix, false);
	_sys.buildAllAtoms();

	CartesianPoint midCA = _sys.getAtom(midpointCAID).getCoor();
	CartesianPoint distance = midCA - O;
	_tr.translate (_sys.getAtomPointers(), distance);

	_sys.duplicateChain("A", "B");

	return 1;
}
