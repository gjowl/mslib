#include <iostream>

#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "RandomNumberGenerator.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "AtomSelection.h"


using namespace MSL;
using namespace std;

string programName = "mutateAndRepackHelixDimer";
string programDescription = "This program builds a sequence onto an ideal symmetric alpha helix dimer, then changes the geometry to a particular configuration.  From there, it performs Monte Carlo backbone movements to produce an optimized structure";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "0.0.0";
string programDate = "12 July 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

struct Options {

	string pdbFile;
	string sequence;
	string rotlibFile;
	string rotLevel;
	
	double xShift;
	double zShift;
	double crossingAngle;
	double axialRotate;

	string outputFile;
	string logFile;
	int startResNum;

	string chain1, chain2;
	
	int MCCycles;
	int MCMaxRejects;
	bool greedy;
	int greedyCycles;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	string configfile; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run
};
/******************************************
 *  
 *  =======  PREDECLARATION OF THE FUNCTIONS  =======
 *
 ******************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);
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
	return "A" + ps + "\nB" + ps;
}
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
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
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized


	trans.rotate(_apvB, m);
	
}
// Just add 10 U(0,1) uniform random variables, offset by 0.5 to make mean = 0 and divide by variance = (10 * var(U(0,1))) 
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}
void computeAllRelevantEnergies(System& _sys, string _tag) {
	/******************************************************************************************
	* For each position, compute
	*	1. Interaction with self and helix backbone
	*	2. Interaction with other sidechains in the same helix
	*	3. Interaction with other helix sidechains
	*	4. Interaction with other helix backbone 
	*******************************************************************************************/
	AtomSelection sel(_sys.getAtomPointers());
	sel.select("bbA,chain A and name N+CA+C+O+HN");
	sel.select("bbB,chain B and name N+CA+C+O+HN");
	sel.select("scA,chain A and not name N+CA+C+O+HN");
	sel.select("scB,chain B and not name N+CA+C+O+HN");


	Chain& chainA = _sys.getChain(0);
	Chain& chainB = _sys.getChain(1);

	for(int i = 0; i < chainA.positionSize(); i++) {
		Position& pos = chainA.getPosition(i);
		string resNum = MslTools::intToString(pos.getResidueNumber());
		string selName = "posA_" + resNum;
		AtomPointerVector self = sel.select(selName + ", chain A and resi " + resNum + " and not name N+HN+CA+C+O");
		double selfE = _sys.calcEnergy(selName,selName);
		double ownChainSC = _sys.calcEnergy(selName,"scA");
		double ownChainBB = _sys.calcEnergy(selName,"bbA");
		double otherChainSC = _sys.calcEnergy(selName,"scB");
		double otherChainBB = _sys.calcEnergy(selName,"bbB");

		cout << _tag << " " << pos.getPositionId() << " " << selfE << " " << ownChainSC << " " << ownChainBB << " " << otherChainSC << " " << otherChainBB << endl;

	}

	for(int i = 0; i < chainB.positionSize(); i++) {
		Position& pos = chainB.getPosition(i);
		string resNum = MslTools::intToString(pos.getResidueNumber());
		string selName = "posB_" + resNum;
		AtomPointerVector self = sel.select(selName + ", chain B and resi " + resNum + " and not name N+HN+CA+C+O");
		double selfE = _sys.calcEnergy(selName,selName);
		double ownChainSC = _sys.calcEnergy(selName,"scB");
		double ownChainBB = _sys.calcEnergy(selName,"bbB");
		double otherChainSC = _sys.calcEnergy(selName,"scA");
		double otherChainBB = _sys.calcEnergy(selName,"bbA");

		cout << _tag << " " << pos.getPositionId() << " " << selfE << " " << ownChainSC << " " << ownChainBB << " " << otherChainSC << " " << otherChainBB << endl;

	}

}

double repackSideChains(System & _sys, SelfPairManager & _spm) {

	_spm.calculateEnergies();
	_spm.runGreedyOptimizer(5);
	return _spm.getMinBound()[0];
}


void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {

	 if (moveType == 0) {
		// Z Shift
		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by

		_trans.translate(_chainA, translateA);

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 1) {
		// Axial Rotation
		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else 	if (moveType == 2) {
		// Crossing Angle 
		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 3) {
		// XShift
		// Helix A interhelical distance
		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);
		_trans.translate(_axisA, translateA);

		// Helix B interhelical distance
		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else {
		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
	}
}




double computeMonomerEnergy(System & _sys, SelfPairManager & _spm) {

	_sys.saveAltCoor("originState");

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector &chainA = _sys.getChain(0).getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES 500 A APART ===
	 ******************************************************************************/
	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor(500.0, 0.0, 0.0);

	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized


	trans.translate(chainA, interDistVect);

	// Repack side chains
	repackSideChains(_sys, _spm);

	_sys.setActiveRotamers(_spm.getMinStates()[0]);
	double monomerEnergy = _spm.getMinBound()[0];
	_sys.calcEnergy();
	_sys.printEnergySummary();
	_sys.writePdb("monomer.pdb");
	//computeAllRelevantEnergies(_sys,"MONOMER");

	_sys.applySavedCoor("originState");

	return monomerEnergy;

}

bool createHelixFromPDB (System &_sys, PDBReader _pdbReader, string _pdbFile, CharmmSystemBuilder &_csb, PolymerSequence &_ps, string _sequence, int _startResNum) {
	string polymerSequence = convertToPolymerSequence (_sequence, _startResNum);
	_ps.setSequence(polymerSequence);
	_csb.buildSystem(_ps);

	_pdbReader.open(_pdbFile);
	if (!_pdbReader.read()) {
		cout << "PDB READER ERROR: CANNOT READ " << _pdbFile << endl;
		exit(0);
	}
	_pdbReader.close();
	_sys.assignCoordinates(_pdbReader.getAtomPointers(), false);
	_sys.buildAllAtoms();
	_sys.writePdb("createdHelix.pdb");

	return 1;

}



void transformCoiledCoil(System &_sys, double _zShift,double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis) {
	AtomPointerVector chainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector chainB = _sys.getChain("B").getAtomPointers();
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(chainA, _axialRotate, _origin, _zAxis);

	//====== Local Crossing Angle ======
	_trans.Xrotate(chainA, (_crossingAngle/2.0));
	//_trans.rotate(_axisA,  (_crossingAngle/2.0), _origin, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor(_xShift/2.0, 0.0, 0.0);
	_trans.translate(chainA, interDistVect);
	//_trans.translate(_axisA, interDistVect);

	c2Symmetry(chainA, chainB);
	//c2Symmetry(_axisA, _axisB);
}

int main(int argc, char *argv[]) {
	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are	
	 *  not always required                            
	 ******************************************************************************/
	Options defaults;                                  

	/******************************************************************************
	 *                             === OPTION PARSING ===
	 *                                                 
	 *  Parse the command line (and possibly input file)
	 *  options. It will also create an input file based
	 *  on the current options in the output directory
	 *  that can be used to re-run the program with the same
	 *  configuration
	 ******************************************************************************/
	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		//printErrors(opt);
		usage();
		exit(1);
	}
	CartesianPoint origin (0.0,0.0,0.0);
	CartesianPoint xAxis (1.0,0.0,0.0);
	CartesianPoint zAxis (0.0,0.0,1.0);

	// Declare System
	System sys;
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	PolymerSequence PS;
	PDBReader pdbRead;

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(0);

	ofstream fout;
	fout.open(opt.logFile.c_str());

	if(!fout.is_open()) {
		cerr << "Cannot open logfile: " << opt.logFile << endl;
		exit(0);
	}



	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, "/library/charmmTopPar/top_all22_prot.inp" , "/library/charmmTopPar/par_all22_prot.inp");

	// Read in template poly-gly PDB File
	//CSB.setBuildTerm("CHARMM_ELEC",false);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW",true);
	CSB.setBuildTerm("SCWRL4_HBOND");


	//Create Helix and duplicate to form ideal dimer
	createHelixFromPDB(sys, pdbRead, opt.pdbFile, CSB, PS, opt.sequence, opt.startResNum);
	
	//sys.duplicateChain("A");
	//CSB.updateNonBonded();



	HydrogenBondBuilder hb(sys,"/data01/sabs/msl_working/mslib/trunk/toppar/scwrl4hb/par_hbond_2.txt");
	if(!hb.buildInteractions(50)) {
		cerr << "Unable to build hbonds " << endl;
		exit(0);
	}

	
	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotlibFile);
	sysRot.defineRotamerSamplingLevels();
	fout << "Loading Rotamers..." << endl;
	/******************************************************************************
	 *                  === LOAD ROTAMERS & SET-UP SPM ===
	 ******************************************************************************/
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);

		// get the active identity in each position and load rotamers for that res
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!sysRot.loadRotamers(&pos, pos.getResidueName(), opt.rotLevel)) { 
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
	fout << "Finished loading rotamers" << endl;

	//Save superimposed helices as template state
	sys.saveAltCoor("Template");
	// Declare SelfPairManager and Set Seed
	SelfPairManager spm;
	spm.seed(RNG1.getSeed()); 
	spm.setSystem(&sys);
	spm.setVerbose(false);


	fout << sys << endl;
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");


	// measure the monomer energy
	double monomerEnergy = computeMonomerEnergy(sys,spm);


	//Transform Dimer to original parameters
	transformCoiledCoil(sys, opt.zShift, opt.crossingAngle, opt.axialRotate, opt.xShift, trans, origin, zAxis, xAxis);

	// measure the initial dimer energy
	double repackEnergy = repackSideChains(sys,spm);
	
	//fout << "RESULT " <<  opt.pdbFile << " ";
	//for(int i = 0; i < opt.posToMutate.size(); i++) {
	//	fout << opt.posToMutate[i] << " ";
	//}
	fout << sys.getEnergySummary() << endl;
	fout << monomerEnergy << " " << repackEnergy << " " << repackEnergy - monomerEnergy << endl;

	sys.setActiveRotamers(spm.getMinStates()[0]);
	//computeAllRelevantEnergies(sys,"DIMER");
	sys.saveAltCoor("savedBestState");
	sys.writePdb("test.pdb");

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	fout << endl << "Monte Carlo repack from best position." << endl;

	double bestEnergy = repackEnergy;
	double currentEnergy = bestEnergy;
	double xShift = opt.xShift;
	double zShift = opt.zShift;
	double crossingAngle = opt.crossingAngle;
	double axialRotate = opt.axialRotate;
	if (opt.MCCycles > 0) {
		MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);

		MCMngr.setEner(bestEnergy);
		
		while(!MCMngr.getComplete()) {

			sys.applySavedCoor("Template");

			//int moveToPreform = RNG1.getRandomInt(3);

			double deltaXShift = getStandardNormal(RNG1) * 0.1;
			double deltaZShift = getStandardNormal(RNG1) * 0.1;
			double deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
			double deltaAxialRotation = getStandardNormal(RNG1) * 1.0; 
			
			//======================================
			//====== Z Shift (Crossing Point) ======
			//======================================
			//if (moveToPreform == 0) {
			//	deltaZShift = getStandardNormal(RNG1) * 0.1;
			//} else if (moveToPreform == 1) {
			////===========================
			////===== Axial Rotation ======
			////===========================
			//	deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
			//} else if (moveToPreform == 2) {
			////==================================
			////====== Local Crossing Angle ======
			////==================================
			//	deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
			//} else if (moveToPreform == 3) {
			////==============================================
			////====== X shift (Interhelical Distance) =======
			////==============================================
			//	deltaXShift = getStandardNormal(RNG1) * 0.1;
			//}
			transformCoiledCoil(sys, zShift+deltaZShift, crossingAngle+deltaCrossingAngle, axialRotate+deltaAxialRotation, xShift+deltaXShift, trans, origin, zAxis, xAxis);
			// Run Optimization
			repackSideChains(sys, spm);
	
			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				fout << "state rejected   energy: " << currentEnergy-monomerEnergy << endl;
			}
			else {
				bestEnergy = currentEnergy;
				sys.saveAltCoor("savedBestState");

				xShift = xShift + deltaXShift;
				crossingAngle = crossingAngle + deltaCrossingAngle;
				axialRotate = axialRotate + deltaAxialRotation;
				zShift = zShift +  deltaZShift;

				fout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotate << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				}
			}

		}
		fout << "Monte Carlo repack complete. " << endl << endl;
	

	
	sys.applySavedCoor("savedBestState");
	double finalEnergy = sys.calcEnergy() - monomerEnergy;
	fout << sys.getEnergySummary() << endl;
	fout << "Final Energy: " << finalEnergy << endl;

	/******************************************************************************
	 *               === END MONTE CARLO REPACK ===
	 ******************************************************************************/


	// Write out pdb file
	if (opt.outputFile != "") {
		fout << "writing " << opt.outputFile << endl;
		if(!sys.writePdb(opt.outputFile)) {
			cerr << "Unable to write " << opt.outputFile << endl;
			exit(0);
		}
	}

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
	opt.required.push_back("rotlibFile");
	opt.required.push_back("rotLevel");
	opt.required.push_back("sequence");
	opt.required.push_back("xShift");
	opt.required.push_back("zShift");
	opt.required.push_back("axialRotate");
	opt.required.push_back("crossingAngle");

	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("chain1");
	opt.allowed.push_back("chain2");
	opt.allowed.push_back("outputFile");
	opt.allowed.push_back("logFile");
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");

	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	opt.allowed.push_back("configfile");


	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
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

	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			return opt;
		}
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
	opt.rotlibFile = OP.getString("rotlibFile");
	if (OP.fail()) {
		opt.errorMessages = "rotlib file not specified";
		opt.errorFlag = true;
	}
	opt.rotLevel = OP.getString("rotLevel");
	if (OP.fail()) {
		opt.errorMessages = "rotLevel not specified";
		opt.errorFlag = true;
	}
	opt.outputFile = OP.getString("outputFile");
	if (OP.fail()) {
		opt.warningMessages = "outputFile not specified - not writing repackedPdb";
		opt.warningFlag = true;
		opt.outputFile = "";
	}
	opt.logFile = OP.getString("logFile");
	if (OP.fail()) {
		opt.warningMessages = "log file not specified - writing to terminal";
		opt.warningFlag = true;
		opt.logFile = "";
	}

	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages = "polymer sequence not specified";
		opt.errorFlag = true;
	}
	opt.startResNum= OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages = "startResNum not specified using 1";
		opt.warningFlag = true;
		opt.startResNum = 1;
	}

	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.errorMessages = "xShift not specified";
		opt.errorFlag = true;
	}


	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.errorMessages = "zShift not specified";
		opt.errorFlag = true;
	}


	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngle not specified";
		opt.errorFlag = true;
	}


	opt.axialRotate = OP.getDouble("axialRotate");
	if (OP.fail()) {
		opt.errorMessages = "axialRotate not specified";
		opt.errorFlag = true;
	}
	
	opt.chain1= OP.getString("chain1");
	if (OP.fail()) {
		opt.warningMessages = "chain1 not specified using A";
		opt.warningFlag = true;
		opt.chain1 = "A";
	}

	opt.chain2= OP.getString("chain2");
	if (OP.fail()) {
		opt.warningMessages = "chain2 not specified using B";
		opt.warningFlag = true;
		opt.chain2 = "B";
	}
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "MCCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.MCCycles = 10;
	}
	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "MCMaxRejects not specified using 2\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 2;
	}
	opt.greedy = OP.getBool("greedyOptimizer");
	if (OP.fail()) {
		opt.warningMessages += "greedyOptimizer not specified using true\n";
		opt.warningFlag = true;
		opt.greedy = true;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
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
	cout << " % mutateAndRepackHelixDimer --pdbFile <name of pdb file> --rotlibFile <name of rotlib file> --rotLevel <rotLevel> --posToMutate <list of positions to mutate>  [--startResNum <default 1> --outputFile <destination of output pdb file>]" << endl;
	cout << endl;
}

