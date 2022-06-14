#include <iostream>

#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "AtomSelection.h"


using namespace MSL;
using namespace std;

string programName = "mutateAndRepackHelixDimer";
string programDescription = "This program reads in a pdb file, makes a mutation, loads rotamers, measures the monomer energy by seperating the helices and then measures the best dimer energy after repack and optionally outputs the repacked pdb file";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "0.0.0";
string programDate = "12 July 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

struct Options {

	string pdbFile;
	string rotlibFile;
	string rotLevel;
	vector<string> posToMutate;

	string outputFile;
	int startResNum;

	string chain1, chain2;


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


void repackSideChains(System & _sys, SelfPairManager & _spm, bool _greedy, int _greedyCycles) {

	_spm.calculateEnergies();

	if(!_greedy) {
		_spm.setMCOptions(1000.0, 0.5, 20000, MonteCarloManager::EXPONENTIAL, 2000, 100, 0.01);
		_spm.setRunDEE(false);
		_spm.setRunEnum(false);
		_spm.setRunSCMF(true);
		_spm.setRunSCMFBiasedMC(true);
		_spm.runOptimizer();
	} else {
		_spm.runGreedyOptimizer(_greedyCycles);
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
	trans.translate(chainA, interDistVect);

	// Repack side chains
	repackSideChains(_sys, _spm);

	_sys.setActiveRotamers(_spm.getMinStates()[0]);
	_sys.calcEnergy();
	_sys.printEnergySummary();

	//computeAllRelevantEnergies(_sys,"MONOMER");

	_sys.applySavedCoor("originState");

	return _spm.getMinBound()[0];

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

	PDBReader pRead(opt.pdbFile);
	pRead.open();
	if(!pRead.read()) {
		cerr << "Unable to open " << opt.pdbFile << endl;
		exit(0);
	}

	map<string,string> mutationMap;
	for (int i=0; i < opt.posToMutate.size(); i++) {
		vector<string> tokens;
		tokens = MslTools::tokenize(opt.posToMutate[i], ":");
		string positionId = tokens[0];
		string newRes = MslTools::toUpper(tokens[1]);
		mutationMap[positionId] = newRes;
	}

	AtomPointerVector atoms = pRead.getAtomPointers();
	string sequence = opt.chain1 + ": {" + MslTools::intToString(opt.startResNum) + "}";
	string oldChainId = opt.chain1;
	// assuming the atoms will be in order
	for(int i = 0; i < atoms.size(); i++) {

		// look at the CA atoms and get the position Id, if there is a match modify sequence accordingly
		if(atoms[i]->getName() != "CA") {
			continue;
		}

		string chainId = atoms[i]->getChainId();
		string positionId = atoms[i]->getPositionId();
		string resName = atoms[i]->getResidueName();

		if(oldChainId != chainId) {
			sequence += "\n" + chainId + ": {" + MslTools::intToString(opt.startResNum) + "}";
			oldChainId = chainId;
		}
		if(mutationMap.find(positionId) == mutationMap.end()) {
			sequence += resName + " " ; 
		} else {
			sequence += mutationMap[positionId] + " " ; 
		}
	}
	cout << sequence << endl;

	PolymerSequence ps;
	ps.setSequence(sequence);


	// Declare System
	System sys;

	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, "/library/charmmTopPar/top_all22_prot.inp" , "/library/charmmTopPar/par_all22_prot.inp");

	// Read in PDB File
	//CSB.setBuildTerm("CHARMM_ELEC",false);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW",true);

	if(!CSB.buildSystem(ps)) {
		cerr << "Unable to build system" << endl;
		exit(0);
	}

	sys.assignCoordinates(atoms,false);
	sys.buildAllAtoms();

	HydrogenBondBuilder hb(sys,"/data01/sabs/msl_working/mslib/trunk/toppar/scwrl4hb/par_hbond_2.txt");
	if(!hb.buildInteractions(50)) {
		cerr << "Unable to build hbonds " << endl;
		exit(0);
	}

	
	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotlibFile);
	sysRot.defineRotamerSamplingLevels();
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

	SelfPairManager spm;
	spm.setSystem(&sys);
	spm.setVerbose(false);

	// measure the monomer energy
	double monomerEnergy = computeMonomerEnergy(sys,spm);

	// measure the repacked energy
	//double repackEnergy = repackSideChains(sys,spm);

	//cout << "RESULT " <<  opt.pdbFile << " ";
	//for(int i = 0; i < opt.posToMutate.size(); i++) {
	//	cout << opt.posToMutate[i] << " ";
	//}
	//cout << monomerEnergy << " " << repackEnergy << " " << repackEnergy - monomerEnergy << endl;

	//sys.setActiveRotamers(spm.getMinStates()[0]);
	//sys.calcEnergy();
	//cout << sys.getEnergySummary() << endl;
	//computeAllRelevantEnergies(sys,"DIMER");

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	//fout << endl << "Monte Carlo repack from best position." << endl;

	xShift = savedXShift;

	if (opt.MCCycles > 0) {
		MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);

		MCMngr.setEner(bestEnergy);
		
		while(!MCMngr.getComplete()) {

			sys.applySavedCoor("savedBestState");
			helicalAxis.applySavedCoor("BestAxis");

			int moveToPreform = RNG1.getRandomInt(3);

			double deltaXShift = 0.0;
			double deltaZShift = 0.0;
			double deltaCrossingAngle = 0.0;
			double deltaAxialRotation = 0.0; 

			//======================================
			//====== Z Shift (Crossing Point) ======
			//======================================
			if (moveToPreform == 0) {
				deltaZShift = getStandardNormal(RNG1) * 0.1;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
			} else if (moveToPreform == 1) {
			//===========================
			//===== Axial Rotation ======
			//===========================
				deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
			} else if (moveToPreform == 2) {
			//==================================
			//====== Local Crossing Angle ======
			//==================================
				deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
			} else if (moveToPreform == 3) {
			//==============================================
			//====== X shift (Interhelical Distance) =======
			//==============================================
				deltaXShift = getStandardNormal(RNG1) * 0.1;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
			}

			// Run Optimization
			repackSideChains(sys, spm, opt.greedy, opt.greedyCycles);
	
			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				//fout << "state rejected   energy: " << currentEnergy << endl;
			}
			else {
				bestEnergy = currentEnergy;
				sys.saveAltCoor("savedBestState");
				helicalAxis.saveAltCoor("BestAxis");

				xShift = xShift + deltaXShift;
				crossingAngle = crossingAngle + deltaCrossingAngle;
				axialRotation = axialRotation + deltaAxialRotation;
				zShift = zShift +  deltaZShift;

				fout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				}
			}

		}
		fout << "Monte Carlo repack complete. " << endl << endl;
	}

	sys.applySavedCoor("savedBestState");
	/******************************************************************************
	 *               === END MONTE CARLO REPACK ===
	 ******************************************************************************/


	// Write out pdb file
	if(opt.outputFile != "") {
		cout << "writing " << opt.outputFile << endl;
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
	opt.required.push_back("posToMutate");

	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("chain1");
	opt.allowed.push_back("chain2");
	opt.allowed.push_back("outputFile");

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
	opt.posToMutate= OP.getStringVector("posToMutate");
	if (OP.fail()) {
		opt.errorMessages = "position to mutate not specified";
		opt.errorFlag = true;
	}

	opt.startResNum= OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages = "startResNum not specified using 1";
		opt.warningFlag = true;
		opt.startResNum = 1;
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

