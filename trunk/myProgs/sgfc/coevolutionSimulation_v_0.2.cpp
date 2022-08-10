// This program will simulate evolution within a membrane helix dimer in 
// order to make hypotheses about how helices coevolve to generate these
// structures.  We can use these predictions to benchmark coevolution 
// algorithms and find coevolution hotspots.
// 
// This program will read in a PDB file of a membrane helix dimer and load
// rotamers for hydrophobic identities.  It will then undergo cycles of 
// "evolution" where a residue is mutated and the resulting structure 
// is optimized.  Using a Monte Carlo procedure, residues will be mutated
// and the PDB files printed out until a "decent" energy is reached.
//
// Output will be a list of sequences that the program sampled along with
// the energy of the puative dimer structure relative to the wt.


#include "System.h"
#include "SelfPairManager.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "OptionParser.h"
#include "SystemRotamerLoader.h"
#include "SysEnv.h"
#include "RandomNumberGenerator.h"
#include "PDBTopologyBuilder.h"
#include "MslTools.h"
#include "Transforms.h"
#include "AtomSelection.h"
#include "BaselineEnergyBuilder.h"
#include "MonteCarloManager.h"


using namespace std;
using namespace MSL;


string programName = "coevolutionSimulation";
string programDescription = "This program mutates positions in a PDB structure for a number of cycles, selecting against unfavorable mutations";
string programAuthor = "Samson Gerald Funakoshi Condon";
string programVersion = "0.1.2";
string programDate = "26 March 2014";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t start_time, finish_time;

static SysEnv SYSENV;

//Options

struct Options{
	// Required
	string pdbFile;


	// Allowed
	vector<int> varPos;	//variable positions to mutate
	vector<string> varID;	//identities to mutate to

	string samplingLevel;	//rotamers to load for each identity


	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;

	int seed;		//Seed for the random number generator.  Time-based seed used if one is not supplied

	string configFile;	//Config file where all the options are listed, one per line.


	// Dependent on the output directory
	string outputDir;	//place to write pdb files for solutions and text files
	string logFile; 	
	string outputFile;	//writes sequence output and energies to a text file
	string fastaFile;	//writes the sequences in fasta format to a text file for alignment
	string rerunConfigFile;	//ConfigFile that would rerun the program with the same options as before


	//MANAGEMENT VARIABLES
	string pwd;	// the present working directory obtained with a getenv
	string host;	// the host name obtained with a getenv
	bool version;	// ask for the program version
	bool help;	// ask for the program help

	bool errorFlag;		// true if there are errors
	bool warningFlag;	// true if there are warnings
	string errorMessages;	// error messages
	string warningMessages;	// warning messages

	vector<string> allowed;			// list of allowed options
	vector<string> required;		// list of required options
	vector<string> dependent;		// list of options that depend on the existence of the first option in the list
	vector< vector <string> > dependentSet;	// entire set of dependent options

	string rerunConfig;	// data for a configuration file that would rerun the job as the current run


	vector<string> disallowed;		// disallowed options that were given
	vector<string> missing;			// required options that were not given
	vector<string> ambiguous;		// required options that were not given
	vector<string> defaultArgs;		// the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent;	// this links short options to long ones (for example -x can be given for --extended)

	string OPerrors;	// the errors from the option parser
};
Options parseOptions(int _argc, char * _argv[], Options defaults);

std::vector< std::vector<bool> > buildResidueMask(std::vector<unsigned int> _residueState, std::vector< std::vector<unsigned int> > _residueRotamers);

int main (int argc, char *argv[]) {
	// store the start time
	time(&start_time);

	Options defaults;
	Options opt = parseOptions(argc, argv, defaults);
	
	if (opt.outputDir != "" && opt.rerunConfigFile != "") {
		string rerunConfig = opt.outputDir + opt.rerunConfigFile;
		ofstream configOut;
		configOut.open(rerunConfig.c_str());
		configOut << opt.rerunConfig << endl;
	}


	Transforms tr;
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);
	MonteCarloManager MC;

	System sys;
	CharmmSystemBuilder CSB;
	HydrogenBondBuilder SCWRL;

	PDBTopologyBuilder PTB (sys, opt.topFile);
	SystemRotamerLoader SRL (sys, opt.rotLibFile);

	CSB.setSystem(sys);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW", true);

	SCWRL.setSystem(sys);
	
	BaselineEnergyBuilder BEB(sys, "/exports/home/scondon/mslib/trunk/myProgs/sgfc/helixBaseline/helix_VdW_SCWRL.bas");

	std::vector <unsigned int> residueState;
	std::vector < string > outputList;
	std::vector < std::vector<bool> > resMask; 

	//==============================================================
	// Set up system and build interactions
	// =============================================================
//
//	if (!sys.readPdb(opt.pdbFile)) {
//		cout << "ERROR_SYS: Did not read PDB File..." << endl;
//		exit(5);
//	}



	if(!CSB.readTopology(opt.topFile)) {
		cout << "ERROR_CSB: Did not read Topology File.." << endl;
		exit (1);
	}

	if(!CSB.readParameters(opt.parFile)) {
		cout << "ERROR_CSB: Did not read Parameter File..." << endl;
		exit (1);
	}
	
	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
		cout << "ERROR_CSB: Did not build system from PDB..." << endl;
		exit (1);
	}

	cout << sys << endl;
	if(!SCWRL.readParameters(opt.hBondFile)) {
		cout << "ERROR_SCWRL: Did not read H-Bond File..." << endl;
		exit(2);
	}


	if(!SCWRL.buildInteractions()) {
		cout << "ERROR_SCWRL: Did not build H-Bond Interactions..." << endl;
		exit(2);
	}
	
	if(!BEB.buildInteractions()) {
		cout << "ERROR_BEB: Did not build baseline interactions..." << endl;
		exit(3);
	}

	cout << sys << endl;
	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;
	cout << opt.topFile << endl;
	cout << opt.parFile << endl;
	cout << opt.hBondFile << endl;


	
	//==============================================================
	// Set the linked positions--assume homomer for now
	//==============================================================
	Chain masterChain = sys.getChain(0);
	for (uint i = 0; i < masterChain.positionSize(); i++) {
		std::vector<string> linked(sys.chainSize());
		for (uint j = 0; j < sys.chainSize(); j++) {
			string posID = sys.getChain(j).getPosition(i).getPositionId();
			linked[j] =posID;
		}
		sys.setLinkedPositions(linked);
	}


	//==============================================================
	// Add identities and load rotamers
	//==============================================================
	for (uint i = 0; i < sys.positionSize(); i++) {
		cout << "\rLoading Identities and Rotamers... " << i+1 << "/" << sys.positionSize() << flush;
		Position &pos = sys.getPosition(i);
		Residue &res = pos.getCurrentIdentity();
		string wtIdentity = res.getResidueName();
		// add rotamers for the current identity
		if (wtIdentity != "ALA" && wtIdentity != "GLY") {
			SRL.loadRotamers(&pos, wtIdentity, opt.samplingLevel);
		}

		// add rotamers for the variable identities
		for (uint j = 0; j < opt.varID.size(); j++) {
			string identity = opt.varID[j];
			if (wtIdentity != identity) {
				CSB.addIdentity(pos, identity, "N CA C O H");
			}
			if (identity != "ALA" && identity != "GLY" && wtIdentity != identity) {
				SRL.loadRotamers(&pos, identity, opt.samplingLevel);
			}
		}


	}
	cout << endl;
	cout << "Rotamers loaded!" << endl;
	cout << "Setting up interactions.  Please wait."  << endl;

	CSB.updateNonBonded(9.0,10.0,11.0);
	SCWRL.update();


	//==============================================================
	// Set up the residue mask
	//==============================================================
	
	//initialize the residue state at 0 (the original identity/rotamer) for each position
	sys.updateVariablePositions();
	residueState = std::vector<unsigned int>(sys.getMasterPositions().size(),0);
	
	cout << "Residue state: " << endl;
	cout << "==================================" << endl;
	for (unsigned int i = 0; i < residueState.size(); i++) {
		cout << residueState[i] << " ";
	}
	cout << endl;


	//for each position, get the number identities and the number of rotamers for each identity
	std::vector < std::vector<unsigned int> > resRots (residueState.size());
	//for each master position
	for (unsigned int i = 0; i < resRots.size(); i++) {
		Position &pos = sys.getPosition(sys.getMasterPositions()[i]);
		resRots[i] = std::vector<unsigned int> (pos.identitySize());
		//cout << pos.getPositionId() << " has " << resRots[i].size() << " identities with total rotamers: "; 

		//for each identity
		for (unsigned int j = 0; j < pos.identitySize(); j++) {
			resRots[i][j] = pos.getTotalNumberOfRotamers(j);
			//cout << resRots[i][j] << " ";
		}
		//cout << endl;
	}
	
	cout << "resRots: " << endl;
	cout << "=================================" << endl;
	for (unsigned int i = 0; i < resRots.size(); i++) {
		for (unsigned int j = 0; j < resRots[i].size(); j++) {
			cout << resRots[i][j] << " ";
		}
		cout << endl;
	}

	//build the residue mask based on the current residue state and the number of rotamers
	resMask = buildResidueMask(residueState, resRots);




	SelfPairManager SPM(&sys);
	SPM.setOnTheFly(true);
	SPM.saveEnergiesByTerm(true);
	SPM.saveInteractionCounts(true);
	SPM.calculateEnergies();


	//Get the energy for the wild type
	cout << "Wild type energy" << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	resMask = buildResidueMask(residueState, resRots);
	SPM.runGreedyOptimizer(3,resMask);
	
	sys.setActiveRotamers(SPM.getMinStates()[0]);
	cout << "Dimer Energy" << endl;
	cout << SPM.getSummary(SPM.getMinStates()[0]) << endl;


	double wtdimerE = SPM.getStateEnergy(SPM.getMinStates()[0]);
/*
	tr.translate(sys.getChain(0).getAtomPointers(), CartesianPoint(500.0, 0.0, 0.0));
	SPM.calculateEnergies();
	SPM.runGreedyOptimizer(3,resMask);

	sys.setActiveRotamers(SPM.getMinStates()[0]);
	cout << "Monomer Energy" << endl;
	
	cout << SPM.getSummary(SPM.getMinStates()[0]) << endl;
	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;


	double wtmonomerE = SPM.getStateEnergy(SPM.getMinStates()[0]);
	string wtSequence;
		for (uint j = 0; j < sys.positionSize(); j++) {
			string threeLetterRes = sys.getPosition(j).getCurrentIdentity().getResidueName();
			string oneLetterRes = MslTools::getOneLetterCode(threeLetterRes);
			wtSequence = wtSequence + oneLetterRes;
		}
	cout << "==================================================" << endl;
	double wteDiff = wtdimerE - wtmonomerE;

	string wtoutput =wtSequence + "," + MslTools::doubleToString(wtmonomerE) + "," + MslTools::doubleToString(wtdimerE) + "," + MslTools::doubleToString(wteDiff);
	cout << wtoutput << endl;
	outputList.push_back(wtoutput);

	sys.applySavedCoor("dimer");
	sys.writePdb("/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/coevolutionDimerReturn.pdb");

	sys.clearSavedCoor("dimer");
*/


	//==============================================================
	// Begin mutating residues
	//==============================================================
	double currentE = wtdimerE;
	unsigned int mutations = 0;
	while (outputList.size() < 100) {
		double dimerE = 0;
//		double monomerE = 0;
//		double eDiff = 0;
		cout << "=======================" << mutations << "=======================" << endl;
		//Mutate a random position to a random available residue

		unsigned int mutatePos = RNG.getRandomInt(residueState.size()-1);
		unsigned int wtRes = residueState[mutatePos];	// we will revert back to this if the mutation is unfavorable
		unsigned int mutateRes = RNG.getRandomInt(resRots[mutatePos].size()-1);

		//If we happen to make a "silent" mutation, try again until we can switch residues
		while (mutateRes == wtRes) {
			mutateRes = RNG.getRandomInt(resRots[mutatePos].size()-1);
		}

		string posID = sys.getPosition(mutatePos).getPositionId();
		string wtResName = sys.getPosition(mutatePos).getResidue(wtRes).getResidueName();
		string mutateResName = sys.getPosition(mutatePos).getResidue(mutateRes).getResidueName();

		residueState[mutatePos] = mutateRes;
		resMask = buildResidueMask(residueState, resRots);
		


		//calculate the dimer energy
//		SPM.calculateEnergies();
		SPM.runGreedyOptimizer(3, resMask);
		dimerE = SPM.getStateEnergy(SPM.getMinStates()[0]);

		
		
		// move the helices away from each other and get the monomer energy
		// NO LONGER NEEDED, CALCULATING BASELINE ENERGIES TO APPROXIMATE MONOMER

//		sys.saveAltCoor("dimer");
//		tr.translate(sys.getChain(0).getAtomPointers(), CartesianPoint(500.0, 0.0, 0.0));
//		SPM.calculateEnergies();
//		SPM.runGreedyOptimizer(3,resMask);
//		monomerE = SPM.getStateEnergy(SPM.getMinStates()[0]);
//		eDiff = dimerE - monomerE;
		
/*
		bool accept = false;
		// Check the energy of the new dimer
		// We will assume that the "wt" dimer is at the baseline of natural selection--more stable
		// variants won't be selected for
		//
		// Also, if a mutation is destabilizing, we should accept any future mutation that brings us closer to 
		// a stable dimer

		if (dimerE <= wtdimerE || dimerE <= currentE) {
			accept = true;
		} else {
			double prob = RNG.getRandomDouble();
			if (prob < 0.1) {
				accept = true;
			}
		}
		//Clash check--ignore results above a certain energy

		if (dimerE >= 1000) {
			accept = false;
		}
*/
		bool accept = MC.accept(dimerE);

		//Now that we have decided to either accept or reject the mutation, we can proceed and
		//actually accept or reject it
		if (accept == true) {
			currentE = dimerE;
			sys.setActiveRotamers(SPM.getMinStates()[0]);
			if(opt.outputDir != "") {
				sys.writePdb(opt.outputDir + "evolutionAccept_" + MslTools::intToString(mutations) + ".pdb");
			}
			//collect output: sequence and energies
			cout << "Mutation accepted!" << endl;
			string sequence;
			for (uint j = 0; j < sys.positionSize(); j++) {
				string threeLetterRes = sys.getPosition(j).getCurrentIdentity().getResidueName();
				string oneLetterRes = MslTools::getOneLetterCode(threeLetterRes);
				//sequence = sequence + "," + oneLetterRes;
				sequence = sequence + oneLetterRes;
			}
			//string output = sequence + "," + MslTools::doubleToString(monomerE) + "," + MslTools::doubleToString(dimerE) + "," + MslTools::doubleToString(eDiff);
			string output = MslTools::intToString(mutations) + "," + posID + "," + wtResName + "," + mutateResName + "," + sequence + "," + MslTools::doubleToString(dimerE);
			outputList.push_back(output);
	
			//revert back to dimer state
//			sys.applySavedCoor("dimer");
//			sys.clearSavedCoor("dimer");


		} else {
			cout << "Mutation rejected..." << endl;
			residueState[mutatePos] = wtRes;
		}

		resMask.clear();
		mutations++;

	}

	if (opt.outputDir != ""  && opt.outputFile != "") {
		string outputFile = opt.outputDir + opt.outputFile;
		ofstream outStream;
		outStream.open(outputFile.c_str());
		for (uint i = 0; i < outputList.size(); i++) {
			outStream << outputList[i] << endl;
		}
		outStream.close();
	} else {
		for (uint i = 0; i < outputList.size(); i++) {
			cout << outputList[i] << endl;
		}
	}
	return 0;
}



Options parseOptions(int _argc, char * _argv[], Options defaults) {
	

	Options opt;



	vector<string> required;
	vector<string> allowed;

	//========================================================
	//
	//================Required Options========================
	//
	//========================================================

	opt.required.push_back("pdbFile");


	//========================================================
	//
	//=================Allowed Options========================
	//
	//========================================================

	opt.allowed.push_back("varPos");
	opt.allowed.push_back("varID");
	opt.allowed.push_back("samplingLevel");


	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("hBondFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("seed");

	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("outputFile");
	opt.allowed.push_back("logFile");
	opt.allowed.push_back("rerunConfigFile");
	opt.allowed.push_back("fastaFile");


	//========================================================
	//
	//===========Dependendent options=========================
	//
	//========================================================

	//These options depend on an output directory	
/*
	opt.dependent.push_back("outputDir");
	opt.dependent.push_back("outputFile");
	opt.dependent.push_back("logFile");
	opt.dependent.push_back("rerunConfigFile");
	opt.dependent.push_back("fastaFile");
	
	opt.dependentSet.push_back(opt.dependent);
	opt.dependent.clear();
*/	
	//Placeholder for more dependent options
	


	//========================================================
	//
	//===========Set the provided options=====================
	//
	//========================================================



	OptionParser OP;

	OP.readArgv(_argc, _argv);

	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.setDependentOn(opt.dependentSet);
	OP.autoExtendOptions();

	//=======================================================
	// CHECK THE GIVEN OPTIONS
	//=======================================================
	
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	//=======================================================
	// CHECK THE GIVEN OPTIONS
	//=======================================================

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages += "pdb file not specified!\n";
		opt.errorFlag = true;
	}
	
	opt.varPos = OP.getIntVector("varPos");
	if (OP.fail()) {
		opt.warningMessages += "variable positions not specified!\n";
		opt.warningFlag = true;
	}

	opt.varID = OP.getStringVector("varID");
	if (OP.fail()) {
		opt.warningMessages += "Identities not specified! Default to all amino acids besides proline.\n";
		opt.warningFlag = true;
		std::vector<string> aaTmp (20,"");
		aaTmp[0] = "ALA";
		aaTmp[1] = "ARG";
		aaTmp[2] = "ASN";
		aaTmp[3] = "ASP";
		aaTmp[4] = "CYS";
		aaTmp[5] = "GLN";
		aaTmp[6] = "GLU";
		aaTmp[7] = "GLY";
		aaTmp[8] = "HSD";
		aaTmp[9] = "HSE";
		aaTmp[10] = "ILE";
		aaTmp[11] = "LEU";
		aaTmp[12] = "LYS";
		aaTmp[13] = "MET";
		aaTmp[14] = "PHE";
		//aaTmp[15] = "PRO";
		aaTmp[15] = "SER";
		aaTmp[16] = "THR";
		aaTmp[17] = "TRP";
		aaTmp[18] = "TYR";
		aaTmp[19] = "VAL";
		
		opt.varID = aaTmp;

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


	opt.samplingLevel = OP.getString("samplingLevel");
	if (OP.fail()) {
		opt.warningMessages += "Sampling Level not specified! Default to SL60.00\n";
		opt.warningFlag = true;
		opt.samplingLevel = "SL60.00";
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
	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_SOLV";
		if(SYSENV.isDefined(envVar)) {
			opt.solvFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n"	;
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
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(SYSENV.isDefined(envVar)) {
			opt.rotLibFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "Seed not specified!  Will start up the RNG using a time based seed.\n";
		opt.warningFlag = true;
		opt.seed = 0;
	}


	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "Output directory for pdbs not specified!  Will not write PDB files.\n";
		opt.warningFlag = true;
	}

	opt.outputFile = OP.getString("outputFile");
	if (OP.fail()) {
		opt.warningMessages += "Output file not specified!  Will print to terminal instead.\n";
		opt.warningFlag = true;
	}

	opt.logFile = OP.getString("logFile");
	if (OP.fail()) {
		opt.warningMessages += "Log file not specified!  Will not print to a log file.\n";
		opt.warningFlag = true;
	}
	
	opt.rerunConfigFile = OP.getString("rerunConfigFile");
	if (OP.fail()) {
		opt.warningMessages += "rerun configuration file not specified!  Will not make a file to rerun this job.\n";
		opt.warningFlag = true;
	}

	opt.fastaFile = OP.getString("fastaFile");
	if(OP.fail()) {
		opt.warningMessages += "Fasta file not specified!  Will not write accepted sequences\n";
		opt.warningFlag = true;
	}


	opt.rerunConfig = "########################################################\n";
	opt.rerunConfig += "#  This configuration file was automatically generated,\n";
	opt.rerunConfig += "#  it will rerun this job with the same options.       \n";
	opt.rerunConfig += "#\n";
	opt.rerunConfig += "#  Run as:\n";
	opt.rerunConfig += "#\n";
	opt.rerunConfig += "#    % " + programName + " --configfile " + opt.rerunConfigFile + "\n";
	opt.rerunConfig += "#\n";
	opt.rerunConfig += "#  Job started on " + (string)ctime(&start_time);
	opt.rerunConfig += "#  on host " + opt.host + ", path " + opt.pwd + "\n";
	if (opt.seed == 0) {
		opt.rerunConfig += "#  seed " + MslTools::intToString(opt.seed) + " (time based)\n";
	} else {
		opt.rerunConfig += "#  seed " + MslTools::intToString(opt.seed) + "\n";
	} 
	opt.rerunConfig += "########################################################\n";
	opt.rerunConfig += "\n";
	opt.rerunConfig += OP.getConfFile();	



	return opt;
}


	
std::vector< std::vector<bool> > buildResidueMask(std::vector<unsigned int> _residueState, std::vector< std::vector<unsigned int> > _residueRotamers) {
	
	std::vector< std::vector<bool> > resMask (_residueState.size());

	//for each position
	for (unsigned int i = 0; i < _residueState.size(); i++) {
		unsigned int currentResidue = _residueState[i];
		if (currentResidue >= _residueRotamers[i].size()) {
			cerr << "ERROR: the current residue number exceeds the number of residues for position " << i << endl; 
			exit(100);
		}
		//for each identity in that position
		for (unsigned int j = 0; j < _residueRotamers[i].size(); j++) {
			if (j==currentResidue) {
				for (unsigned int k = 0; k < _residueRotamers[i][j]; k++) {
					resMask[i].push_back(true);
				}
			} else {
				for (unsigned int k = 0; k < _residueRotamers[i][j]; k++) {
					resMask[i].push_back(false);
				}
			}
		}

		//check to make sure that there are some rotamers set to "true"
		uint counter = 0;
		for (unsigned int j = 0; j < resMask[i].size(); j++) {
			if (resMask[i][j] == 1) {
				counter++;
			}
		}
		if (counter == 0) {
			cout << "ERROR AT POSITION: " << i << endl;
			cout << "Current Residue: " << currentResidue << endl;
			cout << "resRots at this position: " << endl;
			for (uint k = 0; k < _residueRotamers[i].size(); k++) {
				cout << _residueRotamers[i][k] << " ";
			}
			cout << endl;
			cout << "resMask at this position: " << endl;
			for (uint k = 0; k < resMask[i].size(); k++) {
				cout << resMask[i][k] << " ";
			}
			cout << endl;
		}
	}
/*
	cout << "resMask" << endl;
	cout << "================================" << endl;	
	for (unsigned int i = 0; i < resMask.size(); i++) {
		cout << "Position " << i << " ";
		for (unsigned int j = 0; j < resMask[i].size(); j++) {
			cout << resMask[i][j];
		}
		cout << endl;
	}
*/
	return resMask;
}
