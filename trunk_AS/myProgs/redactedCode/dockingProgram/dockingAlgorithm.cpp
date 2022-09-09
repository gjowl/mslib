// This program will accept a template structure of a CATM dimer
// as input and perform Monte Carlo repacks on the helices. Symmetry will be
// enforced between the chains. The templates were built from CATM??
// The distance restraints were added from sort-seq mutational scanning results. 

#include <iostream>
#include <time.h>
#include <unistd.h>

#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "OptionParser.h"		
#include "MslTools.h"			
#include "Transforms.h"			
#include "AtomSelection.h"		
#include "CartesianGeometry.h"		       
#include "System.h"			
#include "SysEnv.h"			
#include "SystemRotamerLoader.h"	
#include "SelfPairManager.h"		
#include "RandomNumberGenerator.h"	
#include "PrincipleComponentAnalysis.h"
#include "MonteCarloManager.h"
#include "AtomBondBuilder.h"
#include "SigmoidInteraction.h"
#include "HelixGenerator.h"

using namespace MSL;
using namespace std;
static SysEnv SYSENV;
string programName = "dockingAlgorithm";
string programDescription = "This program performs a Monte Carlo repack of a candidate geometry for a transmembrane helix bundle with optional sigmoidal distance restraints";
string programAuthor = "Samantha Anderson; Gilbert Loiseau edited this to make a hetero version for this docking algorithm. This version uses many of the defaults for docking designed sequences";
string programVersion = "1.0.3";
string programDate = "29 October 2021";
	//4/3/20 added crossingangle limits
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, mcTime;
double diffTime;

/******************************************
 *  
 *  =======  PREDECLARATION OF THE FUNCTIONS  =======
 *
 ******************************************/
vector < vector<Atom*> > getInterHelicalHbonds(EnergySet* & _ESet);
bool alignOnZ (AtomPointerVector &_apv, Transforms &_trans, CartesianPoint &_atom1, CartesianPoint &_atom1Prime, CartesianPoint &_atom2, CartesianPoint &_atom2Prime);
void applyRotationalSymmetry (System &_sys, Transforms &_tr, CartesianPoint &_axis);
void readGeometryFile(string _filename, vector<string>& _fileVec);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB, Transforms &_trans);
std::vector < std::vector <bool> > getActiveMask(System &_sys);
AtomPointerVector calculateHelicalAxes (AtomPointerVector &_helix);
bool duplicateHelix (System &_sys, int _helices);
double calculateMonomerE (System &_sys, Transforms &_trans, SelfPairManager &_spm);
void deleteTerminalHydrogenBondInteractions(System &_sys, ofstream &_lout, bool &_verbose);
double getStandardNormal(RandomNumberGenerator& RNG) {
	//average 0, std dev ~0.1
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}
double getBoxMuller (RandomNumberGenerator &_RNG) {
	//Box-Muller Transformation of uniform random numbers to standard normal distribution
	//Mean 0, sigma 1
	double U = _RNG.getRandomDouble();
	double V = _RNG.getRandomDouble();
	double stdNormal = sqrt(-2*log(U))*cos(6.28*V);
	return stdNormal;
}
double getRandomNormal(RandomNumberGenerator &_RNG, double _avg, double _stdev) {
	return _avg + _stdev*getBoxMuller(_RNG);
}
void changeRandomRotamer (RandomNumberGenerator &_RNG, std::vector< std::vector <unsigned int> > &_rotIndex, std::vector< unsigned int> &_activeRots);
	
string shiftHelix (AtomPointerVector &_helix, RandomNumberGenerator &_RNG, Transforms &_trans, bool _shiftRotCenter);
string transformTetramer (AtomPointerVector &_ftsb, AtomPointerVector &_ftsbPrime, AtomPointerVector &_bCoor, AtomPointerVector &_ftsl, AtomPointerVector &_ftslPrime, AtomPointerVector &_lCoor, Transforms &_trans, std::vector<double> &_params);
string printParameters (std::vector <double> _params);
bool maskInteractions(std::vector<SigmoidInteraction*>);
//void createHelixBackbone (std::vector<string> &_chainSeqs, std::vector<int> &_numChains, std::vector<int> &_chainStartNum, AtomPointerVector &_apv, HelixGenerator &_hg);
void createHelixBackbone (std::vector<string> &_chainSeqs, int _numChains, std::vector<int> &_chainStartNum, AtomPointerVector &_apv, HelixGenerator &_hg);
void transformHelix (double _axialRotation, double _zShift, double _crossingAngle, double _xShift, AtomPointerVector &_helix, Transforms &_trans);
string convertToThreeLetterSequence(string _seq);
//string convertToFullPolymerSequence(std::vector<string> &_seqs, std::vector<int> &_numChains, std::vector<int> &_chainStartNum);
string convertToFullPolymerSequence(std::vector<string> &_seqs, int _numChains, std::vector<int> &_chainStartNum);
void adjustParameters (RandomNumberGenerator &_RNG, std::vector<double> &_params, std::vector<double> &_steps, double _paramAdjustment, bool _limitCross, double _limitCrossValue, ofstream &_lout);
//string transformBundle (System &_sys, Transforms &_trans, std::vector<int> &_numChains, std::vector<double> &_params);
string transformBundle (System &_sys, Transforms &_trans, int _numChains, std::vector<double> &_params);
string transformBundleToStartingPosition (System &_sys, Transforms &_trans, int _numChains, std::vector<double> &_params);
string transformBundleRandomly (System &_sys, RandomNumberGenerator &_RNG, Transforms &_trans, std::vector<double> &_params, std::vector<double> &_steps, double _paramAdjustment, bool _limitCross, double _limitCrossValue, int _numChains, ofstream &_out);

/******************************************
 *  
 *  =======  OPTIONS  =======
 *
 ******************************************/

struct Options {
	bool verbose;
	bool printIntermediatePdbs;
	bool zipOutputFiles;

	//Parameters required for each set of chains
	std::vector<string> chainSeq;
	int numChains;
	std::vector<int> chainStartNum;
	std::vector<string> chainRotLevels;

	//Output Parameters
	string output;
	string outputName;
	string outputDir;
	int numReplicates;
	int replicateStartNum;
	int numTrajectoryModels;
	
	//Rotamer Parameters
	string rotLevel;

	//MC Parameters
	double startT;
	double endT;
	int mcCycles;
	int mcShape;
	string mcShapeString;
	int mcMaxRejects;
	int convergedSteps;
	double convergedE;
	double probRepack;
	bool shiftRotCenter;
	int numParamMoves;
	int numRotMoves;
	int seed;

	//CSB Parameters
	double ctonnb;
	double ctofnb;
	double cutnb;

	//Orientation Parameter ranges 
	double axialRotStart;
	double axialRotEnd;
	double zShiftStart;
	double zShiftEnd;
	double crossAngleStart;
	double crossAngleEnd;
	double xShiftStart;
	double xShiftEnd;

	//Orientation Parameters for each set of TM helices
	std::vector<double> params;
	double axialRot;
	double zShift;
	double crossingAngle;
	double xShift;
	double designAxialRotation;
	double designZShift;
	double designCrossingAngle;
	double designXShift;

	//Step sizes for each of the above parameters
	double axialRot_step;
	double zShift_step;
	double crossingAngle_step;
	double xShift_step;
	std::vector<double> steps;
	double paramAdjustment;
	bool limitCross;
	double limitCrossValue;
	bool xShiftFirst;

	//Energy Parameters
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	double vdwRadius;
	bool deleteTerminalHbonds;

	//Add restraints to different pairs of atoms
	std::vector<string> restraints;		//Vector of Atom ID pairs to add restraints
	std::vector<double> restraintWeights;	//Vector of restraint strengths
	std::vector<double> restraintSlopes;
	std::vector<double> restraintDists;	//Vector of optimal distances for each restraint
	std::vector<double> restraintIntercepts;
	std::vector<string> restraintGroups;
	double defaultRestraintWeight;
	double defaultRestraintSlope;
	double defaultRestraintDist;
	double defaultRestraintIntercept;
                           
	double repackThreshold;
	//vector<int> rotCount;	
	//Parameter Files
	string topFile;
	string parFile;
	string hbondFile;
	string solvFile;
	string rotLibFile;
	string bbqFile;
	string designPdb;

	//Config Files Options
	vector<int> rotamerSamplingVector;
	
	// input monomerEnergy
	double monomer;
	double monoVdw;
	double monoHbond;
	double monoIMM1;
	
	// input previous dimerEnergy
	double dimer;
	double dimerVdw;
	double dimerHbond;
	double dimerIMM1;

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
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);
void setupOutputDirectory(Options &_opt){
	// for running on the server chtc
	_opt.outputDir = get_current_dir_name();
	_opt.outputDir = _opt.outputDir + "/" + _opt.outputName;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

/******************************************
 *  
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main (int argc, char *argv[]) {
	/******************************************************************************
	 *                  === SET UP OPTIONS AND BOOKKEEPING ===
	 ******************************************************************************/
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
	setupOutputDirectory(opt);
	// Make log file if performing verbose
	string logFile = opt.outputDir + "/docking.log";
	ofstream lout;
	if (opt.verbose){
		lout.open(logFile.c_str());
		if(!lout.is_open()) {
			cerr << "Unable to open " << logFile << endl;
			exit(0);
		}
		lout << opt.rerunConf << endl;
	}
	
	// Make List file
	string listFile = opt.outputDir + "/dockingEnergies.list";
	ofstream eOut;
	eOut.open(listFile.c_str());
	if(!eOut.is_open()) {
		cerr << "unable to open " << listFile << endl;
		exit(0);
	}

	/******************************************************************************
	 *                  === DECLARE SYSTEM ===
	 ******************************************************************************/	
	System sys;
	HelixGenerator hg(opt.bbqFile);

	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);
	if (opt.verbose){lout << "RNG Seed: " << RNG.getSeed() << endl;}
	
	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");//TODO: decide if I need to add in the membrane terms
	CSB.setVdwRescalingFactor(opt.vdwRadius);
	HydrogenBondBuilder hb(sys, opt.hbondFile);

	/******************************************************************************
	 *                  === BUILD TM BUNDLE SYSTEM ===
	 ******************************************************************************/

	//TODO generalize this to make bundles similar to modelCoiledCoils
	AtomPointerVector centered_zAlignedAPV;
	createHelixBackbone (opt.chainSeq, opt.chainSeq.size(), opt.chainStartNum, 
						centered_zAlignedAPV, hg);

	if (opt.verbose){lout << "Backbone Generated!" << endl;}
	
	string fullPS = convertToFullPolymerSequence(opt.chainSeq, opt.chainSeq.size(), 
			opt.chainStartNum);
	cout << fullPS << endl;
	if (opt.verbose){lout << fullPS << endl;}
	PolymerSequence PS;
	PS.setSequence(fullPS);
	CSB.buildSystem(PS);
	if (opt.verbose){lout << "Sequence built! Preparing to assign coordinates..." << endl;}
	sys.wipeAllCoordinates();
	//sys.assignCoordinates(centered_zAlignedAPV, false);
	sys.assignCoordinates(centered_zAlignedAPV, true);
	sys.buildAllAtoms();
	if (opt.verbose){lout << "Coordinates assigned to helix bundle!\n";}
	if(opt.printIntermediatePdbs){
		string centered_zAlignedFile = opt.outputDir + "/centered_zAligned.pdb";
		sys.writePdb(centered_zAlignedFile);
	}
	
	//TODO Make this an option
	if (hb.buildInteractions(30)){
		if (opt.verbose){lout << "Built Interactions" << endl;}
	} else {
		cerr << "Unable to build hydrogen bond interactions." << endl;	
	}

	/******************************************************************************
	 *                  === LOAD ROTAMERS & SET-UP SPM ===
	 ******************************************************************************/
	SelfPairManager spm;

	spm.saveEnergiesByTerm(true);
	spm.setRandomNumberGenerator(&RNG);
	spm.seed(RNG.getSeed());
	spm.setOnTheFly(true);
	spm.setVerbose(false);
	if (opt.verbose){lout << "Loading rotamers..." << endl;}
	/*if (opt.chainRotLevels.size() != 0) {
		//Load position-specific rotamer levels
		if (opt.verbose){lout << "Using position specific rotamers identified chainRotLevels" << endl;}
		unsigned int chainCounter = 0;
		for (uint i = 0; i < opt.chainSeq.size(); i++) {
			std::vector<string> posRotLevels = MslTools::tokenizeAndTrim(opt.chainRotLevels[i]);
			if (posRotLevels.size() != sys.getChain(chainCounter).positionSize()) {
				cerr << "ERROR 4476: For chain bundle " << i << " there was an incorrect number of rotamer positions!" << endl;
				exit(4476);
			}

			for (uint j = 0; j < opt.numChains[i]; j++) {
				if (opt.verbose){lout <<"Bundle " << i << ", chain " << j << endl;}
				for (uint k = 0; k < sys.getChain(chainCounter).positionSize(); k++) {
					Position &pos = sys.getChain(chainCounter).getPosition(k);
					string chainId = sys.getChain(chainCounter).getChainId();
					string posRot = posRotLevels[k];

					if (opt.verbose){lout << pos.getPositionId() << ", " << pos.getResidueName() << ": " << posRot << endl;}
					if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
						if (!sysRot.loadRotamers(&pos, pos.getResidueName(), posRot)) { 
							cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
						}
					}

				}
				chainCounter++;
			}
		}
	} else {
		if (opt.verbose){lout << "Using rotamer level " << opt.rotLevel << " for all positions" << endl;}
		for (uint k=0; k < sys.positionSize(); k++) {
			Position &pos = sys.getPosition(k);
			string posRot = opt.rotLevel;
			if (opt.verbose){lout << pos.getPositionId() << ", " << pos.getResidueName() << ": " << posRot << endl;}
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!sysRot.loadRotamers(&pos, pos.getResidueName(), opt.rotLevel)) { 
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}
	}*///Potentially Redacted: It runs fast enough to not use different rotamer levels for now
	if (opt.verbose){lout << "Using rotamer level " << opt.rotLevel << " for all positions" << endl;}
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);
		string posRot = opt.rotLevel;
		if (opt.verbose){lout << pos.getPositionId() << ", " << pos.getResidueName() << ": " << posRot << endl;}
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!sysRot.loadRotamers(&pos, pos.getResidueName(), opt.rotLevel)) { 
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
	if (opt.verbose){lout << "Rotamers loaded!\n";}
	sys.buildAllAtoms();
	//TODO make these options and set the list option default to be bigger
	CSB.updateNonBonded(opt.ctonnb, opt.ctofnb, opt.cutnb);//change these to the hardcoded values I'm using in design
	sys.updateVariablePositions();
//	//Build a vector of position:rotamer indices for later
//	std::vector < std::vector< unsigned int > > rotamerIndex;
//	for (uint i = 0; i < sys.getVariablePositions().size(); i++) {
//		//Positions
//		unsigned int pos = sys.getVariablePositions()[i];
//		for (uint rot = 0; rot < sys.getPosition(pos).getTotalNumberOfRotamers(); rot++) {
//			//Rotamers
//			std::vector<unsigned int> rotIndex (2,999);
//			rotIndex[0] = i; //Ith variable position
//			rotIndex[1] = rot; //Jth rotamer
//			rotamerIndex.push_back(rotIndex);
//		}
//	}	
	if (opt.verbose){lout << sys << endl;}

	/******************************************************************************
	 * ==== REMOVE THE INTERACTIONS WITHIN THE BACKBONE OF EACH HELIX SINCE THEY DO NOT CHANGE ====
	 ******************************************************************************/
	if (opt.verbose){lout << "Removing interactions within helix backbones..." << endl;}
	EnergySet* pESet = sys.getEnergySet();
	pESet->setWeight("CHARMM_VDW", opt.weight_vdw);
	pESet->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	
	AtomSelection sel(sys.getAllAtomPointers());
	for (unsigned int i = 0; i < sys.chainSize(); i++) {
		string chainId = sys.getChain(i).getChainId();
		string sele = "chain" + chainId;
		string selection = "chain" + chainId + ", chain " + chainId + " and name C+CA+N+O+HA+HA1+HA2+HN+H+HT1+HT2+HT3+OT1+OT2";
		if (opt.verbose){lout << selection << endl;}
		AtomPointerVector chaini = sel.select(selection);
		pESet->deleteInteractionsWithinSelection(sele);
	}

	/******************************************************************************
	*             Delete hydrogen bonds in the 4 (+4)  residues
	* 	      Residue numbers 1,2,3,4 and last,last-1,last-2,last-3
	*******************************************************************************/
	if (opt.deleteTerminalHbonds) {
		if (opt.verbose){lout << "\nDeleting hbonds from the N- and C-terminal extension residues\n" << endl;}
		deleteTerminalHydrogenBondInteractions(sys,lout, opt.verbose);
	}

	/******************************************************************************
	 *               === ADD SIGMOID CONSTRAINTS ===
	 ******************************************************************************/
	std::map<std::string, std::vector<SigmoidInteraction*> > restraintGroups;

	//for (int i = 0; i < opt.restraints.size(); i++) {
	//	std::vector<string> pair = MslTools::tokenize(opt.restraints[i]);
	//	string atom1 = pair[0];
	//	string atom2 = pair[1];
	//	if (sys.atomExists(atom1) && sys.atomExists(atom2)) {
	//		if (opt.restraintGroups.empty()) {
	//			if (opt.verbose){lout << "Adding restraint: " << atom1 << " " << atom2 << endl;}
	//		} else {
	//			if (opt.verbose){lout << "Adding restraint: " << atom1 << " " << atom2 << " to group " << opt.restraintGroups[i] << endl;}
	//		}
	//		if (opt.verbose){lout << opt.restraints[i] << "\t";}
	//		if (opt.verbose){lout << opt.restraintWeights[i] << "\t";}
	//		if (opt.verbose){lout << opt.restraintSlopes[i] << "\t";}
	//		if (opt.verbose){lout << opt.restraintDists[i] << "\t";}
	//		if (opt.verbose){lout << opt.restraintIntercepts[i] << "\t";}
	//		if (opt.verbose){lout << endl;}
	//		SigmoidInteraction *restraint = new SigmoidInteraction(sys.getAtom(atom1), sys.getAtom(atom2), opt.restraintWeights[i], opt.restraintSlopes[i], opt.restraintDists[i], opt.restraintIntercepts[i]);

	//		pESet->addInteraction(restraint);

	//		//Add the interaction pointer to the map for later masking
	//		if (!opt.restraintGroups.empty()){
	//			restraintGroups[opt.restraintGroups[i]].push_back(restraint);
	//		}

	//	} else {
	//		if (opt.verbose){lout << "ERROR 23859: cannot add restraint between " << atom1 << " and " << atom2 << endl;}
	//		cerr << "ERROR 23859: cannot add restraint between " << atom1 << " and " << atom2 << endl;
	//		exit (23859);
	//	}
	//}

	//if (opt.verbose){lout << "Restraints added!\n";}

	sys.saveAltCoor("centered_zAligned");

	/******************************************************************************
	 *                    === PREPARE OUTPUT FILE ===
	 ******************************************************************************/
	string csv = opt.outputDir+"/summary.csv";
	ofstream csvOut;
	csvOut.open(csv.c_str());
	if(!csvOut.is_open()) {
		cerr << "Unable to open " << csv << endl;
		exit(0);
	}
	csvOut << "Sequence,Replicate,Total,Dimer,Monomer,xShift,crossingAngle,axialRotation,zShift" << endl;
	cout << "Sequence,Replicate,Total,Dimer,Monomer,xShift,crossingAngle,axialRotation,zShift" << endl;
	
	/******************************************************************************
	 *               === START LOOP FOR REPLICATES ===
	 ******************************************************************************/
	int bestReplicate = 1;
	double veryBestE = 0;
	for (int rep = 0; rep < opt.numReplicates; rep++){
		// Make Replicate Output File
		stringstream tmp;
		tmp << setw(7) << setfill('0') << (rep + opt.replicateStartNum);
		string repString = tmp.str();
		string repFile = opt.outputDir + "/" + repString + ".out";
		ofstream oRep;
		oRep.open(repFile.c_str());
		if(!oRep.is_open()) {
			cerr << "Unable to open " << repFile << endl;
			exit(0);
		}
		oRep << "Replicate " << rep << endl;
		if (opt.verbose){lout << endl << "Starting Replicate " << rep << endl;}

		/******************************************************************************
		 *               === CALCULATE STARTING POSITION ===
		 ******************************************************************************/
		double axialRot;
		double zShift;
		double crossingAngle;
		double xShift;

		if (opt.axialRot == 1234567890){
			axialRot = RNG.getRandomDouble(opt.axialRotStart, opt.axialRotEnd);
		} else {
			axialRot = opt.axialRot;
		}
		if (opt.zShift == 1234567890){
			zShift = RNG.getRandomDouble(opt.zShiftStart, opt.zShiftEnd);
		} else {
			zShift = opt.zShift;
		}
		if (opt.crossingAngle == 1234567890){
			crossingAngle = RNG.getRandomDouble(opt.crossAngleStart, opt.crossAngleEnd);
		} else {
			crossingAngle = opt.crossingAngle;
		}
		if (opt.xShift == 1234567890){		
			xShift = RNG.getRandomDouble(opt.xShiftStart, opt.xShiftEnd);
		} else {
			xShift = opt.xShift;
		}
		vector<double> params;
		params.push_back(axialRot);
		params.push_back(zShift);
		params.push_back(crossingAngle);
		params.push_back(xShift);
		if (opt.verbose){
			lout << "Starting Orientation Parameters: " << endl;
			lout << printParameters(params)<< endl;
		}
	
		/******************************************************************************
		 *               === MOVE TO STARTING POSITION ===
		 ******************************************************************************/
	
		//Transform the bundle to initial coordinates
		transformBundleToStartingPosition (sys, trans, opt.chainSeq.size(), params);
		if (opt.printIntermediatePdbs){
			string startCoorFile = opt.outputDir + "/startCoor.pdb";
			sys.writePdb(startCoorFile);
		}
		sys.saveAltCoor("startCoor");
		CSB.updateNonBonded(opt.ctonnb, opt.ctofnb, opt.cutnb);	
		sys.buildAllAtoms();
		if (opt.verbose){lout << sys << endl;}
		oRep << sys << endl;
		spm.setSystem(&sys);	
		if (opt.verbose){lout << "Before initial repack:" <<endl;}
		if (opt.verbose){lout << sys.calcEnergy() << endl;}
		if (opt.verbose){lout << sys.getEnergySummary() << endl;}
		double startE = sys.calcEnergy();
		if (opt.verbose){lout << "Calculating monomer energy..." << endl;}
		double monomerE = calculateMonomerE(sys, trans, spm);//this may already work for hetero
		if (opt.printIntermediatePdbs){
			string monomerPdb = opt.outputDir + "/monomer.pdb";
			sys.writePdb(monomerPdb);
		}
		if (opt.verbose){lout << "Monomer Energy: " << monomerE << endl;}
		oRep << "Monomer Energy: " << monomerE << endl;
		string monomerSummary = spm.getSummary(spm.getMinStates()[0]);
		std::vector<unsigned int> repackRots = spm.getMinStates()[0];	
		if (opt.verbose){lout << monomerSummary << endl;}
		oRep << monomerSummary << endl;
		if (opt.verbose){lout << "Repacking TM Bundle..." << endl;}
		sys.applySavedCoor("startCoor");
		
		std::vector <double> initialParams = params;
	
		spm.calculateEnergies();
		spm.runGreedyOptimizer(30);
		repackRots = spm.getMinStates()[0];
		sys.setActiveRotamers(repackRots);	
		double repackE = spm.getMinBound()[0];
	
		if (opt.verbose){lout << "Initial Repack Energy: " << endl;}
		if (opt.verbose){lout << sys.calcEnergy() << endl;}
		if (opt.verbose){lout << sys.getEnergySummary() << endl;}
		
		sys.setActiveRotamers(repackRots);
		sys.saveAltCoor("initialRepack");
		PDBWriter trajectoryWrite;

		if (opt.printIntermediatePdbs){
			string initializePdb = opt.outputDir + "/initialRepack.pdb";	
			sys.writePdb(initializePdb);
			string trajectorypdb = opt.outputDir + "/MC_Trajectory.pdb";
			trajectoryWrite.open(trajectorypdb);
			trajectoryWrite.addRemark("Initial Repack");
			trajectoryWrite.addRemark(printParameters(params));
			trajectoryWrite.writeREMARKS();
			trajectoryWrite.write(sys.getAtomPointers(), true, false, true);
		}
	
		double deltaE = repackE - monomerE;
		if (opt.verbose){lout << "Current Energy: " << deltaE << endl;}
		sys.applySavedCoor("centered_zAligned");
		
		/******************************************************************************
		 * ==== MOVE HELICES CLOSER UNTIL ENERGY GETS WORSE ====
		 ******************************************************************************/
		std::vector<double> xParams = params;
		if(opt.xShiftFirst){
			double bestEnergy = deltaE;
			bool closest = false;
			double bestX = xParams[3];
			lout << "Old X\tNewX\tMonomerE\tCurrentE\toldDE\tNewDE\tPass" << endl;
			while(!closest){
				double newX = xParams[3] - opt.xShift_step;
				lout << xParams[3] << "\t" << newX << "\t";
				xParams[3] = newX;
				sys.applySavedCoor("centered_zAligned");
				transformBundle (sys, trans, opt.chainSeq.size(), xParams);
				sys.setActiveRotamers(repackRots);	
			CSB.updateNonBonded(opt.ctonnb, opt.ctofnb, opt.cutnb);
				double currentE = sys.calcEnergy() - monomerE;
				lout << monomerE << "\t" << sys.calcEnergy() << "\t";
				lout << bestEnergy << "\t" << currentE << "\t";
				//lout << endl << sys.getEnergySummary() << endl;
				if (currentE <= (bestEnergy + 3)){
					lout << "\tBetter" << endl;
					bestEnergy = currentE;
					bestX = newX;
					sys.saveAltCoor("closer");
				} else {
					lout << "\tWorse" << endl;
					closest = true;
				}
			}
			xParams[3] = bestX;
			lout << "Closest x: " << xParams[3] << "\tEnergy: " << bestEnergy << endl;
			sys.applySavedCoor("closer");
			if (opt.printIntermediatePdbs){
				string startCoorFile = opt.outputDir + "/xShifted.pdb";
				sys.writePdb(startCoorFile);
			}

		}
		/******************************************************************************
		 * ==== BEGIN LOCAL MC REPACK OF TM BUNDLE====
		 ******************************************************************************/
		double bestEnergy = deltaE;
		double currentEnergy = deltaE;
		string bestSummary;
		std::vector<unsigned int> bestRotamers = repackRots;
		std::vector<unsigned int> currentRotamers = repackRots;
		std::vector<double> bestParams = xParams;
		MonteCarloManager MC(opt.startT, opt.endT, opt.mcCycles, opt.mcShape, opt.mcMaxRejects, opt.convergedSteps, opt.convergedE);
		MC.setRandomNumberGenerator(&RNG);
		MC.seed(RNG.getSeed());
		MC.setEner(bestEnergy);
	
		if (opt.verbose){lout << "Begin Monte Carlo repacks..." << endl;}
	
		//Stat counters
		double repackTime = 0;
		int repackCounter = 0;
		int acceptCounter = 0;
		int rejectCounter = 0;
		int moveCounter = 0;
		int currentGap = 0;
		int maxGap = 0;
		int avgGap = 0;
		int nthModel = opt.mcCycles / opt.numTrajectoryModels;
		if (nthModel < 1) {
			nthModel = 1;
		}	
		time(&mcTime);
		//TODO updateNonbonded after moves
		while (!MC.getComplete()) {
			moveCounter++;
			std::vector<double> mcParams = bestParams;	
			string moveMade = "";
			//time_t moveStart = time(0);
			sys.applySavedCoor("centered_zAligned");
			//TODO write in limited crossing angles and xshiftfirst
			//adjustParameters(RNG, mcParams, opt.steps, opt.paramAdjustment, opt.limitCross, opt.limitCrossValue, lout);
			moveMade = transformBundleRandomly (sys, RNG, trans, mcParams, opt.steps, opt.paramAdjustment, opt.limitCross, opt.limitCrossValue,  opt.chainSeq.size(), lout);//added adjust parameters to this function to make the changes less symettrical
			if (opt.printIntermediatePdbs){
				trajectoryWrite.write(sys.getAtomPointers(), true, false, true);	
			}

			CSB.updateNonBonded(opt.ctonnb, opt.ctofnb, opt.cutnb);
			spm.setSystem(&sys);
			if (RNG.getRandomDouble() < opt.probRepack) {
				repackCounter++;
				time_t repackStart = time(0);
				if (opt.verbose){lout << "Repacking side chains...";}
				spm.calculateEnergies();
				spm.runGreedyOptimizer(10);
				currentRotamers = spm.getMinStates()[0];
				time_t repackEnd = time(0);
				repackTime += difftime(repackEnd, repackStart);
				if (opt.verbose){lout << difftime(repackEnd, repackStart) << " seconds." << endl;}
			} 	
			for (std::map<string, std::vector<SigmoidInteraction*> >::iterator it = restraintGroups.begin(); it!= restraintGroups.end(); ++it) {
				maskInteractions(it->second);
			}

			sys.setActiveRotamers(currentRotamers);
			currentEnergy = sys.calcEnergy() - monomerE;
	
			//time_t moveEnd = time(0);
			//cout << "Move took " << moveEnd - moveStart << " seconds" << endl;
			if (!MC.accept(currentEnergy)) {
				//lout << "state rejected. Energy : " << currentEnergy << endl;
				if (opt.verbose){
					lout << "Cycle: " << MC.getCurrentCycle() << " Temperature ";
					lout << MC.getCurrentT() << " state rejected! Energy: ";
					lout << currentEnergy << " " << sys.calcEnergy() << endl;
				}
				rejectCounter++;
				currentGap++;
				currentRotamers = bestRotamers;
				sys.applySavedCoor("centered_zAligned");
			} else {
				acceptCounter++;
				if (opt.verbose){lout << "**********************CYCLE " << MC.getCurrentCycle() << "**********************" << endl;}
				if (opt.verbose){lout << moveMade << endl;}
	
				//if (acceptCounter % 10 == 0) {
				//	//Repack side chains every few accepts
				//	spm.calculateEnergies();				
				//	spm.runGreedyOptimizer(10);
				//	double repackEnergy = spm.getMinBound()[0] - monomerE;
				//	if (repackEnergy < currentEnergy) {
				//		currentRotamers = spm.getMinStates()[0];
				//		currentEnergy = repackEnergy;
				//	}
				//}
				avgGap += currentGap;
				if (currentGap > maxGap) {
					maxGap = currentGap;
				}
				currentGap = 0;
				bestEnergy = currentEnergy;
				MC.setEner(currentEnergy);	//We are improving the energy after accepting in some cases, so this will prevent jumping back
				bestParams = mcParams;
				bestRotamers = currentRotamers;
				sys.setActiveRotamers(bestRotamers);
	
				if (opt.verbose){lout << "Cycle: " << MC.getCurrentCycle() << " Temperature " << MC.getCurrentT() << " state accepted! Energy: " << currentEnergy << " " << sys.calcEnergy() << endl;}
				if (opt.verbose){lout << sys.getEnergySummary() << endl;}
				bestSummary = sys.getEnergySummary();
				sys.saveAltCoor("bestState");
				//if (opt.verbose){lout << "==================================SYSTEM==========================" << endl;}
				//if (opt.verbose){lout << sys.getEnergySummary() << endl;}
				if (acceptCounter % (nthModel) == 0) {
					if (opt.verbose){lout << "Parameters: \n";}
					if (opt.verbose){lout << moveMade << endl;}
					if (opt.printIntermediatePdbs){
						if (opt.verbose){lout << "Model written to trajectory file!\n";}				
						trajectoryWrite.clearRemarks();
						trajectoryWrite.addRemark(moveMade);
						trajectoryWrite.addRemark(bestSummary);
						trajectoryWrite.writeREMARKS();
						trajectoryWrite.write(sys.getAtomPointers(), true, false, true);	
					}
				}
			}
		}
	
		if (opt.printIntermediatePdbs){trajectoryWrite.close();}
		if (opt.verbose){lout << "MC Repack Complete!" << endl;}
		if (opt.verbose){lout << MC.getReasonCompleted() << endl;}
		sys.applySavedCoor("bestState");
		sys.calcEnergy();
		if (opt.verbose){lout << sys.getEnergySummary() << endl;}
		//if (opt.verbose){lout << spm.getSummary(bestRotamers) << endl;}
		if (opt.verbose){lout << "Best Energy: " << bestEnergy << endl;}
	
		oRep << "MC Repack Complete!" << endl;
		oRep << MC.getReasonCompleted() << endl;
		oRep << "Final Energy: " << sys.calcEnergy() << endl;
		oRep << sys.getEnergySummary() << endl;
		oRep << "Best Energy: " << bestEnergy << endl;
		oRep << printParameters(bestParams) << endl;

		if (bestEnergy < veryBestE){
			veryBestE = bestEnergy;
			bestReplicate = rep;
		}
	
		string finalRemarks = "Initial State\n";
		finalRemarks += printParameters(initialParams);
		finalRemarks += "Final State\n";
		finalRemarks += printParameters(bestParams);
		finalRemarks += sys.getEnergySummary();
	
		string finalpdb = opt.outputDir + "/" + repString + ".pdb";
		sys.writePdb(finalpdb, finalRemarks);

		/******************************************************************************
		 * ==== Print out useful statistics ====
		 ******************************************************************************/
		double finalXShift = bestParams[3];
		double finalCrossingAngle = bestParams[2];
		double finalZShift = bestParams[1];
		double finalAxialRotation = bestParams[0];
		double dockEnergy = sys.calcEnergy();
		double totalEnergy = dockEnergy-monomerE;
		vector<double> csvOutputs{totalEnergy,dockEnergy,monomerE,finalXShift,finalCrossingAngle,finalAxialRotation,finalZShift};
		csvOut << opt.chainSeq[0] << ',' << repString << ',';
		cout << opt.chainSeq[0] << ',' << repString << ',';
		for (uint i=0; i<csvOutputs.size(); i++){
			if (i < csvOutputs.size()-1){
				csvOut << csvOutputs[i] << ',';
				cout << csvOutputs[i] << ',';
			} else {
				csvOut << csvOutputs[i] << endl;
				cout << csvOutputs[i] << endl;
			}
		}
		
		/******************************************************************************
		 * ==== Print out useful statistics ====
		 ******************************************************************************/
		//TODO Print deltaEs
		if(opt.verbose){
			lout << "Number of moves: " << moveCounter << endl;
			lout << "Number of repacks: " << repackCounter << endl;
			lout << "Number of accepts: " << acceptCounter << endl;
			lout << "Number of rejects: " << rejectCounter << endl;
			if (acceptCounter != 0){
				lout << "Average gap size: " << avgGap/acceptCounter << endl;
			}
			lout << "Max gap between accepts: " << maxGap << endl;
		}

		oRep << "Number of moves: " << moveCounter << endl;
		oRep << "Number of repacks: " << repackCounter << endl;
		oRep << "Number of accepts: " << acceptCounter << endl;
		oRep << "Number of rejects: " << rejectCounter << endl;
		if (acceptCounter != 0){
			oRep << "Average gap size: " << avgGap/acceptCounter << endl;
		}
		oRep << "Max gap between accepts: " << maxGap << endl;
		eOut << repString << "\t" << bestEnergy << endl;		
	}
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	if(opt.verbose){lout << endl << "Total Time: "  << diffTime << " seconds" << endl;}

	// Make Output File
	string outFile = opt.outputDir + "/topModel.out";
	ofstream oOut;
	oOut.open(outFile.c_str());
	if(!oOut.is_open()) {
		cerr << "Unable to open " << outFile << endl;
		exit(0);
	}
	cout << opt.rerunConf << endl;
	oOut << "Best Replicate: " << bestReplicate;
	oOut << "\tEnergy: " << veryBestE << endl;
	oOut << endl << "Total Time: "  << diffTime << " seconds" << endl;
	oOut.close();
	cout << endl << "Total Time: "  << diffTime << " seconds" << endl;

	//zip the folder
	if (opt.zipOutputFiles){
		stringstream tmp;
		tmp << setw(7) << setfill('0') << opt.replicateStartNum;
		string repStartString = tmp.str();
		string tarfile = opt.outputDir + "/" + repStartString + ".tar.gz";
		string zipcommand = "tar -czf " + tarfile + " ./*.out ./*pdb ./*.list";
		system(zipcommand.c_str());
		system("rm *.out");
		system("rm *.pdb");
		system("rm *.list");
	}

	exit(0);
}





/******************************************
 *  
 *  ======= FUNCTIONS =======
 *
 ******************************************/
vector < vector<Atom*> > getInterHelicalHbonds(EnergySet* & _ESet) {
	unsigned int numHbonds = 0;
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];

	vector<string> hbonds;
	vector < vector <Atom*> > hbondapv;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			//cout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << " " << e << endl;
			hbonds.push_back(atoms[0]->getAtomOfIdentityId() + ":" + atoms[2]->getAtomOfIdentityId());
			hbondapv.push_back(atoms);
			numHbonds++;
		}
	}
	// Why are we doing this?
	//_ESet->setAllTermsActive();
	//_ESet->setTermActive("CHARMM_ELEC", false);
	return hbondapv;
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB, Transforms &_trans) {
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
	_trans.rotate(_apvB, m);

}

// Separate chains by 500 Angstroms and repack side chains
// Energy of this system is monomer energy
double calculateMonomerE (System &_sys, Transforms &_trans, SelfPairManager &_spm) {
	for (unsigned int i = 0; i < _sys.chainSize(); i++) {
		double translation = i*500.0;
		_trans.Xtranslate(_sys.getChain(i).getAtomPointers(), translation);
	}
	_spm.calculateEnergies();
	_spm.runGreedyOptimizer(50);
	double monomerE = _spm.getMinBound()[0];
	return monomerE;
}

std::vector < std::vector < bool > > getActiveMask (System &_sys) {
	_sys.updateVariablePositions();
	std::vector <unsigned int> residueState;
	std::vector < std::vector<unsigned int> > resRots(_sys.getMasterPositions().size());
	std::vector < std::vector<bool> > resMask(_sys.getMasterPositions().size());
	//Initialize residue state at the current active identity for each position
	for (unsigned int i = 0; i < _sys.getMasterPositions().size(); i++) {
		Position &pos = _sys.getPosition(_sys.getMasterPositions()[i]);
		unsigned int activeRes = pos.getActiveIdentity();
		residueState.push_back(activeRes);
		
		resRots[i] = std::vector<unsigned int> (pos.identitySize());
		for (unsigned int j = 0; j < pos.identitySize(); j++) {
			resRots[i][j] = pos.getTotalNumberOfRotamers(j);
		}
	}

	for (unsigned int i = 0; i < residueState.size(); i++) {
		unsigned int activeResidue = residueState[i];
		if (activeResidue >= resRots[i].size()) {
			cerr << "ERROR: the current residue number exceeds the number of residues for position " << i << endl;
			exit(100);
		}
		for (unsigned int j = 0; j < resRots[i].size(); j++) {
			if (j==activeResidue) {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(true);
				}
			} else {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(false);
				}
			}
		}
	
		//Sanity check for presence of true rotamers

		bool trueRots = false;
		for (unsigned int j = 0; j < resMask[i].size(); j++) {
			if (resMask[i][j]) {
				trueRots = true;
			}
		}
		if (!trueRots) {
			cerr << "ERROR AT POSITION: " << i << endl;
			cerr << "Current Residue: " << activeResidue << endl;
			cerr << "resRots at this position: " << endl;
			for (uint k = 0; k < resRots[i].size(); k++) {
				cerr << resRots[i][k] << " ";
			}
			cerr << endl;
			cerr << "resMask at this position: " << endl;
			for (uint k = 0; k < resMask[i].size(); k++) {
				cerr << resMask[i][k] << " ";
			}
			cerr << endl;	
			exit(9123);
		}
	}
	return resMask;
}

AtomPointerVector calculateHelicalAxes (AtomPointerVector &_helix) {
	// Uses principal component analysis to determine 
	// axes for an ideal alpha helix
	// Origin: center of mass
	// Z-axis: principal axis 1, extends down length of helix
	// X-axis: principal axis closest to middle alpha carbon
	// Y-axis: the other axis

	PrincipleComponentAnalysis pca;
	pca.computePrincipleComponents(_helix);
	
	Atom bestFitOrigin ("ORI", pca.getLines()[0].getCenter());
	Atom bestFitZ ("ZAX", (pca.getLines()[0].getDirection()+bestFitOrigin.getCoor()));
	Atom bestFitX ("YAX", (pca.getLines()[1].getDirection()+bestFitOrigin.getCoor()));
	Atom bestFitY ("XAX", (pca.getLines()[2].getDirection()+bestFitOrigin.getCoor()));

	AtomPointerVector bestFit;
	bestFit.push_back(new Atom(bestFitOrigin));
	bestFit.push_back(new Atom(bestFitX));
	bestFit.push_back(new Atom(bestFitY));
	bestFit.push_back(new Atom(bestFitZ));
	
	return bestFit;

}

bool duplicateHelix (System &_sys, int _helices=1) {
	AtomBondBuilder Abb;
	for (unsigned int i = 1; i < _helices; i++) {
		_sys.duplicateChain("A");
		Abb.buildConnections(_sys.getChain(i).getAtomPointers());
	}
	
	return 1;
}

void deleteTerminalHydrogenBondInteractions(System &_sys, ofstream &_lout, bool &_verbose) {
	// look at all hbond interactions and remove those
	// remove any interaction that has a donor or acceptor from residues 1 2 3 and n-2 n-1 n on each chain that are not part of the TM

	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize(); // number of positions in chain A = chain B
	// compute the extensions at the beginning and the end
	//int frontExt = _opt.tmStart - _opt.startResNum;
	//int endExt = _opt.endResNum - _opt.tmEnd;
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain& thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int j = 0; j < 4; j++) {
			atoms += positions[j]->getAtomPointers(); 
			if (_verbose){ _lout << "Removing Hbonds from " << positions[j]->getPositionId()  << endl;}
			atoms += positions[positions.size() - 1 - j]->getAtomPointers(); 
			if (_verbose){ _lout << "Removing Hbonds from " << positions[positions.size() - 1 - j]->getPositionId()  << endl;}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");

}

string shiftHelix (AtomPointerVector &_helix, RandomNumberGenerator &_RNG, Transforms &_trans, bool _shiftRotCenter) {
	string moveMade;
	AtomPointerVector helicalAxis = calculateHelicalAxes(_helix);
	double transTemp = 1.0;
	double rotTemp = 10.0;

	double Xshift = getBoxMuller(_RNG) * transTemp;
	double Yshift = getBoxMuller(_RNG) * transTemp;
	double Zshift = getBoxMuller(_RNG) * transTemp;

	double XAngle = getBoxMuller(_RNG) * rotTemp;
	double YAngle = getBoxMuller(_RNG) * rotTemp;
	double ZAngle = getBoxMuller(_RNG) * rotTemp;


	CartesianPoint shiftedPoint (0.0, 0.0, 0.0);
	if (_shiftRotCenter) {
		//Pick random point along the Z-axis to perform the X and Y rotations
		CartesianPoint shiftedRotSlope = helicalAxis[3]->getCoor() - helicalAxis[0]->getCoor();
		
		//An alpha helix rises about 1.5 angstroms per residue
		//shiftedPoint picks any point on the helical axis roughly bounded by the helix itself
		int positionSize = _helix.subdivideByChainAndPosition()[0].size();	
		double shiftLength = positionSize/2*1.5*_RNG.getRandomDouble();
		CartesianPoint shiftedPoint = shiftedRotSlope * shiftLength;
	}
	CartesianPoint shiftedCenter = helicalAxis[0]->getCoor() + shiftedPoint;
	_trans.translate(helicalAxis, shiftedPoint);
	
	//Either _translate a helix or rotate it
	int move = _RNG.getRandomInt(0,5);
	if (move == 0) {
		_trans.Xtranslate(_helix, Xshift);
		moveMade += "translated X: " + MslTools::doubleToString(Xshift);
	} else if (move==1) {
		_trans.Ytranslate(_helix, Yshift);
		moveMade += "translated Y: " + MslTools::doubleToString(Yshift);			
	} else if(move==2) {
		_trans.Ztranslate(_helix, Zshift);
		moveMade += "translated Z: " + MslTools::doubleToString(Zshift);			
	}else if (move==3) {
		_trans.rotate(_helix, XAngle, helicalAxis[1]->getCoor(), helicalAxis[0]->getCoor());
		moveMade += "Rotated X: " + MslTools::doubleToString(XAngle) + " about point " + shiftedCenter.toString();
		
	} else if (move==4) {
		_trans.rotate(_helix, YAngle, helicalAxis[2]->getCoor(), helicalAxis[0]->getCoor());
		moveMade += "Rotated Y: " + MslTools::doubleToString(YAngle) + " about point " + shiftedCenter.toString();
	} else if (move==5) {
		_trans.rotate(_helix, ZAngle, helicalAxis[3]->getCoor(), helicalAxis[0]->getCoor());
		moveMade += "Rotated Z: " + MslTools::doubleToString(ZAngle) + " about point " + shiftedCenter.toString();
	}
	return moveMade;
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

//TODO Force at least one parameter to change
void adjustParameters (RandomNumberGenerator &_RNG, std::vector<double> &_params, std::vector<double> &_steps, double _paramAdjustment, bool _limitCross, double _limitCrossValue, ofstream &_lout) {
	for (unsigned int i = 0; i < _params.size(); i++) {
		if (_steps[i] == 0) {
			continue;
		}
		if (_RNG.getRandomDouble() < _paramAdjustment) {
			if (i == 2 && _limitCross){
				bool withinRange = false;
				double startValue = _params[i];
				double bottomRange = -1 * _limitCrossValue;
				while(!withinRange){
					double xChange = getRandomNormal(_RNG, startValue, _steps[i]);
					if (xChange >= bottomRange && xChange <= _limitCrossValue){
						_params[i] = xChange;
						withinRange = true;
					}
				}
			} else {
				_params[i]= getRandomNormal(_RNG, _params[i], _steps[i]);
			}
		}
	}
	//return _params;
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

/*string transformBundle (System &_sys, Transforms &_trans, std::vector<int> &_numChains, std::vector<double> &_params) {

	unsigned int chainCounter = 0;
	//For each bundle:
	

	for (unsigned int i = 0; i < _numChains.size(); i++) {
		double angle = 0.0;
		//transform the first chain in the bundle
		AtomPointerVector bundleChain = _sys.getChain(chainCounter).getAtomPointers();
		transformHelix(_params[0], _params[1], _params[2], _params[3], bundleChain, _trans);

		//apply rotational symmetry to the rest of the chains in the bundle
		chainCounter++;
		for (unsigned int j = 1; j < _numChains[i]; j++) {

			angle += 360.0/_numChains[i];			
			//Set coordinates of buddy chains to the first chain
			AtomPointerVector chainPrime = _sys.getChain(chainCounter).getAtomPointers();
			for (uint k = 0; k < bundleChain.size(); k++) {
				chainPrime[k]->copyAllCoor(*bundleChain[k]);
			}

			//Apply rotational symmetry to the buddy chain
			_trans.Zrotate(chainPrime, angle);
			chainCounter++;
		}
	}
	string paramString = printParameters(_params);
	return paramString;
}*/

//hetero version
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

string transformBundle (System &_sys, Transforms &_trans, int _numChains, std::vector<double> &_params) {
	//transform the chains in the bundle
	double angle = 0.0;
	for (uint i=0; i<_numChains; i++){
		AtomPointerVector bundleChain = _sys.getChain(i).getAtomPointers();
		//adjust the helix geometry
		transformHelix(_params[0], _params[1], _params[2], _params[3], bundleChain, _trans);
		//give z symmetry to the other chains for each transformation
		if (i>0){
			angle += 360.0/_numChains;			
			AtomPointerVector chainPrime = _sys.getChain(i).getAtomPointers();
			//Apply rotational symmetry to the buddy chain
			_trans.Zrotate(chainPrime, angle);
		}
	}
	string paramString = printParameters(_params);
	return paramString;
}

string transformBundleRandomly (System &_sys, RandomNumberGenerator &_RNG, Transforms &_trans, std::vector<double> &_params, std::vector<double> &_steps, double _paramAdjustment, bool _limitCross, double _limitCrossValue, int _numChains, ofstream &_out) {
	//transform the chains in the bundle
	double angle = 0.0;
	for (uint i=0; i<_numChains; i++){
		AtomPointerVector bundleChain = _sys.getChain(i).getAtomPointers();
		//Randomly change a parameter for transformation
		adjustParameters(_RNG, _params, _steps, _paramAdjustment, _limitCross, _limitCrossValue, _out);
		//adjust the helix geometry using the adjusted parameters
		transformHelix(_params[0], _params[1], _params[2], _params[3], bundleChain, _trans);
		//give z symmetry to the other chains for each transformation
		if (i>0){
			angle += 360.0/_numChains;			
			AtomPointerVector chainPrime = _sys.getChain(i).getAtomPointers();
			//Apply rotational symmetry to the buddy chain
			_trans.Zrotate(chainPrime, angle);
		}
	}
	string paramString = printParameters(_params);
	return paramString;
}

void changeRandomRotamer (RandomNumberGenerator &_RNG, std::vector< std::vector <unsigned int> > &_rotIndex, std::vector< unsigned int> &_activeRots) {
	unsigned int index = _RNG.getRandomInt(0, _rotIndex.size()-1);
	unsigned int pos = _rotIndex[index][0];
	unsigned int rot = _rotIndex[index][1];
	_activeRots.at(pos) = rot;
}

//Function to turn off all but one interaction in a set
//Useful when distinguishing restraints between two chains
bool maskInteractions (std::vector<SigmoidInteraction*> _interactions) {
	//First, unmask all interactions
	std::vector<double> energies;
	for (uint i = 0; i < _interactions.size(); i++) {
		_interactions[i]->setMask(false);
		double energy = _interactions[i]->getEnergy();
		energies.push_back(energy);
	}
	int mindex = distance(energies.begin(), min_element(energies.begin(), energies.end()));
	//Mask all interactions but the lowest energy one
	
	for (uint i = 0; i < _interactions.size(); i++) {
		_interactions[i]->setMask(true);
	}
	_interactions[mindex]->setMask(false);

	return true;
	
}

//Function to initialize placement of chains
/*void createHelixBackbone (std::vector<string> &_chainSeqs, std::vector<int> &_numChains, std::vector<int> &_chainStartNum, AtomPointerVector &_apv, HelixGenerator &_hg) {
	int chainCounter = 0;
	//For each set of chains
	for (unsigned int i = 0; i < _chainSeqs.size(); i++) {
		int chainLength = _chainSeqs[i].length();
		//For the number of chains within the set
		for (unsigned int c = 0; c < _numChains[i]; c++) {

			//create a helix and assign its sequence
			AtomPointerVector helix;
			_hg.generateHelix(helix, chainLength, true, true);

			//Convert APV chain IDs and start numbers to match full polymer sequence
			for (unsigned int r = 0; r < helix.size(); r++) {
				string oldChain = helix[r]->getChainId();
				char newChain = int(oldChain[0]) + (chainCounter);
				string chainID;
				chainID += newChain;
				helix[r]->setChainId(chainID);
				//Figure out how to do chain renumbering from this... 
				helix[r]->setResidueNumber(helix[r]->getResidueNumber() + _chainStartNum[i]);

				Atom *tmpAtom = new Atom(helix(r));
				_apv.push_back(tmpAtom);
				tmpAtom = NULL;
			}
			chainCounter++;
		}
		//apply rotational symmetry to the set
	}
	
}*/

//Hetero version of create helix backbone
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

/*string convertToFullPolymerSequence(std::vector<string> &_seqs, std::vector<int> &_numChains, std::vector<int> &_chainStartNum) {
	string chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	string fullps = "";
	int chainCounter = 0;
	for (unsigned int i = 0; i < _seqs.size(); i++) {
		string threeLetter = convertToThreeLetterSequence(_seqs[i]);

		for (unsigned int j = 0; j < _numChains[i]; j++) {
			string chainID(1, chainIDs[chainCounter]);
			string startPos = "{" + MslTools::intToString(_chainStartNum[i]) + "}";
			fullps += chainID + ": " + startPos + " " + threeLetter + "\n";
			chainCounter++;
		}
	}
	return fullps;
}*/

//hetero version
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



/******************************************
 *  
 *  =======  OPTIONS FUNCTIONS =======
 *
 ******************************************/
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

	//Required
	opt.required.push_back("chainSeq");
	opt.required.push_back("chainStartNum");
	opt.allowed.push_back("startT");
	opt.allowed.push_back("endT");
	opt.allowed.push_back("mcCycles");

	//Geometric Parameterss
	opt.allowed.push_back("axialRotStart");
	opt.allowed.push_back("axialRotEnd");
	opt.allowed.push_back("zShiftStart");
	opt.allowed.push_back("zShiftEnd");
	opt.allowed.push_back("crossAngleStart");
	opt.allowed.push_back("crossAngleEnd");
	opt.allowed.push_back("xShiftStart");
	opt.allowed.push_back("xShiftEnd");
	opt.allowed.push_back("axialRot_step");
	opt.allowed.push_back("zShift_step");
	opt.allowed.push_back("crossingAngle_step");
	opt.allowed.push_back("xShift_step");
	opt.allowed.push_back("axialRot");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("xShift");

	//Output Parameters
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("numReplicates");
	opt.allowed.push_back("replicateStartNum");
	opt.allowed.push_back("printIntermediatePdbs");
	opt.allowed.push_back("zipOutputFiles");
	opt.allowed.push_back("output");
	opt.allowed.push_back("outputName");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("rotLevel");
	opt.allowed.push_back("numTrajectoryModels");

	//Monte Carlo Parameters
	opt.allowed.push_back("mcShape");
	opt.allowed.push_back("mcMaxRejects");
	opt.allowed.push_back("convergedSteps");
	opt.allowed.push_back("convergedE");
	opt.allowed.push_back("probRepack");
	opt.allowed.push_back("shiftRotCenter");
	opt.allowed.push_back("numParamMoves");
	opt.allowed.push_back("numRotMoves");
	opt.allowed.push_back("seed");
	                                         
	opt.allowed.push_back("paramAdjustment");
	opt.allowed.push_back("limitCross");
	opt.allowed.push_back("limitCrossValue");
	opt.allowed.push_back("xShiftFirst");

	//Energy Parameters
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("vdwRadius");
	opt.allowed.push_back("ctonnb");
	opt.allowed.push_back("ctofnb");
	opt.allowed.push_back("cutnb");
	opt.allowed.push_back("chainRotLevels");
	opt.allowed.push_back("deleteTerminalHbonds");

	//Restraint Parameters
	opt.allowed.push_back("restraints");
	opt.allowed.push_back("restraintWeights");
	opt.allowed.push_back("restraintSlopes");
	opt.allowed.push_back("restraintDists");
	opt.allowed.push_back("restraintIntercepts");
	opt.allowed.push_back("restraintGroups");
	opt.allowed.push_back("defaultRestraintWeight");
	opt.allowed.push_back("defaultRestraintDist");
	opt.allowed.push_back("defaultRestraintSlope");
	opt.allowed.push_back("defaultRestraintIntercept");

	//Parameter Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("bbqFile");
	opt.allowed.push_back("configfile");
	
	//Potentiall Redacted
	opt.allowed.push_back("numChains");

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
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/	
	//Required
	//opt.chainSeq = OP.getMultiString("chainSeq");
	opt.chainSeq = OP.getStringVector("chainSeq");//for hetero, get multiple seqs from config
	if (OP.fail()) {
		opt.errorMessages += "chainSeq not specified!\n";
		opt.errorFlag = true;
	}
	opt.numChains = OP.getInt("numChains");
	if (OP.fail()) {
		//opt.errorMessages += "numChains not specified!\n";
		//opt.errorFlag = true;
		opt.warningMessages += "numChains not specified!\n";
		opt.warningFlag = true;
		opt.numChains = 2;
	}
	//opt.chainStartNum = OP.getMultiInt("chainStartNum");
	opt.chainStartNum = OP.getIntVector("chainStartNum");//for hetero, have multiple start numbers (may not need, could just always start them at the same number; I think this is here for her to keep track of her protein start numbers)
	if (OP.fail()) {	
		opt.errorMessages += "chainStartNum not specified!\n";
		opt.errorFlag = true;
	}
	opt.startT = OP.getDouble("startT");
	if (OP.fail()) {
		opt.warningMessages += "start temp not specified, default to 2500!\n";
		opt.warningFlag = true;
		opt.startT = 2500;
	}

	opt.endT = OP.getDouble("endT");
	if (OP.fail()) {
		opt.warningMessages += "end temp not specified, default to 1!\n";
		opt.warningFlag = true;
		opt.endT = 1;
	}                                
	
	opt.mcCycles = OP.getInt("mcCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC cycles not specified!\n";
		opt.warningFlag = true;
		opt.mcCycles = 10000;
	}

	//Geometric Parameters
		//Starting parameters to be generated from these ranges
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
	
		//CHECK THAT NUMBER OF PARAMS ALL MATCH UP!
	opt.params.push_back(opt.axialRot);
	opt.params.push_back(opt.zShift);
	opt.params.push_back(opt.crossingAngle);
	opt.params.push_back(opt.xShift);

	bool equalParamSize = true;
	//if (opt.chainSeq.size() != opt.numChains.size() ) {
	//	equalParamSize = false;
	//}
	//hetero version
	if (opt.chainSeq.size() != opt.numChains ) {
		equalParamSize = false;
	}

	if (opt.chainSeq.size() != opt.chainStartNum.size() ) {
		equalParamSize = false;
	}
	
	if (opt.chainRotLevels.size() != 0 && opt.chainSeq.size() != opt.chainRotLevels.size() ) {
		equalParamSize = false;
	}

//	for (uint i = 0; i < opt.params.size(); i++) {
//		if (opt.chainSeq.size() != opt.params[i].size()) {
//			equalParamSize = false;
//		}
//	}
//	if (!equalParamSize) {
//		opt.errorMessages += "ERROR: NOT A FULL SET OF TMH PARAMS FOR EVERY BUNDLE!\n";
//		opt.errorMessages += "REQUIRED: axialRot, zShift, crossingAngle, xShift\n";
//		opt.errorFlag = true;
//	}
		//Step Sizes
	opt.axialRot_step = OP.getDouble("axialRot_step");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "axialRot_step not specified! Default to 1.0\n";
		opt.axialRot_step = 1.0;
	}

	opt.zShift_step = OP.getDouble("zShift_step");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "zShift_step not specified! Default to 1.0\n";
		opt.zShift_step = 1.0;
	}

	opt.crossingAngle_step = OP.getDouble("crossingAngle_step");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "crossingAngle_step not specified! Default to 1.0\n";
		opt.crossingAngle_step = 1.0;
	}

	opt.xShift_step = OP.getDouble("xShift_step");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "xShift_step not specified! Default to 1.0\n";
		opt.xShift_step = 1.0;
	}

	opt.steps.push_back(opt.axialRot_step);
	opt.steps.push_back(opt.zShift_step);
	opt.steps.push_back(opt.crossingAngle_step);
	opt.steps.push_back(opt.xShift_step);

	//Starting Parameters
	opt.axialRot = OP.getDouble("axialRot");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "axialRot not specified! Setting to 1234567890.\n";
		opt.axialRot = 1234567890;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "zShift not specified! Setting to 1234567890.\n";
		opt.zShift = 1234567890;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "crossingAngle not specified! Setting to 1234567890.\n";
		opt.crossingAngle = 1234567890;
	}
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "xShift not specified! Setting to 1234567890.\n";
		opt.xShift = 1234567890;
	}
	//Design parameters
	opt.designAxialRotation = OP.getDouble("designAxialRotation");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "axialRot not specified! Setting to 1234567890.\n";
		opt.designAxialRotation = 1234567890;
	}
	opt.designZShift = OP.getDouble("designZShift");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "zShift not specified! Setting to 1234567890.\n";
		opt.designZShift = 1234567890;
	}
	opt.designCrossingAngle = OP.getDouble("designCrossingAngle");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "crossingAngle not specified! Setting to 1234567890.\n";
		opt.designCrossingAngle = 1234567890;
	}
	opt.designXShift = OP.getDouble("designXShift");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "designXShift not specified! Setting to 1234567890.\n";
		opt.designXShift = 1234567890;
	}
	opt.paramAdjustment = OP.getDouble("paramAdjustment");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "paramAdjustment not specified! Default to 1.0\n";
		opt.paramAdjustment = 1.0;
	}
	opt.limitCross = OP.getBool("limitCross");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "limitCross not specified! Default to false\n";
		opt.limitCross = false;
	}
	opt.limitCrossValue = OP.getDouble("limitCrossValue");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "limitCrossValue not specified! Default to 60.0\n";
		opt.limitCrossValue = 60.0;
	}
	opt.xShiftFirst = OP.getBool("xShiftFirst");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "xShiftFirst not specified! Default to false\n";
		opt.xShiftFirst = false;
	}


	//Output Parameters
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages+= "verbose not specified, default false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.numReplicates = OP.getInt("numReplicates");
	if (OP.fail()) {
		opt.warningMessages += "numReplicates not specified! Defaulting to 1.\n";
		opt.warningFlag = true;
		opt.numReplicates = 1;
	}
	opt.replicateStartNum = OP.getInt("replicateStartNum");
	if (OP.fail()) {
		opt.warningMessages += "replicateStartNum not specified! Defaulting to 1.\n";
		opt.warningFlag = true;
		opt.replicateStartNum = 1;
	}
	opt.printIntermediatePdbs = OP.getBool("printIntermediatePdbs");
	if (OP.fail()) {
		opt.warningMessages+= "printIntermediatePdbs not specified, default false\n";
		opt.warningFlag = true;
		opt.printIntermediatePdbs = false;
	}
	opt.zipOutputFiles = OP.getBool("zipOutputFiles");
	if (OP.fail()) {
		opt.warningMessages+= "zipOutputFiles not specified, default false\n";
		opt.warningFlag = true;
		opt.zipOutputFiles = false;
	}
	opt.output = OP.getString("output");
	if (OP.fail()) {
		opt.warningMessages += "output file not specified!\n";
		opt.warningFlag = true;
		opt.output = programName + ".log";
	}
	opt.outputName = OP.getString("outputName");
	if (OP.fail()) {
		opt.warningMessages += "output dir not specified!\n";
		opt.warningFlag = true;
		opt.outputName = "./";
	}

	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "output dir not specified! using current directory\n";
		opt.warningFlag = true;
		opt.outputDir = get_current_dir_name();
	}

	opt.numTrajectoryModels = OP.getInt("numTrajectoryModels");
	if (OP.fail()) {
		opt.warningMessages += "numTrajectoryModels not specified, default to 100\n";
		opt.warningFlag = true;
		opt.numTrajectoryModels = 100;
	}

	//Monte Carlo Parameters
	opt.mcShapeString = OP.getString("mcShape");
	if (OP.fail()) {
		opt.warningMessages += "mcShape not specified, default to SIGMOIDAL\n";
		opt.warningFlag = true;
		opt.mcShapeString = "SIGMOIDAL";
	}
	//Check mcShape to make sure it will evaluate correctly
	if (opt.mcShapeString == "CONSTANT") {
		opt.mcShape = 0;
	} else if (opt.mcShapeString == "LINEAR") {
		opt.mcShape = 1;
	} else if (opt.mcShapeString == "EXPONENTIAL") {
		opt.mcShape = 2;
	} else if (opt.mcShapeString == "SIGMOIDAL") {
		opt.mcShape = 3;
	} else if (opt.mcShapeString == "SOFT") {
		opt.mcShape = 4;
	} else {
		opt.errorMessages += "mcShape not CONSTANT, LINEAR, EXPONENTIAL, SIGMOIDAL, or SOFT\n";
		opt.errorFlag = true;
	}
	opt.mcMaxRejects = OP.getInt("mcMaxRejects");
	if (OP.fail()) {
		opt.mcMaxRejects = 5;
		opt.warningMessages += "mcMaxRejects not specified, default to using 5\n";
		opt.warningFlag = true;
	}	

	opt.convergedSteps = OP.getInt("convergedSteps");
	if (OP.fail()) {
		opt.warningMessages += "Max steps without improvement not specified! Default to 250\n";
		opt.warningFlag = true;
		opt.convergedSteps = 250;
	}
	opt.convergedE = OP.getDouble("convergedE");
	if (OP.fail()) {
		opt.warningMessages += "convergedE not specified!\n";
		opt.warningFlag = true;
		opt.convergedE = 0.5;
	}	
	opt.probRepack = OP.getDouble("probRepack");
	if (OP.fail()) {
		opt.warningMessages += "probability to repack side chains not specified, default to 1.0\n";
		opt.warningFlag = true;
		opt.probRepack = 1;
	}

	opt.shiftRotCenter = OP.getBool("shiftRotCenter");
	if (OP.fail()) {
		opt.warningMessages += "shifted center of rotations not specified! Default to false\n";
		opt.warningFlag = true;
		opt.shiftRotCenter = false;
	}
	opt.numParamMoves = OP.getInt("numParamMoves");
	if (OP.fail()) {
		opt.warningMessages += "numParamMoves not specified, default to 5\n";
		opt.warningFlag = true;
		opt.numParamMoves = 5;
	}
	opt.numRotMoves = OP.getInt("numRotMoves");
	if (OP.fail()) {
		opt.warningMessages += "numRotMoves not specified, default to 0\n";
		opt.warningFlag = true;
		opt.numRotMoves = 0;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified, default to 0\n";
		opt.warningFlag = true;
		opt.seed = 0;
	}
	
	//Energy Parameters
	opt.weight_vdw = OP.getDouble("weight_vdw");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_vdw not specified, default 1.0\n";
		opt.weight_vdw = 1.0;
	}
	opt.weight_hbond = OP.getDouble("weight_hbond");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_hbond not specified, default 1.0\n";
		opt.weight_hbond = 1.0;
	}

	opt.vdwRadius = OP.getDouble("vdwRadius");
	if (OP.fail()) {
		opt.warningMessages += "vdwRadius not specified!\n";
		opt.warningFlag = true;
		opt.vdwRadius = 1.0;
	}
	opt.ctonnb = OP.getDouble("ctonnb");
	if (OP.fail()) {
		opt.ctonnb = 9.0;
		opt.warningMessages += "ctonnb not specified! Default to ";
		opt.warningMessages += MslTools::doubleToString(opt.ctonnb);
		opt.warningMessages += "\n";
		opt.warningFlag = true;
	}

	opt.ctofnb = OP.getDouble("ctofnb");
	if (OP.fail()) {
		opt.ctofnb = 12.0;
		opt.warningMessages += "ctofnb not specified! Default to ";
		opt.warningMessages += MslTools::doubleToString(opt.ctofnb);
		opt.warningMessages += "\n";
		opt.warningFlag = true;
	}
	
	opt.cutnb = OP.getDouble("cutnb");
	if (OP.fail()) {
		opt.cutnb = 16.0;
		opt.warningMessages += "cutnb not specified! Default to ";
		opt.warningMessages += MslTools::doubleToString(opt.cutnb);
		opt.warningMessages += "\n";
		opt.warningFlag = true;
	}
	opt.rotLevel = OP.getString("rotLevel");
	if (OP.fail()) {
		opt.warningMessages += "rotamer level not specified\n";
		opt.warningFlag = true;
		opt.rotLevel = "SL80.00";
	}
	opt.chainRotLevels = OP.getMultiString("chainRotLevels");
	if (OP.fail()) {
		opt.warningMessages += "chainRotLevels not specified!\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = false;
	}

	//Restraint Parameters
	opt.restraints = OP.getMultiString("restraints");
	if (OP.fail()) {
		opt.warningMessages+= "No restraints specified!\n";
		opt.warningFlag = true;
	}
	
	opt.defaultRestraintWeight = OP.getDouble("defaultRestraintWeight");
	if (OP.fail()) {
		opt.warningMessages += "No default restraint weight specified, default to 1.0\n";
		opt.warningFlag = true;
		opt.defaultRestraintWeight = 1.0;
	}

	opt.defaultRestraintDist = OP.getDouble("defaultRestraintDist");
	if (OP.fail()) {
		opt.warningMessages += "No default restraint distance specified, default to 5.0\n";
		opt.warningFlag = true;
		opt.defaultRestraintDist = 5.0;
	}
	opt.defaultRestraintSlope = OP.getDouble("defaultRestraintSlope");
	if (OP.fail()) {
		opt.warningMessages += "No default restraint slope specified, default to 1.0\n";
		opt.warningFlag = true;
		opt.defaultRestraintSlope = 1.0;
	}
	opt.defaultRestraintIntercept = OP.getDouble("defaultRestraintIntercept");
	if (OP.fail()) {
		opt.warningMessages += "No default restraint Intercept specified, default to 0.0\n";
		opt.warningFlag = true;
		opt.defaultRestraintIntercept = 0.0;
	}
	opt.restraintWeights = OP.getMultiDouble("restraintWeights");
	if (OP.fail()) {
		opt.warningMessages += "Restraint weights not specified! Default to same restraint weight\n";
		opt.warningFlag = true;
	}

	opt.restraintSlopes = OP.getMultiDouble("restraintSlopes");
	if (OP.fail()) {
		opt.warningMessages += "Restraint slopes not specified! Default to same restraint slopes\n";
		opt.warningFlag = true;
	}

	opt.restraintDists = OP.getMultiDouble("restraintDists");
	if (OP.fail()) {
		opt.warningMessages += "Restraint distances not specified! Default to same restraint distance\n";
		opt.warningFlag = true;
	}
	opt.restraintIntercepts = OP.getMultiDouble("restraintIntercepts");
	if (OP.fail()) {
		opt.warningMessages += "Restraint intercepts not specified! Default to same restraint intercept\n";
		opt.warningFlag = true;
	}

	opt.restraintGroups = OP.getMultiString("restraintGroups");
	if (OP.fail()) {
		opt.warningMessages += "Restraint Groupings not specified! Default to independent restraints\n";
		opt.warningFlag = true;
	}

	//Need to check if the number of restraints matches the number of restraint parameters
	bool restraintFlag = true;

	//Fill out empty restraint params with default values
	if (opt.restraintWeights.empty()) {
		for (int i = 0; i < opt.restraints.size(); i++) {
			opt.restraintWeights.push_back(opt.defaultRestraintWeight);
		}
	}
	if (opt.restraintSlopes.empty()) {
		for (int i = 0; i < opt.restraints.size(); i++) {
			opt.restraintSlopes.push_back(opt.defaultRestraintSlope);
		}
	}
	if (opt.restraintDists.empty()) {
		for (int i = 0; i < opt.restraints.size(); i++) {
			opt.restraintDists.push_back(opt.defaultRestraintDist);
		}
	}
	if (opt.restraintIntercepts.empty()) {
		for (int i = 0; i < opt.restraints.size(); i++) {
			opt.restraintIntercepts.push_back(opt.defaultRestraintIntercept);
		}
	}

	//Check that each restraint now has a set of parameters
	if (opt.restraints.size() != opt.restraintWeights.size()) {
		restraintFlag = false;
	}
	if (opt.restraints.size() != opt.restraintDists.size()) {
		restraintFlag = false;
	}
	if (opt.restraints.size() != opt.restraintSlopes.size()) {
		restraintFlag = false;
	}
	if (opt.restraints.size() != opt.restraintIntercepts.size()) {
		restraintFlag = false;
	}

	//Check if there's either no groups assigned to restraints or they are all assigned to different groups
	if (!opt.restraintGroups.empty() && opt.restraintGroups.size() != opt.restraints.size()) {
		restraintFlag = false;
	}

	if (restraintFlag == false) {
		opt.errorMessages += "ERROR: Number of restraints does not match number of distances or weights!\n";
		opt.errorFlag = true;
	}

	//Parameter Files
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
	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.solvFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(SYSENV.isDefined(envVar)) {
			opt.rotLibFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + ", defaulting to " + SYSENV.getEnv(envVar) + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}
	opt.hbondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hbondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "File not specified using " + opt.hbondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine File - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.bbqFile = OP.getString("bbqFile");
	if (OP.fail()) {
		string envVar = "MSL_BBQ_TABLE";
		if(SYSENV.isDefined(envVar)) {
			opt.bbqFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "bbqFile not specified using " + opt.bbqFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine bbqFile - MSL_BBQ_TABLE - not set\n"	;
			opt.errorFlag = true;
		}
	}

	// return the Options structure
	
	opt.rerunConf = OP.getConfFile();
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
//	cout << "This program runs as:" << endl;
//	cout << " % genHomoUniverse --helixPdbFile <pdbfile> --helicalAxisPdbFile <filename> --output <filename>" << endl;
//	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
//	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
//	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
//	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
//	cout << endl;
}

