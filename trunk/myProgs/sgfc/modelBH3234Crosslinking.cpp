// This program accepts a PDB file of Bax helices 234 and shifts 
// them so that the residues that form crosslinks come into closer 
// contact.  The helices undergo Monte Carlo rearrangement and the 
// energy and correspondence with crosslinking data is optimized.
//
// TODO Figure out strategy
	// Do I work with ideal helices and build off that or generate 
	// template models from poly-Ala? How do I handle 6 chains?
// *TODO Read in PDB File
// *TODO Repack/Relax PDB File
// TODO Add disulfide constraints
// TODO Add loop constraints (same function, but longer distance?)
// TODO Move helices
	// How do we move? MC? All at once? One helix at a time?
		// Translate helix w/ respect to geo. center
		// Rotate helix w/ respect to helical axis and some normal
	// Enforce symmetry?
// TODO Evaluate energetics
	// Combo of good energy and minimized constraints
	// Do I calculate "monomer" energy?
//TODO Print PDBs
//
//


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
#include "SysEnv.h"
#include "CartesianGeometry.h"
#include "MonteCarloManager.h"
#include "SpringInteraction.h"

using namespace MSL;
using namespace std;
static SysEnv ENV;


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
AtomPointerVector calculateHelicalAxes (AtomPointerVector _helixAtoms);
void transformHelix(AtomPointerVector _helix, CartesianPoint _translation);

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

/******************************************
 *  
 *  =======  BEGIN MAIN  =======
 *
 ******************************************/

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
//	Options opt = parseOptions(argc, argv, defaults);
//	if (opt.errorFlag) {
//		cerr << endl;
//		cerr << "The program terminated with errors:" << endl;
//		cerr << endl;
//		cerr << opt.errorMessages << endl;
//		cerr << endl;
//		cerr << opt.OPerrors << endl;
//
//		//printErrors(opt);
//		exit(1);
//	}
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
	RandomNumberGenerator RNG;
	RNG.setSeed(0);

	
	// Set up Charmm System Builder
	string topFile = ENV.getEnv("MSL_CHARMM_TOP");
	string parFile = ENV.getEnv("MSL_CHARMM_PAR");
	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");
	

	CSB.buildSystemFromPDB("/data05/scondon/Projects/bax_oligomer/bh3_234/bh3InGroove_noLoops_polyAla_splitHelices.pdb");
	//CSB.buildSystemFromPDB("/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/31-gly-residue-helix-chainA.pdb");
	sys.buildAllAtoms();



	HydrogenBondBuilder hb(sys,"/data01/sabs/msl_working/mslib/trunk/toppar/scwrl4hb/par_hbond_2.txt");
	if(!hb.buildInteractions(50)) {
		cerr << "Unable to build hbonds " << endl;
		exit(0);
	}
	CSB.updateNonBonded();
	EnergySet* pESet = sys.getEnergySet();
	/******************************************************************************
	 * ==== REMOVE THE INTERACTIONS WITHIN THE BACKBONE OF EACH HELIX SINCE THEY DO NOT CHANGE ====
	 ******************************************************************************/
	AtomSelection sel(sys.getAtomPointers();
	for (unsigned int i = 0; i < sys.chainSize()-1; i++) {
		string sele = "chain" + MslTools::IntToString(i);
		string selection = "chain" + MslTools::IntToString(i) + ", chain " + MslTools::IntToString(i);
		AtomPointerVector chaini = sel.select(selection);
		pESet->deleteInteractionsWithinSelection(sele);
	}

	cout << sys << endl;
	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;

	/******************************************************************************
	 *               === ADD SPRING CONSTRAINTS ===
	 ******************************************************************************/
	


	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	double bestEnergy = sys.calcEnergy();
	double currentEnergy = bestEnergy;
	MonteCarloManager MC(1000.0, 0.5, 1000, MonteCarloManager::EXPONENTIAL, 50);
	MC.setEner(bestEnergy);
	while (!MC.getComplete()) {
		sys.saveAltCoor("previousState");
		unsigned int chain = RNG.getRandomInt(0,sys.chainSize()-1);
		AtomPointerVector helix = sys.getChain(chain).getAtomPointers();
		AtomPointerVector helicalAxis = calculateHelicalAxes(helix);

		double temp = sqrt(MC.getCurrentT());

		double Xshift = getStandardNormal(RNG) * temp;
		double Yshift = getStandardNormal(RNG) * temp;
		double Zshift = getStandardNormal(RNG) * temp;
		CartesianPoint delta (Xshift, Yshift, Zshift);

		double XAngle = getStandardNormal(RNG) * temp;
		double YAngle = getStandardNormal(RNG) * temp;
		double ZAngle = getStandardNormal(RNG) * temp;


		trans.rotate(helix, XAngle, helicalAxis[1]->getCoor(), helicalAxis[0]->getCoor());
		trans.rotate(helicalAxis, XAngle, helicalAxis[1]->getCoor(), helicalAxis[0]->getCoor());
		trans.rotate(helix, YAngle, helicalAxis[2]->getCoor(), helicalAxis[0]->getCoor());
		trans.rotate(helicalAxis, YAngle, helicalAxis[2]->getCoor(), helicalAxis[0]->getCoor());
		trans.rotate(helix, ZAngle, helicalAxis[3]->getCoor(), helicalAxis[0]->getCoor());
		trans.rotate(helicalAxis, ZAngle, helicalAxis[3]->getCoor(), helicalAxis[0]->getCoor());

		cout << "Translating helix " << chain << " by " << delta << endl;		
		transformHelix(helix, delta);
		currentEnergy = sys.calcEnergy();

		if (!MC.accept(currentEnergy)) {
			cout << "state rejected. Energy : " << currentEnergy << endl;
			sys.applySavedCoor("previousState");
		} else {
			bestEnergy = currentEnergy;
			cout << "state accepted! Energy: " << currentEnergy << endl;
		}
	}
	cout << "MC Repack Complete!" << endl;
	cout << sys.getEnergySummary() << endl;
	sys.writePdb("MCtest.pdb");



	exit(1);
}

// Approximate the helical axes for translations and rotations
	// helical origin: CA geometric center
	// helical X: One normal in direction of middle CA
	// helical Y: Other normal is perpendicular to XZ plane
	// helical Z: geometric center of first N/2 residues or first helical turn, whichever is larger

AtomPointerVector calculateHelicalAxes (AtomPointerVector _helixAtoms) {

	AtomSelection sele (_helixAtoms);
	AtomPointerVector backbone = sele.select("backbone, name CA+C+N+O");
	AtomPointerVector CA = sele.select("CA, name CA");
	if (CA.size() % 2 == 0) {
		CA.pop_back();
	}
	CartesianPoint helicalOrigin = CA.getGeometricCenter();
	CartesianPoint helicalX = CA(CA.size()/2).getCoor();


	AtomPointerVector NterminalHelixAtoms;
	if (CA.size() / 2 >= 4) {
		for (AtomPointerVector::const_iterator it = CA.begin(); it != CA.begin() + CA.size()/2; ++it) {
			NterminalHelixAtoms.push_back(*it);
		}
	} else {
		for (AtomPointerVector::const_iterator it = CA.begin(); it!=CA.begin()+3; ++it) {
			NterminalHelixAtoms.push_back(*it);
		}
	}

	CartesianPoint helicalZ = NterminalHelixAtoms.getGeometricCenter();
	CartesianPoint helicalY = helicalOrigin + CartesianGeometry::normalToPlane(helicalZ, helicalOrigin, helicalX);
	Atom* hO = new Atom ("HO", helicalOrigin, "X");
	Atom* hX = new Atom ("HX", helicalX, "X");
	Atom* hY = new Atom ("HY", helicalY, "X");
	Atom* hZ = new Atom ("HZ", helicalZ, "X");

	//cout << "XY PLANE: " << CartesianGeometry::angle(helicalX, helicalOrigin, helicalY) << endl;
	//cout << "XZ PLANE: " << CartesianGeometry::angle(helicalX, helicalOrigin, helicalZ) << endl;
	//cout << "YZ PLANE: " << CartesianGeometry::angle(helicalY, helicalOrigin, helicalZ) << endl;

	AtomPointerVector helicalAxes;
	helicalAxes.push_back(hO);
	helicalAxes.push_back(hX);
	helicalAxes.push_back(hY);
	helicalAxes.push_back(hZ);
	return helicalAxes;
}


void transformHelix(AtomPointerVector _helix, CartesianPoint _translation) {
	AtomPointerVector helicalAxes = calculateHelicalAxes (_helix);

	Transforms trans;
	trans.setTransformAllCoors(true);
	trans.setNaturalMovements(true);

	trans.translate(_helix, _translation);
}
