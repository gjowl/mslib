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
#include "SasaCalculator.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "geomRepack";
string programDescription = "This pairs with seqDesign: takes designed sequences at their given geometry and does a local geometric repack on them";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "19 August 2021";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime, spmStart, spmEnd;
double diffTime, spmTime;

/******************************************
 *        ======= OPTIONS =======
 ******************************************/
struct Options{//TODO: divide these into necessary and non-necessary
	// input files
	string topFile;
	string parFile;
	string solvFile;
	string hbondFile;
	string rotLibFile;
	string designPdb;
	string designCrd;
	string outputDir;

	// Geometry
	double xShift;
	double crossingAngle;
	double axialRotation;
	double zShift;
	int thread;

	// tm
	int tmStart;
	int tmEnd;

	// Design parameters
	string sequence;
	vector<uint> state;
	vector<int> posVector;
	vector<double> burialVector;

	// load rotamers from SASA values (from sgfc)
	bool keepOriginalRotamer;
	std::vector<string> sasaRepackLevel;

	// Monte Carlo Parameters
	int MCCycles;
	int MCMaxRejects;
	double MCStartTemp;
	double MCEndTemp;
	int MCCurve;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	double monomer;
	double monoVdw;
	double monoHbond;
	double monoIMM1;
	
	// input previous dimerEnergy
	double prevDimer;
	double prevDimerVdw;
	double prevDimerHbond;
	double prevDimerIMM1;
	
	// Shift Size
	double deltaX; 
	double deltaCross;
	double deltaAx;
	double deltaZ;

	// Other options
	bool verbose;
	int greedyCycles;
	bool dockHelices;
	string designNumber;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help //figure out how to add stuff to help

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

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
	string runNumber;
	bool useIMM1;
};

//TODO: edit this
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}
/******************************************
 *        ======= FUNCTIONS =======
 ******************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults);
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum);
void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt);
vector<pair<int,double>> calculateResidueBurial (System &_sys);
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<pair<int,double>> &_resiBurial);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType);
void repackSideChains(SelfPairManager & _spm, int _greedyCycles);
double getStandardNormal(RandomNumberGenerator& RNG);
void localMC(System &_sys, Options &_opt, PolymerSequence &_PL, RandomNumberGenerator &_RNG, map<string,map<string,double>> &_seqEnergyMap, ofstream &_sout, ofstream &_err, PDBWriter &_writer){
	System sysDimer;
	CharmmSystemBuilder CSBDimer(sysDimer,_opt.topFile,_opt.parFile,_opt.solvFile);
	CSBDimer.setBuildTerm("CHARMM_ELEC", false);
	CSBDimer.setBuildTerm("CHARMM_ANGL", false);
	CSBDimer.setBuildTerm("CHARMM_BOND", false);
	CSBDimer.setBuildTerm("CHARMM_DIHE", false);
	CSBDimer.setBuildTerm("CHARMM_IMPR", false);
	CSBDimer.setBuildTerm("CHARMM_U-BR", false);
	CSBDimer.setBuildTerm("CHARMM_IMM1REF", true);
	CSBDimer.setBuildTerm("CHARMM_IMM1", true);
	
	CSBDimer.setSolvent("MEMBRANE");
	CSBDimer.setIMM1Params(15, 10);
	
	CSBDimer.setBuildNonBondedInteractions(false);
	if(!CSBDimer.buildSystem(_PL)) {
		cerr << "Unable to build system from " << _PL << endl;
		exit(0);
	}
	
	Chain & chainA = sysDimer.getChain("A");
	Chain & chainB = sysDimer.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sysDimer.assignCoordinates(_sys.getAtomPointers(),false);
	sysDimer.buildAllAtoms();
	
	SystemRotamerLoader sysRot(_sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	HydrogenBondBuilder hb(sysDimer, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct
	
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* EsetDimer = sysDimer.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	EsetDimer->setAllTermsActive();
	EsetDimer->setTermActive("CHARMM_ELEC", false);
	EsetDimer->setTermActive("CHARMM_ANGL", false);
	EsetDimer->setTermActive("CHARMM_BOND", false);
	EsetDimer->setTermActive("CHARMM_DIHE", false);
	EsetDimer->setTermActive("CHARMM_IMPR", false);
	EsetDimer->setTermActive("CHARMM_U-BR", false);
	EsetDimer->setTermActive("CHARMM_IMM1", true);
	EsetDimer->setTermActive("CHARMM_IMM1REF", true);
	EsetDimer->setTermActive("SCWRL4_HBOND", true);
	
	/******************************************************************************
	 *             === SETUP ENERGY SET FOR MONOMER COMPARISON ===
	 ******************************************************************************/
	EsetDimer->setWeight("CHARMM_VDW", 1);
	EsetDimer->setWeight("SCWRL4_HBOND", 1);
	EsetDimer->setWeight("CHARMM_IMM1", 1);
	EsetDimer->setWeight("CHARMM_IMM1REF", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sysDimer,_opt);
	
	/******************************************************************************
	 *                     === SET UP TRANSFORMS OBJECT ===
	 ******************************************************************************/
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	vector<pair <int, double> > resiBurial = calculateResidueBurial(sysDimer);
	sort(resiBurial.begin(), resiBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});
	loadRotamersBySASABurial(sysDimer, sysRot, _opt, resiBurial);

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";
	
	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	//CSBDimer.updateNonBonded(10,12,50);
	CSBDimer.updateNonBonded();
	sysDimer.buildAllAtoms();

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spmDimer;
	spmDimer.seed(_RNG.getSeed());
	spmDimer.setSystem(&sysDimer);
	spmDimer.setVerbose(false);
	spmDimer.getMinStates()[0];
	spmDimer.updateWeights();
	spmDimer.setOnTheFly(true);
	spmDimer.saveEnergiesByTerm(true);
	spmDimer.calculateEnergies();
	
	_sout << "***CALCULATE DIMER ENERGY WITH LOCAL REPACK***" << endl << endl;
	
	sysDimer.saveAltCoor("start");
	helicalAxis.saveAltCoor("start");
	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	//Switch to current sequence
	string sequence = _opt.sequence;
	vector<uint> stateVec = _opt.state;
	sysDimer.setActiveRotamers(stateVec);
	//vector<vector<bool>> mask = getActiveMask(sysDimer);
	cout << spmDimer.getSummary(stateVec) << endl;
	double bestEnergy = spmDimer.getStateEnergy(stateVec);
	
	//Repack dimer
	spmDimer.runGreedyOptimizer(_opt.greedyCycles);
	vector<uint> greedyStateVec = spmDimer.getMinStates()[0];
	cout << spmDimer.getSummary(greedyStateVec) << endl;
	double greedyEnergy = spmDimer.getStateEnergy(greedyStateVec);
	if (greedyEnergy < bestEnergy){
		bestEnergy = greedyEnergy;
		stateVec = greedyStateVec;
	}

	sysDimer.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	cout << spmDimer.getSummary(stateVec) << endl;
	cout << spmDimer.getMinBound()[0] << endl;
	
	double monomerEnergy = _seqEnergyMap[sequence]["Monomer"];
	
	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	double xShift = _opt.xShift;
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;
	
	_sout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	_sout << " xShift: " << xShift << endl << endl;
	
	sysDimer.setActiveRotamers(stateVec);
	double currentEnergy = bestEnergy;
	//double currentEnergy = spmDimer.getMinBound()[0];
	//double bestEnergy = currentEnergy;

	//Redacted for now as of 06_05_2021: docks the helices together
	if (_opt.dockHelices){
		/******************************************************************************
		 *                     === X SHIFT REPACKS ===
		 ******************************************************************************/
		double savedXShift = xShift;
		double previousEnergy = monomerEnergy;
		double deltaXShift = -0.1;
		double globalLowestE = monomerEnergy;
		double xShiftEnd = xShift-0.6;//TODO: probably should change this to an option, or xShift-0.5 or something of the like
	
		
		//if (crossingAngle < 0){
		//	xShiftEnd = 6.4;
		//}
	
		//Redacted for now as of 06_05_2021: docks the helices together
		while (xShift >= xShiftEnd) {
		
			xShift += deltaXShift;
	
			// Move the helix
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );
			
			// Run Optimization
			repackSideChains(spmDimer, _opt.greedyCycles);
		
			vector<unsigned int> MCOFinal;
			MCOFinal = spmDimer.getMinStates()[0];
			sysDimer.setActiveRotamers(MCOFinal);
		
			cout << spmDimer.getSummary(MCOFinal) << endl;
			currentEnergy = spmDimer.getMinBound()[0];
		
			if (currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				savedXShift = xShift;
				sysDimer.saveAltCoor("savedBestState");
				helicalAxis.saveAltCoor("BestAxis");
			}
		
			//_out << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
			cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
	
			// If energy increase twice in a row, and it is above the global lowest energy, quit
			if (currentEnergy < globalLowestE) {
				globalLowestE = currentEnergy;
			}
			if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
				//_out << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
				break;	
			}
			else {
				previousEnergy = currentEnergy;
			}
		} // redacted on 05_27_2021: decided to not do this; if a seqs[i] has to move too much it may not be a good fit for the backbone (and we may be docking later anyways)
		//cout << "Best Energy at x shift: " << bestEnergy-_monomerEnergy << " at " << savedXShift << endl; 
		//xShift = savedXShift;
	}

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	
	if (_opt.verbose){
		cout << "====================================" << endl;
		cout << "Performing Local Monte Carlo Repacks" << endl;
		cout << "====================================" << endl;
	}

	vector<unsigned int> MCOBest = spmDimer.getMinStates()[0];
	
	if (_opt.MCCycles > 0) {
		//MonteCarloManager MCMngr(1000.0, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
		//MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
		MonteCarloManager MCMngr(0, 0, 50, 2, 1);//TODO: try this out to see if it's faster and still converges (trying to see if I'm running too many cycles
	
		MCMngr.setEner(bestEnergy);
		
		unsigned int counter = 0;
		while(!MCMngr.getComplete()) {
	
			sysDimer.applySavedCoor("savedBestState");
			helicalAxis.applySavedCoor("BestAxis");
	
			int moveToPreform = _RNG.getRandomInt(3);
	
			double deltaXShift = 0.0;
			double deltaZShift = 0.0;
			double deltaCrossingAngle = 0.0;
			double deltaAxialRotation = 0.0; 
	
			//======================================
			//====== Z Shift (Crossing Point) ======
			//======================================
			if (moveToPreform == 0) {
				//deltaZShift = getStandardNormal(RNG1) * 0.1;
				deltaZShift = getStandardNormal(_RNG) * _opt.deltaZ;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
			} else if (moveToPreform == 1) {
			//===========================
			//===== Axial Rotation ======
			//===========================
				//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
				deltaAxialRotation = getStandardNormal(_RNG) * _opt.deltaAx;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
			} else if (moveToPreform == 2) {
			//==================================
			//====== Local Crossing Angle ======
			//==================================
				//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
				deltaCrossingAngle = getStandardNormal(_RNG) * _opt.deltaCross;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
			} else if (moveToPreform == 3) {
			//==============================================
			//====== X shift (Interhelical Distance) =======
			//==============================================
				//deltaXShift = getStandardNormal(_RNG1) * 0.1;
				deltaXShift = getStandardNormal(_RNG) * _opt.deltaX;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
			}
	
			// Calculate energy and compare to previous best
			spmDimer.calculateEnergies();
			spmDimer.runGreedyOptimizer(_opt.greedyCycles);
			vector<unsigned int> MCOFinal = spmDimer.getMinStates()[0];
			currentEnergy = spmDimer.getStateEnergy(MCOFinal);
	
			// Run repack every N steps
			//if (counter % 5 == 0) {
			//	spmDimer.calculateEnergies();
			//	spmDimer.runGreedyOptimizer(_opt.greedyCycles, mask);
	
			//	currentEnergy = spmDimer.getMinBound()[0];
			//} else {
			//	currentEnergy = sysDimer.calcEnergy();
			//}

			if (!MCMngr.accept(currentEnergy)) {
				if (_opt.verbose){
					_sout << "state rejected   energy: " << currentEnergy-monomerEnergy << endl;
					cout << "state rejected   energy: " << currentEnergy-monomerEnergy << endl;
				}
			}
			else {
				bestEnergy = currentEnergy;
				sysDimer.saveAltCoor("savedBestState");
				helicalAxis.saveAltCoor("BestAxis");
	
				xShift = xShift + deltaXShift;
				crossingAngle = crossingAngle + deltaCrossingAngle;
				axialRotation = axialRotation + deltaAxialRotation;
				zShift = zShift + deltaZShift;
				MCOBest = MCOFinal;
	
				if (_opt.verbose){
					_sout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
					cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				}
			}
		counter++;
		}
	}
	
	sysDimer.applySavedCoor("savedBestState");

	double dimerEnergy = spmDimer.getStateEnergy(MCOBest);

	double finalEnergy = dimerEnergy-monomerEnergy;

	// Print out info to the summary file
	_sout << "Geometry " << endl;
	_sout << "xShift: " << xShift << endl;
	_sout << "crossingAngle: " << crossingAngle << endl;
	_sout << "axialRotation: " << axialRotation << endl;
	_sout << "zShift: " << zShift << endl << endl;
	_sout << "Energy Summary Below" << endl;
	_sout << "Monomer Energy: " << monomerEnergy << endl;
	_sout << "Dimer Energy: " << sysDimer.calcEnergy() << endl;
	_sout << "Final Energy: " << finalEnergy << endl << endl;
	_sout << EsetDimer->getSummary() << endl << endl;
	cout << EsetDimer->getSummary() << endl << endl;
	
	// Add energies of the sequence to the energy map
	map<string,double> &energyMap = _seqEnergyMap[sequence];
	//outputEnergiesByTerm(spmDimer, MCOBest, energyMap, _opt.energyTermList, "Dimer", 1);
	//saveEnergyDifference(_opt, _seqEnergyMap, sequence);
	_seqEnergyMap[sequence]["Dimer"] = dimerEnergy;
	_seqEnergyMap[sequence]["Total"] = finalEnergy;
	
	// Write pdb into file with all sequences
	_writer.write(sysDimer.getAtomPointers(),true,false,true);
	
	// Write an individual pdb for each sequence
	PDBWriter designWriter;
	designWriter.open(_opt.outputDir + "/design_" + _opt.designNumber + ".pdb");
	designWriter.write(sysDimer.getAtomPointers(), true, false, true);
	designWriter.close();
	
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
}

void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {

	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());
	//fout << x << " " << y << " " << pt << " " << caApV.size() << endl;

	// old code
	//for(int i = 0; i < _apV.size(); i++) {
	//	CartesianPoint& pt = _apV[i]->getCoor();
	//	pt.setZ(pt.getZ() +  zShift);
	//}

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, zShift);
	_trans.translate(_apV, interDistVect);
	_trans.translate(_axis, interDistVect);
	

}

string generatePolymerSequenceFromSequence(string _sequence, int _startResNum) {
	string ps = "";
	for (uint i=0; i<_sequence.length(); i++){
		stringstream tmp;
		tmp << _sequence[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		if(resName == "HIS") {
			resName = "HSE";
		}
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
	
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

/******************************************
 *        ======= BEGIN MAIN =======
 ******************************************/
int main(int argc, char *argv[]){
//TODO: 
//	1. input new config file for a sequence from seqDesign
//	2. transform to the appropriate geometry and state
//	3. check energy to see if it's the same
//	4.  save to the same directory as already originally saved in (make sure I think of a way to keep this continuous. It would be great if I could give my code an option so that it will atuo do this or not; and that this had an option to be connected directly to seqDesgin so that it can be used by lots of people (If I want to be able to setup a server for this, which I think I do:
//		-Make it so that people choose a density type for the kde files, choose the types of AAs they want to use; choose lengths of helices, choose using baseline or not, etc
//		-The tricky part is going to be making this also editable for lots of people if they ever want to edit the code if it breaks on something I'm not thinkin about...add comments in within the next month (it would actually be pretty incredible if I actually am able to say hat I've got a webserver up and running for this by December and also that I have data on the other part too...that would be really clutch and is definitely doable
//	5. after that all works, make a way to take all the config files and generate a submit file so that it's a pretty automatic process to run the seqDesign, then run the repack after (and save in the same folder)
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%Y_%m_%d",timeinfo);
	string date(buffer);
	
	cout << date << endl;

	time(&startTime);
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);

	/******************************************************************************
	 *                       === SETUP OUTPUT FILES ===
	 ******************************************************************************/
	ofstream sout;
	ofstream err;

	string soutfile = opt.outputDir + "/repack.out";
	string errfile  = opt.outputDir + "/errors.out";

	sout.open(soutfile.c_str());
	err.open(errfile.c_str());

	sout << date << endl;
	err << date << endl;

	//TODO: add config file output
	//printOptions(opt, pout);
	
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
	

	/******************************************************************************
	 *                   === INITIALIZE POLYMER SEQUENCE ===
	 ******************************************************************************/
	string polySeq = convertToPolymerSequenceNeutralPatch(opt.sequence, opt.thread);
	//string polySeq = generatePolymerSequenceFromSequence(opt.sequence, opt.thread);
	PolymerSequence PS(polySeq);

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";
	
	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.designCrd); 
	if(!cRead.read()) {
		cerr << "Unable to read " << opt.designCrd << endl;
		exit(0);
	}
	cRead.close();
	AtomPointerVector& designApv = cRead.getAtomPointers();//*/

	System pdb;
	pdb.readPdb(opt.designPdb);//gly69 pdb file; changed from the CRD file during testing to fix a bug but both work and the bug was separate
	
	Chain & chainA1 = pdb.getChain("A");
	Chain & chainB1 = pdb.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA1 = chainA1.getAtomPointers();
	AtomPointerVector & apvChainB1 = chainB1.getAtomPointers();
	
	/******************************************************************************
	 *                      === TRANSFORM TO COORDINATES ===
	 ******************************************************************************/
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// Transformation to zShift, axialRotation, crossingAngle, and xShift
	transformation(apvChainA1, apvChainB1, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
	moveZCenterOfCAMassToOrigin(pdb.getAtomPointers(), helicalAxis.getAtomPointers(), trans);

	/******************************************************************************
	 *             === PREPARE NEW SYSTEM WITH POLYMER SEQUENCE ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);
	
	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);
	
	CSB.setBuildNonBondedInteractions(false);
	//CSB.buildSystemFromPDB(designApv);
	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from pdb" << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                          === SETUP ENERGY SET ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_IMM1", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	
	/******************************************************************************
	 *             === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,opt);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	// Assign number of rotamers by residue burial
	vector<pair <int, double> > resiBurial;
	for (uint i=0; i<opt.posVector.size(); i++){
		int pos = opt.posVector[i];
		double burial = opt.burialVector[i];
		resiBurial.push_back(pair<int,double>(pos, burial));
	}
	loadRotamersBySASABurial(sys, sysRot, opt, resiBurial);

	sys.assignCoordinates(designApv,false);
	sys.buildAllAtoms();
	//slightly convoluted, but something with the helical axis just wouldn't get correct if I didn't tranform it, but if I did it messed with my structure...so I worked around above. I don't think it should affect my results in anyway other than putting my helices in the membrane, and it's consistent with what I did for the results in my design code

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	CSB.updateNonBonded(10,12,50);
	sys.buildAllAtoms();

	RandomNumberGenerator RNG;
	RNG.setTimeBasedSeed();

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	
	cout << "***CALCULATE DIMER ENERGY WITH LOCAL REPACK***" << endl << endl;
	sys.saveAltCoor("start");
	helicalAxis.saveAltCoor("start");
	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	//Switch to current sequence
	string sequence = opt.sequence;
	double bestEnergy = opt.prevDimer;
	
	//Repack dimer
	repackSideChains(spm, 10);
	vector<uint> stateVec = spm.getMinStates()[0];
	cout << spm.getSummary(stateVec) << endl;
	double greedyEnergy = spm.getStateEnergy(stateVec);
	//sys.setActiveRotamers(stateVec);

	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");
	cout << spm.getSummary(stateVec) << endl;
	cout << spm.getMinBound()[0] << endl;
	
	double monomerEnergy = opt.monomer;

	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	double xShift = opt.xShift;
	double crossingAngle = opt.crossingAngle;
	double axialRotation = opt.axialRotation;
	double zShift = opt.zShift;
	
	cout << setiosflags(ios::fixed) << setprecision(3) << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
	cout << " xShift: " << xShift << endl << endl;
	
	double currentEnergy = bestEnergy;
	//Redacted for now as of 06_05_2021: docks the helices together
	if (opt.dockHelices){
		/******************************************************************************
		 *                     === X SHIFT REPACKS ===
		 ******************************************************************************/
		double savedXShift = xShift;
		double previousEnergy = monomerEnergy;
		double deltaXShift = -0.1;
		double globalLowestE = monomerEnergy;
		double xShiftEnd = xShift-0.6;//TODO: probably should change this to an option, or xShift-0.5 or something of the like
	
		//Redacted for now as of 06_05_2021: docks the helices together
		while (xShift >= xShiftEnd) {
		
			xShift += deltaXShift;
	
			// Move the helix
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );
			
			// Run Optimization
			repackSideChains(spm, opt.greedyCycles);
		
			vector<unsigned int> MCOFinal;
			MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
		
			cout << spm.getSummary(MCOFinal) << endl;
			currentEnergy = spm.getMinBound()[0];
		
			if (currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				savedXShift = xShift;
				sys.saveAltCoor("savedBestState");
				helicalAxis.saveAltCoor("BestAxis");
			}
		
			//_out << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
			cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
	
			// If energy increase twice in a row, and it is above the global lowest energy, quit
			if (currentEnergy < globalLowestE) {
				globalLowestE = currentEnergy;
			}
			if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
				//_out << "Energy increasing above global lowest energy... (currently " << globalLowestE-monomerEnergy << ")" << endl;
				break;	
			}
			else {
				previousEnergy = currentEnergy;
			}
		} // redacted on 05_27_2021: decided to not do this; if a seqs[i] has to move too much it may not be a good fit for the backbone (and we may be docking later anyways)
		//cout << "Best Energy at x shift: " << bestEnergy-_monomerEnergy << " at " << savedXShift << endl; 
		//xShift = savedXShift;
	}

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	
	if (opt.verbose){
		cout << "====================================" << endl;
		cout << "Performing Local Monte Carlo Repacks" << endl;
		cout << "====================================" << endl;
	}

	vector<unsigned int> MCOBest = spm.getMinStates()[0];
	
	//MonteCarloManager MCMngr(1000.0, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	//MonteCarloManager MCMngr(opt.MCStartTemp, opt.MCEndTemp, opt.MCCycles, opt.MCCurve, opt.MCMaxRejects);
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);//same amount as in monomer
	
	MCMngr.setEner(bestEnergy);
	cout << "Current Best Energy: " << bestEnergy-monomerEnergy << endl;
	
	unsigned int counter = 0;
	while(!MCMngr.getComplete()) {
	
		sys.applySavedCoor("savedBestState");
		helicalAxis.applySavedCoor("BestAxis");
	
		int moveToPreform = RNG.getRandomInt(3);
	
		double deltaXShift = 0.0;
		double deltaZShift = 0.0;
		double deltaCrossingAngle = 0.0;
		double deltaAxialRotation = 0.0; 
	
		//======================================
		//====== Z Shift (Crossing Point) ======
		//======================================
		if (moveToPreform == 0) {
			//deltaZShift = getStandardNormal(RNG1) * 0.1;
			deltaZShift = getStandardNormal(RNG) * opt.deltaZ;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(RNG) * opt.deltaAx;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(RNG) * opt.deltaCross;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(RNG1) * 0.1;
			deltaXShift = getStandardNormal(RNG) * opt.deltaX;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
		}
	
		// Run Optimization
		repackSideChains(spm, opt.greedyCycles);

		vector<unsigned int> MCOFinal = spm.getMinStates()[0];
		sys.setActiveRotamers(MCOFinal);
		sys.calcEnergy();
		currentEnergy = spm.getMinBound()[0];
	
		if (!MCMngr.accept(currentEnergy)) {
			if (opt.verbose){
				cout << "state rejected   energy: " << currentEnergy-monomerEnergy << endl;
			}
		}
		else {
			bestEnergy = currentEnergy;
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");
	
			xShift = xShift + deltaXShift;
			crossingAngle = crossingAngle + deltaCrossingAngle;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift + deltaZShift;
			MCOBest = MCOFinal;
	
			if (opt.verbose){
				cout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
			}
		}
	counter++;
	}
	
	sys.applySavedCoor("savedBestState");

	double dimerEnergy = spm.getStateEnergy(MCOBest);
	double dimerEnergyDiff = dimerEnergy-opt.prevDimer;
	double finalEnergy = dimerEnergy-monomerEnergy;
	double vdw = spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	double hbond = spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	double imm1 = spm.getStateEnergy(MCOBest, "CHARMM_IMM1")+spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");

	// Print out info to the summary file
	sout << "Geometry " << endl;
	sout << "xShift: " << xShift << endl;
	sout << "crossingAngle: " << crossingAngle << endl;
	sout << "axialRotation: " << axialRotation << endl;
	sout << "zShift: " << zShift << endl << endl;
	sout << "Energy Summary Below" << endl;
	sout << "Monomer Energy: " << monomerEnergy << endl;
	sout << "Energy Before Repack: " << opt.prevDimer << endl << endl;
	sout << "Final Repack Energy: " << finalEnergy << ":" << opt.prevDimerVdw << ":" << vdw << ":" << vdw-opt.prevDimerVdw << ":" << opt.prevDimerHbond << ":" << hbond << ":" << hbond-opt.prevDimerHbond << ":" << opt.prevDimerIMM1 << ":" << imm1 << ":" << imm1-opt.prevDimerIMM1 << ":" << ":" << opt.xShift << ":" << opt.crossingAngle << ":" << opt.axialRotation << ":" << opt.zShift << ":" << xShift << ":" << crossingAngle << ":" << axialRotation << ":" << zShift << opt.thread << ":" << opt.sequence.length() << ":" << endl << endl;

	sout << Eset->getSummary() << endl << endl;
	
	// Add energies of the sequence to the energy map
	//map<string,double> &energyMap = _seqEnergyMap[sequence];
	//outputEnergiesByTerm(spm, MCOBest, energyMap, opt.energyTermList, "", 1);
	//saveEnergyDifference(opt, _seqEnergyMap, sequence);
	//seqEnergyMap[sequence][""] = dimerEnergy;
	//seqEnergyMap[sequence]["Total"] = finalEnergy;
	
	// Write an individual pdb for each sequence
	PDBWriter designWriter;
	designWriter.open(opt.outputDir + "/geomRepack_" + opt.designNumber + ".pdb");
	designWriter.write(sys.getAtomPointers(), true, false, true);
	designWriter.close();
	
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
}

/******************************************
 *  ======= FUNCTIONS DEFINITIONS =======
 ******************************************/
string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		if (it == _seq.begin() || it == _seq.end()-1){
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (it == _seq.begin()){
				if(resName == "HIS") {
					ps = ps + " HSE-ACE";
				} else {
					ps = ps + " " + resName + "-ACE";
				}
			}
			else{
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
		}
		else{
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	int frontExt = _opt.tmStart;
	int endExt = _opt.tmEnd;//TODO: I changed this from endResNum...it should still work but if energies are weird this is likely why
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(frontExt > i) {
				atoms += positions[i]->getAtomPointers();
			}
			if(endExt > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
			}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");
}

vector<pair<int,double>> calculateResidueBurial(System &_sys) {
	/*
	  SASA reference:
	  Protein Engineering vol.15 no.8 pp.659–667, 2002
	  Quantifying the accessible surface area of protein residues in their local environment
	  Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
	*/
	map<string,double> refSasa;
	refSasa["G"] = 83.91;
	refSasa["A"] = 116.40;
	refSasa["S"] = 125.68;
	refSasa["C"] = 141.48;
	refSasa["P"] = 144.80;
	refSasa["T"] = 148.06;
	refSasa["D"] = 155.37;
	refSasa["V"] = 162.24;
	refSasa["N"] = 168.87;
	refSasa["E"] = 187.16;
	refSasa["Q"] = 189.17;
	refSasa["I"] = 189.95;
	refSasa["L"] = 197.99;
	refSasa["H"] = 198.51;
	refSasa["K"] = 207.49;
	refSasa["M"] = 210.55;
	refSasa["F"] = 223.29;
	refSasa["Y"] = 238.30;
	refSasa["R"] = 249.26;
	refSasa["W"] = 265.42;

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator b(_sys.getAtomPointers());
	b.calcSasa();
	for (uint i = 0; i < _sys.positionSize()/2; i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdA = _sys.getPosition(i).getPositionId();
		//string posIdB = _sys.getPosition(i+_sys.positionSize()/2).getPositionId();
		string oneLetter = MslTools::getOneLetterCode(_sys.getPosition(i).getResidueName());
		double resiSasa = b.getResidueSasa(posIdA);
		double burial = resiSasa / refSasa[oneLetter];
		residueBurial.push_back(pair<int,double>(i, burial));
		//residueBurial.push_back(pair<string,double>(posIdB, burial));
		//cout << posId << "\t" << resiSasa << "\t" << burial << endl;
	}

	return residueBurial;
}

void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<pair<int,double>> &_resiBurial){
	//Repack side chains based on sasa scores
	cout << "Loading rotamers based on residue burial in the input structure" << endl;
	int levelCounter = 0;
	for (uint i = 0; i < _resiBurial.size(); i++) {
		double sasaPercentile = double(i) / double(_resiBurial.size());
		if (sasaPercentile > (levelCounter+1)/double(_opt.sasaRepackLevel.size())) {
			levelCounter++;
		}
		Position &posA = _sys.getPosition(_resiBurial[i].first);
		Position &posB = _sys.getPosition(_resiBurial[i].first+_sys.positionSize()/2);
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _opt.sasaRepackLevel[levelCounter];
				cout << posA.getPositionId() << ", " << posA.getResidueName() << "(" << _resiBurial[i].second*100 << "%% exposure): " << posRot << endl;
				cout << posB.getPositionId() << ", " << posB.getResidueName() << "(" << _resiBurial[i].second*100 << "%% exposure): " << posRot << endl;
				if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot, "", _opt.keepOriginalRotamer)) { 
						cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
					}
					if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot, "", _opt.keepOriginalRotamer)) { 
						cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
					}
				}
			}
		}
		else{
			string posRot = _opt.sasaRepackLevel[levelCounter];
			cout << posA.getPositionId() << ", " << posA.getResidueName() << "(" << _resiBurial[i].second*100 << "%% exposure): " << posRot << endl;
			cout << posB.getPositionId() << ", " << posB.getResidueName() << "(" << _resiBurial[i].second*100 << "%% exposure): " << posRot << endl;
			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot, "", _opt.keepOriginalRotamer)) { 
					cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
				}
				if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot, "", _opt.keepOriginalRotamer)) { 
					cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
				}
			}
		}
	}
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	
	/* Faster code 
	for (uint i=0; i < _apvA.size(); i++) {
			_apvB[i]->copyAllCoor(*_apvA[i]);
			vector<CartesianPoint*>& bCoors = _apvB[i]->getAllCoor();

			for(uint j = 0; j < bCoors.size(); j++) {
				bCoors[j]->setX(0 - bCoors[j]->getX());
				bCoors[j]->setY(0 - bCoors[j]->getY());
			}
					
		}
	*/

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
	trans.rotate(_apvB, m);
	
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

void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); 
	_spm.runGreedyOptimizer(_greedyCycles);
}

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
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
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//Required Options
	opt.allowed.push_back("designCrd");
	opt.allowed.push_back("designPdb");
	opt.allowed.push_back("configfile");

	//Design Parameters
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("state");
	opt.allowed.push_back("posVector");
	opt.allowed.push_back("burialVector");

	//Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("thread");
	
	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	
	//input monomerEnergy
	opt.allowed.push_back("monomer");
	opt.allowed.push_back("monoVdw");
	opt.allowed.push_back("monoHbond");
	opt.allowed.push_back("monoIMM1");
	
	//input previous dimerEnergy
	opt.allowed.push_back("prevDimer");
	opt.allowed.push_back("prevDimerVdw");
	opt.allowed.push_back("prevDimerHbond");
	opt.allowed.push_back("prevDimerIMM1");

	//Shift Size
	opt.allowed.push_back("deltaX");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaZ");

	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	//Load Rotamers from SASA values (from sgfc)
	opt.allowed.push_back("sasaRepack");
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("keepOriginalRotamer");
	
	//Other Options
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("dockHelices");
	opt.allowed.push_back("designNumber");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");
	
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

	// Design Parameters
	opt.sequence = OP.getString("sequence");
	if(OP.fail()) {
		opt.errorMessages += "sequence not specified using L\n";
		opt.errorFlag = true;
	}
	opt.state = OP.getUnsignedIntVector("state");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify state, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.posVector = OP.getIntVector("posVector");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify posVector, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.burialVector = OP.getDoubleVector("burialVector");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify burialVector, make sure they are space separated\n";
		opt.errorFlag = true;
	}

	// Geometry
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting to 9.2\n";
		opt.warningFlag = true;
		opt.xShift = 9.2;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to 25\n";
		opt.warningFlag = true;
		opt.crossingAngle = 25;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 40\n";
		opt.warningFlag = true;
		opt.axialRotation = 40;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to 2\n";
		opt.warningFlag = true;
		opt.zShift = 2;
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.thread = 23;
	}
	//TODO: maybe not in this code, but I feel like changing the thread of a sequence to see how it interacts at different threads could be helpful?

	//Load Rotamers using SASA values (from sgfc)
	opt.sasaRepackLevel = OP.getMultiString("sasaRepackLevel");
	if (OP.fail()) {
		opt.warningMessages += "sasaRepacklevel not specified! Default to one level at SL90.00";
		opt.sasaRepackLevel.push_back("SL90.00");
	}
	opt.keepOriginalRotamer = OP.getBool("keepOriginalRotamer");
	if (OP.fail()) {
		opt.warningMessages += "keepOriginalRotamer not specified, default true\n";
		opt.warningFlag = true;
	}

	//Monte Carlo variables
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC cycles not specified, default to 100\n";
		opt.warningFlag = true;
		opt.MCCycles = 100;
	}

	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.MCMaxRejects = 10;
		opt.warningMessages += "Number of MC max rejects not specified, default to using 10\n";
		opt.warningFlag = true;
	}

	opt.MCStartTemp = OP.getDouble("MCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCStartTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.MCStartTemp = 1000.0;
	}
	opt.MCEndTemp = OP.getDouble("MCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.MCEndTemp = 0.5;
	}
	opt.MCCurve = OP.getInt("MCCurve");
	if (OP.fail()) {
		opt.warningMessages += "MCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.MCCurve = 2;
	}

	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaX = 0.1;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 1.0;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 1.0;
	}
	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.1;
	}

	// monomer energies
	opt.monomer = OP.getDouble("monomer");
	if (OP.fail()) {
		opt.errorMessages += "monomer not specified\n";
		opt.errorFlag = true;
	}
	opt.monoVdw = OP.getDouble("monoVdw");
	if (OP.fail()) {
		opt.errorMessages += "monoVdw not specified\n";
		opt.errorFlag = true;
	}
	opt.monoHbond = OP.getDouble("monoHbond");
	if (OP.fail()) {
		opt.errorMessages += "monoHbond not specified\n";
		opt.errorFlag = true;
	}
	opt.monoIMM1 = OP.getDouble("monoIMM1");
	if (OP.fail()) {
		opt.errorMessages += "monoIMM1 not specified\n";
		opt.errorFlag = true;
	}

	// dimer energies
	opt.prevDimer = OP.getDouble("prevDimer");
	if (OP.fail()) {
		opt.errorMessages += "prevDimer not specified\n";
		opt.errorFlag = true;
	}
	opt.prevDimerVdw = OP.getDouble("prevDimerVdw");
	if (OP.fail()) {
		opt.errorMessages += "prevDimerVdw not specified\n";
		opt.errorFlag = true;
	}
	opt.prevDimerHbond = OP.getDouble("prevDimerHbond");
	if (OP.fail()) {
		opt.errorMessages += "prevDimerHbond not specified\n";
		opt.errorFlag = true;
	}
	opt.prevDimerIMM1 = OP.getDouble("prevDimerIMM1");
	if (OP.fail()) {
		opt.errorMessages += "prevDimerIMM1 not specified\n";
		opt.errorFlag = true;
	}

	// Weights
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
	opt.weight_solv = OP.getDouble("weight_solv");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_solv not specified, default 1.0\n";
		opt.weight_solv = 1.0;
	}
	// Input Files
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
		string envVar = "MSL_CHARMM_SOLV";
		if(SYSENV.isDefined(envVar)) {
			opt.solvFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n";
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
			opt.warningMessages += "hbondFile not specified using " + opt.hbondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hbondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine outputDir";
		opt.errorFlag = true;
	}
	opt.designPdb = OP.getString("designPdb");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine designPdb";
		opt.errorFlag = true;
	}
	opt.designCrd = OP.getString("designCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine designCrd";
		opt.errorFlag = true;
	}

	//Other Options
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
	}
	opt.dockHelices = OP.getBool("dockHelices");
	if (OP.fail()) {
		opt.warningMessages += "dockHelices not specified using false\n";
		opt.warningFlag = true;
		opt.dockHelices = false;
	}
	opt.designNumber = OP.getString("designNumber");
	if (OP.fail()) {
		opt.errorMessages += "designNumber not specified \n";
		opt.errorFlag = true;
	}
	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}

	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.sequence.length();
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
	}


	opt.rerunConf = OP.getConfFile();

	return opt;
}
