#include <sstream>
#include <iterator>
#include <unistd.h>
#include "backboneOptimizerFunctions.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;


void deleteTerminalInteractions(System &_sys, BBOptions &_opt, int _firstResiNum, int _lastResiNum){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			// rid of hbonds from first 3 positions
			if(_firstResiNum <= i) {
				atoms += positions[i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			// rid of hbonds from last 3 positions
			if(_lastResiNum > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	for (uint i=0; i<_opt.deleteTerminalInteractions.size(); i++){
		pESet->deleteInteractionsWithAtoms(atoms,_opt.deleteTerminalInteractions[i]);
	}
}

/***********************************
 *output file functions
 ***********************************/
// for running on chtc
void setupOutputDirectory(BBOptions &_opt){
	//_opt.outputDir = string(get_current_dir_name()) + "/" + _opt.sequence;
	_opt.outputDir = string(get_current_dir_name()) + "/" + _opt.uniprotAccession;
	//_opt.outputDir = "/exports/home/gloiseau/mslib/trunk_AS/" + _opt.sequence + "_" + to_string(_opt.seed);
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}


/***********************************
 *geometry setup functions
 ***********************************/

void setGly69ToStartingGeometry(BBOptions &_opt, map<string,double> _geometries, System &_sys, System &_helicalAxis, Transforms &_trans) {
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	_sys.readPdb(_opt.backboneFile);

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// get helical axis atom pointers for setting up...?
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	double xShift = _geometries["xShift"];
	double zShift = _geometries["zShift"];
	double crossingAngle = _geometries["crossingAngle"];
	double axialRotation = _geometries["axialRotation"];

	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

void setGly69ToStartingGeometry(BBOptions &_opt, System &_sys, System &_helicalAxis, Transforms &_trans, double _crossingAngle, double _xShift) {
	/******************************************************************************
	 *         === COPY BACKBONE COORDINATES AND TRANSFORM TO GEOMETRY ===
	 ******************************************************************************/
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	_sys.readPdb(_opt.backboneFile);

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);

	// Set up chain A and chain B atom pointer vectors
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// get helical axis atom pointers for setting up...?
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();
	
	// Transform to chosen geometry
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, _opt.zShift, _opt.axialRotation, _crossingAngle, _xShift, _trans);
	moveZCenterOfCAMassToOrigin(_sys.getAtomPointers(), _helicalAxis.getAtomPointers(), _trans);
}

void prepareSystem(BBOptions &_opt, System &_sys, System &_startGeom, string &_polySeq){	
	// declare system
	CharmmSystemBuilder CSB(_sys,_opt.topFile,_opt.parFile,_opt.solvFile);
	if (_opt.useElec == false){
		CSB.setBuildTerm("CHARMM_ELEC", false);
	} else {
		CSB.setBuildTerm("CHARMM_ELEC", true);
	}
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	//CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	PolymerSequence PL(_polySeq);
	if(!CSB.buildSystem(PL)) {
		cerr << "Unable to build system from " << _polySeq << endl;
		exit(0);
	}

	// get chain A and B from the _system
	Chain & chainA = _sys.getChain("A");
	Chain & chainB = _sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	// assign the coordinates of our system to the given geometry 
	_sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	_sys.buildAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(_sys, _opt.hbondFile);
	hb.buildInteractions(30);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = _sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	//Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	if (_opt.useElec == true){
		Eset->setTermActive("CHARMM_ELEC", true);
		Eset->setWeight("CHARMM_ELEC", _opt.weight_elec);
	}

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	// removes all hydrogen bonding near the termini of our helices
	// (remnant from CATM, but used in the code that was used to get baselines so keeping it to be consistent)
	int firstPos = 0;
    int lastPos = _sys.positionSize();
    deleteTerminalInteractions(_sys,_opt,firstPos,lastPos);

	// Up to here is from readDPBAndCalcEnergy
	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	//CSB.updateNonBonded(10,12,50);

	// TODO: maybe calculate the energy here, then see if there's clashing. If so, then move helices away until no clashing, keeping the same
	// other coordinates? Just for simplicity for now. And if I want to implement this in design, adding this in will likely be a good idea.	
	// as of 2022-7-5: not sure if the above works or needs to be reworked
	PDBWriter writer1;
	writer1.open(_opt.outputDir + "/" + _opt.sequence + "_inputGeometry.pdb");
	writer1.write(_sys.getAtomPointers(), true, false, true);
	writer1.close();
}

/***********************************
 *repack functions
 ***********************************/

void threadDockAndRepack(BBOptions &_opt, System &_helicalAxis, int _thread, int _repackNumber, double _crossingAngle, double _monomerEnergy, 
 map<string,double> _monomerEnergyByTerm, ofstream &_eout){
	// setup the output file
	ofstream sout;  // summary file output
	string soutfile = _opt.outputDir + "/RepackSummaries_" + to_string(_repackNumber) + "/repackSummary_" + to_string(_thread) + "_" + to_string(_crossingAngle) + ".out";
	sout.open(soutfile.c_str());
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// set the helical axis at the origin
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

	// get the starting geometry  
	vector<double> densities;
	map<string,double> geometries;
	// if true, get random axial rotation and z shift from the geometry density file
	if (_opt.getRandomAxAndZ){
		getAxialRotAndZShift(_opt, geometries, densities);
	} else {
		geometries["axialRotation"] = _opt.axialRotation;
		geometries["zShift"] = _opt.zShift;
		densities.push_back(0);
		densities.push_back(0);
	}
	geometries["xShift"] = _opt.xShift;
	geometries["crossingAngle"] = _crossingAngle;

	// Output to summary file
	sout << "***STARTING GEOMETRY:***" << endl;
	for (auto &geometry : geometries){
		sout << geometry.first << ": " << geometry.second << endl;
	}

	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(_opt,geometries,startGeom,helicalAxis,trans);

	/******************************************************************************
	 *                       === SETUP STARTING GEOMETRY ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	string polySeq = convertToPolymerSequence(_opt.sequence, _thread);
	prepareSystem(_opt, sys, startGeom, polySeq);
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?

	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Assign number of rotamers by residue burial
	loadRotamers(sys, sysRot, _opt.SL);
	sout << "SL: " << _opt.SL << endl;
	
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 
	
	// _optimize Initial Starting Position 
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.updateWeights();
	spm.setOnTheFly(true);// changed to make more similar to CATM
	spm.saveEnergiesByTerm(true);
	
	// Repack dimer
	repackSideChains(spm, _opt.greedyCycles);
	vector<uint> startStateVec = spm.getMinStates()[0];
	double currentEnergy = spm.getMinBound()[0];
	sys.setActiveRotamers(startStateVec);

	// calculate the energy of the system
	double startDimer = sys.calcEnergy();
	double totalEnergy = startDimer-_monomerEnergy;
	sout << "Thread " << _thread << endl;
	sout << " -Dimer Energy: " << startDimer << endl;
	sout << " -Monomer Energy: " << _monomerEnergy << endl;
	sout << " -Total Energy: " << totalEnergy << endl;
	sout << spm.getSummary(startStateVec) << endl;

	// get the start xShift and dock helices	
	double savedXShift = _opt.xShift;
	sys.saveAltCoor("savedBestState");
	localXShiftDocking(sys, _opt, totalEnergy, _monomerEnergy, spm, helicalAxis, trans, _thread, savedXShift, sout);

	// set the xShift to the best xShift from docking
	geometries["xShift"] = savedXShift;

	// get geometries from the geometry map
	double xShift = geometries["xShift"];
	double zShift = geometries["zShift"];
	double crossingAngle = geometries["crossingAngle"];
	double axialRotation = geometries["axialRotation"];

	// if the energy is higher than the energy cutoff, don't repack
	if (totalEnergy > _opt.energyCutoff){
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is higher than cutoff " << _opt.energyCutoff << ", skip" << endl;
	} else {
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is lower than cutoff " << _opt.energyCutoff << ", continue" << endl;

		// output starting geometry before repack
		double bestEnergy = startDimer-_monomerEnergy;
		sout << "Starting Geometry" << endl;
		sout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << endl << endl;
		sout << "Current Best Energy: " << bestEnergy << endl;

		// Local backbone monte carlo repacks
		localBackboneRepack(_opt, sys, geometries, densities, spm, helicalAxis, trans, RNG, _monomerEnergyByTerm, _monomerEnergy, _thread, sout, _eout);
	}
	sout.close();
}

// gets a random axial rotation and z shift from the geometry density file
void getAxialRotAndZShift(BBOptions &_opt, map<string,double> &_geometries, vector<double> &_densities){
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed);// not sure if this works without the seed 0 (time based seed)
	// Setup file reader
	Reader reader(_opt.geometryDensityFile);
	reader.open();
	if(!(reader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.geometryDensityFile << endl;
		exit(0);
	}
	vector<string> lines = reader.getAllLines();

	// Extract the geometries from a random line of the geometry file
	int geometryLine = RNG.getRandomInt(1,lines.size()-1);
	vector<string> tokens = MslTools::tokenize(lines[geometryLine], "\t");//xShift, crossingAngle, axialRotation, zShift, angleDistDensity, axialRotationDensity, zShiftDensity
	_geometries["axialRotation"] = MslTools::toDouble(tokens[2]);
	_geometries["zShift"] = MslTools::toDouble(tokens[3]);
	double axialRotationDensity = MslTools::toDouble(tokens[5]);
	double zShiftDensity = MslTools::toDouble(tokens[6]);
	_densities.push_back(axialRotationDensity);
	_densities.push_back(zShiftDensity);
}

void localXShiftDocking(System &_sys, BBOptions &_opt, double &_bestEnergy, double _monomerEnergy, SelfPairManager &_spm, 
 System &_helicalAxis, Transforms &_trans, int _thread, double &_savedXShift, ofstream &_out) {
	// Get helical axis atom pointers 
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();
	
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector & apvChainB = _sys.getChain("B").getAtomPointers();

	vector<unsigned int> MCOBest;
	double deltaXShift = -0.1; // xShift changes
	double xShiftEnd = 6; // xShift to stop at
	double xShift = _opt.xShift; // current xShift
	double previousEnergy = _monomerEnergy; // previous energy to compare to
	double globalLowestE = _monomerEnergy; // Global lowest energy found (if above monomer we won't save anyways)
	double currentEnergy = _bestEnergy; // current energy
	// while loop for x-shifts	
	while (xShift >= xShiftEnd) {
		// add the xShift change to the current xShift
		xShift += deltaXShift;
		// Move the helix
		backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaXShift, 3 );

		// Run Optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal;
		MCOFinal = _spm.getMinStates()[0];
		_sys.setActiveRotamers(MCOFinal);

		// get current energy
		currentEnergy = _spm.getMinBound()[0]-_monomerEnergy;

		// Check if this is the lowest energy
		if (currentEnergy < _bestEnergy) {
			_bestEnergy = currentEnergy;
			_savedXShift = xShift;
			_sys.saveAltCoor("savedBestState");
			MCOBest = MCOFinal;
			//_helicalAxis.saveAltCoor("BestAxis");
		}
		_out << "xShift: " << xShift << " energy: " << currentEnergy << endl;

		// If energy increase twice in a row, and it is above the global lowest energy, quit
		if (currentEnergy < globalLowestE) {
			globalLowestE = currentEnergy;
		}
		if (currentEnergy > (globalLowestE+10.0) && previousEnergy > (globalLowestE+10.0) && currentEnergy > previousEnergy) {
			_out << "Energy increasing above global lowest energy... (currently " << globalLowestE << ")" << endl;
			break;
		} else {
			previousEnergy = currentEnergy;
		}
	}
	_out << "Thread " << _thread << " Best Energy at x shift: " << _bestEnergy << " at " << _savedXShift << endl;
}


void localBackboneRepack(BBOptions &_opt, System &_sys, map<string,double> _geometries, vector<double> _densities, SelfPairManager &_spm,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm,
 double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout){
	_out << "==============================================" << endl;
	_out << " Performing Local Monte Carlo Backbone Repack " << endl;
	_out << "==============================================" << endl;
	// set the system to the best xShift state from docking
	_sys.applySavedCoor("savedBestState");
	// save the current state as the repack state
	_sys.saveAltCoor("savedRepackState");
	// get the best repack energy
	double prevBestEnergy = _sys.calcEnergy();

	// do backbone geometry repacks
	monteCarloRepack(_opt, _sys, _spm, _geometries, _densities, _helicalAxis, _trans, _RNG, prevBestEnergy, 
	 _monomerEnergyByTerm, _monomerEnergy, _thread, _out, _eout);
  
    //outputs a pdb file for the structure 
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.outputDir + "/backboneOptimized_" + to_string(_thread) + ".pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
}

void monteCarloRepack(BBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string,double> _geometries, vector<double> _densities,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 map<string,double> _monomerEnergyByTerm, double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout){
	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);

	// starting geometry
	double xShift = _geometries["xShift"];
	double zShift = _geometries["zShift"];
	double crossingAngle = _geometries["crossingAngle"];
	double axialRotation = _geometries["axialRotation"];

	// Monte Carlo Repack Manager Setup
	MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects);
	MCMngr.setEner(_prevBestEnergy);

	vector<uint> startStateVec = _spm.getMinStates()[0];
	vector<unsigned int> MCOBest = startStateVec;
		
	unsigned int counter = 0;
	double currentEnergy = _prevBestEnergy;
	double startDimer = _prevBestEnergy;
	
	// Get helical axis atom pointers 
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();
	
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector & apvChainB = _sys.getChain("B").getAtomPointers();
	
	PDBWriter writer;

	// setup output file for repack
	ofstream sout;  // summary file output
	string soutfile = _opt.outputDir + "/thread_" + to_string(_thread) + "_" + to_string(crossingAngle) + "_energy.csv";
	sout.open(soutfile.c_str());

	// setup variables for shifts
	bool decreaseMoveSize = _opt.decreaseMoveSize; // if true, decreasing the move size to the minimum move size throughout the repack with accepts
	double deltaX = _opt.deltaX;
	double deltaCross = _opt.deltaCross;
	double deltaAx = _opt.deltaAx;
	double deltaZ = _opt.deltaZ;

	// loop through the MC cycles for backbone repacks
	while(!MCMngr.getComplete()) {
		double startTemp = MCMngr.getCurrentT();

		_sys.applySavedCoor("savedRepackState");
		//_helicalAxis.applySavedCoor("BestRepack");
		
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
			deltaZShift = getStandardNormal(_RNG) * deltaZ;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 1) {
		//===========================
		//===== Axial Rotation ======
		//===========================
			//deltaAxialRotation = getStandardNormal(_RNG1) * 1.0;
			deltaAxialRotation = getStandardNormal(_RNG) * deltaAx;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * deltaCross;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaZShift, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * deltaX;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaZShift, moveToPreform);
		}
		
		// Run Optimization
		repackSideChains(_spm, _opt.greedyCycles);

		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0];
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE
		
		if (!MCMngr.accept(currentEnergy)) {
			sout << "MCReject   xShift: " << xShift+deltaXShift << " crossingAngle: " << crossingAngle+deltaCrossingAngle << " axialRotation: " << axialRotation+deltaAxialRotation << " zShift: " << zShift+deltaZShift << " energy: " << currentEnergy-_monomerEnergy << endl;
		} else {
			_prevBestEnergy = currentEnergy;
			_sys.saveAltCoor("savedRepackState");
			//_helicalAxis.saveAltCoor("BestRepack");
		
			xShift = xShift + deltaXShift;
			crossingAngle = crossingAngle + deltaCrossingAngle;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift + deltaZShift;
			MCOBest = MCOFinal;

			// if accept, decrease the value of the moves by the sigmoid function
			if (decreaseMoveSize == true){
				double endTemp = MCMngr.getCurrentT();
				getCurrentMoveSizes(_opt, startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, decreaseMoveSize);
			}
			sout << "MCAccept " << counter <<  " xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-_monomerEnergy << endl;
			counter++;
			writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	writer.close();
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
		
	_sys.applySavedCoor("savedRepackState");
	currentEnergy = _spm.getMinBound()[0];
	sout << "Best: " << _spm.getSummary(MCOBest) << endl;

	// TODO: make this into a map that saves all of these to be output
	double dimerEnergy = _sys.calcEnergy();
	double finalEnergy = dimerEnergy-_monomerEnergy;
	double vdw = _spm.getStateEnergy(MCOBest, "CHARMM_VDW");
	double hbond = _spm.getStateEnergy(MCOBest, "SCWRL4_HBOND");
	double imm1 = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1");
	double imm1Ref = _spm.getStateEnergy(MCOBest, "CHARMM_IMM1REF");
	double dimerDiff = dimerEnergy-startDimer;

	// get monomer energies from monomerEnergyTerm map
	double monomerVdw = _monomerEnergyByTerm.find("CHARMM_VDW")->second;
	double monomerHbond = _monomerEnergyByTerm.find("SCWRL4_HBOND")->second;
	double monomerImm1 = _monomerEnergyByTerm.find("CHARMM_IMM1")->second;
	double monomerImm1Ref =_monomerEnergyByTerm.find("CHARMM_IMM1REF")->second;

	// calculate the solvent accessible surface area
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	// output the results
	double axialRotDensity = _densities[0];
	double zShiftDensity = _densities[1];
	if (_opt.useElec == false){
		_out << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		_out << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		_out << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		_out << MCMngr.getReasonCompleted() << endl;
		sout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		sout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		sout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		sout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		sout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		sout << MCMngr.getReasonCompleted() << endl;
		_eout << finalEnergy << ',' << _thread << ',' << xShift << ',' << crossingAngle << ',' << axialRotation << ',' << zShift << ',' << axialRotDensity << ',' << zShiftDensity << endl;
	} else {
		double elec = _spm.getStateEnergy(MCOBest, "CHARMM_ELEC");
		double monomerElec =_monomerEnergyByTerm.find("CHARMM_ELEC")->second;
		_out << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, elec, monoElec, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		_out << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		_out << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		_out << MCMngr.getReasonCompleted() << endl;
		sout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, elec, monoElec, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		sout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		sout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		sout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		sout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		sout << MCMngr.getReasonCompleted() << endl;
		_eout << finalEnergy << ',' << _thread << ',' << xShift << ',' << crossingAngle << ',' << axialRotation << ',' << zShift << ',' << axialRotDensity << ',' << zShiftDensity << endl;
	}
	sout << "Thread " << _thread << " finished. Energy: " << finalEnergy << endl;
	sout.close();
}

void getCurrentMoveSizes(BBOptions &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize) {
	double decreaseMultiplier = _endTemp/_currTemp;
	bool decreaseX = true;
	bool decreaseCross = true;
	bool decreaseAx = true;
	bool decreaseZ = true;
	if (_opt.dockHelices == false){
		//if dockHelices false, adjust deltaX to be negative
		_deltaX = increaseMoveSize(_deltaX, _opt.deltaXLimit, decreaseMultiplier, decreaseX);
		if (decreaseX == false){
			_opt.dockHelices = true; // finish bringing helices together
		}
	} else {
		_deltaX = decreaseMoveSize(_deltaX, _opt.deltaXLimit, decreaseMultiplier, decreaseX);
	}
	_deltaCross = decreaseMoveSize(_deltaCross, _opt.deltaCrossLimit, decreaseMultiplier, decreaseCross);
	_deltaAx = decreaseMoveSize(_deltaAx, _opt.deltaAxLimit, decreaseMultiplier, decreaseAx);
	_deltaZ = decreaseMoveSize(_deltaZ, _opt.deltaZLimit, decreaseMultiplier, decreaseZ);	
	if (decreaseX == false && decreaseCross == false && decreaseAx == false && decreaseZ == false){
		_decreaseMoveSize = false;
	}
}

double increaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease) {
	double diffMoveSize = _moveSize - _moveLimit;
	double moveDecrease = diffMoveSize * _decreaseMultiplier;
	double newMoveSize = _moveSize - moveDecrease;
	if (newMoveSize < _moveLimit){
		return newMoveSize;
	} else {
		_decrease = false;
		return _moveSize;
	}
}

double computeMonomerEnergy(System & _sys, BBOptions & _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout){

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	// Objects used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
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
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	if (_opt.useElec){
		CSBMono.setBuildTerm("CHARMM_ELEC", true);
	} else {
		CSBMono.setBuildTerm("CHARMM_ELEC", false);
	}
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());
	

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	int firstPos = 0;
    int lastPos = monoSys.positionSize();
	deleteTerminalInteractions(monoSys,_opt, firstPos, lastPos);

	/******************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	SelfPairManager monoSpm;
	monoSpm.seed(_opt.seed);
	monoSpm.setVerbose(_opt.verbose);

	for (uint k=0; k < monoSys.positionSize(); k++) {
		Position &pos = monoSys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!monoRot.loadRotamers(&pos, pos.getResidueName(), "SL97.00")) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	//_mout << "Monomer - VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;
	if (_opt.useElec == true){
		monoEset->setTermActive("CHARMM_ELEC", true);
		monoEset->setWeight("CHARMM_ELEC", _opt.weight_elec);
	}
	monoSpm.setSystem(&monoSys);
	monoSpm.updateWeights();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	trans.translate(chainA, interDistVect);


	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	trans.translate(chainA, move5Down);
	double bestZ = -5.0;

	monoSys.calcEnergy();

	// Repack side chains
	monoSpm.setOnTheFly(1);
	monoSpm.calculateEnergies();
    monoSpm.runGreedyOptimizer(_opt.greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	helicalAxis.saveAltCoor("BestMonomerAxis");
	//_mout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_mout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		for(int j=0; j<=3; j++) { // test at 4 rotations 0, 90, 180 and 270 degrees
			//==================================
			//====== Axial Rot ======
			//==================================
			monoAxialRotation = j * 90.0;

			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			_mout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			//monoSys.writePdb("mono_" + MslTools::doubleToString(monoTilt) + "_" + MslTools::doubleToString(monoAxialRotation) + ".pdb");

			if(currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				bestTilt = monoTilt;
				bestRotation = monoAxialRotation;
				monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
				monoSys.saveAltCoor("savedBestMonomer");
				helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

		}
	}
	
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	// MonteCarloManager MCMngr(1000.0, 0.5, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		helicalAxis.applySavedCoor("BestMonomerAxis");

		int moveToPreform = _RNG.getRandomInt(2);

		double deltaZShift = 0.0;
		double deltaTilt = 0.0;
		double deltaAxialRotation = 0.0;

		//======================================
		//====== Z Shift ======
		//======================================
		if (moveToPreform == 0) {
			deltaZShift = getStandardNormal(_RNG) * 1.0;
			CartesianPoint translateA = axisA(1).getCoor() - axisA(0).getCoor(); // vector minus helical center
			translateA = translateA.getUnit() * deltaZShift; // unit vector of helical _axis times the amount to shift by
			trans.translate(chainA, translateA);
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_mout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_mout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_mout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_mout << "state rejected   energy: " << currentEnergy << endl;
		}
		else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			//_mout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}

		counter++;
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	//_mout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");
	_mout << monoEset->getSummary();
	_mout << endl;

	// print the monomer
	string monoOutCrdFile  = _opt.outputDir + "/monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}

	string monoOutPdbFile  = _opt.outputDir + "/monomer.pdb";
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	// Store monomer energy by term
	monoSys.calcEnergy();
	_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

	for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
		_mout << it->first << " " << it->second << endl;
	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;
}

/****************************************
 *
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/
BBOptions BBParseOptions(int _argc, char * _argv[]){

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a BBOptions structure
	 *  defined at the head of this file
	 ******************************************/
	BBOptions opt;

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

	//optional
	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	opt.allowed.push_back("weight_elec");
	opt.allowed.push_back("verbose");

	//Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("geometryDensityFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("pdbName");
	opt.allowed.push_back("configfile");

	//
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("seed");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("thread");
	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("greedyCycles");
	
	//Shift Size
	opt.allowed.push_back("deltaX");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaXLimit");
	opt.allowed.push_back("deltaCrossLimit");
	opt.allowed.push_back("deltaAxLimit");
	opt.allowed.push_back("deltaZLimit");
	opt.allowed.push_back("decreaseMoveSize");
	
	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");

	opt.allowed.push_back("negAngle");
	opt.allowed.push_back("negRot");

	//version 2 and up
	opt.allowed.push_back("useElec");
	opt.allowed.push_back("backboneFile");
	opt.allowed.push_back("helicalAxis");
	opt.allowed.push_back("useAlaAtCTerminus");
	opt.allowed.push_back("deleteTerminalInteractions");
	opt.allowed.push_back("uniprotAccession");
	opt.allowed.push_back("dockHelices");
	opt.allowed.push_back("crossAngle");
	opt.allowed.push_back("getRandomAxAndZ");
	opt.allowed.push_back("energyCutoff");
	opt.allowed.push_back("numRepacks");
	
	//Rotamers
	opt.allowed.push_back("SL");

	//Begin Parsing through the options
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	if (OP.countOptions() == 0){
		opt.errorMessages += "No options given!\n";
		return opt;
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	
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

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.negAngle = OP.getBool("negAngle");
	if (OP.fail()) {
		opt.warningMessages += "negAngle not specified using false\n";
		opt.warningFlag = true;
		opt.negAngle = false;
	}
	opt.negRot = OP.getBool("negRot");
	if (OP.fail()) {
		opt.warningMessages += "negRot not specified using false\n";
		opt.warningFlag = true;
		opt.negAngle = false;
	}

	// tm parameters
	opt.sequence = OP.getString("sequence");
	if(OP.fail()) {
		opt.errorMessages += "sequence not specified using L\n";
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
	if (opt.negAngle == true){
		opt.crossingAngle = -opt.crossingAngle;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to 40\n";
		opt.warningFlag = true;
		opt.axialRotation = 40;
	}
	if (opt.negRot == true){
		opt.axialRotation = -opt.axialRotation;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to 2\n";
		opt.warningFlag = true;
		opt.zShift = 2;
	}
	// thread parameters
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 25\n";
		opt.warningFlag = true;
		opt.thread = 25;
	}
	opt.threadStart = OP.getInt("threadStart");
	if (OP.fail()) {
		opt.warningMessages += "threadStart not specified, defaulting to 17\n";
		opt.warningFlag = true;
		opt.threadStart = 17;
	}
	opt.threadEnd = OP.getInt("threadEnd");
	if (OP.fail()) {
		opt.warningMessages += "threadEnd not specified, defaulting to 30\n";
		opt.warningFlag = true;
		opt.threadEnd = 30;
	}
	//TODO: maybe not in this code, but I feel like changing the thread of a sequence to see how it interacts at different threads could be helpful?

	//Monte Carlo variables
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC cycles not specified, default to 100\n";
		opt.warningFlag = true;
		opt.MCCycles = 100;
	}

	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "Number of MC max rejects not specified, default to using 5\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 5;
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
		opt.warningMessages += "MCCurve not specified using SIGMOID(3)\n";
		opt.warningFlag = true;
		opt.MCCurve = 3;
	}
	
    //Weights
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
	opt.weight_elec = OP.getDouble("weight_elec");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_elec not specified, default 1.0\n";
		opt.weight_elec = 1.0;
	}

	//Shift Size
	opt.dockHelices = OP.getBool("dockHelices");
	if (OP.fail()) {
		opt.warningMessages += "dockHelices not specified using true\n";
		opt.warningFlag = true;
		opt.dockHelices = true;
	}
	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaX = 0.5;
	}
	
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 5.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 5.0;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 5.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 5.0;
	}
	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.5;
	}
	opt.deltaXLimit = OP.getDouble("deltaXLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaXLimit not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaXLimit = 0.1;
	}
	opt.deltaCrossLimit = OP.getDouble("deltaCrossLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaCrossLimit not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaCrossLimit = 1.0;
	}
	opt.deltaAxLimit = OP.getDouble("deltaAxLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaAxLimit not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaAxLimit = 1.0;
	}
	opt.deltaZLimit = OP.getDouble("deltaZLimit");
	if (OP.fail()) {
		opt.warningMessages += "deltaZLimit not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaZLimit = 0.1;
	}
	opt.decreaseMoveSize = OP.getBool("decreaseMoveSize");
	if (OP.fail()) {
		opt.warningMessages += "decreaseMoveSize not specified using true\n";
		opt.warningFlag = true;
		opt.decreaseMoveSize = true;
	}
	if (opt.dockHelices == false){
		opt.deltaX = opt.deltaX*(-1);
		opt.deltaXLimit = opt.deltaXLimit*(-1);
	}
    //Parameter files
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

	opt.geometryDensityFile = OP.getString("geometryDensityFile");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine geometryDensityFile, defaulting to original density file\n";
		opt.warningFlag = true;
		opt.geometryDensityFile = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/2021_09_28_geometryDensityFile.txt";
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
		opt.warningMessages += "Unable to determine outputDir, using current directory\n";
		opt.warningFlag = true;
	}

	// Monomer Options
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
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
		opt.MCMaxRejects = 10;
	}

	// version 2 options
	opt.useElec = OP.getBool("useElec");
	if (OP.fail()) {
		opt.warningMessages += "useElec not specified using true\n";
		opt.warningFlag = true;
		opt.useElec = true;
	}
	opt.backboneFile = OP.getString("backboneFile");
	if (OP.fail()) {
		opt.errorMessages += "backboneFile not specified\n";
		opt.errorFlag = true;
	}
	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.errorMessages += "helicalAxis not specified\n";
		opt.errorFlag = true;
	}
	opt.useAlaAtCTerminus = OP.getBool("useAlaAtCTerminus");
	if (OP.fail()) {
		opt.warningMessages += "useAlaAtCTerminus not specified using true\n";
		opt.warningFlag = true;
		opt.useAlaAtCTerminus = true;
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningMessages += "backboneLength not specified using 21\n";
		opt.warningFlag = true;
		opt.backboneLength = 21;
	}
	opt.deleteTerminalInteractions = OP.getMultiString("deleteTerminalInteractions");
	if (OP.fail()) {
		opt.deleteTerminalInteractions.push_back("");
		opt.warningMessages += "deleteTerminalInteractions not specified\n";
		opt.warningFlag = true;
	}
	opt.uniprotAccession = OP.getString("uniprotAccession");
	if (OP.fail()) {
		opt.uniprotAccession = "PROTEIN_UNK";
		opt.warningMessages += "uniprotAccession not specified using " + opt.uniprotAccession + "\n";
		opt.warningFlag = true;
	}
	opt.crossAngle = OP.getMultiDouble("crossAngle");
	if (OP.fail()) {
		opt.crossAngle.push_back(opt.crossingAngle);
		opt.warningMessages += "crossAngle not specified\n";
		opt.warningFlag = true;
	}
	
	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}
	opt.getRandomAxAndZ = OP.getBool("getRandomAxAndZ");
	if (OP.fail()) {
		opt.warningMessages += "getRandomAxAndZ not specified using false\n";
		opt.warningFlag = true;
		opt.getRandomAxAndZ = false;
	}
	opt.energyCutoff = OP.getDouble("energyCutoff");
	if (OP.fail()) {
		opt.warningMessages += "energyCutoff not specified, defaulting to 100\n";
		opt.warningFlag = true;
		opt.energyCutoff = 100;
	}
	opt.numRepacks = OP.getInt("numRepacks");
	if (OP.fail()) {
		opt.warningMessages += "numRepacks not specified, defaulting to 5\n";
		opt.warningFlag = true;
		opt.numRepacks = 5;
	}
	opt.rerunConf = OP.getConfFile();
	return opt;
}
