#include <sstream>
#include <iterator>
#include <unistd.h>
#include "backboneOptimizerFunctions_v2.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

/***********************************
 *version 2 functions
 ***********************************/
string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
	// A:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
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
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
		} else {
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
	return "A" + ps;
}

/***********************************
 *output file functions
 ***********************************/
// for running on chtc
void setupOutputDirectoryChtc(Options &_opt){
	_opt.outputDir = string(get_current_dir_name()) + "/" + _opt.sequence;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

void setupOutputDirectory(Options &_opt){
	_opt.outputDir = _opt.outputDir + "/" + _opt.sequence;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

void deleteTerminalHydrogenBondInteractions(System &_sys, int _firstResiNum, int _lastResiNum){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(_firstResiNum > i) {
				atoms += positions[i]->getAtomPointers();
			}
			if(_lastResiNum > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
			}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");
}

/***********************************
 *functions from geomRepack
 ***********************************/
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt){
	//Repack side chains based on sasa scores
	for (uint i = 0; i < _opt.rotamerSamplingVector.size()/2; i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+_opt.sequence.length());
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _opt.sasaRepackLevel[_opt.rotamerSamplingVector[i]];
				if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot)) { 
						cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
					}
					if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot)) { 
						cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
					}
				}
			}
		} else {
			string posRot = _opt.sasaRepackLevel[_opt.rotamerSamplingVector[i]];
			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&posA, posA.getResidueName(), posRot)) { 
					cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
				}
				if (!_sysRot.loadRotamers(&posB, posB.getResidueName(), posRot)) { 
					cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
				}
			}
		}
	}
}

/***********************************
 *repack functions
 ***********************************/
double computeMonomerEnergy(System & _sys, Options& _opt, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _mout){

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
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
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
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	int firstPos = 0;
    int lastPos = monoSys.positionSize();
	deleteTerminalHydrogenBondInteractions(monoSys, firstPos, lastPos);

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
	_mout << "Monomer - VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;

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
	_mout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		_mout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

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

			_mout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}

		counter++;
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	_mout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;


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
	if(_opt.verbose) {
		monoSys.calcEnergy();
		_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

		for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
			_mout << it->first << " " << it->second << endl;
		}
	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;

}
