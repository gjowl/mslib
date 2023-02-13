#include <iostream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>
#include <chrono>

// MSL Functions
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "CRDReader.h"
#include "SasaCalculator.h"

// My functions
#include "mutateAndBBOptimizeOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "mutateAndBBOptimize";//TODO: better name
string programDescription = "Designs sequences for backbone geometries extracted from the PDB, optimizing specifically for vdW energies";
string programAuthor = "Gilbert Loiseau";
string programVersion = "2";
string programDate = "18 August 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime, spmTime;
auto start = chrono::system_clock::now();

// Functions

/*
    - Read in the pdb file from command line into a system
    - Thread loop:
        - Through backbone optimize procedure for each sequence (mutate one position at a time)
        - Save energies into a map
        - print out the map of energies to compare to
*/
void optimizeBackboneGeometry(MBBOptions &_opt, System &_startGeom, string _sequence, map<string,map<string,double>> &_sequenceEnergyMap, RandomNumberGenerator &_RNG);
double monteCarloRepack(MBBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string, map<string,double>> &_sequenceEnergyMap,
 string _sequence, vector<uint> &_bestState, System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy);
void convertToRelativeAxAndZ(double _axialRot, double _zShift, double &_relativeAx, double &_relativeZ);
void convertToAxAndZForTranformation(MBBOptions &_opt);
void computeMonomerEnergy(System &_sys, System &_helicalAxis, MBBOptions &_opt, Transforms &_trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG) ;
void setupDesignDirectory(MBBOptions &_opt);
void outputGeometry(MBBOptions &_opt, double _xShift, double _crossingAngle, double _axialRotation, double _zShift, ofstream &_out);
void getCurrentMoveSizes(MBBOptions &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);
void deleteTerminalBondInteractions(System &_sys, MBBOptions &_opt, int _firstResiNum, int _lastResiNum);
void prepareSystem(MBBOptions &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS);
void transformHelicalAxis(AtomPointerVector &_axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans);
void outputEnergyFile(MBBOptions &_opt, double _seed, map<string,map<string,double>> &_sequenceEnergyMap);

// help functions
void usage();
void outputErrorMessage(MBBOptions &_opt);
MBBOptions parseMBBOptions(int _argc, char * _argv[]);

int main(int argc, char *argv[]){
	// setup time
	time_t startRealTime = chrono::system_clock::to_time_t(start); 
	cout << "Program Start time: " << ctime(&startRealTime) << endl;
	
	// initialize time variables
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	// start the timer for the program
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	// setup time and date
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);

	//parse command line options
	MBBOptions opt = parseMBBOptions(argc, argv);
	if (opt.errorFlag) {
		outputErrorMessage(opt);
		exit(1);
	} else if (!opt.errorFlag && !opt.warningFlag && opt.errorMessages != ""){
		outputErrorMessage(opt);
		usage();
		exit(0);
	}
	
	// setup output files and directory
	ofstream sout; // summary file output
	ofstream err; // error file output
	ofstream rerun; // rerun config output

	setupDesignDirectory(opt); // makes a directory in the directory that you run from, with design_<runNumber> as the name

	// setup output files
	string soutfile = opt.outputDir + "/summary.out";
	string errfile  = opt.outputDir + "/errors.out";
	string rerunfile = opt.outputDir + "/rerun.config";
	
	// open output files
	sout.open(soutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	// setup the rerun config file
	rerun << opt.rerunConf << endl;
	rerun.close();

	sout << date << endl; // output the date and time to the summary file

	// get the starting geometries; convert to parallelogram axialRot and Z (Mueller 2014; Fig. S1)
	convertToAxAndZForTranformation(opt);
	cout << "***STARTING GEOMETRY:***" << endl;
	cout << "xShift:        " << opt.xShift << endl;
	cout << "crossingAngle: " << opt.crossingAngle << endl;
	cout << "axialRotation: " << opt.axialRotation << endl;
	cout << "zShift:        " << opt.zShift << endl << endl;

	// String for the alternateIds at the interface
	if (opt.xShift <= 7.5){
		opt.Ids.push_back("GLY");
	}

    // Read in the pdb for the starting coordinates for the backbone
	System sys;
	sys.readPdb(opt.inputPdbFile);

    // get the sequence of the pdb
    Chain &chainA = sys.getChain("A");
	string startSequence = convertPolymerSeqToOneLetterSeq(chainA);

	// Initialize RandomNumberGenerator object with seed (time or given seed number)
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed);//defaults to 0, which is a time based seed

    // convert the interface string to an interface vector
    vector<uint> interfacePositions;
    for (uint i=0; i<opt.interface.length(); i++){
        // get interface number
        string posString = opt.interface.substr(i, 1);
        int posInt = MslTools::toInt(posString);
        if (posInt != 0){
            interfacePositions.push_back(i);
        }
    }

    // loop through the interface positions
    // setup the thread vector
    vector<thread> threads;
    map<string, map<string, double>> sequenceEnergyMap;
    cout << startSequence << endl;
    for (uint i=0; i<interfacePositions.size(); i++){
        // get the first interface position
        uint interfacePosition = interfacePositions[i];
        for (uint j=0; j<opt.Ids.size(); j++){
            // get the first id
            string id = opt.Ids[j];
		    string currAA = MslTools::getThreeLetterCode(startSequence.substr(interfacePosition, 1));
		    if (currAA != id){
			    // replace the id at the position in bestSeq with the current id to get current sequence
                string mutantSequence = startSequence;
			    string oneLetterId = MslTools::getOneLetterCode(id);
			    mutantSequence.replace(interfacePosition, 1, oneLetterId);
                cout << mutantSequence << endl;
                // thread loop begin for optimizing backbone geometry for the sequence 
                threads.push_back(thread{optimizeBackboneGeometry, ref(opt), ref(sys), mutantSequence, ref(sequenceEnergyMap), ref(RNG)});
				//optimizeBackboneGeometry(opt, sys, mutantSequence, sequenceEnergyMap, RNG);
            }
        }
    }
    // join all the threads (wait for them all to finish before continuing)
    for (auto& th : threads){
    	th.join();
    }

	// output the energy file
	double seed = RNG.getSeed();
	outputEnergyFile(opt, seed, sequenceEnergyMap);
}

//Functions
void optimizeBackboneGeometry(MBBOptions &_opt, System &_startGeom, string _sequence, map<string,map<string,double>> &_sequenceEnergyMap, RandomNumberGenerator &_RNG){
	
    // Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// set the helical axis at the origin
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

	// 
	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	
	transformHelicalAxis(axisA, axisB, ori, xAxis, zAxis, _opt.zShift, _opt.axialRotation, _opt.crossingAngle, _opt.xShift, trans);

    // prepare polymer sequence 
	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PS(polySeq);
	cout << polySeq << endl;

	// set up the system for the input sequence
	System sys;
	prepareSystem(_opt, sys, _startGeom, PS);
	
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// load the rotamers
	loadRotamers(sys, sysRot, _opt.SL);

	// get chain A and B from the system
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();
	
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	//spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();
	spm.runGreedyOptimizer(_opt.greedyCycles);

    // compute monomer energy
	computeMonomerEnergy(sys, helicalAxis, _opt, trans, _sequenceEnergyMap, _sequence, _RNG);
	double monomerEnergy = _sequenceEnergyMap[_sequence]["Monomer"];
	repackSideChains(spm, _opt.greedyCycles);
	vector<unsigned int> MCOFinal = spm.getMinStates()[0];
	double currentEnergy = spm.getMinBound()[0]-monomerEnergy;
    _sequenceEnergyMap[_sequence]["preBBOptimizeEnergy"] = currentEnergy;
	sys.setActiveRotamers(MCOFinal);
	double finalEnergy = monteCarloRepack(_opt, sys, spm, _sequenceEnergyMap, _sequence, MCOFinal, helicalAxis, trans, _RNG, monomerEnergy);

	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.outputDir + "/" + _sequence + ".pdb");
	writer.write(sys.getAtomPointers(), true, false, false);
	writer.close();
}

void outputEnergyFile(MBBOptions &_opt, double _seed, map<string,map<string,double>> &_sequenceEnergyMap){
	// Setup vector to hold energy file lines
	vector<string> energyLines;
	// get the run parameters
	string t = "\t";
	stringstream enerTerms;
	// For loop to setup the energy file
	uint i = 0;
	for (auto &seq : _sequenceEnergyMap){
		stringstream seqLine;
		string sequence = seq.first;
		seqLine << sequence << t << _opt.interface << t << _seed << t;
		map<string,double> energyMap = _sequenceEnergyMap[sequence];
		// For adding in strings to a line for the energy file; looping through the terms instead of my input terms this time; sort later
		for (auto &ener: energyMap){
			string energyTerm = ener.first;
			double energy = energyMap[energyTerm];
			string term = MslTools::doubleToString(energy)+t;
			seqLine << term;
			if (i == 0){
				enerTerms << energyTerm << t;
			}
		}
		string line = seqLine.str();
		energyLines.push_back(line);
		i++;
	}
	ofstream eout;
	string eoutfile = _opt.outputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	eout << "Sequence" << t <<"Interface" << t << "Seed" << t;
	eout << enerTerms.str() << endl;
	for (uint i=0; i<energyLines.size() ; i++){
		eout << energyLines[i] << endl;
	}
	eout.close();
}

double monteCarloRepack(MBBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string, map<string,double>> &_sequenceEnergyMap,
 string _sequence, vector<uint> &_bestState, System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, double _monomerEnergy){
	// Setup backbone repack file
	ofstream bbout;
	string bboutfile  = _opt.outputDir + "/bbRepack_" + _sequence + ".out";
	bbout.open(bboutfile.c_str());

	// Local Backbone Monte Carlo Repacks Time setup	
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);
	
	// starting geometry
	double xShift = _opt.xShift;	
	double crossingAngle = _opt.crossingAngle;
	double axialRotation = _opt.axialRotation;
	double zShift = _opt.zShift;

	// final geometry
	double finalXShift = _opt.xShift;
	double finalCrossingAngle = _opt.crossingAngle;
	double finalAxialRotation = _opt.axialRotation;
	double finalZShift = _opt.zShift;

	// Get helical axis atom pointers 
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();
	
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector & apvChainB = _sys.getChain("B").getAtomPointers();

	bbout << "***STARTING GEOMETRY***" << endl;
	outputGeometry(_opt, xShift, crossingAngle, axialRotation, zShift, bbout);
	// calculate starting sasa
	SasaCalculator startSasa(_sys.getAtomPointers());
	startSasa.calcSasa();
	double sasa = startSasa.getTotalSasa();
	_sequenceEnergyMap[_sequence]["PreBBOptimizeSasa"] = sasa;
	cout << "Pre BBOptimize SASA: " << sasa << endl;

	// Monte Carlo Repack Manager Setup
	MonteCarloManager MCMngr(_opt.MCStartTemp, _opt.MCEndTemp, _opt.MCCycles, _opt.MCCurve, _opt.MCMaxRejects, _opt.convergedSteps, _opt.convergedE);
	//MonteCarloManager MCMngr(_opt.backboneMCStartTemp, _opt.backboneMCEndTemp, _opt.backboneMCCycles, _opt.backboneMCCurve, _opt.backboneMCMaxRejects);

	vector<unsigned int> MCOBest = _bestState;
	
	unsigned int counter = 0;
	_sys.setActiveRotamers(_bestState);
	double currentEnergy = _spm.getStateEnergy(_bestState)-_monomerEnergy;
	double dimer = _spm.getStateEnergy(_bestState);
	double calcDimer = _sys.calcEnergy();
	cout << "Starting Energy: " << currentEnergy << endl;
	cout << "Starting Dimer Energy: " << dimer << endl;
	cout << "Monomer Energy: " << _monomerEnergy << endl;
	cout << "Calculated Dimer Energy: " << calcDimer << endl;
	_sequenceEnergyMap[_sequence]["DimerPreBBOptimize"] = dimer;
	_sequenceEnergyMap[_sequence]["TotalPreBBOptimize"] = currentEnergy;
	_sequenceEnergyMap[_sequence]["VDWDimerPreBBOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_VDW");
	_sequenceEnergyMap[_sequence]["IMM1DimerPreBBOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1")+_spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	_sequenceEnergyMap[_sequence]["HBONDDimerPreBBOptimize"] = _spm.getStateEnergy(_bestState, "SCWRL4_HBOND");
	
	double bestEnergy = currentEnergy;
	double prevBestEnergy = currentEnergy;
	MCMngr.setEner(currentEnergy);

	// setup variables for shifts: ensures that they start from the proper values for every repack and not just the final value from the initial repack
	bool decreaseMoveSize = _opt.decreaseMoveSize;
	double deltaX = _opt.deltaX;
	double deltaCross = _opt.deltaCross;
	double deltaAx = _opt.deltaAx;
	double deltaZ = _opt.deltaZ;

	bbout << "Starting Repack Cycles" << endl;
	while(!MCMngr.getComplete()) {
		// get the current temperature of the MC
		double startTemp = MCMngr.getCurrentT();
		
		_sys.applySavedCoor("savedRepackState");
		_helicalAxis.applySavedCoor("BestRepack");
		
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
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaAxialRotation, moveToPreform);
		} else if (moveToPreform == 2) {
		//==================================
		//====== Local Crossing Angle ======
		//==================================
			//deltaCrossingAngle = getStandardNormal(_RNG1) * 1.0;
			deltaCrossingAngle = getStandardNormal(_RNG) * deltaCross;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaCrossingAngle, moveToPreform);
		} else if (moveToPreform == 3) {
		//==============================================
		//====== X shift (Interhelical Distance) =======
		//==============================================
			//deltaXShift = getStandardNormal(_RNG1) * 0.1;
			deltaXShift = getStandardNormal(_RNG) * deltaX;
			backboneMovement(apvChainA, apvChainB, axisA, axisB, _trans, deltaXShift, moveToPreform);
		}
		// Run optimization
		repackSideChains(_spm, _opt.greedyCycles);
		vector<unsigned int> MCOFinal = _spm.getMinStates()[0];
		currentEnergy = _spm.getMinBound()[0]-_monomerEnergy;
		_sys.setActiveRotamers(MCOFinal);//THIS WAS NOT HERE BEFORE 2022-8-26 NIGHT! MAKE SURE IT'S IN ALL OTHER CODE, IT'S CRUCIAL TO SAVING THE STATE

		if (counter == 0){
			_sequenceEnergyMap[_sequence]["firstRepackEnergy"] = currentEnergy;
		}

		if (!MCMngr.accept(currentEnergy)) {
			bbout << "MCReject   xShift: " << finalXShift+deltaXShift << " crossingAngle: " << finalCrossingAngle+deltaCrossingAngle << " axialRotation: " << finalAxialRotation+deltaAxialRotation << " zShift: " << finalZShift+deltaZShift << " energy: " << currentEnergy << endl;
		} else {
			bestEnergy = currentEnergy;
			_sys.saveAltCoor("savedRepackState");
			_helicalAxis.saveAltCoor("BestRepack");
		
			finalXShift = finalXShift + deltaXShift;
			finalCrossingAngle = finalCrossingAngle + deltaCrossingAngle;
			finalAxialRotation = finalAxialRotation + deltaAxialRotation;
			finalZShift = finalZShift + deltaZShift;
			MCOBest = MCOFinal;
			
			// if accept, decrease the value of the moves by the sigmoid function
			if (_opt.decreaseMoveSize == true){
				double endTemp = MCMngr.getCurrentT();
				getCurrentMoveSizes(_opt, startTemp, endTemp, deltaX, deltaCross, deltaAx, deltaZ, decreaseMoveSize);
			}
			//cout << "deltaX: " << deltaX << " deltaCross: " << deltaCross << " deltaAx: " << deltaAx << " deltaZ: " << deltaZ << endl;
			bbout << "MCAccept " << counter <<  " xShift: " << finalXShift << " crossingAngle: " << finalCrossingAngle << " axialRotation: " << finalAxialRotation << " zShift: " << finalZShift << " energy: " << currentEnergy << endl;
			counter++;
			//writer.write(_sys.getAtomPointers(), true, false, true);
		}
	}
	bbout << "End Repack Cycles" << endl << endl; 
	//writer.close();
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	_bestState = MCOBest;	
	_sys.applySavedCoor("savedRepackState");
	double dimerEnergy = _spm.getStateEnergy(MCOBest);
	double finalEnergy = dimerEnergy-_monomerEnergy;
	
	// Output change in geometry
	bbout << "***REPACK GEOMETRY***" << endl;
	outputGeometry(_opt, finalXShift, finalCrossingAngle, finalAxialRotation, finalZShift, bbout);
	bbout << "Energy;        Before: " << prevBestEnergy << "; After: " << bestEnergy << endl << endl;

	SasaCalculator endDimerSasa(_sys.getAtomPointers());
	endDimerSasa.calcSasa();
	double endSasa = endDimerSasa.getTotalSasa();
	_sequenceEnergyMap[_sequence]["BBOptimizeSasa"] = endSasa;
	cout << "BBOptimize SASA: " << endSasa << endl;
	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;
	_sequenceEnergyMap[_sequence]["VDWDimerBBOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_VDW");
	_sequenceEnergyMap[_sequence]["IMM1DimerBBOptimize"] = _spm.getStateEnergy(_bestState, "CHARMM_IMM1")+_spm.getStateEnergy(_bestState, "CHARMM_IMM1REF");
	_sequenceEnergyMap[_sequence]["HBONDDimerBBOptimize"] = _spm.getStateEnergy(_bestState, "SCWRL4_HBOND");

	// sets the updated backbone parameters
	_sequenceEnergyMap[_sequence]["endXShift"] = finalXShift;
	_sequenceEnergyMap[_sequence]["endCrossingAngle"] = finalCrossingAngle;
	_sequenceEnergyMap[_sequence]["endAxialRotation"] = finalAxialRotation;
	_sequenceEnergyMap[_sequence]["endZShift"] = finalZShift;
    double relativeAx; 
    double relativeZ; 
	convertToRelativeAxAndZ(_opt.axialRotation, _opt.zShift, relativeAx, relativeZ);
	_sequenceEnergyMap[_sequence]["endAxialRotationPrime"] = relativeAx;
	_sequenceEnergyMap[_sequence]["endZShiftPrime"] = relativeZ;
	bbout << MCMngr.getReasonCompleted() << endl;	
	bbout << "Monte Carlo repack complete. Time: " << diffTimeMC/60 << "min" << endl << endl;
	//TODO: there may be a better way to resolve this, but as of 2022-9-8, I want to get the most data I can before a lab meeting, so putting this here
	if (finalEnergy > 100){
        bbout << "Final energy is " << finalEnergy << " after repack, indicating clashes. Choose a different geometry" << endl;
		exit(0);
	}
	return finalEnergy;
}

void transformHelicalAxis(AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {

	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);

	//===== Axial Rotation ======
	_trans.rotate(_axisA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_axisA, _axisB);
}


void convertToRelativeAxAndZ(double _axialRot, double _zShift, double &_relativeAx, double &_relativeZ){
	// use the positive axial rotation for conversion since our axial rotations are negative
	//double axialRotation = abs(_axialRot);
	double axialRotation = _axialRot;
	// use equations from Mueller 2014; Fig. S1
	_relativeAx = (10*axialRotation/9)+(200*_zShift/27);
	_relativeAx = 100+_relativeAx;
	_relativeZ = (10*_zShift/9)+(0.15*axialRotation/9);
}

void convertToAxAndZForTranformation(MBBOptions &_opt){
	double axialRotation = abs(_opt.axialRotation);
	double zShift = _opt.zShift;
	// convert input axialRotation and zShift to transformation for interfacial parallelogram (Mueller 2014; Fig S1) 
	_opt.axialRotation = axialRotation+(20*zShift/3);
	_opt.axialRotation = -_opt.axialRotation;
	_opt.zShift = zShift+(0.015*axialRotation);
}

void computeMonomerEnergy(System &_sys, System &_helicalAxis, MBBOptions &_opt, Transforms &_trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG) {

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);
	CSBMono.setBuildTerm("CHARMM_IMM1REF", true);
	CSBMono.setBuildTerm("CHARMM_IMM1", true);

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
	monoEset->setAllTermsInactive();
	monoEset->setTermActive("CHARMM_IMM1REF", true);
	monoEset->setTermActive("CHARMM_IMM1", true);
	monoEset->setTermActive("CHARMM_VDW", true);
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	int firstPos = 0;
    int lastPos = monoSys.positionSize();
    deleteTerminalBondInteractions(monoSys,_opt,firstPos,lastPos);

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);
	CSBMono.updateNonBonded(10,12,50);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.updateWeights();
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	_trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), _trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	_trans.translate(chainA, interDistVect);

	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	_trans.translate(chainA, move5Down);
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
	_helicalAxis.saveAltCoor("BestMonomerAxis");
	//_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		//double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_opt.greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		//_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			//_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	_helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		_helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		_trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		_trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		for(int j=0; j<=3; j++) { // test at 4 rotations 0, 90, 180 and 270 degrees
			//==================================
			//====== Axial Rot ======
			//==================================
			monoAxialRotation = j * 90.0;

			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			//_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			//monoSys.writePdb("mono_" + MslTools::doubleToString(monoTilt) + "_" + MslTools::doubleToString(monoAxialRotation) + ".pdb");

			if(currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				bestTilt = monoTilt;
				bestRotation = monoAxialRotation;
				monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
				monoSys.saveAltCoor("savedBestMonomer");
				_helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			_trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

		}
	}

	//MonteCarloManager MCMngr(1000.0, 0.5, _opt.MCCycles, MonteCarloManager::EXPONENTIAL, _opt.MCMaxRejects);
	MonteCarloManager MCMngr(0.5, 0.5, 100, MonteCarloManager::EXPONENTIAL, 5);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		_helicalAxis.applySavedCoor("BestMonomerAxis");

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
			_trans.translate(chainA, translateA);
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			_trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			_trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			_trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_fout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_opt.greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_fout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_fout << "state rejected   energy: " << currentEnergy << endl;
		} else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;
		}
		counter++;
	}

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	//Calculate Monomer energy for output
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getMinBound()[0]*2;

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	monoSasa.calcSasa();
	double monomerSasa = monoSasa.getTotalSasa();
	double totalMonomerSasa = monomerSasa*2;

	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "Monomer", true);
	_sequenceEnergyMap[_seq]["Monomer"] = monomerEnergy;
	_sequenceEnergyMap[_seq]["MonomerSasa"] = totalMonomerSasa;

	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	_helicalAxis.clearSavedCoor("BestMonomerAxis");
	_helicalAxis.clearSavedCoor("bestZ");
}

// help functions
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputErrorMessage(MBBOptions &_opt){
	cout << endl;
	cout << "The program terminated with errors:" << endl;
	cout << endl;
	cout << _opt.errorMessages << endl;
	cout << endl;
	cout << _opt.OPerrors << endl;
	usage();
}


MBBOptions parseMBBOptions(int _argc, char * _argv[]){

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a BBOptions structure
	 *  defined at the head of this file
	 ******************************************/
	MBBOptions opt;

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
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("outputDir");
	opt.allowed.push_back("inputPdbFile");
	opt.allowed.push_back("configfile");

	opt.allowed.push_back("seed");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("negAngle");
	opt.allowed.push_back("negRot");
	opt.allowed.push_back("thread");
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
	opt.allowed.push_back("convergedSteps");
	opt.allowed.push_back("convergedE");

	opt.allowed.push_back("useElec");
	opt.allowed.push_back("helicalAxis");
	opt.allowed.push_back("useAlaAtTermini");
	opt.allowed.push_back("energyTermList");
	opt.allowed.push_back("deleteTerminalInteractions");
	
	//Rotamers
	opt.allowed.push_back("SL");
	opt.allowed.push_back("Ids");
    opt.allowed.push_back("interface");

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
	opt.interface = OP.getString("interface");//TODO: in the future, add in a way to get the interface from SASA or some other method
	if(OP.fail()) {
		opt.errorMessages += "interface not specified\n";
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
	opt.convergedSteps = OP.getInt("convergedSteps");
	if (OP.fail()) {
		opt.warningMessages += "convergedSteps not specified using MCCycles value\n";
		opt.warningFlag = true;
		opt.convergedSteps = opt.MCCycles;
	}
	opt.convergedE = OP.getInt("convergedE");
	if (OP.fail()) {
		opt.warningMessages += "convergedE not specified using 0.001\n";
		opt.warningFlag = true;
		opt.convergedE = 0.001;
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
	
    opt.inputPdbFile = OP.getString("inputPdbFile");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine inputPdbFile\n";
		opt.errorFlag = true;
	}

	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine outputDir, using current directory\n";
		opt.warningFlag = true;
	}

	// Monomer Options
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 0\n";
		opt.warningFlag = true;
		opt.seed = 0;
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
	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.errorMessages += "helicalAxis not specified\n";
		opt.errorFlag = true;
	}
	opt.useAlaAtTermini = OP.getBool("useAlaAtTermini");
	if (OP.fail()) {
		opt.warningMessages += "useAlaAtTermini not specified using true\n";
		opt.warningFlag = true;
		opt.useAlaAtTermini = true;
	}
	opt.deleteTerminalInteractions = OP.getMultiString("deleteTerminalInteractions");
	if (OP.fail()) {
		opt.deleteTerminalInteractions.push_back("");
		opt.warningMessages += "deleteTerminalInteractions not specified\n";
		opt.warningFlag = true;
	}
	opt.energyTermList = OP.getStringVector("energyTermList");
	if (OP.fail()) {
		//This works, but I think if you ever want to output more terms in the future, need to add them to the terms above
		//TODO: write in an error that will tell you if the above is the case
		opt.energyTermList.push_back("CHARMM_VDW");
		opt.energyTermList.push_back("SCWRL4_HBOND");
		opt.energyTermList.push_back("CHARMM_IMM1");
		opt.energyTermList.push_back("CHARMM_IMM1REF");
        opt.warningMessages += "energyTermList not specified, defaulting to CHARMM_VDW, SCWRL4_HBOND, CHARMM_IMM1, CHARMM_IMM1REF\n";
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
	// alternate identities
	opt.Ids = OP.getStringVector("Ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.rerunConf = OP.getConfFile();
	return opt;
}

void setupDesignDirectory(MBBOptions &_opt){
	_opt.outputDir = string(get_current_dir_name()) + "/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}
void outputGeometry(MBBOptions &_opt, double _xShift, double _crossingAngle, double _axialRotation, double _zShift, ofstream &_out){
	_out << "xShift;        Before: " << _opt.xShift << "; After: " << _xShift << endl;
	_out << "crossingAngle; Before: " << _opt.crossingAngle << "; After: " << _crossingAngle << endl;
	_out << "axialRotation; Before: " << _opt.axialRotation << "; After: " << _axialRotation << endl;
	_out << "zShift;        Before: " << _opt.zShift << "; After: " << _zShift << endl << endl;
	double relativeAxBefore;
	double relativeZBefore;
	convertToRelativeAxAndZ(_opt.axialRotation, _opt.zShift, relativeAxBefore, relativeZBefore);
	double relativeAxAfter;
	double relativeZAfter;
	convertToRelativeAxAndZ(_axialRotation, _zShift, relativeAxAfter, relativeZAfter);
	// convert axialRotation and zShift for output (from parallelogram from Ben's paper for interface to square)
	_out << "axialRotationPrime; Before: " << relativeAxBefore << "; After: " << relativeAxAfter << endl;
	_out << "zShiftPrime;        Before: " << relativeZBefore << "; After: " << relativeZAfter << endl << endl;
}

void getCurrentMoveSizes(MBBOptions &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize) {
	double decreaseMultiplier = _endTemp/_currTemp;
	bool decreaseX = true;
	bool decreaseCross = true;
	bool decreaseAx = true;
	bool decreaseZ = true;
	_deltaX = decreaseMoveSize(_deltaX, _opt.deltaXLimit, decreaseMultiplier, decreaseX);
	_deltaCross = decreaseMoveSize(_deltaCross, _opt.deltaCrossLimit, decreaseMultiplier, decreaseCross);
	_deltaAx = decreaseMoveSize(_deltaAx, _opt.deltaAxLimit, decreaseMultiplier, decreaseAx);
	_deltaZ = decreaseMoveSize(_deltaZ, _opt.deltaZLimit, decreaseMultiplier, decreaseZ);
	if (decreaseX == false && decreaseCross == false && decreaseAx == false && decreaseZ == false){
		_decreaseMoveSize = false;
	}
}

void deleteTerminalBondInteractions(System &_sys, MBBOptions &_opt, int _firstResiNum, int _lastResiNum){
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

void prepareSystem(MBBOptions &_opt, System &_sys, System &_startGeom, PolymerSequence &_PS){	
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
	CSB.setBuildTerm("CHARMM_IMM1REF", true);
	CSB.setBuildTerm("CHARMM_IMM1", true);

	// load the membrane as solvent
	CSB.setSolvent("MEMBRANE");
	// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	CSB.setIMM1Params(15, 10);
	// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	CSB.setBuildNonBondedInteractions(false);

	// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	if(!CSB.buildSystem(_PS)) {
		cerr << "Unable to build system from " << _PS << endl;
		exit(0);
	}
	
	// assign the coordinates of our system to the given geometry 
	_sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	_sys.buildAllAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(_sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Initialize EnergySet that contains energies for the chosen terms for our design
	EnergySet* Eset = _sys.getEnergySet();
	// Set all terms inactive and explicitly set the given terms as active
	Eset->setAllTermsInactive();
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
	int firstPos = 0;// should this be based on the thread and the
    int lastPos = _sys.positionSize();
    deleteTerminalBondInteractions(_sys,_opt,firstPos,lastPos);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	CSB.updateNonBonded(10,12,50);
}
