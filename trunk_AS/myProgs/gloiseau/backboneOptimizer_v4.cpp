#include <sstream>
#include <iterator>
#include <unistd.h>
#include <thread>

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
#include "backboneOptimizerFunctions_v2.h"
#include "backboneOptimizerOptions.h"
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "backboneOptimizer_v3";
string programDescription = "This is an updated version of backboneOptimizer: backboneOptimizer does local backbone moves on given backbones\n\
 	from a given input geometry. Functions in this are more improved, and it uses multithreading to try multiple threaded positions.\n\
	It also acts as a forward folding program, where it will try to minimize the energy of a given sequence and geometry.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "4";
string programDate = "30 August 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

// BBOptions Functions
void prepareSystem(BBOptions &_opt, System &_sys, System &_startGeom, string &_polySeq);
void setGly69ToStartingGeometry(BBOptions &_opt, System &_sys, System &_helicalAxis, Transforms &_trans, double _crossingAngle, double _xShift);
void setGly69ToStartingGeometry(BBOptions &_opt, map<string,double> _geometries, System &_sys, System &_helicalAxis, Transforms &_trans);
void threadSearchBackboneOptimization(BBOptions &_opt, System &_helicalAxis, int _thread, map<string,double> _geometries, vector<uint> _stateVector,
 map<string,double> &_monomerEnergyByTerm, double _monomerEnergy, ofstream &_sout, ofstream &_eout);
void defineRotamerLevels(System &_sys, BBOptions &_opt, vector<pair <int, double>> _resiBurial, 
 vector<uint> &_interfacePositions, string &_variablePositionString, string &_rotamerLevels);
void defineInterfaceAndRotamerSampling(BBOptions &_opt, System &_startGeom, PolymerSequence _PS, string &_rotamerLevels,
 string &_polySeq, string &_variablePositionString, string &_rotamerSamplingString, vector<int> &_linkedPositions, 
 vector<uint> &_allInterfacePositions, vector<uint> &_interfacePositions, vector<int> &_rotamerSamplingPerPosition, ofstream &_out);
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, BBOptions &_opt, string _seq);
vector<uint> getAllInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition);
vector<uint> getInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition);
void localXShiftDocking(System &_sys, BBOptions &_opt, double &_bestEnergy, double _monomerEnergy, SelfPairManager &_spm, 
 System &_helicalAxis, Transforms &_trans, int _thread, double &_savedXShift, ofstream &_out);
void localBackboneRepack(BBOptions &_opt, System &_sys, map<string,double> _geometries, SelfPairManager &_spm,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm,
 double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout);
void monteCarloRepack(BBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string,double> _geometries,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, double _prevBestEnergy,
 map<string,double> _monomerEnergyByTerm, double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout);
void checkOptionErrors(BBOptions &_opt);	
void getCurrentMoveSizes(BBOptions &_opt, double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 bool &_decreaseMoveSize);
double increaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);
double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease);
void threadDocking(BBOptions &_opt, System &_helicalAxis, int _thread, double _crossingAngle, double _monomerEnergy, 
map<double, map<double, map<string, double>>> &_threadGeometries, map<double, map<double, vector<uint>>> &_threadStateVector,
map<string,double> _monomerEnergyByTerm, ofstream &_eout);
void getAxialRotAndZShift(BBOptions &_opt, map<string,double> &_geometries);

/***********************************
 *help functions
 ***********************************/
void usage();
void help(BBOptions defaults);
void outputErrorMessage(BBOptions &_opt);
void outputWarningMessage(BBOptions &_opt);

int main(int argc, char *argv[]){
	
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	// Start the timer
	time(&startTime);
	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);

	string date(buffer);
    
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	//Add in some default options that can easily be changed here
	BBOptions opt = BBParseOptions(argc, argv);
	checkOptionErrors(opt);

	/******************************************************************************
	 *                       === SETUP OUTPUTS ===
	 ******************************************************************************/
	ofstream sout;  // summary file output
	ofstream eout;  // summary file output
	ofstream mout;  // monomer file output
	ofstream err;   // error file output
	ofstream rerun; // rerun config output

	// function that defines the output directory
	setupOutputDirectoryChtc(opt);

	// setup the output files
	string soutfile = opt.outputDir + "/summary.out";
	string eoutfile = opt.outputDir + "/allEnergies.csv";
	string moutfile = opt.outputDir + "/monomer.out";
	string errfile  = opt.outputDir + "/errors.out";
	string rerunfile = opt.outputDir + "/rerun.config";

	// open the output files
	sout.open(soutfile.c_str());
	eout.open(eoutfile.c_str());
	mout.open(moutfile.c_str());
	err.open(errfile.c_str());
	rerun.open(rerunfile.c_str());

	// write and close the rerun config file
	rerun << opt.rerunConf << endl;
	rerun.close();

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();
	helicalAxis.saveCoor("originState");
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	// threading for docking at all threads
	map<double, map<double, map<string, double>>> threadGeometries;
	map<double, map<double, vector<uint>>> threadStateVector;
	
	PDBWriter writer;
	/******************************************************************************
	 *                    === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(opt,startGeom,helicalAxis,trans,opt.crossingAngle,opt.xShift);
	helicalAxis.applySavedCoor("originState");
	System sys;
	string polySeq = convertToPolymerSequence(opt.sequence, opt.thread);
	prepareSystem(opt, sys, startGeom, polySeq);
	RandomNumberGenerator RNG;
	RNG.setSeed(opt.seed); 
    map<string,double> monomerEnergyByTerm;
    double monomerEnergy = computeMonomerEnergy(sys, opt, RNG, monomerEnergyByTerm, mout);
    //double monomerEnergy = computeMonomerEnergy1(opt, startGeom, 25, monomerEnergyByTerm, mout);
	for (uint i=opt.threadStart; i<opt.threadEnd; i++){
		vector<thread> dockThreads;
		//cout << "Docking Thread: " << i << endl;
		for (uint j=0; j<opt.crossAngle.size();j++){
			double crossingAngle = opt.crossAngle[j]; 
			//cout << "Angle: " << crossingAngle << endl;
			dockThreads.push_back(thread(threadDocking, ref(opt), ref(helicalAxis), i, crossingAngle, monomerEnergy,
			 ref(threadGeometries), ref(threadStateVector), ref(monomerEnergyByTerm), ref(eout)));
		}
		for (auto &t : dockThreads){
			t.join();
		}
	}
	// use this instead, rid of anything with too high of an energy before here to trim how much of the backbone space to search here
	//for (auto &threads : threadGeometries){
	//	int thread = threads.first;
	//	map<double, map<string, double>> &crossingAngles = threads.second;
	//	for (auto &crossingAngle : crossingAngles){
	//		double angle = crossingAngle.first;
	//		map<string, double> &geometries = crossingAngle.second;
	//		vector<uint> stateVector = threadStateVector[thread][crossingAngle];
	//		cout << v.first << " " << v.second << endl;
	//	}
	//}
	//for (uint i=opt.threadStart; i<opt.threadEnd; i++){
	//	vector<thread> threads;
	//	cout << "Repack Thread: " << i << endl;
	//	for (uint j=0; j<opt.crossAngle.size(); j++){
	//		double crossingAngle = opt.crossAngle[j];
	//		cout << "Angle: " << crossingAngle << endl;
	//		map<string,double> geometries = threadGeometries[i][crossingAngle];
	//		cout << "Geometries: " << endl;
	//		vector<uint> stateVec = threadStateVector[i][crossingAngle];
	//		for (auto &g : geometries){
	//			cout << g.first << ": " << g.second << endl;
	//		}
	//		threads.push_back(thread(threadSearchBackboneOptimization, ref(opt), ref(helicalAxis), i, geometries, stateVec, ref(monomerEnergyByTerm), monomerEnergy, ref(sout), ref(eout)));	
	//	}
	//	for (auto& th: threads){
	//			th.join();
	//	}
	//}
	
	// output the total time
	time(&endTime);
	diffTime = difftime (endTime, startTime);
	sout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
	cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime/60 << " minutes" << endl;
    // close all of the output file writers
    mout.close();
    err.close();
	sout.close();
	eout.close();
}

//Functions
//TODO: to add this into the design code, I think I just need to add in a mask here to the optimizer
// set a number of cycles
void getAxialRotAndZShift(BBOptions &_opt, map<string,double> &_geometries){
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 
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
	//double axialRotationDensity = MslTools::toDouble(tokens[5]);
	//double zShiftDensity = MslTools::toDouble(tokens[6]);
}

void threadDocking(BBOptions &_opt, System &_helicalAxis, int _thread, double _crossingAngle, double _monomerEnergy, 
 map<double, map<double, map<string, double>>> &_threadGeometries, map<double, map<double, vector<uint>>> &_threadStateVector,
 map<string,double> _monomerEnergyByTerm, ofstream &_eout){
	// setup the output file
	ofstream sout;  // summary file output
	string soutfile = _opt.outputDir + "/repackSummary_" + to_string(_thread) + "_" + to_string(_crossingAngle) + ".out";
	sout.open(soutfile.c_str());
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	map<string,double> geometries;

	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);

	if (_opt.getRandomAxAndZ){
		getAxialRotAndZShift(_opt, geometries);
	} else {
		geometries["axialRotation"] = _opt.axialRotation;
		geometries["zShift"] = _opt.zShift;
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

	Transforms trans1;
	trans1.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans1.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	/******************************************************************************
	 *                       === SETUP STARTING GEOMETRY ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	string polySeq = convertToPolymerSequence(_opt.sequence, _thread);
	prepareSystem(_opt, sys, startGeom, polySeq);
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?

	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Assign number of rotamers by residue burial
	loadRotamers(sys, sysRot, _opt.SL);
	//sout << "Number of rotamers: " << sysRot.getNumberOfRotamers() << endl;
	sout << "SL: " << _opt.SL << endl;
	
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 
	
	// _optimize Initial Starting Position 
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	//spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);// changed to make more similar to CATM
	spm.saveEnergiesByTerm(true);
	//spm.calculateEnergies();
	
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
	
	double savedXShift = _opt.xShift;
	localXShiftDocking(sys, _opt, totalEnergy, _monomerEnergy, spm, helicalAxis, trans, _thread, savedXShift, sout);
	//localCrossingAngleSearch(sys, _opt, totalEnergy, _monomerEnergy, spm, _helicalAxis, trans, _thread, savedXShift);
	//localAxialRotSearch(sys, _opt, totalEnergy, _monomerEnergy, spm, _helicalAxis, trans, _thread, savedXShift);
	//localZShiftSearch(sys, _opt, totalEnergy, _monomerEnergy, spm, _helicalAxis, trans, _thread, savedXShift);
	//moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?

	geometries["xShift"] = savedXShift;
	double xShift = geometries["xShift"];
	double zShift = geometries["zShift"];
	double crossingAngle = geometries["crossingAngle"];
	double axialRotation = geometries["axialRotation"];
	if (totalEnergy > _opt.energyCutoff){
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is higher than cutoff (100), skip" << endl;
	} else {
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is lower than cutoff (100), continue" << endl;
		double bestEnergy = startDimer-_monomerEnergy;

		/******************************************************************************
		 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
		 ******************************************************************************/
		sout << "Starting Geometry" << endl;
		sout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << endl << endl;
		sout << "Current Best Energy: " << bestEnergy-_monomerEnergy << endl;

		double bestRepackEnergy = sys.calcEnergy();
		vector<uint> bestRepackState;

		// Local backbone monte carlo repacks
		localBackboneRepack(_opt, sys, geometries, spm, helicalAxis, trans, RNG, _monomerEnergyByTerm, _monomerEnergy, _thread, sout, _eout);
	}
	sout.close();
}

void threadSearchBackboneOptimization(BBOptions &_opt, System &_helicalAxis, int _thread, map<string,double> _geometries, vector<uint> _stateVector,
 map<string,double> &_monomerEnergyByTerm, double _monomerEnergy, ofstream &_sout, ofstream &_eout){
	
	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	System helicalAxis;
	helicalAxis.readPdb(_opt.helicalAxis);
	// get the starting geometry using poly glycine
	System startGeom;
	setGly69ToStartingGeometry(_opt,_geometries,startGeom,helicalAxis,trans);
	
	double xShift = _geometries["xShift"];
	double zShift = _geometries["zShift"];
	double crossingAngle = _geometries["crossingAngle"];
	double axialRotation = _geometries["axialRotation"];
	
	// setup the output file
	ofstream sout;  // summary file output
	string soutfile = _opt.outputDir + "/backboneOptimizerSummary_" + to_string(_thread) + "_" + to_string(crossingAngle) + ".out";
	sout.open(soutfile.c_str());
	sout << "Geometry: " << endl;
	sout << " -xShift: " << xShift << endl;
	sout << " -zShift: " << zShift << endl;
	sout << " -crossingAngle: " << crossingAngle << endl;
	sout << " -axialRotation: " << axialRotation << endl << endl;

	/******************************************************************************
	 *                       === SETUP STARTING GEOMETRY ===
	 ******************************************************************************/
	// set up the system for the input sequence
	System sys;
	string polySeq = convertToPolymerSequence(_opt.sequence, _thread);
	prepareSystem(_opt, sys, startGeom, polySeq);
	moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), _helicalAxis.getAtomPointers(), trans);//compared to CATM, my structures were moved up by like 4 AAs. Could it be because of this?

	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.outputDir + "/geom_" + to_string(_thread) + "_" + to_string(crossingAngle) + "_" + to_string(xShift) + ".pdb");
	writer.write(sys.getAtomPointers(), true, false, true);
	writer.close();

	/******************************************************************************
	 *       === IDENTIFY INTERFACIAL POSITIONS AND GET ROTAMER ASSIGNMENTS ===
	 ******************************************************************************/
	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Assign number of rotamers by residue burial
	loadRotamers(sys, sysRot, _opt.SL);
	
	// setup random number generator object
	RandomNumberGenerator RNG;
	RNG.setSeed(_opt.seed); 

	sys.buildAtoms();
	// _optimize Initial Starting Position 
	SelfPairManager spm;
	spm.seed(RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	//spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);// changed to make more similar to CATM
	spm.saveEnergiesByTerm(true);
	//spm.calculateEnergies();
	
	// Repack dimer
	//repackSideChains(spm, _opt.greedyCycles);
	//vector<uint> startStateVec = spm.getMinStates()[0];
	//double currentEnergy = spm.getMinBound()[0];
	sys.setActiveRotamers(_stateVector);

	// calculate the energy of the system
	double startDimer = sys.calcEnergy();
	double totalEnergy = startDimer-_monomerEnergy;
	sout << "Energies " << _thread << endl;
	sout << " -Dimer Energy: " << startDimer << endl;
	sout << " -Monomer Energy: " << _monomerEnergy << endl;
	sout << " -Total Energy: " << totalEnergy << endl;
	sout << spm.getSummary(_stateVector) << endl;
	if (totalEnergy > _opt.energyCutoff){
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is higher than cutoff (100), skip" << endl;
	} else {
		sout << "Thread " << _thread << " energy " << totalEnergy << " at crossingAngle " << crossingAngle << " and xShift " << xShift << " is lower than cutoff (100), continue" << endl;
	
		sys.saveAltCoor("startingState");
		//_helicalAxis.saveAltCoor("startingAxis");
		sys.saveAltCoor("savedBestState");
		//_helicalAxis.saveAltCoor("BestAxis");
	
		/******************************************************************************
		 *                     === X SHIFT REPACKS ===
		 ******************************************************************************/
		double bestEnergy = startDimer-_monomerEnergy;
		double savedXShift = xShift;

		/******************************************************************************
		 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
		 ******************************************************************************/
		sout << "Starting Geometry" << endl;
		sout << setiosflags(ios::fixed) << setprecision(3) << "xShift: " << savedXShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << endl << endl;
		sout << "Current Best Energy: " << bestEnergy-_monomerEnergy << endl;

		double bestRepackEnergy = sys.calcEnergy();
		vector<uint> bestRepackState;

		// Local backbone monte carlo repacks
		localBackboneRepack(_opt, sys, _geometries, spm, helicalAxis, trans, RNG, _monomerEnergyByTerm, _monomerEnergy, _thread, sout, _eout);
	}
	sout.close();
}

void localBackboneRepack(BBOptions &_opt, System &_sys, map<string,double> _geometries, SelfPairManager &_spm,
 System &_helicalAxis, Transforms &_trans, RandomNumberGenerator &_RNG, map<string,double> _monomerEnergyByTerm,
 double _monomerEnergy, int _thread, ofstream &_out, ofstream &_eout){
	_out << "==============================================" << endl;
	_out << " Performing Local Monte Carlo Backbone Repack " << endl;
	_out << "==============================================" << endl;
	// set the system to the original xShift state
	_sys.applySavedCoor("savedBestState");
	//_helicalAxis.applySavedCoor("BestAxis");
	// save the current state as the repack state
	_sys.saveAltCoor("savedRepackState");
	//_helicalAxis.saveAltCoor("BestRepack");
	// get the best repack energy
	double prevBestEnergy = _sys.calcEnergy();

	// do backbone geometry repacks
	monteCarloRepack(_opt, _sys, _spm, _geometries, _helicalAxis, _trans, _RNG, prevBestEnergy, 
	 _monomerEnergyByTerm, _monomerEnergy, _thread, _out, _eout);
  
    //outputs a pdb file for the structure 
	// Initialize PDBWriter
	PDBWriter writer;
	writer.open(_opt.outputDir + "/backboneOptimized_" + to_string(_thread) + ".pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
}

void monteCarloRepack(BBOptions &_opt, System &_sys, SelfPairManager &_spm, map<string,double> _geometries,
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
	bool decreaseMoveSize = _opt.decreaseMoveSize;
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
	
	if (_opt.useElec == false){
		_out << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		_out << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		_out << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		_out << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		_out << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		_out << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		_out << MCMngr.getReasonCompleted() << endl;
		//_eout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		//_eout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		//_eout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		//_eout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		//_eout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		//_eout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		//_eout << MCMngr.getReasonCompleted() << endl;
		sout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		sout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		sout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',';
		sout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		sout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		sout << MCMngr.getReasonCompleted() << endl;
		_eout << finalEnergy << ',' << _thread << ',' << xShift << ',' << crossingAngle << endl;
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
		//_eout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, elec, monoElec, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		//_eout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		//_eout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		//_eout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		//_eout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		//_eout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		//_eout << MCMngr.getReasonCompleted() << endl;
		sout << "Sequence, thread, Energy, Dimer, Monomer, SASA, vdw, monoVdw, hbond, monoHbond, imm1, monoImm1, imm1Ref, monoImm1Ref, elec, monoElec, startXShift, finalXShift, startCrossingAngle, finalCrossingAngle, startAxialRot, finalAxialRot, startZShift, finalZShift" << endl;
		sout << _opt.sequence << ',' << _thread << ',' << finalEnergy << ',' << dimerEnergy << ',' << _monomerEnergy << ',';
		sout << dimerSasa << ',' << vdw << ',' << monomerVdw << ',' << hbond << ',' << monomerHbond << ',' << imm1 << ',' << monomerImm1 << ',' << imm1Ref << ',' << monomerImm1Ref << ',' << elec << ',' << monomerElec << ',';
		sout << _opt.xShift << ',' << xShift << ',' << _opt.crossingAngle << ',' << crossingAngle << ',';
		sout << _opt.axialRotation << ',' << axialRotation << ',' << _opt.zShift << ',' << zShift << endl;
		sout << "Monte Carlo repack complete. Time: " << diffTimeMC << " seconds" << endl << endl;
		sout << MCMngr.getReasonCompleted() << endl;
		_eout << finalEnergy << ',' << _thread << ',' << xShift << ',' << crossingAngle << endl;
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

double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease) {
	// edited to make sure that the move size is decreasing properly down to the move limit: add in detail here
	double diffMoveSize = _moveSize - _moveLimit;
	double moveDecrease = diffMoveSize * _decreaseMultiplier;
	double newMoveSize = _moveSize - moveDecrease;
	if (newMoveSize > _moveLimit){
		return newMoveSize;
	} else {
		_decrease = false;
		return _moveSize;
	}
}

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

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, BBOptions &_opt, string _seq) {
	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);
	PolymerSequence PS(polySeq);

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
	}

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	System gly69;
	// initialize the gly69 backbone coordinates and transform it to the chosen geometry
	gly69.readPdb(_opt.backboneFile);

	AtomPointerVector& glyAPV = gly69.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator dimerSasa(_sys.getAtomPointers());
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	dimerSasa.calcSasa();
	monoSasa.calcSasa();
	dimerSasa.setTempFactorWithSasa(true);

	for (uint i = 0; i < monoSys.positionSize(); i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdMonomer = monoSys.getPosition(i).getPositionId();
		string posIdDimer = _sys.getPosition(i).getPositionId();
		double resiSasaMonomer = monoSasa.getResidueSasa(posIdMonomer);
		double resiSasaDimer = dimerSasa.getResidueSasa(posIdDimer);
		double burial = resiSasaDimer/resiSasaMonomer;
		residueBurial.push_back(pair<int,double>(i, burial));

		//set sasa for each residue in the b-factor
		AtomSelection selA(_sys.getPosition(i).getAtomPointers());
		AtomSelection selB(_sys.getPosition(i+_opt.sequence.length()).getAtomPointers());
		AtomPointerVector atomsA = selA.select("all");
		AtomPointerVector atomsB = selB.select("all");
		for (AtomPointerVector::iterator k=atomsA.begin(); k!=atomsA.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
		for (AtomPointerVector::iterator k=atomsB.begin(); k!=atomsB.end();k++) {
			// set the residue sasa in the b-factor
			(*k)->setTempFactor(burial);
		}
	}

	PDBWriter writer;
	writer.open(_opt.outputDir+"/interfaceSASA.pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	return residueBurial;
}

vector<uint> getAllInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=0; k<_opt.sequence.length(); k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(BBOptions &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.sequence.length(); k++){
	for (uint k=3; k<_opt.sequence.length()-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
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

		// Check if this is the lowest energvoid threadDocking(BBOptions &_opt, System &_helicalAxis, int _thread, double _crossingAngle, double _monomerEnergy, 
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

/***********************************
 *help functions
 ***********************************/
void usage() {
	cout << endl;
	cout << "Run as :" << endl;
	cout << "   % " << programName << " --configfile <file.config>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;//TODO: add in some help options
}

void outputWarningMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program has the following warning:" << endl;
		cerr << endl;
		cerr << _opt.warningMessages << endl;
		cerr << endl;
}
void outputErrorMessage(BBOptions &_opt){
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << _opt.errorMessages << endl;
		cerr << endl;
		cerr << _opt.OPerrors << endl;
		usage();
}

//TODO: finish writing up this help
void help(BBOptions defaults) {
	cout << "This program runs as:" << endl;
	cout << " % seqDesign " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --numberOfStructuresToMCRepack <int> --energyCutOff <double> --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --MCStartTemp <double> --MCEndTemp <double> --MCCurve <CONSTANT-0, LINEAR-1, EXPONENTIAL-2, SIGMOIDAL-3, SOFT-4>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --thread <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   --weight_hbond <double> --weight_vdw <double> --weight_solv <double> --weight_seqEntropy <double>" << endl;
	cout << "   --sasaRepackLevel <rotLevel> (in format SL95.00; 4 levels used by default) --interfaceLevel <int> " << endl << endl;
	cout << "Template Configuration file (copy and paste the below into a file.config and run code as bin/seqDesign --config file.config" << endl;
	cout << "#Input Files" << endl;
	cout << setw(20) << "topFile " << defaults.topFile << endl;
	cout << setw(20) << "parFile " << defaults.parFile << endl;
	cout << setw(20) << "rotLibFile " << defaults.rotLibFile << endl;
	cout << setw(20) << "solvFile " << defaults.solvFile << endl;
	cout << setw(20) << "hbondFile " << defaults.hbondFile << endl;

	cout << "#Booleans" << endl;
	cout << setw(20) << "verbose " << defaults.verbose << endl;

	cout << endl << "#Energy term weights" << endl;
	cout << setw(20) << "weight_vdw " << defaults.weight_vdw << endl;
	cout << setw(20) << "weight_hbond " << defaults.weight_hbond << endl;
	cout << setw(20) << "weight_solv " << defaults.weight_solv << endl;
	cout << endl;
}

// check through error options and exit if too many
void checkOptionErrors(BBOptions &_opt){	
	if (_opt.errorFlag) {
		outputErrorMessage(_opt);
		exit(1);
	} else if (!_opt.errorFlag && !_opt.warningFlag && _opt.errorMessages != ""){
		outputWarningMessage(_opt);
		usage();
		exit(0);
	}
}