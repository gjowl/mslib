#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
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
#include "BaselineInteraction.h"
#include "BaselineIMM1Interaction.h"
#include "BaselinePairInteraction.h"
#include "BaselineEnergyBuilder.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "baselineIMM1Comparison";
string programDescription = "Updated Version of baselineSelfPairComparison that includes IMM1 baseline: This program generates random monomer sequences and tests the baseline energies
 from generateSelfPairBaseline AND generateBaselineIMM1, aiming to see if there is a good correlation between the calculated energy scores and the scores summed from the baseline energies. 
 This code is test code for the baselines prior to using them for seqDesign. If there is good correlation, these energies can be used to estimate the monomer energy during design, allowing 
 the monomer to be taken into account without separate calculation during the design process.";
string programAuthor = "Gilbert Loiseau";
string programVersion = "6";
string programDate = "1 September 2021";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;
//TODO: edited but haven't tested
/******************************************************************************************************************************************************************************/

struct Options{
	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int thread;
	
	// input files
	string helixGeoFile;
	string topFile;
	string parFile;
	string hbondFile;
	string solvFile;
	string rotLibFile;
	string backboneCrd;	
	string selfEnergyFile;
	string imm1EnergyFile;
	string pairEnergyFile;
	string pdbOutputDir;
	string infile;

	// sequence parameters
	string backboneAA;
	int backboneLength;
	int seqNumber;
	int pairDist;
	
	// side-chain repack variable
	bool verbose;
	int greedyCycles;
	int seed;

	string SL; //number of rotamers

	// optional
	int tmStart;
	int tmEnd;
	bool deleteTerminalHbonds;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	
	// alternate identities and weights
	vector<string> ids;
	vector<int> weights;

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
	string runNumber;

	string configfile;
	string datafile;
};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
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
	return "A" + ps;
}//changed this for only one chain

string generatePolymerSequence(string _backboneAA, int _sequenceLength, int _startResNum) {
	string ps = "";
	string resName = MslTools::getThreeLetterCode(_backboneAA);
	if(resName == "HIS") {
		resName = "HSE";
	}
	for (uint i=0; i<_sequenceLength; i++){
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {

	_spm.setOnTheFly(1);
	_spm.calculateEnergies();
	_spm.runGreedyOptimizer(_greedyCycles);
}

void deleteTerminalHydrogenBondInteractions(System &_sys, Options &_opt){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	int frontExt = _opt.tmStart;
	int endExt = _opt.endResNum;
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

Options parseOptions(int _argc, char * _argv[], Options defaults);

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Datafile: " << _op.datafile << endl << endl;

	_fout << "Options from Datafile" << endl;
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "Warning Messages: " << _op.warningMessages << endl << endl;

	_fout << "Other Parameters" << endl;
	_fout << "backboneCrd " << _op.backboneCrd << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	_fout << "tmStart " << _op.tmStart << endl;
	_fout << "tmEnd " << _op.tmEnd << endl;

	_fout << "helixGeoFile " << _op.helixGeoFile << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "hbondFile " << _op.hbondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "deleteTerminalHbonds " << _op.deleteTerminalHbonds << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.identitySize() > 1){
			for (uint j=0; j < pos.getNumberOfIdentities(); j++){
				pos.setActiveIdentity(j);
				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
					}
				}
				pos.setActiveIdentity(0);
			}
		}
		else{
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}
	}	
}

vector<string> weightAAs(vector<string> &_AAs, vector<int> &_weights, int _seqLength, RandomNumberGenerator &_RNG){
	vector<string> weightedAAs;
	for (uint j=0; j<_AAs.size(); j++){
		for (uint k=0; k<_weights[j]; k++){
			weightedAAs.push_back(_AAs[j]);
		}
	}
	random_shuffle(weightedAAs.begin(), weightedAAs.end(), _RNG);
	return weightedAAs;
}

string randomAASeqDiffTermini(vector<string> _weightedAAs, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	for (uint i=0; i<_seqLength; i++){
		if (i < 4 || i > _seqLength-5){
			if (i ==_seqLength-3 || i ==_seqLength-1){
				seq += "I";
			}
			else{
			    seq += "L";
			}
		}
		else{
			int a=0;
			a = _RNG.getRandomInt(_weightedAAs.size()-1);
			//cout << i << ": " << a << endl;
			if (a == -1){
				while (a == -1){
					a = _RNG.getRandomInt(_weightedAAs.size()-1);
				}
			}//I don't think I need this as long as I change a to an int
			//cout << "AA: " << weightedAAs[a] << endl;
			//randAA.push_back(a);
			seq += _weightedAAs[a];
		}
	}
	return seq;
}

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, int _thread){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=_thread; i<_seqLength+_thread; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi"); 
		ener.push_back(resi);
	}
	sel.clearStoredSelections();
	return ener;
}
				
vector<double> calcPairBaselineEnergies(System &_sys, string _seq, Options _opt, int _seqLength, int _thread){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=_thread; i<_seqLength+_thread; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i);
		sel.select(residue += num1);
		char s1 = _seq[i-_thread];
		for (uint j=i+1; j<_seqLength+_thread; j++){
			int dist = j-i;
			if (dist <= _opt.pairDist){
				char s2 = _seq[j-_thread];
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				//cout << dist << " " << s1 << ":" << s2 << " = " << pair << endl;
				ener.push_back(pair);
			}
			else{
				j = _seqLength+_thread;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

double sumEnergyVector(vector<double> _energies){
	double ener = 0;
	for (uint i=0; i<_energies.size(); i++){
		ener = ener + _energies[i];
	}
	return ener;
}

double computeMonomerEnergy(Options& _opt, string &_seq, double &_monomerHbond, double &_monomerVdw, double &_monomerSelf, double &_monomerPair, RandomNumberGenerator &_RNG) {
	string polySeq = convertToPolymerSequenceNeutralPatch(_seq, _opt.thread);
	PolymerSequence PS(polySeq);
	cout << PS << endl;

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	//CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile);
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
	
	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
	}

	//CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(_opt.backboneCrd); 
	if(!cRead.read()) {
		cerr << "Unable to read " << _opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	monoSys.assignCoordinates(glyAPV,false);
	monoSys.buildAllAtoms();

	CSBMono.updateNonBonded(10,12,50);

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

	/******************************************************************************
	 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
	 ******************************************************************************/
	PDBWriter writer;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	writer.open(dir + "/monomer.pdb");
	for (uint i=0; i<monoSys.atomSize(); i++){
		Atom at = monoSys.getAtom(i);
		if (!at.hasCoor()){
			cout << "Atom " << i << " was not assigned coordinates; program termination";
			break;
		}
		else{
			continue;
		}
	}
	cout << "All atoms have coordinates" << endl;


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
	monoEset->setTermActive("CHARMM_IMM1REF", true);
	monoEset->setTermActive("CHARMM_IMM1", true);

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);
	monoEset->setWeight("CHARMM_IMM1REF", 1);
	monoEset->setWeight("CHARMM_IMM1", 1);
//	if (_useIMM1 == false){
//		monoEset->setTermActive("CHARMM_IMM1REF", false);
//		monoEset->setTermActive("CHARMM_IMM1", false);
//		monoEset->setWeight("CHARMM_IMM1REF", 0);
//		monoEset->setWeight("CHARMM_IMM1", 0);
//	}
//	else{
//		monoEset->setTermActive("CHARMM_IMM1REF", true);
//		monoEset->setTermActive("CHARMM_IMM1", true);
//		monoEset->setWeight("CHARMM_IMM1REF", 1);
//		monoEset->setWeight("CHARMM_IMM1", 1);
//	}

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);
	
	/*****************************************************************************
	 *                 === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, _opt.SL);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(true);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.calculateEnergies();

	repackSideChains(monoSpm, _opt.greedyCycles);

	/*****************************************************************************
	 *            === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ===
	 ******************************************************************************/
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	double monoEner = monoSys.calcEnergy()*2;
	double monomerHbond = monoEset->getTermEnergy("SCWRL4_HBOND")*2; // double the energy for 2 helices
	double monomerVdw = monoEset->getTermEnergy("CHARMM_VDW")*2; // double the energy for 2 helices
	writer.write(monoSys.getAtomPointers(), true, false, true);
	writer.close();
	cout << "Monomer Calc: " << monoEner << endl;
	cout << "Monomer Hbond Calc: " << monomerHbond << endl;
	cout << "Monomer Vdw Calc: " << monomerVdw << endl;
	cout << monoSys.getEnergySummary() << endl;

	_monomerHbond = monomerHbond;
	_monomerVdw = monomerVdw;
	vector<double> self = calcBaselineEnergies(monoSys, 21, 25);
	vector<double> pair = calcPairBaselineEnergies(monoSys, _seq, _opt, 21, 25);
	_monomerSelf = sumEnergyVector(self)*2; // double the energy for 2 helices
	_monomerPair = sumEnergyVector(pair)*2; // double the energy for 2 helices
	double finalEnergy = monoEner; // double the energy for 2 helices
	return finalEnergy;
}

map<string, double> readSingleParameters(string _baselineFile){
	Reader selfReader(_baselineFile);
	selfReader.open();
	map<string, double> selfEnergies;

	if(!(selfReader.is_open())){
		cerr << "WARNING: Unable to open " << _baselineFile << endl;
		exit(0);
	}

	vector<string> lines = selfReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], "\t");
		if (tokens.size() <1){
			continue;
		}
		if (tokens.size() != 2){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) Energy(double)";
			continue;
		}
		selfEnergies[MslTools::toUpper(tokens[0])] = MslTools::toDouble(tokens[1]);
	}

	selfReader.close();
	return selfEnergies;

}

map<string,map<string,map<uint, double>>> readPairParameters(string _baselineFile){
	Reader pairReader(_baselineFile);
	pairReader.open();
	map<string,map<string,map<uint, double>>> pairEnergies;

	if(!(pairReader.is_open())){
		cerr << "WARNING: Unable to open " << _baselineFile << endl;
		exit(0);
	}

	vector<string> lines = pairReader.getAllLines();

	for (int i=0; i<lines.size(); i++){
		vector<string> tokens = MslTools::tokenize(lines[i], " ");
		if (tokens.size() < 1){
			continue;
		}
		if (tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) ResName(string) Distance(uint) Energy(double)";
			continue;
		}
		if (tokens[0].compare(tokens[1]) == 0){
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		}
		else{
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
			pairEnergies[MslTools::toUpper(tokens[1])][MslTools::toUpper(tokens[0])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		}

	}

	pairReader.close();
	return pairEnergies;

}

void buildBaselineIMM1Interactions(System &_sys, map<string, double> &_selfMap, vector<double> &_selfBaselines){
	EnergySet* Eset = _sys.getEnergySet();

	for (uint i=0; i<1; i++){
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for(uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);
				string baseId = res.getResidueName();
				if (p-positions.begin() < 4 || p-positions.begin() > positions.size()-5){
					baseId = baseId.append("-OUT");
				}
				try{
					double ener = _selfMap.at(baseId);
					Atom *a = &res.getAtom("CA");
					_selfBaselines.push_back(ener);
					Eset->addInteraction(new BaselineIMM1Interaction(*a, ener));
				}
				catch (const out_of_range& e){
					continue;
				}
			}
		}
	}
}

void buildSelfInteractions(System &_sys, map<string, double> &_selfMap, vector<double> &_selfBaselines){
	EnergySet* Eset = _sys.getEnergySet();

	for (uint i=0; i<1; i++){
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for(uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);
				string baseId = res.getResidueName();
				if (p-positions.begin() < 4 || p-positions.begin() > positions.size()-5){
					baseId = baseId.append("-OUT");
				}
				try{
					double ener = _selfMap.at(baseId);
					Atom *a = &res.getAtom("CA");
					_selfBaselines.push_back(ener);
					Eset->addInteraction(new BaselineInteraction(*a, ener));
				}
				catch (const out_of_range& e){
					continue;
				}
			}
		}
	}
}

void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>> &_pairMap, vector<double> &_pairBaselines){
	EnergySet* Eset = _sys.getEnergySet();
	ofstream pout;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	string poutName  = dir + "/pairBaseline_2.out";	
	pout.open(poutName.c_str());

	for (uint i=0; i<1; i++){
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();

		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for(uint j=0; j<(*p)->identitySize();j++){
				Residue &res1 = (*p)->getIdentity(j);
				string baseId1 = res1.getResidueName();
				if (p-positions.begin() < 1){
					baseId1 = baseId1.append("-ACE");
				}
				//if (p-positions.begin() > positions.size()-2){
				//	baseId1 = baseId1.append("-CT2");
				//}
				for(vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					if (d <= 10){
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							//if (p2-positions.begin() < 1){
							//	baseId2 = baseId2.append("-ACE");
							//}
							if (p2-positions.begin() > positions.size()-2){
								baseId2 = baseId2.append("-CT2");
							}
							try{
								map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
								map<uint,double> AA2 = AA1.at(baseId2);
								double ener = AA2.at(d);
								//cout << d << " " << baseId1 << ": " << baseId2 << " = " << ener << endl;
								Atom *a = &res1.getAtom("CA");
								Atom *b = &res2.getAtom("CA");
								_pairBaselines.push_back(ener);
								Eset->addInteraction(new BaselinePairInteraction(*a,*b,-1*ener));
							}
							catch (const out_of_range& e){
								continue;
							}
						}
					}
				}
			}
		}
	}
	pout.close();
}

void printEnerFile(string _seq, Options _opt, vector<double> &_ener, vector<double> &_bEner, int _seqNumber, bool _selfEner, ofstream &_out){
	if (_selfEner ==true){
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "AA:     Position:     Energy:     Baseline" << endl;
		}
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << ":    " << i+1 << ":        " << _ener[i] << ":      " << _bEner[i] << endl;
		}
	}
	else{
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "Pair: Distance: Position1: Position2: Energy: Baseline" << endl;//could be interesting to add rotamer number to this
		}
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				int d = j-i;
				if (d <= _opt.pairDist){
					stringstream ss;
					ss <<_seq[i] << _seq[j];
					string p = ss.str();
					pairs.push_back(p);
					dists.push_back(d);
				}
			}
		}
		string sp = ": ";
		int pos1 = 0;
		for (uint i=0; i<pairs.size(); i++){
			if (dists[i] == 1){
				pos1++;
			}
			int pos2 = pos1 + dists[i];
			_out << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << _bEner[i] << endl;
		}
	}
}

/******************************************
 *  
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%Y_%m_%d",timeinfo);
	string date(buffer);
	
	cout << date << endl;

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
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

	/**********************************************************************************
	*
	*    printOutFiles
	*
	**********************************************************************************/
	ofstream bsout;
	ofstream sout;
	ofstream pout;
	
	opt.pdbOutputDir = opt.pdbOutputDir + "/" + date;
	//string dir = opt.pdbOutputDir;
	string dir = "/exports/home/gloiseau/mslib/trunk_AS";
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string bsoutName  = dir + "/baselineEnergySummary_" + opt.runNumber + ".out";
	string soutName   = dir + "/selfEnergyComparison_" + opt.runNumber + ".out";
	string poutName   = dir + "/pairEnergyComparison_" + opt.runNumber + ".out";
	
	bsout.open(bsoutName.c_str());
	sout.open(soutName.c_str());
	pout.open(poutName.c_str());
	
	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd); 
	if(!cRead.read()) {
		cout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	/******************************************************************************
	 *                      === INITIALIZE PDBWRITERS ===
	 ******************************************************************************/
	PDBWriter writer;
	writer.open(dir + "/polyLeu.pdb");
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	//vector<int> highEnergySequences;
	int seqDiscard = 0;
	int seqAccept = 0;
	int seqNumber = opt.seqNumber;
	
	//Random Number Generator
	RandomNumberGenerator RNG;
	//RNG.setTimeBasedSeed();
	RNG.setSeed(opt.seed);
	
	/******************************************************************************
	 *                     === GIVE SEQUENCES WEIGHTS ===
	 ******************************************************************************/
	vector<string> weightedAAs = weightAAs(opt.ids, opt.weights, opt.backboneLength, RNG);
	
	if (opt.backboneLength > weightedAAs.size()){
		cout << "ERROR: Length of protein is too long for total number of AAs given." << endl;
		exit(1);
	}
	
	/******************************************************************************
	 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
	 ******************************************************************************/
	while (seqAccept < seqNumber){
		cout << "Number Accepted: " << seqAccept << endl;
		cout << "Number Discarded: " << seqDiscard << endl;
		
		/******************************************************************************
		 *                     === RANDOMIZE SEQUENCES ===
		 ******************************************************************************/
		//string seq = randomAASeqDiffTermini(weightedAAs, opt.backboneLength, RNG);
		string seq = "LLLLLIAFLLGILLWALLILI";
		cout << seq << endl;
		string polySeq = convertToPolymerSequenceNeutralPatch(seq, opt.thread);
		PolymerSequence PS(polySeq);
		cout << PS << endl;
		
		/******************************************************************************
		 *                     === DECLARE SYSTEM ===
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
		CSB.setIMM1Params(15,10);

		CSB.setBuildNonBondedInteractions(false);

		if(!CSB.buildSystem(PS)) {
			cerr << "Unable to build system from " << polySeq << endl;
			exit(0);
		} else {
			//fout << "CharmmSystem built for sequence" << endl;
		}
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();
		
		CSB.updateNonBonded(10,12,50);	
		
		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
		sysRot.defineRotamerSamplingLevels();
		
		// Add hydrogen bond term
		HydrogenBondBuilder hb(sys, opt.hbondFile);
		hb.buildInteractions(50);
	
		/******************************************************************************
		 *                === CHECK TO SEE IF ALL ATOMS ARE BUILT ===
		 ******************************************************************************/
		for (uint i=0; i<sys.atomSize(); i++){
			Atom at = sys.getAtom(i);
			if (!at.hasCoor()){
				cout << "Atom " << i << " was not assigned coordinates; program termination";
				break;
			}
			else{
				continue;
			}
		}

		/******************************************************************************
		 *                     === INITIAL VARIABLE SET UP ===
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
		Eset->setTermActive("CHARMM_IMM1REF", true);
		Eset->setTermActive("CHARMM_IMM1", true);
		//Eset->setTermActive("CHARMM_IMM1REF", false);
		//Eset->setTermActive("CHARMM_IMM1", false);
		
		// Set weights
		Eset->setWeight("CHARMM_VDW", 1);
		Eset->setWeight("SCWRL4_HBOND", 1);
		Eset->setWeight("CHARMM_IMM1REF", 1);
		Eset->setWeight("CHARMM_IMM1", 1);
		
		/******************************************************************************
		 *           === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
		 ******************************************************************************/
		deleteTerminalHydrogenBondInteractions(sys,opt);

		/******************************************************************************
		 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
		 ******************************************************************************/
		loadRotamers(sys, sysRot, opt.SL);
		
		SelfPairManager spm;
		spm.seed(RNG.getSeed());
		spm.setSystem(&sys);
		spm.setVerbose(true);
		spm.getMinStates()[0];
		spm.updateWeights();
		spm.setOnTheFly(true);
		spm.calculateEnergies();
		
		repackSideChains(spm, 10);
		
		/******************************************************************************
		 *                  === SET SYSTEM TO BEST SPM ROTAMERS ===
		 ******************************************************************************/
		sys.setActiveRotamers(spm.getMinStates()[0]);
		sys.calcEnergy();
		double totalEnergy = sys.calcEnergy();
	
		/******************************************************************************
		 *             === EXTRACT SELF PAIR ENERGIES FOR PRINTING ===
		 ******************************************************************************/
		vector<double> selfEnerVec = calcBaselineEnergies(sys, seq.length(), opt.thread); //reapplying this function from generateSelfPairBaselinevdW_01.cpp to extract self and pair energies for e`ach sequence	
		vector<double> pairEnerVec = calcPairBaselineEnergies(sys, seq, opt, seq.length(), opt.thread); //reapplying this function from generateSelfPairBaselinevdW_01.cpp to extract self and pair energies for e`ach sequence	
		/******************************************************************************
		 *               === CALCULATE ENERGIES FOR EACH POSITION ===
		 ******************************************************************************/
		//set up way to select for the 22 inner residues (selection seems a little weird to me, so I think I'm going to try the alternate way subtracting from the outermeans first)
		totalEnergy = sys.calcEnergy();
		cout << "Total Energy: " << totalEnergy << endl;
		cout << sys.getEnergySummary() << endl;
		writer.write(sys.getAtomPointers(), true, false, true);

		double selfEner = sumEnergyVector(selfEnerVec);
		double pairEner = sumEnergyVector(pairEnerVec);

		double monomerHbond = 0;	
		double monomerVdw = 0;	
		double monomerSelf = 0;	
		double monomerPair = 0;	
		//double monomer = computeMonomerEnergy(opt, seq, monomerHbond, monomerVdw, monomerSelf, monomerPair, RNG);

		//map<string, double> selfMap = readSingleParameters("/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanSelf_par.txt");
		//map<string,map<string,map<uint,double>>> pairMap = readPairParameters("/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanPair_par.txt");
		map<string, double> selfMap = readSingleParameters(opt.selfEnergyFile);
		//map<string, double> imm1Map = readSingleParameters(opt.imm1EnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(opt.pairEnergyFile);
	
		vector<double> selfBaselines;
		//vector<double> imm1Baselines;
		vector<double> pairBaselines;
		buildSelfInteractions(sys, selfMap, selfBaselines);
		//buildBaselineIMM1Interactions(sys, imm1Map, imm1Baselines);
		buildPairInteractions(sys, pairMap, pairBaselines);

		sys.calcEnergy();

		double seqSelfEner = Eset->getTermEnergy("BASELINE");
		//double seqBaselineIMM1Ener = Eset->getTermEnergy("BASELINE_IMM1");
		double seqPairEner = Eset->getTermEnergy("BASELINE_PAIR");
		
		double seqBaselineEnergy = seqSelfEner + seqPairEner;
		double selfPairDif = totalEnergy+seqBaselineEnergy;
		cout << "Total Energy: " << totalEnergy << endl;
		//cout << "Self Energy: " << selfEner << endl;
		//cout << "Pair Energy: " << pairEner << endl;
		cout << "Baseline Self Energy: " << seqSelfEner << endl;
		//cout << "Baseline IMM1 Energy: " << seqBaselineIMM1Ener << endl;
		cout << "Baseline Pair Energy: " << seqPairEner << endl;
		//cout << "SelfPair Energy: " << seqBaselineEnergy << endl;
		
		cout << "Total - Baseline Energy: " << selfPairDif << endl;

		//cout << "Monomer: " << monomer << endl;
		cout << "Monomer Self: " << monomerSelf << "; Baseline Self: " << seqSelfEner << ", " << selfEner <<  endl;
		cout << "Monomer Pair: " << monomerPair << "; Baseline Pair: " << seqPairEner << ", " << pairEner <<  endl;
		cout << "Monomer Hbond: " << monomerHbond << endl;
		cout << "Monomer VDW: " << monomerVdw << endl;
		
		/******************************************************************************
		 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
		 ******************************************************************************/
		if (seqAccept == 0 && opt.seed == 0){
			bsout << "Monomer: SelfPair: Difference:" << endl;
		}
		bsout << totalEnergy << ":    " << seqBaselineEnergy << ":   " << selfPairDif << endl;
		if (pairEnerVec.size() > pairBaselines.size()){
			while (pairEnerVec.size() != pairBaselines.size()){
				pairBaselines.push_back(0);//for some reason my code isn't getting consistent numbers of baselines; this is a cheat until I fix the functions (doesn't work; check YA on line 1570 to be sure if my function fixes work or not
			}
		}
		else{
			if (pairEnerVec.size() < pairBaselines.size()){
				cout << "Pair size is not equal to baseline size: Pair " << pairEnerVec.size() << "vs Baseline " << pairBaselines.size() << "." << endl;
				break;
			}
		}
		printEnerFile(seq, opt, selfEnerVec, selfBaselines, seqAccept, true, sout);
		printEnerFile(seq, opt, pairEnerVec, pairBaselines, seqAccept, false, pout);
		seqAccept++;
		sys.reset();
	}
	time (&endTime);
	diffTime = difftime(endTime, startTime);
	cout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	writer.close();
	bsout.close();
	sout.close();
	pout.close();
}

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
	 *  /exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/helixGenerator.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//opt.required.push_back("");
	//opt.allowed.push_back("");

	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	opt.allowed.push_back("seqNumber");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("thread");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("SL");
	
	// energy weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	// input monomerEnergy
	opt.allowed.push_back("inputMonomerE");
	opt.allowed.push_back("monoE_vdw");
	opt.allowed.push_back("monoE_hbond");
	opt.allowed.push_back("monoE_solv");
	opt.allowed.push_back("monoE_solvRef");
	opt.allowed.push_back("printTermEnergies");

	// input files
	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("selfEnergyFile");
	opt.allowed.push_back("imm1EnergyFile");
	opt.allowed.push_back("pairEnergyFile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("solvFile");
	
	// alternate ids and weights
	opt.allowed.push_back("ids");
	opt.allowed.push_back("weights");
	opt.allowed.push_back("baselines");
	opt.allowed.push_back("baselinesInner");
	opt.allowed.push_back("baselinesvdw");
	opt.allowed.push_back("baselinesvdwInner");
	
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("pairDist");

	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		cerr << "No options given!" << endl;
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
	
	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()){
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}
	opt.pairDist = OP.getInt("pairDist");
	if (OP.fail()){
		opt.warningMessages += "pairDist not specified, using 15\n";
		opt.warningFlag = true;
		opt.pairDist = 15;
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
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}

	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}

	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.tmStart+opt.backboneLength;
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
	}

	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(opt.tmStart) + "\n";
		opt.warningFlag = true;
		opt.startResNum = opt.tmStart;
	}

	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.tmEnd;
	}

	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified using 25\n";
		opt.warningFlag = true;
		opt.thread = 25;
	}

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}

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
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "70.00";
	}
	opt.SL = "SL"+opt.SL;

	opt.backboneAA = OP.getString("backboneAA");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneAA not specified, default to glycine\n";
		opt.backboneAA = "G";
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 35\n";
		opt.backboneLength = 35;
	}
	opt.seqNumber = OP.getInt("seqNumber");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "seqNumber not specified, default to 10\n";
		opt.seqNumber = 10;
	}

	opt.helixGeoFile = OP.getString("helixGeoFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "helixGeoFile not specified, default to /data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
		opt.helixGeoFile = "/data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
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
	opt.selfEnergyFile = OP.getString("selfEnergyFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "selfEnergyFile not specified, default to /exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_selfBaseline.txt";
		opt.selfEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_selfBaseline.txt";
	}
	opt.imm1EnergyFile = OP.getString("imm1EnergyFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "imm1EnergyFile not specified, default to /exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_selfBaseline.txt";
		opt.imm1EnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_selfBaseline.txt";
	}
	opt.pairEnergyFile = OP.getString("pairEnergyFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "pairEnergyFile not specified, default to /exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_pairBaseline.txt";
		opt.pairEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_pairBaseline.txt";
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
	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine backboneCrd";
		opt.errorFlag = true;
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

	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) { 
		opt.warningMessages += "solvFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.solvFile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	opt.infile = OP.getString("infile");
	if (OP.fail()) { 
		opt.warningMessages += "infile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.infile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	
	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine pdbOutputDir";
		opt.errorFlag = true;
	}

	opt.ids = OP.getStringVector("ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}
	opt.weights = OP.getIntVector("weights");
	if (OP.fail()) {
		opt.warningMessages += "Unable to identify alternate AA weights, defaulting each weight to 10\n";
		opt.warningFlag = true;
		for (uint i=0; i<opt.ids.size(); i++){
			opt.weights.push_back(10);
		}
	}
	if (opt.weights.size() != opt.ids.size()){
		opt.errorMessages += "Unable to identify alternate AA weights, make sure to correspond a weight to each AA\n";
		opt.errorFlag = true;
	}
	opt.rerunConf = OP.getConfFile();

	return opt;
}
