#include <iostream>
#include <fstream>
#include <sstream>

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

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "generateBaselinevdW_IMM1";
string programDescription = "This program designs randomized sequences to calculate baselines for future sequence design of vdW dependent sequences; update to include IMM1";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "25 November 2019";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

/******************************************************************************************************************************************************************************/

struct Options{
	string sequence;
	string backboneAA;
	int backboneLength;
	int seqNumber;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int sequenceStart;

	string SL; //number of rotamers

	// input files
	string helixGeoFile;
	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string monoRotLibFile;
	string infile;
	string rulesFile;

	// side-chain repack variable
	int mcCycles;
	int mcMaxRejects;
	double mcStartTemp;
	double mcEndTemp;
	int mcCurve;

	bool verbose;
	int greedyCycles;
	int seed;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	int start;
	int end;

	// alternate identities and weights
	vector<string> ids;
	vector<double> weights;

	// housekeeping variables
	string runNumber;
	bool rejectvdw;

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

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

string convertToPolymerSequence(string _seq, int _startResNum) {
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
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}//changed this for only one chain

string generatePolyLeu(string _backboneAA, int _sequenceLength) {
	string polyLeu = "";
	for (uint i=0; i<_sequenceLength; i++){
		polyLeu = polyLeu + _backboneAA;
	}
	return polyLeu;
}

string generateMultiIDPolymerSequence(string _seq, vector<string> _alternateIds) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		ps = ps + " [";
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
		for (uint i=0; i<_alternateIds.size(); i++){
			if(_alternateIds[i] == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + _alternateIds[i];
			}
		}
		ps = ps + "] ";
		counter++;
	}
	ps = ":{" + MslTools::intToString(0) + "} " + ps;
	//return "A" + ps + "\nB" + ps;
	return "A" + ps;
}

string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _varPos) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if (_varPos[counter] == 1){
			ps = ps + " [";
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
			for (uint i=0; i<_alternateIds.size(); i++){
				if(_alternateIds[i] == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + _alternateIds[i];
				}
			}
			ps = ps + "] ";
		} else {
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		}
		counter++;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

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

//void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable) {
void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {

	_spm.setOnTheFly(1);
	//_spm.recalculateNonSavedEnergies(_savedEnergyFlagTable);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
}
//TODO: figure out how below works and use below instead of my outputs?
map<string,double> getEnergyByTerm(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

map<string,double> getEnergyByTermDoubled(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  2.0* _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

Options parseOptions(int _argc, char * _argv[], Options defaults);

//skeleton code for adding given sequence into different portions of a polyLeu chain
string generateSequence(int _startResNum, int _sequenceLength, int _backboneLength, string _backboneAA, string &_sequence, bool _backbone=true){
	string seq = "";
	for(uint i=0; i<_sequenceLength; i++){
		//somehow multiply the amount of _backboneAA, then add in the _sequence, then add in leftover _backboneAA
		//might be able to switch this to a while loop like below?
		seq += _backboneAA;	
	}
	//seq += _sequence; use at some point if I want to do something like threading through the sequence
	if (_backbone == true){
		while (seq.length() < _backboneLength){
			seq += _backboneAA;
			//cout << seq << endl;
		}
	}
	return seq;
}

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "Warning Messages: " << _op.warningMessages << endl << endl;

	_fout << "Other Parameters" << endl;
	_fout << "backboneCrd " << _op.backboneCrd << endl;
	//_fout << "logFile " << _op.logFile << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	//_fout << "fullSequence " << _op.fullSequence << endl;
	_fout << "tmStart " << _op.tmStart << endl;
	_fout << "tmEnd " << _op.tmEnd << endl;

	_fout << "helixGeoFile " << _op.helixGeoFile << endl;
	_fout << "rulesFile " << _op.rulesFile << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "solvFile " << _op.solvFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;
	_fout << "monoRotLibFile " << _op.monoRotLibFile << endl;

	_fout << "MCCycles " << _op.mcCycles << endl;
	_fout << "MCMaxRejects " << _op.mcMaxRejects << endl;
	_fout << "MCStartTemp " << _op.mcStartTemp << endl;
	_fout << "MCEndTemp " << _op.mcEndTemp << endl;
	_fout << "MCCurve " << _op.mcCurve << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "fullSequenceStart " << _op.sequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;
	_fout << "weight_solv " << _op.weight_solv << endl;

	_fout << "SL " << _op.SL << endl;

	_fout << "ids ";
	for (uint i=0; i<_op.ids.size(); i++){
		_fout << _op.ids[i];
	}
	_fout << endl;

	_fout << "weights ";
	for (uint i=0; i<_op.weights.size(); i++){
		_fout << _op.weights[i];
	}
	_fout << endl;

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

vector<string> weightAAs(vector<string> &_AAs, vector<double> &_weights, int _seqLength, RandomNumberGenerator &_RNG){
	vector<string> weightedAAs;
	for (uint j=0; j<_AAs.size(); j++){
		for (uint k=0; k<_weights[j]; k++){
			weightedAAs.push_back(_AAs[j]);
		}
	}
	random_shuffle(weightedAAs.begin(), weightedAAs.end(), _RNG);
	return weightedAAs;
}

//TODO: update this to get randomAASequence without using integers; I think it works now?
string randomAASequence(vector<string> _ids, vector<double> _weights, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	double sumWeights=0;
	for (uint i=0; i<_weights.size(); i++){
		sumWeights += _weights[i];
	}
	for (uint i=0; i<_seqLength; i++){
		double a=0;
		a = _RNG.getRandomDouble(sumWeights);
		for (uint j=0; j<_weights.size(); j++){
			if (a<_weights[j]){
				seq += _ids[j];
				j=_weights.size();
			}
			a -= _weights[j];
		}
	}
	return seq;
}

string randomAASequence(vector<string> _weightedAAs, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	for (uint i=0; i<_seqLength; i++){
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
	return seq;
}

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, double &_totalProt){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	sel.select("allProt, all");
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+28);//TODO: make this an option based on threading
		sel.select(residue += number);
		double prot = _sys.calcEnergy("resi", "allProt")/2;
		double resi = _sys.calcEnergy("resi")/2; 
		double protEnergy = prot+resi;
		ener.push_back(protEnergy);
		_totalProt += protEnergy;
	}
	return ener;
}
				
void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, double _totalEnergy, double _hbond, double _vdw, double _imm1, int _seqNumber, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	_out << "AA      Position      Energy" << endl;//could be interesting to add rotamer number to this
	for (uint i=0; i<_seq.length(); i++){
		_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
	}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	_out << "Total Energy: " << _totalEnergy << endl;
	_out << "H-bond Energy: " << _hbond << endl;
	_out << "VDW Energy: " << _vdw << endl << endl;
	_out << "IMM1 Energy: " << _imm1 << endl << endl;
}

void printEnerFile(string _seq, vector<double> &_ener, int _seqNumber, ofstream &_out){
	if (_seqNumber == 0){
		_out << "AA      Position      Energy" << endl;//could be interesting to add rotamer number to this
	}
	for (uint i=0; i<_seq.length(); i++){
		_out << _seq[i] << ":     " << i+1 << ":            " << _ener[i] << endl;
	}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
}

//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

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
	
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);
	
	cout << date << endl;

	time(&startTime);
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
	ofstream pout;
	ofstream rout;
	ofstream nout;
	//opt.pdbOutputDir = opt.pdbOutputDir + "/" + date + "/" + opt.datafile;
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string poutName = dir + "/Options_" + opt.runNumber + ".out";
	string routName = dir + "/Residue_Energies_" + opt.runNumber + ".out";
	string noutName = dir + "/Sequences_" + opt.runNumber + ".out";
	
	pout.open(poutName.c_str());
	rout.open(routName.c_str());
	nout.open(noutName.c_str());
	
	/******************************************************************************
	 *                         === INITIALIZE POLYGLY ===
	 ******************************************************************************/
	//string polyGly69 = generateSequence(0, 69, opt.backboneLength, G, opt.sequence, true);
	//cout << polyGly69 << endl;
	
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
	PDBWriter writer2;
	//PDBWriter writer6;
	
	writer.open(dir + "/allSeq_" + opt.runNumber + ".pdb");
	writer2.open(dir + "/acceptedSeq_" + opt.runNumber + ".pdb");
	//writer6.open(dir + "/BestUnbiasedMCState.pdb");
	//TODO: starting here, start a for loop to generate and write energies of multiple sequences
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	//vector<int> highEnergySequences;
	int seqDiscard = 0;
	int seqAccept = 0;
	uint a = 0;
	int seqNumber = opt.seqNumber;

	double bestEnergy = 0;
	double bestvdw = 0;
	double bestimm1 = 0;
	string bestSeq = "";
	string bestvdwSeq = "";
	string bestimm1Seq = "";
	double worstEnergy = -10000000;
	double worstvdw = -10000000;
	double worstimm1 = -10000000;
	string worstSeq = "";
	string worstvdwSeq = "";
	string worstimm1Seq = "";
	//vector<string> completeSequences;
	vector<string> str;
	vector<double> doub;
	
	//Random Number Generator
	RandomNumberGenerator RNG;
	RNG.setTimeBasedSeed();

	//TODO: add in a way to check and make sure I'm not getting the same sequence?
	/******************************************************************************
	 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
	 ******************************************************************************/
	while (seqAccept < seqNumber){
		cout << "Number Accepted: " << seqAccept << endl;
		cout << "Number Discarded: " << seqDiscard << endl;
		a = seqAccept;
		
		/******************************************************************************
		 *                     === RANDOMIZE SEQUENCES ===
		 ******************************************************************************/
		string seq = randomAASequence(opt.ids, opt.weights, opt.backboneLength, RNG);
		cout << seq << endl;
		string polySeq = convertToPolymerSequence(seq, 28);//to keep thread consistent with the design code
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
		
		CSB.setSolvent("MEMBRANE");
		CSB.setIMM1Params(15, 10);
		
		if(!CSB.buildSystem(PS)) {
			cerr << "Unable to build system from " << polySeq << endl;
			exit(0);
		} else {
			//fout << "CharmmSystem built for sequence" << endl;
		}
		
		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
		sysRot.defineRotamerSamplingLevels();
		
		// Add hydrogen bond term
		HydrogenBondBuilder hb(sys, opt.hBondFile);
		hb.buildInteractions(30);
		
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
		Eset->setTermActive("SCWRL4_HBOND", false);
		
		// Set weights
		Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
		//Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
		//Eset->setWeight("CHARMM_IMM1REF", 1);
		//Eset->setWeight("CHARMM_IMM1", 1);
		Eset->setWeight("CHARMM_IMM1REF", 1);
		Eset->setWeight("CHARMM_IMM1", 1);
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);
		sys.buildAllAtoms();
		
		writer.write(sys.getAtomPointers(), true, false, true);
		
		/******************************************************************************
		 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
		 ******************************************************************************/
		CSB.updateNonBonded();
		sys.buildAllAtoms();
		
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
		double vdw = Eset->getTermEnergy("CHARMM_VDW");
		double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
		double imm1 = Eset->getTermEnergy("CHARMM_IMM1REF") + Eset->getTermEnergy("CHARMM_IMM1");
		
		if (opt.rejectvdw){
			if (totalEnergy < 0 && vdw < 0){
				/******************************************************************************
				 *               === CALCULATE ENERGIES FOR EACH POSITION ===
				 ******************************************************************************/
				double totalProt = 0;
				vector<double> AAener = calcBaselineEnergies(sys, seq.length(), totalProt);
				writer2.write(sys.getAtomPointers(), true, false, true);
			
				totalEnergy = sys.calcEnergy();
				cout << "Total Energy: " << totalEnergy << endl;
				cout << sys.getEnergySummary() << endl;
				
				int finalEnergy = totalEnergy;
				int finalTotal = totalProt;
			
				//cout << "Total - ener: " << finalEnergy-finalTotal << endl;//to prove that the energies are the same since it was giving me some very small number that wasn't 0
				
				if (finalEnergy<bestEnergy){
					bestEnergy = finalEnergy;
					bestSeq = seq;
				}
				if (finalEnergy>worstEnergy){
					worstEnergy = finalEnergy;
					worstSeq = seq;
				}
				if (vdw<bestvdw){
					bestvdw = vdw;
					bestvdwSeq = seq;
				}
				if (imm1<bestimm1){
					bestimm1 = imm1;
					bestimm1Seq = seq;
				}
				if (vdw>worstvdw){
					worstvdw = vdw;
					worstvdwSeq = seq;
				}
				if (imm1>worstimm1){
					worstimm1 = imm1;
					worstimm1Seq = seq;
				}
				/******************************************************************************
				 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
				 ******************************************************************************/
				printSeqFile(PS, seq, AAener, totalEnergy, hbond, vdw, imm1, a, nout);
				printEnerFile(seq, AAener, a, rout);
				seqAccept++;
				sys.reset();
			}
			else{
				if (totalEnergy > 0){
					cout << "Total Energy " << totalEnergy << " is too high; go to next random sequence." << endl;
				}
				if (vdw > 0){
					cout << "VDW Energy " << vdw << " is too high; go to next random sequence." << endl;
				}
				seqDiscard++;
				sys.reset();
			}
		}
		else{
			if (totalEnergy < 0){
				/******************************************************************************
				 *               === CALCULATE ENERGIES FOR EACH POSITION ===
				 ******************************************************************************/
				double totalProt = 0;
				vector<double> AAener = calcBaselineEnergies(sys, seq.length(), totalProt);
				writer2.write(sys.getAtomPointers(), true, false, true);
			
				totalEnergy = sys.calcEnergy();
				cout << "Total Energy: " << totalEnergy << endl;
				cout << sys.getEnergySummary() << endl;
				
				int finalEnergy = totalEnergy;
				int finalTotal = totalProt;
			
				//cout << "Total - ener: " << finalEnergy-finalTotal << endl;//to prove that the energies are the same since it was giving me some very small number that wasn't 0
				
				if (finalEnergy<bestEnergy){
					bestEnergy = finalEnergy;
					bestSeq = seq;
				}
				if (finalEnergy>worstEnergy){
					worstEnergy = finalEnergy;
					worstSeq = seq;
				}
				if (vdw<bestvdw){
					bestvdw = vdw;
					bestvdwSeq = seq;
				}
				if (imm1<bestimm1){
					bestimm1 = imm1;
					bestimm1Seq = seq;
				}
				if (vdw>worstvdw){
					worstvdw = vdw;
					worstvdwSeq = seq;
				}
				if (imm1>worstimm1){
					worstimm1 = imm1;
					worstimm1Seq = seq;
				}
				/******************************************************************************
				 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
				 ******************************************************************************/
				printSeqFile(PS, seq, AAener, totalEnergy, hbond, vdw, imm1, a, nout);
				printEnerFile(seq, AAener, a, rout);
				seqAccept++;
				sys.reset();
			}
			else{
				if (totalEnergy > 0){
					cout << "Total Energy " << totalEnergy << " is too high; go to next random sequence." << endl;
				}
				seqDiscard++;
				sys.reset();
			}
		}
	}
	if (seqDiscard > 0){
		cout << "Total sequences not accepted: " << seqDiscard << endl;
		nout << "Total sequences not accepted: " << seqDiscard << endl;
	}
	
	time(&endTime);
	diffTime = difftime (endTime, startTime);

	nout << endl << "Best and Worsts: " << endl << endl;
	
	nout << "Best Energy: " << bestEnergy << ": " << bestSeq << endl;
	nout << "Best vdW: " << bestvdw << ": " << bestvdwSeq << endl;
	nout << "Best IMM1: " << bestimm1 << ": " << bestimm1Seq << endl;
	nout << "Worst Energy: " << worstEnergy << ": " << worstSeq << endl;
	nout << "Worst vdW: " << worstvdw << ": " << worstvdwSeq << endl;
	nout << "Worst IMM1: " << worstimm1 << ": " << worstimm1Seq << endl;
	
	nout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	
	writer.close();
	writer2.close();
	pout.close();
	nout.close();
	rout.close();
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

	//opt.allowed.push_back("");
	
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	opt.allowed.push_back("seqNumber");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("threadBool");
	
	//side-chain repack variable
	opt.allowed.push_back("mcCycles");
	opt.allowed.push_back("mcMaxRejects");
	opt.allowed.push_back("mcStartTemp");
	opt.allowed.push_back("mcEndTemp");
	opt.allowed.push_back("mcCurve");

	opt.allowed.push_back("SL");
	
	// energy weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("start");
	opt.allowed.push_back("end");
	
	opt.allowed.push_back("ener");
	
	opt.allowed.push_back("ivalues");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	// input files
	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("solvFile");
	
	// alternate ids and weights
	opt.allowed.push_back("ids");
	opt.allowed.push_back("weights");

	// housekeeping variables
	opt.allowed.push_back("runNumber");
	opt.allowed.push_back("rejectvdw");

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
	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()) {
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}

	opt.rejectvdw = OP.getBool("rejectvdw");
	if (OP.fail()) {
		opt.warningMessages += "rejectvdw not specified using false\n";
		opt.warningFlag = true;
		opt.rejectvdw = false;
	}

	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages += "sequence (1 letter aa) not specified\n";
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

	opt.rulesFile = OP.getString("rulesFile");
	if (OP.fail()) {
		opt.rulesFile = "/data01/sabs/tmRepacks/GLY_69_Homo_2/tmRules/rules_10kcals_vdw_only/tmRules.out";
		opt.warningMessages += "rulesFile not specified using " + opt.rulesFile + "\n";
		opt.warningFlag = true;
	}

	opt.sequenceStart = OP.getInt("sequenceStart");
	if (OP.fail()) {
		opt.warningMessages += "sequenceStart not specified using 1\n";
		opt.warningFlag = true;
		opt.sequenceStart = 1;
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

	// MC
	opt.mcCycles = OP.getInt("mcCycles");
	if (OP.fail()) {
		opt.errorMessages += "Number of MC cycles not specified!\n";
		opt.errorFlag = true;
	}

	opt.mcMaxRejects = OP.getInt("mcMaxRejects");
	if (OP.fail()) {
		opt.mcMaxRejects = 10;
		opt.warningMessages += "Number of MC max rejects not specified, default to using 10\n";
		opt.warningFlag = true;
	}

	opt.mcStartTemp = OP.getDouble("MCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCStartTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.mcStartTemp = 1000.0;
	}
	opt.mcEndTemp = OP.getDouble("MCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.mcEndTemp = 0.5;
	}
	opt.mcCurve = OP.getInt("MCCurve");
	if (OP.fail()) {
		opt.warningMessages += "MCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.mcCurve = 2;
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
	opt.weight_solv = OP.getDouble("weight_solv");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_solv not specified, default 1.0\n";
		opt.weight_solv = 1.0;
	}

	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	}
	else{
		opt.SL = "SL"+opt.SL;
	}

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

	opt.start = OP.getInt("start");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "Start point not specified, default to 0\n";
		opt.start = 0;
	}
	opt.end = OP.getInt("end");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "End point not specified, default to 60\n";
		opt.end = 60;
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
	opt.monoRotLibFile = OP.getString("monoRotLibFile");
	if (OP.fail()) {
		opt.warningMessages += "monoRotLibFile not specified using " + opt.rotLibFile + "\n";
		opt.warningFlag = true;
		opt.monoRotLibFile = opt.rotLibFile;
	}

	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine backboneCrd";
		opt.errorFlag = true;
	}
	
	opt.hBondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hBondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hbondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hbondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
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
	opt.weights = OP.getDoubleVector("weights");
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
