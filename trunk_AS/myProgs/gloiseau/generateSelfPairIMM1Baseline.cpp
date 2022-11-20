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
#include "functions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "generateSelfPairIMM1Baseline";
string programDescription = "This program generate baselines for my sequences: I was previously using only Leucine ends, but because we use LILI as the C terminal of our helices, I needed to generate new baselines for this termini";
string programAuthor = "Gilbert Loiseau";
string programVersion = "3";
string programDate = "13 October 2022";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;
//TODO: edited but haven't tested yet
/******************************************************************************************************************************************************************************/

struct Options{
	string backboneAA;
	int backboneLength;
	int thread;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	int startResNum;
	int endResNum;
	int sequenceStart;

	bool deleteTerminalHbonds;
	
	string SL; //number of rotamers

	bool transform;
	
	// input files
	string helixGeoFile;
	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string solvFile;
	string hbondFile;
	string rotLibFile;
	string infile;
	string helicalAxis;

	int seqNumber;
	bool verbose;
	int greedyCycles;
	int seed;

	// energy weights
	double weight_vdw;
	double weight_hbond;
	
	// alternate identities and weights
	vector<string> ids;
	vector<int> weights;

	//Pair parameters
	int pairDist;

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

string convertToPolymerSequenceNeutralPatchMono(string _seq, int _startResNum) {
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

void deleteTerminalHydrogenBondInteractions(System &_sys, Options& _opt) {
	// look at all hbond interactions and remove those
	// remove any interaction that has a donor or acceptor from residues 1 2 3 and n-2 n-1 n on each chain that are not part of the TM

	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize(); // number of positions in chain A = chain B
	// compute the extensions at the beginning and the end
	//int frontExt = _opt.tmStart - _opt.startResNum;
	//int endExt = _opt.endResNum - _opt.tmEnd;
	int frontExt = _opt.tmStart;
	int endExt = _opt.endResNum;
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain& thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(frontExt > i) {
				atoms += positions[i]->getAtomPointers(); 
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			if(endExt > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers(); 
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");
}

Options parseBaselineOptions(int _argc, char * _argv[]);

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

void loadRotamers1(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.identitySize() > 1){
			for (uint j=0; j < pos.getNumberOfIdentities(); j++){
				pos.setActiveIdentity(j);
				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "ALA-ACE" && pos.getResidueName() != "ALA-CT2" && pos.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
					}
				}
				pos.setActiveIdentity(0);
			}
		}
		else{
			//if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "ALA-ACE" && pos.getResidueName() != "ALA-CT2" && pos.getResidueName() != "PRO") {
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

string randomAASeqChangedTermini(vector<string> _weightedAAs, int _seqLength, RandomNumberGenerator &_RNG){
	string seq = "";
	//vector<int> randAA;
	for (uint i=0; i<_seqLength; i++){
		if (i < 4 || i > _seqLength-5){
			if (i == _seqLength-3 || i == _seqLength-1){
				seq += "A";
			}
			else{
			   	seq += "A";
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

vector<double> calcBaselinePositionEnergies(System &_sys, int _seqLength, double &_totalProt, int _thread){
	vector <double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	sel.select("allProt, all");
	for (uint i=_thread; i<_seqLength+_thread; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sel.select(residue += number);
		double prot = _sys.calcEnergy("resi", "allProt")/2;//why divide by 2?
		double resi = _sys.calcEnergy("resi")/2; 
		double protEnergy = prot+resi;
		ener.push_back(protEnergy);
		_totalProt += protEnergy;
	}
	sel.clearStoredSelections();
	return ener;
}

vector<double> calcBaselineEnergies(System &_sys, int _seqLength, double &_totalProt, int _thread){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	for (uint i=_thread; i<_seqLength+_thread; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		//cout << i << ": " << resi << endl;
		ener.push_back(resi);
		_totalProt += resi;
	}
	sel.clearStoredSelections();
	return ener;
}

vector<double> calcPairBaselineEnergies(System &_sys, Options _opt, int _seqLength, double &_totalProt, int _thread){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=_thread; i<_seqLength+_thread; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i);
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength+_thread; j++){
			int dist = j-i;
			if (dist <= _opt.pairDist){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
				_totalProt += pair;
				//cout << residue << "; " << resi1 << " : " << pair << endl;
			} else {
				j = _seqLength+_thread;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, double _totalEnergy, double _hbond, double _vdw, int _seqNumber, bool _selfEner, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	if (_selfEner == true){
		_out << "AA\tPosition\tEnergy" << endl;//could be interesting to add rotamer number to this
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << "\t" << i+1 << "\t" << _ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		_out << "Pair\tDistance\tEnergy" << endl;
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				stringstream ss;
				ss <<_seq[i] << _seq[j];
				string p = ss.str();
				int d = j-i;
				pairs.push_back(p);
				dists.push_back(d);
			}
		}
		string spacer = "\t";
		for (uint i=0; i<pairs.size(); i++){
			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << endl;
		}
	}
	_out << "Total Energy: " << _totalEnergy << endl;
	_out << "H-bond Energy: " << _hbond << endl;
	_out << "VDW Energy: " << _vdw << endl << endl;
}

void printSeqFile(PolymerSequence &_PS, string _seq, vector<double> &_ener, vector<double> &_hEner, vector<double> &_vdWEner, double _totalEnergy, double _hbond, double _vdw, int _seqNumber, bool _selfEner, ofstream &_out){
	_out << "Sequence: " << _seqNumber+1 << endl;
	_out << _seq << endl;
	_out << _PS;
	if (_selfEner == true){
		_out << "AA\tPosition\tEnergy\tHBond\tvdW" << endl;//could be interesting to add rotamer number to this
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << "\t" << i+1 << "\t" << _ener[i] << "\t" <<  _hEner[i] << "\t" << _vdWEner[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		_out << "Pair\tDistance\tEnergy\tHBond\tvdW" << endl;
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				stringstream ss;
				ss <<_seq[i] << _seq[j];
				string p = ss.str();
				int d = j-i;
				pairs.push_back(p);
				dists.push_back(d);
			}
		}
		string spacer = "\t";
		for (uint i=0; i<pairs.size(); i++){
			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << spacer << _hEner[i] << spacer << _vdWEner[i] << endl;
		}
	}
	_out << "Total Energy: " << _totalEnergy << endl;
	_out << "H-bond Energy: " << _hbond << endl;
	_out << "VDW Energy: " << _vdw << endl << endl;
}

void printEnerFile(string _seq, vector<double> &_ener, int _seqNumber, bool _selfEner, ofstream &_out){
	if (_selfEner == true){
		if (_seqNumber == 0){
			_out << "AA\tPosition\tEnergy" << endl;//could be interesting to add rotamer number to this
		}
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << "\t" << i+1 << "\t" << _ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		if (_seqNumber == 0){
			_out << "Pair\tDistance\tEnergy" << endl;//could be interesting to add rotamer number to this
		}
		vector<string> pairs;
		vector<double> dists;
		for (uint i=0; i<_seq.length(); i++){
			for (uint j=i+1; j<_seq.length(); j++){
				stringstream ss;
				ss <<_seq[i] << _seq[j];
				string p = ss.str();
				int d = j-i;
				pairs.push_back(p);
				dists.push_back(d);
			}
		}
		string spacer = "\t";
		for (uint i=0; i<pairs.size(); i++){
			_out << pairs[i] << spacer << dists[i] << spacer << _ener[i] << endl;
		}
	}
}

void printEnerFile(string _seq, Options _opt, vector<double> &_ener, vector<double> &_hEner, vector<double> &_vdWEner, vector<double> &_imm1Ener, int _seqNumber, bool _selfEner, ofstream &_out){
	if (_selfEner == true){
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "AA\tPosition\tEnergy\tHbond\tvdW\tIMM1" << endl;//could be interesting to add rotamer number to this
		}
		for (uint i=0; i<_seq.length(); i++){
			_out << _seq[i] << "\t" << i+1 << "\t" << _ener[i] << "\t" <<  _hEner[i] << "\t" << _vdWEner[i] << "\t" << _imm1Ener[i] << endl;
		}//I think the easiest way is to run and append the list for each sequence, then put these all together in an excel spreadsheet and organize there
	}
	else{
		if (_seqNumber == 0 && _opt.seed == 0){
			_out << "Pair\tDistance\tPosition1\tPosition2\tEnergy\tHbond\tvdW\tIMM1" << endl;//could be interesting to add rotamer number to this
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
		string sp = "\t";
		int pos1 = 0;

		for (uint i=0; i<pairs.size(); i++){
			if (dists[i] == 1){
				pos1++;
			}
			int pos2 = pos1 + dists[i];
			_out << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << _hEner[i] << sp << _vdWEner[i] << sp << _imm1Ener[i] << endl;
			//cout << pairs[i] << sp << dists[i] << sp << pos1 << sp << pos2 << sp << _ener[i] << sp << _hEner[i] << sp << _vdWEner[i] << endl;
		}
	}
}
//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

vector<vector<string>> positionToString(System &_sys, vector<int> &_varPos){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_varPos.size()/2; k++){
		//cout << "string" << endl;
		if (_varPos[k] == 1){
			vector<string> tempPos;

			Position &posA = _sys.getPosition(k);
			Position &posB = _sys.getPosition(k+_varPos.size()/2);

			string A = posA.toString();
			string B = posB.toString();

			string delimiter = " ";
			
			size_t p = 0;
			p = A.find(delimiter);

			tempPos.push_back(A.substr(0, p));
			tempPos.push_back(B.substr(0, p));

			stringPositions.push_back(tempPos);
		}
		else{
			continue;
		}
	}
	return stringPositions;
}

vector<double> getSubtractBetweenVectors(vector<double> &_v1, vector<double> &_v2){
	vector<double> v3;
	if (_v1.size() != _v2.size()){
		cout << "Vectors of different sizes; terminate program." << endl;
	}
	else{
		for (uint i=0; i<_v1.size(); i++){
			v3.push_back(_v1[i]-_v2[i]);
		}
	}
	return v3;
}

vector<double> getAdditionBetweenVectors(vector<double> &_v1, vector<double> &_v2){
	vector<double> v3;
	if (_v1.size() != _v2.size()){
		cout << "Vectors of different sizes; terminate program." << endl;
	}
	else{
		for (uint i=0; i<_v1.size(); i++){
			v3.push_back(_v1[i]+_v2[i]);
		}
	}
	return v3;
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

	time (&startTime);
	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options opt = parseBaselineOptions(argc, argv);
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
	ofstream sout;
	ofstream pout;
	ofstream ssout;
	ofstream spout;
	
	string dir = opt.pdbOutputDir + "/" + date;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string soutName  = dir + "/Self_Energies_" + opt.runNumber + ".out";
	string poutName  = dir + "/Pair_Energies_" + opt.runNumber + ".out";
	string ssoutName = dir + "/Sequences_Self_" + opt.runNumber + ".out";
	string spoutName = dir + "/Sequences_Pair_" + opt.runNumber + ".out";
	
	sout.open(soutName.c_str());
	pout.open(poutName.c_str());
	ssout.open(ssoutName.c_str());
	spout.open(spoutName.c_str());
	
	/******************************************************************************
	 *                    === INITIALIZE POLYGLY AS BACKBONE ===
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
	writer.open(dir + "/Run" + opt.runNumber + ".pdb");
	
	/******************************************************************************
	 *                      === HOUSEKEEPING VARIABLES ===
	 ******************************************************************************/
	int seqDiscard = 0;
	int seqAccept = 0;
	int seqNumber = opt.seqNumber;
	
	//Random Number Generator
	RandomNumberGenerator RNG;
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
	 *                         === HELICAL AXIS SET UP ===
	 ******************************************************************************/
	// System for the helical axis that sets protein around the origin (0.0, 0.0, 0.0)
	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxis);

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	// Set up object used for transformations
	Transforms trans;
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// save starting coordinates of the helical axis
	helicalAxis.saveAltCoor("start");
	/******************************************************************************
	 *               === BEGIN LOOP FOR RANDOMIZING SEQUENCES ===
	 ******************************************************************************/
	while (seqAccept < seqNumber){
		cout << "Number Accepted: " << seqAccept << endl;
		cout << "Number Discarded: " << seqDiscard << endl;
		
		/******************************************************************************
		 *                         === RANDOMIZE SEQUENCES ===
		 ******************************************************************************/
		string seq = randomAASeqChangedTermini(weightedAAs, opt.backboneLength, RNG);//To have leucine ends for my sequences, similar to how it will be for my designs (basically getting baselines for those ends)
		cout << seq << endl;
		string polySeq = convertToPolymerSequenceNeutralPatchMono(seq, opt.thread);//To have neutral patches at the ends of my sequences to prevent any hydrogen bonding that could affect energy score
		PolymerSequence PS(polySeq);
		cout << PS << endl;

		/******************************************************************************
		 *                           === DECLARE SYSTEM ===
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
		
		// Reference points for Helices
		CartesianPoint ori(0.0,0.0,0.0);
		CartesianPoint xAxis(1.0,0.0,0.0);
		CartesianPoint zAxis(0.0,0.0,1.0);

		// load starting coordinates of the helical axis
		helicalAxis.applySavedCoor("start");

		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.assignCoordinates(glyAPV,false);
		moveZCenterOfCAMassToOrigin(sys.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
		sys.buildAllAtoms();
		
		CSB.updateNonBonded(10,12,50); //eventually change these to options for cuton, cutoff, and cutnb
		
		SystemRotamerLoader sysRot(sys, opt.rotLibFile);
		sysRot.defineRotamerSamplingLevels();
		
		// Add hydrogen bond term
		HydrogenBondBuilder hb(sys, opt.hbondFile);
		hb.buildInteractions(50);//set this to same as cutnb
		
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
		cout << "All atoms have coordinates" << endl;
		writer.write(sys.getAtomPointers(), true, false, true);

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
		
		// Set weights
		Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
		Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
		Eset->setWeight("CHARMM_IMM1REF", 1);
		Eset->setWeight("CHARMM_IMM1", 1);

		/******************************************************************************
		 *            === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
		 ******************************************************************************/
		deleteTerminalHydrogenBondInteractions(sys,opt);

		/******************************************************************************
		 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
		 ******************************************************************************/
		CSB.updateNonBonded(10,12,50);
		sys.buildAllAtoms();

		loadRotamers1(sys, sysRot, opt.SL);
		
		SelfPairManager spm;
		spm.seed(RNG.getSeed());
		spm.setSystem(&sys);
		spm.setVerbose(false);
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
		cout << Eset->getSummary() << endl;
		double totalEnergy = sys.calcEnergy();
		double vdw = Eset->getTermEnergy("CHARMM_VDW");
		double hbond = Eset->getTermEnergy("SCWRL4_HBOND");
		double imm1 = Eset->getTermEnergy("CHARMM_IMM1REF") + Eset->getTermEnergy("CHARMM_IMM1");
		
		/******************************************************************************
		 *               === CALCULATE ENERGIES FOR EACH POSITION ===
		 ******************************************************************************/
		double totalProt = 0;
		double totalSelfPair = 0;
		cout <<"Sequence Length: " << seq.length() << endl;
		vector<double> AAener = calcBaselineEnergies(sys, seq.length(), totalSelfPair, opt.thread);
		vector<double> AAPosEner = calcBaselinePositionEnergies(sys, seq.length(), totalProt, opt.thread);
		vector<double> pairAAEner = calcPairBaselineEnergies(sys, opt, seq.length(), totalSelfPair, opt.thread);
		
		/******************************************************************************
		 *      === CALCULATE ENERGIES FOR EACH POSITION AFTER ERASING A TERM ===
		 ******************************************************************************/
		Eset->eraseTerm("CHARMM_VDW");
	
		vector<double> AAHbondIMM1 = calcBaselineEnergies(sys, seq.length(), totalProt, opt.thread);
		vector<double> AAPosHbondIMM1 = calcBaselinePositionEnergies(sys, seq.length(), totalProt, opt.thread);
		vector<double> pairAAHbondIMM1 = calcPairBaselineEnergies(sys, opt, seq.length(), totalProt, opt.thread);
		cout << 4 << endl;
		
		/******************************************************************************
		 *                 === CALCULATE ENERGIES FOR ERASED TERM ===
		 ******************************************************************************/
		Eset->eraseTerm("SCWRL4_HBOND");

		vector<double> AAIMM1 = calcBaselineEnergies(sys, seq.length(), totalProt, opt.thread);
		vector<double> AAPosIMM1 = calcBaselinePositionEnergies(sys, seq.length(), totalProt, opt.thread);
		vector<double> pairAAIMM1 = calcPairBaselineEnergies(sys, opt, seq.length(), totalProt, opt.thread);
		cout << 5 << endl;
		
		/******************************************************************************
		 *                 === CALCULATE ENERGIES FOR ERASED TERM ===
		 ******************************************************************************/
		vector<double> AAvdW     = getSubtractBetweenVectors(AAener, AAHbondIMM1);
		vector<double> AAPosvdW  = getSubtractBetweenVectors(AAPosEner, AAPosHbondIMM1);
		vector<double> pairAAvdW = getSubtractBetweenVectors(pairAAEner, pairAAHbondIMM1);

		/******************************************************************************
		 *                 === CALCULATE ENERGIES FOR ERASED TERM ===
		 ******************************************************************************/
		vector<double> AAHbond     = getSubtractBetweenVectors(AAHbondIMM1, AAIMM1);
		vector<double> AAPosHbond  = getSubtractBetweenVectors(AAPosHbondIMM1, AAPosIMM1);
		vector<double> pairAAHbond = getSubtractBetweenVectors(pairAAHbondIMM1, pairAAIMM1);

		/******************************************************************************
		 *        === TROUBLESHOOTING TO SEE IF ENERGIES ADD UP PROPERLY ===
		 ******************************************************************************/
		//vector<double> AAtot     = getSubtractBetweenVectors(AAvdW, AAHbond);
		//vector<double> AAPostot  = getSubtractBetweenVectors(AAPosvdW, AAPosHbond);
		//vector<double> pairAAtot = getSubtractBetweenVectors(pairAAvdW, pairAAHbond);
		
		/******************************************************************************
		 *                 === OUTPUT ENERGIES FOR EACH TERM ===
		 ******************************************************************************/
		double selfEner         = sumEnergyVector(AAener);
		double selfEnerHbond    = sumEnergyVector(AAHbond);
		double selfEnervdW      = sumEnergyVector(AAvdW);
		double selfEnerIMM1     = sumEnergyVector(AAIMM1);
		//double selfPosEner      = sumEnergyVector(AAPosEner);
		//double selfPosEnerHbond = sumEnergyVector(AAPosHbond);
		//double selfPosEnervdW   = sumEnergyVector(AAPosvdW);
		double pairEner         = sumEnergyVector(pairAAEner);
		double pairEnerHbond    = sumEnergyVector(pairAAHbond);
		double pairEnervdW      = sumEnergyVector(pairAAvdW);
		double pairEnerIMM1     = sumEnergyVector(pairAAIMM1);
		
		cout << "Total Energy: " << totalEnergy << endl;
		cout << "Hbond Energy: " << hbond << endl;
		cout << "vdW Energy: " << vdw << endl;
		cout << "IMM1 Energy: " << imm1 << endl;
		cout << "Self Energy: " << selfEner << endl;
		cout << "Self Energy Hbond: " << selfEnerHbond << endl;
		cout << "Self Energy vdW: " << selfEnervdW << endl;
		cout << "Self Energy IMM1: " << selfEnerIMM1 << endl;
		cout << "Pair Energy: " << pairEner << endl;
		cout << "Pair Energy Hbond: " << pairEnerHbond << endl;
		cout << "Pair Energy vdW: " << pairEnervdW << endl;
		cout << "Pair Energy IMM1: " << pairEnerIMM1 << endl;
		//cout << "Self Position Energy: " << selfPosEner << endl;
		cout << "Self Pair Energy: " << selfEner + pairEner << endl;
		cout << "Self Pair Energy: " << totalSelfPair << endl;

		/******************************************************************************
		 *           === PRINT BASELINE AA ENERGIES INTO OUTPUT FILES ===
		 ******************************************************************************/
		printSeqFile(PS, seq, AAener, AAHbond, AAvdW, totalEnergy, hbond, vdw, seqAccept, true, ssout);
		printEnerFile(seq, opt, AAener, AAHbond, AAvdW, AAIMM1, seqAccept, true, sout);
		printSeqFile(PS, seq, pairAAEner, pairAAHbond, pairAAvdW, totalEnergy, hbond, vdw, seqAccept, false, spout);
		printEnerFile(seq, opt, pairAAEner, pairAAHbond, pairAAvdW, pairAAIMM1, seqAccept, false, pout);
		seqAccept++;
		sys.reset();
	}
	time (&endTime);
	diffTime = difftime(endTime, startTime);
	ssout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	spout << "Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << "seconds" << endl;
	writer.close();
	sout.close();
	pout.close();
	ssout.close();
	spout.close();
}

Options parseBaselineOptions(int _argc, char * _argv[]){
	
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
	opt.allowed.push_back("thread");
	opt.allowed.push_back("seqNumber");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("threadBool");
	
	opt.allowed.push_back("SL");
	
	// energy weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	// input files
	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("helicalAxis");
	
	// alternate ids and weights
	opt.allowed.push_back("ids");
	opt.allowed.push_back("weights");

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

	opt.SL = OP.getString("rotLevel");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
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
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "thread not specified, default to 23\n";
		opt.thread = 23;
	}
	opt.seqNumber = OP.getInt("seqNumber");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "seqNumber not specified, default to 10\n";
		opt.seqNumber = 100;
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
	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.warningMessages += "helicalAxis not specified\n";
		opt.warningFlag = true;
		opt.helicalAxis = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/parameterFiles/helicalAxis.pdb";
	}

	opt.rerunConf = OP.getConfFile();

	return opt;
}
