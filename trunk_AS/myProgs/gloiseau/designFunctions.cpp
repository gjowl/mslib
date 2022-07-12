/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/12
 * @Email:  gjowl04@gmail.com
 * @Filename: design.cpp
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/22
 */

#include <sstream>
#include <iterator>
#include <unistd.h>
#include "designFunctions.h"
#include "homodimerFunctions.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

/***********************************
 *string output functions
 ***********************************/
string generateMonomerMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
	// A:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
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
			} else {
				if(resName == "HIS") {
					ps = ps + " HSE-CT2";
				} else {
					ps = ps + " " + resName + "-CT2";
				}
			}
			counter++;
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if (_interfacialPositions[counter] == 1){
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
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}


/***********************************
 *  calculate residue burial
 ***********************************/

//TODO: fix this function and make it less gross

/***********************************
 *output file functions
 ***********************************/
//TODO: make changes so that this can be run locally vs external server
void setupDesignDirectory(Options &_opt, string _date){
	//_opt.pdbOutputDir = _opt.pdbOutputDir + "/" + _date;
	_opt.pdbOutputDir = get_current_dir_name();
	_opt.pdbOutputDir = _opt.pdbOutputDir + "/12_06_2021";//Had problems with dynamic date because of runs starting one day and ending another
	string cmd = "mkdir -p " + _opt.pdbOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	_opt.pdbOutputDir = _opt.pdbOutputDir + "/design_" + _opt.runNumber;
	cmd = "mkdir -p " + _opt.pdbOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

void outputEnergyFile(Options &_opt, string _interface, vector<string> _allDesigns){
	ofstream eout;
	string eoutfile = _opt.pdbOutputDir + "/energyFile.csv";
	eout.open(eoutfile.c_str());
	string tab = "\t";

	// Outputs
	eout << tab << endl;
	for (uint i=0; i<_opt.energyTermsToOutput.size(); i++){
		eout << _opt.energyTermsToOutput[i] << tab;
		cout << _opt.energyTermsToOutput[i] << tab;
	}
	eout << "Starting Sequence" << tab;
	eout << "Baseline" << tab;
	eout << "Sequence" << tab;
	eout << "InterfaceSeq" << tab;
	eout << "xShift" << tab;
	eout << "crossingAngle" << tab;
	eout << "axialRotation" << tab;
	eout << "zShift" << tab;
	eout << "angleDistDensity" << tab;
	eout << "axialRotationDensity" << tab;
	eout << "zShiftDensity" << tab;
	eout << "repackLevels" << tab;
	eout << "interfaceLevels" << tab;
	eout << "backboneLength" << tab;
	eout <<	"PDBPath" << endl;

	cout << "Starting Sequence" << tab;
	cout << "Baseline" << tab;
	cout << "Sequence" << tab;
	cout << "InterfaceSeq" << tab;
	cout << "xShift" << tab;
	cout << "crossingAngle" << tab;
	cout << "axialRotation" << tab;
	cout << "zShift" << tab;
	cout << "angleDistDensity" << tab;
	cout << "axialRotationDensity" << tab;
	cout << "zShiftDensity" << tab;
	cout << "repackLevels" << tab;
	cout << "interfaceLevels" << tab;
	cout << "backboneLength" << tab;
	cout <<	"PDBPath" << endl;

	for (uint i=0; i<_allDesigns.size() ; i++){
		eout << _allDesigns[i] << endl;
		cout << _allDesigns[i] << endl;
	}
	eout.close();
}

//void makeRepackConfig(Options &_opt, string _sequence, vector<uint> _state, string _designNumber, string _pdbPath, ofstream &_globalRepackOut){
void makeRepackConfig(Options &_opt, string _sequence, string _designDir, vector<uint> _state, string _seqNumber, string _pdbPath, string _crdPath, map<string,double> _energyMap, vector<int> _rotamerSampling){
	ofstream dout;
	string doutfile = _opt.pdbOutputDir + "/" + _sequence +"/repack.config";
	dout.open(doutfile.c_str());
	//string designFileDir = "/exports/home/gloiseau/mslib/trunk_AS/designFiles/";//Running on chtc server and needed these to access the proper directories after transferring them back to our server
	//TODO: make sure to add the designFiles directory to this directory
	string designFileDir = "/data02/gloiseau/Sequence_Design_Project/vdwSequenceDesign/sequenceDesign/12_06_2021/";

	//TODO: also when I run this locally, the output for input files is poor. I'll need to fix these outputs in the future (maybe just add an inputfile directory?
	dout << "#Input Files" << endl;
	dout << "topFile            " << designFileDir + _opt.topFile << endl;
	dout << "parFile            " << designFileDir + _opt.parFile << endl;
	dout << "solvFile           " << designFileDir + _opt.solvFile << endl;
	dout << "hbondFile          " << designFileDir + _opt.hbondFile << endl;
	dout << "rotLibFile         " << designFileDir + _opt.rotLibFile << endl;
	dout << "designPdb          " << _pdbPath << endl;
	dout << "designCrd          " << _crdPath << endl;
	dout << "outputDir          " << _designDir << endl << endl;

	dout << "#Geometry" << endl;
	dout << "xShift             " << _opt.xShift << endl;
	dout << "crossingAngle      " << _opt.crossingAngle << endl;
	dout << "axialRotation      " << _opt.axialRotation << endl;
	dout << "zShift             " << _opt.zShift << endl;
	dout << "thread             " << _opt.thread << endl;
	dout << "tmStart            " << _opt.tmStart << endl;
	dout << "tmEnd              " << _opt.tmEnd << endl << endl;

	dout << "#Design Parameters" << endl;
	dout << "sequence              " << _sequence << endl;
	dout << "rotamerSamplingString ";
	for (uint i=0; i<_rotamerSampling.size()/2; i++){
		if (i<(_rotamerSampling.size()/2)-1){
			dout << _rotamerSampling[i];
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << "seed               " << _opt.seed << endl;
	dout << endl;
	dout << "#Rotamer Sampling Vector" << endl;
	dout << "rotamerSampling           ";
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (i<_rotamerSampling.size()-1){
			dout << _rotamerSampling[i] << " ";
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << endl;

	// SASA Rotamer Values
	dout << "#SASA Rotamers" << endl;
	for (uint i=0; i<_opt.sasaRepackLevel.size(); i++){
		dout << "sasaRepackLevel  " << _opt.sasaRepackLevel[i] << endl;
	}
	dout << endl;

	// Weights
	dout << "#Weights" << endl;
	dout << "weight_vdw         " << _opt.weight_vdw << endl;
	dout << "weight_hbond       " << _opt.weight_hbond << endl;
	dout << "weight_solv        " << _opt.weight_solv << endl << endl;

	// Monomer Energies
	dout << "#Monomer Energies" << endl;
	dout << "monomer            " << _energyMap.at("Monomer") << endl;
	dout << "monoVdw            " << _energyMap.at("VDWMonomer") << endl;
	dout << "monoHbond          " << _energyMap.at("HBONDMonomer") << endl;
	dout << "monoIMM1           " << _energyMap.at("IMM1Monomer") << endl << endl;

	// Energies before repack
	dout << "#Dimer Energies" << endl;
	dout << "dimer            " << _energyMap.at("Dimer") << endl;
	dout << "dimerVdw         " << _energyMap.at("VDWDimer") << endl;
	dout << "dimerHbond       " << _energyMap.at("HBONDDimer") << endl;
	dout << "dimerIMM1        " << _energyMap.at("IMM1Dimer") << endl << endl;

	dout << "#Total Energy" << endl;
	dout << "Total            " << _energyMap.at("Total") << endl;
	dout.close();
}

void makeDockingConfig(Options &_opt, string _sequence, string _designDir, string _pdbPath, map<string,double> _energyMap, vector<int> _rotamerSampling){
	ofstream dout;
	string outputDir = _opt.pdbOutputDir + "/" + _sequence;
	string doutfile = outputDir + "/docking.config";
	dout.open(doutfile.c_str());
	//TODO: fix this spacing output so that it's cleaner (I think setw(21)?)
	//TODO: also when I run this locally, the output for input files is poor. I'll need to fix these outputs in the future
	dout << "#Design Parameters" << endl;
	dout << "chainSeq              " << _sequence << " " << _sequence << endl;
	dout << "rotamerSamplingString ";
	for (uint i=0; i<_rotamerSampling.size()/2; i++){
		if (i<(_rotamerSampling.size()/2)-1){
			dout << _rotamerSampling[i];
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << "chainStartNum        " << _opt.thread << endl;
	dout << "seed                 " << _opt.seed << endl << endl;

	string designFileDir = "/data02/gloiseau/Sequence_Design_Project/vdwSequenceDesign/sequenceDesign/12_06_2021/";
	dout << "#Input Files" << endl;
	dout << "topFile              " << designFileDir + _opt.topFile << endl;
	dout << "parFile              " << designFileDir + _opt.parFile << endl;
	dout << "solvFile             " << designFileDir + _opt.solvFile << endl;
	dout << "hbondFile            " << designFileDir + _opt.hbondFile << endl;
	dout << "rotLibFile           " << designFileDir + _opt.rotLibFile << endl;
	dout << "bbqFile              " << designFileDir << "PiscesBBQTable.txt" << endl;
	dout << "designPdb            " << _pdbPath << endl;//In case I ever wnat to do an RMSD calculation against the orignal
	dout << "outputDir            " << _designDir << endl << endl;

	dout << "#Geometry" << endl;
	dout << "xShift               " << _opt.xShift << endl;
	dout << "crossingAngle        " << _opt.crossingAngle << endl;
	dout << "axialRotation        " << _opt.axialRotation << endl;
	dout << "zShift               " << _opt.zShift << endl;
	dout << endl;

	dout << "#Rotamer Sampling Vector" << endl;
	dout << "rotamerSampling      ";
	for (uint i=0; i<_rotamerSampling.size()/2; i++){
		if (i<(_rotamerSampling.size()/2)-1){
			dout << _rotamerSampling[i] << " ";
		} else {
			dout << _rotamerSampling[i] << endl;
		}
	}
	dout << endl;

	dout << "#Weights" << endl;
	dout << "weight_vdw           " << _opt.weight_vdw << endl;
	dout << "weight_hbond         " << _opt.weight_hbond << endl;
	dout << "weight_solv          " << _opt.weight_solv << endl << endl;

	// Monomer Energies
	dout << "#Monomer Energies" << endl;
	dout << "monomer              " << _energyMap.at("Monomer") << endl;
	dout << "monoVdw              " << _energyMap.at("VDWMonomer") << endl;
	dout << "monoHbond            " << _energyMap.at("HBONDMonomer") << endl;
	dout << "monoIMM1             " << _energyMap.at("IMM1Monomer") << endl << endl;

	// Energies before repack
	dout << "#Dimer Energies" << endl;
	dout << "dimer                " << _energyMap.at("Dimer") << endl;
	dout << "dimerVdw             " << _energyMap.at("VDWDimer") << endl;
	dout << "dimerHbond           " << _energyMap.at("HBONDDimer") << endl;
	dout << "dimerIMM1            " << _energyMap.at("IMM1Dimer") << endl << endl;

	dout << "#Total Energy" << endl;
	dout << "Total            " << _energyMap.at("Total") << endl;
	dout.close();
}

void outputDesignFiles(Options &_opt, string _interface, vector<int> _rotamerSamplingPerPosition, vector<pair<string,vector<uint>>> _sequenceStatePair, map<string,map<string,double>> _sequenceEnergyMap, vector<double> _densities){
	// Setup vector to hold energy file lines
	vector<string> energyLines;

	// Setup the parameters for this specific run
	string tab = "\t";
	string xShift = MslTools::doubleToString(_opt.xShift);
	string crossingAngle = MslTools::doubleToString(_opt.crossingAngle);
	string axialRotation = MslTools::doubleToString(_opt.axialRotation);
	string zShift = MslTools::doubleToString(_opt.zShift);
	string angleDistDensity = MslTools::doubleToString(_densities[0]);
	string axialRotationDensity = MslTools::doubleToString(_densities[1]);
	string zShiftDensity = MslTools::doubleToString(_densities[2]);
	string thread = MslTools::intToString(_opt.thread);
	string repackLevels = MslTools::intToString(_opt.sasaRepackLevel.size());
	string interfaceLevel = MslTools::intToString(_opt.interfaceLevel);
	string bbLength = MslTools::intToString(_opt.backboneLength);
	//string runParameters = xShift+tab+crossingAngle+tab+axialRotation+tab+zShift+tab+thread+tab+repackLevels+tab+interfaceLevel;
	string runParameters = xShift+tab+crossingAngle+tab+axialRotation+tab+zShift+tab+angleDistDensity+tab+axialRotationDensity+tab+zShiftDensity+tab+repackLevels+tab+interfaceLevel+tab+bbLength;
	// For loop to setup the energy file
	string startSequence = _sequenceStatePair[0].first;
	for (uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		vector<uint> state = _sequenceStatePair[i].second;
		map<string,double> energyMap = _sequenceEnergyMap.at(sequence);
		// For adding in strings to a line for the energy file
		string tmp = "Sequence Info: ";
		for (uint j=0; j<_opt.energyTermsToOutput.size(); j++){
			string energyTerm = _opt.energyTermsToOutput[j];
			double energy = energyMap.at(energyTerm);
			//cout << sequence << ": " << energyTerm << " = " << energy << endl;
			tmp.append(MslTools::doubleToString(energy));
			tmp.append(tab);
		}

		//Add in path to design PDB and make repack and docking configuration files
		string seqNumber = MslTools::doubleToString(energyMap.at("SequenceNumber"));
		string designDir = "/data02/gloiseau/Sequence_Design_Project/vdwSequenceDesign/sequenceDesign/12_06_2021/design_" + _opt.runNumber + "/"+ sequence;
		string pdbPath = designDir + "/" + sequence + ".pdb";
		string crdPath = designDir + "/" + sequence + ".crd";
		makeRepackConfig(_opt, sequence, designDir, state, seqNumber, pdbPath, crdPath, energyMap, _rotamerSamplingPerPosition);
		makeDockingConfig(_opt, sequence, designDir, pdbPath, energyMap, _rotamerSamplingPerPosition);

		// Append other important features to the end energy files lines
		string interfaceSequence = getInterfaceSequence(_opt,_interface, sequence);
		tmp.append(startSequence);
		tmp.append(tab);
		tmp.append(sequence);
		tmp.append(tab);
		tmp.append(interfaceSequence);
		tmp.append(tab);
		tmp.append(runParameters);
		tmp.append(tab);
		tmp.append(pdbPath);
		energyLines.push_back(tmp);
		tmp.clear();
	}
	outputEnergyFile(_opt, _interface, energyLines);
}








/***********************************
 *baseline energy helper functions
 ***********************************/
//Function to calculate the self energies of a chain
vector<double> calcBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	//cout << "Self Energies" << endl;
	for (uint i=0; i<_seqLength; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i+1);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
		//cout << number << ": " << resi << endl;
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to calculate the pair energies of a chain
vector<double> calcPairBaselineEnergies(System &_sys, int _seqLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=0; i<_seqLength; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i+1);
		sel.select(residue += num1);
		for (uint j=i+1; j<_seqLength;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j+1);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			} else {
				j = _seqLength;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to get the sum of a vector of doubles, typically energies
double sumEnergyVector(vector<double> _energies){
	double ener = 0;
	for (uint i=0; i<_energies.size(); i++){
		ener = ener + _energies[i];
	}
	return ener;
}

/***********************************
 *calculate energies
 ***********************************/
void computeDimerEnergiesLinked(System &_sys, Options &_opt, map<string,map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err) {

	EnergySet* Eset = _sys.getEnergySet();
	_sys.calcEnergy();

	Eset->eraseTerm("BASELINE");
	Eset->eraseTerm("BASELINE_PAIR");
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	for(uint i=0; i<_sequenceStatePair.size(); i++){
		string sequence = _sequenceStatePair[i].first;
		vector<uint> stateVec = _sequenceStatePair[i].second;
		_sys.setActiveRotamers(stateVec);

		//Calculate Dimer energy for output
		double dimerEnergy = _sys.calcEnergy();
		double finalEnergy = dimerEnergy-_sequenceEnergyMap[sequence]["Monomer"];//TODO: make this more elegant by pulling from the map
		_sout << "Sequence: " << sequence << endl;
		_sout << "Dimer Energy: " << dimerEnergy << endl;
		_sout << "Final Energy = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl << endl;

		map<string,double> &energyMap = _sequenceEnergyMap[sequence];
		outputEnergiesByTermLinked(Eset, energyMap, _opt.energyTermList, "Dimer");

		//Setup directory for individual design
		string repackDir = _opt.pdbOutputDir + "/" + MslTools::intToString(i);
		string cmd = "mkdir -p " + repackDir;
		if (system(cmd.c_str())){
			_err << "Unable to make directory" << endl;
			exit(0);
		}
		_sequenceEnergyMap[sequence]["Dimer"] = dimerEnergy;
		_sequenceEnergyMap[sequence]["Total"] = finalEnergy;

		_sys.setActiveRotamers(stateVec);
		PDBWriter designWriter;
		designWriter.open(repackDir + "/" + MslTools::intToString(i) + ".pdb");
		designWriter.write(_sys.getAtomPointers(), true, false, true);
		designWriter.close();
		_writer.write(_sys.getAtomPointers(), true, false, true);
		saveEnergyDifference(_opt, _sequenceEnergyMap, sequence);

		//Write CRD
		CRDWriter writer;
		writer.open(repackDir + "/" + MslTools::intToString(i) + ".crd");
		writer.write(_sys.getAtomPointers(), false);
		writer.close();
	}
}

void computeDimerEnergies(System &_sys, Options &_opt, map<string, map<string,double>> &_sequenceEnergyMap, vector<pair<string,vector<uint>>> &_sequenceStatePair, vector<int> _rotamerSamplingPerPosition, vector<vector<string>> &_linkedPos, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	// Setup timer for LocalMC Cycles
	time_t startTimeLMC, endTimeLMC;
	double diffTimeLMC;
	time(&startTimeLMC);

	// Initialize PDBWriter for design
	PDBWriter writer;
	writer.open(_opt.pdbOutputDir + "/allDesigns.pdb");

	_sout << "Calculating dimer energies..." << endl;
	cout << "Calculating dimer energies..." << endl;
	if (_opt.linkInterfacialPositions){
		computeDimerEnergiesLinked(_sys, _opt, _sequenceEnergyMap, _sequenceStatePair, _rotamerSamplingPerPosition, _linkedPos, _RNG, writer, _sout, _err);
	} else {
		for (uint i=0; i<_sequenceStatePair.size(); i++){
			string sequence = _sequenceStatePair[i].first;
			vector<uint> state = _sequenceStatePair[i].second;
			cout << "Sequence " << i+1 << ": " << sequence << endl;
			_sout << "Sequence " << i+1 << ": " << sequence << endl;

			computeDimerEnergy(_sys, _opt, _sequenceEnergyMap, sequence, state, _rotamerSamplingPerPosition, _linkedPos, i, _RNG, writer, _sout, _err);
			_sequenceStatePair[i].second = state;
		}
	}

	time(&endTimeLMC);
	diffTimeLMC = difftime (endTimeLMC, startTimeLMC);
	writer.close();
	cout << "End dimer Calculations: " << diffTimeLMC << "s" << endl << endl;
	_sout << "End dimer Calculations: " << diffTimeLMC << "s" << endl << endl;
}

void computeDimerEnergy(System &_sys, Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_sequence, vector<uint> &_stateVec, vector<int> &_rotamerSampling, vector<vector<string>> &_linkedPos, int _seqNumber, RandomNumberGenerator &_RNG, PDBWriter &_writer, ofstream &_sout, ofstream &_err) {

	string polySeq = convertToPolymerSequenceNeutralPatch(_sequence, _opt.thread);
	PolymerSequence PL(polySeq);

	// Declare new system
	System sys;
	CharmmSystemBuilder CSB(sys, _opt.topFile, _opt.parFile, _opt.solvFile);
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
	if(!CSB.buildSystem(PL)) {
		_err << "Unable to build system from " << PL << endl;
		exit(0);
	}

	/******************************************************************************
	 *           === TRANSFORM HELICES TO INITIAL STARTING POSITION ===
	 ******************************************************************************/
	sys.assignCoordinates(_sys.getAtomPointers(),false);
	sys.buildAllAtoms();

	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);
	Eset->setTermActive("CHARMM_IMM1REF", true);
	Eset->setTermActive("CHARMM_IMM1", true);

	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);
	Eset->setWeight("CHARMM_IMM1REF", 1);
	Eset->setWeight("CHARMM_IMM1", 1);

	/******************************************************************************
	 *                === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(sys,_opt);

	/******************************************************************************
	 *              === LOAD ROTAMERS AND CHOOSE TO LINK INTERFACE ===
	 ******************************************************************************/
	loadRotamers(sys, sysRot, _opt, _rotamerSampling);
	CSB.updateNonBonded(10,12,50);

	/******************************************************************************
	 *                  === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	// Setup spm and calculate energies
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&sys);
	spm.setVerbose(false);
	spm.getMinStates()[0];
	spm.updateWeights();
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();

	//Calculate Dimer energy for output
	repackSideChains(spm, 10);
	_stateVec = spm.getMinStates()[0];
	sys.setActiveRotamers(_stateVec);

	double dimerEnergy = spm.getStateEnergy(_stateVec);
	double finalEnergy = dimerEnergy-_sequenceEnergyMap[_sequence]["Monomer"];//TODO: make this more elegant by pulling from the map
	_sout << "Dimer Energy: " << dimerEnergy << endl;
	_sout << "Final Energy w/ IMM1 = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;
	cout << "Dimer Energy: " << dimerEnergy << endl;
	cout << "Final Energy = Dimer - Monomer = " << dimerEnergy << "+" << _sequenceEnergyMap[_sequence]["Monomer"]*(-1) << "=" << finalEnergy << endl;

	map<string,double> &energyMap = _sequenceEnergyMap[_sequence];
	outputEnergiesByTerm(spm, _stateVec, energyMap, _opt.energyTermList, "Dimer", true);

	//Setup directory for individual design
	string repackDir = _opt.pdbOutputDir + "/" + MslTools::intToString(_seqNumber);
	string cmd = "mkdir -p " + repackDir;
	if (system(cmd.c_str())){
		_err << "Unable to make directory" << endl;
		exit(0);
	}
	_sequenceEnergyMap[_sequence]["Dimer"] = dimerEnergy;
	_sequenceEnergyMap[_sequence]["Total"] = finalEnergy;

	sys.setActiveRotamers(_stateVec);
	PDBWriter designWriter;
	designWriter.open(repackDir + "/" + MslTools::intToString(_seqNumber) + ".pdb");
	designWriter.write(sys.getAtomPointers(), true, false, true);
	designWriter.close();
	_writer.write(sys.getAtomPointers(), true, false, true);
	saveEnergyDifference(_opt, _sequenceEnergyMap, _sequence);

	//Write CRD
	CRDWriter writer;
	writer.open(repackDir + "/" + MslTools::intToString(_seqNumber) + ".crd");
	writer.write(sys.getAtomPointers(), false);
	writer.close();
}

void computeMonomerEnergyNoIMM1(Options& _opt, map<string,map<string,double>> &_sequenceEnergyMap, string &_seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {

	string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, _opt.thread);
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

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

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

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);

	/*****************************************************************************
	 *                 === GREEDY TO OPTIMIZE ROTAMERS ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, "SL95.00");
	CSBMono.updateNonBonded(10,12,50);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	/*****************************************************************************
	 *            === SET SYSTEM TO BEST SPM ROTAMERS AND OUTPUT ===
	 ******************************************************************************/
	repackSideChains(monoSpm, 10);
	vector<uint> stateVec = monoSpm.getMinStates()[0];

	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getStateEnergy(stateVec)*2;

	// Add energy to sequence energy map
	map<string,double> &energyMap = _sequenceEnergyMap[_seq];
	outputEnergiesByTerm(monoSpm, stateVec, energyMap, _opt.energyTermList, "MonomerNoIMM1", false);
	_sequenceEnergyMap[_seq]["MonomerNoIMM1"] = monomerEnergy;
	vector<double> selfVec = calcBaselineEnergies(monoSys, _opt.backboneLength);
	vector<double> pairVec = calcPairBaselineEnergies(monoSys, _opt.backboneLength);
	//double self = sumEnergyVector(selfVec);
	//double pair = sumEnergyVector(pairVec);
	_sout << "Monomer Energy No IMM1: " << monomerEnergy << endl;
	//cout << "Self Energy:  " << self << endl;
	//cout << "Pair Energy:  " << pair << endl;
	//cout << "Total Energy: " << self+pair << endl;
}

void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {

	//string polySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, _opt.thread);//fixed monomer calculation issue on 05_12_2021
	string polySeq = generateMonomerPolymerSequenceFromSequence(_seq, _opt.thread);
	PolymerSequence PS(polySeq);

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

	CSBMono.setBuildNonBondedInteractions(false);
	if (!CSBMono.buildSystem(PS)){
		cerr << "Unable to build system from " << polySeq << endl;
	}

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

	SystemRotamerLoader monoRot(monoSys, _opt.rotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hbondFile);
	monohb.buildInteractions(50);

	/*****************************************************************************
	 *              === DELETE TERMINAL HYDROGEN BOND INTERACTIONS ===
	 ******************************************************************************/
	deleteTerminalHydrogenBondInteractions(monoSys,_opt);

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
	monoEset->setTermActive("CHARMM_VDW", true);
	monoEset->setTermActive("SCWRL4_HBOND", true);

	monoEset->setWeight("CHARMM_VDW", 1);
	monoEset->setWeight("SCWRL4_HBOND", 1);
	monoEset->setWeight("CHARMM_IMM1REF", 1);
	monoEset->setWeight("CHARMM_IMM1", 1);

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
		cout << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	/*****************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	loadRotamers(monoSys, monoRot, "SL95.00");
	CSBMono.updateNonBonded(10,12,50);

	// Optimize Initial Starting Position (using Baseline to get back to original result)
	SelfPairManager monoSpm;
	monoSpm.seed(_RNG.getSeed());
	monoSpm.setSystem(&monoSys);
	monoSpm.setVerbose(false);
	monoSpm.getMinStates()[0];
	monoSpm.updateWeights();
	monoSpm.setOnTheFly(true);
	monoSpm.saveEnergiesByTerm(true);
	monoSpm.calculateEnergies();

	repackSideChains(monoSpm, _opt.greedyCycles);
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);

	monoSys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");

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
	monoSys.calcEnergy();

	// move center of mass to origin
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
	helicalAxis.saveAltCoor("BestMonomerAxis");
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
				helicalAxis.saveAltCoor("BestMonomerAxis");
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
			helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			//_fout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}
		counter++;
	}

	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	//Calculate Monomer energy for output
	monoSys.applySavedCoor("savedBestMonomer");
	helicalAxis.applySavedCoor("BestMonomerAxis");
        monoSpm.runGreedyOptimizer(_opt.greedyCycles);
	vector<uint> stateVec = monoSpm.getMinStates()[0];
	monoSys.setActiveRotamers(stateVec);
	double monomerEnergy = monoSpm.getStateEnergy(stateVec)*2;
	_sout << "Monomer Energy w/ IMM1: " << monomerEnergy << endl;
	cout << "Monomer Energy w/ IMM1: " << monomerEnergy << endl;

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
	helicalAxis.clearSavedCoor("BestMonomerAxis");
	helicalAxis.clearSavedCoor("bestZ");
}

void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs, RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	cout << "Calculating monomer energies..." << endl;
	_sout << "Calculating monomer energies..." << endl;
	for (uint m=0; m<_seqs.size(); m++){
		cout << "Sequence " << m+1 << ": " << _seqs[m] << endl;
		_sout << "Sequence " << m+1 << ": " << _seqs[m] << endl;
		//TODO: make sure this works
		if (_opt.useIMM1){
			//TODO: how to make this work for heterodimer? another vector? a pair of vector strings? compute for both sequences...
			computeMonomerEnergyNoIMM1(_opt, _sequenceEnergyMap, _seqs[m], _RNG, _sout, _err);
		} else {
			computeMonomerEnergyIMM1(_opt, _trans, _sequenceEnergyMap, _seqs[m], _RNG, _sout, _err);
		}
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	cout << "End monomer calculations: " << diffTimeMono << "s" << endl << endl;
	_sout << "End monomer calculations: " << diffTimeMono << "s" << endl << endl;
}