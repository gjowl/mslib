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
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

/***********************************
 * load rotamers 
 ***********************************/
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling){
	//Repack side chains based on sasa scores
	// get backbone length
	int backboneLength = _sys.getChain("A").positionSize();
	for (uint i = 0; i < _rotamerSampling.size()/2; i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+backboneLength);
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _opt.sasaRepackLevel[_rotamerSampling[i]];
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
			string posRot = _opt.sasaRepackLevel[_rotamerSampling[i]];
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

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<int> &_rotamerSampling){
	// if using the SASA to identify the interface, then load the rotamers by the SASA burial
	if (_opt.useSasaBurial){
		loadRotamersBySASABurial(_sys, _sysRot, _opt, _rotamerSampling);
	} else {
		loadRotamers(_sys, _sysRot, _opt.SL);
	}
}

/***********************************
 *geometry
 ***********************************/
void getGeometry(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_densities, ofstream &_out){
	// Setup file reader
	Reader reader(_opt.geometryDensityFile);
	reader.open();
	if(!(reader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.geometryDensityFile << endl;
		exit(0);
	}
	vector<string> lines = reader.getAllLines();

	// Extract the geometries from a random line of the geometry file
	int geometryLine = _RNG.getRandomInt(1,lines.size()-1);
	vector<string> tokens = MslTools::tokenize(lines[geometryLine], "\t");//xShift, crossingAngle, axialRotation, zShift, angleDistDensity, axialRotationDensity, zShiftDensity
	_opt.xShift = MslTools::toDouble(tokens[0]);
	_opt.crossingAngle = MslTools::toDouble(tokens[1]);
	_opt.axialRotation = MslTools::toDouble(tokens[2]);
	_opt.zShift = MslTools::toDouble(tokens[3]);
	double angleDistDensity = MslTools::toDouble(tokens[4]);
	double axialRotationDensity = MslTools::toDouble(tokens[5]);
	double zShiftDensity = MslTools::toDouble(tokens[6]);
	_densities.push_back(angleDistDensity);
	_densities.push_back(axialRotationDensity);
	_densities.push_back(zShiftDensity);

	// Output to summary file
	_out << "***STARTING GEOMETRY:***" << endl;
	_out << "xShift:        " << _opt.xShift << "\tDensity: " << angleDistDensity << endl;
	_out << "crossingAngle: " << _opt.crossingAngle << "\tDensity: " << angleDistDensity << endl;
	_out << "axialRotation: " << _opt.axialRotation << "\tDensity: " << axialRotationDensity << endl;
	_out << "zShift:        " << _opt.zShift << "\tDensity: " << zShiftDensity << endl << endl;
}

void getAxialRotAndZShift(Options &_opt, RandomNumberGenerator &_RNG, vector<double> &_densities, ofstream &_out){
	// Setup file reader
	Reader reader(_opt.geometryDensityFile);
	reader.open();
	if(!(reader.is_open())){
		cerr << "WARNING: Unable to open " << _opt.geometryDensityFile << endl;
		exit(0);
	}
	vector<string> lines = reader.getAllLines();

	// Extract the geometries from a random line of the geometry file
	int geometryLine = _RNG.getRandomInt(1,lines.size()-1);
	vector<string> tokens = MslTools::tokenize(lines[geometryLine], "\t");//xShift, crossingAngle, axialRotation, zShift, angleDistDensity, axialRotationDensity, zShiftDensity
	_opt.axialRotation = MslTools::toDouble(tokens[2]);
	_opt.zShift = MslTools::toDouble(tokens[3]);
	double angleDistDensity = MslTools::toDouble(tokens[4]);
	double axialRotationDensity = MslTools::toDouble(tokens[5]);
	double zShiftDensity = MslTools::toDouble(tokens[6]);
	_densities.push_back(0);//temp 0 for now; need to determine density of given geometry later
	_densities.push_back(axialRotationDensity);
	_densities.push_back(zShiftDensity);

	// Output to summary file
	_out << "***STARTING GEOMETRY:***" << endl;
	_out << "xShift:        " << _opt.xShift << "\tDensity: " << angleDistDensity << endl;
	_out << "crossingAngle: " << _opt.crossingAngle << "\tDensity: " << angleDistDensity << endl;
	_out << "axialRotation: " << _opt.axialRotation << "\tDensity: " << axialRotationDensity << endl;
	_out << "zShift:        " << _opt.zShift << "\tDensity: " << zShiftDensity << endl << endl;
}

/***************************************
 *define interface and rotamer sampling
 ***************************************/
//Identify which positions are found at the identified interface
//Example: Sequence:  LLLLIGLLIGLLIGLLLL
//         Interface: 000011001100110000
// Positions at interface are 1 and non-interfacial are 0
vector<int> getLinked(vector<int> _rotamerSampling, int _backboneLength, int _interfaceLevel, int _highestRotamerLevel){
	vector<int> linkedPositions;
	for (uint i=0; i<_backboneLength; i++){
		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
			linkedPositions.push_back(i);
		}
	}
	return linkedPositions;
}

// Convert positions to string for setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions) which uses "A,19" "B,19" format!
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<int> &_interfacialPositions, int _backboneLength){
	vector<vector<string>> stringPositions;
	//TODO: if you want to make it possible to do sequences that are not linked, need to change this function! 
	for (uint k=0; k<_backboneLength; k++){
		// search through the interface positions vector to get the interface
		if (_interfacialPositions[k] == 1){
			vector<string> tempPos;

			// get the positions at the backbone position k
			Position &posA = _sys.getPosition(k);
			Position &posB = _sys.getPosition(k+_backboneLength);

			// convert position to string
			string A = posA.toString();
			string B = posB.toString();

			// find the space delimiter	
			string delimiter = " ";
			size_t p = 0;
			p = A.find(delimiter);

			// add the string without anything after the delimiter to the tempPos vector
			tempPos.push_back(A.substr(0, p));
			tempPos.push_back(B.substr(0, p));
			
			// add the tempPos vector to the stringPositions vector
			stringPositions.push_back(tempPos);
		}
	}
	return stringPositions;
}

//makes the best state applicable for unlinked positions by duplicating the rotamer at each interfacial position on the opposite chain
void unlinkBestState(Options &_opt, vector<uint> &_bestState, vector<int> _rotamerSampling, int _backboneLength){
	vector<int> linkedPositions = getLinked(_rotamerSampling, _backboneLength, _opt.interfaceLevel, _opt.sasaRepackLevel.size()-1);

	for (uint i=_backboneLength; i<_backboneLength*2; i++){
		vector<int>::iterator itr;
		itr = find(linkedPositions.begin(), linkedPositions.end(), i-_backboneLength);
		if (itr != linkedPositions.end()){
			vector<uint>::iterator itPos;
			itPos = _bestState.begin()+i;
			_bestState.insert(itPos, _bestState[i-_backboneLength]);
		}
	}
}

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

// converts a polymer sequence to a string for one chain (for homodimer sequences)
string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getCurrentIdentity().getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}

// output energies by term into a referenced energyMap
void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1){
	if (_includeIMM1 == false){//No IMM1 Energy (for the Monte Carlos, both dimer and monomer)
		for (uint i=0; i<_energyTermList.size(); i++){
			string energyTerm = _energyTermList[i]; //CHARMM_ and SCWRL4_ terms
			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;//Removes the CHARMM_ and SCWRL4_ before energyTerm names
			if (energyTerm.find("IMM1") != string::npos){
				continue;
			} else {
				if (_energyDescriptor.find("Monomer") != string::npos){
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
				}
			}
		}
		if (_energyDescriptor.find("Monomer") != string::npos){
			//skip if monomer; could add calc baseline here at some point
		} else {
			_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
			_energyMap["DimerSelfBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE");
			_energyMap["DimerPairBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
		}
	} else if (_includeIMM1 == true){//IMM1 Energies
		for (uint i=0; i<_energyTermList.size(); i++){
			string energyTerm = _energyTermList[i];
			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;
			if (_energyDescriptor.find("Monomer") != string::npos){
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF"))*2;
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
				}
			} else {
				if (energyTerm.find("IMM1") != string::npos){
					_energyMap["IMM1Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF");
				} else {
					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
				}
			}
		}
		//_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
	}
}

// get the string of an interface sequence in 00010001 format where 1 is an interface residue and 0 is a non-interface residue
string getInterfaceSequence(Options &_opt, string _interface, string _sequence){
	string interfaceSequence = "";
	for(string::iterator it = _interface.begin(); it != _interface.end(); it++) {
		stringstream ss;
		ss << *it;
		int pos = it-_interface.begin();
		int tmp = MslTools::toInt(ss.str());
		if (tmp > _opt.interfaceLevel-1){//interfaceLevel counts from 1 but rotamerLevel coutns from 0
			interfaceSequence = interfaceSequence + "-";
		} else {
			interfaceSequence = interfaceSequence + _sequence[pos];
		}
	}
	return interfaceSequence;
}

/***********************************
 *  calculate residue burial
 ***********************************/
// function from Samson that calculates the residue burial of a system by comparing it to the typical burial of an amino acid
std::vector<pair <int, double> > calculateResidueBurial (System &_sys) {
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
	//b.setProbeRadius(3.0);
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

//Calculate Residue Burial and output a PDB that highlights the interface
std::vector<pair <int, double> > calculateResidueBurial (Options &_opt, System &_startGeom, string _seq) {
	// polymer sequences have: chain, starting position of chain residue, three letter AA code
	string polySeq = generatePolymerSequence(_opt.backboneAA, _opt.backboneLength, _opt.thread);
	PolymerSequence PS(polySeq);

	// Declare system for dimer
	System sys;
	CharmmSystemBuilder CSB(sys,_opt.topFile,_opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setBuildNonBondedInteractions(false);
	//CSB.setBuildNoTerms();

	if(!CSB.buildSystem(PS)) {
		cout << "Unable to build system from " << PS << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

	// Set weights
	Eset->setWeight("CHARMM_VDW", 1);
	Eset->setWeight("SCWRL4_HBOND", 1);

	CSB.updateNonBonded(10,12,50);

	// initialize the object for loading rotamers into our _system
	SystemRotamerLoader sysRot(sys, _opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, _opt.hbondFile);
	hb.buildInteractions(50);//when this is here, the HB weight is correct

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(_startGeom.getAtomPointers(),false);
	sys.buildAllAtoms();

	string monoPolySeq = convertToPolymerSequenceNeutralPatchMonomer(_seq, 1);
	PolymerSequence MPS(monoPolySeq);

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
	if (!CSBMono.buildSystem(MPS)){
		cerr << "Unable to build system from " << monoPolySeq << endl;
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

	std::vector<pair <int, double> > residueBurial;
	SasaCalculator dimerSasa(sys.getAtomPointers());
	SasaCalculator monoSasa(monoSys.getAtomPointers());
	dimerSasa.calcSasa();
	monoSasa.calcSasa();
	dimerSasa.setTempFactorWithSasa(true);

	for (uint i = 0; i < monoSys.positionSize(); i++) {//Changed this to account for linked positions in the dimer; gives each AA same number of rotamers as correspnding chain
		string posIdMonomer = monoSys.getPosition(i).getPositionId();
		string posIdDimer = sys.getPosition(i).getPositionId();
		double resiSasaMonomer = monoSasa.getResidueSasa(posIdMonomer);
		double resiSasaDimer = dimerSasa.getResidueSasa(posIdDimer);
		double burial = resiSasaDimer/resiSasaMonomer;
		residueBurial.push_back(pair<int,double>(i, burial));

		//set sasa for each residue in the b-factor
		AtomSelection selA(sys.getPosition(i).getAtomPointers());
		AtomSelection selB(sys.getPosition(i+_opt.backboneLength).getAtomPointers());
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
	// sort in descending order of burial
	sort(residueBurial.begin(), residueBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});
	return residueBurial;
}

// get a vector of all interfacial positions, including the ends
vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition, int _backboneLength){
	vector<uint> variableInterfacePositions;
	for (uint k=0; k<_backboneLength; k++){//TODO: make this not hardcoded to skip RAS
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

// get a vector of interface positions that doesn't include the ends of the sequence
vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition, int _backboneLength){
	vector<uint> variableInterfacePositions;
	for (uint k=3; k<_backboneLength-5; k++){//TODO: make this not hardcoded to skip RAS
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

/***********************************
 *output file functions
 ***********************************/
void setupDesignDirectory(Options &_opt){
	_opt.pdbOutputDir = string(get_current_dir_name()) + "/design_" + _opt.runNumber;
	//_opt.pdbOutputDir = "/exports/home/gloiseau/mslib/trunk_AS/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _opt.pdbOutputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
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

// build in the baseline energies to the system EnergySet
void buildBaselines(System &_sys, Options &_opt){
		map<string, double> selfMap = readSingleParameters(_opt.selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_opt.pairEnergyFile);
		buildSelfInteractions(_sys, selfMap);
		buildPairInteractions(_sys, pairMap);
}
// reads in the self energy baseline file
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
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 2){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) Energy(double)";
			continue;
		}
		selfEnergies[MslTools::toUpper(tokens[0])] = MslTools::toDouble(tokens[1]);
		//cout << tokens[0] << " " << tokens[1] << " " << tokens[2] << " = " << tokens[3] << endl;
	}

	selfReader.close();
	return selfEnergies;
}

// reads in the pair energy baseline file
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
		if(tokens.size() < 1){
			continue;
		}
		if(tokens.size() != 4){
			cerr << "WARNING: Line\"" << lines[i] << "\" is not in FORMAT: ResName(string) ResName(string) Distance(uint) Energy(double)";
			continue;
		}
		if (tokens[0].compare(tokens[1]) == 0){//Added in on 03-18-2021: apparently the code that I was using to flip the AAs in buildPairInteractions isn't good, but this adds the flips to the map which works better and is cleaner
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		} else {
			pairEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
			pairEnergies[MslTools::toUpper(tokens[1])][MslTools::toUpper(tokens[0])][MslTools::toInt(tokens[2])] = MslTools::toDouble(tokens[3]);
		}
	}

	pairReader.close();
	return pairEnergies;
}

// builds the self energy interactions into energy set
void buildSelfInteractions(System &_sys, map<string, double> &_selfMap){
	EnergySet* ESet = _sys.getEnergySet();

	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j<(*p)->identitySize();j++){
				Residue &res = (*p)->getIdentity(j);
				string baseId = res.getResidueName();
				if (p-positions.begin() < 3 || p-positions.begin() > positions.size()-5){//On 03_18_2021 I found this error; position.size() is weird, so need to use 5 instead of 4; on 11_20_2021 saw that a lot of clashing occurs at hte 4th position, so changed this to only use 3
					baseId = baseId.append("-OUT");
				}
				try{
					double ener = _selfMap.at(baseId);
					Atom *a = &res.getAtom("CA");
					ESet->addInteraction(new BaselineInteraction(*a,ener));
				}
				catch (const out_of_range& e){
					continue;
				}
			}
		}
	}
}

// builds the pair energy interactions into energy set
void buildPairInteractions(System &_sys, map<string,map<string,map<uint,double>>>& _pairMap){
	EnergySet* ESet = _sys.getEnergySet();
	for(uint i = 0; i < _sys.chainSize(); i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++){
			for (uint j=0; j < (*p)->identitySize(); j++){
				Residue &res1 = (*p)->getIdentity(j);
				string baseId1 = res1.getResidueName();
				if (p-positions.begin() < 1){
					baseId1 = baseId1.append("-ACE");
				}
				//Changed this on 11_24_2021 for the designFiles/2021_11_22_IMM1Self and Pair baselines
				//if (p-positions.begin() > positions.size()-5){//
				//	baseId1 = baseId1.append("-CT2");
				//}
				//cout << "Identity " << j << ": " << baseId1 << endl;
				for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					//cout << "Position 2: " << p2-positions.begin() << endl;
					if (d <= 10){
						//cout << "Distance: " << d << endl;
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							if (p2-positions.begin() < 3){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-2){
								baseId2 = baseId2.append("-CT2");
							}
							try{
								map<string,map<uint,double>> AA1 = _pairMap.at(baseId1);
								map<uint,double> AA2 = AA1.at(baseId2);
								double ener = AA2.at(d);
								Atom *a = &res1.getAtom("CA");
								Atom *b = &res2.getAtom("CA");
								ESet->addInteraction(new BaselinePairInteraction(*a,*b,-1*ener));//I forgot that this needs to be the opposite sign to actually counteract the energies of vdW and hydrogen bonding; switched it here but should eventually just switch in my baseline parameter file
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
}

/***********************************
 *calculate energies
 ***********************************/
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
	int firstPos = 0;
    int lastPos = monoSys.positionSize();
    deleteTerminalBondInteractions(monoSys,_opt,firstPos,lastPos);

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

void computeMonomerEnergyIMM1(Options& _opt, Transforms & _trans, map<string,map<string,double>> &_sequenceEnergyMap, string _seq,
 RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err) {

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
	int firstPos = 0;
    int lastPos = monoSys.positionSize();
    deleteTerminalBondInteractions(monoSys,_opt,firstPos,lastPos);

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
	moveZCenterOfCAMassToOrigin(chainA, helicalAxis.getAtomPointers(), _trans);
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
			helicalAxis.saveAltCoor("BestMonomerAxis");
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
	double dimerEnergy = _sequenceEnergyMap[_seq]["Dimer"];
	cout << "Dimer Energy: " << _seq << ": " << dimerEnergy << endl;
	double totalEnergy = dimerEnergy-monomerEnergy;
	_sout << "-Dimer - Monomer = " << dimerEnergy << " - " << monomerEnergy << " = " << totalEnergy << endl;
	cout << "-Dimer - Monomer = " << dimerEnergy << " - " << monomerEnergy << " = " << totalEnergy << endl;
	_sequenceEnergyMap[_seq]["Total"] = totalEnergy;

	// Clear saved coordinates
	monoSys.clearSavedCoor("savedBestMonomer");
	monoSys.clearSavedCoor("bestZ");
	helicalAxis.clearSavedCoor("BestMonomerAxis");
	helicalAxis.clearSavedCoor("bestZ");
}

void computeMonomerEnergies(Options &_opt, Transforms &_trans, map<string, map<string,double>> &_sequenceEnergyMap, vector<string> &_seqs,
 RandomNumberGenerator &_RNG, ofstream &_sout, ofstream &_err){
	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	cout << "Calculating monomer energies..." << endl;
	_sout << "Calculating monomer energies..." << endl;
	vector<thread> threads;
	for (auto &seq : _sequenceEnergyMap) {
		string sequence = seq.first;
		threads.push_back(thread(computeMonomerEnergyIMM1, ref(_opt), ref(_trans), ref(_sequenceEnergyMap), sequence, ref(_RNG), ref(_sout), ref(_err)));
	}
	for (auto &t: threads){
		t.join();
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	cout << "End monomer calculations: " << diffTimeMono/60 << "min" << endl << endl;
	_sout << "End monomer calculations: " << diffTimeMono/60 << "min" << endl << endl;
}

void deleteTerminalBondInteractions(System &_sys, Options &_opt, int _firstResiNum, int _lastResiNum){
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

void getSasaForStartingSequence(System &_sys, string _sequence, vector<uint> _state, map<string, map<string,double>> &_sequenceEnergyMap){
	_sys.setActiveRotamers(_state);

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	_sequenceEnergyMap[_sequence]["DimerSasa"] = dimerSasa;
}

/***********************************
* sequence search functions
 ***********************************/
vector<uint> runSCMFToGetStartingSequence(System &_sys, Options &_opt, RandomNumberGenerator &_RNG, string _rotamerSamplingString,
 string _variablePositionString, vector<string> _seqs, vector<uint> _interfacialPositions, map<string, map<string,double>> &_sequenceEnergyMap, 
 map<string, double> _sequenceEntropyMap, ofstream &_out){
	// Setup time variables
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);

	// SelfPairManager setup
	SelfPairManager spm;
	spm.seed(_RNG.getSeed());
	spm.setSystem(&_sys);
	spm.setVerbose(false);
	spm.setRunDEE(_opt.runDEESingles, _opt.runDEEPairs);
	spm.setOnTheFly(false);
	spm.setMCOptions(1000, 0.5, 5000, 3, 10, 1000, 0.01);//changed to sigmoid and added up to 5000
	spm.saveEnergiesByTerm(true); //added back in on 09_21_2021 to get the vdw and hbond energies
	spm.calculateEnergies();

	//Setup running SCMF or UnbiasedMC
	if (_opt.runSCMF == true){
		cout << "Running Self Consistent Mean Field" << endl;
		_out << "Running Self Consistent Mean Field" << endl;
		spm.setRunSCMF(true);
		spm.setRunSCMFBiasedMC(true);
		spm.setRunUnbiasedMC(false);
	} else {
		cout << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		_out << "runSCMF is false; Run Unbiased Monte Carlo " << endl;
		spm.setRunSCMF(false);
		spm.setRunSCMFBiasedMC(false);
		spm.setRunUnbiasedMC(true);
	}

	// run and find a sequence using the chosen parameters (MCOptions, SCMF, DEE, etc.)
	spm.runOptimizer();
	time(&endTime);
	diffTime = difftime (endTime, startTime);

	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = spm.getBestSCMFBiasedMCState();
	_sys.setActiveRotamers(bestState);
	string startSequence = convertPolymerSeqToOneLetterSeq(_sys.getChain("A")); //used for outputting starting sequence
	string interfaceSeq = getInterfaceSequence(_opt, _rotamerSamplingString, startSequence);

	// output spm run optimizer information
	spmRunOptimizerOutput(spm, _sys, interfaceSeq, _variablePositionString, diffTime, _out);
	
	//Add energies for initial sequences into the sequenceEnergyMap
	_seqs.insert(_seqs.begin(), startSequence);
	pair<string,vector<uint>> startSequenceStatePair = make_pair(startSequence, bestState);
	getEnergiesForStartingSequence(_opt, spm, startSequence, bestState, _interfacialPositions, _sequenceEnergyMap, _sequenceEntropyMap);
	getSasaForStartingSequence(_sys, startSequence, bestState, _sequenceEnergyMap);

	return bestState;
}

void getEnergiesForStartingSequence(Options &_opt, SelfPairManager &_spm, string _startSequence,
vector<uint> &_stateVector, vector<uint> _interfacialPositions, map<string, map<string, double>> &_sequenceEnergyMap, map<string, double> &_entropyMap){
	if (_opt.useBaseline){
		double baseline = _spm.getStateEnergy(_stateVector, "BASELINE")+_spm.getStateEnergy(_stateVector, "BASELINE_PAIR");
		_sequenceEnergyMap[_startSequence]["Dimer"] = _spm.getStateEnergy(_stateVector)-baseline;
		_sequenceEnergyMap[_startSequence]["Baseline"] = baseline;
	} else {
		_sequenceEnergyMap[_startSequence]["Dimer"] = _spm.getStateEnergy(_stateVector);
	}
	_sequenceEnergyMap[_startSequence]["VDWDimer"] = _spm.getStateEnergy(_stateVector, "CHARMM_VDW");
	_sequenceEnergyMap[_startSequence]["HBONDDimer"] = _spm.getStateEnergy(_stateVector, "SCWRL4_HBOND");
	_sequenceEnergyMap[_startSequence]["IMM1Dimer"] = _spm.getStateEnergy(_stateVector, "CHARMM_IMM1REF")+_spm.getStateEnergy(_stateVector, "CHARMM_IMM1");

	map<string,int> seqCountMap;
	double numberOfPermutations;
	// Found mistake on 2022-8-18: the sequence entropies for this sequence are wrong; use interfaceAASequenceEntropy instead
	interfaceAASequenceEntropySetup(_startSequence, seqCountMap, numberOfPermutations, _interfacialPositions);
	double SEProb = calculateSequenceProbability(seqCountMap, _entropyMap, numberOfPermutations);

	_sequenceEnergyMap[_startSequence]["SequenceProbability"] = SEProb;
}


void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, string _variablePosString, double _spmTime, ofstream &_out){
	// output this information about the SelfPairManager run below
	// vector for the initial SCMF state
	vector<unsigned int> initialState = _spm.getSCMFstate();
	// vector for the SCMF state after the biased monte carlo
	vector<unsigned int> bestState = _spm.getBestSCMFBiasedMCState();
	
	// checks to see if the sequence changed between the initialState and the bestState
	_sys.setActiveRotamers(initialState);
	string initialSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));
	_sys.setActiveRotamers(bestState);
	string SCMFBestSeq = convertPolymerSeqToOneLetterSeq(_sys.getChain("A"));

	// energy terms from the startingState
	double bestEnergy = _spm.getStateEnergy(bestState);
	double hbondEnergy = _spm.getStateEnergy(bestState, "SCWRL4_HBOND");
	double vdwEnergy = _spm.getStateEnergy(bestState, "CHARMM_VDW");

	// outputs	
	_out << "Initial Sequence:   " << initialSeq << endl;
	_out << "Best Sequence:      " << SCMFBestSeq << endl;
	_out << "Interface Sequence: " << _interfaceSeq << endl;
	_out << "Interface:          " << _variablePosString << endl;
	_out << "Total Energy:       " << bestEnergy << endl;
	_out << "VDW:                " << vdwEnergy << endl;
	_out << "HBOND:              " << hbondEnergy << endl;
	_out << endl << "End SelfPairManager Optimization: " << _spmTime << "s" << endl;
}

/***********************************
 *sequence entropy functions
 ***********************************/
// gets the sequence entropy value for the current sequence and adds it to the energy
void calculateInterfaceSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq,
 map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_prevEntropy,
 double &_currEntropy, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal, vector<uint> _interfacePositionsList){

	_prevSEProb = getInterfaceSequenceEntropyProbability(_opt, _prevSeq, _entropyMap, _interfacePositionsList);
	_currSEProb = getInterfaceSequenceEntropyProbability(_opt, _currSeq, _entropyMap, _interfacePositionsList);

	//Calculate the probability of each sequence compared to the other
	//On 10_26_21: after changing to the 2021_10_25_seqEntropies.txt file with more accurate membrane composition from Liu et al. 2002, I started seeing a lot of GxxxG sequences. I also checked my runs from before this date, and there were qute a bit of GxxxG sequences as well. Uncomment the below if you would like to try running using the internal probabilities, which should give a bit more flexibiilty to choosing AAs with lower likelilhoods of being in the membrane
	double totSEProb = _prevSEProb+_currSEProb;
	double prevSeqProp = _prevSEProb/totSEProb;
	double currSeqProp = _currSEProb/totSEProb;

	//Convert the probability of each sequence to an entropy term that can be added to the original energy to directly compare two sequences
	_prevEntropy = -log(prevSeqProp)*0.592*_opt.weight_seqEntropy;//Multiplies the energy (log proportion times RT in kcal/mol) by the sequence entropy weight (weight of 1 gives entropy same weight as other terms)
	_currEntropy = -log(currSeqProp)*0.592*_opt.weight_seqEntropy;

	//Calculate energy total for best sequence vs current sequence
	//The below includes the baseline energy, which is an estimate of monomer energy
	_bestEnergyTotal = _bestEnergy+_prevEntropy;
	_currEnergyTotal = _currEnergy+_currEntropy;
}

map<string,int> getAACountMap(vector<string> _seq){
	map<string,int> AAcounts;
	for (uint i=0; i<_seq.size(); i++){
		try{
			if (AAcounts.count(_seq[i]) > 0){
				AAcounts.at(_seq[i])++;
			} else {
				AAcounts[_seq[i]] = 1;
			}
		}
		catch(const out_of_range& e){
			continue;
		}
	}
	return AAcounts;
}

double calcNumberOfPermutations(map<string,int> _seqAACounts, int _seqLength){
	//This function calculates number of permutations using following equation: n!/(n-r)! where n is number of positions and r is number of AAs
	double numPermutation = 1;
	double permutationDenominator = 1;
	for(uint i=_seqLength; i>1; i--){
		numPermutation = numPermutation*i;
	}
	map<string,int>::iterator itr;
	for(itr = _seqAACounts.begin(); itr != _seqAACounts.end(); itr++){
		for (uint j=itr->second; j>1; j--){
			permutationDenominator = permutationDenominator*j;
		}
	}
	numPermutation = numPermutation/permutationDenominator;

	return numPermutation;
}

void interfaceAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<uint> _interfacialPositionsList){
	//Get residue name for each interfacial identity
	vector<string> seqVector;
	int numInterfacials = _interfacialPositionsList.size()-1;
	for (uint i=0; i<numInterfacials; i++){
		//Position &pos = _sys.getPosition(_interfacialPositionsList[i]);
		//Residue &resi = pos.getCurrentIdentity();
		//string id = resi.getResidueName();
		stringstream tmp;
		tmp << _seq[_interfacialPositionsList[i]];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
		//cout << _interfacialPositionsList[i] << " : " << resName << endl;
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, numInterfacials);
}

double getInterfaceSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap, vector<uint> _interfacialPositionsList){
	map<string,int> AACountMap;
	double numberOfPermutations;
	interfaceAASequenceEntropySetup(_sequence, AACountMap, numberOfPermutations, _interfacialPositionsList);
	double seqProb = calculateSequenceProbability(AACountMap, _entropyMap, numberOfPermutations);
	return seqProb;
}

double calculateSequenceProbability(map<string,int> &_seqCountMap, map<string,double> &_entropyMap, double _numberOfPermutations){
	//Find AA in count map, get counts of AA, then calculate the sequence probability by multiplying each AAs membrane probability contribution
	double seqProb = 1;
	map<string,int>::iterator it;
	for (it=_seqCountMap.begin(); it != _seqCountMap.end(); it++){
		double memProb = _entropyMap.at(it->first);
		int count = it->second;
		seqProb = seqProb*(pow(memProb, count));
	}
	//cout << "pSeq = " << seqProb << endl;
	seqProb = seqProb*_numberOfPermutations;
	//cout << "comb = " << _numberOfPermutations << endl;
	//cout << "pRes = " << seqProb << endl;
	return seqProb;
}

/****************************************
 *
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/
Options parseOptions(int _argc, char * _argv[]){
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
	 *  /exports/home/gloiseau/mslib/trunk_AS/config/seqDesign.config
	 *
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//opt.required.push_back("");
	//opt.allowed.push_back("");
	
	// Input Files
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("geometryDensityFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("backboneFile");
	opt.allowed.push_back("selfEnergyFile");
	opt.allowed.push_back("pairEnergyFile");
	opt.allowed.push_back("sequenceEntropyFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("helicalAxis");

	// sequence parameters
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	opt.allowed.push_back("interface");

	// booleans
	opt.allowed.push_back("getGeoFromPDBData");
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("useSasaBurial");
	opt.allowed.push_back("useTimeBasedSeed");
	opt.allowed.push_back("deleteTerminalBonds");
	opt.allowed.push_back("linkInterfacialPositions");
	opt.allowed.push_back("useAlaAtCTerminus");
	opt.allowed.push_back("useBaseline");
	opt.allowed.push_back("getRandomAxRotAndZShift");

	// repack parameters
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");


	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("endResNum");

	//Geometry
	opt.allowed.push_back("xShift");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("density");
	opt.allowed.push_back("negAngle");
	opt.allowed.push_back("negRot");
	opt.allowed.push_back("thread");

	//Monte Carlo variables
	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");
	opt.allowed.push_back("MCStartTemp");
	opt.allowed.push_back("MCEndTemp");
	opt.allowed.push_back("MCCurve");
	opt.allowed.push_back("MCConvergedSteps");
	opt.allowed.push_back("MCConvergedE");

	// Backbone Monte Carlo variables
	opt.allowed.push_back("backboneMCCycles");
	opt.allowed.push_back("backboneMCMaxRejects");
	opt.allowed.push_back("backboneMCStartTemp");
	opt.allowed.push_back("backboneMCEndTemp");
	opt.allowed.push_back("backboneMCCurve");
	opt.allowed.push_back("backboneConvergedSteps");
	opt.allowed.push_back("backboneConvergedE");

	// use different energy parameters
	opt.allowed.push_back("useIMM1");
	opt.allowed.push_back("useElec");
	
	//Weights
	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_solv");
	opt.allowed.push_back("weight_elec");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_seqEntropy");

	// Alternate Identities
	opt.allowed.push_back("Ids");

	// MonteCarlo Arguments
	opt.allowed.push_back("numStatesToSave");

	// SelfPairManager Arguments
	opt.allowed.push_back("runDEESingles");
	opt.allowed.push_back("runDEEPairs");
	opt.allowed.push_back("runSCMF");

	// Energy Terms to Output
	opt.allowed.push_back("energyTermList");
	opt.allowed.push_back("deleteTerminalInteractions");

	//Rotamers
	// Load Rotamers from SASA values (from sgfc)
	opt.allowed.push_back("sasaRepackLevel");
	opt.allowed.push_back("interfaceLevel");
	opt.allowed.push_back("SL");
	
	// Shift Size
	opt.allowed.push_back("deltaX");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaXLimit");
	opt.allowed.push_back("deltaCrossLimit");
	opt.allowed.push_back("deltaAxLimit");
	opt.allowed.push_back("deltaZLimit");
	opt.allowed.push_back("decreaseMoveSize");

	opt.allowed.push_back("runNumber");

	// Begin Parsing through the options
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
	opt.backboneFile = OP.getString("backboneFile");
	if (OP.fail()) {
		opt.warningMessages += "backboneFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.backboneFile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	opt.selfEnergyFile = OP.getString("selfEnergyFile");
	if (OP.fail()) {
		opt.warningMessages += "selfEnergyFile not specified, default \n";
		opt.warningFlag = true;
		opt.selfEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanSelf_par.txt";
	}
	opt.pairEnergyFile = OP.getString("pairEnergyFile");
	if (OP.fail()) {
		opt.warningMessages += "pairEnergyFile not specified, default \n";
		opt.warningFlag = true;
		opt.pairEnergyFile = "/exports/home/gloiseau/mslib/trunk_AS/DesignFiles/2020_10_07_meanPair_par.txt";
	}
	opt.sequenceEntropyFile = OP.getString("seqEntropyFile");
	if (OP.fail()) {
		opt.warningMessages += "seqEntropyFile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.sequenceEntropyFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/sequenceEntropies.txt";
	}
	opt.helicalAxis = OP.getString("helicalAxis");
	if (OP.fail()) {
		opt.warningMessages += "helicalAxis not specified\n";
		opt.warningFlag = true;
		opt.helicalAxis = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/parameterFiles/helicalAxis.pdb";
	}
	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine pdbOutputDir, default to current directory\n";
		opt.warningFlag = true;
	}
	
	// sequence parameters
	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.warningMessages += "sequence not specified, defaulting to polyleu\n";
		opt.warningFlag = true;
		opt.sequence = "";
	}
	opt.backboneAA = OP.getString("backboneAA");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneAA not specified, default to V\n";
		opt.backboneAA = "V";
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 21\n";
		opt.backboneLength = 21;
	}
	// override backboneLength if a sequence is specified
	if (opt.sequence != "") {
		opt.backboneLength = opt.sequence.length();
	}
	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(1) + "\n";
		opt.warningFlag = true;
		opt.startResNum = 1;
	}
	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.startResNum+opt.backboneLength) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.startResNum+opt.backboneLength;
	}
	opt.interface = OP.getString("interface");
	if (OP.fail()) {
		opt.warningMessages += "interface not specified\n";
		opt.warningFlag = true;
		opt.interface = "";
	}
	// if sequence is specified, define interface for sequence
	if (opt.sequence != "" && opt.interface == "") {
		opt.interface = "";
		for (uint i=0; i<opt.sequence.length(); i++) {
			// assumes polyLeu backbone
			if (i > 3 || i < opt.sequence.length()-5) {
				// if not Leu, then interface
				if (opt.sequence[i] != 'L') {
					opt.interface += "1";
				} else {
					opt.interface += "0";
				}
			} else {
				opt.interface += "0";
			}
		}
	}
	if (opt.interface != "") {
		if (opt.interface.length() != opt.backboneLength) {
			opt.errorMessages += "interface string and backbone length must be the same length\n";
			opt.errorFlag = true;
		}
	}

	// booleans
	opt.getGeoFromPDBData = OP.getBool("getGeoFromPDBData");
	if (OP.fail()) {
		opt.warningMessages += "getGeoFromPDBData not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = false;
	}
	opt.deleteTerminalBonds = OP.getBool("deleteTerminalBonds");
	if (OP.fail()) {
		opt.deleteTerminalBonds = true;
		opt.warningMessages += "deleteTerminalBonds not specified using true\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalInteractions = OP.getMultiString("deleteTerminalInteractions");
	if (OP.fail()) {
		opt.deleteTerminalInteractions.push_back("SCWRL4_HBOND");
		opt.warningMessages += "deleteTerminalInteractions not specified, defaulting to delete SCWRL4_HBOND\n";
		opt.warningFlag = true;
	}
	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.linkInterfacialPositions = OP.getBool("linkInterfacialPositions");
	if (OP.fail()) {
		opt.warningMessages += "linkInterfacialPositions not specified using true\n";
		opt.warningFlag = true;
		opt.linkInterfacialPositions = true;
	}
	opt.useTimeBasedSeed = OP.getBool("useTimeBasedSeed");
	if (OP.fail()) {
		opt.warningMessages += "useTimeBasedSeed not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useTimeBasedSeed = false;
	}
	opt.useAlaAtCTerminus = OP.getBool("useAlaAtCTerminus");
	if (OP.fail()) {
		opt.warningMessages += "useAlaAtCTerminus not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useAlaAtCTerminus = false;
	}
	opt.useBaseline = OP.getBool("useBaseline");
	if (OP.fail()) {
		opt.warningMessages += "useBaseline not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useBaseline = false;
	}
	opt.getRandomAxRotAndZShift = OP.getBool("getRandomAxRotAndZShift");
	if (OP.fail()) {
		opt.warningMessages += "getRandomAxRotAndZShift not specified, defaulting to false";
		opt.warningFlag = true;
		opt.getRandomAxRotAndZShift = false;
	}
	//use different energy parameters
	opt.useIMM1 = OP.getBool("useIMM1");
	if (OP.fail()) {
		opt.warningMessages += "useIMM1 not specified, defaulting to true\n";
		opt.warningFlag = true;
		opt.useIMM1 = true;
	}
	opt.useElec = OP.getBool("useElec");
	if (OP.fail()) {
		opt.warningMessages += "useElec not specified using false\n";
		opt.warningFlag = true;
		opt.useElec = false;
	}
	
	// starting geometry
	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting getGeoFromPDBData true\n";
		opt.warningFlag = true;
		opt.getGeoFromPDBData = true;
	}
	opt.density = OP.getDouble("density");
	if (OP.fail()) {
		opt.warningMessages += "density not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.density = 0;
	}
	opt.negAngle = OP.getBool("negAngle");
	if (OP.fail()) {
		opt.warningMessages += "negAngle not specified using false\n";
		opt.warningFlag = true;
		opt.negAngle = false;
	}
	if (opt.negAngle == true){
		opt.crossingAngle = -opt.crossingAngle;
	}
	opt.negRot = OP.getBool("negRot");
	if (OP.fail()) {
		opt.warningMessages += "negRot not specified using false\n";
		opt.warningFlag = true;
		opt.negRot = false;
	}
	if (opt.negRot == true){
		opt.axialRotation = opt.axialRotation-100;
	}
	opt.thread = OP.getInt("thread");
	if (OP.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.thread = 0;
	}

	//Load Rotamers using SASA values (from sgfc)
	opt.useSasaBurial = OP.getBool("useSasaBurial");
	if (OP.fail()) {
		opt.warningMessages += "useSasaBurial not specified, default true\n";
		opt.warningFlag = true;
	}
	opt.sasaRepackLevel = OP.getMultiString("sasaRepackLevel");
	if (OP.fail()) {
		opt.warningMessages += "sasaRepacklevel not specified! Default to one level at " + opt.SL;
		opt.sasaRepackLevel.push_back(opt.SL);
	}
	opt.interfaceLevel = OP.getInt("interfaceLevel");
	if (OP.fail()) {
		opt.warningMessages += "interfaceLevel not specified using 1\n";
		opt.warningFlag = true;
		opt.interfaceLevel = 1;
	}
	//rotlevel
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL95.00\n";
		opt.SL = "SL95.00";
	} else {
		opt.SL = "SL"+opt.SL;
	}

	//Monte Carlo parameters
	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.errorMessages += "Number of MC cycles not specified!\n";
		opt.errorFlag = true;
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
	opt.MCConvergedSteps = OP.getInt("MCConvergedSteps");
	if (OP.fail()) {
		opt.warningMessages += "MCConvergedSteps not specified using 10\n";
		opt.warningFlag = true;
		opt.MCConvergedSteps = 10;
	}
	//opt.MCConvergedE = OP.getDouble("MCConvergedE");
	//if (OP.fail()) {
	//	opt.warningMessages += "MCConvergedE not specified using 0.01\n";
	//	opt.warningFlag = true;
	//	opt.MCConvergedE = 0.01;
	//}

	// Backbone Monte Carlo parameters
	opt.backboneMCCycles = OP.getInt("backboneMCCycles");
	if (OP.fail()) {
		opt.backboneMCCycles = 100;
		opt.warningMessages += "Number of backboneMC cycles not specified, default to 100\n";
		opt.warningFlag = true;
	}
	opt.backboneMCMaxRejects = OP.getInt("backboneMCMaxRejects");
	if (OP.fail()) {
		opt.backboneMCMaxRejects = 5;
		opt.warningMessages += "Number of backboneMC max rejects not specified, default to using 5\n";
		opt.warningFlag = true;
	}
	opt.backboneMCStartTemp = OP.getDouble("backboneMCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCStartTemp not specified using 100.0\n";
		opt.warningFlag = true;
		opt.backboneMCStartTemp = 100.0;
	}
	opt.backboneMCEndTemp = OP.getDouble("backboneMCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.backboneMCEndTemp = 0.5;
	}
	opt.backboneMCCurve = OP.getInt("backboneMCCurve");
	if (OP.fail()) {
		opt.warningMessages += "backboneMCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.backboneMCCurve = 2;
	}
	opt.backboneConvergedSteps = OP.getInt("backboneConvergedSteps");
	if (OP.fail()) {
		opt.backboneConvergedSteps = opt.backboneMCCycles/2;
		opt.warningMessages += "backboneConvergedSteps not specified using half of given cycles (" + to_string(opt.backboneConvergedSteps) + "\n";
		opt.warningFlag = true;
	}
	opt.backboneConvergedE = OP.getDouble("backboneConvergedE");
	if (OP.fail()) {
		opt.warningMessages += "backboneConvergedE not specified using 0.001\n";
		opt.warningFlag = true;
		opt.backboneConvergedE = 0.001;
	}

	// repack parameters	
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.greedyCycles = 10;
	}

	// energy weights
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
	opt.weight_seqEntropy = OP.getDouble("weight_seqEntropy");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_seqEntropy not specified, default 1.0\n";
		opt.weight_seqEntropy = 1.0;
	}
	opt.weight_elec = OP.getDouble("weight_elec");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_elec not specified, default 1.0\n";
		opt.weight_elec = 1.0;
	}

	// alternate identities
	opt.Ids = OP.getStringVector("Ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}

	//SelfPairManager Optimization Options
	opt.runDEESingles = OP.getBool("runDEESingles");
	if (OP.fail()) {
		opt.warningMessages += "runDEESingles not specified, defaulting to false";
		opt.warningFlag = true;
		opt.runDEESingles = false;
	}
	opt.runDEEPairs = OP.getBool("runDEEPairs");
	if (OP.fail()) {
		opt.warningMessages += "runDEEPairs not specified, defaulting to false";
		opt.warningFlag = true;
		opt.runDEEPairs = false;
	}
	opt.runSCMF = OP.getBool("runSCMF");
	if (OP.fail()) {
		opt.warningMessages += "runSCMF not specified, defaulting to true";
		opt.warningFlag = true;
		opt.runSCMF = true;
	}

	//Energy Terms to Output
	opt.energyTermList = OP.getStringVector("energyTermList");
	if (OP.fail()) {
		//This works, but I think if you ever want to output more terms in the future, need to add them to the terms above
		//TODO: write in an error that will tell you if the above is the case
		opt.energyTermList.push_back("CHARMM_VDW");
		opt.energyTermList.push_back("SCWRL4_HBOND");
		opt.energyTermList.push_back("CHARMM_IMM1");
		opt.energyTermList.push_back("CHARMM_IMM1REF");
	}
	
	// backbone repack variables
	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.5\n";
		opt.warningFlag = true;
		opt.deltaX = 0.5;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 3.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 3.0;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 3.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 3.0;
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
	
	// run parameters	
	opt.runNumber = OP.getString("runNumber");
	if (OP.fail()) {
		opt.warningMessages += "runNumber not specified, using 1\n";
		opt.warningFlag = true;
		opt.runNumber = MslTools::intToString(1);
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}
	opt.rerunConf = OP.getConfFile();

	return opt;
}