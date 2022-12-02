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
#include "multiCodeFunctions.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

/***********************************
 * load rotamers 
 ***********************************/
void loadRotamersBySASABurial(System &_sys, SystemRotamerLoader &_sysRot, vector<string> _repackLevels, vector<uint> _rotamerSampling){
	//Repack side chains based on sasa scores
	// get backbone length
	int backboneLength = _sys.getChain("A").positionSize();
	for (uint i = 0; i < _rotamerSampling.size(); i++) {
		Position &posA = _sys.getPosition(i);
		Position &posB = _sys.getPosition(i+backboneLength);
		if (posA.identitySize() > 1){
			for (uint j=0; j < posA.getNumberOfIdentities(); j++){
				posA.setActiveIdentity(j);
				posB.setActiveIdentity(j);
				string posRot = _repackLevels[_rotamerSampling[i]];
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
			string posRot = _repackLevels[_rotamerSampling[i]];
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

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, Options &_opt, vector<uint> &_rotamerSampling){
	// if using the SASA to identify the interface, then load the rotamers by the SASA burial
	if (_opt.useSasaBurial){
		loadRotamersBySASABurial(_sys, _sysRot, _opt.sasaRepackLevel, _rotamerSampling);
	} else {
		loadRotamers(_sys, _sysRot, _opt.SL);
	}
}


/***************************************
 *define interface and rotamer sampling
 ***************************************/
//Identify which positions are found at the identified interface
//Example: Sequence:  LLLLIGLLIGLLIGLLLL
//         Interface: 000011001100110000
// Positions at interface are 1 and non-interfacial are 0
// Convert positions to string for setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions) which uses "A,19" "B,19" format!
vector<vector<string>> convertToLinkedFormat(System &_sys, vector<uint> &_interfacePositions, int _backboneLength){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_interfacePositions.size(); k++){
		int interfacePosition = _interfacePositions[k];
		// search through the interface positions vector to get the interface
		vector<string> tempPos;

		// get the positions at the backbone position k
		Position &posA = _sys.getPosition(interfacePosition);
		Position &posB = _sys.getPosition(interfacePosition+_backboneLength);

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
	return stringPositions;
}

//makes the best state applicable for unlinked positions by duplicating the rotamer at each interfacial position on the opposite chain
vector<uint> unlinkBestState(vector<uint> _bestState, vector<uint> _interfacePositions, int _backboneLength){
	int originalSize = _bestState.size();
	for (uint i=0; i<originalSize; i++){
		_bestState.push_back(_bestState[i]);
	}
	return _bestState;
}

/***********************************
 *string output functions
 ***********************************/
// get the string of an interface sequence in 00010001 format where 1 is an interface residue and 0 is a non-interface residue
string getInterfaceSequence(int _interfaceLevelLimit, string _interface, string _sequence){
	string interfaceSequence = "";
	for(string::iterator it = _interface.begin(); it != _interface.end(); it++) {
		stringstream ss;
		ss << *it;
		int pos = it-_interface.begin();
		int tmp = MslTools::toInt(ss.str());
		if (tmp > _interfaceLevelLimit-1){//interfaceLevel counts from 1 but rotamerLevel coutns from 0
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
	string polySeq = convertToPolymerSequence(_seq, _opt.thread);
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
	// sort in descending order of burial (most buried first)
	sort(residueBurial.begin(), residueBurial.end(), [](auto &left, auto &right) {
			return left.second < right.second;
	});
	return residueBurial;
}

// get a vector of all interfacial positions, including the ends
vector<uint> getAllInterfacePositions(Options &_opt, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength){
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
vector<uint> getInterfacePositions(Options &_opt, vector<uint> &_rotamerSamplingPerPosition, int _backboneLength){
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
	_opt.outputDir = string(get_current_dir_name()) + "/design_" + _opt.runNumber;
	//_opt.outputDir = "/exports/home/gloiseau/mslib/trunk_AS/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _opt.outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}

/***********************************
 *baseline energy helper functions
 ***********************************/
//Function to calculate the non-termini self energies of a chain
vector<double> calcBaselineEnergies(System &_sys, int _thread, int _backboneLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());
	//cout << "Self Energies" << endl;
	for (uint i=_thread+3; i<_thread+_backboneLength-5; i++){
		string residue = "resi, chain A and resi ";
		string number = to_string(i);
		sel.select(residue += number);
		double resi = _sys.calcEnergy("resi");
		ener.push_back(resi);
		//cout << number << ": " << resi << endl;
	}
	sel.clearStoredSelections();
	return ener;
}

//Function to calculate the non-termini pair energies of a chain
vector<double> calcPairBaselineEnergies(System &_sys, int _thread, int _backboneLength){
	vector<double> ener;
	AtomSelection sel(_sys.getAtomPointers());

	for (uint i=_thread+3; i<_thread+_backboneLength-5; i++){
		string residue = "resi1, chain A and resi ";
		string num1 = to_string(i);
		sel.select(residue += num1);
		for (uint j=i+1; j<_thread+_backboneLength-5;j++){
			int dist = j-i;
			if (dist <= 10){
				string resi1 = "resi2, chain A and resi ";
				string num2 = to_string(j);
				sel.select(resi1 += num2);
				double pair = _sys.calcEnergy("resi1", "resi2");
				ener.push_back(pair);
			} else {
				j = _thread+_backboneLength-5;
			}
		}
	}
	sel.clearStoredSelections();
	return ener;
}

// build in the baseline energies to the system EnergySet
void buildBaselines(System &_sys, string _selfEnergyFile, string _pairEnergyFile){
		map<string, double> selfMap = readSingleParameters(_selfEnergyFile);
		map<string,map<string,map<uint,double>>> pairMap = readPairParameters(_pairEnergyFile);
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
				//cout << "Adding self energy for " << baseId << " Position: " << p-positions.begin() << endl;
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
				if (p-positions.begin() < 3){
					baseId1 = baseId1.append("-ACE");
				}
				//Changed this on 11_24_2021 for the designFiles/2021_11_22_IMM1Self and Pair baselines
				if (p-positions.begin() > positions.size()-5){//
					baseId1 = baseId1.append("-CT2");
				}
				for (vector<Position*>::iterator p2 = p+1; p2 != positions.end(); p2++){
					uint d = p2-p;
					if (d <= 10){
						//cout << "Distance: " << d << endl;
						for (uint k=0; k < (*p2)->identitySize(); k++){
							Residue &res2 = (*p2)->getIdentity(k);
							string baseId2 = res2.getResidueName();
							if (p2-positions.begin() < 3){
								baseId2 = baseId2.append("-ACE");
							}
							if (p2-positions.begin() > positions.size()-5){
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
void getSasaForStartingSequence(System &_sys, string _sequence, vector<uint> _state, map<string, map<string,double>> &_sequenceEnergyMap){
	_sys.setActiveRotamers(_state);

	//Setup SasaCalculator to calculate the monomer SASA
	SasaCalculator sasa(_sys.getAtomPointers());
	sasa.calcSasa();
	double dimerSasa = sasa.getTotalSasa();

	_sequenceEnergyMap[_sequence]["DimerSasa"] = dimerSasa;
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
	internalAASequenceEntropySetup(_startSequence, seqCountMap, numberOfPermutations, _interfacialPositions);
	double SEProb = calculateSequenceProbability(seqCountMap, _entropyMap, numberOfPermutations);

	_sequenceEnergyMap[_startSequence]["SequenceProbability"] = SEProb;
}

void spmRunOptimizerOutput(SelfPairManager &_spm, System &_sys, string _interfaceSeq, ofstream &_out){
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
	_out << "Total Energy:       " << bestEnergy << endl;
	_out << "VDW:                " << vdwEnergy << endl;
	_out << "HBOND:              " << hbondEnergy << endl;
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

void calculateInternalSequenceEntropy(Options &_opt, string _prevSeq, string _currSeq,
 map<string,double> _entropyMap, double &_prevSEProb, double &_currSEProb, double &_prevEntropy,
 double &_currEntropy, double _bestEnergy, double _currEnergy, double &_bestEnergyTotal, double &_currEnergyTotal, vector<uint> _interfacePositionsList){

	_prevSEProb = getInternalSequenceEntropyProbability(_opt, _prevSeq, _entropyMap, _interfacePositionsList);
	_currSEProb = getInternalSequenceEntropyProbability(_opt, _currSeq, _entropyMap, _interfacePositionsList);

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
	//cout << "Number Positions = " << _seqLength << endl;
	for(uint i=_seqLength; i>1; i--){
		numPermutation = numPermutation*i;
	}
	// was running the below before 2022-10-19; realized it doesn't give me the n!/(n-r)!, so not sure why I thought it was
	// as of 2022-11-3: running this again per Alessandro; it should be n!/(numAA1!*numAA2!*...)
	map<string,int>::iterator itr;
	for(itr = _seqAACounts.begin(); itr != _seqAACounts.end(); itr++){
		//cout << itr->first << ": " << itr->second << endl;
		for (uint j=itr->second; j>1; j--){
			permutationDenominator = permutationDenominator*j;
		}
	}
	//double denomFactorial = _seqLength-_seqAACounts.size();
	//cout << "Denominator Factorial = " << denomFactorial << "!" << endl;
	//cout << "Num unique AAs = " << _seqAACounts.size();
	//for (uint i=denomFactorial; i>1; i--){
	//	permutationDenominator = permutationDenominator*i;
	//}
	//cout << "Number Permutations: " << numPermutation << endl;
	//cout << "Permutation Denominator: " << permutationDenominator << endl;
	numPermutation = numPermutation/permutationDenominator;
	//cout << "Number Permutations/Denominator: " << numPermutation << endl;

	return numPermutation;
}

void interfaceAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<uint> _interfacialPositionsList){
	//Get residue name for each interfacial identity
	vector<string> seqVector;
	int numInterfacials = _interfacialPositionsList.size();
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

void internalAASequenceEntropySetup(string _seq, map<string,int> &_seqCountMap, double &_numberOfPermutations, vector<uint> _interfacialPositionsList){
	//Get residue name for each interfacial identity
	vector<string> seqVector;
	int numInterfacials = _interfacialPositionsList.size();
	for (uint i=3; i<_seq.length()-4; i++){
		//Position &pos = _sys.getPosition(_interfacialPositionsList[i]);
		//Residue &resi = pos.getCurrentIdentity();
		//string id = resi.getResidueName();
		stringstream tmp;
		tmp << _seq[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		seqVector.push_back(resName);
		//cout << _interfacialPositionsList[i] << " : " << resName << endl;
	}
	_seqCountMap = getAACountMap(seqVector);
	_numberOfPermutations = calcNumberOfPermutations(_seqCountMap, seqVector.size());
}

double getInternalSequenceEntropyProbability(Options &_opt, string _sequence, map<string,double> &_entropyMap, vector<uint> _interfacialPositionsList){
	map<string,int> AACountMap;
	double numberOfPermutations;
	internalAASequenceEntropySetup(_sequence, AACountMap, numberOfPermutations, _interfacialPositionsList);
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
	opt.allowed.push_back("outputDir");
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
	opt.allowed.push_back("numberOfSequencesToSave");

	// booleans
	opt.allowed.push_back("getGeoFromPDBData");
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("useSasaBurial");
	opt.allowed.push_back("useTimeBasedSeed");
	opt.allowed.push_back("deleteTerminalBonds");
	opt.allowed.push_back("linkInterfacialPositions");
	opt.allowed.push_back("useAlaAtTermini");
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
	opt.allowed.push_back("MCResetTemp");
	opt.allowed.push_back("MCResetCycles");

	// Backbone Monte Carlo variables
	opt.allowed.push_back("backboneMCCycles");
	opt.allowed.push_back("backboneMCMaxRejects");
	opt.allowed.push_back("backboneMCStartTemp");
	opt.allowed.push_back("backboneMCEndTemp");
	opt.allowed.push_back("backboneMCCurve");
	opt.allowed.push_back("backboneConvergedSteps");
	opt.allowed.push_back("backboneConvergedE");
	opt.allowed.push_back("backboneSearchCycles");

	// use different energy parameters
	opt.allowed.push_back("useIMM1");
	opt.allowed.push_back("useElec");
	opt.allowed.push_back("compareSasa");
	
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
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages += "Unable to determine outputDir, default to current directory\n";
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
	opt.numberOfSequencesToSave = OP.getInt("numberOfSequencesToSave");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "numberOfSequencesToSave not specified, default to 5\n";
		opt.numberOfSequencesToSave = 5;
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
	opt.useAlaAtTermini = OP.getBool("useAlaAtTermini");
	if (OP.fail()) {
		opt.warningMessages += "useAlaAtTermini not specified, defaulting to false";
		opt.warningFlag = true;
		opt.useAlaAtTermini = false;
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
	opt.compareSasa = OP.getBool("compareSasa");
	if (OP.fail()) {
		opt.warningMessages += "compareSasa not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.compareSasa = false;
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
		//opt.axialRotation = -opt.axialRotation;//for CATM geometries?
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
		opt.useSasaBurial = true;
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
		opt.warningMessages += "Number of MCResetCycles not specified, using 100\n";
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
	opt.MCResetTemp = OP.getDouble("MCResetTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCResetTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.MCResetTemp = 3649; // 50% likelihood of accepting a solution within 5kcals
	}
	opt.MCResetCycles = OP.getInt("MCResetCycles");
	if (OP.fail()) {
		opt.warningMessages += "Number of MCResetCycles not specified, using 100\n";
		opt.warningFlag = true;
		opt.MCResetCycles = 100;
	}

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
	opt.backboneSearchCycles = OP.getInt("backboneSearchCycles");
	if (OP.fail()) {
		opt.backboneSearchCycles = 5;
		opt.warningMessages += "Number of backbone search cycles not specified, default to 5\n";
		opt.warningFlag = true;
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
	// String for the alternateIds at the interface
	if (opt.xShift <= 7.5){
		opt.Ids.push_back("GLY");
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

///***********************************
// *geometry
// ***********************************/
//void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {
//	AtomSelection sel(_apV);
//	AtomPointerVector & caApV = sel.select("name CA");
//	double zShift = 0.0;
//	for(int i = 0; i < caApV.size(); i++) {
//		zShift += (caApV[i]->getCoor()).getZ();
//	}
//	zShift = -1.0 * zShift/double(caApV.size());
//
//	CartesianPoint interDistVect;
//	interDistVect.setCoor(0.0, 0.0, zShift);
//	_trans.translate(_apV, interDistVect);
//	_trans.translate(_axis, interDistVect);
//}
//
//void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
//	// Set coordinates of chain A to chain B
//	for (uint i=0; i < _apvA.size(); i++) {
//		_apvB[i]->copyAllCoor(*_apvA[i]);
//	}
//
//	// Rotation matrix for 180 degrees
//	// flips the sign on the x and y coordinates
//	Matrix m(3,3,0.0);
//	m[0][0] = -1.0;
//	m[0][1] = 0.0;
//	m[0][2] = 0.0;
//	m[1][0] = 0.0;
//	m[1][1] = -1.0;
//	m[1][2] = 0.0;
//	m[2][0] = 0.0;
//	m[2][1] = 0.0;
//	m[2][2] = 1.0;
//
//	// Rotate chain B around Z axis
//	Transforms trans;
//	trans.rotate(_apvB, m);
//}
//
//void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
//
//	//====== Z Shift (Crossing Point) ======
//	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
//	_trans.translate(_chainA, zShiftCP);
//
//	//===== Axial Rotation ======
//	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);
//
//	//====== Local Crossing Angle ======
//	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
//	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);
//
//	//====== X shift (Interhelical Distance) =======
//	CartesianPoint interDistVect;
//	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
//	_trans.translate(_chainA, interDistVect);
//	_trans.translate(_axisA, interDistVect);
//
//	c2Symmetry(_chainA, _chainB);
//	c2Symmetry(_axisA, _axisB);
//}
//
//void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {
//
//	 if (moveType == 0) {
//		// Z Shift
//		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center
//		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by
//
//		_trans.translate(_chainA, translateA);
//
//		c2Symmetry(_chainA, _chainB);
//		c2Symmetry(_axisA, _axisB);
//
//	} else if (moveType == 1) {
//		// Axial Rotation
//		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());
//
//		c2Symmetry(_chainA, _chainB);
//		c2Symmetry(_axisA, _axisB);
//
//	} else 	if (moveType == 2) {
//		// Crossing Angle
//		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
//		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
//
//		c2Symmetry(_chainA, _chainB);
//		c2Symmetry(_axisA, _axisB);
//
//	} else if (moveType == 3) {
//		// XShift
//		// Helix A interhelical distance
//		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center
//		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by
//
//		_trans.translate(_chainA, translateA);
//		_trans.translate(_axisA, translateA);
//
//		// Helix B interhelical distance
//		c2Symmetry(_chainA, _chainB);
//		c2Symmetry(_axisA, _axisB);
//
//	} else {
//		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
//	}
//}
//
//map<string,double> getEnergyByTerm(EnergySet* _eSet) {
//	// get all terms
//	map<string,double> eByTerm;
//	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
//	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
//		if(_eSet->isTermActive(it->first)) {
//			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
//		}
//	}
//	return eByTerm;
//}
//
//map<string,double> getEnergyByTermDoubled(EnergySet* _eSet) {
//	// get all terms
//	map<string,double> eByTerm;
//	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
//	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
//		if(_eSet->isTermActive(it->first)) {
//			eByTerm[it->first] =  2.0* _eSet->getTermEnergy(it->first);
//		}
//	}
//	return eByTerm;
//}
//
//void checkIfAtomsAreBuilt(System &_sys, ofstream &_err){
//	for (uint i=0; i<_sys.atomSize(); i++){
//		Atom atom = _sys.getAtom(i);
//		if (!atom.hasCoor()){
//			_err << "Atom " << i << " was not assigned coordinates; program termination";
//			cout << "Atom " << i << " was not assigned coordinates; program termination";
//			break;
//		} else {
//			continue;
//		}
//	}
//}
//
//string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum) {
//	string ps = "";
//	for (uint i=0; i<_sequence.length(); i++){
//		stringstream tmp;
//		tmp << _sequence[i];
//		string aa = tmp.str();
//		string resName = MslTools::getThreeLetterCode(aa);
//		if(resName == "HIS") {
//			resName = "HSE";
//		}
//		ps = ps + " " + resName;
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps;
//}
//string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum) {
//	string ps = "";
//	string resName = MslTools::getThreeLetterCode(_backboneAA);
//	if(resName == "HIS") {
//		resName = "HSE";
//	}
//	for (uint i=0; i<_backboneLength; i++){
//		ps = ps + " " + resName;
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps + "\nB" + ps;
//}
//
//
///***********************************
// *load rotamer functions
// ***********************************/
//void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot){
//	for (uint k=0; k<_sys.positionSize(); k++) {
//		Position &pos = _sys.getPosition(k);
//		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
//			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), "SL90.00")) {//lower rotamer level because I did baselines at this level
//				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
//			}
//		}
//	}
//}
//void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
//	for (uint k=0; k < _sys.positionSize(); k++) {
//		Position &pos = _sys.getPosition(k);
//
//		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
//			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
//				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
//			}
//		}
//	}
//}
//
//// converts a polymer sequence to a string for one chain (for homodimer sequences)
//string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
//	string seq = "";
//	for (uint i=0; i<_chain.positionSize(); i++){
//		string resName = _chain.getPosition(i).getCurrentIdentity().getResidueName();
//		string resID = MslTools::getOneLetterCode(resName);
//		seq += resID;
//	}
//	return seq;
//}
//
//// output energies by term into a referenced energyMap
//void outputEnergiesByTerm(SelfPairManager &_spm, vector<uint> _stateVec, map<string,double> &_energyMap,
//vector<string> _energyTermList, string _energyDescriptor, bool _includeIMM1){
//	if (_includeIMM1 == false){//No IMM1 Energy (for the Monte Carlos, both dimer and monomer)
//		for (uint i=0; i<_energyTermList.size(); i++){
//			string energyTerm = _energyTermList[i]; //CHARMM_ and SCWRL4_ terms
//			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;//Removes the CHARMM_ and SCWRL4_ before energyTerm names
//			if (energyTerm.find("IMM1") != string::npos){
//				continue;
//			} else {
//				if (_energyDescriptor.find("Monomer") != string::npos){
//					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
//				} else {
//					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
//				}
//			}
//		}
//		if (_energyDescriptor.find("Monomer") != string::npos){
//			//skip if monomer; could add calc baseline here at some point
//		} else {
//			_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
//			_energyMap["DimerSelfBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE");
//			_energyMap["DimerPairBaseline"] = _spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
//		}
//	} else if (_includeIMM1 == true){//IMM1 Energies
//		for (uint i=0; i<_energyTermList.size(); i++){
//			string energyTerm = _energyTermList[i];
//			string energyLabel = energyTerm.substr(7,energyTerm.length())+_energyDescriptor;
//			if (_energyDescriptor.find("Monomer") != string::npos){
//				if (energyTerm.find("IMM1") != string::npos){
//					_energyMap["IMM1Monomer"] = (_spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF"))*2;
//				} else {
//					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm)*2;
//				}
//			} else {
//				if (energyTerm.find("IMM1") != string::npos){
//					_energyMap["IMM1Dimer"] = _spm.getStateEnergy(_stateVec,"CHARMM_IMM1")+_spm.getStateEnergy(_stateVec,"CHARMM_IMM1REF");
//				} else {
//					_energyMap[energyLabel] = _spm.getStateEnergy(_stateVec, energyTerm);
//				}
//			}
//		}
//		//_energyMap["Baseline"] = _spm.getStateEnergy(_stateVec,"BASELINE")+_spm.getStateEnergy(_stateVec,"BASELINE_PAIR");
//	}
//}
//
////void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
////	for (uint k=0; k<_sys.positionSize(); k++) {
////		Position &pos = _sys.getPosition(k);
////		if (pos.identitySize() > 1){
////			for (uint j=0; j < pos.getNumberOfIdentities(); j++){
////				pos.setActiveIdentity(j);
////				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
////					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
////						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
////					}
////				}
////			}
////			pos.setActiveIdentity(0);
////		} else {
////			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
////				if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
////					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
////				}
////			}
////		}
////	}
////}
//
////below function only loads rotamers onto the interfacial positions by interfacialPositions (01 where 0 = non-interfacial and 1 = interfacial)
//void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<int> _interface){
//	for (uint k=0; k<_interface.size(); k++) {
//		if (_interface[k] < _numRotamerLevels){
//			Position &pos = _sys.getPosition(k);
//			if (pos.identitySize() > 1){
//				for (uint j=0; j < pos.getNumberOfIdentities(); j++){
//					pos.setActiveIdentity(j);
//					if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
//						if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
//							cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
//						}
//					}
//					pos.setActiveIdentity(0);
//				}
//			} else {
//				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
//					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
//						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
//					}
//				}
//			}
//		}
//	}
//}
//
//void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
//	_spm.setOnTheFly(1);
//	_spm.calculateEnergies(); // CHANGE BACK!!!
//	_spm.runGreedyOptimizer(_greedyCycles);
//}
//
///***********************************
// *functions from designFunctions
// ***********************************/
//string getAlternateIdString(vector<string> _alternateIds){
//	string alternateIdsString = "";
//	for (uint i=0; i<_alternateIds.size(); i++){
//		if (i == _alternateIds.size()-1){
//			alternateIdsString += _alternateIds[i];
//		} else {
//			alternateIdsString += _alternateIds[i] += " ";
//		}
//	}
//	return alternateIdsString;
//}
//
//
//string convertVectorUintToString(vector<uint> _inputVector){
//	string outputString = "";
//	for (uint i=0; i<_inputVector.size(); i++){
//		outputString += MslTools::intToString(_inputVector[i]);
//	}
//	return outputString;
//}
//
//// get the positions that will be linked on the interface (will have same AA identity and rotamer for self consistent mean field)
//vector<uint> getLinkedPositions(vector<uint> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel){
//	vector<uint> positionsToLink;
//	for (uint i=0; i<_rotamerSampling.size(); i++){
//		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
//			positionsToLink.push_back(1);
//		} else {
//			positionsToLink.push_back(0);
//		}
//	}
//	return positionsToLink;
//}
//
//// define the rotamer level for each position in the backbone
//vector<uint> convertStringToVectorUint(string _inputString){
//	vector<uint> outputVec;
//	for (uint i=0; i<_inputString.size(); i++){
//		stringstream ss;
//		ss << _inputString[i];
//		uint stringToInt = MslTools::toUnsignedInt(ss.str());
//		outputVec.push_back(stringToInt);
//	}
//	return outputVec;
//}
//
//// get a backbone sequence with an alanine cap at the beginning and end as an option
//string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) {
//	// initial start of sequence
//	string str = "";
//	//2021-09-21: add in an alanine cap to allow for more variable positions at the leucine region
//	for (uint i=0; i<_length-3; i++){
//		if (i<3){
//			if (_useAlaCap == true){
//				str = str + "A";
//			} else {
//				str = str + _backboneAA;
//			}
//		} else {
//			str = str + _backboneAA;
//		}
//	}
//	// Adds in the LILI at the end of the sequence which is necessary for our TOXCAT plasmids
//	if (_useAlaCap == true){
//		str = str + "AAA";
//	} else {
//		str = str + "ILI";
//	}
//	return str;
//}
//
//// generate string for backbone sequence
//string generateString(string _backbone, int _length) {
//	string str = "";
//	for (uint i=0; i<_length; i++){
//		str = str + _backbone;
//	}
//	return str;
//}
//
//string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions) {
//	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
//	// A:{32}ALA ILE GLY GLY GLY
//	// B:{32}ALA ILE GLY GLY GLY
//	string ps = "";
//	int counter = 0;
//	int startPos = _startResNum;
//	int endPos = _startResNum+_seq.length();
//	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
//		int pos = it-_seq.begin()+_startResNum;
//		if (it == _seq.begin() || it == _seq.end()-1){
//		//if (it == _seq.begin()){
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if (it == _seq.begin()){
//				if(resName == "HIS") {
//					ps = ps + " HSE-ACE";
//				} else {
//					ps = ps + " " + resName + "-ACE";
//				}
//			} else {
//				if(resName == "HIS") {
//					ps = ps + " HSE-CT2";
//				} else {
//					ps = ps + " " + resName + "-CT2";
//				}
//			}
//			counter++;
//		} else if (pos < startPos+3 || pos > endPos-5){
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if(resName == "HIS") {
//				ps = ps + " HSE";
//			} else {
//				ps = ps + " " + resName;
//			}
//		} else {
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			//cout << pos << endl;
//			if (find(_interfacialPositions.begin(), _interfacialPositions.end(), pos) != _interfacialPositions.end()){
//				ps = ps + " [";
//				for (uint i=0; i<_alternateIds.size(); i++){
//					if(_alternateIds[i] == "HIS") {
//						ps = ps + " HSE";
//					} else {
//						ps = ps + " " + _alternateIds[i];
//					}
//				}
//				ps = ps + "] ";
//			} else {
//				if(resName == "HIS") {
//					ps = ps + " HSE";
//				} else {
//					ps = ps + " " + resName;
//				}
//			}
//			counter++;
//		}
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps + "\nB" + ps;
//}
//
////Code Samson made a while back that should get each active ID and set a mask for anything that isn't active
//std::vector < std::vector < bool > > getActiveMask (System &_sys) {
//	_sys.updateVariablePositions();
//	std::vector <unsigned int> residueState;
//	std::vector < std::vector<unsigned int> > resRots(_sys.getMasterPositions().size());
//	std::vector < std::vector<bool> > resMask(_sys.getMasterPositions().size());
//	//Initialize residue state at the current active identity for each position
//	for (unsigned int i = 0; i < _sys.getMasterPositions().size(); i++) {
//		Position &pos = _sys.getPosition(_sys.getMasterPositions()[i]);
//		unsigned int activeRes = pos.getActiveIdentity();
//		residueState.push_back(activeRes);
//
//		resRots[i] = std::vector<unsigned int> (pos.identitySize());
//		for (unsigned int j = 0; j < pos.identitySize(); j++) {
//			resRots[i][j] = pos.getTotalNumberOfRotamers(j);
//		}
//	}
//
//	for (unsigned int i = 0; i < residueState.size(); i++) {
//		unsigned int activeResidue = residueState[i];
//		if (activeResidue >= resRots[i].size()) {
//			cerr << "ERROR: the current residue number exceeds the number of residues for position " << i << endl;
//			exit(100);
//		}
//		for (unsigned int j = 0; j < resRots[i].size(); j++) {
//			if (j==activeResidue) {
//				for (unsigned int k = 0; k < resRots[i][j]; k++) {
//					resMask[i].push_back(true);
//				}
//			} else {
//				for (unsigned int k = 0; k < resRots[i][j]; k++) {
//					resMask[i].push_back(false);
//				}
//			}
//		}
//
//		//Sanity check for presence of true rotamers
//
//		bool trueRots = false;
//		for (unsigned int j = 0; j < resMask[i].size(); j++) {
//			if (resMask[i][j]) {
//				trueRots = true;
//			}
//		}
//		if (!trueRots) {
//			cerr << "ERROR AT POSITION: " << i << endl;
//			cerr << "Current Residue: " << activeResidue << endl;
//			cerr << "resRots at this position: " << endl;
//			for (uint k = 0; k < resRots[i].size(); k++) {
//				cerr << resRots[i][k] << " ";
//			}
//			cerr << endl;
//			cerr << "resMask at this position: " << endl;
//			for (uint k = 0; k < resMask[i].size(); k++) {
//				cerr << resMask[i][k] << " ";
//			}
//			cerr << endl;
//			exit(9123);
//		}
//	}
//	return resMask;
//}
//
//string convertToPolymerSequence(string _seq, int _startResNum) {
//	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
//	// A:{32}ALA ILE GLY GLY GLY
//	// B:{32}ALA ILE GLY GLY GLY
//	string ps = "";
//	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
//		stringstream ss;
//		ss << *it;
//		string resName = MslTools::getThreeLetterCode(ss.str());
//		if(resName == "HIS") {
//			ps = ps + " HSE";
//		} else {
//			ps = ps + " " + resName;
//		}
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps + "\nB" + ps;
//}
//
//string convertToPolymerSequenceNeutralPatch(string _seq, int _startResNum) {
//	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
//	// A:{32}ALA ILE GLY GLY GLY
//	// B:{32}ALA ILE GLY GLY GLY
//	string ps = "";
//	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
//		if (it == _seq.begin() || it == _seq.end()-1){
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if (it == _seq.begin()){
//				if(resName == "HIS") {
//					ps = ps + " HSE-ACE";
//				} else {
//					ps = ps + " " + resName + "-ACE";
//				}
//			} else {
//				if(resName == "HIS") {
//					ps = ps + " HSE-CT2";
//				} else {
//					ps = ps + " " + resName + "-CT2";
//				}
//			}
//		} else {
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if(resName == "HIS") {
//				ps = ps + " HSE";
//			} else {
//				ps = ps + " " + resName;
//			}
//		}
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps + "\nB" + ps;
//}
//
//string convertToPolymerSequenceNeutralPatchMonomer(string _seq, int _startResNum) {
//	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
//	// A:{32}ALA ILE GLY GLY GLY
//	string ps = "";
//	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
//		if (it == _seq.begin() || it == _seq.end()-1){
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if (it == _seq.begin()){
//				if(resName == "HIS") {
//					ps = ps + " HSE-ACE";
//				} else {
//					ps = ps + " " + resName + "-ACE";
//				}
//			} else {
//				if(resName == "HIS") {
//					ps = ps + " HSE-CT2";
//				} else {
//					ps = ps + " " + resName + "-CT2";
//				}
//			}
//		} else {
//			stringstream ss;
//			ss << *it;
//			string resName = MslTools::getThreeLetterCode(ss.str());
//			if(resName == "HIS") {
//				ps = ps + " HSE";
//			} else {
//				ps = ps + " " + resName;
//			}
//		}
//	}
//	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
//	return "A" + ps;
//}
//
////Function to get the sum of a vector of doubles, typically energies
//double sumEnergyVector(vector<double> _energies){
//	double ener = 0;
//	for (uint i=0; i<_energies.size(); i++){
//		ener = ener + _energies[i];
//	}
//	return ener;
//}
//
//void resetEnergySet(System &_sys, vector<string> _energyTermList){
//	for (uint i=0; i<_energyTermList.size(); i++){
//		string energyTerm = _energyTermList[i];
//		_sys.getEnergySet()->eraseTerm(energyTerm);
//	}
//}
//
//void writePdb(System &_sys, string _outputDir, string _pdbName){
//	PDBWriter writer;
//	writer.open(_outputDir + "/" + _pdbName + ".pdb");
//	writer.write(_sys.getAtomPointers(), true, false, false);
//	writer.close();
//}
//