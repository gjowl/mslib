#include <sstream>
#include <iterator>
#include <unistd.h>
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv ENV;

// General functions that many of my programs use

// Just add 10 U(0,1) uniform random variables, offset by 0.5 to make mean = 0 and divide by variance = (10 * var(U(0,1)))
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

/***********************************
 *geometry
 ***********************************/
void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {
	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, zShift);
	_trans.translate(_apV, interDistVect);
	_trans.translate(_axis, interDistVect);
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
	Matrix m(3,3,0.0);
	m[0][0] = -1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = -1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;

	// Rotate chain B around Z axis
	Transforms trans;
	trans.rotate(_apvB, m);
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {

	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {

	 if (moveType == 0) {
		// Z Shift
		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by

		_trans.translate(_chainA, translateA);

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 1) {
		// Axial Rotation
		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else 	if (moveType == 2) {
		// Crossing Angle
		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 3) {
		// XShift
		// Helix A interhelical distance
		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center
		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);
		_trans.translate(_axisA, translateA);

		// Helix B interhelical distance
		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else {
		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
	}
}

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

void checkIfAtomsAreBuilt(System &_sys, ofstream &_err){
	for (uint i=0; i<_sys.atomSize(); i++){
		Atom atom = _sys.getAtom(i);
		if (!atom.hasCoor()){
			_err << "Atom " << i << " was not assigned coordinates; program termination";
			cout << "Atom " << i << " was not assigned coordinates; program termination";
			break;
		} else {
			continue;
		}
	}
}
string generateMonomerPolymerSequenceFromSequence(string _sequence, int _startResNum) {
	string ps = "";
	for (uint i=0; i<_sequence.length(); i++){
		stringstream tmp;
		tmp << _sequence[i];
		string aa = tmp.str();
		string resName = MslTools::getThreeLetterCode(aa);
		if(resName == "HIS") {
			resName = "HSE";
		}
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}
string generatePolymerSequence(string _backboneAA, int _backboneLength, int _startResNum) {
	string ps = "";
	string resName = MslTools::getThreeLetterCode(_backboneAA);
	if(resName == "HIS") {
		resName = "HSE";
	}
	for (uint i=0; i<_backboneLength; i++){
		ps = ps + " " + resName;
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}


/***********************************
 *load rotamer functions
 ***********************************/
void loadMonomerRotamers(System &_sys, SystemRotamerLoader &_sysRot){
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), "SL90.00")) {//lower rotamer level because I did baselines at this level
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
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
		} else {
			if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
				if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
					cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
				}
			}
		}
	}
}

//below function only loads rotamers onto the interfacial positions by interfacialPositions (01 where 0 = non-interfacial and 1 = interfacial)
void loadInterfacialRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL, int _numRotamerLevels, vector<int> _interface){
	for (uint k=0; k<_interface.size(); k++) {
		if (_interface[k] < _numRotamerLevels){
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
			} else {
				if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
					if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
						cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
					}
				}
			}
		}
	}
}

void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
}

/***********************************
 *functions from designFunctions
 ***********************************/
string getAlternateIdString(vector<string> _alternateIds){
	string alternateIdsString = "";
	for (uint i=0; i<_alternateIds.size(); i++){
		if (i == _alternateIds.size()-1){
			alternateIdsString += _alternateIds[i];
		} else {
			alternateIdsString += _alternateIds[i] += " ";
		}
	}
	return alternateIdsString;
}

vector<uint> getAllInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=0; k<_opt.backboneLength; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

vector<uint> getInterfacePositions(Options &_opt, vector<int> &_rotamerSamplingPerPosition){
	vector<uint> variableInterfacePositions;
	//TODO: make this variable in case I eventually decide that I actually want to mutate everything but the final Leu or something
	//for (uint k=0; k<_opt.backboneLength; k++){
	for (uint k=3; k<_opt.backboneLength-5; k++){
		if (_rotamerSamplingPerPosition[k] < _opt.interfaceLevel){
			variableInterfacePositions.push_back(k);
		} else {
			continue;
		}
	}
	return variableInterfacePositions;
}

string getInterfaceString(vector<int> _interface, int _seqLength){
	string interfaceString = "";
	for (uint i=0; i<_interface.size(); i++){
		if (i == _seqLength){
			i = _interface.size();
		} else {
			interfaceString += MslTools::intToString(_interface[i]);
		}
	}
	return interfaceString;
}

// get the positions that will be linked on the interface (will have same AA identity and rotamer for self consistent mean field)
vector<int> getLinkedPositions(vector<int> _rotamerSampling, int _interfaceLevel, int _highestRotamerLevel){
	vector<int> positionsToLink;
	for (uint i=0; i<_rotamerSampling.size(); i++){
		if (_rotamerSampling[i] < _interfaceLevel || _rotamerSampling[i] == _highestRotamerLevel){
			positionsToLink.push_back(1);
		} else {
			positionsToLink.push_back(0);
		}
	}
	return positionsToLink;
}

// define the rotamer level for each position in the backbone
vector<int> getRotamerSampling(string _rotamerLevels){
	vector<int> rotamerSampling;
	for (uint n=0; n<2; n++){
		for (uint i=0; i<_rotamerLevels.size(); i++){
			stringstream ss;
			ss << _rotamerLevels[i];
			rotamerSampling.push_back(MslTools::toInt(ss.str()));
		}
	}
	return rotamerSampling;
}

// get a backbone sequence with an alanine cap at the beginning and end as an option
string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) {
	string str = "";
	//2021-09-21: add in an alanine cap to allow for more variable positions at the leucine region
	for (uint i=0; i<_length-4; i++){
		if (i<4){
			if (_useAlaCap == true){
				str = str + "A";
			} else {
				str = str + _backboneAA;
			}
		} else {
			str = str + _backboneAA;
		}
	}
	// Adds in the LILI at the end of the sequence which is necessary for our TOXCAT plasmids
	str = str + "LILI";
	return str;
}

// calculate the SASA burial for each residue
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
std::vector<pair <int, double> > calculateResidueBurial (System &_sys, Options &_opt, string _seq) {
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
		AtomSelection selB(_sys.getPosition(i+_opt.backboneLength).getAtomPointers());
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
	writer.open(_opt.pdbOutputDir+"/interfaceSASA.pdb");
	writer.write(_sys.getAtomPointers(), true, false, true);
	writer.close();
	return residueBurial;
}

// generate string for backbone sequence
string generateString(string _backbone, int _length) {
	string str = "";
	for (uint i=0; i<_length; i++){
		str = str + _backbone;
	}
	return str;
}
