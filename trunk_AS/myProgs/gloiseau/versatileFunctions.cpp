#include <sstream>
#include <iterator>
#include <unistd.h>
#include "versatileFunctions.h"

using namespace std;
using namespace MSL;

/***********************************
 * System Functions
 ***********************************/
void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	for (uint k=0; k < _sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
        // check if there is more than one identity for this position
        if (pos.identitySize() > 1) {
            for (uint i=0; i<pos.identitySize(); i++){
                // set the active identity
			    pos.setActiveIdentity(i);
		        if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
		        	if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
		        		cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
		        	}
		        }
            }
        } else if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(),_SL)) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}
}

void deleteTerminalBondInteractions(System &_sys, vector<string> &_deleteTerminalInteractionList){
	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize();
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain & thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		int firstResidueNumber = 0; // hardcoded, this is done stupidly but works
		int lastResidueNumber = positions[positions.size()-1]->getResidueNumber();
		for(int i = firstResidueNumber; i < firstResidueNumber+3; i++) {
			// rid of hbonds from first 3 positions
			if(firstResidueNumber <= i) {
				atoms += positions[i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			// rid of hbonds from last 3 positions
			if(lastResidueNumber > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	for (uint i=0; i<_deleteTerminalInteractionList.size(); i++){
		pESet->deleteInteractionsWithAtoms(atoms,_deleteTerminalInteractionList[i]);
	}
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

// set the active identity for each position to the identity in the given sequence
void setActiveSequence(System &_sys, string _sequence){
	// loop through the sequence
	for (uint i=0; i<_sequence.size(); i++){
		// get the ith residue in the sequence
		string res = _sequence.substr(i,1);
		// loop through all of the chains in the system
		for (uint j=0; j<_sys.chainSize(); j++){
			Chain &chain = _sys.getChain(j);
			// get the ith position in the system
			Position &pos = chain.getPosition(i);
			// get the position id for the ith position
			string posId = pos.getPositionId();
			// convert the residue id to three letter code
			string aa = MslTools::getThreeLetterCode(res);
			// set active identity
			_sys.setActiveIdentity(posId, aa);
		}
	}
}

// Code Samson made a while back that should get each active ID and set a mask for anything that isn't active
std::vector < std::vector < bool > > getActiveMask (System &_sys) {
	_sys.updateVariablePositions();
	std::vector <unsigned int> residueState;
	std::vector < std::vector<unsigned int> > resRots(_sys.getMasterPositions().size());
	std::vector < std::vector<bool> > resMask(_sys.getMasterPositions().size());
	//Initialize residue state at the current active identity for each position
	for (unsigned int i = 0; i < _sys.getMasterPositions().size(); i++) {
		Position &pos = _sys.getPosition(_sys.getMasterPositions()[i]);
		unsigned int activeRes = pos.getActiveIdentity();
		residueState.push_back(activeRes);

		resRots[i] = std::vector<unsigned int> (pos.identitySize());
		for (unsigned int j = 0; j < pos.identitySize(); j++) {
			resRots[i][j] = pos.getTotalNumberOfRotamers(j);
		}
	}

	for (unsigned int i = 0; i < residueState.size(); i++) {
		unsigned int activeResidue = residueState[i];
		if (activeResidue >= resRots[i].size()) {
			cerr << "ERROR: the current residue number exceeds the number of residues for position " << i << endl;
			exit(100);
		}
		for (unsigned int j = 0; j < resRots[i].size(); j++) {
			if (j==activeResidue) {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(true);
				}
			} else {
				for (unsigned int k = 0; k < resRots[i][j]; k++) {
					resMask[i].push_back(false);
				}
			}
		}

		//Sanity check for presence of true rotamers

		bool trueRots = false;
		for (unsigned int j = 0; j < resMask[i].size(); j++) {
			if (resMask[i][j]) {
				trueRots = true;
			}
		}
		if (!trueRots) {
			cerr << "ERROR AT POSITION: " << i << endl;
			cerr << "Current Residue: " << activeResidue << endl;
			cerr << "resRots at this position: " << endl;
			for (uint k = 0; k < resRots[i].size(); k++) {
				cerr << resRots[i][k] << " ";
			}
			cerr << endl;
			cerr << "resMask at this position: " << endl;
			for (uint k = 0; k < resMask[i].size(); k++) {
				cerr << resMask[i][k] << " ";
			}
			cerr << endl;
			exit(9123);
		}
	}
	return resMask;
}


/***********************************
 * EnergySet Functions
 ***********************************/
void resetEnergySet(System &_sys, vector<string> _energyTermList){
	for (uint i=0; i<_energyTermList.size(); i++){
		string energyTerm = _energyTermList[i];
		_sys.getEnergySet()->eraseTerm(energyTerm);
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

// function for writing a pdb
void writePdb(System &_sys, string _outputDir, string _pdbName){
	PDBWriter writer;
	writer.open(_outputDir + "/" + _pdbName + ".pdb");
	writer.write(_sys.getAtomPointers(), true, false, false);
	writer.close();
}

// function to get the sum of a vector of doubles, typically energies
double sumDoubleVector(vector<double> _vector){
	double sum = 0;
	for (uint i=0; i<_vector.size(); i++){
		sum = sum + _vector[i];
	}
	return sum;
}

void repackSideChains(SelfPairManager & _spm, int _greedyCycles) {
	_spm.setOnTheFly(1);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
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
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical axis times the amount to shift by

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
// Just add 10 U(0,1) uniform random variables, offset by 0.5 to make mean = 0 and divide by variance = (10 * var(U(0,1)))
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

// adjust the move sizes based on the monte carlo function (decrease by multiplying by the change in temperature)
void getCurrentMoveSizes(double &_currTemp, double &_endTemp, double &_deltaX, double &_deltaCross, double &_deltaAx, double &_deltaZ,
 double _deltaXLimit, double _deltaCrossLimit, double _deltaAxLimit, double _deltaZLimit, bool &_decreaseMoveSize, int _moveToPerform) {
	// define the temperature change
	double decreaseMultiplier = _endTemp/_currTemp;
	// setup the decrease booleans; if the move size has reached its minimum, don't decrease it further and set to false
	bool decreaseZ = false;
	bool decreaseAx = false;
	bool decreaseCross = false;
	bool decreaseX = false;

	// move to perform is the move that was just performed from the backboneOptimizer function
	if (_moveToPerform == 0){
		decreaseZ = true;
		_deltaZ = decreaseMoveSize(_deltaZ, _deltaZLimit, decreaseMultiplier, decreaseZ);
	} else if (_moveToPerform == 1) {
		decreaseAx = true;
		_deltaAx = decreaseMoveSize(_deltaAx, _deltaAxLimit, decreaseMultiplier, decreaseAx);
	} else if (_moveToPerform == 2) {
		decreaseCross = true;
		_deltaCross = decreaseMoveSize(_deltaCross, _deltaCrossLimit, decreaseMultiplier, decreaseCross);
	} else if (_moveToPerform == 3) {
		decreaseX = true;
		_deltaX = decreaseMoveSize(_deltaX, _deltaXLimit, decreaseMultiplier, decreaseX);
	}

	// if all of the move sizes have reached their minimum, set the decrease move size boolean to false	
	if (decreaseX == false && decreaseCross == false && decreaseAx == false && decreaseZ == false){
		_decreaseMoveSize = false;
	}
}

double decreaseMoveSize(double _moveSize, double _moveLimit, double _decreaseMultiplier, bool &_decrease) {
	// edited to make sure that the move size is decreasing properly down to the move limit: add in detail here
	double diffMoveSize = _moveSize-_moveLimit;
	double moveDecrease = diffMoveSize-(diffMoveSize*_decreaseMultiplier);// edited on 2022-10-6: it now works properly; before it decreased to the move limit
	double newMoveSize = _moveSize - moveDecrease;
	if (newMoveSize > _moveLimit){
		return newMoveSize;
	} else {
		_decrease = false;
		return _moveSize;
	}
}

/***********************************
* general functions
 ***********************************/
// set the active identity for each position to the identity in the given sequence (only for homodimers)
string extractSequence(System &_sys){
	// initialize the sequence string
	string sequence = "";
	// get the first chain from the system
	Chain &chain = _sys.getChain(0);
	// loop through the chain
	for (uint i=0; i<chain.positionSize(); i++){
		// get the ith position in the system
		Position &pos = chain.getPosition(i);
		// get the residue name of the ith position
		string res = pos.getResidueName();
		// convert the residue name to one letter code
		string aa = MslTools::getOneLetterCode(res);
		// add the one letter code to the sequence string
		sequence += aa;
	}
	return sequence;
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
	return "A" + ps + "\nB" + ps;
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
	return "A" + ps + "\nB" + ps;
}

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

string convertVectorUintToString(vector<uint> _inputVector){
	string outputString = "";
	for (uint i=0; i<_inputVector.size(); i++){
		outputString += MslTools::intToString(_inputVector[i]);
	}
	return outputString;
}

// get a backbone sequence with an alanine cap at the beginning and end as an option
string generateBackboneSequence(string _backboneAA, int _length, bool _useAlaCap) {
	// initial start of sequence
	string str = "";
	//2021-09-21: add in an alanine cap to allow for more variable positions at the leucine region
	for (uint i=0; i<_length-3; i++){
		if (i<3){
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
	if (_useAlaCap == true){
		str = str + "AAA";
	} else {
		str = str + "ILI";
	}
	return str;
}

// generate string for backbone sequence
string generateString(string _backbone, int _length) {
	string str = "";
	for (uint i=0; i<_length; i++){
		str = str + _backbone;
	}
	return str;
}

string generateMultiIDPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _interfacialPositions) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	int startPos = _startResNum;
	int endPos = _startResNum+_seq.length();
	for(string::iterator it = _seq.begin(); it != _seq.end(); it++) {
		int pos = it-_seq.begin()+_startResNum;
		if (it == _seq.begin() || it == _seq.end()-1){
		//if (it == _seq.begin()){
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
		} else if (pos < startPos+3 || pos > endPos-5){
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		} else {
			stringstream ss;
			ss << *it;
			string resName = MslTools::getThreeLetterCode(ss.str());
			//cout << pos << endl;
			if (find(_interfacialPositions.begin(), _interfacialPositions.end(), pos) != _interfacialPositions.end()){
				ps = ps + " [";
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
	return "A" + ps + "\nB" + ps;
}

string convertPolymerSeqToOneLetterSeq(Chain &_chain) {
	string seq = "";
	for (uint i=0; i<_chain.positionSize(); i++){
		string resName = _chain.getPosition(i).getCurrentIdentity().getResidueName();
		string resID = MslTools::getOneLetterCode(resName);
		seq += resID;
	}
	return seq;
}

// define the rotamer level for each position in the backbone
vector<uint> convertStringToVectorUint(string _inputString){
	vector<uint> outputVec;
	for (uint i=0; i<_inputString.size(); i++){
		stringstream ss;
		ss << _inputString[i];
		uint stringToInt = MslTools::toUnsignedInt(ss.str());
		outputVec.push_back(stringToInt);
	}
	return outputVec;
}
