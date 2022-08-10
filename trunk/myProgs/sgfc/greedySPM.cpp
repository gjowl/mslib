#include "System.h"
#include "RandomNumberGenerator.h"
#include "SystemRotamerLoader.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "SysEnv.h"
#include "FormatConverter.h"
#include "SelfPairManager.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

vector <unsigned int> greedy (System &_sys, unsigned int _maxCycles = 3);
std::vector<unsigned int> greedySPM (System &_sys, SelfPairManager &_spm, RandomNumberGenerator &_rng);

//PUT GETRANDOMORDER INTO RANDOMNUMBERGENERATOR.H SOMETIME
std::vector <unsigned int> getRandomOrder (unsigned int _size);
std::vector <unsigned int> getRandomOrder (unsigned int _start, unsigned int _end);
std::vector <unsigned int> getState (MSL::System &_Sys);


int main () {

//	Build system and Interactions
	RandomNumberGenerator RNG;

	string pathName = "/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/";
	string fileName = "O15162_01-crd";
	string extension = ".pdb";
	string fullName = pathName + fileName + extension;

	System sys;
	if (sys.readPdb(fullName)) {
		cout << "Pdb read by sys successfully!" << endl;
		cout << sys << endl;
	} else {
		cout << "Did not read Pdb file..." << endl;
		exit(0);
	}
	

//	CharmmSystemBuilder CSB(sys, SYSENV.getEnv("MSL_CHARMM_TOP"), SYSENV.getEnv("MSL_CHARMM_PAR"));
//	HydrogenBondBuilder SCWRL(sys, SYSENV.getEnv("MSL_HBOND_PAR"));
//	SystemRotamerLoader SRL(sys, SYSENV.getEnv("MSL_ROTLIB"));
	CharmmSystemBuilder CSB;
	HydrogenBondBuilder SCWRL;
	SystemRotamerLoader SRL;

	CSB.setSystem(sys);
	CSB.readTopology(SYSENV.getEnv("MSL_CHARMM_TOP"));
	CSB.readParameters(SYSENV.getEnv("MSL_CHARMM_PAR"));

	if (CSB.buildSystemFromPDB(fullName)) {
		cout << "Pdb " << fullName << " read and built successfully!" << endl;
		cout << sys << endl;
	} else {
		cout << "ERROR: Pdb " << fullName << " not read..." << endl;
		exit(0);
	}
	
	SCWRL.setSystem(sys);
	SCWRL.readParameters(SYSENV.getEnv("MSL_HBOND_PAR"));

	if (SCWRL.buildInteractions()) {
		cout << "Hydrogen bond interactions built successfully!" << endl;
	} else {
		cout << "ERROR: Hydrogen bonds not built..." << endl;
		exit(0);
	}

	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;

//	SelfPairManager spm(&sys);
//	spm.calculateEnergies();
//	cout << spm.getFixedEnergy() << endl;

//	Load rotamers

	cout << "Loading rotamers..." << endl;

	//for each position in the chain
	for (uint i = 0; i < sys.positionSize(); i++) {
		Position &pos = sys.getPosition(i);
		string posId = pos.getPositionId();
		Residue res = pos.getCurrentIdentity();
		string resName = res.getResidueName();
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if(!SRL.loadRotamers(&pos, resName, "SL60.00")) {
			cout << "ERROR: Rotamers for position " << posId << ", residue " << resName << " NOT LOADED..." << endl;
			}
		}
	}

	SelfPairManager spm(&sys);
	spm.setOnTheFly(true);
	spm.calculateEnergies();

	spm.runGreedyOptimizer(5);
/*
	vector <unsigned int> greedy = greedySPM(sys, spm, RNG);
	sys.setActiveRotamers(greedy);
	if (sys.writePdb(pathName + fileName + "-greedySPM" + extension)) {
		cout << "Pdb written successfully!" << endl;
	} else {
		cout << "Did not write pdb..." << endl;
	}
*/
			
		

//	Begin greedy algorithms

//	cout << sys.getEnergySummary() << endl;
//	sys.setActiveRotamer("A,286",9);
//	cout << sys.getEnergySummary() << endl;
//	spm.setSystem(&sys);
//	spm.calculateEnergies();
//
//	vector <unsigned int> varPos = sys.getVariablePositions();
//	vector <unsigned int> state (sys.getVariablePositions().size(),0);
//	vector <unsigned int> rotState (sys.getVariablePositions().size(),1);
//	cout << spm.getSummary(state) << endl;
//	state[0] = 9;
//	cout << spm.getSummary(state) << endl;
//	cout << spm.getSummary(rotState) << endl;

	


/*
	vector < vector<unsigned int> > greedyVectors;
	vector<unsigned int> greedyCycle;
	for (uint i = 0; i < 5; i++) {
		greedyVectors.push_back(greedy(sys,3));
		cout << spm.getSummary(greedyVectors[i]) << endl;
		if (i!=0 && greedyVectors[i] == greedyVectors[i-1]) {
			cout << "Greedy cycles converged at cycle " << i << endl;
			if (sys.writePdb(pathName + fileName + "-greedy" + extension)) {
				cout << "Pdb file " << pathName + fileName + "-greedy" + extension << "written!" << endl;
			}
			break;
		}
	}
*/
	return 0;
}

//greedySPM

std::vector<unsigned int> greedySPM (System &_sys, SelfPairManager &_spm, RandomNumberGenerator &_rng) {
	cout << "====================================" << endl;
	cout << "Begin SPM version of Greedy" << endl;
	cout << "====================================" << endl;

	vector<unsigned int> varPos = _sys.getVariablePositions();
	vector<unsigned int> state(_sys.getVariablePositions().size(),0);
	vector<unsigned int> cycleState( _sys.getVariablePositions().size() );
	vector<unsigned int> rots = _spm.getNumberOfRotamers();
	double bestE = _spm.getStateEnergy(state);
	
	
	for (uint cycles = 0; cycles < 5; cycles++) {
		cout << "====================================" << endl;
		cout << "Cycle " << cycles << endl;
		cout << "====================================" << endl;

		double cycleE = _spm.getStateEnergy(cycleState);
		vector<unsigned int> order = _rng.getRandomOrder(varPos.size()); //the order that we will go through the varPos vector
		for (uint i = 0; i < order.size(); i++) {
			uint bestRot = 0;
			for (uint j = 0; j < rots[order[i]]; j++) {
				cycleState[order[i]] = j;
				double currentE = _spm.getStateEnergy(cycleState);
				if (currentE <= cycleE) {
					cycleE = currentE;
					bestRot = j;
				}
	
			}
			cycleState[order[i]] = bestRot;
			cout << "Best rotamer for position " << varPos[order[i]] << " is " << cycleState[order[i]] << " with current energy: " << cycleE << endl;
		}

		cout << "Cycle " << cycles <<" state: " << endl;
		for (uint i = 0; i < cycleState.size(); i++) {	
			cout << cycleState[i] << " ";
		}
		cout << endl;

		//check energy for this cycle: is it better than the previous best?
		if (state == cycleState && cycles > 0) {
			cout << "System converged with energy " << bestE << " at cycle " << cycles << endl;
			state = cycleState;

			break;
		}
		if (cycleE < bestE) {
			bestE = cycleE;
			state = cycleState;
		}



	}
	cout << _spm.getStateEnergy(state)<< endl;
	cout << _spm.getSummary(state) << endl;
	return state;
}

/*
vector<unsigned int> getState (System * _sys) {
	vector <unsigned int> varPos = _sys.getVariablePositions();
	vector <unsigned int> state (_sys.getVariablePositions().size());
	for (uint i = 0; i < varPos.size(); i++) {
		state[i] = _sys.getPosition(varPos[i]).getActiveIdentity();
		cout << "Position " << i << "; Rotamer " << state[i] << endl;
	}
	return state;
}
*/
/*
//	Begin greedy random
vector<unsigned int> greedy (System & _sys, unsigned int _maxCycles) {
	cout << "====================================" << endl;
	cout << "Begin Random Greedy" << endl;
	cout << "====================================" << endl;
	RandomNumberGenerator RNG;
	double bestE = _sys.calcEnergy();
	vector <unsigned int> bestRots;
	for (uint cycles = 0; cycles < _maxCycles; cycles++) {
		cout << "====================================" << endl;
		cout << "Cycle " << cycles << endl;
		cout << "====================================" << endl;
		
		vector<unsigned int> bestRandomRots(_sys.getVariablePositions().size());
		double minE = _sys.calcEnergy();
		vector<unsigned int> order = _sys.getVariablePositions();	
		for (vector<unsigned int>::iterator i = order.begin(); i < order.end(); ++i) {
			int currentPosition = *i;
			string posID = _sys.getPosition(currentPosition).getPositionId();
			uint bestPos = 0;
			for (uint j = 0; j < _sys.getPosition(currentPosition).getTotalNumberOfRotamers(); j++) {
				_sys.setActiveRotamer(posID, j);
				double currentE = _sys.calcEnergy();
				if (currentE <= minE) {
					bestPos = j;
					minE = currentE;
				}
			}
			//bestRandomRots[currentPosition] = bestPos;
			_sys.setActiveRotamer(posID, bestPos);
			cout << "Best rotamer for position " << currentPosition << ": " << bestPos << " with current energy: " << minE << endl;
		}
	

		//check energy for this cycle: is it better than the previous best?
		//if it's the same, exit
		if (abs(minE - bestE) < 0.00000001) {
			cout << "System converged with energy " << bestE << " at cycle " << cycles << endl;
			//bestRots = bestRandomRots;
			break;
		}
		if (minE < bestE) {
			bestE = minE;
			//bestRots = bestRandomRots;
		}
	}

	cout << "Best rotamers returned as vector:" << endl;
	cout << "====================================" << endl;
	for (std::vector<unsigned int>::iterator iterate = getState(_sys).begin(); iterate < getState(_sys).end(); ++iterate) {
		cout << *iterate << " ";
	}
	cout << endl << "====================================" << endl;
	
	return getState(_sys);
}
*/




