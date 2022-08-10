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


std::vector <unsigned int> getRandomOrder (unsigned int _size);
std::vector <unsigned int> getRandomOrder (unsigned int _start, unsigned int _end);
std::vector <unsigned int> getState (MSL::System &_Sys);


int main (int argc, char* argv[]) {

//	Build system and Interactions
	RandomNumberGenerator RNG;

//	This program takes a CRD-formatted protein structure file, not PDB
//	/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/O15162_01-crd.pdb
	string fullName = argv[1];
	System sys;
	if (sys.readPdb(fullName)) {
		cout << "Pdb read by sys successfully!" << endl;
		cout << sys << endl;
	} else {
		cout << "Did not read Pdb file..." << endl;
		exit(0);
	}
	


	SystemRotamerLoader SRL(sys, SYSENV.getEnv("MSL_ROTLIB"));
	CharmmSystemBuilder CSB;
	HydrogenBondBuilder SCWRL;


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
	
	cout <<"========================================="<<endl;
	cout <<"      Energy of structure                "<<endl;
	cout <<"========================================="<<endl;

	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;

	//load rotamers onto each position
	cout << "Loading rotamers..." << endl;

	//for each position in the chain
	for (uint i = 0; i < sys.positionSize(); i++) {
		Position &pos = sys.getPosition(i);
		string posId = pos.getPositionId();
		Residue res = pos.getCurrentIdentity();
		string resName = res.getResidueName();
		//do not load rotamers for glycine, alanine, and proline
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if(!SRL.loadRotamers(&pos, resName, "SL80.00")) {
			cout << "ERROR: Rotamers for position " << posId << ", residue " << resName << " NOT LOADED..." << endl;
			}
		}
	}
	//Calculate the Self Pair Manager Table of Interactions
	SelfPairManager spm(&sys);
	spm.setOnTheFly(true);
	spm.saveEnergiesByTerm(true);
	spm.calculateEnergies();


	vector <unsigned int> greedyState = greedySPM(sys, spm, RNG);


	//Load the rotamers found by greedySPM onto the system and calculate the system energy
	//Check to make sure SPM and sys give the same energies
	cout << "===============================" << endl;
	cout << "System Energy Calculation" << endl;
	cout << "===============================" << endl;
	sys.setActiveRotamers(greedyState);
	cout <<sys.calcEnergy() << endl;
	cout <<sys.getEnergySummary() << endl;

			
		

	return 0;
}

//greedySPM

std::vector<unsigned int> greedySPM (System &_sys, SelfPairManager &_spm, RandomNumberGenerator &_rng) {
	cout << "====================================" << endl;
	cout << "Begin SPM version of Greedy" << endl;
	cout << "====================================" << endl;

	vector<unsigned int> varPos = _sys.getVariablePositions();
	//The state is a vector where each position in the vector corresponds to a variable position
	//and the value of that position corresponds to the particular rotamer being queried
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
	cout << "========================================" << endl;
	cout << _spm.getStateEnergy(state)<< endl;
	cout << _spm.getSummary(state) << endl;
	return state;
}





