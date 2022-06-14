#include "System.h"
#include "RandomNumberGenerator.h"
#include "SystemRotamerLoader.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "SysEnv.h"
#include "FormatConverter.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

vector <unsigned int> greedy (System &_sys, unsigned int _maxCycles = 3);

//PUT GETRANDOMORDER INTO RANDOMNUMBERGENERATOR.H SOMETIME
std::vector <unsigned int> getRandomOrder (unsigned int _size);
std::vector <unsigned int> getRandomOrder (unsigned int _start, unsigned int _end);



int main () {

//	Build system and Interactions

	string pathName = "/exports/home/scondon/mslib/trunk/myProgs/sgfc/greedySimple/";
	string fileName = "O15162_01-crd";
	string extension = ".pdb";
	string fullName = pathName + fileName + extension;

	System sys;
	if (sys.readPdb(fullName)) {
		cout << "Pdb read by sys successfully!" << endl;
		cout << sys << endl;
	} else {
		exit(0);
	}


//	FormatConverter format("PDB3", "CHARMM22");
//	format.convert (sys.getAtomPointers());
//	sys.writePdb(pathName + fileName + "-converted" + extension);

//	CharmmSystemBuilder CSB(sys, SYSENV.getEnv("MSL_CHARMM_TOP"), SYSENV.getEnv("MSL_CHARMM_PAR"));
	CharmmSystemBuilder CSB(sys, "/exports/home/scondon/mslib/trunk/toppar/charmm22.top", "/exports/home/scondon/mslib/trunk/toppar/charmm22.par");
	CSB.setBuildTerm("CHARMM_ELEC", false);
//	CSB.setBuildTerm("CHARMM_VDW", false);
	HydrogenBondBuilder SCWRL(sys, SYSENV.getEnv("MSL_HBOND_PAR"));
	SystemRotamerLoader SRL(sys, SYSENV.getEnv("MSL_ROTLIB"));

	if (CSB.buildSystemFromPDB(fullName)) {
		cout << "Pdb " << fullName << " read and built successfully!" << endl;
		cout << sys << endl;
	} else {
		cout << "ERROR: Pdb " << fullName << " not read..." << endl;
		exit(0);
	}
	


	if (SCWRL.buildInteractions()) {
		cout << "Hydrogen bond interactions built successfully!" << endl;
	} else {
		cout << "ERROR: Hydrogen bonds not built..." << endl;
		exit(0);
	}

	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary() << endl;
//	Load rotamers

	cout << "Loading rotamers..." << endl;
//	for (uint i = 0; i < sys.chainSize(); i++) {
//		Chain chainNum = sys.getChain(i);
		for (uint j = 0; j < sys.positionSize(); j++) {
	
			Position &pos = sys.getPosition(j);
			string posId = pos.getPositionId();
			Residue res = pos.getCurrentIdentity();
			string resName = res.getResidueName();
			if(SRL.loadRotamers(&pos, pos.getResidueName(), "SL90.00")) {
				cout << "Rotamers for position " << pos.getPositionId() << ", residue " << resName << " loaded!" << endl;
			} else {
			cout << "ERROR: Rotamers for position " << posId << ", residue " << resName << " NOT LOADED..." << endl;
			}
		}
//	}

//	Begin greedy algorithms
	vector < vector<unsigned int> > greedyVectors;
	vector<unsigned int> greedyCycle;
	for (uint i = 0; i < 5; i++) {
		greedyVectors.push_back(greedy(sys,3));
		if (i!=0 && greedyVectors[i] == greedyVectors[i-1]) {
			cout << "Greedy cycles converged at cycle " << i << endl;
			if (sys.writePdb(pathName + fileName + "-greedy" + extension)) {
				cout << "Pdb file " << pathName + fileName + "-greedy" + extension << "written!" << endl;
			}
			break;
		}
	}
	return 0;
}





//getRandomOrder
/*
std::vector <unsigned int> getRandomOrder (unsigned int _size) {
	unsigned int _start = 0;
	return getRandomOrder (_start, _size); 
}


std::vector<unsigned int> getRandomOrder (unsigned int _start, unsigned int _end) {

//	Debug comments to terminal
//	cout << "==============================" << endl;
//	cout << "Ordered vector" << endl;
//	cout << "==============================" << endl;

	RandomNumberGenerator RNG;
	std::vector <unsigned int> ordered;
	std::vector <unsigned int> random;
	for (unsigned int i = _start; i < _end; i++) {
		ordered.push_back (i);
//		cout << ordered.back() << endl;
	}

//	Debug comments to terminal
//	cout << "==============================" << endl;
//	cout << "Random vector" << endl;
//	cout << "==============================" << endl;
	

	while (!ordered.empty()) {
		double randInit = RNG.getRandomDouble();
		double randNum = randInit * ordered.size();
		int randInt = randNum;
		random.push_back (ordered[randInt]);
		ordered.erase(ordered.begin() + randInt);
//		cout << random.back() << endl;
	}

	return random;

}
*/
//greedy
/*
bool greedy(System _sys) {
	for (uint i = 0; i < _sys.positionSize(); i++) {
		string posID = _sys.getPosition(i).getPositionId();
		for (uint j = 0; j < _sys.getPosition(i).getTotalNumberOfRotamers(); j++) {
			_sys.setActiveRotamer(posID, j);

		}
	}

	return true;
}
*/



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
		
		vector<unsigned int> bestRandomRots(_sys.positionSize());
		double minE = _sys.calcEnergy();
		vector<unsigned int> order = RNG.getRandomOrder(_sys.positionSize());
	
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
			bestRandomRots[currentPosition] = bestPos;
			_sys.setActiveRotamer(posID, bestPos);
			cout << "Best rotamer for position " << currentPosition << ": " << bestPos << " with current energy: " << minE << endl;
		}
	

		//check energy for this cycle: is it better than the previous best?
		//if it's the same, exit
		if (abs(minE - bestE) < 0.00000001) {
			cout << "System converged with energy " << bestE << " at cycle " << cycles << endl;
			bestRots = bestRandomRots;
			break;
		}
		if (minE < bestE) {
			bestE = minE;
			bestRots = bestRandomRots;
		}
	}

	cout << "Best rotamers returned as vector:" << endl;
	cout << "====================================" << endl;
	for (std::vector<unsigned int>::iterator iterate = bestRots.begin(); iterate < bestRots.end(); ++iterate) {
		cout << *iterate << " ";
	}
	cout << endl << "====================================" << endl;
	
	return bestRots;
}






/*
	cout << "====================================" << endl;
	cout << "Begin Simple Greedy" << endl;
	cout << "====================================" << endl;
	
	vector<unsigned int> bestRotamers;
	double minE = sys.calcEnergy();

	for (uint i = 0; i < sys.positionSize(); i++) {
		string posID = sys.getPosition(i).getPositionId();
		uint bestPos = 0;
		for (uint j = 0; j < sys.getPosition(i).getTotalNumberOfRotamers(); j++) {
			sys.setActiveRotamer(posID, j);
			double currentE = sys.calcEnergy();
			if (currentE <= minE) {
				bestPos = j;
				minE = currentE;
			}
		}
		bestRotamers.push_back(bestPos);
		sys.setActiveRotamer(posID, bestPos);
		cout << "Best rotamer for position " << i << ": " << bestPos << endl;
		
		
	}
	
	cout << sys.getEnergySummary() << endl;
	sys.writePdb(pathName + "bestGreedy.pdb");
*/

