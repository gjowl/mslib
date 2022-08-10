//This function will take a protein sequence with multiple possible residues at a position and make a mask so that only the current residue is examined.  This will be useful for mutational analysis where we will want to first find the energy of a wild type sequence, then mutate a position.  We can load the rotamers of all identities all at once and then not have to loop through the entire set of rotamers.
//
//
#include "System.h"
#include "SystemRotamerLoader.h"
#include "SysEnv.h"
#include "PDBTopologyBuilder.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;


std::vector< std::vector<bool> > getResidueMask (System &_sys);

int main () {

	System sys;
	std::vector< std::vector<bool> > resMask; 

	if (sys.readPdb("/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/4HRO-crd_greedy.pdb")) {
		cout << "PDB read successfully!" << endl;
	} else {
		cout << "Did not read the file...try again?" << endl;
		exit(1);
	}
	sys.seed();
	sys.buildAtoms();
	SystemRotamerLoader SRL(sys, SYSENV.getEnv("MSL_ROTLIB"));
	PDBTopologyBuilder PTB(sys, SYSENV.getEnv("MSL_CHARMM_TOP"));

	cout << "=========================================" << endl;
	cout << "No rotamers loaded:" << endl;
	cout << "=========================================" << endl;
	resMask = getResidueMask(sys);
	for (uint i = 0; i < resMask.size(); i++) {
		cout << "Position " << i << ": ";
		for (uint j = 0; j < resMask[i].size(); j++) {
			cout << resMask[i][j] << " ";
		}
		cout << endl;
	}

	resMask.clear();


	cout << "=========================================" << endl;
	cout << "Loading Rotamers..." << endl;
	cout << "=========================================" << endl;


	//	Load rotamers


	//for each position in the chain
	for (uint i = 0; i < sys.positionSize(); i++) {
		Position &pos = sys.getPosition(i);
		string posId = pos.getPositionId();
		Residue res = pos.getCurrentIdentity();
		string resName = res.getResidueName();
		if(!PTB.addIdentity(pos, "ARG", "N CA C O H")) {
			cout << "Did not add ARG to the list of identities..." << endl;
		}
		if(!PTB.addIdentity(pos, "TRP", "N CA C O H")) {
			cout << "Did not add ARG to the list of identities..." << endl;
		}
		if(!SRL.loadRotamers(&pos, "ARG", 3)) {
			cout << "ARG rotamers not loaded..." << endl;
		}
		if(!SRL.loadRotamers(&pos, "TRP", 5)){
			cout << "TRP rotamers not loaded..." << endl;
		}

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if(SRL.loadRotamers(&pos, resName, "SL60.00")) {
				cout << "Rotamers for position " << posId << ", residue " << resName << " loaded!" << endl;
			} else {
			cout << "ERROR: Rotamers for position " << posId << ", residue " << resName << " NOT LOADED..." << endl;
			}
		}

		if (i == 6) {
			cout << "=========================================" << endl;
			cout << "Rotamers loaded for first seven residues" << endl;
			cout << "=========================================" << endl;
			
			resMask = getResidueMask(sys);
			for (uint i = 0; i < resMask.size(); i++) {
				cout << "Position " << i << ": ";
				for (uint j = 0; j < resMask[i].size(); j++) {
					cout << resMask[i][j] << " ";
				}
				cout << endl;
			}
		
			resMask.clear();
		}




	}



	cout << "=========================================" << endl;
	cout << "SL60.00, 3 ARG and 5 TRP rotamers added to each position" << endl;
	cout << "=========================================" << endl;

	resMask = getResidueMask(sys);
	for (uint i = 0; i < resMask.size(); i++) {
		cout << "Position " << i << ": ";
		for (uint j = 0; j < resMask[i].size(); j++) {
			cout << resMask[i][j] << " ";
		}
		cout << endl;
	}

	resMask.clear();

	cout << "=========================================" << endl;
	cout << "Change A,65 to ARG identity" << endl;
	cout << "=========================================" << endl;

	sys.setActiveRotamer("A,65,ARG", 2);

	resMask = getResidueMask(sys);
	for (uint i = 0; i < resMask.size(); i++) {
		cout << sys.getPosition(i).getPositionId()<< ": ";
		for (uint j = 0; j < resMask[i].size(); j++) {
			cout << resMask[i][j] << " ";
		}
		cout << endl;
	}
	resMask.clear();
	
	return 0;
}

std::vector< std::vector<bool> > getResidueMask (System &_sys) {
	std::vector< std::vector<bool> > mask(_sys.positionSize());
	
	//for each position
	for (unsigned int i = 0; i < _sys.positionSize(); i++) {
		string posID = _sys.getPosition(i).getPositionId();
		string identity = _sys.getPosition(i).getCurrentIdentity().getIdentityId();
		//cout << "Identity at position " << posID << ": " << identity << endl;

		//set up the mask for this particular position
		mask[i] = std::vector<bool> (_sys.getTotalNumberOfRotamers(i));

		//for each residue
		for (unsigned int j = 0; j < _sys.getTotalNumberOfRotamers(i); j++) {
			//for some reason sys.setActiveRotamer(i) isn't working but when I get the position first it works...
			_sys.getPosition(i).setActiveRotamer(j);
			//cout << _sys.getPosition(i).getCurrentIdentity().getIdentityId() << " " << _sys.getPosition(i).getRotamerId() << endl;
			if (identity == _sys.getPosition(i).getCurrentIdentity().getIdentityId()) {
				mask[i][j] = true;
			} else {
				mask[i][j] = false;
			}
		}
		_sys.getPosition(i).setActiveRotamer(0); //reset the rotamer to its original position
	}
	return mask;	


}

