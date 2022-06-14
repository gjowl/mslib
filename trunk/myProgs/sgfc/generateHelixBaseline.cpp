//This program will make a baseline file based on a poly-alanine helix.
//A poly A helix will be generated and the energy calculated.
//The middle of the helix will be changed to each residue and the energy for each
//residue calculated.  The difference between the poly A and the mutant will
//be the baseline delta delta G energy for that residue.
//

#include "System.h"
#include "SelfPairManager.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "SystemRotamerLoader.h"
#include "PDBTopologyBuilder.h"
#include "SystemRotamerLoader.h"
#include "MslTools.h"
#include "HelixGenerator.h"
#include "SysEnv.h"
#include "MslTools.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

int main(int argc, char **argv) {
	System sys;
	CharmmSystemBuilder CSB;
	HydrogenBondBuilder SCWRL;
	PDBTopologyBuilder PDB (sys, SYSENV.getEnv("MSL_CHARMM_TOP"));
	SystemRotamerLoader SRL (sys, SYSENV.getEnv("MSL_ROTLIB"));


	CSB.setSystem(sys);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW", true);
	CSB.setBuildTerm("CHARMM_IMM1",true);


	SCWRL.setSystem(sys);
	string directory = "/exports/home/scondon/mslib/trunk/myProgs/sgfc/PDBs/helixBaseline/";
	string filename = "polyAla_monomer.pdb";
//	string directory = "/data06/scraven/mslib/trunk/myProgs/scraven/pdbs/";
//	string filename = "ftsb_TM_alignedZ-crd.pdb";
/*
	if (!sys.readPdb(directory + filename)) {
		cout << "PDB not read!" << endl;
		exit(5);
	}
*/
	

	if(!CSB.readTopology(SYSENV.getEnv("MSL_CHARMM_TOP"))) {

		cout << "ERROR_CSB: Did not read Topology File.." << endl;
		exit (1);
	}

	if(!CSB.readParameters(SYSENV.getEnv("MSL_CHARMM_PAR"))) {

		cout << "ERROR_CSB: Did not read Parameter File..." << endl;
		exit (1);
	}

	if(!CSB.readSolvation(SYSENV.getEnv("MSL_CHARMM_SOLV"))) {

		cout << "ERROR_CSB: Did not read Solvation file..." << endl;
		exit(1);
	}
	
	if(!CSB.buildSystemFromPDB(directory + filename)) {
		cout << "ERROR_CSB: Did not build system from PDB..." << endl;
		exit (1);
	}

	

	if(!SCWRL.readParameters(SYSENV.getEnv("MSL_HBOND_CA_PAR"))) {
		cout << "ERROR_SCWRL: Did not read H-Bond File..." << endl;
		exit(2);
	}
	

	if(!SCWRL.buildInteractions()) {
		cout << "ERROR_SCWRL: Did not build H-Bond Interactions..." << endl;
		exit(2);
	}	
	

	cout << sys << endl;
	cout << sys.calcEnergy() << endl;
	sys.printEnergySummary();

	std::vector<std::string> aaTmp(20,"");
	aaTmp[0] = "ALA";
	aaTmp[1] = "ARG";
	aaTmp[2] = "ASN";
	aaTmp[3] = "ASP";
	aaTmp[4] = "CYS";
	aaTmp[5] = "GLN";
	aaTmp[6] = "GLU";
	aaTmp[7] = "GLY";
	aaTmp[8] = "HSD";
	aaTmp[9] = "HSE";
	aaTmp[10] = "ILE";
	aaTmp[11] = "LEU";
	aaTmp[12] = "LYS";
	aaTmp[13] = "MET";
	aaTmp[14] = "PHE";
	//aaTmp[15] = "PRO";
	aaTmp[15] = "SER";
	aaTmp[16] = "THR";
	aaTmp[17] = "TRP";
	aaTmp[18] = "TYR";
	aaTmp[19] = "VAL";

	double wt = sys.calcEnergy();
	double mutant;
	double ddG;
	std::vector<string> output;

	for (uint i = 0; i < aaTmp.size(); i++) {
		//mutate Position 15 (the middle of the helix) to each residue and get the energy
		if (aaTmp[i] != "ALA") {
			Position &pos = sys.getPosition(14);
			string resName = pos.getResidue(0).getResidueName();
			CSB.addIdentity(pos, aaTmp[i]);
			CSB.removeIdentity(pos, resName);
			if(!SRL.loadRotamers(&pos, aaTmp[i], "SL90.00")) {
				cout << "ERROR: DID NOT LOAD ROTAMERS FOR " << aaTmp[i] << endl;
			}
		}
		CSB.updateNonBonded();
		SCWRL.buildInteractions();

		SelfPairManager SPM(&sys);
		SPM.calculateEnergies();
		SPM.runGreedyOptimizer(3);
		mutant = SPM.getStateEnergy(SPM.getMinStates()[0]);
		ddG = mutant - wt;
		
		string baseline = aaTmp[i] + " CA " + MslTools::doubleToString(ddG);
		output.push_back(baseline);

		sys.setActiveRotamers(SPM.getMinStates()[0]);
		sys.writePdb(directory+"baseline"+aaTmp[i]+".pdb");
		



	}

	for (uint i = 0; i < output.size(); i++) {
		cout << output[i] << endl;
	}
	/*
	CSB.updateNonBonded();
	sys.updateVariablePositions();
	cout << sys << endl;


	SelfPairManager SPM (&sys);
	SPM.setOnTheFly(true);

	SPM.calculateEnergies();


	for (uint i = 0; i < aaTmp.size(); i++) {
		Position &pos = sys.getPosition(14);
		pos.setActiveRotamer(i);
		sys.writePdb(directory+"baseline"+aaTmp[i]+".pdb");
		cout << sys.calcEnergy() << endl;
	}


	for (uint i = 0; i < aaTmp.size(); i++) {
		std::vector <uint> state(sys.getVariablePositions().size());
		state[0] = i;
		cout << SPM.getSummary(state) << endl;
	}
*/
	return 0;
}
