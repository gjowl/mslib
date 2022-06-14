#include <iostream>
#include <fstream>
#include <sstream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"


using namespace MSL;
using namespace std;

struct Options{
	string sequence;
	string helixGeoFile;

	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string AA;

	int seed;
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	int MCCycles;
        int MCMaxRejects;	
	int AAi;
	int AA2;

	vector<string> allowed, required;

	string configfile;
};

void replaceAA(CharmmSystemBuilder &_CSB, System &_sys, string _AA, int _AA1, int _AA2){
	//add in if statement for if AA ints not within sys
	if (_AA1 > _sys.getChain("A").positionSize()){
		cerr << "Amino Acid Position " << _AA1 << "is not in System." << endl;
		exit(0);
	}

	/*for (uint i=0; i<_sys.getChain("A").positionSize(); i++){	
		Position &pos1 = _sys.getChain("A").getPosition(i);
		Position &pos2 = _sys.getChain("B").getPosition(i);
		if (pos1.getActiveIdentity() != 0){
			_CSB.removeIdentity(pos1, "HSZ");
		}
		if (pos2.getActiveIdentity() != 0){
			_CSB.removeIdentity(pos2, "HSZ");
		}
	}*/	

	Position &pos1 = _sys.getChain("A").getPosition(_AA1);
	Position &pos2 = _sys.getChain("B").getPosition(_AA1);
	_CSB.addIdentity(pos1, _AA);//maybe I should add ID to all, and then just switch?
	_CSB.addIdentity(pos2, _AA);
	pos1.setActiveIdentity(1, false);
	pos2.setActiveIdentity(1, false);
	
	if (_AA2 != 1 || _AA2 != 3 || _AA2 != 4){//this works right but is it right...? need to find a way to do it for all...? maybe a for loop with int _AA2 being i<=4 (2 just shouldn't happen)
		Position &pos3 = _sys.getChain("A").getPosition(_AA1+_AA2);
		Position &pos4 = _sys.getChain("B").getPosition(_AA1+_AA2);
		_CSB.addIdentity(pos3, _AA);
		_CSB.addIdentity(pos4, _AA);
		pos3.setActiveIdentity(1, false);
		pos4.setActiveIdentity(1, false);
	}
	//_sys.buildAllAtoms();
}

int main(int argc, char *argv[]){	

	Options defaults;

	defaults.sequence = "GGGGGGGGGGGGGGG";
	defaults.helixGeoFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/CENTROIDS_0.5_cutoff-6_byHbond";
	defaults.pdbOutputDir = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau";
	defaults.topFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/toppar_HSZ/top_all22_prot.inp";
	defaults.parFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/toppar_HSZ/par_all22_prot.inp";
	defaults.solvFile = "";
	defaults.hBondFile = "/data01/sabs/tmRepacks/charmm/par_hbond_CA_2.txt";
	defaults.rotLibFile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/EBL_11-2011_CHARMM22_Zinc.txt";
	defaults.seed = 1;
	defaults.weight_vdw = 1;
	defaults.weight_hbond = 1;
	defaults.weight_solv = 1;
	defaults.backboneCrd = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd";
	defaults.MCCycles = 1;
	defaults.MCMaxRejects = 10;
	defaults.AA = "MEZ";
	defaults.AAi = 4;
	defaults.AA2 = 4;
	
	//string infile = "/exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/inputPDB_tomin.pdb";
	string infile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	Options opt = defaults;

	/****************Open Output File************************/
	ofstream pout;
	string poutName = opt.pdbOutputDir + "/CATM_generated.out";
	
	pout.open(poutName.c_str());


	/*****************Declare System*************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);
	
	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	CSB.buildSystemFromPDB(infile);
	PDBWriter writer;
	PDBWriter writer2;
	sys.saveAltCoor("start");
	
	writer.open("/exports/home/gloiseau/MEZ_rotamers.pdb");
	writer.write(sys.getAtomPointers(), true, false, true);

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	vector<Position*>& positionsChainA = sys.getChain("A").getPositions();
	vector<Position*>& positionsChainB = sys.getChain("B").getPositions();

	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd); 
	if(!cRead.read()) {
		cout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();
	AtomPointerVector& glyAPV = cRead.getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(opt.seed);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	//sys.wipeAllCoordinates();
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();
	
	writer.write(sys.getAtomPointers(), true, false, true);
	//writer.write(sys.getChain("B").getAtomPointers(), true, false, true);
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	hb.buildInteractions(30);//maybe the +charges?

	replaceAA(CSB, sys, opt.AA, 1, 1);
	
	/******************************************************************************
	 *                     === LOAD ROTAMERS ===
	 ******************************************************************************/
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!sysRot.loadRotamers(&pos, pos.getResidueName(),"SL95.00")) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

	int rotTotal = positionsChainA[1]->getTotalNumberOfRotamers();

	//stringstream ss;
	for (uint i=0; i<rotTotal; i++){
		stringstream ss;
	        ss << setfill('0') << setw(2) << i;
		//ss << i;
		positionsChainA[1]->setActiveRotamer("MEZ",i);
		writer.write(sys.getAtomPointers(), true, false, true);
		//string l = "/exports/home/gloiseau/Generated_PDBs/model" + ss.str() + ".pdb"; 
		//cout << "L: " << l << endl;
		//writer2.open();
		//writer2.write(sys.getAtomPointers(), true, false, true);
		//writer2.close();
	}

	pout.close();
	writer.close();
}
