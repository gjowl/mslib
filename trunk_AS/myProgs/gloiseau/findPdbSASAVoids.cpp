#include <sstream>
#include <iterator>
#include <unistd.h>

// MSL Functions
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "SasaCalculator.h"
#include "DeadEndElimination.h"
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
#include "ResidueSelection.h"
#include "versatileFunctions.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = ""; 
string programDescription = "Reads a PDB into MSL and then calculates the energy; made for reading rosetta PDBs and calculating the energy";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "13 Feb 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

double getSasaAtPosition(vector<string> _chainIds, int _position, SasaCalculator _sasa){
    double positionSasa = 0;
    for (uint i=0; i<_chainIds.size()-1; i++){
        // get the residue on the chain
        string chainResi = _chainIds[i] + ',' + MslTools::intToString(_position);
        // calculate the SASA of the mutated pdb
        double resiSasa = _sasa.getResidueSasa(chainResi);
        positionSasa += resiSasa;
    }
    return positionSasa;
}

	
int main(int argc, char *argv[]){
    // get the current working directory
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    string outputDir = string(cwd) + "/";

    // parse through the command line arguments
    OptionParser OP;
	OP.readArgv(argc, argv);
	OP.autoExtendOptions();
	// check if there is a config file
	if(OP.getString("config") != "NA"){
		OP.readFile(OP.getString("config"));
	}
    int thread = OP.getInt("thread"); // thread of the helix
	string pdbFile = OP.getString("pdbFile"); // input pdb file
    string topFile = OP.getString("topFile"); // topology file
    string parFile = OP.getString("parFile"); // parameter file
    string solvFile = OP.getString("solvFile"); // solvent file
    string backboneFile = OP.getString("backboneFile"); // backbone file
    string helicalAxisFile = OP.getString("helicalAxisFile"); // helical axis file

	// output the command line arguments
    cout << "Command line arguments:" << endl;
	cout << "thread: " << thread << endl;
	cout << "pdbFile: " << pdbFile << endl;
	cout << "topFile: " << topFile << endl;
	cout << "parFile: " << parFile << endl;
	cout << "solvFile: " << solvFile << endl;
	cout << "backboneFile: " << backboneFile << endl;
	cout << "helicalAxisFile: " << helicalAxisFile << endl;

	// read in the input pdb file
	System pdb;
	pdb.readPdb(pdbFile,true);
    pdb.buildAllAtoms();
	pdb.saveAltCoor("startState");

    // initialize the SASA calculator
	SasaCalculator sasa(pdb.getAtomPointers());
	sasa.calcSasa();

    // initialize the residue to be added
    Residue resi("ALA", 26);
    // TODO:
    // 1. loop through the pdb positions
    // 2. mutate each position on both helices to ala
    // 3. calculate the SASA of the mutated pdb; SASA of the pdb at that position before and after that change
    // 4. save that value to a map with the sequence as a key
    // 5. output the map to a csv file

    // loop through the pdb positions
    vector<Position*> positions = pdb.getPositions();
    for (uint i=0; i<positions.size(); i++){
        // get the position
        Position* pos = positions[i];
        // add the alanine identity to the position
        resi.setResidueNumber(i+23);
        pos->addIdentity(resi);
    }
    pdb.buildAllAtoms();
    cout << "Identity added to all positions" << endl;

    // get the chains from the pdb
    vector<Chain*> chains = pdb.getChains();
    // get a list of the chain ids
    vector<string> chainIds;
    for (uint i=0; i<chains.size(); i++){
        chainIds.push_back(chains[i]->getChainId());
    }
    map<string, map<string, double>> sequenceSasaMap;
    for (uint i=0; i<chains[0]->positionSize(); i++){
        // get the previous identity of the position
        Residue prevResi = positions[i]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        // check if the position is an alanine
        if (resi == "ALA"){
            cout << "Position " << i << " is an alanine" << endl;
            continue;
        }
        cout << "Position " << i << " is a " << resi << endl;
        map<string,double> sasaMap;
        double startSasa = getSasaAtPosition(chainIds, i+23, sasa);
        cout << "SASA: " << startSasa << endl;
        // loop through the chains
        for (uint j=0; j<chains.size(); j++){
            // get the position on the chain
            Position& pos = chains[j]->getPosition(i);
            cout << pos.toString() << endl;
            // switch the identity to alanine
            pos.setActiveIdentity(1);
            cout << pos.identitySize() << endl;
            cout << pos.toString() << endl;
        }
        // initialize the sasa for the position; make this a function and add this to before switching the aa
        double currentSasa = getSasaAtPosition(chainIds, i+23, sasa);

	    string currentSequence = extractSequence(pdb);
        sasaMap["Start"] = startSasa;
        sasaMap["Current"] = currentSasa;
        
        cout << "SASA: " << currentSasa << endl;
        sequenceSasaMap[currentSequence] = sasaMap;
        // get the SASA of the pdb at that position before and after that change
	    // write the pdb
	    writePdb(pdb, outputDir, "dimer");
        exit(0);

        // set the identity back to the previous identity
        //pos->setActiveIdentity(prevIdentity);
    }

	//string sequence = extractSequence(pdb);
 //   // make the polymer sequence
	//string polySeq = convertToPolymerSequence(sequence, thread);
 //   PolymerSequence PS(polySeq);
	//cout << PS << endl;

 //   // initialize system
 //   System sys;
	//CharmmSystemBuilder CSB(sys,topFile,parFile,solvFile);
	//// load the membrane as solvent
	//CSB.setSolvent("MEMBRANE");
	//// set the midpoint length of the membrane and the exponential factor for the membrane (src/CharmmEnergy.cpp: IMM1ZtransFunction) 
	//CSB.setIMM1Params(15, 10);
	//// sets all nonbonded interactions to 0, excluding interactions between far atoms (src/CharmmSystemBuilder.cpp: updateNonbonded)
	//CSB.setBuildNonBondedInteractions(false);
	//// Setup polymer sequence and build the sequence using CharmmSystemBuilder
	//if(!CSB.buildSystem(PS)) {
	//	cerr << "Unable to build system from " << PS << endl;
	//	exit(0);
	//}

	//// assign the coordinates of our system to the given geometry 
	//sys.assignCoordinates(pdb.getAtomPointers(),false);
	//sys.buildAllAtoms();

	// write the pdb
	writePdb(pdb, outputDir, "dimer");
}