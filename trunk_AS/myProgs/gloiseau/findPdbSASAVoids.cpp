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

double getSasaAtPosition(System &_pdb, vector<string> _chainIds, int _position){
    double positionSasa = 0;
	SasaCalculator sasa(_pdb.getAtomPointers());
	sasa.calcSasa();
    for (uint i=0; i<_chainIds.size(); i++){
        // get the residue on the chain
        string chainResi = _chainIds[i] + ',' + MslTools::intToString(_position);
        // calculate the SASA of the mutated pdb
        double resiSasa = sasa.getResidueSasa(chainResi);
        positionSasa += resiSasa;
        cout << "chainResi: " << chainResi << "; resiSasa: " << resiSasa << endl;
        cout << sasa.getResidueSasaTable() << endl;
    }
    return positionSasa;
}

void setAminoAcidAtPosition(System &_pdb, vector<Chain*> _chains, int _position, int _chainPosition, string _aa){
    for (uint j=0; j<_chains.size(); j++){
        // get the position on the chain
        Position& pos = _chains[j]->getPosition(_position);
        string chain = _chains[j]->getChainId();
        string posId = chain+','+MslTools::intToString(_chainPosition);
        cout << "posId: " << posId << endl;
        // switch the identity to alanine
        _pdb.setActiveIdentity(posId,_aa);
    }
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
    vector<string> ids = OP.getStringVector("ids"); // chain ids
    vector<int> interfacePositions = OP.getIntVector("interfacePositions"); // interface positions

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
    CharmmSystemBuilder CSB(pdb,topFile,parFile,solvFile);
    CSB.setSolvent("MEMBRANE");
    CSB.setIMM1Params(15, 10);
    CSB.setBuildNonBondedInteractions(false);
    CSB.buildSystemFromPDB(pdbFile);

	PDBWriter writer;
	writer.open(outputDir + "/" +  "test.pdb");
	writer.write(pdb.getAtomPointers(), true, false, true);
    // save the starting state of the pdb
    pdb.saveAltCoor("start");
    // initialize the SASA calculator
	SasaCalculator sasa(pdb.getAtomPointers());
	sasa.calcSasa();
    
    // get the chains from the pdb
    vector<Chain*> chains = pdb.getChains();
    vector<Position*> positions = pdb.getPositions();

    string startPosition = positions[0]->getPositionId(1);
    // convert start position to a number
    int startPos = MslTools::toInt(startPosition);

    // get the chains from the pdb
    for (uint i=0; i<positions.size(); i++){
        // get the position
        Position* pos = positions[i];
        string posId = pos->getPositionId();
        // add identity to the position
        Residue prevResi = positions[i]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        //double resiSasa = sasa.getResidueSasa(posId);
        cout << "Position " << posId << " is " << resi << endl;
        CSB.addIdentity(posId,"ALA");
        pdb.setActiveIdentity(posId,"ALA");
        Residue currResi = positions[i]->getCurrentIdentity();
        string resi1 = currResi.getResidueName();
        cout << "Position " << posId << " is " << resi1 << endl;
        //resiSasa = sasa.getResidueSasa(posId);
        pdb.setActiveIdentity(posId,resi);
    }
    pdb.buildAllAtoms();

	
    // get the chain ids
    vector<string> chainIds;
    for (uint i=0; i<chains.size(); i++){
        chainIds.push_back(chains[i]->getChainId());
    }

    // get the SASA of mutated pdbs
    map<string, map<string, double>> sequenceSasaMap;
    for (uint i=0; i<chains[0]->positionSize(); i++){
        int pos = i;
        int chainPos = i+startPos;
        // get the previous identity of the position
        Residue prevResi = positions[pos]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        // check if the position is an alanine
        if (resi == "ALA"){
            cout << "Position " << pos << " is ALA" << endl;
            continue;
        }
        cout << "Position " << pos << " is " << resi << endl;
        cout << sasa.getResidueSasaTable() << endl;
        // initialize the sasa map
        map<string,double> sasaMap;
        // get the sasa of the position
        double startSasa = getSasaAtPosition(pdb, chainIds, chainPos);
        cout << "Start SASA: " << startSasa << endl;
	    sasa.calcSasa();
        //cout << pdb.getAtomPointers().toString() << endl;
	    double totalSasa = sasa.getTotalSasa();
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, "ALA");
        Residue currResi = positions[pos]->getCurrentIdentity();
        string resi1 = currResi.getResidueName();
        cout << "Position " << pos << " is " << resi1 << endl;
        pdb.buildAllAtoms();
	    sasa.calcSasa();
	    totalSasa = sasa.getTotalSasa();
        //cout << totalSasa << endl;
        //cout << sasa.getResidueSasaTable() << endl;
        double posSasa = positions[pos]->getSasa();
        //cout << "Position " << pos << " SASA: " << posSasa << endl;
        
        // initialize the sasa for the position; make this a function and add this to before switching the aa
        double currentSasa = getSasaAtPosition(pdb, chainIds, chainPos);
        cout << "Current SASA: " << currentSasa << endl;
        exit(0);

	    string currentSequence = extractSequence(pdb);
        sasaMap["Start"] = startSasa;
        sasaMap["Current"] = currentSasa;
        sasaMap["SasaDifference"] = startSasa - currentSasa;
        sequenceSasaMap[currentSequence] = sasaMap;

	    // write the pdb
	    writer.write(pdb.getAtomPointers(), true, false, true);
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, resi);
        pdb.applySavedCoor("start");
    }
    cout << sasa.getResidueSasaTable() << endl;
    exit(0);
    
    // TODO:
    // 1. loop through the pdb positions
    // 2. mutate each position on both helices to ala
    // 3. calculate the SASA of the mutated pdb; SASA of the pdb at that position before and after that change
    // 4. save that value to a map with the sequence as a key
    // 5. output the map to a csv file
    string sequence = extractSequence(pdb);
    // make the polymer sequence
    string polySeq = generateMultiIDPolymerSequence(sequence, thread, ids, interfacePositions);
    PolymerSequence PS(polySeq);
    System sys;
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

	writer.close();
}