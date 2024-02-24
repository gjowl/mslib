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

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "getPolyAlaSasa"; 
string programDescription = "Reads a PDB into MSL and then calculates the energy; made for reading rosetta PDBs and calculating the energy";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "13 Feb 2023";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

string extractSequence(System &_sys);
double getSasaAtPosition(System &_pdb, vector<string> _chainIds, int _position, double &_totalSasa);
map<string, map<string, double>> getMonomerSasa(System &_pdb, string _topFile, string _parFile, string _solvFile);
void setAminoAcidAtPosition(System &_pdb, vector<Chain*> _chains, int _position, int _chainPosition, string _aa);
void setupDirectory(string &_outputDir);

int main(int argc, char *argv[]){
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
    string outputDir = OP.getString("outputDir"); // output directory
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
	
    // setup the output directory
    setupDirectory(outputDir);
    
    // read in the input pdb file
    System startGeom;
    startGeom.readPdb(pdbFile);

	// read in the input pdb file
    System pdb;
    CharmmSystemBuilder CSB(pdb,topFile,parFile,solvFile);
    CSB.setSolvent("MEMBRANE");
    CSB.setIMM1Params(15, 10);
    CSB.setBuildNonBondedInteractions(false);
    CSB.buildSystemFromPDB(pdbFile);

    // save the starting state of the pdb (already repacked, don't need to load energy terms for another repack)
    pdb.saveAltCoor("start");

    // get the chains from the pdb
    vector<Chain*> chains = pdb.getChains();
    vector<Position*> positions = pdb.getPositions();

    // get the start position and convert to a number
    string startPosition = positions[0]->getPositionId(1);
    int startPos = MslTools::toInt(startPosition);

    // get the chains from the pdb
    for (uint i=0; i<positions.size(); i++){
        // get the position
        Position* pos = positions[i];
        string posId = pos->getPositionId();
        // add identity to the position
        Residue prevResi = positions[i]->getCurrentIdentity();
        string resi = prevResi.getResidueName();
        
        // if the first or last position of a chain, remove the identity (saw issues with unbuilt atoms in the center of the pdb at 0,0,0)
        if (i == 0 || i == chains[0]->positionSize()-1 || i == chains[0]->positionSize() || i == positions.size()-1){
            CSB.removeIdentity(posId,resi);
            CSB.addIdentity(posId,"ALA");
            pdb.setActiveIdentity(posId,"ALA");
        } else {
            CSB.addIdentity(posId,"ALA");
            pdb.setActiveIdentity(posId,"ALA");
        }
    }
    pdb.assignCoordinates(startGeom.getAtomPointers(), true);
    pdb.buildAllAtoms();
	string startSequence = extractSequence(pdb);
	PDBWriter writer;
	writer.open(outputDir + "/" + startSequence + "_voids.pdb");
	writer.write(pdb.getAtomPointers(), true, false, true);

    // get the chain ids
    vector<string> chainIds;
    for (uint i=0; i<chains.size(); i++){
        chainIds.push_back(chains[i]->getChainId());
    }

    // get the SASA of mutated pdbs
    map<string, map<string, double>> monomerSasas = getMonomerSasa(pdb, topFile, parFile, solvFile);

    // write the sasa map to a file
    ofstream sasaFile;
    sasaFile.open(outputDir + "/" + "sasaMap.txt");
    sasaFile << "" << endl;
    for (auto it=monomerSasas.begin(); it!=monomerSasas.end(); it++){
        sasaFile << it->first << ",";
        for (auto it2=it->second.begin(); it2!=it->second.end(); it2++){
            sasaFile << it2->second << ",";
        }
        sasaFile << endl;
    }
    sasaFile.close();

    // close the pdb writer
	writer.close();
}

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

map<string, map<string, double>> getMonomerSasa(System &_pdb, string _topFile, string _parFile, string _solvFile){
	Chain & inputChain = _pdb.getChain(0);

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _topFile, _parFile, _solvFile);
	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());
    
    for (uint i=0; i<inputChain.positionSize(); i++){
        // get the position
        Position& pos = inputChain.getPosition(i);
        string posId = pos.getPositionId();
        // add identity to the position
        Residue prevResi = pos.getCurrentIdentity();
        string resi = prevResi.getResidueName();
        if (resi == "ALA"){
            continue;
        } else {
            CSBMono.addIdentity(posId,"ALA");
            monoSys.setActiveIdentity(posId,"ALA");
        }
    }
    monoSys.assignCoordinates(inputChain.getAtomPointers());
    monoSys.buildAllAtoms();

    // save the coordinates of the monomer
    monoSys.saveAltCoor("start");

    // get the initial sasa of the monomer
    SasaCalculator startSasa(monoSys.getAtomPointers());
    startSasa.calcSasa();
    double startTotalSasa = startSasa.getTotalSasa()*2;

    // get the sequence 
    string startSeq = extractSequence(monoSys);

    // initialize the map to store the SASA
    map<string, map<string, double>> sasaMap;
    vector<string> chainIds;
    // get the chain from the system
    Chain & monoChain = monoSys.getChain(0);
    chainIds.push_back(monoChain.getChainId());
    // loop through the identities at each position and get the SASA
    for (uint i=0; i<monoChain.positionSize(); i++){
        // get the position on the chain
        Position& pos = monoChain.getPosition(i);
        Residue currResi = pos.getCurrentIdentity();
        string resi = currResi.getResidueName();
        string posId = pos.getPositionId(1);

        int chainPos = MslTools::toInt(posId);
        // get the sasa of the position
        double startTotalSasa = 0;
        double sasa = getSasaAtPosition(monoSys, chainIds, chainPos, startTotalSasa);
        sasaMap[startSeq]["Pos_"+posId] = sasa;
    }
	PDBWriter writer;
	writer.open("voids.pdb");
	writer.write(monoSys.getAtomPointers(), true, false, true);
    writer.close();
    return sasaMap;
}

double getSasaAtPosition(System &_pdb, vector<string> _chainIds, int _position, double &_totalSasa){
    double positionSasa = 0;
	SasaCalculator sasa(_pdb.getAtomPointers());
	sasa.calcSasa();
    for (uint i=0; i<_chainIds.size(); i++){
        // get the residue on the chain
        string chainResi = _chainIds[i] + ',' + MslTools::intToString(_position);
        // calculate the SASA of the mutated pdb
        double resiSasa = sasa.getResidueSasa(chainResi);
        positionSasa += resiSasa;
    }
    // calculate the total SASA of the mutated pdb
    _totalSasa = sasa.getTotalSasa();
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

void setupDirectory(string &_outputDir){
	_outputDir = string(get_current_dir_name()) + "/" + _outputDir;
	//_opt.outputDir = "/exports/home/gloiseau/mslib/trunk_AS/design_" + _opt.runNumber;
	string cmd = "mkdir -p " + _outputDir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}