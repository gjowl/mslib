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
        }
    }

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
    // loop through the identities at each position and get the SASA
    for (uint i=0; i< inputChain.positionSize(); i++){
        // get the position on the chain
        Position& pos = inputChain.getPosition(i);
        Residue currResi = pos.getCurrentIdentity();
        string resi = currResi.getResidueName();
        string posId = pos.getPositionId();

        // check if the position is an alanine
        if (resi == "ALA"){
            continue;
        }

        // set the identity to alanine
        monoSys.setActiveIdentity(posId,"ALA");
        cout << posId << endl;
        
        string mutantSeq = extractSequence(monoSys);
    
	    //Setup SasaCalculator to calculate the monomer SASA
	    SasaCalculator monoSasa(monoSys.getAtomPointers());
	    monoSasa.calcSasa();
	    double monomerSasa = monoSasa.getTotalSasa();
	    double totalMonomerSasa = monomerSasa*2;

        // add values to the sasaMap
        sasaMap[mutantSeq]["Mut_MonomerSasa"] = totalMonomerSasa;
        sasaMap[mutantSeq]["WT_MonomerSasa"] = startTotalSasa;

        // reset the monomer to the original identity and coordinates
        monoSys.setActiveIdentity(posId,resi);
        monoSys.applySavedCoor("start");
    }
    return sasaMap;
}

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
    System pdb;
    CharmmSystemBuilder CSB(pdb,topFile,parFile,solvFile);
    CSB.setSolvent("MEMBRANE");
    CSB.setIMM1Params(15, 10);
    CSB.setBuildNonBondedInteractions(false);
    CSB.buildSystemFromPDB(pdbFile);

	string startSequence = extractSequence(pdb);
	PDBWriter writer;
	writer.open(outputDir + "/" + startSequence + "_voids.pdb");
	writer.write(pdb.getAtomPointers(), true, false, true);

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
        CSB.addIdentity(posId,"ALA");
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
        // initialize the sasa map
        map<string,double> sasaMap;
        // get the sasa of the position
        double startTotalSasa = 0;
        double startSasa = getSasaAtPosition(pdb, chainIds, chainPos, startTotalSasa);
        cout << "Start SASA: " << startSasa << endl;

        // switch the identity to alanine
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, "ALA");
        Residue currResi = positions[pos]->getCurrentIdentity();
        string resi1 = currResi.getResidueName();
        cout << "Position " << pos << " is " << resi1 << endl;
        double posSasa = positions[pos]->getSasa();
        
        // initialize the sasa for the position; make this a function and add this to before switching the aa
        double mutantTotalSasa = 0;
        double currentSasa = getSasaAtPosition(pdb, chainIds, chainPos, mutantTotalSasa);
        cout << "Mutant SASA: " << currentSasa << endl;

	    string currentSequence = extractSequence(pdb);
        sasaMap["Start"] = startSasa;
        sasaMap["Total"] = startTotalSasa;
        sasaMap["Mutant"] = currentSasa;
        sasaMap["TotalMutant"] = mutantTotalSasa;
        sequenceSasaMap[currentSequence] = sasaMap;

	    // write the pdb
	    writer.write(pdb.getAtomPointers(), true, false, true);
        // set the amino acid back to the original
        setAminoAcidAtPosition(pdb, chains, pos, chainPos, resi);
        // reset the pdb
        pdb.applySavedCoor("start");
    }
    map<string, map<string, double>> monomerSasas = getMonomerSasa(pdb, topFile, parFile, solvFile);
    // append the monomer sasas to the sequence sasa map
    for (auto it=monomerSasas.begin(); it!=monomerSasas.end(); it++){
        sequenceSasaMap[it->first]["WT_MonomerSasa"] = it->second["WT_MonomerSasa"];
        sequenceSasaMap[it->first]["Mut_MonomerSasa"] = it->second["Mut_MonomerSasa"];
    }

    // write the sasa map to a file
    ofstream sasaFile;
    sasaFile.open(outputDir + "/" + "sasaMap.txt");
    sasaFile << "Sequence,Mutant_MonomerSasa,Mutant_AA,SasaDifference,WT_AA,WT_Sasa,Mutant_Sasa,WT_MonomerSasa" << endl;
    for (auto it=sequenceSasaMap.begin(); it!=sequenceSasaMap.end(); it++){
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