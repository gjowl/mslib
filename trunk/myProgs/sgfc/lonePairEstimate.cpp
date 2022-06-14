/* 
This program will identify the putative location of h-bonding lone pair electrons by reading in a system containing sp2 and sp3 oxygens and generating an Energy Set of Hydrogen bond interactions.  For these interactions, the lone pair electrons will be estimated using cartesian points which will be written to a pdb file.  

If the Scwrl 4 Hbond interaction is working correctly, the lone pair electrons of the sp2 oxygen should lie in the plane of the double bond and the sp3 lone pairs should be arranged in a tetrahedral alignment.  If it is not working, then the sp3 lone pairs will lie in a plane defined by the two atoms bonded to the oxygen just like the sp2 oxygen.

*/

#include "EnergySet.h"
#include "System.h"
#include "CartesianGeometry.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "PolymerSequence.h"
#include "SysEnv.h"
#include "Transforms.h"
#include "AtomSelection.h"

using namespace std;
using namespace MSL;


int main () {
	//SysEnv ENV;

	string topFile = "/exports/home/scondon/mslib/trunk/toppar/charmm22.top";
	string parFile = "/exports/home/scondon/mslib/trunk/toppar/charmm22.par";
	string solvFile = "/library/charmmTopPar/solvpar22.inp";
	string hBondFile = "/exports/home/scondon/mslib/trunk/toppar/scwrl4hb/par_hbond_CA_2.txt";



	//create a system
	Transforms tr;
	System sys;
	PolymerSequence seq ("A: SER THR B: THR SER");
	CharmmSystemBuilder charmm(sys, topFile, parFile);
	charmm.buildSystem(seq);
	HydrogenBondBuilder scwrl (sys, hBondFile);

	vector<Atom> lonePairs;

//	sys.printIcTable();
	
	//provide a starting seed to build the protein from its internal coordinates

	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
	cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}

	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}
	sys.buildAtoms();


	AtomSelection sele (sys.getAllAtomPointers());
	sele.select("chainB, chain B");

	tr.translate(sele.getSelection("chainB"), CartesianPoint(8.5,0,0));
	

	scwrl.buildInteractions(4.0);
	sys.writePdb("lonePairInitialSystem.pdb");

	


	//get energy for the system
	
	EnergySet *energy = sys.getEnergySet();
	energy->calcEnergy();
	energy->printSummary();


	//examine the H-bond energies
	energy->setAllTermsInactive();
	energy->setTermActive("SCWRL4_HBOND");
	energy->calcEnergy();
	energy->printSummary();

	//get each of the H-bond interactions


	std::map<std::string, std::vector<Interaction*> > *energyMap = energy->getEnergyTerms();
//	std::map<std::string, std::vector<Interaction*> >::iterator it;

	for (std::map<std::string, std::vector<Interaction*> >::iterator it=energyMap->begin(); it!=energyMap->end(); ++it) {
		if (it->first != "SCWRL4_HBOND") {
			// only look at HBond interactions
			continue;
		}


		for(vector<Interaction*>::const_iterator l=it->second.begin(); l!=it->second.end(); ++l) {
			vector<Atom*> &HBondPointers = (*l)->getAtomPointers();
			vector<double> HBondParams = (*l)->getParams();


			CartesianPoint e1 (0,0,0);
			Atom lonePair1 ("C,1,LNP,X"); //set arbitrary fake atom info to make finding them in pymol easy

			CartesianPoint e2 (0,0,0);
			Atom lonePair2 ("C,1,LNP,X");
			
			for (uint i = 0; i < HBondPointers.size(); i++) {
				cout << *HBondPointers[i] << endl;

			}
			for (uint i = 0; i < HBondParams.size(); i++) {
				cout << HBondParams[i] << endl;
			}
			
			e1 = (CartesianGeometry::buildRadians(HBondPointers[2]->getCoor(),HBondPointers[3]->getCoor(),HBondPointers[4]->getCoor(),HBondParams[0],HBondParams[1],HBondParams[2]));
			lonePair1.setCoor(e1);
			cout << lonePair1 << endl;

			e2 = (CartesianGeometry::buildRadians(HBondPointers[2]->getCoor(),HBondPointers[3]->getCoor(),HBondPointers[4]->getCoor(),HBondParams[0],HBondParams[1],HBondParams[3]));
			lonePair2.setCoor(e2);
			cout << lonePair2 << endl;

//add the acceptors as CartesianPoints to the system

			lonePairs.push_back(lonePair1);
			lonePairs.push_back(lonePair2);


			cout << "Energy: " << (*l)->getEnergy() << endl;
		}
	}



	//write the system to a PDB file
	for (uint i = 0; i < lonePairs.size(); i++) {
		cout << lonePairs[i].getCoor() << endl;
		lonePairs[i].setResidueNumber(i);//make each lone pair a different residue to aid pymol lookup
		sys.addAtom(lonePairs[i]);
	}

	sys.writePdb("lonePairsAdded.pdb");

	return 0;
}
