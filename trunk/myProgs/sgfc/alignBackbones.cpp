#include <iostream>

#include "System.h"
#include "Transforms.h"
#include "AtomSelection.h"


using namespace std;
using namespace MSL;


/* 
This program will take two pdb files with identical alpha carbon backbones and align them.
*/

int main () {

	//Read pdb files into system

	string inputPdb1 = "myProgs/1C3W.pdb";
	string inputPdb2 = "myProgs/1FBB.pdb";


	System pdb1;
	System pdb2;

	Transforms tr;

	if (pdb1.readPdb(inputPdb1) && pdb2.readPdb(inputPdb2)) {
		cout << "Both pdb files read successfully" << endl;
		cout << "Pdb1: " << pdb1 << endl;
		cout << "Pdb2: " << pdb2 << endl;
	} else {
		cout << "There was a problem in reading one of the pdb files." << endl;
	return 0;
	}

	//Generate a selection of C-alpha carbons from each pdb

	AtomPointerVector pdb1Pointers = pdb1.getAtomPointers();
	AtomPointerVector pdb2Pointers = pdb2.getAtomPointers();

	AtomSelection pdb1Sel (pdb1Pointers);
	AtomSelection pdb2Sel (pdb2Pointers);

	AtomPointerVector pdb1Alpha = pdb1Sel.select("alpha_C, name ca and resi 5-156");
	AtomPointerVector pdb2Alpha = pdb2Sel.select("alpha_C, name ca and resi 5-156");

	cout << "PDB #1: " << endl << pdb1Pointers << endl;
	cout << "PDB #1 Alpha Carbons: " << endl << pdb1Alpha << endl;

	cout << "PDB #2: " << endl << pdb2Pointers << endl;
	cout << "PDB #2 Alpha Carbons: " << endl << pdb2Alpha << endl;

	cout << "Pre-align RMSD: " << pdb1Alpha.rmsd(pdb2Alpha) << endl;

	if ( tr.rmsdAlignment (pdb2Alpha, pdb1Alpha, pdb2Pointers) ){
		cout << "Pdb files 1 and 2 aligned!" << endl;
		cout << "Post-align RMSD: " << pdb1Alpha.rmsd(pdb2Alpha) << endl;
		Matrix rotMatrix = tr.getLastRotationMatrix();
		cout << "Rotation Matrix: " << rotMatrix << endl;
		CartesianPoint translationPoint = tr.getLastTranslation();
		cout << "Translation Point: " << translationPoint << endl;

		/*	
		tr.revertRmsdAlignment (pdb1Alpha, pdb2Alpha, pdb1Alpha);


		tr.rotate ( pdb2Pointers, rotMatrix, translationPoint );
		tr.translate ( pdb2Pointers, translationPoint );
		*/

	
		pdb1.writePdb("myProgs/pdb1-Align.pdb");
		pdb2.writePdb("myProgs/pdb2-Align.pdb");

	} else {
		cout << "Did not align..." << endl;
	}


	return 0;
}
