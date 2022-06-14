
#include <iostream>

#include "System.h"
#include "Transforms.h"
#include "AtomSelection.h"


using namespace std;
using namespace MSL;


/* 
This program will read a pdb file and copy it into two systems, original and transform.  The "transform" system will undergo some rotations and translations which will be written to a pdb file.  The transform system will then be realigned to the original system and a new pdb file will be written.
*/

int main () {

	//Read pdb files into system

	string inputPdb = "myProgs/1C3W.pdb";
	System pdbOriginal;
	System pdbTransform;
	Transforms tr;
	if (pdbOriginal.readPdb(inputPdb)) {
		cout << "Pdb file read successfully!" << endl;

		 //Make a copy of the system for manipulation
		pdbTransform = pdbOriginal;
		cout << pdbOriginal << endl;
		cout << pdbTransform << endl;
	} else {
		cout << "Did not read pdb file..." << endl;
		return 0;
	}

	AtomPointerVector pOriginal = pdbOriginal.getAtomPointers();
	AtomPointerVector pTransform = pdbTransform.getAtomPointers();

	tr.rotate(pTransform, 50, CartesianPoint(1.1, 3.2, -1.5), CartesianPoint(2.4, 5.5, 0.0));
	tr.translate(pTransform, CartesianPoint(50.0, 100.9, -100.5));
	pdbTransform.writePdb("myProgs/1C3W-trans-system.pdb"); 
	cout << "Rotated file written" << endl;


	cout << "Transformed RMSD: " << pTransform.rmsd(pOriginal) << endl;
	if (tr.rmsdAlignment(pTransform, pOriginal)) {
		cout << "Transformed pdb file realigned" << endl;
		cout << "New RMSD: " <<pTransform.rmsd(pOriginal) << endl;
		pdbTransform.writePdb("myProgs/1C3W-align-system.pdb");
	}

	return 0;
}
