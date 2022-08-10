/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#include <iostream>

#include "Transforms.h"
#include "AtomContainer.h"

using namespace std;

using namespace MSL;


int main() {





	string pdbtext = "\
ATOM      1 1H   ALA A   1      -0.079  -2.341  -3.811  1.00  0.00           H  \n\
ATOM      2 2H   ALA A   1       0.350  -0.712  -3.601  1.00  0.00           H  \n\
ATOM      3 3H   ALA A   1      -0.809  -1.422  -2.584  1.00  0.00           H  \n\
ATOM      4  N   ALA A   1       0.072  -1.586  -3.112  1.00  0.00           N  \n\
ATOM      5  CA  ALA A   1       1.121  -1.981  -2.192  1.00  0.00           C  \n\
ATOM      6  HA  ALA A   1       0.834  -2.871  -1.644  1.00  0.00           H  \n\
ATOM      7  CB  ALA A   1       2.417  -2.319  -2.965  1.00  0.00           C  \n\
ATOM      8 1HB  ALA A   1       2.218  -3.144  -3.682  1.00  0.00           H  \n\
ATOM      9 2HB  ALA A   1       2.775  -1.441  -3.546  1.00  0.00           H  \n\
ATOM     10 3HB  ALA A   1       3.228  -2.646  -2.281  1.00  0.00           H  \n\
ATOM     11  C   ALA A   1       1.370  -0.901  -1.154  1.00  0.00           C  \n\
ATOM     12  O   ALA A   1       1.549  -1.196   0.026  1.00  0.00           O  \n\
TER      13                                                                     \n\
END                                                                             \n";

	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Write a pdb file and read it into an atom vector     |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	
	// write the input test PDB file
	ofstream pdb_fs;
	pdb_fs.open("/tmp/testPdb.pdb");
	pdb_fs << pdbtext;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file /tmp/testPdb.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file /tmp/testPdb.pdb" << endl;
	}

	pdb_fs.close();

	Transforms tr;
	tr.setTransformAllCoors(true);
	

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                X-Rotate atom by atom with the             |" << endl;
	cout << "|                         old function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont1X;
	cont1X.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av1X = cont1X.getAtomPointers();
	cout << "Read atom vector with size " << av1X.size() << endl;
	if (av1X.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av1X.begin(); k != av1X.end() ; k++){
		cout << *(*k) << endl;
		tr.Xrotate(*(*k), 180);
	}
	cont1X.writePdb("/tmp/rotated1X.pdb");
	cout << "After 180 X rotation" << endl;
	for (AtomPointerVector::iterator k = av1X.begin(); k != av1X.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                X-Rotate atom by atom with the             |" << endl;
	cout << "|                         new function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont2X;
	cont2X.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av2X = cont2X.getAtomPointers();
	cout << "Read atom vector with size " << av2X.size() << endl;
	if (av2X.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av2X.begin(); k != av2X.end() ; k++){
		cout << *(*k) << endl;
		tr.Xrotate180(*(*k));
	}
	cont2X.writePdb("/tmp/rotated2X.pdb");
	cout << "After 180 X rotation" << endl;
	for (AtomPointerVector::iterator k = av2X.begin(); k != av2X.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                X-Rotate atom by atom vector               |" << endl;
	cout << "|                    with the new function                  |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont3X;
	cont3X.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av3X = cont3X.getAtomPointers();
	cout << "Read atom vector with size " << av3X.size() << endl;
	if (av3X.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av3X.begin(); k != av3X.end() ; k++){
		cout << *(*k) << endl;
	}
	tr.Xrotate180(av3X);
	cont3X.writePdb("/tmp/rotated3X.pdb");
	cout << "After 180 X rotation" << endl;
	for (AtomPointerVector::iterator k = av3X.begin(); k != av3X.end() ; k++){
		cout << *(*k) << endl;
	}












	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Y-Rotate atom by atom with the             |" << endl;
	cout << "|                         old function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont1Y;
	cont1Y.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av1Y = cont1Y.getAtomPointers();
	cout << "Read atom vector with size " << av1Y.size() << endl;
	if (av1Y.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av1Y.begin(); k != av1Y.end() ; k++){
		cout << *(*k) << endl;
		tr.Yrotate(*(*k), 180);
	}
	cont1Y.writePdb("/tmp/rotated1Y.pdb");
	cout << "After 180 Y rotation" << endl;
	for (AtomPointerVector::iterator k = av1Y.begin(); k != av1Y.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Y-Rotate atom by atom with the             |" << endl;
	cout << "|                         new function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont2Y;
	cont2Y.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av2Y = cont2Y.getAtomPointers();
	cout << "Read atom vector with size " << av2Y.size() << endl;
	if (av2Y.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av2Y.begin(); k != av2Y.end() ; k++){
		cout << *(*k) << endl;
		tr.Yrotate180(*(*k));
	}
	cont2Y.writePdb("/tmp/rotated2Y.pdb");
	cout << "After 180 Y rotation" << endl;
	for (AtomPointerVector::iterator k = av2Y.begin(); k != av2Y.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Y-Rotate atom by atom vector               |" << endl;
	cout << "|                    with the new function                  |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont3Y;
	cont3Y.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av3Y = cont3Y.getAtomPointers();
	cout << "Read atom vector with size " << av3Y.size() << endl;
	if (av3Y.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av3Y.begin(); k != av3Y.end() ; k++){
		cout << *(*k) << endl;
	}
	tr.Yrotate180(av3Y);
	cont3Y.writePdb("/tmp/rotated3Y.pdb");
	cout << "After 180 Y rotation" << endl;
	for (AtomPointerVector::iterator k = av3Y.begin(); k != av3Y.end() ; k++){
		cout << *(*k) << endl;
	}








	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Z-Rotate atom by atom with the             |" << endl;
	cout << "|                         old function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont1Z;
	cont1Z.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av1Z = cont1Z.getAtomPointers();
	cout << "Read atom vector with size " << av1Z.size() << endl;
	if (av1Z.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av1Z.begin(); k != av1Z.end() ; k++){
		cout << *(*k) << endl;
		tr.Zrotate(*(*k), 180);
	}
	cont1Z.writePdb("/tmp/rotated1Z.pdb");
	cout << "After 180 Z rotation" << endl;
	for (AtomPointerVector::iterator k = av1Z.begin(); k != av1Z.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Z-Rotate atom by atom with the             |" << endl;
	cout << "|                         new function                      |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont2Z;
	cont2Z.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av2Z = cont2Z.getAtomPointers();
	cout << "Read atom vector with size " << av2Z.size() << endl;
	if (av2Z.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av2Z.begin(); k != av2Z.end() ; k++){
		cout << *(*k) << endl;
		tr.Zrotate180(*(*k));
	}
	cont2Z.writePdb("/tmp/rotated2Z.pdb");
	cout << "After 180 Z rotation" << endl;
	for (AtomPointerVector::iterator k = av2Z.begin(); k != av2Z.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Z-Rotate atom by atom vector               |" << endl;
	cout << "|                    with the new function                  |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	AtomContainer cont3Z;
	cont3Z.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av3Z = cont3Z.getAtomPointers();
	cout << "Read atom vector with size " << av3Z.size() << endl;
	if (av3Z.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	for (AtomPointerVector::iterator k = av3Z.begin(); k != av3Z.end() ; k++){
		cout << *(*k) << endl;
	}
	tr.Zrotate180(av3Z);
	cont3Z.writePdb("/tmp/rotated3Z.pdb");
	cout << "After 180 Z rotation" << endl;
	for (AtomPointerVector::iterator k = av3Z.begin(); k != av3Z.end() ; k++){
		cout << *(*k) << endl;
	}



	tr.setStoreTransformHistory(true);
	CartesianPoint trans(6.4, 2.3, -1.2);
	CartesianPoint trans2(-3.1, 0.5, -3.4);
	CartesianPoint center(2.3, -1.7, 3.4);

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Apply history to a single atom             |" << endl;
	cout << "|                    with the old function                  |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	cout << endl;
	AtomContainer cont3H;
	cont3H.readPdb("/tmp/testPdb.pdb");
	Atom copy3H(cont3H[0]); // create a copy
	cout << cont3H[0] << endl;
	tr.Yrotate(cont3H[0], 180);
	tr.translate(cont3H[0], trans);
	tr.Xrotate(cont3H[0], 180);
	tr.Xrotate(cont3H[0], 73.4);
	tr.rotate(cont3H[0], 63.2, center);
	tr.Zrotate(cont3H[0], 180);
	tr.translate(cont3H[0], trans2);
	cout << "After transformations: original" << endl;
	cout << cont3H[0] << endl;
	tr.applyHistory(copy3H);
	cout << "Copy" << endl;
	cout << copy3H << endl;
	
	tr.resetHistory();


	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                Apply history to an atom vector            |" << endl;
	cout << "|                    with the old function                  |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	cout << endl;
	AtomContainer cont4H;
	cont4H.readPdb("/tmp/testPdb.pdb");
	AtomPointerVector av4H = cont3Z.getAtomPointers();
	for (AtomPointerVector::iterator k = av4H.begin(); k != av4H.end() ; k++){
		cout << *(*k) << endl;
	}
	AtomPointerVector copy4H;
	for (unsigned int i=0; i<av4H.size(); i++) {
		copy4H.push_back(new Atom(*av4H[i]));
	}
	tr.Yrotate(av4H, 180);
	tr.translate(av4H, trans);
	tr.Xrotate(av4H, 180);
	tr.Xrotate(av4H, 73.4);
	tr.rotate(av4H, 63.2, center);
	tr.Zrotate(av4H, 180);
	tr.translate(av4H, trans2);
	cout << "After transformations: original" << endl;
	for (AtomPointerVector::iterator k = av4H.begin(); k != av4H.end() ; k++){
		cout << *(*k) << endl;
	}
	tr.applyHistory(copy4H);
	cout << "Copy" << endl;
	for (AtomPointerVector::iterator k = copy4H.begin(); k != copy4H.end() ; k++){
		cout << *(*k) << endl;
		delete *k;
	}




	return 0;



}
