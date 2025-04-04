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
#include <string>

#include "AtomSelection.h"
#include "IcEntry.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "System.h"

using namespace MSL;
using namespace std;


int main(){


	string pdbtext = "\
ATOM      1  N   ALA A   1      -0.557   1.474 -12.560  1.00  0.00              \n\
ATOM      2  CA  ALA A   1      -1.686   1.490 -11.633  1.00  0.00              \n\
ATOM      3  CB  ALA A   1      -3.019   1.370 -12.406  1.00  0.00              \n\
ATOM      4  C   ALA A   1      -1.553   0.374 -10.582  1.00  0.00              \n\
ATOM      5  O   ALA A   1      -0.603  -0.407 -10.563  1.00  0.00              \n\
ATOM      6  N   ALA A   2      -2.542   0.281  -9.659  1.00  0.00              \n\
ATOM      7  CA  ALA A   2      -2.634  -0.689  -8.571  1.00  0.00              \n\
ATOM      8  CB  ALA A   2      -1.341  -0.604  -7.727  1.00  0.00              \n\
ATOM      9  C   ALA A   2      -3.884  -0.536  -7.686  1.00  0.00              \n\
TER      10      ALA A   2                                                      \n\
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

	PDBReader rAv;
	//rAv.open(argv[1]);
	rAv.open("/tmp/testPdb.pdb");
	rAv.read();
	AtomPointerVector av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();
	System sys(av);
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}


	/********************************************************
	 *  Create a moledule, a ALA ALA dipeptide, param 19 atoms
	 *
	 *   RESI ALA          0.00
	 *   GROUP   
	 *   ATOM N    NH1    -0.47  !     |
	 *   ATOM HN   H       0.31  !  HN-N
	 *   ATOM CA   CT1     0.07  !     |     HB1
	 *   ATOM HA   HB      0.09  !     |    /
	 *   GROUP                   !  HA-CA--CB-HB2
	 *   ATOM CB   CT3    -0.27  !     |    \
	 *   ATOM HB1  HA      0.09  !     |     HB3
	 *   ATOM HB2  HA      0.09  !   O=C
	 *   ATOM HB3  HA      0.09  !     |
	 *   GROUP                   !
	 *   ATOM C    C       0.51
	 *   ATOM O    O      -0.51
	 *   BOND CB CA  N  HN  N  CA  
	 *   BOND C  CA  C  +N  CA HA  CB HB1  CB HB2  CB HB3 
	 *   DOUBLE O  C 
	 *   IMPR N -C CA HN  C CA +N O   
	 *   DONOR HN N   
	 *   ACCEPTOR O C   
	 *   IC -C   CA   *N   HN    1.3551 126.4900  180.0000 115.4200  0.9996
	 *   IC -C   N    CA   C     1.3551 126.4900  180.0000 114.4400  1.5390
	 *   IC N    CA   C    +N    1.4592 114.4400  180.0000 116.8400  1.3558
	 *   IC +N   CA   *C   O     1.3558 116.8400  180.0000 122.5200  1.2297
	 *   IC CA   C    +N   +CA   1.5390 116.8400  180.0000 126.7700  1.4613
	 *   IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
	 *   IC N    C    *CA  HA    1.4592 114.4400 -120.4500 106.3900  1.0840
	 *   IC C    CA   CB   HB1   1.5390 111.0900  177.2500 109.6000  1.1109
	 *   IC HB1  CA   *CB  HB2   1.1109 109.6000  119.1300 111.0500  1.1119
	 *   IC HB1  CA   *CB  HB3   1.1109 109.6000 -119.5800 111.6100  1.1114
	 *
	 ********************************************************/

	sys.addIcEntry("A 0 C",  "A 1 N",  "A 1 CA", "A 1 C",  1.3551, 126.4900, 180.0000, 114.4400, 1.5390);
	sys.addIcEntry("A 1 N",  "A 1 CA", "A 1 C",  "A 2 N",  1.4592, 114.4400, 180.0000, 116.8400, 1.3558);
	sys.addIcEntry("A 2 N",  "A 1 CA", "A 1 C",  "A 1 O",  1.3558, 116.8400, 180.0000, 122.5200, 1.2297, true);
	sys.addIcEntry("A 1 CA", "A 1 C",  "A 2 N",  "A 2 CA", 1.5390, 116.8400, 180.0000, 126.7700, 1.4613);
	sys.addIcEntry("A 1 N",  "A 1 C",  "A 1 CA", "A 1 CB", 1.4592, 114.4400, 123.2300, 111.0900, 1.5461, true);
	sys.addIcEntry("A 1 C",  "A 2 N",  "A 2 CA", "A 2 C",  1.3551, 126.4900, 180.0000, 114.4400, 1.5390);
	sys.addIcEntry("A 2 N",  "A 2 CA", "A 2 C",  "A 3 N",  1.4592, 114.4400, 180.0000, 116.8400, 1.3558);
	sys.addIcEntry("A 2 N",  "A 2 C",  "A 2 CA", "A 2 CB", 1.4592, 114.4400, 123.2300, 111.0900, 1.5461, true);

	sys.fillIcFromCoor();
	sys.printIcTable();
	sys.wipeAllCoordinates();

	cout << "Test seeding with A 1 C, A 1 CA, A 1 N" << endl;
	sys.seed("A 1 C", "A 1 CA", "A 1 N");
	sys.buildAllAtoms();
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}
	string filename = "/tmp/builtAtoms-A_1_C-A_1_CA-A_1_N.pdb";
	PDBWriter writer(filename);
    writer.open();
	writer.write(sys.getAtomPointers());
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;


	cout << "Test seeding with A 1 N, A 1 CA, A 1 C" << endl;
	sys.wipeAllCoordinates();
	sys.seed("A 1 N", "A 1 CA", "A 1 C");
	sys.buildAllAtoms();
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}
	filename = "/tmp/builtAtoms-A_1_N-A_1_CA-A_1_C.pdb";
	writer.open(filename);
	writer.write(sys.getAtomPointers());
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;

	cout << "Test seeding with A 2 N, A 1 C, A 1 CA" << endl;
	sys.wipeAllCoordinates();
	sys.seed("A 2 N", "A 1 C", "A 1 CA");
	sys.buildAllAtoms();
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}
	filename = "/tmp/builtAtoms-A_2_N-A_1_C-A_1_CA.pdb";
	writer.open(filename);
	writer.write(sys.getAtomPointers());
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;

	cout << "Test seeding with A 1 C, A 1 CA, A 1 CB" << endl;
	sys.wipeAllCoordinates();
	sys.seed("A 1 C", "A 1 CA", "A 1 CB");
	sys.buildAllAtoms();
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}
	filename = "/tmp/builtAtoms-A_1_C-A_1_CA-A_1_CB.pdb";
	writer.open(filename);
	writer.write(sys.getAtomPointers());
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;



/*
	Atom N1("N" ,  -0.557,   1.474, -12.560);
	Atom CA1("CA",  -1.686,   1.490, -11.633);
	Atom CB1("CB",   0.000,   0.000,   0.000);
	Atom C1("C" ,  -1.553,   0.374, -10.582);
	Atom O1("O" ,   0.000,   0.000,   0.000);
	Atom N2("N" ,   0.000,   0.000,   0.000);
	Atom CA2("CA",   0.000,   0.000,   0.000);
	Atom CB2("CB",   0.000,   0.000,   0.000);
	Atom C2("C" ,   0.000,   0.000,   0.000);

	CB1.wipeCoordinates();
	O1.wipeCoordinates();
	N2.wipeCoordinates();
	CA2.wipeCoordinates();
	CB2.wipeCoordinates();
	C2.wipeCoordinates();

	N1.setResidueName("ALA");
        CA1.setResidueName("ALA");
        CB1.setResidueName("ALA");
        C1.setResidueName("ALA");
        O1.setResidueName("ALA");
        N2.setResidueName("ALA");
        CA2.setResidueName("ALA");
        CB2.setResidueName("ALA");
	C2.setResidueName("ALA");

	N1.setResidueNumber(1);
        CA1.setResidueNumber(1);
        CB1.setResidueNumber(1);
        C1.setResidueNumber(1);
        O1.setResidueNumber(1);
        N2.setResidueNumber(2);
        CA2.setResidueNumber(2);
        CB2.setResidueNumber(2);
	C2.setResidueNumber(2);

	vector<IcEntry*> IC;
	//   IC N    CA   C    +N    1.4592 114.4400  180.0000 116.8400  1.3558
	IC.push_back(new IcEntry(N1,   CA1,  C1 ,   N2 ,   1.4592, 114.4400,  180.0000, 116.8400,  1.3558));
	N1.addIcEntry(IC[0]);
	N2.addIcEntry(IC[0]);
	//   IC +N   CA   *C   O     1.3558 116.8400  180.0000 122.5200  1.2297
	IC.push_back(new IcEntry(N2,   CA1,  C1 ,   O1 ,   1.3558, 116.8400,  180.0000, 122.5200,  1.2297, true));
	N2.addIcEntry(IC[1]);
	O1.addIcEntry(IC[1]);
	//   IC CA   C    +N   +CA   1.5390 116.8400  180.0000 126.7700  1.4613
	IC.push_back(new IcEntry(CA1,  C1 ,  N2 ,   CA2,   1.5390, 116.8400,  180.0000, 126.7700,  1.4613));
	CA1.addIcEntry(IC[2]);
	CA2.addIcEntry(IC[2]);
	//   IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
	IC.push_back(new IcEntry(N1,   C1 ,  CA1,   CB1,   1.4592, 114.4400,  123.2300, 111.0900,  1.5461, true));
	N1.addIcEntry(IC[3]);
	CB1.addIcEntry(IC[3]);
	//   IC -C   N    CA   C     1.3551 126.4900  180.0000 114.4400  1.5390
	IC.push_back(new IcEntry(C1,   N2 ,  CA2,   C2 ,   1.3551, 126.4900,  180.0000, 114.4400,  1.5390));
	C1.addIcEntry(IC[4]);
	C2.addIcEntry(IC[4]);
	//   IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
	IC.push_back(new IcEntry(N2,   C2 ,  CA2,   CB2,   1.4592, 114.4400,  123.2300, 111.0900,  1.5461, true));
	N2.addIcEntry(IC[5]);
	CB2.addIcEntry(IC[5]);

	cout << "N1   " << N1.getCoor() << endl;
	cout << "CA1  " << CA1.getCoor() << endl;
	cout << "CB1  " << CB1.getCoor() << endl;
	cout << "C1   " << C1.getCoor() << endl;
	cout << "O1   " << O1.getCoor() << endl;
	cout << "N2   " << N2.getCoor() << endl;
	cout << "CA2  " << CA2.getCoor() << endl;
	cout << "CB2  " << CB2.getCoor() << endl;
	cout << "C2   " << C2.getCoor() << endl;

	cout << endl;
	cout << "=========================" << endl;
	cout << "Build CB2" << endl;
	cout << endl;

	CB2.buildFromIc();

	cout << "N1   " << N1.getCoor() << endl;
	cout << "CA1  " << CA1.getCoor() << endl;
	cout << "CB1  " << CB1.getCoor() << endl;
	cout << "C1   " << C1.getCoor() << endl;
	cout << "O1   " << O1.getCoor() << endl;
	cout << "N2   " << N2.getCoor() << endl;
	cout << "CA2  " << CA2.getCoor() << endl;
	cout << "CB2  " << CB2.getCoor() << endl;
	cout << "C2   " << C2.getCoor() << endl;

	cout << endl;
	cout << "=========================" << endl;
	cout << "Build CB1" << endl;
	cout << endl;

	CB1.buildFromIc();

	cout << "N1   " << N1.getCoor() << endl;
	cout << "CA1  " << CA1.getCoor() << endl;
	cout << "CB1  " << CB1.getCoor() << endl;
	cout << "C1   " << C1.getCoor() << endl;
	cout << "O1   " << O1.getCoor() << endl;
	cout << "N2   " << N2.getCoor() << endl;
	cout << "CA2  " << CA2.getCoor() << endl;
	cout << "CB2  " << CB2.getCoor() << endl;
	cout << "C2   " << C2.getCoor() << endl;

	cout << endl;
	cout << "=========================" << endl;
	cout << "Build O1" << endl;
	cout << endl;

	O1.buildFromIc();

	cout << "N1   " << N1.getCoor() << endl;
	cout << "CA1  " << CA1.getCoor() << endl;
	cout << "CB1  " << CB1.getCoor() << endl;
	cout << "C1   " << C1.getCoor() << endl;
	cout << "O1   " << O1.getCoor() << endl;
	cout << "N2   " << N2.getCoor() << endl;
	cout << "CA2  " << CA2.getCoor() << endl;
	cout << "CB2  " << CB2.getCoor() << endl;
	cout << "C2   " << C2.getCoor() << endl;

	cout << "=========================" << endl;
	cout << "Check consistency:" << endl;
	double val = N2.distance(C1);
	cout << "N2/C1 distance: " << val;
	if (val - 1.3558 < -0.0001 || val - 1.3558 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = N2.angle(C1, CA1);
	cout << "N2/C1/CA1 angle: " << val;
	if (val - 116.8400 < -0.0001 || val - 116.8400 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = N2.dihedral(C1, CA1, N1);
	while (val <= -180.0) {
		val += 360.0;
	}
	while (val > 180.0) {
		val += 360.0;
	}
	cout << "N2/C1/CA1/N1 dihedral: " << val;
	if (val - 180.0000 < -0.0001 || val - 180.0000 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}

	val = O1.distance(C1);
	cout << "O1/C1 distance: " << val;
	if (val - 1.2297 < -0.0001 || val - 1.2297 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = O1.angle(C1, CA1);
	cout << "O1/C1/CA1 angle: " << val;
	if (val - 122.5200 < -0.0001 || val - 122.5200 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = O1.dihedral(C1, CA1, N2);
	while (val <= -180.0) {
		val += 360.0;
	}
	while (val > 180.0) {
		val += 360.0;
	}
	cout << "O1/C1/CA1/N2 dihedral: " << val;
	if (val - 180.0000 < -0.0001 || val - 180.0000 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}

	val = CA2.distance(N2);
	cout << "CA2/N2 distance: " << val;
	if (val - 1.4613 < -0.0001 || val - 1.4613 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = CA2.angle(N2, C1);
	cout << "CA2/N2/C1 angle: " << val;
	if (val - 126.7700 < -0.0001 || val - 126.7700 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}
	val = CA2.dihedral(N2, C1, CA1);
	while (val <= -180.0) {
		val += 360.0;
	}
	while (val > 180.0) {
		val += 360.0;
	}
	cout << "CA2/N2/C1/CA1 dihedral: " << val;
	if (val - 180.0000 < -0.0001 || val - 180.0000 > 0.0001) {
		cout << " NOT OK" << endl;
	} else {
		cout << " OK" << endl;
	}


	AtomPointerVector av;	
	av.push_back(&N1);
	av.push_back(&CA1);
	av.push_back(&CB1);
	av.push_back(&C1);
	av.push_back(&O1);
	av.push_back(&N2);
	av.push_back(&CA2);
	av.push_back(&CB2);
	av.push_back(&C2);

	string filename = "/tmp/builtAtoms.pdb";
	PDBWriter writer(filename);
    writer.open();
	writer.write(av);
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;
*/
}
