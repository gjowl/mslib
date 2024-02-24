#include <iostream>
#include "System.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;

int main(int argc, char *argv[])
{
    // Read in a PDB file

	string pdbFile = argv[1];
	Transforms tr;
	CartesianPoint Zaxis(0.0,0.0,1.0);

	System sys;
	sys.readPdb(pdbFile);

	for (uint i = 0; i < sys.getChain("A").atomSize(); i++) {
		Atom& modifiedAtom = sys.getChain("A").getAtom(i);
		tr.Xtranslate(modifiedAtom, 12.0);
		tr.rotate(modifiedAtom, 180.0, Zaxis);
	}

    sys.writePdb("monomer_b.pdb");
}
