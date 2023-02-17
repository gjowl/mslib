#include <sstream>
#include <iterator>
#include <unistd.h>
#include "HelixDimer.h"

using namespace std;
using namespace MSL;

HelixDimer::HelixDimer(string _name,double _energy, int _thread) : AtomContainer() {
	name = _name;
	energy = _energy;
	thread = _thread;
}

void HelixDimer::print() {
	cout << name << " " << energy << " " << hbonds.size() << " "  << thread << endl;
}

void HelixDimer::setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds, int _centroidIDNum) {
	xshift = _x;
	zshift = _z;
	axialrot = _a;
	crossang = _c;
	interface = _interface;
	centroidIDNum = _centroidIDNum;

	prolineMask = _prolineMask;
	hbonds = _hbonds;
	relAxialrot = 10.0/9.0 * axialrot + 200.0/27.0 * zshift;
	relZshift = 10.0/9.0 * zshift + axialrot/60.0;

	// set the 35 - thread - {0,1,4,5} as 2 - to mark the trapezoid
	int trapCorner = 35-thread;
	if(relAxialrot > 0) {
		trapCorner -= 1;
	} else if (relAxialrot < -100) {
		trapCorner += 1;
	}
	if(relZshift > 6) {
		trapCorner -= 4;
	} else if (relZshift < 0) {
		trapCorner += 4;
	}
	//cout << "TRAPCORNER " << trapCorner << endl;
	//cout << "interface " << interface << endl;
	if(trapCorner >= 0) {
		interface.replace(trapCorner,1,"2");
	}
	if(trapCorner >= 1) {
		interface.replace(trapCorner-1,1,"2");
	}
	if(trapCorner >= 4) {
		interface.replace(trapCorner-4,1,"2");
	}
	if(trapCorner >= 5) {
		interface.replace(trapCorner-5,1,"2");
	}

}

void HelixDimer::printHelixDimerDetails(ofstream & _fout) {
	_fout << "interfacial " << interface << endl;
	_fout << "prolineMask " << prolineMask << endl;
	_fout << "name " << name << endl;
	_fout << "thread " << thread << endl;
	_fout << "Xshift " << xshift << endl;
	_fout << "Zshift " << zshift << endl;
	_fout << "relZshift " << relZshift << endl;
	_fout << "axialRot " << axialrot << endl;
	_fout << "relAxialRot " << relAxialrot << endl;
	_fout << "crossAng " << crossang << endl;
	_fout << "energy " << energy << endl;
	string hbondsList = "";
	for(int h = 0; h < hbonds.size(); h++) {
		if(h == 0) {
			hbondsList += hbonds[h];
		} else {
			hbondsList += ";" + hbonds[h];
		}
	}
	_fout << "hBonds " << hbondsList << endl;
}