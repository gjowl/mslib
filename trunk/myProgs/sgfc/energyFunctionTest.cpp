#include <CharmmAngleInteraction.h>
#include <CharmmBondInteraction.h>
#include <CharmmElectrostaticInteraction.h>
#include <CharmmDihedralInteraction.h>
#include <CharmmVdwInteraction.h>
#include <Scwrl4HBondInteraction.h>
#include <Atom.h>
#include <MslTools.h>
#include <Transforms.h>
#include <math.h>

using namespace MSL;
using namespace std;

int main () {	
	Transforms tr;
	//set up atoms

	Atom atomA;
	atomA.setCoor (-3.906, -3.504, -12.660);

	Atom atomB;
	atomB.setCoor (-4.107, -3.524, -11.681);
	
	Atom atomC;
	atomC.setCoor (-3.571, -4.269, -14.486);

	Atom atomD;
	atomD.setCoor (-3.724, -4.228, -15.705);

	Atom atomE;
	atomE.setCoor (-3.595, -3.085, -16.409);



	//set up parameters
	//Bond parameters
	double kConstant = 20;
	double dConstant = 5;

	//vdw Parameters
	double rmin = 10;
	double eMin = 20;

	//angle Parameters
	double kTheta = 25;
	double rad = 0;

	//dihedral Parameters
	double kChi = 20;
	double mult = 2;
	
	//Electrostatic parameters
	atomA.setCharge(1.0);
	atomB.setCharge(-1.0);

	double dielectric = 1.0;
	double rescale = 1.0;
	bool Rdielectric = false;

	//Scwrl4 parameters
	double dist = atomA.distance(atomC);
	double ang = 3.14;
	double e1_dihe = -120*6.28/360;
	double e2_dihe = 120*6.28/360;
	double d0 = 3;
	double sigD = 2;
	double B = 20;
	double alphaMax = 2.0;
	double betaMax = 2.0;
	double scaling = 1.0;
	atomC.setCharge(0.55);
	atomD.setCharge(-0.60);
	atomE.setCharge(-0.1);



	//Bond Energy
	for (unsigned int i = 0; i < 50; i++) {
		atomB.setCoor (0,0,i);
		CharmmBondInteraction bond (atomA, atomB, kConstant, dConstant);

		cout << "Bond," << atomA.distance(atomB) << "," << bond.getEnergy() << endl;
	}

	//Vdw Energy
	for (unsigned int i = 0; i < 50; i++) {
		atomB.setCoor (0,0,i);
		CharmmVdwInteraction vdw (atomA, atomB, rmin, eMin);

		cout << "Vdw," << atomA.distance(atomB) << "," << vdw.getEnergy() << endl;
	}

	//Angle Energy
	for (unsigned int i = 0; i < 360; i++) {
		tr.setBondAngle (atomA, atomB, atomC, i);
		CharmmAngleInteraction angle (atomA, atomB, atomC, kTheta, rad);

		cout << "angle," << atomA.angle(atomB, atomC) << "," <<angle.getEnergy() << endl;
	}

	//Dihedral Energy
	for (unsigned int i = 0; i < 360; i++) {
		tr.setDihedral (atomA, atomB, atomC, atomD, i, true) ;
		CharmmDihedralInteraction dihedral (atomA, atomB, atomC, atomD, kChi, mult, rad);
/*
		cout << "Intended angle: " << i << endl;
		cout <<"A: "<<atomA.getCoor() << endl;
		cout <<"B: "<<atomB.getCoor() << endl;
		cout <<"C: "<<atomC.getCoor() << endl;
		cout <<"D: "<<atomD.getCoor() << endl;
		cout << "Dihedral: " << atomA.dihedral(atomB, atomC, atomD) << endl;

		cout << "==========================================" << endl;
*/	

		cout << "Dihedral," << atomA.dihedral(atomB, atomC, atomD) << "," << dihedral.getEnergy() << endl;

	}
	//Electrostatic Energy
	for (uint i = 0; i < 30; i++) {
		atomB.setCoor(1,1,i);
		CharmmElectrostaticInteraction electric (atomA, atomB, dielectric, rescale, Rdielectric);

		cout << "Electrostatic," << atomA.distance(atomB) << "," << electric.getEnergy() << endl;
	}
		





	//SCWRL Energy
	for (uint i = 0; i < 360; i ++) {
		tr.setBondAngle (atomB, atomC, atomD, i) ;
	
		Scwrl4HBondInteraction hBond ( atomA, atomB, atomC, atomD, atomE, dist, ang, e1_dihe, e2_dihe, d0, sigD, B, alphaMax, betaMax, scaling );
	
		cout << "HBond," << atomB.angle(atomC, atomD) << "," << hBond.getEnergy() << endl;

/*
//check Alpha
		CartesianPoint e0 = atomA.getCoor() - atomB.getCoor();
		CartesianPoint n = atomA.getCoor() - atomC.getCoor();
		double cosAlpha = cos(CartesianGeometry::angleRadians((n * -1.0), e0));

//check Beta

		CartesianPoint eA = (CartesianGeometry::buildRadians(atomC.getCoor(), atomD.getCoor(), atomE.getCoor(), dist, ang, e1_dihe)) - atomB.getCoor();
		double cosBeta_eA = cos(CartesianGeometry::angleRadians(n, eA));

		hBond.printParameters();
		cout << "Beta: " << cosBeta_eA << endl;
		cout << "Alpha: " << cosAlpha << endl;
*/
	}



		return 0;
}
