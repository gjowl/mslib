#include <iostream>
#include <typeinfo>
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "Transforms.h"
using namespace MSL;
using namespace std;

int main(){
        /*      
        System sys;
        string topFile = "/data03/CATM/topparsolvhbond/top_all22_prot.inp";
        string parFile = "/data03/CATM/topparsolvhbond/par_all22_prot.inp";
        string pdbFile = "/exports/home/yudong/2kpf-hse.pdb";

        CharmmSystemBuilder CSB(sys,topFile,parFile);
        // PDBReader PR;
        // if(!PR.open(pdbFile) || !PR.read()){
        //        cerr << "Errrrrrrrrrrror~1" << endl;
        // }
        
        if(!CSB.buildSystemFromPDB(pdbFile)){
                cerr << "Errrrrrrrrrrror~2" << endl;
                exit(0);
        }
        
        EnergySet * Eset = sys.getEnergySet();
        double energy = Eset->calcEnergy();
        cout << "Energy of the System is" << energy << endl;
        return 0;
        */

        System sys;
        sys.readPdb("/exports/home/yudong/2kpf.pdb");
        Residue ala = sys.getResidue(4);
        cout << "get Residue ALA "<< ala << endl;
        AtomPointerVector apv = ala.getAtomPointers();
        cout << "APV for this ALA is\n" << apv ;
        apv[0]->setType("NH1");
        apv[1]->setType("CT1");
        apv[2]->setType("C");
        apv[3]->setType("O");
        apv[4]->setType("CT3");
        apv[5]->setType("H"); apv[5]->setName("HN");
        apv[6]->setType("HB");
        apv[7]->setType("HA");
        apv[8]->setType("HA");
        apv[9]->setType("HA");
        cout << "APV after modify\n";

        for(int i=0;i<10;i++)
                cout << apv[i]->getName() << " "  << apv[i]->getType() << endl;
        
        CharmmSystemBuilder CSB(sys,"/exports/home/yudong/mslib/trunk/toppar/charmm22.top","/exports/home/yudong/mslib/trunk/toppar/charmm22.par");
        CSB.buildSystemFromPDB(apv);
        sys.buildAllAtoms();
        apv = sys.getAtomPointers();
        cout << "After Build CS singleResSys has the following atoms:\n" << sys.getAtomPointers() ;
        EnergySet * Eset = sys.getEnergySet();
        double energy = Eset->calcEnergy();
        cout << "Energy " << energy << endl;


        Transforms tr;
        //cout << *apv[7] << *apv[6] << *apv[1] << *apv[0];
        cout << "dihedral between HB1,CB-CA-N is " << apv[7]->dihedral(*(apv[6]),*(apv[1]),*(apv[0])) << endl;
        tr.setDihedral(*apv[7],*apv[6],*apv[1],*apv[0],(apv[7]->dihedral(*(apv[6]),*(apv[1]),*(apv[0]))+30));
         cout << "dihedral after transform between HB1,CB-CA-N is " << apv[7]->dihedral(*(apv[6]),*(apv[1]),*(apv[0])) << endl;

        //cout << "distance between C and CA is: " << apv[1]->distance(*(apv[10])) << endl;
        //tr.setBondDistance(*(apv[1]),*(apv[10]),(apv[1]->distance(*(apv[10]))+0.2));

        //cout << "new distance between C and CA is: " << apv[1]->distance(*(apv[10])) << endl;
        cout << "--------CORD AFTER SET NEW BOND LENGTH------------\n";
        cout << apv << sys.getAtomPointers();
        cout << "new energy after set the new C-CA bond length is: " << Eset->calcEnergy() << endl;
}
