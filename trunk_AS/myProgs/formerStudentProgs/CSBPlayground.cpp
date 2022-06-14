#include <iostream>
#include "System.h"
#include "CharmmSystemBuilder.h"
using namespace MSL;
using namespace std;

int main(){
        System sys;
        sys.readPdb("/exports/home/yudong/2kpf.pdb");
        sys.writePdb("/exports/home/yudong/2kpf_msl.pdb");
        
        CharmmSystemBuilder CSB(sys,"/exports/home/yudong/mslib/trunk/toppar/charmm22.top","/exports/home/yudong/mslib/trunk/toppar/charmm22.par");
        
        PolymerSequence seq (sys.getAtomPointers());
        cout << seq << endl;
        seq.setPositionIdentity(0,5,0,"HSE");
        seq.setPositionIdentity(0,6,0,"HSE");
        seq.setPositionIdentity(1,5,0,"HSE");
        seq.setPositionIdentity(1,6,0,"HSE");
        cout << seq << endl;
         
        System singularitySys;
        CharmmSystemBuilder csb(singularitySys,"/exports/home/yudong/mslib/trunk/toppar/charmm22.top","/exports/home/yudong/mslib/trunk/toppar/charmm22.par");
        csb.buildSystem(seq);
        singularitySys.writePdb("/exports/home/yudong/2kpf_msl_big_band.pdb");
        
        singularitySys.seed();
        singularitySys.writePdb("/exports/home/yudong/2kpf_msl_after_seed.pdb");

        singularitySys.buildAllAtoms();
        singularitySys.writePdb("/exports/home/yudong/2kpf_msl_build_from_top.pdb");

        EnergySet *eset = singularitySys.getEnergySet();
        cout << "Energy From calcEnergies:EnergySet: " << eset->calcEnergy() << endl;
        cout << singularitySys.getEnergySummary();

        return 0;
} 
