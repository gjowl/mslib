#include <iostream>
#include "Atom.h"
#include "CharmmBondInteraction.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"

using namespace std;
using namespace MSL;

int main(){
        Atom a1("a1",3.14,6.62,6.02);
        Atom a2("a2",3.14,6.62,6.02);
        CharmmBondInteraction cbi(a1,a2,10.5,3.103);
        RandomNumberGenerator RNG1;
        RandomNumberGenerator RNG2;
        Transforms tr;
        double energy;
        double distance = 0;
        bool reach_left_bound = 0;
        bool reach_right_bound = 0;
        //double distance = (a1.getCoor() - a2.getCoor()).length();
        cout << "CharmmBondInteraction set with k = 3.21 d0 = 3.103" << endl;

        while(!(reach_left_bound && reach_right_bound)){
                double move = RNG1.getRandomDouble(0.1,0.3); 
                if((RNG2.getRandomDouble() < 0.5 && !reach_left_bound) || reach_right_bound)
                        move = -move;
                //cout << "the random move is:" << move << endl;
                tr.Xtranslate(a1,move);
                distance += move;
                if(distance > 15.1)
                        reach_right_bound = 1;
                if(distance < -15.1)
                        reach_left_bound = 1;

                cout << "The current distance between two atoms is: " << distance << endl;
                energy = cbi.getEnergy();

                cout << "The current energy of this interaction is: " << energy << endl ; 
        }
        return 0;
}
