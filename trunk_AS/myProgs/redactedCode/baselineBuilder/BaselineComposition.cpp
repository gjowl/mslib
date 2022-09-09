//Modified version of Rosetta AACompositionEnergy found in /exports/home/gloiseau/rosetta_src/main/source/src/core/energy_methods/AAComposition.cc
//

#include "BaselineComposition.h"

using namespace MSL;
using namespace std;

const string BaselineComposition::typeName = "AAComposition";



AACompositionEnergy::~AACompositionEnergy() = default;

void AACompositionEnergy::finalizeTotalEnergy(System &_sys){
//TODO:
// 1. initialize a vector of residue pointers
// 2. populate a vector of residue pointers for the system
// 3. write in a totalEnergy to set the AAComposition to using the calculateEnergy function below
}

double AACompositionEnergy::calculateEnergy(){
//TODO:
// 1. Get expected number residues for each type by multiplying frac by size of internal AAs
// 2. Loop through residues and count residues in each type set
// 3. loop through counts and accumulate the penalty
//
// The penalty reader comes from the AACompositionEnergySetup.cc file, where it actually reads the values in the comp file. I'm assuming that they did this because of the way that Rosetta already works (accumulates a bunch of functional scoring term file values as pointers that can be accessed easily (?), and this is why it's so messy looking. I should just be able to make my own .comp like files and read through to get the vector for this; long term, it would be good to add this type of functionality to this or another file, but for now, I'm just going to write it into my design code like I did with baselinePairInteractions, THEN call it here as a vector of vectors to be looped through
}


//There's a helper boject for this object in Rosetta that I'm not sure I need. It stores the following:
// - names of residues types to count
// - names of residues not to count
// - properties (?) necessary to be counted or not counted
//
//The first two I think are easy: read through a system, get the identities, and count only identities found
//And it looks like properties are types of residues, which so long as I'm getting identities it shouldn't matter because I'm planning on counting all of them
// - in the future this could be helpful if say I only want to count aromatics or things that can Hbond, but for now I don't think the properties are important
//
//So seems like I can rid of the helper function as well and just get to writing the rest of the calculate energy without those
//I'll have to learn how to run glycine in their AAComposition function or something because it may be difficult to properly compare this energy without some of these functionsi
//
//After another look through I DO need the types and maybe properties: This is what allows me to set caps on how many of say ALA can be found
//It calculates the difference between desired and observed number, then imposes a penalty as a function of that difference
// - so if outside of ideal range of AAs in the membrane by x AAs, penalty score will increase sharply as x increases
//
// I think this is way easier to code than I'm making it out to be. Now I just need to make sure that the vdW energy is also optimized for my sequences (it may have been, but the DEE part that I just got rid of was what made it so poor)
// All I have to do is make a penalty file that has TYP, start ot end, Penalties 100 0 100, Fraction for amount AA, and then run
//
// From what I can tell the person who coded this did a really poor job: GREAT Commenting, HORRIBLE execution of the code
//   - there are a lot of copies of the same function in the AACompositionEnergySetup.hh file that makes it really difficult to interpret what is happening
//   - I'm pretty sure I understand what they did (create a vector for each AA or type of penalty 0 and penalty 1 values, then have a SEPARATE method function for each)
//   - I think I'm going to write it in as an option to change the function to whatever from low to high value of penalty to make it simple, probably take values from whatever is used in the MC
