//This program is designed to compare the energies of a wild-type PDB 
//structure and a mutant in a fixed-backbone simulation.  It will load
//a PDB file and calculate its energy, then add a new identity at a
//certain position, removing the old identity.  Next, the structure will
//be repacked and the lowest energy obtained.  The best structure will have
//its energy reported and the structure will be written to a PDB file.
//
//
