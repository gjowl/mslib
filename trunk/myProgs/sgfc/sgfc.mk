# THIS MAKE FILE IS USED TO MAKE PROGRAMS OR OBJECT THAT DO NOT BELONG 
# TO THE CENTRAL REPOSITIORY WITHOUT ALTERING THE MASTER MAKEFILE
#
# 1) Rename this file to myProgs.mk
#
# 2) Place you own programs (i.e. prog1.cpp prog2.cpp) and objects (Obj1.h, Obj1.cpp) 
#    in the myProgs sub-directory 
#
# 3) Add the program next to the MYPROG definition and the objects next to MYSOURCE
#
# Example:
# for objects Obj1.h/Obj1.cpp and Obj2.h/Obj2.cpp
#
# MYSOURCE   = Obj1 Obj2
#
# for programs prog1.cpp prog2.cpp
#
# MYPROGS    = prog1 prog2
#
# To compile (from the above directory):
# Objects:
#   make objs/Obj1.o
#   make objs/Obj2.o
#
# Programs:
#   make bin/prog1
#   make bin/prog2

MYPROGS = testGetRandomOrder greedySimple greedySPM greedySPM_new randomTest residueMask coevolutionSimulation coevolutionSimulation_v_0.2 generateHelixBaseline mutateAndRepackHelixDimer genHomoUniverse genHomoUniverse_v_1.0.6 genHomoUniverse_v_1.0.7 calculateDimerEnergy mutateAndMCRepackHelixDimer addConservationBFactor modelBH3234Crosslinking CATM_v22
MYSOURCE = 
MYDIR = myProgs/sgfc

