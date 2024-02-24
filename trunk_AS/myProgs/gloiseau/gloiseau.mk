# THIS MAKE FILE IS USED TO MAKE PROGRAMS OR OBJECT THAT DO NOT BELONG 
# TO THE CENTRAL REPOSITIORY WITHOUT ALTERING THE MASTER MAKEFILE
#
# PROTOCOL FOR myProgs/USERNAME usage
#
# 1) Rename this file to myProgs.mk
#
# 2) Create a USERNAME directory
#
# 3) Place your own programs (i.e. prog1.cpp prog2.cpp) and objects (Obj1.h, Obj1.cpp) 
#    in the myProgs/USERNAME sub-directory 
#
# 4) Create USERNAME.mk in myProgs/USERNAME/
#
# 5) Add the program next to the MYPROG definition and the objects next to MYSOURCE in myProgs/USERNAME.mk
#
# 6) Add "MYDIR = myProgs/USERNAME" in myProgs/USERNAME.mk
#
# 7) Add "include myProgs/USERNAME/USERNAME.mk" to the end of this file
#
# To compile (from the MSL/trunk directory):
# Objects:
#   make objs/Obj1.o
#   make objs/Obj2.o
#
# Programs:
#   make bin/prog1
#   make bin/prog2

# filtering in make files: https://stackoverflow.com/questions/5674816/makefile-excluding-files

# OBJECTS
# The below gets the list of objects in the directory MYDIR (h and cpp files with the same name) and puts them in the variable MYSOURCE
# get all h files
objFileList = $(wildcard $(MYDIR)/*.h)
# convert the list of file paths to filenames without extension
objs = $(basename $(objFileList))
# add to global makefile variable MYSOURCE
MYSOURCE = $(notdir $(objs))
#
# PROGRAMS
# The below gets the list of programs in the directory MYDIR (files that don't match with h file names) and puts them in the variable MYPROGS
# get list of files in dir and get all files that don't match up with the file .h or .cpp, add to MYPROGS
progFileList = $(wildcard $(MYDIR)/*.cpp)
# convert the list of file paths to filenames without extension
progFiles = $(basename $(progFileList))
# filter the cpp file list for anything that doesn't have a .h file with the same name
progs = $(filter-out $(objs),$(progFiles))
#prog = $(filter-out $(addsuffix /%,$(MYSOURCE)),$(progFiles))
# add to global makefile variable MYPROGS
MYPROGS = $(notdir $(progs))
MYHEADERS = 