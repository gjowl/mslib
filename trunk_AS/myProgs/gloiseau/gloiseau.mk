# THIS MAKE FILE IS USED TO MAKE PROGRAMS OR OBJECT THAT DO NOT BELONG 
# TO THE CENTRAL REPOSITIORY WITHOUT ALTERING THE MASTER MAKEFILE
#
# UPDATED PROTOCOL FOR myProgs/USERNAME usage
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
# Example:
# for objects Obj1.h/Obj1.cpp and Obj2.h/Obj2.cpp
#
# MYSOURCE   = Obj1 Obj2
#
# for programs prog1.cpp prog2.cpp
#
# MYPROGS    = prog1 prog2
#
# To compile (from the MSL/trunk directory):
# Objects:
#   make objs/Obj1.o
#   make objs/Obj2.o
#
# Programs:
#   make bin/prog1
#   make bin/prog2

#TODO: add in a description of what I did in this make file
# filtering in make files: https://stackoverflow.com/questions/5674816/makefile-excluding-files

# TODO: maybe to make more elegant, have files for each of these programs that when running the input,
# it will go to those files and get the objects and programs from there?

# The below works for now; now need to see if it works when I have other programs in here
# cannot have cpp and h files with the same name as the program!
# get list of files in dir, get all h files, then get all cpp files with same name, then add to MYSOURCE
# get all h files
objFileList = $(wildcard $(MYDIR)/*.h)
# convert the list of file paths to filenames without extension
objs = $(basename $(objFileList))
# add to global makefile variable MYSOURCE
MYSOURCE = $(notdir $(objs))

# get list of files in dir and get all files that don't match up with the file .h or .cpp, add to MYPROGS
# output list of files in the objDir to the terminal
progFileList = $(wildcard $(MYDIR)/*.cpp)
# convert the list of file paths to filenames without extension
progFiles = $(basename $(progFileList))
# filter the cpp file list for anything that doesn't have a .h file with the same name
progs = $(filter-out $(objs),$(progFiles))
#prog = $(filter-out $(addsuffix /%,$(MYSOURCE)),$(progFiles))
# add to global makefile variable MYPROGS
MYPROGS = $(notdir $(progs))
MYHEADERS = 

#2022-06-17: Added in input directories for objects and programs for some separation
# didn't want to have to change the includes for many files and objects, so I'm
# switching back to the orignal method of including the make files, but
# just going through all files in the directory rather than typing them in
#MYOBJDIR  = $(MYDIR)/objs
#MYPROGDIR = $(MYDIR)/progs

# To output to the command line for testing:
#$(info $$myprogs is [${MYPROGS}])
#$(info $$myprogs is [${MYSOURCE}])

## get the list of files in the objDir
#objs = $(wildcard $(MYOBJDIR)/*)
## convert the list of file paths to filenames without extension
#objFiles = $(basename $(objs))
## add to global makefile variable MYSOURCE
#MYSOURCE = $(notdir $(objFiles))
#$(info $$myprogs is [${MYSOURCE}])
#$(info $$myprogs is [${MYOBJDIR}])