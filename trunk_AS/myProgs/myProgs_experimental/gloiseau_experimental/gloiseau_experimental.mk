# THIS MAKE FILE IS USED TO MAKE PROGRAMS OR OBJECT THAT DO NOT BELONG 
# TO THE CENTRAL REPOSITIORY WITHOUT ALTERING THE MASTER MAKEFILE
#
# PROTOCOL FOR myProgs/USERNAME_experimental usage: this is a copy of the original myProgs/USERNAME.mk file
# that I was playing with for a bit of my grad career. It worked well, but unfortunately only if you had a single
# program to build. It builds everything in the myProgs/USERNAME directory, but it compiles them all as a single
# program. I hated having to add things to a make file whenever I made a new program, and this works for that. But
# unfortunately I couldn't get it to build only a single file at a time. This led to issues with having the same 
# functions in multiple programs, and I ended up making a bunch of alternate versions of code that will be found in
# this myProgs/USERNAME_experimental directory. I'm leaving it here for posterity, but I transitioned my work to the myProgs/USERNAME.mk. 
# Files in this directory are accurate as of edits on 2024-2-23.
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
# 

# OBJECTS
# The below gets the list of objects in the directory MYDIR (h and cpp files with the same name) and puts them in the variable MYSOURCE
# get all h files
objFileList = $(wildcard $(MYDIR)/*.h)
# convert the list of file paths to filenames without extension
objs = $(basename $(objFileList))
# add to global makefile variable MYSOURCE
MYSOURCE = $(notdir $(objs))

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