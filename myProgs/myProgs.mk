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

#MYSOURCE  = 
#MYPROGS   = moveProteinCenter
#MYHEADERS = 

include myProgs/yudong/yudong.mk
