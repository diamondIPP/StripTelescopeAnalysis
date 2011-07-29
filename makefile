# define compile command

CC = g++

CFLAGS  = -g -Wno-deprecated
LDFLAGS = -L/usr/local/lib

OBJ 	= diamondAnalysis.cpp
HEAD    = 

SVNDEV := -D'SVN_REV="$(shell svnversion -n .)"'

ROOTCFLAGS    = -D_REENTRANT -pthread  $(root-config --cflags) -I$(ROOTSYS)/include
ROOTLIBS      = $(root-config --libs) -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl
ROOTGLIBS     = $(root-config --glibs) -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lGui -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl

CFLAGS       += $(SVNDEV) $(ROOTCFLAGS)


all: $(OBJ) $(HEAD) makefile
	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJ) -o diamondAnalysis 

clean:	
	rm *.o

