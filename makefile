# define compile command

CC = g++ -arch x86_64

CFLAGS  = -g -Wno-deprecated
LDFLAGS = -L/usr/local/lib

OBJ 	= main.cpp
HEAD    = 

ROOTCFLAGS    = -D_REENTRANT -pthread -I$(ROOTSYS)/include
ROOTLIBS      = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl
ROOTGLIBS     = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lGui -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl

CFLAGS       += $(ROOTCFLAGS)


all: $(OBJ) $(HEAD) makefile
	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJ) -o main 

clean:	
	rm *.o
