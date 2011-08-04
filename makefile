
# define compile command
ROOTCFLAGS  	:= -D_REENTRANT -pthread  $(shell root-config --cflags) 
#-I$(ROOTSYS)/include
ROOTLIBS    	:= $(shell root-config --libs)
#-L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl
ROOTGLIBS   	:= $(shell root-config --glibs)
# -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lGui -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl

SVNDEV 			:= -D'SVN_REV="$(shell svnversion -n .)"'

CC 				:= g++
SRCDIR          := src
INCLUDEDIR      := include

CFLAGS  		:= -g -Wno-deprecated -I$(INCLUDEDIR)
CFLAGS      	+= $(SVNDEV) $(ROOTCFLAGS)
LDFLAGS 		:= -L/usr/local/lib

OBJ 			:= diamondAnalysis.cpp
HEAD    		:= 





LD              := g++
LDFLAGS         := -g $(LLABLDFLAGS) $(ROOTLIBS)

LIBFILES		:=	HistogrammSaver.class.o ChannelScreen.o TDetectorPlane.o TDiamondTrack.o TDetectorAlignment.o 
LIBFILES		+=  FidCutRegion.o Cluster.class.o ClusteredEvent.class.o Clustering.class.o TDetector_Data.o TTrigger_Event.o
LIBFILES		+=  TPed_and_RMS.o TEvent_Array.o SlidingPedestal.class.o PSDetector.class.o PSEvent.class.o
LIBFILES		+=	RawEvent.class.o RawDetector.class.o Track.class.o AlignmentClass.o TADCEventReader.o
LIBFILES		+=	TSettings.class.o

PROGS			:= diamondAnalysis


all: diamondAnalysis

#all: $(OBJ) $(HEAD) makefile
#	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJ) -o diamondAnalysis 
	
	
	
diamondAnalysis: $(LIBFILES)
        
$(PROGS):
        #
        # linking $@
        #
		$(LD) $^ $(LDFLAGS)  $(ROOTGLIBS) $(OBJ) $(CFLAGS) -o $@

%.o: $(SRCDIR)/%.cpp $(INCLUDEDIR)/%.hh
        #
        # compiling $@
        #
        #(cd $(SRCDIR); g++ $(CPPFLAGS) -c $< )
		g++ $(CFLAGS) -c $<
        # DONE
        #


clean:	
	rm -fv *.o diamondAnalysis

