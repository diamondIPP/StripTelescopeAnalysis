
# define compile command
ROOTCFLAGS  	:=  $(shell root-config --cflags) 
#-I$(ROOTSYS)/include
ROOTLIBS    	:= $(shell root-config --libs) -lMinuit
#-L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl
ROOTGLIBS   	:= $(shell root-config --glibs) -lMinuit
# -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lGui -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl

ROOTCINT		:= $(ROOTSYS)/bin/rootcint
SVNDEV 			:= -D'SVN_REV="$(shell svnversion -n .)"'

CC 				:= g++
SRCDIR          := src
INCLUDEDIR      := include
OBJDIR			:= obj
LIBDIR			:= ~/lib


OPTIMAZATIONFLAG :=	-O0 -g
CFLAGS  		:= -g -Wall -I$(INCLUDEDIR) -D_REENTRANT 
CFLAGS      	+= $(SVNDEV) $(ROOTCFLAGS) -fPIC $(OPTIMAZATIONFLAG)


OBJ 			:= diamondAnalysis.cpp
HEAD    		:= 





LD              := g++

LDFLAGS 		:= -L/usr/local/lib  $(ROOTGLIBS) -g $(LLABLDFLAGS) -fPIC -Wall -m64 $(OPTIMAZATIONFLAG)

LIBFILES		:=	HistogrammSaver.class.o  TDetectorPlane.o TDiamondTrack.o TPlaneProperties.o
LIBFILES		+=  FidCutRegion.o Cluster.class.o ClusteredEvent.class.o TDetector_Data.o TTrigger_Event.o
LIBFILES		+=  TPed_and_RMS.o TEvent_Array.o SlidingPedestal.class.o PSDetector.class.o PSEvent.class.o
LIBFILES		+=	RawEvent.class.o RawDetector.class.o Track.class.o TADCEventReader.o
LIBFILES		+=	TRawEventReader.o  TRawEventSaver.o TPedestalCalculation.o
LIBFILES		+=	TAnalysisOfClustering.o TAnalysisOfPedestal.o TTracking.o 
LIBFILES		+=	TTransparentAnalysis.o TAnalysisOfAlignment.o
#LIBFILES		+=  TTransparentClustering.o TTransparentAnalysis.o AlignmentClass.o Clustering.class.o 
LIBFILES		+=  TSelectionClass.o TPositionPrediction.o  THTMLGenerator.o 
LIBFILES		+=  TAlignment.o TClustering.o TTrack.o TResidual.o
LIBFILES 		+=  TChannelMapping.o TSettings.class.o ChannelScreen.o LandauGaussFit.o
LIBFILES		+=	libTEvent.so


ROOTLIBFILES	:=	TEventDict.o TEvent.o  TPlane.o  TCluster.o TDetectorAlignment.o
PROGS			:= diamondAnalysis


all: rootclean diamondAnalysis

#all: $(OBJ) $(HEAD) makefile
#	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJ) -o diamondAnalysis 
	
	
	
diamondAnalysis: $(LIBFILES)
		#
		# Please do: export LD_LIBRARY_PATH+=$$LD_LIBRARY_PATH:~/lib
        #
        
$(PROGS):
        #
        # linking $@
        #
		$(LD) $^ $(LDFLAGS)  $(ROOTGLIBS) $(OBJ) $(CFLAGS) -o $@
		@echo  "\n\nPlease do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib"


libTEvent.so: $(ROOTLIBFILES)	
		#
		# Creating Shared ROOT Lib
		#
		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
		#
		$(LD) -m64 -shared $(LDFLAGS) -o $@ $^
		cp -rfv libTEvent.so $(LIBDIR) 
 		#
 		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
 		#
 
TEventDict.cpp: $(INCLUDEDIR)/TEvent.hh $(INCLUDEDIR)/TEventLinkDef.h
		#
		# compiling $@
		#
		#echo $(ROOTCINT) -v $(CFLAGS) -f TClusterDict.cpp -c $(INCLUDEDIR)/TCluster.hh $(INCLUDEDIR)/LinkDef.h
		$(ROOTCINT) -v  -f TEventDict.cpp  -c -p -I$(INCLUDEDIR) TCluster.hh TPlane.hh TDetectorAlignment.hh TEvent.hh TEventLinkDef.h

TEventDict.o: TEventDict.cpp
		#
		#
		#
		$(CC) $(CFLAGS) -fPIC -c -m64 -o $@ $<
				
%.o: $(SRCDIR)/%.cpp $(INCLUDEDIR)/%.hh
        #
        # compiling $@
        #
        #(cd $(SRCDIR); g++ $(CPPFLAGS) -c $< )
		$(CC) $(CFLAGS) -c $<
        # DONE
        #

clean:	
	rm -fv *.o diamondAnalysis
	
rootclean:
	rm -fv *Dict.* && rm -fv *.so

