
# define compile command
ROOTCFLAGS  	:=  $(shell root-config --cflags) 
#-I$(ROOTSYS)/include
ROOTLIBS    	:= $(shell root-config --libs)
#-L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl
ROOTGLIBS   	:= $(shell root-config --glibs)
# -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lGui -pthread -Wl,-rpath,$(ROOTSYS)/lib -lm -ldl

SVNDEV 			:= -D'SVN_REV="$(shell svnversion -n .)"'

CC 				:= g++
SRCDIR          := src
INCLUDEDIR      := include

CFLAGS  		:= -g -Wall -I$(INCLUDEDIR) -D_REENTRANT 
CFLAGS      	+= $(SVNDEV) $(ROOTCFLAGS) -fPIC
LDFLAGS 		:= -L/usr/local/lib $(ROOTLIBS) $(ROOTGLIBS) -fPIC


OBJ 			:= diamondAnalysis.cpp
HEAD    		:= 





LD              := g++
LDFLAGS         := -g $(LLABLDFLAGS) $(ROOTLIBS)

LIBFILES		:=	HistogrammSaver.class.o ChannelScreen.o TDetectorPlane.o TDiamondTrack.o TDetectorAlignment.o 
LIBFILES		+=  FidCutRegion.o Cluster.class.o ClusteredEvent.class.o Clustering.class.o TDetector_Data.o TTrigger_Event.o
LIBFILES		+=  TPed_and_RMS.o TEvent_Array.o SlidingPedestal.class.o PSDetector.class.o PSEvent.class.o
LIBFILES		+=	RawEvent.class.o RawDetector.class.o Track.class.o AlignmentClass.o TADCEventReader.o
LIBFILES		+=	TSettings.class.o TRawEventReader.o TTransparentClustering.o TRawEventSaver.o TPedestalCalculation.o
LIBFILES		+=	TAnalysisOfClustering.o TAnalysisOfPedestal.o TTransparentAnalysis.o 
LIBFILES		+=  TSelectionClass.o
LIBFILES		+=  TAlignment.o TClustering.o
#LIBFILES 		+=	libTCluster.so 
LIBFILES		+=	libTPlane.so #libTEvent.so

PROGS			:= diamondAnalysis


all: rootclean diamondAnalysis

#all: $(OBJ) $(HEAD) makefile
#	$(CC) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJ) -o diamondAnalysis 
	
	
	
diamondAnalysis: $(LIBFILES)
		#
		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
        #
        
$(PROGS):
        #
        # linking $@
        #
		$(LD) $^ $(LDFLAGS)  $(ROOTGLIBS) $(OBJ) $(CFLAGS) -o $@
		@echo  "\n\nPlease do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib"

#libTCluster.so: TClusterDict.o TCluster.o 	
#		#
#		# Creating Shared ROOT Lib
#		#
#		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
#		#
#		g++  -g -fPIC -Wall -m64 -shared $(LDFLAGS) -o $@ $^
#		cp -rfv libTCluster.so ~/lib/ 
# 		#
# 		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
# 		#
 
TClusterDict.cpp: $(INCLUDEDIR)/TCluster.hh $(INCLUDEDIR)/TClusterLinkDef.h
		#
		# compiling $@
		#
		#echo $(ROOTSYS)/bin/rootcint -v $(CFLAGS) -f TClusterDict.cpp -c $(INCLUDEDIR)/TCluster.hh $(INCLUDEDIR)/LinkDef.h
		$(ROOTSYS)/bin/rootcint -v  -f TClusterDict.cpp  -c $(INCLUDEDIR)/TCluster.hh $(INCLUDEDIR)/TClusterLinkDef.h

TClusterDict.o: TClusterDict.cpp
		#
		#
		#
		g++ $(CFLAGS) -fPIC -c -m64 -o $@ $<

libTPlane.so: TPlaneDict.o TPlane.o TCluster.o TClusterDict.o
		#TClusterDict.o TCluster.o	
		#
		# Creating Shared ROOT Lib
		#
		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
		#
		g++  -g -fPIC -Wall -m64 -shared $(LDFLAGS) -o $@ $^
		cp -rfv libTPlane.so ~/lib/ 
 		#
 		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
 		#
 
TPlaneDict.cpp: $(INCLUDEDIR)/TPlane.hh $(INCLUDEDIR)/TPlaneLinkDef.h
		#
		# compiling $@
		#
		#echo $(ROOTSYS)/bin/rootcint -v  -f TPlaneDict.cpp  -c $(INCLUDEDIR)/TPlane.hh $(INCLUDEDIR)/TPlaneLinkDef.h
		$(ROOTSYS)/bin/rootcint -v -f TPlaneDict.cpp -c -p -I$(INCLUDEDIR) TCluster.hh TPlane.hh TPlaneLinkDef.h

TPlaneDict.o: TPlaneDict.cpp
		#
		#
		#
		g++ $(CFLAGS) -fPIC -c -m64 -o $@ $<
		

libTEvent.so: TEventDict.o TEvent.o 	
		#
		# Creating Shared ROOT Lib
		#
		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
		#
		g++  -g -fPIC -Wall -m64 -shared $(LDFLAGS) -o $@ $^
		cp -rfv libTEvent.so ~/lib/ 
 		#
 		# Please do: export LD_LIBRARY_PATH+=$LD_LIBRARY_PATH:~/lib
 		#
 
TEventDict.cpp: $(INCLUDEDIR)/TEvent.hh $(INCLUDEDIR)/TEventLinkDef.h
		#
		# compiling $@
		#
		#echo $(ROOTSYS)/bin/rootcint -v $(CFLAGS) -f TClusterDict.cpp -c $(INCLUDEDIR)/TCluster.hh $(INCLUDEDIR)/LinkDef.h
		$(ROOTSYS)/bin/rootcint -v  -f TEventDict.cpp  -c $(INCLUDEDIR)/TEvent.hh $(INCLUDEDIR)/TEventLinkDef.h

TEventDict.o: TEventDict.cpp
		#
		#
		#
		g++ $(CFLAGS) -fPIC -c -m64 -o $@ $<
				
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
	
rootclean:
	rm -fv *Dict.* && rm -fv *.so

