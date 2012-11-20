//organize pedestal analysis code in a class structure (SlidingPedestal.class.cpp) and executable function (PedestalAnalyze.cpp)
//2010-07-11 Function wrapper for class SlidingPedestal
//2010-07-21 Compiles now but runs at 1/10 speed with command: g++ -o PedestalAnalyze.exe PedestalAnalyze.cpp `root-config --cflags --glibs`

#include "SlidingPedestal.class.hh"
#include <sstream>
#include <iostream>

// argv[]={PedestalAnalyze, int RunNumber, int NEvents, int NSkip}
int main(int argc, char* argv[]) {

   int RunNumber = (int)strtod(argv[1],0);
   int NEvents = (int)strtod(argv[2],0);
   int NSkip = (int)strtod(argv[3],0);
   
   std::cout<<"RunNumber = "<<RunNumber<<"\tNEvents = "<<NEvents<<"\tNSkip = "<<NSkip<<std::endl;

   SlidingPedestal sp(RunNumber);
   sp.Slide(NEvents,NSkip);
   return 0;
}
