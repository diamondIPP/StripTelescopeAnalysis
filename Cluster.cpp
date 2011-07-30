//organize pedestal analysis code in a class structure (SlidingPedestal.class.cpp) and executable function (PedestalAnalyze.cpp)
//2010-07-11 Function wrapper for class SlidingPedestal

#include "Clustering.class.hh"

int main() {
   Clustering* cl = new Clustering();
   cl->ClusterRun();
   delete cl;
   return 0;
}
