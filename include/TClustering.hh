/*
 * TClustering.hh
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#ifndef TCLUSTERING_HH_
#define TCLUSTERING_HH_

using namespace std;

#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"

#include "TADCEventReader.hh"
#include "TCluster.hh"
#define N_DET_CHANNELS 256


class TClustering {

public:
	typedef  vector< vector<TCluster*> > vecvecTCluster;
	TClustering(int runNumber,int seedSigma=10,int hitSigma=7);
	virtual ~TClustering();
	void ClusterEvents(int nEvents);
private:
	void clusterEvent();
	void clusterPlane(int det);
	int combineCluster(int det,int ch);
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    vector<TCluster*> vecCluster;
    vecvecTCluster vecvecCluster;
    int nEvent;
    int seedSigma;
    int hitSigma;
    int verbosity;
    bool createClusterTree(int nEvents);
    void setBranchAdresses();
	stringstream  filepath;
    TTree *clusterTree;
    TFile *clusterFile;
    int runNumber;
    UInt_t nClusters[9];
};

#endif /* TCLUSTERING_HH_ */
