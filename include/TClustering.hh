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
#include "TSettings.class.hh"
#define N_DET_CHANNELS 256
#define N_DIA_CHANNELS 128

class TClustering {

public:
	TClustering(int runNumber,int seedDetSigma=10,int hitDetSigma=7,int seedDiaSigma=5, int hitDiaSigma=3);
	virtual ~TClustering();
	void ClusterEvents(int nEvents);
	void setSettings(TSettings* settings);
private:
	void clusterEvent();
	void clusterPlane(int det);
	int combineCluster(int det,int ch,int maxAdcValue=255);
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings *settings;
    vector<TCluster> vecCluster[9];
    UInt_t clusterRev;
    TCluster::vecvecTCluster vecvecCluster;
    TCluster::vecvecTCluster* pVecvecCluster;
    int nEvent;
    int seedDetSigma;
    int hitDetSigma;
    int seedDiaSigma;
    int hitDiaSigma;
    int verbosity;
    bool createClusterTree(int nEvents);
    void setBranchAdresses();
	stringstream  filepath;
	stringstream rawFilePath;
    TTree *clusterTree;
    TFile *clusterFile;
    int runNumber;
    UInt_t nClusters[9];
    UShort_t maxDetAdcValue;
    UShort_t maxDiaAdcValue;

};

#endif /* TCLUSTERING_HH_ */
