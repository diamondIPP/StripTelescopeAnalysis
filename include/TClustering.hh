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
#include "TEvent.hh"
#include "TPlane.hh"
#include "TSettings.class.hh"
//#define N_DET_CHANNELS 256
//#define N_DIA_CHANNELS 128

class TClustering {

public:
	TClustering(TSettings *settings);//int runNumber,int seedDetSigma=10,int hitDetSigma=7,int seedDiaSigma=5, int hitDiaSigma=3);
	virtual ~TClustering();
	void ClusterEvents(UInt_t nEvents);
	void setSettings(TSettings* settings);
	static TH1F* createEtaIntegral(TH1F* histo, std::string histName);
private:
	void addToEtaDistributions();
	void saveEtaCorrections();
	void clusterEvent();
	void clusterDetector(UInt_t det);
	int combineCluster(UInt_t det,UInt_t ch);
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings *settings;
    vector<TCluster> vecCluster[9];
    UInt_t clusterRev;
    UInt_t nEvent;
//    int seedDetSigma;
//    int hitDetSigma;
//    int seedDiaSigma;
//    int hitDiaSigma;
    int verbosity;
    bool createdTree;
    bool createClusterTree(int nEvents);
    void setBranchAdresses();
	stringstream  filepath;
	stringstream rawFilePath;
    TTree *clusterTree;
    TFile *clusterFile;
    UInt_t runNumber;
    UInt_t nClusters[9];
    UShort_t maxDetAdcValue;
    UShort_t maxDiaAdcValue;
    TEvent *pEvent;
    TH1F* hEtaDistribution[9];
    UInt_t nInvalidReadout;

};

#endif /* TCLUSTERING_HH_ */
