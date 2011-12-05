//
//  TAnalysisOfPedestal.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 30.11.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#ifndef TANALYSISOFPEDESTAL_HH_
#define TANALYSISOFPEDESTAL_HH_

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

#define N_DET_CHANNELS 256
using namespace std;

class TAnalysisOfPedestal {
public:
	TAnalysisOfPedestal(int runnumber,int seedSigma=10,int hitSigma=7);
	virtual ~TAnalysisOfPedestal();
	void	doAnalysis(int nEvents=0);
private:
	void saveHistos();
	void checkForDeadChannels();
	void analyseForSeeds();
	void getBiggestHit();
	void initialiseHistos();
	void checkForSaturatedChannels();
	void analyseCluster();
    void analyseBiggestHit();
	TH1F *hSaturatedChannels[8];
	TH1F *hSeedMap[8];
	TH1F *hSeedMap2[9];
	TH1F *hClusterMap[8];
	TH1F* hNumberOfSeeds[8];
	TH1F* hChannelBiggestHit[8];
	TH1F* hPulsHeightBiggestHit[8];
	TH1F* hPulsHeightNextBiggestHit[8];
	TH1F* hNumberOfClusters[9];
	TH1F* hClusterSize[9];
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    int nEvent;
    int seedSigma;
    int hitSigma;
    TH1F *histo_pulseheight_sigma[8];
	TH1F *histo_pulseheight_sigma_second[8];
	TH1F *histo_pulseheight_sigma125[8];
	TH1F *histo_second_biggest_hit_direction[8];
	TH1F *histo_pulseheight_sigma_second_left[8];
	TH1F *histo_pulseheight_sigma_second_right[8];
	TH1F *histo_biggest_hit_map[8];
	TH1F *histo_pulseheight_left_sigma[8];
	TH1F *histo_pulseheight_left_sigma_second[8];
	TH1F *histo_pulseheight_right_sigma[8];
	TH1F *histo_pulseheight_right_sigma_second[8];
};

#endif /* TANALYSISOFPEDESTAL_HH_ */

