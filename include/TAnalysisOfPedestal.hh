/**
 *  TAnalysisOfPedestal.hh
 *  Diamond Analysis
 *
 *
 *  Created by Lukas Baeni on 30.11.11.
 *
 *
 *  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
 *
 */


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
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTMLPedestal.hh"

#include "TADCEventReader.hh"
#include "TSettings.class.hh"

using namespace std;

class TAnalysisOfPedestal {
public:
	TAnalysisOfPedestal(TSettings* settings);
	virtual ~TAnalysisOfPedestal();
	void	doAnalysis(UInt_t nEvents=0);
private:
	void updateMeanCalulation();
	void createPedestalMeanHistos();
	void saveHistos();
	void checkForDeadChannels();
	void analyseForSeeds();
	void getBiggestHit();
	void initialiseHistos();
	void checkForSaturatedChannels();
	void analyseCluster();
    void analyseBiggestHit();
	TH1F *hSaturatedChannels[9];
	TH1F *hSeedMap[9];
	TH1F *hSeedMap2[9];
	TH1F *hClusterMap[9];
	TH1F* hNumberOfSeeds[9];
	TH1F* hChannelBiggestHit[9];
	TH1F* hPulsHeightBiggestHit[9];
	TH1F* hPulsHeightNextBiggestHit[9];
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings *settings;
    int nEvent;
    int seedSigma;
    int hitSigma;
    TH1F *histo_pulseheight_sigma[9];
	TH1F *histo_pulseheight_sigma_second[9];
	TH1F *histo_pulseheight_sigma125[9];
	TH1F *histo_second_biggest_hit_direction[9];
	TH1F *histo_pulseheight_sigma_second_left[9];
	TH1F *histo_pulseheight_sigma_second_right[9];
	TH1F *histo_biggest_hit_map[9];
	TH1F *histo_pulseheight_left_sigma[9];
	TH1F *histo_pulseheight_left_sigma_second[9];
	TH1F *histo_pulseheight_right_sigma[9];
	TH1F *histo_pulseheight_right_sigma_second[9];
	TH1F *hAllAdcNoise[9];
private:
	std::vector< std::vector<Float_t> > pedestalMeanValue,pedestalSigmaValue;
	std::vector< std::vector<UInt_t> > nPedestalHits;
	std::vector< std::vector<UInt_t> > diaRawADCvalues; //vector of vector of adc Value (ch, eventNo)
	THTMLPedestal *htmlPedestal;
};

#endif /* TANALYSISOFPEDESTAL_HH_ */

