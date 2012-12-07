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
#include "TResults.hh"

using namespace std;

#define N_DETECTORS 9
class TAnalysisOfPedestal {
public:
	TAnalysisOfPedestal(TSettings* settings);
	virtual ~TAnalysisOfPedestal();
	void setResults(TResults* results){this->res=results;};
	void	doAnalysis(UInt_t nEvents=0);
private:
	void analyseEvent();
	TResults *res;
	void updateMeanCalulation();
	void createPedestalMeanHistos();
	void saveHistos();
	void savePHinSigmaHistos();
	void findYPlotRangeForPHHisto(TH1F* histo, Float_t hitCut);
	void checkForDeadChannels();
	void analyseForSeeds();
	void getBiggestHit();
	void initialiseHistos();
	void checkForSaturatedChannels();
//	void analyseCluster();
    void analyseBiggestHit(bool CMN_corrected=false);
	TH1F *hSaturatedChannels[9];
	TH1F *hSeedMap[9];
	TH1F *hSeedMap2[9];
	TH1F *hClusterMap[9];
	TH1F* hNumberOfSeeds[9];
	TH1F* hChannelBiggestSignal[9];
	TH1F* hSNR_BiggestSignal[9];
	TH1F* hSNR_BiggestAdjacent[9];
	TH1F* hCMNoiseDistribution;
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings *settings;
	UInt_t verbosity;
    UInt_t nEvent;
    int seedSigma;
    int hitSigma;
    TH1F *hBiggestSignalInSigma[9];
    TH1F *hBiggestSignalInSigmaCMN[9];
	TH1F *hBiggestAdjacentSignalInSigma[9];
	TH1F *hBiggestAdjacentSignalInSigmaCMN[9];
	TH1F *histo_pulseheight_sigma125[9];
	TH1F *hHitOrderMap[9];
	TH1F *histo_pulseheight_sigma_second_left[9];
	TH1F *histo_pulseheight_sigma_second_right[9];
	TH1F *hBiggestHitChannelMap[9];
	TH1F *hBiggestHitChannelMapCMN[9];
	TH1F *histo_pulseheight_left_sigma[9];
	TH1F *histo_pulseheight_left_sigma_second[9];
	TH1F *histo_pulseheight_right_sigma[9];
	TH1F *histo_pulseheight_right_sigma_second[9];
	TH1F *hAllAdcNoise[9];
	TH1F *hDiaAllAdcNoise;
	TH1F *hDiaAllAdcNoiseCMN;

private:
	std::vector< std::vector<Float_t> > pedestalMeanValue,pedestalSigmaValue;
	std::vector< std::vector<UInt_t> > nPedestalHits;
	std::vector< Float_t > vecAvrgPed;
	std::vector< Float_t > vecAvrgPedCMN;
	std::vector< Float_t > vecAvrgSigma;
	std::vector< Float_t > vecAvrgSigmaCMN;
	std::vector< Float_t > vecCMNoise;
	std::vector< Float_t > vecEventNo;

	std::vector<Float_t>vecSaturatedChannels[N_DETECTORS];

	std::vector<Float_t>vecBiggestSignalInSigma[N_DETECTORS];
	std::vector<Float_t>vecBiggestSignal[N_DETECTORS];
	std::vector<Int_t> vecBiggestHitChannel[N_DETECTORS];
	std::vector<Float_t>vecBiggestAdjacentSignalInSigma[N_DETECTORS];
	std::vector<Float_t>vecBiggestAdjacentSignal[N_DETECTORS];
	std::vector<Int_t> vecBiggestAdjacentHitChannel[N_DETECTORS];

	std::vector<Float_t>vecBiggestSignalInSigmaCMN[N_DETECTORS];
	std::vector<Int_t> vecBiggestHitChannelCMN[N_DETECTORS];
	std::vector<Float_t>vecBiggestAdjacentSignalCMN[N_DETECTORS];
	std::vector<Int_t> vecBiggestAdjacentHitChannelCMN[N_DETECTORS];

	std::vector<Int_t> vecHitOrder[N_DETECTORS];
	std::vector< std::vector<UInt_t> > diaRawADCvalues; //vector of vector of adc Value (ch, eventNo)
	THTMLPedestal *htmlPedestal;
};

#endif /* TANALYSISOFPEDESTAL_HH_ */

