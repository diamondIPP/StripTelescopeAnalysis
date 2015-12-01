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
#include <algorithm>
#include <set>
#include <deque>

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "TMultiGraph.h"
#include "TPolyMarker.h"

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

	void initialiseHistos();
	void saveHistos();
	void saveAdcVsEventProfiles();
	TH1F* doSlidingWindowAnalysis(TH1D* histo,Int_t nAvrg,bool absolute=true);
	void savePHinSigmaHistos();
	Float_t findYPlotRangeForPHHisto(TH1F* histo, Float_t hitCut);
	void createPedestalMeanHistos();
private:
	void analyseEvent();
	void updateMeanCalulation(UInt_t det,UInt_t ch);
	void checkForDeadChannels(UInt_t det,UInt_t ch);
	void checkCommonModeNoise();
    void analyseBiggestHit(UInt_t det,bool CMN_corrected=false);
    void findBiggestSignalInDet(UInt_t det, UInt_t ch);
    void findBiggestSignalInDia(UInt_t pattern);//, UInt_t ch);
    void SetYRangeForSignalInSigmaPlot(TH1F* histo);
private:
    map<Int_t,TH1F*> hBiggestSignalInSigmaDiaPattern;
    map<Int_t,TH1F*> hBiggestSignalInSigmaDiaPatternCMN;
    map<Int_t,TH1F*> hBiggestAdjacentSignalInSigmaDiaPattern;
    map<Int_t,TH1F*> hBiggestAdjacentSignalInSigmaDiaPatternCMN;
	TResults *res;
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
	TH1F* hRelCmnUncertainty;
	TH1F* hCmnNUsedChannels;
	TH1F* hCmnUsedChannels;
	TH1F* hNewComonModeNoise;
	TH1F* hCmnChannelWeight;
    TH2F* hCmnChannelWeightVsChannel;
    TH2F* hCmnFractionVsChannel;
	TH2F* hCmnNewVsNUsedChannels;
	TH2F* hCmnVsNewCmn;
	TH2F* hNewCmnVsEventNo;

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
	TH2F *hDiaAllAdcNoiseChannel;
	TH2F *hDiaAllAdcNoiseCMNChannel;
	TH2F *hLeftVsRightSignal[9];
	TH2F *hLeftVsRightSignalCMN[9];
	vector <Float_t> adcValues;
	vector <Float_t> pedestalValues;
	vector <Float_t> upperHitCutValues;
	vector <Float_t> lowerHitCutValues;
	vector <Float_t> upperSeedCutValues;
	vector <Float_t> lowerSeedCutValues;
	vector <Float_t> pedestalValuesCMN;
	vector <Float_t> upperHitCutValuesCMN;
	vector <Float_t> lowerHitCutValuesCMN;
	vector <Float_t> upperSeedCutValuesCMN;
	vector <Float_t> lowerSeedCutValuesCMN;
	vector <Float_t> eventNumbers;
	std::map<TString,TH1*> hHistoMap;
private:
	Float_t numberOfSeeds;
	Float_t sumPed;
	Float_t sumPedCMN;
	Float_t sumNoise;
	Float_t sumNoiseCMN;
	int nSumPed;
	int nSumPedCMN;
	int nSumNoiseCMN;
	int nSumNoise;

	Float_t biggestSignal;
	UInt_t biggestHitChannel;
	Float_t biggestSignalCMN;
	UInt_t biggestHitChannelCMN;

	Float_t cmNoise;

	UInt_t adc;
	Float_t snr;
	Float_t sigma;
	Float_t sigmaCMN;

	Float_t pedestal;
	Float_t noise;
	Float_t signal;

	Float_t signalCMN;
	Float_t pedestalCMN ;
	Float_t noiseCMN;
	Float_t cmn;
	bool isSaturated;
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
	std::vector<Float_t>vecBiggestAdjacentSignalInSigmaCMN[N_DETECTORS];
	std::vector<Int_t> vecBiggestAdjacentHitChannelCMN[N_DETECTORS];

	std::vector<Int_t> vecHitOrder[N_DETECTORS];
	std::vector< std::vector<UInt_t> > diaRawADCvalues; //vector of vector of adc Value (ch, eventNo)
	std::vector< std::set <UInt_t> > invalidBiggestHitChannels;
	THTMLPedestal *htmlPedestal;
};

#endif /* TANALYSISOFPEDESTAL_HH_ */

