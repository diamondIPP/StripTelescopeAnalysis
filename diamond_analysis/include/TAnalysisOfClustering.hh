/*
 * TDeadChannels.hh
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#ifndef TDEADCHANNELS_HH_
#define TDEADCHANNELS_HH_

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
#include "TSettings.class.hh"
#include "LandauGaussFit.hh"
#include "TProfile.h"

#include "TADCEventReader.hh"
#include "THTMLCluster.hh"
#include "TClustering.hh"

#define N_DET_CHANNELS 256
#define N_DIA_CHANNELS 128
using namespace std;

class TAnalysisOfClustering {
public:
	TAnalysisOfClustering(TSettings *settings);
	virtual ~TAnalysisOfClustering();
	void	doAnalysis(int nEvents=0);
	void setSettings(TSettings* settings);
private:
	void analyseEvent();
	void saveHistos();
	void saveEtaIntegrals();
	void savePHHistos();

	void initialiseHistos();
	void checkForDeadChannels();
	void compareCentroid_ChargeWeightedMean();
	void analyseForSeeds();
	void analyse2ndHighestHit();
	void analyseClusterPosition();
	void checkForSaturatedChannels();
	void analyseCluster();
	void createPHDistribution();

	void etaInvestigation();
	void saveEtaInvestigationHistos();
	void saveEtaDividedHistos(TH3F* h3DLeft,TH3F* h3DRight, TH2F* h2DLeft, TH2F* h2DRight,string name_comparision, Float_t etaWidth=.2);
	TH1F *hSaturatedChannels[9];
	TH1F *hSeedMap[9];
	TH1F *hSeedMap2[9];
	TH1F *hClusterMap[9];
	TH1F* hNumberOfSeeds[9];
	TH1F* hChannelBiggestHit[9];
	TH1F* hPulsHeightBiggestHit[9];
	TH1F* hPulsHeightNextBiggestHit[9];
	TH1F* hNumberOfClusters[9];
	TH1F* hClusterSize[9];
	TH1F* hClusterSeedSize[9];
	TH2F* hBiggestHitVsClusterSize[9];
	vector <double>vecClusterSizes,vecClusterSeedSizes,vecNumberOfClusters;
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings* settings;
    UInt_t verbosity;
    int nEvent;
    int seedSigma;
    int hitSigma;
    vector<double>vecPHMeans;
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
	TH1F *histo_H2C_biggestHit;
	TH2F *histo_CWM_biggestHit;
private:
	THTMLCluster *htmlClus;
	TH1F *h2ndBiggestHitSignal[9];
	TH1F *h2ndBiggestHitOverCharge[9];
	TH1F *h2ndBiggestHitPosition[9];
	TH1F *hLeftHitOverLeftAndRight[9];
	TH1F *hDeltaLeftRightHitOverLeftAndRight[9];
	TH1F *hSignal2ndHighestOverSignalHighestRatio[9];
	TH2F *hRelativeClusterPositionCWM[9];
	TH2F *hRelativeClusterPositionCorEta[9];
	TH2F *hRelativeClusterPositionEta[9];
	TH1F *hClusterPosition[9];
	TH1F *hEtaDistribution[9];
	TH1F *hEtaDistributionCMN[9];
	TH2F* hEtaDistributionVsLeftChannel[9];
	TH2F* hEtaDistributionVsClusterSize[9];
	TH2F* hEtaDistributionVsCharge[9];
	TH1F *hEtaDistribution5Percent[9];
	TH2F *hEtaDistributionVsSignalRight[9];
	TH2F *hEtaDistributionVsSignalLeft[9];
	TH2F *hEtaDistributionVsSignalSum[9];
	TH2F *hSignalLeftVsSignalRight[9];
	TH2F *hPHDistribution[9];
private:
	vector < vector <Float_t> > vecvecSignalLeftLeft;
	vector < vector <Float_t> > vecvecSignalRightRight;
	vector < vector <Float_t> > vecvecLeftEtaSignal;
	vector < vector <Float_t> > vecvecRightEtaSignal;
	vector < vector <Float_t> > vecvecSignalLeftOfHighest;
	vector < vector <Float_t> > vecvecSignalRightOfHighest;
	vector < vector <Float_t> > vecvecSignalHighest;
	vector < vector <Float_t> > vecvecEta;
};
#endif /* TDEADCHANNELS_HH_ */

