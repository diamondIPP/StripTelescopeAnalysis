//
//  TTransparentAnalysis.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 05.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TTRANSPARENTANALYSIS_HH_
#define TTRANSPARENTANALYSIS_HH_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstring>
#include <deque>


//ROOT Class Headers
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h" // for adding your own classes to ROOT's library
#include "TStyle.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "FidCutRegion.hh"

#include "TSettings.class.hh"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"

using namespace std;

class TTransparentAnalysis {
public:
	TTransparentAnalysis(int runNumber);
	virtual ~TTransparentAnalysis();
	void	doAnalysis(int nEvents=0);
	
private:
	void initHistograms();
	void saveHistograms();
	void fitTrack();
	void analyze(int nEvents);
	
	
    TSystem* sys;
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
	
	
	// histograms
	TH1F* histo_transparentclustering_landau[10];
    TH1F* histo_transparentclustering_landau_mean;
    TH1F* histo_transparentclustering_eta;
   	TH1F* histo_transparentclustering_hitdiff;
   	TH2F* histo_transparentclustering_hitdiff_scatter;
   	TH1F* histo_transparentclustering_2Channel_PulseHeight;
   	TH1F* histo_transparentclustering_residuals[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
   	TH2F* histo_transparentclustering_residuals_scatter[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
   	TH1F* histo_transparentclustering_residuals_largest_hit[10];
   	TH2F* histo_transparentclustering_residuals_largest_hit_scatter[10];
   	TH1F* histo_transparentclustering_residuals_2largest_hits;
   	TH2F* histo_transparentclustering_residuals_2largest_hits_scatter;
   	TH1F* histo_transparentclustering_SNR_vs_channel;
   	TH1F* histo_transparentclustering_chi2X;
   	TH1F* histo_transparentclustering_chi2Y;
	
	
	
//	void saveHistos();
//	void checkForDeadChannels();
//	void analyseForSeeds();
//	void getBiggestHit();
//	void initialiseHistos();
//	void checkForSaturatedChannels();
//	void analyseCluster();
//    void analyseBiggestHit();
//	TH1F *hSaturatedChannels[8];
//	TH1F *hSeedMap[8];
//	TH1F *hSeedMap2[9];
//	TH1F *hClusterMap[8];
//	TH1F* hNumberOfSeeds[8];
//	TH1F* hChannelBiggestHit[8];
//	TH1F* hPulsHeightBiggestHit[8];
//	TH1F* hPulsHeightNextBiggestHit[8];
//	TH1F* hNumberOfClusters[9];
//	TH1F* hClusterSize[9];
//	TADCEventReader* eventReader;
//	HistogrammSaver* histSaver;
//    TSystem* sys;
//    int nEvent;
//    int seedSigma;
//    int hitSigma;
//    TH1F *histo_pulseheight_sigma[8];
//	TH1F *histo_pulseheight_sigma_second[8];
//	TH1F *histo_pulseheight_sigma125[8];
//	TH1F *histo_second_biggest_hit_direction[8];
//	TH1F *histo_pulseheight_sigma_second_left[8];
//	TH1F *histo_pulseheight_sigma_second_right[8];
//	TH1F *histo_biggest_hit_map[8];
//	TH1F *histo_pulseheight_left_sigma[8];
//	TH1F *histo_pulseheight_left_sigma_second[8];
//	TH1F *histo_pulseheight_right_sigma[8];
//	TH1F *histo_pulseheight_right_sigma_second[8];
};


#endif /* TTRANSPARENTANALYSIS_HH_ */
