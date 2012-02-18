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
#include "TTrack.hh"
#include "TTracking.hh"
#include "TRawEventSaver.hh"
#include "TCluster.hh"

using namespace std;

class TTransparentAnalysis {
public:
	TTransparentAnalysis(int runNumber, TSettings settings);
	virtual ~TTransparentAnalysis();
	void	doAnalysis(int nEvents=0);
	void analyze(int nEvents, int startEvent);
	void setSettings(TSettings* settings);
	
private:
	void initHistograms();
	void fillHistograms();
	void saveHistograms();
	void deleteHistograms();
	void printCutFlow();
	void fitTrack();
	void analyzeTrack(TTrack track);
	TCluster makeTransparentCluster(UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	bool checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	
	UInt_t subjectDetector;
	vector<UInt_t> refPlanes;
	TCluster::calculationMode_t clusterCalcMode;
	
	TPositionPrediction* positionPrediction;
	vector<TCluster> transparentClusters;
	Float_t predXPosition, predYPosition;
	Float_t positionInDetSystem, predPerpPosition, predPosition;
	
	
    TSystem* sys;
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSettings* settings;
	TTracking* tracking;
	
	// cut flow
	UInt_t nAnalyzedEvents;
	UInt_t regionNotOnPlane;
	UInt_t saturatedChannel;
	UInt_t screenedChannel;
	UInt_t noValidTrack;
	UInt_t noFidCutRegion;
	
	
	// histograms
	vector<TH1F*> hLaundau;
	vector<TH1F*> hEta;
	vector<TH1F*> hResidual;
	
	vector<TH1F*> hLaundau2Hightest;
	vector<TH1F*> hEta2Hightest;
	vector<TH1F*> hResidual2Hightest;
	
	
	
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
	
	
	

};


#endif /* TTRANSPARENTANALYSIS_HH_ */
