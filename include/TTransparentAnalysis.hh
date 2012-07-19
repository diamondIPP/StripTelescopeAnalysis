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
//#include "FidCutRegion.hh"

#include "TSettings.class.hh"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TTrack.hh"
#include "TTracking.hh"
#include "TRawEventSaver.hh"
#include "TCluster.hh"
#include "THTMLTransparentAnalysis.hh"
#include "LandauGaussFit.hh"

using namespace std;

class TTransparentAnalysis {
public:
	TTransparentAnalysis(TSettings* settings);
	virtual ~TTransparentAnalysis();
	void	doAnalysis(int nEvents=0);
	void analyze(UInt_t nEvents, UInt_t startEvent);
	void setSettings(TSettings* settings);
	
private:
	void initHistograms();
	void fillHistograms();
	TF1* doGaussFit(TH1F *histo);
	void fitHistograms();
	void saveHistograms();
	void deleteHistograms();
	void printCutFlow();
	void fitTrack();
	void analyzeTrack(TTrack track);
	TCluster makeTransparentCluster(UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	bool checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	void printEvent();
	void printCluster(TCluster cluster);
	Float_t getResidual(TCluster cluster, TCluster::calculationMode_t clusterCalculationMode);
	
	UInt_t subjectDetector, subjectPlane;
	TPlaneProperties::enumCoordinate subjectDetectorCoordinate;
	vector<UInt_t> refPlanes;
	TCluster::calculationMode_t clusterCalcMode;
	UInt_t verbosity;
	UInt_t nEvent;
	
	TPositionPrediction* positionPrediction;
	vector<TCluster> transparentClusters;
	Float_t predXPosition, predYPosition;
	Float_t positionInDetSystem, predPerpPosition, predPosition;
	
	
    TSystem* sys;
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSettings* settings;
	TTracking* tracking;
	THTMLTransparentAnalysis* htmlTransAna;
	LandauGaussFit* landauGauss;
	
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
	vector<TH1F*> hResidualChargeWeighted;
	
	vector<TH1F*> hLaundau2Highest;
//	vector<TH1F*> hEta2Hightest;
	vector<TH1F*> hResidualHighest2Centroid;
	
	
	
//	TH1F* histo_transparentclustering_landau[10];
//    TH1F* histo_transparentclustering_landau_mean;
//    TH1F* histo_transparentclustering_eta;
//   	TH1F* histo_transparentclustering_hitdiff;
//   	TH2F* histo_transparentclustering_hitdiff_scatter;
//   	TH1F* histo_transparentclustering_2Channel_PulseHeight;
//   	TH1F* histo_transparentclustering_residuals[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
//   	TH2F* histo_transparentclustering_residuals_scatter[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
//   	TH1F* histo_transparentclustering_residuals_largest_hit[10];
//   	TH2F* histo_transparentclustering_residuals_largest_hit_scatter[10];
//   	TH1F* histo_transparentclustering_residuals_2largest_hits;
//   	TH2F* histo_transparentclustering_residuals_2largest_hits_scatter;
//   	TH1F* histo_transparentclustering_SNR_vs_channel;
//   	TH1F* histo_transparentclustering_chi2X;
//   	TH1F* histo_transparentclustering_chi2Y;
	
	// results
	vector<Float_t> vecMPLandau;
	vector<Float_t> vecMPLandau2Highest;
	vector<Float_t> vecMeanLandau;
	vector<Float_t> vecMeanLandau2Highest;
	vector< pair <Float_t,Float_t> > vecResidualChargeWeighted;
	vector< pair <Float_t,Float_t> > vecResidualHighest2Centroid;
	

};


#endif /* TTRANSPARENTANALYSIS_HH_ */
