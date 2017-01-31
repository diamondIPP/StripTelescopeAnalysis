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
#include <map>


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
#include "TProfile2D.h"
#include "TRandom.h"
#include "TFitResult.h"//TFitResultPtr.h"
//#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TList.h"
#include "TPolyMarker.h"
#include "TObjArray.h"
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
#include "TClustering.hh"
#include "TResults.hh"
#include "TAnalysisOfSelection.hh"
#include "TAnalysisOfClustering.hh"
#include "TSpectrum.h"


using namespace std;

class TTransparentAnalysis {
public:
	TTransparentAnalysis(TSettings* settings, TSettings::alignmentMode mode = TSettings::normalMode);
	virtual ~TTransparentAnalysis();
//	void	doAnalysis(int nEvents=0);
	void analyze(UInt_t nEvents, UInt_t startEvent);
	void calcEtaCorrectedResiduals();
	void setSettings(TSettings* settings);
	static TCluster makeTransparentCluster(TTracking *reader,TSettings* set, UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	void setResults(TResults* res){cout<<"Setting results!"<<endl;results=res;};
private:
	void clearEventVector();
	void createEventVector(Int_t startEvent = 0);
	void analyseNonHitEvents();
	void initHistograms();
	void initHistograms2();
	void initPedestalAndNoiseHistos(UInt_t maxEvents=1e6);
	void initClusteredHistos(UInt_t startEvent = 1e5, UInt_t maxEvents=1e6);
	void fillHistograms();
	void fillClusteredHistos();
	void fillPedestalsAndNoiseHistos();
	void saveClusteredHistos();
	void savePedestalHistos();
    void saveNoiseHistos();
    Float_t GetFractionOutsideNSigma(TH1* hRes,Float_t mean, Float_t sigma,Int_t nSigma=2);
	TF1* doGaussFit(TH1F *histo,Float_t xmin =0,Float_t xmax=0);
	TF1* doDoubleGaussFit(TH1F *histo);
	TF1* doFixedDoubleGaussFit(TH1F* histo);
	TF1* doGaussPlusStepFunction(TH1F* histo);
	void createEtaIntegrals();
	void fitHistograms();
	void createEfficiencyPlots(TH1F* hLandauDistribution);
	void analyseEtaDistributions();
	void analyseEtaDistribution(TH1F* hEtaDist);
	void saveHistograms();
	void saveLandausVsPositionPlots(UInt_t clusterSize);
	void SaveLandauVsEventNoPlots(UInt_t clusterSize);
	void AnalyzeLandauVsEventNo(TH2* hLandauVsEventNo);
    void AnalyzeLandauVsEventNoMaxBin(TH2* hLandauVsEventNo);
    void AnalyzeLandauVsEventNoFitSlices(TH2* hLandauVsEventNo);
	void deleteHistograms();
	void deleteFits();
	void printCutFlow();
//	void fitTrack();
//	void analyzeTrack(TTrack track);
	bool predictPositions(bool savePrediction = true);
	bool checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize);
	static int getSignedChannelNumber(Float_t position);
	void printEvent();
	void printCluster(TCluster cluster);
	Float_t getResidual(TCluster cluster,bool cmnCorrected, TCluster::calculationMode_t clusterCalculationMode, TH1F* hEtaInt=0);
	std::pair<Float_t,Float_t >  getFWCrossingPoint(TH1F* hRes,Float_t crossingPoint);
	void saveResolutionPlot(TH1F* hRes, UInt_t clusterSize,TString additionalInfo="");
	// run variables
	UInt_t subjectDetector, subjectPlane;
	TPlaneProperties::enumCoordinate subjectDetectorCoordinate;
	vector<UInt_t> refPlanes;
	TCluster::calculationMode_t clusterCalcMode;
	UInt_t verbosity;
	UInt_t nEvent;
	
	// event variables
	TPositionPrediction* positionPrediction;
	TCluster transparentClusters;
	vector<TCluster> noHitClusters;
	vector<TEvent* > vecEvents;
	Float_t predXPosition, predYPosition;
	Float_t positionInDetSystemMetric,positionInDetSystemChannelSpace, positionInDetSystemMetricY, predPerpPosition, predPosition;
	
	// sys variables
    TSystem* sys;
	HistogrammSaver* histSaver;
    TSettings* settings;
    TResults* results;
	TTracking* eventReader;
	THTMLTransparentAnalysis* htmlTransAna;
	LandauGaussFit* landauGauss;
	
	// cut flow
	UInt_t nEvents;
	UInt_t nAnalyzedEvents;
	UInt_t regionNotOnPlane;
	UInt_t saturatedChannel;
	UInt_t screenedChannel;
	UInt_t noValidTrack;
	UInt_t noFidCutRegion;
	UInt_t usedForAlignment;
	UInt_t highChi2;
	//	UInt_t usedForSiliconAlignment;
	// data for Histos

	vector<Float_t> vecEta;
	vector<Float_t> vecSignalLeftOfEta;
	vector<Float_t> vecSignalRightOfEta;
	vector<Float_t> vecSignalLeftOfHighest;
	vector<Float_t> vecSignalRightOfHighest;
	vector<Float_t> vecClusterCharge;
	vector<Float_t> vecHighestSignal;

	vector< vector < Float_t> > vecvecResXChargeWeighted;
	vector< vector < Float_t> > vecvecResXEtaCorrected;
	vector< vector < Float_t> > vecvecResXHighest2Centroid;
	vector< vector < Float_t> > vecvecResXHighestHit;
	vector< vector<Float_t> > vecvecRelPos;
	vector< vector<Float_t> > vecvecRelPos2;
	vector< vector<Float_t> > vecvecEta;
	vector< vector<Float_t> > vecvecEtaCMNcorrected;
	vector< Float_t> vecDeltaEta;
	vector< Float_t> vecRelatedEta2;
	vector< Float_t> vecRelatedEta10;
	vector< Float_t> vecRelatedResXEtaCorrected;

	// transparent results tree
	static const int trTree_nMaxStrips = 10;
	void bookTransparentTree();
	void resetTransparentTree();
	void fillTransparentTree();
	void writeTransparentTree();
	TFile* transparentTreeFile;
	TTree* transparentTree;
	Bool_t  trTree_ValidTrack;
	Bool_t  trTree_AlignmentTrack;
	Bool_t  trTree_inFiducialRegion;
	Bool_t  trTree_ValidPredRegion;
	Bool_t  trTree_ValidChi2;
	Bool_t  trTree_ScreenedCluster;
	Float_t trTree_Chi2;
	Int_t   trTree_RunNumber;
	Float_t trTree_tseed;
	Float_t trTree_thit;
	Int_t   trTree_Event;
	Int_t   trTree_Direction;
	Int_t   trTree_NStrips;
	Int_t   trTree_CenterStrip;
	Float_t trTree_ADC                 [10];
	Float_t trTree_Pedestal            [10];
	Float_t trTree_PedestalSigma       [10];
	Float_t trTree_Signal              [10];
	Float_t trTree_CMN;
	Float_t trTree_PedestalCMNcorr     [10];
	Float_t trTree_PedestalSigmaCMNcorr[10];
	Float_t trTree_SignalCMNcorr       [10];
	Float_t trTree_HighestStripADC                 [10];
	Float_t trTree_HighestStripPedestal            [10];
	Float_t trTree_HighestStripPedestalSigma       [10];
	Float_t trTree_HighestStripSignal              [10];
	Float_t trTree_HighestStripSNR                 [10];
	Float_t trTree_HighestStripADCCMNcorr          [10];
	Float_t trTree_HighestStripPedestalCMNcorr     [10];
	Float_t trTree_HighestStripPedestalSigmaCMNcorr[10];
	Float_t trTree_HighestStripSignalCMNcorr       [10];
	Float_t trTree_HighestStripSNRCMNcorr          [10];
	Float_t trTree_StripADC                        [128];
	Float_t trTree_StripPedestal                   [128];
	Float_t trTree_StripPedestalSigma              [128];
	Float_t trTree_StripSignal                     [128];
	Float_t trTree_StripSNR                        [128];
	Float_t trTree_StripADCCMC                     [128];
	Float_t trTree_StripPedestalCMC                [128];
	Float_t trTree_StripPedestalSigmaCMC           [128];
	Float_t trTree_StripSignalCMC                  [128];
	Float_t trTree_StripSNRCMC                     [128];

	Float_t inf;
	// histograms
	vector<TH1F*> hLandau;
	vector<TProfile2D*> hLandau2HighestProfile2D;
	vector< TH2F*> hLandau2HighestFidCutX;
	vector< TH2F*> hLandau2HighestFidCutY;
	vector< TH2F*> hLandau2HighestPredX;
	vector< TH2F*> hLandau2HighestPredY;
	vector< vector< Float_t> > vecVecLandau;
	vector< Float_t> vecPredictedChannel;
	vector< Float_t> vecPredictedDetectorPositionY;
	Int_t predChannel;
	vector< Float_t> vecVecFidCutX;
	vector< Float_t> vecVecFidCutY;
	vector< Float_t> vecPredX;
	vector<Float_t> vecPredY;
	vector<TH1F*> hEta;
	vector<TH1F*> hEtaCMNcorrected;
	vector< vector<Float_t> > vecVecEta;

	TH1F* hLandauMean;
	TH1F* hLandauMP;
	TH1F* hPredictedPositionInStrip;

	TH2F* hSignalInSNR_Centroid2Strips_Dia;
	TH2F* hSignalInSNR_Centroid2Strips_Dia_cmnCor;
	TH2F* hSignalInSNR_Centroid2StripsSorted_Dia;
	TH2F* hSignalInSNR_Centroid2StripsSorted_Dia_cmnCor;
	TH2F* hSignalInSNR_TrackStripLeft_Dia;
	TH2F* hSignalInSNR_TrackStripLeft_Dia_cmnCor;
	TH2F* hSignalInSNR_TrackStripRight_Dia;
	TH2F* hSignalInSNR_TrackStripRight_Dia_cmnCor;
	
	vector<TH1F*> hLandau2Highest;
	vector<TH1F*> hLandau2Highest_nonCMC;
	vector<TH1F*> hLandau1Highest;
	vector<TH1F*> hLandau2HighestFixedNoise;
	vector <TH1F*> hLandauFixedNoise;
//	vector<TH1F*> hEta2Hightest;
	vector<TH1F*> hResidualHighest2Centroid;
	vector<TH1F*> hResidualHighestHit;
	vector<TH1F*> hResidualEtaCorrected;
	vector<TH1F*> hResidualChargeWeighted;
	vector<TH2F*> hResidualVsHitPositionChargeWeighted;
	vector<TH2F*> hResidualVsHitPositionHigehest2Centroid;
	vector<TH2F*> hResidualVsHitPositionEtaCorrected;
	TH1F* hLandau2HighestMean;
	TH1F* hLandau2HighestMP;
	TH2F* hSelectedTracksAvrgSiliconHitPos;
	vector<TH1F*> hEtaIntegrals;
	vector<TH1F*> hEtaCMNcorrectedIntegrals;
	TH1F* hEtaIntegral_clustered;
//	TH2F* hResidualEtaVsEstimatedHitPosition,hResidualChargeWeightedVsEstimatedHitPosition,hResidualHighest2CentroidVsEstimatedHitPosition;
	
	// fits
	vector<TF1*> fitLandau;
	vector<TF1*> fitLandauFixedNoise;
	vector<TF1*> fitLandau2Highest;
	vector<TF1*> fitLandau2HighestFixedNoise;
	vector<TF1*> fitResidualChargeWeighted;
	vector<TF1*> fitResidualHighest2Centroid;
	vector<TF1*> fitResidualEtaCorrected;
	TSettings::alignmentMode alignMode;
	
	
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
	vector<UInt_t> eventNumbers;
	vector<TCluster> vecTransparentClusters;
	vector<Float_t> vecMPLandau;
	vector<Float_t> vecMPLandau2Highest;
	vector<Float_t> vecMeanLandau;
	vector<Float_t> vecMeanLandau2Highest;
	vector< pair <Float_t,Float_t> > vecResidualChargeWeighted;
	vector< pair <Float_t,Float_t> > vecResidualHighest2Centroid;
	vector< pair <Float_t,Float_t> > vecResidualEtaCorrected;
	vector< pair <Float_t,Float_t> > vecResidualEtaCorrected_2ndGaus;
	vector<Float_t> vectorEventNo;
	vector<Float_t> vectorCMN;
	vector< vector<Float_t> > vecVecPh2Highest;
	
	vector<Float_t> vecPredictedPosition, vecRelPredictedPosition;
	vector<Float_t> vecChi2;
	Float_t predXMin, predXMax, predYMin, predYMax;
	bool cmCorrected;
	vector<Float_t> noiseWidths;
	vector<Float_t> noiseWidthsCMN;
	vector<Float_t> noiseWidths2OutOfX;
	vector<Float_t> noiseWidths2OutOfXCMN;

	std::map< UInt_t, TProfile* > hPedestalVsEvenNo;
    std::map< UInt_t, TProfile* > hNoiseVsEvenNo;

	TProfile* hCmnVsEventNo;

	vector<TProfile2D*> vecPHVsEventNo_Areas;
	vector<TProfile2D*> vecPH2HighestVsEventNo_Areas;


	TH2F* hLandauVsEventNo_Clustered;
	TH1F* hClusterSize_Clustered;
	TH2F* hLandauVsClusterSize_Clustered;
	TH1F* hNClusteres_Clustered;
	TH2F* hEtaVsClusterSize_clustered;
	TH2F* hEtaCMNcorrectedVsClusterSize_clustered;
	TH1F* hEta_clustered;
	TH1F* hResidualHighest2Centroid_clustered;
	TH1F* hResidualHighestHit_clustered;
	TH1F* hResidualEtaCorrected_clustered;
	TH1F* hResidualChargeWeighted_clustered;
	TH2F* hResidualHighest2CentroidVsClusterSize_clustered;
	TH2F* hResidualHighestHitVsClusterSize_clustered;
	TH2F* hResidualEtaCorrectedVsClusterSize_clustered;
	TH2F* hResidualChargeWeightedVsClusterSize_clustered;
	TH2F* hResidualTranspEtaCorrectedVsClusterSize_clustered;
	TH2F* hResidualEtaCorrectedOnlyHitsVsClusterSize_clustered;
	TH2F* hResidualTranspEtaCorrectedOnlyHitsVsClusterSize_clustered;

	Float_t minX,maxX,minY,maxY;

public:
	static UInt_t GetHitArea(TSettings* set, Float_t xVal, Float_t yVal,UInt_t xDivisions,UInt_t yDivisions);
private:
	void initDividedAreaAxis(TAxis *axis);
public:
	static TString GetNameOfArea(Int_t x,Int_t y);
private:
	void initPHvsEventNoAreaPlots(UInt_t nStart = 0, UInt_t nEnd = 1e6);
	void fillPHvsEventNoAreaPlots(UInt_t area, UInt_t clusterSize, UInt_t charge, UInt_t chargeOfTwo);
	TH2D* hPHVsEventNo;
	TH2D* hPH2OutOf10VsEventNo;
	void savePHvsEventNoAreaPlots();
	UInt_t xDivisions;
	UInt_t yDivisions;
};


#endif /* TTRANSPARENTANALYSIS_HH_ */
