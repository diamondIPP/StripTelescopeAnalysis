/*
 * TAnalysisOfSelection.hh
 *
 *  Created on: May 18, 2012
 *      Author: bachmair
 */

#ifndef TANALYSISOFSELECTION_HH_
#define TANALYSISOFSELECTION_HH_
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
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTMLPedestal.hh"
#include "LandauGaussFit.hh"
#include "THTMLLandaus.hh"
#include "THTMLSelectionAnalysis.hh"
#include "TTransparentAnalysis.hh"

#include "TADCEventReader.hh"
#include "TSettings.class.hh"

class TAnalysisOfSelection {
public:
	TAnalysisOfSelection(TSettings *settings);
	virtual ~TAnalysisOfSelection();
	void	doAnalysis(UInt_t nEvents=0);
    static void savePHvsEventNoAreaPlots(HistogrammSaver* histSaver, TSettings* settings, TProfile2D* prof2D,UInt_t xDivisions, UInt_t yDivisions);
private:
	void analyseEvent();
	void initialiseHistos();
	void saveHistos();
	void saveFidCutHistos();
	void saveDiamondAreaHistos();
	void initPHvsEventNoAreaPlots(UInt_t nStart, UInt_t nEnd);
	void initDividedAreaAxis(TAxis* axis);
	void fillPHvsEventNoAreaPlots(UInt_t area, UInt_t charge, UInt_t chargeOfTwo);
private:
	Int_t verbosity;
	TSettings *settings;
	HistogrammSaver *histSaver;
	TADCEventReader* eventReader;
	TSystem *sys;
	UInt_t nEvent;
	TH2F* histoLandauDistribution;
	TH2F* histoLandauDistribution2D;
	TH2F* histoLandauDistribution2DNoBorderSeed;
	TH2F* histoLandauDistribution2DNoBorderHit;
	TH2F* histoLandauDistribution2D_unmasked;
	TH2F* histoLandauDistribution2DNoBorderSeed_unmasked;
	TH2F* histoLandauDistribution2DNoBorderHit_unmasked;
	vector<TH2F*> hChargeVsFidX;
	vector<TH2F*> hChargeVsFidY;

	TH2F* hEtaVsLeftChannelNo;
	TH2F* hEtaVsRelPos;
	TH2F* hEtaCMNcorrectedVsLeftChannelNo;
	TH1F* hClusterPosition;
	TH1F* h3dDiamond;
	TH1F* hNoDiamond;
	TH1F* h3dDiamond_hit;
	TH1F* hNoDiamond_hit;
	TH2F* hFidCut;
	TH2F* hFidCutOneDiamondCluster;
	TH2F* hValidSiliconAndDiamondHit;
	TH2F* hValidSiliconAndOneDiamondHit;
	TH2F* hValidSiliconAndOneDiamondHitNotMasked;
	TH2F* hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels;
	TH2F* hValidSiliconAndOneDiamondHitInOneArea;
	TH2F* hValidSiliconAndOneDiamondHitInSameAreaAndFidCut;
	TH2F* hTwoClustersArea;
	TH1F* hNDiaClusters;
	TH3F* hChargeVsFidCut;
	TH2F* hFidCutXvsChannelPos;
	TH2F* hClusterSizeVsChannelPos;
	TH1F* hOneClusterHitChannel;
	TH2F* hOneClusterHitChannelVsArea;
	TH2F* hOneClusterHitChannelVsFiducialArea;
	TH2F* hOneClusterHitChannelAreaVsFiducialArea;
	THTMLLandaus *htmlLandau;
	THTMLSelectionAnalysis *htmlSelection;
	vector<Float_t> vecEta;
	vector<Float_t> vecSignalLeftOfEta;
	vector<Float_t> vecSignalRightOfEta;
	vector<Float_t> vecSignalLeftOfHighest;
	vector<Float_t> vecSignalRightOfHighest;
	vector<Float_t> vecHighestSignal;
	vector<Float_t> vecClusterCharge;
	vector<Float_t> vecTest;


    UInt_t xDivisions;
    UInt_t yDivisions;
    TProfile2D* hPHVsEventNo_Areas;
    TProfile2D* hPH2HighestVsEventNo_Areas;
};

#endif /* TANALYSISOFSELECTION_HH_ */
