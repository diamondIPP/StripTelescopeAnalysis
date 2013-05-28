/*
 * TAnalysisOf3dDiamonds.hh
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair, iain
 */
#ifndef TANALYSISOF3DDIAMONDS_HH_
#define TANALYSISOF3DDIAMONDS_HH_
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
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTMLPedestal.hh"
#include "LandauGaussFit.hh"

#include "TADCEventReader.hh"
#include "TTracking.hh"
#include "TSettings.class.hh"
#include "TCluster.hh"

using namespace std;

class TAnalysisOf3dDiamonds {
public:
	TAnalysisOf3dDiamonds(TSettings *settings);
	virtual ~TAnalysisOf3dDiamonds();
	void	doAnalysis(UInt_t nEvents=0);

private:
	void initialiseHistos();
	void saveHistos();
	void analyseEvent();

	//YAlignment
	void initialiseYAlignmentHistos();
	void saveYAlignmentHistos();
	void YAlignment();
	void DrawYAlignmentFidCutRegions();
	int YAlignmentFiducialCut(int nRegion);
	void DrawMetallisationGrid();
	Float_t GetYPositionInDetSystem();

	//new functions
	void HitandSeedCount(TCluster* nCluster, int ni);
	void ClusterShape(TCluster* nCluster);
	void RemoveLumpyClusters(TCluster* nCluster);
	void RemoveEdgeHits(TCluster* nCluster, int* nEdgeChannel);
	int FiducialCut(int Detector);
	void DrawFidCutRegions();
	void PredictChannelHit(TCluster* nCluster);


private:

	string FileNameEnd;
	vector<TH2F*> hPHvsPredictedChannel, hPHvsChannel, hPHvsPredictedXPos, hPredictedPositionDiamondHit, hHitandSeedCount, hChi2XChi2Y, hFidCutXvsFidCutY;
	vector<TH1F*>hEventsvsChannel;
	TH1F* hNumberofClusters;
	TH1F* hDoubleClusterPos;
	TCanvas* DoubleCluster;
	TH1F* hDoubleClusterPos0;
	TH1F* hDoubleClusterPos1;
	TH1F* hLandauCluster1;
	TH1F* hLandauCluster2;
	TH1F* hLandauDoubleCombined;

	//vector<TH3F*> hFidCutXvsFidCutYvsCharge;
	//For hFidCutXvsFidCutYvsMeanCharge
	vector<TCanvas*> ptrCanvas, ptrCanvasEvents, ptrCanvasMean;
	vector<TH2D*> hFidCutXvsFidCutYvsCharge, hFidCutXvsFidCutYvsEvents, hFidCutXvsFidCutYvsMeanCharge;
	TH2D* hFidCutXvsFidCutYvsMeanChargeAllDetectors;
	//For XPosvsYPosvsMeanCharge
	TH2D* hXPosvsYPosvsCharge;
	TH2D* hXPosvsYPosvsEvents;
	TH2D* hXPosvsYPosvsMeanCharge;
	TH2D* hFidCutXvsFidCutYvsPredictedEvents;
	TH2D* hFidCutXvsFidCutYvsSeenEvents;
	TH2D* hFidCutXvsFidCutY0Clusters;
	TH2D* hFidCutXvsFidCutY1_1Clusters;
	TH2D* hFidCutXvsFidCutY1Clusters;
	TH2D* hFidCutXvsFidCutY2Clusters;
	TH2D* hFidCutXvsFidCutY2Clusters0;
	TH2D* hFidCutXvsFidCutY2Clusters1;
	TH2D* hFidCutXvsFidCutY3Clusters;
	TH2D* hEfficiency;
	//TH3F* hFidCutXvsFidCutYvsCharge;
	vector<TH1F*> hLandau;
	TH2F* hPHvsChannelStrip;
	TH2F* hPHvsPredictedXPosStrip;
	TH2F* hPHvsChannel3dNoHoles;
	TH2F* hPHvsPredictedXPos3dNoHoles;
	TH2F* hPHvsChannel3dWithHoles;
	TH2F* hPHvsPredictedXPos3dWithHoles;

	//For Fiducial Cut Boundaries Histograms
	vector <TBox*> FidCutChannelTBox;
	TCanvas* hCombinedMeanCharge;    //TCanvas for MeanCharge All detectors
	TCanvas* hPredictedEvents;
	TCanvas* hSeenEvents;
	TCanvas* h0Clusters;
	TCanvas* h1_1Clusters;
	TCanvas* h1Clusters;
	TCanvas* h2Clusters;
	TCanvas* h2Clusters0;
	TCanvas* h2Clusters1;
	TCanvas* h3Clusters;
	float FiducialMetric;			// Which Fiducial cut to apply
	float FiducialChannel;
	TCanvas* FidCudBoundMetricCanvas;
	TCanvas* FidCudBoundChannelCanvas;
	TH2F* FidCudBoundMetric;
	TH2F* FidCudBoundChannel;

	/////////////
	//For YAlignment Histos
	/////////////
	TCanvas* hCombinedMeanChargeYAlignment;
	TH2D* hFidCutXvsFidCutYvsChargeYAlignment;
	TH2D* hFidCutXvsFidCutYvsEventsYAlignment;
	TH2D* hFidCutXvsFidCutYvsMeanChargeYAlignment;
	TCanvas* h3DdetMeanCharge;
	TH2D* hDetXvsDetY3D;
	TH2D* hDetXvsDetY3DvsEvents;
	TH2D* hDetXvsDetY3DMeanCharge;
	TCanvas* h3DdetMeanChargeRebinned;
	TH2D* hDetXvsDetY3DRebinned;
	TH2D* hDetXvsDetY3DvsEventsRebinned;
	TH2D* hDetXvsDetY3DMeanChargeRebinned;
	vector<TH2D*> hCellsCharge;
	vector<TH2D*> hCellsEvents;
	TCanvas* hCellsOverlayedCanvas;
	TH2D* hCellsOverlayedCharge;
	TH2D* hCellsOverlayedEvents;
	TH2D* hCellsOverlayedMeanCharge;
	TH1F* hBinnedMeanCharge;
	//xEdge
	TH1F* hEdgeCharge;
	TH1F* hEdgeChargeEvents;
	TCanvas* hEdgeMeanChargeCanvas;
	TH1F* hEdgeMeanCharge;
	//yEdge
	TH1F* hyEdgeCharge;
	TH1F* hyEdgeChargeEvents;
	TCanvas* hyEdgeMeanChargeCanvas;
	TH1F* hyEdgeMeanCharge;
	//DeadCellProfile
	TCanvas* hDeadCellMeanChargeCanvas;
	TH1F* hDeadCell;
	TH1F* hDeadCellEvents;
	TH1F* hDeadCellMeanCharge;

	int* DeadCellsArrayPointer;
	TH1F* hDeadCellsProfileCharge;
	TH1F* hDeadCellsProfileEvents;
	TH1F* hDeadCellsProfileMeanCharge;
	vector<TH1F*> hDeadCellsCharge;
	vector<TH1F*> hDeadCellsEvents;
	//Fiducial Cut
	vector <TBox*> FidCutChannelYAlignmentTBox;
	vector<int*> FidCutYAlignment;
	//vector<float> FidYAlignment;


	TSettings *settings;
	HistogrammSaver *histSaver;
	TTracking* eventReader;
	UInt_t nEvent;
	vector<Float_t> vecXPredictedDiamondHit,vecYPredictedDiamondHit,vecPHDiamondHitStrip,vecPHDiamondHit3dNoHoles,vecPHDiamondHit3dWithHoles,vecClusterSize,vecChi2Y,vecChi2X,vecClusterSeedSize,vecXPredictedStrip,vecYPredictedStrip;
	vector< vector<Float_t>* > vecPHDiamondHit, vecXPredicted, vecYPredicted, vecClusterFrequency;
	int HitCount, SeedCount;

	vector<int*> Detector;
	vector<int*> FidCut;

};

#endif /* TANALYSISOF3DDIAMONDS_HH_ */
