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
#include "TCellAnalysisClass.hh"

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
	void initialise3DGridReference();
	void initialise3DYAlignmentHistos();
	void initialise3DOverviewHistos();
	void initialiseYAlignmentHistos();
	void saveYAlignmentHistos();
	void YAlignment();
	void DrawYAlignmentFidCutRegions();
	Int_t* CellToChannel(int nCell,float nXcell);
	int YAlignmentFiducialCut(int nRegion);
	void DrawMetallisationGrid(TCanvas* nCanvas);
	Float_t GetYPositionInDetSystem();
	float* SortArrayBtoS(float* nArray, int nSize);
	void printArray(float* nArray, int nSize, const std::string& space);

	//new functions
	void HitandSeedCount(TCluster* nCluster, int ni);
	void ClusterShape(TCluster* nCluster);
	void RemoveLumpyClusters(TCluster* nCluster);
	void RemoveEdgeHits(TCluster* nCluster, int* nEdgeChannel);
	void ClusterChannels(TCluster* nCluster);
	int RemoveClustersWithHitOutside3D(TCluster* nCluster);
	int RemoveEdge3DClusters(TCluster* nCluster);
	int FiducialCut(int Detector);
	void DrawFidCutRegions();
	void PredictChannelHit(TCluster* nCluster);

	pair<int,int> getCellNo(Float_t xPos, Float_t yPos);


private:
	TCellAnalysisClass* clusteredAnalysis;
	float SortArrayPointer[4];

	string FileNameEnd;
	vector<TH2F*> hPHvsPredictedChannel, hPHvsChannel, hPHvsPredictedXPos, hPredictedPositionDiamondHit, hHitandSeedCount, hChi2XChi2Y, hFidCutXvsFidCutY;
	vector<TH1F*>hEventsvsChannel;
	TH1F* hNumberofClusters;
	TH1F* hDoubleClusterPos;
	TCanvas* cDoubleCluster;
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
	TCanvas* cPredictedEvents;
	TCanvas* cSeenEvents;
	TCanvas* c0Clusters;
	TCanvas* c1_1Clusters;
	TCanvas* c1Clusters;
	TCanvas* c2Clusters;
	TCanvas* c2Clusters0;
	TCanvas* c2Clusters1;
	TCanvas* c3Clusters;
	float FiducialMetric;			// Which Fiducial cut to apply
	float FiducialChannel;
	TCanvas* FidCudBoundMetricCanvas;
	TCanvas* cFidCudBoundChannel;
	TH2F* FidCudBoundMetric;
	TH2F* FidCudBoundChannel;

	/////////////
	//For YAlignment Histos
	/////////////
	TCanvas* hGridReferenceCanvas;
	TH2D* hGridReference;
	TCanvas* cCombinedMeanChargeYAlignment;
	TH2D* hFidCutXvsFidCutYvsChargeYAlignment;
	TH2D* hFidCutXvsFidCutYvsEventsYAlignment;
	TH2D* hFidCutXvsFidCutYvsMeanChargeYAlignment;
	TCanvas* c3DdetMeanCharge;
	TH2D* hDetXvsDetY3D;
	TH2D* hDetXvsDetY3DvsEvents;
	TH2D* hDetXvsDetY3DMeanCharge;
	TCanvas* c3DdetMeanChargeWithMeanClusterSize;
	TH2D* hCellsMeanClusteSize;
	TCanvas* c3DdetMeanChargeRebinned;
	TCanvas* h3DdetEventsRebinned;
	TH2D* hDetXvsDetY3DRebinned;
	TH2D* hDetXvsDetY3DvsEventsRebinned;
	TH2D* hDetXvsDetY3DMeanChargeRebinned;
	TCanvas* cRebinnedQuarterCellFails;
	TH2D* RebinnedQuarterCellFails;
	TCanvas* cDetXvsDetY3DRebinnedRMS;
	TH2D* hDetXvsDetY3DRebinnedRMS;
	vector<TH1F*> hQuaterCellsLandau;
	vector<TH1F*> hQuarterCellsClusterSize;
	TCanvas* c2DClusterSizeQuarterCell;
	TH2D* h2DClusterSizeQuarterCell;
	TH2D* h2DClusterSizeQuarterCellClone;
	TH2D* h2DClusterSizeClone;
	TH2D* h2DClusterSizeQuarterCellClone1;
	TH2D* h2DClusterSizeClone1;
	TCanvas* c2DMeanClusterSizeQuarterCell;
	TH2D* h2DMeanClusterSizeQuarterCell;
	TH2D* h2DMeanClusterSizeQuarterCellTotal;
	TH2D* h2DMeanClusterSizeQuarterCellEvents;
	TH2D* h2DClusterSizeXAxis;
	TCanvas* cDetXvsDetY3DMeanChargeRebinnedQuarterCell;
	TH2D* hDetXvsDetY3DMeanChargeRebinnedQuarterCell;
	TCanvas* cDetXvsDetY3DRebinnedQuarterCellRMSCanvas;
	TH2D* hDetXvsDetY3DRebinnedQuarterCellRMS;
	vector<TCanvas*> cDetXvsDetY3DMeanChargeQuarterCellGrading;
	vector<TH2D*> hDetXvsDetY3DMeanChargeQuarterCellGrading;
	vector<TH1F*> hDetXvsDetY3DMeanChargeQuarterCellGradingLandau;
	TCanvas* cDetXvsDetY3DMeanChargeHighlightedQuarters;
	TCanvas* c3DdetQuarterCellClusterSize;
	TH2D* hQuarterCellsMeanClusterSize;
	TCanvas* hOverview;
	TH2D* hDetXvsDetY3DOverview;
	TCanvas* cCellNumbering;
	TH2D* hCellNumbering;
	TCanvas* c3DdetDeltaXChannel;
	TH2D* h3DdetDeltaXChannel;
	vector<TH1F*> hCellsDeltaX;
	vector<TH1F*> hCellsDeltaXQuarterCellGrading;
	TCanvas* c3DdetQuarterCellFluctuation1;
	TH2D* h3DdetQuarterCellFluctuation;
	TCanvas* c3DdetQuarterCellFluctuation;
	TH2D* h3DdetQuarterCellFluctuation1;
	TCanvas* c3DdetDeltaXChannelAbove1000;
	TH2D* h3DdetDeltaXChannelAbove1000;
	vector<TH2D*> hCellsCharge;
	vector<TH2D*> hCellsEvents;
	vector<TH2D*> hCellsEventsCheck;
	TH1F* hTransparentCharge3D;
	vector<TH1F*> hCellTransparentLandau;
	vector<TH1F*> hQuarterCellGradedTransparentLandau;
	vector<TCanvas*> cCellsTransparentHitPositionCellGraded;
	vector<TH2D*> hCellsTransparentHitPosition;
	vector<TH2D*> hCellsTransparentHitPositionCellGraded;
	vector<TH2D*> hCellsChargeBinAlignment[9];
	vector<TH2D*> hCellsEventsBinAlignment[9];
	TCanvas* cCellsOverlayed;
	TH2D* hCellsOverlayedCharge;
	TH2D* hCellsOverlayedEvents;
	TH2D* hCellsOverlayedMeanCharge;
	vector<TCanvas*> hCellsOverlayedBinAlignmentCanvas;
	vector<TH2D*> hCellsOverlayedChargeBinAlignment;
	vector<TH2D*> hCellsOverlayedEventsBinAlignment;
	vector<TH2D*> hCellsOverlayedMeanChargeBinAlignment;
	vector<TCanvas*> cCellsOverlayedBinAlignment;
	vector<TH2D*> hCellsOverlayedChargeBinAlignment1;
	vector<TH2D*> hCellsOverlayedEventsBinAlignment1;
	vector<TH2D*> hCellsOverlayedMeanChargeBinAlignment1;
	vector<TH1F*> hCellsLandau;
	vector<TH1F*> hCellsClusteSize;
	TCanvas* c2DClusterSizeCanvas;
	TH2D* h2DClusterSize;
	vector<TH1F*> hCellsLandauGraded;
	//hCellsGoodandBad
	TH1F* hCellsHarris18Good;
	TH1F* hCellsHarris10Bad;
	TH1F* hCellsAlexAllQuarters;
	TCanvas* cCellsLandau2D;
	TCanvas* cCellsLandau2DHighlightedQuarters;
	TH2D* hCellsLandau2D;
	TH2D* hCellsLandau2DQuarterFail;
	//Removed Columns
	TH2D* hCellsColumnCheck55;
	TH2D* hCellsColumnCheck1010;
	vector<TH1F*> hCellsOverlayBinSpec55;
	vector<TH1F*> hCellsOverlayBinSpec1010;
	TH2D* hCellsOverlayed55RMS;
	TH2D* hCellsOverlayed1010RMS;
	TH2D* hCellsOverlayed1010Significance;
	TF1* Landau;
	TCanvas* cCellsEventsNoColumn;
	TH2D* hCellsOverlayedEventsNoColumns;
	TH1F* hCellsOverlayedEntriesNoColumns;
	TH1F* hCellsOverlayedLandauNoColumn;
	TH1F* hCellsOverlayedColumnLandau;
	vector<TH2D*> hCellsEventsNoColumn;
	vector<TH1F*> hCellsLandauNoColumn;
	vector<TH1F*> hCellsLandauGradedNoColumn;
	TH1F* hBinnedMeanCharge;
	//xEdge
	TH1F* hEdgeCharge;
	TH1F* hEdgeChargeEvents;
	TCanvas* cEdgeMeanCharge;
	TH1F* hEdgeMeanCharge;
	//yEdge
	TH1F* hyEdgeCharge;
	TH1F* hyEdgeChargeEvents;
	TCanvas* cyEdgeMeanCharge;
	TH1F* hyEdgeMeanCharge;
	//DeadCellProfile
	TCanvas* cDeadCellMeanCharge;
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
	vector< vector<Float_t>* > vecPHDiamondHit, vecXPredicted, vecYPredicted, vecClusterFrequency,vecXPreditedDetector,vecYPreditedDetector;
	int HitCount, SeedCount;

	vector<int*> Detector;
	vector<TFiducialCut*> FidCut;
	TFidCutRegions fidCuts;
	TPositionPrediction *predictedPosition;

	Int_t verbosity;


private:
};

#endif /* TANALYSISOF3DDIAMONDS_HH_ */
