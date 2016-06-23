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
#include <algorithm>    // std::sort
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include <algorithm>    // std::min_element, std::max_element

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTML3DAnalysis.hh"
#include "TAnalysisOfAnalysisDifferences.hh"
#include "TAnalysisOf3DResolution.hh"
#include "LandauGaussFit.hh"

#include "TADCEventReader.hh"
#include "TTracking.hh"
#include "TSettings.class.hh"
#include "TCluster.hh"
#include "TCellAnalysisClass.hh"
#include "TTransparentAnalysis.hh"
#include "THStack.h"
#include "TPaletteAxis.h"

using namespace std;

class TAnalysisOf3dDiamonds {
public:
	TAnalysisOf3dDiamonds(TSettings *settings);
	virtual ~TAnalysisOf3dDiamonds();
	void	doAnalysis(UInt_t nEvents=0);

private:
	void PrintPositions();
	bool eventValid();
	void createTreeTestHistos();
	void fillClusterDistributions();

	void initialiseHistos();
	void initialiseLongAnalysisHistos();
	void initialiseShortAnalysisHistos();
	void initialiseTransparentAnalysisHistos();
	void InitialiseStripAnalysisHistos();

	void saveHistos();
	void saveGlobalHistos();
	void SaveShortAnalysisHistos();
	void SaveLongAnalysisHistos();
	void SaveLongAnalysisHistos2();
	void saveTransparentAnalysisHistos();
	void SaveStripAnalysisHistos();

	void ShortAnalysis();
	void LongAnalysis();
	bool TransparentAnalysis();
	void StripAnalysis();

    TCluster *diamondCluster;
    TCluster transparentCluster;
    TCluster clusteredCluster;
    bool validClusteredAnalysis;
    bool validTransparentAnalysis;
    void LongAnalysis_checkClusteredAnalysis();
    void LongAnalysis_checkTransparentAnalysis();

	//YAlignment
	void initialise3DYAlignmentHistos();
	void initialise3DOverviewHistos();
	void initialise3DCellOverlayHistos();
	void initialise3DCellCentralColumnOverlayHistos();
	void initialise3DCellBiasColumnOverlayHistos();
	void initialise3DCellOverlayIndividualBinHistos();
	void initialise3DOffsetOverlayHistos();
	void initialise3DOffsetAlignmentOverlayHistos();
	void initialiseEdgeFreeHistos();
	//new functions
	void HitandSeedCount(TCluster* nCluster);
	void ClusterPlots(int nClusters, float nfiducialValueX, float nfiducialValueY);
	void RemoveLumpyClusters(TCluster* nCluster);
	int RemoveEdgeHits(TCluster* nCluster, pair<int,int> nDetector);
//	int RemoveClustersWithHitOutside3D(TCluster* nCluster);
//	int RemoveEdgeClusters(TCluster* nCluster,  int nDetector);
	vector<Float_t> LongAnalysis_GradeCellByQuarters(int quarterFailCriteriaTyp, vector<TH1F*> hQuarterLandaus);
private:
private:
	void ShortAnalysis_FillEdgeAlignmentHistos();
	void ShortAnalysis_Analyse1Cluster(UInt_t clusterNo=0);
	void ShortAnalysis_Analyse2Cluster();
	void ShortAnalysis_Save2ClusterPlots();
	void ShortAnalysis_FillEdgeDistributions(Float_t clusterCharge);
	void ShortAnalysis_SaveEdgeDistributions();
	void ShortAnalysis_FillMeanChargeVector(Float_t clusterCharge);
	void ShortAnalysis_SaveMeanChargeVector();

private:
	void MakeGhostCluster(TCluster *diamondCluster,Int_t clusterSize);
	void LongAnalysisSaveCellAndQuaterNumbering();
	void LongAnalysis_FillOverlayedHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,Float_t clusterCharge,Float_t ClusterSize);
	void LongAnalysis_Fill3DCellOverlayIndividualBinHistos(Float_t xRelPosDet,Float_t yRelPosDet,Float_t clusterCharge, Float_t ClusterSize);
	void LongAnalysis_FillOverlayCentralColumnHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,Float_t clusterCharge, Float_t ClusterSize, TCluster* diamondCluster);
	void LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet, Float_t ClusterSize, TCluster* diamondCluster);
	void LongAnalysis_FillOverlayBiasColumnHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,Float_t clusterCharge, Float_t ClusterSize, TCluster* diamondCluster);
	void LongAnalysis_Fill3DOffsetOverlayBiasColumnAlignment(Float_t xRelPosDet,Float_t yRelPosDet, Float_t clusterCharge, Float_t ClusterSize);
	void LongAnalysis_FillEdgeFreeHistos(Float_t xPredDet, Float_t yPredDet, Float_t charge);
	TAnalysisOf3DResolution* resolutionAnalysis = 0;// DA: initialize with 0
	void LongAnalysis_FillResolutionPlots();
	void LongAnalysis_InitResolutionPlots();
	void LongAnalysis_CreateResolutionPlots();
	void LongAnalysis_CreateResolutionPlots(vector<TH1F*> *vec,TString kind);
	void LongAnalysis_CreateTH2_CellPlots(vector<TH2F*> *vec,TString kind,TString prefix= "hResolution");
	void LongAnalysis_InitChargeSharingPlots();
	void LongAnalysis_FillChargeSharingPlots();
	void LongAnalysis_SaveChargeSharingPlots();
	void LongAnalysis_SaveRawPulseHeightPlots();
	void LongAnalysis_SaveSNRPerCell();
	void LongAnalysis_SaveGoodAndBadCellLandaus();
	void LongAnalysis_InitGoodCellsLandaus();
	void LongAnalysis_FillGoodCellsLandaus(Float_t charge);
	void LongAnalysis_SaveGoodCellsLandaus();
	void LongAnalysis_SaveDeadCellProfile();
	void LongAnalysis_SaveCellsOverlayMeanCharge();
	void LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge();
	void LongAnalysis_SaveCellsBiasColumnOverlayMeanCharge();
	void LongAnalysis_SaveCellsOverlayBiasColumnAndCentralColumnStack();
	void LongAnalysis_Save3DCellOverlayIndividualBinHistos();
	void LongAnalysis_Save3D3DOffsetOverlayBiasColumnAlignment();
	Float_t LongAnalysis_CalculateRMS(vector<Float_t>* nVector);
	void LongAnalysis_Fill2DCellHitsBelowCutRelative(TH1F* nHisto, Int_t Alignment);
	void LongAnalysis_SaveCellsOverlayOffsetMeanCharge();
	void LongAnalysis_SaveEdgeFreeHistos();
	void LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors(int quarterFailCriteriaTyp = 0);
	void LongAnalysis_SaveFailedQuarters();
	void LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail();
	void LongAnalysis_SaveCellsClusterSize2DVsGrading();
	void LongAnalysis_SaveQuarterCellsClusterSize2DVsGrading();
	void LongAnalysis_SaveMeanChargePlots();
	bool LongAnalysis_IsDeadCell(vector<TH1F*> nhQuarterCellsLandau, Float_t nThreshold);
	void LongAnalysis_FillRelativeAddedTransparentCharge();

	void LongAnalysis_InitialiseRelativeAddedTransparentCharge();
	void LongAnalysis_SaveRelativeAddedTransparentCharge();
	void LongAnalysis_CreateRelativeAddedTransparentChargeComparisonPlots();

	void LongAnalysis_CompareTransparentAndClusteredAnalysis_Maps();

	void DoMonteCarloOfAvrgChargePerBinInOverlay(TProfile2D* profOverlay,TH1F* hLandauOfOverlay);
private:
	bool isTransparentCluster;
	bool useCMN;
	TString appendix;
private:
    TCluster GhostCluster;
	TH1F* hLandauStrip;
	TH2F* hLandauStripNegativeCharges;
	TH1F* hLandauStripNegativeChargesFraction;

	TH2F* hLandauStripNegativeChargesClPos;
	TH2F* hLandauStripNegativeChargePosition;
	TH1F* hLandau3DWithColumns;
	TH1F* hLandau3DPhantom;
    TH1F* hLandau3DPhantomCentral;
	TProfile2D* hLandauStripFidCutXvsFidCutY;
	TH2F* hLandauStripFiducialPosition;
	TH2F* hLandau3DWithColumnsFidCutXvsFidCutY;
	TH2F* hLandau3DPhantomFidCutXvsFidCutY;
	TProfile2D* hNegativeChargeFieldWireFraction;
	TH1F* hNegativeChargeFieldWireSpectrumAll;
	TH1F* hNegativeChargeFieldWireSpectrumAllButBad;
	TH1F* hNegativeChargeFieldWireSpectrumGood;
    TH1F* hNegativeChargeFullCellSpectrumAll;
    TH1F* hNegativeChargeFullCellSpectrumAllButBad;
    TH1F* hNegativeChargeFullCellSpectrumGood;

	TH2D* hNegativeChargeFieldWirePositions;
	TH2F* hNegativeChargeFieldWirePositionsOverlay;
	TCellAnalysisClass* clusteredAnalysis;
	vector<float> SortArrayPointer;
	Float_t fiducialValueX, fiducialValueY, chi2x,chi2y,xPredicted,yPredicted,xPredDet,yPredDet;
	Int_t PulseHeightBins, PulseHeightMin, PulseHeightMax,PulseHeightMaxMeanCharge,PulseHeightMinMeanCharge;

	string FileNameEnd;
	vector<TH2F*> hPHvsPredictedChannel, hPHvsChannel, hPHvsPredictedXPos, hPredictedPositionDiamondHit, hHitandSeedCount, hChi2XChi2Y, hFidCutXvsFidCutY;
	vector<TH1F*>hEventsvsChannel;
	TH1F* hEventsvsChannelCombined;
	TH1F* hNumberofClusters;
	TH1F* hDoubleClusterPos;
	TCanvas* cDoubleCluster;
	TH1F* hDoubleClusterPos0;
	TH1F* hDoubleClusterPos1;
	TH1F* hLandauCluster1;
	TH1F* hLandauCluster2;
	TH1F* hLandauDoubleCombined;
	TProfile2D* hFidCutsVsMeanCharge;

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
	vector<TH2D*> hFidCutXvsFidCutYClusters;
	TH2D* hEfficiency;
	//TH3F* hFidCutXvsFidCutYvsCharge;
	vector<TH1F*> hLandau;
	vector<TH1F*> hLandauTransparent;
	vector<TH1F*> hLandauTransparentBadCellsRemoved;
	TH2D* hCellOverlayvsCharge;
	TH2D* hCellOverlayvsEvents;
	TH2D* hCellOverlayvsMeanCharge;
	TCanvas* cCellOverlayvsMeanCharge;

	TH2F* hPHvsChannelStrip;
	TH2F* hPHvsPredictedXPosStrip;
	TH2F* hPHvsChannel3dNoHoles;
	TH2F* hPHvsPredictedXPos3dNoHoles;
	TH2F* hPHvsChannel3dWithHoles;
	TH2F* hPHvsPredictedXPos3dWithHoles;

	//For Fiducial Cut Boundaries Histograms
	vector <TBox*> FidCutChannelTBox;
	TCanvas* cPredictedEvents;
	TCanvas* cSeenEvents;
	vector<TCanvas*> cClusters;
	float FiducialMetric;			// Which Fiducial cut to apply
	float FiducialChannel;
	TCanvas* FidCudBoundMetricCanvas;
	TCanvas* cFidCudBoundChannel;
	TH2F* FidCudBoundMetric;
	TH2F* FidCudBoundChannel;

	/////////////
	//For YAlignment Histos
	/////////////
	TCanvas* cCombinedMeanChargeYAlignment;
	TCanvas* c3DdetMeanCharge;

	TCanvas* c3DdetMeanChargeWithMeanClusterSize;
	TH2D* hCellsMeanClusteSize;
	TCanvas* c3DdetMeanChargeRebinned;

	TCanvas* cRebinnedQuarterCellFails;
	TH2D* RebinnedQuarterCellFails;

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

	vector<TCanvas*> cDetXvsDetY3DMeanChargeQuarterCellGrading;
	vector<TH2D*> hDetXvsDetY3DMeanChargeQuarterCellGrading;
	vector<TH1F*> hDetXvsDetY3DMeanChargeQuarterCellGradingLandau;
	TCanvas* cDetXvsDetY3DMeanChargeHighlightedQuarters;
	TCanvas* c3DdetQuarterCellClusterSize;
	TH2D* hQuarterCellsMeanClusterSize;
	TCanvas* hOverview;
	TH2D* hDetXvsDetY3DOverview;
	TCanvas* cCellNumbering;
//	TH2D* hCellNumbering;
	TCanvas* c3DdetDeltaXChannel;
	TH2D* h3DdetDeltaXChannel;
	vector<TH1F*> hCellsDeltaX;
	vector<TH1F*> hCellsDeltaXQuarterCellGrading;
	TCanvas* c3DdetQuarterCellFluctuation1;
	TH2D* h3DdetQuarterCellFluctuation;
	TCanvas* c3DdetQuarterCellFluctuation;
	TH2D* h3DdetQuarterCellFluctuation1;
	//TCanvas* c3DdetDeltaXChannelAbove1000;
	//TH2D* h3DdetDeltaXChannelAbove1000;
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

	vector<TCanvas*> hCellsOverlayedBinAlignmentCanvas;
	vector<TH2D*> hCellsOverlayedChargeBinAlignment;
	vector<TH2D*> hCellsOverlayedEventsBinAlignment;
	vector<TH2D*> hCellsOverlayedMeanChargeBinAlignment;
	vector<TCanvas*> cCellsOverlayedBinAlignment;
	vector<TH2D*> hCellsOverlayedChargeBinAlignment1;
	vector<TH2D*> hCellsOverlayedEventsBinAlignment1;
	vector<TH2D*> hCellsOverlayedMeanChargeBinAlignment1;

	//LongAnalysis

	TCanvas* cDetXvsDetY3DMeanCharge;
	vector<TProfile2D*> hPulseHeightVsDetectorHitPostionXY_trans;
	TProfile2D* hPulseHeightVsDetectorHitPostionXY;
	TProfile2D* hPulseHeightVsDetectorHitPostionXYGoodCells;
	TProfile2D* hPulseHeightVsCell;
	TH1F* hLandauGoodCellsWithoutEdges;
	TH1F* hLandauGoodCellsWithoutColumns;
	TH1F* hLandauGoodCells;
//	TH2D* hDetXvsDetY3DEvents;
	TH1F* hCellMeanCharge;
	TH2F* hValidEventsFiducialSpace;
	TH2F* hValidEventsDetSpace;
	TH2F* hClusterEventsDetSpace;

	TH2D* hNegativeChargePosition;
	TH1F* hNegativeChargeFraction;
	TH2D* hNegativeChargeRatio;
	TH2D* hNegativeChargeRatioMax;
	TH2D* hNegativeChargeRatioAbs;
	TProfile2D* hNegativeChargeRatioOverlay;

	//vector<TCanvas*> cDeadCellMeanCharge;
	vector<TProfile*> hDeadCellCharge;
	vector<TH2D*> hDeadCellPositions;

	vector<TH1F*> hCellsLandau;
	vector<TH1F*> hCellsClusteSize;

	vector< vector<TH1F*> > hQuarterCellsLandau;
	vector< vector<TH1F*> > hQuarterCellsClusterSize;
	vector< vector<Int_t> > vecQuarterCellsPassFail;
	vector< vector<Int_t> > vecQuarterCellsFluctuation;
	vector<Int_t> CellGrading;
	vector<Int_t> HighlightGrading;
	Float_t hLandauGoodCellsMean;

	//TCanvas* cCellsOverlayMeanCharge;
	vector<TProfile2D*> hCellsOverlayAvrgCharge;
	vector<TProfile2D*> hCellsOverlayAvrgChargeMinusBadCells;
    vector<TH1F*> hCellsLandauMinusBadCells;
	vector<TProfile2D*> hCellsOverlayAvrgChargeGoodCells;
	vector<TProfile2D*> hCellsOverlayAvrgChargeBadCells;
	vector<TProfile2D*> hCellsOverlayAvrgChargeNoColumnHit;
	vector<TProfile2D*> hCellsCentralColumnOverlayAvrgCharge;
	vector<TProfile2D*> hCellsCentralColumnOverlayAvrgChargeMinusBadCells;
	TProfile2D* hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells;
	vector<TH2F*> hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents;
	vector<TProfile2D*> hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut;
	vector<TProfile2D*> hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis;
	vector<TH1F* > hLandauMinusBadCellsOffsetAnalysis;
	vector<TProfile2D*> hCellsCentralColumnOverlayAvrgChargeGoodCells;
	vector<TProfile2D*> hCellsOffsetOverlayAvrgCharge;
	vector<TProfile2D*> hCellsOffsetOverlayAvrgChargeMinusBadCells;
	vector<TProfile2D*> hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment;
	vector<TProfile2D*> hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut;
	TProfile2D* hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells;
	vector<TProfile2D*> hCellsOffsetOverlayAvrgChargeGoodCells;
	TProfile2D* hPulseHeigthCentralRegion;
	TProfile2D* hPulseHeigthEdgeRegion;
	TH2F* hEventsCentralRegion;
	TH2F* hEventsEdgeRegion;

	TProfile2D* hCellsBiasColumnOverlayAvrgCharge;
	TProfile2D* hCellsBiasColumnOverlayAvrgChargeMinusBadCells;
	TProfile2D* hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells;
	TH1F* hCellsBiasColumnOverlayLandauMinusBadCells;
	TH2F* hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents;
	TProfile2D* hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut;

	vector<TH1F*> hCellOverlayWithColumnLandau;
	vector<TH1F*> hCellOverlayNoColumnLandau;
	vector<TH1F*> hCellOverlayColumnLandau;

	vector<TH1F*> hCellsCentralColumnOverlayLandau;
	vector<TH1F*> hCellsCentralColumnOverlayLandauMinusBadCells;
	vector<TH1F*> hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis;
	vector<TH1F*> hCellsCentralColumnOverlayLandauGoodCells;

	//TH1F* hCellsHarrisGood;
	//TH1F* hCellsHarrisBad;

	TCanvas* cCellsLandau2D;
	TCanvas* cCellsLandau2DHighlightedQuarters;
	TH2D* hCellsLandau2D;
	TH2D* hCellsLandau2DQuarterFail;

	//LongAnalysis End

	TCanvas* c2DClusterSizeCanvas;
	TH2D* h2DClusterSize;
	vector<TH1F*> hCellsLandauGraded;
	//hCellsGoodandBad

	TH1F* hCellsAlexAllQuarters;
	TCanvas* cHarrisGoodandStripNormailsed;
	TCanvas* cGoodGradedandStripNormailsed;
	TCanvas* cStripFidCutXFidCutYvsMeanCharge;
	vector<TH2D*> hXdetvsYdetvsCharge;
	vector<TH2D*> hXdetvsYdetvsEvents;
	vector<TH2D*> hXdetvsYdetvsMeanCharge;
	vector<TCanvas*> ptrCanvasXdetvsYdetMeanCharge;
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

	int* DeadCellsArrayPointer;
	TH1F* hDeadCellsProfileCharge;
	TH1F* hDeadCellsProfileEvents;
	TH1F* hDeadCellsProfileMeanCharge;
	vector<TH1F*> hDeadCellsCharge;
	vector<TH1F*> hDeadCellsEvents;
	//Fiducial Cut
	vector <TBox*> FidCutChannelYAlignmentTBox;
	//vector<float> FidYAlignment;


	TSettings *settings;
	THTML3DAnalysis* html3D;
	HistogrammSaver *histSaver;
	TTracking* eventReader;
	UInt_t nEvent;
	vector<Float_t> vecXPredictedDiamondHit,vecYPredictedDiamondHit,vecPHDiamondHitStrip,vecPHDiamondHit3dNoHoles,vecPHDiamondHit3dWithHoles,vecClusterSize,vecChi2Y,vecChi2X,vecClusterSeedSize,vecXPredictedStrip,vecYPredictedStrip;
	vector< vector<Float_t>* > vecPHDiamondHit, vecXPredicted, vecYPredicted, vecClusterFrequency,vecXPreditedDetector,vecYPreditedDetector;
	vector< vector <Float_t> > vecEdgePredX,vecEdgePredY,vecEdgePulseHeight;
	vector<TH2F*> vecHChargeSharing;
	vector<TH1F*> vecHResolutionPerCell_maxValue;
	vector<TH1F*> vecHResolutionPerCell_chargeWeighted;
	vector<TH1F*> vecHResolutionPerCell_highest2Centroid;
    vector<TH1F*> vecHResolutionPerCell_h2C_WithCut;
	vector<TH2F*> vecHResolutionPerCell_maxValue_vs_SNR;
	vector<TH2F*> vecHResolutionPerCell_chargeWeighted_vs_SNR;
	vector<TH2F*> vecHResolutionPerCell_highest2Centroid_vs_SNR;
    vector<TH2F*> vecHResolutionPerCell_h2C_WithCut_vs_SNR;
    vector<TH2F*> vecHResolutionPerCell_maxValue_vs_PredHit;
    vector<TH2F*> vecHResolutionPerCell_maxValue_vs_PredHitY;
    vector<TH2F*> vecHResolutionPerCell_chargeWeighted_vs_PredHit;
    vector<TH2F*> vecHResolutionPerCell_chargeWeighted_vs_PredHitY;
    vector<TH2F*> vecHResolutionPerCell_highest2Centroid_vs_PredHit;
    vector<TH2F*> vecHResolutionPerCell_highest2Centroid_vs_PredHitY;
    vector<TH2F*> vecHResolutionPerCell_h2C_WithCut_vs_PredHit;
    vector<TH2F*> vecHResolutionPerCell_h2C_WithCut_vs_PredHitY;
    TH2F* hAdjacentChannels_SNR;
    TH2F* hAdjacentSNR_vs_cellNo;
    TH2F* hAdjacentChannels_Signal;
private:
	TH2F* hLongAnalysisInvalidCellNo;
	TH2F* hLongAnalysisInvalidCluster;
	TH2D* hLongAnalysisQuarterFluctuations, MeanOfLandauGoodCells;
	TH2F* hShortAnalysis2ClusterHitPattern_1stCluster;
	TH2F* hShortAnalysis2ClusterHitPattern_2ndCluster;

	TH2F* hRelativeChargeTwoClustersX;
	TH2F* hRelativeChargeTwoClustersY;
	TProfile2D* hRelativeChargeTwoClustersXY;
	TProfile2D* hShortAnalysis2TotalChargeXY;
	TProfile* hRelatviveNumberOfMultipleClusterEvents;
	TProfile* hRelatviveNumberOfMultipleClusterEventsSamePattern;
	TH2F* hHitPositionNoCluster;
	TH2F* hHitPositionOneCluster;
	TH2F* hHitPositionMultiCluster;
	TH2F* hHitPositionTwoCluster;
	TH1F* hNClusters;

	TProfile2D* hTotalAvrgChargeXY;
	TProfile2D* hTotalAvrgChargeXY_electrons;
	vector<Int_t>vecDeadCells;
private:
	vector <Float_t> vecPredDetX_ShortAna,vecPredDetY_ShortAna,vecPulseHeight_ShortAna;
	vector <Float_t> vecPH_Cluster1_ShortAna,vecPH_Cluster2_ShortAna,vecCh_Cluster2_ShortAna,vecCh_Cluster1_ShortAna;
private:
	int HitCount, SeedCount;

	TPositionPrediction *predictedPosition;

	Int_t verbosity;
	vector<UInt_t> vecSilPlanes;
	UInt_t subjectDetector;
	UInt_t subjectPlane;

private:
	//TransparentAnalysis
	TH2F* hTransparentAnalysisInvalidCluster;
	TH2F* hTransparentAnalysisValidCluster;
	TH2F* hTransparentAnalysisValidClusterFidCutXvsFidCutY;

	vector<TH1F*> hTransparentAnalysisTransparentCharge;
	vector<TH1F*> hTransparentAnalysisTransparentChargeGoodCells;
	vector<TH1F*> hTransparentAnalysisTransparentChargeBadCells;
	vector<TProfile2D*> hTransparentAnalysisTransparentChargeProfile;

	vector<TH1F*> hTransparentAnalysisTransparentAddedCharge;
	vector<TH1F*> hTransparentAnalysisTransparentAddedChargeGoodCells;
	vector<TH1F*> hTransparentAnalysisTransparentAddedChargeBadCells;
	vector<TProfile2D*> hTransparentAnalysisTransparentAddedChargeProfile;

	vector<TH1F*> hTransparentAnalysisRelativeAddedCharge;
	vector<TH1F*> hTransparentAnalysisRelativeAddedChargeGoodCells;
	vector<TH1F*> hTransparentAnalysisRelativeAddedChargeBadCells;
	vector<TProfile2D*> hTransparentAnalysisRelativeAddedChargeProfile;

    vector<TH1F*> hTransparentAnalysisRelativeCharge;
    vector<TH1F*> hTransparentAnalysisRelativeChargeGoodCells;
    vector<TH1F*> hTransparentAnalysisRelativeChargeBadCells;
    vector<TProfile2D*> hTransparentAnalysisRelativeChargeProfile;

    vector<TH1F*> hTransparentAnalysisTransparentChargeWithoutEdge;
    vector<TH1F*> hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge;
    vector<TH1F*> hTransparentAnalysisTransparentChargeBadCellsWithoutEdge;
    vector<TProfile2D*> hTransparentAnalysisTransparentChargeProfileWithoutEdge;
    TH2D* hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells;
    TH2D* hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells;

    vector<TH1F*> hTransparentAnalysisTransparentAddedChargeWithoutEdge;
    vector<TH1F*> hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge;
    vector<TH1F*> hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge;
    vector<TProfile2D*> hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge;

    vector<TH1F*> hTransparentAnalysisRelativeAddedChargeWithoutEdge;
    vector<TH1F*> hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge;
    vector<TH1F*> hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge;
    vector<TProfile2D*> hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge;

    vector<TH1F*> hTransparentAnalysisTransparentRelativeHitChannel;
    vector<TH1F*> hTransparentAnalysisTransparentRelativeHitChannelGoodCells;
    vector<TH1F*> hTransparentAnalysisTransparentRelativeHitChannelBadCells;
    vector<TProfile2D*> hTransparentAnalysisTransparentRelativeHitChannelProfile;

    vector<TProfile2D*> VecOverlayCellBinHistos;
    vector<TH1F*> VecOverlayCellBinLandaus;
    TH1F* hOverlayCellBinHits;
    TH1F* hOverlayCellBinHitsBelowCut;
    TH1F* hOverlayCellOffsetBinHits;
    TH1F* hOverlayCellOffsetBinHitsBelowCut;
    vector<TH1F*> hOverlayCellOffsetAlignmentBinHits;
    vector<TH1F*> hOverlayCellOffsetAlignmentBinHitsBelowCut;
    TH1F* hOverlayCellUnEvenBinningBinHits;
    TH1F* hOverlayCellUnEvenBinningBinHitsBelowCut;

    vector<Float_t> ShiftX;
    vector<Float_t> ShiftY;


    UInt_t maxClusterSize3d;

    map<Int_t, TCluster> mapClusteredAnalysisGoodCells;
    map<Int_t, TCluster> mapTransparentAnalysisGoodCells;
    map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsGoodCells;


    map<Int_t, TCluster> mapClusteredAnalysisAllCells;
    map<Int_t, TCluster> mapTransparentAnalysisAllCells;
    map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsAllCells;

    map<Int_t, TCluster> mapClusteredAnalysisAllButBadCells;
    map<Int_t, TCluster> mapTransparentAnalysisAllButBadCells;
    map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsAllButBadCells;

    Float_t maxsnr;
};

#endif /* TANALYSISOF3DDIAMONDS_HH_ */
