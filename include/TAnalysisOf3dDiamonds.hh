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
#include "TCutG.h"
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
	bool eventValid();
	void initialiseShortAnalysisHistos();
	void createTreeTestHistos();
	void saveShortAnalysisHistos();
	void ShortAnalysis();

	void initialiseTransparentAnalysisHistos();
	void TransparentAnalysis();
	void saveTransparentAnalysisHistos();
	Float_t getTransparentCharge(Int_t nDiamondPattern, Int_t nChannelHit);

	float* VectorToArray(vector<float> nvector);

	//YAlignment
	void initialise3DGridReference();
	void initialise3DYAlignmentHistos();
	void initialise3DOverviewHistos();
	void initialise3D2DLandauAndClustersizeHistos();
	void initialise3DCellOverlayHistos();
	void saveLongAnalysisHistos();
	void LongAnalysis();
	float* SortArrayBtoS(float* nArray, int nSize);
	void printArray(float* nArray, int nSize, const std::string& space);

	//new functions
	void HitandSeedCount(TCluster* nCluster);
	void ClusterPlots(int nClusters, float nfiducialValueX, float nfiducialValueY);
	void RemoveLumpyClusters(TCluster* nCluster);
	int RemoveEdgeHits(TCluster* nCluster, pair<int,int> nDetector);
	int RemoveClustersWithHitOutside3D(TCluster* nCluster);
	int RemoveEdgeClusters(TCluster* nCluster,  int nDetector);

private:
	TCellAnalysisClass* clusteredAnalysis;
	vector<float> SortArrayPointer;
	Float_t fiducialValueX, fiducialValueY, chi2x,chi2y,xPredicted,yPredicted,xPredDet,yPredDet;

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
	TCanvas* hCombinedMeanCharge;    //TCanvas for MeanCharge All detectors
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
	TCanvas* hGridReferenceCanvas;
	TH2D* hGridReferenceDetSpace;
	TH2D* hGridReferenceCellSpace;
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
	TCanvas* cHarrisGoodandStripNormailsed;
	TCanvas* cGoodGradedandStripNormailsed;
	TCanvas* cStripFidCutXFidCutYvsMeanCharge;
	TH2D* hStripFidCutXFidCutYvsCharge;
	TH2D* hStripFidCutXFidCutYvsEvents;
	TH2D* hStripFidCutXFidCutYvsMeanCharge;
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
	TH1F* hBinnedMeanCharge;

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
	//vector<float> FidYAlignment;


	TSettings *settings;
	HistogrammSaver *histSaver;
	TTracking* eventReader;
	UInt_t nEvent;
	vector<Float_t> vecXPredictedDiamondHit,vecYPredictedDiamondHit,vecPHDiamondHitStrip,vecPHDiamondHit3dNoHoles,vecPHDiamondHit3dWithHoles,vecClusterSize,vecChi2Y,vecChi2X,vecClusterSeedSize,vecXPredictedStrip,vecYPredictedStrip;
	vector< vector<Float_t>* > vecPHDiamondHit, vecXPredicted, vecYPredicted, vecClusterFrequency,vecXPreditedDetector,vecYPreditedDetector;
	vector< vector <Float_t> > vecEdgePredX,vecEdgePredY,vecEdgePulseHeight;
	int HitCount, SeedCount;

	TPositionPrediction *predictedPosition;

	Int_t verbosity;
	vector<UInt_t> vecSilPlanes;
	UInt_t subjectDetector;
	UInt_t subjectPlane;

private:
};

#endif /* TANALYSISOF3DDIAMONDS_HH_ */