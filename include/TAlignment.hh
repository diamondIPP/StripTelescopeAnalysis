/*
 * TAlignment.hh
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#ifndef TALIGNMENT_HH_
#define TALIGNMENT_HH_

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
#include "TH2F.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "THStack.h"


#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "TCluster.hh"
#include "TADCEventReader.hh"
#include "TDiamondTrack.hh"
#include "TSettings.class.hh"
#include "TDetectorAlignment.hh"
#include "TTrack.hh"
#include "TPlane.hh"
#include "TResults.hh"
#include "TResidual.hh"
#include "TTracking.hh"
#include "TTransparentAnalysis.hh"

#include "THTMLAlignment.hh"


/**
 * @brief alignment of all subdetectors
 *
 * @author Felix Bachmair
 * @date 12.1.2012 15:30		Last Change - Felix Bachmair
 *
 * Alignment of all stripe detectors
 * main function is Align()
 * it creates first a vector of events to make the whole alignment faster
 * afterwards it calculates the residuals before alignment
 * last step alignment of the planes.
 *
 *
 */
class TAlignment {
private:
	enum resCalcMode  { normalCalcMode,chi2CalcMode,resolutionCalcMode};
	enum resolutionUpdateMode {normalMode=1,resolutionMode=2,noUpdate=0};
	enum alignmentMode {singleStrip,doubleStrip};
public:
	enum enumDetectorsToAlign {diaAlignment=0,silAlignment=1,bothAlignment=2};
public:
	TAlignment(TSettings* settings,TSettings::alignmentMode mode = TSettings::normalMode);
	virtual ~TAlignment();
	int Align(UInt_t nEvents=0,UInt_t startEvent=0,enumDetectorsToAlign detAlign=bothAlignment);
	int AlignSilicon(UInt_t nEvents=0,UInt_t startEvent=0){ return Align(nEvents,startEvent,silAlignment);}
	int AlignDiamond(UInt_t nEvents=0,UInt_t startEvent=0){ return Align(nEvents,startEvent,diaAlignment);}
	void createEventVectors(UInt_t nEvents=0, UInt_t startEvent=0,enumDetectorsToAlign detAlign=bothAlignment);
	void createTransparentEventVectors(UInt_t nEvents=0, UInt_t startEvent=0);
	void setSettings(TSettings* settings);
	void PrintEvents(UInt_t maxEvent=0,UInt_t startEvent=0);
	void setVerbosity(UInt_t verb){verbosity=verb;};
	void setResults(TResults *res){results=res;}
private:
	void clearMeasuredVectors();
	void initialiseDetectorAlignment(TSettings::alignmentMode mode =TSettings::normalMode);
	void loadDetectorAlignment(TSettings::alignmentMode mode =TSettings::normalMode);
	void AlignSiliconPlanes();
	void doPreAlignment();
	bool siliconAlignmentStep(bool bPrint=false,bool bUpdateAlignment = true);
	void AlignDiamondPlane();
	bool checkPullDistributions();
	bool checkMaxDeviationDetectorPair(TPlaneProperties::enumCoordinate cor, UInt_t pl1, UInt_t pl2, Float_t maxDeviation = .1);
	bool checkDeviationInDifferentDetectorPairs(TPlaneProperties::enumCoordinate cor,Float_t expctedFactor=2./3.,  Float_t maxDeviation=.5);
	bool comparePullDistributions(TPlaneProperties::enumCoordinate cor, Float_t maxDeviation = .2);
	void saveAlignment(TSettings::alignmentMode mode = TSettings::normalMode);
	void getChi2Distribution(Float_t maxChi2=1000);
	void AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2){alignDetector(TPlaneProperties::XY_COR,subjectPlane,refPlane1,refPlane2);};
	void AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2){alignDetector(TPlaneProperties::X_COR,subjectPlane,refPlane1,refPlane2);};
	void AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2){alignDetector(TPlaneProperties::Y_COR,subjectPlane,refPlane1,refPlane2);};
	void DoEtaCorrectionSilicon(UInt_t correctionStep=0);
	void getFinalSiliconAlignmentResuluts();
	void UpdateResolutions(vector<Float_t> residuals, vector<Float_t> resolutions);
	void setSiliconDetectorResolution(Float_t maxChi2);

	void CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, string refPlaneString, bool bPlot, bool bUpdateResolution, bool bChi2=false);
//	void CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,string refPlaneString,bool bPlot=true, resolutionUpdateMode bUpdateResolution = noUpdate,bool bChi2=false);
	void SetResolutionsWithUserInput();
	void LoadResolutionFromSettingsFile();
	void inputResolution(UInt_t plane, TPlaneProperties::enumCoordinate cor);

	TResidual alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t>vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));
	TResidual alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=false,TResidual res=TResidual(true));
	TResidual alignStripDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t>vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));

	TResidual CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=true,TResidual  res=TResidual(true));
	TResidual CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t> vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));
	TResidual CheckStripDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t> vecRefPlanes, bool bAlign=false,bool bPlot=false,TResidual res=TResidual(true));
	TResidual CheckStripDetectorAlignmentChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2);

	TString GetReferencePlaneString( vector<UInt_t> *vecRefPlanes);
	TResidual Residual(alignmentMode aligning, TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bAlign=false, bool plot=false,TResidual res=TResidual(true),TCluster::calculationMode_t mode=TCluster::highest2Centroid,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10);
	TResidual getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=false,TResidual res=TResidual(true),TCluster::calculationMode_t mode=TCluster::highest2Centroid,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10);
	TResidual getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot=false,TResidual resOld=TResidual(true),TCluster::calculationMode_t mode=TCluster::highest2Centroid,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10){
        return Residual(doubleStrip,cor,subjectPlane,vecRefPlanes,false,bPlot,resOld,mode,calcMode,maxChi2);
    };
	/**
	 * @brief creates element TResidual to adjust the alignment
	 *
	 * creates a vector of pedicted X positions, predicted Y positions, delta X and delta Y
	 * and use the function calculateResidual to get the residual with this vectors
	 * @param   cor coordindate for which the residual is calculated
	 * @param   subjectPlane plane for which the residual should be calculated
	 * @param   refPlane1 first reference plane
	 * @param   refPlane2 second reference plane
	 * @param   bPlot   variable to create plots or not
	 *
	 * @return
	 */
	TResidual getStripResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,
	        bool bAlign=false,bool bPlot=false,TResidual resOld=TResidual(true),
	        TCluster::calculationMode_t mode=TCluster::maxValue,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10){
	    return Residual(singleStrip,cor,subjectPlane,vecRefPlanes,bAlign,bPlot,resOld,mode,calcMode,maxChi2);
	}


	TTracking* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TRandom rand;
    UInt_t nEvent;
    Float_t alignmentPercentage;
    UInt_t runNumber;
    TSettings *settings;

    TDetectorAlignment* align;
    Double_t detectorD0Z; // by definition
    Double_t detectorD1Z; // by definition
    Double_t detectorD2Z; // by definition
    Double_t detectorD3Z; // by definition
    Double_t detectorDiaZ; // by definition
    int verbosity;
    TTrack* myTrack;
    Float_t res_keep_factor;
	
	vector<TEvent> events;
	vector <Float_t> fiducialValueX;
    vector <Float_t> fiducialValueY;
	vector<bool> telescopeAlignmentEvent;
	Int_t nAlignmentStep;
    Int_t nAlignSteps;

    bool silicon_chi2_cut;
    Int_t nDiaAlignmentStep;
    Int_t nDiaAlignSteps;
    TCluster::calculationMode_t diaCalcMode;
    TCluster::calculationMode_t silCalcMode;
    THTMLAlignment *htmlAlign;
    bool bPlotAll;
    TSettings::alignmentMode mode;
private:
    TString GetPlotPreName(UInt_t subjectPlane);
    TString GetPlotPostName(bool bChi2);

    Float_t CreateSigmaOfPredictionXPlots(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot);
    void CreateDistributionPlotDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, Float_t xPredictionSigma);
    void CreateDistributionPlotDeltaY(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, Float_t yPredictionSigma);

    void CreateScatterPlotObsXvsObsY(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateAngularDistributionPlot(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);

    void CreateScatterPlotPredXvsDeltaY(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotPredYvsDeltaY(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotPredYvsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotPredXvsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterChi2vsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);

    void CreateScatterPlotMeasXvsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotPredXDetvsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotPredXChvsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateFidValueXVsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateFidValueYVsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    Float_t CreateSigmaOfPredictionYPlots(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotEtaVsDeltaX(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot);

    void CreateRelHitPosXPredDetMetricVsUseEventPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot);
    void CreateRelHitPosXMeasDetMetricVsUseEventPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot);
    void CreateRelHitPosMeasXPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot);
    void CreateRelHitPosPredXPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot);
    void CreateResidualVsAngle(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot);

    void CreateScatterPlotDeltaXvsChi2X(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);
    void CreateScatterPlotDeltaYvsChi2Y(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString,bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment);

    void CreateRelHitPosVsChi2Plots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString);
    void CreateChi2DistributionPlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString);


    TResidual resPlane1,resPlane2,resPlane3;
    TResults* results;
	std::vector<TResidual> vecRes103;
	vector<Float_t> vecXLabPredMetric;
	vector<Float_t> vecYLabPredMetric;
	vector<Float_t> vecXDetPredMetric,vecYDetPredMetric;
	vector<Float_t> vecXDetPredChannel;
	vector<Float_t> vecXLabMeasMetric;
	vector<Float_t> vecYLabMeasMetric;
	vector<Float_t> vecXLabDeltaMetric;
	vector<Float_t> vecYLabDeltaMetric;
	vector <Float_t> vecXPullDist;
	vector <Float_t> vecYPullDist;
	vector<Float_t> vecXDetMeasMetric;
	vector<Float_t> vecYDetMeasMetric;
	vector<Float_t> vecXChi2;
	vector<Float_t> vecYChi2;
	vector<Float_t> vecXPhi;
	vector<Float_t> vecYPhi;
	vector<Float_t> vecEta;
	vector<Float_t> vecXResPrediction;
	vector<Float_t> vecYResPrediction;
	vector<Float_t> vecClusterSize;
	vector<Float_t> vecXFidValue;
	vector<Float_t> vecYFidValue;
	vector<pair<Float_t,Float_t> > pullValuesX;
	vector<pair<Float_t,Float_t> > pullValuesY;
	vector<Float_t > trackResValuesX;
	vector<Float_t > trackResValuesY;
	vector<pair<Float_t,Float_t> > gausFitValuesX;
	vector<pair<Float_t,Float_t> > gausFitValuesY;
	vector<Float_t> vecXDetRelHitPosPredMetricAll;
	vector<Float_t> vecXDetRelHitPosMeasMetricAll;
    vector<Float_t> vecDeltaXMetricAll;
    vector<Float_t> vecDeltaYMetricAll;
	vector<Float_t> vecUsedEventAll;
	TCluster::calculationMode_t getClusterCalcMode(Int_t plane){return TPlaneProperties::isDiamondPlane(plane)?diaCalcMode:silCalcMode;}
};

#endif /* TALIGNMENT_HH_ */
