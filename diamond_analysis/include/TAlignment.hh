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
public:
	enum enumDetectorsToAlign {diaAlignment,silAlignment,bothAlignment};
public:
	TAlignment(TSettings* settings,TSettings::alignmentMode mode = TSettings::normalMode);
	virtual ~TAlignment();
	int Align(UInt_t nEvents=0,UInt_t startEvent=0,enumDetectorsToAlign detAlign=bothAlignment);
	int AlignSilicon(UInt_t nEvents=0,UInt_t startEvent=0);
	int AlignDiamond(UInt_t nEvents=0,UInt_t startEvent=0);
	void createEventVectors(UInt_t nEvents=0, UInt_t startEvent=0);
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
	void UpdateResolutions();
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

	TResidual getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool plot=false,TResidual res=TResidual(true),TCluster::calculationMode_t mode=TCluster::highest2Centroid,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10);
	TResidual getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool plot=false,TResidual res=TResidual(true),TCluster::calculationMode_t mode=TCluster::highest2Centroid,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10);
	TResidual getStripResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bAlign=false,bool plot=false,TResidual res=TResidual(true),TCluster::calculationMode_t mode=TCluster::maxValue,resCalcMode calcMode = normalCalcMode,Float_t maxChi2=10);
//	TResidual getStripResidualChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bAlign=false,bool plot=false,Float_t maxChi2=1000,TCluster::calculationMode_t mode=TCluster::maxValue);

//	TResidual calculateResidual(TPlaneProperties::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>* yPred,vector<Float_t>* deltaY) {return calculateResidual(cor,xPred,deltaX,yPred,deltaY,TResidual(true));};
//	TResidual calculateResidual(TPlaneProperties::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>* yPred,vector<Float_t>* deltaY,TResidual res);

//	TResidual calculateResidualWithChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> refPlane,Float_t maxChi2=10,bool bAlign=false,bool plot=false);
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
	Int_t nAlignmentStep;
    Int_t nAlignSteps;

    Int_t nDiaAlignmentStep;
    Int_t nDiaAlignSteps;
    TCluster::calculationMode_t diaCalcMode;
    TCluster::calculationMode_t silCalcMode;
    THTMLAlignment *htmlAlign;
    bool bPlotAll;
    TSettings::alignmentMode mode;
private:
    TResidual resPlane1,resPlane2,resPlane3;
    TResults* results;
	std::vector<TResidual> vecRes103;
	vector<Float_t> vecXPred;
	vector<Float_t> vecYPred;
	vector<Float_t> vecXObs;
	vector<Float_t> vecYObs;
	vector<Float_t> vecXDelta;
	vector<Float_t> vecYDelta;
	vector <Float_t> vecXPullDist;
	vector <Float_t> vecYPullDist;
	vector<Float_t> vecXMeasured;
	vector<Float_t> vecYMeasured;
	vector<Float_t> vecXChi2;
	vector<Float_t> vecYChi2;
	vector<Float_t> vecXPhi;
	vector<Float_t> vecYPhi;
	vector<Float_t> vecXResPrediction;
	vector<Float_t> vecYResPrediction;
	vector<Float_t> vecClusterSize;
	vector<pair<Float_t,Float_t> > pullValuesX;
	vector<pair<Float_t,Float_t> > pullValuesY;
	vector<Float_t > trackResValuesX;
	vector<Float_t > trackResValuesY;
	vector<pair<Float_t,Float_t> > gausFitValuesX;
	vector<pair<Float_t,Float_t> > gausFitValuesY;
	TCluster::calculationMode_t getClusterCalcMode(Int_t plane){return TPlaneProperties::isDiamondPlane(plane)?diaCalcMode:silCalcMode;}
};

#endif /* TALIGNMENT_HH_ */
