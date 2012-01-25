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


#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "Track.class.hh"
#include "Cluster.class.hh"
#include "TADCEventReader.hh"
#include "TDiamondTrack.hh"
#include "AlignmentClass.hh"
#include "TSettings.class.hh"
#include "TDetectorAlignment.hh"
#include "TTrack.hh"
#include "TPlane.hh"
#include "TResidual.hh"


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
public:
	TAlignment(TSettings* settings);
	virtual ~TAlignment();
	int Align(UInt_t nEvents=0,UInt_t startEvent=0);
	void createEventVectors(UInt_t nEvents=0, UInt_t startEvent=0);
	void setSettings(TSettings* settings);
	void PrintEvents(UInt_t maxEvent=0,UInt_t startEvent=0);
private:
	void AlignSiliconPlanes();
	void AlignDiamondPlane();
	void saveAlignment();
	void getChi2Distribution(Float_t maxChi2=1000);
	void AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void getFinalSiliconAlignmentResuluts();
	void setSiliconDetectorResolution(Float_t maxChi2);
	void CreatePlots(TPlane::enumCoordinate cor, UInt_t subjectPlane,string refPlaneString,bool bPlot=true, bool bUpdateAlignment=false);
	TResidual AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t>vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));
	TResidual AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=false,TResidual res=TResidual(true));
	TResidual AlignStripDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t>vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));

	TResidual CheckDetectorAlignment(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=true,TResidual  res=TResidual(true));
	TResidual CheckDetectorAlignment(TPlane::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t> vecRefPlanes,bool bPlot=false,TResidual res=TResidual(true));
	TResidual CheckStripDetectorAlignment(TPlane::enumCoordinate cor, UInt_t subjectPlane,vector<UInt_t> vecRefPlanes, bool bAlign=false,bool bPlot=false,TResidual res=TResidual(true));

	TResidual getResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool plot=false,TResidual res=TResidual(true));
	TResidual getResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool plot=false,TResidual res=TResidual(true));
	TResidual getStripResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bAlign=false,bool plot=false,TResidual res=TResidual(true));

	TResidual calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>* yPred,vector<Float_t>* deltaY);
	TResidual calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>* yPred,vector<Float_t>* deltaY,TResidual res);

	TResidual calculateResidualWithChi2(TPlane::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> refPlane,Float_t maxChi2=10,bool bAlign=false,bool plot=false);
	TADCEventReader* eventReader;
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
private:
	std::vector<TResidual> vecRes103;
	vector<Float_t> vecXPred;
	vector<Float_t> vecYPred;
	vector<Float_t> vecXObs;
	vector<Float_t> vecYObs;
	vector<Float_t> vecDeltaX;
	vector<Float_t> vecDeltaY;
	vector<Float_t> vecXChi2;
	vector<Float_t> vecYChi2;

};

#endif /* TALIGNMENT_HH_ */
