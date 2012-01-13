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
 * @bug The alignment of the y Coordinate is not working 100%
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
	void AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	TResidual AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=false,TResidual res=TResidual(true));
	TResidual CheckDetectorAlignment(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot=true,TResidual  res=TResidual(true));
	TResidual getResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool plot=false,TResidual res=TResidual(true));
	TResidual calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY);
	TResidual calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY,TResidual res);
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TRandom rand;
    vector<TDiamondTrack> tracks;
    vector<bool> tracks_masked;
    vector<TDiamondTrack> tracks_fidcut;
    vector<bool> tracks_masked_fidcut;
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
    int nAlignSteps;
    TTrack* myTrack;
    Float_t res_keep_factor;
	
	vector<TEvent> events;
	Int_t alignmentSteps;
	Int_t nAlignmentStep;

private:
	std::vector<TResidual> vecRes103;
};

#endif /* TALIGNMENT_HH_ */
