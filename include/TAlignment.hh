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

	struct TResiduum{ Float_t resXMean, resXSigma,resYMean,resYSigma;Float_t sumRx;
	Float_t sumRy;
	Float_t sumVx;
	Float_t sumVy;
	Float_t sumV2x;
	Float_t sumV2y;
	Float_t sumVRx;
	Float_t sumVRy;
	UInt_t nUsedTracks;};

class TAlignment {
public:
	TAlignment(TSettings* settings);
	virtual ~TAlignment();
	int Align();
	void createVectors(UInt_t nEvents);
	void setSettings(TSettings* settings);
private:
	void createVectors();
	void initialiseHistos();
	void saveHistos();
	void addEventToTracks();
	void doDetAlignmentStep();
	void doDiaAlignmentStep();
	void AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	void AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	TResiduum getResiduum(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2);
	TResiduum calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY);
	TResiduum calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY,TResiduum res);
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
private:

	TH2F* hScatterPosition[4];
	TH1F* hXPositionDifference[3];
	TH2F* hXYPositionDifference[3];
	TH2F* hXXPositionDifference[3];
	TH1F* hYPositionDifference[3];
	TH2F* hYYPositionDifference[3];
	TH2F* hYXPositionDifference[3];
};

#endif /* TALIGNMENT_HH_ */
