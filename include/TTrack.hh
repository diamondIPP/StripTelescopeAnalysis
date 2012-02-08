//
//  TTrack.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TTRACK_HH_
#define TTRACK_HH_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstring>
#include <deque>


//ROOT Class Headers
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h" // for adding your own classes to ROOT's library
#include "TStyle.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLinearFitter.h"
//#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "FidCutRegion.hh"

#include "TSettings.class.hh"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TDetectorAlignment.hh"
#include "TEvent.hh"
#include "TPositionPrediction.hh"

class TTrack {
public:

	TTrack(TDetectorAlignment *alignment);
	virtual ~TTrack();
	UInt_t getNClusters(int det);
//	bool  isValidSiliconTrack(){event.isValidSiliconTrack();};
	void setEvent(TEvent *newEvent){this->event = newEvent;};
	Float_t getStripXPositionOfCluster(UInt_t plane,TCluster xCluster, Float_t yPred,TCluster::calculationMode_t mode);
	Float_t getStripXPosition(UInt_t plane,Float_t yPred,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getPositionOfCluster(TPlane::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster, TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getPosition(TPlane::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getXPosition(UInt_t plane);
	Float_t getYPosition(UInt_t plane);
	Float_t getZPosition(UInt_t plane);
	Float_t getXMeasured(UInt_t plane);
	Float_t getYMeasured(UInt_t plane);
	Float_t getMeasured(TPlane::enumCoordinate cor,UInt_t plane);
	TPositionPrediction* predictPosition(UInt_t subjectPlane,vector<UInt_t> vecRefPlanes,bool bPrint=false);
	vector<Float_t> getSiXPositions();
	vector<Float_t> getSiYPositions();
	TEvent* getEvent(){return event;}
	void setDetectorAlignment(TDetectorAlignment *alignment);
    UInt_t getVerbosity() const;
    void setVerbosity(UInt_t verbosity = 0){this->verbosity = verbosity;}
	UInt_t getRawChannelNumber(UInt_t det, Float_t xPred, Float_t yPred); // returns the raw channel number for a x,y position in lab system
	Float_t getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred);
	Float_t getXPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred);
	Float_t getYPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred);

    ;
private:
    //	void setPositions(TEvent event);
    Float_t getXOffset(int plane)
    {
        return this->alignment->GetXOffset(plane);
    }

    ; // TODO: get offsets
    Float_t getYOffset(int plane)
    {
        return this->alignment->GetYOffset(plane);
    }

    ;
    Float_t getPhiXOffset(int plane)
    {
        return this->alignment->GetPhiXOffset(plane);
    }

    ;
    Float_t getPhiYOffset(int plane)
    {
        return this->alignment->GetPhiYOffset(plane);
    }

    ;
    TEvent *event;
    TDetectorAlignment *alignment;
    TPositionPrediction* predictedPosition;
    UInt_t verbosity;
    TLinearFitter *linFitX;
    TLinearFitter *linFitY;
};

#endif
