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
	void setEvent(TEvent *event){this->event = event;};
	Float_t getPosition(TPlane::enumCoordinate cor,UInt_t plane);
	Float_t getXPosition(UInt_t plane);
	Float_t getYPosition(UInt_t plane);
	Float_t getZPosition(UInt_t plane);
	TPositionPrediction predictPosition(UInt_t subjectPlane,vector<UInt_t> vecRefPlanes);
	vector<Float_t> getSiXPositions();
	vector<Float_t> getSiYPositions();
	TEvent* getEvent(){return event;}
	void setDetectorAlignment(TDetectorAlignment *alignment);
	
private:
//	void setPositions(TEvent event);
	Float_t getXOffset(int plane){return this->alignment->GetXOffset(plane);}; // TODO: get offsets
	Float_t getYOffset(int plane){return this->alignment->GetYOffset(plane);};
	Float_t getPhiXOffset(int plane){return this->alignment->GetPhiXOffset(plane);};
	Float_t getPhiYOffset(int plane){return this->alignment->GetPhiYOffset(plane);};
	
	TEvent *event;
	TDetectorAlignment *alignment;
	UInt_t verbosity;
	
};

#endif
