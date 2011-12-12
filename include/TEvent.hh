//
//  TEvent.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TEVENT_HH_
#define TEVENT_HH_

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
#include "TObject.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "FidCutRegion.hh"

#include "TSettings.class.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TDetectorAlignment.hh"
#include "TPlane.hh"

class TEvent:public TObject {
public:
	TEvent(UInt_t nEvent=0);
	virtual ~TEvent();
	TPlane getPlane(int plane){return planes[plane];};
	void addPlane(TPlane plane,Int_t pos=-1);
	UInt_t getNXClusters(UInt_t plane){if(planes.size()>plane)return planes.at(plane).getNXClusters();else return 0;}
	UInt_t getNYClusters(UInt_t plane){if(planes.size()>plane)return planes.at(plane).getNXClusters();else return 0;}
	UInt_t getNPlanes();
	bool  isValidSiliconEvent();
	UInt_t getEventNumber(){return eventNumber;};
//	void setPositions(TDetectorAlignment align);
	
private:
	
	vector<TPlane> planes;
	UInt_t eventNumber;

    ClassDef(TEvent,1);
public:
	
};

#endif //TEVENT_HH_
