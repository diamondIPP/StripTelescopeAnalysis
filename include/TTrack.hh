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

class TTrack {
public:
	TTrack(TEvent event);
	virtual ~TTrack();
	UInt_t getNClusters(int det);
	bool  isValidSiliconTrack(){event.isValidSiliconTrack();};
	
private:
	void setPositions(TEvent event);
	Float_t getXOffset(int plane){return 0;}; // TODO: get offsets
	Float_t getYOffset(int plane){return 0;};
	Float_t getPhiXOffset(int plane){return 0;};
	Float_t getPhiYOffset(int plane){return 0;};
	
	
	
	vector<Float_t> xPositions, yPositions, zPositions;
	
public:
	
};

#endif
