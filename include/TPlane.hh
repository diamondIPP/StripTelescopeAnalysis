//
//  TPlane.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TPlane_hh
#define TPlane_hh

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
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TDetectorAlignment.hh"
#include "TCluster.hh"

//class TCluster;
class TPlane:public TObject {
public:
	TPlane();
	TPlane(vector<TCluster> xClusters, vector<TCluster> yClusters);
	virtual ~TPlane();
	Float_t getXPosition(int cl){return xClusters[cl].getPosition();};
	Float_t getYPosition(int cl){return yClusters[cl].getPosition();};
private:
	
	
	vector<TCluster> xClusters, yClusters;
    ClassDef(TPlane,1);
};

#endif // TPlane_hh
