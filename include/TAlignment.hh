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

class TAlignment {
public:
	TAlignment(int runNumber);
	virtual ~TAlignment();
	int Align();
	void createVectors(UInt_t nEvents);
private:
	void createVectors();
	void initialiseHistos();
	void saveHistos();
	void addEventToTracks();
	void doDetAlignmentStep();
	void doDiaAlignmentStep();
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
