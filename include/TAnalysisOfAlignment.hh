/*
 * TAnalysisOfAlignment.hh
 *
 *  Created on: Mar 2, 2012
 *      Author: bachmair
 */

#ifndef TANALYSISOFALIGNMENT_HH_
#define TANALYSISOFALIGNMENT_HH_
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
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "TTracking.hh"
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"

#include "TADCEventReader.hh"
class TAnalysisOfAlignment {
public:
	TAnalysisOfAlignment(TSettings *settings);
	virtual ~TAnalysisOfAlignment();
	void setSettings(TSettings *settings);
	void doAnalysis(UInt_t nEvents=0);
private:
	UInt_t nEvent;
	void DoEtaCorrection(UInt_t correctionStep);
	void initialiseHistos();
	void chi2Distribution();
	void saveHistos();
	TTracking* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    TSettings* settings;
    UInt_t verbosity;
    TH1F *hChi2DistributionX;
    TH1F *hChi2DistributionY;
    TH1F *hChi2DistributionXY;
    TH2F *hChi2DistributionXY2D;
    TH2F *hAngleDistribution;
    TFile* correctedEtaFile;
};

#endif /* TANALYSISOFALIGNMENT_HH_ */
