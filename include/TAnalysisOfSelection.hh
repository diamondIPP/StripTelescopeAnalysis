/*
 * TAnalysisOfSelection.hh
 *
 *  Created on: May 18, 2012
 *      Author: bachmair
 */

#ifndef TANALYSISOFSELECTION_HH_
#define TANALYSISOFSELECTION_HH_
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
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTMLPedestal.hh"
#include "LandauGaussFit.hh"
#include "THTMLLandaus.hh"

#include "TADCEventReader.hh"
#include "TSettings.class.hh"

class TAnalysisOfSelection {
public:
	TAnalysisOfSelection(TSettings *settings);
	virtual ~TAnalysisOfSelection();
	void	doAnalysis(UInt_t nEvents=0);
private:
	void initialiseHistos();
	void saveHistos();
	void analyseEvent();
private:
	TSettings *settings;
	HistogrammSaver *histSaver;
	TADCEventReader* eventReader;
	TSystem *sys;
	UInt_t nEvent;
	TH2F* histoLandauDistribution;
	THTMLLandaus *htmlLandau;
};

#endif /* TANALYSISOFSELECTION_HH_ */
