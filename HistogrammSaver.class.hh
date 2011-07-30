/*
 * HistogrammSaver.class.hh
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 */

#ifndef HISTOGRAMMSAVER_CLASS_HH_
#define HISTOGRAMMSAVER_CLASS_HH_
//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TFile.h"

class HistogrammSaver {
public:
	HistogrammSaver(int verbosity=0);
	virtual ~HistogrammSaver();
    void SaveHistogram(TH1F* histo);
    void SaveHistogram(TH2F* histo);
    void SaveHistogramPNG(TH1F* histo);
    void SaveHistogramPNG(TH2F* histo);
    void SaveHistogramROOT(TH1F* histo);
    void SaveHistogramROOT(TH2F* histo);
    void SaveHistogramPDF(TH1F* histo);
    void SaveHistogramPDF(TH2F* histo);
    void SetVerbosity(unsigned int i);
    void SetRunNumber(unsigned int runNumber);
    void SetNumberOfEvents(unsigned int nEvents);
    void SetPlotsPath(string Path);
    void SetStyle(TStyle newStyle);
private:
    unsigned int verbosity;
    TPaveText *pt;
    TDatime dateandtime;
    string plots_path;
    unsigned int runNumber;
    unsigned int nEvents;
    void UpdatePaveText();
    TStyle currentStyle;
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
