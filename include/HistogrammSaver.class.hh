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
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TFile.h"
#include "TImage.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TSystem.h"

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
    void SetPlotsPath(std::string Path);
    std::string GetPlotsPath(){return plots_path;}
    void SetStyle(TStyle newStyle);
    void SetDuckStyle();

    static void SaveCanvasPNG(TCanvas *canvas, std::string location, std::string file_name);
    static void SaveCanvasC(TCanvas *canvas, std::string location, std::string file_name);
    static void SaveCanvasRoot(TCanvas *canvas, std::string location, std::string file_name);
private:
    unsigned int verbosity;
    TPaveText *pt;
    TDatime dateandtime;
    std::string plots_path;
    unsigned int runNumber;
    unsigned int nEvents;
    void UpdatePaveText();
    TStyle *currentStyle;
    TSystem *sys;
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
