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
#include "TGraph.h"
#include "TPaveText.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TFile.h"
#include "TImage.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TF1.h"
#include <sys/dir.h>
#include <stdio.h>
#include <stdlib.h>

//#include <sys/dirent.h>
#include <sys/stat.h>

class HistogrammSaver {
public:
	enum EnumAxisRange{
		maxWidth,fiveSigma,threeSigma,positiveArea,positiveSigma,manual
	};
	HistogrammSaver(int verbosity=0);
	virtual ~HistogrammSaver();
    void SaveHistogram(TH1F* histo, bool fitGauss = 0);
    void SaveHistogram(TH2F* histo);
    void SaveGraph(TGraph* graph,std::string name,std::string option="AP");
    void SaveHistogramPNG(TH1F* histo);
    void SaveHistogramPNG(TH2F* histo);
    void SaveGraphPNG(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramFitGaussPNG(TH1F* histo);
    void SaveHistogramROOT(TH1F* histo);
    void SaveHistogramROOT(TH2F* histo);
    void SaveGraphROOT(TGraph* graph,std::string name,std::string option="AP");
    void SaveHistogramPDF(TH1F* histo);
    void SaveHistogramPDF(TH2F* histo);
    void SetVerbosity(unsigned int i);
    void SetRunNumber(unsigned int runNumber);
    void SetNumberOfEvents(unsigned int nEvents);
    void SetPlotsPath(std::string Path);
    std::string GetPlotsPath(){return plots_path;}
    void SetStyle(TStyle newStyle);
    void SetDuckStyle();
    void SetRange(Float_t min,Float_t max);

    static TH2F CreateScatterHisto(std::string name,std::vector<Float_t> posX, std::vector<Float_t> posY,UInt_t nBins=4096);
    static TGraph CreateDipendencyGraph(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos);
    static TH2F CreateDipendencyHisto(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos,UInt_t nBins=4096);
    static TH1F CreateDistributionHisto(std::string name, std::vector<Float_t> vec,UInt_t nBins=4096,EnumAxisRange range=maxWidth,Float_t xmin=-1,Float_t xmax=1);
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
    TFile* histoFile;
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
