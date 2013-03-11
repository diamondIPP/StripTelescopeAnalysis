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
#include "TH3F.h"
#include "TH3.h"
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
#include "TLegend.h"
#include "TGaxis.h"
#include <sys/dir.h>
#include <stdio.h>
#include <stdlib.h>
#include "TGraphErrors.h"
#include "TString.h"
#include <limits>

//#include <sys/dirent.h>
#include <sys/stat.h>

class HistogrammSaver {
private:
	unsigned int verbosity;
public:

	enum EnumAxisRange{
		maxWidth,fiveSigma,threeSigma,positiveArea,positiveSigma,manual
	};
	HistogrammSaver(int verbosity=0);
	virtual ~HistogrammSaver();
	void SetOptStat(std::string optStat){gStyle->SetOptStat(optStat.c_str());}
	void SetOptStat(Int_t stat){gStyle->SetOptStat(stat);}
	void SetOptFit(Int_t fitOpt){gStyle->SetOptFit(fitOpt);}
	//void SetOptFit(std::string fitOpt){gStyle->SetOptFit(fitOpt.c_str());}
	void SaveCanvas(TCanvas* canvas);
	void SaveCanvasROOT(TCanvas* canvas);
	void SaveCanvasPNG(TCanvas* canvas);
	void SaveTwoHistos(std::string canvasName,TH1F* histo1,TH1F* histo2,double refactorSecond=1, UInt_t verbosity=0);
	void SaveHistogram(TH1* histo, bool fitGauss = 0,bool adjustRange =0);
	void SaveHistogramWithFit(TH1F* histo, TF1* fit, UInt_t verbosity=0);
	void SaveHistogramWithCutLine(TH1F *histo,Float_t cutValue);
	void SaveHistogramLogZ(TH2F* histo);
	void SaveHistogram(TH2F* histo,bool drawStatBox=true);
	void SaveGraph(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramPNG(TH1* histo);
	void SaveHistogramPNG(TH2F* histo);
	void SaveGraphPNG(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramFitGaussPNG(TH1* histo);
	void SaveHistogramROOT(TH1* histo);
	void SaveHistogramROOT(TH2F* histo);
	void SaveHistogramROOT(TH3F* histo);
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
	void SaveStringToFile(std::string name,std::string data);
	void SetRange(Float_t min,Float_t max);
	static std::pair<Float_t, Float_t> OptimizeXRange(TH1F* histo);
	static void OptimizeXRange(TH2F* histo);
	static void OptimizeYRange(TH2F* histo);
	static void OptimizeXYRange(TH2F* histo);


	static TH2F* CreateScatterHisto(std::string name,std::vector<Float_t> posY, std::vector<Float_t> posX,UInt_t nBins=512,Float_t minRange =-1.*std::numeric_limits<float>::infinity(), Float_t maxRange=std::numeric_limits<float>::infinity());
	static TGraph CreateDipendencyGraph(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos);
	static TH2F* CreateDipendencyHisto(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos,UInt_t nBins=512);
	static TH1F* CreateDistributionHisto(std::string name, std::vector<Float_t> vec,UInt_t nBins=4096,EnumAxisRange range=maxWidth,Float_t xmin=-1*std::numeric_limits<float>::infinity(),Float_t xmax = std::numeric_limits<float>::infinity());
	static Float_t GetMean(std::vector<Float_t> vec);
	static void SaveCanvasPNG(TCanvas *canvas, std::string location, std::string file_name);
	static void SaveCanvasC(TCanvas *canvas, std::string location, std::string file_name);
	static void SaveCanvasRoot(TCanvas *canvas, std::string location, std::string file_name);
	static TGraphErrors CreateErrorGraph(std::string name,std::vector<Float_t> x, std::vector<Float_t> y, std::vector<Float_t> ex, std::vector<Float_t> ey);
	static void CopyAxisRangesToHisto(TH1F* changingHisto,TH1F* axisInputHisto);
private:
	Float_t xRangeMin,xRangeMax;

    TPaveText *pt;
    TDatime dateandtime;
    std::string plots_path;
    unsigned int runNumber;
    unsigned int nEvents;
    void UpdatePaveText();
    TStyle *currentStyle;
    TStyle *currentStyle2D;
    TSystem *sys;
    TFile* histoFile;
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
