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
#include "TProfile2D.h"
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
#include "TCutG.h"
#include "TPaveStats.h"
#include "TSettings.class.hh"
#include "TText.h"
#include "TFidCutRegions.hh"
#include "TFiducialCut.hh"

//#include <sys/dirent.h>
#include <sys/stat.h>

class HistogrammSaver {
private:
	unsigned int verbosity;
public:

	enum EnumAxisRange{
		maxWidth,fiveSigma,threeSigma,positiveArea,positiveSigma,manual
	};
	HistogrammSaver(TSettings* settings,int verbosity=0);
	virtual ~HistogrammSaver();
	void InitializeGridReferenceDetSpace();
	void SetOptStat(std::string optStat){gStyle->SetOptStat(optStat.c_str());}
	void SetOptStat(Int_t stat){gStyle->SetOptStat(stat);}
	void SetOptFit(Int_t fitOpt){gStyle->SetOptFit(fitOpt);}
	TPaveText* updateMean(TH1F* histo,Float_t minX = (-1)*std::numeric_limits<float>::infinity(),
								Float_t maxX =      std::numeric_limits<float>::infinity() );
	TPaveText* GetUpdatedLandauMeans(TH1F* histo,Float_t mpv);
	//void SetOptFit(std::string fitOpt){gStyle->SetOptFit(fitOpt.c_str());}
	void SaveCanvas(TCanvas* canvas);
	void SaveCanvasROOT(TCanvas* canvas);
	void SaveCanvasPNG(TCanvas* canvas);
	void SaveTwoHistos(std::string canvasName,TH1F* histo1,TH1F* histo2,double refactorSecond=1, UInt_t verbosity=0);
	void SaveTwoHistosNormalized(std::string canvasName,TH1* histo1,TH1* histo2,double refactorSecond=1, UInt_t verbosity=0);
	void SaveHistogramLandau(TH1F* histo);
	void SaveHistogram(TH2* histo,bool drawStatBox=true);
	void SaveHistogram(TH1* histo, bool fitGauss = 0,bool adjustRange =0,bool drawStatsBox = true);
	void SaveHistogramWithFit(TH1F* histo, TF1* fit,UInt_t verbosity=0);
	void SaveHistogramWithFit(TH1F* histo, TF1* fit, Float_t xmin, Float_t xmax, UInt_t verbosity=0);
	void SaveHistogramWithCutLine(TH1F *histo,Float_t cutValue);
	void SaveHistogramLogZ(TH2* histo);
	void SaveGraph(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramPNG(TH1* histo);
	void SaveHistogramPNG(TH2* histo);
	void SaveGraphPNG(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramFitGaussPNG(TH1* histo);
	void SaveHistogramROOT(TH1* histo);
	void SaveHistogramROOT(TH2* histo);
	void SaveHistogramROOT(TH3F* histo);
	void SaveGraphROOT(TGraph* graph,std::string name,std::string option="AP");
	void SaveHistogramPDF(TH1F* histo);
	void SaveHistogramPDF(TH2* histo);
	void DrawFailedQuarters(vector< pair<Int_t,Int_t> > failedQuarters,TCanvas*c1);
//	void SaveHistogramWithCellGrid(TH2D* histo){SaveHistogramWithCellGrid(TH2*)histo)};
//	void SaveHistogramWithCellGrid(TH2* histo){SaveHistogramWithCellGrid(TH2*)histo)}
	TCanvas* DrawHistogramWithCellGrid(TH2* histo,TH2* histo2=0);
	void SaveHistogramWithCellGrid(TH2* histo, TH2* histo2=0);
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
	static void OptimizeXRange(TH2* histo);
	static void OptimizeYRange(TH2* histo);
	static void OptimizeXYRange(TH2* histo);

	static TProfile2D* CreateProfile2D(std::string name, std::vector<Float_t> posX, std::vector<Float_t> posY, std::vector<Float_t> posZ,
			UInt_t nBinsX=128, UInt_t nBinsY=128,
			Float_t minRangeX = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeX= std::numeric_limits<float>::infinity(),
			Float_t minRangeY = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeY= std::numeric_limits<float>::infinity(),
			Float_t minRangeZ = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeZ= std::numeric_limits<float>::infinity(),
			Float_t factor = 0.05);

	static TH3F* Create3DHisto(std::string name, std::vector<Float_t> posX, std::vector<Float_t> posY, std::vector<Float_t> posZ,
			UInt_t nBinsX=128, UInt_t nBinsY=128,UInt_t nBinsZ=128,
			Float_t minRangeX = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeX= std::numeric_limits<float>::infinity(),
			Float_t minRangeY = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeY= std::numeric_limits<float>::infinity(),
			Float_t minRangeZ = (-1)*std::numeric_limits<float>::infinity(),Float_t maxRangeZ= std::numeric_limits<float>::infinity(),
			Float_t factor = 0.05);
	static TH2F* CreateScatterHisto(std::string name,std::vector<Float_t> posY, std::vector<Float_t> posX,UInt_t nBinsX=512,UInt_t nBinsY=512,
			Float_t minRangeX = (-1)*std::numeric_limits<float>::infinity(), Float_t maxRangeX = std::numeric_limits<float>::infinity(),
			Float_t minRangeY = (-1)*std::numeric_limits<float>::infinity(), Float_t maxRangeY = std::numeric_limits<float>::infinity(),
			Float_t factor = 0.05);
	static TGraph CreateDipendencyGraph(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos, ULong_t maxSize=-1);
	static TH2F* CreateDipendencyHisto(std::string name,std::vector<Float_t> Delta, std::vector<Float_t> pos,UInt_t nBins=512);
	static TH1F* CreateDistributionHisto(std::string name, std::vector<Float_t> vec,UInt_t nBins=4096,EnumAxisRange range=maxWidth,
			Float_t xmin = -1*std::numeric_limits<float>::infinity(),Float_t xmax = std::numeric_limits<float>::infinity(),
			Float_t factor = 0.05);
	static Float_t GetMean(std::vector<Float_t> vec);
	static void SaveCanvasPNG(TCanvas *canvas, std::string location, std::string file_name);
	static void SaveCanvasC(TCanvas *canvas, std::string location, std::string file_name);
	static void SaveCanvasRoot(TCanvas *canvas, std::string location, std::string file_name);
	static TGraphErrors CreateErrorGraph(std::string name,std::vector<Float_t> x, std::vector<Float_t> y, std::vector<Float_t> ex, std::vector<Float_t> ey);
	static void CopyAxisRangesToHisto(TH1F* changingHisto,TH1F* axisInputHisto);
	TH2D* GetHistoBinedInQuarters(TString name);
	TH2D* GetHistoBinedInCells(TString name,Int_t binsPerCellAxis=1);
	TProfile2D* GetProfile2dBinedInCells(TString name,Int_t binsPerCellAxis=1);
	TH3D* Get3dHistoBinedInCells(TString name,UInt_t binsz, Float_t minz,Float_t maxz, Int_t binsPerCellAxis=1);
private:
	Float_t xRangeMin,xRangeMax;
    TPaveText *pt;
    TH2D* hGridReferenceDetSpace;
    TH2D* hGridReferenceCellSpace;
    TDatime dateandtime;
    std::string plots_path;
    unsigned int runNumber;
    unsigned int nEvents;
    void UpdatePaveText();
    TStyle *currentStyle;
    TStyle *currentStyle2D;
    TSystem *sys;
    TFile* histoFile;
    TSettings* settings;
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
