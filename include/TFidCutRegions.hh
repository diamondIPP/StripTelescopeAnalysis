/*
 * TFidCutRegions.hh
 *
 *  Created on: Jul 13, 2012
 *      Author: bachmair
 */

#ifndef TFIDCUTREGIONS_HH_
#define TFIDCUTREGIONS_HH_
#include <vector>
//#include <pair>
#include "TFiducialCut.hh"
#include "TROOT.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "TH2F.h"
#include "TPlaneProperties.hh"
#include "TObject.h"
#include "TMath.h"
#include "TCluster.hh"
#include "TCutG.h"
#include "TLegend.h"
class TFidCutRegions:public TObject {
public:
	TFidCutRegions();
	TFidCutRegions(std::vector< std::pair < Float_t, Float_t> > xInt, std::vector< std::pair < Float_t, Float_t> > yInt, UInt_t nDiamonds);
	TFidCutRegions(Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh,UInt_t nDiamonds);
	TFidCutRegions(TH2F* histo, int nDiamonds,Float_t fidCutPercentage);
	virtual ~TFidCutRegions();
	void SetName(std::string newName){this->name=newName;};
	TCanvas* getFiducialCutCanvas(TPlaneProperties::enumCoordinate cor);
	TFiducialCut* getFidCut(std::string);
	TFiducialCut* getFidCut(UInt_t index);
	void addFiducialCut(Float_t xLow,Float_t xHigh,Float_t yLow, Float_t yHigh);
	void addFiducialCut(TFiducialCut* fidCut);
	Float_t getXLow(UInt_t i=0);
	Float_t getYLow(UInt_t i=0);
	Float_t getXHigh(UInt_t i=0);
	Float_t getYHigh(UInt_t i=0);
	Float_t getHigh(TPlaneProperties::enumCoordinate cor, UInt_t i);
	Float_t getLow(TPlaneProperties::enumCoordinate cor, UInt_t i);

	void Print(int intend = 0);
	void setRunDescription(std::string runDes);
	Int_t getFiducialCutIndex(Float_t xVal, Float_t yVal);
	bool IsInFiducialCut(Float_t xVal,Float_t yVal);
	int getFidCutRegion(Float_t xVal,Float_t yVal);
	TCanvas* getAllFiducialCutsCanvas(TH2F* hScatterPlot=0, bool optimizeAxisRange = false);
	void setHistogramm(TH2F* hEventScatterPlot){this->hEventScatterPlot=hEventScatterPlot;}
	UInt_t getNFidCuts(){return fidCuts.size();}
	Float_t getMaxFiducialX(UInt_t index = 0);
	Float_t getMinFiducialX(UInt_t index = 0);
	Float_t getMaxFiducialY(UInt_t index = 0);
	Float_t getMinFiducialY(UInt_t index = 0);
	void DrawFiducialCutsToCanvas(TCanvas* c1,bool DrawLegend=false);
	TCutG* getFiducialAreaCut(UInt_t nFidCut);
	TCutG* getCutG(TString name, Float_t xLow,Float_t yLow, Float_t xHigh, Float_t yHigh);
	void Reset();
	UInt_t size(){return getNFidCuts();}
	Float_t getAddionalCutXHigh() const {return addionalCut_xHigh;}
	void setAddionalCutXHigh(Float_t addionalCutXHigh) { addionalCut_xHigh = addionalCutXHigh; }
	Float_t getAddionalCutXLow() const { return addionalCut_xLow; }
	void setAddionalCutXLow(Float_t addionalCutXLow) { addionalCut_xLow = addionalCutXLow; }
	Float_t getAddionalCutYHigh() const { return addionalCut_yHigh; }
	void setAddionalCutYHigh(Float_t addionalCutYHigh) { addionalCut_yHigh = addionalCutYHigh; }
	Float_t getAddionalCutYLow() const { return addionalCut_yLow; }
	void setAddionalCutYLow(Float_t addionalCutYLow) { addionalCut_yLow = addionalCutYLow; }
	UInt_t getActiveIndex(){return index;}
private:
	std::string name;
	void initVariables();
	std::vector<std::pair<Float_t,Float_t> >findFiducialCutIntervall(TH1D* hProj,Float_t fidCutPercentage);
	TCanvas*  getFiducialCutProjectionCanvas(TH1D* hProj,std::vector< std::pair<Float_t,Float_t> > intervals);

	TPaveText* getFiducialAreaPaveText(UInt_t nFidCut);
	UInt_t index;
	void createFidCuts();
	std::vector<std::pair<Float_t,Float_t> > xInt,yInt;
	UInt_t nDiamonds;
	UInt_t nFidCuts;
	std::vector<TFiducialCut*> fidCuts;
	std::string runDescription;
	TH2F* hEventScatterPlot;
	UInt_t verbosity;
	Float_t addionalCut_xLow;
	Float_t addionalCut_xHigh;
	Float_t addionalCut_yLow;
	Float_t addionalCut_yHigh;

    ClassDef(TFidCutRegions,1)
};

#endif /* TFIDCUTREGIONS_HH_ */
