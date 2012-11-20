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
class TFidCutRegions {
public:
  TFidCutRegions(std::vector<std::pair <Float_t, Float_t> > xInt,std::vector<std::pair <Float_t, Float_t> >yInt,UInt_t nDiamonds);
  TFidCutRegions(Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh,UInt_t nDiamonds);
  TFidCutRegions(TH2F* histo, int nDiamonds,Float_t fidCutPercentage);
  virtual ~TFidCutRegions();
  TCanvas* getFiducialCutCanvas(TPlaneProperties::enumCoordinate cor);
  TFiducialCut* getFidCut(std::string);
  Float_t getXLow(UInt_t i);
  Float_t getYLow(UInt_t i);
  Float_t getXHigh(UInt_t i);
  Float_t getYHigh(UInt_t i);
  void Print(int intend);
  void setRunDescription(std::string runDes);
  bool isInFiducialCut(Float_t xVal,Float_t yVal);
  int getFidCutRegion(Float_t xVal,Float_t yVal);
  TCanvas* getAllFiducialCutsCanvas(TH2F* hScatterPlot=0);
  void setHistogramm(TH2F* hEventScatterPlot){this->hEventScatterPlot=hEventScatterPlot;}
private:
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
};

#endif /* TFIDCUTREGIONS_HH_ */