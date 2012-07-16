/*
 * TFidCutRegions.hh
 *
 *  Created on: Jul 13, 2012
 *      Author: bachmair
 */

#ifndef TFIDCUTREGIONS_HH_
#define TFIDCUTREGIONS_HH_
#include <vector>
#include "TFiducialCut.hh"
#include "TROOT.h"
#include "TH1D.h"
class TFidCutRegions {
public:
  TFidCutRegions(std::vector<std::pair <Float_t, Float_t> > xInt,std::vector<std::pair <Float_t, Float_t> >yInt,UInt_t nDiamonds);
  TFidCutRegions(TH1D* histo, int nDiamonds);
  virtual ~TFidCutRegions();
  TFiducialCut* getFidCut(std::string);
  Float_t getXLow(UInt_t i);
  Float_t getYLow(UInt_t i);
  Float_t getXHigh(UInt_t i);
  Float_t getYHigh(UInt_t i);
  void Print(int intend);
private:
  void createFidCuts();
  std::vector<std::pair<Float_t,Float_t> > xInt,yInt;
  UInt_t nDiamonds;
  UInt_t nFidCuts;
  std::vector<TFiducialCut*> fidCuts;
};

#endif /* TFIDCUTREGIONS_HH_ */
