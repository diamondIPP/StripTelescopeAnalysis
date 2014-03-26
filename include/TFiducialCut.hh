//Classes and functions for AutoFidCut()
//is loaded by Clustering.class.cpp
//
//2010-11-22 started (max)
//
#ifndef TFIDUCIALCUT_HH
#define TFIDUCIALCUT_HH
#include <iostream>
#include <iomanip>
#include "TROOT.h"
#include "TObject.h"
#include "TCluster.hh"
#include "TCutG.h"
#include "TCanvas.h"
class TFiducialCut:public TObject {
	int index;
	bool active;
	Float_t x_low, x_high, y_low, y_high;
	TString name;
public:
    TCutG* GetFiducialAreaCut(bool bEmphasis=false);
	TFiducialCut(int i,Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh);
	TFiducialCut(int i = 0);
	bool IsInFiducialCut(Float_t xVal, Float_t yVal)const {
	  return ((x_low<xVal)&&(xVal<x_high)&&(y_low<yVal)&&(yVal<y_high));
	}
	void SetName(TString name){this->name=name;}
	TString GetName(){return name;};
	void SetAllValuesZero();
	void SetXLow(Float_t xl) {x_low = xl;};
	void SetXHigh(Float_t xh) {x_high = xh;};
	void SetYLow(Float_t yl) {y_low = yl;};
	void SetYHigh(Float_t yh) {y_high = yh;};
	void Print(UInt_t i=0);
	void SetActive(bool i) {active = i;};
	bool GetActive() {return active;};
	Float_t GetXLow() {return x_low;};
	Float_t GetXHigh() {return x_high;};
	Float_t GetYLow() {return y_low;};
	Float_t GetYHigh() {return y_high;};
	void SetAllValues(Float_t xLow,Float_t xHigh, Float_t yLow, Float_t yHigh);
	void DrawFiducialCutToCanvas(TCanvas* c1, bool bEmphasis);

	ClassDef(TFiducialCut,1);

};
#endif /*TFIDUCIALCUT_HH*/

