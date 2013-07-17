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
class TFiducialCut:public TObject {
	int index;
	bool active;
	int x_low, x_high, y_low, y_high;
public:
	TFiducialCut(int i,Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh);
	TFiducialCut(int i = 0);
	bool isInFiducialCut(Float_t xVal, Float_t yVal)const {
	  return ((x_low<xVal)&&(xVal<x_high)&&(y_low<yVal)&&(yVal<y_high));
	}
	void SetAllValuesZero();
	void SetXLow(int xl) {x_low = xl;};
	void SetXHigh(int xh) {x_high = xh;};
	void SetYLow(int yl) {y_low = yl;};
	void SetYHigh(int yh) {y_high = yh;};
	void Print(UInt_t i=0);
	void SetActive(bool i) {active = i;};
	bool GetActive() {return active;};
	int GetXLow() {return x_low;};
	int GetXHigh() {return x_high;};
	int GetYLow() {return y_low;};
	int GetYHigh() {return y_high;};

	ClassDef(TFiducialCut,1);

};
#endif /*TFIDUCIALCUT_HH*/

