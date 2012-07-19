//Classes and functions for AutoFidCut()
//is loaded by Clustering.class.cpp
//
//2010-11-22 started (max)
//
#ifndef TFIDUCIALCUT_HH
#define TFIDUCIALCUT_HH
#include <iostream>
#include "TROOT.h"
class TFiducialCut {
	int index;
	bool active;
	int x_low, x_high, y_low, y_high;
public:
	TFiducialCut(int i,Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh);
	TFiducialCut(int i);
	bool isInFiducialCut(Float_t xVal, Float_t yVal)const {
	  return xVal<x_high&&xVal>x_low&&yVal>y_low&&yVal<y_high;
	}
	void SetAllValuesZero();
	void SetXLow(int xl) {x_low = xl;};
	void SetXHigh(int xh) {x_high = xh;};
	void SetYLow(int yl) {y_low = yl;};
	void SetYHigh(int yh) {y_high = yh;};
	void Print();
	void SetActive(bool i) {active = i;};
	bool GetActive() {return active;};
	int GetXLow() {return x_low;};
	int GetXHigh() {return x_high;};
	int GetYLow() {return y_low;};
	int GetYHigh() {return y_high;};
};
#endif /*TFIDUCIALCUT_HH*/

