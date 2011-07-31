//Classes and functions for AutoFidCut()
//is loaded by Clustering.class.cpp
//
//2010-11-22 started (max)
//
#ifndef FITCUTREGION_HH
#define FITCUTREGION_HH
#include <iostream>
class FidCutRegion {
	int index;
	bool active;
	int x_low, x_high, y_low, y_high;
public:
	FidCutRegion(int i);
	void SetAllValuesZero();
	void SetValueXLow(int xl) {x_low = xl;};
	void SetValueXHigh(int xh) {x_high = xh;};
	void SetValueYLow(int yl) {y_low = yl;};
	void SetValueYHigh(int yh) {y_high = yh;};
	void GetAllValues();
	void SetActive(bool i) {active = i;};
	bool GetActive() {return active;};
	int GetValueXLow() {return x_low;};
	int GetValueXHigh() {return x_high;};
	int GetValueYLow() {return y_low;};
	int GetValueYHigh() {return y_high;};
};
#endif /*CHANNELSCREEN_H*/

