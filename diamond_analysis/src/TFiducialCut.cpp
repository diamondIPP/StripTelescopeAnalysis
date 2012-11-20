/*
 * FidCutRegion.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */
#include "TFiducialCut.hh"


TFiducialCut::TFiducialCut(int i,Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh){
	SetAllValuesZero();
	index=i;
	SetXLow(xLow);
	SetXHigh(xHigh);
	SetYLow(yLow);
	SetYHigh(yHigh);
}
TFiducialCut::TFiducialCut(int i) {
	SetAllValuesZero();
	index = i;
}

void TFiducialCut::SetAllValuesZero() {
	active = 0;
	x_low = 0;
	x_high = 0;
	y_low = 0;
	y_high = 0;
}

void TFiducialCut::Print () {
	std::cout << "FidCutRegion #" << index << ":\t XLow:\t" << x_low << "\t XHigh:\t" << x_high << "\t YLow:\t" << y_low << "\t YHigh:\t" << y_high << "\n"<<std::flush;
}
