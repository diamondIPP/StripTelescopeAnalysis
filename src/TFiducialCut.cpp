/*
 * FidCutRegion.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */
#include "TFiducialCut.hh"

ClassImp(TFiducialCut);
using namespace std;
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
	std::cout << "FidCutRegion #" << index << ":"<<
			" X: " <<setw(6) <<std::right<< x_low << " - " <<setw(6)<<std::left<< x_high <<
			" Y: " <<setw(6) <<std::right<< y_low << " - " <<setw(6)<<std::left<<y_high <<std::right<< "\n"<<std::flush;
}
