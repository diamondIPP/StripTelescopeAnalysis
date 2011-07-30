/*
 * FidCutRegion.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */
#include "FidCutRegion.hh"

FidCutRegion::FidCutRegion(int i) {
	index = i;
}

void FidCutRegion::SetAllValuesZero() {
	active = 0;
	x_low = 0;
	x_high = 0;
	y_low = 0;
	y_high = 0;
}

void FidCutRegion::GetAllValues () {
	std::cout << "FidCutRegion #:\t" << index << "\t XLow:\t" << x_low << "\t XHigh:\t" << x_high << "\t YLow:\t" << y_low << "\t YHigh:\t" << y_high << "\n"<<std::flush;
}
