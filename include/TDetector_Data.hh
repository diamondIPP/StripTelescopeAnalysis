/*
 * TDetector_Data.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TDETECTOR_DATA_HH_
#define TDETECTOR_DATA_HH_
#include "TMath.h"
#include <iostream>

using namespace std;
class TDetector_Data {
public:
	TDetector_Data();
	virtual ~TDetector_Data();
	int GetADC_value(Int_t index);
	void SetADC_value(Int_t index, Int_t value);
	int ADC_values[256];
};

#endif /* TDETECTOR_DATA_HH_ */
