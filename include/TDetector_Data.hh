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

class TDetector_Data {
public:
	TDetector_Data();
	virtual ~TDetector_Data();
	uint16_t GetADC_value(Int_t index);
	void SetADC_value(Int_t index, UShort_t value);

	uint16_t ADC_values[256];
};

#endif /* TDETECTOR_DATA_HH_ */
