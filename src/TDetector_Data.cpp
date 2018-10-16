/*
 * TDetector_Data.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TDetector_Data.hh"
using namespace std;

TDetector_Data::TDetector_Data()
{
	for(Int_t i=0; i<256; i++)
	{
		ADC_values[i]=0;
	}
}

TDetector_Data::~TDetector_Data() {
}



uint16_t TDetector_Data::GetADC_value(Int_t index)
{
	return ADC_values[index];
}

void TDetector_Data::SetADC_value(Int_t index, UShort_t value)
{
	ADC_values[index] = value;
}
