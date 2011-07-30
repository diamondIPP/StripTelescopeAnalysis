/*
 * TEvent_Array.h
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TEVENT_ARRAY_H_
#define TEVENT_ARRAY_H_
#include "TMath.h"
#include <iostream>
using namespace std;
class TEvent_Array {
public:
	TEvent_Array(Int_t size);
	virtual ~TEvent_Array();
    Int_t GetChannelValue(Int_t index);
    Int_t GetSize();
    void SetChannelValue(Int_t value, Int_t index);
    void SetChannelValueSquared(Int_t value, Int_t index);
    Float_t CalculateSigma(Float_t &new_mean);
    Float_t CalculateSigma();

 private:
    Int_t Iteration_Size;
    vector<Int_t> Channel_ADC_values;
    //Int_t Channel_ADC_values[1000];
    //Int_t Channel_ADC_values_squared[1000]; //Taylor
};

#endif /* TEVENT_ARRAY_H_ */
