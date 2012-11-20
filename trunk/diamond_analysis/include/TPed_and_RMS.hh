/*
 * TPed_and_RMS.h
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TPED_AND_RMS_H_
#define TPED_AND_RMS_H_
#include "TMath.h"
#include <iostream>
class TPed_and_RMS {
public:
	TPed_and_RMS();
	virtual ~TPed_and_RMS();
    void SetPedValues(Int_t ped_index, Float_t ped_value);
    void SetRMSValues(Int_t RMS_index, Float_t RMS_value);
    Float_t GetPedValues(Int_t ped_index){return Pedestal_values[ped_index];}
    Float_t GetRMSValues(Int_t RMS_index){return RMS_values[RMS_index];}

 protected:
    Float_t Pedestal_values[256];
    Float_t RMS_values[256];
};

#endif /* TPED_AND_RMS_H_ */
