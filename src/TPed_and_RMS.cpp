/*
 * TPed_and_RMS.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TPed_and_RMS.hh"
TPed_and_RMS::TPed_and_RMS()
{
   for(Int_t i=0; i<256; i++)
   {
      Pedestal_values[i] = 0;
      RMS_values[i] = 1000;
   }
}

TPed_and_RMS::~TPed_and_RMS() {
	// TODO Auto-generated destructor stub
}


void TPed_and_RMS::SetPedValues(Int_t ped_index, Float_t ped_value)
{
   Pedestal_values[ped_index]=ped_value;
}
void TPed_and_RMS::SetRMSValues(Int_t RMS_index, Float_t RMS_value)
{
   RMS_values[RMS_index]=RMS_value;
}
