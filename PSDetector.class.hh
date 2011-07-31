/*
 * PSDetector.class.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef PSDETECTOR_CLASS_HH_
#define PSDETECTOR_CLASS_HH_


#include <iostream>
#include "TObject.h"

class PSDetector {//: public TObject {
   public:
      //functions
  //    PSDetector();
      PSDetector(int NChannels = 256);
      virtual ~PSDetector();
      void StoreHit(int Channel, short RawADC, float PedestalMean, float PedestalWidth);
      int GetNChannels() {return nchannels;};
      //void SetNChannels(int NChannels) {nchannels = NChannels;};
      float GetSNR(int index) {return (adc[index]-pedmean[index])/pedwidth[index];};
      float GetPSADC(int index) {return adc[index]-pedmean[index];};
      float GetPedestalMean(int index) {return pedmean[index];};
      float GetPedestalWidth(int index) {return pedwidth[index];};
      short GetADC(int index) {return adc[index];};
      int GetChannel(int index) {return channel[index];};

   //protected:
      unsigned char channel[256];
      short adc[256];
      float pedmean[256];
      float pedwidth[256];
      unsigned char itterator;
      unsigned short nchannels;

   //ClassDef(PSDetector,1);
};


#endif /* PSDETECTOR_CLASS_HH_ */
