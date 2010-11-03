//Class for storing pedestal subtracted data (zero suppressed to exclude non-hits) for a detector plane
//2010-07-11 Static arrays for data (not writing this class to disk anymore)

#ifndef __PSDetector__
#define __PSDetector__

#include <iostream>
#include "TObject.h"

class PSDetector {//: public TObject {
   public:
      //functions
  //    PSDetector();
      PSDetector(int NChannels = 256);
      ~PSDetector();
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

#endif

PSDetector::PSDetector(int NChannels) {
   itterator = 0;
   nchannels = NChannels;
   //channel = new unsigned char[NChannels];
   //adc = new short[NChannels];
   //pedmean = new float[NChannels];
   //pedwidth = new float[NChannels];
   //std::cout << "nchannels = " << (int) nchannels << "\tNChannels = " << NChannels << std::endl;
}
PSDetector::~PSDetector() {}

void PSDetector::StoreHit(int Channel, short RawADC, float PedestalMean, float PedestalWidth) {
   if(itterator>=nchannels) {
      std::cout << "PSDetector::StoreHit: Error! Tried to store more channels than allocated memory!\nnchannels = " << int(nchannels) << "\titterator = " << int(itterator) << std::endl;
      return;
   }
   channel[itterator] = (unsigned char) Channel;
   adc[itterator] = (RawADC);
   pedmean[itterator] = (PedestalMean);
   pedwidth[itterator] = (PedestalWidth);
   itterator++;
}

//ClassImp(PSDetector);


//#ifdef __CINT__
//#pragma link C++ class PSDetector+;
//#endif
