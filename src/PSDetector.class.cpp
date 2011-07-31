//Class for storing pedestal subtracted data (zero suppressed to exclude non-hits) for a detector plane
//2010-07-11 Static arrays for data (not writing this class to disk anymore)

#include "PSDetector.class.hh"

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
