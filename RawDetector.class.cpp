//Diamond data ranges from 0 to 2^12=4096 while silicon data ranges from 0 to 2^8=256 so use short (16 bits)
#include <iostream>

class RawDetector {
   public:
      //functions
      RawDetector();
      ~RawDetector();
      void SetADC(int channel, unsigned short adc);
      unsigned short GetADC(Int_t channel);
      
   protected:
      unsigned short raw_detector_data[256];
};

RawDetector::RawDetector() {}
RawDetector::~RawDetector() {}

void RawDetector::SetADC(int channel, unsigned short adc) {
   raw_detector_data[channel] = adc;
}

unsigned short RawDetector::GetADC(Int_t channel) {
   if(channel>=0 && channel<=255) return raw_detector_data[channel];
   else {
      std::cout<< "RawDetector::GetADC : channel " << channel << " is out of bounds!"<<std::endl;
      return -1;
   }
}
