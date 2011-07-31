//Class for storing pedestal subtracted adc and noise for each channel

#include "PSEvent.class.hh"

PSEvent::PSEvent() {
   //clear flags
   CMN3Flag = false;
   CMN5Flag = false;
   ErrorFlag = false;
   store_threshold = 0;
}

PSEvent::~PSEvent() {
}

PSDetector PSEvent::GetDetector(Int_t det) {
   if(det>=0 && det<9) return pedsub_detector_data[det];
   else {
      std::cout<< "PSEvent::SetDetector : Detector " << det << " is out of bounds. There are only 9 detectors."<<std::endl;
      return PSDetector();
   }
}

void PSEvent::SetDetector(Int_t det, TDetector_Data Detector, TPed_and_RMS *Pedestal) {
   int numberofchannelstosave;
   if(det>=0 && det<9) {
      //figure out how many channels to save
      numberofchannelstosave = 0;
      for(Int_t i=0; i<256; i++)
         if(TMath::Abs(Detector.ADC_values[i]-Pedestal->GetPedValues(i))>store_threshold*Pedestal->GetRMSValues(i))
            numberofchannelstosave++;
//      std::cout << "number of channels = " << numberofchannelstosave <<std::endl;
//      //create new PSDetector with enough room to save channels
//      pedsub_detector_data[det] = new PSDetector(numberofchannelstosave);
      //reset detector channel counters for storing new event
      pedsub_detector_data[det].nchannels = numberofchannelstosave;
      pedsub_detector_data[det].itterator = 0;
      //save channels
      for(Int_t i=0; i<256; i++)
         if(TMath::Abs(Detector.ADC_values[i]-Pedestal->GetPedValues(i))>store_threshold*Pedestal->GetRMSValues(i))
            pedsub_detector_data[det].StoreHit(i, Detector.ADC_values[i], Pedestal->GetPedValues(i), Pedestal->GetRMSValues(i));
   }
   else std::cout<< "PSEvent::SetDetector : Detector " << det << " is out of bounds. There are only 9 detectors."<<std::endl;
}

//ClassImp(PSEvent);


//#ifdef __CINT__
//#pragma link C++ class PSEvent+;
//#endif
