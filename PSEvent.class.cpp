//Class for storing pedestal subtracted adc and noise for each channel

#ifndef __PSEvent__
#define __PSEvent__

#include "PSDetector.class.cpp"
#include "Event_Classes.h" //Data Storage and Processing Events
#include "TMath.h"
#include "TObject.h"

class PSEvent {//: public TObject {
   public:
      //functions
      PSEvent();
      ~PSEvent();
      void SetDetector(Int_t det, TDetector_Data Detector, TPed_and_RMS *Pedestal);
      PSDetector GetDetector(Int_t det);
      bool CMN3Flag, CMN5Flag, ErrorFlag;
      void SetStoreThreshold(float StoreThreshold) {store_threshold = StoreThreshold;};
      void SetEventNumber(int EventNumber) {event_number = EventNumber;};
      int GetEventNumber() {return event_number;};
      
   //protected:
      float store_threshold;
      int event_number;
      PSDetector pedsub_detector_data[9]; // scheme for storing and retrieving detectors: (2*n+i) where n={0,1,2,3,4} is detector and i={0,1} is x or y component
      
      //ClassDef(PSEvent,1);
};

#endif

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
