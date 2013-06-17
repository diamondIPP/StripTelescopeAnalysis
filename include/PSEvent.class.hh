/*
 * PSEvent.class.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef PSEVENT_CLASS_HH_
#define PSEVENT_CLASS_HH_

#include "PSDetector.class.hh"
//#include "Event_Classes.h" //Data Storage and Processing Events
#include "TMath.h"
#include "TObject.h"
#include "TDetector_Data.hh"
#include "TPed_and_RMS.hh"
class PSEvent {//: public TObject {
   public:
      //functions
      PSEvent();
      virtual ~PSEvent();
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


#endif /* PSEVENT_CLASS_HH_ */
