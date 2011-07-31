/*
 * RawEvent.class.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef RAWEVENT_CLASS_HH_
#define RAWEVENT_CLASS_HH_

#include "RawDetector.class.hh"
#include "RZEvent.struct.cpp"

class RawEvent {
   public:
      //functions
      RawEvent() {};
      RawEvent(Int_t run_number, RZEvent rzevent);
      virtual ~RawEvent();
      RawDetector GetDetector(Int_t det);
      unsigned short* GetUnflippedDetector(Int_t det);

   protected:
      //Event Header Variables
      unsigned int EvTrig;
      unsigned int EvNo;
      unsigned int EvPos;
      int EvTag;
      int EvDate;
      int EvTime;
      unsigned int TrigCnt;
      unsigned int EvVmeTime;
      int VFasCnt[8]; //Assuming NB_MAX_VFAS = 8
      int VFasReg[8];
      int EvNetTime;
      short int MeasNo;
      short int EvInMeasNo;
      int Reserved[2]; //Assuming EVENT_HEADER_RESERVED_ESZ-2 = 2
      int Eor; //Event Trailer Variable
      int RunNumber;
      RawDetector raw_detector_data[10]; // scheme for storing and retrieving detectors: (2*n+i) where n={0,1,2,3,4} is detector and i={0,1} is x or y component
      //RawSiliconDetector raw_si_detector_data[8]; // scheme for storing and retrieving detectors: (2*n+i) where n={0,1,2,3} is detector and i={0,1} is x or y component
      //RawDiamondDetector raw_di_detector_data[2]; // scheme for storing and retrieving detectors: (2*n+i) where n={0,1,2,3} is detector and i={0,1} is x or y component
      unsigned short raw_detector_data_chained_unflipped[4][512];
};
#endif /* RAWEVENT_CLASS_HH_ */
