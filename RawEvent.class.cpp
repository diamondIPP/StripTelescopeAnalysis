//Class for storing raw data in a more standard format
//2010-07-28 Realize that due to scrambling of the data (example 0x3615abcd written in header as 0x1536cdab) and realizing that both silicon and diamond data is written down as 4-byte words, detectors should be swapped as follows: 0->3, 1->2, 2->1, 0->3, 4->5, 5->4; the last two say that dia0->dia1 and dia1->dia0

#include "RawDetector.class.cpp"

   
class RawEvent {
   public: 
      //functions
      RawEvent() {};
      RawEvent(Int_t run_number, RZEvent rzevent);
      ~RawEvent();
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

RawEvent::RawEvent(Int_t run_number, RZEvent rzevent) {
   
   RunNumber = run_number;
   EvTrig = rzevent.EvTrig;
   EvNo = rzevent.EvNo;
   EvPos = rzevent.EvPos;
   EvTag = rzevent.EvTag;
   EvDate = rzevent.EvDate;
   EvTime = rzevent.EvTime;
   TrigCnt = rzevent.TrigCnt;
   EvVmeTime = rzevent.EvVmeTime;
   EvNetTime = rzevent.EvNetTime;
   MeasNo = rzevent.MeasNo;
   EvInMeasNo = rzevent.EvInMeasNo;
   Eor = rzevent.Eor;

   for(int i=0; i<8; i++) {
      VFasCnt[i] = rzevent.VFasCnt[i];
      VFasReg[i] = rzevent.VFasReg[i];
   }
   
   for(int i=0; i<2; i++) {
      Reserved[i] = rzevent.Reserved[i];
   }
   
   //store silicon data
   for(int detector=0; detector<4; detector++) {
      for(int orientation=0; orientation<2; orientation++) {
         for(int channel=0; channel<256; channel++) {
            //The 3X and 1X detectors are actually physically flipped so we need to reverse the channels
            if((detector==1||detector==3)&&orientation==0) raw_detector_data[2*detector+orientation].SetADC(255-channel,rzevent.Input[4*channel+(3-detector)+1024*orientation]);
            //The other detectors are fine
            else raw_detector_data[2*detector+orientation].SetADC(channel,rzevent.Input[4*channel+(3-detector)+1024*orientation]);
         }
      }
      for(int channel=0; channel<512; channel++) 
         raw_detector_data_chained_unflipped[detector][channel] = rzevent.Input[4*channel+(3-detector)];
   }
   
   //store diamond data
   for(int channel=0; channel<256; channel++) {
      raw_detector_data[2*4+1].SetADC(channel,rzevent.RD42[2*channel]);
      //cout << "dia0\t" << channel << "\t" << rzevent.RD42[2*channel] << "\t" << raw_detector_data[2*4+0].GetADC(channel) <<endl;
   }
   for(int channel=0; channel<256; channel++) {
      raw_detector_data[2*4+0].SetADC(channel,rzevent.RD42[2*channel+1]);
      //cout << "dia1\t" << channel << "\t" << rzevent.RD42[2*channel+1] << "\t" << raw_detector_data[2*4+1].GetADC(channel) <<endl;
   }
}

RawEvent::~RawEvent() {}

RawDetector RawEvent::GetDetector(Int_t det) {
   return raw_detector_data[det];
}

unsigned short* RawEvent::GetUnflippedDetector(Int_t det) {
   return &raw_detector_data_chained_unflipped[det][0];
}
