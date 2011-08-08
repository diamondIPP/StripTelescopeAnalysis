//*********************Diamond and Silicon Detectors Structure Header*****************//
#include <string>

//***Variables and Constants***//
const int DIAMOND_MEM = 256; //263 appears to be the necessary number of 4-byte words that should be
                             //allotted for the diamond data so that the last 4-byte words is the Trailer.
                             //Or, it is 526 2-byte (16-bit) words.

//UPDATE -- 526 is the number that should be used for October 2006 Test Beam Data.
//          256 is the number that should be used for August 2007 Test Beam Data.

//Hardware Channel Offset--the number of empty diamond channels that are run through the analysis (determined by hardware setup)
//For October 2006 Test beam OFFSET = 7, for August 2007, OFFSET = 0.

//const Int_t DIA_OFFSET = 0; // 0 is usual but 7 is needed for 2006 data
const int DIA_OFFSET = 0; // 0 is usual but 7 is needed for 2006 data

const int SILICON_MEM = 2048; //Number of bytes of memory for each 8-bit ADC input for the telescope detectors (4 detectors *  2 layers * 256 channels)


//const int Sil_Chan_Num = 256; //deprecated: Number of channels (and also ADC inputs) in each silicon telescope detector

//const int Dia_Chan_Num = 263; //deprecated: Number of ADC inputs read in for the diamond detector. The number of actual diamond channels is 128 (see note for DIAMOND_MEM)
                              //USe 263 for October 2006 Data
                            
//Int_t Number_events = 0; //deprecated: Initialization of Number of events read in variable

//string di_name = "Dia1"; //deprecated
//string di_name = "Dia0"; //Only used for Oct 2006 Runs

//***Structures***//
struct RZEvent
{
      //The Event structure is the structure that holds the memory sizes for each type of data that is read in from the raw binary files. It determines
      //how many bits are to be allotted for each data value. An entire event includes header data, the actual ADC values for both the silicon detectors
      //and the diamond detectors and a trailer data value to signal the end of the event. The names for the header and trailer data are taken straight 
      //from the original DAQ file.

      //int space[804*1]; 
      //space holding variable. Assuming that the total event size is 3216 bytes, I can use
      //this variable to push down the outputted data into a later event.
      

      //start of actually Event Structure data

      //Event Header Variables
      unsigned int EvTrig;
      unsigned int EvNo;
      unsigned int EvPos;

      int EvTag;
      int EvDate;
      int EvTime;

      unsigned int TrigCnt;
      unsigned int EvVmeTime;

      //Assuming NB_MAX_VFAS = 8
      int VFasCnt[8];
      int VFasReg[8];

      int EvNetTime;
      short int MeasNo;
      short int EvInMeasNo;

      //Assuming EVENT_HEADER_RESERVED_ESZ-2 = 2
      int Reserved[2];


      

      //Event Data Variables
      unsigned char Input[SILICON_MEM];  //The silicon reference telescope detectors read out an 8-bit ADC value.
      unsigned short int RD42[DIAMOND_MEM]; //The diamond detectors read out a 16-bit ADC value although only the first 12 bits are significant.



      //Event Trailer Variable
      unsigned int Eor;

};

//end of structure declaration
