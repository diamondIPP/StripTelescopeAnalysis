 /*******************************Raw Event Routine**********************************/

/*----------------------------------------------------------------------------------------------------------------
As arguements, Eventread_final takes in the Run name in the form "RUN_100XY" and the integer number of data files 
in the entire run.

The Event read in routine is a root program that takes the raw binary data from the October 2006 CERN Test Beam run,
reads it in, and then sorts the bits into various sized variables as described by the Event structure which is listed
in the Diamondstruct.h header file. Each data file contains up to 10000 events where each event has header and trailer
data and ADC values for 4 silicon detectors (the reference telescope) and a diamond detector. An issue with the raw data
is that it was originally written using one type of endianness, but is being read in and analyzed by a computer that the
opposite endianess. Thus, after reading in the bits, the byte order is swapped before being read out as a data value. 
A number of endian swap routines are listed before the main program for various data types. The silicon ADC values are 
only one byte in size and do not need to be endian swapped. Anything larger than one byte however is swapped. The data
values are then sorted and stored into arrays depending on which detector or layer (X or Y) they correspond to. The original
DAQ readme listed the order in which the 2048 bytes were read in. Once sorted, the arrays are stored as branches on a TTree 
and then written to a root file. This root file is then opened by other programs including Eventanalyze which displays
the ADC values on an event by event basis. 

Note on Variable names: The 4 silicon detectors have 2 layers each and are listed as detectors 0 through 4 and layers
X and Y. Thus the first detector, X layer is D0X and so on. The diamond detector is ususally referred to as Dia or Diamond.
--------------------------------------------------------------------------------------------------------------------*/


#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "RZEvent.struct.cpp"
#include "RawEvent.class.cpp"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;



/**** Endian Functions ****/


//This function swaps the endianess of the read in data for a signed 32-bit integer.
void endian_swap(int& x)
{
   x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24);
}

//This function is overloaded to swap the endianness of the read in data for an unsigned
//32-bit integer.
void uendian_swap(unsigned int& x)
{
   x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24);
}

//This function swaps the endianess of the read in data for a signed 16-bit integer.
void short_endian_swap(short int& x)
{
   x = (x>>8) | (x<<8);
}

//This function swaps the endianness of the read in data for an unsigned 16-bit integer.
void ushort_endian_swap(unsigned short int& x)
{
   x = (x>>8) | (x<<8);
}

//Note: Char variables do not need an endian swap routine because endianness is only effected
//at the byte level, and since the char is only 1 byte, a swaping routine does nothing.



/********************Main Routine ************************/
/**** Reading in Data, Swapping Endianness and Output ****/

RawEvent GetRawEvent(int RunNumber, int EventNumber, bool verbose = 0) { 

   //Declarations
   RZEvent TEvent;
   Int_t EventsPerFile = 10000;

   //Check that TEvent size is consistant -- TEvent is declared in Diamondstruct.h


   //Filename to lookup event
   std::ostringstream filename;
   filename << "RUN_" << RunNumber << "_" << EventNumber/EventsPerFile << ".rz";

   //Opening the file that we specified to be read
   ifstream file(filename.str().c_str(),ios::in | ios::binary);  //The .c_str() must be added for ifstream to be able to read in the file name string. 

   if(verbose) {
      cout << "**************************************" << endl;
      cout << "Opening " << filename.str() << " for read" << endl;
      cout << "**************************************" << endl;
   }

   if (!file) {
      cout << "File open error: " << filename.str() << " not found" << endl;
      return RawEvent();
   }

   //Read in event data (Header, Data, and Trailer) using Event structure as specified in Diamondstuct.h
   file.seekg(EventNumber%EventsPerFile * sizeof(TEvent),ios::beg);
   file.read(reinterpret_cast<char*>(&TEvent),sizeof(TEvent));
   file.close(); // don't forget to close the file

   //Changing Endianness of Event Header Read in Data
   uendian_swap(TEvent.EvTrig);
   uendian_swap(TEvent.EvNo);
   uendian_swap(TEvent.EvPos);
   endian_swap(TEvent.EvTag);
   endian_swap(TEvent.EvDate);
   endian_swap(TEvent.EvTime);
   uendian_swap(TEvent.TrigCnt);
   uendian_swap(TEvent.EvVmeTime);
   for (int i=0; i<8; i++)
   {
      endian_swap(TEvent.VFasCnt[i]);
      endian_swap(TEvent.VFasReg[i]);
   }
   endian_swap(TEvent.EvNetTime);
   short_endian_swap(TEvent.MeasNo);
   short_endian_swap(TEvent.EvInMeasNo);
   endian_swap(TEvent.Reserved[0]);
   endian_swap(TEvent.Reserved[1]);

   //Endian Swap for Diamond Data and Outputing the test values for the data following the Telescope Reference Detectors
   for (int i=0; i<DIAMOND_MEM; i++)
      ushort_endian_swap(TEvent.RD42[i]);

   //Swaping Endianness and then Outputing the Event Trailer data
   //uendian_swap(TEvent.Eor);

   //Reading out Event Header Data to Screen
   if(verbose) {
      cout << "Header dump:" << endl;
      cout << "EvTrig: " << TEvent.EvTrig << endl;
      cout << "EvNo: " << TEvent.EvNo << endl;
      cout << "EvPos: " << TEvent.EvPos << endl;
      cout << "EvTag: " << TEvent.EvTag << endl;
      cout << "EvDate: " << TEvent.EvDate << endl;
      cout << "EvTime: " << TEvent.EvTime << endl;
      cout << "TrigCnt: " << TEvent.TrigCnt << endl;
      cout << "EvVmeTime: " << TEvent.EvVmeTime << endl;
      for (int j=0; j<8; j++) cout << "VFasCnt[" << j << "]: " << TEvent.VFasCnt[j] << endl;
      for (int j=0; j<8; j++) cout << "VFasReg[" << j << "]: " << TEvent.VFasReg[j] << endl;
      cout << "EvNetTime: " << TEvent.EvNetTime << endl;
      cout << "MeasNo: " << TEvent.MeasNo << endl;
      cout << "EvInMeasNo: " << TEvent.EvInMeasNo << endl;
      cout << "Reserved[EVENT_HEADER_RESERVED_ESZ-2]: " << TEvent.Reserved[0] << " and " << TEvent.Reserved[1] << endl;
      cout << "Eor: " << TEvent.Eor << endl;
      cout<<endl;
   }
   
   cout<<"TEvent.RD42[55] = "<<TEvent.RD42[55]<<" (should correspond to dia0, channel"<<55/2<<")"<<endl;
   cout<<"TEvent.RD42[161] = "<<TEvent.RD42[161]<<" (should correspond to dia0, channel"<<161/2<<")"<<endl;
   
   /*
   //print out addresses of RZEvent members to make sure that we got it right
   cout<<"&(TEvent) = "<<&(TEvent)<<endl;
   cout<<"&(TEvent.EvTrig) = "<<&(TEvent.EvTrig)<<endl;
   cout<<"&(TEvent.EvNo) = "<<&(TEvent.EvNo)<<endl;
   cout<<"&(TEvent.EvPos) = "<<&(TEvent.EvPos)<<endl;
   cout<<"&(TEvent.EvTag) = "<<&(TEvent.EvTag)<<endl;
   cout<<"&(TEvent.EvDate) = "<<&(TEvent.EvDate)<<endl;
   cout<<"&(TEvent.EvTime) = "<<&(TEvent.EvTime)<<endl;
   cout<<"&(TEvent.TrigCnt) = "<<&(TEvent.TrigCnt)<<endl;
   cout<<"&(TEvent.EvVmeTime) = "<<&(TEvent.EvVmeTime)<<endl;
   cout<<"&(TEvent.VFasCnt[0]) = "<<&(TEvent.VFasCnt[0])<<endl;
   cout<<"&(TEvent.VFasCnt[7]) = "<<&(TEvent.VFasCnt[7])<<endl;
   cout<<"&(TEvent.VFasReg[0]) = "<<&(TEvent.VFasReg[0])<<endl;
   cout<<"&(TEvent.VFasReg[7]) = "<<&(TEvent.VFasReg[7])<<endl;
   cout<<"&(TEvent.EvNetTime) = "<<&(TEvent.EvNetTime)<<endl;
   cout<<"&(TEvent.MeasNo) = "<<&(TEvent.MeasNo)<<endl;
   cout<<"&(TEvent.EvInMeasNo) = "<<&(TEvent.EvInMeasNo)<<endl;
   cout<<"&(TEvent.Reserved[0]) = "<<&(TEvent.Reserved[0])<<endl;
   cout<<"&(TEvent.Reserved[1]) = "<<&(TEvent.Reserved[1])<<endl;
   cout<<"&(TEvent.Input[0]) = "<<(uint*)&(TEvent.Input[0])<<endl;
   cout<<"&(TEvent.Input["<<SILICON_MEM-1<<"]) = "<<(uint*)&(TEvent.Input[SILICON_MEM-1])<<endl;
   cout<<"&(TEvent.RD42[0]) = "<<&(TEvent.RD42[0])<<endl;
   cout<<"&(TEvent.RD42["<<DIAMOND_MEM-1<<"]) = "<<&(TEvent.RD42[DIAMOND_MEM-1])<<endl;
   cout<<"&(TEvent.Eor) = "<<(void*)&(TEvent.Eor)<<endl;
   //cout<<"&(TEvent.Eor)-&(TEvent) = "<<(uint)&(TEvent.Eor)-(uint)&(TEvent)<<endl;
   //cout<<"&(TEvent.Eor)-&(TEvent.EvTrig) = "<<(void*)&(TEvent.Eor)-(void*)&(TEvent.EvTrig)<<endl;
   */
   
   //store the data in a RawEvent class object
   RawEvent rawevent(RunNumber,TEvent);
   
   return rawevent;
   
} //end of binary data read-in while loop



