/*
 * SlidingPedestal.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef SLIDINGPEDESTAL_HH_
#define SLIDINGPEDESTAL_HH_


//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>

//ROOT Class Headers
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TMath.h"
using namespace TMath;
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TImage.h"
//#include "TROOT.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TText.h"

//Taylor
#include "TMultiGraph.h"
#include "TH2I.h"
//AYSHA CHECK
//Class and stuct definitions
#include "TEvent_Array.hh"
#include "TPed_and_RMS.hh"
#include "TDetector_Data.hh"
#include "TTrigger_Event.hh"
//#include "Event_Classes.h" //Data Storage and Processing Events
//#include "PSEvent.class.hh" //Pedestal Subtracted Data Saved in PSEvent Class
#include "ChannelScreen.hh" //Channel Screen Class
//#include "SaveToFile.h" //Functions to save plots as .png, .C or .root files
#include "HistogrammSaver.class.hh"
#include "RZEvent.struct.hh" //the header file that is connected to the Diamond/telescope data
#include "TSettings.class.hh"
#include "TRawEventReader.hh"

class SlidingPedestal {
   public:
      //functions
      SlidingPedestal(unsigned int RunNumber, std::string RunDescription = "");
      ~SlidingPedestal();
      void SetDetector(Int_t det, TDetector_Data Detector, TPed_and_RMS *Pedestal);
      void PedIteration(TH1F *hist, TEvent_Array *Event, TPed_and_RMS *array, Int_t channel, Float_t hit_factor);
      void BufferFill(TDetector_Data &Anyf, TPed_and_RMS *ped, std::vector< std::deque<Int_t> > *buffer_deque, Int_t channel_number, Int_t deque_size, Float_t Hit, Int_t initial_event);
      void RunningPedestal(TDetector_Data detector_buff, TPed_and_RMS initial, TPed_and_RMS *store, std::vector< std::deque<Int_t> > &buffer_deque, Int_t channel_number, Int_t deque_size, Float_t threshold_factor, Int_t *zeroRMS, Int_t event);
     // void RunningCommonMode(std::deque<Double_t> &CMN_deque, Double_t &CMN_Mean, Double_t &CMN_RMS, Double_t new_ave, Int_t deque_size, Int_t event, Int_t *zeroRMS);
      void RunningPedestal_CMN(TDetector_Data detector_buff, TPed_and_RMS initial, TPed_and_RMS *store, std::vector< std::deque<Float_t> > &buffer_deque, Int_t channel_number, Int_t deque_size, Float_t threshold_factor, Int_t *zeroRMS, Int_t event, Float_t correction);
      void Hit_Occupancy(ChannelScreen screen, TH1F *occup, TDetector_Data detector_buffer, TPed_and_RMS *ped_store, Float_t RMS_factor, Int_t chan_begin, Int_t chan_end, Int_t const dia_offset, Float_t correction);
      void PedRMSCalcFromBuffer(std::vector< std::deque<Int_t> > &buffer_deque, TPed_and_RMS *ped_store);
      void Slide(Int_t NEvents, Int_t Initial_Event = 1000, Int_t hit_occupancy = 0);

//      /**** Endian functions for interpreting the rz data ****/
//      //This function swaps the endianess of the read in data for a signed 32-bit integer.
//      void endian_swap(int& x) { x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24); }
//      //This function is overloaded to swap the endianness of the read in data for an unsigned 32-bit integer.
//      void uendian_swap(unsigned int& x) { x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24); }
//      //This function swaps the endianess of the read in data for a signed 16-bit integer.
//      void short_endian_swap(short int& x) { x = (x>>8) | (x<<8); }
//      //This function swaps the endianness of the read in data for an unsigned 16-bit integer.
//      void ushort_endian_swap(unsigned short int& x) { x = (x>>8) | (x<<8); }
//      //Note: Char variables do not need an endian swap routine because endianness is only effected at the byte level, and since the char is only 1 byte, a swaping routine does nothing.
//      /********************Main Routine for Reading in Data, Swapping Endianness and Output ****/
//      int ReadRawEvent(int EventNumber, bool verbose = 0);

private:
      void initialiseDiamondHistogramms(Int_t EventNumber);
      TMultiGraph *DiamondChannelM[128];
      TGraph *DiamondChannelADC[128];
      TGraph *DiamondChannelPedestal[128];
      TGraph *DiamondChannelPedUp[128];
      TGraph *DiamondChannelPedUp2[128];
      TGraph *DiamondChannelPedUp3[128];
      TGraph *DiamondChannelPedUp5[128];
      TGraph *DiamondChannelPedDown[128];
      TGraph *DiamondChannelPedDown2[128];
      TGraph *DiamondChannelPedDown3[128];
      TGraph *DiamondChannelPedDown5[128];
      void initialisetSingleChannel2000plots(Int_t EventNumber);
      TMultiGraph *SingleChannelM[8][256];
      TGraph *SingleChannelADC[8][256];
      TGraph *SingleChannelPedestal[8][256];
      TGraph *SingleChannelPedUp[8][256];
      TGraph *SingleChannelPedUp2[8][256];
      TGraph *SingleChannelPedUp3[8][256];
      TGraph *SingleChannelPedUp5[8][256];
      TGraph *SingleChannelPedDown[8][256];
      TGraph *SingleChannelPedDown2[8][256];
      TGraph *SingleChannelPedDown3[8][256];
      TGraph *SingleChannelPedDown5[8][256];

      void createPlotTag(Int_t NEvents);
      TDatime dateandtime;
      TPaveText *pt;
   protected:

      int run_number, event_number;

      //Code revision declarations (to avoid editing ClusterVar.h needlessly):
   //   float fix_dia_noise; // fix_dia_noise<0 disables diamond noise-fixing
      float store_threshold; // zero suppression to reduce amount of data stored

   private:
      TRawEventReader* rawEventReader;
      TSettings *settings;

      ChannelScreen Det_channel_screen[9];

      //paths
      std::string plots_path;
      std::string png_file_char;
      std::string C_file_char;
      std::string root_file_char;
      TSystem* sys;
      std::string settings_file;
      std::stringstream plotspath;
      std::stringstream settingspath;
      std::stringstream pedfilepath;

      //event storage; written to pedtree
      UInt_t Det_NChannels[9];
      UChar_t Det_Channels[9][256];
      UChar_t Det_ADC[8][256];
      UShort_t Dia_ADC[256];
      Float_t Det_PedMean[9][256];
      Float_t Det_PedWidth[9][256];
//      bool CMNEvent_flag;
      bool ZeroDivisorEvent_flag;

      TH1F *hRMSDifference;
      Int_t numberPlottedBufferNoiseHistos;
      std::vector<TH1F*> hBufferNoise;
      std::vector<int> plottedBufferEvents;

      Int_t nBufferPlots; //counter
      Float_t rms_sigma_difference_cut;

      TH1F *hPullDist[128];

      //end Taylor's stuff
};



#endif /* SLIDINGPEDESTAL_HH_ */
