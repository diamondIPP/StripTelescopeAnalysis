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
using namespace std;

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
//#include "PSEvent.class.cpp" //Pedestal Subtracted Data Saved in PSEvent Class
#include "ChannelScreen.hh" //Channel Screen Class
//#include "SaveToFile.h" //Functions to save plots as .png, .C or .root files
#include "HistogrammSaver.class.hh"
#include "RZEvent.struct.cpp" //the header file that is connected to the Diamond/telescope data


class SlidingPedestal {
   public:
      //functions
      SlidingPedestal(unsigned int RunNumber, string RunDescription = "");
      ~SlidingPedestal();
      void LoadSettings();
      void ParseIntArray(string value, vector<int> &vec);
      void SetDetector(Int_t det, TDetector_Data Detector, TPed_and_RMS *Pedestal);
      void PedIteration(TH1F *hist, TEvent_Array *Event, TPed_and_RMS *array, Int_t channel, Float_t hit_factor);
      void BufferFill(TDetector_Data &Anyf, TPed_and_RMS *ped, vector< deque<Int_t> > *buffer_deque, Int_t channel_number, Int_t deque_size, Float_t Hit, Int_t initial_event);
      void RunningPedestal(TDetector_Data detector_buff, TPed_and_RMS initial, TPed_and_RMS *store, vector< deque<Int_t> > &buffer_deque, Int_t channel_number, Int_t deque_size, Float_t threshold_factor, Int_t *zeroRMS, Int_t event);
      void RunningCommonMode(deque<Double_t> &CMN_deque, Double_t &CMN_Mean, Double_t &CMN_RMS, Double_t new_ave, Int_t deque_size, Int_t event, Int_t *zeroRMS);
      void Hit_Occupancy(ChannelScreen screen, TH1F *occup, TDetector_Data detector_buffer, TPed_and_RMS *ped_store, Float_t RMS_factor, Int_t chan_begin, Int_t chan_end, Int_t const dia_offset);
      void PedRMSCalcFromBuffer(vector< deque<Int_t> > &buffer_deque, TPed_and_RMS *ped_store);
      void Slide(Int_t NEvents, Int_t Initial_Event = 1000, Int_t hit_occupancy = 0);

      /**** Endian functions for interpreting the rz data ****/
      //This function swaps the endianess of the read in data for a signed 32-bit integer.
      void endian_swap(int& x) { x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24); }
      //This function is overloaded to swap the endianness of the read in data for an unsigned 32-bit integer.
      void uendian_swap(unsigned int& x) { x = (x >> 24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24); }
      //This function swaps the endianess of the read in data for a signed 16-bit integer.
      void short_endian_swap(short int& x) { x = (x>>8) | (x<<8); }
      //This function swaps the endianness of the read in data for an unsigned 16-bit integer.
      void ushort_endian_swap(unsigned short int& x) { x = (x>>8) | (x<<8); }
      //Note: Char variables do not need an endian swap routine because endianness is only effected at the byte level, and since the char is only 1 byte, a swaping routine does nothing.
      /********************Main Routine for Reading in Data, Swapping Endianness and Output ****/
      int ReadRawEvent(int EventNumber, bool verbose = 0);


   protected:
      string current_rz_filename;
      ifstream current_rz_file;
      RZEvent TEvent;
      TDetector_Data D0X, D0Y, D1X, D1Y, D2X, D2Y, D3X, D3Y, Dia0, Dia1;
      //PSEvent PedSubEvent;
      int run_number, event_number;

      //Code revision declarations (to avoid editing ClusterVar.h needlessly):
      float fix_dia_noise; // fix_dia_noise<0 disables diamond noise-fixing
      float store_threshold; // zero suppression to reduce amount of data stored

   private:
      //settings
      Int_t SaveAllFilesSwitch; //1 for save files, 0 for don't
      Int_t ClosePlotsOnSave;
      Int_t IndexProduceSwitch;
      bool single_channel_analysis_enable; // enable channel noise analysis
      Int_t single_channel_analysis_eventwindow; // Number of events to put in each histogram
      Int_t dia_input; // 0 for 2006 and 1 for the rest
      Float_t Si_Pedestal_Hit_Factor;
      Float_t Di_Pedestal_Hit_Factor;
      /***Common Mode Noise Constraints***/
      //Range of Corrected Events from
      Int_t CMN_corr_low;
      //to
      Int_t CMN_corr_high;
      // i.e. any event in between 4sigma and 7sigma will be corrected
      //CMN Cut Factor
      Int_t CMN_cut;  //Should be less than or equal to CMN_coor_high
      Int_t Iter_Size; //buffer size
      Int_t Taylor_speed_throttle; //# of events to recalculate RMS the old way; set to 1 to disable

      //Channels to Screen
      vector<int> single_channel_analysis_channels;

      //Channels to Screen
      vector<int> Det_channel_screen_channels[9];
      vector<int> Det_channel_screen_regions[9];
      ChannelScreen Det_channel_screen[9];

      //paths
      string plots_path;
      string png_file_char;
      string C_file_char;
      string root_file_char;
      TSystem* sys;
      string settings_file;
      ostringstream plotspath;
      ostringstream settingspath;
      ostringstream pedfilepath;

      //event storage; written to pedtree
      UInt_t Det_NChannels[9];
      UChar_t Det_Channels[9][256];
      UChar_t Det_ADC[8][256];
      UShort_t Dia_ADC[256];
      Float_t Det_PedMean[9][256];
      Float_t Det_PedWidth[9][256];
      bool CMNEvent_flag;
      bool ZeroDivisorEvent_flag;

      //Taylor's stuff
      Int_t plotChannel_on; //make RMS Difference plot for all detectors, and Buffer Noise plots for D0X
      Int_t plotDiamond; //make Buffer Noise plots for the diamond instead
      Int_t makeBufferPlots; //make Buffer Plot whenever sigma and rms differ by rms_sigma_difference_cut
      //NOTE: only works if plotChannel_on = 1 and plottedChannel < 256
      Int_t SingleChannel2000plots; //make SC_Pedestal plots for all silicon detectors and channels
      Int_t makeDiamondPlots; //make DC_Pedestal plots for all diamond channels
      Int_t makeHits2D; //make 2D histogram of hits and seeds
      Int_t makeNoise2D; //make 2D histogram of noise per channel
      Int_t makePullDist; //make pull distribution
      Int_t makePedRMSTree; //make .root file of pedestal and rms values
      Int_t eventPrintHex; //print hex (should match .rz data)

      UInt_t plottedChannel; //256 = enter channel on run. also, set to 256 and type 256 to turn off buffer noise plots
      TH1F *hRMSDifference;
      //const Int_t numberPlottedBufferNoiseHistos;
      Int_t numberPlottedBufferNoiseHistos;
      vector<TH1F*> hBufferNoise;
      vector<int> plottedBufferEvents;
      //Taylor, change #6 7/6

      Int_t maxBufferPlots;
      Int_t nBufferPlots; //counter
      Float_t rms_sigma_difference_cut;
      Int_t high_rms_cut; //cut on absolute rms value instead of comparing to Gaussian
      Float_t rms_cut; //value to use if high_rms_cut
      //deque<TH1F*> bufferPlotsDeque;

      Int_t zoomDiamondPlots; //zoom in on DC_Pedestal (100 event / window)

      Int_t singleTrack2D; //plot single tracks only in 2D hits histogram
      Int_t singleTrack2DmaxClusterSize; //max size of clusters in silicon track (cluster = Di_Hit_Factor hits; no check for seeds/shoulders)

      Float_t maxNoise2D; //highest noise value plotted in 2D noise histogram

      TH1F *hPullDist[128];
      //end Taylor's stuff
};



#endif /* SLIDINGPEDESTAL_HH_ */
