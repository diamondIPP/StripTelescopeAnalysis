//organize pedestal analysis code in a class structure (SlidingPedestal.class.cpp) and executable function (PedestalAnalyze.cpp)
//2010-07-02 Based on Duris' mod of Schaffner's master-2010-06-03/SlidingPedestalV4.cpp 
//           Organized original code in a class structure
//           Replaced TTree pedtree lookup (telescope.root) with GetRawEvent.cpp code (take from .rz files) 
//           Future plan: modular design...keep most data as class members, create in SlidingPedestal() and destroy in ~SlidingPedestal()
//2010-07-06 Based on dev-2010-07-02/SlidingPedestal.class.cpp 
//           Rewriting output tree
//2010-07-08 Finished writing new class structure for output events using zero suppression, but each event is really big (10's to 100's of kb/event)!
//2010-07-11 Now organizing data in static sized PSEvent, but writing variable sized arrays to disk (significant reduction in memory usage)
//2010-07-20 Reading settings from flat file now
//2010-07-21 Compiles now but runs at 1/10 speed with command: g++ -o PedestalAnalyze.exe PedestalAnalyze.cpp `root-config --cflags --glibs`
//2010-07-22 Renamed TChannel_Screen to ChannelScreen and Ped_Subtracted.root to PedestalSubtracted.root
//2010-07-27 Renamed settings file from Setting.ini to Setting.XXXXX.ini, where XXXXX=run number
//2010-07-28 Added Taylor's speed tweak to SlidingPedestal.class.cpp
//           Added Taylor's variable buffer size to SlidingPedestal.class.cpp
//2010-08-14 Removed CMN and "Bad" (aka zero divisor) event cuts; now flag, count, and save those events
//           NOTE: Realize that due to scrambling of the data (example 0x3615abcd written in header as 0x1536cdab) and realizing that both silicon and diamond data is written down as 4-byte words, detectors should be swapped as follows: 0->3, 1->2, 2->1, 0->3, 4->5, 5->4; the last two say that dia0->dia1 and dia1->dia0
//           Reversed silicon plane order and exchanged diamond inputs
//2010-08-30 Flagged CMN and ZeroDivisor events for exclusion from Clustering analysis

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

//Class and stuct definitions
//#include "Event_Classes.h" //Data Storage and Processing Events
//#include "PSEvent.class.cpp" //Pedestal Subtracted Data Saved in PSEvent Class
#include "ChannelScreen.hh" //Channel Screen Class
#include "HistogrammSaver.class.hh" //Functions to save plots as .png, .C or .root files
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
      void Slide(Int_t Event_Number, Int_t Initial_Event = 1000, Int_t hit_occupancy = 0);
      
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
      void ReadRawEvent(int EventNumber, bool verbose = 0);
      
      
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
      Float_t Si_Hit_Factor;
      Float_t Di_Hit_Factor;
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
};


SlidingPedestal::SlidingPedestal(unsigned int RunNumber, string RunDescription) {
   //Code revision declarations (to avoid editing ClusterVar.h needlessly):
   run_number = RunNumber;
   
   //default settings
   SaveAllFilesSwitch = 1; //1 for save files, 0 for don't
   ClosePlotsOnSave = 1;
   IndexProduceSwitch = 1;
   fix_dia_noise = -1;//7.7; // fix_dia_noise<0 disables diamond noise-fixing
   store_threshold = 2;
   //PedSubEvent.SetStoreThreshold(store_threshold);
   single_channel_analysis_enable = false; // enable channel noise analysis
   single_channel_analysis_eventwindow = 5000; // Number of events to put in each histogram 
   dia_input = 0; // 1 for 2006 and 0 for the rest
   Si_Hit_Factor = 5;
   Di_Hit_Factor = 5;
   /***Common Mode Noise Constraints***/
   //Range of Corrected Events from
   CMN_corr_low = 3;
   //to
   CMN_corr_high = 7;
   // i.e. any event in between 4sigma and 7sigma will be corrected
   //CMN Cut Factor
   CMN_cut = 4;  //Should be less than or equal to CMN_coor_high
   Iter_Size = 500; //buffer size
   Taylor_speed_throttle = 1000; //# of events to recalculate RMS the old way; set to 1 to disable


   //default paths
   sys = gSystem;
   
   plotspath << sys->pwd() << "/plots-" << RunNumber;
   if(RunDescription=="") plotspath << "/";
   else plotspath << "-" << RunDescription << "/";
   plots_path = plotspath.str();
   png_file_char = plotspath.str();
   C_file_char = plotspath.str();
   root_file_char = plotspath.str();
   //make plots dir
   if(sys->mkdir(plots_path.c_str())==0) cout<<"Made directory "<<plots_path<<endl;
   
   settingspath << sys->pwd() << "/Settings." << RunNumber;
   if(RunDescription=="") settingspath << ".ini";
   else settingspath << "-" << RunDescription << ".ini";
   settings_file = settingspath.str();
   
   pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
   if(RunDescription=="") pedfilepath << ".root";
   else pedfilepath << "-" << RunDescription << ".root";
   
   
   LoadSettings();
   
   //screen channels
   for(int det=0; det<9; det++) {
      Det_channel_screen[det].ScreenChannels(Det_channel_screen_channels[det]);
      Det_channel_screen[det].ScreenRegions(Det_channel_screen_regions[det]);
      cout<<"Detector "<<det<<" screened channels: ";
      Det_channel_screen[det].PrintScreenedChannels();
      cout<<endl;
   }
}


SlidingPedestal::~SlidingPedestal() {
   current_rz_file.close(); // don't forget to close the file when deleting the class
}

void SlidingPedestal::LoadSettings() {
   
   cout<<endl<<"Overriding default settings with settings in Settings.ini"<<endl<<endl;
   
   ifstream file(settings_file.c_str());
   if(!file) {
      cout << "An error has encountered while trying to open file " << settings_file << endl;
      cout << "Keeping default settings; no channels will be screened." << endl;
      return;
   }
   else cout << settings_file << " successfully opened." << endl << endl;
   

   while(!file.eof()) {
      
      //get next line
      string line;
      getline(file,line);
      
      //check if comment or empty line
      if ((line.substr(0, 1) == ";") || (line.substr(0, 1) == "#") || (line.substr(0, 1) == "/") || line.empty()) {
         continue;
      }
      
      //find the index of first '=' character on the line
      string::size_type offsetl = line.find_first_of('=');
      string::size_type offsetr = line.find_first_of('=');

      //extract the key (LHS of the ini line)
      string key = line.substr(0, offsetl);
      
      //trim spaces from key
      while(line.at(offsetl-1)==' ') {
         offsetl--;
      }
      key = line.substr(0, offsetl);
      
      //extract the value (RHS of the ini line)
      string value = line.substr(offsetr+1, line.length()-(offsetr+1));
      
      //trim spaces from value
      while(line.at(offsetr+1)==' ') {
         offsetr++;
      }
      value = line.substr(offsetr+1, line.length()-(offsetr+1));
      
      //trim end ';' from end of key if found
      if(value.find_first_of(';')!=string::npos) {
         value = line.substr(offsetr+1, value.find_first_of(';'));//line.length()-(offsetr+1)-1);
      }
      
      //cant switch on strings so use if statements
      if(key=="SaveAllFilesSwitch") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         SaveAllFilesSwitch = (int)strtod(value.c_str(),0);
      }
      if(key=="ClosePlotsOnSave") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ClosePlotsOnSave = (int)strtod(value.c_str(),0);
      }
      if(key=="IndexProduceSwitch") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         IndexProduceSwitch = (int)strtod(value.c_str(),0);
      }
      if(key=="single_channel_analysis_enable") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         single_channel_analysis_enable = (int)strtod(value.c_str(),0);
      }
      if(key=="single_channel_analysis_eventwindow") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         single_channel_analysis_eventwindow = (int)strtod(value.c_str(),0);
      }
      if(key=="dia_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         dia_input = (int)strtod(value.c_str(),0);
      }
      if(key=="fix_dia_noise") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         fix_dia_noise = (int)strtod(value.c_str(),0);
      }
      if(key=="store_threshold") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         store_threshold = (float)strtod(value.c_str(),0);
      }
      if(key=="Si_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Si_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Di_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Di_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="CMN_corr_low") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         CMN_corr_low = (int)strtod(value.c_str(),0);
      }
      if(key=="CMN_corr_high") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         CMN_corr_high = (int)strtod(value.c_str(),0);
      }
      if(key=="CMN_cut") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         CMN_cut = (int)strtod(value.c_str(),0);
      }
      if(key=="Iter_Size") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Iter_Size = (int)strtod(value.c_str(),0);
      }
      if(key=="Taylor_speed_throttle") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Taylor_speed_throttle = (int)strtod(value.c_str(),0);
      }
      if(key=="single_channel_analysis_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,single_channel_analysis_channels);
      }
      if(key=="D0X_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[0]);
      }
      if(key=="D0Y_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[1]);
      }
      if(key=="D1X_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[2]);
      }
      if(key=="D1Y_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[3]);
      }
      if(key=="D2X_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[4]);
      }
      if(key=="D2Y_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[5]);
      }
      if(key=="D3X_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[6]);
      }
      if(key=="D3Y_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[7]);
      }
      if(key=="Dia_channel_screen_channels") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_channels[8]);
      }
      if(key=="D0X_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[0]);
      }
      if(key=="D0Y_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[1]);
      }
      if(key=="D1X_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[2]);
      }
      if(key=="D1Y_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[3]);
      }
      if(key=="D2X_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[4]);
      }
      if(key=="D2Y_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[5]);
      }
      if(key=="D3X_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[6]);
      }
      if(key=="D3Y_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[7]);
      }
      if(key=="Dia_channel_screen_regions") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_regions[8]);
      }
         
      
   }
   
   file.close();
   cout<<endl<<"Finished importing settings from Settings.ini"<<endl<<endl;
}

void SlidingPedestal::ParseIntArray(string value, vector<int> &vec) {
   int index=0;
   string::size_type offset1 = value.find_first_of('{')+1;
   string::size_type offset2 = value.find_first_of(',');
   
   vec.push_back((int)strtod(value.substr(offset1,offset2-offset1).c_str(),0));
   
   //cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(offset1,offset2-offset1)="<<value.substr(offset1,offset2-offset1)<<"\tstrtod(value.substr(offset1,offset2-offset1),0)="<<strtod(value.substr(offset1,offset2-offset1).c_str(),0)<<endl;//check
   
   value = value.substr(offset2+1,value.length()-(offset2+1));
   
   while(value.length()>2) {
      offset2 = TMath::Min(value.find_first_of(','),value.find_first_of('}'));
      vec.push_back((int)strtod(value.substr(0,offset2).c_str(),0));
      index++;
      //cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(0,offset2)="<<value.substr(0,offset2)<<"\tstrtod(value.substr(0,offset2),0)="<<strtod(value.substr(0,offset2).c_str(),0)<<endl;//check
      if(value.find_first_of(';')<2) break;
      value = value.substr(offset2+1,value.length()-(offset2+1));
   }
}

void SlidingPedestal::SetDetector(Int_t det, TDetector_Data Detector, TPed_and_RMS *Pedestal) {
   int numberofchannelstosave, itterator;
   if(det<8) {
      //figure out how many channels to save
      numberofchannelstosave = 0;
      for(Int_t i=0; i<256; i++)
         if(TMath::Abs(Detector.ADC_values[i]-Pedestal->GetPedValues(i))>store_threshold*Pedestal->GetRMSValues(i))
            numberofchannelstosave++;
//      std::cout << "number of channels = " << numberofchannelstosave <<std::endl;

      //reset detector channel counters for storing new event
      Det_NChannels[det] = numberofchannelstosave;
      itterator = 0;
      //save channels
      for(Int_t i=0; i<256; i++)
         if(TMath::Abs(Detector.ADC_values[i]-Pedestal->GetPedValues(i))>store_threshold*Pedestal->GetRMSValues(i)) {
            Det_Channels[det][itterator] = (unsigned char) i; //channel number
            Det_ADC[det][itterator] = (unsigned char) Detector.ADC_values[i];
            Det_PedMean[det][itterator] = Pedestal->GetPedValues(i);
            Det_PedWidth[det][itterator] = Pedestal->GetRMSValues(i);
            itterator++;
         }
      
   }
   if(det==8) {
      //how many channels to save (always save all of diamond data)
      numberofchannelstosave = 128;

      //reset detector channel counters for storing new event
      Det_NChannels[det] = numberofchannelstosave;
      //save channels
      for(Int_t i=0; i<128; i++) {
         Det_Channels[det][i] = (unsigned char) i; //channel number
         Dia_ADC[i] = (unsigned short) Detector.ADC_values[i];
         Det_PedMean[det][i] = Pedestal->GetPedValues(i);
         Det_PedWidth[det][i] = Pedestal->GetRMSValues(i);
      }
      
   }
}


void SlidingPedestal::PedIteration(TH1F *hist, TEvent_Array *Event, TPed_and_RMS *array, Int_t channel, Float_t hit_factor) {
   
   Int_t Size = Event->GetSize();
   for(Int_t i=0; i<Size; i++)
   {
      if(Abs((Event->GetChannelValue(i))-(array->GetPedValues(channel)))<(hit_factor*(array->GetRMSValues(channel))))
      {
         hist->Fill(Event->GetChannelValue(i));
      }
   }
   array->SetPedValues(channel,hist->GetMean());
   array->SetRMSValues(channel,hist->GetRMS());
}


void SlidingPedestal::BufferFill(TDetector_Data &Anyf, TPed_and_RMS *ped, vector< deque<Int_t> > *buffer_deque, Int_t channel_number, Int_t deque_size, Float_t Hit, Int_t initial_event) {
   Int_t buffer_marker = 0;
   Int_t event_counter = initial_event; //This should start at the same event (or about) as the initialization pedestal as the difference is what led to the Channel 57 problem
   Float_t deque_sizef = (Float_t) deque_size;
   Float_t Hit_Factor = Hit;
   vector< deque<Int_t> > buffer(256);
   
   Int_t size_scale[channel_number];
   for(Int_t i=0; i<channel_number; i++) {size_scale[i]=0;}

   //while(buffer_marker != channel_number && event_counter < (initial_event+500))
   while (buffer_marker != channel_number && event_counter < initial_event + Iter_Size + 400)
   {
      buffer_marker=0;
      ReadRawEvent(event_counter);
      TTrigger_Event *Event = new TTrigger_Event;
      
      Event->SetAny(Anyf);

      for(Int_t s=0; s<256; s++)
      {
         TDetector_Data *Data = new TDetector_Data;
         *Data = Event->GetAny();
         if(Abs(Data->GetADC_value(s)-ped->GetPedValues(s))>(Hit_Factor*ped->GetRMSValues(s)))
         {
            continue;
         }
         else
         {
            if(buffer[s].size() != deque_sizef)
            {
               if(s==64)
               {
                  //cout << "Buffer Input: " << Data->GetADC_value(s) << endl;
               }
               buffer[s].push_back(Data->GetADC_value(s));
               size_scale[s] = buffer[s].size();
            }
         }
         delete Data;
         Data = 0;
      }


      for(Int_t i=0; i<256; i++)
      {
         if(size_scale[i]==deque_size)
         {
            buffer_marker++;
            //cout << "Buffer Marker " << buffer_marker << endl; 
         }
         else
         {
            //break;
         }
      }
      delete Event;
      Event = 0;
      event_counter++;
      //cout << "Event: " << event_counter << endl;
      if(buffer_marker==256)
      {
         cout << "One detector all done at Event: " << event_counter << endl;
      }
      //if(event_counter==(initial_event+499))
      if (event_counter == initial_event + Iter_Size + 399)
      {
         cout << "One detector overflow" << " buffer marker: " << buffer_marker << endl;
         for(Int_t s=0; s<256; s++)
         {
            while(buffer[s].size()!=deque_sizef)
            {
               Int_t last = buffer[s].size();
               buffer[s].push_back((last!=0) ? (buffer[s][last-1]) : 0);
               size_scale[s] = buffer[s].size();
               //cout << "buffer " << s << " growing to size " << buffer[s].size() << endl;
            }
         }
      }

   }
   *buffer_deque = buffer;
}


void SlidingPedestal::RunningPedestal(TDetector_Data detector_buff, TPed_and_RMS initial, TPed_and_RMS *store, vector< deque<Int_t> > &buffer_deque, Int_t channel_number, Int_t deque_size, Float_t threshold_factor, Int_t *zeroRMS, Int_t event)
{
   Float_t divisor = (Float_t) deque_size;
   TPed_and_RMS Ped_buffer;
   //TEvent_Array RMS_buffer(deque_size);
   Ped_buffer = *store;
   /* for(Int_t s=0; s<channel_number; s++)
   {
      if(Abs(detector_buff.GetADC_value(s)-store->GetPedValues(s))>(threshold_factor*store->GetRMSValues(s)))
      {
         continue;
      }
      Float_t Remove = buffer_deque[s][0]/divisor;
      Float_t Add = (detector_buff.GetADC_value(s)/divisor);
      if(s==64)  //(57+7))
      {
         //cout << "s is: " << s << " ADC Value: " << detector_buff.GetADC_value(s) << " Ped Before: " << Ped_buffer.GetPedValues(s) << " and after: " << Ped_buffer.GetPedValues(s)+Add-Remove << " Add: " << Add << " Remove: " << Remove << endl;
      }
      store->SetPedValues(s,(Ped_buffer.GetPedValues(s))+Add-Remove);
      buffer_deque[s].pop_front();
      buffer_deque[s].push_back(detector_buff.GetADC_value(s));
      for(Int_t y=0; y<deque_size; y++)
      {
         RMS_buffer.SetChannelValue(buffer_deque[s][y],y);
         RMS_buffer.SetChannelValueSquared(buffer_deque[s][y],y);
      }
   store->SetRMSValues(s,RMS_buffer.CalculateSigma()); */ //old code

   Int_t detector_buff_ADC_value, buffer_front;
   Float_t new_mean, new_sigma, store_Ped_Value, store_RMS_Value;
   Double_t Remove, Add;
   for (Int_t s = 0; s < channel_number; s++) {
     detector_buff_ADC_value = detector_buff.GetADC_value(s);
     store_Ped_Value = store->GetPedValues(s);
     store_RMS_Value = store->GetRMSValues(s);

     if (Abs(detector_buff_ADC_value - store_Ped_Value) > threshold_factor * store_RMS_Value)
       continue;

     buffer_front = buffer_deque[s][0];
     buffer_deque[s].pop_front();
     buffer_deque[s].push_back(detector_buff_ADC_value);

     if (event % Taylor_speed_throttle) {
       Remove = buffer_front / divisor;
       Add = detector_buff_ADC_value / divisor;

       new_mean = store_Ped_Value + Add - Remove;

       new_sigma = Sqrt(store_RMS_Value * store_RMS_Value * divisor * divisor + (divisor - 1.) * (detector_buff_ADC_value + buffer_front) * (detector_buff_ADC_value - buffer_front) - 2 * (divisor * store_Ped_Value - buffer_front) * (detector_buff_ADC_value - buffer_front)) / divisor;
       //new_sigma = Sqrt(store_RMS_Value * store_RMS_Value + (divisor - 1.) * (Add * Add - Remove * Remove) - 2 * (store_Ped_Value - Remove) * (Add - Remove)); //Taylor: less precise value (maximum difference about 0.3%), but accurate and much faster

       //compare sigma values
       /* if (s == (Int_t)plottedChannel) {
         TEvent_Array RMS_buffer(deque_size);

         for (Int_t y = 0; y < deque_size; y++)
           RMS_buffer.SetChannelValue(buffer_deque[s][y], y);

         Float_t old_sigma = RMS_buffer.CalculateSigma();

         hRMSDifference->Fill((new_sigma - old_sigma) / old_sigma);
       } // */ //Taylor
     } else {
       TEvent_Array RMS_buffer(deque_size);

       for (Int_t y = 0; y < deque_size; y++)
         RMS_buffer.SetChannelValue(buffer_deque[s][y], y);

       new_sigma = RMS_buffer.CalculateSigma(new_mean);
     } //Taylor: calculate old way to make sure value doesn't drift away

     store->SetPedValues(s, new_mean);
     store->SetRMSValues(s, new_sigma);
     //new code


      if(store->GetRMSValues(s)< 0.1*(initial.GetRMSValues(s)) || store->GetRMSValues(s) == 0)  //0.1 is the scale factor for the initial value (brings it closer to zero, but makes sure it never hits zero)
      {
         store->SetRMSValues(s,0.1*(initial.GetRMSValues(s)));
         if(store->GetRMSValues(s)==0)
         {
            store->SetRMSValues(s,0.1);
         }
      }
      if(store->GetRMSValues(s)==0)
      {
         *zeroRMS = 1;
         cout << "RMS zeroed!!!" << " at event: " << event << " and channel " << s << endl;
         break;
      }
         
   }
}


void SlidingPedestal::RunningCommonMode(deque<Double_t> &CMN_deque, Double_t &CMN_Mean, Double_t &CMN_RMS, Double_t new_ave, Int_t deque_size, Int_t event, Int_t *zeroRMS)
{
   Double_t divisor = (Double_t) deque_size;
   Double_t RMS_buffer2 = 0;
   Double_t mean_buffer = CMN_Mean;


   // Change Mean
   Double_t Remove = CMN_deque[0]/divisor;
   Double_t Add = new_ave/divisor;
   CMN_Mean = mean_buffer+Add-Remove;
   CMN_deque.pop_front();
   CMN_deque.push_back(new_ave);

   // Change RMS
   for(Int_t i=0; i<deque_size; i++)
   {
      RMS_buffer2 = RMS_buffer2 + ((CMN_deque[i])*(CMN_deque[i]));
   }
   CMN_RMS = (Float_t) Sqrt((RMS_buffer2/divisor)-(CMN_Mean*CMN_Mean));
   if(CMN_RMS < 0.1)
   {
      CMN_RMS = 0.1;
   }
   if(CMN_RMS == 0)
   {
      *zeroRMS = 1;
      cout << "CMN RMS has zeroed at event " << event << endl;
   }
   //if(CMN_RMS !<>= 0.0)
   if(isnan(CMN_RMS))
   {
      *zeroRMS = 1;
      cout << "CMN RMS has diverged." << endl;
      cout << "Stats: " << endl;
      cout << "Mean: " << CMN_Mean  << endl;
      cout << "Divisor: " << divisor << endl;
      cout << "RMS_buffer2/divisor: " << RMS_buffer2/divisor << endl;
      cout << "Mean^2 : " << CMN_Mean*CMN_Mean << endl;
   }
}


void SlidingPedestal::Hit_Occupancy(ChannelScreen screen, TH1F *occup, TDetector_Data detector_buffer, TPed_and_RMS *ped_store, Float_t RMS_factor, Int_t chan_begin, Int_t chan_end, Int_t const dia_offset) {
   for(Int_t i = chan_begin+dia_offset; i<chan_end+dia_offset; i++)
   {
      if(screen.CheckChannel(i-dia_offset)==1)
      {
         if(TMath::Abs(detector_buffer.GetADC_value(i)-ped_store->GetPedValues(i)) > RMS_factor*ped_store->GetRMSValues(i))
         {
            occup->Fill(i-dia_offset);
         }
      }
   }
}


void SlidingPedestal::PedRMSCalcFromBuffer(vector< deque<Int_t> > &buffer_deque, TPed_and_RMS *ped_store) {
   Float_t Pedestal_check[256];
   Float_t Pedestal_check_squared[256];
   Float_t Pedestal_check_RMS[256];
   //Float_t Ped_check_sum = 0;
   //Float_t Ped_check_sum_squared = 0;
   Int_t Ped_check_sum = 0;
   Int_t Ped_check_sum_squared = 0;
   Int_t buffer_ij;
   for(Int_t i=0; i<256; i++)
   {
      Ped_check_sum = 0;
      Ped_check_sum_squared = 0;
      /* for(Int_t j=0; j<100; j++)
      {
         Ped_check_sum = Ped_check_sum + buffer_deque[i][j];
         Ped_check_sum_squared = Ped_check_sum_squared + (buffer_deque[i][j]*buffer_deque[i][j]);
      }
      Pedestal_check[i] = Ped_check_sum/100.0;
      Pedestal_check_squared[i] = Ped_check_sum_squared/100.0; // */ //Taylor

      for (Int_t j = 0; j < Iter_Size; j++) {
        buffer_ij = buffer_deque[i][j];
        Ped_check_sum += buffer_ij;
        Ped_check_sum_squared += buffer_ij * buffer_ij;
      }

      Pedestal_check[i] = Ped_check_sum / ((Float_t)Iter_Size);
      Pedestal_check_squared[i] = Ped_check_sum_squared / ((Float_t)Iter_Size);

      //RMS
      Pedestal_check_RMS[i] = Sqrt(Pedestal_check_squared[i]-(Pedestal_check[i]*Pedestal_check[i]));

      ped_store->SetPedValues(i,Pedestal_check[i]);
      ped_store->SetRMSValues(i,Pedestal_check_RMS[i]);
   }
}


void SlidingPedestal::Slide(Int_t Event_Number, Int_t Initial_Event, Int_t hit_occupancy) {
   //Setting Style Attributes for Plots (Starting with Plain style and then making custom adjustments)
   gROOT->SetStyle("Plain"); //General style (see TStyle)
   gStyle->SetOptStat(1110); //Stat options to be displayed
   gStyle->SetOptFit(1111);  //Fit options to be displayed
   gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadTopMargin(0.15);
   //gStyle->SetTitleColor(19,"");
   gStyle->SetStatH(0.12); //Sets Height of Stats Box
   gStyle->SetStatW(0.15); //Sets Width of Stats Box

   //Start Timer//
   TStopwatch Watch;
   Watch.Start();
   TDatime dateandtime;
   ///////////////
   

   //Plot Tag attached to all plots: Gives Run Number, Cut Threshold, Initial Number of Events, and Date/Time of Plot Creation
   char thresh1 [50];
   //char thresh2 [50];
   char datasetsize [50];
   char *pthresh1 = &thresh1[0];
   //char *pthresh2 = &thresh2[0];
   char *pdatasetsize = &datasetsize[0];
   sprintf(thresh1,"%i.%i", (int)Di_Hit_Factor,(int)(Di_Hit_Factor*10)%10);
   //sprintf(thresh2,"%i", CLUS_FACTOR);
   sprintf(datasetsize,"%i", Event_Number);
   //char dash[10] = "-";
   char space[10] = " ";
   char eventsinset[50] = " Events in Data Set";
   //char cut[10] = " cut | ";
   char sigma[50] = " Sig Cut";
   //char *pdash = &dash[0];
   char *pspace = &space[0];
   char *peventsinset = &eventsinset[0];
   //char *pcut = &cut[0];
   char *psigma = &sigma[0];
   strcat(pthresh1,psigma);
   //strcat(pthresh1,pthresh2);  //char of theshold i.e. "5-3"
   //strcat(pthresh1,pcut);
   strcat(pthresh1,pspace);
   strcat(pthresh1,pdatasetsize);
   strcat(pthresh1,peventsinset);
   //char run[10] = "Run ";
   //char *prun = &run[0];
   //char *run_number_label = &run_number[0];
   //strcat(prun,run_number_label);
   std::ostringstream run_number_label;
   run_number_label << "Run " << run_number;
   TPaveText *pt = new TPaveText(0.07,0,0.22,0.10,"NDC");  //Normalized CoordinateSystem: Define with x1,y1 is left bottom of box text, x2,y2 is upper right of text box. Goes from 0,0 at bottom left corner of pad to 1,1 of upper right corner
   pt->SetTextSize(0.0250);
   pt->AddText(run_number_label.str().c_str());
   pt->AddText(pthresh1);
   pt->AddText(dateandtime.AsSQLString());
   pt->SetBorderSize(0); //Set Border to Zero
   pt->SetFillColor(0); //Set Fill to White



   //Data is read in with function ReadRawEvent and stored as class member data




   //**********BEGINNING OF PEDESTAL CALC*******************************************//

   deque<TTrigger_Event> Events_deque;
   //Int_t const Iter_Size = 100; //because of def of TEvent_Array in Event_Classes.h this cannot be more than 300
   Int_t const Channel_Number = 256;
   Int_t const Iteration_Number = 6;
   //Detector_Data test_data;


//Channel Screen
   ChannelScreen screen;

   Int_t screen_array_size = Det_channel_screen_channels[8].size();//sizeof(Det_channel_screen_channels[8])/sizeof(Det_channel_screen_channels[8][0]);
   cout << "Number of Channels Screened is " << screen_array_size << endl;
   for(Int_t i=0; i<screen_array_size; i++)
   {
      screen.ChannelOff(Det_channel_screen_channels[8][i]);
   }

   for(Int_t i=Initial_Event; i<(Initial_Event+Iter_Size); i++)  //Initialzation begins at the first event
   {
      ReadRawEvent(i);
      TTrigger_Event Event;
      Event.SetD0X(D0X);
      Event.SetD0Y(D0Y);
      Event.SetD1X(D1X);
      Event.SetD1Y(D1Y);
      Event.SetD2X(D2X);
      Event.SetD2Y(D2Y);
      Event.SetD3X(D3X);
      Event.SetD3Y(D3Y);
      Event.SetDia0(Dia0);
      Event.SetDia1(Dia1);
      Events_deque.push_back(Event);
      //cout << i << "\t" << test_data.GetADC_value(161) << endl;
   }
   
   TPed_and_RMS *Initial_D0X = new TPed_and_RMS;
   TPed_and_RMS *Initial_D0Y = new TPed_and_RMS;
   TPed_and_RMS *Initial_D1X = new TPed_and_RMS;
   TPed_and_RMS *Initial_D1Y = new TPed_and_RMS;
   TPed_and_RMS *Initial_D2X = new TPed_and_RMS;
   TPed_and_RMS *Initial_D2Y = new TPed_and_RMS;
   TPed_and_RMS *Initial_D3X = new TPed_and_RMS;
   TPed_and_RMS *Initial_D3Y = new TPed_and_RMS;
   TPed_and_RMS *Initial_Dia0 = new TPed_and_RMS;
   TPed_and_RMS *Initial_Dia1 = new TPed_and_RMS;

   for(Int_t s=0; s<Channel_Number; s++) //loop over  # of channels in each detector
   {
         
      TEvent_Array *EA_D0X = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D0Y = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D1X = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D1Y = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D2X = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D2Y = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D3X = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_D3Y = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_Dia0 = new TEvent_Array(Iter_Size);
      TEvent_Array *EA_Dia1 = new TEvent_Array(Iter_Size);

      for(Int_t n=0; n<Iter_Size; n++)
      {
         int secretchan = 254;
         if(s==secretchan) cout<<n<<"\t";
         TTrigger_Event Event;
         Event = Events_deque[n];
         TDetector_Data Detector;
         Detector = Event.GetD0X();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D0X->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD0Y();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D0Y->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD1X();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D1X->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD1Y();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D1Y->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD2X();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D2X->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD2Y();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D2Y->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD3X();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D3X->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetD3Y();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_D3Y->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetDia0();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<"\t";
         EA_Dia0->SetChannelValue(Detector.GetADC_value(s),n);
         Detector = Event.GetDia1();
         if(s==secretchan) cout<<Detector.GetADC_value(s)<<endl;
         EA_Dia1->SetChannelValue(Detector.GetADC_value(s),n);
            
      }
      for(Int_t j=0; j<Iteration_Number; j++)
      {
         /*
         TH1F *D0X_pedestal = new TH1F("D0X_pedestal", "DOX Pedestal value for each channel",5000,0,5000);
         TH1F *D0Y_pedestal = new TH1F("D0Y_pedestal", "D0Y Pedestal value for each channel",5000,0,5000);
         TH1F *D1X_pedestal = new TH1F("D1X_pedestal", "D1X Pedestal value for each channel",5000,0,5000);
         TH1F *D1Y_pedestal = new TH1F("D1Y_pedestal", "D1Y Pedestal value for each channel",5000,0,5000);
         TH1F *D2X_pedestal = new TH1F("D2X_pedestal", "D2X Pedestal value for each channel",5000,0,5000);
         TH1F *D2Y_pedestal = new TH1F("D2Y_pedestal", "D2Y Pedestal value for each channel",5000,0,5000);
         TH1F *D3X_pedestal = new TH1F("D3X_pedestal", "D3X Pedestal value for each channel",5000,0,5000);
         TH1F *D3Y_pedestal = new TH1F("D3Y_pedestal", "D3Y Pedestal value for each channel",5000,0,5000);
         TH1F *Dia0_pedestal = new TH1F("Dia0_pedestal", "0 Pedestal value for each channel",5000,0,5000);
         TH1F *Dia1_pedestal = new TH1F("Dia1_pedestal", "1 Pedestal value for each channel",5000,0,5000);*/
         int silicon_adc_res = 255;
         int diamond_adc_res = 4095;
         TH1F *D0X_pedestal = new TH1F("D0X_pedestal", "DOX Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D0Y_pedestal = new TH1F("D0Y_pedestal", "D0Y Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D1X_pedestal = new TH1F("D1X_pedestal", "D1X Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D1Y_pedestal = new TH1F("D1Y_pedestal", "D1Y Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D2X_pedestal = new TH1F("D2X_pedestal", "D2X Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D2Y_pedestal = new TH1F("D2Y_pedestal", "D2Y Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D3X_pedestal = new TH1F("D3X_pedestal", "D3X Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *D3Y_pedestal = new TH1F("D3Y_pedestal", "D3Y Pedestal value for each channel",silicon_adc_res+1,0,silicon_adc_res);
         TH1F *Dia0_pedestal = new TH1F("Dia0_pedestal", "Dia0 Pedestal value for each channel",diamond_adc_res+1,0,diamond_adc_res);
         TH1F *Dia1_pedestal = new TH1F("Dia1_pedestal", "Dia1 Pedestal value for each channel",diamond_adc_res+1,0,diamond_adc_res);
            
         //cout << "For Channel " << s << " RMS: " << Initial_D0X->GetRMSValues(s) << " Ped: " << Initial_D0X->GetPedValues(s) << endl;
         PedIteration(D0X_pedestal, EA_D0X, Initial_D0X, s, Si_Hit_Factor);
         PedIteration(D0Y_pedestal, EA_D0Y, Initial_D0Y, s, Si_Hit_Factor);
         PedIteration(D1X_pedestal, EA_D1X, Initial_D1X, s, Si_Hit_Factor);
         PedIteration(D1Y_pedestal, EA_D1Y, Initial_D1Y, s, Si_Hit_Factor);
         PedIteration(D2X_pedestal, EA_D2X, Initial_D2X, s, Si_Hit_Factor);
         PedIteration(D2Y_pedestal, EA_D2Y, Initial_D2Y, s, Si_Hit_Factor);
         PedIteration(D3X_pedestal, EA_D3X, Initial_D3X, s, Si_Hit_Factor);
         PedIteration(D3Y_pedestal, EA_D3Y, Initial_D3Y, s, Si_Hit_Factor);
         PedIteration(Dia0_pedestal, EA_Dia0, Initial_Dia0, s, Di_Hit_Factor);
         PedIteration(Dia1_pedestal, EA_Dia1, Initial_Dia1, s, Di_Hit_Factor);

         delete D0X_pedestal;
         delete D0Y_pedestal;
         delete D1X_pedestal;
         delete D1Y_pedestal;
         delete D2X_pedestal;
         delete D2Y_pedestal;
         delete D3X_pedestal;
         delete D3Y_pedestal;
         delete Dia0_pedestal;
         delete Dia1_pedestal;
         D0X_pedestal = 0;
         D0Y_pedestal = 0;
         D1X_pedestal = 0;
         D1Y_pedestal = 0;
         D2X_pedestal = 0;
         D2Y_pedestal = 0;
         D3X_pedestal = 0;
         D3Y_pedestal = 0;
         Dia0_pedestal = 0;
         Dia1_pedestal = 0;
      }
      
      delete EA_D0X;
      delete EA_D0Y;
      delete EA_D1X;
      delete EA_D1Y;
      delete EA_D2X;
      delete EA_D2Y;
      delete EA_D3X;
      delete EA_D3Y;
      delete EA_Dia0;
      delete EA_Dia1;

      EA_D0X = 0;
      EA_D0Y = 0;
      EA_D1X = 0;
      EA_D1Y = 0;
      EA_D2X = 0;
      EA_D2Y = 0;
      EA_D3X = 0;
      EA_D3Y = 0;
      EA_Dia0 = 0;
      EA_Dia1 = 0;

   }
   for(Int_t i=0; i<128; i++)
   {
      //cout << "Channel: " << i << " Pedestal: " << Initial_Dia0->GetPedValues(i) << " RMS: " << Initial_Dia0->GetRMSValues(i) << endl;
   }
   cout << endl << "Finished Initializing Pedestal and RMS..." << endl;

   Int_t size_scale[256];
   for(Int_t i=0; i<256; i++) {size_scale[i] = 0;}
   //Int_t buffer_marker = 0;
   //Int_t event_counter = Iter_Size;
   vector< deque<Int_t> > D0Xbuffer(256);
//   vector< deque<Int_t> > *pD0Xbuffer(256) = &D0Xbuffer;
   vector< deque<Int_t> > D0Ybuffer(256);
   vector< deque<Int_t> > D1Xbuffer(256);
   vector< deque<Int_t> > D1Ybuffer(256);
   vector< deque<Int_t> > D2Xbuffer(256);
   vector< deque<Int_t> > D2Ybuffer(256);
   vector< deque<Int_t> > D3Xbuffer(256);
   vector< deque<Int_t> > D3Ybuffer(256);
   vector< deque<Int_t> > Dia0buffer(256);
   vector< deque<Int_t> > Dia1buffer(256);

   //Buffer Filling should begin at the same starting event as the Initialization of the Pedestal
   
   BufferFill(D0X, Initial_D0X, &D0Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event); 
   BufferFill(D0Y, Initial_D0Y, &D0Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D1X, Initial_D1X, &D1Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D1Y, Initial_D1Y, &D1Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D2X, Initial_D2X, &D2Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D2Y, Initial_D2Y, &D2Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D3X, Initial_D3X, &D3Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(D3Y, Initial_D3Y, &D3Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, Initial_Event);
   BufferFill(Dia0, Initial_Dia0, &Dia0buffer, Channel_Number, Iter_Size, Di_Hit_Factor, Initial_Event);
   BufferFill(Dia1, Initial_Dia1, &Dia1buffer, Channel_Number, Iter_Size, Di_Hit_Factor, Initial_Event);

   cout << endl << "Finished Filling Deques for Pedestal and RMS..." << endl;

   PedRMSCalcFromBuffer(D0Xbuffer, Initial_D0X);
   PedRMSCalcFromBuffer(D0Ybuffer, Initial_D0Y);
   PedRMSCalcFromBuffer(D1Xbuffer, Initial_D1X);
   PedRMSCalcFromBuffer(D1Ybuffer, Initial_D1Y);
   PedRMSCalcFromBuffer(D2Xbuffer, Initial_D2X);
   PedRMSCalcFromBuffer(D2Ybuffer, Initial_D2Y);
   PedRMSCalcFromBuffer(D3Xbuffer, Initial_D3X);
   PedRMSCalcFromBuffer(D3Ybuffer, Initial_D3Y);
   PedRMSCalcFromBuffer(Dia0buffer, Initial_Dia0);
   PedRMSCalcFromBuffer(Dia1buffer, Initial_Dia1);


   //Initialize Common Mode Noise Cut -Take average pedestal value over diamond channels and then average these over a sample of events (typically 100) to get an RMS for a cutoff value
   deque<Double_t> Common_Mode_Dia0_deque;
   deque<Double_t> Common_Mode_Dia1_deque;
   deque<Double_t> Common_Mode_Dia0_deque2;
   deque<Double_t> Common_Mode_Dia1_deque2;
   Double_t Common_Mode_Dia0_Mean;
   Double_t Common_Mode_Dia1_Mean;
   Double_t Common_Mode_Dia0_RMS;
   Double_t Common_Mode_Dia1_RMS;

   //Placeholders for analysis
   Double_t Common_Mode_Mean = 0;
   Double_t Common_Mode_RMS = 0;


   for(Int_t i=0; i<Iter_Size; i++)
   {
      Double_t CMN_buffer1=0;
      Double_t CMN_buffer2=0;
      Double_t divisor1 = (Double_t) 128;
      Double_t divisor2 = (Double_t) 128;
      for(Int_t s=0+DIA_OFFSET; s<(128+DIA_OFFSET); s++)
      {
         if(screen.CheckChannel(s)==0)
         {
            divisor1--;
            continue;
         }
         CMN_buffer1 = CMN_buffer1 + Dia0buffer[s][i];
      }
      for(Int_t s=0+DIA_OFFSET; s<(128+DIA_OFFSET); s++)
      {
         if(screen.CheckChannel(s)==0)
         {
            divisor2--;
            continue;
         }
         CMN_buffer2 = CMN_buffer2 + Dia1buffer[s][i];
      }
      Common_Mode_Dia0_deque.push_back(CMN_buffer1/divisor1);
      Common_Mode_Dia1_deque.push_back(CMN_buffer2/divisor2);
      //cout << endl << endl << endl << "Test" << endl;
      if(dia_input == 0){
         cout << "Diamond Input 0 Pedestal Average = " << CMN_buffer1/divisor1 << " or " << CMN_buffer1 << "/" << divisor1 << endl;}
         if(dia_input == 1){
            cout << "Diamond Input 1 Pedestal Average = " << CMN_buffer2/divisor2 << " or " << CMN_buffer2 << "/" << divisor2 << endl;}
            Common_Mode_Dia0_deque2.push_back((CMN_buffer1/divisor1)*(CMN_buffer1/divisor1));
            Common_Mode_Dia1_deque2.push_back((CMN_buffer2/divisor2)*(CMN_buffer2/divisor2));
   }
   Double_t Sum0 = 0;
   Double_t Sum1 = 0;
   Double_t Sum02 = 0;
   Double_t Sum12 = 0;
   for(Int_t i=0; i<Iter_Size; i++)
   {
      
      Sum0 = Sum0 + Common_Mode_Dia0_deque[i];
      Sum1 = Sum1 + Common_Mode_Dia1_deque[i];
      Sum02 = Sum02 + Common_Mode_Dia0_deque2[i];
      Sum12 = Sum12 + Common_Mode_Dia1_deque2[i];
   }
   Double_t div = (Double_t) Iter_Size;
   Common_Mode_Dia0_Mean = Sum0/div;
   Common_Mode_Dia1_Mean = Sum1/div;
   Common_Mode_Dia0_RMS = Sqrt((Sum02/div)-((Sum0/div)*(Sum0/div)));
   Common_Mode_Dia1_RMS = Sqrt((Sum12/div)-((Sum1/div)*(Sum1/div)));

   
   
   if(dia_input == 0)
   {
      cout << endl << "Common Mode Noise Mean: " << Common_Mode_Dia0_Mean << endl;
      cout << "Common Mode Noise RMS: " << Common_Mode_Dia0_RMS << endl;
      Common_Mode_Mean = Common_Mode_Dia0_Mean;
      Common_Mode_RMS = Common_Mode_Dia0_RMS;
   }
   if(dia_input == 1)
   {
      cout << endl << "Common Mode Noise Mean: " << Common_Mode_Dia1_Mean << endl;
      cout << "Common Mode Noise RMS: " << Common_Mode_Dia1_RMS << endl;
      Common_Mode_Mean = Common_Mode_Dia1_Mean;
      Common_Mode_RMS = Common_Mode_Dia1_RMS;
   }

      
   cout << endl << "Finished Initializing Common Mode Noise..." << endl;

   TPed_and_RMS *Ped_Store_D0X = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D0Y = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D1X = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D1Y = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D2X = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D2Y = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D3X = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_D3Y = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_Dia0 = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_Dia1 = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store = new TPed_and_RMS;
   TPed_and_RMS *Ped_Store_fixed_noise = new TPed_and_RMS;
   TPed_and_RMS Initialized_Values_D0X;
   TPed_and_RMS Initialized_Values_D0Y;
   TPed_and_RMS Initialized_Values_D1X;
   TPed_and_RMS Initialized_Values_D1Y;
   TPed_and_RMS Initialized_Values_D2X;
   TPed_and_RMS Initialized_Values_D2Y;
   TPed_and_RMS Initialized_Values_D3X;
   TPed_and_RMS Initialized_Values_D3Y;
   TPed_and_RMS Initialized_Values_Dia0;
   TPed_and_RMS Initialized_Values_Dia1;
   Initialized_Values_D0X = *Initial_D0X;
   Initialized_Values_D0Y = *Initial_D0Y;
   Initialized_Values_D1X = *Initial_D1X;
   Initialized_Values_D1Y = *Initial_D1Y;
   Initialized_Values_D2X = *Initial_D2X;
   Initialized_Values_D2Y = *Initial_D2Y;
   Initialized_Values_D3X = *Initial_D3X;
   Initialized_Values_D3Y = *Initial_D3Y;
   Initialized_Values_Dia0 = *Initial_Dia0;
   Initialized_Values_Dia1 = *Initial_Dia1;
   Ped_Store_D0X = Initial_D0X;
   Ped_Store_D0Y = Initial_D0Y;
   Ped_Store_D1X = Initial_D1X;
   Ped_Store_D1Y = Initial_D1Y;
   Ped_Store_D2X = Initial_D2X;
   Ped_Store_D2Y = Initial_D2Y;
   Ped_Store_D3X = Initial_D3X;
   Ped_Store_D3Y = Initial_D3Y;
   Ped_Store_Dia0 = Initial_Dia0;
   Ped_Store_Dia1 = Initial_Dia1;
   if(dia_input == 0)
   {
      Ped_Store = Initial_Dia0;
   }
   if(dia_input == 1)
   {
      Ped_Store = Initial_Dia1;
   }
   TTrigger_Event Event_buffer;
   TDetector_Data detector_buffer;
   TPed_and_RMS Ped_Buffer;
   TEvent_Array RMS_buffer(Iter_Size);
   TTrigger_Event New_Event;
   Int_t rmscount = 0;
   Int_t *prmscount = &rmscount;

   //store output pedestal subtracted data
   TFile *PedFile = new TFile(pedfilepath.str().c_str(),"recreate");
   TTree *PedTree = new TTree("PedTree","A Tree for Pedestal Subtracted Data");
//   PedTree->SetMaxTreeSize(8000000000); 
//   PedTree->SetMaxTreeSize(Long64_t(TMath::Power(2,63))); // http://lists.healthgrid.org/archives/gate-users/2004-December/000396.html
   PedTree->SetMaxTreeSize(Long64_t(TMath::Power(2,45)));  // Gives ~1TB (actually 4TB) max file size
   //PedTree->Branch("EventBranch","PSEvent",&store_event,24000,0); //address must be address of a pointer to a stored event
   
   //Event Header Branches
   PedTree->Branch("RunNumber",&run_number,"RunNumber/i");
   PedTree->Branch("EventNumber",&event_number,"EventNumber/i");
   PedTree->Branch("StoreThreshold",&store_threshold,"StoreThreshold/F");
   PedTree->Branch("CMNEvent_flag",&CMNEvent_flag,"CMNEvent_flag/O");
   PedTree->Branch("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag,"ZeroDivisorEvent_flag/O");
   
   //PedTree->Branch("EventNumber",&PedSubEvent.event_number,"EventNumber/i");
   //PedTree->Branch("StoreThreshold",&PedSubEvent.store_threshold,"StoreThreshold/F");

   //Telescope Data Branches
   PedTree->Branch("D0X_NChannels",&Det_NChannels[0],"D0X_NChannels/i");
   PedTree->Branch("D0Y_NChannels",&Det_NChannels[1],"D0Y_NChannels/i");
   PedTree->Branch("D1X_NChannels",&Det_NChannels[2],"D1X_NChannels/i");
   PedTree->Branch("D1Y_NChannels",&Det_NChannels[3],"D1Y_NChannels/i");
   PedTree->Branch("D2X_NChannels",&Det_NChannels[4],"D2X_NChannels/i");
   PedTree->Branch("D2Y_NChannels",&Det_NChannels[5],"D2Y_NChannels/i");
   PedTree->Branch("D3X_NChannels",&Det_NChannels[6],"D3X_NChannels/i");
   PedTree->Branch("D3Y_NChannels",&Det_NChannels[7],"D3Y_NChannels/i");
   PedTree->Branch("Dia_NChannels",&Det_NChannels[8],"Dia_NChannels/i");
   PedTree->Branch("D0X_Channels",&Det_Channels[0],"D0X_Channels[D0X_NChannels]/b");
   PedTree->Branch("D0Y_Channels",&Det_Channels[1],"D0Y_Channels[D0Y_NChannels]/b");
   PedTree->Branch("D1X_Channels",&Det_Channels[2],"D1X_Channels[D1X_NChannels]/b");
   PedTree->Branch("D1Y_Channels",&Det_Channels[3],"D1Y_Channels[D1Y_NChannels]/b");
   PedTree->Branch("D2X_Channels",&Det_Channels[4],"D2X_Channels[D2X_NChannels]/b");
   PedTree->Branch("D2Y_Channels",&Det_Channels[5],"D2Y_Channels[D2Y_NChannels]/b");
   PedTree->Branch("D3X_Channels",&Det_Channels[6],"D3X_Channels[D3X_NChannels]/b");
   PedTree->Branch("D3Y_Channels",&Det_Channels[7],"D3Y_Channels[D3Y_NChannels]/b");
   PedTree->Branch("Dia_Channels",&Det_Channels[8],"Dia_Channels[Dia_NChannels]/b");
   PedTree->Branch("D0X_ADC",&Det_ADC[0],"D0X_ADC[D0X_NChannels]/b");
   PedTree->Branch("D0Y_ADC",&Det_ADC[1],"D0Y_ADC[D0Y_NChannels]/b");
   PedTree->Branch("D1X_ADC",&Det_ADC[2],"D1X_ADC[D1X_NChannels]/b");
   PedTree->Branch("D1Y_ADC",&Det_ADC[3],"D1Y_ADC[D1Y_NChannels]/b");
   PedTree->Branch("D2X_ADC",&Det_ADC[4],"D2X_ADC[D2X_NChannels]/b");
   PedTree->Branch("D2Y_ADC",&Det_ADC[5],"D2Y_ADC[D2Y_NChannels]/b");
   PedTree->Branch("D3X_ADC",&Det_ADC[6],"D3X_ADC[D3X_NChannels]/b");
   PedTree->Branch("D3Y_ADC",&Det_ADC[7],"D3Y_ADC[D3Y_NChannels]/b");
   PedTree->Branch("Dia_ADC",&Dia_ADC,"Dia_ADC[Dia_NChannels]/s");
   PedTree->Branch("D0X_PedMean",&Det_PedMean[0],"D0X_PedMean[D0X_NChannels]/F");
   PedTree->Branch("D0Y_PedMean",&Det_PedMean[1],"D0Y_PedMean[D0Y_NChannels]/F");
   PedTree->Branch("D1X_PedMean",&Det_PedMean[2],"D1X_PedMean[D1X_NChannels]/F");
   PedTree->Branch("D1Y_PedMean",&Det_PedMean[3],"D1Y_PedMean[D1Y_NChannels]/F");
   PedTree->Branch("D2X_PedMean",&Det_PedMean[4],"D2X_PedMean[D2X_NChannels]/F");
   PedTree->Branch("D2Y_PedMean",&Det_PedMean[5],"D2Y_PedMean[D2Y_NChannels]/F");
   PedTree->Branch("D3X_PedMean",&Det_PedMean[6],"D3X_PedMean[D3X_NChannels]/F");
   PedTree->Branch("D3Y_PedMean",&Det_PedMean[7],"D3Y_PedMean[D3Y_NChannels]/F");
   PedTree->Branch("Dia_PedMean",&Det_PedMean[8],"Dia_PedMean[Dia_NChannels]/F");
   PedTree->Branch("D0X_PedWidth",&Det_PedWidth[0],"D0X_PedWidth[D0X_NChannels]/F");
   PedTree->Branch("D0Y_PedWidth",&Det_PedWidth[1],"D0Y_PedWidth[D0Y_NChannels]/F");
   PedTree->Branch("D1X_PedWidth",&Det_PedWidth[2],"D1X_PedWidth[D1X_NChannels]/F");
   PedTree->Branch("D1Y_PedWidth",&Det_PedWidth[3],"D1Y_PedWidth[D1Y_NChannels]/F");
   PedTree->Branch("D2X_PedWidth",&Det_PedWidth[4],"D2X_PedWidth[D2X_NChannels]/F");
   PedTree->Branch("D2Y_PedWidth",&Det_PedWidth[5],"D2Y_PedWidth[D2Y_NChannels]/F");
   PedTree->Branch("D3X_PedWidth",&Det_PedWidth[6],"D3X_PedWidth[D3X_NChannels]/F");
   PedTree->Branch("D3Y_PedWidth",&Det_PedWidth[7],"D3Y_PedWidth[D3Y_NChannels]/F");
   PedTree->Branch("Dia_PedWidth",&Det_PedWidth[8],"Dia_PedWidth[Dia_NChannels]/F");
   /*
   PedTree->Branch("D0X_NChannels",&PedSubEvent.pedsub_detector_data[0].nchannels,"D0X_NChannels/i");
   PedTree->Branch("D0Y_NChannels",&PedSubEvent.pedsub_detector_data[1].nchannels,"D0Y_NChannels/i");
   PedTree->Branch("D1X_NChannels",&PedSubEvent.pedsub_detector_data[2].nchannels,"D1X_NChannels/i");
   PedTree->Branch("D1Y_NChannels",&PedSubEvent.pedsub_detector_data[3].nchannels,"D1Y_NChannels/i");
   PedTree->Branch("D2X_NChannels",&PedSubEvent.pedsub_detector_data[4].nchannels,"D2X_NChannels/i");
   PedTree->Branch("D2Y_NChannels",&PedSubEvent.pedsub_detector_data[5].nchannels,"D2Y_NChannels/i");
   PedTree->Branch("D3X_NChannels",&PedSubEvent.pedsub_detector_data[6].nchannels,"D3X_NChannels/i");
   PedTree->Branch("D3Y_NChannels",&PedSubEvent.pedsub_detector_data[7].nchannels,"D3Y_NChannels/i");
   PedTree->Branch("Dia_NChannels",&PedSubEvent.pedsub_detector_data[8].nchannels,"Dia_NChannels/i");
   PedTree->Branch("D0X_Channels",&PedSubEvent.pedsub_detector_data[0].channel,"D0X_Channels[D0X_NChannels]/b");
   PedTree->Branch("D0Y_Channels",&PedSubEvent.pedsub_detector_data[1].channel,"D0Y_Channels[D0Y_NChannels]/b");
   PedTree->Branch("D1X_Channels",&PedSubEvent.pedsub_detector_data[2].channel,"D1X_Channels[D1X_NChannels]/b");
   PedTree->Branch("D1Y_Channels",&PedSubEvent.pedsub_detector_data[3].channel,"D1Y_Channels[D1Y_NChannels]/b");
   PedTree->Branch("D2X_Channels",&PedSubEvent.pedsub_detector_data[4].channel,"D2X_Channels[D2X_NChannels]/b");
   PedTree->Branch("D2Y_Channels",&PedSubEvent.pedsub_detector_data[5].channel,"D2Y_Channels[D2Y_NChannels]/b");
   PedTree->Branch("D3X_Channels",&PedSubEvent.pedsub_detector_data[6].channel,"D3X_Channels[D3X_NChannels]/b");
   PedTree->Branch("D3Y_Channels",&PedSubEvent.pedsub_detector_data[7].channel,"D3Y_Channels[D3Y_NChannels]/b");
   PedTree->Branch("Dia_Channels",&PedSubEvent.pedsub_detector_data[8].channel,"Dia_Channels[Dia_NChannels]/b");
   PedTree->Branch("D0X_ADC",&PedSubEvent.pedsub_detector_data[0].adc,"D0X_ADC[D0X_NChannels]/b");
   PedTree->Branch("D0Y_ADC",&PedSubEvent.pedsub_detector_data[1].adc,"D0Y_ADC[D0Y_NChannels]/b");
   PedTree->Branch("D1X_ADC",&PedSubEvent.pedsub_detector_data[2].adc,"D1X_ADC[D1X_NChannels]/b");
   PedTree->Branch("D1Y_ADC",&PedSubEvent.pedsub_detector_data[3].adc,"D1Y_ADC[D1Y_NChannels]/b");
   PedTree->Branch("D2X_ADC",&PedSubEvent.pedsub_detector_data[4].adc,"D2X_ADC[D2X_NChannels]/b");
   PedTree->Branch("D2Y_ADC",&PedSubEvent.pedsub_detector_data[5].adc,"D2Y_ADC[D2Y_NChannels]/b");
   PedTree->Branch("D3X_ADC",&PedSubEvent.pedsub_detector_data[6].adc,"D3X_ADC[D3X_NChannels]/b");
   PedTree->Branch("D3Y_ADC",&PedSubEvent.pedsub_detector_data[7].adc,"D3Y_ADC[D3Y_NChannels]/b");
   PedTree->Branch("Dia_ADC",&PedSubEvent.pedsub_detector_data[8].adc,"Dia_ADC[Dia_NChannels]/s");
   PedTree->Branch("D0X_PedMean",&PedSubEvent.pedsub_detector_data[0].pedmean,"D0X_PedMean[D0X_NChannels]/F");
   PedTree->Branch("D0Y_PedMean",&PedSubEvent.pedsub_detector_data[1].pedmean,"D0Y_PedMean[D0Y_NChannels]/F");
   PedTree->Branch("D1X_PedMean",&PedSubEvent.pedsub_detector_data[2].pedmean,"D1X_PedMean[D1X_NChannels]/F");
   PedTree->Branch("D1Y_PedMean",&PedSubEvent.pedsub_detector_data[3].pedmean,"D1Y_PedMean[D1Y_NChannels]/F");
   PedTree->Branch("D2X_PedMean",&PedSubEvent.pedsub_detector_data[4].pedmean,"D2X_PedMean[D2X_NChannels]/F");
   PedTree->Branch("D2Y_PedMean",&PedSubEvent.pedsub_detector_data[5].pedmean,"D2Y_PedMean[D2Y_NChannels]/F");
   PedTree->Branch("D3X_PedMean",&PedSubEvent.pedsub_detector_data[6].pedmean,"D3X_PedMean[D3X_NChannels]/F");
   PedTree->Branch("D3Y_PedMean",&PedSubEvent.pedsub_detector_data[7].pedmean,"D3Y_PedMean[D3Y_NChannels]/F");
   PedTree->Branch("Dia_PedMean",&PedSubEvent.pedsub_detector_data[8].pedmean,"Dia_PedMean[Dia_NChannels]/F");
   PedTree->Branch("D0X_PedWidth",&PedSubEvent.pedsub_detector_data[0].pedwidth,"D0X_PedWidth[D0X_NChannels]/F");
   PedTree->Branch("D0Y_PedWidth",&PedSubEvent.pedsub_detector_data[1].pedwidth,"D0Y_PedWidth[D0Y_NChannels]/F");
   PedTree->Branch("D1X_PedWidth",&PedSubEvent.pedsub_detector_data[2].pedwidth,"D1X_PedWidth[D1X_NChannels]/F");
   PedTree->Branch("D1Y_PedWidth",&PedSubEvent.pedsub_detector_data[3].pedwidth,"D1Y_PedWidth[D1Y_NChannels]/F");
   PedTree->Branch("D2X_PedWidth",&PedSubEvent.pedsub_detector_data[4].pedwidth,"D2X_PedWidth[D2X_NChannels]/F");
   PedTree->Branch("D2Y_PedWidth",&PedSubEvent.pedsub_detector_data[5].pedwidth,"D2Y_PedWidth[D2Y_NChannels]/F");
   PedTree->Branch("D3X_PedWidth",&PedSubEvent.pedsub_detector_data[6].pedwidth,"D3X_PedWidth[D3X_NChannels]/F");
   PedTree->Branch("D3Y_PedWidth",&PedSubEvent.pedsub_detector_data[7].pedwidth,"D3Y_PedWidth[D3Y_NChannels]/F");
   PedTree->Branch("Dia_PedWidth",&PedSubEvent.pedsub_detector_data[8].pedwidth,"Dia_PedWidth[Dia_NChannels]/F");
    */
   
   ofstream pedestal_file("running_ped.txt");
   ofstream average_ped_file("average_ped.txt");
   ofstream ped_subtracted_file("ped_subtracted.txt");
   ofstream runstats_file("run_statistics.txt");
   Int_t Output_Events = 0;
 
   cout << endl; 
   cout << "Beginning Event Input - Number of Read-in Events: " << Event_Number << endl;
   cout << endl;

   Int_t Triggered_Events = 0; //Number of Triggers
   Int_t Bad_Events = 0; //Number of events which have a zero divisor
   Int_t ZeroDivisor_Events = 0; //Number of events which have a zero divisor
   Int_t CMN_Events = 0; // number of cmn events
   Int_t Events_Sent_Into_Corrector = 0; //Number of events which are remaining after removal of bad events
   Int_t Corrected_Events = 0; //Number of Events that get corrected
   Int_t PostCorr_Bad_Events = 0; //Number of events which have a zero divisor after the CMN correction (should remain zero in general)
   Int_t PostCorr_Events_Outside_CMN_win = 0; //Events that are cut after the correction
   Int_t Events_Saved = 0; //Final Number of Events saved 

   //Int_t Events_Outside_Large_Win = 0; //Number of Events whose Pedestal Average fell outside the large window
   //Int_t Events_Inside_Large_Win = 0; //Number of Events whose Pedestal Average was inside large window and was used to calc running Pedestal


   Int_t CMN_upper_range = (int) (/*Common_Mode_Dia1_Mean*/0 + 100);
   Int_t CMN_lower_range = (int) (/*Common_Mode_Dia1_Mean*/0 - 100);

   TH1F *CMN_noise = new TH1F("CMN_noise", "CMN_Noise", 400,CMN_lower_range,CMN_upper_range);
   CMN_noise->SetDirectory(0);

   TH1F *CMN_noise_saved = new TH1F("CMN_noise_saved", "CMN_Noise Saved", 400,CMN_lower_range,CMN_upper_range);
   CMN_noise_saved->SetDirectory(0);

   TH1F *noise = new TH1F("noise", "Noise", 200,-100,100);
   noise->SetDirectory(0);

   TH1F *corr_dist = new TH1F("corr_dist", "Distribution of Corrections", 200,-100,100);
   corr_dist->SetDirectory(0);

   TH1F *CMN_RMS_dist = new TH1F("CMN_RMS_dist", "CMN RMS Distribution", 100,0,40);
   CMN_RMS_dist->SetDirectory(0);


   TH1F *hit_occup = new TH1F("hit_occup","Diamond Hit Occupancy",128,-0.5,127.5);
   hit_occup->SetDirectory(0);

   TH2F *raw_ADC_by_event = new TH2F("raw_ADC_by_event","Raw ADC By Event",Event_Number,Initial_Event,Initial_Event+Event_Number,400,CMN_lower_range,CMN_upper_range);
   raw_ADC_by_event->SetDirectory(0);

   TGraph *raw_ADC_by_event_graph = new TGraph(Event_Number+1);
   TGraph *RMS_threshold_by_event_graph_up = new TGraph(Event_Number+1);
   TGraph *RMS_threshold_by_event_graph_down = new TGraph(Event_Number+1);

   TGraph *raw_ADC_by_event_CMN_cut_graph = new TGraph(Event_Number+1);

   TGraph *PS_ADC_by_event_graph = new TGraph(Event_Number+1);
   TGraph *PS_RMS_threshold_by_event_graph_up = new TGraph(Event_Number+1);
   TGraph *PS_RMS_threshold_by_event_graph_down = new TGraph(Event_Number+1);

   TGraph *PS_ADC_by_event_CMN_cut_graph = new TGraph(Event_Number+1);
   
   Int_t event_index = (Initial_Event-1); //index for TGraphs
   
   //-----------------------------------------
   // Initializing the single channel analysis
   Int_t NumberSingleChannelAnalysisChannels = single_channel_analysis_channels.size();
   Int_t NumberSingleChannelAnalysisWindows = Event_Number / single_channel_analysis_eventwindow;
   Int_t NumberSingleChannelAnalysisHistos = NumberSingleChannelAnalysisChannels * NumberSingleChannelAnalysisWindows;
   
   TH1F* SingleChannelAnalysisNoiseHistos[NumberSingleChannelAnalysisHistos];

   for(int s=0; s<NumberSingleChannelAnalysisChannels; s++)
   {
      for(int t=0; t<NumberSingleChannelAnalysisWindows; t++)
      {
         string temporaryname = "Channel_Noise_ch"; // creates a string called temporaryname
         stringstream out;
         out << single_channel_analysis_channels[s] << "_window" << t;
         temporaryname += out.str();
         SingleChannelAnalysisNoiseHistos[s + t*NumberSingleChannelAnalysisChannels] = new TH1F(temporaryname.c_str(),temporaryname.c_str(), 120, -60, 60);  // creates a histogram and returns the location in memory of the histogram and assigns it to the sth pointer to a TH1F
      }
   }

   TH1F* SingleChannelAnalysisPulseHeightHistos[NumberSingleChannelAnalysisHistos];

   for(int s=0; s<NumberSingleChannelAnalysisChannels; s++)
   {
      for(int t=0; t<NumberSingleChannelAnalysisWindows; t++)
      {
         string temporaryname = "Channel_PulseHeight_ch"; // creates a string called temporaryname
         stringstream out;
         out << single_channel_analysis_channels[s] << "_window" << t;
         temporaryname += out.str();
         SingleChannelAnalysisPulseHeightHistos[s + t*NumberSingleChannelAnalysisChannels] = new TH1F(temporaryname.c_str(),temporaryname.c_str(), 120, 0, 1000);  // creates a histogram and returns the location in memory of the histogram and assigns it to the sth pointer to a TH1F
      }
   }
   
   // Finished initializing the single channel analysis
   //--------------------------------------------------
   
   
   //-----------------------------------------------
   // Starting the main sliding pedestal calculation
   for(Int_t e=Initial_Event+Iter_Size; e<(Initial_Event+Iter_Size+Event_Number); e++)  //The actual pedestal calculation should begin at the event just after that used for the Initialization and Buffer Fill (some overlap is tolerated due to the low frequency of hits
   {
      
      event_index++;
      Triggered_Events++;
      if(e%10000==0)
      {
         cout << "on loop " << e << "..." << endl;
         //Watch.Print("u");
      }
      ReadRawEvent(e);
      
      //Ped_Buffer = *Ped_Store_D0X;
      //TTrigger_Event *Old_Event = new TTrigger_Event;
      //*Old_Event = Events_deque[0];
      //TTrigger_Event *New_Event = new TTrigger_Event;
      New_Event.SetD0X(D0X);
      New_Event.SetD0Y(D0Y);
      New_Event.SetD1X(D1X);
      New_Event.SetD1Y(D1Y);
      New_Event.SetD2X(D2X);
      New_Event.SetD2Y(D2Y);
      New_Event.SetD3X(D3X);
      New_Event.SetD3Y(D3Y);
      New_Event.SetDia0(Dia0);
      New_Event.SetDia1(Dia1);

      // Common Mode Noise Test
      if(dia_input == 0)
      {
         detector_buffer = New_Event.GetDia0();
      }
      if(dia_input == 1)
      {
         detector_buffer = New_Event.GetDia1();
      }
      
      Double_t sum_buffer = 0;
      Double_t ps_sum_buffer = 0;
      //Double_t ps_divisor = 128;
      Double_t divisor = (Double_t) 128;
      Float_t correction = 0;
      
      for(Int_t s=0+DIA_OFFSET; s<(128+DIA_OFFSET); s++)
      {
	  
         if(Abs(detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s))>(Di_Hit_Factor*Ped_Store->GetRMSValues(s)) || screen.CheckChannel(s)==0)
         {
            divisor--;
            continue;
         }
         sum_buffer = sum_buffer + detector_buffer.GetADC_value(s);
         ps_sum_buffer = ps_sum_buffer + (detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s)); //Sum of Pedestal Subtracted channels
      }

      if(e==500 || e==3203 || e==4336 || e==4335 || e==4337)
      {
         cout << "for event " << e << " sum buffer/divisor is " << sum_buffer/divisor << " and ps sum buffer is " << ps_sum_buffer/divisor << endl;
      }
      

      CMNEvent_flag = 0;
      ZeroDivisorEvent_flag = 0;
            
      //Remove any bad events that have a divisor of 0
      if(divisor==0)
      {
         ZeroDivisorEvent_flag = 1;
         Bad_Events++;
         ZeroDivisor_Events++;
         event_index--; //keeps index of TGraph from skipping a value
         continue; //don't calculate CMN or save event because of a bad event and go onto next event
      }
     
      //Fill Pre-CMN cut Plots and files
      
      //Plots and Files
      if(divisor!=0) {
         CMN_noise->Fill(ps_sum_buffer/divisor);
         CMN_RMS_dist->Fill(Common_Mode_RMS);
         average_ped_file << sum_buffer/divisor << endl;
         //raw_ADC_by_event->Fill(e,sum_buffer/divisor);
         ped_subtracted_file << ps_sum_buffer/divisor << endl;
         //Graphs
         Double_t doub_event = (Double_t) e; //converts e into double for TGraph argument requirements
         raw_ADC_by_event_graph->SetPoint(event_index,doub_event,sum_buffer/divisor);
         RMS_threshold_by_event_graph_up->SetPoint(event_index,doub_event,Common_Mode_Mean+(Common_Mode_RMS*CMN_cut));
         RMS_threshold_by_event_graph_down->SetPoint(event_index,doub_event,Common_Mode_Mean-(Common_Mode_RMS*CMN_cut));
         PS_ADC_by_event_graph->SetPoint(event_index,doub_event,ps_sum_buffer/divisor);
         PS_RMS_threshold_by_event_graph_up->SetPoint(event_index,doub_event,(0+Common_Mode_RMS*CMN_cut));
         PS_RMS_threshold_by_event_graph_down->SetPoint(event_index,doub_event,0-(Common_Mode_RMS*CMN_cut));
         Events_Sent_Into_Corrector++;


         //Begin Correction of some events (i.e. correct events approx above 3-4sigma and below 6-7sigma)
         //if(Abs((sum_buffer/divisor)-(Common_Mode_Mean)) > ((CMN_corr_low)*(Common_Mode_RMS) && Abs((sum_buffer/divisor)-(Common_Mode_Mean) < (CMN_corr_high)*(Common_Mode_RMS)))) //&& Abs((sum_buffer/divisor)-Common_Mode_Dia1_Mean)<)
         if(Abs((sum_buffer/divisor)-Common_Mode_Mean) > ((CMN_corr_low)*Common_Mode_RMS) && Abs((sum_buffer/divisor)-Common_Mode_Mean) < (CMN_corr_high)*Common_Mode_RMS) //&& Abs((sum_buffer/divisor)-Common_Mode_Dia1_Mean)<)
         {
            //2010-08-08: just store as is and flag as CMN event.
            CMNEvent_flag = 1;
            CMN_Events++;
            
            
            Corrected_Events++;
            correction = (sum_buffer/divisor)-Common_Mode_Mean;  //difference between average pedestal of channels and average of this average
            corr_dist->Fill(correction);  //fill plot showing distribution of corrections

            /*for(Int_t i=0; i<256; i++)
            {
               //cout << detector_buffer.GetADC_value(i)-Ped_Store_Dia1->GetPedValues(i) << endl;
            }
            */
            continue; //go to next event (don't store)

            //Correct each non-noisy channel of the event with the above calculated correction by subtracting it from each channel
            for(Int_t s=0+DIA_OFFSET; s<(128+DIA_OFFSET); s++)
            {

               if(screen.CheckChannel(s)==1)
               {
                  detector_buffer.SetADC_value(s,detector_buffer.GetADC_value(s)-correction); //correction is subtracted for each channel
                  
               }
               
               
            }
            
         }
      }
      
      //Once the events needing correction have been corrected, we repeat the process
      sum_buffer = 0; //reset raw sum value
      ps_sum_buffer = 0; //reset ped_subtracted sum value
      divisor = 128; //reset divisor value

      //scanning through all channels
      for(Int_t s=0+DIA_OFFSET; s<(128+DIA_OFFSET); s++)
      {
	
         if(Abs(detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s))>(Di_Hit_Factor*Ped_Store->GetRMSValues(s)) || screen.CheckChannel(s)==0)
         {
            if(single_channel_analysis_enable){
		//we have a hit so fill the pulse height histo
               for(Int_t t = 0; t < NumberSingleChannelAnalysisChannels; t++)
               {
                  if(s-DIA_OFFSET == single_channel_analysis_channels[t])
                  {
                     SingleChannelAnalysisPulseHeightHistos[t + NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow) ]->Fill(detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s));
                  }
               }
            }

            divisor--;
            continue;
         }	  
         
         
         if(single_channel_analysis_enable){
		// we have no hit so fill the noise histo
            for(Int_t t = 0; t < NumberSingleChannelAnalysisChannels; t++)
            {
	//            if(e%100==0) cout<<"s-DIA_OFFSET = "<<s-DIA_OFFSET<<"\tt = "<<t<<"\tsingle_channel_analysis_channels[t] = "<<single_channel_analysis_channels[t]<<endl;
               if(s-DIA_OFFSET == single_channel_analysis_channels[t])
               {
	//               cout<<"t = "<<t<<"\tNumberSingleChannelAnalysisChannels = "<<NumberSingleChannelAnalysisChannels<<endl;
	//               cout<<"e = "<<e<<"\tInitial_Event = "<<Initial_Event<<"\tIter_Size = "<<Iter_Size<<"\tsingle_channel_analysis_eventwindow = "<<single_channel_analysis_eventwindow<<endl;
	//               cout<<"Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow) = "<<Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow)<<endl;
	//               cout<<"NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow) = "<<NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow)<<endl;
	//               cout<<"t + NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow) = "<<t + NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow)<<endl;
	
                  SingleChannelAnalysisNoiseHistos[t + NumberSingleChannelAnalysisChannels * Int_t((e-Initial_Event-Iter_Size) / single_channel_analysis_eventwindow) ]->Fill(detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s));
	//               if(e%100==0) cout<<"Fill: detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s) = "<<detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s)<<endl;
               }
            }
         }
         
         
         sum_buffer = sum_buffer + detector_buffer.GetADC_value(s);
         ps_sum_buffer = ps_sum_buffer + (detector_buffer.GetADC_value(s)-Ped_Store->GetPedValues(s)); //Sum of Pedestal Subtracted channels 
         
      }
      
      
      //2010-08-13: flag as bad event.
      //remove any remaining bad events (bad defined as no unmasked channels without hits
      if(divisor==0)
      {
         if(!ZeroDivisorEvent_flag) {
            ZeroDivisorEvent_flag = 1;
            ZeroDivisor_Events++;
         }
         
         PostCorr_Bad_Events++;
         event_index--; //keeps index of TGraph from skipping a value
         continue; //don't calculate CMN or save event because of a bad event and go onto next event
         
      }
      
      //2010-08-13: flag as CMN event.
      //remove any events now with a average pedestal greater than some sigma which should be equal to or less than the correction window (i.e. about 5 sigma)
      if(Abs((sum_buffer/divisor)-Common_Mode_Mean)>((CMN_cut)*Common_Mode_RMS))
      {
         if(!CMNEvent_flag) {
            CMNEvent_flag = 1;
            CMN_Events++;
         }
         
         
         PostCorr_Events_Outside_CMN_win++;
         event_index--; //keeps index of TGraph from skipping a value
         continue;
         
      }
      //cout << " and still is: " << detector_buffer.GetADC_value(57) << endl;
      //cout << " and sum_buffer/divisor after is: " << sum_buffer/divisor << " (" << sum_buffer << "/" << divisor << ")" << endl;

      //Fill plots and files and graphs
      
      CMN_noise_saved->Fill(ps_sum_buffer/divisor);
      //raw_ADC_by_event_CMN_cut_graph->SetPoint(event_index,doub_event,sum_buffer/divisor);
      //PS_ADC_by_event_CMN_cut_graph->SetPoint(event_index,doub_event,ps_sum_buffer/divisor);
      
      Events_Saved++;

      //if(sum_buffer/divisor<1640)
      //{
	  //cout << "Low Average Event at event: " << e << " with average: " << sum_buffer/divisor << " and mean and RMS at: " << Common_Mode_Dia1_Mean << " and " << Common_Mode_Dia1_RMS << endl;
      //}
 

      
      
      //detector_buffer = New_Event.GetDia1();
      
      //input << "Input ADC: " << "\t" << detector_buffer.GetADC_value(57+DIA_OFFSET) << endl;
      
      
      //fill various plots and files
      if(e%1==0)
      {
         Hit_Occupancy(screen, hit_occup, detector_buffer, Ped_Store, Di_Hit_Factor, 0, 128, DIA_OFFSET);
         for(Int_t i=(0+DIA_OFFSET); i<(128+DIA_OFFSET); i++)
         {
            if(TMath::Abs(detector_buffer.GetADC_value(i)-Ped_Store->GetPedValues(i) <= Di_Hit_Factor*Ped_Store->GetRMSValues(i)))
            {
               noise->Fill(detector_buffer.GetADC_value(i)-Ped_Store->GetPedValues(i));
            }
         }
         for(Int_t i=(75+DIA_OFFSET); i<(76+DIA_OFFSET); i++)
         {
            //if(detector_buffer.GetADC_value(i)-Ped_Store_Dia0->GetPedValues(i)>(Hit_Factor*Ped_Store_Dia0->GetRMSValues(i)))
            //{
            //   break;
            // }
            
            char channel_name[10];
            Int_t channel_number = 0;
            channel_number = sprintf(channel_name,"%i",i-DIA_OFFSET);
            //pedestal_file << channel_name << "\t" << detector_buffer.GetADC_value(i)-Ped_Store_Dia1->GetPedValues(i) << "\t" << detector_buffer.GetADC_value(i)-Initialized_Values_Dia1.GetPedValues(i)<< "\t" << Ped_Store_Dia1->GetPedValues(i) << "\t" << Initialized_Values_Dia1.GetPedValues(i) <<"\t" <<Ped_Store_Dia1->GetRMSValues(i) <<"\t" <<(2+Hit_Factor)*Ped_Store_Dia1->GetRMSValues(i) << "\t" << Initialized_Values_Dia1.GetRMSValues(i) << "\t" << Hit_Factor*Initialized_Values_Dia1.GetRMSValues(i) << "\t" << detector_buffer.GetADC_value(i) << "\t";
            
         }
         //pedestal_file << sum_buffer/divisor << endl;
      }
      for(Int_t i=61; i<62; i++)
      {
         //cout << "Event: " << e << " Channel: " << i << "Raw ADC is :" << detector_buffer.GetADC_value(i) << " Pedestal: " << Ped_Store_Dia1->GetPedValues(i) << " RMS: " << Ped_Store_Dia1->GetRMSValues(i) << endl;
      }

      
      //Store all current values for the event and write the the TTree
      SetDetector(0,D0X,Ped_Store_D0X);
      SetDetector(1,D0Y,Ped_Store_D0Y);
      SetDetector(2,D1X,Ped_Store_D1X);
      SetDetector(3,D1Y,Ped_Store_D1Y);
      SetDetector(4,D2X,Ped_Store_D2X);
      SetDetector(5,D2Y,Ped_Store_D2Y);
      SetDetector(6,D3X,Ped_Store_D3X);
      SetDetector(7,D3Y,Ped_Store_D3Y);
      if(fix_dia_noise>0) for(int ii=0; ii<256; ii++) {
         Ped_Store_fixed_noise->SetPedValues(ii,Ped_Store->GetPedValues(ii));
         Ped_Store_fixed_noise->SetRMSValues(ii,fix_dia_noise);
      }
      if(fix_dia_noise>0) SetDetector(8, detector_buffer, Ped_Store_fixed_noise);
      else SetDetector(8,detector_buffer,Ped_Store);
      event_number = e;
      //PedSubEvent.SetEventNumber(e);
      //PedSubEvent->SetTriggerAmount(Event_Number); //total number of events analyzed
      //cout << "Store Event " << PedSubEvent->GetTriggerAmount() << endl;
      if(dia_input==1 && Dia1.GetADC_value(84)!=detector_buffer.GetADC_value(84)) {
         cout<<"The detector_buffer has changed for event "<<e<<"; it's not longer the raw ADC data anymore :("<<endl;
         cout<<"Original data Dia1.GetADC_value(84) = "<<Dia1.GetADC_value(84)<<endl;
         cout<<"Changed data detector_buffer.GetADC_value(84) = "<<detector_buffer.GetADC_value(84)<<endl;
      }
      PedTree->Fill();
      Output_Events++;

      //Recalculate the running pedestal and running average Pedestal (CMN) for the next event

      
      //Running Pedestal for each channel
      if(dia_input == 0)
      {
         RunningPedestal(detector_buffer, Initialized_Values_Dia0, Ped_Store, Dia0buffer, Channel_Number, Iter_Size, Di_Hit_Factor, prmscount, e);
      }
      if(dia_input == 1)
      {
         RunningPedestal(detector_buffer, Initialized_Values_Dia1, Ped_Store, Dia1buffer, Channel_Number, Iter_Size, Di_Hit_Factor, prmscount, e);
      }

      //if running pedestal RMS of any channel hits zero, the loops is broken and code stopped. Something's wrong.
      if(rmscount==1)
      {
         break;
      }

      

      //Calculate Running Common Mode Noise (running average Pedestal)
      if(dia_input == 0)
      {
         RunningCommonMode(Common_Mode_Dia0_deque, Common_Mode_Mean, Common_Mode_RMS, sum_buffer/divisor, Iter_Size, e, prmscount);
         //Common_Mode_Mean = Common_Mode_Dia0_Mean;
         //Common_Mode_RMS = Common_Mode_Dia0_RMS;
      }
      if(dia_input == 1)
      {
         RunningCommonMode(Common_Mode_Dia1_deque, Common_Mode_Mean, Common_Mode_RMS, sum_buffer/divisor, Iter_Size, e, prmscount);
         //Common_Mode_Mean = Common_Mode_Dia1_Mean;
         //Common_Mode_RMS = Common_Mode_Dia1_RMS;
      }


      //Calculate running peds of rest of the detectors for the next event.
      detector_buffer = New_Event.GetD0X();
      RunningPedestal(detector_buffer, Initialized_Values_D0X, Ped_Store_D0X, D0Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD0Y();
      RunningPedestal(detector_buffer, Initialized_Values_D0Y, Ped_Store_D0Y, D0Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD1X();
      RunningPedestal(detector_buffer, Initialized_Values_D1X, Ped_Store_D1X, D1Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD1Y();
      RunningPedestal(detector_buffer, Initialized_Values_D1Y, Ped_Store_D1Y, D1Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD2X();
      RunningPedestal(detector_buffer, Initialized_Values_D2X, Ped_Store_D2X, D2Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD2Y();
      RunningPedestal(detector_buffer, Initialized_Values_D2Y, Ped_Store_D2Y, D2Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD3X();
      RunningPedestal(detector_buffer, Initialized_Values_D3X, Ped_Store_D3X, D3Xbuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);
      detector_buffer = New_Event.GetD3Y();
      RunningPedestal(detector_buffer, Initialized_Values_D3Y, Ped_Store_D3Y, D3Ybuffer, Channel_Number, Iter_Size, Si_Hit_Factor, prmscount, e);




   }//End of loop over events
   //store_event->SetTriggerAmount(Event_Number);
   //Int_t num = (Int_t) (Event_Number - (PedTree->GetEntries()));
   //cout << "Event Number is " << Event_Number << " and Entries is " << PedTree->GetEntries() << endl;
   //cout << "Num is " << num << endl;
   //store_event->SetNumRemovedEvents(num);

   //PedTree->Fill();

   //cout <<"Number events removed is " << store_event->GetNumRemovedEvents() << endl;

   if(single_channel_analysis_enable){
      for(int s=0; s<NumberSingleChannelAnalysisHistos; s++)
      {
         TCanvas *SingleChannelAnalysisCanvas = new TCanvas("SingleChannelAnalysisCanvas", "SingleChannelAnalysisCanvas", 200, 400, 600, 800);
	
         SingleChannelAnalysisNoiseHistos[s]->Draw();
         SingleChannelAnalysisNoiseHistos[s]->Fit("gaus", "Q");
	
         gSystem->ProcessEvents();
		
         HistogrammSaver::SaveCanvasPNG(SingleChannelAnalysisCanvas, png_file_char, (char*)SingleChannelAnalysisNoiseHistos[s]->GetName());
         delete SingleChannelAnalysisNoiseHistos[s];
         delete SingleChannelAnalysisCanvas;
	
      }
	
      for(int s=0; s<NumberSingleChannelAnalysisHistos; s++)
      {
         TCanvas *SingleChannelAnalysisCanvas = new TCanvas("SingleChannelAnalysisCanvas", "SingleChannelAnalysisCanvas", 200, 400, 600, 800);
	
         SingleChannelAnalysisPulseHeightHistos[s]->Draw();
	
         gSystem->ProcessEvents();
		
         HistogrammSaver::SaveCanvasPNG(SingleChannelAnalysisCanvas, png_file_char, (char*)SingleChannelAnalysisPulseHeightHistos[s]->GetName());
         delete SingleChannelAnalysisPulseHeightHistos[s];
         delete SingleChannelAnalysisCanvas;
	
      }
   }
   
   
   //Initialized_Values_D0X, Ped_Store_D0X

   //Ploting Initial RMS Values for Silicon versus Channel Number
   Float_t Si_RMS_array[8][256];
   Float_t Si_Ped_array[8][256];
   Float_t Si_Ped_array_max[8];
   for(Int_t det = 0; det<8; det++)
      Si_Ped_array_max[det] = 0;
   Float_t Si_channel_index[256];
   for(Int_t i=0; i<256; i++)
      Si_channel_index[i] = i;
   TPed_and_RMS* pedestal_data_pointer = 0;

   for(Int_t det = 0; det<8; det++) {
      switch(det) {
         case 0:
            pedestal_data_pointer = &Initialized_Values_D0X;
            break;
         case 1:
            pedestal_data_pointer = &Initialized_Values_D1X;
            break;
         case 2:
            pedestal_data_pointer = &Initialized_Values_D2X;
            break;
         case 3:
            pedestal_data_pointer = &Initialized_Values_D3X;
            break;
         case 4:
            pedestal_data_pointer = &Initialized_Values_D0Y;
            break;
         case 5:
            pedestal_data_pointer = &Initialized_Values_D1Y;
            break;
         case 6:
            pedestal_data_pointer = &Initialized_Values_D2Y;
            break;
         case 7:
            pedestal_data_pointer = &Initialized_Values_D3Y;
            break;
      }
      
      for(Int_t i=0; i<256; i++) {
         Si_Ped_array[det][i] = pedestal_data_pointer->GetPedValues(i);
         if(Si_Ped_array[det][i]>Si_Ped_array_max[det])
            Si_Ped_array_max[det] = Si_Ped_array[det][i];
         Si_RMS_array[det][i] = pedestal_data_pointer->GetRMSValues(i);
      }
   }
   
   

   //Ploting Initial RMS Values for Diamond versus Channel Number
   Float_t RMS_array[128];
   Float_t RMS_array_final[128];
   Float_t Ped_array[128];
   Float_t Ped_array_max = 0;
   Float_t channel_index[128];
   for(Int_t i=0; i<128; i++)
   {
      if(dia_input == 0)
      {
         Ped_array[i] = Initialized_Values_Dia0.GetPedValues(i+DIA_OFFSET);
      }
      if(dia_input == 1)
      {
         Ped_array[i] = Initialized_Values_Dia1.GetPedValues(i+DIA_OFFSET);
      }
      if(Ped_array[i]>Ped_array_max)
      {
         Ped_array_max = Ped_array[i];
      }
      if(screen.CheckChannel(i)==1)
      {
         if(dia_input == 0)
         {
            RMS_array[i] = Initialized_Values_Dia0.GetRMSValues(i+DIA_OFFSET);
         }
         if(dia_input == 1)
         {
            RMS_array[i] = Initialized_Values_Dia1.GetRMSValues(i+DIA_OFFSET);
         }
         if(fix_dia_noise>0) RMS_array_final[i] = Ped_Store_fixed_noise->GetRMSValues(i+DIA_OFFSET);
         else RMS_array_final[i] = Ped_Store->GetRMSValues(i+DIA_OFFSET);
      }
      else
      {
         RMS_array[i] = 0;
         RMS_array_final[i] = 0;
      }
      channel_index[i] = i; 
   }

   cout << endl << "Printing out final channel noise:" <<endl;
   for (int ch=0; ch<128; ch++) cout <<ch<<"\t"<<Ped_Store->GetRMSValues(ch+DIA_OFFSET);
   cout << endl;

   TCanvas *channel_noise_can = new TCanvas("channel_noise_can","Channel Noise",200,400,800,600);
   channel_noise_can->cd(1);
   TGraph *channel_noise = new TGraph(128, channel_index, RMS_array);
   channel_noise->Draw("AB");
   channel_noise->GetXaxis()->SetTitle("Channel");
   channel_noise->GetYaxis()->SetTitle("Noise in ADC Counts");
   channel_noise->SetTitle("Noise Per Channel");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(channel_noise_can, png_file_char, "Channel_Noise_Initial");
      HistogrammSaver::SaveCanvasC(channel_noise_can, C_file_char, "Channel_Noise_Initial");
      HistogrammSaver::SaveCanvasRoot(channel_noise_can, root_file_char, "Channel_Noise_Initial");
      if(ClosePlotsOnSave == 1)
      {
         delete channel_noise_can;
      }
   }

   TCanvas *channel_noise_fin_can = new TCanvas("channel_noise_fin_can","Final Channel Noise",200,400,800,600);
   channel_noise_fin_can->cd(1);
   TGraph *channel_noise_fin = new TGraph(128, channel_index, RMS_array_final);
   channel_noise_fin->Draw("AB");
   channel_noise_fin->GetXaxis()->SetTitle("Channel");
   channel_noise_fin->GetYaxis()->SetTitle("Noise in ADC Counts");
   channel_noise_fin->SetTitle("Noise Per Channel");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(channel_noise_fin_can, png_file_char, "Channel_Noise_Final");
      HistogrammSaver::SaveCanvasC(channel_noise_fin_can, C_file_char, "Channel_Noise_Final");
      HistogrammSaver::SaveCanvasRoot(channel_noise_fin_can, root_file_char, "Channel_Noise_Final");
      if(ClosePlotsOnSave == 1)
      {
         delete channel_noise_fin_can;
      }
   }

   TCanvas *channel_pedestal_values_can = new TCanvas("channel_pedestal_values_can","Pedestal Values by Channel",200,400,800,600);
   channel_pedestal_values_can->cd(1);
   TGraph *channel_pedestal_values = new TGraph(128, channel_index, Ped_array);
   channel_pedestal_values->Draw("AB");
   channel_pedestal_values->GetXaxis()->SetTitle("Channel");
   channel_pedestal_values->GetYaxis()->SetRangeUser(0,Ped_array_max+100);
   //channel_pedestal_values->GetYaxis()->SetTitle("Pedestal Value");
   channel_pedestal_values->SetTitle("Pedestal Value by Channel");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(channel_pedestal_values_can, png_file_char, "Pedestal_Values");
      HistogrammSaver::SaveCanvasC(channel_pedestal_values_can, C_file_char, "Pedestal_Values");
      HistogrammSaver::SaveCanvasRoot(channel_pedestal_values_can, root_file_char, "Pedestal_Values");
      if(ClosePlotsOnSave == 1)
      {
         delete channel_pedestal_values_can;
      }
   }
   
   TGraph *si_channel_pedestals = 0;
   string det_name, ped_title;
   for(Int_t det = 0; det<8; det++) {
      switch(det) {
         case 0:
            det_name = "D0X";
            break;
         case 1:
            det_name = "D1X";
            break;
         case 2:
            det_name = "D2X";
            break;
         case 3:
            det_name = "D3X";
            break;
         case 4:
            det_name = "D0Y";
            break;
         case 5:
            det_name = "D1Y";
            break;
         case 6:
            det_name = "D2Y";
            break;
         case 7:
            det_name = "D3Y";
            break;
      }
      
      ped_title = "Pedestal_mean_vs_channel_" + det_name;
      channel_pedestal_values_can = new TCanvas("channel_pedestal_values_can","Pedestal Values by Channel",200,400,800,600);
      channel_pedestal_values_can->cd(1);
      si_channel_pedestals = new TGraph(256, Si_channel_index, Si_Ped_array[det]);
      si_channel_pedestals->Draw("AB");
      si_channel_pedestals->GetXaxis()->SetTitle("Channel");
      si_channel_pedestals->GetYaxis()->SetRangeUser(0,Si_Ped_array_max[det]*1.1);
      si_channel_pedestals->GetYaxis()->SetTitle("Pedestal Mean");
      si_channel_pedestals->SetTitle(ped_title.c_str());
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(channel_pedestal_values_can, png_file_char, (char*) ped_title.c_str());
         SaveCanvasC(channel_pedestal_values_can, C_file_char, (char*) ped_title.c_str());
         SaveCanvasRoot(channel_pedestal_values_can, root_file_char, (char*) ped_title.c_str());
         if(ClosePlotsOnSave == 1)
         {
            delete channel_pedestal_values_can;
         }
      }
      delete si_channel_pedestals;
      
      ped_title = "Pedestal_rms_vs_channel_" + det_name;
      channel_noise_can = new TCanvas("channel_noise_can","Channel Noise",200,400,800,600);
      channel_noise_can->cd(1);
      si_channel_pedestals = new TGraph(256, Si_channel_index, Si_RMS_array[det]);
      si_channel_pedestals->Draw("AB");
      si_channel_pedestals->GetXaxis()->SetTitle("Channel");
      si_channel_pedestals->GetYaxis()->SetTitle("Pedestal RMS");
      si_channel_pedestals->SetTitle(ped_title.c_str());
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         HistogrammSaver::SaveCanvasPNG(channel_noise_can, png_file_char, (char*) ped_title.c_str());
         HistogrammSaver::SaveCanvasC(channel_noise_can, C_file_char, (char*) ped_title.c_str());
         HistogrammSaver::SaveCanvasRoot(channel_noise_can, root_file_char, (char*) ped_title.c_str());
         if(ClosePlotsOnSave == 1)
         {
            delete channel_noise_can;
         }
      }
      delete si_channel_pedestals;
      
   }


   TCanvas *noise_can = new TCanvas("noise_can","Noise Canvas",200,400,800,600);
   noise_can->cd(1);
   noise->Draw();
   noise->GetXaxis()->SetTitle("Non-Hit ADC of Single Channel");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(noise_can, png_file_char, "noise_can");
      HistogrammSaver::SaveCanvasC(noise_can, C_file_char, "noise_can");
      HistogrammSaver::SaveCanvasRoot(noise_can, root_file_char, "noise_can");
      if(ClosePlotsOnSave == 1)
      {
         delete noise_can;
      }
   }

   TCanvas *CMN_noise_can = new TCanvas("CMN_noise_can","CMN Noise Canvas",200,400,800,600);
   CMN_noise_can->cd(1);
   CMN_noise->Draw();
   CMN_noise->GetXaxis()->SetTitle("Average Value of Non-Hit Channel ADC");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(CMN_noise_can, png_file_char, "CMN_noise_can");
      HistogrammSaver::SaveCanvasC(CMN_noise_can, C_file_char, "CMN_noise_can");
      HistogrammSaver::SaveCanvasRoot(CMN_noise_can, root_file_char, "CMN_noise_can");
      if(ClosePlotsOnSave == 1)
      {
         delete CMN_noise_can;
      }
   }

   TCanvas *CMN_noise_saved_can = new TCanvas("CMN_noise_saved_can","CMN Noise Saved Canvas",200,400,800,600);
   CMN_noise_saved_can->cd(1);
   CMN_noise_saved->Draw();
   CMN_noise_saved->GetXaxis()->SetTitle("Average Value of Non-Hit Channel ADC");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(CMN_noise_saved_can, png_file_char, "CMN_noise_saved_can");
      HistogrammSaver::SaveCanvasC(CMN_noise_saved_can, C_file_char, "CMN_noise_saved_can");
      HistogrammSaver::SaveCanvasRoot(CMN_noise_saved_can, root_file_char, "CMN_noise_saved_can");
      if(ClosePlotsOnSave == 1)
      {
         delete CMN_noise_saved_can;
      }
   }

   TCanvas *hit_occup_can = new TCanvas("hit_occup_can","Hit Occupancy Canvas",200,400,800,600);
   hit_occup_can->cd(1);
   hit_occup->Draw();
   hit_occup->GetXaxis()->SetTitle("Channel");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(hit_occup_can, png_file_char, "hit_occup_can_dia");
      HistogrammSaver::SaveCanvasC(hit_occup_can, C_file_char, "hit_occup_can_dia");
      HistogrammSaver::SaveCanvasRoot(hit_occup_can, root_file_char, "hit_occup_can_dia");
      if(ClosePlotsOnSave == 1)
      {
         delete hit_occup_can;
      }
   }


   TCanvas *raw_ADC_can = new TCanvas("raw_ADC_can","Raw ADC Values by Event",200,400,1400,600);
   raw_ADC_can->cd(1);
   raw_ADC_by_event->Draw("arr");
   raw_ADC_by_event->GetXaxis()->SetTitle("Event");
   raw_ADC_by_event->GetYaxis()->SetTitle("ADC Value");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(raw_ADC_can, png_file_char, "raw_ADC_by_event");
      HistogrammSaver::SaveCanvasC(raw_ADC_can, C_file_char, "raw_ADC_by_event");
      HistogrammSaver::SaveCanvasRoot(raw_ADC_can, root_file_char, "raw_ADC_by_event");
      if(ClosePlotsOnSave == 1)
      {
         delete raw_ADC_can;
      }
   }


   TCanvas *raw_ADC_graph_can = new TCanvas("raw_ADC_graph_can","Raw ADC Values by Event Graph",200,400,1400,600);
   raw_ADC_graph_can->cd(1);
   raw_ADC_by_event_graph->Draw("AL");
   raw_ADC_by_event_graph->GetXaxis()->SetTitle("Event");
   raw_ADC_by_event_graph->GetYaxis()->SetTitle("ADC Value");
   raw_ADC_by_event_graph->GetYaxis()->SetRangeUser(Common_Mode_Dia1_Mean-100,Common_Mode_Dia1_Mean+100);
   pt->Draw();
   RMS_threshold_by_event_graph_up->SetLineColor(kRed);
   RMS_threshold_by_event_graph_up->Draw("sameL");
   RMS_threshold_by_event_graph_down->SetLineColor(kRed);
   RMS_threshold_by_event_graph_down->Draw("sameL");
   if(SaveAllFilesSwitch == 1)
   {

      SaveCanvasC(raw_ADC_graph_can, C_file_char, "Raw_ADC_vs_Event");
      SaveCanvasRoot(raw_ADC_graph_can, root_file_char, "Raw_ADC_vs_Event");
   }
   raw_ADC_by_event_graph->GetXaxis()->SetRangeUser(Initial_Event+500,Initial_Event+10500);
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      SaveCanvasPNG(raw_ADC_graph_can, png_file_char, "Raw_ADC_vs_Event_zoom");
      if(ClosePlotsOnSave == 1)
      {
         delete raw_ADC_graph_can;
      }
   }

   TCanvas *raw_ADC_graph_CMN_cut_can = new TCanvas("raw_ADC_graph_CMN_cut_can","Raw ADC Values by Event Graph with CMN Cut",200,400,1400,600);
   raw_ADC_graph_CMN_cut_can->cd(1);
   raw_ADC_by_event_CMN_cut_graph->Draw("AL");
   raw_ADC_by_event_CMN_cut_graph->GetXaxis()->SetTitle("Event");
   raw_ADC_by_event_CMN_cut_graph->GetYaxis()->SetTitle("ADC Value");
   raw_ADC_by_event_CMN_cut_graph->GetYaxis()->SetRangeUser(Common_Mode_Dia1_Mean-100,Common_Mode_Dia1_Mean+100);
   pt->Draw();
   RMS_threshold_by_event_graph_up->SetLineColor(kRed);
   RMS_threshold_by_event_graph_up->Draw("sameL");
   RMS_threshold_by_event_graph_down->SetLineColor(kRed);
   RMS_threshold_by_event_graph_down->Draw("sameL");
   if(SaveAllFilesSwitch == 1)
   {
      SaveCanvasC(raw_ADC_graph_CMN_cut_can, C_file_char, "Raw_ADC_vs_Event_CMN_cut");
      SaveCanvasRoot(raw_ADC_graph_CMN_cut_can, root_file_char, "Raw_ADC_vs_Event_CMN_cut");
   }
   raw_ADC_by_event_CMN_cut_graph->GetXaxis()->SetRangeUser(Initial_Event+500,Initial_Event+10500);
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(raw_ADC_graph_CMN_cut_can, png_file_char, "Raw_ADC_vs_Event_CMN_cut_zoom");
      if(ClosePlotsOnSave == 1)
      {
         delete raw_ADC_graph_CMN_cut_can;
      }
   }

   TCanvas *PS_ADC_graph_can = new TCanvas("PS_ADC_graph_can","Ped Subtracted ADC Values by Event Graph",200,400,1400,600);
   PS_ADC_graph_can->cd(1);
   PS_ADC_by_event_graph->Draw("AL");
   PS_ADC_by_event_graph->GetXaxis()->SetTitle("Event");
   PS_ADC_by_event_graph->GetYaxis()->SetRangeUser(-100,100);
   PS_ADC_by_event_graph->GetYaxis()->SetTitle("ADC Value");
   pt->Draw();
   PS_RMS_threshold_by_event_graph_up->SetLineColor(kRed);
   PS_RMS_threshold_by_event_graph_up->Draw("sameL");
   PS_RMS_threshold_by_event_graph_down->SetLineColor(kRed);
   PS_RMS_threshold_by_event_graph_down->Draw("sameL");
   if(SaveAllFilesSwitch == 1)
   {
      SaveCanvasC(PS_ADC_graph_can, C_file_char, "PS_ADC_vs_Event");
      SaveCanvasRoot(PS_ADC_graph_can, root_file_char, "PS_ADC_vs_Event");
   }
   PS_ADC_by_event_graph->GetXaxis()->SetRangeUser(Initial_Event+500,Initial_Event+10500);
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(PS_ADC_graph_can, png_file_char, "PS_ADC_vs_Event_zoom");
      if(ClosePlotsOnSave == 1)
      {
         delete PS_ADC_graph_can;
      }
   }

   TCanvas *PS_ADC_graph_CMNcut_can = new TCanvas("PS_ADC_graph_CMNcut_can","Ped Subtracted ADC Values by Event With CMN Cut Graph",200,400,1400,600);
   PS_ADC_graph_CMNcut_can->cd(1);
   PS_ADC_by_event_CMN_cut_graph->Draw("AL");
   PS_ADC_by_event_CMN_cut_graph->GetXaxis()->SetTitle("Event");
   PS_ADC_by_event_CMN_cut_graph->GetYaxis()->SetRangeUser(-100,100);
   PS_ADC_by_event_CMN_cut_graph->GetYaxis()->SetTitle("ADC Value");
   pt->Draw();
   PS_RMS_threshold_by_event_graph_up->SetLineColor(kRed);
   PS_RMS_threshold_by_event_graph_up->Draw("sameL");
   PS_RMS_threshold_by_event_graph_down->SetLineColor(kRed);
   PS_RMS_threshold_by_event_graph_down->Draw("sameL");
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(PS_ADC_graph_CMNcut_can, png_file_char, "PS_ADC_graph_CMNcut_can");
      HistogrammSaver::SaveCanvasC(PS_ADC_graph_CMNcut_can, C_file_char, "PS_ADC_graph_CMNcut_can");
      HistogrammSaver::SaveCanvasRoot(PS_ADC_graph_CMNcut_can, root_file_char, "PS_ADC_graph_CMNcut_can");
      if(ClosePlotsOnSave == 1)
      {
         delete PS_ADC_graph_CMNcut_can;
      }
   }

   TCanvas *corr_dist_can = new TCanvas("corr_dist_can","Corrected Distribution Canvas",200,400,800,600);
   corr_dist_can->cd(1);
   corr_dist->Draw();
   corr_dist->GetXaxis()->SetTitle("ADC of Correction");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(corr_dist_can, png_file_char, "corr_dist_can");
      HistogrammSaver::SaveCanvasC(corr_dist_can, C_file_char, "corr_dist_can");
      HistogrammSaver::SaveCanvasRoot(corr_dist_can, root_file_char, "corr_dist_can");
      if(ClosePlotsOnSave == 1)
      {
         delete corr_dist_can;
      }
   }

   TCanvas *CMN_RMS_dist_can = new TCanvas("CMN_RMS_dist_can","CMN RMS Distribution Canvas",200,400,800,600);
   CMN_RMS_dist_can->cd(1);
   CMN_RMS_dist->Draw();
   CMN_RMS_dist->GetXaxis()->SetTitle("CMN RMS in ADC");
   pt->Draw();
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      HistogrammSaver::SaveCanvasPNG(CMN_RMS_dist_can, png_file_char, "CMN_RMS_dist_can");
      HistogrammSaver::SaveCanvasC(CMN_RMS_dist_can, C_file_char, "CMN_RMS_dist_can");
      HistogrammSaver::SaveCanvasRoot(CMN_RMS_dist_can, root_file_char, "CMN_RMS_dist_can");
      if(ClosePlotsOnSave == 1)
      {
         delete CMN_RMS_dist_can;
      }
   }

   PedTree->Write();
   PedFile->Write();
   PedTree->Print();
   cout << "Number of Events in PedTree: " << PedTree->GetEntries() << endl;
   PedFile->Close();

   cout << "Number of Events Stored: " << Output_Events << endl;


   //Event Stats
   cout << "Triggered Events: " << Triggered_Events << endl;
   cout << "Bad Events: " << Bad_Events << endl;
   cout << "Events flagged with ZeroDivisorEvent_flag: " << ZeroDivisor_Events << endl;
   cout << "Events flagged with CMNEvent_flag: " << CMN_Events << endl;
   cout << "Events Sent into Corrector: " << Events_Sent_Into_Corrector << endl;
   cout << "Corrected Events: " << Corrected_Events << endl;
   cout << "Post Correction Bad Events: " << PostCorr_Bad_Events << endl;
   cout << "Post Correct CMN cut events: " << PostCorr_Events_Outside_CMN_win << endl;
   cout << "Final Number of Events Saved: " << Events_Saved << endl;

   runstats_file << "Triggered Events: " << Triggered_Events << endl;
   runstats_file << "Bad Events: " << Bad_Events << endl;
   runstats_file << "Events flagged with ZeroDivisorEvent_flag: " << ZeroDivisor_Events << endl;
   runstats_file << "Events flagged with CMNEvent_flag: " << CMN_Events << endl;
   runstats_file << "Events Sent into Corrector: " << Events_Sent_Into_Corrector << endl;
   runstats_file << "Corrected Events: " << Corrected_Events << endl;
   runstats_file << "Post Correction Bad Events: " << PostCorr_Bad_Events << endl;
   runstats_file << "Post Correct CMN cut events: " << PostCorr_Events_Outside_CMN_win << endl;
   runstats_file << "Final Number of Events Saved: " << Events_Saved << endl;




   
   /*
   if(hit_occupancy == 1)
   {
      TFile *pedestalfile = new TFile("PedestalSubtracted.root");
      TTree *pedtree = (TTree*)pedestalfile->Get("PedTree");
      if (!pedtree)
      {
         cerr << "Tree not found!" << endl;
         return;
      }
      TStored_Event_Ped_Sub ReadIn_Event;
      TStored_Event_Ped_Sub* pReadIn_Event = &ReadIn_Event;
       
      cout << "Number of Events in Still in pedtree: " << pedtree->GetEntries() << endl;

      pedtree->SetBranchAddress("EventBranch",&pReadIn_Event);
      TPSDetector_Data iD0X, iD0Y, iD1X, iD1Y, iD2X, iD2Y, iD3X, iD3Y, iDia0, iDia1, iDia;
      TCanvas *chan_freq_D0X = new TCanvas("chan_freq_D0X","D0X Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D0X->Divide(1,0);
      TCanvas *chan_freq_D0Y = new TCanvas("chan_freq_D0Y","D0Y Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D0Y->Divide(1,0);
      TCanvas *chan_freq_D1X = new TCanvas("chan_freq_D1X","D1X Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D1X->Divide(1,0);
      TCanvas *chan_freq_D1Y = new TCanvas("chan_freq_D1Y","D1Y Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D1Y->Divide(1,0);
      TCanvas *chan_freq_D2X = new TCanvas("chan_freq_D2X","D2X Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D2X->Divide(1,0);
      TCanvas *chan_freq_D2Y = new TCanvas("chan_freq_D2Y","D2Y Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D2Y->Divide(1,0);
      TCanvas *chan_freq_D3X = new TCanvas("chan_freq_D3X","D3X Number of Hits Per Channel", 200,400,800,600);
      chan_freq_D3X->Divide(1,0);
      TCanvas *chan_freq_D3Y = new TCanvas("chan_freq_D3Y","D3YNumber of Hits Per Channel", 200,400,800,600);
      chan_freq_D3Y->Divide(1,0);
      TCanvas *chan_freq_Dia0 = new TCanvas("chan_freq_Dia0","Diamond_Input_0 Number of Hits Per Channel", 200,400,800,600);
      chan_freq_Dia0->Divide(1,0);
      TCanvas *chan_freq_Dia1 = new TCanvas("chan_freq_Dia1","Diamond_Input_1 Number of Hits Per Channel", 200,400,800,600);
      chan_freq_Dia1->Divide(1,0);//

      //Histograms for Hit Occupancy / Channel Frequency
      TH1F *Dia0_Channel_Frequency = new TH1F("Dia0ChanFreq","Dia0 Number of Hits per Channel",256,-0.5,255.5);
      TH1F *Dia1_Channel_Frequency = new TH1F("Dia1ChanFreq","Dia1 Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D0X_Channel_Frequency = new TH1F("D0XChanFreq","D0X Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D0Y_Channel_Frequency = new TH1F("D0YChanFreq","D0Y Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D1X_Channel_Frequency = new TH1F("D1XChanFreq","D1X Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D1Y_Channel_Frequency = new TH1F("D1YChanFreq","D1Y Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D2X_Channel_Frequency = new TH1F("D2XChanFreq","D2X Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D2Y_Channel_Frequency = new TH1F("D2YChanFreq","D2Y Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D3X_Channel_Frequency = new TH1F("D3XChanFreq","D3X Number of Hits per Channel",256,-0.5,255.5);
      TH1F *D3Y_Channel_Frequency = new TH1F("D3YChanFreq","D3Y Number of Hits per Channel",256,-0.5,255.5);
      Dia0_Channel_Frequency->SetDirectory(0);
      Dia1_Channel_Frequency->SetDirectory(0);
      D0X_Channel_Frequency->SetDirectory(0);
      D0Y_Channel_Frequency->SetDirectory(0);
      D1X_Channel_Frequency->SetDirectory(0);
      D1Y_Channel_Frequency->SetDirectory(0);
      D2X_Channel_Frequency->SetDirectory(0);
      D2Y_Channel_Frequency->SetDirectory(0);
      D3X_Channel_Frequency->SetDirectory(0);
      D3Y_Channel_Frequency->SetDirectory(0);

      Int_t badchannel = 0;         //skip filling bad channels
      Int_t const FACTOR = 5;
      Int_t amount = (Int_t) pedtree->GetEntries();
      for(Int_t i=0; i<amount; i++)
      {
         pedtree->GetEntry(i);
         if(i==9637)
         {
            cout << "Number of events removed is " << ReadIn_Event.GetNumRemovedEvents() << endl;
         }

         iD0X = ReadIn_Event.GetPSD0X();
         iD0Y = ReadIn_Event.GetPSD0Y();
         iD1X = ReadIn_Event.GetPSD1X();
         iD1Y = ReadIn_Event.GetPSD1Y();
         iD2X = ReadIn_Event.GetPSD2X();
         iD2Y = ReadIn_Event.GetPSD2Y();
         iD3X = ReadIn_Event.GetPSD3X();
         iD3Y = ReadIn_Event.GetPSD3Y();
         iDia0 = ReadIn_Event.GetPSDia0();
         iDia1 = ReadIn_Event.GetPSDia1();

         pedestal_file << "Entry " << i << " has Event NUmber " << ReadIn_Event.GetEventNumber() << " difference: " << ReadIn_Event.GetEventNumber()-i << endl;
         for(Int_t x=0; x<256; x++)
         {
         //skip filling bad channels
            badchannel = 0;
         for(unsigned s=0; s<Det_channel_screen_channels[0].size(); s++)
            {
               if(Det_channel_screen_channels[0][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD0X.GetPSADCValue(x) > FACTOR*iD0X.GetPSRMSValue(x))
               {
                  D0X_Channel_Frequency->Fill(x);
               }
            }
              
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[1].size(); s++)
            {
               if(Det_channel_screen_channels[1][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD0Y.GetPSADCValue(x) > FACTOR*iD0Y.GetPSRMSValue(x))
               {
                  D0Y_Channel_Frequency->Fill(x);
               }
            }
               
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[2].size(); s++)
            {
               if(Det_channel_screen_channels[2][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD1X.GetPSADCValue(x) > FACTOR*iD1X.GetPSRMSValue(x))
               {
                  D1X_Channel_Frequency->Fill(x);
               }
            }
               
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[3].size(); s++)
            {
               if(Det_channel_screen_channels[3][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD1Y.GetPSADCValue(x) > FACTOR*iD1Y.GetPSRMSValue(x))
               {
                  D1Y_Channel_Frequency->Fill(x);
               }
            }
               
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[4].size(); s++)
            {
               if(Det_channel_screen_channels[4][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD2X.GetPSADCValue(x) > FACTOR*iD2X.GetPSRMSValue(x))
               {
                  D2X_Channel_Frequency->Fill(x);
               }
            }
               
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[5].size(); s++)
            {
               if(Det_channel_screen_channels[5][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD2Y.GetPSADCValue(x) > FACTOR*iD2Y.GetPSRMSValue(x))
               {
                  D2Y_Channel_Frequency->Fill(x);
               }
            }
               
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[6].size(); s++)
            {
               if(Det_channel_screen_channels[6][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD3X.GetPSADCValue(x) > FACTOR*iD3X.GetPSRMSValue(x))
               {
                  D3X_Channel_Frequency->Fill(x);
               }
            }
           
               //skip filling bad channels
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[7].size(); s++)
            {
               if(Det_channel_screen_channels[7][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iD3Y.GetPSADCValue(x) > FACTOR*iD3Y.GetPSRMSValue(x))
               {
                  D3Y_Channel_Frequency->Fill(x);
               }
            }
         }
         for(Int_t x=0+DIA_OFFSET; x<256+DIA_OFFSET; x++)
         {
              //skip filling bad channels   /// might have problem w/ this DIA_OFFSET in 2006 data
            badchannel = 0;
            for(unsigned s=0; s<Det_channel_screen_channels[8].size(); s++)
            {
               if(Det_channel_screen_channels[8][s] == x) 
                  badchannel = 1;
            }
            if(badchannel==0)
            {
               if(iDia0.GetPSADCValue(x) > FACTOR*iDia0.GetPSRMSValue(x))
               {
                  Dia0_Channel_Frequency->Fill(x-DIA_OFFSET);
               }
               if(iDia1.GetPSADCValue(x) > FACTOR*iDia1.GetPSRMSValue(x))
               {
                  Dia1_Channel_Frequency->Fill(x-DIA_OFFSET);
               }
            }
         }

      
      }

      
      chan_freq_D0X->cd(1);
      D0X_Channel_Frequency->Draw();
      D0X_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D0X_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D0X, png_file_char, "chan_freq_D0X");
         SaveCanvasC(chan_freq_D0X, C_file_char, "chan_freq_D0X");
         SaveCanvasRoot(chan_freq_D0X, root_file_char, "chan_freq_D0X");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D0X;
         }
      }

      chan_freq_D0Y->cd(1);
      D0Y_Channel_Frequency->Draw();
      D0Y_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D0Y_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D0Y, png_file_char, "chan_freq_D0Y");
         SaveCanvasC(chan_freq_D0Y, C_file_char, "chan_freq_D0Y");
         SaveCanvasRoot(chan_freq_D0Y, root_file_char, "chan_freq_D0Y");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D0Y;
         }
      }

      chan_freq_D1X->cd(1);
      D1X_Channel_Frequency->Draw();
      D1X_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D1X_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D1X, png_file_char, "chan_freq_D1X");
         SaveCanvasC(chan_freq_D1X, C_file_char, "chan_freq_D1X");
         SaveCanvasRoot(chan_freq_D1X, root_file_char, "chan_freq_D1X");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D1X;
         }
      }

      chan_freq_D1Y->cd(1);
      D1Y_Channel_Frequency->Draw();
      D1Y_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D1Y_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D1Y, png_file_char, "chan_freq_D1Y");
         SaveCanvasC(chan_freq_D1Y, C_file_char, "chan_freq_D1Y");
         SaveCanvasRoot(chan_freq_D1Y, root_file_char, "chan_freq_D1Y");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D1Y;
         }
      }

      chan_freq_D2X->cd(1);
      D2X_Channel_Frequency->Draw();
      D2X_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D2X_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D2X, png_file_char, "chan_freq_D2X");
         SaveCanvasC(chan_freq_D2X, C_file_char, "chan_freq_D2X");
         SaveCanvasRoot(chan_freq_D2X, root_file_char, "chan_freq_D2X");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D2X;
         }
      }

      chan_freq_D2Y->cd(1);
      D2Y_Channel_Frequency->Draw();
      D2Y_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D2Y_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D2Y, png_file_char, "chan_freq_D2Y");
         SaveCanvasC(chan_freq_D2Y, C_file_char, "chan_freq_D2Y"); 
       */
   /*
         SaveCanvasRoot(chan_freq_D2Y, root_file_char, "chan_freq_D2Y");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D2Y;
         }
      }

      chan_freq_D3X->cd(1);
      D3X_Channel_Frequency->Draw();
      D3X_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D3X_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D3X, png_file_char, "chan_freq_D3X");
         SaveCanvasC(chan_freq_D3X, C_file_char, "chan_freq_D3X");
         SaveCanvasRoot(chan_freq_D3X, root_file_char, "chan_freq_D3X");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D3X;
         }
      }

      chan_freq_D3Y->cd(1);
      D3Y_Channel_Frequency->Draw();
      D3Y_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      D3Y_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_D3Y, png_file_char, "chan_freq_D3Y");
         SaveCanvasC(chan_freq_D3Y, C_file_char, "chan_freq_D3Y");
         SaveCanvasRoot(chan_freq_D3Y, root_file_char, "chan_freq_D3Y");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_D3Y;
         }
      }

      chan_freq_Dia0->cd(1);
      Dia0_Channel_Frequency->Draw();
      Dia0_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      Dia0_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_Dia0, png_file_char, "chan_freq_Dia0_nocoinc");
         SaveCanvasC(chan_freq_Dia0, C_file_char, "chan_freq_Dia0_nocoinc");
         SaveCanvasRoot(chan_freq_Dia0, root_file_char, "chan_freq_Dia0_nocoinc");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_Dia0;
         }
      }

      chan_freq_Dia1->cd(1);
      Dia1_Channel_Frequency->Draw();
      Dia1_Channel_Frequency->GetXaxis()->SetTitle("Channel");
      Dia1_Channel_Frequency->SetTitle("Frequency of Hits Per Channel");
      pt->Draw();
      if(SaveAllFilesSwitch == 1)
      {
         gSystem->ProcessEvents();
         SaveCanvasPNG(chan_freq_Dia1, png_file_char, "chan_freq_Dia1_nocoinc");
         SaveCanvasC(chan_freq_Dia1, C_file_char, "chan_freq_Dia1_nocoinc");
         SaveCanvasRoot(chan_freq_Dia1, root_file_char, "chan_freq_Dia1_nocoinc");
         if(ClosePlotsOnSave == 1)
         {
            delete chan_freq_Dia1;
         }
      }

   }
   
   */

   if(IndexProduceSwitch==1 && single_channel_analysis_enable)
   {
       
       //Produce singlechannelanal.html for easy browsing
  
       //diamond page
      char index_html[] = "singlechannelanal.html";
      char loc[500];
      memcpy(loc,png_file_char.c_str(),strlen(png_file_char.c_str())+1);
      strcat(loc,index_html);
      ofstream index_file(loc);
      cout << index_html << " created at " << loc << endl;
   
      index_file << "<HTML>" << endl; 
      index_file << "<Title>Diamond Analysis Results Package Logbook</title>" << endl;
      index_file << "<body bgcolor=\"LightGreen\">" << endl;
      index_file << "<a name=\"top\"></a>" << endl;
//      index_file << "<a href=\"index.html\">Diamond</a><t><a href=\"silicon.html\">Silicon</a><t><a href=\"channelnoise.html\">SingleChannelAnalysisAnalysis</a><p>" << endl;
      index_file << "<CENTER>" << endl;
      index_file << "<Font color=\"330099\">" << endl;
      index_file << "<H1>Diamond Analysis Results Package Logbook</H1>" << endl;
      index_file << "<H3>Single Channel Analysis</H3>" << endl;
      index_file << "</font>" << endl;
      index_file << "</Center>" << endl;
      index_file << "<HR SIZE=\"10\" Color=\"#FFFF33\">" << endl;
      index_file << "<h3>" << run_number << " Series Results Package</h3>" << endl;


      index_file << "<p>***Runs Statistics*** (posted " << dateandtime.GetMonth() << "/ " << dateandtime.GetDay() << "/" << dateandtime.GetYear() << " at " << dateandtime.GetHour() << ":" << dateandtime.GetMinute() << ":" << dateandtime.GetSecond() << ")</p>" << endl;
      index_file << "<p>The following channels were selected for single channel analysis with a window of " << single_channel_analysis_eventwindow << " events: " ;
      for(unsigned i=0; i<single_channel_analysis_channels.size(); i++)
      {
         index_file << single_channel_analysis_channels[i] << ", ";
      }
      index_file << "</p>" << endl;
      index_file << "<hr size=\"5\">" << endl;
      
      
       // Plot index
      index_file << "<a name=\"plot_index\"></a>" << endl;
//      index_file << "<font color=\"f00000\"><a href=\"index.html\">For diamond plots click here.</a></font><p>" << endl;
//      index_file << "<font color=\"f00000\"><a href=\"silicon.html\">For silicon plots click here.</a></font><p>" << endl;
      index_file << "<p><h2>Plot Index</h2></p>" << endl;
      for(unsigned i=0; i<single_channel_analysis_channels.size(); i++)
      {
         index_file << "<p><a href=\"#Channel_"<< single_channel_analysis_channels[i] <<"\">Channel "<< single_channel_analysis_channels[i] <<" analysis </a></p>" << endl;
      }
      index_file << "<hr size=\"5\">" << endl;
      
      
      for(unsigned i=0; i<single_channel_analysis_channels.size(); i++)
      {
         index_file << "<a name=\"Channel_"<< single_channel_analysis_channels[i] <<"\"></a>" << endl;
         index_file << "<h1>Channel "<< single_channel_analysis_channels[i] <<" analysis</h1>" << endl;
         for(Int_t j=0; j<Int_t(Event_Number / single_channel_analysis_eventwindow); j++)
         {
            //title
            index_file << "<p>Pedestal subtracted ADC value for non-hit events (left) and for hit events (right) between events "<< j * single_channel_analysis_eventwindow + Initial_Event << " and "<< (j+1) * single_channel_analysis_eventwindow + Initial_Event <<"</p>" << endl;
            //plot
            index_file << "<img src=\"Channel_Noise_ch"<<single_channel_analysis_channels[i]<<"_window"<<j<<".png\">" << endl;
            index_file << "<img src=\"Channel_PulseHeight_ch"<<single_channel_analysis_channels[i]<<"_window"<<j<<".png\">" << endl;
         }
         index_file << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
         index_file << "<hr size=\"5\">" << endl;
      }
      
      index_file << "</HTML>" << endl;
      index_file.close();
   }

   pedestal_file.close();
   runstats_file.close();
   average_ped_file.close();
   ped_subtracted_file.close();
   
   Watch.Stop();
   Watch.Print("u");
} //End of PedestalAnalyze


void SlidingPedestal::ReadRawEvent(int EventNumber, bool verbose) { 

   //Declarations
   Int_t EventsPerFile = 10000;

   //Check that TEvent size is consistant -- TEvent is declared in Diamondstruct.h


   //Filename to lookup event
   std::ostringstream filename;
   filename << "RUN_" << run_number << "_" << EventNumber/EventsPerFile << ".rz";

   //Open the desired rz file if not open already
   if(current_rz_filename!=filename.str()) {
      current_rz_file.close(); // desired filename is different so close old file
      current_rz_filename=filename.str();
      if(verbose) {
         cout << "**************************************" << endl;
         cout << "Opening " << filename.str() << " for read" << endl;
         cout << "**************************************" << endl;
      }
      current_rz_file.open(filename.str().c_str(),ios::in | ios::binary);  //The .c_str() must be added for ifstream to be able to read in the file name string. 
      if (!current_rz_file) {
         cout << "File open error: " << filename.str() << " not found" << endl;
         return;
      }
   }
   
   //Read in event data (Header, Data, and Trailer) using Event structure as specified in Diamondstuct.h
   current_rz_file.seekg(EventNumber%EventsPerFile * sizeof(TEvent),ios::beg);
   current_rz_file.read(reinterpret_cast<char*>(&TEvent),sizeof(TEvent));

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
   uendian_swap(TEvent.Eor);

   //Reading out Event Header Data to Screen
   if(verbose) {
      cout << "Header dump:" << endl;
      cout << "EvTrig: " << TEvent.EvTrig << endl;
      cout << "EvNo: " << TEvent.EvNo << endl;
      cout << "EvPos: " << TEvent.EvPos << endl;
      cout << "EvTag: " << TEvent.EvTag << endl;
      cout << "EvDate: " << TEvent.EvTime << endl;
      cout << "TrigCnt: " << TEvent.TrigCnt << endl;
      cout << "EvVmeTime: " << TEvent.EvVmeTime << endl;
      for (int j=0; j<8; j++) cout << "VFasCnt[" << j << "]: " << TEvent.VFasCnt[j] << endl;
      for (int j=0; j<8; j++) cout << "VFasReg[" << j << "]: " << TEvent.VFasReg[j] << endl;
      cout << "EvNetTime: " << TEvent.EvNetTime << endl;
      cout << "MeasNo: " << TEvent.MeasNo << endl;
      cout << "EvInMeasNo: " << TEvent.EvInMeasNo << endl;
      cout << "Reserved[EVENT_HEADER_RESERVED_ESZ-2]: " << TEvent.Reserved[0] << " and " << TEvent.Reserved[1] << endl;
      cout << "Eor: " << TEvent.Eor<< endl;
   }

   
   //Sorting the values of intoutput into different columns for the x and y strips of the 4 detectors
   
   /* ---------------------------------------------------------------------------------------------------------------------
   As described in the original DAQ readme, the 2048 bytes for each event of the silicon telescope data is read in the 
   following order: To start, the first channel of the X layer ADC values are read in for each of the 4 detectors and 
   then so on for each channel until 256. Then the Y layer ADC values are read in a similar manner. So explicitly, channel 0 
   of the D0X is read in first, then channel 0 of D1X is next, then channel 0 of D2X, then channel 0 of D3X, then channel 1
   of D0X and so on. This explains the the way the values are then sorted into the 256 byte arrays for each Detector layer.
   ------------------------------------------------------------------------------------------------------------------------ */
   
   //store the raw adc values in class member storage
   for (int i=0; i<256; i++)
   {
      //NOTE: Realized that due to scrambling of the data (example 0x3615abcd written in header as 0x1536cdab) 
      //and realizing that both silicon and diamond data is written down as 4-byte words, detectors should be 
      //swapped as follows: 0->3, 1->2, 2->1, 0->3, 4->5, 5->4; the last two say that dia0->dia1 and dia1->dia0

      /*
      //Old way of mapping detectors
      D0X.ADC_values[255-i]=TEvent.Input[4*i]; //The 0X and 2X detectors are actually physically flipped so we need to reverse the values with the program
      D1X.ADC_values[i]=TEvent.Input[4*i+1];
      D2X.ADC_values[255-i]=TEvent.Input[4*i+2];
      D3X.ADC_values[i]=TEvent.Input[4*i+3];
      
      D0Y.ADC_values[i]=TEvent.Input[4*i+1024];
      D1Y.ADC_values[i]=TEvent.Input[4*i+1+1024];
      D2Y.ADC_values[i]=TEvent.Input[4*i+2+1024];
      D3Y.ADC_values[i]=TEvent.Input[4*i+3+1024];
      
      Dia0.ADC_values[i]=TEvent.RD42[i*2];
      Dia1.ADC_values[i]=TEvent.RD42[i*2+1];
      */
      
      //New way of mapping detectors
      D0X.ADC_values[255-i]=TEvent.Input[4*i+3]; //The 0X and 2X detectors are actually physically flipped so we need to reverse the values with the program
      D1X.ADC_values[i]=TEvent.Input[4*i+2];
      D2X.ADC_values[255-i]=TEvent.Input[4*i+1];
      D3X.ADC_values[i]=TEvent.Input[4*i];
      
      D0Y.ADC_values[i]=TEvent.Input[4*i+3+1024];
      D1Y.ADC_values[i]=TEvent.Input[4*i+2+1024];
      D2Y.ADC_values[i]=TEvent.Input[4*i+1+1024];
      D3Y.ADC_values[i]=TEvent.Input[4*i+1024];
      
      Dia0.ADC_values[i]=TEvent.RD42[i*2+1];
      Dia1.ADC_values[i]=TEvent.RD42[i*2];
   }
         
   //Memory Consistancy check. If The amounof Diamond memory is not 263 bytes worth, then this will output a non sensical number.
   if (EventNumber%500 == 0)
   {
      cout << "For requested event " << EventNumber << ", the current value of EvTrig and EvPos is: " << TEvent.EvTrig << " and " << TEvent.EvPos << endl;
   }
   
   return;
   
} //end of binary data read-in while loop


