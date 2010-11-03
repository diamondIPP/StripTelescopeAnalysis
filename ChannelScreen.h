//#ifndef ChannelScreen
//#define ChannelScreen

//Channel Screen Cluster
//Each Channel is given a 1 or a 0 depending on whether it should be counted or not

#include <iostream>
#include "TMath.h"

class ChannelScreen {

   public:
      ChannelScreen();
      ~ChannelScreen();
      void ScreenChannels(vector<int> channels_to_screen);
      void ScreenRegions(vector<int> regions_to_screen);
      void ChannelOff(Int_t index);
      void PrintScreenedChannels();
      Int_t CheckChannel(Int_t index);
      
   private:
      Int_t channel_switch[256];

};

ChannelScreen::ChannelScreen()
{
   for(Int_t i=0; i<256; i++)
   {
      channel_switch[i] = 1;
   }
}
ChannelScreen::~ChannelScreen() {}

void ChannelScreen::ScreenChannels(vector<int> channels_to_screen) {
   for(uint i=0; i<channels_to_screen.size(); i++) ChannelOff(channels_to_screen[i]);
}

void ChannelScreen::ScreenRegions(vector<int> regions_to_screen) {
   if(regions_to_screen.size()%2) {
      std::cout<<"ChannelScreen::ScreenRegions: Invalid region(s) to screen"<<std::endl;
      return;
   }

   int ch1, ch2;
   
   for(uint i=0; i<regions_to_screen.size()/2; i++) {
      ch1=regions_to_screen[2*i];
      ch2=regions_to_screen[2*i+1];
      for(int j=TMath::Min(ch1,ch2); j<TMath::Max(ch1,ch2); j++)
         ChannelOff(j);
   }
}

void ChannelScreen::ChannelOff(Int_t index)
{
   if(index>255||index<0) {
      std::cout<<"ChannelScreen::ChannelOff: "<<index<<" is an invalid channel to screen"<<std::endl;
      return;
   }
   else channel_switch[index] = 0;
}

Int_t ChannelScreen::CheckChannel(Int_t index)
{
   if(index>255||index<0) {
      std::cout<<"ChannelScreen::CheckChannel: "<<index<<" is an invalid channel to check"<<std::endl;
      return 1;
   }
   else return channel_switch[index];
}

void ChannelScreen::PrintScreenedChannels() {
   for(int i=0; i<256; i++)
      if(channel_switch[i]==0) std::cout<<i<<", ";
}


//#endif
