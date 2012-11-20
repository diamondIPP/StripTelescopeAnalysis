/*
 * ChannelScreen.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "ChannelScreen.hh"

ClassImp(ChannelScreen);
ChannelScreen::ChannelScreen()
{
   for(Int_t i=0; i<256; i++)
   {
      channel_switch[i] = 1;
   }
}
ChannelScreen::~ChannelScreen() {}

void ChannelScreen::ScreenChannels(std::vector<int> channels_to_screen) {
	std::cout<<"Screen Channels: "<<channels_to_screen.size()<<" "<<std::endl;
   for(uint i=0; i<channels_to_screen.size(); i++) {
	   ChannelOff(channels_to_screen[i]);
   }
   this->PrintScreenedChannels();
}

void ChannelScreen::ScreenRegions(std::vector<int> regions_to_screen) {
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
	   std::cout<<index<<":";
   if(index>255||index<0) {
      std::cout<<"ChannelScreen::ChannelOff: "<<index<<" is an invalid channel to screen"<<std::endl;
      return;
   }
   else channel_switch[index] = 0;
   std::cout<<channel_switch[index]<<" ";
}

Int_t ChannelScreen::CheckChannel(Int_t index)
{
   if(index>255||index<0) {
      std::cout<<"ChannelScreen::CheckChannel: "<<index<<" is an invalid channel to check"<<std::endl;
      return 1;
   }
   else return channel_switch[index];
}

bool ChannelScreen::isScreened(UInt_t channel){
	return (CheckChannel(channel)==0);
}

void ChannelScreen::PrintScreenedChannels() {
	std::cout<<std::endl;
   for(int i=0; i<256; i++){
	  // std::cout<<CheckChannel(i);
      if(channel_switch[i]==0) std::cout<<i<<", ";
   }
   std::cout<<std::endl;
}
