//#ifndef ChannelScreen
//#define ChannelScreen

//Channel Screen Cluster
//Each Channel is given a 1 or a 0 depending on whether it should be counted or not
#ifndef CHANNELSCREEN_H
#define CHANNELSCREEN_H
#include <iostream>
#include "TMath.h"
#include "TObject.h"
#include <vector>
#include "TPlaneProperties.hh"

typedef unsigned int uint;

class ChannelScreen :public TObject {

   public:
      ChannelScreen(UInt_t det=0);
      ~ChannelScreen();
      void setDetectorNumber(UInt_t det){this->det=det;}
      void ScreenChannels(std::vector<int> channels_to_screen);
      void ScreenRegions(std::vector<int> regions_to_screen);
      void ChannelOff(Int_t index);
      void PrintScreenedChannels();
      Int_t CheckChannel(Int_t index);
      bool isScreened(UInt_t channel);
      
   private:
      Int_t channel_switch[256];
      UInt_t det;

      ClassDef(ChannelScreen,1);
};
#endif/*CHANNELSCREEN_H*/


//#endif /*CHANNELSCREEN_H*/
