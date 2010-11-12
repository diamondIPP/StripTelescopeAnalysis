//Class for storing clustered data
//2010-08-24 Added Lumpy and Saturated cluster flags
//           GetTotalADC function for raw landaus

#pragma once //limit each file from being included more than once

#include <vector>
#include <iostream>
#include "TMath.h"

typedef unsigned int uint;

class Cluster {
   public:
      Cluster(Float_t HitFactor = 3, Float_t SeedFactor = 5);
      ~Cluster();
      void Clear();
      void AddHit(unsigned char Channel, short ADC, float PedMean, float PedWidth);
      unsigned char GetChannel(int index);
      short GetADC(int index);
      float GetPedMean(int index);
      float GetPedWidth(int index);
      float GetPSADC(int index);
      float GetCharge();
      float GetTotalADC();
      float GetSNR(bool conservative = 0);
      int GetPeakHit();
      void SetEta(float Eta);
      float GetEta();
      int GetNHits();
      int GetNSeeds();
      float Get1stMoment(bool highest2 = 0);
      float Get2ndMoment();
      float GetRMS();
      float GetStandardDeviation();
      float GetPositionLocal();
      float GetPositionGlobal();
      void FlagGoldenGateCluster(bool flag = 1);
      void FlagBadChannelCluster(bool flag = 1);
      void FlagLumpyCluster(bool flag = 1);
      void FlagSaturatedCluster(bool flag = 1);
      bool IsGoldenGateCluster();
      bool IsBadChannelCluster();
      bool IsLumpyCluster();
      bool IsSaturatedCluster();
      
      
   public:
      Float_t hit_threshold, seed_threshold;
      //data relevant to clustering
      /*
      unsigned char nchannels;
      unsigned char channel[10];
      short adc[10];
      float pedmean[10];
      float pedwidth[10];
      */
      vector<unsigned char> channel;
      vector<short> adc;
      vector<float> pedmean;
      vector<float> pedwidth;
      float highest2_eta, highest2_centroid;
      bool flag_golden_gate_cluster;
      bool flag_bad_channel_cluster;
      bool flag_lumpy_cluster;
      bool flag_saturated_cluster;
};

Cluster::Cluster(Float_t HitFactor, Float_t SeedFactor) {
   hit_threshold = HitFactor;
   seed_threshold = SeedFactor;
   flag_golden_gate_cluster = 0;
   highest2_eta = -1;
   highest2_centroid = -1;
}

Cluster::~Cluster() {}

void Cluster::Clear() {
   channel.clear();
   adc.clear();
   pedmean.clear();
   pedwidth.clear();
}

void Cluster::AddHit(unsigned char Channel, short ADC, float PedMean, float PedWidth) {
   //store channel
   channel.push_back(Channel);
   adc.push_back(ADC);
   pedmean.push_back(PedMean);
   pedwidth.push_back(PedWidth);
}
unsigned char Cluster::GetChannel(int index) {
   if((uint)index<channel.size()) return channel[index];
   else return 0;
}

short Cluster::GetADC(int index) {
   if((uint)index<adc.size()) return adc[index];
   else return -1;
}

float Cluster::GetPedMean(int index) {
   if((uint)index<pedmean.size()) return pedmean[index];
   else return -1;
}

float Cluster::GetPedWidth(int index) {
   if((uint)index<pedwidth.size()) return pedwidth[index];
   else return -1;
}

float Cluster::GetPSADC(int index) {
   if((uint)index<adc.size()) return adc[index]-pedmean[index];
   else return -1;
}

float Cluster::GetCharge() {
   double totalcharge = 0, pedsubadc;
   for(uint i=0; i<channel.size(); i++) {
      pedsubadc=adc[i]-pedmean[i];
      if(pedsubadc>hit_threshold*pedwidth[i]) totalcharge += pedsubadc;
   }
   return totalcharge;
}

float Cluster::GetTotalADC() {
   double totalcharge = 0, pedsubadc;
   for(uint i=0; i<channel.size(); i++) {
      pedsubadc=adc[i];
      if(adc[i]-pedmean[i]>hit_threshold*pedwidth[i]) totalcharge += pedsubadc;
   }
   return totalcharge;
}

float Cluster::GetSNR(bool conservative) {
   
   if(conservative) {
      double pedwidthavg=0;
      int nhits=GetNHits();
      for(int ch=0; ch<nhits; ch++)
         pedwidthavg += GetPedWidth(ch);
      pedwidthavg = pedwidthavg/nhits;
      return GetCharge()/pedwidthavg;
   }
   
   else {
      double totalsnr=0, snr;
      for(uint i=0; i<channel.size(); i++) {
         snr=(adc[i]-pedmean[i])/pedwidth[i];
         if(snr>hit_threshold) totalsnr += snr;
      }
      return totalsnr;
   }
}

int Cluster::GetPeakHit() {
   int clus_peak = -1;
   int clus_peak_adc = -1;
   for(int h=0; h<GetNHits(); h++) {
      if(GetADC(h)>clus_peak_adc) {
         clus_peak = h;
         clus_peak_adc = GetADC(h);
      }
   }
   
   return clus_peak;
}

void Cluster::SetEta(float Eta) {
   highest2_eta = Eta;
}

float Cluster::GetEta() {
   return highest2_eta;
}

int Cluster::GetNHits() {
   int nhits = 0;
   for(uint i=0; i<channel.size(); i++) {
      if(adc[i]-pedmean[i]>hit_threshold*pedwidth[i]) nhits++;
   }
   return nhits;
}

int Cluster::GetNSeeds() {
   int nseeds = 0;
   for(uint i=0; i<channel.size(); i++) {
      if(adc[i]-pedmean[i]>seed_threshold*pedwidth[i]) nseeds++;
   }
   return nseeds;
}

float Cluster::Get1stMoment(bool highest2) {
   if(highest2) 
      return highest2_centroid;
   else {
      double firstmoment = 0;
      double psadc = 0;
      for(uint i=0; i<channel.size(); i++) {
         psadc = adc[i]-pedmean[i];
         if(psadc>hit_threshold*pedwidth[i]) firstmoment += psadc*channel[i];
      }
      return firstmoment/GetCharge();
   }
}

float Cluster::Get2ndMoment() {
   double secondmoment = 0;
   double psadc = 0;
   for(uint i=0; i<channel.size(); i++) {
      psadc = adc[i]-pedmean[i];
      if(psadc>hit_threshold*pedwidth[i]) secondmoment += psadc*channel[i]*channel[i];
   }
   return secondmoment/GetCharge();
}

float Cluster::GetRMS() {
   double mean = Get1stMoment();
   double meansquaredev = 0;
   double psadc, deviation;
   for(uint i=0; i<channel.size(); i++) {
      psadc = adc[i]-pedmean[i];
      deviation = channel[i] - mean;
      if(psadc>hit_threshold*pedwidth[i]) meansquaredev += psadc*deviation*deviation;
   }
   return TMath::Sqrt(meansquaredev/GetCharge());
}

float Cluster::GetStandardDeviation() {
   double firstmom = Get1stMoment();
   return TMath::Sqrt(Get2ndMoment()-firstmom*firstmom);
}

float Cluster::GetPositionLocal() {
   return 0;
}

float Cluster::GetPositionGlobal() {
   return 0;
}

void Cluster::FlagGoldenGateCluster(bool flag) {
   flag_golden_gate_cluster = flag;
}

void Cluster::FlagBadChannelCluster(bool flag) {
   flag_bad_channel_cluster = flag;
}

void Cluster::FlagLumpyCluster(bool flag) {
   flag_lumpy_cluster = flag;
}

void Cluster::FlagSaturatedCluster(bool flag) {
   flag_saturated_cluster = flag;
}

bool Cluster::IsGoldenGateCluster() {
   return flag_golden_gate_cluster;
}

bool Cluster::IsBadChannelCluster() {
   return flag_bad_channel_cluster;
}

bool Cluster::IsLumpyCluster() {
   return flag_lumpy_cluster;
}

bool Cluster::IsSaturatedCluster() {
   return flag_saturated_cluster;
}
