/*
 * Cluster.class.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef CLUSTER_CLASS_HH_
#define CLUSTER_CLASS_HH_
//#pragma once //limit each file from being included more than once

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
      std::vector<unsigned char> channel;
      std::vector<short> adc;
      std::vector<float> pedmean;
      std::vector<float> pedwidth;
      float highest2_eta, highest2_centroid;
      bool flag_golden_gate_cluster;
      bool flag_bad_channel_cluster;
      bool flag_lumpy_cluster;
      bool flag_saturated_cluster;
};

#endif /* CLUSTER_CLASS_HH_ */
