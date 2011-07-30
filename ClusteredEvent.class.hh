/*
 * ClusteredEvent.class.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef CLUSTEREDEVENT_CLASS_HH_
#define CLUSTEREDEVENT_CLASS_HH_

#include <vector>
#include <iostream>

#include "Cluster.class.hh"

typedef unsigned int uint;
using namespace std;

class ClusteredEvent {
   public:
      ClusteredEvent(int EventNumber = -1, float SiHitFactor = 3, float SiSeedFactor = 5, float DiHitFactor = 3, float DiSeedFactor = 5);
      ~ClusteredEvent();
      Cluster* AddCluster(int det, bool cut_cluster = 0);
      Cluster* GetCluster(int det, int index, bool cut_cluster = 0);
      Cluster* GetCutCluster(int det, int index);
      unsigned int GetNClusters(int det, bool cut_cluster = 0);
      unsigned int GetNCutClusters(int det);
      unsigned int GetNGoodClusters(int det);
      bool HasGoldenGateCluster(int det = -1 /*det=-1 searches all detectors*/);
      bool HasBadChannelCluster(int det = -1 /*det=-1 searches all silicon detectors*/);
      bool HasLumpyCluster(int det = -1 /*det=-1 searches all silicon detectors*/);
      bool HasSaturatedCluster(int det = -1 /*det=-1 searches all silicon detectors*/);
      bool HasOneSiliconTrack();
      bool HasCoincidentClustersXY(int det /*det=0,1,2,3*/);
      void SetEventNumber(int EventNumber);
      void SetThresholds(float SiHitFactor = 3, float SiSeedFactor = 5, float DiHitFactor = 3, float DiSeedFactor = 5);
      void Clear();
      void Print();

   public:
      vector<Cluster> detector_clusters[9]; //array (index for each detector) of vectors of Clusters
      vector<Cluster> detector_clusters_cut[9]; //array (index for each detector) of vectors of Clusters
      int event_number;
      Float_t Si_Cluster_Hit_Factor;
      Float_t Si_Cluster_Seed_Factor;
      Float_t Di_Cluster_Hit_Factor;
      Float_t Di_Cluster_Seed_Factor;
};

#endif /* CLUSTEREDEVENT_CLASS_HH_ */
