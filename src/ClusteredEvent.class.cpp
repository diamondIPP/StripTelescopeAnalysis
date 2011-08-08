//class for organizing clusters in an event
//2010-08-04 HasBadChannelCluster() now ignores 1-hit clusters in a masked channel
//           HasOneSiliconTrack() now ignores 1-hit clusters in a masked channel
//2010-08-19 Cut bad clusters
//           HasNGoodClusters now deprecated

//#pragma once //limit each file from being included more than once
#include "ClusteredEvent.class.hh"

using namespace std;

ClusteredEvent::ClusteredEvent(int EventNumber, float SiHitFactor, float SiSeedFactor, float DiHitFactor, float DiSeedFactor) {
   event_number = EventNumber;
   Si_Cluster_Hit_Factor = SiHitFactor;
   Si_Cluster_Seed_Factor = SiSeedFactor;
   Di_Cluster_Hit_Factor = DiHitFactor;
   Di_Cluster_Seed_Factor = DiSeedFactor;
}

ClusteredEvent::~ClusteredEvent() {}

Cluster* ClusteredEvent::AddCluster(int det, bool cut_cluster) {
   if(cut_cluster) {
      if(det<8) detector_clusters_cut[det].push_back(Cluster(Si_Cluster_Hit_Factor,Si_Cluster_Seed_Factor));
      if(det==8) detector_clusters_cut[det].push_back(Cluster(Di_Cluster_Hit_Factor,Di_Cluster_Seed_Factor));
      return &detector_clusters_cut[det][detector_clusters_cut[det].size()-1];
   }
   else {
      if(det<8) detector_clusters[det].push_back(Cluster(Si_Cluster_Hit_Factor,Si_Cluster_Seed_Factor));
      if(det==8) detector_clusters[det].push_back(Cluster(Di_Cluster_Hit_Factor,Di_Cluster_Seed_Factor));
      return &detector_clusters[det][detector_clusters[det].size()-1];
   }
}

Cluster* ClusteredEvent::GetCluster(int det, int index, bool cut_cluster) {
   if(cut_cluster) {
      if((int)detector_clusters_cut[det].size()>index) return &detector_clusters_cut[det][index];
      else {
         cout<<"ClusteredEvent::GetCluster: Index "<<index<<" is out of array bounds for detector "<<det<<"; there are only "<<GetNCutClusters(det)<<" cut clusters"<<endl;
         return (Cluster*)0;
      }
   }
   else {
      if((int)detector_clusters[det].size()>index) return &detector_clusters[det][index];
      else {
         cout<<"ClusteredEvent::GetCluster: Index "<<index<<" is out of array bounds for detector "<<det<<"; there are only "<<GetNClusters(det)<<" clusters"<<endl;
         return (Cluster*)0;
      }
   }
}

Cluster* ClusteredEvent::GetCutCluster(int det, int index) {
   return GetCluster(det,index,1);
}

unsigned int ClusteredEvent::GetNClusters(int det, bool cut_cluster) {
   if(cut_cluster) return detector_clusters_cut[det].size();
   else return detector_clusters[det].size();
}

unsigned int ClusteredEvent::GetNCutClusters(int det) {
   return detector_clusters_cut[det].size();
}

//deprecated
unsigned int ClusteredEvent::GetNGoodClusters(int det) {
   unsigned int ngoodclusters = 0;
   for(uint clus=0; clus<detector_clusters[det].size(); clus++) {
      if(GetCluster(det,clus)->IsBadChannelCluster() && GetCluster(det,clus)->GetNHits()==1) continue;
      else ngoodclusters++;
   }
   return ngoodclusters;
}

bool ClusteredEvent::HasGoldenGateCluster(int detector) {
   bool goldengateflag = 0;
   for(uint cutclus=0; cutclus<2; cutclus++) {
      if(detector==-1||detector==10) {
         for(int det=0; det<9; det++) {
            for(uint clus=0; clus<GetNClusters(det,cutclus); clus++) {
               if(GetCluster(det,clus,cutclus)->IsGoldenGateCluster()) goldengateflag=1;
            }
         }
      }
      else { 
         for(uint clus=0; clus<GetNClusters(detector,cutclus); clus++) {
            if(GetCluster(detector,clus,cutclus)->IsGoldenGateCluster()) goldengateflag=1;
         }
      }
   }
   return goldengateflag;
}

bool ClusteredEvent::HasBadChannelCluster(int detector) {
   bool badchannelflag = 0;
   for(uint cutclus=0; cutclus<2; cutclus++) {
      if(detector==-1||detector==10) {
         for(int det=0; det<9; det++) {
            for(uint clus=0; clus<GetNClusters(det,cutclus); clus++){
               if(GetCluster(det,clus,cutclus)->IsBadChannelCluster() && GetCluster(det,clus,1)->GetNHits()>1) badchannelflag=1;
            }
         }
      }
      else {
         for(uint clus=0; clus<GetNClusters(detector,cutclus); clus++) {
            if(GetCluster(detector,clus,cutclus)->IsBadChannelCluster() && GetCluster(detector,clus,1)->GetNHits()>1) badchannelflag=1;
         }
      }
   }
   return badchannelflag;
}

bool ClusteredEvent::HasLumpyCluster(int detector) {
   bool lumpyflag = 0;
   for(uint cutclus=0; cutclus<2; cutclus++) {
      if(detector==-1||detector==10) {
         for(int det=0; det<9; det++) {
            for(uint clus=0; clus<GetNClusters(det,cutclus); clus++){
               if(GetCluster(det,clus,cutclus)->IsLumpyCluster()) lumpyflag=1;
            }
         }
      }
      else {
         for(uint clus=0; clus<GetNClusters(detector,cutclus); clus++) {
            if(GetCluster(detector,clus,cutclus)->IsLumpyCluster()) lumpyflag=1;
         }
      }
   }
   return lumpyflag;
}

bool ClusteredEvent::HasSaturatedCluster(int detector) {
   bool saturatedflag = 0;
   for(uint cutclus=0; cutclus<2; cutclus++) {
      if(detector==-1||detector==10) {
         for(int det=0; det<9; det++) {
            for(uint clus=0; clus<GetNClusters(det,cutclus); clus++){
               if(GetCluster(det,clus,cutclus)->IsSaturatedCluster()) saturatedflag=1;
            }
         }
      }
      else {
         for(uint clus=0; clus<GetNClusters(detector,cutclus); clus++) {
            if(GetCluster(detector,clus,cutclus)->IsSaturatedCluster()) saturatedflag=1;
         }
      }
   }
   return saturatedflag;
}

bool ClusteredEvent::HasOneSiliconTrack() {
   bool one_and_only_one = 1;
   for(int det=0; det<8; det++) 
      if(GetNClusters(det)!=1) one_and_only_one = 0;
   return one_and_only_one;
}

bool ClusteredEvent::HasCoincidentClustersXY(int det /*det=0,1,2,3*/) {
   bool one_and_only_one = 1;
   for(int d=2*det; d<2*det+2; d++) 
      if(GetNClusters(det)!=1) one_and_only_one = 0;
   return one_and_only_one;
}

void ClusteredEvent::SetEventNumber(int EventNumber) {
   event_number = EventNumber;
}

void ClusteredEvent::SetThresholds(float SiHitFactor, float SiSeedFactor, float DiHitFactor, float DiSeedFactor) {
   Si_Cluster_Hit_Factor = SiHitFactor;
   Si_Cluster_Seed_Factor = SiSeedFactor;
   Di_Cluster_Hit_Factor = DiHitFactor;
   Di_Cluster_Seed_Factor = DiSeedFactor;
}

void ClusteredEvent::Clear() {
   for(int det=0; det<9; det++) detector_clusters[det].clear();
   for(int det=0; det<9; det++) detector_clusters_cut[det].clear();
}

void ClusteredEvent::Print() {
   cout<<endl<<"--------------------"<<endl;
   cout<<"event_number="<<event_number<<endl;
   cout<<"Si_Cluster_Hit_Factor="<<Si_Cluster_Hit_Factor<<endl;
   cout<<"Si_Cluster_Seed_Factor="<<Si_Cluster_Seed_Factor<<endl;
   cout<<"Di_Cluster_Hit_Factor="<<Di_Cluster_Hit_Factor<<endl;
   cout<<"Di_Cluster_Seed_Factor="<<Di_Cluster_Seed_Factor<<endl;
   if(HasOneSiliconTrack()) cout<<"HasOneSiliconTrack"<<endl;
   for(int det=0; det<4; det++)
      cout<<"D"<<det<<" HasCoincidentClustersXY"<<endl;
   for(uint det=0; det<9; det++) {
      cout<<endl<<"Detector "<<det<<" has "<<detector_clusters[det].size()<<" clusters ("<<GetNGoodClusters(det)<<" of which are good clusters):"<<endl;
      for(uint clus=0; clus<detector_clusters[det].size(); clus++) {
         cout<<endl;
         cout<<"   Cluster "<<clus<<" has "<<detector_clusters[det][clus].GetNHits()<<" hits and ";
         cout<<detector_clusters[det][clus].GetNSeeds()<<" seeds"<<endl;
         if(detector_clusters[det][clus].IsGoldenGateCluster()) cout<<"\tFlagged as goldengate cluster!"<<endl;
         if(detector_clusters[det][clus].IsBadChannelCluster()) cout<<"\tFlagged as bad channel cluster!"<<endl;
         if(detector_clusters[det][clus].IsLumpyCluster()) cout<<"\tFlagged as lumpy cluster!"<<endl;
         if(detector_clusters[det][clus].IsSaturatedCluster()) cout<<"\tFlagged as saturated cluster!"<<endl;
         cout<<"\tTotal charge is "<<detector_clusters[det][clus].GetCharge()<<endl;
         cout<<"\t1st moment is "<<detector_clusters[det][clus].Get1stMoment()<<endl;
         cout<<"\tStandard deviation is "<<detector_clusters[det][clus].GetStandardDeviation()<<endl;
         cout<<"\tAssuming gaussian cluster: 68% charge in "<<2*detector_clusters[det][clus].GetStandardDeviation()
               <<" strips, 95% of charge in "<<4*detector_clusters[det][clus].GetStandardDeviation()
               <<" strips, and 99.7% of charge in "<<6*detector_clusters[det][clus].GetStandardDeviation()<<" strips"<<endl;
         cout<<"\tindex\tchan\tadc\tpedmean\tpedrms\tpsadc\tsnr"<<endl;
         for(int in=0; in<detector_clusters[det][clus].GetNHits(); in++) {
            cout<<"\t"<<in<<"\t"<<(int)this->GetCluster(det,clus)->GetChannel(in)<<"\t"
                  <<(int)this->GetCluster(det,clus)->GetADC(in)<<"\t"<<this->GetCluster(det,clus)->GetPedMean(in)
                  <<"\t"<<this->GetCluster(det,clus)->GetPedWidth(in)<<"\t"
                  <<this->GetCluster(det,clus)->GetADC(in)-this->GetCluster(det,clus)->GetPedMean(in)
                  <<"\t"<<(this->GetCluster(det,clus)->GetADC(in)-this->GetCluster(det,clus)->GetPedMean(in))/this->GetCluster(det,clus)->GetPedWidth(in)<<endl;
         }
      }
      cout<<endl<<"Detector "<<det<<" has "<<detector_clusters_cut[det].size()<<" CUT CLUSTERS:"<<endl;
      for(uint clus=0; clus<detector_clusters_cut[det].size(); clus++) {
         cout<<endl;
         cout<<"   Cluster "<<clus<<" has "<<detector_clusters_cut[det][clus].GetNHits()<<" hits and ";
         cout<<detector_clusters_cut[det][clus].GetNSeeds()<<" seeds"<<endl;
         if(detector_clusters_cut[det][clus].IsGoldenGateCluster()) cout<<"\tFlagged as goldengate cluster!"<<endl;
         if(detector_clusters_cut[det][clus].IsBadChannelCluster()) cout<<"\tFlagged as bad channel cluster!"<<endl;
         if(detector_clusters_cut[det][clus].IsLumpyCluster()) cout<<"\tFlagged as lumpy cluster!"<<endl;
         if(detector_clusters_cut[det][clus].IsSaturatedCluster()) cout<<"\tFlagged as saturated cluster!"<<endl;
         cout<<"\tTotal charge is "<<detector_clusters_cut[det][clus].GetCharge()<<endl;
         cout<<"\t1st moment is "<<detector_clusters_cut[det][clus].Get1stMoment()<<endl;
         cout<<"\tStandard deviation is "<<detector_clusters_cut[det][clus].GetStandardDeviation()<<endl;
         cout<<"\tAssuming gaussian cluster: 68% charge in "<<2*detector_clusters_cut[det][clus].GetStandardDeviation()
               <<" strips, 95% of charge in "<<4*detector_clusters_cut[det][clus].GetStandardDeviation()
               <<" strips, and 99.7% of charge in "<<6*detector_clusters_cut[det][clus].GetStandardDeviation()<<" strips"<<endl;
         cout<<"\tindex\tchan\tadc\tpedmean\tpedrms\tpsadc\tsnr"<<endl;
         for(int in=0; in<detector_clusters_cut[det][clus].GetNHits(); in++) {
            cout<<"\t"<<in<<"\t"<<(int)this->GetCluster(det,clus,1)->GetChannel(in)<<"\t"
                  <<(int)this->GetCluster(det,clus,1)->GetADC(in)<<"\t"<<this->GetCluster(det,clus,1)->GetPedMean(in)
                  <<"\t"<<this->GetCluster(det,clus,1)->GetPedWidth(in)<<"\t"
                  <<this->GetCluster(det,clus,1)->GetADC(in)-this->GetCluster(det,clus,1)->GetPedMean(in)
                  <<"\t"<<(this->GetCluster(det,clus,1)->GetADC(in)-this->GetCluster(det,clus,1)->GetPedMean(in))/this->GetCluster(det,clus,1)->GetPedWidth(in)<<endl;
         }
      }
   }
   cout<<"--------------------"<<endl<<endl;
}

