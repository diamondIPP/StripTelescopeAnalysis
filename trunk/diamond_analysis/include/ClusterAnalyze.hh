/*
 * ClusterAnalyze.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef CLUSTERANALYZE_HH_
#define CLUSTERANALYZE_HH_


class ClusterAnalyze {
   public:
      ClusterAnalyze(unsigned int RunNumber, string RunDescription = ""); //open files
      ~ClusterAnalyze(); //close files
      void LoadSettings();
      void ParseIntArray(string value, vector<int> &vec);
      void Cluster_Clear(int det);
      void ClearClusters();
      void BookHistograms();

   private:
      //settings
      Int_t SaveAllFilesSwitch; //1 for save files, 0 for don't
      Int_t ClosePlotsOnSave;
      Int_t IndexProduceSwitch;

      Float_t Si_Cluster_Seed_Factor;
      Float_t Si_Cluster_Hit_Factor;
      Float_t Di_Cluster_Seed_Factor;
      Float_t Di_Cluster_Hit_Factor;

      vector<int> single_channel_analysis_channels;

      //Channels to Screen
      vector<int> Det_channel_screen_input[9];
      ChannelScreen Det_channel_screen[9];

      //paths
      string png_file_char;
      string C_file_char;
      string root_file_char;
      TSystem* sys;
      string settings_file;

      //processed event storage; read from pedtree
      UInt_t run_number;
      UInt_t event_number;
      Float_t store_threshold;
      UInt_t Det_NChannels[9];
      UChar_t Det_Channels[9][256];
      UChar_t Det_ADC[8][256];
      UShort_t Dia_ADC[256];
      Float_t Det_PedMean[9][256];
      Float_t Det_PedWidth[9][256];

      //output clusters
      UInt_t Det_NClusters[9];
      //Cluster Det_Cluster[9][10]; //how many max clusters per detector should we allow?
      Int_t Det_Most_Signif_Cluster[9];
      Float_t Det_Most_Signif_Cluster_Seed_PSADC[9];

      UInt_t Det_Cluster_NChannels[9];
      UChar_t Det_Cluster_Channels[9][256];
      UChar_t Det_Cluster_ADC[8][256];
      UShort_t Dia_Cluster_ADC[256];
      Float_t Det_Cluster_PedMean[9][256];
      Float_t Det_Cluster_PedWidth[9][256];

      //io
      UInt_t current_event;
      TFile *ClusterFile;
      TTree *ClusterTree;

      //histograms
      TH1F* histo_hitocc[9]; //if there's a hit, fill the bin
      TH1F* histo_noise[9][2]; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
      TH1F* histo_trackocc[9][11][2]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut
      TH1F* histo_landau[9][11][2][2]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut; fourth index: psadc, snr
      TH1F* histo_clustersize[9][2]; //second index: no fidcut, fidcut
      TH2F* histo_scatter[6][2]; //first index: d0,d1,d2,d3,<d1&d2>,<all>; second index: no fidcut, fidcut
      TH1F* histo_eta[9][4][2]; //second index: transparent,2hit,low,high; third index: no fidcut, fidcut
      TH1F* histo_trackfreq_vs_time; //show how many tracks at a given time
};

#endif /* CLUSTERANALYZE_HH_ */
