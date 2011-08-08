/*
 * Clustering.class.hh
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 */

#ifndef CLUSTERING_CLASS_HH_
#define CLUSTERING_CLASS_HH_


//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


//ROOT Class Headers
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h" // for adding your own classes to ROOT's library
#include "TStyle.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "FidCutRegion.hh"


using namespace TMath;
#include "Cluster.class.hh"
#include "ClusteredEvent.class.hh"
#include "ChannelScreen.hh"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "AlignmentClass.hh"
#include "TSettings.class.hh"
typedef unsigned int uint;

class Clustering {
   public:
      Clustering(unsigned int RunNumber, std::string RunDescription = ""); //open files
      ~Clustering(); //close files
      void ClusterEvent(bool verbose = 0);
	void ClusterEventSeeds(bool verbose = 0);
      void BookHistograms();
      void InitializeHistograms();
      void DeleteHistograms();
      void DrawHistograms();
      void GenerateHTML();
      void ClusterRun(bool plots = 1);

      void Alignment(bool plots, bool CutFakeTracksOn);
	  void HistCleaner(int regid, TH2F* histo);
      void AutoFidCut();
//	void LinTrackFit(std::vector<Float_t> x_positions, std::vector<Float_t> y_positions, std::vector<Float_t> &par);
	void EventMonitor(int CurrentEvent = 0);
	void SetRunParameters(int reg, FidCutRegion current_region, bool MultiRegions = false);
   public://See if needded
	bool getUseAutoFidCut(){return settings->getUseAutoFidCut();};
	void setAlternativeClustering(bool value){settings->setAlternativeClustering(value);}
	void setUseAutoFidCut(bool value){settings->setUseAutoFidCut(value);}
//
//	bool UseAutoFidCut;
//	bool AlternativeClustering;

   private:
	TSettings *settings;
    AlignmentClass* alignment;


      std::vector<int> single_channel_analysis_channels;


      //Alignment
      std::vector<TDiamondTrack> tracks, tracks_fidcut;
      std::vector<bool> tracks_mask, tracks_fidcut_mask;
	std::vector<Float_t> dia_offset;

      //Telescope geometry
      Double_t detectorD0Z;
      Double_t detectorD1Z;
      Double_t detectorD2Z;
      Double_t detectorD3Z;
      Double_t detectorDiaZ;

      //Is the diamond aligned to Silicon x coordinates?


      //How should charge interpolation be done for two hit clusters?



      //paths
      std::string plots_path;
      std::string png_file_char;
      std::string C_file_char;
      std::string root_file_char;
      TSystem* sys;
      std::string settings_file;
	std::string pedfile_path;

	TADCEventReader *eventReader;
      //processed event storage; read from pedtree
      //UInt_t run_number;

      //output clusters
      ClusteredEvent clustered_event;

      //output tracks
      Cluster track_clusters[9];
      Float_t track_fit_xint, track_fit_xint_err, track_fit_yint, track_fit_yint_err, track_fit_xslope, track_fit_xslope_err, track_fit_yslope, track_fit_yslope_err, track_fit_rmsdev;

      //io
      UInt_t current_event;
      TDatime dateandtime;
      TRandom3 rand;
      HistogrammSaver *histSaver;


   public:

      //plots
      TPaveText *pt; //plot tag attached to all plots
      TH1F* histo_dianoise[2]; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
    TH1F* histo_detnoise[8];
      TH1F* histo_clusterocc[9][11][3]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
      TH1F* histo_landau[9][11][3][2]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
      TH1F* histo_clusterfreq[9][2]; //second index: nclusters, ngoodclusters
      TH1F* histo_clustersizefreq[9][3][3]; //second index: nhits, nseeds, standard deviation; third index: no fidcut, fidcut, no track req
      TH2F* histo_scatter[5][3]; //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut,
      TH1F* histo_eta[9][4+2][3][2]; //second index: transparent,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req; fourth index: first readout chip, 2nd readout chip
      TH1F* histo_etaintegral[9][3][2]; //second index: left to right, right to left, average; third index: first readout chip, 2nd readout chip
      TH1F* histo_trackfreq_vs_time; //show how many tracks at a given time
      TH2F* histo_eta_vs_Q[9][2][2]; //second index: transparent,2hit; third index: first readout chip, 2nd readout chip
      TH1F* histo_hitocc_saturated[9][2]; //second index: number saturated, fraction saturated
      TH1F* histo_clusters_average[9][3]; //second index: no fidcut, fidcut, no track req

	  TH2F* histo_afc_clone;
	  TH1F* histo_afc_x; //defined for AutoFidCut(), sum of hits per x bin, (max, 19.11.2010)
	  TH1F* histo_afc_x_nzv; //definded for AutoFidCut(), # of non-zero values in x bins (max, 19.11.2010)
	  TH1F* histo_afc_x_mean; //defined for AutoFidCut(), mean for x bins (max, 19.11.2010)
	  TH1F* histo_afc_y; //defined for AutoFidCut(), sum of hits per y bin (max, 19.11.2010)
	  TH1F* histo_afc_y_nzv; //defined for AutoFidCut(), # of non-zero values in y bins (max, 19.11.2010)
	  TH1F* histo_afc_y_mean; //defined for AutoFidCut(), mean for y bins (max, 19.11.2010)

      TH1F* histo_afc_unit_histo_1f;
      TH2F* histo_afc_unit_histo_2f;
      TH2F* histo_afc_region_1_mask;

	  TH1F* histo_afc_x_cut;
	  TH1F* histo_afc_y_cut;

	TH1F* histo_transparentclustering_landau[10]; // index: n channels
    TH1F* histo_transparentclustering_landau_mean;
	TH1F* histo_transparentclustering_eta;
	TH1F* histo_transparentclustering_hitdiff;
    TH2F* histo_transparentclustering_hitdiff_scatter;
    TH1F* histo_transparentclustering_2Channel_PulseHeight;

	TH2F* histo_scatter_autofidcut;

	  FidCutRegion* FCR[4];

      //added for handling of overall verbosity
      int verbosity;

   private:

      //plots toggle switches
      bool histoswitch_dianoise[2]; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
      bool histoswitch_clusterocc[9][11][3]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
      bool histoswitch_landau[9][11][3][2]; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
      bool histoswitch_clusterfreq[9][2]; //second index: nclusters, ngoodclusters
      bool histoswitch_clustersizefreq[9][3][3]; //second index: nhits, nseeds, standard deviation; third index: no fidcut, fidcut, no track req
      bool histoswitch_scatter[5][3]; //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut,
      bool histoswitch_eta[9][4+2][3][2]; //second index: 2highest,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req; fourth index: first readout chip, 2nd readout chip
      bool histoswitch_etaintegral[9][3][2]; //second index: left to right, right to left, average; third index: first readout chip, 2nd readout chip
      bool histoswitch_trackfreq_vs_time; //show how many tracks at a given time
      bool histoswitch_eta_vs_Q[9][2][2]; //second index: transparent,2hit; third index: first readout chip, 2nd readout chip
      bool histoswitch_hitocc_saturated[9][2]; //second index: number saturated, fraction saturated
      bool histoswitch_clusters_average[9][3]; //second index: no fidcut, fidcut, no track req

	bool histoswitch_transparentclustering_landau[5]; // index: n channels
	bool histoswitch_transparentclustering_eta[5]; // index: n channels

      //fits
      TF1* histofit_dianoise[2];

      //cut flow counters
      Int_t total_events, goldengatecluster_events[10], badchannelcluster_events[11], lumpycluster_events[10], saturatedcluster_events[10], totalsurviving_events, singlesitrack_events, singlesitrack_1diamondclus_events, singlesitrack_fidcut_events, singlesitrack_fidcut_1diamondclus_events, detectorxycluster_events[4], CMNEvents, ZeroDivisorEvents;
      //alignment counters
      int counter_alignment_tracks;
      int counter_alignment_fidcut_tracks;
      int counter_alignment_only_tracks;
      int counter_alignment_only_fidcut_tracks;
      int counter_alignment_tracks_zero_suppressed;

};

#endif /* CLUSTERING_CLASS_HH_ */
