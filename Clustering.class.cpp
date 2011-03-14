//machinery for clustering
//2010-07-11 Goal: would like to multithread by splitting a run into several windows, then on each window calc the ped, cluster, track, and at each stage combine windows to fill histograms. Can have histograms talk to each other by making trees of data to histogram or writing histograms.root and appending there. 
//2010-07-23 Import settings from Setting.ini
//2010-07-27 Clustering analysis done for the most part
//           Decided to save clusters as arrays of basic datatypes (can't save classes to trees)
//           Putting off data reduction to PedestalAnalyze and ClusterAnalyze 
//2010-08-01 Started major histogramming efforts
//2010-08-02 deleting histograms CRASHES root on second call to the constructor so no more deleting
//           draws histograms and prints cutflow stats
//2010-08-04 Fixed settings array parsing bug (first element in array always was zero)
//           Better channel screening
//2010-08-05 Added run stats to histograms
//2010-08-14 Skip analysis for CMN and ZeroDivisor events (flagged in Sliding Pedestal)
//2010-08-15 Created HTML summary page
//           Fixed a slew of problems with the HTML summary page
//           Fixed some cuts (now skip bad cluster events in silicon and clusters w/ bad seeds in diamond)
//           Various other fixes
//           TODO: figure out why D1 has no coincident clusters for runs 14109 and 14110 while it works for 14107
//2010-08-19 Created LumpyCluster flag (if cluster is not monotonically decreasing from the peak)
//           TODO: create plots for cut clusters
//           TODO: create plots for high/low etas
//           TODO: write cut cluster event numbers to csv files
//2010-09-02 Eta vs Q plots and saturated cluster removal and sat clus hit occ plots
//2010-09-18 TODO: Plan for tracking is as follows
//                    -> align adhoc without tracking
//                    -> track with residuals' widths as uncertainties
//                    -> align with tracking (this makes sure tracks are really tracks; thereby improving residuals)
//2010-11-12 Blinded the alignment resolution determination
//2010-12-01 Created AutoFidCut() to Clustering.class.cpp
//2011-03-07 XCode 4.0 Test max AND second Xcode 4 Test!

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

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
#include "FidCutRegion.h"


using namespace TMath;

//Class and stuct definitions
//#include "ChannelScreen.h" //Channel Screen Class
#include "Cluster.class.cpp"
#include "ClusteredEvent.class.cpp"
#include "AlignmentClasses.h"

typedef unsigned int uint;

class Clustering {
   public:
      Clustering(unsigned int RunNumber, string RunDescription = ""); //open files
      ~Clustering(); //close files
      void LoadSettings();
      void ParseIntArray(string value, vector<int> &vec);
      void ParseFloatArray(string value, vector<float> &vec);
      void ClusterEvent(bool verbose = 0);
	void ClusterEventSeeds(bool verbose = 0);
      void BookHistograms();
      void SaveHistogram(TH1F* histo);
      void SaveHistogram(TH2F* histo);
      void SaveHistogramPNG(TH1F* histo);
      void SaveHistogramPNG(TH2F* histo);
      void SaveHistogramROOT(TH1F* histo);
      void SaveHistogramROOT(TH2F* histo);
      void SaveHistogramPDF(TH1F* histo);
      void SaveHistogramPDF(TH2F* histo);
      void InitializeHistograms();
      void DeleteHistograms();
      void DrawHistograms();
      void GenerateHTML();
      void ClusterRun(bool plots = 1);
      void Align(bool plots = 1, bool CutFakeTracksOn = false);
	  void HistCleaner(int regid, TH2F* histo);
      void AutoFidCut();
	void TransparentClustering(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);
//	void LinTrackFit(vector<Float_t> x_positions, vector<Float_t> y_positions, vector<Float_t> &par);
	void EventMonitor(int CurrentEvent = 0);
	void SetRunParameters(int reg, FidCutRegion current_region, bool MultiRegions = false);
	
	bool UseAutoFidCut;
	bool AlternativeClustering;
      
   private:
      //general settings
      Int_t SaveAllFilesSwitch; //1 for save files, 0 for don't
      Int_t ClosePlotsOnSave;
      Int_t IndexProduceSwitch;
      
      //pedestal settings
      float fix_dia_noise; // fix_dia_noise<0 disables diamond noise-fixing
      Int_t Iter_Size; //buffer size
      Int_t Taylor_speed_throttle; //# of events to recalculate RMS the old way; set to 1 to disable
      Int_t dia_input; // 1 for 2006 and 0 for the rest
      Float_t Si_Pedestal_Hit_Factor;
      Float_t Di_Pedestal_Hit_Factor;
      Int_t CMN_cut;  //Should be less than or equal to CMN_coor_high

      //clustering settings
      Float_t Si_Cluster_Seed_Factor;
      Float_t Si_Cluster_Hit_Factor;
      Float_t Di_Cluster_Seed_Factor;
      Float_t Di_Cluster_Hit_Factor;
      Float_t si_avg_fidcut_xlow;
      Float_t si_avg_fidcut_xhigh;
      Float_t si_avg_fidcut_ylow;
      Float_t si_avg_fidcut_yhigh;
      Int_t pulse_height_num_bins;
      Float_t pulse_height_si_max;
      Float_t pulse_height_di_max;
      Float_t snr_distribution_si_max;
      Float_t snr_distribution_di_max;
            
      //Hi/low eta slices
      Float_t eta_lowq_slice_low;
      Float_t eta_lowq_slice_hi;
      Float_t eta_hiq_slice_low;
      Float_t eta_hiq_slice_hi;

      //Number of slices (<1 to disable)
      Int_t etavsq_n_landau_slices;
      
      Int_t snr_plots_enable;
      vector<int> single_channel_analysis_channels;
      
      //Channels to Screen
      vector<int> Det_channel_screen_channels[9];
      vector<int> Det_channel_screen_regions[9];
      ChannelScreen Det_channel_screen[9];
         
         
      //Alignment
      vector<TDiamondTrack> tracks, tracks_fidcut;
      vector<bool> tracks_mask, tracks_fidcut_mask;
	Float_t alignment_chi2;
	vector<Float_t> dia_offset;

      //Telescope geometry
      Double_t detectorD0Z;
      Double_t detectorD1Z;
      Double_t detectorD2Z;
      Double_t detectorD3Z;
      Double_t detectorDiaZ;

      //Is the diamond aligned to Silicon x coordinates?
      bool dia_x_aligned;

      //How should charge interpolation be done for two hit clusters?
      bool eta_correction;

      //Filter tracks not in good fiducial region w/o bad strips
      Int_t align_sil_fid_xlow;
      Int_t align_sil_fid_xhi;
      Int_t align_sil_fid_ylow;
      Int_t align_sil_fid_yhi;
      
      //Alignment constants for tracking
      vector<Float_t> alignment_x_offsets;
      vector<Float_t> alignment_y_offsets;
      vector<Float_t> alignment_phi_offsets;
      vector<Float_t> alignment_z_offsets;
      
      //fraction of tracks used to derive alignment constants
      Float_t alignment_training_track_fraction;
      
      //paths
      string plots_path;
      string png_file_char;
      string C_file_char;
      string root_file_char;
      TSystem* sys;
      string settings_file;
	string pedfile_path;
      
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
      bool CMNEvent_flag;
      bool ZeroDivisorEvent_flag;
      
      //output clusters
      ClusteredEvent clustered_event;
      
      //output tracks
      Cluster track_clusters[9];
      Float_t track_fit_xint, track_fit_xint_err, track_fit_yint, track_fit_yint_err, track_fit_xslope, track_fit_xslope_err, track_fit_yslope, track_fit_yslope_err, track_fit_rmsdev;
      
      //io
      UInt_t current_event;
      TFile *PedFile;
      TTree *PedTree;
      TDatime dateandtime;
      TRandom3 rand;
      
   public:
      
      //plots
      TPaveText *pt; //plot tag attached to all plots
      TH1F* histo_dianoise[2]; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
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
	TH1F* histo_transparentclustering_eta;
	TH1F* histo_transparentclustering_hitdiff;
	
	TH2F* histo_scatter_autofidcut;
	
	  FidCutRegion* FCR[4];
      
      //added for handling of overall verbosity
      int Verbosity;
	
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

Clustering::Clustering(unsigned int RunNumber, string RunDescription) {
   
   //default general settings
   SaveAllFilesSwitch = 1; //1 for save files, 0 for don't
   ClosePlotsOnSave = 1;
   IndexProduceSwitch = 1;
   
   //default pedestal settings
   fix_dia_noise = -1;//7.7; // fix_dia_noise<0 disables diamond noise-fixing
   store_threshold = 2;
   dia_input = 0; // 1 for 2006 and 0 for the rest
   Si_Pedestal_Hit_Factor = 5;
   Di_Pedestal_Hit_Factor = 5;
   CMN_cut = 4;  //Should be less than or equal to CMN_coor_high
   Iter_Size = 500; //buffer size
   Taylor_speed_throttle = 1000; //# of events to recalculate RMS the old way; set to 1 to disable
   
   //default clustering settings
   snr_plots_enable = 0;
   Si_Cluster_Seed_Factor = 5;
   Si_Cluster_Hit_Factor = 3;
   Di_Cluster_Seed_Factor = 5;
   Di_Cluster_Hit_Factor = 3;
   si_avg_fidcut_xlow = 90;
   si_avg_fidcut_xhigh = 165;
   si_avg_fidcut_ylow = 80;
   si_avg_fidcut_yhigh = 160;
   pulse_height_num_bins = 300;
   pulse_height_si_max = 300;
   pulse_height_di_max = 3000;
   snr_distribution_si_max = 2500;
   snr_distribution_di_max = 2500;
	UseAutoFidCut = 0;
	AlternativeClustering = 0;
   
   //Hi/low eta slices
   eta_lowq_slice_low = 600;
   eta_lowq_slice_hi = 700;
   eta_hiq_slice_low = 1200;
   eta_hiq_slice_hi = 1500;
   
   //Number of slices (<1 to disable)
   etavsq_n_landau_slices = 0;
   
   
   //Alignment

   //Telescope geometry (looks like this is the same geom for all runs)
   //Note: if alignment gives large residuals, try switching wide to compact geometry or vs
   //wide geometry; edges: 0, 2.40, 9, 18, 20.40 (si modules 2.40cm wide, x/y planes spaced 2mm, D0/D1 interspacing 9mm, dia module 1.9cm wide)
   detectorD0Z = 0.725; // by definition
   detectorD1Z = 1.625; // by definition
   detectorD2Z = 18.725; // by definition
   detectorD3Z = 19.625; // by definition
   detectorDiaZ = 10.2; // by definition
   //compact geometry; edges: 0, 2.40, 6, 12, 14.40 (si modules 2.40cm wide, x/y planes spaced 2mm, D0/D1 interspacing 9mm, dia module 1.9cm wide)
   //Double_t detectorD0Z = 0.725; // by definition
   //Double_t detectorD1Z = 1.625; // by definition
   //Double_t detectorD2Z = 12.725; // by definition
   //Double_t detectorD3Z = 13.625; // by definition
   //Double_t detectorDiaZ = 7.2; // by definition
	
	//cut tracks with chi2 > alignment_chi2
	alignment_chi2 = 9999.;

   //Is the diamond aligned to Silicon x coordinates?
   dia_x_aligned = true;

   //How should charge interpolation be done for two hit clusters?
   eta_correction = false;

   //Filter tracks not in good fiducial region w/o bad strips (defaults to silicon fiducial region)
   align_sil_fid_xlow = si_avg_fidcut_xlow;
   align_sil_fid_xhi = si_avg_fidcut_xhigh;
   align_sil_fid_ylow = si_avg_fidcut_ylow;
   align_sil_fid_yhi = si_avg_fidcut_yhigh;
   
   //fraction of tracks used to derive alignment constants
   alignment_training_track_fraction = 0.25;
   
   //default paths
   sys = gSystem;
   
   ostringstream plotspath;
   plotspath << sys->pwd() << "/plots-" << RunNumber;
   if(RunDescription=="") plotspath << "/";
   else plotspath << "-" << RunDescription << "/";
   plots_path = plotspath.str();
   png_file_char = plotspath.str();
   C_file_char = plotspath.str();
   root_file_char = plotspath.str();
   //make plots dir
   sys->mkdir(plots_path.c_str());
   
   
   ostringstream settingspath;
   settingspath << sys->pwd() << "/Settings." << RunNumber;
   if(RunDescription=="") settingspath << ".ini";
   else settingspath << "-" << RunDescription << ".ini";
   settings_file = settingspath.str();
   
   ostringstream pedfilepath;
   pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
   if(RunDescription=="") pedfilepath << ".root";
   else pedfilepath << "-" << RunDescription << ".root";
	pedfile_path = pedfilepath.str();
	
	// create file histograms.root
	ostringstream root_histo_file_path;
	root_histo_file_path << plotspath.str().c_str() << "histograms.root";
	cout << "creating " << root_histo_file_path.str().c_str() << " .." << endl;
	TFile fplots(root_histo_file_path.str().c_str(),"RECREATE");
   
   LoadSettings();
   
   //screen channels
   for(int det=0; det<9; det++) {
      Det_channel_screen[det].ScreenChannels(Det_channel_screen_channels[det]);
      Det_channel_screen[det].ScreenRegions(Det_channel_screen_regions[det]);
      cout<<"Detector "<<det<<" screened channels: ";
      Det_channel_screen[det].PrintScreenedChannels();
      cout<<endl;
   }
   
   //Get pedestal subtracted event
   
   PedFile = new TFile(pedfilepath.str().c_str());
   PedTree = (TTree*)PedFile->Get("PedTree");
   if (!PedTree)
   {
      cerr << "PedTree not found!" << endl;
   }
      
   cout << PedTree->GetEntries() << " events in PedTree " << endl;
      
   //Event Header Branches
   PedTree->SetBranchAddress("RunNumber",&run_number);
   PedTree->SetBranchAddress("EventNumber",&event_number);
   PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
   PedTree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
   PedTree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);
   
   //Telescope Data Branches
   PedTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
   PedTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
   PedTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
   PedTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
   PedTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
   PedTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
   PedTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
   PedTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
   PedTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
   PedTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
   PedTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
   PedTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
   PedTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
   PedTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
   PedTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
   PedTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
   PedTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
   PedTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
   PedTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
   PedTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
   PedTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
   PedTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
   PedTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
   PedTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
   PedTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
   PedTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
   PedTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
   PedTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
   PedTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
   PedTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
   PedTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
   PedTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
   PedTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
   PedTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
   PedTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
   PedTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
   PedTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
   PedTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
   PedTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
   PedTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
   PedTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
   PedTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
   PedTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
   PedTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
   PedTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
   
   current_event = 0;
   PedTree->GetEvent(current_event);
   cout<< "Loaded first event in PedTree: "<<event_number<<endl;
   cout<< "RunNumber is: "<<run_number<<endl;
   cout<< "StoreThreshold is: "<<store_threshold<<endl;
   
   //Determine which histos to make
   for(int cut=0; cut<2; cut++) histoswitch_dianoise[cut]=1; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
   for(int det=0; det<9; det++) for(int hits=0; hits<11; hits++) for(int cut=0; cut<3; cut++) {
      if(det==8 && (cut==0||cut==1||(cut==2&&hits==10))) histoswitch_clusterocc[det][hits][cut] = 1; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
      else if(det<8 && (cut==0||cut==1||cut==2) && (hits==0||hits==1||hits==2||hits==3||hits==4||hits==10)) histoswitch_clusterocc[det][hits][cut] = 1;
      else histoswitch_clusterocc[det][hits][cut] = 0;
   }
   for(int det=0; det<9; det++) for(int hits=0; hits<11; hits++) for(int cut=0; cut<3; cut++) for(int snr=0; snr<2; snr++) {
      if(det==8 && (cut==0||cut==1) && snr==0) histoswitch_landau[det][hits][cut][snr] = 1; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
      else if(det<8 && (cut==1||cut==2) && snr==0 && (hits==0||hits==1||hits==2||hits==3||hits==4||hits==10)) histoswitch_landau[det][hits][cut][snr] = 1;
      else histoswitch_landau[det][hits][cut][snr] = 0;
   }
   for(int det=0; det<9; det++) for(int cut=0; cut<2; cut++) histoswitch_clusterfreq[det][cut] = 1;
   for(int det=0; det<9; det++) for(int metric=0; metric<3; metric++) for(int cut=0; cut<3; cut++) {
      if(metric!=2 && det==8 && (cut==0||cut==1)) histoswitch_clustersizefreq[det][metric][cut] = 1; //second index: nhits, nseeds, standard deviation; third index: no fidcut, fidcut, no track req
      else if(metric!=2 && det<8 && cut==0) histoswitch_clustersizefreq[det][metric][cut] = 1;
      else histoswitch_clustersizefreq[det][metric][cut] = 0;
   }
   for(int det=0; det<5; det++) for(int cut=0; cut<3; cut++) histoswitch_scatter[det][cut] = 1; //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut
   for(int det=0; det<9; det++) for(int type=0; type<4; type++) for(int cut=0; cut<3; cut++) for(int chip=0; chip<2; chip++) {
      if(det==8 && (type==0||type==1) && (cut==0||cut==1)) histoswitch_eta[det][type][cut][chip] = 1; //second index: transparent,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req; fourth index: 1st readout chip, 2nd readout chip
      else if(det<8 && (type==0||type==1) && (cut==0||cut==1)) histoswitch_eta[det][type][cut][chip] = 1;
      else histoswitch_eta[det][type][cut][chip] = 0;
   }
   for(int det=0; det<9; det++) for(int type=0; type<3; type++) for(int chip=0; chip<2; chip++) {
      if(1) histoswitch_etaintegral[det][type][chip] = 1; //second index: transparent,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req
      else histoswitch_etaintegral[det][type][chip] = 0;
   }
   histoswitch_trackfreq_vs_time = 0; //show how many tracks at a given time
   for(int det=0; det<9; det++) for(int etatype=0; etatype<2; etatype++) for(int chip=0; chip<2; chip++) histoswitch_eta_vs_Q[det][etatype][chip] = 1;
   for(int det=0; det<9; det++) for(int frac=0; frac<2; frac++) {
      if(det<8 && frac==0) histoswitch_hitocc_saturated[det][frac] = 1;
      else if(det==8) histoswitch_hitocc_saturated[det][frac] = 1;
      else histoswitch_hitocc_saturated[det][frac] = 0;
   }
   for(int det=0; det<9; det++) for(int trackcut=0; trackcut<3; trackcut++) { //second index: no fidcut, fidcut, no track req
      if(trackcut==0) histoswitch_clusters_average[det][trackcut] = 1;
      else histoswitch_clusters_average[det][trackcut] = 0;
   }
   
   //initialize histograms
   InitializeHistograms();
   string histo_name, det_name, trackcut_name, nhitsopt_name, etaopt_name, clussizeopt_name, chip_name;
   
   histo_dianoise[0] = new TH1F("noise_diamond","noise_diamond",80,-40,40);
   histo_dianoise[1] = new TH1F("noise_diamond_fidcutevents","noise_diamond_fidcutevents",80,-40,40);
   
   for(int det=0; det<5; det++) {
      
      switch(det) {
      case 0: det_name = "_D0"; break;
      case 1: det_name = "_D1"; break;
      case 2: det_name = "_D2"; break;
      case 3: det_name = "_D3"; break;
      }
      
      for(int trackcut=0; trackcut<3; trackcut++) {
         if(det==4) histo_name = "Silicon8HitsScatter_Average";
         else histo_name = "Silicon8HitsScatter" + det_name;
         if(trackcut==1) histo_name += "_1DiamondCluster";
         if(trackcut==2) histo_name += "_Fidcut";
         histo_scatter[det][trackcut] = new TH2F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5,256,-0.5,255.5); //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut
      }
   }
	
	histo_scatter_autofidcut = new TH2F("AutoFidCut_Silicon8HitsScatter_Average_1DiamondCluster","AutoFidCut_Silicon8HitsScatter_Average_1DiamondCluster",256,-0.5,255.5,256,-0.5,255.5);
   
   for(int det=0; det<9; det++) {
      
      switch(det) {
         case 0: det_name = "_D0X"; break;
         case 1: det_name = "_D0Y"; break;
         case 2: det_name = "_D1X"; break;
         case 3: det_name = "_D1Y"; break;
         case 4: det_name = "_D2X"; break;
         case 5: det_name = "_D2Y"; break;
         case 6: det_name = "_D3X"; break;
         case 7: det_name = "_D3Y"; break;
         case 8: det_name = "_Dia"; break;
      }
      
      
      for(int chip=0; chip<2; chip++) {
         if(chip) chip_name = "_RightChip";
         else chip_name = "_LeftChip";
         
         histo_name = "EtaVsQ_TransparentEta_8HitsFidcut" + det_name + chip_name;
         if(det<8) histo_eta_vs_Q[det][0][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,0,pulse_height_si_max,101,-0.005,1.005);
         if(det==8) histo_eta_vs_Q[det][0][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,0,pulse_height_di_max,101,-0.005,1.005);
         
         histo_name = "EtaVsQ_2HitEta_8HitsFidcut" + det_name + chip_name;
         if(det<8) histo_eta_vs_Q[det][1][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,0,pulse_height_si_max,101,-0.005,1.005);
         if(det==8) histo_eta_vs_Q[det][1][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,0,pulse_height_di_max,101,-0.005,1.005);
         
      }
      for(int frac=0; frac<2; frac++) {
         histo_name = "HitOccupancy_SaturatedChannels";
         if(frac==1) histo_name += "_RelativeFrequency";
         histo_name += det_name;
         if(det<8) histo_hitocc_saturated[det][frac] = new TH1F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5);
         if(det==8) histo_hitocc_saturated[det][frac] = new TH1F(histo_name.c_str(),histo_name.c_str(),128,-0.5,127.5);
      }
      
      for(int cut=0; cut<2; cut++) {
         histo_name = "ClusterFrequency" + det_name;
         if(cut) histo_name += "_CutClusters";
         histo_clusterfreq[det][cut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-0.5,10.5);
      }
      
      for(int trackcut=0; trackcut<3; trackcut++) {
         
         switch(trackcut) {
            case 0: trackcut_name = "_8Hits"; break;
            case 1: trackcut_name = "_8HitsFidcut"; break;
            case 2: trackcut_name = "_No8Hits"; break;
         }
         
         histo_name = "ClustersAverage" + det_name + trackcut_name;
         histo_clusters_average[det][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-5.5,5.5);
         
         for(int clussizeopt=0; clussizeopt<3; clussizeopt++) {
            switch(clussizeopt) {
               case 0: clussizeopt_name = "_NHits"; break;
               case 1: clussizeopt_name = "_NSeeds"; break;
               case 2: clussizeopt_name = "_StdDev"; break;
            }
            histo_name = "ClusterSizeFrequency" + det_name + clussizeopt_name + trackcut_name;
            histo_clustersizefreq[det][clussizeopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-0.5,10.5); //second index: nhits, nseeds, 2nd moment; third index: no fidcut, fidcut, no track req
         }
         
         for(int etaopt=0; etaopt<6; etaopt++) {
            switch(etaopt) {
               case 0: etaopt_name = "_2Largest"; break;
               case 1: etaopt_name = "_2HitCluster"; break;
               case 2: etaopt_name = "_LowCharge"; break;
               case 3: etaopt_name = "_HighCharge"; break;
               case 4: etaopt_name = "_LowCharge_Landau"; break;
               case 5: etaopt_name = "_HighCharge_Landau"; break;
            }
            
            for(int chip=0; chip<2; chip++) {
               if(chip) chip_name = "_RightChip";
               else chip_name = "_LeftChip";
               histo_name = "Eta" + det_name + chip_name + etaopt_name + trackcut_name;
               histo_eta[det][etaopt][trackcut][chip] = new TH1F(histo_name.c_str(),histo_name.c_str(),1001,-0.0005,1.0005); //second index: 2highest,2hit,low,high; third index: track but no fidcut, track and fidcut, no track req
            }
         }
         
         if(trackcut==0) for(int etaopt=0; etaopt<3; etaopt++) {
            switch(etaopt) {
               case 0: etaopt_name = "_2Largest_LeftToRight"; break;
               case 1: etaopt_name = "_2Largest_RightToLeft"; break;
               case 2: etaopt_name = "_2Largest_Average"; break;
            }
            
            for(int chip=0; chip<2; chip++) {
               if(chip) chip_name = "_RightChip";
               else chip_name = "_LeftChip";
               histo_name = "EtaIntegral" + det_name + chip_name + etaopt_name + trackcut_name;
               histo_etaintegral[det][etaopt][chip] = new TH1F(histo_name.c_str(),histo_name.c_str(),1001,-0.0005,1.0005); //second index: left to right, right to left, average
            }
         }
         
         for(int nhitsopt=0; nhitsopt<11; nhitsopt++) {
            
            switch(nhitsopt) {
               case 0: nhitsopt_name = "_1HitClusters"; break;
               case 1: nhitsopt_name = "_2HitClusters"; break;
               case 2: nhitsopt_name = "_3HitClusters"; break;
               case 3: nhitsopt_name = "_4HitClusters"; break;
               case 4: nhitsopt_name = "_5orMoreHitClusters"; break;
               case 5: nhitsopt_name = "_1-2HitClusters"; break;
               case 6: nhitsopt_name = "_1-3HitClusters"; break;
               case 7: nhitsopt_name = "_1-4HitClusters"; break;
               case 8: nhitsopt_name = "_2and3HitClusters"; break;
               case 9: nhitsopt_name = "_3and4HitClusters"; break;
               case 10: nhitsopt_name = "_AllHitClusters"; break;
            }
            
            histo_name = "ClusterCentroidOccupancy" + det_name + nhitsopt_name + trackcut_name;
            if(det==8) histo_clusterocc[det][nhitsopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),128,-0.5,127.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
            else histo_clusterocc[det][nhitsopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5);//-0.5/5,255.5/5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
            
            for(int snr=0; snr<2; snr++) {
               if(snr) {
                  histo_name = "SNRDistribution" + det_name + nhitsopt_name + trackcut_name;
                  if(det==8) histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,-0.5,snr_distribution_di_max+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
                  else histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,-0.5,snr_distribution_si_max+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
               }
               else {
                  histo_name = "PulseHeight" + det_name + nhitsopt_name + trackcut_name;
                  if(det==8) histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,-0.5,pulse_height_di_max+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
                  else histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),pulse_height_num_bins,-0.5,pulse_height_si_max+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
               }
            }
         }
      }
   }
/*	cout << "init new histos..";
	for (int i = 0; i < 5; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_Dia_" << (i+1) << "HitTransparClusters";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),2000,0.,2000.);
		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
		histo_transparentclustering_eta[i] = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),2000,0.,2000.);
	}
	cout << " done." << endl;*/

   
	histo_afc_x = new TH1F("histo_afc_x","histo_afc_x",256,0,255); // for afc (max, 19.11.2010)
	histo_afc_y = new TH1F("histo_afc_y","histo_afc_y",256,0,255); // for afc (max, 19.11.2010)
	histo_afc_clone = new TH2F("histo_afc_clone","histo_afc_clone",256,0,255,256,0,255);
	
	histo_afc_x_cut = new TH1F("histo_afc_x_cut","histo_afc_x_cut",256,0,255);
	histo_afc_y_cut = new TH1F("histo_afc_y_cut","histo_afc_y_cut",256,0,255);

    histo_afc_unit_histo_1f = new TH1F("histo_afc_unit_histo_1f","histo_afc_unit_histo_1f",256,0,255);
	histo_afc_unit_histo_2f = new TH2F("histo_afc_unit_histo_2f","histo_afc_unit_histo_2f",256,0,255,256,0,255);
	histo_afc_region_1_mask = new TH2F("histo_afc_region_1_mask","histo_afc_region_1_mask",256,0,255,256,0,255);
    
	for(int i=0;i<4;i++) {
	 FCR[i] = new FidCutRegion(i);
	 FCR[i]->SetAllValuesZero();
	}
	
   //initialize counters
   total_events = 0, totalsurviving_events = 0, singlesitrack_events = 0, singlesitrack_1diamondclus_events = 0, singlesitrack_fidcut_events = 0, singlesitrack_fidcut_1diamondclus_events = 0, CMNEvents = 0, ZeroDivisorEvents = 0;
   for(int i=0; i<4; i++) detectorxycluster_events[i] = 0;
   for(int d=0; d<10; d++) {
      goldengatecluster_events[d] = 0;
      badchannelcluster_events[d] = 0;
      lumpycluster_events[d] = 0;
      saturatedcluster_events[d] = 0;
   }
   //alignment counters
   counter_alignment_tracks = 0;
   counter_alignment_fidcut_tracks = 0;
   counter_alignment_only_tracks = 0;
   counter_alignment_only_fidcut_tracks = 0;
   counter_alignment_tracks_zero_suppressed = 0;
   
   //Setting Style Attributes for Plots (Starting with Plain style and then making custom adjustments)
   gROOT->SetStyle("Plain"); //General style (see TStyle)
   gStyle->SetOptStat(111110); //Stat options to be displayed			without under- and overflow use gStyle->SetOptStat(1110);
   gStyle->SetOptFit(1111);  //Fit options to be displayed
   gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadTopMargin(0.15);
   //gStyle->SetTitleColor(19,"");
   gStyle->SetStatH(0.12); //Sets Height of Stats Box
   gStyle->SetStatW(0.15); //Sets Width of Stats Box
   gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
   
   //Plot Tag attached to all plots: Gives Run Number, Cut Threshold, Initial Number of Events, and Date/Time of Plot Creation
   char thresh1 [50];
   char datasetsize [50];
   char *pthresh1 = &thresh1[0];
   char *pdatasetsize = &datasetsize[0];
   sprintf(thresh1,"%i.%i", (int)Di_Cluster_Hit_Factor,(int)(Di_Cluster_Hit_Factor*10)%10);
   sprintf(datasetsize,"%i", (int)PedTree->GetEntries());
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
   strcat(pthresh1,pspace);
   strcat(pthresh1,pdatasetsize);
   strcat(pthresh1,peventsinset);
   std::ostringstream run_number_label;
   run_number_label << "Run " << run_number;
   pt = new TPaveText(0.07,0,0.22,0.10,"NDC");  //Normalized CoordinateSystem: Define with x1,y1 is left bottom of box text, x2,y2 is upper right of text box. Goes from 0,0 at bottom left corner of pad to 1,1 of upper right corner
   pt->SetTextSize(0.0250);
   pt->AddText(run_number_label.str().c_str());
   pt->AddText(pthresh1);
   pt->AddText(dateandtime.AsSQLString());
   pt->SetBorderSize(0); //Set Border to Zero
   pt->SetFillColor(0); //Set Fill to White
   
}

Clustering::~Clustering() {
   PedFile->Close();
   //delete PedTree;
   //delete PedFile; // delete this line if root segfaults
   
   //DeleteHistograms();
}


void Clustering::LoadSettings() {

   cout << endl << "Overriding default settings with settings in Settings.ini" << endl << endl;
   
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
      if(key=="snr_plots_enable") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         snr_plots_enable = (int)strtod(value.c_str(),0);
      }
      if(key=="fix_dia_noise") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         fix_dia_noise = (int)strtod(value.c_str(),0);
      }
      if(key=="store_threshold") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         store_threshold = (float)strtod(value.c_str(),0);
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
      if(key=="dia_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         dia_input = (int)strtod(value.c_str(),0);
      }
      if(key=="Si_Pedestal_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Si_Pedestal_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Di_Pedestal_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Di_Pedestal_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Si_Cluster_Seed_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Si_Cluster_Seed_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Di_Cluster_Seed_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Di_Cluster_Seed_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Si_Cluster_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Si_Cluster_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="Di_Cluster_Hit_Factor") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         Di_Cluster_Hit_Factor = (float)strtod(value.c_str(),0);
      }
      if(key=="eta_lowq_slice_low") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         eta_lowq_slice_low = (float)strtod(value.c_str(),0);
      }
      if(key=="eta_lowq_slice_hi") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         eta_lowq_slice_hi = (float)strtod(value.c_str(),0);
      }
      if(key=="eta_hiq_slice_low") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         eta_hiq_slice_low = (float)strtod(value.c_str(),0);
      }
      if(key=="eta_hiq_slice_hi") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         eta_hiq_slice_hi = (float)strtod(value.c_str(),0);
      }
      if(key=="etavsq_n_landau_slices") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         etavsq_n_landau_slices = (float)strtod(value.c_str(),0);
      }
      if(key=="alignment_x_offsets") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseFloatArray(value,alignment_x_offsets);
      }
      if(key=="alignment_y_offsets") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseFloatArray(value,alignment_y_offsets);
      }
      if(key=="alignment_phi_offsets") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseFloatArray(value,alignment_phi_offsets);
      }
      if(key=="alignment_z_offsets") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseFloatArray(value,alignment_z_offsets);
      }
      if(key=="alignment_training_track_fraction") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         alignment_training_track_fraction = (float)strtod(value.c_str(),0);
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
      if(key=="si_avg_fidcut_xlow") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         si_avg_fidcut_xlow = (int)strtod(value.c_str(),0);
      }
      if(key=="si_avg_fidcut_xhigh") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         si_avg_fidcut_xhigh = (int)strtod(value.c_str(),0);
      }
      if(key=="si_avg_fidcut_ylow") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         si_avg_fidcut_ylow = (int)strtod(value.c_str(),0);
      }
      if(key=="si_avg_fidcut_yhigh") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         si_avg_fidcut_yhigh = (int)strtod(value.c_str(),0);
      }
      if(key=="pulse_height_num_bins") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         pulse_height_num_bins = (int)strtod(value.c_str(),0);
      }
      if(key=="pulse_height_si_max") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         pulse_height_si_max = (int)strtod(value.c_str(),0);
      }
      if(key=="pulse_height_di_max") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         pulse_height_di_max = (int)strtod(value.c_str(),0);
      }
      if(key=="snr_distribution_si_max") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         snr_distribution_si_max = (int)strtod(value.c_str(),0);
      }
      if(key=="snr_distribution_di_max") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         snr_distribution_di_max = (int)strtod(value.c_str(),0);
      }
	   if (key == "alignment_chi2") {
		   cout << key.c_str() << " = " << value.c_str() << endl;
		   alignment_chi2 = (Float_t)strtod(value.c_str(),0);
	   }
	   if (key == "UseAutoFidCut") {
		   cout << key.c_str() << " = " << value.c_str() << endl;
		   UseAutoFidCut = (bool)strtod(value.c_str(),0);
	   }
	   if (key == "AlternativeClustering") {
		   cout << key.c_str() << " = " << value.c_str() << endl;
		   AlternativeClustering = (bool)strtod(value.c_str(),0);
	   }
   }
   
   file.close();
   cout<<endl<<"Finished importing settings from "<<settings_file<<endl<<endl;
}

void Clustering::ParseIntArray(string value, vector<int> &vec) {
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

void Clustering::ParseFloatArray(string value, vector<float> &vec) {
   int index=0;
   string::size_type offset1 = value.find_first_of('{')+1;
   string::size_type offset2 = value.find_first_of(',');
   
   vec.push_back((float)strtod(value.substr(offset1,offset2-offset1).c_str(),0));
   
   //cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(offset1,offset2-offset1)="<<value.substr(offset1,offset2-offset1)<<"\tstrtod(value.substr(offset1,offset2-offset1),0)="<<strtod(value.substr(offset1,offset2-offset1).c_str(),0)<<endl;//check
   
   value = value.substr(offset2+1,value.length()-(offset2+1));
   
   while(value.length()>2) {
      offset2 = TMath::Min(value.find_first_of(','),value.find_first_of('}'));
      vec.push_back((float)strtod(value.substr(0,offset2).c_str(),0));
      index++;
      //cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(0,offset2)="<<value.substr(0,offset2)<<"\tstrtod(value.substr(0,offset2),0)="<<strtod(value.substr(0,offset2).c_str(),0)<<endl;//check
      if(value.find_first_of(';')<2) break;
      value = value.substr(offset2+1,value.length()-(offset2+1));
   }
}
      
//take current event and cluster
void Clustering::ClusterEvent(bool verbose) {
   PedTree->GetEvent(current_event);
   if(verbose) cout<<endl<<endl;
   if(current_event%1000==0) cout<<"Clustering::ClusterEvent(): current_event = "<<current_event<<endl;
   
   clustered_event.Clear();
   clustered_event.SetEventNumber(current_event);
   
   vector<int> hits, cluster, badchannelclusterflags, goldengateclusterflags, lumpyclusterflags, saturatedclusterflags;
   vector< vector<int> > clusters;
   int previouschan, currentchan, hasaseed, hasmasked, previousseed, isgoldengate, islumpy, peakchan, hassaturated;
   float peakchan_psadc, currentchan_psadc, previouschan_psadc;
   
   for(int det=0; det<9; det++) {
      hits.clear();
      cluster.clear();
      clusters.clear();
      
      //look for hits
      if(verbose) cout<<endl<<endl<<"Detector "<<det<<" hits: ";
      for(int i=0; i<(int)Det_NChannels[det]; i++) {
         if(det<8) if(Det_ADC[det][i]-Det_PedMean[det][i] > Si_Cluster_Hit_Factor*Det_PedWidth[det][i]) {
            hits.push_back(i);
            if(verbose) {
               cout<<(int)Det_Channels[det][i];
               if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
               cout<<", ";
            }
         }
         if(det==8) if(Dia_ADC[i]-Det_PedMean[det][i] > Di_Cluster_Hit_Factor*Det_PedWidth[det][i]) {
            hits.push_back(i);
            if(verbose) {
               cout<<(int)Det_Channels[det][i];
               if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
               cout<<", ";
            }
         }
      }
      if(verbose) {
         cout<<endl<<"Channels screened: ";
         for(int i=0; i<256; i++) if(!Det_channel_screen[det].CheckChannel(i)) cout<<i<<", ";
         cout<<endl<<hits.size()<<" hits found in "<<(int)Det_NChannels[det]<<" saved channels"<<endl;
      }
      if(hits.size()==0) {
         if(verbose) cout<<"No hits found so skipping to next detector."<<endl;
         continue;
      }
      hits.push_back(-1);
      if(verbose) cout<<"pushed back -1 as hit for determining end of hits"<<endl;
      
      //now look for contiguous regions in hits and save clustered hits with seeds
      if(verbose) cout<<"before clear(): cluster.size()="<<cluster.size()<<" and cluster.size()="<<cluster.size()<<endl;
      cluster.clear();
      clusters.clear();
      goldengateclusterflags.clear();
      badchannelclusterflags.clear();
      lumpyclusterflags.clear();
      saturatedclusterflags.clear();
      if(verbose) cout<<"after clear(): cluster.size()="<<cluster.size()<<" and cluster.size()="<<cluster.size()<<endl;
      previouschan=-1;
      for(uint i=0; i<hits.size(); i++) {
         currentchan = Det_Channels[det][hits[i]];
         if(verbose) {
            if(hits[i]==-1) cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or end of hits"<<endl;
            else {
               cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or channel "<<currentchan<<" (";
               if(det==8) currentchan_psadc = Dia_ADC[hits[i]]-Det_PedMean[det][hits[i]];
               if(det<8) currentchan_psadc = Det_ADC[det][hits[i]]-Det_PedMean[det][hits[i]];
               cout<<currentchan_psadc<<" psadc, "<<Det_PedWidth[det][hits[i]]<<" pedrms, "<<currentchan_psadc/Det_PedWidth[det][hits[i]]<<" snr)"<<endl;
            }
         }
         //build a cluster of hits
         if((previouschan==-1 || currentchan==previouschan+1) && hits[i]!=-1) {
            if(verbose) cout<<"adding channel to cluster"<<endl;
            cluster.push_back(hits[i]);
            previouschan = currentchan;
         }
         //found end of cluster so search current cluster for a seed
         else {
            if(hits[i]!=-1) i--;
            if(verbose) cout<<"found end of cluster; looking for seed:"<<endl;
            hasaseed=0;
            hasmasked=0;
            previousseed=-1;
            isgoldengate=0;
            islumpy=0;
            hassaturated=0;
            peakchan=-1;
            peakchan_psadc=-5000;
            
            //require no masked channels adjacent to the cluster
            if((int)Det_Channels[det][cluster[0]]==0 || (int)Det_Channels[det][cluster[cluster.size()-1]]==255) {
               hasmasked = 1;
               if(verbose) cout<< "Cluster is up against edge of detector; flagging cluster as bad channel cluster." << endl;
            }
            else if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][cluster[0]]-1)
               || !Det_channel_screen[det].CheckChannel((int)Det_Channels[det][cluster[cluster.size()-1]]+1)) {
               hasmasked = 1; 
               if(verbose) cout<< "Channel(s) adjacent to the cluster is masked; flagging cluster as bad channel cluster." << endl;
            }
            
            for(uint j=0; j<cluster.size(); j++) {
               currentchan = cluster[j];
               if(verbose) cout<<"cluster["<<j<<"]="<<cluster[j]<<" or channel "<<(int)Det_Channels[det][currentchan];
               if(det<8) {
                  currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
                  //check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
                  if(cluster.size()>2) if(currentchan_psadc > peakchan_psadc) {
                     peakchan_psadc = currentchan_psadc;
                     peakchan = j;
                  }
                  //check for seed
                  if(currentchan_psadc > Si_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
                     hasaseed=1; 
                     if(verbose) cout<<" is a seed";
                     if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
                        isgoldengate=1;
                        if(verbose) cout<<" (goldengate cluster)";
                     }
                     else previousseed=Det_Channels[det][currentchan];
                     //check to see if the hit saturated the adc
                     if(Det_ADC[det][currentchan]>254) {
                        hassaturated=1;
                        if(verbose) cout<<" (saturated adc)";
                     }
                  }
                  //require no masked hits in silicon
                  if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked hit)";}
               }
               if(det==8) {
                  currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
                  //check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
                  if(cluster.size()>2) if(currentchan_psadc > peakchan_psadc) {
                     peakchan_psadc = currentchan_psadc;
                     peakchan = j;
                  }
                  //check for seed
                  if(currentchan_psadc > Di_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
                     hasaseed=1; 
                     if(verbose) cout<<" is a seed";
                     if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
                        isgoldengate=1;
                        if(verbose) cout<<" (goldengate cluster)";
                     }
                     else previousseed=Det_Channels[det][currentchan];
                     //check to see if the hit saturated the adc
                     if(Dia_ADC[currentchan]>4094) {
                        hassaturated=1;
                        if(verbose) cout<<" (saturated adc)";
                     }
                     //require no masked seeds in diamond
                     if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked seed)";}
                  }
               }
               if(verbose) cout<<endl;
            }
            //Lumpy cluster check (can't have lumpy cluster with less than 3 hits)
            if(cluster.size()>2) {
               if(verbose) cout<<"Found peak at cluster hit "<<peakchan<<" or channel index "<<cluster[peakchan]<<" with PSADC of "<<peakchan_psadc<<endl;
               //now scan right of peak to see if monotonically decreasing
               previouschan_psadc = peakchan_psadc;
               for(uint j=peakchan+1; j<cluster.size(); j++) {
                  currentchan = cluster[j];
                  if(det<8) currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
                  if(det==8) currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
                  if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
                  if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)Det_Channels[det][currentchan]<<")"<<endl;}
                  previouschan_psadc = currentchan_psadc;
               }
               //now scan left of peak to see if monotonically decreasing
               previouschan_psadc = peakchan_psadc;
               for(int j=peakchan-1; j>=0; j--) {
                  currentchan = cluster[j];
                  if(det<8) currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
                  if(det==8) currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
                  if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
                  if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)Det_Channels[det][currentchan]<<")"<<endl;}
                  previouschan_psadc = currentchan_psadc;
               }
            }
            //if there's a seed in the cluster, save it
            if(hasaseed) {
               clusters.push_back(cluster); 
               badchannelclusterflags.push_back(hasmasked); 
               goldengateclusterflags.push_back(isgoldengate); 
               lumpyclusterflags.push_back(islumpy);
               saturatedclusterflags.push_back(hassaturated);
               if(verbose) {
                  cout<<"storing cluster"<<endl;
                  if(hasmasked) cout<<"Flagged as BadChannelCluster!"<<endl;
                  if(isgoldengate) cout<<"Flagged as GoldenGateCluster!"<<endl;
                  if(islumpy) cout<<"Flagged as LumpyCluster!"<<endl;
                  if(hassaturated) cout<<"Flagged as SaturatedCluster!"<<endl;
               }
            }
            cluster.clear(); //start storing a new cluster
            previouschan=-1; //save next channel
            //move on to the next cluster
         }//end else (found end of candidate cluster)
      }//end loop over hits in detector
      
      if(clusters.size()==0) {
         if(verbose) cout<<"No clusters found so skipping to next detector."<<endl;
         continue;
      }
      
      //now that we have lists of channels belonging to clusters, let's create Cluster objects to store the data
      Cluster* current_cluster = 0;
      int highest_index, nexthighest_index;
      float highest_psadc, nexthighest_psadc, current_psadc;
      if(clusters.size()>0) {
         if(verbose) cout<<"det="<<det<<"\tclusters.size()="<<clusters.size()<<endl;
         for(uint i=0; i<clusters.size(); i++) {
            //add cluster to a different list in current event depending on flags
            current_cluster = clustered_event.AddCluster(det,badchannelclusterflags[i]); //||goldengateclusterflags[i]||lumpyclusterflags[i]||saturatedclusterflags[i]);
            //flag cluster
            current_cluster->FlagBadChannelCluster(badchannelclusterflags[i]);
            current_cluster->FlagGoldenGateCluster(goldengateclusterflags[i]);
            current_cluster->FlagLumpyCluster(lumpyclusterflags[i]);
            current_cluster->FlagSaturatedCluster(saturatedclusterflags[i]);
            //calculate transparent eta while saving each channel to cluster (note that due to "zero supression" in the saved pedestal info, not all clusters will have an eta and will be reported as -1)
            highest_index = -1; nexthighest_index = -1;
            highest_psadc = 0; nexthighest_psadc = 0; current_psadc = 0;
            for(uint j=0; j<clusters[i].size(); j++) { //save each cluster one channel at a time
               currentchan = clusters[i][j];
               if(det<8) {
                  current_cluster->AddHit(Det_Channels[det][currentchan], Det_ADC[det][currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
                  current_psadc = Det_ADC[det][currentchan] - Det_PedMean[det][currentchan];
               }
               if(det==8) {
                  current_cluster->AddHit(Det_Channels[det][currentchan], Dia_ADC[currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
                  current_psadc = Dia_ADC[currentchan] - Det_PedMean[det][currentchan];
               }
               if(verbose) cout<<"Detector "<<det<<": cluster_saved/all_clusters="<<i+1<<"/"<<clusters.size()<<"\tchan_to_save="<<j+1<<"/"<<clusters[i].size()<<"\tcurrentchan="<<(int)Det_Channels[det][currentchan]<<endl;
               //locate highest seed
               if(current_psadc>highest_psadc) {
                  highest_index = currentchan;
                  highest_psadc = current_psadc;
               }
            }//end loop over channels in a cluster to save
            //check which channel has next highest psadc
            if(highest_index<(int)Det_NChannels[det]-1) if(Det_Channels[det][highest_index+1]==Det_Channels[det][highest_index]+1 && Det_channel_screen[det].CheckChannel(Det_Channels[det][highest_index+1])) {
               //first if the next channel is available and an ok channel, just assume it has the nexthighest_psadc
               nexthighest_index = highest_index+1;
               if(det==8) nexthighest_psadc = Dia_ADC[nexthighest_index] - Det_PedMean[det][nexthighest_index];
               else nexthighest_psadc = Det_ADC[det][nexthighest_index] - Det_PedMean[det][nexthighest_index]; 
            }
            if(highest_index>0) if(Det_Channels[det][highest_index-1]==Det_Channels[det][highest_index]-1 && Det_channel_screen[det].CheckChannel(Det_Channels[det][highest_index-1])) {
               //now if the previous channel is available and an ok channel, check whether it has a higher psadc than the next one
               if(det==8) current_psadc = Dia_ADC[highest_index-1] - Det_PedMean[det][highest_index-1];
               else current_psadc = Det_ADC[det][highest_index-1] - Det_PedMean[det][highest_index-1];
               if(current_psadc>nexthighest_psadc) {
                  nexthighest_index = highest_index-1;
                  nexthighest_psadc = current_psadc;
               }
            }
            //calculate eta and centroid for highest 2 channels (used later for alignment and tracking)
            if(nexthighest_index>-1) {
               //eta
               if(highest_index>nexthighest_index) current_cluster->SetEta(highest_psadc/(highest_psadc+nexthighest_psadc));
               else current_cluster->SetEta(nexthighest_psadc/(highest_psadc+nexthighest_psadc));
               
               //centroid
               current_cluster->highest2_centroid=(Det_Channels[det][highest_index]*highest_psadc+Det_Channels[det][nexthighest_index]*nexthighest_psadc)/(highest_psadc+nexthighest_psadc);
               
               /*
               if(totalsurviving_events%1000==0) {
                  cout<<"Event "<<&current_event-1<<"; detector "<<det<<":"<<endl; 
                  cout<<"highest_index = "<<highest_index<<"\thighest_psadc="<<highest_psadc<<endl;
                  cout<<"nexthighest_index = "<<nexthighest_index<<"\tnexthighest_psadc="<<nexthighest_psadc<<endl;
                  cout<<"highest_psadc/(highest_psadc+nexthighest_psadc) = "<<highest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
                  cout<<"nexthighest_psadc/(highest_psadc+nexthighest_psadc) = "<<nexthighest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
               }
               */
            }
         }//end loop over clusters
      }
   }//end loop over detectors
   
   //get ready for next event
   current_event++;
   
}

// alternative function to ClusterEvent
// channels are searched for a seed. clusters are build around located seeds.
void Clustering::ClusterEventSeeds(bool verbose) {
	PedTree->GetEvent(current_event);
	if(verbose) cout<<endl<<endl;
	if(current_event%1000==0) cout<<"Clustering::ClusterEvent(): current_event = "<<current_event<<endl;
	
	clustered_event.Clear();
	clustered_event.SetEventNumber(current_event);
	
	vector<int> hits, cluster, badchannelclusterflags, goldengateclusterflags, lumpyclusterflags, saturatedclusterflags, seeds, cluster_channels;
	vector< vector<int> > clusters;
	int previouschan, currentchan, hasaseed, hasmasked, previousseed, isgoldengate, islumpy, peakchan, hassaturated;
	float peakchan_psadc, currentchan_psadc, previouschan_psadc;
	
	// -- loop over detectors
	for(int det=0; det<9; det++) {
		hits.clear();
		cluster.clear();
		clusters.clear();
		seeds.clear();
		
		// -- loop over channels and look for seeds
		if (verbose) cout << endl << endl << "Detector " << det << " seeds: ";
		for (int i = 0; i < (int)Det_NChannels[det]; i++) {
			if (det < 8 && Det_ADC[det][i]-Det_PedMean[det][i] > Si_Cluster_Seed_Factor*Det_PedWidth[det][i]) {
				seeds.push_back(i);
				if (verbose) {
					cout<<(int)Det_Channels[det][i];
					if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
					cout<<", ";
				}
			}
			if (det == 8 && Dia_ADC[i]-Det_PedMean[det][i] > Di_Cluster_Seed_Factor*Det_PedWidth[det][i]) {
				seeds.push_back(i);
				if (verbose) {
					cout<<(int)Det_Channels[det][i];
					if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
					cout<<", ";
				}
			}
		}	// end loop over channels
		
		goldengateclusterflags.clear();
		badchannelclusterflags.clear();
		lumpyclusterflags.clear();
		saturatedclusterflags.clear();
		cluster_channels.clear();
		
		/*		for (int i=0; i < seeds.size()-1; i++) {
		 if (TMath::Abs(seeds[i]-seeds[i+1]) > 1) {
		 goldengateclusterflags.push_back(seeds[i]);
		 }
		 }*/
		
		// detect cluster for every seed
		for (int i = 0; i < seeds.size(); i++) {
			bool channel_already_used = false;
			int first_hit = 0;
			int last_hit = 0;
			hasaseed=1;
			hasmasked=0;
			previousseed=-1;
			isgoldengate=0;
			islumpy=0;
			hassaturated=0;
			peakchan=-1;
			peakchan_psadc=-5000;
			cluster.clear();
			
			// check if channel is already used in another cluster
			for (int j = 0; j < cluster_channels.size(); j++) {
				if (seeds[i] == cluster_channels[j]) channel_already_used = true;
			}
			if (channel_already_used) continue;
			
			if (det < 8) {
				// look for first hit of cluster
				first_hit = seeds[i];
				for (int j = seeds[i]-1; j >= 0 && Det_ADC[det][j]-Det_PedMean[det][j] > Si_Cluster_Hit_Factor*Det_PedWidth[det][j]; j--) {
					first_hit = j;
				}
				
				// look for last hit of cluster and write hits to cluster
				for (int j = first_hit; j < (int)Det_NChannels[det] && Det_ADC[det][j]-Det_PedMean[det][j] > Si_Cluster_Hit_Factor*Det_PedWidth[det][j]; j++) {
					cluster.push_back(j);
					cluster_channels.push_back(j);
				}
				
				
/*				for (int j = seeds[i]-1; j >= 0 && Det_ADC[det][j]-Det_PedMean[det][j] > Si_Cluster_Hit_Factor*Det_PedWidth[det][j]; j--) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (Det_ADC[det][j]-Det_PedMean[det][j] > Det_ADC[det][j+1]-Det_PedMean[det][j+1])
							islumpy = 1;
					}
				}
				for (int j = seeds[i]+1; j >= 0 && Det_ADC[det][j]-Det_PedMean[det][j] > Si_Cluster_Hit_Factor*Det_PedWidth[det][j]; j++) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (Det_ADC[det][j]-Det_PedMean[det][j] > Det_ADC[det][j-1]-Det_PedMean[det][j-1])
							islumpy = 1;
					}
				}*/
			}
			if (det == 8) {
				// look for first hit of the cluster
				first_hit = seeds[i];
				for (int j = seeds[i]-1; j >= 0 && Dia_ADC[j]-Det_PedMean[det][j] > Di_Cluster_Hit_Factor*Det_PedWidth[det][j]; j--) {
					first_hit = j;
				}
				
				// look for last hit of cluster and write hits to cluster
				for (int j = first_hit; j < (int)Det_NChannels[det] && Dia_ADC[j]-Det_PedMean[det][j] > Di_Cluster_Hit_Factor*Det_PedWidth[det][j]; j++) {
					cluster.push_back(j);
					cluster_channels.push_back(j);
				}
				
				
/*				for (int j = seeds[i]-1; j >= 0 && Dia_ADC[j]-Det_PedMean[det][j] > Di_Cluster_Hit_Factor*Det_PedWidth[det][j]; j--) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (Dia_ADC[j]-Det_PedMean[det][j] > Dia_ADC[j+1]-Det_PedMean[det][j+1])
							islumpy = 1;
					}
				}
				for (int j = seeds[i]+1; j >= 0 && Dia_ADC[j]-Det_PedMean[det][j] > Di_Cluster_Hit_Factor*Det_PedWidth[det][j]; j++) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (Dia_ADC[j]-Det_PedMean[det][j] > Dia_ADC[j-1]-Det_PedMean[det][j-1])
							islumpy = 1;
					}
				}*/
			}
			if (cluster.size() > 0) {
				
				//require no masked channels adjacent to the cluster
				if((int)Det_Channels[det][cluster[0]]==0 || (int)Det_Channels[det][cluster[cluster.size()-1]]==255) {
					hasmasked = 1;
					if(verbose) cout<< "Cluster is up against edge of detector; flagging cluster as bad channel cluster." << endl;
				}
				else if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][cluster[0]]-1)
						|| !Det_channel_screen[det].CheckChannel((int)Det_Channels[det][cluster[cluster.size()-1]]+1)) {
					hasmasked = 1; 
					if(verbose) cout<< "Channel(s) adjacent to the cluster is masked; flagging cluster as bad channel cluster." << endl;
				}
				
				// -- look for peak, ckeck golden gate, check saturation and require no masked hits in silicon
				for(uint j=0; j<cluster.size(); j++) {
					currentchan = cluster[j];
					if(verbose) cout<<"cluster["<<j<<"]="<<cluster[j]<<" or channel "<<(int)Det_Channels[det][currentchan];
					if(det<8) {
						currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
						
						//check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
//						if(cluster.size()>2) {
							if(currentchan_psadc > peakchan_psadc) {
								peakchan_psadc = currentchan_psadc;
								peakchan = j;
							}
//						}
						
						// check for golden gate
						if(currentchan_psadc > Si_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
							if(verbose) cout<<" is a seed";
							if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
								isgoldengate=1;
								if(verbose) cout<<" (goldengate cluster)";
							}
							previousseed=Det_Channels[det][currentchan];
						}
						
						//check to see if the hit saturated the adc
						if(Det_ADC[det][currentchan]>254) {
							hassaturated=1;
							if(verbose) cout<<" (saturated adc)";
						}
						
						//require no masked hits in silicon
						if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked hit)";}
					}
					if(det==8) {
						currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
						
						//check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
//						if(cluster.size()>2)
						if(currentchan_psadc > peakchan_psadc) {
							peakchan_psadc = currentchan_psadc;
							peakchan = j;
						}
						
						//check for golden gate
						if(currentchan_psadc > Di_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
							if(verbose) cout<<" is a seed";
							if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
								isgoldengate=1;
								if(verbose) cout<<" (goldengate cluster)";
							}
							previousseed=Det_Channels[det][currentchan];
							//require no masked seeds in diamond
							if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked seed)";}
						}
						
						//check to see if the hit saturated the adc
						if(Dia_ADC[currentchan]>4094) {
							hassaturated=1;
							if(verbose) cout<<" (saturated adc)";
						}
					}
					if(verbose) cout<<endl;
				}
				
				//Lumpy cluster check (can't have lumpy cluster with less than 3 hits)
				if(cluster.size()>2) {
					if(verbose) cout<<"Found peak at cluster hit "<<peakchan<<" or channel index "<<cluster[peakchan]<<" with PSADC of "<<peakchan_psadc<<endl;
					//now scan right of peak to see if monotonically decreasing
					previouschan_psadc = peakchan_psadc;
					for(uint j=peakchan+1; j<cluster.size(); j++) {
						currentchan = cluster[j];
						if(det<8) currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
						if(det==8) currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
						if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
						if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)Det_Channels[det][currentchan]<<")"<<endl;}
						previouschan_psadc = currentchan_psadc;
					}
					//now scan left of peak to see if monotonically decreasing
					previouschan_psadc = peakchan_psadc;
					for(int j=peakchan-1; j>=0; j--) {
						currentchan = cluster[j];
						if(det<8) currentchan_psadc = Det_ADC[det][currentchan]-Det_PedMean[det][currentchan];
						if(det==8) currentchan_psadc = Dia_ADC[currentchan]-Det_PedMean[det][currentchan];
						if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
						if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)Det_Channels[det][currentchan]<<")"<<endl;}
						previouschan_psadc = currentchan_psadc;
					}
				}
				
				//if there's a seed in the cluster, save it
				if(hasaseed) {
					clusters.push_back(cluster); 
					badchannelclusterflags.push_back(hasmasked); 
					goldengateclusterflags.push_back(isgoldengate); 
					lumpyclusterflags.push_back(islumpy);
					saturatedclusterflags.push_back(hassaturated);
					if(verbose) {
						cout<<"storing cluster"<<endl;
						if(hasmasked) cout<<"Flagged as BadChannelCluster!"<<endl;
						if(isgoldengate) cout<<"Flagged as GoldenGateCluster!"<<endl;
						if(islumpy) cout<<"Flagged as LumpyCluster!"<<endl;
						if(hassaturated) cout<<"Flagged as SaturatedCluster!"<<endl;
					}
				}
			}
		} // end loop over seeds
		
		//now that we have lists of channels belonging to clusters, let's create Cluster objects to store the data
		Cluster* current_cluster = 0;
		int highest_index, nexthighest_index;
		float highest_psadc, nexthighest_psadc, current_psadc;
		if(clusters.size()>0) {
			if(verbose) cout<<"det="<<det<<"\tclusters.size()="<<clusters.size()<<endl;
			for(uint i=0; i<clusters.size(); i++) {
				//add cluster to a different list in current event depending on flags
				current_cluster = clustered_event.AddCluster(det,badchannelclusterflags[i]); //||goldengateclusterflags[i]||lumpyclusterflags[i]||saturatedclusterflags[i]);
				//flag cluster
				current_cluster->FlagBadChannelCluster(badchannelclusterflags[i]);
				current_cluster->FlagGoldenGateCluster(goldengateclusterflags[i]);
				current_cluster->FlagLumpyCluster(lumpyclusterflags[i]);
				current_cluster->FlagSaturatedCluster(saturatedclusterflags[i]);
				//calculate transparent eta while saving each channel to cluster (note that due to "zero supression" in the saved pedestal info, not all clusters will have an eta and will be reported as -1)
				highest_index = -1; nexthighest_index = -1;
				highest_psadc = 0; nexthighest_psadc = 0; current_psadc = 0;
				for(uint j=0; j<clusters[i].size(); j++) { //save each cluster one channel at a time
					currentchan = clusters[i][j];
					if(det<8) {
						current_cluster->AddHit(Det_Channels[det][currentchan], Det_ADC[det][currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
						current_psadc = Det_ADC[det][currentchan] - Det_PedMean[det][currentchan];
					}
					if(det==8) {
						current_cluster->AddHit(Det_Channels[det][currentchan], Dia_ADC[currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
						current_psadc = Dia_ADC[currentchan] - Det_PedMean[det][currentchan];
					}
					if(verbose) cout<<"Detector "<<det<<": cluster_saved/all_clusters="<<i+1<<"/"<<clusters.size()<<"\tchan_to_save="<<j+1<<"/"<<clusters[i].size()<<"\tcurrentchan="<<(int)Det_Channels[det][currentchan]<<endl;
					//locate highest seed
					if(current_psadc>highest_psadc) {
						highest_index = currentchan;
						highest_psadc = current_psadc;
					}
				}//end loop over channels in a cluster to save
				//check which channel has next highest psadc
				if(highest_index<(int)Det_NChannels[det]-1) if(Det_Channels[det][highest_index+1]==Det_Channels[det][highest_index]+1 && Det_channel_screen[det].CheckChannel(Det_Channels[det][highest_index+1])) {
					//first if the next channel is available and an ok channel, just assume it has the nexthighest_psadc
					nexthighest_index = highest_index+1;
					if(det==8) nexthighest_psadc = Dia_ADC[nexthighest_index] - Det_PedMean[det][nexthighest_index];
					else nexthighest_psadc = Det_ADC[det][nexthighest_index] - Det_PedMean[det][nexthighest_index]; 
				}
				if(highest_index>0) if(Det_Channels[det][highest_index-1]==Det_Channels[det][highest_index]-1 && Det_channel_screen[det].CheckChannel(Det_Channels[det][highest_index-1])) {
					//now if the previous channel is available and an ok channel, check whether it has a higher psadc than the next one
					if(det==8) current_psadc = Dia_ADC[highest_index-1] - Det_PedMean[det][highest_index-1];
					else current_psadc = Det_ADC[det][highest_index-1] - Det_PedMean[det][highest_index-1];
					if(current_psadc>nexthighest_psadc) {
						nexthighest_index = highest_index-1;
						nexthighest_psadc = current_psadc;
					}
				}
				//calculate eta and centroid for highest 2 channels (used later for alignment and tracking)
				if(nexthighest_index>-1) {
					//eta
					if(highest_index>nexthighest_index) current_cluster->SetEta(highest_psadc/(highest_psadc+nexthighest_psadc));
					else current_cluster->SetEta(nexthighest_psadc/(highest_psadc+nexthighest_psadc));
					
					//centroid
					current_cluster->highest2_centroid=(Det_Channels[det][highest_index]*highest_psadc+Det_Channels[det][nexthighest_index]*nexthighest_psadc)/(highest_psadc+nexthighest_psadc);
					
					/*
					 if(totalsurviving_events%1000==0) {
					 cout<<"Event "<<&current_event-1<<"; detector "<<det<<":"<<endl; 
					 cout<<"highest_index = "<<highest_index<<"\thighest_psadc="<<highest_psadc<<endl;
					 cout<<"nexthighest_index = "<<nexthighest_index<<"\tnexthighest_psadc="<<nexthighest_psadc<<endl;
					 cout<<"highest_psadc/(highest_psadc+nexthighest_psadc) = "<<highest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
					 cout<<"nexthighest_psadc/(highest_psadc+nexthighest_psadc) = "<<nexthighest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
					 }
					 */
				}
			}//end loop over clusters
		}
		
	}	// end loop over detectors
	current_event++;
	// -- end of lukas' code
	
}

void SetDuckStyle() {
	TStyle* DuckStyle = new TStyle("DuckStyle", "Famous Duck Style");
	DuckStyle->SetOptStat(1110); //Stat options to be displayed
	//	DuckStyle->SetOptFit(1111);  //Fit options to be displayed
	DuckStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
	DuckStyle->SetPadRightMargin(0.15);
	DuckStyle->SetPadTopMargin(0.15);
	DuckStyle->SetTitleColor(kBlack);
	DuckStyle->SetFrameLineWidth(0.05);
	DuckStyle->SetStatFont(42);
	DuckStyle->SetStatBorderSize(0);
	DuckStyle->SetPadBorderSize(0);
	DuckStyle->SetPadTopMargin(1);
	DuckStyle->SetLegendBorderSize(0);
	DuckStyle->SetStatFontSize(0.02);
	DuckStyle->SetStatStyle(0);
	DuckStyle->SetStatH(0.12); //Sets Height of Stats Box
	DuckStyle->SetStatW(0.15); //Sets Width of Stats Box
	DuckStyle->SetStatX(0.9);	
	DuckStyle->SetStatY(0.97);
	DuckStyle->SetTitleOffset(1.0,"Y");
	DuckStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)	
	DuckStyle->SetCanvasBorderMode(0);
	DuckStyle->SetTitleFont(42,"XYZ");
	DuckStyle->SetTitleFontSize(0.038);
	//	DuckStyle->SetTitleTextSize(0.03);
	DuckStyle->SetTitleTextColor(kBlack);
	DuckStyle->SetFrameLineStyle(0);
	DuckStyle->SetGridStyle(0);
	DuckStyle->SetHatchesLineWidth(0);
	//	DuckStyle->SetOptTitle(0);
	DuckStyle->SetPadTickX(0);
	DuckStyle->SetPadTickY(0);
	DuckStyle->SetTitleX(0.07);
	DuckStyle->SetTitleY(0.925);
	
	DuckStyle->SetTitleSize(0.02,"XYZ");
	DuckStyle->SetLineWidth(0.05);
	DuckStyle->SetHistLineWidth(1);
	//	DuckStyle->SetTitleStyle(0);
	DuckStyle->SetTitleBorderSize(0);
	DuckStyle->SetTitleFillColor(0);
	DuckStyle->SetHistLineWidth(0.05);
	
	//	DuckStyle->SetTickLength(0,"XY");
	//	gStyle->SetBarOffset(0.5);
	//	gStyle->SetStatFont(42);
	//	gStyle->SetTextSize(0.01);
	DuckStyle->SetLabelFont(42,"XYZ");
	DuckStyle->SetLabelColor(kBlack,"XYZ");
	DuckStyle->SetLabelSize(0.025,"XYZ");
	//DuckStyle->SetTitleOffset(1.8, "Y"); // Another way to set the Offset
	//	gStyle->SetTitleOffset(1.2, "X"); // Another way to set the Offset
	DuckStyle->SetTitleOffset(1.2,"X");
	
	DuckStyle->cd();
	
	cout << "Using DuckStyle" << endl;
}

void Clustering::SaveHistogram(TH1F* histo) {
   SaveHistogramPNG(histo);
   SaveHistogramROOT(histo);
}

void Clustering::SaveHistogram(TH2F* histo) {
   SaveHistogramPNG(histo);
   SaveHistogramROOT(histo);
}

void Clustering::SaveHistogramPDF(TH1F* histo) {
	TCanvas plots_canvas("plots_canvas","plots_canvas");
	plots_canvas.cd();
	histo->Draw();
	pt->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas.Print(plot_filename.str().c_str());
}

void Clustering::SaveHistogramPDF(TH2F* histo) {
	TCanvas plots_canvas("plots_canvas","plots_canvas");
	//plots_canvas.cd();
//	SetDuckStyle();
	plots_canvas.cd();
	cout << "Using SaveHistogrammPDF on TH2F histogram " << histo->GetName() << endl;
	//histo->Draw();
	gStyle->SetTitleFont(42);
	gStyle->SetMarkerSize(0);
	pt->SetTextSize(0.0250);
	pt->SetTextColor(kBlack);
	histo->SetTitleFont(42);
	histo->UseCurrentStyle();
	histo->Draw("colz");
	pt->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas.Print(plot_filename.str().c_str());	
	//pt->SetTextSize(0.1);
}

void Clustering::SaveHistogramPNG(TH1F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".png";
   plots_canvas.Print(plot_filename.str().c_str());
}

void Clustering::SaveHistogramROOT(TH1F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   pt->Draw();
   ostringstream plot_filename;
	ostringstream histo_filename;
//	histo_filename << plots_path << histo->GetName() << "_histo.root";
	histo_filename << plots_path << "histograms.root";
   plot_filename << plots_path << histo->GetName() << ".root";
   plots_canvas.Print(plot_filename.str().c_str());
	TFile f(histo_filename.str().c_str(),"UPDATE");
	histo->Write();
}

void Clustering::SaveHistogramPNG(TH2F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   histo->Draw("colz");
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".png";
   plots_canvas.Print(plot_filename.str().c_str());
}

void Clustering::SaveHistogramROOT(TH2F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   histo->Draw("colz");
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".root";
   plots_canvas.Print(plot_filename.str().c_str());
}

void Clustering::InitializeHistograms() {
   for(int i=0; i<2; i++) histo_dianoise[i] = 0;
   for(int i=0; i<5; i++) for(int j=0; j<3; j++) histo_scatter[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) histo_clusterocc[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) histo_landau[i][j][k][l] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) histo_clusterfreq[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) histo_clustersizefreq[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<6; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) histo_eta[i][j][k][l] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<2; k++) histo_etaintegral[i][j][k] = 0;
   histo_trackfreq_vs_time = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) histo_eta_vs_Q[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) histo_hitocc_saturated[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) histo_clusters_average[i][j] = 0;
	histo_scatter_autofidcut = 0;
}

void Clustering::DeleteHistograms() {
   for(int i=0; i<2; i++) delete histo_dianoise[i];
   for(int i=0; i<2; i++) delete histofit_dianoise[i];
   for(int i=0; i<5; i++) for(int j=0; j<3; j++) delete histo_scatter[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) delete histo_clusterocc[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) delete histo_landau[i][j][k][l];
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) delete histo_clusterfreq[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) delete histo_clustersizefreq[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<6; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) delete histo_eta[i][j][k][l];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<2; k++) delete histo_etaintegral[i][j][k];
   delete histo_trackfreq_vs_time;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) delete histo_eta_vs_Q[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) delete histo_hitocc_saturated[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) delete histo_clusters_average[i][j];
	delete histo_scatter_autofidcut;
}

void Clustering::DrawHistograms() {
   
   string tempname = "tempname";
   for(int cut=0; cut<2; cut++) if(histoswitch_dianoise[cut]) {
      tempname += "0";
      histofit_dianoise[cut] = new TF1(tempname.c_str(),"gaus",histo_dianoise[cut]->GetMean()-histo_dianoise[cut]->GetRMS(),histo_dianoise[cut]->GetMean()+histo_dianoise[cut]->GetRMS());
      histofit_dianoise[cut]->SetLineColor(kBlue);
      histo_dianoise[cut]->Fit(tempname.c_str(),"r"); // fit option "r" restricts the range of the fit
      SaveHistogramPNG(histo_dianoise[cut]);
   }
   
   for(int det=0; det<5; det++)
      for(int trackcut=0; trackcut<3; trackcut++)
         if(histoswitch_scatter[det][trackcut]) SaveHistogramPNG(histo_scatter[det][trackcut]);
   
   gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
   for(int det=0; det<9; det++) {
      
      for(int frac=0; frac<2; frac++)
         if(histoswitch_hitocc_saturated[det][frac]) SaveHistogramPNG(histo_hitocc_saturated[det][frac]);
      
      for(int goodcut=0; goodcut<2; goodcut++)
         if(histoswitch_clusterfreq[det][goodcut]) SaveHistogramPNG(histo_clusterfreq[det][goodcut]);
         
      histo_clusters_average[det][0]->SetNormFactor(histo_clusters_average[det][0]->Integral(-100,100)/singlesitrack_events);
      
      for(int chip=0; chip<2; chip++) {
         for(int etatype=0; etatype<2; etatype++) {
            if(histoswitch_eta_vs_Q[det][etatype][chip]) SaveHistogramPNG(histo_eta_vs_Q[det][etatype][chip]);
         }
         
         //agregate eta plots
         TCanvas etacanvas("etacanvas");
         if(histoswitch_etaintegral[det][2][chip]) histo_etaintegral[det][2][chip]->Draw();
         if(histoswitch_etaintegral[det][0][chip]) histo_etaintegral[det][0][chip]->Draw("same");
         if(histoswitch_etaintegral[det][1][chip]) histo_etaintegral[det][1][chip]->Draw("same");
         pt->Draw();
         ostringstream plot_filename;
         if(chip) plot_filename << plots_path << "EtaIntegrals_D" << det << "_RightChip.png";
         else plot_filename << plots_path << "EtaIntegrals_D" << det << "_LeftChip.png";
         etacanvas.Print(plot_filename.str().c_str());
         
         //individual eta plots_canvas
         if(histo_etaintegral[det][0][chip]) SaveHistogramPNG(histo_etaintegral[det][0][chip]);
         if(histo_etaintegral[det][2][chip]) SaveHistogramPNG(histo_etaintegral[det][2][chip]);
            
      }
      
      for(int trackcut=0; trackcut<3; trackcut++) {

         if(histoswitch_clusters_average[det][trackcut]) SaveHistogramPNG(histo_clusters_average[det][trackcut]);
         
         for(int clussizeopt=0; clussizeopt<3; clussizeopt++)
            if(histoswitch_clustersizefreq[det][clussizeopt][trackcut]) SaveHistogramPNG(histo_clustersizefreq[det][clussizeopt][trackcut]);
         
         for(int etaopt=0; etaopt<6; etaopt++) {
            if(det==8) if(histoswitch_eta[det][etaopt][trackcut][0]) SaveHistogramPNG(histo_eta[det][etaopt][trackcut][0]);
            if(det<8) for(int chip=0; chip<2; chip++)
               if(histoswitch_eta[det][etaopt][trackcut][chip]) SaveHistogramPNG(histo_eta[det][etaopt][trackcut][chip]);
         }
         
         for(int nhitsopt=0; nhitsopt<11; nhitsopt++) {
            //hitoccs
            if(histoswitch_clusterocc[det][nhitsopt][trackcut]) 
               SaveHistogramPNG(histo_clusterocc[det][nhitsopt][trackcut]);
            
            //pulse heights
            for(int snr=0; snr<2; snr++) 
               if(histoswitch_landau[det][nhitsopt][trackcut][snr]) SaveHistogram(histo_landau[det][nhitsopt][trackcut][snr]);
         }
      }
   }
}


void Clustering::BookHistograms() {
   
   TDetectorPlane D0, D1, D2, D3, Dia;
   D0.SetZ(detectorD0Z);
   D1.SetZ(detectorD1Z);
   D2.SetZ(detectorD2Z);
   D3.SetZ(detectorD3Z);
   Dia.SetZ(detectorDiaZ);
   bool mask;
   
   //cut flow counters
   //Int_t total_events, goldengatecluster_events, badchannelcluster_events, singlesitrack_events, singlesitrack_1diamondclus_events, singlesitrack_fidcut_events, singlesitrack_fidcut_1diamondclus_events;
   
   //count total events
   total_events++;
   
   //count bad events
   //if(clustered_event.HasSaturatedCluster()) saturatedcluster_events[9]++;
   
   
   //cut events containing the following clusters
   //silicon: saturated, lumpy, goldengate
   //diamond: saturated,
   if(CMNEvent_flag || ZeroDivisorEvent_flag || clustered_event.HasSaturatedCluster() || clustered_event.HasSaturatedCluster(8) || clustered_event.HasLumpyCluster() || clustered_event.HasGoldenGateCluster() || clustered_event.HasBadChannelCluster()) {//|| !clustered_event.HasOneSiliconTrack()) {//clustered_event.HasBadChannelCluster()) {
      if(CMNEvent_flag) CMNEvents++;
      if(ZeroDivisorEvent_flag) ZeroDivisorEvents++;
      
      //count cut types of events here
      
      return;
   }
   
   //count various types of events surviving the above cuts
   if(clustered_event.HasGoldenGateCluster()) goldengatecluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasGoldenGateCluster(det)) goldengatecluster_events[det]++;
   if(clustered_event.HasBadChannelCluster()) badchannelcluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasBadChannelCluster(det)) badchannelcluster_events[det]++;
   if(clustered_event.HasLumpyCluster()) lumpycluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasLumpyCluster(det)) lumpycluster_events[det]++;
   if(clustered_event.HasSaturatedCluster()) saturatedcluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasSaturatedCluster(det)) {
      saturatedcluster_events[det]++;
      for(uint clus=0; clus<clustered_event.GetNCutClusters(det); clus++)
         for(int frac=0; frac<2; frac++)
            if(clustered_event.GetCutCluster(det,clus)->IsSaturatedCluster()) for(int hit=0; hit<clustered_event.GetCutCluster(det,clus)->GetNHits(); hit++) {
               if(det<8) if(clustered_event.GetCutCluster(det,clus)->GetADC(hit)>4094 && Det_channel_screen[det].CheckChannel(clustered_event.GetCutCluster(det,clus)->GetChannel(hit))) histo_hitocc_saturated[det][frac]->Fill(clustered_event.GetCutCluster(det,clus)->GetChannel(hit));
               if(det==8) if(clustered_event.GetCutCluster(det,clus)->GetADC(hit)>4094 && Det_channel_screen[det].CheckChannel(clustered_event.GetCutCluster(det,clus)->GetChannel(hit)) && clustered_event.GetCutCluster(det,clus)->GetChannel(hit)>63) histo_hitocc_saturated[det][frac]->Fill(clustered_event.GetCutCluster(det,clus)->GetChannel(hit)-63);
            }
   }
   totalsurviving_events++;
   
   Float_t centroid, charge, snr, stddev, eta, eta2ch = -1;
   Int_t nhits, nseeds, clus_peak;
   Cluster* current_cluster = 0;
   
   //simple tracking: require one and only one cluster in each silicon plane
   bool one_and_only_one = clustered_event.HasOneSiliconTrack();
   if(one_and_only_one) singlesitrack_events++;
   //else if((current_event-1)%1000==0) clustered_event.Print();
   
   //look to see if there are coincident hits in the silicon planes
   for(int d=0; d<4; d++)
      if(clustered_event.HasCoincidentClustersXY(d)) detectorxycluster_events[d]++;
   
   //if a track, check whether it's in the silicon fiducial region
   bool fiducial_track = 0;
   Float_t si_avg_x=0, si_avg_y=0, det_centroid[8];
   if(one_and_only_one) {
      for(int det=0; det<4; det++) {
         si_avg_x += clustered_event.GetCluster(2*det,0)->Get1stMoment();
         si_avg_y += clustered_event.GetCluster(2*det+1,0)->Get1stMoment();
      }
      si_avg_x = si_avg_x/4;
      si_avg_y = si_avg_y/4;
      
      if(si_avg_x>si_avg_fidcut_xlow && si_avg_x<si_avg_fidcut_xhigh && si_avg_y>si_avg_fidcut_ylow && si_avg_y<si_avg_fidcut_yhigh) {
         fiducial_track=1;
         singlesitrack_fidcut_events++;
      }
      
      //silicon track clusters' centroids
      for(int det=0; det<8; det++)
         det_centroid[det] = clustered_event.GetCluster(det,0)->Get1stMoment();
      
      //silicon track scatter for each plane
      for(int det=0; det<4; det++)
         histo_scatter[det][0]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
      //silicon track average scatter
      histo_scatter[4][0]->Fill(si_avg_x,si_avg_y);
      
      //scatter of tracks in the diamond
      if(clustered_event.GetNClusters(8)==1) {
         //scatter of tracks in diamond in each plane
         for(int det=0; det<4; det++)
            histo_scatter[det][1]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
         //scatter of tracks in diamond averaged over silicon planes
         histo_scatter[4][1]->Fill(si_avg_x,si_avg_y);
		  //check usual events
//		  if (si_avg_x < 80) {
//			  EventMonitor(current_event-1);
//		  }
         //count diamond track events
         singlesitrack_1diamondclus_events++;
         //count fidcut diamond track events
         if(fiducial_track) singlesitrack_fidcut_1diamondclus_events++;
      }
      
      //scatter of diamond tracks in the silicon fiducial region
      if(fiducial_track && clustered_event.GetNClusters(8)==1) {
         //scatter of diamond tracks in the silicon fiducial region
         for(int det=0; det<4; det++) 
            histo_scatter[det][2]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
         //scatter of diamond tracks in the silicon fiducial region averaged over silicon planes
         histo_scatter[4][2]->Fill(si_avg_x,si_avg_y);
      }
      
      //diamond noise plots
      float psadc;
      for(uint i=0; i<Det_NChannels[8]; i++) {
         psadc=Dia_ADC[i]-Det_PedMean[8][i];
         if(Det_channel_screen[8].CheckChannel(Det_Channels[8][i]) && psadc<Di_Cluster_Hit_Factor*Det_PedWidth[8][i]) {
            histo_dianoise[0]->Fill(psadc);
            if(fiducial_track) histo_dianoise[1]->Fill(psadc);
         }
      }
      
      //save tracks for alignment
      
      bool zero_suppressed_track = 0;
      for(int det=0; det<8; det++)
         if(clustered_event.GetCluster(det,0)->highest2_centroid==-1)
            zero_suppressed_track=1;
         
      if(zero_suppressed_track) //skip events where pedestal calc zero-suppression prevents eta correction
         counter_alignment_tracks_zero_suppressed++;
      
      else { //if track is ok, then add to list of alignment tracks
         
         //now tracking with centroid of the seed and the highest neighbor (for eta correction)
         D0.SetX(clustered_event.GetCluster(0,0)->highest2_centroid);
         D1.SetX(clustered_event.GetCluster(2,0)->highest2_centroid);
         D2.SetX(clustered_event.GetCluster(4,0)->highest2_centroid);
         D3.SetX(clustered_event.GetCluster(6,0)->highest2_centroid);
         
         D0.SetY(clustered_event.GetCluster(1,0)->highest2_centroid);
         D1.SetY(clustered_event.GetCluster(3,0)->highest2_centroid);
         D2.SetY(clustered_event.GetCluster(5,0)->highest2_centroid);
         D3.SetY(clustered_event.GetCluster(7,0)->highest2_centroid);

         //now tracks in diamond
         if(clustered_event.GetNClusters(8)==1 && fiducial_track) {
            if(clustered_event.GetCluster(8,0)->highest2_centroid!=-1) {
               Dia.SetX(clustered_event.GetCluster(8,0)->highest2_centroid);
               TDiamondTrack track_fidcut = TDiamondTrack(clustered_event.event_number,D0,D1,D2,D3,Dia);
               tracks_fidcut.push_back(track_fidcut);
               counter_alignment_fidcut_tracks++;
               mask = rand.Uniform()<alignment_training_track_fraction;
               tracks_fidcut_mask.push_back(mask);
               counter_alignment_only_fidcut_tracks += mask;
            }
         }
         else
            Dia.SetX(-1);
         
         TDiamondTrack track = TDiamondTrack(clustered_event.event_number,D0,D1,D2,D3);
         //track.SetEventNumber(Silicon_tracks[t].Event_number);
         tracks.push_back(track);
         counter_alignment_tracks++;
         mask = rand.Uniform()<alignment_training_track_fraction;
         tracks_mask.push_back(mask);
         counter_alignment_only_tracks += mask;
         
      }//end else
      
   }//end if one_and_only_one
   
   uint nclusters;
   //loop over detectors
   for(int det=0; det<9; det++) {
      
      //cluster frequency
      nclusters=clustered_event.GetNClusters(det);
      histo_clusterfreq[det][0]->Fill(nclusters);
      histo_clusterfreq[det][1]->Fill(clustered_event.GetNCutClusters(det));
      
      //loop over all *good* clusters
      for(uint clus=0; clus<nclusters; clus++) {
         
         current_cluster = clustered_event.GetCluster(det,clus);
         
         //if(current_cluster->IsBadChannelCluster()) continue; // skip bad clusters
         
         nhits = current_cluster->GetNHits();
         nseeds = current_cluster->GetNSeeds();
         stddev = current_cluster->GetStandardDeviation();
         centroid = current_cluster->Get1stMoment();
         charge = current_cluster->GetCharge();
         snr = current_cluster->GetSNR();
         eta = current_cluster->GetEta();
         //if(totalsurviving_events%1000==0) cout<<"Event "<<current_event-1<<": Transparent eta for detector "<<det<<" is "<<eta<<endl;
         if(nhits==2) eta2ch = (current_cluster->GetADC(1)-current_cluster->GetPedMean(1))/current_cluster->GetCharge();
         clus_peak = current_cluster->GetPeakHit();
         
         //all good clusters
         //------------
         
         //cluster size frequency
         histo_clustersizefreq[det][0][2]->Fill(nhits);
         histo_clustersizefreq[det][1][2]->Fill(nseeds);
         histo_clustersizefreq[det][2][2]->Fill(stddev);
         
         //transparent eta
         if(eta!=-1) {
            if(centroid<128) histo_eta[det][0][2][0]->Fill(eta); //first chip
            else histo_eta[det][0][2][1]->Fill(eta); //second chip
         }
         
         //cluster property frequency plots
         if(nhits==1) {
            histo_clusterocc[det][0][2]->Fill(centroid); // fill 1 cluster histo
            histo_landau[det][0][2][0]->Fill(charge);
            histo_landau[det][0][2][1]->Fill(snr);
            histo_clusterocc[det][5][2]->Fill(centroid); // fill 1&2 cluster histo
            histo_landau[det][5][2][0]->Fill(charge);
            histo_landau[det][5][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==2) {
            if(centroid<128) histo_eta[det][1][2][0]->Fill(eta2ch); // fill 2hit cluster eta for first chip
            else histo_eta[det][1][2][1]->Fill(eta2ch); // fill 2hit cluster eta for second chip
            histo_clusterocc[det][1][2]->Fill(centroid); // fill 2 cluster histo
            histo_landau[det][1][2][0]->Fill(charge);
            histo_landau[det][1][2][1]->Fill(snr);
            histo_clusterocc[det][5][2]->Fill(centroid); // fill 1&2 cluster histo
            histo_landau[det][5][2][0]->Fill(charge);
            histo_landau[det][5][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][8][2]->Fill(centroid); // fill 2&3 cluster histo
            histo_landau[det][8][2][0]->Fill(charge);
            histo_landau[det][8][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==3) {
            histo_clusterocc[det][2][2]->Fill(centroid); // fill 3 cluster histo
            histo_landau[det][2][2][0]->Fill(charge);
            histo_landau[det][2][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][8][2]->Fill(centroid); // fill 2&3 cluster histo
            histo_landau[det][8][2][0]->Fill(charge);
            histo_landau[det][8][2][1]->Fill(snr);
            histo_clusterocc[det][9][2]->Fill(centroid); // fill 3&4 cluster histo
            histo_landau[det][9][2][0]->Fill(charge);
            histo_landau[det][9][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==4) {
            histo_clusterocc[det][3][2]->Fill(centroid); // fill 4 cluster histo
            histo_landau[det][3][2][0]->Fill(charge);
            histo_landau[det][3][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][9][2]->Fill(centroid); // fill 3&4 cluster histo
            histo_landau[det][9][2][0]->Fill(charge);
            histo_landau[det][9][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits>4) {
            histo_clusterocc[det][4][2]->Fill(centroid); // fill >=5 cluster histo
            histo_landau[det][4][2][0]->Fill(charge);
            histo_landau[det][4][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         
         //cut on clean tracks
         //-------------------
         
         if(one_and_only_one) {
            
            //cluster size frequency
            histo_clustersizefreq[det][0][0]->Fill(nhits);
            histo_clustersizefreq[det][1][0]->Fill(nseeds);
            histo_clustersizefreq[det][2][0]->Fill(stddev);
            
            //transparent eta
            if(eta!=-1) {
               if(centroid<128) histo_eta[det][0][0][0]->Fill(eta); //first chip
               else histo_eta[det][0][0][1]->Fill(eta); //second chip
            }
            
            //stack clusters
            for(int h=0; h<nhits; h++)
               histo_clusters_average[det][0]->Fill(h-clus_peak,current_cluster->GetPSADC(h));
            
            //cluster property frequency plots
            if(nhits==1) {
               histo_clusterocc[det][0][0]->Fill(centroid); // fill 1 cluster histo
               histo_landau[det][0][0][0]->Fill(charge);
               histo_landau[det][0][0][1]->Fill(snr);
               histo_clusterocc[det][5][0]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][0][0]->Fill(charge);
               histo_landau[det][5][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==2) {
               if(centroid<128) histo_eta[det][1][0][0]->Fill(eta2ch); // fill 2hit cluster eta
               else histo_eta[det][1][0][1]->Fill(eta2ch); // fill 2hit cluster eta
               histo_clusterocc[det][1][0]->Fill(centroid); // fill 2 cluster histo
               histo_landau[det][1][0][0]->Fill(charge);
               histo_landau[det][1][0][1]->Fill(snr);
               histo_clusterocc[det][5][0]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][0][0]->Fill(charge);
               histo_landau[det][5][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][8][0]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][0][0]->Fill(charge);
               histo_landau[det][8][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==3) {
               histo_clusterocc[det][2][0]->Fill(centroid); // fill 3 cluster histo
               histo_landau[det][2][0][0]->Fill(charge);
               histo_landau[det][2][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][8][0]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][0][0]->Fill(charge);
               histo_landau[det][8][0][1]->Fill(snr);
               histo_clusterocc[det][9][0]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][0][0]->Fill(charge);
               histo_landau[det][9][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==4) {
               histo_clusterocc[det][3][0]->Fill(centroid); // fill 4 cluster histo
               histo_landau[det][3][0][0]->Fill(charge);
               histo_landau[det][3][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][9][0]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][0][0]->Fill(charge);
               histo_landau[det][9][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits>4) {
               histo_clusterocc[det][4][0]->Fill(centroid); // fill >=5 cluster histo
               histo_landau[det][4][0][0]->Fill(charge);
               histo_landau[det][4][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
         }
         
         //cut fiducial tracks
         //-------------------
         
         if(fiducial_track && (det<8 || (det==8 && nclusters==1))) {
            
            //cluster size frequency
            histo_clustersizefreq[det][0][1]->Fill(nhits);
            histo_clustersizefreq[det][1][1]->Fill(nseeds);
            histo_clustersizefreq[det][2][1]->Fill(stddev);
            
            //transparent eta
            if(eta!=-1) {
               if(centroid<128) {
                  histo_eta[det][0][1][0]->Fill(eta);
                  histo_eta_vs_Q[det][0][0]->Fill(charge,eta); //fill eta vs charge plot
               }
               else {
                  histo_eta[det][0][1][1]->Fill(eta);
                  histo_eta_vs_Q[det][0][1]->Fill(charge,eta); //fill eta vs charge plot
               }
            }
            
            
            //cluster property frequency plots
            if(nhits==1) {
               histo_clusterocc[det][0][1]->Fill(centroid); // fill 1 cluster histo
               histo_landau[det][0][1][0]->Fill(charge);
               //histo_landau[det][0][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][0][1][1]->Fill(snr);
               histo_clusterocc[det][5][1]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][1][0]->Fill(charge);
               //histo_landau[det][5][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][5][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==2) {
               if(centroid<128) {
                  histo_eta[det][1][1][0]->Fill(eta2ch); // fill 2hit cluster eta
                  histo_eta_vs_Q[det][1][0]->Fill(charge,eta2ch); //fill eta vs charge plot
               }
               else {
                  histo_eta[det][1][1][1]->Fill(eta2ch); // fill 2hit cluster eta
                  histo_eta_vs_Q[det][1][1]->Fill(charge,eta2ch); //fill eta vs charge plot
               }
               histo_clusterocc[det][1][1]->Fill(centroid); // fill 2 cluster histo
               histo_landau[det][1][1][0]->Fill(charge);
               histo_landau[det][1][1][1]->Fill(snr);
               histo_clusterocc[det][5][1]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][1][0]->Fill(charge);
               //histo_landau[det][5][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][5][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][8][1]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][1][0]->Fill(charge);
               histo_landau[det][8][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==3) {
               histo_clusterocc[det][2][1]->Fill(centroid); // fill 3 cluster histo
               histo_landau[det][2][1][0]->Fill(charge);
               histo_landau[det][2][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][8][1]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][1][0]->Fill(charge);
               histo_landau[det][8][1][1]->Fill(snr);
               histo_clusterocc[det][9][1]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][1][0]->Fill(charge);
               histo_landau[det][9][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==4) {
               histo_clusterocc[det][3][1]->Fill(centroid); // fill 4 cluster histo
               histo_landau[det][3][1][0]->Fill(charge);
               histo_landau[det][3][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][9][1]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][1][0]->Fill(charge);
               histo_landau[det][9][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits>4) {
               histo_clusterocc[det][4][1]->Fill(centroid); // fill >=5 cluster histo
               histo_landau[det][4][1][0]->Fill(charge);
               histo_landau[det][4][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
         }
      }//end loop over clusters
   }//end loop over detectors
}

void Clustering::GenerateHTML() {
   
   int section, subsection;
   
   //summary page
   string html_summary_path = plots_path + "/index.html";
   ofstream html_summary(html_summary_path.c_str());
   cout << "Run summary HTML created at " << html_summary_path << endl;
   
   html_summary << "<html>" << endl; 
   html_summary << "<title>Run "<<run_number<<" analysis results - Summary</title>" << endl;
   html_summary << "<body bgcolor=\"White\">" << endl;
   html_summary << "<font face=\"Arial,Sans\" size=2>"<<endl;
   html_summary << "<a name=\"top\"></a>" << endl;
   html_summary << "<b>Summary</b>||"
                << "<a href=\"d8.html\">Diamond</a>||"
                << "<a href=\"d0.html\">D0X</a>||"
                << "<a href=\"d1.html\">D0Y</a>||"
                << "<a href=\"d2.html\">D1X</a>||"
                << "<a href=\"d3.html\">D1Y</a>||"
                << "<a href=\"d4.html\">D2X</a>||"
                << "<a href=\"d5.html\">D2Y</a>||"
                << "<a href=\"d6.html\">D3X</a>||"
                << "<a href=\"d7.html\">D3Y</a><p>"<<endl;
   html_summary << "<center><font color=\"#000000\"><h1>Run "<<run_number<<" analysis results - Summary</h1></font></center>"<<endl;
   html_summary << "<hr size=\"10\" Color=\"#ffff33\">" << endl; 
   html_summary << "Results from " << dateandtime.GetMonth() << "/ " << dateandtime.GetDay() << "/" << dateandtime.GetYear() << " at " << dateandtime.GetHour() << ":" << dateandtime.GetMinute() << ":" << dateandtime.GetSecond() << ")<p>" << endl;
   
   html_summary<<"<h3>Settings for pedestal calculation</h3>"<<endl;
   if(dia_input==0) html_summary<<"Dia0 or sirocco input 4 selected for analysis<br>"<<endl;
   if(dia_input==1) html_summary<<"Dia1 or sirocco input 5 selected for analysis<br>"<<endl;
   html_summary<<"Buffer size of "<<Iter_Size<<" events used for pedestal analysis.<br>"<<endl;
   if(fix_dia_noise>=0) html_summary<<"Noise for all diamond channels FIXED at "<<fix_dia_noise<<"<br>"<<endl;
   html_summary<<"Silicon channel adc and pedestal info saved for psadc fluctuations >"<<store_threshold<<"sigma<br>"<<endl;
   html_summary<<"For pedestal analysis, a silicon hit is >"<< Si_Pedestal_Hit_Factor <<"sigma and a diamond hit is >"<< Di_Pedestal_Hit_Factor <<"sigma<br>"<<endl;
   html_summary<<"Events with diamond common mode noise fluctations >"<<CMN_cut<<"sigma flagged for removal from subsequent clustering analysis<br>"<<endl;
   if(Taylor_speed_throttle==1) html_summary<<"Taylor's RMS speed tweak disabled<br>"<<endl;
   else html_summary<<"Taylor's RMS speed tweak enabled; RMS recalculated normally every "<<Taylor_speed_throttle<<" events<br>"<<endl;
   
   html_summary<<"<h3>Settings for clustering</h3>"<<endl;
   html_summary<<"For silicon cluster analysis, a hit is >"<< Si_Cluster_Hit_Factor <<"sigma and a seed is >"<< Si_Cluster_Seed_Factor <<"sigma<br>"<<endl;
   html_summary<<"For diamond cluster analysis, a hit is >"<< Di_Cluster_Hit_Factor <<"sigma and a seed is >"<< Di_Cluster_Seed_Factor <<"sigma<br>"<<endl;
   html_summary<<"The following diamond channels were screened in the clustering analysis: ";
   for(int i=0; i<128; i++) {
      if(!Det_channel_screen[8].CheckChannel(i)) html_summary<<i<<",";
   }
   html_summary<<"<br>"<<endl;
   
   html_summary<<"<h3>Cut flow</h3>"<<endl;
   //html_summary<<"<ul>"<<endl;//1
   html_summary<<"<b>Events slated for analysis = "<<total_events<<"</b>"<<endl;
   html_summary<<"<ul>"<<endl;//2
   html_summary<<"<li><b>Events cut from analysis= "<<total_events-totalsurviving_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//3a
   html_summary<<"<li>CMNEvents = "<<CMNEvents<<"</li>"<<endl;
   html_summary<<"<li>ZeroDivisorEvents = "<<ZeroDivisorEvents<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//3a
   html_summary<<"<li><b>Events surviving above cuts analyzed = "<<totalsurviving_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//3b
   html_summary<<"<li>Golden gate cluster events = "<<goldengatecluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4a
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" golden gate cluster events = "<<goldengatecluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia golden gate cluster events = "<<goldengatecluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4a
   html_summary<<"<li>Bad channel cluster events = "<<badchannelcluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4b
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" bad channel cluster events = "<<badchannelcluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia bad channel cluster events = "<<badchannelcluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4b
   html_summary<<"<li>Lumpy cluster events = "<<lumpycluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4c
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" lumpy cluster events = "<<lumpycluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia lumpy cluster events = "<<lumpycluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4c
   html_summary<<"<li>Saturated cluster events = "<<saturatedcluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4d
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" saturated cluster events = "<<saturatedcluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia saturated cluster events = "<<saturatedcluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4d
   for(int d=0; d<4; d++)
      html_summary<<"<li>detectorxycluster_events["<<d<<"] = "<<detectorxycluster_events[d]<<"</li>"<<endl;
   html_summary<<"<li><b>Events with a simple <u>track</u> or <i>one and only one cluster in each silicon</i> = "<<singlesitrack_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//4c
   html_summary<<"<li>Events with single track in silicon and single cluster in diamond = "<<singlesitrack_1diamondclus_events<<"</li>"<<endl;
   html_summary<<"<li><b>Events with single track in silicon fiducial region = "<<singlesitrack_fidcut_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//5
   html_summary<<"<li><b>Events with single track in silicon fiducial region and single cluster in diamond = "<<singlesitrack_fidcut_1diamondclus_events<<"</b></li>"<<endl;
   if(singlesitrack_fidcut_1diamondclus_events!=histo_clusterocc[8][10][1]->GetEntries()) cout<<"Clustering::GenerateHTML: We've got a problem here. "<<singlesitrack_fidcut_1diamondclus_events<<" is not equal to "<<histo_clusterocc[8][10][1]->GetEntries()<<endl;
   html_summary<<"<ul>"<<endl;//6
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 1-hit cluster in diamond = "<<(float)histo_clusterocc[8][0][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 2-hit cluster in diamond = "<<(float)histo_clusterocc[8][1][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 3-hit cluster in diamond = "<<(float)histo_clusterocc[8][2][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 4-hit cluster in diamond = "<<(float)histo_clusterocc[8][3][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single >4-hit cluster in diamond = "<<(float)histo_clusterocc[8][4][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//6
   html_summary<<"</ul>"<<endl;//5
   html_summary<<"</ul>"<<endl;//4c
   html_summary<<"</ul>"<<endl;//3b
   html_summary<<"</ul>"<<endl;//2
   //html_summary<<"</ul>"<<endl;//1
   html_summary<<"<hr size=\"5\">" << endl;
   
   
   // Plot index
   section=0;
   subsection=1;
   html_summary << "<a name=\"plot_index\"></a>" << endl;
   //html_summary << "<font color=\"f00000\"><a href=\"silicon.html\">For silicon plots click here.</a></font><p>" << endl;
   html_summary << "<p><h2>Plot Index</h2></p>" << endl;
   html_summary << "<p><a href=\"#Pedestal_Data_jump\">("<<subsection++<<") Diamond pedestal plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Hitocc_Data_jump\">("<<subsection++<<") Diamond track occupancy plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Scatter_Data_jump\">("<<subsection++<<") Silicon scatter plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Pulse_Height_jump\">("<<subsection++<<") Diamond Pulse Height Landau Plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Pulse_Height_fid_jump\">("<<subsection++<<") Diamond Pulse Height Landau Plots with Silicon Fiducial Cut </a></p>" << endl;
   html_summary << "<p><a href=\"#Charge_Dist_jump\">("<<subsection++<<") Adjacent Channel Cluster Asymmetry Plots </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   
   // Diamond pedestal plots 
   section++;
   subsection=1;
   html_summary << "<a name=\"Pedestal_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Diamond pedestal plots</h1>" << endl;
   html_summary << "<p><h3>("<<section<<"."<<subsection++<<") Diamond Pedestal Mean</h3> Diamond pedestal means by Channel</p>" << endl;
   html_summary << "<p><img SRC=\"Pedestal_Values.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<") Diamond Pedestal RMS (Initial)</h3> A diamond pedestal RMS by channel for initial values" << endl;
   html_summary << "<p><img SRC=\"Channel_Noise_Initial.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<") Diamond Pedestal RMS (Final)</h3> A diamond pedestal RMS by channel for final calculation (from running calculation)" << endl;
   html_summary << "<p><img SRC=\"Channel_Noise_Final.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the raw ADC values verses event number (a zoomed in section)." << endl;
   html_summary << "<p><img SRC=\"Raw_ADC_vs_Event_zoom.png\" BORDER=\"1\"></p>";
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the raw ADC values versus event number with CMN events removed." << endl;
   html_summary << "<p><img SRC=\"Raw_ADC_vs_Event_CMN_cut_zoom.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the pedestal subtracted ADC values vesus event number (zoom)" << endl;
   html_summary << "<p><img SRC=\"PS_ADC_vs_Event_zoom.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   
   // Diamond track occupancy plots
   section++;
   subsection=1;
   html_summary << "<a name=\"Hitocc_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Diamond occupancy plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the number of raw hits (>"<<Di_Pedestal_Hit_Factor<<"sigma) in each channel that the pedestal analysis found." << endl;
   html_summary << "<p><img SRC=\"hit_occup_can_dia.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all events in the diamond." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all single cluster events in the diamond corresponding to <b>tracks</b> in the silicon." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all single cluster events in the diamond corresponding to tracks in the silicon <b>fiducial</b> region." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   // Silicon scatter plots
   section++;
   subsection=1;
   html_summary << "<a name=\"Scatter_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Silicon scatter plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A Scatter Plot of the Average X values for the 1 and 2 silicon planes versus an average of the y values for the 1 and 2 silicon planes for events that go through all 8 silicon planes" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the above but restricted to events where the event had a single diamond cluster." << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the above but restricted to events where the event had a single diamond cluster and is in the fiducial cut on the diamond.  A cut on the average value of the X and Y channels of the silicon is made for tracks through all 8 planes and the diamond.</p>" << endl;
   html_summary << "<p>The range in silicon space is: Channels " << si_avg_fidcut_xlow << " to " << si_avg_fidcut_xhigh << " on the x-axis and " << si_avg_fidcut_ylow << " to " << si_avg_fidcut_yhigh << " on the y-axis</p>" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
    
   // Full diamond landaus
   section++;
   subsection=1;
   html_summary << "<a name=\"Pulse_Height_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Full Diamond Pulse Height Landau Plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of diamond cluster size or number of hits (>"<<Di_Cluster_Hit_Factor<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of number of seeds (>"<<Di_Cluster_Seed_Factor<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>all channel size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][10][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][0][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][1][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][2][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>4 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][3][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>>4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][4][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>both 1 and 2 hit clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][5][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, and 3 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][6][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, 3 and 4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][7][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 and 3 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][8][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 and 4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][9][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
     
   // Fidcut diamond landaus
   section++;
   subsection=1;
   html_summary << "<a name=\"Pulse_Height_fid_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Fiducial Cut Pulse Height and Noise</h1>" << endl;
   html_summary << "<p>The following plots are taken from a cut in the silicon planes. A cut on the average value of the X and Y channels of the silicon is made for tracks through all 8 planes and the diamond.</p>" << endl;
   html_summary << "<p>The range in silicon space is: Channels " << si_avg_fidcut_xlow << " to " << si_avg_fidcut_xhigh << " on the x-axis and " << si_avg_fidcut_ylow << " to " << si_avg_fidcut_yhigh << " on the y-axis</p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<p><h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the X and Y channels of the silicon surrounding the diamond with the fiducial cut.</p>" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Projection of the silicon fiducial region in the diamond." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][1]->GetTitle()<<".png\" BORDER=\"1\"></p>";
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of diamond cluster size or number of hits (>"<<Di_Cluster_Hit_Factor<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][0][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of number of seeds (>"<<Di_Cluster_Seed_Factor<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][1][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>all channel size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][10][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][0][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][1][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][2][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>4 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][3][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>>4 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][4][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>both 1 and 2 hit clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][5][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, and 3 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][6][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, 3 and 4 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][7][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 and 3 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][8][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 and 4 hit size clusters for the fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][9][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels when there was a track in the silicon fiducial region.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
      
   // Charge sharing plots
   section++;
   subsection=1;
   html_summary << "<a name=\"Charge_Dist_jump\"></a>" << endl;
   html_summary << "<p><h1>("<<section<<") Adjacent Channel Cluster Asymmetry Analysis</h1></p>" << endl;
   html_summary << "<p>This section shows the left to right distribution of charge for one, two, and three channel clusters and for the transparency hits. The main analysis is looking at the distribution of charge for two channel clusters, that is, which channel tends to have more charge deposited than the other. When ordered 0-128, the channels to the left are the lower number indexed channels and the higher are the right. The plot is made by taking the amount of charge in the left channel and dividing by the total. The histogram fills in bins less than 0.5 for left-biased clusters and greater than 0.5 for right-biased clusters. For one channel clusters, the two channels used are the actual cluster, and then the next highest adjacent channel. For three channel clusters, the two largest adjacent charge values are used. For the transparency, any hit above 5 sigma is used as a seed and the next highest adjacent channel is also take.</p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Two Channel Cluster Distribution" << endl;
   html_summary << "<p><img SRC=\""<<histo_eta[8][1][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Transparency Distribution" << endl;
   html_summary << "<p><img SRC=\""<<histo_eta[8][0][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   /*
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> One Channel Cluster and Next Highest Channel Distribution" << endl;
   html_summary << "<p><img SRC=\"One_Channel_Cluster_and_next_highest_channel_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster Two Highest Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_two_highest_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster High to Low Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_hitolow_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster Two Lowest Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_twolow_Distribution.png\" BORDER=\"1\"></p>" << endl;
   */
   html_summary << "<hr size=\"5\">" << endl;

   html_summary << "</font>"<<endl;
   html_summary << "</HTML>" << endl;
   html_summary.close();
   
}



/******************88
Note that the following ClusterRun is a flexible script;
for example, we could draw histograms and generate html
for the first 100k events, then change plot directory
and plot the second 100k events, 
********************/

//sequentially cluster all events
void Clustering::ClusterRun(bool plots) {
   
   //Start Timer
   TStopwatch Watch;
   Watch.Start();
   
   //record large clusters
   string large_clusters_filename = plots_path + "large_clusters.txt";
   ofstream large_clusters_file(large_clusters_filename.c_str());
   large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
   uint nclusters, nhits;
   Cluster* current_cluster;
   
   //histograms of hit positions to find out if the eta correction helped
   TH1F* histo_hitpos[9][2];
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         ostringstream histotitle;
         histotitle << "Hit_Position_Interstrip_D" << det << "_chip" << chip;
         histo_hitpos[det][chip] = new TH1F(histotitle.str().c_str(), histotitle.str().c_str(), 101, -0.005, 1.005);
      }
   }
	
	
	current_event = 0;
	if (false) { // TODO: remove this section (moved to AutoFidCut())
		// produce necessary plots to detect fidcut region
		cout << endl << endl << "-- produce plot for AutoFidCut().." << endl;
		cout << "running over " << PedTree->GetEntries() << " events.." << endl;
		for (uint e=0; e<PedTree->GetEntries(); e++) { // maybe it's not necessary to run over all events for the AutoFidCut?!
			if (!AlternativeClustering) ClusterEvent();
			else ClusterEventSeeds();
			if (e%10000==0) clustered_event.Print();
			
//			current_cluster = 0;
			
			// -- produce scatter plot for AutoFidCut
			bool one_and_only_one = clustered_event.HasOneSiliconTrack();
			Float_t si_avg_x=0, si_avg_y=0;
			if (CMNEvent_flag || ZeroDivisorEvent_flag || clustered_event.HasSaturatedCluster() || clustered_event.HasSaturatedCluster(8) || clustered_event.HasLumpyCluster() || clustered_event.HasGoldenGateCluster() || clustered_event.HasBadChannelCluster())
				continue;
			if (one_and_only_one) {
				for (int det=0; det<4; det++) {
					si_avg_x += clustered_event.GetCluster(2*det,0)->Get1stMoment();
					si_avg_y += clustered_event.GetCluster(2*det+1,0)->Get1stMoment();
				}
				si_avg_x = si_avg_x/4;
				si_avg_y = si_avg_y/4;
				
				if (clustered_event.GetNClusters(8)==1)
					histo_scatter_autofidcut->Fill(si_avg_x,si_avg_y);
			}
		}
		cout << "plot production for AutoFidCut: done." << endl << endl;
		SaveHistogram(histo_scatter_autofidcut);
	}
	current_event = 0;

   //loop over events
   for(uint e=0; e<PedTree->GetEntries(); e++) {
   //for(uint e=0; e<10000; e++) {
//	   bool AlternativeClustering = false;
	   if (!AlternativeClustering) ClusterEvent();
	   else ClusterEventSeeds();
      if(e%10000==0) clustered_event.Print();
      BookHistograms();
      
      //record large clusters
      //if a track, check whether it's in the silicon fiducial region
      if(clustered_event.HasOneSiliconTrack()) {
         
         //loop over detectors to record large clusters to a txt file
         for(int det=0; det<9; det++) {
            nclusters = clustered_event.GetNClusters(det);
            //loop over all *good* clusters
            for(uint clus=0; clus<nclusters; clus++) {
               current_cluster = clustered_event.GetCluster(det,clus);
               nhits = current_cluster->GetNHits();
               if(nhits>8) {
                  large_clusters_file << event_number << "\t" << det << "\t" << clus << "/" << nclusters << "\t" << nhits << "\t" << current_cluster->GetNSeeds();
                  if(current_cluster->IsBadChannelCluster()) large_clusters_file << " (bad chan cluster)"; // skip bad clusters
                  large_clusters_file << endl;
               }
            }
         }//end large cluster record
         
      }
   }
   
   //make the hitoccupancy plots for the saturated hits
   for(int det=0; det<9; det++) {
      for(int bin=0; bin<=histo_hitocc_saturated[det][1]->GetNbinsX(); bin++) {
         histo_hitocc_saturated[det][1]->SetBinContent(bin,(float)histo_hitocc_saturated[det][1]->GetBinContent(bin)/total_events);
      }
      for(int frac=0; frac<2; frac++)
         if(det==8) histo_hitocc_saturated[det][frac]->GetXaxis()->SetRangeUser(1,63);
   }
   
   //generate the integrated etas
   float leftintegral, rightintegral;
   int nbins;
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         nbins = histo_eta[det][0][0][chip]->GetNbinsX();
         for(int bin=1; bin<=nbins; bin++) {
            leftintegral = histo_eta[det][0][0][chip]->Integral(1,bin)/histo_eta[det][0][0][chip]->Integral(1,nbins);
            rightintegral = histo_eta[det][0][0][chip]->Integral(nbins-bin+1,nbins)/histo_eta[det][0][0][chip]->Integral(1,nbins);
            histo_etaintegral[det][0][chip]->SetBinContent(bin,leftintegral);
            histo_etaintegral[det][1][chip]->SetBinContent(bin,rightintegral);
            histo_etaintegral[det][2][chip]->SetBinContent(bin,(leftintegral+rightintegral)/2.);
         }
      }
   }
   
   //eta correct hits
   double hit_position;
   for(uint t=0; t<tracks.size(); t++) {
      for(int det=0; det<8; det++) {
         hit_position = tracks[t].GetDetectorHitPosition(det);
         if(0&&t%1000==0) {
            cout<<histo_etaintegral[det][2][0]->GetTitle()<<endl;
            cout<<hit_position<<endl;
            cout<<int(hit_position)<<endl;
            cout<<int(100*(hit_position-int(hit_position)))<<endl;
            cout<<int(100*(hit_position-int(hit_position))+0.5)<<endl;
            cout<<histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
            cout<<histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
         }
         //left to right integrated eta
         if(hit_position<128) {
            nbins = histo_etaintegral[det][0][0]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            histo_hitpos[det][0]->Fill(hit_position-int(hit_position));
         }
         else {
            nbins = histo_etaintegral[det][0][1]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            histo_hitpos[det][1]->Fill(hit_position-int(hit_position));
         }
         //symmetrzed
         //if(hit_position<128) hit_position = int(hit_position) + histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //else hit_position = int(hit_position) + histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         tracks[t].SetDetectorHitPosition(det, hit_position);
      }
   }
   for(uint t=0; t<tracks_fidcut.size(); t++) {
      for(int det=0; det<9; det++) {
         hit_position = tracks_fidcut[t].GetDetectorHitPosition(det);
         if(0&&t%1000==0) {
            cout<<histo_etaintegral[det][2][0]->GetTitle()<<endl;
            cout<<hit_position<<endl;
            cout<<int(hit_position)<<endl;
            cout<<int(100*(hit_position-int(hit_position)))<<endl;
            cout<<int(100*(hit_position-int(hit_position))+0.5)<<endl;
            cout<<histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
            cout<<histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
         }
         //left to right integrated eta
         //if(hit_position<128) hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //else hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //symmetrzed
         if(hit_position<128) {
            nbins = histo_etaintegral[det][0][0]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            if(det==8) histo_hitpos[det][0]->Fill(hit_position-int(hit_position));
         }
         else {
            nbins = histo_etaintegral[det][0][1]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            if(det==8) histo_hitpos[det][1]->Fill(hit_position-int(hit_position));
         }
         tracks_fidcut[t].SetDetectorHitPosition(det, hit_position);
      }
   }
   
	AutoFidCut(); //added 2010-12-01 (max)
	
   //Draw plots and generate HTML
   if(plots) {
      DrawHistograms();
      
      
      GenerateHTML();
      ////DeleteHistograms(); //this causes segfaults for some unknown reason
   }
   TCanvas tempcanv("tempcanv");
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         ostringstream tempsstream;
         tempsstream << plots_path << histo_hitpos[det][chip]->GetTitle() << ".png";
         histo_hitpos[det][chip]->Draw();
         tempcanv.Print(tempsstream.str().c_str());
      }
   }
   
   large_clusters_file.close();
   PedFile->Close(); // somehow histo_etaintegral[0][2][0] is deleted by this line?!
   
   
   cout<<"Cut flow is as follows:"<<endl;
   cout<<"total_events = "<<total_events<<endl;
   cout<<"[ ]CMNEvents = "<<CMNEvents<<endl;
   cout<<"[ ]ZeroDivisorEvents = "<<ZeroDivisorEvents<<endl;
   cout<<"[a]totalsurviving_events = "<<totalsurviving_events<<endl;
   cout<<"   [ ]goldengatecluster_events[9] = "<<goldengatecluster_events[9]<<endl;
   cout<<"   [ ]badchannelcluster_events[9] = "<<badchannelcluster_events[9]<<endl;
   cout<<"   [ ]lumpycluster_events[9] = "<<lumpycluster_events[9]<<endl;
   cout<<"   [ ]saturatedcluster_events[9] = "<<saturatedcluster_events[9]<<endl;
   for(int d=0; d<4; d++)
      cout<<"   [ ]detectorxycluster_events["<<d<<"] = "<<detectorxycluster_events[d]<<endl;
   cout<<"   [b]singlesitrack_events = "<<singlesitrack_events<<endl;
   cout<<"      [ ]singlesitrack_1diamondclus_events = "<<singlesitrack_1diamondclus_events<<endl;
   cout<<"      [c]singlesitrack_fidcut_events = "<<singlesitrack_fidcut_events<<endl;
   cout<<"         [ ]singlesitrack_fidcut_1diamondclus_events = "<<singlesitrack_fidcut_1diamondclus_events<<endl;
         
   Watch.Stop();
   Watch.Print("u");
   cout<<"tracks.size() = "<<tracks.size()<<endl;
   cout<<"tracks_fidcut.size() = "<<tracks_fidcut.size()<<endl;
   cout<<"counter_alignment_tracks_zero_suppressed = "<<counter_alignment_tracks_zero_suppressed<<endl;
   cout<<"counter_alignment_tracks = "<<counter_alignment_tracks<<endl;
   cout<<"counter_alignment_only_tracks = "<<counter_alignment_only_tracks<<endl;
   cout<<"counter_alignment_only_tracks/float(counter_alignment_tracks) = "<<counter_alignment_only_tracks/float(counter_alignment_tracks)<<endl;
   cout<<"counter_alignment_fidcut_tracks = "<<counter_alignment_fidcut_tracks<<endl;
   cout<<"counter_alignment_only_fidcut_tracks = "<<counter_alignment_only_fidcut_tracks<<endl;
   cout<<"counter_alignment_only_fidcut_tracks/float(counter_alignment_fidcut_tracks) = "<<counter_alignment_only_fidcut_tracks/float(counter_alignment_fidcut_tracks)<<endl;
   cout<<"singlesitrack_fidcut_events = "<<singlesitrack_fidcut_events<<endl;
   
}

void Clustering::HistCleaner(int regid, TH2F* histo) {
    
    
	cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
	cout << endl << "HistCleaner called on histogram " << histo->GetName() << "." << endl;
	
	int afc_current_histo_x_colsums[256];
	int afc_current_histo_y_colsums[256];
	
	TH1F* afc_current_histo_x = new TH1F("afc_current_histo_x","afc_current_histo_x",256,0,255);
	TH1F* afc_current_histo_y = new TH1F("afc_current_histo_y","afc_current_histo_y",256,0,255);
	
	for(int i=0;i<256;i++) {
		afc_current_histo_x_colsums[i] = 0;
		for(int j=0;j<256;j++) {
			afc_current_histo_x_colsums[i] = afc_current_histo_x_colsums[i] + histo->GetBinContent(i,j);
			afc_current_histo_x->SetBinContent(i,afc_current_histo_x_colsums[i]);
		}
		
	}
	
	//calculate x bin values for each y row and count # of non-zero bins and set the corresponding histogram bins
	for(int j=0;j<256;j++) {
		afc_current_histo_y_colsums[j] = 0;
		for(int i=0;i<256;i++) {
			afc_current_histo_y_colsums[j] = afc_current_histo_y_colsums[j] + histo->GetBinContent(i,j);
			afc_current_histo_y->SetBinContent(j,afc_current_histo_y_colsums[j]);
		}
		
	}
	
	int xpos = 2;
	int ypos = 2;
    
	cout << endl;
	cout << "Results for histo " << histo->GetName() << " with region id " << regid << " are:" << endl << endl;
    
	
	for(int i=xpos;i<256;i++) {
		if( (afc_current_histo_x->GetBinContent(i-1) == 0) && (afc_current_histo_x->GetBinContent(i) > 0) ) {
			//FCR[z]->SetValueXLow(i+5);
			cout << "Setting " << histo->GetName() << ".ValueXLow to: " << i+5 << endl;
			for(int j=i+5;j<256;j++) {
				if ( (afc_current_histo_x->GetBinContent(j-1) > 0) && (afc_current_histo_x->GetBinContent(j) == 0) ) {
                    //					FCR[z]->SetValueXHigh(j-5);
					cout << "Setting " << histo->GetName() << ".ValueXHigh to: " << j-5 << endl;					
					xpos=j+5;
					break;
				}} break;
			//break;
		}}
    
	for(int j=ypos;j<256;j++) {
		if( (afc_current_histo_y->GetBinContent(j-1) == 0) && (afc_current_histo_y->GetBinContent(j) > 0) ) {
			//FCR[z]->SetValueXLow(i+5);
			cout << "Setting " << histo->GetName() << ".ValueYLow to: " << j+5 << endl;
			for(int k=j+5;k<256;k++) {
				if ( (afc_current_histo_y->GetBinContent(k-1) > 0) && (afc_current_histo_y->GetBinContent(k) == 0) ) {
					//					FCR[z]->SetValueXHigh(j-5);
					cout << "Setting " << histo->GetName() << ".ValueYHigh to: " << k-5 << endl;					
					ypos=k+5;
					break;
				}} break;
			//break;
		}}
	
	
	afc_current_histo_x->Delete();
	afc_current_histo_y->Delete();
    
	cout << endl << "Finished HistCleaner on " << histo->GetName() << "." << endl;  
	cout << "HistCleaner exiting..." << endl << endl;
}

void Clustering::AutoFidCut() {
		
	histo_afc_unit_histo_1f->Reset(); // reset 1f unit histogram (all bins zero)
	histo_afc_unit_histo_2f->Reset(); // reset 2f unit histogram (all bins zero)
    
	//SaveHistogramPDF(histo_afc_unit_histo_1f);
	//SaveHistogramPDF(histo_afc_unit_histo_2f);
	
	//define sum and counter vectors for x columns and y rows
	int histo_afc_col_sums_x[256];
	int histo_afc_col_sums_y[256];
	
	
    if(Verbosity>0) {
    cout << "\n";
	cout << "AutoFidCut: I'm the AutoFidCut function.\n\n";
	cout << " START AutoFidCut \n\n";	
    }    

    if(Verbosity>0) {    
	// produce necessary plots to detect fidcut region
	cout << endl << endl << "-- produce plot for AutoFidCut().." << endl;
	cout << "running over " << PedTree->GetEntries() << " events.." << endl;
    }
    
	//record large clusters
	string large_clusters_filename = plots_path + "large_clusters.txt";
	ofstream large_clusters_file(large_clusters_filename.c_str());
	large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
	uint nclusters, nhits;
	Cluster* current_cluster;
	
	current_event = 0;
	for (uint e=0; e<PedTree->GetEntries(); e++) { // maybe it's not necessary to run over all events for the AutoFidCut?!
		if (!AlternativeClustering) ClusterEvent();
		else ClusterEventSeeds();
		if (e%10000==0) clustered_event.Print();
		
		// current_cluster = 0;
		
		// -- produce scatter plot for AutoFidCut
		bool one_and_only_one = clustered_event.HasOneSiliconTrack();
		Float_t si_avg_x=0, si_avg_y=0;
		if (CMNEvent_flag || ZeroDivisorEvent_flag || clustered_event.HasSaturatedCluster() || clustered_event.HasSaturatedCluster(8) || clustered_event.HasLumpyCluster() || clustered_event.HasGoldenGateCluster() || clustered_event.HasBadChannelCluster())
			continue;
		if (one_and_only_one) {
			for (int det=0; det<4; det++) {
				si_avg_x += clustered_event.GetCluster(2*det,0)->Get1stMoment();
				si_avg_y += clustered_event.GetCluster(2*det+1,0)->Get1stMoment();
			}
			si_avg_x = si_avg_x/4;
			si_avg_y = si_avg_y/4;
			
			if (clustered_event.GetNClusters(8)==1)
				histo_scatter_autofidcut->Fill(si_avg_x,si_avg_y);
		}
	}
    
    if(Verbosity>0) {
	cout << "plot production for AutoFidCut: done." << endl << endl;
    }
	
    SaveHistogram(histo_scatter_autofidcut);
	
	// -- end of scatter plot production
    
    //calculate y bin values for each x column and count # of non-zero bins and set the corresponding histogram bins
    for(int i=0;i<256;i++) {
        histo_afc_col_sums_x[i] = 0;
        for(int j=0;j<256;j++) {
            histo_afc_col_sums_x[i] = histo_afc_col_sums_x[i] + histo_scatter_autofidcut->GetBinContent(i,j);
            histo_afc_x->SetBinContent(i,histo_afc_col_sums_x[i]);
        }
        
    }
    
    //calculate x bin values for each y row and count # of non-zero bins and set the corresponding histogram bins
    for(int j=0;j<256;j++) {
        histo_afc_col_sums_y[j] = 0;
        for(int i=0;i<256;i++) {
            histo_afc_col_sums_y[j] = histo_afc_col_sums_y[j] + histo_scatter_autofidcut->GetBinContent(i,j);
            histo_afc_y->SetBinContent(j,histo_afc_col_sums_y[j]);
        }
        
    }
    
    //algorithm for finding numeric differences between maxima of the x and y and a certain cut factor (afc_cut_factor)
    
    //get maximums for x and y histogram
    int histo_afc_x_max = (int)histo_afc_x->GetMaximum();
    int histo_afc_y_max = (int)histo_afc_y->GetMaximum();
    
    //set cut factor for pre-fidcut
    float afc_cut_factor = 0.4;
    
    //write maxima of x and y histograms to stdout
    if(Verbosity>0) {
    cout << "histo_afc_x->GetMaximum() is:\t" << (int)histo_afc_x->GetMaximum() << "\n";
    cout << "histo_afc_y->GetMaximum() is:\t" << (int)histo_afc_y->GetMaximum() << "\n";
    }
        
    //set x and y fidcut histogram bin contents according to afc_cut_factor and maxima (set to zero if values are less than afc_cut_factor*maximum, keep if values are greater)
    for(int i=0;i<256;i++) {
        if( (histo_afc_x->GetBinContent(i)) < (afc_cut_factor*histo_afc_x_max)) {
            histo_afc_x_cut->SetBinContent(i,0);
        } else {
            histo_afc_x_cut->SetBinContent(i,histo_afc_x->GetBinContent(i));
        }
    }
    
    for(int j=0;j<256;j++) {
        if( (histo_afc_y->GetBinContent(j)) < (afc_cut_factor*histo_afc_y_max)) {
            histo_afc_y_cut->SetBinContent(j,0);
        } else {
            histo_afc_y_cut->SetBinContent(j,histo_afc_y->GetBinContent(j));
        }
        
    }	
    
    
    //TH2F histo_afc_scatter_firstfidcut = new TH2F("histo_afc_scatter_firstfidcut",256,0,255,256,0,255);
    TH2F *histo_afc_scatter_firstfidcut = (TH2F*)histo_scatter_autofidcut->Clone();
    histo_afc_scatter_firstfidcut->SetName("histo_afc_scatter_firstfidcut");
    
    int afc_width_min = 20; //variable for defining minimum block width to be recognized as partial fid cut region
    
    int afc_width_counter = 0; //temporary counter for the width of the current block
    int afc_last_start = 0; //temporary variable to store last start position in block finder algorithm
    int afc_last_end = 0; //temporary variable to store last end position in block finder algorithn
    
    int afc_block_count_x = 0;
    int afc_block_count_y = 0;
    
    //loop for y direction to count projected consecutive blocks with a block width over afc_width_min
    for(int j=1;j<256;j++) {
        if ( (histo_afc_x_cut->GetBinContent(j) > 0) && (histo_afc_x_cut->GetBinContent(j-1) == 0) ) {
            afc_last_start = j;
            afc_width_counter++;
            afc_block_count_x++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 1 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        } else if ( (histo_afc_x_cut->GetBinContent(j) > 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) ) {
            afc_width_counter++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 2 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) && (afc_width_counter < afc_width_min) ) {
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 3 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            for(int k=afc_last_start;k<=j;k++) {
                histo_afc_x_cut->SetBinContent(k,0);
            }
            afc_block_count_x--;
            afc_width_counter = 0;
            continue;
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) && (afc_width_counter > afc_width_min) ) {
            afc_width_counter = 0;
            afc_block_count_x++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 5 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) ) {
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 0 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        }
        else {
            afc_width_counter++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 4 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        }
    }
    
    //loop for y direction to count projected consecutive blocks with a block width over afc_width_min    
    for(int j=1;j<256;j++) {
        if ( (histo_afc_y_cut->GetBinContent(j) > 0) && (histo_afc_y_cut->GetBinContent(j-1) == 0) ) {
            afc_last_start = j;
            afc_width_counter++;
            afc_block_count_y++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 1 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        } else if ( (histo_afc_y_cut->GetBinContent(j) > 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) ) {
            afc_width_counter++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 2 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) && (afc_width_counter < afc_width_min) ) {
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 3 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            for(int k=afc_last_start;k<=j;k++) {
                histo_afc_y_cut->SetBinContent(k,0);
            }
            afc_block_count_y--;
            afc_width_counter = 0;
            continue;
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) && (afc_width_counter > afc_width_min) ) {
            afc_width_counter = 0;
            afc_block_count_y++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 5 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) ) {
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 0 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        }
        else {
            afc_width_counter++;
            if(Verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 4 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        }
    }
    
	
	
    //divide the detector plane in a maximum of four quadrants
    
    if(Verbosity>=1) {
    cout << "#-#-#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#" << endl;
    cout << "#" << endl;
    cout << "# afc region recognition started..." << endl;
    cout << "#" << endl;
	}
	
    int afc_x_divider = 0;
    int afc_y_divider = 0;
    bool one_and_only_afc_region_x = 1;
    bool one_and_only_afc_region_y = 1;
    bool one_and_only_afc_region = 0;
	
    if(Verbosity>=2) cout << "# x divider function called." << endl;
    for(int i=2;i<256;i++) {			
        if( (histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0) ) {
            for(int j=i;j<256;j++) {
                if ( (histo_afc_x_cut->GetBinContent(j-1) > 0) && (histo_afc_x_cut->GetBinContent(j) == 0) ) {
                    afc_x_divider = j+5;
                    break;
                }
            }
            break;
        }
    }
    
    for(int k=afc_x_divider-4;k<256;k++) {
        if(histo_afc_x_cut->GetBinContent(k)!=0) {
            one_and_only_afc_region_x = 0;
        }
    }
	
    if (one_and_only_afc_region_x) {
        afc_x_divider = 0;
        if(Verbosity>=1) cout << "# this run has one and only one afc region in x direction. afc_x_divider set to 0." << endl;
    } else {
        if(Verbosity>=1) cout << "# this run has more than one afc region in x direction." << endl;
    }
	
    
    if(Verbosity>=1) cout << "#" << endl << "# y divider function called." << endl;
    for(int i=2;i<256;i++) {
        if( (histo_afc_y_cut->GetBinContent(i-1) == 0) && (histo_afc_y_cut->GetBinContent(i) > 0) ) {
            for(int j=i;j<256;j++) {
                if ( (histo_afc_y_cut->GetBinContent(j-1) > 0) && (histo_afc_y_cut->GetBinContent(j) == 0) ) {
                    afc_y_divider = j+5;
                    break;
                }
            }
            break;
        }
    }
    
    for(int k=afc_y_divider-4;k<256;k++) {
        if(histo_afc_y_cut->GetBinContent(k)!=0) {
            one_and_only_afc_region_y = 0;
        }
    }
	
    if (one_and_only_afc_region_y) {
        afc_y_divider = 0;
        if(Verbosity>=1) cout << "# this run has one and only one afc region in y direction. afc_y_divider set to 0." << endl;
    } else {
        if(Verbosity>=1) cout << "# this run has more than one afc region in y direction." << endl;
    }
	
    if( one_and_only_afc_region_x && one_and_only_afc_region_y ) {
        one_and_only_afc_region = 1;
        if(Verbosity>=1) cout << "#" << endl << "#" << endl << "# this run has one and only one afc region." << endl;
    } else {
        if(Verbosity>=1) cout << "#" << endl << "#" << endl << "# this run has more than one afc region." << endl;
    }
	
    if(Verbosity>=1) cout << "#" << endl << "# afc region recognition ended..." << endl << "#" << endl << "#-#-#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#Ð#" << endl;
	
    //end divider code
	
	if(Verbosity>=1) { 
    cout << endl;
    cout << "# # # # # # # # # # # # # # # # # # # # # # #" << endl;
    cout << "#" << endl;
    cout << "# Famous afc_x_divider has value: " << afc_x_divider << endl;
    cout << "#" << endl;
    cout << "# Famous afc_y_divider has value: " << afc_y_divider << endl;
    cout << "#" << endl;
    cout << "# # # # # # # # # # # # # # # # # # # # # # #" << endl;	
    cout << endl;
    }
	
    for(int z=1;z<2;z++) {
        for(int u=0;u<afc_x_divider;u++) {
            for(int l=1;l<afc_y_divider;l++) {
                histo_afc_unit_histo_2f->SetBinContent(u,l,1);
                
            }
        }
    }
    
    //		for(int z=1;z<2;z++) {
    //			for(int u=afc_x_divider;u<256;u++) {
    //				for(int l=1;l<256;l++) {
    //					histo_afc_unit_histo_2f->SetBinContent(u,l,0);
    //				
    //				}
    //			}
    //		}
	
	
	
    histo_afc_unit_histo_2f->SetName("histo_afc_unit_histo_2f_r1");
    
    SaveHistogramPDF(histo_afc_unit_histo_2f);
	
    (*histo_afc_region_1_mask) = (*histo_afc_unit_histo_2f)*(*histo_afc_scatter_firstfidcut);
	
    histo_afc_region_1_mask->SetName("histo_afc_region_1_mask");
	
    SaveHistogramPDF(histo_afc_region_1_mask);
    if(Verbosity>=2) cout << endl << endl << "afc_region_1 saved" << endl << endl;
	
    // run HistCleaner on 
    if(Verbosity>=2) cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
    HistCleaner(9999, histo_afc_region_1_mask);
    if(Verbosity>=2) cout << endl << "Finished HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
	
    //		SaveHistogramPDF(histo_afc_region_1_mask);
    //		cout << endl << endl << "afc_region_1 saved" << endl << endl;
	
    //end divider stuff
	
    // set histogram bin value for histo_afc_scatter_firstfidcut
    for(int i=0;i<256;i++) {
        if(histo_afc_x_cut->GetBinContent(i) == 0) {
            for(int k=0;k<256;k++) {
                histo_afc_scatter_firstfidcut->SetBinContent(i,k,0);
            }
        }
    }
    
    for(int j=0;j<256;j++) {
        if(histo_afc_y_cut->GetBinContent(j) == 0) {
            for(int k=0;k<256;k++) {
                histo_afc_scatter_firstfidcut->SetBinContent(k,j,0);
            }
        }
    }
    
	//define counter variables for region counters
    int afc_region_counter_x = 0;
    int afc_region_counter_y = 0;
	
    //count number of consecutive regions in x direction
    for(int i=2;i<256;i++) {
        if( ((histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0)) || ((histo_afc_x_cut->GetBinContent(i-1) > 0) && (histo_afc_x_cut->GetBinContent(i) == 0))) {	
            afc_region_counter_x++;
        }
    }
	
    //count number of consecutive regions in y direction
    for(int i=2;i<256;i++) {
        if( ((histo_afc_y_cut->GetBinContent(i-1) == 0) && (histo_afc_y_cut->GetBinContent(i) > 0)) || ((histo_afc_y_cut->GetBinContent(i-1) > 0) && (histo_afc_y_cut->GetBinContent(i) == 0))) {	
			afc_region_counter_y++;
        }
    }
	
    if(Verbosity>=1) {
    cout << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    cout << "It is most likely that this detector setup has " << (afc_region_counter_x/2)*(afc_region_counter_y/2) << " active diamond regions, i guess." << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    }
	
    int xpos = 2;
    int ypos = 2;
	
    for(int z=1;z<((afc_region_counter_x/2)*(afc_region_counter_y/2))+1;z++) {
        
		for(int i=xpos;i<256;i++) {
			if( (histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0) ) {
				FCR[z]->SetValueXLow(i+5);
				if(Verbosity>=2) cout << "Setting FCR[" << z << "].ValueXLow to: " << i+5 << endl;
				for(int j=i+5;j<256;j++) {
					if ( (histo_afc_x_cut->GetBinContent(j-1) > 0) && (histo_afc_x_cut->GetBinContent(j) == 0) ) {
                        FCR[z]->SetValueXHigh(j-5);
                        if(Verbosity>=2) cout << "Setting FCR[" << z << "].ValueXHigh to: " << j-5 << endl;					
                        xpos=j+5;
                        break;
					}} break;
                //break;
			}}
        
		
    }
	
    if(Verbosity>=2) {
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    cout << "FCR[1]->ValueXLow is: " << FCR[1]->GetValueXLow() << endl;
    cout << "FCR[1]->ValueXHigh is: " << FCR[1]->GetValueXHigh() << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    }
    
    int afc_regions_count = (afc_block_count_x/2)*(afc_block_count_y/2);
    int afc_regions_safety_delta = 5;
    
    if(Verbosity>=2) {
    cout << "\nafc_block_count_x is:\t" << afc_block_count_x << "\n";
    cout << "afc_block_count_y is:\t" << afc_block_count_y << "\n";
    cout << "That means, we have " << afc_block_count_x/2 << " block(s) in x- and " << afc_block_count_y/2 << " block(s) in y-direction. so: " << afc_regions_count << " in total.\n\n";
    }
    
    
    //take afc_block_counts_x and afc_block_counts_y and populate FCRs
    
    
    for(int i=0;i<afc_regions_count;i++) {
        FCR[i]->SetAllValuesZero();
        
        
        //		FCR[i]->SetValueXLow(some_var+afc_regions_safety_delta);
        //		FCR[i]->SetValueXHigh(some_var-afc_regions_safety_delta);
        //		FCR[i]->SetValueYLow(some_var+afc_regions_safety_delta);
        //		FCR[i]->SetValueYHigh(some_var-afc_regions_safety_delta);
    }
    
    histo_afc_scatter_firstfidcut->SetTitle("First Fidcut on Scatter Plot");
    histo_afc_scatter_firstfidcut->SetTitleFont(42);
    
    histo_afc_region_1_mask->Reset();
	
    (*histo_afc_region_1_mask) = (*histo_afc_unit_histo_2f)*(*histo_afc_scatter_firstfidcut);
	
    histo_afc_region_1_mask->SetName("histo_afc_region_1_mask_LATEST");
    
    int afc_region_start = 0;
	
    if(afc_y_divider == 0) {
        afc_region_start = 1;
    }
	
	for(int i=1;i<5;i++) {
		
		ostringstream phname;
		phname << "histo_current_2f_" << i;
		
		TH2F* afc_current_histo_2f = new TH2F(phname.str().c_str(),phname.str().c_str(),256,0,255,256,0,255);
		
		histo_afc_unit_histo_2f->Reset();
		if(i==1) {
			for(int a=2;a<afc_x_divider;a++) {
				for(int b=1;b<afc_y_divider;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
        
		if(i==2) {
			for(int a=afc_x_divider;a<256;a++) {
				for(int b=1;b<afc_y_divider;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		if(i==3) {
			for(int a=2;a<afc_x_divider;a++) {
				for(int b=afc_y_divider;b<256;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		if(i==4) {
			for(int a=afc_x_divider;a<256;a++) {
				for(int b=afc_y_divider;b<256;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		
		(*afc_current_histo_2f) = (*histo_afc_scatter_firstfidcut)*(*histo_afc_unit_histo_2f); 
		
		afc_current_histo_2f->SetName(phname.str().c_str());
		
        if(Verbosity>=1) {
		cout << "TENGTARATENG" << endl;
		cout << "FOR LOOP FOR FIDCUT REGIONS. NAME OF CURRENT HISTO IS: " << afc_current_histo_2f->GetName() << "." << endl;
        }
		HistCleaner(i, afc_current_histo_2f);
		
		
		SaveHistogramPDF(afc_current_histo_2f);
		
		afc_current_histo_2f->Delete();
		
	}
	
	
    
    //		cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
    HistCleaner(1, histo_afc_region_1_mask);
    //		 cout << endl << "Finished HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;  
	
	
    SaveHistogramPDF(histo_afc_region_1_mask);
    
    HistCleaner(1, histo_afc_scatter_firstfidcut);
	
    FCR[0]->GetAllValues();
    
    FCR[3]->SetAllValuesZero();
    FCR[3]->SetValueXLow(1);
    FCR[3]->SetValueYHigh(2);
    FCR[3]->GetAllValues();
    
    SaveHistogramPDF(histo_afc_x);
    SaveHistogramPDF(histo_afc_y);
    SaveHistogramPDF(histo_afc_x_cut);
    SaveHistogramPDF(histo_afc_y_cut);
    SaveHistogramPDF(histo_afc_scatter_firstfidcut);
    
    
    SaveHistogramPNG(histo_afc_x);
    SaveHistogramPNG(histo_afc_y);
    SaveHistogramPNG(histo_afc_x_cut);
    SaveHistogramPNG(histo_afc_y_cut);
    SaveHistogramPNG(histo_afc_scatter_firstfidcut);
    
    if(Verbosity>=0) {
    cout << "\n";
    cout << "FINISHED AutoFidCut \n\n\n";
    }
}	

void Clustering::Align(bool plots, bool CutFakeTracksOn) {
	
	if(tracks.size()==0) {
		cout<<"Clustering::Align: No tracks found; calling Clustering::ClusterRun first..."<<endl;
		ClusterRun(plots); // doesn't use alternative clustering
	}
	
	if (tracks.size() == 0) {
		cout << "No tracks available. Alignment not possible. (tracks.size() = " << tracks.size() << ")" << endl;
		return;
	}
	
/*	int count_true = 0;
	int count_false = 0;
	for (int i = 0; i < tracks_mask.size(); i++) {
		cout << "track " << i << " has mask " << tracks_mask[i] << endl;
		if (tracks_mask[i]) count_true++;
		if (!tracks_mask[i]) count_false++;
	}
	cout << count_true << " tracks are masked as true and " << count_false << " tracks as false." << endl;
	return;*/
	
	string plots_path_save = plots_path;
	vector<TDiamondTrack> alignment_tracks = tracks;
	vector<bool> alignment_tracks_mask = tracks_mask;
	vector<TDiamondTrack> alignment_tracks_fidcut = tracks_fidcut;
	vector<bool> alignment_tracks_fidcut_mask = tracks_fidcut_mask;
	ostringstream plots_path_alignment, plots_path_alignment_CutFakeTracks;
	plots_path_alignment << plots_path << "/alignment/";
	plots_path_alignment_CutFakeTracks << plots_path << "/alignment_CutFakeTracks/";
	plots_path = plots_path_alignment.str();
	
	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
	for (int alignStep = 0; alignStep < 2; alignStep++) {
		sys->mkdir(plots_path.c_str());
//		TDetectorAlignment* align = new TDetectorAlignment(plots_path, tracks, tracks_mask);
		TDetectorAlignment* align = new TDetectorAlignment(plots_path, alignment_tracks, alignment_tracks_mask);
		
		Int_t nPasses = 10;
		Double_t plot_width_factor = 3; // scales the widths of the plots; range is a 3*width of distribution centered on mean
		
		align->PlotAngularDistribution(); //look at angular distribution of tracks
		align->PlotCartesianDistribution(); //look at slope distribution of tracks
		
		string prename = "alignment_PrealignmentResiduals";
		string postname = "alignment_PostalignmentResiduals";
		
		// generate residuals before alignment
		align->CheckDetectorAlignmentXYPlots(0, 1, 3, prename);
		align->CheckDetectorAlignmentXYPlots(1, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(2, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(3, 0, 2, prename);
		
		// itterative alignment loop
		for(int i=0; i<nPasses; i++) {
			cout << "\n\nPass " << i+1 << endl<< endl;
			//xy alignment
			align->CheckDetectorAlignmentXY(0, 1, 3);
			align->AlignDetectorXY(1, 0, 3);
			align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);
			//align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);
			//align->AlignDetectorXY(2, 0, 3);
			align->AlignDetectorXY(3, 0, 2);
			//align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);
			
			//phi alignment: this chunk of code causes seg faulting in code at top of loop!
			//align->AlignDetectorPhi(1, 0, 3, false, false);
			//align->AlignDetectorPhi(2, 0, 3, false, false);
			//align->AlignDetectorPhi(3, 0, 2, false, false);
			
			//phi alignment: this chunk of code causes seg faulting in code at top of loop!
			//align->AlignDetectorZ(1, 0, 3, false, false);
			//align->AlignDetectorZ(2, 0, 3, false, false);
			//align->AlignDetectorZ(3, 0, 2, false, false);
		}
		
		cout<<endl;
		cout<<endl;
		cout<<"Checking final residuals"<<endl;
		cout<<endl;
		
		// generate residuals after alignment
		align->CheckDetectorAlignmentXYPlots(0, 1, 3, postname);
		align->CheckDetectorAlignmentXYPlots(1, 0, 3, postname);
		align->CheckDetectorAlignmentXYPlots(2, 0, 3, postname);
		align->CheckDetectorAlignmentXYPlots(3, 0, 2, postname);
		
		cout<<endl;
		
		
		//Now align the diamond
		
		//load fidcut tracks w/ 1 diamond cluster
//		align->LoadTracks(tracks_fidcut, tracks_fidcut_mask);
		align->LoadTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask);
		
		//check that the silicon is still aligned for these tracks_fidcut
		cout<<"Check that the telescope alignment still holds for fidcut tracks w/ single diamond cluster"<<endl;
		align->CheckDetectorAlignmentXY(0, 1, 3);
		align->CheckDetectorAlignmentXY(1, 0, 3);
		align->CheckDetectorAlignmentXY(2, 0, 3);
		align->CheckDetectorAlignmentXY(3, 0, 2);
		
		// generate residuals before alignment
		align->CheckDetectorAlignmentXYPlots(4, 1, 2, prename);
		
		// itterative alignment loop
		for(int i=0; i<5; i++) {
			cout << "\n\nPass " << i+1 << endl<< endl;
			//xy alignment
			align->AlignDetectorXY(4, 1, 2);
			//align->AlignDetectorXY(4, 0, 3);
			//align->AlignDetectorXY(4, 1, 3);
			//align->AlignDetectorXY(4, 0, 2);
			//align->AlignDetectorXY(4, 1, 2);
		}
		
		cout<<endl;
		cout<<endl;
		cout<<"Checking final diamond residuals"<<endl;
		cout<<endl;
		
		// generate residuals after alignment
		align->CheckDetectorAlignmentXYPlots(4, 1, 2, postname);
		
		
		//report results in a file
		
		ostringstream alignment_summary_path;
		alignment_summary_path << plots_path << "Alignment_Summary.txt";
		ofstream alignment_summary(alignment_summary_path.str().c_str());
		
		alignment_summary << "Alignment summary " << endl;
		alignment_summary << "----------------- " << endl << endl;
		
		alignment_summary << "Offsets (multiples of 50um):" << endl << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_x_offset[" << det << "] = " << align->det_x_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_y_offset[" << det << "] = " << align->det_y_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_z_offset[" << det << "] = " << align->det_z_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_phix_offset[" << det << "] = " << align->det_phix_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_phiy_offset[" << det << "] = " << align->det_phiy_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		
		alignment_summary << "Resolutions (multiples of 50um):" << endl << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_x_resolution[" << det << "] = " << align->det_x_resolution[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << "det_y_resolution[" << det << "] = " << align->det_y_resolution[det] <<endl; 
		}
		
		
		alignment_summary << endl << endl;
		alignment_summary << "Alignment summary (for pasting into a spread sheet) " << endl;
		alignment_summary << "--------------------------------------------------- " << endl << endl;
		
		alignment_summary << "Offsets (multiples of 50um):" << endl << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_x_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_y_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_z_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_phix_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_phiy_offset[det] <<endl; 
		}
		alignment_summary << endl; 
		
		alignment_summary << "Resolutions (multiples of 50um):" << endl << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_x_resolution[det] <<endl; 
		}
		alignment_summary << endl; 
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_y_resolution[det] <<endl; 
		}
		alignment_summary.close();
		
		cout << "Intrinsic silicon resolution " << align->GetSiResolution() << " strips or " << align->GetSiResolution() * 50 << "um" << endl;
		
		if (CutFakeTracksOn && alignStep < 1) {
			align->CutFakeTracks(alignment_tracks, alignment_tracks_mask, alignment_chi2, CutFakeTracksOn, true);
			align->CutFakeTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, alignment_chi2, CutFakeTracksOn, true);
			plots_path = plots_path_alignment_CutFakeTracks.str();
		}
		if (!CutFakeTracksOn || alignStep == 1) {
			dia_offset.clear();
//			dia_offset.push_back(align->GetXOffset(5));
//			dia_offset.push_back(align->GetYOffset(5));
//			dia_offset.push_back(align->GetZOffset(5));
//			cout << "align->GetXOffset(4) = " << align->GetXOffset(4) << endl;
//			cout << "align->GetXOffset(5) = " << align->GetXOffset(5) << endl;
//			TransparentClustering(alignment_tracks, alignment_tracks_mask, align);
			TransparentClustering(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, align);
			break;
		}
	} // end alignment loop
   /*
   //Plot out the offsets
   for(Int_t plane=1; plane<4; plane++) {
      align->PlotXOffsetHistory(plane);
      align->PlotYOffsetHistory(plane);
      align->PlotZOffsetHistory(plane);
      align->PlotPhiOffsetHistory(plane);
   }
   
   //Print out the offsets
   for(Int_t plane=1; plane<4; plane++) {
      cout<<"Detector D"<<plane<<":\t";
      cout<<align->GetXOffset(plane)<<"\t";
      cout<<align->GetYOffset(plane)<<"\t";
      cout<<align->GetZOffset(plane)<<"\t";
      cout<<align->GetPhiOffset(plane)<<endl;
   }
   */
	plots_path = plots_path_save;
}

// -- take track after alignment, point into diamond and produce a cluster with the surrounding channels
void Clustering::TransparentClustering(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose) {
	cout << "Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
	int event;
	vector<Float_t> x_positions, y_positions, z_positions, par;
	Float_t diamond_hit_position = 0;
	int diamond_hit_channel = 0;
	int diamond_secondhit_channel = 0;
	int current_channel, current_sign;
	Float_t cluster_adc = 0;
	Float_t transp_eta = 0;
	Float_t firstchannel_adc, secondchannel_adc;
	vector<int> event_numbers;
	Double_t diamond_x_offset = align->GetXOffset(4);
	Double_t diamond_z_position = 10.2; // TODO: is the z position of the diamond always 10.2??
	Float_t eff_diamond_hit_position = 0;
	int eff_diamond_hit_channel = 0;
	Float_t hit_diff = 0;
	
	cout << "diamond x offset: " << diamond_x_offset << endl;
	cout << "diamond hit factor: " << Di_Cluster_Hit_Factor << endl;
	
	cout << "init histograms for transparent clustering.." << endl;
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_Dia_" << (i+1) << "HitTransparClusters";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),200,0.,2000.);
//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	}
	ostringstream histoname_eta;
	histoname_eta << "Eta_Dia_2CentroidHits_TransparClusters";
	histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
	histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_Dia_TransparClusters","DiffEstEffHit_Dia_TransparClusters", 100, 0.,5.);
	cout << " done." << endl;
	
	verbose = true;
	
	PedFile = new TFile(pedfile_path.c_str());
	PedTree = (TTree*)PedFile->Get("PedTree");
	if (!PedTree)
	{
		cerr << "PedTree not found!" << endl;
	}
	
	cout << PedTree->GetEntries() << " events in PedTree " << endl;
	
	//Event Header Branches
	PedTree->SetBranchAddress("RunNumber",&run_number);
	PedTree->SetBranchAddress("EventNumber",&event_number);
	PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
	PedTree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
	PedTree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);
	
	//Telescope Data Branches
	PedTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
	PedTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
	PedTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
	PedTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
	PedTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
	PedTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
	PedTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
	PedTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
	PedTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
	PedTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
	PedTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
	PedTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
	PedTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
	PedTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
	PedTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
	PedTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
	PedTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
	PedTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
	PedTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
	PedTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
	PedTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
	PedTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
	PedTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
	PedTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
	PedTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
	PedTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
	PedTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
	PedTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
	PedTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
	PedTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
	PedTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
	PedTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
	PedTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
	PedTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
	PedTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
	PedTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
	PedTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
	PedTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
	PedTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
	PedTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
	PedTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
	PedTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
	PedTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
	PedTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
	PedTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
	
	PedTree->GetEvent(0);
	cout<< "Loaded first event in PedTree: "<<event_number<<endl;
	cout<< "RunNumber is: "<<run_number<<endl;
	cout<< "StoreThreshold is: "<<store_threshold<<endl;
	
	event_numbers.clear();
	for (int j = 0; j < PedTree->GetEntries(); j++) {
		PedTree->GetEvent(j);
//		cout << endl << endl << endl;
//		cout << "event " << j << " in PedTree has event_number: " << event_number << endl;
		event_numbers.push_back(event_number);
//		for (int blablabla = 0; blablabla < 254; blablabla++) {
//			cout << "Dia_ADC = " << Dia_ADC[blablabla] << ",\tDet_PedMean = " << Det_PedMean[8][blablabla] << endl;
//			cout << "Collected charge in channel " << blablabla << " of diamond: " << Dia_ADC[blablabla]-Det_PedMean[8][blablabla] << endl;
//		}
	}
//	return;
	
	// loop over tracks
	for (int i = 0; i < tracks.size(); i++) {
		if (verbose) cout << " -- starting transparent clustering for track " << i << endl;
		
		// check if track is masked
		if (tracks[i].FakeTrack) {
			if (verbose) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
			continue;
		}
		
		// get event number for track
		if (verbose) cout << "Getting event number.. ";
		event = tracks[i].GetEventNumber();
		if (verbose) cout << " -> track " << i << " corresponds to event " << event << endl;
		
		// check if event number is valid
		if (event < 0) {
			cout << "Track " << i << " has no event number. Skipping this track.." << endl;
			continue;
		}
		
		// load data (apply offset)
		align->LoadData(tracks[i]);
		
		// read out x, y, z positions
		x_positions.clear();
		y_positions.clear();
		z_positions.clear();
		for (int det = 0; det < 4; det++) {
			x_positions.push_back(align->track_holder.GetD(det).GetX());
			y_positions.push_back(align->track_holder.GetD(det).GetY());
			z_positions.push_back(align->track_holder.GetD(det).GetZ());
		}
		
		// read out effictive diamond hit position
		eff_diamond_hit_position = align->track_holder.GetD(4).GetX();
		
		// fit track
		par.clear();
		align->LinTrackFit(z_positions, x_positions, par);
		if (verbose) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;
		
		// estimate hit position in diamond
//		diamond_z_position = align->track_holder.GetD(4).GetZ();
		diamond_hit_position = par[0] + par[1] * diamond_z_position;
		diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
		diamond_hit_channel = (int)diamond_hit_position;
		if (verbose) cout << "z position of diamond is " << diamond_z_position << endl;
		
		// difference between estimated and effective hit in diamond
		hit_diff = TMath::Abs(eff_diamond_hit_position - diamond_hit_position);
//		cout << "effective hit position in diamond:\t" << eff_diamond_hit_position << "\testimated position in diamond:\t" << diamond_hit_position << endl;
		
		//get event
		if (verbose) cout << "getting event " << event << ".." << endl;
		PedTree->GetEvent(event);
		if (verbose) cout << "event_number = " << event_number << endl;
		
		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int i = 0; i < 128; i++) {
			if (Dia_ADC[i]-Det_PedMean[8][i] > (Dia_ADC[eff_diamond_hit_channel]-Det_PedMean[8][eff_diamond_hit_channel])) {
				eff_diamond_hit_channel = i;
			}
		}
		histo_transparentclustering_hitdiff->Fill(TMath::Abs(eff_diamond_hit_channel - diamond_hit_position));
		cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;
		
		// cluster diamond channels around estimated hit position
//		cout << " --" << endl;
		if (verbose) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbose) cout << "Dia_ADC = " << Dia_ADC[diamond_hit_channel] << ",\tDet_PedMean = " << Det_PedMean[8][diamond_hit_channel] << endl;
		if (verbose) cout << "Collected charge in channel " << (int)diamond_hit_position << " of diamond: " << Dia_ADC[(int)diamond_hit_position]-Det_PedMean[8][(int)diamond_hit_position] << endl;
		
		int sign;
		
		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;
		
		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = Dia_ADC[diamond_hit_channel]-Det_PedMean[8][diamond_hit_channel];
		secondchannel_adc = Dia_ADC[diamond_secondhit_channel]-Det_PedMean[8][diamond_secondhit_channel];
		if (sign == 1) transp_eta = firstchannel_adc / (firstchannel_adc + secondchannel_adc);
		else transp_eta = secondchannel_adc / (firstchannel_adc + secondchannel_adc);
		histo_transparentclustering_eta->Fill(transp_eta);
		
		if (verbose) cout << "clusters for track " << i << ":" << endl;
		// loop over different cluster sizes
		for (int j = 0; j < 10; j++) {
			cluster_adc = 0;
			current_channel = diamond_hit_channel;
			if (verbose) cout << "selected channels for " << j+1 << " hit transparent cluster: ";
			current_sign = sign;
			// sum adc for n channel cluster
			for (int channel = 0; channel <= j; channel++) {
				current_channel = current_channel + current_sign * channel;
				current_sign = (-1) * current_sign;
				if (verbose) cout << current_channel;
				if (verbose) if (channel < j) cout << ", ";
				if (current_channel > 0 && current_channel < 128 && Dia_ADC[current_channel]-Det_PedMean[8][current_channel] > Di_Cluster_Hit_Factor*Det_PedWidth[8][current_channel])
					cluster_adc = cluster_adc + Dia_ADC[current_channel]-Det_PedMean[8][current_channel];
			}
			if (verbose) cout << endl;
			if (verbose) cout << "total charge of cluster: " << cluster_adc << endl;
			if (verbose) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
			histo_transparentclustering_landau[j]->Fill(cluster_adc);
		} // end loop over cluster sizes
	} // end loop over tracks
	
	// save histograms
	for (int i = 0; i < 10; i++) {
		SaveHistogram(histo_transparentclustering_landau[i]);
	}
	SaveHistogram(histo_transparentclustering_eta);
	SaveHistogram(histo_transparentclustering_hitdiff);
	
	PedFile->Close();
}

// --
void Clustering::EventMonitor(int CurrentEvent) {
	current_event = current_event - 1;
	
	//record large clusters
	string large_clusters_filename = plots_path + "large_clusters.txt";
	ofstream large_clusters_file(large_clusters_filename.c_str());
	large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
	uint nclusters, nhits;
	Cluster* current_cluster;
	
	//histograms of hit positions to find out if the eta correction helped
	TH1F* histo_hitpos[9][2];
	for(int det=0; det<9; det++) {
		for(int chip=0; chip<2; chip++) {
			ostringstream histotitle;
			histotitle << "Hit_Position_Interstrip_D" << det << "_chip" << chip;
			histo_hitpos[det][chip] = new TH1F(histotitle.str().c_str(), histotitle.str().c_str(), 101, -0.005, 1.005);
		}
	}
	
	if (!AlternativeClustering) ClusterEvent();
	else ClusterEventSeeds();
//	if(e%10000==0) clustered_event.Print();
//	BookHistograms();
	
	//record large clusters
	//if a track, check whether it's in the silicon fiducial region
	if(clustered_event.HasOneSiliconTrack()) {
		
		//loop over detectors to record large clusters to a txt file
		for(int det=0; det<9; det++) {
            nclusters = clustered_event.GetNClusters(det);
            //loop over all *good* clusters
            for(uint clus=0; clus<nclusters; clus++) {
				current_cluster = clustered_event.GetCluster(det,clus);
				nhits = current_cluster->GetNHits();
				if(nhits>8) {
					large_clusters_file << event_number << "\t" << det << "\t" << clus << "/" << nclusters << "\t" << nhits << "\t" << current_cluster->GetNSeeds();
					if(current_cluster->IsBadChannelCluster()) large_clusters_file << " (bad chan cluster)"; // skip bad clusters
					large_clusters_file << endl;
				}
            }
		}//end large cluster record
		
	}
	
	ostringstream histo_diamond_name;
	histo_diamond_name << "Event" << current_event << "DiaClusters";
	TH2F* histo_clusters[5];
	TH1F* histo_diamond;
	histo_diamond = new TH1F(histo_diamond_name.str().c_str(),histo_diamond_name.str().c_str(),256,-0.5,255.5);
	Float_t si_avg_x=0, si_avg_y=0;
	
	for (int det = 0; det < 5; det++) {
		ostringstream histo_clusters_name;
		histo_clusters_name << "Event" << current_event << "Det" << det << "Clusters";
		histo_clusters[det] = new TH2F(histo_clusters_name.str().c_str(),histo_clusters_name.str().c_str(),256,-0.5,255.5,256,-0.5,255.5);
		Float_t si_x=0, si_y=0;
		Float_t cluster_size = 0;
		//		int bin_x = 0, bin_y = 0, bin = 0;
		si_x = clustered_event.GetCluster(2*det, 0)->Get1stMoment();
		if (det < 4) {
			si_y = clustered_event.GetCluster(2*det+1, 0)->Get1stMoment();
			cluster_size = (clustered_event.GetCluster(2*det, 0)->GetNHits() + clustered_event.GetCluster(2*det+1, 0)->GetNHits()) / 2.;
			si_avg_x += si_x;
			si_avg_y += si_y;
		}
//		bin_x = histo_clusters[det]->
//		bin = histo_clusters[det]->GetBin(
		if (det < 4) {
			histo_clusters[det]->SetBinContent(si_x+1,si_y+1,cluster_size);
		}
		else {
			si_avg_x = si_avg_x / 4;
			si_avg_y = si_avg_y / 4;
			histo_clusters[det]->SetBinContent(si_x+1,si_avg_y+1,clustered_event.GetCluster(2*det, 0)->GetNHits());
			
			for (int j = 0; j < Det_NChannels[8]; j++) {
				if (Dia_ADC[j]-Det_PedMean[8][j] > Di_Cluster_Hit_Factor*Det_PedWidth[8][j] || true) {
//					histo_diamond->SetBinContent(j+1,Dia_ADC[j]-Det_PedMean[8][j]);
					histo_diamond->SetBinContent((int)Det_Channels[8][j],Dia_ADC[j]-Det_PedMean[8][j]);
				}
			}
			SaveHistogram(histo_diamond);
//			histo_diamond->SetBinContent(si_x+1,clustered_event.GetCluster(2*det, 0)->GetNHits());
		}

		SaveHistogram(histo_clusters[det]);
	}
	
	current_event++;
}

void Clustering::SetRunParameters(int reg, FidCutRegion current_region, bool MultiRegions) {
	si_avg_fidcut_xlow = current_region.GetValueXLow();
	si_avg_fidcut_xhigh = current_region.GetValueXHigh();
	si_avg_fidcut_ylow = current_region.GetValueYLow();
	si_avg_fidcut_yhigh = current_region.GetValueYHigh();
	
	// TODO: alternative plots path only for more than 1 region!!
	
	if (MultiRegions) {
		
		ostringstream plotspath;
		plotspath << plots_path.c_str() << "-FidCutRegion" << reg;
		
		//	plotspath << sys->pwd() << "/plots-" << RunNumber;
		//	if(RunDescription=="") plotspath << "/";
		//	else plotspath << "-" << RunDescription << "/";
		
		plotspath << "-FidCutRegion" << reg;
		
		plots_path = plotspath.str();
		png_file_char = plotspath.str();
		C_file_char = plotspath.str();
		root_file_char = plotspath.str();
		//make plots dir
		sys->mkdir(plots_path.c_str());
	}
}
