/*
 * AlignmentClass.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef ALIGNMENTCLASS_HH_
#define ALIGNMENTCLASS_HH_

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

#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TSettings.class.hh"
#include "TTransparentClustering.hh"

class AlignmentClass {
public:
	AlignmentClass(std::string fileName,UInt_t nEventNumber=0);
	virtual ~AlignmentClass();
	bool GetEvent(UInt_t nEvent);

	int Align(bool plots = 1, bool CutFakeTracksOn = false);
    void TransparentAnalysis(std::vector<TDiamondTrack> &tracks, std::vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);
	void TransparentClustering(std::vector<TDiamondTrack> &tracks, std::vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);
	void TransparentClustering2(std::vector<TDiamondTrack> &tracks, std::vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);

	void SetSettings(TSettings* settings);
	void SetPlotsPath(std::string plotspath);
	void SetTracks(std::vector<TDiamondTrack> tracks);
	void SetAlignment_tracks_mask(std::vector<bool> tracks_mask);
	void SetAlignment_tracks_fidcut_mask(std::vector<bool> alignment_tracks_fidcut_mask);
	void SetAlignment_tracks_fidcut(std::vector<TDiamondTrack> alignment_tracks_fidcut);
private:
	void createAlignmentSummary();

private:/**needed variables*/
	TDetectorAlignment* align;
	std::string plotsPath;
	std::vector<TDiamondTrack> alignment_tracks;
	std::vector<bool> alignment_tracks_mask;
	std::vector<TDiamondTrack> alignment_tracks_fidcut;
	std::vector<bool> alignment_tracks_fidcut_mask;

	//Event SETTINGS
	TSettings *settings;
	TADCEventReader *eventReader;

//	//clustering settings
//	Float_t Di_Cluster_Hit_Factor;
//    Float_t pulse_height_di_max;
//    Int_t pulse_height_num_bins;

	//histos
    TH1F* histo_transparentclustering_landau[10];
    TH1F* histo_transparentclustering_landau_mean;
    TH1F* histo_transparentclustering_eta;
	TH1F* histo_transparentclustering_hitdiff;
	TH2F* histo_transparentclustering_hitdiff_scatter;
	TH1F* histo_transparentclustering_2Channel_PulseHeight;
	TH1F* histo_transparentclustering_residuals[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
	TH2F* histo_transparentclustering_residuals_scatter[10];	// index: 0 distance to center of hit channel, 1 distance to charge weighted mean of closest two channels, 2 distance to charge weighted mean of closest three channels, ..
	TH1F* histo_transparentclustering_residuals_largest_hit[10];
	TH2F* histo_transparentclustering_residuals_largest_hit_scatter[10];
	TH1F* histo_transparentclustering_residuals_2largest_hits;
	TH2F* histo_transparentclustering_residuals_2largest_hits_scatter;
	TH1F* histo_transparentclustering_SNR_vs_channel;
	TH1F* histo_transparentclustering_chi2X;
	TH1F* histo_transparentclustering_chi2Y;
	
	Float_t dianoise_sigma[2];	//TODO: get from clustering tree!!

private: /* not used at the moment*/
	std::vector<Float_t> dia_offset;


private:/*constants*/
	int nAlignSteps;

private:/*see if needed*/
	HistogrammSaver *histSaver;
	HistogrammSaver *histSaverAlignment;
	HistogrammSaver *histSaverAlignmentFakedTracks;

private:
	std::string PedFileName;
	int verbosity;
};

#endif /* ALIGNMENTCLASS_HH_ */
