/*
 * AlignmentClass.hh
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#ifndef ALIGNMENTCLASS_HH_
#define ALIGNMENTCLASS_HH_

#include <vector>
#include <string>
#include "TTree.h"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
using namespace std;

class AlignmentClass {
public:
	AlignmentClass();
	virtual ~AlignmentClass();
	int Align(bool plots = 1, bool CutFakeTracksOn = false);
    void TransparentAnalysis(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);
	void TransparentClustering(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);


private:/**needed variables*/
	vector<TDiamondTrack> tracks;

	vector<bool> alignment_tracks_mask;
	vector<TDiamondTrack> alignment_tracks_fidcut;
	vector<bool> alignment_tracks_fidcut_mask;

	//Event SETTINGS
	//TODO: Put these variables in an RawEventClass
	UInt_t run_number;
	UInt_t event_number;
	Float_t store_threshold;
	bool CMNEvent_flag;
	bool ZeroDivisorEvent_flag;
	UInt_t Det_NChannels[9];
	UChar_t Det_Channels[9][256];
	UChar_t Det_ADC[8][256];
	UShort_t Dia_ADC[256];
	Float_t Det_PedMean[9][256];
	Float_t Det_PedWidth[9][256];

	//pathes
	string plots_path;
	string pedfile_path;
	//clustering settings
	Float_t Di_Cluster_Hit_Factor;
    Float_t pulse_height_di_max;
    Int_t pulse_height_num_bins;

	//histos
    TH1F* histo_transparentclustering_landau[10];
    TH1F* histo_transparentclustering_landau_mean;
    TH1F* histo_transparentclustering_eta;
	TH1F* histo_transparentclustering_hitdiff;
	TH2F* histo_transparentclustering_hitdiff_scatter;
	TH1F* histo_transparentclustering_2Channel_PulseHeight;

private: /* not used at the moment*/
	vector<Float_t> dia_offset;

private:/*variable from settings:*/
	Float_t alignment_chi2;

private:/*constants*/
	int nAlignSteps;

private:/*see if needed*/
	TSystem* sys;
	TFile *PedFile;
	TTree *PedTree;
	HistogrammSaver *histSaver;
};

#endif /* ALIGNMENTCLASS_HH_ */
