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
#include "TADCEventReader.hh"
#include "TSettings.class.hh"
using namespace std;

class AlignmentClass {
public:
	AlignmentClass(string fileName,UInt_t nEventNumber=0);
	virtual ~AlignmentClass();
	bool GetEvent(UInt_t nEvent);

	int Align(bool plots = 1, bool CutFakeTracksOn = false);
    void TransparentAnalysis(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);
	void TransparentClustering(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose = false);

	void SetSettings(TSettings* settings);
	void SetPlotsPath(string plotspath);
	void SetTracks(vector<TDiamondTrack> tracks);
	void SetAlignment_tracks_mask(vector<bool> tracks_mask);
	void SetAlignment_tracks_fidcut_mask(vector<bool> alignment_tracks_fidcut_mask);
	void SetAlignment_tracks_fidcut(vector<TDiamondTrack> alignment_tracks_fidcut);
private:
	void createAlignmentSummary();

private:/**needed variables*/
	TDetectorAlignment* align;

	vector<TDiamondTrack> alignment_tracks;
	vector<bool> alignment_tracks_mask;
	vector<TDiamondTrack> alignment_tracks_fidcut;
	vector<bool> alignment_tracks_fidcut_mask;

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

private: /* not used at the moment*/
	vector<Float_t> dia_offset;


private:/*constants*/
	int nAlignSteps;

private:/*see if needed*/
	HistogrammSaver *histSaver;
	HistogrammSaver *histSaverAlignment;
	HistogrammSaver *histSaverAlignmentFakedTracks;

private:
	string PedFileName;
	int verbosity;
};

#endif /* ALIGNMENTCLASS_HH_ */
