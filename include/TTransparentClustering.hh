/*
 * TTransparentClustering.hh
 *
 *  Created on: 04.11.2011
 *      Author: bachmair
 */

#ifndef TTRANSPARENTCLUSTERING_HH_
#define TTRANSPARENTCLUSTERING_HH_

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

#include "TSettings.class.hh"
#include "TDiamondTrack.hh"
#include "TDetectorAlignment.hh"
#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"

using namespace std;

class TTransparentClustering {
public:
	TTransparentClustering(string PedFileName);
	virtual ~TTransparentClustering();
	void MakeTransparentClustering();
    Double_t getDiamondPhiOffset() const;
    Double_t getDiamondPhiYOffset() const;
    Double_t getDiamondXOffset() const;
    Double_t getDiamondYOffset() const;
    Double_t getDiamondZPosition() const;
    void setDiamondPhiOffset(Double_t diamondPhiOffset);
    void setDiamondPhiYOffset(Double_t diamondPhiYOffset);
    void setDiamondXOffset(Double_t diamondXOffset);
    void setDiamondYOffset(Double_t diamondYOffset);
    void setDiamondZPosition(Double_t diamondZPosition);
    void setTracks( vector<TDiamondTrack> tracks);
    Float_t getSiRes() const;
    void setSiRes(Float_t siRes);
    TDetectorAlignment *getCurrentAlignment() const;
    void setCurrentAlignment(TDetectorAlignment *currentAlignment);


	void SetSettings(TSettings* settings);
	void SetPlotsPath(std::string plotspath);
    vector<TDiamondTrack> getTracks() const;
private:
    void initHistogramms();
private:
    TSettings *settings;
	HistogrammSaver *histSaver;
	TADCEventReader *eventReader;
    Double_t diamond_z_position;
    Double_t diamond_phi_y_offset;
    Double_t diamond_phi_offset;
    Double_t diamond_y_offset;
    Double_t diamond_x_offset;
    Float_t SiRes;
    TDetectorAlignment *currentAlignment;
    vector<TDiamondTrack> tracks;
	Float_t dianoise_sigma[2];

private:
    Int_t verbosity;

private:
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
};

#endif /* TTRANSPARENTCLUSTERING_HH_ */
