/*
 * TDetectorAlignment.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TDETECTORALIGNMENT_HH_
#define TDETECTORALIGNMENT_HH_
//C++ Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
//#include <deque>
#include <ctime> // seed the random number generator
#include <cstdlib> // random number generator
#include <sstream>

#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"

#include "TDiamondTrack.hh"
#include "TDetectorPlane.hh"


class TDetectorAlignment{

   public:
      TDetectorAlignment(string plots_path_string);
      TDetectorAlignment(string plots_path_string, vector<TDiamondTrack> &input_tracks, vector<bool> &input_tracks_mask);
      ~TDetectorAlignment() {};
      void SaveCanvas(TCanvas* canv, string filename);
      Double_t GetXOffset(Int_t plane) {return det_x_offset[plane];};
      Double_t GetYOffset(Int_t plane) {return det_y_offset[plane];};
      Double_t GetZOffset(Int_t plane) {return det_z_offset[plane];};
      Double_t GetPhiXOffset(Int_t plane) {return det_phix_offset[plane];};
      Double_t GetPhiYOffset(Int_t plane) {return det_phiy_offset[plane];};
      Double_t GetXResolution(Int_t plane) {return det_x_resolution[plane];};
      Double_t GetYResolution(Int_t plane) {return det_y_resolution[plane];};
      Double_t GetSiResolution();
      vector<Double_t> GetXOffsetHistory(Int_t plane) {return det_x_offset_history[plane];};
      vector<Double_t> GetYOffsetHistory(Int_t plane) {return det_y_offset_history[plane];};
      vector<Double_t> GetZOffsetHistory(Int_t plane) {return det_z_offset_history[plane];};
      vector<Double_t> GetPhiOffsetHistory(Int_t plane) {return det_phi_offset_history[plane];};
      void PlotXOffsetHistory(Int_t plane);
      void PlotYOffsetHistory(Int_t plane);
      void PlotZOffsetHistory(Int_t plane);
      void PlotPhiOffsetHistory(Int_t plane);
      void PositionPredictor(int subject_detector, int ref_detector1, int ref_detector2);
      void PositionPredictor(int subject_detector);
      void TrackFit(Double_t &predicted_x, Double_t &predicted_y);
      void TrackFit3(Double_t &predicted_x, Double_t &predicted_y);
      Double_t GetPredictedX() {return predicted_x;}
      Double_t GetPredictedY() {return predicted_y;}
      void LoadTracks(vector<TDiamondTrack> &input_tracks, vector<bool> &input_tracks_mask);
      void LoadData(TDiamondTrack track);
      void PlotAngularDistribution();
      void PlotCartesianDistribution();
      void AlignDetectorXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose = true, bool justcheck = false);
      void AlignDetectorXY(int subject_detector, bool verbose = true, bool justcheck = false);
      void AlignDetectorZ(Int_t subject_detector, bool graphs = false, bool verbose = true);
      void AlignDetectorZ(Int_t subject_detector, Int_t ref_detector1, Int_t ref_detector2, bool graphs = false, bool verbose = true);
      void CheckDetectorAlignmentXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose = true);
      void CheckDetectorAlignmentXY(int subject_detector, bool verbose = true);
      void CheckDetectorAlignmentXYPlots(int subject_detector, int ref_detector1, int ref_detector2, string &histo_title);
	void CheckDetectorAlignmentXYPlots(int subject_detector, string &histo_title);
	void CutFakeTracks(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, Float_t alignment_chi2 = 9999., bool CutFakeTracksOn = false, bool verbose = false);
	Float_t LinTrackFit(vector<Float_t> X, vector<Float_t> Y, vector<Float_t> &par, Float_t res = 9999.);

   protected:
      TDetectorPlane D0;
      TDetectorPlane D1;
      TDetectorPlane D2;
      TDetectorPlane D3;

   public:
      //store tracks and masks for determining alignment constants
      vector<TDiamondTrack> track_storage;
      vector<bool> track_mask_storage;
      //TODO: instead of having multiple copies of tracks, why not store pointer to list in stored in Clustering::tracks
      //TODO: instead of having tracks and tracks_fidcut, store boolean fidcut mask

   public:
      //temp data for calculations
      TDiamondTrack track_holder;
      Double_t x_array[3];
      Double_t y_array[3];
      Double_t z_array[3];
      Double_t subject_z;
      Double_t predicted_x;
      Double_t predicted_y;
      Int_t nDetectors;


      //store global offsets here
      Double_t det_x_offset[6];
      Double_t det_y_offset[6];
      Double_t det_z_offset[6];
      Double_t det_phix_offset[6];
      Double_t det_phiy_offset[6];

      //store resolutions here
      Double_t det_x_resolution[6];
      Double_t det_y_resolution[6];

      //store reconstructed offsets here
      vector<Double_t> det_x_offset_history[6];
      vector<Double_t> det_y_offset_history[6];
      vector<Double_t> det_z_offset_history[6];
      vector<Double_t> det_phi_offset_history[6];

   public:
      /*
      TH1F residualsX;
      TH1F residualsY;
      TH2F residualsXY;
      TH2F residualsXvsY;
      TH2F residualsYvsX;
      */
      Double_t residualsXmean, residualsYmean, residualsXrms, residualsYrms;
      string plots_path;
      int SaveAllFilesSwitch, ClosePlotsOnSave, SaveAllRootFilesSwitch;

   private:
      int verbosity;

};
#endif /* TDETECTORALIGNMENT_HH_ */
