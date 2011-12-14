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
#include <iomanip>
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
	TDetectorAlignment();
      TDetectorAlignment(std::string plots_path_string);
      TDetectorAlignment(std::string plots_path_string, std::vector<TDiamondTrack> &input_tracks, std::vector<bool> &input_tracks_mask);
      ~TDetectorAlignment() {};
      void PrintDiamondTracks();
      void SaveCanvas(TCanvas* canv, std::string filename);
      Double_t GetXOffset(Int_t plane) {return det_x_offset[plane];};
      Double_t GetYOffset(Int_t plane) {return det_y_offset[plane];};
      Double_t GetZOffset(Int_t plane) {return det_z_offset[plane];};
      void SetZOffset(Int_t plane,Float_t zOffset) {if(plane<6)det_z_offset[plane]=zOffset;};
      void AddToXOffset(UInt_t plane, Float_t addXOffset){if(plane<6)det_x_offset[plane]+=addXOffset;}
      void AddToYOffset(UInt_t plane, Float_t addYOffset){if(plane<6)det_y_offset[plane]+=addYOffset;}
      Double_t GetPhiXOffset(Int_t plane) {return det_phix_offset[plane];};
      void AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset){if(plane<6)det_phix_offset[plane]+=addPhiXOffset;}
      Double_t GetPhiYOffset(Int_t plane) {return det_phiy_offset[plane];};
      void AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset){if(plane<6)det_phiy_offset[plane]+=addPhiYOffset;}
      Double_t GetXResolution(Int_t plane) {return det_x_resolution[plane];};
      Double_t GetYResolution(Int_t plane) {return det_y_resolution[plane];};
      Double_t GetSiResolution();
      std::vector<Double_t> GetXOffsetHistory(Int_t plane) {return det_x_offset_history[plane];};
      std::vector<Double_t> GetYOffsetHistory(Int_t plane) {return det_y_offset_history[plane];};
      std::vector<Double_t> GetZOffsetHistory(Int_t plane) {return det_z_offset_history[plane];};
      std::vector<Double_t> GetPhiOffsetHistory(Int_t plane) {return det_phi_offset_history[plane];};
      void PositionPredictor(int subject_detector, int ref_detector1, int ref_detector2);
      void PositionPredictor(int subject_detector);
      void TrackFit(Double_t &predicted_x, Double_t &predicted_y);
      void TrackFit3(Double_t &predicted_x, Double_t &predicted_y);
      Double_t GetPredictedX() {return predicted_x;}
      Double_t GetPredictedY() {return predicted_y;}
      void LoadTracks(std::vector<TDiamondTrack> &input_tracks, std::vector<bool> &input_tracks_mask);
      void LoadData(TDiamondTrack track);
      void PlotAngularDistribution();
      void PlotCartesianDistribution();
      void AlignDetectorXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose = true, bool justcheck = false);
      void AlignDetectorXY(int subject_detector, bool verbose = true, bool justcheck = false);
      void AlignDetectorZ(Int_t subject_detector, bool graphs = false, bool verbose = true);
      void AlignDetectorZ(Int_t subject_detector, Int_t ref_detector1, Int_t ref_detector2, bool graphs = false, bool verbose = true);
      void CheckDetectorAlignmentXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose = true);
      void CheckDetectorAlignmentXY(int subject_detector, bool verbose = true);
      void CheckDetectorAlignmentXYPlots(int subject_detector, int ref_detector1, int ref_detector2, std::string &histo_title);
	void CheckDetectorAlignmentXYPlots(int subject_detector, std::string &histo_title);
	void CutFakeTracks(std::vector<TDiamondTrack> &tracks, std::vector<bool> &tracks_mask, Float_t alignment_chi2 = 9999., bool CutFakeTracksOn = false, bool verbose = false);
	Float_t LinTrackFit(std::vector<Float_t> X, std::vector<Float_t> Y, std::vector<Float_t> &par, Float_t res = 9999.);
    int getVerbosity() const;
    void setVerbosity(int verbosity);

   protected:
      TDetectorPlane D0;
      TDetectorPlane D1;
      TDetectorPlane D2;
      TDetectorPlane D3;

   public:
      //store tracks and masks for determining alignment constants
      std::vector<TDiamondTrack> track_storage;
      std::vector<bool> track_mask_storage;
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
      std:: vector<Double_t> det_x_offset_history[6];
      std::vector<Double_t> det_y_offset_history[6];
      std::vector<Double_t> det_z_offset_history[6];
      std::vector<Double_t> det_phi_offset_history[6];

   public:
      /*
      TH1F residualsX;
      TH1F residualsY;
      TH2F residualsXY;
      TH2F residualsXvsY;
      TH2F residualsYvsX;
      */
      Double_t residualsXmean, residualsYmean, residualsXrms, residualsYrms;
      std::string plots_path;
      int SaveAllFilesSwitch, ClosePlotsOnSave, SaveAllRootFilesSwitch;

   private:
      int verbosity;

};
#endif /* TDETECTORALIGNMENT_HH_ */
