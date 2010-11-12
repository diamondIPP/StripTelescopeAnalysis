//Classes used for silicon and diamond alignment
//based on Alignment_Monte_Carlo_j15a.cpp
//crapload of changes to class structures
//fixed segafaulting in phi and z alignment by moving plot generation to its own member fcn
//2010-09-18 Plucked from dev-2010-04-30-eta_correction for adapatation in dev-2010-09-18-alignment
//2010-09
//           TODO: Several improvements would be tracking and aligning x and y planes separately and of course tearing the damn telescope apart to pin down the geometry

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

typedef unsigned int uint;

using namespace std;
//using namespace TMath;


Int_t RoundValue(Double_t value)
{
   Double_t diff = value - (int)value;
   if(value>=0)
   {
      if(diff>=0.5)
      {
         value = (int)value + 1;
      }
      if(diff<0.5)
      {
         value = (int)value;
      }
   }
   if(value<0)
   {
      if(fabs(diff)>=0.5)
      {
         value = (int)value - 1;
      }
      if(fabs(diff)<0.5)
      {
         value = (int)value;
      }
   }
   
   return (Int_t)value;
}

//_________________________________________________________________________________________________
class TDetectorPlane {

   public:
      TDetectorPlane() {};
      ~TDetectorPlane() {};
      void SetX(Float_t x) {X_position = x;};
      void SetY(Float_t y) {Y_position = y;};
      void SetZ(Float_t z) {Z_position = z;};
      Float_t GetX() const {return X_position;}
      Float_t GetY() const {return Y_position;}
      Float_t GetZ() const {return Z_position;}

   private:
      Float_t X_position;
      Float_t Y_position;
      Float_t Z_position;

};

//_________________________________________________________________________________________________
class TDiamondTrack{

   public:
      TDiamondTrack() {};
      TDiamondTrack(TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3);
      TDiamondTrack(TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3, TDetectorPlane Dia);
      ~TDiamondTrack();
      TDetectorPlane const GetD0() { return D0;};
      TDetectorPlane const GetD1() { return D1;};
      TDetectorPlane const GetD2() { return D2;};
      TDetectorPlane const GetD3() { return D3;};
      TDetectorPlane const GetD4() { return D4;};
      TDetectorPlane const GetDia() { return D4;};
      TDetectorPlane const GetD(Int_t detector);
      Float_t GetDetectorHitPosition(Int_t det /*det = 0 to 8*/);
      void SetDetectorHitPosition(Int_t det /*det = 0 to 8*/, Float_t pos);

   protected:
      TDetectorPlane D0;
      TDetectorPlane D1;
      TDetectorPlane D2;
      TDetectorPlane D3;
      TDetectorPlane D4;
};

TDiamondTrack::TDiamondTrack(TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3)
{
   D0 = Det0;
   D1 = Det1;
   D2 = Det2;
   D3 = Det3;
}

TDiamondTrack::TDiamondTrack(TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3, TDetectorPlane Dia)
{
   D0 = Det0;
   D1 = Det1;
   D2 = Det2;
   D3 = Det3;
   D4 = Dia;
}


TDiamondTrack::~TDiamondTrack() {}

TDetectorPlane const TDiamondTrack::GetD(Int_t detector)
{
   if(detector == 0) return D0;
   if(detector == 1) return D1;
   if(detector == 2) return D2;
   if(detector == 3) return D3;
   if(detector == 4) return D4;
   else {
      cout<<"TDiamondTrack::GetD("<<detector<<"): detector "<<detector<<" is not a valid detector"<<endl;
      return TDetectorPlane();
   }
}

Float_t TDiamondTrack::GetDetectorHitPosition(Int_t det /*det = 0 to 8*/) {
   switch(det) {
      case 0: return D0.GetX(); break;
      case 1: return D0.GetY(); break;
      case 2: return D1.GetX(); break;
      case 3: return D1.GetY(); break;
      case 4: return D2.GetX(); break;
      case 5: return D2.GetY(); break;
      case 6: return D3.GetX(); break;
      case 7: return D3.GetY(); break;
      case 8: return D4.GetX(); break;
      case 9: return D4.GetY(); break;
   }
   return -1;
}

void TDiamondTrack::SetDetectorHitPosition(Int_t det /*det = 0 to 8*/, Float_t pos) {
   switch(det) {
      case 0: D0.SetX(pos); break;
      case 1: D0.SetY(pos); break;
      case 2: D1.SetX(pos); break;
      case 3: D1.SetY(pos); break;
      case 4: D2.SetX(pos); break;
      case 5: D2.SetY(pos); break;
      case 6: D3.SetX(pos); break;
      case 7: D3.SetY(pos); break;
      case 8: D4.SetX(pos); break;
      case 9: D4.SetY(pos); break;
      }
}

//_________________________________________________________________________________________________
class TDetectorAlignment{

   public:
      TDetectorAlignment(string plots_path_string);
      TDetectorAlignment(string plots_path_string, vector<TDiamondTrack> &input_tracks);
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
      void LoadTracks(vector<TDiamondTrack> &input_tracks) {track_storage.clear(); track_storage = input_tracks; cout<<"TDetectAlignment::LoadTracks: "<<track_storage.size()<<" tracks loaded"<<endl;};
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
	void CutFakeTracks();
	Float_t *LinTrackFit(vector<Float_t> X, vector<Float_t> Y, Float_t res);

   protected:
      TDetectorPlane D0;
      TDetectorPlane D1;
      TDetectorPlane D2;
      TDetectorPlane D3;

   public:
      //store tracks
      vector<TDiamondTrack> track_storage;
      
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
      int SaveAllFilesSwitch, ClosePlotsOnSave;
      
};

TDetectorAlignment::TDetectorAlignment(string plots_path_string) {
   
   nDetectors = 5;
   /*residualsX = TH1F("residualsX","residualsX",10000,-100,100);
   residualsY = TH1F("residualsY","residualsY",10000,-100,100);
   residualsXY = TH2F("residualsXY","residualsXY",10000,-100,100,10000,-100,100);
   residualsXvsY = TH2F("residualsXvsY","residualsXvsY",256,-0.5,255.5,10000,-100,100);
   residualsYvsX = TH2F("residualsYvsX","residualsYvsX",256,-0.5,255.5,10000,-100,100);
   */
   
   for(Int_t i=0; i<5; i++) {
      det_x_offset[i] = 0;
      det_y_offset[i] = 0;
      det_z_offset[i] = 0;
      det_phix_offset[i] = 0;
      det_phiy_offset[i] = 0;
      
      det_x_resolution[i] = -1;
      det_y_resolution[i] = -1;
   }
   
   plots_path = plots_path_string;
   SaveAllFilesSwitch = 1;
   ClosePlotsOnSave = 1;
}

TDetectorAlignment::TDetectorAlignment(string plots_path_string, vector<TDiamondTrack> &input_tracks) {
   
   nDetectors = 5;
   /*
   residualsX = TH1F("residualsX","residualsX",10000,-100,100);
   residualsY = TH1F("residualsY","residualsY",10000,-100,100);
   residualsXY = TH2F("residualsXY","residualsXY",10000,-100,100,10000,-100,100);
   residualsXvsY = TH2F("residualsXvsY","residualsXvsY",256,-0.5,255.5,10000,-100,100);
   residualsYvsX = TH2F("residualsYvsX","residualsYvsX",256,-0.5,255.5,10000,-100,100);
   */
   for(Int_t i=0; i<5; i++) {
      det_x_offset[i] = 0;
      det_y_offset[i] = 0;
      det_z_offset[i] = 0;
      det_phix_offset[i] = 0;
      det_phiy_offset[i] = 0;
      
      det_x_resolution[i] = -1;
      det_y_resolution[i] = -1;
   }
   
   LoadTracks(input_tracks);
   
   plots_path = plots_path_string;
   SaveAllFilesSwitch = 1;
   ClosePlotsOnSave = 1;
}

void TDetectorAlignment::SaveCanvas(TCanvas* canv, string filename) {
   ostringstream plot_filename;
   plot_filename << plots_path << filename;
   canv->Print(plot_filename.str().c_str());
}

Double_t TDetectorAlignment::GetSiResolution() {
   Int_t deta[] = {1,0,0,0}, detb[] = {0,1,2,3}, detc[] = {3,3,3,2};
   Double_t zab, zbc, ref_residual_width, numerator=0, denominator=0;
   for(int det=0; det<8; det++) {
      if(det%2) ref_residual_width = det_y_resolution[det/2];
      else ref_residual_width = det_x_resolution[det/2];
      zab = TMath::Abs(track_storage[0].GetD(deta[det/2]).GetZ() - track_storage[0].GetD(detb[det/2]).GetZ());
      zbc = TMath::Abs(track_storage[0].GetD(detb[det/2]).GetZ() - track_storage[0].GetD(detc[det/2]).GetZ());
      numerator += ref_residual_width*ref_residual_width * (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc));
      denominator += (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc)) * (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc));
   }
   
   return TMath::Sqrt(numerator/denominator);
}

void TDetectorAlignment::LoadData(TDiamondTrack track)
{
   TDetectorPlane det_buffer[5];
   
   for(Int_t det=0; det<5; det++) det_buffer[det] = track.GetD(det);
   
   /*
   //rotate planes; the 128's shift the axis of the rotation to go through the center of each plane
   for(Int_t det=0; det<nDetectors; det++) {
      det_buffer[det].SetX((det_buffer[det].GetX()-128)*TMath::Cos(det_phi_offset[det]) -
            (det_buffer[det].GetY()-128)*TMath::Sin(-det_phi_offset[det]) + 128);
      det_buffer[det].SetY((det_buffer[det].GetX()-128)*TMath::Sin(-det_phi_offset[det]) +
            (det_buffer[det].GetY()-128)*TMath::Cos(det_phi_offset[det]) + 128);
   }
   */
   
   //translate planes
   for(Int_t det=0; det<5; det++) {
      det_buffer[det].SetX(det_buffer[det].GetX()-det_x_offset[det]);
      det_buffer[det].SetY(det_buffer[det].GetY()-det_y_offset[det]);
      det_buffer[det].SetZ(det_buffer[det].GetZ()-det_z_offset[det]);
   }
   
   //rotate planes; the 128's shift the axis of the rotation to go through the center of each plane
   for(Int_t det=0; det<4; det++) {
      det_buffer[det].SetX((det_buffer[det].GetX()-128)*TMath::Cos(det_phix_offset[det]) -
            (det_buffer[det].GetY()-128)*TMath::Sin(-det_phix_offset[det]) + 128);
      det_buffer[det].SetY((det_buffer[det].GetX()-128)*TMath::Sin(-det_phiy_offset[det]) +
            (det_buffer[det].GetY()-128)*TMath::Cos(det_phiy_offset[det]) + 128);
   }
   
   //rotate diamond differently since no y data and different cenetr; the 64's shift the axis of the rotation to go through the center of the diamond plane
   /*for(Int_t det=4; det<5; det++) {
      det_buffer[det].SetX((det_buffer[det].GetX()-64)*TMath::Cos(det_phix_offset[det]) + 64);
   }*/
   
   //rotate diamond differently since no y data and different cenetr; the 64's shift the axis of the rotation to go through the center of the diamond plane
   for(Int_t det=4; det<5; det++) {
      det_buffer[det].SetX((det_buffer[det].GetX()-64)*TMath::Cos(det_phix_offset[det]) -
            (det_buffer[det].GetY()-64)*TMath::Sin(-det_phix_offset[det]) + 64);
      det_buffer[det].SetY((det_buffer[det].GetX()-64)*TMath::Sin(-det_phiy_offset[det]) +
            (det_buffer[det].GetY()-64)*TMath::Cos(det_phiy_offset[det]) + 64);
   }
   
   //load track
   TDiamondTrack temptrack(det_buffer[0],det_buffer[1],det_buffer[2],det_buffer[3],det_buffer[4]);
   track_holder = temptrack;
   
   /* //diagnostics
   cout<<det_x_offset[1]<<"\t"<<det_y_offset[1]<<"\t"<<det_z_offset[1]<<"\t"<<det_phi_offset[1]<<endl;
   cout<<track.GetD(0).GetX()<<"\t"<<track.GetD(0).GetY()<<"\t"<<track.GetD(0).GetZ()<<endl;
   cout<<track.GetD(1).GetX()<<"\t"<<track.GetD(1).GetY()<<"\t"<<track.GetD(1).GetZ()<<endl;
   cout<<det_buffer[0].GetX()<<"\t"<<det_buffer[0].GetY()<<"\t"<<det_buffer[0].GetZ()<<endl;
   cout<<det_buffer[1].GetX()<<"\t"<<det_buffer[1].GetY()<<"\t"<<det_buffer[1].GetZ()<<endl;
   cout<<temptrack.GetD(0).GetX()<<"\t"<<temptrack.GetD(0).GetY()<<"\t"<<temptrack.GetD(0).GetZ()<<endl;
   cout<<temptrack.GetD(1).GetX()<<"\t"<<temptrack.GetD(1).GetY()<<"\t"<<temptrack.GetD(1).GetZ()<<endl;
   cout<<track_holder.GetD(0).GetX()<<"\t"<<track_holder.GetD(0).GetY()<<"\t"<<track_holder.GetD(0).GetZ()<<endl;
   cout<<track_holder.GetD(1).GetX()<<"\t"<<track_holder.GetD(1).GetY()<<"\t"<<track_holder.GetD(1).GetZ()<<endl;
   */
}

void TDetectorAlignment::PositionPredictor(int subject_detector, int ref_detector1, int ref_detector2)
{
   x_array[0] = track_holder.GetD(ref_detector1).GetX();
   x_array[1] = track_holder.GetD(ref_detector2).GetX();
   
   y_array[0] = track_holder.GetD(ref_detector1).GetY();
   y_array[1] = track_holder.GetD(ref_detector2).GetY();
   
   z_array[0] = track_holder.GetD(ref_detector1).GetZ();
   z_array[1] = track_holder.GetD(ref_detector2).GetZ();

   subject_z = track_holder.GetD(subject_detector).GetZ();
   
   TrackFit(predicted_x, predicted_y);
   
}

void TDetectorAlignment::PositionPredictor(int subject_detector)
{
   int index=0;
   for(int det=0; det<4; det++){
      if(det!=subject_detector&&det!=4){
         x_array[index] = track_holder.GetD(det).GetX();
         y_array[index] = track_holder.GetD(det).GetY();
         z_array[index] = track_holder.GetD(det).GetZ();
         index++;
      }
   }

   subject_z = track_holder.GetD(subject_detector).GetZ();
   
   TrackFit3(predicted_x, predicted_y);
}

void TDetectorAlignment::TrackFit(Double_t &predicted_x, Double_t &predicted_y)
{
   Double_t slope = (x_array[1]-x_array[0])/(z_array[1]-z_array[0]);
   Double_t intercept = x_array[0]-slope*z_array[0];
   predicted_x = (subject_z*slope)+intercept;
   
   slope = (y_array[1]-y_array[0])/(z_array[1]-z_array[0]);
   intercept = y_array[0]-slope*z_array[0];
   predicted_y = (subject_z*slope)+intercept;
}

void TDetectorAlignment::TrackFit3(Double_t &predicted_x, Double_t &predicted_y)
{
   TF1 *linear = new TF1("linear","([0]*x)+[1]",-40,40); // range is guess at how badly aligned we can go
   TGraph *trackx = new TGraph(3,z_array,x_array);
   trackx->Fit("linear","Q");
   Double_t slope = linear->GetParameter(0);
   //cout << "xslope = " << slope << endl;
   Double_t intercept = linear->GetParameter(1);
   //cout << "xintercept = " << intercept << endl;
   predicted_x = (subject_z*slope)+intercept;

   //i don't think this works...
   //TCanvas *can= new TCanvas("can", "can",800,600,600,600);
   //trackx->Draw("A*");
   
   TGraph *tracky = new TGraph(3,z_array,y_array);
   tracky->Fit("linear","Q");
   slope = linear->GetParameter(0);
   //cout << "yslope = " << slope << endl;
   intercept= linear->GetParameter(1);
   //cout << "yintercept = " << intercept << endl;
   predicted_y = (subject_z*slope)+intercept;

   delete trackx;
   delete tracky;
   delete linear;
}

void TDetectorAlignment::PlotAngularDistribution()
{
   Double_t deltaX = 0;
   Double_t deltaY = 0;
   Double_t deltaZ = track_storage[0].GetD3().GetZ() - track_storage[0].GetD0().GetZ();
   cout << "For first track, deltaZ = " << deltaZ << endl;
   
   Double_t pitch = 0.005; // in cm
   
   Double_t pol_low = 0;
   Double_t pol_hi = TMath::ATan(pitch*512/deltaZ);
   Double_t pol_bin_width = 0.00002;
   
   Double_t azi_low = -TMath::ASin(1)*2;
   Double_t azi_hi = TMath::ASin(1)*2;
   Double_t azi_bin_width = 0.005;
   
   TH1F polar_spread("polar_spread", "polar_spread", (Int_t)(TMath::Abs(pol_hi-pol_low)/pol_bin_width), pol_low, pol_hi);
   TH1F azimuthal_spread("azimuthal_spread", "azimuthal_spread", (Int_t)(TMath::Abs(azi_hi-azi_low)/azi_bin_width), azi_low, azi_hi);
   
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
   {
      TDiamondTrack track = track_storage[t]; // note: (*track_storage)[t] is equivalent to track_storage[0][t]
      TDetectorPlane D0_buf = track.GetD0();
      TDetectorPlane D3_buf = track.GetD3();

      //calculate azimuthal and polar angles of track (using points from D0 and D3 instead of slopes)
      deltaX = D3_buf.GetX() - D0_buf.GetX();
      deltaY = D3_buf.GetY() - D0_buf.GetY();
      deltaZ = D3_buf.GetZ() - D0_buf.GetZ();
      
      polar_spread.Fill(TMath::ATan(pitch*TMath::Sqrt(deltaX*deltaX + deltaY*deltaY)/deltaZ));
      azimuthal_spread.Fill(TMath::ATan(deltaY/deltaX));
   }
   
   cout << "Tracks have mean polar angle " << polar_spread.GetMean() << " and RMS of " << polar_spread.GetRMS() << endl;
   cout << "Tracks have mean azimuthal angle " << azimuthal_spread.GetMean() << " and RMS of " << azimuthal_spread.GetRMS() << endl;
   
   TCanvas *tempcanangle = new TCanvas("anglestempcanv","anglestempcanv",800,600);
   polar_spread.GetXaxis()->SetRangeUser(polar_spread.GetMean()-3*polar_spread.GetRMS(),polar_spread.GetMean()+3*polar_spread.GetRMS());
   polar_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracksangles_polar.png");
   //SaveCanvas(tempcanangle, "alignment_tracksangles_polar.C"); //for many tracks this takes a long time to plot
   delete tempcanangle;
   
   tempcanangle = new TCanvas("anglestempcanv","anglestempcanv",800,600);
   azimuthal_spread.GetXaxis()->SetRangeUser(azimuthal_spread.GetMean()-3*azimuthal_spread.GetRMS(),azimuthal_spread.GetMean()+3*azimuthal_spread.GetRMS());
   azimuthal_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracksangles_azimuth.png");
   delete tempcanangle;
}

void TDetectorAlignment::PlotCartesianDistribution()
{
   Double_t deltaX = 0;
   Double_t deltaY = 0;
   Double_t deltaZ = track_storage[0].GetD3().GetZ() - track_storage[0].GetD0().GetZ();
   cout << "For first track, deltaZ = " << deltaZ << endl;
   
   Double_t xint_low = -128;
   Double_t xint_hi = 128;
   Double_t xint_bin_width = 1/3.;//0.1;
   Double_t yint_low = -128;
   Double_t yint_hi = 128;
   Double_t yint_bin_width = 1/3.;//0.1;
   
   Double_t xslope_low = -1;
   Double_t xslope_hi = 1;
   Double_t xslope_bin_width = 0.01;
   Double_t yslope_low = -1;
   Double_t yslope_hi = 1;
   Double_t yslope_bin_width = 0.01;
   
   
   TH1F xint_spread("xintercept_spread", "xintercept_spread", (Int_t)(TMath::Abs(xint_hi-xint_low)/xint_bin_width), xint_low, xint_hi);
   TH1F yint_spread("yintercept_spread", "yintercept_spread", (Int_t)(TMath::Abs(yint_hi-yint_low)/yint_bin_width), yint_low, yint_hi);
   TH1F xslope_spread("xslope_spread", "xslope_spread", (Int_t)(TMath::Abs(xslope_hi-xslope_low)/xslope_bin_width), xslope_low, xslope_hi);
   TH1F yslope_spread("yslope_spread", "yslope_spread", (Int_t)(TMath::Abs(yslope_hi-yslope_low)/yslope_bin_width), yslope_low, yslope_hi);
   
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
   {
      TDiamondTrack track = track_storage[t]; // note: (*track_storage)[t] is equivalent to track_storage[0][t]
      TDetectorPlane D0_buf = track.GetD0();
      TDetectorPlane D3_buf = track.GetD3();

      //calculate slopes of track (using points from D0 and D3 instead of slopes)
      deltaX = D3_buf.GetX() - D0_buf.GetX();
      deltaY = D3_buf.GetY() - D0_buf.GetY();
      deltaZ = D3_buf.GetZ() - D0_buf.GetZ();
      
      xint_spread.Fill(D0_buf.GetX());
      yint_spread.Fill(D0_buf.GetY());
      xslope_spread.Fill(deltaX/deltaZ);
      yslope_spread.Fill(deltaY/deltaZ);
   }
   
   cout << "Tracks have mean x-intercept of " << xint_spread.GetMean() << " and RMS of " << xint_spread.GetRMS() << endl;
   cout << "Tracks have mean y-intercept of " << yint_spread.GetMean() << " and RMS of " << yint_spread.GetRMS() << endl;
   cout << "Tracks have mean x-slope of " << xslope_spread.GetMean() << " and RMS of " << xslope_spread.GetRMS() << endl;
   cout << "Tracks have mean y-slope of " << yslope_spread.GetMean() << " and RMS of " << yslope_spread.GetRMS() << endl;
   
   TCanvas* tempcanangle = new TCanvas("trackdistributiontempcanv","trackdistributiontempcanv",800,600);
   xint_spread.GetXaxis()->SetRangeUser(xint_spread.GetMean()-3*xint_spread.GetRMS(),xint_spread.GetMean()+3*xint_spread.GetRMS());
   xint_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracks_xintercept.png");
   //SaveCanvas(tempcanangle, "alignment_tracks_xintercept.C");
   delete tempcanangle;
   
   tempcanangle = new TCanvas("trackdistributiontempcanv","trackdistributiontempcanv",800,600);
   yint_spread.GetXaxis()->SetRangeUser(yint_spread.GetMean()-3*yint_spread.GetRMS(),yint_spread.GetMean()+3*yint_spread.GetRMS());
   yint_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracks_yintercept.png");
   //SaveCanvas(tempcanangle, "alignment_tracks_yintercept.C");
   delete tempcanangle;
   
   tempcanangle = new TCanvas("trackdistributiontempcanv","trackdistributiontempcanv",800,600);
   xslope_spread.GetXaxis()->SetRangeUser(xslope_spread.GetMean()-3*xslope_spread.GetRMS(),xslope_spread.GetMean()+3*xslope_spread.GetRMS());
   xslope_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracks_xslope.png");
   //SaveCanvas(tempcanangle, "alignment_tracks_xslope.C");
   delete tempcanangle;
   
   tempcanangle = new TCanvas("trackdistributiontempcanv","trackdistributiontempcanv",800,600);
   yslope_spread.GetXaxis()->SetRangeUser(yslope_spread.GetMean()-3*yslope_spread.GetRMS(),yslope_spread.GetMean()+3*yslope_spread.GetRMS());
   yslope_spread.Draw();
   gSystem->ProcessEvents();
   SaveCanvas(tempcanangle, "alignment_tracks_yslope.png");
   //SaveCanvas(tempcanangle, "alignment_tracks_yslope.C");
   delete tempcanangle;
}

void TDetectorAlignment::AlignDetectorXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose, bool justcheck) {
        
   std::ostringstream detector_name;
   detector_name << "D" << subject_detector;
   
   string titleresx = "residuals" + detector_name.str() + "X";
   string titleresy = "residuals" + detector_name.str() + "Y";
   string titleresxy = "residuals" + detector_name.str() + "XY";
   string titleresxvsy = "residuals" + detector_name.str() + "XvsY";
   string titleresyvsx = "residuals" + detector_name.str() + "YvsX";
   
   TH1F residualsX(titleresx.c_str(),titleresx.c_str(),10000,-100,100);
   TH1F residualsY(titleresy.c_str(),titleresy.c_str(),10000,-100,100);
   TH2F residualsXY(titleresxy.c_str(),titleresxy.c_str(),10000,-100,100,10000,-100,100);
   TH2F residualsXvsY = TH2F(titleresxvsy.c_str(),titleresxvsy.c_str(),256,-0.5,255.5,10000,-100,100);
   TH2F residualsYvsX = TH2F(titleresyvsx.c_str(),titleresyvsx.c_str(),256,-0.5,255.5,10000,-100,100);
   
   vector <float> predictedX, predictedY, observedX, observedY;
   
   Double_t x_offset;
   Double_t y_offset;
   Double_t phix_offset;
   Double_t phiy_offset;
   
   Double_t resxtest, resytest, resxmean, resxrms, resymean, resyrms, res_keep_factor = 2;
   Double_t predx, predy, obsvx, obsvy;
   
   residualsX.Reset();
   residualsY.Reset();
   residualsXY.Reset();
   residualsXvsY.Reset();
   residualsYvsX.Reset();
   
   
   //Align X
   
   //first estimate residuals widths
   resxmean = 0;
   resxrms = 0;
   resymean = 0;
   resyrms = 0;
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
   {
      LoadData(track_storage[t]);
      if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
         PositionPredictor(subject_detector);
      else 
         PositionPredictor(subject_detector, ref_detector1, ref_detector2);
      if(subject_detector==4) track_storage[t].SetDetectorHitPosition(9,GetPredictedY());
      resxtest = track_holder.GetD(subject_detector).GetX()-GetPredictedX();
      resytest = track_holder.GetD(subject_detector).GetY()-GetPredictedY();
      resxmean += resxtest;
      resymean += resytest;
      resxrms += resxtest * resxtest;
      resyrms += resytest * resytest;
   }
   resxmean = resxmean / Double_t(track_storage.size());
   resymean = resymean / Double_t(track_storage.size());
   resxrms = TMath::Sqrt(resxrms / Double_t(track_storage.size()) - resxmean*resxmean);
   resyrms = TMath::Sqrt(resyrms / Double_t(track_storage.size()) - resymean*resymean);
   
   //now select tracks with reasonably small residuals
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
   {
      LoadData(track_storage[t]);
      if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
         PositionPredictor(subject_detector);
      else 
         PositionPredictor(subject_detector, ref_detector1, ref_detector2);
      
      //try guessing the yposition of the diamond tracks
      if(subject_detector==4) track_storage[t].SetDetectorHitPosition(9,GetPredictedY());
      
      predx=GetPredictedX();
      predy=GetPredictedY();
      obsvx=track_holder.GetD(subject_detector).GetX();
      obsvy=track_holder.GetD(subject_detector).GetY();
      
      resxtest=TMath::Abs(obsvx-predx-resxmean)/resxrms/res_keep_factor;
      resytest=TMath::Abs(obsvy-predy-resymean)/resyrms/res_keep_factor;
      
      if(subject_detector<4) {
         if(resxtest<1 && resytest<1) {
            predictedX.push_back(predx);
            predictedY.push_back(predy);
            observedX.push_back(obsvx);
            observedY.push_back(obsvy);
         }
      }
      
      else if(subject_detector==4) {
         if(resxtest<1) {
            predictedX.push_back(predx);
            predictedY.push_back(predy);
            observedX.push_back(obsvx);
            observedY.push_back(predy);
         }
      }
      
      if(0&&t%track_storage.size()/5==0) {
         cout<<"----\nobsvx = "<<obsvx<<"\tpredx = "<<predx<<"\tobsvx-predx = "<<obsvx-predx<<endl;
         cout<<"resxtest = "<<resxtest<<"\tresxmean = "<<resxmean<<"\tresxrms = "<<resxrms<<endl;
         
         if(1&&subject_detector==4) {
            cout<<"resxtest = "<<resxtest<<"\tresxtest<1 = "<<int(resxtest<1)<<endl;
            cout<<"resytest = "<<resytest<<"\tresytest<1 = "<<int(resytest<1)<<"\t(resxtest<1 && resytest<1) = "<<int(resxtest<1 && resytest<1)<<endl;
         }
      }
   }
   cout<<"D"<<subject_detector<<"X: "<<observedX.size()<<" / "<<track_storage.size()<<" = "<<float(observedX.size())/track_storage.size()<<" of tracks survived a "<<res_keep_factor<<" sigma residual cut"<<endl;
   
   // Calculate offsets
   double sumr=0, sumv=0, sumv2=0, sumvr=0;
   for(uint i=0; i<observedX.size(); i++) {
      sumr += observedX[i] - predictedX[i];
      sumv += predictedY[i];
      sumv2 += predictedY[i] * predictedY[i];
      sumvr += predictedY[i] * (observedX[i] - predictedX[i]);
   }
	//TODO: take arctan of that?!
   
   //update offsets and resolutions
   if(!justcheck) {
      x_offset = (sumr * sumv2 - sumvr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
      phix_offset = -(observedX.size() * sumvr - sumr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
      //update x-offsets
      det_x_offset[subject_detector] += x_offset;
      //if(subject_detector)
         det_phix_offset[subject_detector] += phix_offset;
      //if(subject_detector==4) det_phix_offset[subject_detector] = 0.04;
      //update offset history
      det_x_offset_history[subject_detector].push_back(det_x_offset[subject_detector]);
   }
   
   predictedX.clear();
   predictedY.clear();
   observedX.clear();
   observedY.clear();
   residualsX.Reset();
   residualsY.Reset();
   
   
   //Align Y (for silicon only)

   if(subject_detector<4) {
      
      //first estimate residuals widths
      resxmean = 0;
      resxrms = 0;
      resymean = 0;
      resyrms = 0;
      for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
      {
         LoadData(track_storage[t]);
         if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
            PositionPredictor(subject_detector);
         else 
            PositionPredictor(subject_detector, ref_detector1, ref_detector2);
         if(subject_detector==4) track_storage[t].SetDetectorHitPosition(9,GetPredictedY());
         resxtest = track_holder.GetD(subject_detector).GetX()-GetPredictedX();
         resytest = track_holder.GetD(subject_detector).GetY()-GetPredictedY();
         resxmean += resxtest;
         resymean += resytest;
         resxrms += resxtest * resxtest;
         resyrms += resytest * resytest;
      }
      resxmean = resxmean / Double_t(track_storage.size());
      resymean = resymean / Double_t(track_storage.size());
      resxrms = TMath::Sqrt(resxrms / Double_t(track_storage.size()) - resxmean*resxmean);
      resyrms = TMath::Sqrt(resyrms / Double_t(track_storage.size()) - resymean*resymean);
      
      //now select tracks with reasonably small residuals
      for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
      {
         LoadData(track_storage[t]);
         if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
            PositionPredictor(subject_detector);
         else 
            PositionPredictor(subject_detector, ref_detector1, ref_detector2);
         
         //try guessing the yposition of the diamond tracks
         if(subject_detector==4) track_storage[t].SetDetectorHitPosition(9,GetPredictedY());
         
         predx=GetPredictedX();
         predy=GetPredictedY();
         obsvx=track_holder.GetD(subject_detector).GetX();
         obsvy=track_holder.GetD(subject_detector).GetY();
         
         resxtest=TMath::Abs(obsvx-predx-resxmean)/resxrms/res_keep_factor;
         resytest=TMath::Abs(obsvy-predy-resymean)/resyrms/res_keep_factor;
         
         if(0||(resxtest<1 && resytest<1)) {
            predictedX.push_back(predx);
            predictedY.push_back(predy);
            observedX.push_back(obsvx);
            observedY.push_back(obsvy);
         }
         
         if(0&&t%track_storage.size()/5==0) {
            cout<<"----\nobsvx = "<<obsvx<<"\tpredx = "<<predx<<"\tobsvx-predx = "<<obsvx-predx<<endl;
            cout<<"resxtest = "<<resxtest<<"\tresxmean = "<<resxmean<<"\tresxrms = "<<resxrms<<endl;
            cout<<"----\nobsvy = "<<obsvy<<"\tpredy = "<<predy<<"\tobsvy-predy = "<<obsvy-predy<<endl;
            cout<<"resytest = "<<resytest<<"\tresymean = "<<resymean<<"\tresyrms = "<<resyrms<<endl;
         }
      }
      cout<<"D"<<subject_detector<<"Y: "<<observedX.size()<<" / "<<track_storage.size()<<" = "<<float(observedX.size())/track_storage.size()<<" of tracks survived a "<<res_keep_factor<<" sigma residual cut"<<endl;
      
      sumr=0, sumv=0, sumv2=0, sumvr=0;
      for(uint i=0; i<observedX.size(); i++) {
         sumr += observedY[i] - predictedY[i];
         sumv += predictedX[i];
         sumv2 += predictedX[i] * predictedX[i];
         sumvr += predictedX[i] * (observedY[i] - predictedY[i]);
      }
      
      //update offsets and resolutions
      if(!justcheck) {
         y_offset = (sumr * sumv2 - sumvr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
         phiy_offset = (observedX.size() * sumvr - sumr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
         //update y-offsets
         det_y_offset[subject_detector] += y_offset;
         det_phiy_offset[subject_detector] += phiy_offset;
         //update offset history
         det_y_offset_history[subject_detector].push_back(det_y_offset[subject_detector]);
      }
   }
   
   residualsX.Reset();
   residualsY.Reset();
      
      
   //generate residuals to get the new resolutions
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++) {
      LoadData(track_storage[t]);
      if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
         PositionPredictor(subject_detector);
      else 
         PositionPredictor(subject_detector, ref_detector1, ref_detector2);
      residualsX.Fill(track_holder.GetD(subject_detector).GetX()-GetPredictedX());
      residualsY.Fill(track_holder.GetD(subject_detector).GetY()-GetPredictedY());
   }
   
   //zoom residuals to get the right rms (wtf root)
   residualsX.GetXaxis()->SetRangeUser(residualsX.GetMean()-3*residualsX.GetRMS(),residualsX.GetMean()+3*residualsX.GetRMS());
   residualsX.GetXaxis()->SetRangeUser(residualsX.GetMean()-3*residualsX.GetRMS(),residualsX.GetMean()+3*residualsX.GetRMS());
   residualsY.GetXaxis()->SetRangeUser(residualsY.GetMean()-3*residualsY.GetRMS(),residualsY.GetMean()+3*residualsY.GetRMS());
   residualsY.GetXaxis()->SetRangeUser(residualsY.GetMean()-3*residualsY.GetRMS(),residualsY.GetMean()+3*residualsY.GetRMS());
   
   //fit gaussians to center of residuals
   TF1 histofitx("histofitx","gaus",residualsX.GetMean()-2*residualsX.GetRMS(),residualsX.GetMean()+2*residualsX.GetRMS());
   TF1 histofity("histofity","gaus",residualsY.GetMean()-2*residualsY.GetRMS(),residualsY.GetMean()+2*residualsY.GetRMS());
   histofitx.SetLineColor(kBlue);
   histofity.SetLineColor(kBlue);
   residualsX.Fit(&histofitx,"rq"); // fit option "r" restricts the range of the fit; "q" quiets output
   residualsY.Fit(&histofity,"rq");
   
   /*
   //just for investigation purposes
   //it looks like the rms changes depending on how far zoomed in you look
   TCanvas canv("asdfa");
   residualsX.GetXaxis()->SetRangeUser(residualsX.GetMean()-3*residualsX.GetRMS(),residualsX.GetMean()+3*residualsX.GetRMS());
   residualsX.GetXaxis()->SetRangeUser(residualsX.GetMean()-3*residualsX.GetRMS(),residualsX.GetMean()+3*residualsX.GetRMS());
   residualsX.Draw();
   canv.Print("/home/jduris/diamond/runs/newcode_sandbox/residualsX.png");
   */
   
   //update resolution
   det_x_resolution[subject_detector] = histofitx.GetParameter(2);
   det_y_resolution[subject_detector] = histofity.GetParameter(2);
   //det_x_resolution[subject_detector] = residualsX.GetRMS();
   //det_y_resolution[subject_detector] = residualsY.GetRMS();
   
   //cout << "For " << detector_name << "X mean is " << x_offset << " and RMS is " << residualsX->GetRMS() << endl;
   //cout << "For " << detector_name << "Y mean is " << y_offset << " and RMS is " << residualsY->GetRMS() << endl;
   
   if(verbose) {
      cout << "For " << detector_name.str() << "X change in offset is " << x_offset
            << ", new offset is " << det_x_offset[subject_detector]
            << ", phiX change in offset is " << phix_offset
            << ", new offset is " << det_phix_offset[subject_detector]
            << ", resolution is " << det_x_resolution[subject_detector] << endl;
      if(1||subject_detector<4) {
         cout << "For " << detector_name.str() << "Y change in offset is " << y_offset
               << ", new offset is " << det_y_offset[subject_detector]
               << ", phiY change in offset is " << phiy_offset
               << ", new offset is " << det_phiy_offset[subject_detector]
               << ", resolution is " << det_y_resolution[subject_detector] << endl;
      }
   }
}

void TDetectorAlignment::AlignDetectorXY(int subject_detector, bool verbose, bool justcheck) {
   AlignDetectorXY(subject_detector, -1, -1, verbose, justcheck);
}

void TDetectorAlignment::CheckDetectorAlignmentXY(int subject_detector, int ref_detector1, int ref_detector2, bool verbose) {
   AlignDetectorXY(subject_detector, ref_detector1, ref_detector2, verbose, true);
}

void TDetectorAlignment::CheckDetectorAlignmentXY(int subject_detector, bool verbose) {
   AlignDetectorXY(subject_detector, verbose, true);
}

void TDetectorAlignment::CheckDetectorAlignmentXYPlots(int subject_detector, int ref_detector1, int ref_detector2, string& histo_title) {
   TF1 residualsXfit, residualsYfit; //to be fixed
         
   std::ostringstream detector_name;
   detector_name << "D" << subject_detector;
   
   string titleresx = histo_title + detector_name.str() + "X";
   string titleresy = histo_title + detector_name.str() + "Y";
   string titleresxy = histo_title + detector_name.str() + "XY";
   string titleresxvsy = histo_title + detector_name.str() + "XvsY";
   string titleresyvsx = histo_title + detector_name.str() + "YvsX";
   
   TH1F plotresidualsX(titleresx.c_str(),titleresx.c_str(),10000,-100,100);
   TH1F plotresidualsY(titleresy.c_str(),titleresy.c_str(),10000,-100,100);
   TH2F plotresidualsXY(titleresxy.c_str(),titleresxy.c_str(),10000,-100,100,10000,-100,100);
   TH2F plotresidualsXvsY = TH2F(titleresxvsy.c_str(),titleresxvsy.c_str(),256,-0.5,255.5,10000,-100,100);
   TH2F plotresidualsYvsX = TH2F(titleresyvsx.c_str(),titleresyvsx.c_str(),256,-0.5,255.5,10000,-100,100);

   Double_t resxmean = 0, resxrms = 0, resymean = 0, resyrms = 0, resxtest, resytest;
   
   for(Int_t t=0; t<(Int_t)track_storage.size(); t++)
   {
      LoadData(track_storage[t]);
      if(ref_detector1<0||ref_detector1>5||ref_detector2<0||ref_detector2>5) 
         PositionPredictor(subject_detector);
      else 
         PositionPredictor(subject_detector, ref_detector1, ref_detector2);
      plotresidualsX.Fill(track_holder.GetD(subject_detector).GetX()-GetPredictedX());
      plotresidualsY.Fill(track_holder.GetD(subject_detector).GetY()-GetPredictedY());
      plotresidualsXY.Fill(track_holder.GetD(subject_detector).GetX()-GetPredictedX(),track_holder.GetD(subject_detector).GetY()-GetPredictedY());
      //plotresidualsXvsY.Fill(track_holder.GetD(subject_detector).GetY(),track_holder.GetD(subject_detector).GetX()-GetPredictedX());
      //plotresidualsYvsX.Fill(track_holder.GetD(subject_detector).GetX(),track_holder.GetD(subject_detector).GetY()-GetPredictedY());
      plotresidualsXvsY.Fill(GetPredictedY(),track_holder.GetD(subject_detector).GetX()-GetPredictedX());
      plotresidualsYvsX.Fill(GetPredictedX(),track_holder.GetD(subject_detector).GetY()-GetPredictedY());
      
      /*
      resxtest = track_holder.GetD(subject_detector).GetX()-GetPredictedX();
      resytest = track_holder.GetD(subject_detector).GetY()-GetPredictedY();
      resxmean += resxtest;
      resymean += resytest;
      resxrms += resxtest * resxtest;
      resyrms += resytest * resytest;
      */
   }
   /*
   resxmean = resxmean / Double_t(track_storage.size());
   resymean = resymean / Double_t(track_storage.size());
   resxrms = TMath::Sqrt(resxrms / Double_t(track_storage.size()) - resxmean*resxmean);
   resyrms = TMath::Sqrt(resyrms / Double_t(track_storage.size()) - resymean*resymean);
  */
   
   //plot
   Double_t plot_width_factor = 7;
   Double_t plot_fit_factor = 3/4.;
   
   TCanvas *tempcan = new TCanvas("residualstempcanv","residualstempcanv",800,600);
   //plotresidualsX.GetXaxis()->SetRangeUser(resxmean-plot_width_factor*resxrms,resxmean+plot_width_factor*resxrms);
   //TF1 histofitx("histofitx","gaus",resxmean-plot_fit_factor*resxrms,resxmean+plot_fit_factor*resxrms);
   plotresidualsX.GetXaxis()->SetRangeUser(plotresidualsX.GetMean()-plot_width_factor*plotresidualsX.GetRMS(),plotresidualsX.GetMean()+plot_width_factor*plotresidualsX.GetRMS());
   plotresidualsX.GetXaxis()->SetRangeUser(plotresidualsX.GetMean()-plot_width_factor*plotresidualsX.GetRMS(),plotresidualsX.GetMean()+plot_width_factor*plotresidualsX.GetRMS());
   TF1 histofitx("histofitx","gaus",plotresidualsX.GetMean()-plot_fit_factor*plotresidualsX.GetRMS(),plotresidualsX.GetMean()+plot_fit_factor*plotresidualsX.GetRMS());
   histofitx.SetLineColor(kBlue);
   plotresidualsX.Fit(&histofitx,"rq");
   plotresidualsX.Draw();
   //cout<<"plotresidualsX.GetMean()="<<plotresidualsX.GetMean()<<"\tplotresidualsX.GetRMS()="<<plotresidualsX.GetRMS()<<endl;
   cout<<"histofitx.GetParameter(1)="<<histofitx.GetParameter(1)<<"\thistofitx.GetParameter(2)="<<histofitx.GetParameter(2)<<endl;
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      ostringstream titleresxoss;
      titleresxoss << titleresx << ".png";
      SaveCanvas(tempcan, titleresxoss.str().c_str());
      if(ClosePlotsOnSave == 1)
      {
         delete tempcan;
      }
   }
   
   tempcan = new TCanvas("residualstempcanv","residualstempcanv",800,600);
   //plotresidualsY.GetXaxis()->SetRangeUser(resymean-plot_width_factor*resyrms,resymean+plot_width_factor*resyrms);
   //TF1 histofity("histofity","gaus",resymean-plot_fit_factor*resyrms,resymean+plot_fit_factor*resyrms);
   plotresidualsY.GetXaxis()->SetRangeUser(plotresidualsY.GetMean()-plot_width_factor*plotresidualsY.GetRMS(),plotresidualsY.GetMean()+plot_width_factor*plotresidualsY.GetRMS());
   plotresidualsY.GetXaxis()->SetRangeUser(plotresidualsY.GetMean()-plot_width_factor*plotresidualsY.GetRMS(),plotresidualsY.GetMean()+plot_width_factor*plotresidualsY.GetRMS());
   TF1 histofity("histofity","gaus",plotresidualsY.GetMean()-plot_fit_factor*plotresidualsY.GetRMS(),plotresidualsY.GetMean()+plot_fit_factor*plotresidualsY.GetRMS());
   histofity.SetLineColor(kBlue);
   plotresidualsY.Fit(&histofity,"rq"); //r restricts fit range while q quiets output
   plotresidualsY.Draw();
   //cout<<"plotresidualsY.GetMean()="<<plotresidualsY.GetMean()<<"\tplotresidualsY.GetRMS()="<<plotresidualsY.GetRMS()<<endl;
   cout<<"histofity.GetParameter(1)="<<histofity.GetParameter(1)<<"\thistofity.GetParameter(2)="<<histofity.GetParameter(2)<<endl;
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      ostringstream titleresyoss;
      titleresyoss << titleresy << ".png";
      SaveCanvas(tempcan, titleresyoss.str().c_str());
      if(ClosePlotsOnSave == 1)
      {
         delete tempcan;
      }
   }
   
   gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
   tempcan = new TCanvas("residualstempcanv","residualstempcanv",800,600);
   Double_t xtoyratio = tempcan->GetWindowWidth()/Double_t(tempcan->GetWindowHeight());
   //plotresidualsXY.GetXaxis()->SetRangeUser(resxmean-plot_width_factor*resxrms,resxmean+plot_width_factor*resxrms);
   //plotresidualsXY.GetYaxis()->SetRangeUser(resymean-plot_width_factor*resyrms,resymean+plot_width_factor*resyrms);
   plotresidualsXY.GetXaxis()->SetRangeUser(plotresidualsXY.GetMean(1)-2*plot_width_factor*plotresidualsXY.GetRMS(1),plotresidualsXY.GetMean(1)+2*plot_width_factor*plotresidualsXY.GetRMS(1));
   plotresidualsXY.GetYaxis()->SetRangeUser(plotresidualsXY.GetMean(2)-2*plot_width_factor*plotresidualsXY.GetRMS(2)*xtoyratio,plotresidualsXY.GetMean(2)+2*plot_width_factor*plotresidualsXY.GetRMS(2)*xtoyratio);
   plotresidualsXY.Draw("colz");
   //cout<<"plotresidualsXY.GetMean()="<<plotresidualsXY.GetMean()<<"\tplotresidualsXY.GetRMS()="<<plotresidualsXY.GetRMS()<<endl;
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      ostringstream titleresxyoss;
      titleresxyoss << titleresxy << ".png";
      SaveCanvas(tempcan, titleresxyoss.str().c_str());
      if(ClosePlotsOnSave == 1)
      {
         delete tempcan;
      }
   }
   
   tempcan = new TCanvas("residualstempcanv","residualstempcanv",800,600);
   //plotresidualsXvsY.GetXaxis()->SetRangeUser(plotresidualsXvsY.GetMean(1)-plot_width_factor*plotresidualsXY.GetRMS(1),plotresidualsXY.GetMean(1)+plot_width_factor*plotresidualsXY.GetRMS(1));
   plotresidualsXvsY.GetYaxis()->SetRangeUser(plotresidualsXvsY.GetMean(2)-plot_width_factor*plotresidualsXvsY.GetRMS(2),plotresidualsXvsY.GetMean(2)+plot_width_factor*plotresidualsXvsY.GetRMS(2));
   //plotresidualsXvsY.GetYaxis()->SetRangeUser(resxmean-plot_width_factor*resxrms,resxmean+plot_width_factor*resxrms);
   plotresidualsXvsY.Draw("colz");
   //cout<<"plotresidualsXvsY.GetMean()="<<plotresidualsXvsY.GetMean()<<"\tplotresidualsXvsY.GetRMS()="<<plotresidualsXvsY.GetRMS()<<endl;
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      ostringstream titleresxvsyoss;
      titleresxvsyoss << titleresxvsy << ".png";
      SaveCanvas(tempcan, titleresxvsyoss.str().c_str());
      if(ClosePlotsOnSave == 1)
      {
         delete tempcan;
      }
   }
   
   tempcan = new TCanvas("residualstempcanv","residualstempcanv",800,600);
   //plotresidualsXY.GetXaxis()->SetRangeUser(plotresidualsXY.GetMean(1)-plot_width_factor*plotresidualsXY.GetRMS(1),plotresidualsXY.GetMean(1)+plot_width_factor*plotresidualsXY.GetRMS(1));
   plotresidualsYvsX.GetYaxis()->SetRangeUser(plotresidualsYvsX.GetMean(2)-plot_width_factor*plotresidualsYvsX.GetRMS(2),plotresidualsYvsX.GetMean(2)+plot_width_factor*plotresidualsYvsX.GetRMS(2));
   //plotresidualsYvsX.GetYaxis()->SetRangeUser(resymean-plot_width_factor*resyrms,resymean+plot_width_factor*resyrms);
   plotresidualsYvsX.Draw("colz");
   //cout<<"plotresidualsYvsX.GetMean()="<<plotresidualsYvsX.GetMean()<<"\tplotresidualsYvsX.GetRMS()="<<plotresidualsYvsX.GetRMS()<<endl;
   if(SaveAllFilesSwitch == 1)
   {
      gSystem->ProcessEvents();
      ostringstream titleresyvsxoss;
      titleresyvsxoss << titleresyvsx << ".png";
      SaveCanvas(tempcan, titleresyvsxoss.str().c_str());
      if(ClosePlotsOnSave == 1)
      {
         delete tempcan;
      }
   }
}

void TDetectorAlignment::CheckDetectorAlignmentXYPlots(int subject_detector, string &histo_title) {
   CheckDetectorAlignmentXYPlots(subject_detector, -1, -1, histo_title);
}

void TDetectorAlignment::AlignDetectorZ(Int_t subject_detector, bool graphs, bool verbose) {
   AlignDetectorZ(subject_detector, -1, -1, graphs, verbose);
}

void TDetectorAlignment::AlignDetectorZ(Int_t subject_detector, Int_t ref_detector1, Int_t ref_detector2, bool graphs, bool verbose) {
   Double_t report_zoffset_x, report_zoffset_y; //to be deleted
   
   std::ostringstream detector_name;
   detector_name << "D" << subject_detector;
   
   string tempname = detector_name.str() + "ZAlignmentResiduals";
   
   TH1F residualsX, residualsY;
   //TH1F *middleresX = new TH1F("middleresX", "middleresX", 10000,-100,100); 
   //TH1F *middleresY = new TH1F("middleresY", "middleresY", 10000,-100,100); 
   
   TF1 fitdummy;
   vector<TDiamondTrack> temp_tracks;
   vector<Double_t> zvals, resxmean, resymean, resxrms, resyrms, dresxrms, dresyrms, dresxyrms, dresxyrmsscaled;
   Double_t bestz = 0;
   //Double_t old_det_z_offset = det_z_offset[subject_detector];
   
   Int_t track_number = track_storage.size();
   //Float_t initialz = 0;
   Float_t initialz = track_storage[0].GetD(subject_detector).GetZ();
   if(subject_detector==1) initialz = (track_storage[0].GetD(2).GetZ() + track_storage[0].GetD(0).GetZ())/2.;
   if(subject_detector==2) initialz = (track_storage[0].GetD(3).GetZ() + track_storage[0].GetD(1).GetZ())/2.;
   if(subject_detector==3) initialz = (track_storage[0].GetD(2).GetZ() + track_storage[0].GetD(1).GetZ())/2.;
   cout<<"initial z for detector "<<subject_detector<<" is "<<initialz<<endl;
   Float_t deltaz = 0.01; // in cm 
   Float_t zrange = TMath::Max(1.0, 3.*TMath::Abs(initialz - track_storage[0].GetD(subject_detector).GetZ())); // in cm (we'll crawl forward and backward by this amount)
   
   if(verbose) cout<<"z"<<"\t"<<"resXMean"<<"\t"<<"resYMean"<<"\t"<<"resXRMS"<<"\t"<<"resYRMS"<<endl;//<<"\t"<<"resXEntries"<<"\t"<<"resYEntries"<<endl;
   
   //loop over possible z
   for(Double_t z=initialz-zrange; z<initialz+zrange; z+=deltaz)
   //for(Double_t z=initialz; z==initialz; z+=deltaz)
   {
      //First modify z of tracks
      temp_tracks.clear();
      for(Int_t t=0; t<track_number; t++)
      {
         TDiamondTrack track = track_storage[t];
         TDetectorPlane D0_buf = track.GetD0();
         TDetectorPlane D1_buf = track.GetD1();
         TDetectorPlane D2_buf = track.GetD2();
         TDetectorPlane D3_buf = track.GetD3();
      
         if(subject_detector == 0) D0_buf.SetZ(z);
         if(subject_detector == 1) D1_buf.SetZ(z);
         if(subject_detector == 2) D2_buf.SetZ(z);
         if(subject_detector == 3) D3_buf.SetZ(z);
     
         TDiamondTrack new_track = TDiamondTrack(D0_buf, D1_buf, D2_buf, D3_buf);
         temp_tracks.push_back(new_track);
      }
      det_z_offset[subject_detector] = z;
      //Next calculate the residuals for this z
      residualsX.Reset();
      residualsY.Reset();
      CheckDetectorAlignmentXY(subject_detector, ref_detector1, ref_detector2, false);

      /*
      if((z-initialz)*(z-initialz)<deltaz/10.){
      middleresX->Add(&residualsX);
      middleresY->Add(&residualsY);
   }*/
      
      //Tabulate the residuals' means and RMS's
      zvals.push_back(z);
      resxmean.push_back(residualsX.GetMean());
      resymean.push_back(residualsY.GetMean());
      resxrms.push_back(residualsX.GetRMS());
      resyrms.push_back(residualsY.GetRMS());
      
      if(verbose) cout<<z<<"\t"<<residualsX.GetMean()<<"\t"<<residualsY.GetMean()<<"\t"<<residualsX.GetRMS()<<"\t"<<residualsY.GetRMS()<<endl;
            
      //tempname += "_";
      residualsX.Reset();
      residualsY.Reset();
   }
   
   //figure out derivatives
   dresxrms.push_back(TMath::Abs((-3.*resxrms[0] + 4.*resxrms[1] - resxrms[2])/2./deltaz)); //3point end derivative
   for(Int_t i=1; i<(int)resxrms.size()-1; i++) 
      dresxrms.push_back(TMath::Abs((resxrms[i+1] - resxrms[i-1])/2./deltaz));
   dresxrms.push_back(TMath::Abs((3.*resxrms[resxrms.size()-1] - 4.*resxrms[resxrms.size()-2] + resxrms[resxrms.size()-3])/2./deltaz));
   
   dresyrms.push_back(TMath::Abs((-3.*resyrms[0] + 4.*resyrms[1] - resyrms[2])/2./deltaz));
   for(Int_t i=1; i<(int)resyrms.size()-1; i++) 
      dresyrms.push_back(TMath::Abs((resyrms[i+1] - resyrms[i-1])/2./deltaz));
   dresyrms.push_back(TMath::Abs((3.*resyrms[resyrms.size()-1] - 4.*resyrms[resyrms.size()-2] + resyrms[resyrms.size()-3])/2./deltaz));
   
   for(Int_t i=0; i<(int)resxrms.size(); i++) dresxyrms.push_back(TMath::Sqrt(dresxrms[i]*dresxrms[i]));// + dresyrms[i]*dresyrms[i]));
   
   /*
   //print out diagnostics
   cout<<"phi"<<"\t"<<"resXRMS"<<"\t"<<"resYRMS"<<"\t"<<"dresXRMS"<<"\t"<<"dresYRMS"<<"\t"<<"dresXYRMS"<<endl;
   for(Int_t i=0; i<resxrms.size(); i++){
   cout<<phivals[i]<<"\t"<<"\t"<<resxrms[i]<<"\t"<<resyrms[i]<<"\t"<<dresxrms[i]<<"\t"<<dresyrms[i]<<"\t"<<dresxyrms[i]<<endl;
}
   */
   /*
   //figure out best z and rms there
   for(Int_t i=0; i<(int)dresxrms.size(); i++){
   if(TMath::Abs(dresxrms[i])<TMath::Abs(dresxrms[(int)bestz])) bestz=i;
}
   Float_t report_zresolution_x = resxrms[(int)bestz];
   bestz=zvals[(int)bestz];
   cout<<"Minimum for "<<detector_name.str()<<"X at: "<<bestz<<endl;
   report_zoffset_x = bestz;
   bestz=0;
   for(Int_t i=0; i<(int)dresyrms.size(); i++){
   if(TMath::Abs(dresyrms[i])<TMath::Abs(dresyrms[(int)bestz])) bestz=i;
}
   Float_t report_zresolution_y = resyrms[(int)bestz];
   bestz=zvals[(int)bestz];
   cout<<"Minimum for "<<detector_name.str()<<"Y at: "<<bestz<<endl;
   report_zoffset_y = bestz;
   */
   
   //figure out best z and rms there
   for(Int_t i=0; i<(int)resxrms.size(); i++){
      if(TMath::Abs(resxrms[i])<TMath::Abs(resxrms[(int)bestz])) bestz=i;
   }
   Float_t report_zresolution_x = resxrms[(int)bestz];
   bestz=zvals[(int)bestz];
   cout<<"Minimum for "<<detector_name.str()<<"X at: "<<bestz<<endl;
   report_zoffset_x = bestz;
   bestz=0;
   for(Int_t i=0; i<(int)resyrms.size(); i++){
      if(TMath::Abs(resyrms[i])<TMath::Abs(resyrms[(int)bestz])) bestz=i;
   }
   Float_t report_zresolution_y = resyrms[(int)bestz];
   bestz=zvals[(int)bestz];
   cout<<"Minimum for "<<detector_name.str()<<"Y at: "<<bestz<<endl;
   report_zoffset_y = bestz;
   
   //update current offsets
   det_z_offset[subject_detector] += report_zoffset_x;
   
   //update offset history
   det_z_offset_history[subject_detector].push_back(det_z_offset[subject_detector]);
   
   //update resolution
   //det_z_resolution[subject_detector] = report_zresolution_x;
   
   /*
   if(subject_detector==1) {
   report_offset_z = (track_storage[0].GetD(2).GetZ() + track_storage[0].GetD(0).GetZ())/2 - bestz;
   cout<<"Best z for "<<detector_name.str()<<" is: "<<report_offset_z + track_storage[0].GetD(1).GetZ()<<endl;
}
   if(subject_detector==2) {
   report_offset_z = (track_storage[0].GetD(3).GetZ() + track_storage[0].GetD(1).GetZ())/2 - bestz;
   cout<<"Best z for "<<detector_name.str()<<" is: "<<report_offset_z + track_storage[0].GetD(2).GetZ()<<endl;
}
   if(subject_detector==3) {
   report_offset_z = (track_storage[0].GetD(2).GetZ() + track_storage[0].GetD(1).GetZ())/2 - bestz;
   cout<<"Best z for "<<detector_name.str()<<" is: "<<report_offset_z + track_storage[0].GetD(3).GetZ();
}
   cout<<"Best change in z for "<<detector_name.str()<<" is: "<<report_offset_z<<endl;
   */
   
   //Plots
   string titleresx = detector_name.str() + "X residuals rms vs z";
   string titleresy = detector_name.str() + "Y residuals rms vs z";
   string titleresxy = detector_name.str() + "X (red*) and Y (blue*) residuals rms and derivatives (dots) vs z";
   
   string centralxresname = "alignment_z_" + detector_name.str() + "centralXresiduals";
   string centralyresname = "alignment_z_" + detector_name.str() + "centralYresiduals";
   string filenameresx = "alignment_z_" + detector_name.str() + "XresidualRMSvsZ";
   string filenameresy = "alignment_z_" + detector_name.str() + "YresidualRMSvsZ";
   string filenameresxy = "alignment_z_" + detector_name.str() + "XandYresidualsRMSvsZ";
   
   if(graphs) {
   //graphs
      TF1 fitx("fitx","[2]*(x-[0])*(x-[0])+[1]",report_zoffset_x-0.5,report_zoffset_x+0.5);
      fitx.SetParameter(0,report_zoffset_x);
      fitx.SetParameter(1,report_zresolution_x);
      fitx.SetParameter(2,0.01);
      fitx.SetParLimits(0,initialz-zrange,initialz+zrange);
      fitx.SetParLimits(1,0,10*report_zresolution_x);
      fitx.SetParLimits(2,0,1);
      TF1 fity("fity","[2]*(x-[0])*(x-[0])+[1]",report_zoffset_y-0.5,report_zoffset_y+0.5);
      fity.SetParameter(0,report_zoffset_y);
      fity.SetParameter(1,report_zresolution_y);
      fity.SetParameter(2,0.01);
      fitx.SetParLimits(0,initialz-zrange,initialz+zrange);
      fity.SetParLimits(1,0,10*report_zresolution_y);
      fity.SetParLimits(2,0,1);
      
      TGraph *residualRMSXgraph = new TGraph(zvals.size(), (const Double_t*) &zvals[0], (const Double_t*) &resxrms[0]);
      residualRMSXgraph->SetTitle("x residuals rms vs z");
      TCanvas *tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      residualRMSXgraph->Fit("fitx","r"); // fit option "r" restricts the fit range
      residualRMSXgraph->Draw("A*");
      gSystem->ProcessEvents();
      SaveCanvas(tempcan, "alignment_z_residualXRMSvsZ.png");
      delete tempcan;
      
      if(verbose)cout<<detector_name.str()<<"Z-offset is "<<fitx.GetParameter(0)<<"cm with uncertainty "<<fitx.GetParError(0)<<"cm.  ";
      if(verbose)cout<<detector_name.str()<<"X-resolution is "<<fitx.GetParameter(1)<<" units of pitch with uncertainty "<<fitx.GetParError(1)<<" units of pitch."<<endl;
   
      TGraph *residualRMSYgraph = new TGraph(zvals.size(), (const Double_t*) &zvals[0], (const Double_t*) &resyrms[0]);
      residualRMSYgraph->SetTitle("y residuals rms vs z");
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      residualRMSYgraph->Fit("fity","r"); // fit option "r" restricts the fit range
      residualRMSYgraph->Draw("A*");
      gSystem->ProcessEvents();
      SaveCanvas(tempcan, "alignment_z_residualYRMSvsZ.png");
      delete tempcan;
   
      if(verbose)cout<<detector_name.str()<<"Z-offset is "<<fity.GetParameter(0)<<"cm with uncertainty "<<fity.GetParError(0)<<"cm.  ";
      if(verbose)cout<<detector_name.str()<<"Y-resolution is "<<fity.GetParameter(1)<<" units of pitch with uncertainty "<<fity.GetParError(1)<<" units of pitch."<<endl;
      
      TGraph *residualmeanXgraph = new TGraph(zvals.size(), (const Double_t*) &zvals[0], (const Double_t*) &resxmean[0]);
      residualmeanXgraph->SetTitle("x residuals mean vs z");
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      residualmeanXgraph->Draw("A*");
      gSystem->ProcessEvents();
      SaveCanvas(tempcan, "alignment_z_residualXmeanvsZ.png");
      delete tempcan;
   
      TGraph *residualmeanYgraph = new TGraph(zvals.size(), (const Double_t*) &zvals[0], (const Double_t*) &resymean[0]);
      residualmeanYgraph->SetTitle("y residuals mean vs z");
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      residualmeanYgraph->Draw("A*");
      gSystem->ProcessEvents();
      SaveCanvas(tempcan, "alignment_z_residualYmeanvsZ.png");
      delete tempcan;
      /*
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      middleresX->GetXaxis()->SetRangeUser(middleresX->GetMean()-3*middleresX->GetRMS(),middleresX->GetMean()+3*middleresX->GetRMS());
      middleresX->Draw();
      gSystem->ProcessEvents();
      SaveCanvasPNG(tempcan, "", (char*) centralxresname.c_str());
      delete tempcan;
   
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 400);
      tempcan->cd(1);
      middleresY->GetXaxis()->SetRangeUser(middleresY->GetMean()-3*middleresY->GetRMS(),middleresY->GetMean()+3*middleresY->GetRMS());
      middleresY->Draw();
      gSystem->ProcessEvents();
      SaveCanvasPNG(tempcan, "", (char*) centralyresname.c_str());
      delete tempcan;
      */
   //plot residuals and their derivatives on one graph for comparison
      for(Int_t i=0; i<(int)resxrms.size(); i++) dresxyrmsscaled.push_back(dresxyrms[i]*TMath::Sqrt(resxrms[1]*resxrms[1]+resyrms[1]*resyrms[1])/dresxyrms[1]); //scale graph
      TGraph *dresidualRMSXYgraph = new TGraph(zvals.size(), (const Double_t*) &zvals[0], (const Double_t*) &dresxyrmsscaled[0]);
      residualRMSXgraph->SetTitle(titleresxy.c_str());
      tempcan = new TCanvas("tempcanvas", "tempcanvas", 200, 10, 600, 600);
      tempcan->cd(1);
   //scale graphs
      Double_t graphmin = resxrms[0];
      Double_t graphmax = resxrms[0];
      for(Int_t i=0; i<(int)resxrms.size(); i++){
         if(resxrms[i]>graphmax) graphmax=resxrms[i];
         if(resxrms[i]<graphmin) graphmin=resxrms[i];
         if(resyrms[i]>graphmax) graphmax=resyrms[i];
         if(resyrms[i]<graphmin) graphmin=resyrms[i];
         if(dresxyrmsscaled[i]>graphmax) graphmax=dresxyrmsscaled[i];
         if(dresxyrmsscaled[i]<graphmin) graphmin=dresxyrmsscaled[i];
      }
      graphmin -= 0.1*(graphmax-graphmin);
      graphmax += 0.1*(graphmax-graphmin);
   //cout<<"graphmin="<<graphmin<<"\tgraphmax="<<graphmax<<endl;
      residualRMSXgraph->SetMinimum(graphmin);
      residualRMSYgraph->SetMinimum(graphmin);
      dresidualRMSXYgraph->SetMinimum(graphmin);
      residualRMSXgraph->SetMaximum(graphmax);
      residualRMSYgraph->SetMaximum(graphmax);
      dresidualRMSXYgraph->SetMaximum(graphmax);
   //draw graphs
      residualRMSXgraph->SetMarkerColor(kRed);
      residualRMSXgraph->Draw("A*");
      residualRMSYgraph->SetMarkerColor(kBlue);
      residualRMSYgraph->Draw("*same");
      dresidualRMSXYgraph->SetMarkerStyle(20);
      dresidualRMSXYgraph->SetMarkerSize(0.5);
      dresidualRMSXYgraph->Draw("*same");
      gSystem->ProcessEvents();
      ostringstream filenameresxyoss;
      filenameresxyoss << filenameresxy << ".png";
      SaveCanvas(tempcan, filenameresxyoss.str().c_str());
      delete tempcan;
   }
   
   residualsX.Reset();
   residualsY.Reset();
   //middleresX->Clear();
   //middleresY->Clear();
   
   //delete middleresX;
   //delete middleresY;
   
}

void TDetectorAlignment::CutFakeTracks() {
	vector<Float_t> X_strips_X_positions, Y_strips_Y_positions, X_strips_Z_positions, Y_strips_Z_positions;
	Float_t X_Mean, Y_Mean, Z_Mean;
	Float_t res = GetSiResolution();
	//TH1F *middleresX = new TH1F("middleresX", "middleresX", 10000,-100,100);
	TH1F *histo_alignmentfitchi2_XStrips = new TH1F("histo_alignmentfitchi2_XStrips","histo_alignmentfitchi2_XStrips",100,0.,20.);
	TH1F *histo_alignmentfitchi2_YStrips = new TH1F("histo_alignmentfitchi2_YStrips","histo_alignmentfitchi2_YStrips",100,0.,20.);
	for (int i = 0; i < track_storage.size(); i++) {
		LoadData(track_storage[i]);
		X_strips_X_positions.clear();
		Y_strips_Y_positions.clear();
		X_strips_Z_positions.clear();
		Y_strips_Z_positions.clear();
		//		track_storage_aligned.push_back(track_holder);
		for (int det = 0; det < 8; det++) {
			if (det%2 == 0) { // x strips
				X_strips_X_positions.push_back(track_holder.GetD(det).GetX());
				X_strips_Z_positions.push_back(track_holder.GetD(det).GetZ());
			}
			else { // y strips
				Y_strips_Y_positions.push_back(track_holder.GetD(det).GetY());
				Y_strips_Z_positions.push_back(track_holder.GetD(det).GetZ());
			}
		}
		Float_t *X_strips_par[3];
		Float_t *Y_strips_par[3];
		X_strips_par = LinTrackFit(X_strips_Z_positions, X_strips_X_positions, res);
		Y_strips_par = LinTrackFit(Y_strips_Z_positions, Y_strips_Y_positions, res);
		histo_alignmentfitchi2_XStrips->Fill(X_strips_par[2]);
		histo_alignmentfitchi2_YStrips->Fill(Y_strips_par[2]);
	}
	tempcan = new TCanvas("tempcanvas");
	histo_alignmentfitchi2_XStrips->Draw();
	SaveCanvas(tempcan, "alignmentfitchi2_XStrips.png");
	histo_alignmentfitchi2_YStrips->Draw();
	SaveCanvas(tempcan, "alignmentfitchi2_YStrips.png");
	delete tmpcan;
}

Float_t *TDetectorAlignment::LinTrackFit(vector<Float_t> X, vector<Float_t> Y, Float_t res) {
	// -- fits Y = par[0] + par[1] * x
	// -- returns par, with par[2] = chi2
	Float_t X_Mean = 0;
	Float_t Y_Mean = 0;
	Float_t tmp1 = 0;
	Float_t tmp2 = 0;
	Float_t tmp3 = 0;
	Float_t *par[3];
	if (X.size() != Y.size()) return -1;
	for (int i = 0; i < X.size(); i++) {
		X_Mean = X_Mean + X[i];
		Y_Mean = Y_Mean + Y[i];
	}
	X_Mean = X_Mean / X.size();
	Y_Mean = Y_Mean / Y.size();
	for (int i = 0; i < X.size(); i++) {
		tmp1 = tmp1 + ((X[i] - X_Mean) * (Y[i] - Y_Mean));
		tmp2 = tmp2 + ((X[i] - X_Mean) * (X[i] - X_Mean));
	}
	par[1] = tmp1 / tmp2;
	par[0] = Y_Mean - par[1] * X_Mean;
	for (int i = 0; i < X.size(); i++) {
		tmp1 = par[0] + par[1] * X[i];
		tmp3 = tmp3 + (tmp1 - X[i]) * (tmp1 - X[i]);
	}
	par[2] = tmp3 / (res * res);
	return &par[0];
}
