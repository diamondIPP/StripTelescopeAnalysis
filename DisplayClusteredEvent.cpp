//Display pedestal subtracted events from PedestalSubtracted.root (need to run PedestalAnalyze first)
//2010-07-22 First verison based on DisplayRawEvent.cpp
//2010-07-27 Now based on dynamic arrays version of Clustering.class.cpp
//           Now based on ClusteredEvent class and GetClusteredEvent function 
//2010-09-16 Fixed event lookup

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>

#include "Cluster.class.hh"
#include "ClusteredEvent.class.hh"
#include "GetClusteredEvent.cpp"

using namespace std;

// argv[]={DisplayClusteredEvent, int EventNumber}
//int main(int argc, char* argv[]) {
int DisplayClusteredEvent(unsigned int RunNumber, unsigned int EventNumber, string RunDescription = "", bool verbose = 0) {
   
   //if(argc==2 || argc==3) {
   if(1) {
      //int EventNumber = (int)strtod(argv[1],0);
      //bool verbose = 0;
      //if(argc==3) verbose = (int)strtod(argv[3],0);
      
      
      //clusters for event
      ClusteredEvent clustered_event = GetClusteredEvent(RunNumber, EventNumber, RunDescription, verbose);
      clustered_event.Print();
      
      //default paths
      
      TSystem* sys = gSystem;
      string settings_file;
   
      ostringstream settingspath;
      settingspath << sys->pwd() << "/Settings." << RunNumber;
      if(RunDescription=="") settingspath << ".ini";
      else settingspath << "-" << RunDescription << ".ini";
      settings_file = settingspath.str();
   
      ostringstream pedfilepath;
      pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
      if(RunDescription=="") pedfilepath << ".root";
      else pedfilepath << "-" << RunDescription << ".root";
      
      cout<<"EventNumber = "<<EventNumber<<endl;
      
      //Get pedestal subtracted event
      
      TFile *PedestalFile = new TFile(pedfilepath.str().c_str());
      TTree *PedestalTree = (TTree*)PedestalFile->Get("PedTree");
      if (!PedestalTree)
      {
         cerr << "PedTree not found!" << endl;
         return -1;
      }
      
      //event storage; for events in pedtree
      UInt_t run_number;
      UInt_t event_number;
      Float_t store_threshold;
      UInt_t Det_NChannels[9];
      UChar_t Det_Channels[9][256];
      UChar_t Det_ADC[8][256];
      UShort_t Dia_ADC[256];
      Float_t Det_PedMean[9][256];
      Float_t Det_PedWidth[9][256];
      
      //Event Header Branches
//      PedestalTree->SetBranchAddress("RunNumber",&run_number);
      PedestalTree->SetBranchAddress("EventNumber",&event_number);
      PedestalTree->SetBranchAddress("StoreThreshold",&store_threshold);
   
      //Telescope Data Branches
      PedestalTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
      PedestalTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
      PedestalTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
      PedestalTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
      PedestalTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
      PedestalTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
      PedestalTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
      PedestalTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
      PedestalTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
      PedestalTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
      PedestalTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
      PedestalTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
      PedestalTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
      PedestalTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
      PedestalTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
      PedestalTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
      PedestalTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
      PedestalTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
      PedestalTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
      PedestalTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
      PedestalTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
      PedestalTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
      PedestalTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
      PedestalTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
      PedestalTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
      PedestalTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
      PedestalTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
      PedestalTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
      PedestalTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
      PedestalTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
      PedestalTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
      PedestalTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
      PedestalTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
      PedestalTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
      PedestalTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
      PedestalTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
      PedestalTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
      PedestalTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
      PedestalTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
      PedestalTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
      PedestalTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
      PedestalTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
      PedestalTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
      PedestalTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
      PedestalTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
      
      //print stats
      cout << PedestalTree->GetEntries() << " events in PedTree " << endl;
      PedestalTree->GetEvent(0);
      cout<<"First event in PedTree is "<<event_number<<endl;
      PedestalTree->GetEvent(PedestalTree->GetEntries()-1);
      cout<<"Last event in PedTree is "<<event_number<<endl;
      
      //try to locate event
      PedestalTree->GetEvent(0);
      for(int guessed_event = EventNumber-event_number; guessed_event>0; guessed_event--) {
         PedestalTree->GetEvent(guessed_event);
         //cout<<"event_number="<<event_number<<endl;
         if(event_number==EventNumber) break;
      }
      cout<<"Event "<<event_number<<" loaded"<<endl;
      
      PedestalFile->Close();
      
      //TApplication theApp("App", &argc, argv);
      
      TCanvas *PedSubADCCanvas = new TCanvas("PedSubADCCanvas", "Telescope Reference Detectors and Diamond Detector PedSub ADC Output",200, 400, 1200, 900);
      PedSubADCCanvas->Divide(2,5);
      //PedSubADCCanvas->Connect("Closed()", "TApplication", &theApp, "Terminate()");
      
      //For Telescope Reference Detectors Plotting a Graph of ADC_value versus Strip_number
      
      
      TH1F *PedSubADCHisto[10];
      TH1F *Ped3WidthHisto[10];
      TH1F *Ped5WidthHisto[10];
      vector<TH1F*> ClusterPSADCHisto[9];
      string detname[10];
      detname[0] = "D0X";
      detname[1] = "D0Y";
      detname[2] = "D1X";
      detname[3] = "D1Y";
      detname[4] = "D2X";
      detname[5] = "D2Y";
      detname[6] = "D3X";
      detname[7] = "D3Y";
      detname[8] = "Dia0";
      detname[9] = "Dia1";
      string detname2[10];
      detname2[0] = "2D0X";
      detname2[1] = "2D0Y";
      detname2[2] = "2D1X";
      detname2[3] = "2D1Y";
      detname2[4] = "2D2X";
      detname2[5] = "2D2Y";
      detname2[6] = "2D3X";
      detname2[7] = "2D3Y";
      detname2[8] = "2Dia0";
      detname2[9] = "2Dia1";
      string detname3[10];
      detname3[0] = "3D0X";
      detname3[1] = "3D0Y";
      detname3[2] = "3D1X";
      detname3[3] = "3D1Y";
      detname3[4] = "3D2X";
      detname3[5] = "3D2Y";
      detname3[6] = "3D3X";
      detname3[7] = "3D3Y";
      detname3[8] = "3Dia0";
      detname3[9] = "3Dia1";
      
      
      for(int det=0; det<9; det++) {
         
         //pedsub data
         PedSubADCCanvas->cd(det+1);
         if(det==8) {
            PedSubADCHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),128,-0.5,127.5);
            Ped3WidthHisto[det] = new TH1F(detname2[det].c_str(), detname2[det].c_str(), 128, -0.5, 127.5);
            Ped5WidthHisto[det] = new TH1F(detname3[det].c_str(), detname3[det].c_str(), 128, -0.5, 127.5);
            for(int i=0; i<128; i++) 
               PedSubADCHisto[det]->SetBinContent(i+1,0);
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(uint i=0; i<Det_NChannels[det]; i++) {
               PedSubADCHisto[det]->SetBinContent((int)Det_Channels[det][i]+1,float(Dia_ADC[i]-Det_PedMean[det][i]));
               Ped3WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 3 * float(Det_PedWidth[det][i]));
               Ped5WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 5 * float(Det_PedWidth[det][i]));
               if(verbose) cout<<i<<"\t"<<(int)Det_Channels[det][i]<<"\t"<<Dia_ADC[i]<<"\t"<<Det_PedMean[det][i]<<"\t"<<Det_PedWidth[det][i]<<"\t"<<endl;
            }
         }
         else {
            PedSubADCHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),256,-0.5,255.5);
            Ped3WidthHisto[det] = new TH1F(detname2[det].c_str(), detname2[det].c_str(), 256, -0.5, 255.5);
            Ped5WidthHisto[det] = new TH1F(detname3[det].c_str(), detname3[det].c_str(), 256, -0.5, 255.5);
            for(int i=0; i<256; i++) 
               PedSubADCHisto[det]->SetBinContent(i+1,0);
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(uint i=0; i<Det_NChannels[det]; i++) {
               PedSubADCHisto[det]->SetBinContent((int)Det_Channels[det][i]+1,float(Det_ADC[det][i]-Det_PedMean[det][i]));
               Ped3WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 3 * float(Det_PedWidth[det][i]));
               Ped5WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 5 * float(Det_PedWidth[det][i]));
               if(verbose) cout<<i<<"\t"<<(int)Det_Channels[det][i]<<"\t"<<(int)Det_ADC[det][i]<<"\t"<<Det_PedMean[det][i]<<"\t"<<Det_PedWidth[det][i]<<"\t"<<endl;
            }
         }
         
         //clustered data
         ostringstream tempdetname;
         for(uint n=0; n<clustered_event.GetNClusters(det); n++) {
            tempdetname << detname[det] << "_cluster_" << n;
            ClusterPSADCHisto[det].push_back(new TH1F(tempdetname.str().c_str(),tempdetname.str().c_str(),256,-0.5,255.5));
            tempdetname.clear();
            ClusterPSADCHisto[det][n]->SetLineColor(n+2);
            ClusterPSADCHisto[det][n]->SetFillColor(n+2);
            ClusterPSADCHisto[det][n]->SetLineStyle(3);
            ClusterPSADCHisto[det][n]->SetFillStyle(3001);
            if(clustered_event.GetCluster(det,n)->IsGoldenGateCluster() || clustered_event.GetCluster(det,n)->IsBadChannelCluster()) {
               //ClusterPSADCHisto[det][n]->SetLineStyle(2);
               ClusterPSADCHisto[det][n]->SetFillStyle(0);
            }
            
            for(int i=0; i<256; i++) 
               ClusterPSADCHisto[det][n]->SetBinContent(i+1,0);
            
            for(int i=0; i<clustered_event.GetCluster(det,n)->GetNHits(); i++) {
               ClusterPSADCHisto[det][n]->SetBinContent(clustered_event.GetCluster(det,n)->GetChannel(i)+1,float(clustered_event.GetCluster(det,n)->GetADC(i)-clustered_event.GetCluster(det,n)->GetPedMean(i)));
            }
         }
         
         PedSubADCHisto[det]->Draw();
         Ped3WidthHisto[det]->SetLineColor(4);
         Ped5WidthHisto[det]->SetLineColor(2);
         Ped3WidthHisto[det]->Draw("same");
         Ped5WidthHisto[det]->Draw("same");
         PedSubADCHisto[det]->Draw("same");
         PedSubADCHisto[det]->GetXaxis()->SetTitle("Channel");
         PedSubADCHisto[det]->GetXaxis()->CenterTitle();
         PedSubADCHisto[det]->GetYaxis()->SetTitle("PedSub ADC Value");
         PedSubADCHisto[det]->GetYaxis()->CenterTitle();
         detname[det] += " PedSub ADC";
         PedSubADCHisto[det]->SetTitle(detname[det].c_str());
         for(uint n=0; n<clustered_event.GetNClusters(det); n++) 
            ClusterPSADCHisto[det][n]->Draw("same");
      }
      PedSubADCCanvas->Update();
      //theApp.Run();
      return 0;
   }
   
   else {
      //cout<<"Usage:\n"<<argv[0]<<" RunNumber EventNumber (verbose)"<<endl<<endl;
      return -1;
   }
}//End of Eventanalyze_final


