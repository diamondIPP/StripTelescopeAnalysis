//Display pedestal subtracted events from PedestalSubtracted.root (need to run PedestalAnalyze first)
//Compile: g++ DisplaySNREvent.cpp -o DisplaySNREvent `root-config --cflags --glibs`
//2010-07-22 First verison based on DisplayRawEvent.cpp
//2010-08-01 Taylor added error line
//2010-09-16 Fixed event lookup

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TSystem.h"
#include <iostream>
#include <cstdlib>
#include <sstream>

using namespace std;


// argv[]={DisplayPedMean, int RunNumber, int EventNumber}
int main(int argc, char* argv[]) {
   
   if(argc==3 || argc==4) {
      int RunNumber = (int)strtod(argv[1],0);
      int EventNumber = (int)strtod(argv[2],0);
      bool verbose = 0;
      if(argc==4) verbose = (int)strtod(argv[3],0);
      
      cout<<"RunNumber = "<<RunNumber<<"\tEventNumber = "<<EventNumber<<endl;

//void DisplayRawEvent(int RunNumber, int EventNumber, bool show_unflipped_plots=0) {
   
      //Get pedestal subtracted event
         
      TSystem* sys = gSystem;
      ostringstream pedfilepath;
      pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
      //if(RunDescription=="") pedfilepath << ".root";
      //else pedfilepath << "-" << RunDescription << ".root";
      pedfilepath << ".root";
      //TFile *PedestalFile = new TFile("PedestalSubtracted.root");
      TFile *PedestalFile = new TFile(pedfilepath.str().c_str());
      TTree *PedTree = (TTree*)PedestalFile->Get("PedTree");
      if (!PedTree)
      {
         cerr << "Tree not found!" << endl;
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
//      PedTree->SetBranchAddress("RunNumber",&run_number);
      PedTree->SetBranchAddress("EventNumber",&event_number);
      PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
   
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
      
      //print stats
      cout << PedTree->GetEntries() << " events in PedTree " << endl;
      PedTree->GetEvent(0);
      cout<<"First event in PedTree is "<<event_number<<endl;
      PedTree->GetEvent(PedTree->GetEntries()-1);
      cout<<"Last event in PedTree is "<<event_number<<endl;
      
      //try to locate event
      PedTree->GetEvent(0);
      for(int guessed_event = EventNumber-event_number; guessed_event>0; guessed_event--) {
         PedTree->GetEvent(guessed_event);
         //cout<<"event_number="<<event_number<<endl;
         if(event_number==EventNumber) break;
      }
      cout<<"Event "<<event_number<<" loaded"<<endl;
      
      PedestalFile->Close();
      
      TApplication theApp("App", &argc, argv);
      
      TCanvas *PedSubADCCanvas = new TCanvas("PedSubADCCanvas", "Telescope Reference Detectors and Diamond Detector PedSub ADC Output",200, 400, 1200, 900);
      PedSubADCCanvas->Divide(2,5);
      PedSubADCCanvas->Connect("Closed()", "TApplication", &theApp, "Terminate()");
      
      //For Telescope Reference Detectors Plotting a Graph of ADC_value versus Strip_number
      
      TH1F *PedSubADCHisto[10];
      TH1F *Ped3WidthHisto[10];
      TH1F *Ped5WidthHisto[10];
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
         PedSubADCCanvas->cd(det+1);
         if(det==8) {
            PedSubADCHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),128,-0.5,127.5);
            Ped3WidthHisto[det] = new TH1F(detname2[det].c_str(), detname2[det].c_str(), 128, -0.5, 127.5);
            Ped5WidthHisto[det] = new TH1F(detname3[det].c_str(), detname3[det].c_str(), 128, -0.5, 127.5);
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(int i=0; i<Det_NChannels[det]; i++) {
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
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(int i=0; i<Det_NChannels[det]; i++) {
               PedSubADCHisto[det]->SetBinContent((int)Det_Channels[det][i]+1,float(Det_ADC[det][i]-Det_PedMean[det][i]));
               Ped3WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 3 * float(Det_PedWidth[det][i]));
               Ped5WidthHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, 5 * float(Det_PedWidth[det][i]));
               if(verbose) cout<<i<<"\t"<<(int)Det_Channels[det][i]<<"\t"<<(int)Det_ADC[det][i]<<"\t"<<Det_PedMean[det][i]<<"\t"<<Det_PedWidth[det][i]<<"\t"<<endl;
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
         PedSubADCCanvas->cd(det+1);
      }
      PedSubADCCanvas->Update();
      theApp.Run();
   
      return 0;
   }
   
   else {
      cout<<"Usage:\n"<<argv[0]<<" RunNumber EventNumber (verbose)"<<endl<<endl;
      return -1;
   }
}//End of Eventanalyze_final


