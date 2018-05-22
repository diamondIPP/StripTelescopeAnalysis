//Display pedestal subtracted events from PedestalSubtracted.root (need to run PedestalAnalyze first)
//Compile: g++ DisplaySNREvent.cpp -o DisplaySNREvent `root-config --cflags --glibs`
//2010-07-22 First verison based on DisplayRawEvent.cpp
//2010-08-01 Taylor fixed SetBinContent channel offset problem 
//2010-09-16 Fixed event lookup

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <cstdlib>
using namespace std;

// argv[]={DisplayPedMean, int RunNumber, int EventNumber}
int main(int argc, char* argv[]) {
   
   if(argc==3 || argc==4) {
      int RunNumber = (int)strtod(argv[1],0);
      int EventNumber = (int)strtod(argv[2],0);
      bool verbose = 0;
      if(argc==4) verbose = (int)strtod(argv[3],0);
      
      cout<<"RunNumber = "<<RunNumber<<"\tEventNumber = "<<EventNumber<<endl;
      
      //Get pedestal subtracted event
      
      TFile *PedestalFile = new TFile("PedestalSubtracted.root");
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
      
      TCanvas *PedMeanCanvas = new TCanvas("PedMeanCanvas", "Telescope Reference Detectors and Diamond Detector Pedestal Means",200, 400, 1200, 900);
      PedMeanCanvas->Divide(2,5);
      PedMeanCanvas->Connect("Closed()", "TApplication", &theApp, "Terminate()");
      
      //For Telescope Reference Detectors Plotting a Graph of ADC_value versus Strip_number
      
      TH1F *PedMeanHisto[10];
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
      
      for(int det=0; det<9; det++) {
         PedMeanCanvas->cd(det+1);
         if(det==8) {
            PedMeanHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),128,-0.5,127.5);
            for(int i=0; i<128; i++) 
               PedMeanHisto[det]->SetBinContent(i + 1, 0);
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(int i=0; i<Det_NChannels[det]; i++) {
               PedMeanHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, float(Det_PedMean[det][i]));
               if(verbose) cout<<i<<"\t"<<(int)Det_Channels[det][i]<<"\t"<<Dia_ADC[i]<<"\t"<<Det_PedMean[det][i]<<"\t"<<Det_PedWidth[det][i]<<"\t"<<endl;
            }
         }
         else {
            PedMeanHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),256,-0.5,255.5);
            for(int i=0; i<256; i++) 
               PedMeanHisto[det]->SetBinContent(i + 1, 0);
            if(verbose) {
               cout<<endl<<"Detector "<<det<<endl<<endl;
               cout<<"index\tchan\tadc\tpedmean\tpedrms"<<endl;
            }
            for(int i=0; i<Det_NChannels[det]; i++) {
               PedMeanHisto[det]->SetBinContent((int)Det_Channels[det][i] + 1, float(Det_PedMean[det][i]));
               if(verbose) cout<<i<<"\t"<<(int)Det_Channels[det][i]<<"\t"<<(int)Det_ADC[det][i]<<"\t"<<Det_PedMean[det][i]<<"\t"<<Det_PedWidth[det][i]<<"\t"<<endl;
            }
         }
         PedMeanHisto[det]->Draw();
         PedMeanHisto[det]->GetXaxis()->SetTitle("Channel");
         PedMeanHisto[det]->GetXaxis()->CenterTitle();
         PedMeanHisto[det]->GetYaxis()->SetTitle("Pedestal Mean");
         PedMeanHisto[det]->GetYaxis()->CenterTitle();
         detname[det] += " PedMean";
         PedMeanHisto[det]->SetTitle(detname[det].c_str());
         PedMeanCanvas->cd(det+1);
      }
      PedMeanCanvas->Update();
      theApp.Run();
   
      return 0;
   }
   
   else {
      cout<<"Usage:\n"<<argv[0]<<" RunNumber EventNumber (verbose)"<<endl<<endl;
      return -1;
   }
}//End of Eventanalyze_final


