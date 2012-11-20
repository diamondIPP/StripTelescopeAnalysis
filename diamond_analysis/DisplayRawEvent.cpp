//Display raw events from RZ files
//2010-07-22 Now compiles and runs as stand alone program :)
//2010-08-01 Taylor fixed SetBinContent channel offset problem 

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "GetRawEvent.cpp" //the header file that is connected to the Diamond/telescope data
#include <iostream>
#include <cstdlib>
using namespace std;


/*
// argv[]={DisplayRawEvent, int RunNumber, int EventNumber}
int main(int argc, char* argv[]) {

   int RunNumber = (int)strtod(argv[1],0);
   int EventNumber = (int)strtod(argv[2],0);
   bool show_unflipped_plots = 0;
   
   cout<<"RunNumber = "<<RunNumber<<"\tEventNumber = "<<EventNumber<<endl;
   
   TApplication theApp("App", &argc, argv);
   //TApplication theApp2("App", &argc, argv);
*/
void DisplayRawEvent(int RunNumber, int EventNumber, bool show_unflipped_plots=0) {
   
   RawEvent rawevent = GetRawEvent(RunNumber, EventNumber, true);
   
   TCanvas *RawADCCanvas = new TCanvas("RawADCCanvas", "Telescope Reference Detectors and Diamond Detector Raw ADC Output",200, 400, 1200, 900);
   TCanvas *RawADCCanvasUnflipped = 0;
   RawADCCanvas->Divide(2,5);
   //RawADCCanvas->Connect("Closed()", "TApplication", &theApp, "Terminate()");
   
   //For Telescope Reference Detectors Plotting a Graph of ADC_value versus Strip_number
   
   TH1F *RawADCHisto[15];
   string detname[15];
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
   detname[10] = "D0";
   detname[11] = "D1";
   detname[12] = "D2";
   detname[13] = "D3";
   detname[14] = "Diamond";
   
   for(int det=0; det<10; det++) {
      RawADCCanvas->cd(det+1);
      if(det==8||det==9) {
         RawADCHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),128,-0.5,127.5);
         for(int i=0; i<128; i++) 
            RawADCHisto[det]->SetBinContent(i + 1, rawevent.GetDetector(det).GetADC(i));
      }
      else {
         RawADCHisto[det] = new TH1F(detname[det].c_str(),detname[det].c_str(),256,-0.5,255.5);
         for(int i=0; i<256; i++) 
            RawADCHisto[det]->SetBinContent(i + 1, rawevent.GetDetector(det).GetADC(i));
      }
      RawADCHisto[det]->Draw();
      RawADCHisto[det]->GetXaxis()->SetTitle("Channel");
      RawADCHisto[det]->GetXaxis()->CenterTitle();
      RawADCHisto[det]->GetYaxis()->SetTitle("Raw ADC Value");
      RawADCHisto[det]->GetYaxis()->CenterTitle();
      detname[det] += " Raw ADC";
      RawADCHisto[det]->SetTitle(detname[det].c_str());
      RawADCCanvas->cd(det+1);
   }
   RawADCCanvas->Update();
   //theApp.Run();
   
   if(show_unflipped_plots) {
      RawADCCanvasUnflipped = new TCanvas("RawADCCanvasUnflipped", "Not Flipped Telescope Reference Detectors and Diamond Detector Raw ADC Output",200, 400, 1200, 900);
      RawADCCanvasUnflipped->Divide(2,2);//1,4);
      //RawADCCanvasUnflipped->Connect("Closed()", "TApplication", &theApp, "Terminate()");
      for(int det=0; det<4; det++) {
         RawADCCanvasUnflipped->cd(det+1);
         RawADCHisto[det+10] = new TH1F(detname[det+10].c_str(),detname[det+10].c_str(),512,-0.5,511.5);
         for(int i=0; i<512; i++) 
            RawADCHisto[det+10]->SetBinContent(i+1, rawevent.GetUnflippedDetector(det)[i]);
         RawADCHisto[det+10]->Draw();
         RawADCHisto[det+10]->GetXaxis()->SetTitle("Channel");
         RawADCHisto[det+10]->GetXaxis()->CenterTitle();
         RawADCHisto[det+10]->GetYaxis()->SetTitle("Raw ADC Value");
         RawADCHisto[det+10]->GetYaxis()->CenterTitle();
         RawADCHisto[det+10]->SetMinimum(0);
         detname[det+10] += " Raw ADC (Not Flipped)";
         RawADCHisto[det+10]->SetTitle(detname[det+10].c_str());
      }
   RawADCCanvasUnflipped->Update();
   }
   //theApp.Run();

   //return 0;
}//End of Eventanalyze_final


