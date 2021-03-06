/*
 * HistogrammSaver.class.cpp
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 *
 */

#include "HistogrammSaver.class.hh"

using namespace std;

HistogrammSaver::HistogrammSaver(TSettings * newSettings,int verbosity) {
    if(!newSettings){
        cerr<<"[HistogrammSaver::HistogrammSaver]: settings == NULL "<<endl;
        exit(-1);
    }

    this->settings = newSettings;
    sys=NULL;
    pt=NULL;
    this->verbosity=verbosity;
    runNumber=0;
    nEvents=0;
    hGoodCells = 0;
    hBadCells = 0;
    hAllCells = 0;
    hAllButBadCells = 0;
    plots_path=".";
    pt = new TPaveText(0.07,0,0.22,0.10,"NDC");
    UpdatePaveText();
    this->SetNumberOfEvents(newSettings->getNEvents());
    this->SetRunNumber(settings->getRunNumber());
    if(verbosity)cout<<"HistogrammSaver::HistogrammSaver:: get new TSystem"<<endl;
    //	sys=new TSystem();
    sys=gSystem;
    if(verbosity)cout<<"HistogrammSaver::HistogrammSaver:: Set Style"<<endl;
    currentStyle=gROOT->GetStyle("Plain_RD42");
    if(currentStyle!=0)
        currentStyle->cd();
    bool paperMode = settings->IsPaperMode();
    if(gStyle!=0)
        if(!gStyle->IsZombie()){

            if((string)gStyle->GetName()!="Plain_RD42"){
                gROOT->SetStyle("Plain"); //General style (see TStyle)

                //	    gStyle->SetOptStat(221111111); //Stat options to be displayed			without under- and overflow use gStyle->SetOptStat(1110);
                if(gStyle->GetOptStat()!=221111111){
                    if(paperMode)
                        gStyle->SetOptStat("emr");
                    else
                        gStyle->SetOptStat("nemr");//KSiou
                }

                if(gStyle->GetOptFit()!=1111){
                    gStyle->SetOptFit(1111);  //Fit options to be displayed
                    gStyle->SetStatH(0.12); //Sets Height of Stats Box
                    gStyle->SetStatW(0.15); //Sets Width of Stats Box
                }
                if(gStyle->GetPadBottomMargin()!=0.15) gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
                //	    gStyle->SetPadRightMargin(0.15);
                if(gStyle->GetPadTopMargin()!=0.15) gStyle->SetPadTopMargin(0.15);
                //gStyle->SetTitleColor(19,"");
                gStyle->SetPalette(53);
                currentStyle= (TStyle*)gStyle->Clone("Plain_RD42");
                currentStyle->SetPalette(53);
                currentStyle2D= (TStyle*)currentStyle->Clone("Plain_RD42_2D");
                currentStyle2D->SetOptStat("ne");
                currentStyle2D ->SetPalette(53);
                currentStyle->cd();
                //				gROOT->SetStyle("Plain_RD42");
            }
        }
//    if (paperMode) gStyle->SetOptTitle(false);
    gStyle->SetPalette(53); //
    SetPaperPlotStyle();
    if(verbosity)cout<<"HistogrammSaver::HistogrammSaver::Created instance of HistogrammSaver"<<endl;
    gErrorIgnoreLevel=3001;
    InitializeGridReferenceDetSpace();

}

HistogrammSaver::~HistogrammSaver() {
    DeleteCellHistos();
    //	TString string1 = sys->GetFromPipe(".! mkdir root-Files");
    //	cout<<string1<<endl;
    MoveRootFilesInSubfolder();
    //	string1 = sys->GetFromPipe(".!mv -v *.root root-Files");
    //	cout<<string1<<endl;
    if (pt) this->pt->Delete();
}


void HistogrammSaver::MoveRootFilesInSubfolder(){
    stringstream test;
    test << "find "<<plots_path<<" -maxdepth 1 -iname \"*.root\" -exec mv -f {} "<<plots_path<<"/root/  \\;";
    TString toPath = plots_path +"/root/";
    TString fromPath = plots_path +"/*.root";
    TString executeString = "mkdir -p "+toPath;
//    cout<<"MoveRootFilesInSubfolder: Execute: \""<<executeString<<"\", "<<flush;
    system(executeString);
    executeString ="mv  -f -u "+fromPath+" "+toPath;
//    cout<<"Execute: \""<<executeString<<"\",..."<<flush;
    system(executeString);//t.str();//<<"\""<<endl;
//    cout<<"DONE"<<endl;
}

void HistogrammSaver::SetPaperPlotStyle(){
    if (!settings->IsPaperMode())
        return;
    cout<<"SET PAPER PLOT STYLE!"<<endl;
//    if (!gStyle)
//      gStyle = new TStyle("gStyle","Style for P-TDR");
      // For the canvas:
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetCanvasColor(kWhite);
      gStyle->SetCanvasDefH(1000); //Height of canvas
      gStyle->SetCanvasDefW(1000); //Width of canvas
      gStyle->SetCanvasDefX(0);   //POsition on screen
      gStyle->SetCanvasDefY(0);

      // For the Pad:
      gStyle->SetPadBorderMode(0);
      // gStyle->SetPadBorderSize(Width_t size = 1);
      gStyle->SetPadColor(kWhite);
      gStyle->SetPadGridX(false);
      gStyle->SetPadGridY(false);
      gStyle->SetGridColor(0);
      gStyle->SetGridStyle(3);
      gStyle->SetGridWidth(1);

      // For the frame:
      gStyle->SetFrameBorderMode(0);
      gStyle->SetFrameBorderSize(1);
      gStyle->SetFrameFillColor(0);
      gStyle->SetFrameFillStyle(0);
      gStyle->SetFrameLineColor(1);
      gStyle->SetFrameLineStyle(1);
      gStyle->SetFrameLineWidth(1);

      // For the histo:
//      gStyle->SetHistFillColor(63);
      gStyle->SetHistFillStyle(0);
      gStyle->SetHistLineColor(1);
      gStyle->SetHistLineStyle(0);
      gStyle->SetHistLineWidth(2);
      // gStyle->SetLegoInnerR(Float_t rad = 0.5);
      // gStyle->SetNumberContours(Int_t number = 20);

    //  gStyle->SetEndErrorSize(0);
      gStyle->SetErrorX(0.);
    //  gStyle->SetErrorMarker(20);

      gStyle->SetMarkerStyle(20);

      //For the fit/function:
      gStyle->SetOptFit(1);
      gStyle->SetFitFormat("5.4g");
      gStyle->SetFuncColor(2);
      gStyle->SetFuncStyle(1);
      gStyle->SetFuncWidth(1);

      //For the date:
      gStyle->SetOptDate(0);
      // gStyle->SetDateX(Float_t x = 0.01);
      // gStyle->SetDateY(Float_t y = 0.01);

      // For the statistics box:
      gStyle->SetOptFile(0);
      gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
      gStyle->SetStatColor(kWhite);
      gStyle->SetStatFont(42);
      gStyle->SetStatFontSize(0.025);
      gStyle->SetStatTextColor(1);
      gStyle->SetStatFormat("6.4g");
      gStyle->SetStatBorderSize(1);
      gStyle->SetStatH(0.1);
      gStyle->SetStatW(0.15);
      // gStyle->SetStatStyle(Style_t style = 1001);
      // gStyle->SetStatX(Float_t x = 0);
      // gStyle->SetStatY(Float_t y = 0);

      // Margins:
      gStyle->SetPadTopMargin(0.13);
      gStyle->SetPadBottomMargin(0.2);
      gStyle->SetPadLeftMargin(0.13);
      gStyle->SetPadRightMargin(0.2);

      // For the Global title:

      //  gStyle->SetOptTitle(0);
      gStyle->SetTitleFont(42);
      gStyle->SetTitleColor(kBlack);
      gStyle->SetTitleTextColor(1);
      gStyle->SetTitleFillColor(10);
      gStyle->SetTitleFontSize(0.05);
      gStyle->SetOptTitle(0);
      // gStyle->SetTitleH(0); // Set the height of the title box
      // gStyle->SetTitleW(0); // Set the width of the title box
      // gStyle->SetTitleX(0); // Set the position of the title box
      // gStyle->SetTitleY(0.985); // Set the position of the title box
      // gStyle->SetTitleStyle(Style_t style = 1001);
      // gStyle->SetTitleBorderSize(2);

      // For the axis titles:

      gStyle->SetTitleColor(kBlack, "XYZ");
      gStyle->SetTitleFont(42, "XYZ");
      gStyle->SetTitleSize(0.04, "XYZ");
      // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
      // gStyle->SetTitleYSize(Float_t size = 0.02);
      gStyle->SetTitleXOffset(1.1);
      gStyle->SetTitleYOffset(1.3);
      gStyle->SetTitleOffset(1.4, "Z"); // Another way to set the Offset

      // For the axis labels:

      gStyle->SetLabelColor(kBlack, "XYZ");
      gStyle->SetLabelFont(42, "XYZ");
      gStyle->SetLabelOffset(0.010, "XYZ");
      gStyle->SetLabelOffset(0.015, "Z");
      gStyle->SetLabelSize(0.045, "XYZ");

      // For the axis:

      gStyle->SetAxisColor(1, "XYZ");
      gStyle->SetStripDecimals(kTRUE);
      gStyle->SetTickLength(0.03, "XYZ");
      gStyle->SetNdivisions(510, "XYZ");
      gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
      gStyle->SetPadTickY(1);

      // Change for log plots:
      gStyle->SetOptLogx(0);
      gStyle->SetOptLogy(0);
      gStyle->SetOptLogz(0);

      //Legend
      gStyle->SetLegendFont(42);

      // Postscript options:
      // gStyle->SetPaperSize(15.,15.);
      // gStyle->SetLineScalePS(Float_t scale = 3);
      // gStyle->SetLineStyleString(Int_t i, const char* text);
      // gStyle->SetHeaderPS(const char* header);
      // gStyle->SetTitlePS(const char* pstitle);

      // gStyle->SetBarOffset(Float_t baroff = 0.5);
      // gStyle->SetBarWidth(Float_t barwidth = 0.5);
      // gStyle->SetPaintTextFormat(const char* format = "g");
      // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
      // gStyle->SetTimeOffset(Double_t toffset);
      // gStyle->SetHistMinimumZero(kTRUE);

      gStyle->cd();
}

void HistogrammSaver::InitializeGridReferenceDetSpace(){
    TString nameDet = "hGridRefenrenceDetSpace";
    TString nameCell = "hGridRefenrenceCellSpace";
    Float_t xBins = settings->getNColumns3d();
    TFidCutRegions* metallisationFidCuts = settings->get3dMetallisationFidCuts();
    if(settings->is3dDiamond())metallisationFidCuts->Print(1);
    TFiducialCut* fidCut = metallisationFidCuts->getFidCut((UInt_t)3);
    Float_t xLow = 0;
    Float_t xHigh = 0;
    Float_t yBins = 0;
    Float_t yLow = 0;
    Float_t yHigh = 0;
    if(fidCut){
        xLow = fidCut->GetXLow();//getXMetalisationStart3d;
        xHigh = fidCut->GetXHigh();//getXMetalisationEnd3d;
        yBins = settings->getNRows3d();
        yLow = fidCut->GetYLow();
        yHigh = fidCut->GetYHigh();//getYMetalisationEnd3d;
    }
    //	cout<<"nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh"<<endl;
    //	cout<<nameDet<<" "<<nameDet<<" "<<xBins<<" "<<xLow<<" "<<xHigh<<" "<<yBins<<" "<<yLow<<" "<<yHigh<<endl;
    hGridReferenceDetSpace = new TH2D(nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh);
    hGridReferenceCellSpace = new TH2D(nameCell,nameCell,xBins,0,xBins,yBins,0,yBins);

    for(UInt_t i=0;i<settings->getNColumns3d();i++){
        hGridReferenceDetSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
        hGridReferenceCellSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
    }
    for(UInt_t j=0;j<settings->getNRows3d();j++){
        hGridReferenceDetSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
        hGridReferenceCellSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
    }
    hGridReferenceDetSpace->SetStats(kFALSE);
    hGridReferenceDetSpace->SetTickLength(0.0, "X");
    hGridReferenceDetSpace->SetTickLength(0.0, "Y");
    hGridReferenceDetSpace->GetXaxis()->SetTitle("column of cell");
    hGridReferenceDetSpace->GetYaxis()->SetTitle("row of cell");
    hGridReferenceCellSpace->SetStats(kFALSE);
    hGridReferenceCellSpace->SetTickLength(0.0, "X");
    hGridReferenceCellSpace->SetTickLength(0.0, "Y");
}


void HistogrammSaver::SaveTwoHistosNormalized(TString canvasName, TH1 *histo1, TH1 *histo2,double refactorSecond, TString position, UInt_t verbosity){
//    cout<<"Save2HistosNormalized: "<<histo1<<" "<<histo2<<endl;
    bool internalVerbose = false;
    if(!histo1&&!histo2)return;
    if(!histo1||!histo2){
        if (histo1) SaveHistogram(histo1);
        else SaveHistogram(histo2);
        return;
    }
    cout<<"Save2HistosNormalized: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
    TCanvas *c1 = new TCanvas(canvasName,canvasName);
    cout<<"new Canvas: "<<c1->GetName()<<endl;
    c1->cd();
    c1->SetObjectStat(false);
    Float_t min1 = histo1->GetMinimum()/histo1->Integral();;
    Float_t min2 = histo2->GetMinimum()/histo2->Integral();;
    Float_t min = TMath::Min(min1,min2);
    Float_t max1 = histo1->GetMaximum()/histo1->Integral();
    Float_t max2 = histo2->GetMaximum()/histo2->Integral();;

    //	Float_t range1 = max1-min1;
    //	Float_t range2 = max2-min2;
    Float_t max = TMath::Max(max1,max2);
    Float_t range = max - min;
    Float_t middle = (max+min)/2.;
    if(min>=0&&(middle - range/2.*1.1)<0)
        min =0;
    else
        min = middle - range/2.*1.1;
    max = middle + range/2.*1.4;
    //	int stat = gStyle->GetOptStat();
    if(max2*refactorSecond>max1)
        refactorSecond=max2/max1*0.5;
    if(refactorSecond!=1)histo2->Scale(refactorSecond);
    if (verbosity>2||internalVerbose) cout<<"min: "<<min<<" max: "<<max;
    if (verbosity>2||internalVerbose) cout<<" refactorSecond:"<<refactorSecond<<"\thisto1:"<<max1<<"\thisto2:"<<max2<<flush;
    if (verbosity>2||internalVerbose) cout<<endl<<"Nhisto1: "<<histo1->GetEntries()<<" Nhisto2:"<<histo2->GetEntries()<<flush;
    histo1->SetStats(false);
    histo2->SetStats(false);
    TH1F* histo1Normalized;
    TH1F* histo2Normalized;
    if(max1>max2){
        if (verbosity>2||internalVerbose) cout<<"\tdraw1-"<<flush;
        histo1Normalized = (TH1F*)((TH1F*)(histo1->Clone()))->DrawNormalized("");
        if(histo1Normalized)
            histo1Normalized->GetYaxis()->SetRangeUser(min,max);
        else
            histo1->GetYaxis()->SetRangeUser(min,max);
        if (verbosity>2||internalVerbose) cout<<"draw2 "<<flush;
        histo2Normalized = (TH1F*)((TH1F*)(histo2->Clone()))->DrawNormalized("same");
        //		histo2->GetYaxis()->SetRangeUser(min,max);
    }
    else{
        if (verbosity>2||internalVerbose) cout<<"\tdraw2-"<<flush;
        histo2Normalized = (TH1F*)((TH1F*)(histo2->Clone()))->DrawNormalized("");
        if (histo2Normalized)
            histo2Normalized->GetYaxis()->SetRangeUser(min,max);
        else
            histo2->GetYaxis()->SetRangeUser(min,max);
        if (verbosity>2||internalVerbose) cout<<"draw1 "<<flush;
        histo1Normalized =  (TH1F*)((TH1F*)(histo1->Clone()))->DrawNormalized("same");
        //		histo1->GetYaxis()->SetRangeUser(min,max);
    }
    c1->Update();
    //	TVirtualPad *pad =c1->GetPad(0);
    if (verbosity>2||internalVerbose) cout<<"MIN: "<<min<<"-->";
    min=(double)(min/refactorSecond);
    if (verbosity>2||internalVerbose) cout<<min<<"\t\tMAX: "<<max<<"--->";
    max = (double)(max/refactorSecond);
    if (verbosity>2||internalVerbose) cout<<max<<endl;
    c1->Update();
    TLegend *leg;
    position.ToLower();
    if (position== "right")
        leg = new TLegend(0.52,0.75,0.9,0.95);
    else
        leg =new TLegend(0.1,0.75,0.4,0.9);
    leg->SetFillColor(kWhite);
    leg->SetHeader("Legend");
    //	if(histo1Normalized)
    //		leg->AddEntry(histo1Normalized,histo1Normalized->GetName());
    //	else
    if(histo1&&!histo1->IsZombie())
        leg->AddEntry(histo1,histo1->GetName());
    //	if(histo2Normalized)
    //		leg->AddEntry(histo2Normalized,histo2Normalized->GetName());
    //	else
    if(histo2&&!histo2->IsZombie())
        leg->AddEntry(histo2,histo2->GetName());
    leg->Draw("same");
    //	TPaveText* pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvasName.Data()));
    //	if(!settings->IsPaperMode()) pt2->Draw("same");
    c1->Update();
    cout<<"Save Canvas "<< c1<<": "<<c1->GetName()<<endl;
    SaveCanvas(c1);
}

void HistogrammSaver::SaveTwoHistosScaled(TString canvasName, TH1 *histo1, TH1 *histo2, TString position, TString drawOption, UInt_t verbosity){
    if(!histo1||!histo2){
        if (histo1) SaveHistogram(histo1);
        if (histo2) SaveHistogram(histo2);
        return;
    }
    Float_t scale1 = 1./histo1->GetBinContent(histo1->GetMaximumBin());
    Float_t scale2 = 1./histo2->GetBinContent(histo2->GetMaximumBin());
    this->SaveTwoHistosComparison(canvasName,histo1,histo2,scale1,scale2,position,verbosity);

}
void HistogrammSaver::SaveTwoHistosComparison(TString canvasName, TH1 *histo1, TH1 *histo2,double scale1, double scale2, TString position, TString drawOption, UInt_t verbosity){
    if(!histo1&&!histo2)return;
      if(!histo1||!histo2){
          if (histo1) SaveHistogram(histo1);
          else SaveHistogram(histo2);
          return;
      }
      histo1 = (TH1*)histo1->Clone();
      histo2 = (TH1*)histo2->Clone();
      histo1->Scale(scale1);
      histo2->Scale(scale2);
      if(histo1->GetLineColor() == histo2->GetLineColor())
          histo2->SetLineColor(histo1->GetLineColor()+1);
      cout<<"SaveTwoHistosComparison: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
      TCanvas *c1 = new TCanvas(canvasName,canvasName);
      c1->cd();
      c1->SetObjectStat(true);
      TString title = canvasName;
      title+= ";";
      title+= histo1->GetXaxis()->GetTitle();
      title+= ";";
      title+= histo1->GetYaxis()->GetTitle();
      THStack* stack = new THStack(canvasName,title);
      stack->Add(histo1);
      stack->Add(histo2);
      histo1->Draw("goff");
      histo2->Draw("goff");
      if (verbosity>2) cout<<endl<<"Nhisto1: "<<histo1->GetEntries()<<" Nhisto2:"<<histo2->GetEntries()<<flush;
      histo1->SetStats(false);
      histo2->SetStats(false);
      stack->Draw("nostack");
      c1->Update();
      TLegend *leg;
      position.ToLower();
      if (position== "right")
          leg = new TLegend(0.52,0.75,0.9,0.95);
      else
          leg = new TLegend(0.1,0.75,0.48,0.9);
      leg->SetFillColor(kWhite);
//      leg->SetHeader("Legend");
      leg->AddEntry(histo1,histo1->GetName(),"LP");
      leg->AddEntry(histo2,histo2->GetName(),"LP");
      leg->Draw("same");
      if (drawOption.Contains("logy"))
          c1->SetLogy();
      if (drawOption.Contains("logx"))
                c1->SetLogx();
      c1->Update();
      SaveCanvas(c1);
}
void HistogrammSaver::SaveTwoHistos(TString canvasName, TH1 *histo1, TH1 *histo2,double refactorSecond, TString position, UInt_t verbosity)
{
//    cout<<"Save2Histos: "<<histo1<<" "<<histo2<<endl;
    if(!histo1&&!histo2)return;
    if(!histo1||!histo2){
        if (histo1) SaveHistogram(histo1);
        else SaveHistogram(histo2);
        return;
    }
    histo1 = (TH1*)histo1->Clone();
    histo2 = (TH1*)histo2->Clone();
    if(histo1->GetLineColor() == histo2->GetLineColor())
        histo2->SetLineColor(histo1->GetLineColor()+1);
    cout<<"Save2Histos: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
    if (verbosity>2) cout<<"Save2Histos: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
    TCanvas *c1 = new TCanvas(canvasName,canvasName);
    c1->cd();
    c1->SetObjectStat(true);
    Float_t min1 = histo1->GetBinContent(histo1->GetMinimumBin());//GetMinimum();
    Float_t min2 = histo2->GetBinContent(histo2->GetMinimumBin());
    Float_t min = TMath::Min(min1,min2);
    Float_t max1 =  histo1->GetBinContent(histo1->GetMaximumBin());//GetMinimum();
    Float_t max2 = histo2->GetBinContent(histo2->GetMaximumBin());
    Float_t max = TMath::Max(max1,max2);
    Float_t factor = 1.1;
    Float_t range = max - min;
    Float_t middle = (max+min)/2.;
    if(min>=0&&(middle - range/2.*factor)<0)
        min =0;
    else
        min = middle - range/2.*factor;
    max = middle + range/2.*factor;
    //	int stat = gStyle->GetOptStat();
    if(refactorSecond!=1&&histo2->GetMaximum()*refactorSecond>histo1->GetMaximum())
        refactorSecond=histo2->GetMaximum()/histo1->GetMaximum()*0.5;
    histo1->Draw("goff");
    histo2->Draw("goff");
    Float_t xmin1 = histo1->GetXaxis()->GetBinLowEdge(histo1->GetXaxis()->GetFirst());
    Float_t xmin2 = histo2->GetXaxis()->GetBinLowEdge(histo2->GetXaxis()->GetFirst());
    Float_t xmax1 = histo1->GetXaxis()->GetBinLowEdge(histo1->GetXaxis()->GetLast());
    Float_t xmax2 = histo2->GetXaxis()->GetBinLowEdge(histo2->GetXaxis()->GetLast());
    Float_t xmin = TMath::Min(xmin1,xmin2);
    Float_t xmax = TMath::Max(xmax1,xmax2);

    histo1->GetXaxis()->SetRangeUser(xmin,xmax);
    histo2->GetXaxis()->SetRangeUser(xmin,xmax);
    if(refactorSecond!=1)histo2->Scale(refactorSecond);
    if (verbosity>2) cout<<"min: "<<min<<" max: "<<max;
    if (verbosity>2) cout<<" refactorSecond:"<<refactorSecond<<"\thisto1:"<<histo1->GetMaximum()<<"\thisto2:"<<histo2->GetMaximum()<<flush;
    if (verbosity>2) cout<<endl<<"Nhisto1: "<<histo1->GetEntries()<<" Nhisto2:"<<histo2->GetEntries()<<flush;
    histo1->SetStats(false);
    histo2->SetStats(false);
    Float_t ymin = TMath::Min(histo1->GetBinContent(histo1->GetMinimumBin()),
            histo2->GetBinContent(histo2->GetMinimumBin()));
    ymin = TMath::Min(ymin,min);
    Float_t ymax = TMath::Max(histo1->GetBinContent(histo1->GetMaximumBin()),
            histo2->GetBinContent(histo2->GetMaximumBin()));
    TH1F* frame = c1->DrawFrame(xmin,xmax,ymin,ymax);
    frame->SetTitle("frame");
    frame->GetXaxis()->SetTitle(histo1->GetXaxis()->GetTitle());
    frame->GetYaxis()->SetTitle(histo1->GetYaxis()->GetTitle());
    frame->Draw();
    histo1->Draw("same");
    histo2->Draw("same");
    if(histo1->GetMaximum()>histo2->GetMaximum()){
        if (verbosity>2) cout<<"\tdraw1-"<<flush;
        histo1->Draw("");
        histo1->GetYaxis()->SetRangeUser(min,max);
        if (verbosity>2) cout<<"draw2 "<<flush;
        histo2->Draw("same");
        //		histo2->GetYaxis()->SetRangeUser(min,max);
    }
    else{
        if (verbosity>2) cout<<"\tdraw2-"<<flush;
        histo2->Draw("");
        histo2->GetYaxis()->SetRangeUser(min,max);
        if (verbosity>2) cout<<"draw1 "<<flush;
        histo1->Draw("same");
        //		histo1->GetYaxis()->SetRangeUser(min,max);
    }
    c1->Update();
    TVirtualPad *pad =c1->GetPad(0);
    if (verbosity>2) cout<<"MIN: "<<min<<"-->";
    min=(double)(min/refactorSecond);
    if (verbosity>2) cout<<min<<"\t\tMAX: "<<max<<"--->";
    max = (double)(max/refactorSecond);
    if (verbosity>2) cout<<max<<endl;
    TGaxis *axis = new TGaxis(pad->GetUxmax(),pad->GetUymin(),pad->GetUxmax(), pad->GetUymax(),min,max,510,"+L");
    axis->SetLineColor(histo2->GetLineColor());
    axis->SetLabelColor(histo2->GetLineColor());
    axis->SetTextColor(histo2->GetLineColor());
    axis->SetTitle(histo2->GetYaxis()->GetTitle());
    axis->Draw("same");
    c1->Update();
    TLegend *leg;

    position.ToLower();
    if (position== "right")
        leg = new TLegend(0.52,0.75,0.9,0.95);
    else
        leg =new TLegend(0.1,0.75,0.48,0.9);
    leg->SetFillColor(kWhite);
    leg->SetHeader("Legend");
    leg->AddEntry(histo1,histo1->GetName());
    leg->AddEntry(histo2,histo2->GetName());
    leg->Draw("same");
    if(settings->IsPaperMode()){

    }
    else{
        TPaveText* pt2 = 0;
        if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvasName.Data()));
        if(pt2 && !settings->IsPaperMode()) pt2->Draw("same");
    }
    c1->Update();
    SaveCanvas(c1);
}


void HistogrammSaver::SaveStringToFile(string name, string data)
{
    std::ofstream file;
    stringstream outputFileName;
    cout<<"FILE: "<<name<<endl;
    //cout<<"PATH: "<<sys->pwd()<<endl;
    outputFileName<<this->GetPlotsPath()<<"/"<<name;
    cout<<"create String to file: \""<<outputFileName.str()<<"\""<<endl;
    file.open(outputFileName.str().c_str());
    file<<data;
    file.close();
}

/**
 *
 * @param histo
 * @param minX
 * @param maxX
 * @return
 * @todo add sigma value to pt
 */
TPaveText* HistogrammSaver::updateMean(TH1F* histo, Float_t minX, Float_t maxX) {
    Int_t minBin = histo->FindBin(minX);
    Int_t maxBin = histo->FindBin(maxX);
    Float_t mean = 0;
    Float_t sigma = 0;
    Float_t nEntries = 0;
    for (Int_t bin = minBin; bin<maxBin+1;bin++){
        Float_t weighted = histo->GetBinContent(bin)*histo->GetBinCenter(bin);
        nEntries+=histo->GetBinContent(bin);
        mean += weighted;
        sigma+= weighted*weighted;
    }
    mean = mean/nEntries;
    //	cout<<"new calculation of mean for range ["<<minX<<","<<maxX<<"]"<<endl;
    TCanvas *c1 = new TCanvas();
    histo->Draw();
    c1->Update();
    Float_t maxY = histo->GetBinContent(histo->GetMaximumBin());
    Float_t maxYPos = histo->GetBinCenter(histo->GetMaximumBin());
    Int_t bin1 = histo->FindFirstBinAbove(maxY/2);
    Int_t bin2 = histo->FindLastBinAbove(maxY/2);
    Float_t fwhm = histo->GetBinCenter(bin2) - histo->GetBinCenter(bin1);
    TPaveStats* hstat = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
    if (verbosity)
        cout<<"Got stats" <<hstat<<endl;
    if (!hstat)
        histo->GetListOfFunctions()->Print();

    TF1* fit = 0;
    fit = (TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("Fitfcn_%s",histo->GetName()));
    if(!fit)
        fit = (TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("fLangauFixedNoise_%s",histo->GetName()));

    TPaveText* hstat2 = 0;
    if (hstat) hstat2 = (TPaveText*) hstat->Clone();
    else
        histo->GetListOfFunctions()->Print();
    if(hstat2){
        TText * text = hstat2->AddText("");
        text->SetTextSize(0);
        text = hstat2->AddText(TString::Format("Mean_{> %.1f}  =   %.1f",minX,mean));
        text->SetTextSize(0);
        text = hstat2->AddText(TString::Format("Mean_{all}  =   %.1f",histo->GetMean()));
        text->SetTextSize(0);
        text = hstat2->AddText("");
        text->SetTextSize(0);
        text = hstat2->AddText(TString::Format("FWHM  =   %.1f",fwhm));
        text->SetTextSize(0);
        text = hstat2->AddText(TString::Format("MP_{histo}  =   %.1f",maxYPos));
        text->SetTextSize(0);
        text = hstat2->AddText(TString::Format("FWHM/MP_{histo}  =   %.3f",fwhm/maxYPos));
        text->SetTextSize(0);
        text = hstat2->AddText("");
        text->SetTextSize(0);
        if(fit){
            Float_t width = fit->GetParameter(0);
            Float_t gsigma = fit->GetParameter(3);
            text = hstat2->AddText(TString::Format("Width/GSigma  =   %.3f",width/gsigma));
            text->SetTextSize(0);
            text = hstat2->AddText("");
            text->SetTextSize(0);
        }
        Float_t yNDC = 0.5;
        hstat2->SetY1NDC(yNDC);
    }
    else{
        cout<<"something is bad..."<<endl;
    }
    return hstat2;

}


TPaveText* HistogrammSaver::GetUpdatedLandauMeans(TH1F* histo,Float_t mpv,Float_t gSigma){
    if (!histo)
        return 0;
    Float_t minX,maxX;
    minX = (-1.) *   std::numeric_limits<float>::infinity();
    maxX =  std::numeric_limits<float>::infinity();
    //Find good mean calculation Area
    Int_t startBin = histo->FindBin(mpv);
    //	cout<<"Start Bin: " <<startBin<<endl;
    Float_t max = histo->GetBinContent(startBin);
    Int_t bin;
    for(bin = startBin;bin>0;bin--){
        if(histo->GetBinContent(bin)<.05*max)
            break;
    }
    Int_t deltaBins = startBin - bin;
    bin = startBin-deltaBins*1.5;
    if(bin>0)
        minX = histo->GetBinLowEdge(bin);
    else
        minX = mpv*.5;
    //Add a "fit" to histo
    TPaveText* ptMean = updateMean(histo,minX,maxX);
    if(gSigma>0&&ptMean)
        ptMean->AddText(TString::Format("GSigma  =   %.1f",gSigma));
    maxX = histo->GetBinLowEdge(histo->GetNbinsX());
    TF1* fMeanCalculationArea = new TF1("fMeanCalculationArea","pol0",minX,maxX);
    fMeanCalculationArea->SetLineColor(kGreen);
    fMeanCalculationArea->FixParameter(0,0);
    fMeanCalculationArea->SetLineWidth(5);
    histo->Fit(fMeanCalculationArea,"Q+","",minX,maxX);
    return ptMean;
}

void HistogrammSaver::AddMarkedCells(TCanvas *c1){
    if (!c1) return;
    c1->cd();
    vector<TCutG*> goodCells;
    TCutG *cell;
    for (UInt_t i = 0; i < settings->GetNCells3d();i++){
        int cellType = 0;
        if (settings->IsGoodCell(3,i))
            cellType += 1;
        if (settings->isBadCell(3,i))
            cellType += 2;
        if (settings->isDeadCell(3,i))
            cellType += 4;
        if (cellType==0)
            continue;
        int column = settings->getColumnOfCell(i);
        int row = settings->getRowOfCell(i);
        int cellno = settings->get3DCellNo(column,row);
        //cout<<"Cell: "<<i<<"|"<<cellno<<": "<<column<<"/"<<row<<": "<<cellType<<endl;
        cell = settings->GetCell(i,c1->GetName());
        cell->SetLineWidth(2);
        switch (cellType){
            case 1: cell->SetLineColor(kGreen); break;
            case 2: cell->SetLineColor(kRed); break;
            case 4: cell->SetLineColor(kPink); break;
            case 6: cell->SetLineColor(kBlue); break;
            default: cell->SetLineWidth(4);break;
        }
        cell->Draw("same");
        if (cellType == 1)
            goodCells.push_back(cell);
    }
    for (UInt_t i =0; i<goodCells.size();i++)
        goodCells.at(i)->Draw("same");
}
TCanvas* HistogrammSaver::DrawHistogramWithCellGridAndMarkedCells(TH2* histo,TH2* histo2){
    TCanvas * c1 = DrawHistogramWithCellGrid(histo,histo2);
    AddMarkedCells(c1);
    return c1;
}
TCanvas* HistogrammSaver::DrawHistogramWithCellGrid(TH2* histo,TH2* histo2){
    if (verbosity) cout<<"HistogrammSaver::DrawHistogramWithCellGrid: "<<flush;
    if (verbosity) cout<<"H1: ";
    if (histo)
        if (verbosity) cout<<histo->GetName();
    else
        if (verbosity) cout<<histo;
    if (verbosity) cout<<"\tH2: ";
    if (histo2)
        if (verbosity) cout<<histo2->GetName();
    else
        if (verbosity) cout<<histo2;
    if (verbosity) cout<<flush;

    TString name = histo->GetName();
    if (name.BeginsWith("h"))
        name.Replace(0,1,"c");
    else
        name.Insert(0,"c_");
    TCanvas* c1 = new TCanvas(name,name);
    c1->cd();
    c1->SetRightMargin(.15);
    hGridReferenceDetSpace->SetTitle(histo->GetTitle());		//Set title to require
    hGridReferenceDetSpace->Draw("COL");
    histo->GetZaxis()->SetTitleOffset(1.3);
    histo->GetZaxis()->SetLabelOffset(0);
    if (histo){
        histo->SetContour(100);
        histo->Draw("sameCOLZ");
        hGridReferenceDetSpace->Draw("sameCOL");
    }
    if (histo2){
        histo2->Draw("sameTEXTAH");
    }
    settings->DrawMetallisationGrid(c1, 3);
    name.Replace(0,1,"f");
    TFile *f = GetFilePointer(GetROOTFileName(name),"RECREATE");
    f->cd();
    histo->Clone()->Write();
    hGridReferenceDetSpace->Clone()->Write();
    if(histo2)
        histo2->Clone()->Write();
    delete f;
    if (verbosity) cout<<"-Done"<<endl;
    return c1;
}

void HistogrammSaver::SaveHistogramWithCellGridAndMarkedCells(TH2* histo,TH2* histo2){
    if (!histo)
        return;
    gStyle->SetPaintTextFormat("4.2f");
    TCanvas *c1 = DrawHistogramWithCellGridAndMarkedCells(histo,histo2);
    this->SaveCanvas(c1);
}

void HistogrammSaver::SaveHistogramWithCellGrid(TH2* histo,TH2* histo2) {
    //	cout<<"[HistogrammSaver::SaveHistogramWithCellGrid]\t"<<flush;
    if (!histo)
        return;
    gStyle->SetPaintTextFormat("4.2f");
    TCanvas *c1 = DrawHistogramWithCellGrid(histo,histo2);
    this->SaveCanvas(c1);
}

void HistogrammSaver::DrawFailedQuarters(
        vector<pair<Int_t, Int_t> > failedQuarters, TCanvas* c1) {
    if (!c1){
        cerr<<"DrawFailedQuarters - Canvas =0"<<endl;
        return;
    }
    if(verbosity)cout<<"[HistogrammSaver::DrawFailedQuarters]"<<endl;
    UInt_t DiamondPattern = 3;
    Float_t xStart = settings->get3dMetallisationFidCuts()->getXLow(DiamondPattern);
    Float_t yStart =settings->get3dMetallisationFidCuts()->getYLow(DiamondPattern);
    UInt_t det = TPlaneProperties::getDetDiamond();
    Float_t cellwidth = settings->GetCellWidth(det,DiamondPattern-1);
    Float_t cellheight = settings->GetCellHeight();
    TCutG * failedQuarter;
    int i =0;
    for (vector<pair<Int_t, Int_t> >::iterator quarter = failedQuarters.begin();
            quarter != failedQuarters.end(); ++quarter){
        i++;
        Int_t ncell = (*quarter).first;
        Int_t nquarter = (*quarter).second;
        if(!settings->isValidCellNo((UInt_t)ncell)){
            cerr<< "Invalid Cell No: " << ncell<<endl;
            continue;
        }
        int column = settings->getColumnOfCell(ncell);
        int row = settings->getRowOfCell(ncell);
        float xLow = xStart + (column+.5*(nquarter/2))*cellwidth;
        float yLow = yStart + (row+.5*(nquarter%2))*cellheight;
        float xHigh = xLow+cellwidth/2;
        float yHigh = yLow+cellheight/2;
        if(verbosity)cout<<ncell<<"."<<nquarter;
        if(verbosity)cout<<": X:" <<xLow<<"-"<<xHigh;
        if(verbosity)cout<<"; Y:" <<yLow<<"-"<<yHigh<<"; "<<flush;
        TString name = c1->GetName();
        name.Append(TString::Format("_FailedQuarter_%dOf%d",i,(int)failedQuarters.size()));
        if(verbosity)cout<<" DRAW: "<< name<<flush;;
        failedQuarter = new TCutG(name,5);
        failedQuarter->SetPoint(0,xLow,yLow);
        failedQuarter->SetPoint(1,xLow,yHigh);
        failedQuarter->SetPoint(2,xHigh,yHigh);
        failedQuarter->SetPoint(3,xHigh,yLow);
        failedQuarter->SetPoint(4,xLow,yLow);
        failedQuarter->SetFillStyle(3001);
        failedQuarter->SetLineWidth(0);
        failedQuarter->SetFillColor(kRed);
        failedQuarter->Draw("sameF");
        if(verbosity)cout<<" - Done"<<endl;
    }
}

TH2D* HistogrammSaver::GetHistoBinedInQuarters(TString name) {
    return GetHistoBinedInCells(name,2);
}

TProfile2D* HistogrammSaver::GetProfileBinnedAroundFieldWires(TString name, Int_t binsPerCellAxis) {
    cout<<"[HistogrammSaver::GetHistoBinnedAroundFieldWires] create "<<name<<endl;
    TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
    Int_t xbins = settings->getNColumns3d();
    Int_t ybins = settings->getNRows3d();
    Float_t xlow = diaMetFidCut->GetXLow();
    Float_t xup = diaMetFidCut->GetXHigh();
    Float_t dx = (xup-xlow)/xbins;
    Float_t ylow = diaMetFidCut->GetYLow();
    Float_t yup = diaMetFidCut->GetYHigh();
    Float_t dy = (yup-ylow)/ybins;
    xlow -=dx/2.;
    xup  +=dx/2.;
    ylow -=dy/2.;
    yup  +=dy/2.;
    xbins = (xbins+1)*binsPerCellAxis;
    ybins = (ybins+1)*binsPerCellAxis;
    TProfile2D* histo = new TProfile2D(name,name,xbins,xlow,xup,ybins,ylow,yup);
    return histo;
}


TH2D* HistogrammSaver::GetHistoBinnedAroundFieldWires(TString name, Int_t binsPerCellAxis) {
    cout<<"[HistogrammSaver::GetHistoBinnedAroundFieldWires] create "<<name<<endl;
    TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
    Int_t xbins = settings->getNColumns3d();
    Int_t ybins = settings->getNRows3d();
    Float_t xlow = diaMetFidCut->GetXLow();
    Float_t xup = diaMetFidCut->GetXHigh();
    Float_t dx = (xup-xlow)/xbins;
    Float_t ylow = diaMetFidCut->GetYLow();
    Float_t yup = diaMetFidCut->GetYHigh();
    Float_t dy = (yup-ylow)/ybins;
    xlow -=dx/2.;
    xup  +=dx/2.;
    ylow -=dy/2.;
    yup  +=dy/2.;
    xbins = (xbins+1)*binsPerCellAxis;
    ybins = (ybins+1)*binsPerCellAxis;
    TH2D* histo = new TH2D(name,name,xbins,xlow,xup,ybins,ylow,yup);
    return histo;
}

TH2D* HistogrammSaver::GetHistoBinedInCells(TString name, Int_t binsPerCellAxis) {
    if (verbosity) cout<<"create "<<name<<endl;
    TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
    TH2D* histo = new TH2D(name,name,
            settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
            settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh());
    return histo;
}

TH3D* HistogrammSaver::Get3dHistoBinedInCells(TString name, UInt_t binsz,
        Float_t minz, Float_t maxz, Int_t binsPerCellAxis) {
    if (verbosity) cout<<"create "<<name<<endl;
    TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
    TH3D* histo = new TH3D(name,name,
            settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
            settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh(),
            binsz, minz,maxz);
    return histo;
}

TProfile2D* HistogrammSaver::GetProfile2dBinedInCells(TString name,
        Int_t binsPerCellAxis) {
    if (verbosity) cout<<"create "<<name<<endl;
    TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
    if (verbosity) cout<<" HistogrammSaver::GetProfile2dBinedInCells: "<<name<<"  binsPerCellAxis: "<<binsPerCellAxis<<endl;
    diaMetFidCut->Print();
    if (verbosity) cout<<"Columns: "<<settings->getNColumns3d()<<"/ Rows:"<<settings->getNRows3d()<<endl;
    TProfile2D* histo = new TProfile2D(name,name,
            settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
            settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh());
    return histo;
}

TProfile2D* HistogrammSaver::CreateProfile2D(TString name,
        std::vector<Float_t> posX, std::vector<Float_t> posY,
        std::vector<Float_t> posZ, UInt_t nBinsX, UInt_t nBinsY,
        Float_t minRangeX, Float_t maxRangeX, Float_t minRangeY,
        Float_t maxRangeY, Float_t minRangeZ, Float_t maxRangeZ,
        Float_t factor) {
    if(posX.size()!=posY.size()||posX.size()==0) {
        cerr<<"ERROR HistogrammSaver::CreateScatterHisto vectors have different size "<<posX.size()<<" "<<posY.size()<<" "<<name<<endl;
        return new TProfile2D(name,name,2,0,1,2,0,1);
    }
    cout<<"CreateProfile2D "<<name<<endl;
    cout<<TString::Format("maxRange:\nX: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minRangeX,maxRangeX,minRangeY,maxRangeY,minRangeZ,maxRangeZ)<<endl;
    Float_t maxX = posX.at(0);
    Float_t maxY = posY.at(0);
    Float_t maxZ = posZ.at(0);
    Float_t minX = posX.at(0);
    Float_t minY = posY.at(0);
    Float_t minZ = posZ.at(0);
    //  cout<<" Create Histo: '"<<name<<"' - Range ("<<minRangeX<<"-"<<maxRangeX<<"),  ("
    //          <<minRangeY<<"-"<<maxRangeY<<"), ("<<minRangeZ<<"-"<<maxRangeZ<<")"<<endl;
    for(UInt_t i=0;i<posX.size();i++){
        if (posX.at(i)<minRangeX||posX.at(i)>maxRangeX)
            continue;
        if (posY.at(i)<minRangeY||posY.at(i)>maxRangeY)
            continue;
        if (posZ.at(i)<minRangeZ||posZ.at(i)>maxRangeZ)
            continue;
        if(posX.at(i)>maxX)maxX=posX.at(i);
        else if(posX.at(i)<minX)minX=posX.at(i);
        if(posY.at(i)>maxY)maxY=posY.at(i);
        else if(posY.at(i)<minY)minY=posY.at(i);
        if(posZ.at(i)>maxZ)maxZ=posZ.at(i);
        else if(posZ.at(i)<minZ)minZ=posZ.at(i);
    }
    //  cout<<TString::Format("X: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minX,maxX,minY,maxY,minZ,maxZ)<<endl;
    Float_t factorX = factor;
    Float_t factorY = factor;
    Float_t factorZ = factor;
    Float_t deltaXMax = maxRangeX - minRangeX;
    Float_t deltaYMax = maxRangeY - minRangeY;
    Float_t deltaZMax = maxRangeZ - minRangeZ;
    Float_t maxDiff = 0.02;
    if ( TMath::Abs(maxRangeX-maxX)/deltaXMax <= maxDiff && TMath::Abs(minRangeX - minX)/deltaXMax <= maxDiff ) {
        factorX = 0;
        maxX = maxRangeX;
        minX = minRangeX;
    }
    if ( TMath::Abs(maxRangeY-maxY)/deltaYMax <= maxDiff && TMath::Abs(minRangeY - minY)/deltaYMax <= maxDiff ) {
        factorY = 0;
        minY = minRangeY;
        maxY = maxRangeY;
    }
    if ( TMath::Abs(maxRangeZ-maxZ)/deltaZMax <= maxDiff && TMath::Abs(minRangeZ - minZ)/deltaZMax <= maxDiff ) {
        factorZ = 0;
        minZ = minRangeZ;
        maxZ = maxRangeZ;
    }
    Float_t deltaX=maxX-minX;
    Float_t deltaY=maxY-minY;
    Float_t deltaZ=maxZ-minZ;
    //  cout<<"\t"<<deltaX<<" "<<deltaY<<" "<<deltaZ<<endl;
    //  cout<<"\t"<<factorX<<" "<<factorY<<" "<<factorZ<<endl;
    minX = minX-factorX*deltaX;
    maxX = maxX+factorX*deltaX;
    minY = minY-factorY*deltaY;
    maxY = maxY+factorY*deltaY;
    minZ = minZ-factorZ*deltaZ;
    maxZ = maxZ+factorZ*deltaZ;
    cout<<TString::Format("X: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minX,maxX,minY,maxY,minZ,maxZ)<<endl;
    //  char t; cin>>t;
    TProfile2D* histo = new TProfile2D(name,name,
            nBinsX,minX,maxX,
            nBinsY,minY,maxY,
            minZ,maxZ);
    for(UInt_t i=0;i<posX.size();i++){
        if (posX.at(i)<minRangeX||posX.at(i)>maxRangeX)
            continue;
        if (posY.at(i)<minRangeY||posY.at(i)>maxRangeY)
            continue;
        if (posZ.at(i)<minRangeZ||posZ.at(i)>maxRangeZ)
            continue;
        histo->Fill(posX.at(i),posY.at(i),posZ.at(i));
    }
    return histo;
}

TString HistogrammSaver::GetFileName(TString name, TString extension){
    if (name.Contains("/"))
        name = name(name.Last('/')+1,name.Length());
    TString path = (TString)plots_path;
    if (!path.EndsWith("/"))
        path += "/";
    if (name.EndsWith(".root"))
        name.ReplaceAll(".root","");
    if (name.EndsWith(".png"))
        name.ReplaceAll(".eps","");
    path+=extension+"/";
    if (name == "")
        name="histograms";
    path+=(TString)name;
    TString appendix = TString::Format(".%d.",runNumber);
    appendix+=extension;
//    if (!path.EndsWith(appendix))
    path.Append(appendix);
//    path.ReplaceAll("//","/");
    return path;
}

TFile* HistogrammSaver::GetFilePointer(TString name, TString option) {
    TFile* file = new TFile(GetROOTFileName(name),option);
    if (!file and option.Contains("UPDATE"))
        file = new TFile(GetROOTFileName(name),"RECREATE");
    return file;
}

void HistogrammSaver::DrawGoodCellsRegion(TCanvas* c1) {
    if (!c1)
        return;
    cout<<"HistogrammSaver::DrawGoodCellsRegion"<<endl;
    std::vector< std::vector<Int_t> > goodCellRegions = settings->getGoodCellRegions3d();
    Float_t xmin = 1e9;
    Float_t xmax = -1e9;
    Float_t ymin = 1e9;
    Float_t ymax = -1e9;
    for (UInt_t region = 0; region < goodCellRegions.size();region++){
        for (UInt_t cell =0; cell< goodCellRegions.at(region).size();cell++){
            Int_t cellNo = goodCellRegions.at(region).at(cell);
            Int_t row = settings->getRowOfCell(cellNo);
            Int_t column = settings->getColumnOfCell(cellNo);
            TString name = TString::Format("gGoodCell_%d_%d",region,cellNo);
            Int_t oldVerbose = verbosity;
            verbosity=10;
            cout<<region<<" "<<cell<<" "<<cellNo<<": "<<row<<" "<<column<<endl;
            Float_t xx1[]={ hGridReferenceDetSpace->GetXaxis()->GetBinLowEdge(column+1),hGridReferenceDetSpace->GetXaxis()->GetBinUpEdge(column+1)};
            Float_t yy1[]={ hGridReferenceDetSpace->GetYaxis()->GetBinLowEdge(row+1),hGridReferenceDetSpace->GetYaxis()->GetBinUpEdge(row+1)};
            if (xx1[0]<xmin) xmin=xx1[0];
            if (xx1[1]>xmax) xmax=xx1[1];
            if (yy1[0]<ymin) ymin=yy1[0];
            if (yy1[1]>ymax) ymax=yy1[1];
//            Float_t x = hGridReferenceDetSpace->GetXaxis()->GetBinCenter(column+1);
//            Float_t y = hGridReferenceDetSpace->GetYaxis()->GetBinCenter(row+1);
//            TCutG* gCell = this->GetCutGofBin(name,hGridReferenceDetSpace,x,y);
//            verbosity = oldVerbose;
//            gCell->SetLineColor(kCyan+2);
//            gCell->SetLineWidth(4);
//            for (UInt_t i =0; i< gCell->GetN();i++){
//                cout<<"\t"<<gCell->GetX()[i]<<" / "<<gCell->GetY()[i]<<endl;
//                if(gCell->GetX()[i]<xmin) xmin = gCell->GetX()[i];
//                if(gCell->GetX()[i]>xmax) xmax = gCell->GetX()[i];
//                if(gCell->GetY()[i]<ymin) ymin = gCell->GetY()[i];
//                if(gCell->GetY()[i]>ymax) ymax = gCell->GetY()[i];
//
//            }
////            gCell->Draw("same");
//            delete gCell;
        }
//        ymax = ymax-.03*(ymax-ymin);
        Float_t xx[] = {xmin,xmax,xmax,xmin,xmin};
        Float_t yy[] = {ymin,ymin,ymax,ymax,ymin};
        TString name = TString::Format("gGoodCellRegion_%d",region);
        TCutG * cut = new TCutG(name,5,xx,yy);
        cut->SetLineColor(kBlack);
        cut->SetLineWidth(7);
        cut->Draw("same");
    }
}

TH1F* HistogrammSaver::GetBinsInHistogram(TH2* histo2, UInt_t nbins) {
    if (!histo2)
        return 0;
    Float_t min = 1e9;
    Float_t max = -1e9;
    Float_t content;
    for (UInt_t xbin = 1; xbin<= histo2->GetNbinsX();xbin++)
    for (UInt_t ybin = 1; ybin<= histo2->GetNbinsY();xbin++){
        content = histo2->GetBinContent(xbin,ybin);
        if (content<min)
            min=content;
        if (content>max)
            max = content;
    }
    Float_t delta = max-min;
    min -= .1*delta;
    max += .1*delta;
    TString name = histo2->GetName()+(TString)"_binned";
    TString title = histo2->GetTitle();
    TH1F* histo = new TH1F(name,title,nbins,min,max);
    histo->GetXaxis()->SetTitle(histo2->GetZaxis()->GetTitle());
    for (UInt_t xbin = 1; xbin<= histo2->GetNbinsX();xbin++)
       for (UInt_t ybin = 1; ybin<= histo2->GetNbinsY();xbin++){
           histo->Fill(histo2->GetBinContent(xbin,ybin));
       }
    return histo;
}

TH1F* HistogrammSaver::SaveProjectionZ(TH2* histo, bool isOverlay,bool doFit,Float_t minX,Float_t maxX,Int_t nbins) {
    if (!histo)
        return 0;
    if (minX>maxX){
        minX =  histo->GetBinContent(histo->GetMinimumBin());
        maxX = histo->GetBinContent(histo->GetMaximumBin());
        Float_t deltaX = maxX-minX;
        maxX += deltaX*0.2;
        minX -= deltaX*0.1;
    }
//    cout<<"HistogrammSaver::SaveOverlayDistribution"<<histo<<" range: "<<minX<<"-"<<maxX<<flush;
    TString name = histo->GetName()+TString("_1D");
    TH1F *histo_1D = new TH1F(name,"all other bins",nbins,minX,maxX);
    name = histo->GetName()+TString("_1D_middle");
    TH1F *histo_1D_middle = new TH1F("hOverlay_1D_middle","readout column bin",nbins,minX,maxX);
    name = histo->GetName()+TString("_1D_corners");
    TH1F *histo_1D_corner = new TH1F("hOverlay_1D_corners","bias column bins",nbins,minX,maxX);
    Float_t content,xlow,ylow,xhigh,yhigh;
//    cout<<" fill content"<<flush;
    for (UInt_t bx=1;bx<= histo->GetNbinsX();bx++){
        for (UInt_t by = 1; by <= histo->GetNbinsY();by++){
            content = histo->GetBinContent(bx,by);
            xlow = histo->GetXaxis()->GetBinLowEdge(bx);
            xhigh = histo->GetXaxis()->GetBinUpEdge(bx);
            ylow = histo->GetYaxis()->GetBinLowEdge(by);
            yhigh = histo->GetYaxis()->GetBinUpEdge(by);
            if ( isOverlay && (xlow < 75 && 75 < xhigh ) && ( ylow < 75 && 75 < yhigh))
                histo_1D_middle->Fill(content);
            else if ( isOverlay && ( xlow <=0 && 0 <= xhigh) && (ylow<=0 && 0<=yhigh))
                histo_1D_corner->Fill(content);
            else if( isOverlay && ( xlow<=150 && 150 <=xhigh) && ( ylow<=0 && 0<=yhigh))
                histo_1D_corner->Fill(content);
            else if( isOverlay && ( xlow<=150 && 150 <=xhigh) && (ylow<=150 && 150<=yhigh))
                histo_1D_corner->Fill(content);
            else if( isOverlay && ( xlow<=0 && 0<=xhigh) && (ylow<=150 && 150<=yhigh))
                histo_1D_corner->Fill(content);
            else
                histo_1D->Fill(content);
        }
    }
//    cout<<"-done "<<histo_1D->GetEntries()<<"/"<<histo_1D_corner->GetEntries()<<histo_1D_middle->GetEntries()<<flush;
    histo_1D->SetLineColor(kBlack);
//    cout<<" fitGaus"<<flush;
    TF1* fGaus = new TF1("fgaus","gaus",minX,maxX);
    fGaus->SetLineStyle(2);
    fGaus->SetLineColor(kGreen);
    if (doFit)
        histo_1D->Fit(fGaus,"Q");

    histo_1D_middle->SetLineColor(kRed);
    histo_1D_middle->SetLineStyle(2);
    histo_1D_middle->SetFillColor(kRed);
    histo_1D_middle->SetFillStyle(3013);
    histo_1D_corner->SetLineColor(kBlue);
    histo_1D_corner->SetLineStyle(3);
    histo_1D_corner->SetFillColor(kBlue);
    histo_1D_corner->SetFillStyle(3021);
    name= TString("stack_")+histo->GetName()+TString("_1D");
//    cout<<" saveStack"<<flush;
    TString title;
    if (isOverlay)
        title = TString("Overlay distribution;");
    else
        title = histo->GetTitle()+(TString)" pz;";
    title+=histo->GetZaxis()->GetTitle()+TString(";number of entries");
    THStack* hstack = new THStack(name,title);
    hstack->Add(histo_1D);
    if (isOverlay){
        hstack->Add(histo_1D_middle);
        hstack->Add(histo_1D_corner);
        this->SaveStack(hstack,"",true);
    }
    else
        this->SaveStack(hstack,"",false);
    delete hstack;
    delete histo_1D_middle;
    delete histo_1D_corner;
    return histo_1D;
//    cout<<" DONE"<<endl;
}

void HistogrammSaver::SaveNegativeChargeOverlay(TH2* prof) {
    if (!prof)
        return;
    TString name = prof->GetName();
    prof->SetName("cc_"+name);
    Float_t max = prof->GetBinContent(prof->GetMaximumBin());
    Float_t min = prof->GetBinContent(prof->GetMinimumBin());
    Float_t range = TMath::Max(TMath::Abs(min),TMath::Abs(max));
    int  NCont = 999;
    UInt_t NRGBs = 3;
    gStyle->SetNumberContours(NCont);
    Double_t stops[] = { 0.00, 0.50, 1.00};
    Double_t red[]   = { 0.00, 1.00, 1.0};
    Double_t green[] = { 0.00, 1.00, 0.00};
    Double_t blue[]  = { 1.00, 1.00, 0.00};

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    prof->SetContour(999);
    min = prof->GetMinimum();
    max = prof->GetMaximum();
    prof->SetMinimum(range*-1.05);
    prof->SetMinimum(range*1.05);
    prof->GetZaxis()->SetRangeUser(-1.05*range,1.05*range);
    this->SaveOverlay(prof);
    prof->SetName("ch_"+name);
    this->SaveHistogram(prof,true,false);
    gStyle->SetPalette(53); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
    prof->SetName(name);
    prof->SetMinimum(min);
    prof->SetMaximum(max);
}

void HistogrammSaver::SaveProjectionX(TH2* histo) {
    if (!histo)
        return;
    TString name = histo->GetName();
    name += "_px";
    TH1D *px = histo->ProjectionX(name);
    this->SaveHistogram(px);
    delete px;
}

void HistogrammSaver::SaveBinnedProjectionX(TH2* histo) {
    if (!histo)
        return;
    TString name = histo->GetName();
    name += "_BinnedPx";
    cout<<"[HistogrammSaver::SaveBinnedProjectionX]: "<<name<<endl;
    Int_t nBinsY = histo->GetNbinsY();
    TString yTitle = histo->GetYaxis()->GetTitle();
    vector<TH1D*> projections;
    TString title = name;
    title+=";";
    title+=histo->GetXaxis()->GetTitle();
    title+=";number of entries";
    THStack *stack = new THStack(name,name);
    Int_t colors[] = {kBlack,kGreen,kBlue,kOrange};
    Int_t styles[] = {3,1,2,4};
    Int_t nStyles = 4;
    Int_t n=0;
    cout<<"\tCreate stack out of "<<nBinsY<<"Projections."<<endl;
    for (UInt_t i = 1; i <= nBinsY;i++){
        TString hname = name;
        hname+=TString::Format("_bin%d",i);
        TH1D *px = histo->ProjectionX(hname,i,i);
        px->SetTitle(yTitle+TString::Format(", Bin %d",i));
        cout<<"\t  * "<<i<<": "<<hname<<" "<<px->GetEntries()<<endl;
        if (px->GetEntries()){
            px->SetLineStyle(styles[n%nStyles]);
            px->SetLineColor(styles[n%nStyles]);
            projections.push_back(px);
            stack->Add(px);
            n+=1;
        }
        else
            delete px;
    }
    cout<<"\tSave Stack \""<<stack->GetName()<<"\" with  "<<projections.size()<<" Projections."<<endl;
    this->SaveStack(stack,"nostack",true,false);
    delete stack;
    for (UInt_t i = 0; i < projections.size();i++)
        delete projections[i];
    cout<<"[HistogrammSaver::SaveBinnedProjectionX] Done."<<endl;
}

void HistogrammSaver::SaveProjectionY(TH2* histo) {
    if (!histo)
        return;
    TString name = histo->GetName();
    name += "_py";
    TH1D *py = histo->ProjectionY(name);
    this->SaveHistogram(py);
    delete py;
}

void HistogrammSaver::SaveIntegral(TH1* histo,bool bRelative) {
    if (!histo) return;
    TString name = histo->GetName();
    name += "_Integral";
    if (bRelative)
        name+="Rel";
    Float_t factor = 1;
    Int_t xbins = histo->GetNbinsX();
    if (bRelative)
        factor = histo->Integral(0,xbins+1);
    if (factor==0)
        factor =1;
    Float_t xlow = histo->GetXaxis()->GetXmin();
    Float_t xup = histo->GetXaxis()->GetXmax();
    TH1D* hIntegral = new TH1D(name,name,xbins,xlow,xup);
    for (Int_t i = 1; i <= xbins; i++)
        hIntegral->SetBinContent(i,histo->Integral(0,i)/factor);
    hIntegral->SetEntries(histo->GetEntries());
    hIntegral->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    hIntegral->GetYaxis()->SetTitle("integrated no of events up to thrs.");
    if (bRelative){
        hIntegral->SetMinimum(0);
        hIntegral->SetMaximum(1.05);
    }
    this->SaveHistogram(hIntegral,false,false,true);
    delete hIntegral;
}

void HistogrammSaver::UpdatePaveText(){
    if (!pt)
        return;
    if (settings->IsPaperMode()){
        delete pt;
        pt =0;
        return;
    }
        ;
    pt->Clear();
    if (settings->IsPaperMode())
        return;
    pt->SetTextSize(0.0250);
    std::ostringstream svnRev_label;
    svnRev_label <<"Rev: "<<GIT_REV;
    pt->AddText(svnRev_label.str().c_str());
    std::ostringstream run_number_label;
    run_number_label << "Run " <<runNumber;
    pt->AddText(run_number_label.str().c_str());
    std::ostringstream pthresh2;
    pthresh2 <<"with "<< nEvents << " Events";
    pt->AddText(pthresh2.str().c_str());
    pt->AddText(dateandtime.AsSQLString());
    pt->SetBorderSize(0); //Set Border to Zero
    pt->SetFillColor(0); //Set Fill to White
}

void HistogrammSaver::SetRunNumber(unsigned int newRunNumber){
    runNumber=newRunNumber;
    if(verbosity)cout<<"HistogrammSaver: Set RunNumber="<<runNumber<<endl;
    UpdatePaveText();
}

void HistogrammSaver::SetNumberOfEvents(unsigned int nNewEvents){
    nEvents=nNewEvents;
    if(verbosity)cout<<"HistogrammSaver: Set Number of Events ="<<nEvents<<endl;
    UpdatePaveText();
}
void HistogrammSaver::SetPlotsPath(string path){
    if (plots_path.size()>1)
        MoveRootFilesInSubfolder();
    plots_path.assign(path);
    if(verbosity)cout<<"HistogrammSaver::Set Plotspath: \""<<plots_path<<"\""<<endl;
    int isNotCreated=sys->mkdir(plots_path.c_str(),true);
    if (isNotCreated!=0){
        //		cout<<"***************************************************\n";
        //		cout<<"********** Directory not created ******************\n";
        //		cout<<"***************************************************\n";
        if(verbosity)cout<<plots_path<<endl;
    }
    sys->mkdir(plots_path.c_str(),true);
    int stat = mkdir(plots_path.c_str(),0777);//0777(S_IRWXO||S_IRWXG||S_IRWXU));// S_IRWXU|S_IRGRP|S_IXGRP||S_IRWXU||S_IRWXG||S_IRWXO);
    if(!stat)cout<<"Verzeichnis angelegt: \""<<plots_path<<"\""<<endl;
    //	else cout<<"Verzeichnis konnte nicht angelegt werden..."<<endl;

    stringstream rootPath;
    rootPath<<plots_path<<"/root";
    sys->mkdir(rootPath.str().c_str(),true);
    stat = mkdir(rootPath.str().c_str(),0777);//0777(S_IRWXO||S_IRWXG||S_IRWXU));// S_IRWXU|S_IRGRP|S_IXGRP||S_IRWXU||S_IRWXG||S_IRWXO);
    if(!stat)cout<<"Verzeichnis angelegt: \""<<rootPath.str().c_str()<<"\""<<endl;
    stringstream epsPath;
    epsPath<<plots_path<<"/eps";
    sys->mkdir(epsPath.str().c_str(),true);
    stat = mkdir(epsPath.str().c_str(),0777);//0777(S_IRWXO||S_IRWXG||S_IRWXU));// S_IRWXU|S_IRGRP|S_IXGRP||S_IRWXU||S_IRWXG||S_IRWXO);
    if(!stat)cout<<"Verzeichnis angelegt: \""<<rootPath.str().c_str()<<"\""<<endl;

    //	else cout<<"Verzeichnis konnte nicht angelegt werden..."<<endl;

}

void HistogrammSaver::SetStyle(TStyle newStyle){
    //	currentStyle.TStyle(&newStyle);//.Clone());
    delete currentStyle;
    currentStyle = new TStyle(newStyle);
    currentStyle->cd();
    SetPaperPlotStyle();

}


void HistogrammSaver::SaveHistogramLandau(TH1F* histo){
    if(histo==0)return;
    cout<<"Save "<<histo->GetName()<<" "<<histo->GetEntries()<<endl;
    if(histo->GetEntries()==0)return;

    bool fixedNoise = false;
    TF1* fit = (TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("Fitfcn_%s",histo->GetName()));
    if(!fit){
        fit = (TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("fLangauFixedNoise_%s",histo->GetName()));
        fixedNoise = true;
    }
    Float_t mpv = histo->GetMean();
    Float_t gSigma = -1;
    if(fit){
        mpv = fit->GetParameter(1);
        if(fixedNoise)
            gSigma = fit->GetParameter(3);
    }
    //	cout<<"MPV: "<<mpv<<" mean: "<<histo->GetMean()<<" "<<fit->GetName()<<endl;
    TPaveText* stats = (TPaveText*) this->GetUpdatedLandauMeans(histo,mpv,gSigma);
    TString name = histo->GetName();
    if (name.First('h')==0)
        name.Replace(0,1,"c");
    else
        name = TString::Format("c_%s",histo->GetName());
    TCanvas *c1 = new TCanvas(name,name);
    c1->cd();
    histo->Draw();
    if(stats)stats->Draw();
    //	cout<<"Saving: "<<c1->GetName()<<endl;
    SaveCanvas(c1);
    //create ROOT
    SaveHistogramROOT(histo);
}

/**
 * *********************************************************
 * *********************************************************
 */
void HistogrammSaver::SaveHistogram(TH1* histo, bool fitGauss,bool adjustRange,bool drawStatBox,TString drawOption) {
    if(histo==0)return;
    if(histo->GetEntries()==0)return;
    if (!drawStatBox)
        histo->SetStats(false);
    if(adjustRange){
        int binxMin=0;
        for(binxMin=0;binxMin<histo->GetNbinsX();binxMin++)if(histo->GetBinContent(binxMin))break;
        int binxMax;
        for(binxMax=histo->GetNbinsX();binxMax>0;binxMax--)if(histo->GetBinContent(binxMax))break;
        histo->GetXaxis()->SetRangeUser(histo->GetBinLowEdge(binxMin),histo->GetBinLowEdge(binxMax+1));
    }
    //create PNG
    if (fitGauss) SaveHistogramFitGaussPNG(histo);
    else SaveHistogramPNG(histo,drawOption);
    //create ROOT
    SaveHistogramROOT(histo);
}

void HistogrammSaver::SaveHistogramWithExtendedFit(TH1* histo,TF1* fit, Float_t minX,Float_t maxX){
    if(histo==0)return;
    if(histo->GetEntries()==0)return;
    if (fit==0) return SaveHistogram(histo);
    TF1* fitExtended = (TF1*)fit->Clone();
    if(minX>maxX){
        Float_t temp = minX;
        minX = maxX;
        maxX = temp;
    }
    fitExtended->SetRange(minX,maxX);
    if(fitExtended->GetLineWidth()>1)
        fitExtended->SetLineWidth(fitExtended->GetLineWidth()-1);
    fitExtended->SetLineStyle(2);
    TString name = histo->GetName();
    name.Replace(0,1,'c');
    TCanvas* c1 = new TCanvas(name,name);
    c1->cd();
    TH1F *htemp = (TH1F*)histo->Clone();
    htemp->Draw();
    fitExtended->Draw("same");
    SaveCanvas(c1);
    delete c1;
}

void HistogrammSaver::SaveHistogramWithFit(TH1* histo,TF1* fit, UInt_t verbosity){
    return SaveHistogramWithFit(histo,fit,1,-1,verbosity);
    if(histo==0)return;
    if(histo->GetEntries()==0)return;
    if(fit==0) SaveHistogram(histo);
    if (verbosity>0) cout<<"Save Histogram With Fit:"<<histo->GetTitle()<<endl;
    TString name = histo->GetName();
    if (name.First('h') == 0)
        name.Replace(0,1,"c");
    else
        name = (TString)"c_"+name;
    TCanvas *plots_canvas =  new TCanvas(name,name );
    plots_canvas->Clear();
    plots_canvas->cd();
    TH1F *htemp = (TH1F*)histo->Clone();
    TF1* fittemp = (TF1*)fit->Clone();
    TPaveText * pt2 = 0;
    if(pt) pt2 =(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));

    htemp->Draw();
    fittemp->SetLineColor(kRed);
    fittemp->Draw("same");
    if(pt2 && ! settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plots_canvas->Print(GetROOTFileName(histo->GetName()));
    TString rootFileName = GetROOTFileName();
    TFile *f = GetFilePointer(rootFileName,"UPDATE");
    f->cd();
    TH1F* h_clone = (TH1F*)histo->Clone();
    h_clone->SetDirectory(f);
    h_clone->Write();
    TF1* f_clone = (TF1*)fit->Clone();
    f_clone->Write();
    plots_canvas->Write();
    SaveCanvas(plots_canvas);
    plot_filename << plots_path << histo->GetName() << ".png";
    plots_canvas->Print(plot_filename.str().c_str());
    plots_canvas->Print(GetEpsFileName(histo->GetName()));
    f->Close();
    //	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramWithFit(TH1* histo,TF1* fit,Float_t xmin,Float_t xmax, UInt_t verbosity){
    if(histo==0)return;
    if(histo->GetEntries()==0)return;
    if(fit==0) SaveHistogram(histo);
    if (verbosity>0) cout<<"Save Histogram With Fit:"<<histo->GetTitle()<<endl;
    TString name = histo->GetName();
    if (name.First('h') == 0)
        name.Replace(0,1,"c");
    else
        name = (TString)"c_"+name;
    TCanvas *plots_canvas =  new TCanvas(name,name );
    plots_canvas->Clear();
    plots_canvas->cd();
    //	TH1F *htemp = (TH1F*)histo->Clone();
    TPaveText * pt2 = 0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
    if (xmin<xmax){
        if(verbosity){
            cout<<"Fitting: "<<fit->GetName()<<" "<<xmin<<" - "<<xmax<<endl;
            histo->Fit(fit,"","",xmin,xmax);
        }
        else
            histo->Fit(fit,"Q+","",xmin,xmax);
    }
    else{//xmin>xmax ==> fitting wiout range
        if(verbosity){
            cout<<"Fitting: "<<fit->GetName()<<endl;
            histo->Fit(fit);
        }
        else
            histo->Fit(fit,"Q+");
    }
    histo->Draw();
    //	TF1* fittemp = (TF1*)fit->Clone();
    //	fittemp->SetLineColor(kRed);
    fit->SetLineStyle(3);
    fit->Draw("same");
    if(pt2 && !settings->IsPaperMode()) pt2->Draw();

    TString path = GetROOTFileName(histo->GetName());
    plots_canvas->Print(path);
    TString rootFileName = GetROOTFileName();
    TFile* f = this->GetFilePointer(rootFileName,"UPDATE");
    f->cd();
    TH1F* h_clone =(TH1F*)histo->Clone();
    h_clone->SetDirectory(f);
    h_clone->Write();
    ((TF1*)fit->Clone())->Write();
    plots_canvas->Write();
    path = (TString)plots_path + histo->GetName() +(TString)".png";
    plots_canvas->Print(path);
    plots_canvas->Print(GetEpsFileName(histo->GetName()));
    f->Close();
    delete f;
    //	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramWithCutLine(TH1F *histo,Float_t cutValue){
    if(histo==0)return;
    TCanvas *c2 = new TCanvas(TString::Format("c%s",histo->GetName()),histo->GetTitle());
    c2->cd();
    histo->Draw();
    double xCor[] = {cutValue,cutValue};
    double yCor[] = {0,histo->GetMaximum()*2};
    TGraph* line = new TGraph(2,xCor,yCor);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("Lsame");
    this->SaveCanvas(c2);
    delete c2;
}
void HistogrammSaver::SaveHistogramLogZ(TH2* histo){
    if(histo==0)return;
    TString canvasName = "c_";
    canvasName +=histo->GetName();
    TCanvas *c1 = new TCanvas(canvasName,canvasName);
    c1->cd();
    c1->SetLogz();
    TH2* htemp = (TH2*) histo->Clone();
    htemp->Draw("colz");
    this->SaveCanvas(c1);
    delete htemp;
    delete c1;
}

/**
 * The function Get ProfileX and GetProfileY for TProfile2D objects are only taken from class TH2D
 * This is not what one wants. This is an implementation based on
 *  http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=12103
 * To copy a bin content one needs to:
 * - sum of bin y values
 * - sum of bin y*y values
 * - sum of bin weight = number of bin entries if profile is filled with weights =1
 * - sum of bin weight ^2 if the profile is filled with weights different than 1
 *
 * @param prof2d
 * @param firstybin
 * @param lastybin
 * @return
 */
TProfile* HistogrammSaver::GetProfileX(TProfile2D* prof2d, TString name, Int_t firstybin, Int_t lastybin){
    if(!prof2d) return 0;
    if(lastybin==-1) lastybin = prof2d->GetYaxis()->GetNbins();
    cout<<"[HistogrammSaver::GetProfileX] "<<prof2d->GetName()<<" "<<firstybin<<"-"<<lastybin<<endl;
    if(name=="_pfx"||name =="")
        name = prof2d->GetName()+(TString)"_pfx";
    TString title = prof2d->GetTitle() + (TString)" profile X";
    UInt_t binsX = prof2d->GetXaxis()->GetNbins();
    Float_t xLow = prof2d->GetXaxis()->GetXmin();
    Float_t xHigh = prof2d->GetXaxis()->GetXmax();
    TProfile* prof = new TProfile(name,title,binsX,xLow,xHigh);
    cout<<"new Prof: "<<prof->GetName()<<endl;
    prof->Draw("goff");
    if(prof->GetXaxis())
        prof->GetXaxis()->SetTitle(prof2d->GetXaxis()->GetTitle());
    if (prof2d->GetZaxis()) title= prof2d->GetZaxis()->GetTitle();
    else title = "";
    if(firstybin!=1||lastybin!= prof2d->GetYaxis()->GetNbins())
        title.Append(TString::Format(" yrange: %f - %f",prof2d->GetYaxis()->GetBinLowEdge(firstybin),prof2d->GetYaxis()->GetBinLowEdge(lastybin)));
    if(prof->GetZaxis()) prof->GetYaxis()->SetTitle(title);
    cout<<"xaxis: "<<prof->GetXaxis()->GetTitle()<<endl;
    cout<<"yaxis: "<<prof->GetYaxis()->GetTitle()<<endl;
    Int_t firstxbin = 1;
    Int_t lastxbin = prof2d->GetXaxis()->GetNbins();
    for(Int_t biny = firstybin; biny <= lastybin; biny++){
        name = prof->GetName()+TString::Format("_biny%d:",biny);
        TProfile* profyBin = new TProfile(name,title,binsX,xLow,xHigh);
        cout<<"new ybin: "<<biny<<" "<<prof->GetName()<<": "<<endl;
        profyBin->Reset();
        Int_t nEntries =0;
        for(Int_t binx = firstxbin; binx <= lastxbin; binx++){
            Int_t bin2d = prof2d->GetBin(binx,biny);
            (*profyBin)[binx] = (*prof2d)[bin2d];    // copy bin y values
            (*profyBin->GetSumw2())[binx] =  (*prof2d->GetSumw2())[bin2d];   // copy bin y*y values
            profyBin->SetBinEntries(binx, prof2d->GetBinEntries(bin2d) );    // copy bin entries
            nEntries+=prof2d->GetBinEntries(bin2d);
            // copy (if needed) bin sum of weight square
            if ( prof2d->GetBinSumw2()->fN > bin2d ) {
                profyBin->Sumw2();
                (*profyBin->GetBinSumw2())[binx] = (*prof2d->GetBinSumw2())[bin2d];
            }
        }
        profyBin->SetEntries(nEntries);
        prof->Add(profyBin);
        cout<<""<<profyBin->GetEntries()<<" ---> "<<prof->GetEntries()<<endl;
        delete profyBin;
    }
    cout<<" final profile: "<<prof->GetEntries()<<"/"<<prof2d->GetEntries()<<endl;
    //    char t; cin>>t;
    return prof;
}

TProfile* HistogrammSaver::GetProfileY(TProfile2D* prof2d,TString name,Int_t firstxbin, Int_t lastxbin){
    if(!prof2d) return 0;
    if(lastxbin==-1) lastxbin = prof2d->GetXaxis()->GetNbins();
    cout<<"[HistogrammSaver::GetProfileY] "<<prof2d->GetName()<<" "<<firstxbin<<"-"<<lastxbin<<endl;
    if(name=="_pfy"||name =="")
        name = prof2d->GetName()+(TString)"_pfy";
    TString title = prof2d->GetTitle() + (TString)" profile Y";
    UInt_t binsX = prof2d->GetYaxis()->GetNbins();
    Float_t xLow = prof2d->GetYaxis()->GetXmin();
    Float_t xHigh = prof2d->GetYaxis()->GetXmax();
    TProfile* prof = new TProfile(name,title,binsX,xLow,xHigh);
    cout<<"new Prof: "<<prof->GetName()<<endl;
    prof->Draw("goff");
    if(prof->GetXaxis())
        prof->GetXaxis()->SetTitle(prof2d->GetYaxis()->GetTitle());
    if (prof2d->GetZaxis()) title= prof2d->GetZaxis()->GetTitle();
    else title = "";
    if(firstxbin!=1||lastxbin!= prof2d->GetXaxis()->GetNbins())
        title.Append(TString::Format(" yrange: %f - %f",prof2d->GetYaxis()->GetBinLowEdge(firstxbin),prof2d->GetYaxis()->GetBinLowEdge(lastxbin)));
    if(prof->GetZaxis()) prof->GetYaxis()->SetTitle(title);
    cout<<"xaxis: "<<prof->GetXaxis()->GetTitle()<<endl;
    cout<<"yaxis: "<<prof->GetYaxis()->GetTitle()<<endl;
    Int_t firstybin = 1;
    Int_t lastybin = prof2d->GetXaxis()->GetNbins();
    for(Int_t binx = firstybin; binx <= lastybin; binx++){
        name = prof->GetName()+TString::Format("_binx%d:",binx);
        TProfile* profyBin = new TProfile(name,title,binsX,xLow,xHigh);
        cout<<"new ybin: "<<binx<<" "<<prof->GetName()<<": "<<endl;
        profyBin->Reset();
        Int_t nEntries =0;
        for(Int_t biny = firstxbin; biny <= lastxbin; biny++){
            Int_t bin2d = prof2d->GetBin(binx,biny);
            (*profyBin)[biny] = (*prof2d)[bin2d];    // copy bin y values
            (*profyBin->GetSumw2())[biny] =  (*prof2d->GetSumw2())[bin2d];   // copy bin y*y values
            profyBin->SetBinEntries(biny, prof2d->GetBinEntries(bin2d) );    // copy bin entries
            nEntries+=prof2d->GetBinEntries(bin2d);
            // copy (if needed) bin sum of weight square
            if ( prof2d->GetBinSumw2()->fN > bin2d ) {
                profyBin->Sumw2();
                (*profyBin->GetBinSumw2())[biny] = (*prof2d->GetBinSumw2())[bin2d];
            }
        }
        profyBin->SetEntries(nEntries);
        prof->Add(profyBin);
        cout<<""<<profyBin->GetEntries()<<" ---> "<<prof->GetEntries()<<endl;
        delete profyBin;
    }
    cout<<" final profile: "<<prof->GetEntries()<<"/"<<prof2d->GetEntries()<<endl;
    //    char t; cin>>t;
    return prof;
}

void HistogrammSaver::SaveProfile2DWithEntriesAsText(TProfile2D* prof, bool drawStatBox){
    TString name = prof->GetName();
    if (name.First('h')==0)name[0]='c';
    TCanvas *c1 = new TCanvas(name);
    if (!drawStatBox)
        c1->SetObjectStat(false);
    prof->SetContour(100);
    prof->Draw("colz");
    TH2D* histo = prof->ProjectionXY(prof->GetName()+(TString)"_binEntries","B");
    histo->Draw("TEXTsame");
    SaveCanvas(c1);
}
void HistogrammSaver::SaveHistogram(TH2* histo, bool drawStatBox,bool optimizeRange,TString drawOption) {
    if (!histo)return;
    if(histo->GetEntries()==0)return;
    if (!drawStatBox)
        histo->SetStats(false);
    //	histo->SetStats(false);
    SaveHistogramPNG(histo,optimizeRange,drawOption);
    SaveHistogramROOT(histo,optimizeRange,drawOption);
}

void HistogrammSaver::SaveOverlay(TH2* histo,TString drawOption) {
    std::pair<Float_t,Float_t> readoutColumn = make_pair(75.,75.);
    std::pair<Float_t,Float_t> biasColumn = make_pair(0,0);
    if (!histo)return;
    if(histo->GetEntries()==0)return;
    TH1F* h = this->SaveProjectionZ(histo);
    delete h;
//    TH1F *h = this->GetBinsInHistogram(histo);
//    this->SaveHistogram(h);
//    delete h;
    Int_t style = gStyle->GetOptStat();
    if (settings->IsPaperMode())
        gStyle->SetOptStat("e");
    else
        gStyle->SetOptStat("ne");
    TString name = histo->GetName();
    if(name.First('h')==0)
        name.Replace(0,1,'c');
    else
        name.Insert(0,"c_");
    TCanvas *c1 = new TCanvas(name);
    c1->SetRightMargin(.15);
    histo->SetContour(100);
    histo->Draw("colz");
    histo->GetZaxis()->SetTitleOffset(1.2);
    TH1F* frame = c1->DrawFrame(0,0,150,150,histo->GetTitle());
    frame->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    frame->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    histo->Draw(drawOption+" same");
    vector<TMarker *> markers;
    vector<TCutG*> cells;


    markers.push_back(new TMarker(readoutColumn.first,readoutColumn.second,20));
    markers.back()->SetMarkerColor(kPink);
    markers.back()->SetMarkerSize(1.5);
    markers.back()->Draw();
    cells.push_back(GetCutGofBin("readoutBin",histo,readoutColumn.first,readoutColumn.second));
    if(cells.back()) cells.back()->SetLineColor(kBlack);
    if(cells.back()) cells.back()->SetLineWidth(7);
    if(cells.back()) cells.back()->Draw();
    for (Int_t i = 0; i < 4; i++){
        Float_t x = biasColumn.first;
        Float_t y = biasColumn.second;
        x+= settings->GetCellHeight()*(i%2);
        y+= settings->GetCellHeight()*(i/2);
        markers.push_back(new TMarker(x,y,20));
        markers.back()->SetMarkerColor(kBlack);
        markers.back()->SetMarkerSize(1.5);
        markers.back()->Draw();
        cells.push_back(GetCutGofBin(TString::Format("biasBin_%d",i),histo,x,y));
        if(cells.back()) cells.back()->SetLineColor(kBlack);
        if(cells.back()) cells.back()->SetLineWidth(7);
        if(cells.back()) cells.back()->Draw();
    }

    c1->Update();
    SaveCanvas(c1);
    delete c1;
}

TCutG* HistogrammSaver::GetCutGofBin(TString name, TH2* histo,Float_t x,Float_t y){
    if (!histo)return 0;

    Int_t bin = histo->FindBin(x,y);
    Int_t binx,biny,binz;
    histo->GetBinXYZ(bin,binx,biny,binz);
    if(verbosity>5)
        cout<<"HistogramSaver::GetCutGofBin \""<<name<<"\" "<<x<<" "<<y<<"->"<<bin<<": "<<binx<<"/"<<biny<<"/"<<binz<<endl;
    Float_t xlow = histo->GetXaxis()->GetBinLowEdge(binx);
    Float_t xup = histo->GetXaxis()->GetBinUpEdge(binx);
    Float_t ylow = histo->GetYaxis()->GetBinLowEdge(biny);
    Float_t yup = histo->GetYaxis()->GetBinUpEdge(biny);
    if(verbosity>5)
        cout<<"X:"<<xlow<<"-"<<xup<<"\tY:"<<ylow<<"-"<<yup<<endl;
    while (xup > histo->GetXaxis()->GetXmax()){
        Float_t delta = xup-xlow;
        xup = xlow;
        xlow = xup -delta;
    }
    while (yup > histo->GetYaxis()->GetXmax()){
        Float_t delta = yup-ylow;
        yup = ylow;
        ylow = yup -delta;
    }
    if(verbosity>5)
        cout<<"X:"<<xlow<<"-"<<xup<<"\tY:"<<ylow<<"-"<<yup<<endl;

    Float_t xx[] = {xlow,xup,xup,xlow,xlow};
    Float_t yy[] = {ylow,ylow,yup,yup,ylow};
    TCutG * cut = new TCutG(name,5,xx,yy);
    cut->SetLineWidth(3);
    return cut;
}

void HistogrammSaver:: Save1DProfileYWithFitAndInfluence(TH2* histo, TString function, bool drawStatbox){
    TString name = "fit_" + (TString)histo->GetName();
    TF1* fit = new TF1(name,function);
    return Save1DProfileYWithFitAndInfluence(histo,fit,drawStatbox);
}

void HistogrammSaver::Save1DProfileYWithFitAndInfluence(TH2* htemp,TF1* fit, bool drawStatbox){
    if(!fit)
        fit = new TF1("fit","pol1");
    TProfile *prof=0;
    if(!htemp)
        return;
    prof = htemp->ProfileY();
    Save1DProfileWithFitAndInfluence(prof,fit,drawStatbox);
    if(prof) delete prof;
}


TProfile* HistogrammSaver:: CreateAndSave1DProfileXWithFitAndInfluence(TH2* histo, TString function, bool drawStatbox){
    TString name = "fit_" + (TString)histo->GetName();
    if (function=="")
        function="pol1";
    TF1* fit = new TF1(name,function);
    return CreateAndSave1DProfileXWithFitAndInfluence(histo,fit,drawStatbox);
}

TProfile* HistogrammSaver::CreateAndSave1DProfileXWithFitAndInfluence(TH2* htemp,TF1* fit, bool drawStatbox){
    if(!fit)
        fit = new TF1("fit","pol1");
    TProfile *prof=0;
    if(!htemp)
        return prof;
    prof = htemp->ProfileX();
    if(prof){
        Save1DProfileWithFitAndInfluence(prof,fit,drawStatbox);
    }
    return prof;
}

void HistogrammSaver::Save1DProfileWithFitAndInfluence(TProfile *prof, TF1* fit, bool drawStatBox){
    if(prof==0) return;
    if (!prof) return;
    TCanvas* c1 = new TCanvas( (TString)("c_"+(TString)prof->GetName()) );
    c1->cd();
    if (!drawStatBox)
        prof->SetStats(false);
    TPaveText *text = 0;
    prof->Draw("");
    if(fit && !fit->IsZombie()&& prof){
        prof->Fit(fit,"Q");
        if (prof->GetXaxis() && prof->GetYaxis()){
            Float_t xmin = prof->GetXaxis()->GetXmin();
            Float_t xmax = prof->GetXaxis()->GetXmax();
            Float_t ymin = fit->GetMinimum(xmin,xmax);
            Float_t ymax = fit->GetMaximum(xmin,xmax);
            text = new TPaveText(.2,.2,.5,.3,"brNDC");
            text->SetFillColor(0);
            text->AddText(TString::Format("relative Influence: #frac{#Delta_{x}}{x_{max}} = %2.2f %%",(ymax-ymin)/ymax*100));
        }
    }
    prof->Draw("");
    if ( prof->GetXaxis() &&  prof->GetYaxis() && prof){
        Float_t xmin = prof->GetXaxis()->GetXmin();
        Float_t xmax = prof->GetXaxis()->GetXmax();
        Float_t ymin = prof->GetBinContent(prof->GetMinimumBin())-1.5*prof->GetBinError(prof->GetMinimumBin());
        Double_t ymax = prof->GetBinContent(prof->GetMaximumBin())+1.5*prof->GetBinError(prof->GetMaximumBin());
        ymax = TMath::Max(prof->GetYaxis()->GetXmax(),ymax);
        Float_t delta = ymax-ymin;
        ymax = (delta)*1.35+ymin;
        ymin = ymin - .05*delta;
        if(xmax<xmin){
            Float_t b = xmax;
            xmax = xmin;
            xmin = b;
        }
        if(ymax<ymin){
            Float_t b = ymax;
            ymax = ymin;
            ymin = b;
        }
        //    cout<<xmin<<"-"<<xmax<<" "<<ymin<<"-"<<ymax<<endl;
        TH1 *frame1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
        frame1->SetTitle(prof->GetTitle());
        frame1->GetXaxis()->SetTitle(prof->GetXaxis()->GetTitle());
        frame1->GetYaxis()->SetTitle(prof->GetYaxis()->GetTitle());
        frame1->Draw();
        prof->Draw("sames");
    }
    else if (prof)
        prof->Draw("");
    if(fit)fit->Draw("same");
    if(text)    text->Draw("same");
    SaveCanvas(c1);
}

void HistogrammSaver::SaveCanvas(TCanvas *canvas,TString name)
{
    if(canvas==0)
        return;
    SaveCanvasPNG(canvas,name);
    SaveCanvasROOT(canvas,name);
}

void HistogrammSaver::SaveGraph(TGraph* graph,std::string name,std::string option){
    if(graph==0)return;
    if(graph->GetN()==0)return;
    SaveGraphPNG(graph,name,option);
    SaveGraphROOT(graph,name,option);
}

void HistogrammSaver::SaveHistogramPDF(TH1F* histo) {
    if(!histo){
        cerr<<"HistogrammSaver::SaveHistogramPDF(TH1F*) \t histo == 0"<<endl;
        return;
    }
    if(histo->GetEntries()==0)return;
    TCanvas *plots_canvas = new TCanvas(TString::Format("cPdf_%s",histo->GetName()),TString::Format("c_%s",histo->GetName()));
    plots_canvas->cd();
    UInt_t maxBinX =histo->GetNbinsX();
    for(UInt_t i=histo->GetNbinsX();i>0;i--)
        if(histo->GetBinContent(i)==0)maxBinX=i;
    UInt_t minBinX =0;
    for(Int_t i=0;i<histo->GetNbinsX();i++)
        if(histo->GetBinContent(i)==0)minBinX=i;
    Float_t xmin = histo->GetXaxis()->GetBinLowEdge(minBinX);
    Float_t xmax = histo->GetXaxis()->GetBinLowEdge(maxBinX+1);
    histo->GetXaxis()->SetRangeUser(xmin,xmax);
    histo->Draw();
    TPaveText *pt2=0;
    if(pt)pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << histo->GetName() << ".pdf";
    plots_canvas->Print(plot_filename.str().c_str());
    //	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramPDF(TH2* histo) {
    if(histo==0)return;
    if(histo->GetEntries()==0)return;
    TCanvas *plots_canvas = new TCanvas(TString::Format("cPdf_%s",histo->GetName()),TString::Format("c_%s",histo->GetName()));
    plots_canvas->cd();
    //plots_canvas.cd();
    //	SetDuckStyle();
    if(verbosity)cout << "Using SaveHistogrammPDF on TH2 histogram " << histo->GetName() << endl;
    //histo->Draw();
    TPaveText *pt2=0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
    gStyle->SetTitleFont(42);
    gStyle->SetMarkerSize(0);
    if (pt2) pt2->SetTextSize(0.0250);
    if (pt2) pt2->SetTextColor(kBlack);
    histo->SetTitleFont(42);
    histo->UseCurrentStyle();
    histo->SetContour(100);
    histo->Draw("colz");
    if(pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << histo->GetName() << ".pdf";
    plots_canvas->Print(plot_filename.str().c_str());
    if(plots_canvas)delete plots_canvas;
    //pt->SetTextSize(runNumber0.1);
}

void HistogrammSaver::SaveHistogramPNG(TH1* histo,TString drawOption,bool optimizedRange) {
    if(histo==0)return;
    if(!histo){
        cout<<"Histogram is not existing..."<<endl;
        return;
    }
    if(histo->GetEntries()==0){
        if(verbosity)cout<<"Histogram "<<histo->GetName()<<" has no entries..."<<endl;
        return;
    }
    stringstream histoName;
    histoName<<histo->GetName()<<"_Clone";
    TH1* htemp=(TH1*)histo->Clone(histoName.str().c_str());
    if(htemp==0)return;
    TCanvas *plots_canvas = new TCanvas(TString::Format("cPng_%s",histo->GetName()),TString::Format("c_%s",histo->GetName()));
    plots_canvas->cd();
    drawOption.ToLower();
    if(drawOption.Contains("logy")){
        htemp->SetMinimum(1e-5);
        plots_canvas->SetLogy();
        cout<<"Draw " <<histo->GetName()<<" logy"<<endl;
        drawOption.Remove(drawOption.Index("logy"),4);
    }
    if(drawOption.Contains("logx")){
        cout<<"Draw " <<histo->GetName()<<" logx"<<endl;
        plots_canvas->SetLogx();
        drawOption.Remove(drawOption.Index("logx"),4);
    }
    if (optimizedRange)    htemp->SetMinimum(0);
    htemp->Draw(drawOption);
    TPaveText *pt2=0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
    if(pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << histo->GetName() << ".png";
    plots_canvas->Print(plot_filename.str().c_str());
    plots_canvas->Print(GetEpsFileName(histo->GetName()));
    //	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveCanvasROOT(TCanvas *canvas,TString name)
{
    if(!canvas)
        return;
    if (name=="")
        name = canvas->GetName();
    TCanvas* plots_canvas=(TCanvas*)canvas->Clone();
    plots_canvas->cd();

    TPaveText *pt2=0;
    if(pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvas->GetName()));
    if(pt2 && !settings->IsPaperMode()) pt2->Draw();

    TFile f(GetROOTFileName(name),"UPDATE");
    canvas->Write();
}

void HistogrammSaver::SaveCanvasPNG(TCanvas *canvas, TString name)
{
    if(canvas==0)
        return;
    if (name=="")
        name = canvas->GetName();
    canvas->cd();
    TPaveText *pt2=0;
    if(pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvas->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << name<<".png";
    canvas->Print(plot_filename.str().c_str());
    canvas->Print(GetEpsFileName(name));
}

void HistogrammSaver::SaveGraphPNG(TGraph* graph,string name,string option){
    if(!graph)
        return;
    if(graph->GetN()==0)return;
    TCanvas plots_canvas(TString::Format("c_%s",name.c_str()),TString::Format("c_%s",name.c_str()));
    plots_canvas.cd();
    graph->Draw(option.c_str());

    TPaveText *pt2 = 0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",graph->GetName()));
    if(pt2 && !settings->IsPaperMode()) pt2->Draw();

    ostringstream plot_filename;
    plot_filename << plots_path << name << ".png";
    plots_canvas.Print(plot_filename.str().c_str());
    plots_canvas.Print(GetEpsFileName(name));
}

void HistogrammSaver::SaveHistogramFitGaussPNG(TH1* htemp) {
    if(!htemp)
        return;
    TH1* histo = (TH1*)htemp->Clone(TString::Format("%s_Clone",htemp->GetName()));
    if(histo->GetEntries()==0)return;

    TF1 histofitx("histofitx","gaus",histo->GetMean()-2*histo->GetRMS(),histo->GetMean()+2*histo->GetRMS());
    histofitx.SetLineColor(kBlue);
    histo->Fit(&histofitx,"rq");

    TCanvas plots_canvas("plots_canvas","plots_canvas");
    plots_canvas.cd();
    histo->Draw();
    histofitx.Draw("same");
    TPaveText *pt2=0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",htemp->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();

    ostringstream plot_filename;
    plot_filename << plots_path << htemp->GetName() << ".png";
    plots_canvas.Print(plot_filename.str().c_str());
    plots_canvas.Print(GetEpsFileName(htemp->GetName()));
}

void HistogrammSaver::SaveHistogramROOT(TH1* htemp) {
    if(!htemp)return;
    //	if(htemp->GetEntries()==0)return;

//    ostringstream plots_filename;
//    plots_filename << plots_path<<"/" << htemp->GetName() << "."<<runNumber<<".root";
    TCanvas *plots_canvas =  new TCanvas(TString::Format("c_%s", htemp->GetName()), TString::Format("croot_%s", htemp->GetName()));
    plots_canvas->cd();
    TH1* histo = (TH1*)htemp->Clone();
    if(!histo)
        return;

    TPaveText *pt2 = 0;
    if(pt) pt2 = (TPaveText*)pt->Clone(TString::Format("ptcRoot_%s",htemp->GetName()));

    plots_canvas->Clear();
    plots_canvas->cd();
    histo->Draw();
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    plots_canvas->Draw();
    histo->Draw();

    //write to own root File
    plots_canvas->Write(GetROOTFileName(htemp->GetName()));
//    plots_canvas->Write(plots_filename.str().c_str());
    TFile *f = GetFilePointer(GetROOTFileName());
    if(!f)
        return;
    f->cd();
    TCanvas *plots_canvas2 = (TCanvas*) plots_canvas->Clone(TString::Format("cc_%s",htemp->GetName()));
    if(!plots_canvas2)
        return;
    //add to histograms.root
    f->cd();
    plots_canvas2->Write();
    TH1* h_clone = (TH1*)htemp->Clone();
    h_clone->SetDirectory(f);
    h_clone->Write();

    if (f && !f->IsZombie() && f->IsOpen())
        f->Close();
    //	if(plots_canvas)delete plots_canvas;

}

void HistogrammSaver::SaveHistogramPNG(TH2* histo,bool optimizeRange,TString drawOption) {

    if(!histo){
        cerr<<"HistogrammSaver::SaveHistogramPNG(TH2*), histogram ==0"<<endl;
        return;
    }
    if(histo->GetEntries()==0)return;
    //	gROOT->SetStyle("Plain_RD42_2D");
    //	gROOT->ForceStyle(true);
    TCanvas *plots_canvas =  new TCanvas(TString::Format("cPng_%s", histo->GetName()), TString::Format("c_%s", histo->GetName()));
    plots_canvas->Clear();
    plots_canvas->cd();
    TH2* htemp = (TH2*)histo->Clone();
    if(optimizeRange) HistogrammSaver::OptimizeXYRange(htemp);
    htemp->Draw(drawOption);

    TPaveText *pt2=0;
    if (pt) pt2 = (TPaveText*)pt->Clone(TString::Format("ptPng_%s",histo->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << histo->GetName() << ".png";
    plots_canvas->Print(plot_filename.str().c_str());
    plots_canvas->Print(GetEpsFileName(histo->GetName()));

    //	gROOT->SetStyle("Plain_RD42");
    //currentStyle->cd();
    //	if(plots_canvas)delete plots_canvas;
    if (htemp) delete htemp;
    if (plots_canvas) delete plots_canvas;
}

void HistogrammSaver::SaveHistogramROOT(TH2* histo,bool optimizeRange,TString drawOption) {
    if(!histo){
        cerr<<"HistogrammSaver::SaveHistogramROOT(TH2*) histogram == 0"<<endl;
        return;
    }
    if(histo->GetEntries()==0)return;
    TCanvas *plots_canvas =  new TCanvas(TString::Format("cRoot_%s", histo->GetName()), TString::Format("c_%s", histo->GetName()));
    plots_canvas->Clear();

    plots_canvas->cd();
    TH2* htemp = (TH2*)histo->Clone();
    if(htemp==0)
        return;
    htemp->Draw();
    if(optimizeRange)
        HistogrammSaver::OptimizeXYRange(htemp);
    htemp->Draw(drawOption);

    TPaveText *pt2 = 0;
    if (pt) (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    ostringstream plot_filename;
    plot_filename << plots_path << histo->GetName() << "."<<runNumber<<".root";
    plots_canvas->Print(GetROOTFileName(histo->GetName()));
    TString rootFileName = GetROOTFileName();
    TFile *f = this->GetFilePointer(rootFileName,"UPDATE");
    f->cd();
    plots_canvas->Write();
    f->Close();
    if (htemp) delete htemp;
    if (plots_canvas) delete plots_canvas;

    //	if (plots_canvas) delete plots_canvas;
}


void HistogrammSaver::SaveHistogramROOT(TH3F* histo){
    if(!histo){
        cerr<<"HistogrammSaver::SaveHistogramROOT(TH2*) histogram == 0"<<endl;
        return;
    }
    if(histo->GetEntries()==0)return;
    TH3F* htemp = (TH3F*)histo->Clone();
    if(htemp==0)
        return;
    htemp->Draw();
    TString fileName = GetROOTFileName(histo->GetName());
    htemp->Write(fileName);
    htemp->Write();
    htemp->Print(GetROOTFileName(histo->GetName()));
    stringstream histo_filename;
    TString rootFileName = GetROOTFileName();
    TFile *f = GetFilePointer(rootFileName,"UPDATE");
    f->cd();
    htemp->Write();
    f->Close();
    if (htemp)
        delete htemp;
}


void HistogrammSaver::SaveGraphROOT(TGraph* graph,std::string name,std::string option){
    if(!graph) {
        cerr<<"HistogrammSaver::SaveGraphROOT(TGraph* ) graph == 0"<<endl;
        return;
    }
    if(graph->GetN()==0)return;
    //	TCanvas *plots_canvas = ((TCanvas *)(gROOT->GetListOfCanvases()->FindObject("plots_canvas")));
    //	if (plots_canvas) plots_canvas->Clear();
    //	else plots_canvas = new TCanvas("plots_canvas", "plots_canvas");
    TString cName = name;//histo->GetName();
    //        if (cName.First('h') == 0)
    //            cName.Replace(0,1,"c");
    //        else
    //            cName = (TString)"c_"+cName;
    TCanvas *plots_canvas =  new TCanvas(cName,cName );
    plots_canvas->Clear();

    plots_canvas->cd();
    TGraph* gTemp = (TGraph*)graph->Clone();
    gTemp->Draw(option.c_str());

    TPaveText *pt2=0;
    if(pt) pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",graph->GetName()));
    if (pt2 && !settings->IsPaperMode()) pt2->Draw();
    plots_canvas->Print(GetROOTFileName(name));
    //	if(plots_canvas)	delete plots_canvas;
}

void HistogrammSaver::SetVerbosity(unsigned int i)
{
    this->verbosity=i;
    if(verbosity) cout<<"HistogrammSaver::Set Verbosity ON"<<endl;
}

//
//void HistogrammSaver::SaveCanvasRoot(TCanvas *canvas, string location, string file_name)
//{
//    if(canvas==0)
//        return;
//    char loc[500];
//    memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//    char rt[] = ".root";
//    char *rtloc = loc;
//    //Saving .root file
//    strcat(rtloc,file_name.c_str());
//    strcat(rtloc,rt);
//    char const *rt_file = &rtloc[0];
//    TObjArray list(0);
//    list.Add(canvas);
//    TFile f(rt_file,"recreate");
//    list.Write();
//    f.Close();
//    cout << ".root file was created at: " << rt_file << endl;
//}

//void SaveCanvasC(TCanvas *canvas, char* location, char* file_name);
//void SaveCanvasC(TCanvas *canvas, string location, string file_name)
//{
//    if(canvas==0 )return;
//    char loc[500];
//    memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//    char cmac[] = ".C";
//    char *cmacloc = loc;
//
//    //Saving .C macro
//    strcat(cmacloc,file_name.c_str());
//    strcat(cmacloc,cmac);
//    char const *cmac_file = &cmacloc[0];
//    canvas->SaveSource(cmac_file);
//    cout << ".c macro was created at: " << cmac_file << endl;
//}
//
//void HistogrammSaver::SaveCanvasPNG(TCanvas *canvas, string location, string file_name)
//{
//    if(canvas==0)return;
//    char loc[500];
//    memcpy(loc,location.c_str(),strlen(location.c_str())+1);
//    char png[] = ".png";
//    char *pngloc = loc;
//
//    //Save .png file
//    strcat(pngloc,file_name.c_str());
//    strcat(pngloc,png);
//    char const *png_file = &pngloc[0]; //assigns address of first element of the file string to the char_file pointer
//    TImage *img = TImage::Create();
//    img->FromPad(canvas);
//    img->WriteImage(png_file);
//    cout << ".png file was created at: " << png_file << endl;
//    delete img;
//}

void HistogrammSaver::SetDuckStyle() {
    TStyle* DuckStyle = new TStyle("DuckStyle", "Famous Duck Style");
    //	DuckStyle->SetOptFit(1111);  //Fit options to be displayed
    DuckStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
    DuckStyle->SetPadRightMargin(0.15);
    DuckStyle->SetPadTopMargin(0.15);
    DuckStyle->SetTitleColor(kBlack);
    DuckStyle->SetFrameLineWidth(1);
    DuckStyle->SetStatFont(42);
    DuckStyle->SetStatBorderSize(0);
    DuckStyle->SetPadBorderSize(0);
    DuckStyle->SetPadTopMargin(1);
    DuckStyle->SetLegendBorderSize(0);
    DuckStyle->SetStatFontSize(0.02);
    DuckStyle->SetStatStyle(0);
    DuckStyle->SetStatH(0.12); //Sets Height of Stats Box
    DuckStyle->SetStatW(0.15); //Sets Width of Stats Box
    DuckStyle->SetStatX(0.9);
    DuckStyle->SetStatY(0.97);
    DuckStyle->SetTitleOffset(1.0,"Y");
    DuckStyle->SetPalette(53); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
    DuckStyle->SetCanvasBorderMode(0);
    DuckStyle->SetTitleFont(42,"XYZ");
    DuckStyle->SetTitleFontSize(0.038);
    //	DuckStyle->SetTitleTextSize(0.03);
    DuckStyle->SetTitleTextColor(kBlack);
    DuckStyle->SetFrameLineStyle(0);
    DuckStyle->SetGridStyle(0);
    DuckStyle->SetHatchesLineWidth(0);
    //	DuckStyle->SetOptTitle(0);
    DuckStyle->SetPadTickX(0);
    DuckStyle->SetPadTickY(0);
    DuckStyle->SetTitleX(0.07);
    DuckStyle->SetTitleY(0.925);

    DuckStyle->SetTitleSize(0.02,"XYZ");
    DuckStyle->SetLineWidth(1);
    DuckStyle->SetHistLineWidth(1);
    //	DuckStyle->SetTitleStyle(0);
    DuckStyle->SetTitleBorderSize(0);
    DuckStyle->SetTitleFillColor(0);
    DuckStyle->SetHistLineWidth(1);

    //	DuckStyle->SetTickLength(0,"XY");
    //	gStyle->SetBarOffset(0.5);
    //	gStyle->SetStatFont(42);
    //	gStyle->SetTextSize(0.01);
    DuckStyle->SetLabelFont(42,"XYZ");
    DuckStyle->SetLabelColor(kBlack,"XYZ");
    DuckStyle->SetLabelSize(0.04,"XYZ");
    //DuckStyle->SetTitleOffset(1.8, "Y"); // Another way to set the Offset
    //	gStyle->SetTitleOffset(1.2, "X"); // Another way to set the Offset
    DuckStyle->SetTitleOffset(1.2,"X");
    DuckStyle->cd();

    cout << "Using DuckStyle" << endl;
}


/**
 * @brief creates a scatter histogram with posX_vs_posY as an input
 *
 * @return TH3F histogram
 */
TH3F* HistogrammSaver::Create3DHisto(std::string name, std::vector<Float_t> posX, std::vector<Float_t> posY, std::vector<Float_t> posZ,UInt_t nBinsX, UInt_t nBinsY,UInt_t nBinsZ,
        Float_t minRangeX,Float_t maxRangeX,Float_t minRangeY,Float_t maxRangeY,Float_t minRangeZ,Float_t maxRangeZ, Float_t factor)
{
    if(posX.size()!=posY.size()||posX.size()!=posZ.size()||posX.size()==0) {
        cerr<<"ERROR HistogrammSaver::CreateScatterHisto vectors have different size "<<posX.size()<<" "<<posY.size()<<" "<<name<<endl;
        return new TH3F();
    }
    cout<<"Creating 3dHisto: "<<name<<endl;
    cout<<TString::Format("maxRange:\nX: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minRangeX,maxRangeX,minRangeY,maxRangeY,minRangeZ,maxRangeZ)<<endl;
    Float_t maxX = posX.at(0);
    Float_t maxY = posY.at(0);
    Float_t maxZ = posZ.at(0);
    Float_t minX = posY.at(0);
    Float_t minY = posY.at(0);
    Float_t minZ = posZ.at(0);
    //	cout<<" Create Histo: '"<<name<<"' - Range ("<<minRangeX<<"-"<<maxRangeX<<"),  ("
    //			<<minRangeY<<"-"<<maxRangeY<<"), ("<<minRangeZ<<"-"<<maxRangeZ<<")"<<endl;
    for(UInt_t i=0;i<posX.size();i++){
        if (posX.at(i)<minRangeX||posX.at(i)>maxRangeX)
            continue;
        if (posY.at(i)<minRangeY||posY.at(i)>maxRangeY)
            continue;
        if (posZ.at(i)<minRangeZ||posZ.at(i)>maxRangeZ)
            continue;
        if(posX.at(i)>maxX)maxX=posX.at(i);
        else if(posX.at(i)<minX)minX=posX.at(i);
        if(posY.at(i)>maxY)maxY=posY.at(i);
        else if(posY.at(i)<minY)minY=posY.at(i);
        if(posZ.at(i)>maxZ)maxZ=posZ.at(i);
        else if(posZ.at(i)<minZ)minZ=posZ.at(i);
    }
    //	cout<<TString::Format("X: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minX,maxX,minY,maxY,minZ,maxZ)<<endl;
    Float_t factorX = factor;
    Float_t factorY = factor;
    Float_t factorZ = factor;
    Float_t deltaXMax = maxRangeX - minRangeX;
    Float_t deltaYMax = maxRangeY - minRangeY;
    Float_t deltaZMax = maxRangeZ - minRangeZ;
    Float_t maxDiff = 0.02;
    if ( TMath::Abs(maxRangeX-maxX)/deltaXMax <= maxDiff && TMath::Abs(minRangeX - minX)/deltaXMax <= maxDiff ) {
        factorX = 0;
        maxX = maxRangeX;
        minX = minRangeX;
    }
    if ( TMath::Abs(maxRangeY-maxY)/deltaYMax <= maxDiff && TMath::Abs(minRangeY - minY)/deltaYMax <= maxDiff ) {
        factorY = 0;
        minY = minRangeY;
        maxY = maxRangeY;
    }
    if ( TMath::Abs(maxRangeZ-maxZ)/deltaZMax <= maxDiff && TMath::Abs(minRangeZ - minZ)/deltaZMax <= maxDiff ) {
        factorZ = 0;
        minZ = minRangeZ;
        maxZ = maxRangeZ;
    }
    Float_t deltaX=maxX-minX;
    Float_t deltaY=maxY-minY;
    Float_t deltaZ=maxZ-minZ;
    //	cout<<"\t"<<deltaX<<" "<<deltaY<<" "<<deltaZ<<endl;
    //	cout<<"\t"<<factorX<<" "<<factorY<<" "<<factorZ<<endl;
    minX = minX-factorX*deltaX;
    maxX = maxX+factorX*deltaX;
    minY = minY-factorY*deltaY;
    maxY = maxY+factorY*deltaY;
    minZ = minZ-factorZ*deltaZ;
    maxZ = maxZ+factorZ*deltaZ;
    cout<<TString::Format("X: [%f,%f],\tY: [%f,%f],\tZ: [%f,%f]",minX,maxX,minY,maxY,minZ,maxZ)<<endl;
    //	char t; cin>>t;
    TH3F* histo = new TH3F(name.c_str(),name.c_str(),
            nBinsX,minX,maxX,
            nBinsY,minY,maxY,
            nBinsZ,minZ,maxZ);
    for(UInt_t i=0;i<posX.size();i++){
        if (posX.at(i)<minRangeX||posX.at(i)>maxRangeX)
            continue;
        if (posY.at(i)<minRangeY||posY.at(i)>maxRangeY)
            continue;
        if (posZ.at(i)<minRangeZ||posZ.at(i)>maxRangeZ)
            continue;
        histo->Fill(posX.at(i),posY.at(i),posZ.at(i));
    }
    histo->GetXaxis()->SetTitle("X-Position");
    histo->GetYaxis()->SetTitle("Y-Position");
    histo->GetZaxis()->SetTitle("Z-Position");
    histo->Draw();
    int minXbin = histo->GetXaxis()->FindBin(minX);
    int maxXbin = histo->GetXaxis()->FindBin(maxX);
    histo->GetXaxis()->SetRange(minXbin,maxXbin);
    int minYbin = histo->GetYaxis()->FindBin(minY);
    int maxYbin = histo->GetYaxis()->FindBin(maxY);
    histo->GetYaxis()->SetRange(minYbin,maxYbin);
    int minZbin = histo->GetZaxis()->FindBin(minZ);
    int maxZbin = histo->GetZaxis()->FindBin(maxZ);
    histo->GetZaxis()->SetRange(minZbin,maxZbin);

    return histo;
}
/**
 * @brief creates a scatter histogram with posX_vs_posY as an input
 *
 * @return TH2F histogram
 */
TH2F* HistogrammSaver::CreateScatterHisto(std::string name, std::vector<Float_t> posY, std::vector<Float_t> posX,
        UInt_t nBinsX, UInt_t nBinsY, Float_t minRangeX,Float_t maxRangeX, Float_t minRangeY, Float_t maxRangeY,Float_t factor, bool auto_limits)
{
    //	Float_t factor = 0.05;//5% bigger INtervall...
    if(posX.size()!=posY.size()||posX.size()==0) {
        cerr<<"ERROR HistogrammSaver::CreateScatterHisto vectors have different size "<<posX.size()<<" "<<posY.size()<<" "<<name<<endl;
        return new TH2F();
    }
    TString hName = name;
    bool  verb = hName.Contains("_ScatterPlot_RelHitPosXPredDet_vs_UseEvent_Plane_");
    if (hName.Contains("hCalcModeCutProblem")) verb = true;
    if (verb)    cout<<"[HistogrammSaver::CreateScatterHisto]" <<name <<endl;;
    UInt_t entries = posX.size();
    if (verb) cout<<" entries: "<<entries;
    std::vector<Float_t> posX2 = posX;
    std::vector<Float_t> posY2 = posY;
    std::sort(posX2.begin(),posX2.end());
    std::sort(posY2.begin(),posY2.end());
    UInt_t nLow = 0;
    UInt_t nUp  = entries-1;
    Float_t minX,maxX,minY,maxY,deltaX,deltaY;
    Int_t calcMethod = 0;
    if (maxRangeX == (+1) * std::numeric_limits<float>::infinity() ||
            minRangeX == (-1) * std::numeric_limits<float>::infinity() ||
            maxRangeY == (+1) * std::numeric_limits<float>::infinity() ||
            minRangeY == (-1) * std::numeric_limits<float>::infinity() )
        calcMethod = 1;
    if (verb) cout<<" CalcMethod = "<<calcMethod<<endl;
    if (calcMethod == 1){
        nLow = .05 * entries;
        nUp = .95 * entries;
        if (nUp >= entries)
            nUp = entries-1;
    }
    if (verb) cout<<" nLow "<<nLow<<"/"<<posX2.at(nLow)<<"/"<<posY2.at(nLow)<<endl;
    if (verb) cout<<" nUp  "<<nUp<<"/" <<posX2.at(nUp)<<"/"<<posY2.at(nUp)<<endl;
    while (posX2.at(nLow)!=posX2.at(nLow) && nLow < entries-1)
        nLow++;
    while (posX2.at(nUp)!=posX2.at(nUp) && nUp>0 )
        nUp--;
    if (verb) cout<<" nLow "<<nLow<<"/"<<posX2.at(nLow)<<"/"<<posY2.at(nLow)<<endl;
    if (verb) cout<<" nUp  "<<nUp<<"/" <<posX2.at(nUp)<<"/"<<posY2.at(nUp)<<endl;
    if (calcMethod == 1){
        minX = posX2.at(nLow);
        maxX = posX2.at(nUp);
        deltaX = maxX-minX;
        minX = minX-.1*deltaX;
        maxX = maxX+.1*deltaX;
        minY = posY2.at(nLow);
        maxY = posY2.at(nUp);
        deltaY = maxY-minY;
        minY = minY-.1*deltaY;
        maxY = maxY+.1*deltaY;
    }
    else{
        minX = posX2.at(nLow);
        maxX = posX2.at(nUp);
        deltaX = maxX-minX;
        minY = posY2.at(nLow);
        maxY = posY2.at(nUp);
        deltaY = maxY-minY;
    }

    if (verb) cout <<" X:"<<minX<<"-"<<maxX<<"\tY:"<<minY<<"-"<<maxY<<" --- "<<posX2.at(0)<<"-"<<posX2.at(entries-1)<<"/"<<posY2.at(0)<<"-"<<posY2.at(entries-1)<<endl;
    if (minX < minRangeX) minX = minRangeX;
    if (maxX > maxRangeX) maxX = maxRangeX;
    if (minY < minRangeY) minY = minRangeY;
    if (maxY > maxRangeY) maxY = maxRangeY;
    if (verb) cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
    deltaX=maxX-minX;
    deltaY=maxY-minY;
	if (!auto_limits) {
		minX = minRangeX;
		maxX = maxRangeX;
		minY = minRangeY;
		maxY = maxRangeY;
		deltaX = 0.;
		deltaY = 0.;
	}
    TH2F* histo = new TH2F(name.c_str(),name.c_str(),nBinsX,minX-factor*deltaX,maxX+factor*deltaX,nBinsY,minY-factor*deltaY,maxY+factor*deltaY);
    for(UInt_t i=0;i<posX.size();i++){
        if (verb){
            if ((abs(posX.at(i)) >50) ||(abs(posX.at(i)) >50) )
                cout <<"Filling: " << i<< "/"<<posX.at(i)<<" "<<posY.at(i)<<endl;;
        }
        if (posX.at(i) < minRangeX || posX.at(i) > maxRangeX)
            continue;
        if (posY.at(i) < minRangeY || posY.at(i) > maxRangeY)
            continue;
        histo->Fill(posX.at(i),posY.at(i));
    }
    histo->GetXaxis()->SetTitle("X-Position");
    histo->GetYaxis()->SetTitle("Y-Position");

    //set xrange
    TH1D* hProj=histo->ProjectionX(histo->GetName()+(TString)"__px");
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;

    //set yRange
    hProj=histo->ProjectionY(histo->GetName()+(TString)"__py");
    int binyMin=0;
    for(binyMin=0;binyMin<hProj->GetNbinsX();binyMin++)if(hProj->GetBinContent(binyMin))break;
    int binyMax;
    for(binyMax=hProj->GetNbinsX();binyMax>0;binyMax--)if(hProj->GetBinContent(binyMax))break;

    histo->GetYaxis()->SetRangeUser(hProj->GetBinLowEdge(binyMin-1),hProj->GetBinLowEdge(binyMax+1));
    delete hProj;

    return histo;
}

TGraph HistogrammSaver::CreateDipendencyGraph(std::string name, std::vector<Float_t> vecY, std::vector<Float_t> vecX,ULong_t maxSize)
{
    if(vecY.size()!=vecX.size()||vecX.size()==0) {
        cerr<<"ERROR HistogrammSaver::CreateDipendencyGraph vectors have different size "<<vecY.size()<<" "<<vecX.size()<<": "<<name<<endl;
        return TGraph();
    }
    //cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
    if (maxSize==-1)
        maxSize = vecY.size();
    ULong_t size = TMath::Min(maxSize,(ULong_t)vecY.size());
    TGraph hGraph = TGraph(size,&vecX.at(0),&vecY.at(0));
    hGraph.Draw("AP");
    hGraph.GetXaxis()->SetName("X axis");
    hGraph.GetYaxis()->SetName("y axis");
    hGraph.SetTitle(name.c_str());
    return hGraph;
}

void HistogrammSaver::CopyAxisRangesToHisto(TH1F* changingHisto,TH1F* axisInputHisto){
    if(axisInputHisto&&changingHisto){
        changingHisto->Draw("goff");
        axisInputHisto->Draw("goff");
        Float_t xmin = axisInputHisto->GetXaxis()->GetXmin();
        Float_t xmax = axisInputHisto->GetXaxis()->GetXmax();
        Float_t ymin = axisInputHisto->GetYaxis()->GetXmin();
        Float_t ymax = axisInputHisto->GetYaxis()->GetXmax();
        if(ymax==1)
            ymax= axisInputHisto->GetBinContent(axisInputHisto->GetMaximumBin());
        changingHisto->GetXaxis()->SetRangeUser(xmin,xmax);
        changingHisto->GetYaxis()->SetRangeUser(ymin,ymax);
        cout<<"copyAxisRangeToHisto: x: "<<xmin<<"-"<<xmax<<"\ty:"<<ymin<<"-"<<ymax<<endl;
    }
    else
        cerr<<"HistogrammSaver::CopyAxisRangesToHisto::\tOne of the histogram is a pointer to Null: "<<changingHisto<<" "<<axisInputHisto<<endl;
}

TGraphErrors HistogrammSaver::CreateErrorGraph(std::string name, std::vector<Float_t> x, std::vector<Float_t> y, std::vector<Float_t> ex, std::vector<Float_t> ey)
{
    if(x.size()!=y.size()||x.size()!=ex.size()||ex.size()!=ey.size()||x.size()==0) {
        cerr<<"ERROR HistogrammSaver::CreateErrorGraph vectors have different size "<<x.size()<<" "<<y.size()<<endl;
        return TGraphErrors();
    }

    cout<<"HistogrammSaver::CREATE CreateErrorGraph:\""<<name<<"\" with "<<x.size()<<" Entries"<<endl;
    TGraphErrors hGraph = TGraphErrors(x.size(),&x.at(0),&y.at(0),&ex.at(0),&ey.at(0));
    hGraph.SetTitle(name.c_str());
    return hGraph;
}

TH2F* HistogrammSaver::CreateDipendencyHisto(std::string name, std::vector<Float_t> Delta, std::vector<Float_t> pos, UInt_t nBins)
{
    TH2F *histo = CreateScatterHisto(name,pos,Delta,nBins);
    histo->GetXaxis()->SetTitle("Position");
    histo->GetYaxis()->SetTitle("Difference");
    return histo;
}

void HistogrammSaver::SetRange(Float_t min,Float_t max){
    if (min<max){
        this->xRangeMin=min;
        this->xRangeMax=max;
    }
}


Float_t HistogrammSaver::GetMean(std::vector<Float_t> vec){
    Float_t mean = 0;
    Float_t mean2 = 0;
    Float_t nEntries = vec.size();
    for(UInt_t i=0;i<vec.size();i++){
        mean+=vec.at(i);
        mean2+=vec.at(i)*vec.at(i);
    }
    mean=mean/nEntries;
    mean2=mean2/nEntries;
    Float_t sigma = TMath::Sqrt(mean2-mean*mean);
    cout<<"Mean: "<<mean*100<<" +/- " <<sigma*100<<"\t"<<vec.size() << mean<<"/"<<mean2<<endl;
    return mean;
}
TH1F* HistogrammSaver::CreateDistributionHisto(std::string name, std::vector<Float_t> vec, UInt_t nBins,EnumAxisRange range,Float_t xmin,Float_t xmax, Float_t factor)
{
    TString hName = name;
    bool verb = hName.BeginsWith("hSilicon_PostAlignment_Distribution_DeltaY_Plane_0");
    verb = verb || hName.BeginsWith("hDiamond_PostAlignment_ScatterPlot_RelHitPosXPredDet_vs_UseEvent_Plane_");
    verb = verb || hName.BeginsWith("hDiamond_PostAlignment_Distribution_DeltaX_Plane_");
    if (verb) cout<<"\n[CreateDistributionHisto]"<<hName<<endl;
    int verbosity = verb*6;
    //	Float_t factor = 0.05;//5% bigger INtervall...
    if(vec.size()==0)
        return new TH1F(name.c_str(),name.c_str(),nBins,0.,1.);
    std::vector<Float_t>vec2 = vec;
    std::sort (vec2.begin(), vec2.end());
    int entries = vec2.size();
    int low = entries *.05;
    int up = entries *.95;
    if (verbosity>3)
        cout<<"Kicking out max and mins: "<<entries<<" "<<low<<"-"<<up<<"\tLOW: "<<
            vec2.at(0)<<"->"<<vec2.at(low)<<"\tUP: "<<vec2.at(up)<<"<--"<<vec2.at(entries-1)<<endl;
    Float_t max = vec.at(low);
    Float_t min = vec.at(low);
    if(verbosity>3)cout<<"Create Histo "<<name<<", mode "<<range<<" "<<flush;
    if (range==maxWidth){
        for(UInt_t i=0;i<vec.size();i++){
            if (max<vec.at(i))max=vec.at(i);
            if (min>vec.at(i))min=vec.at(i);
        }
        Float_t delta = max-min;
        min =min-delta*factor;
        max=max+delta*factor;
        if(min-max==0){
            min-=0.5*min;
            max+=0.5*min;
        }
        if(verbosity>3)cout<<" maxWidth "<<min <<"-"<<max<<endl;
    }
    else if(range==fiveSigma||range==threeSigma){
        Float_t  mean2 =0;
        Float_t sigma2 = 0;
        int n=0;
        for(UInt_t i=low;i<=up && i < entries;i++){
            Float_t x = vec2.at(i);
            if(x<xmin||x>xmax)
                continue;
            mean2+=x;
            sigma2+=x*x;
            n++;
        }
        mean2/=(Float_t)n;
        sigma2/=(Float_t)n;
        if(verbosity>3)cout<<"mean2: "<<mean2<<"\tsigma2:"<<sigma2<<endl;

        Float_t mean=0;
        Float_t sigma=0;
        UInt_t nEvents=0;
        Float_t nsigma = range==fiveSigma?5:3;
        for(UInt_t i=0;i<vec.size();i++){
            Float_t x = vec.at(i);
            if(x<xmin||x>xmax)
                continue;
            if( ((x-mean2)/sigma2)>nsigma)
                continue;
            mean+=x;
            sigma+=x*x;
            nEvents++;
        }
        mean/=(Float_t)nEvents;
        sigma/=(Float_t)nEvents;
        sigma = sigma -mean*mean;
        sigma=TMath::Sqrt((Double_t)sigma);
        if(sigma==0)
            sigma = .5 *mean;
        UInt_t nSigma = (range==fiveSigma)? 5:3;
        max=mean+nSigma*sigma;
        min=mean-nSigma*sigma;
        if(verbosity>3)cout<<""<<nSigma<<"Sigma: "<<mean<<"+/-"<<sigma<<" ==> "<<min <<"-"<<max<<endl;
    }
    else if(range==positiveArea){
        min=0;
        for(UInt_t i=0;i<vec.size();i++)
            if (max<vec.at(i))max=vec.at(i);
        max*=(1+factor);
        if(verbosity>3)cout<<" positiveArea: 0 -"<<max;
    }
    else if(range==positiveSigma){
        min=0;
        Float_t mean=0;
        Float_t sigma=0;
        for(UInt_t i=0;i<vec.size();i++){
            mean+=vec.at(i);
            sigma+=vec.at(i)*vec.at(i);
        }
        mean/=(Float_t)vec.size();
        sigma/=(Float_t)vec.size();
        sigma = sigma -mean*mean;
        sigma=TMath::Sqrt(sigma);
        UInt_t nSigma = 3;
        max=mean+nSigma*sigma;
        if(verbosity>3)cout<<" positiveSigma: 0 - "<<max<<endl;
    }
    else if(range==manual){
        max = xmax;
        min = xmin;
        if(verbosity>3)cout<<" manual: "<<min<<" - " << max <<endl;
    }

    TH1F* histo = new TH1F(name.c_str(),name.c_str(),nBins,min,max);
    for(UInt_t i=0;i<vec.size();i++){
        Float_t x = vec.at(i);
//        if(x<xmin||x>xmax)
//            continue;
        histo->Fill(vec.at(i));
    }
    int ntries=0;
    entries = histo->GetEntries();
    while (true){
        Double_t max2 = histo->GetBinContent(histo->GetMaximumBin());
        Double_t fraction = max2 / (Double_t)entries;
        if (histo->GetNbinsX() < 20)
            break;
        if (max2>0.05)
            break;
        if(ntries>=3)//todo change hardcoding
            break;
        histo->Rebin();
        ntries++;
    }
    histo->Draw("goff");
    Int_t binsx = histo->GetNbinsX();
    if (histo->GetXaxis()->GetBinCenter(1) > min)
        min = histo->GetXaxis()->GetBinLowEdge(1);
    if (histo->GetXaxis()->GetBinCenter(binsx) <= max)
        max = histo->GetXaxis()->GetBinUpEdge(binsx);
    if(verbosity>3)cout<<"Update Xaxis Range: "<<min<<"-"<<max<<endl;
    histo->GetXaxis()->SetRangeUser(min,max);
    histo->GetYaxis()->SetTitle("number of entries #");
    return histo;
}

std::pair<Float_t, Float_t> HistogrammSaver::OptimizeXRange(TH1F* histo){
    histo->Draw();
    Float_t xmax = histo->GetXaxis()->GetXmax();
    Int_t maxBin = 0;
    Int_t minBin = 0;
    Int_t i=0;
    for( i=0;i<histo->GetNbinsX()&&histo->GetBinContent(i)==0;i++)
        minBin=i;
    for(i=0;i<histo->GetNbinsX();i++)
        if(histo->GetBinContent(i)>0)
            maxBin = i;
    if(xmax>histo->GetXaxis()->GetBinCenter(maxBin))
        xmax = histo->GetXaxis()->GetBinCenter(maxBin+1);
    Float_t xmin = histo->GetXaxis()->GetBinCenter(minBin-1);
    histo->GetXaxis()->SetRangeUser(xmin,xmax);
    histo->Draw();
    return std::make_pair(xmin,xmax);
}

void HistogrammSaver::OptimizeXRange(TH2* histo){
    histo->Draw();
    TH1F* htemp = (TH1F*)histo->ProjectionX("htemp");
    std::pair<Float_t,Float_t> range = HistogrammSaver::OptimizeXRange(htemp);
    Float_t xmin = range.first;
    Float_t xmax = range.second;
    delete htemp;
    histo->GetXaxis()->SetRangeUser(xmin,xmax);
}


void HistogrammSaver::OptimizeYRange(TH2* histo){
    histo->Draw();
    TH1F* htemp = (TH1F*)histo->ProjectionY("htemp");
    std::pair<Float_t,Float_t> range = HistogrammSaver::OptimizeXRange(htemp);
    Float_t xmin = range.first;
    Float_t xmax = range.second;
    delete htemp;
    histo->GetYaxis()->SetRangeUser(xmin,xmax);
}

void HistogrammSaver::OptimizeXYRange(TH2* histo){
    HistogrammSaver::OptimizeXRange(histo);
    HistogrammSaver::OptimizeYRange(histo);
}

void HistogrammSaver::SaveStack(THStack* stack, TString drawOption,bool bDrawLegend,bool bDrawOnCanvas,TString xTitle,TString yTitle) {
    if (!stack) return;
    TString name = stack->GetName();
    if(name.First('h')==0)
        name.Replace(0,1,'c');
    else
        name.Insert(0,"c_");
    TCanvas *c1 = new TCanvas(name);
    c1->cd();
    c1->SetObjectStat(0);
    stack->SetObjectStat(false);
    c1->Clear();
    stack->Draw();
    if (xTitle=="")
        if(stack->GetXaxis()) xTitle = stack->GetXaxis()->GetTitle();
    if (yTitle=="")
        if(stack->GetYaxis()) yTitle = stack->GetYaxis()->GetTitle();
    if(verbosity>6)cout<<"xTitle: '"<<xTitle<<"'  yTitle: '"<<yTitle<<"'"<<endl;

    if(verbosity>6)cout<<"SaveStack: Xaxis: ";
    stack->GetXaxis()->SetTitle(xTitle);
    if(verbosity>6)cout<<stack->GetXaxis()->GetTitle();
    if(verbosity>6)cout<<endl;
    if(verbosity>6)cout<<"SaveStack: Yaxis: ";
    stack->GetYaxis()->SetTitle(yTitle);
    if(verbosity>6)cout<<stack->GetYaxis()->GetTitle();
    if(verbosity>6)cout<<endl;
    if(drawOption=="")
        stack->Draw();
    else
        stack->Draw(drawOption);
    c1->Draw();
    Double_t xmin = stack->GetXaxis()->GetXmin();
    Double_t xmax = stack->GetXaxis()->GetXmax();
    Double_t ymin = 0;
    Double_t ymax = stack->GetMaximum();
    if(stack->GetYaxis()){
        if(verbosity>6)cout<<"Setting y axis range"<<stack->GetYaxis()<<endl;
        Double_t min = stack->GetYaxis()->GetXmin();
        Double_t max = stack->GetYaxis()->GetXmax();
        if(verbosity>6)cout<<min<<" "<<max<<endl;
        min = TMath::Max(min,stack->GetMinimum(drawOption));
        max = TMath::Max(max,stack->GetMaximum(drawOption));
        if(verbosity>6)cout<<min<<" "<<max<<endl;
        Float_t delta = max-min;

        if(verbosity>6)cout<<delta<<endl;
        stack->GetYaxis()->SetRange(min - .1*delta, max +.2*delta );
        stack->GetYaxis()->SetTitle(yTitle);
        if(verbosity>6)cout<<"set"<<endl;
        c1->Update();
        ymin = min;
        ymax = max;
    }
    if(verbosity>6)   cout<<"draw"<<endl;
    Double_t delta = ymax-ymin;
    ymin = ymin - .1 *(delta);
    ymax = ymax + .3*delta;
    if(bDrawOnCanvas){
        TH1F* frame = c1->DrawFrame(xmin,ymin,xmax,ymax,stack->GetTitle());
        frame->GetXaxis()->SetTitle(xTitle);
        frame->GetYaxis()->SetTitle(yTitle);
        c1->Update();
        gPad->Update();
        if(drawOption=="")
            stack->Draw("same");
        else{
            drawOption=(TString)"same "+drawOption;
            stack->Draw(drawOption);
        }
        c1->Update();
    }
    else{
        if(drawOption=="")
            stack->Draw("");
        else
            stack->Draw(drawOption);
    }
    if (xTitle!="")
        stack->GetXaxis()->SetTitle(xTitle);
    if (yTitle!="")
        stack->GetYaxis()->SetTitle(yTitle);
    c1->Update();
    if (bDrawLegend){
        TLegend* leg = c1->BuildLegend();
        leg->SetX1NDC(0.45);
        leg->SetX2NDC(0.9);
        leg->SetY1NDC(0.74);
        leg->SetY2NDC(0.999);
        //        leg->AddEntry(stack);
        if(leg) leg->SetFillColor(kWhite);
        if(leg) leg->Draw();
    }
    /*
    ostringstream plot_filename;
    plot_filename << plots_path << canvas->GetName()<<".root";
    //*/
    stack->SaveAs(GetROOTFileName(stack->GetName()));
    stack->Write();
    c1->Modified();
    SaveCanvas(c1);
}

TH2D* HistogrammSaver::GetBinContentHisto(TProfile2D* prof) {
    if(!prof) return 0;
    TString name = prof->GetName() + (TString)"_binEntries";
    TH2D* histo = prof->ProjectionXY(name,"B");
    histo->GetZaxis()->SetTitle("number of entries #");
    return histo;
}

void HistogrammSaver::CreateResolutionPlots(vector<TH1*>*vec,TString kind, Int_t subjectDetector, TString appendix){
    UInt_t nBins = 128;
    Float_t minX = -1*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX =settings->GetCellWidth(subjectDetector,2);
    TString name = "hResolutionGoodCells_"+kind+appendix;
    TH1F* hResolutionGoodCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionGoodCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionBadCells_"+kind+appendix;
    TH1F* hResolutionBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionBadCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllCells_"+kind+appendix;
    TH1F* hResolutionAllCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllButBadCells_"+kind+appendix;
    TH1F* hResolutionAllButBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllButBadCells->GetXaxis()->SetTitle("Residual / #mum");
    string plots_path = this->GetPlotsPath();
    HistogrammSaver newHistSaver(settings);
    newHistSaver.SetPlotsPath(plots_path+(string)"/hResolution/"+(string)kind+(string)"/");
    for(UInt_t cell=0;cell< vec->size();cell++){
        TH1F* histo = (TH1F*)vec->at(cell);
        if (!histo)
            continue;
        if (settings->IsGoodCell(3,cell))
            hResolutionGoodCells->Add(histo);
        if (settings->isBadCell(3,cell))
            hResolutionBadCells->Add(histo);
        else
            hResolutionAllButBadCells->Add(histo);
        hResolutionAllCells->Add(histo);
        newHistSaver.SaveHistogram(histo);
        vec->at(cell)= 0;
        delete histo;
    }
    Float_t xmin = -100;
    Float_t xmax = +100;
    TF1* fitX = new TF1("fit","[0]*TMath::Sqrt(TMath::Pi()/2)*[1]*(TMath::Erf(([2]+[3]-x)/TMath::Sqrt(2)/[1])+TMath::Erf(([3]-[2]+x)/TMath::Sqrt(2)/[1]))",xmin,xmax);
    fitX->FixParameter(3,settings->GetCellWidth(subjectDetector,2)/2);//TODO
    fitX->SetParLimits(1,0,40);
    fitX->SetParNames("Integral","sigma of Gaus","position");
    fitX->SetParameter(2,0.);
    fitX->SetParameter(1,10);
    Int_t statOpt = gStyle->GetOptStat();
    gStyle->SetOptStat(1111);
    hResolutionGoodCells->GetYaxis()->SetTitle("number of entries #");
    this->SaveHistogram(hResolutionGoodCells,false,false,true);
    hResolutionBadCells->GetYaxis()->SetTitle("number of entries #");
    this->SaveHistogram(hResolutionBadCells,false,false,true);
    hResolutionAllCells->GetYaxis()->SetTitle("number of entries #");
    this->SaveHistogram(hResolutionAllCells);
    hResolutionAllButBadCells->GetYaxis()->SetTitle("number of entries #");
    this->SaveHistogram(hResolutionAllButBadCells,false,false,true);
    gStyle->SetOptStat(statOpt);
    delete hResolutionGoodCells;
    delete hResolutionBadCells;
    delete hResolutionAllCells;
    delete hResolutionAllButBadCells;
}

TH2D* HistogrammSaver::GetTH2DOfCut(TH2* histo, Float_t cut, EnumDirection dir) {
    if (!histo) return 0;

    TString name = histo->GetName();
    name.Append("_cut");
    cout<<"[HistogrammSaver::GetTH2DOfCut]: "<<name<<" x ";
    if (dir == positive)
        cout<<" > "<<cut<<endl;
    else
        cout<<" < "<<cut<<endl;
    TString title = histo->GetTitle();
    if (dir == positive)
        title.Append(TString::Format(" cut: x > %f",cut));
    else
        title.Append(TString::Format(" cut: x < %f",cut));
    TH2D* h = (TH2D*) histo->Clone(name);
    h->SetTitle(title);
    int n = 0;
    for (UInt_t xbin = 1; xbin <= h->GetNbinsX(); xbin++)
    for (UInt_t ybin = 1; ybin <= h->GetNbinsY(); ybin++){
        Float_t val = h->GetBinContent(xbin,ybin);
        Int_t set_val  = 0;
        if ( (dir == positive && val > cut) ||
             (dir == negative && val < cut) )
            set_val = 1;
        if (set_val)
            n++;
        h->SetBinContent(xbin,ybin,set_val);
        if (dir == positive)
            cout<<TString::Format(" * %2d/%2d: %.3f > %.3f --> %d",xbin,ybin,val,cut,set_val)<<endl;
        else
            cout<<TString::Format(" * %2d/%2d: %.3f < %.3f --> %d",xbin,ybin,val,cut,set_val)<<endl;
    }
    cout<<" Found "<<n<<" out of "<<h->GetNbinsX()*h->GetNbinsY()<<" bins fullfilling the cut."<<endl;
    return h;
}

void HistogrammSaver::DeleteCellHistos(){
    if (hGoodCells){
        delete hGoodCells;
        hGoodCells = 0;
    }
    if (hBadCells){
        delete hBadCells;
        hBadCells = 0;
    }
    if (hAllCells){
        delete hAllCells;
        hAllCells = 0;
    }
    if (hAllButBadCells){
         delete hAllButBadCells;
         hAllButBadCells = 0;
    }
}

void HistogrammSaver::CreateTH2_CellPlots(vector<TH1*> *vec,TString kind,TString prefix, TString appendix){
    if (!vec) return;
     if (vec->size() == 0) return;
     TH2* histo = (TH2*)vec->at(0);
     if (!histo) return;
     cout<<"[HistogrammSaver::CreateTH2_CellPlots] "<<prefix<<" "<<kind<<" "<<appendix<<endl;
     if (kind != "" && !kind.BeginsWith("_"))
         kind.Prepend("_");
     DeleteCellHistos();
     TString name = prefix+"GoodCells"+kind+appendix;
     hGoodCells = (TH2F*)histo->Clone(name);
     hGoodCells->Reset();
     hGoodCells->SetTitle(prefix + " Good Cells "+ kind +appendix);

     name = prefix+"BadCells"+kind+appendix;
     hBadCells =  (TH2F*)hGoodCells->Clone(name);
     hBadCells->Reset();
     hBadCells->SetTitle(prefix + " Bad Cells "+ kind +appendix);

     name = prefix+"AllButBadCells"+kind+appendix;
     hAllButBadCells =  (TH2F*)hGoodCells->Clone(name);
     hAllButBadCells->Reset();
     hAllButBadCells->SetTitle(prefix + " AllButBad Cells "+ kind +appendix);

     name = prefix+"AllCells"+kind+appendix;
     hAllCells =  (TH2F*)hGoodCells->Clone(name);
     hAllCells->Reset();
     hAllCells->SetTitle(prefix + " All Cells "+ kind +appendix);

     string plots_path = this->GetPlotsPath();
     string new_plots_path = plots_path;
     new_plots_path+=(string)prefix;
     TString subfolder = kind(1,kind.Length());
     new_plots_path+="/"+(string)subfolder+"/";
     HistogrammSaver newHistSaver(settings);
     newHistSaver.SetPlotsPath(new_plots_path);
     //TH2* histo;

     for(UInt_t cell=0;cell< vec->size();cell++){
         TRawEventSaver::showStatusBar(cell,vec->size(),1);
         histo = (TH2*)vec->at(cell);
         if (settings->isBadCell(3,cell))
              hBadCells->Add(histo);
         else if (settings->IsGoodCell(3,cell))
             hGoodCells->Add(histo);
         else
             hAllButBadCells->Add(histo);
         newHistSaver.SaveHistogram(histo,true,false);
         vec->at(cell)= 0;
         delete histo;
     }
     cout<<"#"<<flush;
     hAllButBadCells->Add(hGoodCells);
     hAllCells->Add(hAllButBadCells);
     hAllCells->Add(hBadCells);
     cout<<"DONE"<<endl;
     ((TH2F*)hGoodCells)->GetZaxis()->SetTitle("number of entries #");
     this->SaveHistogram((TH2F*)hGoodCells,true,false);
     TH1D* proj = ((TH2F*)hGoodCells)->ProjectionX();
     this->SaveHistogram(proj);
     delete proj;
     ((TH2F*)hBadCells)->GetZaxis()->SetTitle("number of entries #");
     this->SaveHistogram((TH2F*)hBadCells,true,false);
     proj = ((TH2F*)hBadCells)->ProjectionX();
     this->SaveHistogram(proj);
     delete proj;
     ((TH2F*)hAllCells)->GetZaxis()->SetTitle("number of entries #");
     this->SaveHistogram((TH2F*)hAllCells,true,false);
     proj = ((TH2F*)hAllCells)->ProjectionX();
     this->SaveHistogram(proj);
     delete proj;
     ((TH2F*)hAllButBadCells)->GetZaxis()->SetTitle("number of entries #");
     this->SaveHistogram((TH2F*)hAllButBadCells,true,false);
     proj = ((TH2F*)hAllButBadCells)->ProjectionX();
     this->SaveHistogram(proj);
     delete proj;
     cout<<"Test: "<<((TH2F*)hAllButBadCells)->GetName()<<" with "<<((TH2F*)hAllButBadCells)->GetEntries()<<endl;
}
