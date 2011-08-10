/*
 * HistogrammSaver.class.cpp
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 */

#include "HistogrammSaver.class.hh"

using namespace std;

HistogrammSaver::HistogrammSaver(int verbosity) {
	sys=NULL;
	pt=NULL;
	this->verbosity=verbosity;
	runNumber=0;
	nEvents=0;
	plots_path=".";
	pt = new TPaveText(0.07,0,0.22,0.10,"NDC");
	UpdatePaveText();
	if(verbosity)cout<<"HistogrammSaver::HistogrammSaver:: get new TSystem"<<endl;
	sys=new TSystem();
	if(verbosity)cout<<"HistogrammSaver::HistogrammSaver:: Set Style"<<endl;
	currentStyle= new TStyle("HistSaverStyle","HistSaverStyle");
	currentStyle->SetPalette(1);
	gROOT->SetStyle("Plain"); //General style (see TStyle)
	gStyle->SetOptStat(111110); //Stat options to be displayed			without under- and overflow use gStyle->SetOptStat(1110);
	gStyle->SetOptFit(1111);  //Fit options to be displayed
	gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadTopMargin(0.15);
	//gStyle->SetTitleColor(19,"");
	gStyle->SetStatH(0.12); //Sets Height of Stats Box
	gStyle->SetStatW(0.15); //Sets Width of Stats Box
	gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
	if(verbosity)cout<<"HistogrammSaver::HistogrammSaver::Created instance of HistogrammSaver"<<endl;
}

HistogrammSaver::~HistogrammSaver() {
	this->pt->Delete();
	currentStyle->Delete();
	sys->Delete();
}

void HistogrammSaver::UpdatePaveText(){
	pt->Clear();
	pt->SetTextSize(0.0250);
	std::ostringstream svnRev_label;
	svnRev_label<<"SVN-Rev: "<<SVN_REV;
	pt->AddText(svnRev_label.str().c_str());
	std::ostringstream run_number_label;
	run_number_label << "Run " <<runNumber;
	pt->AddText(run_number_label.str().c_str());
	std::ostringstream pthresh2;
	pthresh2 << nEvents << " Events in Data Set";
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
	plots_path.assign(path);
	if(verbosity)cout<<"HistogrammSaver::Set Plotspath: \""<<plots_path<<"\""<<endl;
	sys->mkdir(plots_path.c_str());
}

void HistogrammSaver::SetStyle(TStyle newStyle){
	//	currentStyle.TStyle(&newStyle);//.Clone());
	delete currentStyle;
	currentStyle = new TStyle(newStyle);
	currentStyle->cd();
}

/**
 * *********************************************************
 * *********************************************************
 */
void HistogrammSaver::SaveHistogram(TH1F* histo) {
   SaveHistogramPNG(histo);
   SaveHistogramROOT(histo);
}

void HistogrammSaver::SaveHistogram(TH2F* histo) {
   SaveHistogramPNG(histo);
   SaveHistogramROOT(histo);
}

void HistogrammSaver::SaveHistogramPDF(TH1F* histo) {
	TCanvas plots_canvas("plots_canvas","plots_canvas");
	plots_canvas.cd();
	histo->Draw();
	pt->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas.Print(plot_filename.str().c_str());
}

void HistogrammSaver::SaveHistogramPDF(TH2F* histo) {
	TCanvas plots_canvas("plots_canvas","plots_canvas");
	//plots_canvas.cd();
//	SetDuckStyle();
	plots_canvas.cd();
	if(verbosity)cout << "Using SaveHistogrammPDF on TH2F histogram " << histo->GetName() << endl;
	//histo->Draw();
	gStyle->SetTitleFont(42);
	gStyle->SetMarkerSize(0);
	pt->SetTextSize(0.0250);
	pt->SetTextColor(kBlack);
	histo->SetTitleFont(42);
	histo->UseCurrentStyle();
	histo->Draw("colz");
	pt->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas.Print(plot_filename.str().c_str());
	//pt->SetTextSize(runNumber0.1);
}

void HistogrammSaver::SaveHistogramPNG(TH1F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".png";
   plots_canvas.Print(plot_filename.str().c_str());
}

void HistogrammSaver::SaveHistogramROOT(TH1F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   pt->Draw();
   ostringstream plot_filename;
	ostringstream histo_filename;
//	histo_filename << plots_path << histo->GetName() << "_histo.root";
	histo_filename << plots_path << "histograms.root";
   plot_filename << plots_path << histo->GetName() << ".root";
   plots_canvas.Print(plot_filename.str().c_str());
	TFile f(histo_filename.str().c_str(),"UPDATE");
	histo->Write();
}

void HistogrammSaver::SaveHistogramPNG(TH2F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   histo->Draw("colz");
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".png";
   plots_canvas.Print(plot_filename.str().c_str());
}

void HistogrammSaver::SaveHistogramROOT(TH2F* histo) {
   TCanvas plots_canvas("plots_canvas","plots_canvas");
   plots_canvas.cd();
   histo->Draw();
   histo->Draw("colz");
   pt->Draw();
   ostringstream plot_filename;
   plot_filename << plots_path << histo->GetName() << ".root";
   plots_canvas.Print(plot_filename.str().c_str());
}

void HistogrammSaver::SetVerbosity(unsigned int i)
{
	this->verbosity=i;
	if(verbosity) cout<<"HistogrammSaver::Set Verbosity ON"<<endl;
}


void HistogrammSaver::SaveCanvasRoot(TCanvas *canvas, string location, string file_name)
{
    char loc[500];
    memcpy(loc,location.c_str(),strlen(location.c_str())+1);
    char rt[] = ".root";
    char *rtloc = loc;
    //Saving .root file
    strcat(rtloc,file_name.c_str());
    strcat(rtloc,rt);
    char const *rt_file = &rtloc[0];
    TObjArray list(0);
    list.Add(canvas);
    TFile f(rt_file,"recreate");
    list.Write();
    f.Close();
    cout << ".root file was created at: " << rt_file << endl;
}

//void SaveCanvasC(TCanvas *canvas, char* location, char* file_name);
void SaveCanvasC(TCanvas *canvas, string location, string file_name)
{
   char loc[500];
   memcpy(loc,location.c_str(),strlen(location.c_str())+1);
   char cmac[] = ".C";
   char *cmacloc = loc;

   //Saving .C macro
   strcat(cmacloc,file_name.c_str());
   strcat(cmacloc,cmac);
   char const *cmac_file = &cmacloc[0];
   canvas->SaveSource(cmac_file);
   cout << ".c macro was created at: " << cmac_file << endl;
}

void HistogrammSaver::SaveCanvasPNG(TCanvas *canvas, string location, string file_name)
{
   char loc[500];
   memcpy(loc,location.c_str(),strlen(location.c_str())+1);
   char png[] = ".png";
   char *pngloc = loc;

   //Save .png file
   strcat(pngloc,file_name.c_str());
   strcat(pngloc,png);
   char const *png_file = &pngloc[0]; //assigns address of first element of the file string to the char_file pointer
   TImage *img = TImage::Create();
   img->FromPad(canvas);
   img->WriteImage(png_file);
   cout << ".png file was created at: " << png_file << endl;
   delete img;
}

void HistogrammSaver::SetDuckStyle() {
	TStyle* DuckStyle = new TStyle("DuckStyle", "Famous Duck Style");
	DuckStyle->SetOptStat(1110); //Stat options to be displayed
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
	DuckStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
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
	DuckStyle->SetLabelSize(0.025,"XYZ");
	//DuckStyle->SetTitleOffset(1.8, "Y"); // Another way to set the Offset
	//	gStyle->SetTitleOffset(1.2, "X"); // Another way to set the Offset
	DuckStyle->SetTitleOffset(1.2,"X");
	DuckStyle->cd();

	cout << "Using DuckStyle" << endl;
}
