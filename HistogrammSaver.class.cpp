/*
 * HistogrammSaver.class.cpp
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 */

#include "HistogrammSaver.class.hh"

HistogrammSaver::HistogrammSaver(int verbose) {
	// TODO Auto-generated constructor stub
	verbosity=verbose;
	runNumber=0;
	nEvents=0;
	plots_path=".";
	pt = new TPaveText(0.07,0,0.22,0.10,"NDC");
	UpdatePaveText();
	if(verbosity)cout<<"Created instance of HistogrammSaver"<<endl;
}

HistogrammSaver::~HistogrammSaver() {
	// TODO Auto-generated destructor stub
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
	if(verbosity)cout<<"Set Plotspath: \""<<plots_path<<"\""<<endl;
}

void HistogrammSaver::SetStyle(TStyle newStyle){
	//	currentStyle.TStyle(&newStyle);//.Clone());
	newStyle.cd();
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


