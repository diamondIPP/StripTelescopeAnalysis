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
	//	sys=new TSystem();
	sys=gSystem;
	if(verbosity)cout<<"HistogrammSaver::HistogrammSaver:: Set Style"<<endl;
	currentStyle=gROOT->GetStyle("Plain_RD42");
	if(currentStyle!=0)
		currentStyle->cd();
	else if(gStyle!=0)
		if(!gStyle->IsZombie()){
			if((string)gStyle->GetName()!="Plain_RD42"){
				gROOT->SetStyle("Plain"); //General style (see TStyle)
				//	    gStyle->SetOptStat(221111111); //Stat options to be displayed			without under- and overflow use gStyle->SetOptStat(1110);
				if(gStyle->GetOptStat()!=221111111)
					gStyle->SetOptStat("nemr");//KSiou
				if(gStyle->GetOptFit()!=1111){
					gStyle->SetOptFit(1111);  //Fit options to be displayed
					gStyle->SetStatH(0.12); //Sets Height of Stats Box
					gStyle->SetStatW(0.15); //Sets Width of Stats Box
				}
				if(gStyle->GetPadBottomMargin()!=0.15) gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
				//	    gStyle->SetPadRightMargin(0.15);
				if(gStyle->GetPadTopMargin()!=0.15) gStyle->SetPadTopMargin(0.15);
				//gStyle->SetTitleColor(19,"");
				gStyle->SetPalette(1); //
				gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
				currentStyle= (TStyle*)gStyle->Clone("Plain_RD42");
				currentStyle->SetPalette(1);
				currentStyle2D= (TStyle*)currentStyle->Clone("Plain_RD42_2D");
				currentStyle2D->SetOptStat("ne");
				currentStyle2D ->SetPalette(1);
				currentStyle->cd();
				gROOT->SetStyle("Plain_RD42");
			}
		}

	gStyle->SetPalette(1); //
	if(verbosity)cout<<"HistogrammSaver::HistogrammSaver::Created instance of HistogrammSaver"<<endl;
	gErrorIgnoreLevel=3001;

}

HistogrammSaver::~HistogrammSaver() {

	//	TString string1 = sys->GetFromPipe(".! mkdir root-Files");
	//	cout<<string1<<endl;
	stringstream test;
	test<< "mv -f "<<plots_path<<"/*.root "<<plots_path<<"/root/";
	//cout<<"\""<<test.str()<<"\""<<endl;
	system(test.str().c_str());//t.str();//<<"\""<<endl;
	//	string1 = sys->GetFromPipe(".!mv -v *.root root-Files");
	//	cout<<string1<<endl;
	this->pt->Delete();
}


void HistogrammSaver::SaveTwoHistos(std::string canvasName, TH1F *histo1, TH1F *histo2,double refactorSecond, UInt_t verbosity)
{
	if (verbosity>2) cout<<"Save2Histos: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
	TCanvas *c1 = new TCanvas(canvasName.c_str(),canvasName.c_str());
	c1->cd();
	Float_t min1 = histo1->GetMinimum();
	Float_t min2 = histo2->GetMinimum();
	Float_t min = TMath::Min(min1,min2);
	Float_t max1 = histo1->GetMaximum();
	Float_t max2 = histo2->GetMaximum();
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
	if(histo2->GetMaximum()*refactorSecond>histo1->GetMaximum())
		refactorSecond=histo2->GetMaximum()/histo1->GetMaximum()*0.5;
	if(refactorSecond!=1)histo2->Scale(refactorSecond);
	if (verbosity>2) cout<<"min: "<<min<<" max: "<<max;
	if (verbosity>2) cout<<" refactorSecond:"<<refactorSecond<<"\thisto1:"<<histo1->GetMaximum()<<"\thisto2:"<<histo2->GetMaximum()<<flush;
	if (verbosity>2) cout<<endl<<"Nhisto1: "<<histo1->GetEntries()<<" Nhisto2:"<<histo2->GetEntries()<<flush;
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
	TLegend *leg =new TLegend(0.1,0.75,0.48,0.9);
	leg->SetFillColor(kWhite);
	leg->SetHeader("Legend");
	leg->AddEntry(histo1,histo1->GetName());
	leg->AddEntry(histo2,histo2->GetName());
	leg->Draw("same");
	TPaveText* pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvasName.c_str()));
	pt2->Draw("same");
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
	int isNotCreated=sys->mkdir(plots_path.c_str(),true);
	if (isNotCreated!=0){
		//		cout<<"***************************************************\n";
		//		cout<<"********** Directory not created ******************\n";
		//		cout<<"***************************************************\n";
		cout<<plots_path<<endl;
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
	//	else cout<<"Verzeichnis konnte nicht angelegt werden..."<<endl;

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
void HistogrammSaver::SaveHistogram(TH1* histo, bool fitGauss,bool adjustRange) {
	if(histo->GetEntries()==0)return;
	if(adjustRange){
		int binxMin=0;
		for(binxMin=0;binxMin<histo->GetNbinsX();binxMin++)if(histo->GetBinContent(binxMin))break;
		int binxMax;
		for(binxMax=histo->GetNbinsX();binxMax>0;binxMax--)if(histo->GetBinContent(binxMax))break;
		histo->GetXaxis()->SetRangeUser(histo->GetBinLowEdge(binxMin),histo->GetBinLowEdge(binxMax+1));
	}
	//create PNG
	if (fitGauss) SaveHistogramFitGaussPNG(histo);
	else SaveHistogramPNG(histo);
	//create ROOT
	SaveHistogramROOT(histo);
}
void HistogrammSaver::SaveHistogramWithFit(TH1F* histo,TF1* fit, UInt_t verbosity){
	if(histo==0)return;
	if(histo->GetEntries()==0)return;
	if(fit==0) SaveHistogram(histo);
	if (verbosity>0) cout<<"Save Histogram With Fit:"<<histo->GetTitle()<<endl;
	TCanvas *plots_canvas =  new TCanvas( TString::Format("c_%s", histo->GetName() ) , TString::Format("c_%s", histo->GetName() ) );
	plots_canvas->Clear();
	plots_canvas->cd();
	TH1F *htemp = (TH1F*)histo->Clone();
	TF1* fittemp = (TF1*)fit->Clone();
	TPaveText * pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));

	htemp->Draw();
	fittemp->SetLineColor(kRed);
	fittemp->Draw("same");
	pt2->Draw();
	ostringstream plot_filename;
	ostringstream histo_filename;
	histo_filename << plots_path << "histograms.root";
	plot_filename << plots_path << histo->GetName() << ".root";
	plots_canvas->Print(plot_filename.str().c_str());
	TFile f(histo_filename.str().c_str(),"UPDATE");
	f.cd();
	((TH1F*)histo->Clone())->Write();
	((TF1*)fit->Clone())->Write();
	plots_canvas->Write();
	plot_filename.clear();
	plot_filename.str("");
	plot_filename.clear();
	plot_filename << plots_path << histo->GetName() << ".png";
	plots_canvas->Print(plot_filename.str().c_str());
	f.Close();
//	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramWithCutLine(TH1F *histo,Float_t cutValue){
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
void HistogrammSaver::SaveHistogramLogZ(TH2F* histo){
	TString canvasName = "c_";
	canvasName +=histo->GetName();
	TCanvas *c1 = new TCanvas(canvasName,canvasName);
	c1->cd();
	c1->SetLogz();
	TH2F* htemp = (TH2F*) histo->Clone();
	htemp->Draw("colz");
	this->SaveCanvas(c1);
	delete htemp;
	delete c1;
}

void HistogrammSaver::SaveHistogram(TH2F* histo, bool drawStatBox) {
	if (!histo)return;
	if(histo->GetEntries()==0)return;
	if (!drawStatBox)
			histo->SetStats(false);
//	histo->SetStats(false);
	SaveHistogramPNG(histo);
	SaveHistogramROOT(histo);
}

void HistogrammSaver::SaveCanvas(TCanvas *canvas)
{
	if(canvas==0)
		return;
	SaveCanvasPNG(canvas);
	SaveCanvasROOT(canvas);
}
void HistogrammSaver::SaveGraph(TGraph* graph,std::string name,std::string option){
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
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas->Print(plot_filename.str().c_str());
//	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramPDF(TH2F* histo) {
	if(histo->GetEntries()==0)return;
	TCanvas *plots_canvas = new TCanvas(TString::Format("cPdf_%s",histo->GetName()),TString::Format("c_%s",histo->GetName()));
	plots_canvas->cd();
	//plots_canvas.cd();
	//	SetDuckStyle();
	if(verbosity)cout << "Using SaveHistogrammPDF on TH2F histogram " << histo->GetName() << endl;
	//histo->Draw();
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	gStyle->SetTitleFont(42);
	gStyle->SetMarkerSize(0);
	pt2->SetTextSize(0.0250);
	pt2->SetTextColor(kBlack);
	histo->SetTitleFont(42);
	histo->UseCurrentStyle();
	histo->Draw("colz");
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".pdf";
	plots_canvas->Print(plot_filename.str().c_str());
	if(plots_canvas)delete plots_canvas;
	//pt->SetTextSize(runNumber0.1);
}

void HistogrammSaver::SaveHistogramPNG(TH1* histo) {
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

	htemp->SetMinimum(0.);
	htemp->Draw();
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".png";
	plots_canvas->Print(plot_filename.str().c_str());
//	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveCanvasROOT(TCanvas *canvas)
{
	if(!canvas)
		return;
	ostringstream plot_filename;
	plot_filename << plots_path << canvas->GetName()<<".root";
	TCanvas* plots_canvas=(TCanvas*)canvas->Clone();
	plots_canvas->cd();

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",canvas->GetName()));
	pt2->Draw();

	TFile f(plot_filename.str().c_str(),"UPDATE");
	canvas->Write();
}

void HistogrammSaver::SaveCanvasPNG(TCanvas *canvas)
{
	if(canvas==0)
		return;
	canvas->cd();
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",canvas->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << canvas->GetName()<<".png";
	canvas->Print(plot_filename.str().c_str());
}

void HistogrammSaver::SaveGraphPNG(TGraph* graph,string name,string option){
	if(!graph)
		return;
	if(graph->GetN()==0)return;
	TCanvas plots_canvas(TString::Format("c_%s",name.c_str()),TString::Format("c_%s",name.c_str()));
	plots_canvas.cd();
	graph->Draw(option.c_str());

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",graph->GetName()));
	pt2->Draw();

	ostringstream plot_filename;
	plot_filename << plots_path << name << ".png";
	plots_canvas.Print(plot_filename.str().c_str());
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
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",htemp->GetName()));
	pt2->Draw();

	ostringstream plot_filename;
	plot_filename << plots_path << htemp->GetName() << ".png";
	plots_canvas.Print(plot_filename.str().c_str());
}

void HistogrammSaver::SaveHistogramROOT(TH1* htemp) {
	if(!htemp)return;
//	if(htemp->GetEntries()==0)return;

	ostringstream plots_filename;
	ostringstream histo_filename;
	plots_filename << plots_path<<"/" << htemp->GetName() << ".root";
	histo_filename << plots_path << "histograms.root";
	TCanvas *plots_canvas =  new TCanvas(TString::Format("c_%s", htemp->GetName()), TString::Format("c_%s", htemp->GetName()));
	plots_canvas->cd();
	TH1* histo = (TH1*)htemp->Clone();

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("ptRoot_%s",htemp->GetName()));

	plots_canvas->Clear();
	plots_canvas->cd();
	histo->Draw();
	pt2->Draw();
	plots_canvas->Draw();
	histo->Draw();

	//write to own root File
	plots_canvas->Write(plots_filename.str().c_str());
	plots_canvas->Write(plots_filename.str().c_str());
	TFile *f = new TFile(histo_filename.str().c_str(),"UPDATE");
	TCanvas *plots_canvas2 = (TCanvas*) plots_canvas->Clone(TString::Format("cc_%s",htemp->GetName()));
	//add to histograms.root
	f->cd();
	plots_canvas2->Write();
	f->Close();
//	if(plots_canvas)delete plots_canvas;

}

void HistogrammSaver::SaveHistogramPNG(TH2F* histo) {

	if(!histo){
		cerr<<"HistogrammSaver::SaveHistogramPNG(TH2F*), histogram ==0"<<endl;
				return;
	}
	if(histo->GetEntries()==0)return;
	gROOT->SetStyle("Plain_RD42_2D");
	gROOT->ForceStyle(true);
	//currentStyle2D->cd();
	//	cout<<"Save: \""<<histo->GetName()<<"\""<<endl;
//	cout<<"\tgStyle2D: "<<gStyle<<"-->"<<flush;
	TCanvas *plots_canvas =  new TCanvas(TString::Format("cPng_%s", histo->GetName()), TString::Format("c_%s", histo->GetName()));
	plots_canvas->Clear();
	plots_canvas->cd();
	TH2F* htemp = (TH2F*)histo->Clone();
	HistogrammSaver::OptimizeXYRange(htemp);
	htemp->Draw("colz");

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("ptPng_%s",histo->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".png";
	plots_canvas->Print(plot_filename.str().c_str());
	//currentStyle->cd();
	gROOT->SetStyle("Plain_RD42");
	cout<<"gStyle: "<<gStyle<<endl;
	//currentStyle->cd();
	//	if(plots_canvas)delete plots_canvas;
}

void HistogrammSaver::SaveHistogramROOT(TH2F* histo) {
	if(!histo){
		cerr<<"HistogrammSaver::SaveHistogramROOT(TH2F*) histogram == 0"<<endl;
		return;
	}
	if(histo->GetEntries()==0)return;
	TCanvas *plots_canvas =  new TCanvas(TString::Format("cRoot_%s", histo->GetName()), TString::Format("c_%s", histo->GetName()));
	plots_canvas->Clear();

	plots_canvas->cd();
	TH2F* htemp = (TH2F*)histo->Clone();
	if(htemp==0)
		return;
	htemp->Draw();
	HistogrammSaver::OptimizeXYRange(htemp);
	htemp->Draw("colz");

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".root";
	plots_canvas->Print(plot_filename.str().c_str());

	stringstream histo_filename;
	histo_filename << plots_path << "histograms.root";
	TFile *f = new TFile(histo_filename.str().c_str(),"UPDATE");
	f->cd();
	plots_canvas->Write();
	f->Close();
//	if (plots_canvas) delete plots_canvas;
}


void HistogrammSaver::SaveHistogramROOT(TH3F* histo){
	if(!histo){
			cerr<<"HistogrammSaver::SaveHistogramROOT(TH2F*) histogram == 0"<<endl;
			return;
		}
	if(histo->GetEntries()==0)return;
	TH3F* htemp = (TH3F*)histo->Clone();
	if(htemp==0)
		return;
	htemp->Draw();
	TString name = histo->GetName();
//	cout<<"Histo Name: "<<name<<endl;
//	cout<<"plots path: "<<plots_path<<endl;
	string fileName = plots_path.c_str();
	fileName.append(name);
	fileName.append(".root");
//	cout<<"pwd: "<<gSystem->pwd()<<endl;;
//	cout<<"Write 3d histo: " <<fileName<<flush;char t; cin>>t;
	htemp->Write(fileName.c_str());
	htemp->Write();
	stringstream plot_filename;
	plot_filename<< plots_path << histo->GetName() << ".root";
	htemp->Print(plot_filename.str().c_str());
	stringstream histo_filename;
	histo_filename << plots_path << "histograms.root";
	TFile *f = new TFile(histo_filename.str().c_str(),"UPDATE");
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

	TCanvas *plots_canvas =  new TCanvas(TString::Format("c_%s", name.c_str()), TString::Format("c_%s", name.c_str()));
	plots_canvas->Clear();

	plots_canvas->cd();
	TGraph* gTemp = (TGraph*)graph->Clone();
	gTemp->Draw(option.c_str());

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",graph->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << name<< ".root";
	plots_canvas->Print(plot_filename.str().c_str());
//	if(plots_canvas)	delete plots_canvas;
}

void HistogrammSaver::SetVerbosity(unsigned int i)
{
	this->verbosity=i;
	if(verbosity) cout<<"HistogrammSaver::Set Verbosity ON"<<endl;
}


void HistogrammSaver::SaveCanvasRoot(TCanvas *canvas, string location, string file_name)
{
	if(canvas==0)
		return;
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
	if(canvas==0 )return;
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
	if(canvas==0)return;
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


/**
 * @brief creates a scatter histogram with posX_vs_posY as an input
 *
 * @return TH2F histogram
 */
TH2F* HistogrammSaver::CreateScatterHisto(std::string name, std::vector<Float_t> posY, std::vector<Float_t> posX, UInt_t nBins)
{
	Float_t factor = 0.05;//5% bigger INtervall...
	if(posX.size()!=posY.size()||posX.size()==0) {
		cerr<<"ERROR HistogrammSaver::CreateScatterHisto vectors have different size "<<posX.size()<<" "<<posY.size()<<" "<<name<<endl;
		return new TH2F();
	}
	Float_t maxX = posX.at(0);
	Float_t maxY = posY.at(0);
	Float_t minX = posY.at(0);
	Float_t minY = posY.at(0);
	for(UInt_t i=0;i<posX.size();i++){
		if(posX.at(i)>maxX)maxX=posX.at(i);
		else if(posX.at(i)<minX)minX=posX.at(i);
		if(posY.at(i)>maxY)maxY=posY.at(i);
		else if(posY.at(i)<minY)minY=posY.at(i);
	}
	//cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
	Float_t deltaX=maxX-minX;
	Float_t deltaY=maxY-minY;
	TH2F* histo = new TH2F(name.c_str(),name.c_str(),nBins,minX-factor*deltaX,maxX+factor*deltaX,nBins,minY-factor*deltaY,maxY+factor*deltaY);
	for(UInt_t i=0;i<posX.size();i++)
		histo->Fill(posX.at(i),posY.at(i));
	histo->GetXaxis()->SetTitle("X-Position");
	histo->GetYaxis()->SetTitle("Y-Position");

	//set xrange
	TH1D* hProj=histo->ProjectionX();
	int binxMin=0;
	for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
	int binxMax;
	for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
	histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
	delete hProj;

	//set yRange
	hProj=histo->ProjectionY();
	int binyMin=0;
	for(binyMin=0;binyMin<hProj->GetNbinsX();binyMin++)if(hProj->GetBinContent(binyMin))break;
	int binyMax;
	for(binyMax=hProj->GetNbinsX();binyMax>0;binyMax--)if(hProj->GetBinContent(binyMax))break;

	histo->GetYaxis()->SetRangeUser(hProj->GetBinLowEdge(binyMin-1),hProj->GetBinLowEdge(binyMax+1));
	delete hProj;

	return histo;
}

TGraph HistogrammSaver::CreateDipendencyGraph(std::string name, std::vector<Float_t> vecY, std::vector<Float_t> vecX)
{
	if(vecY.size()!=vecX.size()||vecX.size()==0) {
		cerr<<"ERROR HistogrammSaver::CreateDipendencyGraph vectors have different size "<<vecY.size()<<" "<<vecX.size()<<": "<<name<<endl;
		return TGraph();
	}
	//cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
	TGraph hGraph = TGraph(vecY.size(),&vecX.at(0),&vecY.at(0));
	hGraph.GetXaxis()->SetName("PredictedPosition");
	hGraph.GetYaxis()->SetName("Delta");
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
TH1F* HistogrammSaver::CreateDistributionHisto(std::string name, std::vector<Float_t> vec, UInt_t nBins,EnumAxisRange range,Float_t xmin,Float_t xmax)
{
	Float_t factor = 0.05;//5% bigger INtervall...
	if(vec.size()==0)
		return new TH1F(name.c_str(),name.c_str(),nBins,0.,1.);
	Float_t max = vec.at(0);
	Float_t min = vec.at(0);
	cout<<"Create Histo: "<<range<<" "<<flush;
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
		cout<<" maxWidth "<<min <<"-"<<max<<endl;
	}
	else if(range==fiveSigma||range==threeSigma){
		Float_t  mean2 =0;
		Float_t sigma2 = 0;
		for(UInt_t i=0;i<vec.size();i++){
			mean2+=vec.at(i);
			sigma2+=vec.at(i)*vec.at(i);
		}
		mean2/=(Float_t)vec.size();
		sigma2/=(Float_t)vec.size();

		Float_t mean=0;
		Float_t sigma=0;
		UInt_t nEvents=0;
		for(UInt_t i=0;i<vec.size();i++){
			if( (vec.at(i)-mean2)/sigma2<3.){
			mean+=vec.at(i);
			sigma+=vec.at(i)*vec.at(i);
			nEvents++;
			}
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
		//cout<<" nSigma:"<<nSigma<<" "<<mean<<"+/-"<<sigma<<" ==> "<<min <<"-"<<max<<endl;
	}
	else if(range==positiveArea){
		min=0;
		for(UInt_t i=0;i<vec.size();i++)
			if (max<vec.at(i))max=vec.at(i);
		max*=(1+factor);
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
	}
	else if(range==manual){
		max =xmax;
		min=xmin;
	}

	TH1F* histo = new TH1F(name.c_str(),name.c_str(),nBins,min,max);
	for(UInt_t i=0;i<vec.size();i++){
		histo->Fill(vec.at(i));
	}
	int ntries=0;
	while ((histo->GetBinContent(histo->GetMaximumBin())/histo->GetEntries())<0.05&&ntries<3){//todo change hardcoding
		histo->Rebin();ntries++;
	}
	histo->GetXaxis()->SetRangeUser(min,max);
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

void HistogrammSaver::OptimizeXRange(TH2F* histo){
	histo->Draw();
	TH1F* htemp = (TH1F*)histo->ProjectionX("htemp");
	std::pair<Float_t,Float_t> range = HistogrammSaver::OptimizeXRange(htemp);
	Float_t xmin = range.first;
	Float_t xmax = range.second;
	delete htemp;
	histo->GetXaxis()->SetRangeUser(xmin,xmax);
}


void HistogrammSaver::OptimizeYRange(TH2F* histo){
	histo->Draw();
	TH1F* htemp = (TH1F*)histo->ProjectionY("htemp");
	std::pair<Float_t,Float_t> range = HistogrammSaver::OptimizeXRange(htemp);
	Float_t xmin = range.first;
	Float_t xmax = range.second;
	delete htemp;
	histo->GetYaxis()->SetRangeUser(xmin,xmax);
}

void HistogrammSaver::OptimizeXYRange(TH2F* histo){
	HistogrammSaver::OptimizeXRange(histo);
	HistogrammSaver::OptimizeYRange(histo);
}
