/*
 * HistogrammSaver.class.cpp
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
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
	InitializeGridReferenceDetSpace();

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

	for(UInt_t i=0;i<settings->getNRows3d();i++){
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
	hGridReferenceDetSpace->GetXaxis()->SetTitle("Row of Cell");
	hGridReferenceDetSpace->GetYaxis()->SetTitle("Column of Cell");
	hGridReferenceCellSpace->SetStats(kFALSE);
	hGridReferenceCellSpace->SetTickLength(0.0, "X");
	hGridReferenceCellSpace->SetTickLength(0.0, "Y");
}


void HistogrammSaver::SaveTwoHistosNormalized(TString canvasName, TH1 *histo1, TH1 *histo2,double refactorSecond, UInt_t verbosity){
	cout<<"Save2HistosNormalized: "<<histo1<<" "<<histo2<<endl;
	if(!histo1&&!histo2)return;
	if(!histo1||!histo2){
		if (histo1) SaveHistogram(histo1);
		else SaveHistogram(histo2);
		return;
	}
	cout<<"Save2HistosNormalized: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
	TCanvas *c1 = new TCanvas(canvasName,canvasName);
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
	if (verbosity>2) cout<<"min: "<<min<<" max: "<<max;
	if (verbosity>2) cout<<" refactorSecond:"<<refactorSecond<<"\thisto1:"<<max1<<"\thisto2:"<<max2<<flush;
	if (verbosity>2) cout<<endl<<"Nhisto1: "<<histo1->GetEntries()<<" Nhisto2:"<<histo2->GetEntries()<<flush;
	histo1->SetStats(false);
	histo2->SetStats(false);
	TH1F* histo1Normalized;
	TH1F* histo2Normalized;
	if(max1>max2){
		if (verbosity>2) cout<<"\tdraw1-"<<flush;
		histo1Normalized = (TH1F*)((TH1F*)(histo1->Clone()))->DrawNormalized("");
		if(histo1Normalized)
			histo1Normalized->GetYaxis()->SetRangeUser(min,max);
		else
			histo1->GetYaxis()->SetRangeUser(min,max);
		if (verbosity>2) cout<<"draw2 "<<flush;
		histo2Normalized = (TH1F*)((TH1F*)(histo2->Clone()))->DrawNormalized("");
		//		histo2->GetYaxis()->SetRangeUser(min,max);
	}
	else{
		if (verbosity>2) cout<<"\tdraw2-"<<flush;
		histo2Normalized = (TH1F*)((TH1F*)(histo2->Clone()))->DrawNormalized("");
		if (histo2Normalized)
			histo2Normalized->GetYaxis()->SetRangeUser(min,max);
		else
			histo2->GetYaxis()->SetRangeUser(min,max);
		if (verbosity>2) cout<<"draw1 "<<flush;
		histo1Normalized =  (TH1F*)((TH1F*)(histo1->Clone()))->DrawNormalized("same");
		//		histo1->GetYaxis()->SetRangeUser(min,max);
	}
	c1->Update();
//	TVirtualPad *pad =c1->GetPad(0);
	if (verbosity>2) cout<<"MIN: "<<min<<"-->";
	min=(double)(min/refactorSecond);
	if (verbosity>2) cout<<min<<"\t\tMAX: "<<max<<"--->";
	max = (double)(max/refactorSecond);
	if (verbosity>2) cout<<max<<endl;
	c1->Update();
	TLegend *leg =new TLegend(0.52,0.75,0.9,0.9);
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
	TPaveText* pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvasName.Data()));
	pt2->Draw("same");
	c1->Update();
	SaveCanvas(c1);
}
void HistogrammSaver::SaveTwoHistos(TString canvasName, TH1 *histo1, TH1 *histo2,double refactorSecond, UInt_t verbosity)
{
	cout<<"Save2Histos: "<<histo1<<" "<<histo2<<endl;
	if(!histo1&&!histo2)return;
	if(!histo1||!histo2){
		if (histo1) SaveHistogram(histo1);
		else SaveHistogram(histo2);
		return;
	}
	if(histo1->GetLineColor() == histo2->GetLineColor())
		histo2->SetLineColor(histo1->GetLineColor()+1);
	cout<<"Save2Histos: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
	if (verbosity>2) cout<<"Save2Histos: "<<histo1->GetName()<<" "<<histo2->GetName()<<" to "<<canvasName<<endl;
	TCanvas *c1 = new TCanvas(canvasName,canvasName);
	c1->cd();
	c1->SetObjectStat(false);
	Float_t min1 = histo1->GetBinContent(histo1->GetMinimumBin());//GetMinimum();
	Float_t min2 = histo2->GetBinContent(histo2->GetMinimumBin());
	Float_t min = TMath::Min(min1,min2);
	Float_t max1 =  histo1->GetBinContent(histo1->GetMaximumBin());//GetMinimum();
	Float_t max2 = histo2->GetBinContent(histo2->GetMaximumBin());
	//	Float_t range1 = max1-min1;
	//	Float_t range2 = max2-min2;
	Float_t factor = 1.1;
	Float_t max = TMath::Max(max1,max2);
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
	TPaveText* pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",canvasName.Data()));
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

	TF1* fit = 0;
	fit = (TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("Fitfcn_%s",histo->GetName()));
	if(!fit)
		(TF1*)histo->GetListOfFunctions()->FindObject(TString::Format("fLangauFixedNoise_%s",histo->GetName()));

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
	TPaveText* pt = updateMean(histo,minX,maxX);
	if(gSigma>0&&pt)
		pt->AddText(TString::Format("GSigma  =   %.1f",gSigma));
	maxX = histo->GetBinLowEdge(histo->GetNbinsX());
	TF1* fMeanCalculationArea = new TF1("fMeanCalculationArea","pol0",minX,maxX);
	fMeanCalculationArea->SetLineColor(kGreen);
	fMeanCalculationArea->FixParameter(0,0);
	fMeanCalculationArea->SetLineWidth(5);
	histo->Fit(fMeanCalculationArea,"Q+","",minX,maxX);
	return pt;
}


TCanvas* HistogrammSaver::DrawHistogramWithCellGrid(TH2* histo,TH2* histo2){
	TString name = histo->GetName();
		if (name.BeginsWith("h"))
			name.Replace(0,1,"c");
		else
			name.Insert(0,"c_");
		TCanvas* c1 = new TCanvas(name,name);
		c1->cd();
		hGridReferenceDetSpace->SetTitle(histo->GetTitle());		//Set title to require
		hGridReferenceDetSpace->Draw("COL");
	    if (histo)
		    histo->Draw("sameCOLZAH");
	//	TLegend* leg = 0;
		if (histo2){
			histo2->Draw("sameTEXTAH");
	//		if(histo2!=histo){
	//			leg = c1->BuildLegend();
	//			leg->Clear();
	//			leg->AddEntry(histo);
	//			leg->AddEntry(histo2);
	//		}
		}
		//hGridReference->Draw("COL");
		settings->DrawMetallisationGrid(c1, 3);
	//	cout<<c1->GetName()<<endl;
	//	if (leg)
	//		leg->Draw();
	return c1;
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
		if(!settings->isValidCellNo((*quarter).first)){
			cerr<< "Invalid Cell No: " << (*quarter).first<<endl;
			continue;
		}

		int column = settings->getColumnOfCell((*quarter).first);
		int row = settings->getRowOfCell((*quarter).first);
		float xLow = xStart + (column+.5*((*quarter).second/2))*cellwidth;
		float yLow = yStart + (row+.5*((*quarter).second%2))*cellheight;
		float xHigh = xLow+cellwidth/2;
		float yHigh = yLow+cellheight/2;
		TString name = c1->GetName();
		name.Append(TString::Format("_FailedQuarter_%dOf%d",i,(int)failedQuarters.size()));
//		cout<<" DRAW: "<< name<<endl;
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
	}
}

TH2D* HistogrammSaver::GetHistoBinedInQuarters(TString name) {
	return GetHistoBinedInCells(name,2);
}

TH2D* HistogrammSaver::GetHistoBinedInCells(TString name, Int_t binsPerCellAxis) {
	cout<<"create "<<name<<endl;
	TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
	TH2D* histo = new TH2D(name,name,
			settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
			settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh());
	return histo;
}

TH3D* HistogrammSaver::Get3dHistoBinedInCells(TString name, UInt_t binsz,
		Float_t minz, Float_t maxz, Int_t binsPerCellAxis) {
	cout<<"create "<<name<<endl;
		TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
		TH3D* histo = new TH3D(name,name,
				settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
				settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh(),
				binsz, minz,maxz);
		return histo;
}

TProfile2D* HistogrammSaver::GetProfile2dBinedInCells(TString name,
		Int_t binsPerCellAxis) {
	cout<<"create "<<name<<endl;
	TFiducialCut* diaMetFidCut = settings->get3dMetallisationFidCuts()->getFidCut(3);
	TProfile2D* histo = new TProfile2D(name,name,
			settings->getNColumns3d()*binsPerCellAxis,diaMetFidCut->GetXLow(),diaMetFidCut->GetXHigh(),
			settings->getNRows3d()*binsPerCellAxis,diaMetFidCut->GetYLow(),diaMetFidCut->GetYHigh());
	return histo;
}

TProfile2D* HistogrammSaver::CreateProfile2D(std::string name,
		std::vector<Float_t> posX, std::vector<Float_t> posY,
		std::vector<Float_t> posZ, UInt_t nBinsX, UInt_t nBinsY,
		Float_t minRangeX, Float_t maxRangeX, Float_t minRangeY,
		Float_t maxRangeY, Float_t minRangeZ, Float_t maxRangeZ,
		Float_t factor) {
    return new TProfile2D();//todp
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
	//	else cout<<"Verzeichnis konnte nicht angelegt werden..."<<endl;

}

void HistogrammSaver::SetStyle(TStyle newStyle){
	//	currentStyle.TStyle(&newStyle);//.Clone());
	delete currentStyle;
	currentStyle = new TStyle(newStyle);
	currentStyle->cd();
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
	TString name = TString::Format("c_%s",histo->GetName());
	TCanvas *c1 = new TCanvas(name,name);
	c1->cd();
	histo->Draw();
	stats->Draw();
//	cout<<"Saving: "<<c1->GetName()<<endl;
	SaveCanvas(c1);
	//create ROOT
	SaveHistogramROOT(histo);
}

/**
 * *********************************************************
 * *********************************************************
 */
void HistogrammSaver::SaveHistogram(TH1* histo, bool fitGauss,bool adjustRange,bool drawStatBox) {
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

void HistogrammSaver::SaveHistogramWithFit(TH1F* histo,TF1* fit,Float_t xmin,Float_t xmax, UInt_t verbosity){
	if(histo==0)return;
	if(histo->GetEntries()==0)return;
	if(fit==0) SaveHistogram(histo);
	if (verbosity>0) cout<<"Save Histogram With Fit:"<<histo->GetTitle()<<endl;
	TCanvas *plots_canvas =  new TCanvas( TString::Format("c_%s", histo->GetName() ) , TString::Format("c_%s", histo->GetName() ) );
	plots_canvas->Clear();
	plots_canvas->cd();
//	TH1F *htemp = (TH1F*)histo->Clone();
	TPaveText * pt2 = (TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	if(verbosity){
		cout<<"Fitting: "<<fit->GetName()<<" "<<xmin<<" - "<<xmax<<endl;
		histo->Fit(fit,"","",xmin,xmax);
	}
	else
		histo->Fit(fit,"Q+","",xmin,xmax);
	histo->Draw();
//	TF1* fittemp = (TF1*)fit->Clone();
//	fittemp->SetLineColor(kRed);
	fit->SetLineStyle(3);
	fit->Draw("same");
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

void HistogrammSaver::SaveHistogram(TH2* histo, bool drawStatBox) {
	if (!histo)return;
	if(histo->GetEntries()==0)return;
	if (!drawStatBox)
			histo->SetStats(false);
//	histo->SetStats(false);
	SaveHistogramPNG(histo);
	SaveHistogramROOT(histo);
}

void HistogrammSaver:: Save1DProfileYWithFitAndInfluence(TH2* histo, TString function){
	TString name = "fit_" + (TString)histo->GetName();
	TF1* fit = new TF1(name,function);
	return Save1DProfileYWithFitAndInfluence(histo,fit);
}

void HistogrammSaver::Save1DProfileYWithFitAndInfluence(TH2* htemp,TF1* fit){
	if(!fit)
		fit = new TF1("fit","pol1");
	TProfile *prof=0;
	if(!htemp)
		return;
	prof = htemp->ProfileY();
	if(prof){
		TCanvas* c1 = new TCanvas( (TString)("c_"+(TString)prof->GetName()) );
		prof->Fit(fit);
		Float_t xmin = prof->GetXaxis()->GetXmin();
		Float_t xmax = prof->GetXaxis()->GetXmax();
		Float_t ymin = fit->GetMinimum(xmin,xmax);
		Float_t ymax = fit->GetMaximum(xmin,xmax);
		TPaveText *text = new TPaveText(.2,.2,.5,.3,"brNDC");
		text->SetFillColor(0);
		text->AddText(TString::Format("relative Influence: #frac{#Delta_{y}}{y_{max}} = %2.2f %%",(ymax-ymin)/ymax*100));
		text->Draw("same");
		SaveCanvas(c1);
		delete prof;
		prof = 0;
	}

}

void HistogrammSaver::SaveCanvas(TCanvas *canvas)
{
	if(canvas==0)
		return;
	SaveCanvasPNG(canvas);
	SaveCanvasROOT(canvas);
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
	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("pt_%s",histo->GetName()));
	pt2->Draw();
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

void HistogrammSaver::SaveHistogramPNG(TH2* histo) {

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
	HistogrammSaver::OptimizeXYRange(htemp);
	htemp->Draw("colz");

	TPaveText *pt2=(TPaveText*)pt->Clone(TString::Format("ptPng_%s",histo->GetName()));
	pt2->Draw();
	ostringstream plot_filename;
	plot_filename << plots_path << histo->GetName() << ".png";
	plots_canvas->Print(plot_filename.str().c_str());
//	gROOT->SetStyle("Plain_RD42");
	//currentStyle->cd();
	//	if(plots_canvas)delete plots_canvas;
	if (htemp) delete htemp;
	if (plots_canvas) delete plots_canvas;
}

void HistogrammSaver::SaveHistogramROOT(TH2* histo) {
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
	TString name = histo->GetName();
	string fileName = plots_path.c_str();
	fileName.append(name);
	fileName.append(".root");
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
		UInt_t nBinsX, UInt_t nBinsY, Float_t minRangeX,Float_t maxRangeX, Float_t minRangeY, Float_t maxRangeY,Float_t factor)
{
//	Float_t factor = 0.05;//5% bigger INtervall...
	if(posX.size()!=posY.size()||posX.size()==0) {
		cerr<<"ERROR HistogrammSaver::CreateScatterHisto vectors have different size "<<posX.size()<<" "<<posY.size()<<" "<<name<<endl;
		return new TH2F();
	}
	Float_t maxX = posX.at(0);
	Float_t maxY = posY.at(0);
	Float_t minX = posY.at(0);
	Float_t minY = posY.at(0);
	for(UInt_t i=0;i<posX.size();i++){
		if (posX.at(i)<minRangeX||posX.at(i)>maxRangeX)
			continue;
		if (posY.at(i)<minRangeY||posY.at(i)>maxRangeY)
					continue;
		if(posX.at(i)>maxX)maxX=posX.at(i);
		else if(posX.at(i)<minX)minX=posX.at(i);
		if(posY.at(i)>maxY)maxY=posY.at(i);
		else if(posY.at(i)<minY)minY=posY.at(i);
	}
	//cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
	Float_t deltaX=maxX-minX;
	Float_t deltaY=maxY-minY;
	TH2F* histo = new TH2F(name.c_str(),name.c_str(),nBinsX,minX-factor*deltaX,maxX+factor*deltaX,nBinsY,minY-factor*deltaY,maxY+factor*deltaY);
	for(UInt_t i=0;i<posX.size();i++){
		if (posX.at(i) < minRangeX || posX.at(i) > maxRangeX)
			continue;
		if (posY.at(i) < minRangeY || posY.at(i) > maxRangeY)
					continue;
		histo->Fill(posX.at(i),posY.at(i));
	}
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

TGraph HistogrammSaver::CreateDipendencyGraph(std::string name, std::vector<Float_t> vecY, std::vector<Float_t> vecX,ULong_t maxSize)
{
	if(vecY.size()!=vecX.size()||vecX.size()==0) {
		cerr<<"ERROR HistogrammSaver::CreateDipendencyGraph vectors have different size "<<vecY.size()<<" "<<vecX.size()<<": "<<name<<endl;
		return TGraph();
	}
	//cout<<"HistogrammSaver::CREATE Scatterplot:\""<<name<<"\" with "<<posX.size()<<" Entries"<<endl;
	ULong_t size = TMath::Min(maxSize,(ULong_t)vecY.size());
	TGraph hGraph = TGraph(size,&vecX.at(0),&vecY.at(0));
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
TH1F* HistogrammSaver::CreateDistributionHisto(std::string name, std::vector<Float_t> vec, UInt_t nBins,EnumAxisRange range,Float_t xmin,Float_t xmax, Float_t factor)
{
	int verbosity = 0;
//	Float_t factor = 0.05;//5% bigger INtervall...
	if(vec.size()==0)
		return new TH1F(name.c_str(),name.c_str(),nBins,0.,1.);
	Float_t max = vec.at(0);
	Float_t min = vec.at(0);
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
		for(UInt_t i=0;i<vec.size();i++){
			Float_t x = vec.at(i);
			if(x<xmin||x>xmax)
				continue;
			mean2+=x;
			sigma2+=x*x;
			n++;
		}
		mean2/=(Float_t)n;
		sigma2/=(Float_t)n;

		Float_t mean=0;
		Float_t sigma=0;
		UInt_t nEvents=0;
		for(UInt_t i=0;i<vec.size();i++){
			Float_t x = vec.at(i);
			if(x<xmin||x>xmax)
				continue;
			if( (x-mean2)/sigma2>3.)
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
		max =xmax;
		min=xmin;
		if(verbosity>3)cout<<" manual: "<<min<<" - " << max <<endl;
	}

	TH1F* histo = new TH1F(name.c_str(),name.c_str(),nBins,min,max);
	for(UInt_t i=0;i<vec.size();i++){
		Float_t x = vec.at(i);
		if(x<xmin||x>xmax)
			continue;
		histo->Fill(vec.at(i));
	}
	int ntries=0;
	while ((histo->GetBinContent(histo->GetMaximumBin())/histo->GetEntries())<0.05&&ntries<3){//todo change hardcoding
		histo->Rebin();ntries++;
	}
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
