/*
 * TAnalysisOfSelection.cpp
 *
 *  Created on: May 18, 2012
 *      Author: bachmair
 */

#include "../include/TAnalysisOfSelection.hh"

TAnalysisOfSelection::TAnalysisOfSelection(TSettings *settings) {
	if(settings!=0)
		this->settings=settings;
	else exit(0);

	sys = gSystem;
	UInt_t runNumber=settings->getRunNumber();

  sys->MakeDirectory(settings->getAbsoluteOuputPath().c_str());
	htmlLandau=new THTMLLandaus(settings);
  sys->cd(settings->getAbsoluteOuputPath().c_str());

	stringstream  filepath;
	filepath.str("");
	filepath<<"selectionData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
//	htmlPedestal->setMainPath((string)(sys->pwd()));
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str(),settings->getRunNumber());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("selectionAnalysis");
	sys->cd("selectionAnalysis");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
//	htmlPedestal->setSubdirPath("selectionAnalysis");
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	htmlLandau->setFileGeneratingPath(plotsPath.str());
	sys->cd("..");
	initialiseHistos();

	cout<<"end initialise"<<endl;
}

TAnalysisOfSelection::~TAnalysisOfSelection() {
	htmlLandau->generateHTMLFile();
	delete eventReader;
	delete histSaver;
	delete htmlLandau;
	sys->cd("..");
}

void TAnalysisOfSelection::doAnalysis(UInt_t nEvents)
{
	cout<<"analyze selection data..."<<endl;
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
		analyseEvent();
	}
	saveHistos();

}

void TAnalysisOfSelection::initialiseHistos()
{
	histoLandauDistribution= new TH2F("hLandauDiamond_OneCluster","hLandauDiamond_OneCluster",512,0,4096,8,0.5,8.5);
	histoLandauDistribution->GetXaxis()->SetTitle("Charge in ADC counts");
	histoLandauDistribution->GetYaxis()->SetTitle("ClusterSize");
	hFidCut= new TH2F("hFidCut","hFidCut",256,0,255,256,0,255);
	hFidCut->GetXaxis()->SetTitle("FidCutValue in X");
	hFidCut->GetYaxis()->SetTitle("FidCutValue in Y");
	hClusterPosition=new TH1F("hClusterPositionDia","Events which have a valid Silicon Track",128,0,127);
	hClusterPosition->GetXaxis()->SetTitle("highes Cluster Channel Position");
	hClusterPosition->GetYaxis()->SetTitle("number of Events #");
	h3dDiamond = new TH1F("h3dDiamond","Sum of Charge for all 18 3d-channels",4096,0,4095);
	h3dDiamond_hit = new TH1F("h3dDiamond_hit","Sum of Charge for all 18 3d-channels with a Hit",4096,0,4095);
	hNoDiamond = new TH1F("hNoDiamond","Sum of Charge for all 18 no-channels",4096,0,4095);
	hNoDiamond_hit = new TH1F("hNoDiamond_hit","Sum of Charge for all 18 no-channels with a Hit",4096,0,4095);
}

void TAnalysisOfSelection::saveHistos()
{
	LandauGaussFit landauGauss;
	histSaver->SaveHistogram(histoLandauDistribution);
	vector <Float_t> vecMP;
	vector <Float_t> vecClusSize;
	vector <Float_t> vecWidth;
	vector <Float_t> vecXError;
	vector <Float_t> vecHistoMax;
	vector <Float_t> vecHistoMean;
	vector <Float_t> vecHistoMeanGaus;
	vector <Float_t> vecHistoMeanLandau;
	TH1F* histoClusSize = (TH1F*)histoLandauDistribution->ProjectionY("ClusterSizeDiamond",0,4096);
	TH1F *histo = (TH1F*)histoLandauDistribution->ProjectionX("hPulseHeightDiamondAll",0,8);
	histo->GetYaxis()->SetTitle("number of Entries #");
	Float_t histoMean,histoMax,histoRMS,histoMeanGausFit;
	Double_t xmin,xmax;
	TF1* fit=0;
	TF1* gausFit=0;
	histoMean = histo->GetMean();
	histoMax = histo->GetBinCenter(histo->GetMaximumBin());
	histoRMS = histo->GetRMS();
	xmin=histoMax-histoRMS;
	xmax=histoMax+histoRMS;
	gausFit = new TF1("gausFit","gaus",xmin,xmax);
//	cout<<"gausFit: "<<gausFit<<endl;
	histo->Fit(gausFit,"","same+",xmin,xmax);
	fit = landauGauss.doLandauGaussFit(histo);
//	cout <<"gausFit:"<<gausFit->GetTitle()<<" is a:"<< gausFit->ClassName()<<" "<<gausFit->GetNpar()<<endl;
	histoMeanGausFit = gausFit->GetParameter(1);
	vecWidth.push_back(fit->GetParameter(0));
	vecHistoMax.push_back(histoMax);
	vecHistoMean.push_back(histoMean);
	vecHistoMeanGaus.push_back(histoMeanGausFit);
	vecHistoMeanLandau.push_back(fit->GetParameter(1));

	histSaver->SaveHistogram(histo);
	Float_t width=fit->GetParameter(0);
	Float_t MP = fit->GetParameter(1);
	Float_t area = fit->GetParameter(2);
	Float_t gWidth = fit->GetParameter(3);
	//	vecMP.push_back(MP);
	//	vecWidth.push_back(width);
	//vecClusSize.push_back(0);
	//vecXError.push_back(0);
	for(UInt_t clusSize=1;clusSize<8;clusSize++){
		stringstream name;
		name<< "hPulseHeigthDiamond_"<<clusSize<<"_ClusterSize";
		TH1F* histo = (TH1F*)histoLandauDistribution->ProjectionX(name.str().c_str(),clusSize,clusSize);

		histo->SetTitle(name.str().c_str());
		histo->GetYaxis()->SetTitle("number of Entries #");
		TF1* fitCS=0;
		if(clusSize<5){
		  int nTries=0;
		  while(histo->GetMaximum()<histo->GetEntries()*0.1&&nTries<5)
		    histo->Rebin(),nTries++;
		  histoMean = histo->GetMean();
		  histoMax = histo->GetBinCenter(histo->GetMaximumBin());
		  histoRMS = histo->GetRMS();
		  xmin=histoMax-histoRMS, xmax=histoMax+histoRMS;
		  if(gausFit!=0)delete gausFit;
		  gausFit = new TF1("gausFit","gaus",xmin,xmax);
		  histo->Fit(gausFit,"","sames+",xmin,xmax);
		  histoMeanGausFit = gausFit->GetParameter(1);
		  if(fitCS!=0)delete fitCS;
			fitCS = landauGauss.doLandauGaussFit(histo);
			vecMP.push_back(fitCS->GetParameter(1));
			vecClusSize.push_back(clusSize);
			vecXError.push_back(.5);
			vecWidth.push_back(fitCS->GetParameter(0));
			vecHistoMax.push_back(histoMax);
			vecHistoMean.push_back(histoMean);
			vecHistoMeanGaus.push_back(histoMeanGausFit);
			vecHistoMeanLandau.push_back(fitCS->GetParameter(1));
		}
		histSaver->SaveHistogram(histo);
		delete histo;
	}
	cout<<"Create ErrorGraph"<<endl;
	TGraphErrors* graph = new TGraphErrors(vecMP.size(),&vecClusSize.at(0),&vecMP.at(0),&vecXError.at(0),&vecWidth.at(0));
	stringstream name;
	name<<"MPV of Landau for one ClusterSizes";
	graph->SetTitle(name.str().c_str());
	name.clear();name.str("");name.clear();
	name<<"hMPV_Landau_diff_ClusterSizes";
	graph->SetName(name.str().c_str());
	graph->Draw("APLE1 goff");
	graph->GetXaxis()->SetTitle("Cluster Size");
	graph->GetYaxis()->SetTitle("Most Probable Value of Landau");
	graph->SetMarkerColor(kGreen);
	graph->SetMarkerStyle(22);
	graph->SetFillColor(kWhite);
	graph->SetLineWidth(2);
	cout<<"Create Canvas"<<endl;
	TCanvas *c1= new TCanvas("cMVP_Landau_vs_ClusterSize","cMVP_Landau_vs_ClusterSize",800,600);
	c1->cd();
	Float_t xVal[] = {0,5};
	Float_t exVal[] = {0.5,0.5};
	Float_t yVal[] = {MP,MP};
	Float_t eyVal[]= {width,width};
	cout<<"Create ErrorGraph MEAN"<<endl;
	TGraphErrors *gMVP = new TGraphErrors(2,xVal,yVal,exVal,eyVal);
	gMVP->SetName("gMPV_ALL");
	gMVP->SetTitle("MVP of all Clusters");
	gMVP->SetFillColor(kRed);
	gMVP->SetFillStyle(3002);
	gMVP->SetLineColor(kBlue);
	cout<<"Create MultiGraph"<<endl;
	TMultiGraph *mg = new TMultiGraph("mgMVP_ClusterSize","MVP of Landau vs. ClusterSize");
	mg->Add(gMVP,"3L");
	mg->Add(graph,"PLE1");
	cout<<"Draw Canvas"<<endl;
	mg->Draw("a");
	mg->GetXaxis()->SetTitle("Cluster Size of Diamond");
	mg->GetYaxis()->SetTitle("MPV of Landau  ");
	mg->GetXaxis()->SetRangeUser(0.5,4.5);
	TLegend *leg = c1->BuildLegend(0.15,0.55,0.6,0.8);
	leg->SetFillColor(kWhite);
	cout<<"Save Canvas"<<endl;
	histSaver->SaveCanvas(c1);

//	TLine *lMVP = new TLine(graph->GetXaxis()->GetXmin(),MP,graph->GetXaxis()->GetXmax(),MP);
//	TLine *lMVPplus = new TLine(graph->GetXaxis()->GetXmin(),MP+width,graph->GetXaxis()->GetXmax(),MP+width);
//	TLine *lMVPminus = new TLine(graph->GetXaxis()->GetXmin(),MP-width,graph->GetXaxis()->GetXmax(),MP-width);
	histSaver->SaveGraph(graph,name.str(),"APLE1");
	htmlLandau->addLandauDiamond(width,MP,area,gWidth);
	htmlLandau->addLandauDiamondTable(vecHistoMean,vecHistoMax,vecHistoMeanGaus,vecHistoMeanLandau);


	histoClusSize->SetTitle("ClusterSize Diamond");
	histoClusSize->GetXaxis()->SetTitle("ClusterSize");
	histoClusSize->GetYaxis()->SetTitle("Number of Entries #");
	histSaver->SaveHistogram(histoClusSize);

	htmlLandau->addSection("ClusterSize Diamond",htmlLandau->putImageOfPath("ClusterSizeDiamond","png",50));
	if(fit!=0)delete fit;
	delete histo;
	delete histoClusSize;
	delete histoLandauDistribution;
	delete mg;
	delete c1;

	histSaver->SaveHistogram(hFidCut);
	delete hFidCut;
	histSaver->SaveHistogram(hClusterPosition,0,1);
	delete hClusterPosition;
	cout<<"Save histos "<<endl;
	cout<<h3dDiamond->GetEntries()<<endl;
	histSaver->SaveHistogram(h3dDiamond,0,1);
	cout<<hNoDiamond->GetEntries()<<endl;
	histSaver->SaveHistogram(hNoDiamond,0,1);
	cout<<h3dDiamond_hit->GetEntries()<<endl;
  histSaver->SaveHistogram(h3dDiamond_hit,0,1);
  cout<<hNoDiamond_hit->GetEntries()<<endl;
  histSaver->SaveHistogram(hNoDiamond_hit,0,1);
  cout<<"Save canvas"<<endl;
  TCanvas *c2= new TCanvas("c2","c2",1024,800);
  c2->cd();
  h3dDiamond_hit->Draw();
  hNoDiamond_hit->SetLineColor(kBlue);
  hNoDiamond_hit->Draw("same");
  histSaver->SaveCanvas(c2);
  delete c2;
	delete h3dDiamond;
	delete hNoDiamond;
}

void TAnalysisOfSelection::analyseEvent()
{
  Float_t fiducialValueX=0;
  Float_t fiducialValueY=0;

	if(eventReader->isValidTrack()){//
	//if(eventReader->useForAnalysis()||eventReader->useForAlignment()){

    for(UInt_t plane=0;plane<4;plane++){
      fiducialValueX+=eventReader->getCluster(plane,TPlaneProperties::X_COR,0).getPosition();
      fiducialValueY+=eventReader->getCluster(plane,TPlaneProperties::Y_COR,0).getPosition();
    }
    fiducialValueX/=4.;
    fiducialValueY/=4.;
    Float_t charge3d=0;
    for(UInt_t ch=0;ch<18;ch++){
      Float_t rawSignal=eventReader->getRawSignal(TPlaneProperties::getDetDiamond(),ch);
      if(charge3d<rawSignal)charge3d=rawSignal;
    }
    h3dDiamond->Fill(charge3d);
    Float_t chargeNo=0;
    for(UInt_t ch=36;ch<54;ch++){
      Float_t rawSignal=eventReader->getRawSignal(TPlaneProperties::getDetDiamond(),ch);
      if(chargeNo<rawSignal)charge3d=rawSignal;
    }
    hNoDiamond->Fill(chargeNo);
    if(eventReader->getNClusters(TPlaneProperties::getDetDiamond())<=0)
        return;
    hNoDiamond_hit->Fill(chargeNo);
    h3dDiamond_hit->Fill(charge3d);

    hFidCut->Fill(fiducialValueX,fiducialValueY);
		TCluster cluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
		Float_t charge = cluster.getCharge(false);
		UInt_t clustSize = cluster.size();
		if(clustSize>8) clustSize=8;
//		cout<<nEvent<<":\t"<<charge<<endl;
		histoLandauDistribution->Fill(charge,clustSize);
		Float_t pos = cluster.getPosition(TCluster::maxValue,0);
		hClusterPosition->Fill(pos);
	}
}









