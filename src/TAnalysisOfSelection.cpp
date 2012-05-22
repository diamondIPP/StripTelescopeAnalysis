/*
 * TAnalysisOfSelection.cpp
 *
 *  Created on: May 18, 2012
 *      Author: bachmair
 */

#include "../include/TAnalysisOfSelection.hh"

TAnalysisOfSelection::TAnalysisOfSelection(TSettings *settings) {
	// TODO Auto-generated constructor stub
	if(settings!=0)
		this->settings=settings;
	else exit(0);

	sys = gSystem;
	stringstream  runString;
	UInt_t runNumber=settings->getRunNumber();
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());
	htmlLandau=new THTMLLandaus(settings);
	sys->cd(runString.str().c_str());
	htmlLandau->setMainPath(sys->pwd());
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
	sys->cd("..");
	initialiseHistos();

	cout<<"end initialise"<<endl;
}

TAnalysisOfSelection::~TAnalysisOfSelection() {
	// TODO Auto-generated destructor stub
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
}

void TAnalysisOfSelection::saveHistos()
{
	LandauGaussFit landauGauss;
	histSaver->SaveHistogram(histoLandauDistribution);
	vector <Float_t> vecMP;
	vector <Float_t> vecClusSize;
	vector <Float_t> vecWidth;
	vector <Float_t> vecXError;
	TH1F* histoClusSize = (TH1F*)histoLandauDistribution->ProjectionY("ClusterSizeDiamond",0,4096);
	TH1F *histo = (TH1F*)histoLandauDistribution->ProjectionX("hPulseHeightDiamondAll",0,8);
	histo->GetYaxis()->SetTitle("number of Entries #");
	TF1* fit = landauGauss.doLandauGaussFit(histo);
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
		TF1* fit;
		if(clusSize<5){
			fit = landauGauss.doLandauGaussFit(histo);
			vecMP.push_back(fit->GetParameter(1));
			vecClusSize.push_back(clusSize);
			vecXError.push_back(.5);
			vecWidth.push_back(fit->GetParameter(0));
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
	histoClusSize->SetTitle("ClusterSize Diamond");
	histoClusSize->GetXaxis()->SetTitle("ClusterSize");
	histoClusSize->GetYaxis()->SetTitle("Number of Entries #");
	histSaver->SaveHistogram(histoClusSize);
	htmlLandau->addSection("ClusterSize Diamond",htmlLandau->putImageOfPath("ClusterSizeDiamond","png",50));
	delete histo;
	delete histoClusSize;
	delete histoLandauDistribution;
	delete mg;
	delete c1;
}

void TAnalysisOfSelection::analyseEvent()
{
	if(eventReader->useForAnalysis()||eventReader->useForAlignment()){
		TCluster cluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
		Float_t charge = cluster.getCharge(false);
		UInt_t clustSize = cluster.size();
		if(clustSize>8)clustSize=8;
//		cout<<nEvent<<":\t"<<charge<<endl;
		histoLandauDistribution->Fill(charge,clustSize);
	}
}









