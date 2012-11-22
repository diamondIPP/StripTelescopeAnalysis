/*
 * TAnalysisOf3dDiamonds.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair,iain
 */

#include "../include/TAnalysisOf3dDiamonds.hh"

TAnalysisOf3dDiamonds::TAnalysisOf3dDiamonds(TSettings *newSettings) {
	if(newSettings!=0)
		this->settings=newSettings;
	else exit(0);

	UInt_t runNumber=settings->getRunNumber();

	//htmlLandau=new THTMLLandaus(settings);

	settings->goTo3dDiamondTreeDir();
	eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),runNumber);
	histSaver=new HistogrammSaver();
	settings->goTo3dDiamondAnalysisDir();

	histSaver->SetPlotsPath(settings->get3dDiamondAnalysisPath());
	histSaver->SetRunNumber(runNumber);
	//htmlLandau->setFileGeneratingPath(settings->getSelectionAnalysisPath());//todo Write html3dDiamond
	settings->goTo3dDiamondTreeDir();
	initialiseHistos();

	cout<<"end initialise"<<endl;
}

TAnalysisOf3dDiamonds::~TAnalysisOf3dDiamonds() {
	//htmlLandau->generateHTMLFile();
	if(eventReader!=0) delete eventReader;
	if(histSaver!=0)   delete histSaver;
	//if(htmlLandau!=0)  delete htmlLandau;
	settings->goToOutputDir();
}

void TAnalysisOf3dDiamonds::doAnalysis(UInt_t nEvents) {
	cout<<"analyze selection data..."<<endl;
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		eventReader->LoadEvent(nEvent);
		analyseEvent();
	}
	saveHistos();
}

void TAnalysisOf3dDiamonds::initialiseHistos() {
}

void TAnalysisOf3dDiamonds::saveHistos() {
	TH2F *hPredictedPositionDiamond = HistogrammSaver::CreateScatterHisto("hPredictedPositionDiamond",vecYPredicted,vecXPredicted,512);
	hPredictedPositionDiamond->SetTitle("Predicted Position of Silicon Track into the diamond Plane");
	hPredictedPositionDiamond->GetXaxis()->SetTitle("predicted x Position");
	hPredictedPositionDiamond->GetYaxis()->SetTitle("predicted y Position");
	histSaver->SaveHistogram(hPredictedPositionDiamond);

	TH2F *hPredictedPositionDiamondHit = HistogrammSaver::CreateScatterHisto("hPredictedPositionDiamondHit",vecYPredictedDiamondHit,vecXPredictedDiamondHit,512);
	hPredictedPositionDiamondHit->SetTitle("Predicted Position of Silicon Track into the diamond Plane with one diamond Cluster");
	hPredictedPositionDiamondHit->GetXaxis()->SetTitle("predicted x Position");
	hPredictedPositionDiamondHit->GetYaxis()->SetTitle("predicted y Position");
	histSaver->SaveHistogram(hPredictedPositionDiamondHit);

	TH1F* hLandauBeforeCut = HistogrammSaver::CreateDistributionHisto("hLandauBeforeCut",vecPHDiamondHit,256);
	hLandauBeforeCut->SetTitle("Landau before a Lab Fiducial Cut");
	hLandauBeforeCut->GetXaxis()->SetTitle("PH of diamond cluster");
	hLandauBeforeCut->GetYaxis()->SetTitle("number of entries #");
	histSaver->SaveHistogram(hLandauBeforeCut);

	vector<Float_t> vecPHdiamond_LAB_fiducut;
	Float_t yLabFidCutMin = 75;
	Float_t yLabFidCutMax = 130;
	Float_t xLabFidCutMin = 125;
	Float_t xLabFidCutMax = 160;

	for(UInt_t i=0;i<vecPHDiamondHit.size();i++){
		Float_t xPos = vecXPredictedDiamondHit.at(i);
		Float_t yPos = vecYPredictedDiamondHit.at(i);
		Float_t PH = vecPHDiamondHit.at(i);
		if ( ( xLabFidCutMin <xPos && xPos< xLabFidCutMax) && (yLabFidCutMin <yPos && yPos<yLabFidCutMax))
			vecPHdiamond_LAB_fiducut.push_back(PH);
	}
	TH1F* hLandauAfterCut = HistogrammSaver::CreateDistributionHisto("hLandauAfterCut",vecPHdiamond_LAB_fiducut,256);
	hLandauAfterCut->SetTitle("Landau before a Lab Fiducial Cut");
	hLandauAfterCut->GetXaxis()->SetTitle("PH of diamond cluster");
	hLandauAfterCut->GetYaxis()->SetTitle("number of entries #");
	histSaver->SaveHistogram(hLandauAfterCut);
	TCanvas* c1 = new TCanvas("cPulseHeights");
	c1->cd();
	hLandauBeforeCut->Draw();
	hLandauAfterCut->SetLineColor(kBlue);
	hLandauAfterCut->Draw("same");
	histSaver->SaveCanvas(c1);
	TH2F*hPHvsPredictedXPos = HistogrammSaver::CreateScatterHisto("hPHvsPredictedXPos",vecXPredictedDiamondHit,vecPHDiamondHit,512);
	hPHvsPredictedXPos->SetTitle("predicted pos of hit in x vs Pulse Height of Diamond Cluster");
	hPHvsPredictedXPos->GetYaxis()->SetTitle("predicted x Position");
	hPHvsPredictedXPos->GetXaxis()->SetTitle("PH of Diamond Cluster in ADC counts");
	histSaver->SaveHistogram(hPHvsPredictedXPos);

	TH2F*hPHvsPredictedYPos = HistogrammSaver::CreateScatterHisto("hPHvsPredictedYPos",vecYPredictedDiamondHit,vecPHDiamondHit,512);
	hPHvsPredictedYPos->SetTitle("predicted pos of hit in y vs Pulse Height of Diamond Cluster");
	hPHvsPredictedYPos->GetYaxis()->SetTitle("predicted y Position");
	hPHvsPredictedYPos->GetXaxis()->SetTitle("PH of Diamond Cluster in ADC counts");
	histSaver->SaveHistogram(hPHvsPredictedYPos);
}

void TAnalysisOf3dDiamonds::analyseEvent() {
	if(!eventReader->isValidTrack())
		return;
	vector<UInt_t> vecSilPlanes;
	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++)
		vecSilPlanes.push_back(pl);
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	TPositionPrediction *predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	Float_t xPos = predictedPosition->getPositionX();
	Float_t yPos = predictedPosition->getPositionY();
	if(!(predictedPosition->getChi2X()<2 && predictedPosition->getChi2Y()<2))
		return;
	if(xPos<-500||yPos<-500)
		return;
	vecXPredicted.push_back(xPos);
	vecYPredicted.push_back(yPos);
	if(!eventReader->getNDiamondClusters()==1)
		return;
	vecXPredictedDiamondHit.push_back(xPos);
	vecYPredictedDiamondHit.push_back(yPos);
	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	vecPHDiamondHit.push_back(diamondCluster.getCharge());
}
