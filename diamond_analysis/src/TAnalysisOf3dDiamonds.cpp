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
	verbosity = settings->getVerbosity();
	settings->goTo3dDiamondTreeDir();
	eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
	histSaver=new HistogrammSaver();
	settings->goTo3dDiamondAnalysisDir();

	histSaver->SetPlotsPath(settings->get3dDiamondAnalysisPath());
	histSaver->SetRunNumber(runNumber);
	//htmlLandau->setFileGeneratingPath(settings->getSelectionAnalysisPath());//todo Write html3dDiamond
	settings->goTo3dDiamondTreeDir();
	initialiseHistos();
	vecSilPlanes.clear();
	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++)
		vecSilPlanes.push_back(pl);
	subjectPlane = TPlaneProperties::getDiamondPlane();
	predictedPosition=0;
	events=0;
	noValidTrack=0;
	tooHighChi2=0;
	invalidPositionPrediction=0;
	notOneDiaCluster=0;
	notInFidCut=0;
	saturatedCluster=0;
	areaAndFidCutDoNotAgree=0;
	hitNotInOneArea=0;
	validEvents=0;

	vecChannelMetricPos.clear();
	vecChannel.clear();
	cout<<"end initialise"<<endl;
}

TAnalysisOf3dDiamonds::~TAnalysisOf3dDiamonds() {
	Int_t n = TMath::Log10(events)+2;
	cout<<"\n\nCounting: "<<endl;
	cout<<"\tnoValidTrack:             \t"<<setw(n)<<right<<noValidTrack<<"\t"<<right<<setw(6)<<(Float_t)noValidTrack/events*100.<<" %\n";;
	cout<<"\tnoDiaCluster:             \t"<<setw(n)<<right<<noDiaCluster<<"\t"<<right<<setw(6)<<(Float_t)noDiaCluster/events*100.<<" %\n";;
	cout<<"\tnotOneDiaCluster:         \t"<<setw(n)<<right<<notOneDiaCluster<<"\t"<<right<<setw(6)<<(Float_t)notOneDiaCluster/events*100.<<" %\n";;
	cout<<"\tsaturatedCluster:         \t"<<setw(n)<<right<<saturatedCluster<<"\t"<<right<<setw(6)<<(Float_t)saturatedCluster/events*100.<<" %\n";;
	cout<<"\tmaskedClusters:		   \t"<<setw(n)<<right<<maskedClusters<<"\t"<<right<<setw(6)<<(Float_t)maskedClusters/events*100.<<" %\n";
	cout<<"\thitNotInOneArea:          \t"<<setw(n)<<right<<hitNotInOneArea<<"\t"<<right<<setw(6)<<(Float_t)hitNotInOneArea/events*100.<<" %\n";;
	cout<<"\tnotInFidCut:              \t"<<setw(n)<<right<<notInFidCut<<"\t"<<right<<setw(6)<<(Float_t)notInFidCut/events*100.<<" %\n";;
	cout<<"\ttooHighChi2:              \t"<<setw(n)<<right<<tooHighChi2<<"\t"<<right<<setw(6)<<(Float_t)tooHighChi2/events*100.<<" %\n";;
	cout<<"\tareaAndFidCutDoNotAgree:  \t"<<setw(n)<<right<<areaAndFidCutDoNotAgree<<"\t"<<right<<setw(6)<<(Float_t)areaAndFidCutDoNotAgree/events*100.<<" %\n";;
	cout<<"\tinvalidPositionPrediction:\t"<<setw(n)<<right<<invalidPositionPrediction<<"\t"<<right<<setw(6)<<(Float_t)invalidPositionPrediction/events*100.<<" %\n";;
	cout<<"\t                          \t"<<setw(n)<<setfill('-')<<"-"<<setfill(' ')<<"\n";
	cout<<"\tValid Events:             \t"<<setw(n)<<right<<validEvents<<"\t"<<right<<setw(6)<<(Float_t)validEvents/events*100.<<" %\n\n"<<endl;;
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
	hChi2X = new TH1F("hChi2X","hChi2X",512,0,20);
	hChi2X->GetXaxis()->SetTitle("#chi^{2}_{X}");
	hChi2X->GetYaxis()->SetTitle("number of entries");
	hChi2Y = new TH1F("hChi2Y","hChi2Y",512,0,20);
	hChi2Y->GetXaxis()->SetTitle("#chi^{2}_{Y}");
	hChi2Y->GetYaxis()->SetTitle("number of entries");
	hChi2XY = new TH2F("hChi2XY","hChi2XY",512,0,20,512,0,20);
	hChi2XY->GetXaxis()->SetTitle("#chi^{2}_{X}");
	hChi2XY->GetYaxis()->SetTitle("#chi^{2}_{Y}");
	hChi2XY->GetZaxis()->SetTitle("number of entries");

	hChi2 = new TH1F("hChi2","Max of #chi^{2} in x and y",512,0,20);
	hChi2->GetXaxis()->SetTitle("Max(#chi^{2}_{X},#chi^{2}_{Y})");
	hChi2->GetYaxis()->SetTitle("number of entries");
	hValidSiliconAndOneDiamondHit = new TH2F("hValidSiliconAndOneDiamondHit","hValidSiliconAndOneDiamondHit",512,5500,10000,512,3800,5000);
	hValidSiliconAndOneDiamondHit->GetXaxis()->SetTitle("Predicted X Position / #mum");
	hValidSiliconAndOneDiamondHit->GetYaxis()->SetTitle("Predicted Y Position / #mum");
	hValidSiliconAndOneDiamondHitNotMasked = new TH2F("hValidSiliconAndOneDiamondHitNotMasked","hValidSiliconAndOneDiamondHitNotMasked",512,5500,10000,512,3800,5000);
	hValidSiliconAndOneDiamondHitNotMasked->GetXaxis()->SetTitle("Predicted X Position / #mum");
	hValidSiliconAndOneDiamondHitNotMasked->GetYaxis()->SetTitle("Predicted Y Position / #mum");
	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels = new TH2F("hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels","hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels",512,5500,10000,512,3800,5000);
	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels->GetXaxis()->SetTitle("Predicted X Position / #mum");
	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels->GetYaxis()->SetTitle("Predicted Y Position / #mum");
	hValidSiliconAndOneDiamondHitInOneArea = new TH2F("hValidSiliconAndOneDiamondHitInOneArea","hValidSiliconAndOneDiamondHitInOneArea",512,5500,10000,512,3800,5000);
	hValidSiliconAndOneDiamondHitInOneArea->GetXaxis()->SetTitle("Predicted X Position / #mum");
	hValidSiliconAndOneDiamondHitInOneArea->GetYaxis()->SetTitle("Predicted Y Position / #mum");
	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut = new TH2F("hValidSiliconAndOneDiamondHitInSameAreaAndFidCut","hValidSiliconAndOneDiamondHitInSameAreaAndFidCut",512,5500,10000,512,3800,5000);
	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->GetXaxis()->SetTitle("Predicted X Position / #mum");
	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->GetYaxis()->SetTitle("Predicted Y Position / #mum");
	hAreaVsFidCut = new TH2F("hAreaVsFidCut","hAreaVsFidCut",6,-1.5,4.5,6,-1.5,4.5);
	hAreaVsFidCut->GetYaxis()->SetTitle("Diamond Area");
	hAreaVsFidCut->GetXaxis()->SetTitle("Fiducial Cut");
	for (UInt_t i=1;i<=settings->get3dFidCuts()->getNFidCuts();i++){
		Float_t xMin = settings->get3dFidCuts()->getMinFiducialX(i);
		Float_t xMax = settings->get3dFidCuts()->getMaxFiducialX(i);
		cout<<"3d - fidCut no "<<i<<" of " <<settings->get3dFidCuts()->getNFidCuts();
		cout<< "\t"<<xMin<<"-"<<xMax<<endl;
		UInt_t xBins = (Int_t)(2*(xMax-xMin));
		TString name = TString::Format("hChargeVsFidX_HitInFidCutNo%d",i);
	}

}

void TAnalysisOf3dDiamonds::saveHistos() {
	histSaver->SaveHistogram(hChi2X);
	histSaver->SaveHistogram(hChi2Y);
	histSaver->SaveHistogramWithCutLine(hChi2,settings->getChi2Cut3D());
	histSaver->SaveHistogram(hChi2XY);
	delete hChi2X;
	delete hChi2Y;
	delete hChi2;
	delete hChi2XY;
	histSaver->SaveHistogram(hValidSiliconAndOneDiamondHit);
	histSaver->SaveHistogram(hValidSiliconAndOneDiamondHitNotMasked);
	histSaver->SaveHistogram(hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels);
	histSaver->SaveHistogram(hValidSiliconAndOneDiamondHitInOneArea);
	histSaver->SaveHistogram(hValidSiliconAndOneDiamondHitInSameAreaAndFidCut);
	delete hValidSiliconAndOneDiamondHit;
	delete hValidSiliconAndOneDiamondHitNotMasked;
	delete hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels;
	delete hValidSiliconAndOneDiamondHitInOneArea;
	delete hValidSiliconAndOneDiamondHitInSameAreaAndFidCut;

	TH2F* hChannelVsPosition = HistogrammSaver::CreateScatterHisto("hChannelVsPosition",vecChannelMetricPos,vecChannel,512);
	hChannelVsPosition->SetTitle("highest Channel vs predicted Postion");
	hChannelVsPosition->GetXaxis()->SetTitle("highest Channel /ch");
	hChannelVsPosition->GetYaxis()->SetTitle("predicted Position /#mum");
	histSaver->SaveHistogram(hChannelVsPosition);
	delete hChannelVsPosition;

	TH1F *hSpecialChannelMetricPos = HistogrammSaver::CreateDistributionHisto("hSpecialChannelMetricPos",vecSpecialChannelMetricPos);
	hSpecialChannelMetricPos->SetTitle("metric xPos(lab) of channels 29, 88 or 59");
	hSpecialChannelMetricPos->GetXaxis()->SetTitle("metric Pos / #mum");
	hSpecialChannelMetricPos->GetYaxis()->SetTitle("number of entries");
	histSaver->SaveHistogram(hSpecialChannelMetricPos);
	delete hSpecialChannelMetricPos;
	cout<<vecXFidCut.size()<<"-"<<vecXPredicted.size()<<endl;

	TH2F *hPredictionVsFidCutX = HistogrammSaver::CreateScatterHisto("hPredictionVsFidCutX",vecXFidCut,vecXPredicted,512);
	hPredictionVsFidCutX->SetTitle("Predicted Position of Silicon Track into the diamond Plane");
	hPredictionVsFidCutX->GetXaxis()->SetTitle("predicted x Position /#mum");
	hPredictionVsFidCutX->GetYaxis()->SetTitle("FidCut X /ch");
	histSaver->SaveHistogram(hPredictionVsFidCutX);
	delete hPredictionVsFidCutX;

	TH2F *hPredictionVsFidCutY = HistogrammSaver::CreateScatterHisto("hPredictionVsFidCutY",vecYFidCut,vecYPredicted,512);
	hPredictionVsFidCutY->SetTitle("Predicted Position of Silicon Track into the diamond Plane");
	hPredictionVsFidCutY->GetXaxis()->SetTitle("predicted y Position /#mum");
	hPredictionVsFidCutY->GetYaxis()->SetTitle("FidCut Y /ch");
	histSaver->SaveHistogram(hPredictionVsFidCutY);

	TH2F *hPredictedPositionDiamond = HistogrammSaver::CreateScatterHisto("hPredictedPositionDiamond",vecYPredicted,vecXPredicted,512);
	hPredictedPositionDiamond->SetTitle("Predicted Position of Silicon Track into the diamond Plane");
	hPredictedPositionDiamond->GetXaxis()->SetTitle("predicted x Position");
	hPredictedPositionDiamond->GetYaxis()->SetTitle("predicted y Position");
	histSaver->SaveHistogram(hPredictedPositionDiamond);
	delete hPredictedPositionDiamond;

	TH2F *hPredictedPositionDiamondHit = HistogrammSaver::CreateScatterHisto("hPredictedPositionDiamondHit",vecYPredictedDiamondHit,vecXPredictedDiamondHit,512);
	hPredictedPositionDiamondHit->SetTitle("Predicted Position of Silicon Track into the diamond Plane with one diamond Cluster");
	hPredictedPositionDiamondHit->GetXaxis()->SetTitle("predicted x Position");
	hPredictedPositionDiamondHit->GetYaxis()->SetTitle("predicted y Position");
	histSaver->SaveHistogram(hPredictedPositionDiamondHit);
	delete hPredictedPositionDiamondHit;


	TH2F *hPredictedPositionDiamondHitFidCut = HistogrammSaver::CreateScatterHisto("hPredictedPositionDiamondHitFidCut",vecYPredictedDiamondHit,vecXPredictedDiamondHit,512);
	hPredictedPositionDiamondHitFidCut->SetTitle("Predicted Position of Silicon Track into the diamond Plane with one diamond Cluster and in FidCut");
	hPredictedPositionDiamondHitFidCut->GetXaxis()->SetTitle("predicted x Position");
	hPredictedPositionDiamondHitFidCut->GetYaxis()->SetTitle("predicted y Position");
	histSaver->SaveHistogram(hPredictedPositionDiamondHitFidCut);
	delete hPredictedPositionDiamondHitFidCut;

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
		Float_t xPos = vecXPredictedDiamondHitFidCut.at(i);
		Float_t yPos = vecYPredictedDiamondHitFidCut.at(i);
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
	events++;
	if(!eventReader->isValidTrack()){
		noValidTrack++;
		return;
	}
	Int_t nDiaClusters = eventReader->getNDiamondClusters();
	if(nDiaClusters==0){
		noDiaCluster++;
		return;
	}
	if(nDiaClusters>1){
		notOneDiaCluster++;
		return;
	}
	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster()){
		saturatedCluster++;
		return;
	}
	Float_t pos = diamondCluster.getPosition(TCluster::maxValue,0);
	Int_t area = settings->getDiaDetectorAreaOfChannel(pos);
	bool isInOneArea = !(area==-1);
	Float_t fiducialValueX = eventReader->getFiducialValueX();
	Float_t fiducialValueY = eventReader->getFiducialValueY();
	bool isInFidCut = settings->getSelectionFidCuts()->isInFiducialCut(fiducialValueX,fiducialValueY);
	Int_t fidRegionIndex = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY)-1;

	bool isMasked = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),diamondCluster,false);
	bool isMaskedAdjacentChannels = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),diamondCluster,true);

	if(isMasked){
		maskedClusters++;
		return;
	}
	if(!isInOneArea){
		cout<<nEvent<<" - "<<fiducialValueX<<"/"<<fiducialValueY<<" = "<<isInFidCut<<" \t"<< pos <<" -> "<<area<<" "<<isInOneArea<<endl;
		hitNotInOneArea++;
		return;
	}
	if(!isInFidCut){
		notInFidCut++;
		return;
	}
	hAreaVsFidCut ->Fill(area,fidRegionIndex);
	if(area!=fidRegionIndex){
		areaAndFidCutDoNotAgree++;
		return;
	}
	if(predictedPosition)
		delete predictedPosition;
	predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	Float_t xPos = predictedPosition->getPositionX();
	Float_t yPos = predictedPosition->getPositionY();
	Float_t chi2X = predictedPosition->getChi2X();
	Float_t chi2Y = predictedPosition->getChi2Y();
	Float_t chi2 = TMath::Max(chi2X,chi2Y);

	hChi2X->Fill(chi2X);
	hChi2Y->Fill(chi2Y);
	hChi2XY->Fill(chi2X,chi2Y);
	hChi2->Fill(chi2);
	if(chi2>settings->getChi2Cut3D()){
		tooHighChi2++;
		return;
	}
	if(xPos<-500||yPos<-500){
		invalidPositionPrediction++;
		return;
	}
	vecXPredicted.push_back(xPos);
	vecYPredicted.push_back(yPos);
	vecXFidCut.push_back(fiducialValueX);
	vecYFidCut.push_back(fiducialValueY);
	if(verbosity>5)cout<<nEvent<<":  "<<xPos<<"/"<<yPos<<" with "<<chi2<<" chi2.\t"<<nDiaClusters<<" "<<isInFidCut<<endl;

	if (nDiaClusters==1){
		hValidSiliconAndOneDiamondHit ->Fill(xPos,yPos);
		if(!isMasked)
			hValidSiliconAndOneDiamondHitNotMasked -> Fill(xPos,yPos);
		if(!isMaskedAdjacentChannels)
			hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels -> Fill(xPos,yPos);
		if(isInOneArea)
			hValidSiliconAndOneDiamondHitInOneArea-> Fill(xPos,yPos);
		if (area==fidRegionIndex)
			hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->Fill(xPos,yPos);
	}
	vecXPredictedDiamondHit.push_back(xPos);
	vecYPredictedDiamondHit.push_back(yPos);

	vecXPredictedDiamondHitFidCut.push_back(xPos);
	vecYPredictedDiamondHitFidCut.push_back(yPos);

	vecPHDiamondHit.push_back(diamondCluster.getCharge());
	//59,88,28
	Int_t highCh = diamondCluster.getHighestSignalChannel();
	if(highCh != pos)
		cout<<"High Ch != pos: "<<highCh<<", "<<pos<<endl;
	if(highCh == 59 ||highCh==88||highCh==28)
		vecSpecialChannelMetricPos.push_back(xPos);
			//if(eventReader->useForAnalysis()||eventReader->useForAlignment()){

	UInt_t clustSize = diamondCluster.size();

	if(clustSize>8) clustSize=8;
	vecChannelMetricPos.push_back(xPos);
	vecChannel.push_back(highCh);
	if(clustSize<=2&&nDiaClusters==1&&area==fidRegionIndex){
		//		if(fidRegionIndex<hChargeVsFidX.size()&&fidRegionIndex>=0){
//			hChargeVsFidX[fidRegionIndex]->Fill(charge,fiducialValueX);
//			hChargeVsFidY[fidRegionIndex]->Fill(charge,fiducialValueY);
//		}
//		else{
//			cout<<"fidRegion not valid: "<<fidRegionIndex<<" "<<fiducialValueX<<"/"<<fiducialValueY<<endl;
//			settings->getSelectionFidCuts()->Print();
//		}
//		hChargeVsFidCut->Fill(fiducialValueX,fiducialValueY,charge);
//		hFidCutXvsChannelPos->Fill(fiducialValueX,pos);
//		histoLandauDistribution2D_unmasked->Fill(charge,pos);
//		bool isBorderSeedCluster = settings->hasBorderSeed(TPlaneProperties::getDetDiamond(),cluster);
//		bool isBorderHitCluster = settings->hasBorderHit(TPlaneProperties::getDetDiamond(),cluster);
//		if (!isBorderSeedCluster)
//			histoLandauDistribution2DNoBorderSeed_unmasked->Fill(charge,pos);
//		if (!isBorderHitCluster)
//			histoLandauDistribution2DNoBorderHit_unmasked->Fill(charge,pos);
//		bool isMaskedCluster = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),cluster,false);
//		if(!isMaskedCluster){
//			histoLandauDistribution2D->Fill(charge,pos);
//			if (!isBorderSeedCluster)
//				histoLandauDistribution2DNoBorderSeed->Fill(charge,pos);
//			if (!isBorderHitCluster)
//				histoLandauDistribution2DNoBorderHit->Fill(charge,pos);
//		}
	}
	validEvents++;
}
