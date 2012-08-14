//
//  TTransparentAnalysis.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 05.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TTransparentAnalysis.hh"

TTransparentAnalysis::TTransparentAnalysis(TSettings* settings) {
	cout<<"**********************************************************"<<endl;
	cout<<"********TTransparentAnalysis::TTransparentAnalysis********"<<endl;
	cout<<"**********************************************************"<<endl;
	// TODO Auto-generated constructor stub
	sys = gSystem;
	setSettings(settings);
	UInt_t runNumber =settings->getRunNumber();
	
	settings->goToAlignmentRootDir();
	eventReader = new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),runNumber);
	// TODO: load settings!!!
	
	histSaver = new HistogrammSaver();
//	settings->goToTransparentAnalysisDir();
	histSaver->SetPlotsPath(settings->getTransparentAnalysisDir());
	histSaver->SetRunNumber(settings->getRunNumber());
	htmlTransAna = new THTMLTransparentAnalysis(settings);
	htmlTransAna->setFileGeneratingPath(settings->getTransparentAnalysisDir());
	clustering = new TClustering(settings);
	
	
	// TODO: move these setting to the proper place
	subjectDetector = 8;
	subjectPlane = subjectDetector/2;
	if (subjectDetector%2 == 0) {
		subjectDetectorCoordinate = TPlaneProperties::X_COR;
	}
	else {
		subjectDetectorCoordinate = TPlaneProperties::Y_COR;
	}
	for (int i = 0; i < 4; i++) {
		refPlanes.push_back(i);
	}
	clusterCalcMode = TCluster::highest2Centroid;
	verbosity = 0;
	
//  settings->goToAlignmentRootDir();
	initHistograms();
//	this->seedSigma=seedSigma;
//	this->hitSigma=hitSigma;
	cout<<"end initialise"<<endl;
	
}

TTransparentAnalysis::~TTransparentAnalysis() {
	// TODO Auto-generated destructor stub
	cout<<"\n\nClosing TTransparentAnalysis"<<endl;
	fitHistograms();
	saveHistograms();
	deleteHistograms();
	deleteFits();
	
	vector<vector <Float_t> > meanPulseHeights;
	vector<vector <Float_t> > mpPulseHeights;
	vector<vector <pair <Float_t,Float_t> > > resolutions;
	
	mpPulseHeights.push_back(vecMPLandau);
	mpPulseHeights.push_back(vecMPLandau2Highest);
	meanPulseHeights.push_back(vecMeanLandau);
	meanPulseHeights.push_back(vecMeanLandau2Highest);
	resolutions.push_back(vecResidualChargeWeighted);
	resolutions.push_back(vecResidualHighest2Centroid);
	resolutions.push_back(vecResidualEtaCorrected);
	
	htmlTransAna->createPulseHeightPlots(meanPulseHeights, mpPulseHeights);
	htmlTransAna->createResolutionPlots(resolutions);
	htmlTransAna->createEtaPlots();
	htmlTransAna->createEtaIntegrals();
	htmlTransAna->generateHTMLFile();
	if (eventReader!=0) delete eventReader;
	if (histSaver!=0) delete histSaver;
	if (htmlTransAna!=0) delete htmlTransAna;
	if (clustering!=0) delete clustering;
	settings->goToOutputDir();
}

void TTransparentAnalysis::analyze(UInt_t nEvents, UInt_t startEvent) {
	cout<<"\n\n******************************************\n";
	cout<<    "******Start Transparent Analysis...*******\n";
	cout<<"******************************************\n\n"<<endl;
	nAnalyzedEvents = 0;
	regionNotOnPlane = 0;
	saturatedChannel = 0;
	screenedChannel = 0;
	noValidTrack = 0;
	cout<<"Current Dir: "<<sys->pwd()<<endl;
	if (nEvents+startEvent >= eventReader->GetEntries()) {
		cout << "only "<<eventReader->GetEntries()<<" in tree!\n";
		nEvents = eventReader->GetEntries()-startEvent;
	}
	for (nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
//		if (verbosity > 4) cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << eventReader<<endl;
		eventReader->LoadEvent(nEvent);
//		if (eventReader->isValidTrack() == 0) {
		if (eventReader->useForAnalysis() == 0) {
			if (verbosity > 6) printEvent();
			noValidTrack++;
			continue;
		}
		if (eventReader->isInFiducialCut() == 0) {
			noFidCutRegion++;
			continue;
		}
		transparentClusters.clear();
		
		this->predictPositions();
		
		
		/*
		positionPrediction = eventReader->predictPosition(subjectPlane,refPlanes,false);
		this->predXPosition = positionPrediction->getPositionX();
		this->predYPosition = positionPrediction->getPositionY();
//		if (verbosity > 4) cout << "predicted x position:\t" << this->predXPosition << "\ty position:\t" << this->predYPosition << endl;
		if (subjectDetector%2 == 0) {
			this->predPerpPosition = this->predYPosition;
			this->predPosition = this->predXPosition;
		}
		else {
			this->predPerpPosition = this->predXPosition;
			this->predPosition = this->predYPosition;
		}
		// TODO: position in det system
		this->positionInDetSystem = eventReader->getPositionInDetSystem(subjectDetector, this->predXPosition, this->predYPosition);
//		if (verbosity > 4) cout << "position in det system:\t" << this->positionInDetSystem << endl;
//		if (verbosity > 4)
//			cout << "clustered analysis strip position:\t" << eventReader->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
		*/
		
		
		if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) continue;
		for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(subjectDetector, this->positionInDetSystem, clusterSize));
		}
		nAnalyzedEvents++;
		this->fillHistograms();
		if (verbosity > 4) printEvent();
		
		
		// save clusters for eta corrected analysis
		vecTransparentClusters.push_back(transparentClusters);
		eventNumbers.push_back(nEvent);
	}
	this->printCutFlow();
	createEtaIntegrals();
	calcEtaCorrectedResiduals();
}

void TTransparentAnalysis::calcEtaCorrectedResiduals() {
	cout << "\n\ncalculating eta corrected residuals..\n";
	if (eventNumbers.size() != vecTransparentClusters.size()) {
		cout << "now we are in deep trouble!! size of eventNumbers and transparentClusters do not match!" << endl;
		return;
	}
	if (vecTransparentClusters.size()==0 || eventNumbers.size()==0) {
		cout << "oh boy.. you didn't run the analysis!" << endl;
	}
	for (UInt_t iEvent = 0; iEvent < eventNumbers.size(); iEvent++) {
		TRawEventSaver::showStatusBar(iEvent,eventNumbers.size(),100);
		nEvent = eventNumbers.at(iEvent);
		eventReader->LoadEvent(nEvent);
		this->predictPositions();
		for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
			// TODO: TCluster::corEta or TCluster::eta???
			hResidualEtaCorrected[clusterSize]->Fill(getResidual(vecTransparentClusters.at(iEvent).at(clusterSize),TCluster::corEta,hEtaIntegrals[clusterSize]));
//			if (clusterSize == 1) printCluster(vecTransparentClusters.at(iEvent).at(clusterSize));
		}
	}
}

void TTransparentAnalysis::predictPositions() {
	positionPrediction = eventReader->predictPosition(subjectPlane,refPlanes,false);
	this->predXPosition = positionPrediction->getPositionX();
	this->predYPosition = positionPrediction->getPositionY();
	//		if (verbosity > 4) cout << "predicted x position:\t" << this->predXPosition << "\ty position:\t" << this->predYPosition << endl;
	if (subjectDetector%2 == 0) {
		this->predPerpPosition = this->predYPosition;
		this->predPosition = this->predXPosition;
	}
	else {
		this->predPerpPosition = this->predXPosition;
		this->predPosition = this->predYPosition;
	}
	// TODO: position in det system
	this->positionInDetSystem = eventReader->getPositionInDetSystem(subjectDetector, this->predXPosition, this->predYPosition);
	//		if (verbosity > 4) cout << "position in det system:\t" << this->positionInDetSystem << endl;
	//		if (verbosity > 4)
	//			cout << "clustered analysis strip position:\t" << eventReader->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
}

bool TTransparentAnalysis::checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// get channel and direction for clustering
	UInt_t centerChannel;
	int direction;
	direction = getSignedChannelNumber(centerPosition);
	centerChannel = TMath::Abs(direction);
	if (direction < 0) direction = -1;
	else direction = 1;
	
	// check predicted cluster channels
	int currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		if (currentChannel < 0) {
			if (verbosity > 5) cout << "channel " << currentChannel << " is not on this detector.." << endl;
			regionNotOnPlane++;
			return false;
		}
		if (currentChannel > TPlaneProperties::getNChannels(det)-1) {
			if (verbosity > 5) cout << "channel " << currentChannel << " is not on this detector.." << endl;
			regionNotOnPlane++;
			return false;
		}
		if (this->settings->getDet_channel_screen(det).isScreened(currentChannel) == true) {
			if (verbosity > 5) cout << "channel " << currentChannel << " is screened.." << endl;
			screenedChannel++;
			return false;
		}
		if (eventReader->isSaturated(det, currentChannel) == true) {
			if (verbosity > 5) cout << "channel " << currentChannel << " has saturated.." << endl;
			saturatedChannel++;
			return false;
		}
	}
	return true;
}

// TODO: avoid wrong channel numbers (>128, <0)
TCluster TTransparentAnalysis::makeTransparentCluster(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// get channel and direction for clustering
	UInt_t centerChannel;
	int direction;
	direction = getSignedChannelNumber(centerPosition);
//	cout << "centerPosition: " << centerPosition << "\tdirection: " << direction << endl;
	centerChannel = TMath::Abs(direction);
	if (direction < 0) direction = -1;
	else direction = 1;
	
	// make cluster
	TCluster transparentCluster = TCluster(eventReader->getEvent_number(), det, -99, -99, TPlaneProperties::getNChannels(det));
	int currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		transparentCluster.addChannel(currentChannel, eventReader->getRawSignal(det,currentChannel), eventReader->getRawSignalInSigma(det,currentChannel), eventReader->getAdcValue(det,currentChannel), eventReader->isSaturated(det,currentChannel), settings->isDet_channel_screened(det,currentChannel));
	}
	return transparentCluster;
}

int TTransparentAnalysis::getSignedChannelNumber(Float_t position) {
	if (position < 0) return -9999;
	UInt_t channel = 0;
	int direction;
	if (position-(int)position < 0.5) {
		channel = (UInt_t)position;
		direction = 1;
	}
	else {
		channel = (UInt_t)position+1;
		direction = -1;
	}
	return direction * channel;
}

void TTransparentAnalysis::setSettings(TSettings* settings) {
	this->settings=settings;
}

void TTransparentAnalysis::initHistograms() {
	UInt_t bins=100;
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		// TODO: take care of histogram names and bins!!
		stringstream histNameLaundau, histNameLaundau2Highest, histNameEta, histNameResidualChargeWeighted, histNameResidualHighest2Centroid, histNameResidualEtaCorrected;
		// TODO: histogram naming!!
		histNameLaundau << "hDiaTranspAnaPulseHeightOf" << clusterSize+1 << "Strips";
		histNameLaundau2Highest << "hDiaTranspAnaPulseHeightOf2HighestIn" << clusterSize+1 << "Strips";
		histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize+1 << "Strips";
		histNameResidualChargeWeighted << "hDiaTranspAnaResidualChargeWeightedIn" << clusterSize+1 << "StripsMinusPred";
		histNameResidualHighest2Centroid << "hDiaTranspAnaResidualHighest2CentroidIn" << clusterSize+1 << "StripsMinusPred";
		histNameResidualEtaCorrected << "hDiaTranspAnaResidualEtaCorrectedIn" << clusterSize+1 << "StripsMinusPred";
		hLaundau.push_back(new TH1F(histNameLaundau.str().c_str(),histNameLaundau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hLaundau2Highest.push_back(new TH1F(histNameLaundau2Highest.str().c_str(),histNameLaundau2Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hEta.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
		hResidualChargeWeighted.push_back(new TH1F(histNameResidualChargeWeighted.str().c_str(),histNameResidualChargeWeighted.str().c_str(),bins,-5.,5.));
		hResidualHighest2Centroid.push_back(new TH1F(histNameResidualHighest2Centroid.str().c_str(),histNameResidualHighest2Centroid.str().c_str(),bins,-5.,5.));
		hResidualEtaCorrected.push_back(new TH1F(histNameResidualEtaCorrected.str().c_str(),histNameResidualEtaCorrected.str().c_str(),bins,-5.,5.));
	}
	hLaundauMean = new TH1F("hDiaTranspAnaPulseHeightMean","hDiaTranspAnaPulseHeightMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLaundauMP = new TH1F("hDiaTranspAnaPulseHeightMP","hDiaTranspAnaPulseHeightMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLaundau2HighestMean = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMean","hDiaTranspAnaPulseHeightOf2HighestMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLaundau2HighestMP = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMP","hDiaTranspAnaPulseHeightOf2HighestMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hPredictedPositionInStrip = new TH1F("hPredictedPositionInStrip","hPredictedPositionInStrip",2,-1.5,1.5);
}

void TTransparentAnalysis::fillHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		hLaundau[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge());
		hLaundau2Highest[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge(2,false));
		hEta[clusterSize]->Fill(this->transparentClusters[clusterSize].getEta());
//		if (clusterSize == 1 /* && this->transparentClusters[clusterSize].getCharge() != this->transparentClusters[clusterSize].getCharge(2,false)*/) printCluster(this->transparentClusters[clusterSize]);
//		if (clusterSize > 0 && this->transparentClusters[clusterSize].getEta() < 0) printCluster(this->transparentClusters[clusterSize]);
		// TODO: why is the eta distribution for 2 channel clusters more symmetric than for 3 and more channel clusters?
//		if (clusterSize == 2 && this->transparentClusters[clusterSize-1].getEta() != this->transparentClusters[clusterSize].getEta()) {
//			if (this->transparentClusters[clusterSize-1].getHighestSignalChannel()!=this->transparentClusters[clusterSize].getHighestSignalChannel()
//				&&
//				this->transparentClusters[clusterSize-1].getHighest2Centroid()!=this->transparentClusters[clusterSize].getHighest2Centroid()) {
//			printCluster(this->transparentClusters[clusterSize-1]);
//			printCluster(this->transparentClusters[clusterSize]);
//			}
//		}
		if (clusterSize+1 != transparentClusters[clusterSize].getClusterSize()) {
			cout << "wrong cluster size!" << endl;
			cout << "clusterSize+1 = " << clusterSize+1 << "\ttransparentClusters[clusterSize].getClusterSize() = " << transparentClusters[clusterSize].getClusterSize() << endl;
		}
		hResidualChargeWeighted[clusterSize]->Fill(this->getResidual(this->transparentClusters[clusterSize],TCluster::chargeWeighted));
		hResidualHighest2Centroid[clusterSize]->Fill(this->getResidual(this->transparentClusters[clusterSize],TCluster::highest2Centroid));
	}
//	hPredictedPositionInStrip->Fill();
}

TF1* TTransparentAnalysis::doGaussFit(TH1F *histo) {
//	TH1* histo = (TH1*)htemp->Clone();
	if (histo->GetEntries()==0) return 0;
	TF1* histofitx = new TF1("histofitx","gaus",histo->GetMean()-2*histo->GetRMS(),histo->GetMean()+2*histo->GetRMS());
	histofitx->SetLineColor(kBlue);
	histo->Fit(histofitx,"rq");
	return histofitx;
}

void TTransparentAnalysis::createEtaIntegrals() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		stringstream histName;
		histName << "hDiaTranspAnaEtaIntegral2HighestIn"<<clusterSize+1<<"Strips";
		hEtaIntegrals.push_back(clustering->createEtaIntegral(hEta[clusterSize], histName.str()));
	}
}

void TTransparentAnalysis::fitHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		// fit histograms
		fitLandau.push_back(landauGauss->doLandauGaussFit(hLaundau[clusterSize]));
		fitLandau2Highest.push_back(landauGauss->doLandauGaussFit(hLaundau2Highest[clusterSize]));
		fitResidualChargeWeighted.push_back(doGaussFit(hResidualChargeWeighted[clusterSize]));
		fitResidualHighest2Centroid.push_back(doGaussFit(hResidualHighest2Centroid[clusterSize]));
		fitResidualEtaCorrected.push_back(doGaussFit(hResidualEtaCorrected[clusterSize]));
		
		// save fit parameters
		vecMPLandau.push_back(fitLandau[clusterSize]->GetParameter(1));
		vecMPLandau2Highest.push_back(fitLandau2Highest[clusterSize]->GetParameter(1));
		hLaundauMP->SetBinContent(clusterSize+1,fitLandau[clusterSize]->GetParameter(1));
		hLaundau2HighestMP->SetBinContent(clusterSize+1,fitLandau2Highest[clusterSize]->GetParameter(1));
		vecMeanLandau.push_back(hLaundau[clusterSize]->GetMean());
		vecMeanLandau2Highest.push_back(hLaundau2Highest[clusterSize]->GetMean());
		hLaundauMean->SetBinContent(clusterSize+1,hLaundau[clusterSize]->GetMean());
		hLaundau2HighestMean->SetBinContent(clusterSize+1,hLaundau2Highest[clusterSize]->GetMean());
		pair <Float_t,Float_t> tempPair;
		if (fitResidualChargeWeighted[clusterSize]!=0) {
			tempPair.first = fitResidualChargeWeighted[clusterSize]->GetParameter(1);
			tempPair.second = fitResidualChargeWeighted[clusterSize]->GetParameter(2);
		}
		else {
			tempPair.first = 0;
			tempPair.second = 0;
		}
		vecResidualChargeWeighted.push_back(tempPair);
		if (fitResidualHighest2Centroid[clusterSize]!=0) {
			tempPair.first = fitResidualHighest2Centroid[clusterSize]->GetParameter(1);
			tempPair.second = fitResidualHighest2Centroid[clusterSize]->GetParameter(2);
		}
		else {
			tempPair.first = 0;
			tempPair.second = 0;
		}
		vecResidualHighest2Centroid.push_back(tempPair);
		if (fitResidualEtaCorrected[clusterSize]!=0) {
			tempPair.first = fitResidualEtaCorrected[clusterSize]->GetParameter(1);
			tempPair.second = fitResidualEtaCorrected[clusterSize]->GetParameter(2);
		}
		else {
			tempPair.first = 0;
			tempPair.second = 0;
		}
		vecResidualEtaCorrected.push_back(tempPair);
	}
	hLaundauMean->Scale(1./hLaundauMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLaundauMP->Scale(1./hLaundauMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLaundau2HighestMean->Scale(1./hLaundau2HighestMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLaundau2HighestMP->Scale(1./hLaundau2HighestMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
}

void TTransparentAnalysis::saveHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		histSaver->SaveHistogramWithFit(hLaundau[clusterSize],fitLandau[clusterSize]);
		histSaver->SaveHistogramWithFit(hLaundau2Highest[clusterSize],fitLandau2Highest[clusterSize]);
		histSaver->SaveHistogram(hEta[clusterSize],0);
		histSaver->SaveHistogramWithFit(hResidualChargeWeighted[clusterSize],fitResidualChargeWeighted[clusterSize]);
		histSaver->SaveHistogramWithFit(hResidualHighest2Centroid[clusterSize],fitResidualHighest2Centroid[clusterSize]);
		histSaver->SaveHistogramWithFit(hResidualEtaCorrected[clusterSize],fitResidualEtaCorrected[clusterSize]);
		histSaver->SaveHistogram(hEtaIntegrals[clusterSize],0);
	}
	histSaver->SaveHistogram(hLaundauMean);
	histSaver->SaveHistogram(hLaundauMP);
	histSaver->SaveHistogram(hLaundau2HighestMean);
	histSaver->SaveHistogram(hLaundau2HighestMP);
}

void TTransparentAnalysis::deleteHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		delete hLaundau[clusterSize];
		delete hLaundau2Highest[clusterSize];
		delete hEta[clusterSize];
		delete hResidualChargeWeighted[clusterSize];
		delete hResidualHighest2Centroid[clusterSize];
	}
	delete hLaundauMean;
	delete hLaundauMP;
	delete hLaundau2HighestMean;
	delete hLaundau2HighestMP;
}

void TTransparentAnalysis::deleteFits() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		delete fitLandau[clusterSize];
		delete fitLandau2Highest[clusterSize];
		delete fitResidualChargeWeighted[clusterSize];
		delete fitResidualHighest2Centroid[clusterSize];
	}
}

void TTransparentAnalysis::printCutFlow() {
	cout << "\n\n\n";
	cout << "TTransparentAnalysis has analyzed " << nAnalyzedEvents << " events." << endl;
	cout << "region not on plane\t" << regionNotOnPlane << endl;
	cout << "saturated channel\t" << saturatedChannel << endl;
	cout << "screened channel\t" << screenedChannel << endl;
	cout << "track not valid\t" << noValidTrack << endl;
	cout << "track not in fidutial cut region\t" << noFidCutRegion << endl;
}

void TTransparentAnalysis::printEvent() {
	cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << endl;
	if (eventReader->useForAnalysis() == 0) {
		cout << "this track is not used for the analysis.." << endl;
		return;
	}
	cout << "predicted pos in lab system:\t" << this->predPosition << "\tpredicted perp position:\t" << this->predPerpPosition << endl;
	cout << "predicted pos in det system:\t" << this->positionInDetSystem << endl;
	cout << "clustered analysis position in lab system:\t" << eventReader->getStripXPosition(subjectPlane,this->predPerpPosition,clusterCalcMode) << endl;
	cout << "clustered analysis position in det system:\t" << eventReader->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
	if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) {
		cout << "this track did not pass the check.." << endl;
		return;
	}
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		cout << "transparent cluster of size " << clusterSize+1 << ":" << endl;
		cout << "\tpulse height:\t" << this->transparentClusters[clusterSize].getCharge() << endl;
		cout << "\teta:\t" << this->transparentClusters[clusterSize].getEta() << endl;
		cout << "\tresidual:\t" << this->getResidual(this->transparentClusters[clusterSize],this->clusterCalcMode) << endl;
		cout << "\tcluster pos in det system:\t" << this->transparentClusters[clusterSize].getPosition(this->clusterCalcMode) << endl;
		cout << "\tcluster pos in lab system:\t" << eventReader->getPositionOfCluster(subjectDetector, this->transparentClusters[clusterSize], this->predPerpPosition, this->clusterCalcMode) << endl;
	}
	return;
}

void TTransparentAnalysis::printCluster(TCluster cluster) {
	cout << "\n--- event " << nEvent;
	cout << "\n\tcluster size: " << cluster.getClusterSize();
	cout << "\n\tcharge: " << cluster.getCharge(false);
	cout << "\n\tcharge of 2 highest centroid: " << cluster.getCharge(2,false);
	cout << "\n\thighest channel: " << cluster.getHighestSignalChannel();
	cout << "\n\thighest 2 centroid: " << cluster.getHighest2Centroid();
	cout << "\n\tcluster position of highest channel: " << cluster.getClusterPosition(cluster.getHighestSignalChannel());
	cout << "\n\thighest channel is seed? " << cluster.isSeed(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
	cout << "\n\thighest channel is hit? " << cluster.isHit(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
	cout << "\n\tseed sigma: " << cluster.getSeedSigma();
	cout << "\n\thit sigma: " << cluster.getHitSigma();
	cout << "\n\tpredicted channel: " << positionInDetSystem;
	cout << "\n\tpredicted position: " << predPosition;
	cout << "\n\tcharge weighted position (TCluster::chargeWeighted): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid position (TCluster::highest2Centroid): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::highest2Centroid);
	cout << "\n\tcharge weighted residual: " << getResidual(cluster,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid residual: " << getResidual(cluster,TCluster::highest2Centroid);
	cout << "\n\teta: " << cluster.getEta();
	if (hEtaIntegrals.size() != 0) {
		cout << "\n\teta corrected position (TCluster::corEta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::corEta,hEtaIntegrals[cluster.getClusterSize()]);
		cout << "\n\teta corrected position (TCluster::eta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::eta,hEtaIntegrals[cluster.getClusterSize()]);
	}
	cout << "\n\t";
	cluster.Print();
}

Float_t TTransparentAnalysis::getResidual(TCluster cluster, TCluster::calculationMode_t clusterCalculationMode, TH1F* hEtaInt) {
	return eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,clusterCalculationMode, hEtaInt)-this->predPosition;
}

