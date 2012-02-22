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
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());
	
	sys->cd(runString.str().c_str());
	stringstream  filepath, alignmentFileName;
	filepath.str("");
	alignmentFileName.str("");
	filepath<<"selectionData."<<runNumber<<".root";
	alignmentFileName << "alignment." << runNumber << ".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	
	tracking = new TTracking(filepath.str(),alignmentFileName.str(),runNumber);
	
	// TODO: load settings!!!
	
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("transparentAnalysis");
	sys->cd("transparentAnalysis");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	sys->cd("..");
	
	
	// TODO: move these setting to the proper place
	subjectDetector = 8;
	subjectPlane = subjectDetector/2;
	if (subjectDetector%2 == 0) {
		subjectDetectorCoordinate = TPlane::X_COR;
	}
	else {
		subjectDetectorCoordinate = TPlane::Y_COR;
	}
	for (int i = 0; i < 4; i++) {
		refPlanes.push_back(i);
	}
	clusterCalcMode = TCluster::highest2Centroid;
	verbosity = 0;
	
	
	initHistograms();
//	this->seedSigma=seedSigma;
//	this->hitSigma=hitSigma;
	cout<<"end initialise"<<endl;
	
	
	
	
}

TTransparentAnalysis::~TTransparentAnalysis() {
	// TODO Auto-generated destructor stub
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
	for (nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
//		if (verbosity > 4) cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << endl;
		tracking->LoadEvent(nEvent);
//		if (tracking->isValidTrack() == 0) {
		if (tracking->useForAnalysis() == 0) {
			if (verbosity > 6) printEvent();
			noValidTrack++;
			continue;
		}
		if (tracking->isInFiducialCut() == 0) {
			noFidCutRegion++;
			continue;
		}
		transparentClusters.clear();
		positionPrediction = tracking->predictPosition(subjectPlane,refPlanes,false);
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
		this->positionInDetSystem = tracking->getPositionInDetSystem(subjectDetector, this->predXPosition, this->predYPosition);
//		if (verbosity > 4) cout << "position in det system:\t" << this->positionInDetSystem << endl;
//		if (verbosity > 4)
//			cout << "clustered analysis strip position:\t" << tracking->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
		if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) continue;
		for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(subjectDetector, this->positionInDetSystem, clusterSize));
		}
		nAnalyzedEvents++;
		this->fillHistograms();
		if (verbosity > 4) printEvent();
	}
	this->saveHistograms();
	this->printCutFlow();
}

bool TTransparentAnalysis::checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// TODO: cut flow!!
	UInt_t centerChannel = (UInt_t)centerPosition;
	int direction = 1;
	if (centerPosition-(int)centerPosition<0.5) {
		direction = -1;
	}
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
		if (tracking->isSaturated(det, currentChannel) == true) {
			if (verbosity > 5) cout << "channel " << currentChannel << " has saturated.." << endl;
			saturatedChannel++;
			return false;
		}
	}
	return true;
}

// TODO: avoid wrong channel numbers (>128, <0)
TCluster TTransparentAnalysis::makeTransparentCluster(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	UInt_t centerChannel = (UInt_t)centerPosition;
	int direction = 1;
	if (centerPosition-(int)centerPosition<0.5) {
		direction = -1;
	}
	TCluster transparentCluster = TCluster(tracking->getEvent_number(), det, 0, 0, TPlaneProperties::getNChannels(det));
	int currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		transparentCluster.addChannel(currentChannel, tracking->getSignal(det,currentChannel), tracking->getSignalInSigma(det,currentChannel), tracking->getAdcValue(det,currentChannel), tracking->isSaturated(det,currentChannel), settings->isDet_channel_screened(det,currentChannel));
	}
	return transparentCluster;
}

void TTransparentAnalysis::setSettings(TSettings* settings){
	this->settings=settings;
}

void TTransparentAnalysis::initHistograms() {
	UInt_t bins=100;
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		// TODO: take care of histogram names and bins!!
		stringstream histNameLaundau, histNameEta, histNameResidual;
		// TODO: histogram naming!!
		histNameLaundau << "hDiaTranspAnaPulseHightOf" << clusterSize+1 << "Strips";
		histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize+1 << "Strips";
		histNameResidual << "hDiaTranspAnaResidual2HighestIn" << clusterSize+1 << "StripsMinusPred";
		hLaundau.push_back(new TH1F(histNameLaundau.str().c_str(),histNameLaundau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hEta.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
		hResidual.push_back(new TH1F(histNameResidual.str().c_str(),histNameResidual.str().c_str(),bins,-5.,5.));
	}
}

void TTransparentAnalysis::fillHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		hLaundau[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge());
		hEta[clusterSize]->Fill(this->transparentClusters[clusterSize].getEta());
		hResidual[clusterSize]->Fill(this->getResidual(this->transparentClusters[clusterSize]));
	}
}

void TTransparentAnalysis::saveHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		histSaver->SaveHistogram(hLaundau[clusterSize],0);
		histSaver->SaveHistogram(hEta[clusterSize],0);
		histSaver->SaveHistogram(hResidual[clusterSize],1);
	}
}

// TODO: call TTransparentAnalysis::deleteHistograms
void TTransparentAnalysis::deleteHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		delete hLaundau[clusterSize];
		delete hEta[clusterSize];
		delete hResidual[clusterSize];
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
	if (tracking->useForAnalysis() == 0) {
		cout << "this track is not used for the analysis.." << endl;
		return;
	}
	cout << "predicted pos in lab system:\t" << this->predPosition << "\tpredicted perp position:\t" << this->predPerpPosition << endl;
	cout << "predicted pos in det system:\t" << this->positionInDetSystem << endl;
	cout << "clustered analysis position in lab system:\t" << tracking->getStripXPosition(subjectPlane,this->predPerpPosition,clusterCalcMode) << endl;
	cout << "clustered analysis position in det system:\t" << tracking->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
	if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) {
		cout << "this track did not pass the check.." << endl;
		return;
	}
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		cout << "transparent cluster of size " << clusterSize+1 << ":" << endl;
		cout << "\tpulse height:\t" << this->transparentClusters[clusterSize].getCharge() << endl;
		cout << "\teta:\t" << this->transparentClusters[clusterSize].getEta() << endl;
		cout << "\tresidual:\t" << this->getResidual(this->transparentClusters[clusterSize]) << endl;
		cout << "\tcluster pos in det system:\t" << this->transparentClusters[clusterSize].getPosition(this->clusterCalcMode) << endl;
		cout << "\tcluster pos in lab system:\t" << tracking->getPositionOfCluster(subjectDetector, this->transparentClusters[clusterSize], this->predPerpPosition, this->clusterCalcMode) << endl;
	}
	return;
}

Float_t TTransparentAnalysis::getResidual(TCluster cluster) {
	return tracking->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,this->clusterCalcMode)-this->predPosition;
}

