//
//  TTransparentAnalysis.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 05.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TTransparentAnalysis.hh"

TTransparentAnalysis::TTransparentAnalysis(int runNumber, TSettings settings) {
	cout<<"**********************************************************"<<endl;
	cout<<"********TTransparentAnalysis::TTransparentAnalysis********"<<endl;
	cout<<"**********************************************************"<<endl;
	// TODO Auto-generated constructor stub
	sys = gSystem;
	setSettings(&settings);
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());
	
	sys->cd(runString.str().c_str());
	stringstream  filepath, alignmentFileName;
	filepath.str("");
	alignmentFileName.str("");
	filepath<<"clusterData."<<runNumber<<".root";
	alignmentFileName << "alignment." << runNumber << ".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	
	tracking = new TTracking(filepath.str(),alignmentFileName.str());
	
	// TODO: load settings!!!
	
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("transparentAnalysis");
	sys->cd("transparentAnalysis");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	initHistograms();
//	this->seedSigma=seedSigma;
//	this->hitSigma=hitSigma;
	cout<<"end initialise"<<endl;
	
	
	
	
}

TTransparentAnalysis::~TTransparentAnalysis() {
	// TODO Auto-generated destructor stub
}

void TTransparentAnalysis::analyze(int nEvents, int startEvent) {
	cout<<"\n\n******************************************\n";
	cout<<    "******Start Transparent Analysis...*******\n";
	cout<<"******************************************\n\n"<<endl;
	
	
	// TODO: move these setting to the proper place
	subjectDetector = 8;
	for (int i = 0; i < 8; i++) {
		refPlanes.push_back(i);
	}
	transparentMaxClusterSize = 10;
	
	
	
	nAnalyzedEvents = 0;
	regionNotOnPlane = 0;
	saturatedChannel = 0;
	screenedChannel = 0;
	noValidTrack = 0;
	
	for (int nEvent = startEvent; nEvent < nEvents; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
		tracking->LoadEvent(nEvent);
		if (tracking->isValidTrack() == 0) {
			noValidTrack++;
			continue;
		}
		transparentClusters.clear();
		positionPrediction = tracking->predictPosition(subjectDetector,refPlanes);
		this->predXPosition = positionPrediction->getPositionX();
		this->predYPosition = positionPrediction->getPositionY();
		// TODO: position in det system
//		this->positionInDetSystem = getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred)
		if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, transparentMaxClusterSize) == false) continue;
		for (UInt_t clusterSize = 1; clusterSize < transparentMaxClusterSize+1; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(subjectDetector, this->positionInDetSystem, clusterSize));
		}
		nAnalyzedEvents++;
		this->fillHistograms();
	}
}

bool TTransparentAnalysis::checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// TODO: cut flow!!
	UInt_t centerChannel = (UInt_t)centerPosition;
	int direction = 1;
	if (centerPosition-(int)centerPosition<0.5) {
		direction = -1;
	}
	int currentChannel = centerChannel;
	for (int iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		if (currentChannel < 0) {
			regionNotOnPlane++;
			return false;
		}
		if (currentChannel > TPlaneProperties::getNChannels(det)-1) {
			regionNotOnPlane++;
			return false;
		}
		if (this->settings->getDet_channel_screen(det).isScreened(currentChannel) == true) {
			screenedChannel++;
			return false;
		}
		if (tracking->isSaturated(det, currentChannel) == true) {
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
	
}

void TTransparentAnalysis::fillHistograms() {
	
}
