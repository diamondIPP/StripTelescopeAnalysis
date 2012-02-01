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
	setSettings(settings);
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
	for (int i = 0; i < 8; i++) {
		siliconPlanes.push_back(i);
	}
	
	for (int nEvent = startEvent; nEvent < nEvents; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
		tracking->LoadEvent(nEvent);
		if (tracking->isValidTrack() == 0) continue;	// TODO: cut flow
		transparentClusters.clear();
		positionPrediction = tracking->predictPosition(8,siliconPlanes);
		for (UInt_t clusterSize = 1; clusterSize < 11; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(8, positionPrediction->getPositionX(), clusterSize, 1)); // TODO: check centerChannel and direction!!
		}
	}
}

// TODO: avoid wrong channel numbers (>128, <0)
TCluster TTransparentAnalysis::makeTransparentCluster(UInt_t det, UInt_t centerChannel, UInt_t clusterSize, int direction) {
	direction = direction / TMath::Abs(direction);
	TCluster transparentCluster = TCluster(tracking->getEvent_number(), det, 0, 0, 128); // TODO: modify for to use for telescope planes
	int currentChannel = centerChannel;
	for (UInt_t iChannel; iChannel < clusterSize; iChannel++) {
		currentChannel -= direction * iChannel;
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