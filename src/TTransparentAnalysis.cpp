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
	
	
	// TODO: move these setting to the proper place
	subjectDetector = 8;
	for (int i = 0; i < 8; i++) {
		refPlanes.push_back(i);
	}
	clusterCalcMode = TCluster::highest2Centroid;
	
	
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
	nAnalyzedEvents = 0;
	regionNotOnPlane = 0;
	saturatedChannel = 0;
	screenedChannel = 0;
	noValidTrack = 0;
	for (int nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
		tracking->LoadEvent(nEvent);
		if (tracking->isValidTrack() == 0) {
			noValidTrack++;
			continue;
		}
		if (tracking->isInFiducialCut() == 0) {
			noFidCutRegion++;
			continue;
		}
		transparentClusters.clear();
		positionPrediction = tracking->predictPosition(subjectDetector,refPlanes);
		this->predXPosition = positionPrediction->getPositionX();
		this->predYPosition = positionPrediction->getPositionY();
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
		if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystem, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) continue;
		for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(subjectDetector, this->positionInDetSystem, clusterSize));
		}
		nAnalyzedEvents++;
		this->fillHistograms();
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
	UInt_t bins=512;
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		// TODO: take care of histogram names and bins!!
		stringstream histNameLaundau, histNameEta, histNameResidual;
		// TODO: histogram naming!!
		histNameLaundau << "hDiaTranspAnaPulseHightOf" << clusterSize << "Strips";
		histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize << "Strips";
		histNameResidual << "hDiaTranspAnaResidual2HighestIn" << clusterSize << "StripsMinusPred";
		hLaundau.push_back(new TH1F(histNameLaundau.str().c_str(),histNameLaundau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hEta.push_back(new TH1F(histNameLaundau.str().c_str(),histNameLaundau.str().c_str(),bins,0,1));
		hResidual.push_back(new TH1F("","",bins,-5.,5.));
	}
}

void TTransparentAnalysis::fillHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		hLaundau[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge());
		hEta[clusterSize]->Fill(this->transparentClusters[clusterSize].getEta());
		Float_t clusterLabPos = tracking->getPositionOfCluster(subjectDetector, this->transparentClusters[clusterSize], this->predPerpPosition, this->clusterCalcMode);
		hResidual[clusterSize]->Fill(clusterLabPos-this->predPosition);
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
