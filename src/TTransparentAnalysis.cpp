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
	sys->MakeDirectory(settings->getRelativePath().c_str());;
	sys->cd(settings->getRelativePath().c_str());
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
	htmlTransAna= new THTMLTransparentAnalysis(settings);
	htmlTransAna->setFileGeneratingPath(sys->pwd());
	sys->cd("..");
	
	
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
	
	
	initHistograms();
//	this->seedSigma=seedSigma;
//	this->hitSigma=hitSigma;
	cout<<"end initialise"<<endl;
	
	
	
	
}

TTransparentAnalysis::~TTransparentAnalysis() {
	// TODO Auto-generated destructor stub
	cout<<"\n\nClosing TTransparentAnalysis"<<endl;
	saveHistograms();
	deleteHistograms();
	
	// TODO: replace this!
	vector<vector <Float_t> > meanPulseHeights;
	vector<vector <Float_t> > mpPulseHeights;
	vector<vector <pair <Float_t,Float_t> > > resolutions;
//	meanPulseHeights.resize(2);
	resolutions.resize(2);
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
//		meanPulseHeights.at(0).push_back(100.*clusterSize);
//		meanPulseHeights.at(1).push_back(200.*clusterSize);
		pair <Float_t,Float_t> resolution;
		resolution.first = 0.;
		resolution.second = 50. * 0.1381;
		resolutions.at(0).push_back(resolution);
		resolutions.at(1).push_back(resolution);
	}
	
	mpPulseHeights.push_back(vecMPLandau);
	mpPulseHeights.push_back(vecMPLandau2Highest);
	meanPulseHeights.push_back(vecMeanLandau);
	meanPulseHeights.push_back(vecMeanLandau2Highest);
	
	htmlTransAna->createPulseHeightPlots(meanPulseHeights);
	htmlTransAna->createResolutionPlots(resolutions);
	htmlTransAna->generateHTMLFile();
	if(eventReader!=0)delete eventReader;
	if(histSaver!=0)delete histSaver;
	if(htmlTransAna)delete htmlTransAna;
	if(tracking!=0) delete tracking;
	sys->cd("..");
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
	TCluster transparentCluster = TCluster(tracking->getEvent_number(), det, -99, -99, TPlaneProperties::getNChannels(det));
	int currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		transparentCluster.addChannel(currentChannel, tracking->getRawSignal(det,currentChannel), tracking->getRawSignalInSigma(det,currentChannel), tracking->getAdcValue(det,currentChannel), tracking->isSaturated(det,currentChannel), settings->isDet_channel_screened(det,currentChannel));
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
		stringstream histNameLaundau, histNameLaundau2Highest, histNameEta, histNameResidualChargeWeighted, histNameResidualHighest2Centroid;
		// TODO: histogram naming!!
		histNameLaundau << "hDiaTranspAnaPulseHightOf" << clusterSize+1 << "Strips";
		histNameLaundau2Highest << "hDiaTranspAnaPulseHightOf2HighestIn" << clusterSize+1 << "Strips";
		histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize+1 << "Strips";
		histNameResidualChargeWeighted << "hDiaTranspAnaResidualChargeWeightedIn" << clusterSize+1 << "StripsMinusPred";
		histNameResidualHighest2Centroid << "hDiaTranspAnaResidualHighest2CentroidIn" << clusterSize+1 << "StripsMinusPred";
		hLaundau.push_back(new TH1F(histNameLaundau.str().c_str(),histNameLaundau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hLaundau2Highest.push_back(new TH1F(histNameLaundau2Highest.str().c_str(),histNameLaundau2Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hEta.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
		hResidualChargeWeighted.push_back(new TH1F(histNameResidualChargeWeighted.str().c_str(),histNameResidualChargeWeighted.str().c_str(),bins,-5.,5.));
		hResidualHighest2Centroid.push_back(new TH1F(histNameResidualHighest2Centroid.str().c_str(),histNameResidualHighest2Centroid.str().c_str(),bins,-5.,5.));
	}
}

void TTransparentAnalysis::fillHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		hLaundau[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge());
		hLaundau2Highest[clusterSize]->Fill(this->transparentClusters[clusterSize].getCharge(2,false));
		hEta[clusterSize]->Fill(this->transparentClusters[clusterSize].getEta());
//		if (clusterSize == 1 /*&& this->transparentClusters[clusterSize].getCharge() != this->transparentClusters[clusterSize].getCharge(2,false)*/) printEvent(this->transparentClusters[clusterSize]);
		if (clusterSize+1 != transparentClusters[clusterSize].getClusterSize()) {
			cout << "wrong cluster size!" << endl;
			cout << "clusterSize+1 = " << clusterSize+1 << "\ttransparentClusters[clusterSize].getClusterSize() = " << transparentClusters[clusterSize].getClusterSize() << endl;
		}
		hResidualChargeWeighted[clusterSize]->Fill(this->getResidual(this->transparentClusters[clusterSize],TCluster::chargeWeighted));
		hResidualHighest2Centroid[clusterSize]->Fill(this->getResidual(this->transparentClusters[clusterSize],TCluster::highest2Centroid));
	}
}

TF1* TTransparentAnalysis::doGaussFit(TH1F *histo) {
//	TH1* histo = (TH1*)htemp->Clone();
	if (histo->GetEntries()==0) return 0;
	TF1* histofitx = new TF1("histofitx","gaus",histo->GetMean()-2*histo->GetRMS(),histo->GetMean()+2*histo->GetRMS());
	histofitx->SetLineColor(kBlue);
	histo->Fit(histofitx,"rq");
//	histo->Draw();
	
//	TCanvas plots_canvas("plots_canvas","plots_canvas");
//	plots_canvas.cd();
//	histo->Draw();
//	pt->Draw();
//	ostringstream plot_filename;
//	plot_filename << plots_path << histo->GetName() << ".png";
//	plots_canvas.Print(plot_filename.str().c_str());
//	if (histo!=0) delete histo;
	return histofitx;
}

void TTransparentAnalysis::saveHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		TF1 *fitLandau=0, *fitLandau2Highest=0, *fitResidualChargeWeighted=0, *fitResidualHighest2Centroid=0;
		fitLandau = landauGauss->doLandauGaussFit(hLaundau[clusterSize]);
		fitLandau2Highest = landauGauss->doLandauGaussFit(hLaundau2Highest[clusterSize]);
		fitResidualChargeWeighted = doGaussFit(hResidualChargeWeighted[clusterSize]);
		fitResidualHighest2Centroid = doGaussFit(hResidualHighest2Centroid[clusterSize]);
		vecMPLandau.push_back(fitLandau->GetParameter(1));
		vecMPLandau2Highest.push_back(fitLandau2Highest->GetParameter(1));
		vecMeanLandau.push_back(hLaundau[clusterSize]->GetMean());
		vecMeanLandau2Highest.push_back(hLaundau2Highest[clusterSize]->GetMean());
		histSaver->SaveHistogramWithFit(hLaundau[clusterSize],fitLandau);
		histSaver->SaveHistogramWithFit(hLaundau2Highest[clusterSize],fitLandau2Highest);
		histSaver->SaveHistogram(hEta[clusterSize],0);
		histSaver->SaveHistogramWithFit(hResidualChargeWeighted[clusterSize],fitResidualChargeWeighted);
		histSaver->SaveHistogramWithFit(hResidualHighest2Centroid[clusterSize],fitResidualHighest2Centroid);
	}
}

// TODO: call TTransparentAnalysis::deleteHistograms
void TTransparentAnalysis::deleteHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		delete hLaundau[clusterSize];
		delete hLaundau2Highest[clusterSize];
		delete hEta[clusterSize];
		delete hResidualChargeWeighted[clusterSize];
		delete hResidualHighest2Centroid[clusterSize];
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
		cout << "\tresidual:\t" << this->getResidual(this->transparentClusters[clusterSize],this->clusterCalcMode) << endl;
		cout << "\tcluster pos in det system:\t" << this->transparentClusters[clusterSize].getPosition(this->clusterCalcMode) << endl;
		cout << "\tcluster pos in lab system:\t" << tracking->getPositionOfCluster(subjectDetector, this->transparentClusters[clusterSize], this->predPerpPosition, this->clusterCalcMode) << endl;
	}
	return;
}

void TTransparentAnalysis::printEvent(TCluster cluster) {
	cout << "\n--- event " << nEvent;
	cout << "\n\tcluster size: " << cluster.getClusterSize();
	cout << "\n\tcharge: " << cluster.getCharge(false);
	cout << "\n\tcharge of 2 highest centroid: " << cluster.getCharge(2,false);
	cout << "\n\thighest channel: " << cluster.getHighestSignalChannel();
	cout << "\n\tcluster position of highest channel: " << cluster.getClusterPosition(cluster.getHighestSignalChannel());
	cout << "\n\thighest channel is seed? " << cluster.isSeed(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
	cout << "\n\thighest channel is hit? " << cluster.isHit(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
	cout << "\n\tseed sigma: " << cluster.getSeedSigma();
	cout << "\n\thit sigma: " << cluster.getHitSigma();
	cout << "\n\tpredicted channel: " << positionInDetSystem;
	cout << "\n\tpredicted position: " << predPosition;
	cout << "\n\tcharge weighted position: " << tracking->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid position: " << tracking->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::highest2Centroid);
	cout << "\n\tcharge weighted residual: " << getResidual(cluster,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid residual: " << getResidual(cluster,TCluster::highest2Centroid);
	cout << endl;
	cluster.Print();
}

Float_t TTransparentAnalysis::getResidual(TCluster cluster, TCluster::calculationMode_t clusterCalculationMode) {
	return tracking->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,clusterCalculationMode)-this->predPosition;
}

