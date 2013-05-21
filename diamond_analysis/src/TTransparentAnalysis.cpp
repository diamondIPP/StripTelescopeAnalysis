//
//  TTransparentAnalysis.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 05.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TTransparentAnalysis.hh"

TTransparentAnalysis::TTransparentAnalysis(TSettings* settings, TSettings::alignmentMode mode) {
	cout<<"**********************************************************"<<endl;
	cout<<"********TTransparentAnalysis::TTransparentAnalysis********"<<endl;
	cout<<"**********************************************************"<<endl;
	// TODO Auto-generated constructor stub
	if(settings==0){
		cerr<<"Settings invalid:"<<settings<<endl;
		exit(-1);
	}
	sys = gSystem;
	setSettings(settings);
	
	settings->goToAlignmentRootDir();
	eventReader = new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(mode),settings->getEtaDistributionPath(),settings);
	// TODO: load settings!!!
	
	histSaver = new HistogrammSaver();
//	settings->goToTransparentAnalysisDir();
	histSaver->SetPlotsPath(settings->getTransparentAnalysisDir(mode));
	histSaver->SetRunNumber(settings->getRunNumber());
	htmlTransAna = new THTMLTransparentAnalysis(settings);
	htmlTransAna->setFileGeneratingPath(settings->getTransparentAnalysisDir(mode));
	
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
	verbosity = settings->getVerbosity();
	
	initHistograms();
	cout<<"end initialise"<<endl;
	positionPrediction = 0;
	inf  = std::numeric_limits<float>::infinity();

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
	resolutions.push_back(vecResidualEtaCorrected_2ndGaus);
	
	htmlTransAna->createPulseHeightPlots(meanPulseHeights, mpPulseHeights);
	htmlTransAna->createResolutionPlots(resolutions);
	htmlTransAna->createEtaPlots();
	htmlTransAna->createEtaIntegrals();
	htmlTransAna->generateHTMLFile();
	if (eventReader!=0) delete eventReader;
	if (histSaver!=0) delete histSaver;
	if (htmlTransAna!=0) delete htmlTransAna;
	settings->goToOutputDir();
}

void TTransparentAnalysis::analyze(UInt_t nEvents, UInt_t startEvent) {
	cout<<"\n\n******************************************\n";
	cout<<    "******Start Transparent Analysis...*******\n";
	cout<<"******************************************\n\n"<<endl;
	if(verbosity>10&&verbosity%2==1){
		cout<< "Press a Key and enter to continue..."<<flush;
		char t;
		cin >>t;
	}
//	usedForSiliconAlignment = 0;
	if(verbosity>6)cout<<"Current Dir: "<<sys->pwd()<<endl;
	if (nEvents+startEvent > eventReader->GetEntries()) {
		cout << "only "<<eventReader->GetEntries()<<" in tree!\n";
		nEvents = eventReader->GetEntries()-startEvent;
	}
	this->nEvents = nEvents;
	createEventVector(startEvent);
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
	for(UInt_t clusterSize = 0; clusterSize <TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);clusterSize++){
		vecvecResXChargeWeighted[clusterSize].clear();
		vecvecRelPos[clusterSize].clear();
		vecvecRelPos2[clusterSize].clear();
		vecvecEta[clusterSize].clear();
		vecvecEtaCMNcorrected[clusterSize].clear();
		vecvecResXHighest2Centroid[clusterSize].clear();
		vecvecResXEtaCorrected[clusterSize].clear();
		vecvecResXHighestHit[clusterSize].clear();
	}
//	vecPredictedPosition.clear();
//	vecRelPredictedPosition.clear();
//	vecPredictedChannel.clear();
//	vecChi2.clear();
	UInt_t maxSize = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);
	for (UInt_t iEvent = 0; iEvent < eventNumbers.size(); iEvent++) {
		TRawEventSaver::showStatusBar(iEvent,eventNumbers.size(),100);
		nEvent = eventNumbers.at(iEvent);
		eventReader->LoadEvent(nEvent);
		if(!this->predictPositions(false))
			continue;
		Float_t etaClusSizeOf2 = -1;
		for (UInt_t clusterSize = 0; clusterSize < maxSize; clusterSize++) {
			if (clusterSize == 2 && false) {
				cout << "using " << hEtaIntegrals[clusterSize]->GetName() << " to fill " << hResidualEtaCorrected[clusterSize]->GetName() << endl;
				printCluster(vecTransparentClusters.at(iEvent).at(clusterSize));
			}
			hResidualEtaCorrected[clusterSize]->Fill(getResidual(vecTransparentClusters.at(iEvent).at(clusterSize),TCluster::corEta,hEtaIntegrals[clusterSize]));
//			if (clusterSize == 1) printCluster(vecTransparentClusters.at(iEvent).at(clusterSize));
			Float_t metricPosInDetSystem = eventReader->getPositionInDetSystem(subjectDetector,predXPosition,predYPosition);
			Float_t channelPosInDetSystem = settings->convertMetricToChannelSpace(subjectDetector,metricPosInDetSystem);
			Float_t resXChargeWeighted = this->getResidual(this->vecTransparentClusters.at(iEvent)[clusterSize],TCluster::chargeWeighted,hEtaIntegrals[clusterSize]);
			Float_t resXEtaCorrected = this->getResidual(this->vecTransparentClusters.at(iEvent)[clusterSize],TCluster::corEta,hEtaIntegrals[clusterSize]);
			Float_t resXHighest2Centroid = this->getResidual(this->vecTransparentClusters.at(iEvent)[clusterSize],TCluster::highest2Centroid,hEtaIntegrals[clusterSize]);
			Float_t relChannelPos = channelPosInDetSystem - (int)(channelPosInDetSystem+.5);
			Float_t resXHighestHit= this->getResidual(this->vecTransparentClusters.at(iEvent)[clusterSize],TCluster::maxValue,hEtaIntegrals[clusterSize]);

			//			Float_t relHitPos = this->predPosition- (int)(predPosition+.5);
			Int_t leftChannel=-1;
			Float_t eta = vecTransparentClusters[iEvent][clusterSize].getEta(leftChannel);
			Float_t etaCMNcorrected = vecTransparentClusters[iEvent][clusterSize].getEta(true);
			if (clusterSize == 2){
				etaClusSizeOf2 = eta;
			}
			if(verbosity>4)
				cout<<nEvent<<": "<<clusterSize<<"clusterSize: "<<channelPosInDetSystem<<"-->"<<relChannelPos<<" <-> "<<resXChargeWeighted<<", "<<resXEtaCorrected<<", "<<resXHighest2Centroid<<endl;
			vecvecRelPos[clusterSize].push_back(relChannelPos);
			vecvecRelPos2[clusterSize].push_back(relChannelPos+.5);
			vecvecEta[clusterSize].push_back(eta);
			vecvecEtaCMNcorrected[clusterSize].push_back(etaCMNcorrected);
			//			if (resXChargeWeighted > -6000)
			vecvecResXChargeWeighted[clusterSize].push_back(resXChargeWeighted);
			//			if (resXHighest2Centroid > -6000)
			vecvecResXHighest2Centroid[clusterSize].push_back(resXHighest2Centroid);
			vecvecResXHighestHit[clusterSize].push_back(resXHighestHit);
			//			if (resXEtaCorrected > -6000)
			vecvecResXEtaCorrected[clusterSize].push_back(resXEtaCorrected);
			if(clusterSize==maxSize-1){
				Float_t deltaEta = eta-etaClusSizeOf2;
				vecDeltaEta.push_back(deltaEta);
				vecRelatedEta2.push_back(etaClusSizeOf2);
				vecRelatedEta10.push_back(eta);
				vecRelatedResXEtaCorrected.push_back(resXEtaCorrected);
				Float_t signalLeftOfEta = vecTransparentClusters[iEvent][clusterSize].getSignalOfChannel(leftChannel-1);
				Float_t signalRightOfEta = vecTransparentClusters[iEvent][clusterSize].getSignalOfChannel(leftChannel+2);
				Int_t highestClusterPos =  vecTransparentClusters[iEvent][clusterSize].getHighestHitClusterPosition();
				Float_t leftOfHighestSignal =  vecTransparentClusters[iEvent][clusterSize].getSignal(highestClusterPos-1);
				Float_t rightOfHighestSignal =  vecTransparentClusters[iEvent][clusterSize].getSignal(highestClusterPos+1);
				Float_t charge = vecTransparentClusters[iEvent][clusterSize].getCharge(2,false);
				Float_t highestSignal = vecTransparentClusters[iEvent][clusterSize].getHighestSignal();
				this->vecSignalLeftOfEta.push_back(signalLeftOfEta);
				this->vecSignalRightOfEta.push_back(signalRightOfEta);
				this->vecSignalLeftOfHighest.push_back(leftOfHighestSignal);
				this->vecSignalRightOfHighest.push_back(rightOfHighestSignal);
				this->vecClusterCharge.push_back(charge);
				this->vecHighestSignal.push_back(highestSignal);
				this->vecEta.push_back(eta);

			}
		}
	}
}

bool TTransparentAnalysis::predictPositions(bool savePrediction) {
	if (positionPrediction) delete positionPrediction;
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
	this->positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, this->predXPosition, this->predYPosition);
	this->positionInDetSystemChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,positionInDetSystemMetric);
	if(verbosity>5)cout<<"\nEventNo: "<<nEvent<<":\t"<<predPosition<<"/"<<predPerpPosition<<"--->"<<positionInDetSystemMetric<<" mum --> "<<positionInDetSystemChannelSpace<<" ch"<<flush;
	//		if (verbosity > 4) cout << "position in det system:\t" << this->positionInDetSystem << endl;
	//		if (verbosity > 4)
	//			cout << "clustered analysis strip position:\t" << eventReader->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
	Float_t chi2 = positionPrediction->getChi2();
	if(savePrediction){
	vecPredictedPosition.push_back(positionInDetSystemChannelSpace);
	vecRelPredictedPosition.push_back(positionInDetSystemChannelSpace-(int)(positionInDetSystemChannelSpace));
	vecChi2.push_back(chi2);}

	if(chi2>settings->getTransparentChi2())
		return false;
	return true;
}

bool TTransparentAnalysis::checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// get channel and direction for clustering
	UInt_t centerChannel;
	if(verbosity>3)
		cout<<"\ncheck Pred Region: "<<nEvent<< " "<<det<<" "<<centerPosition<<" "<<clusterSize<<endl;
	int direction;
	direction = getSignedChannelNumber(centerPosition);
	centerChannel = TMath::Abs(direction);
	if (direction < 0) direction = -1;
	else direction = 1;
//	cout<<"\t"<<direction<<" x "<< centerChannel<<endl;
	
	// check predicted cluster channels
	UInt_t currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		direction *= -1;
		currentChannel += direction * iChannel;
		if (currentChannel < 0) {
			if (verbosity > 5) cout << "\tchannel " << currentChannel << " is not on this detector.." << endl;
			regionNotOnPlane++;
			return false;
		}
		if (currentChannel > TPlaneProperties::getNChannels(det)-1) {
			if (verbosity > 5) cout << "\tchannel " << currentChannel << " is not on this detector.." << endl;
			regionNotOnPlane++;
			return false;
		}
		if (this->settings->getDet_channel_screen(det).isScreened(currentChannel) == true) {
			if (verbosity > 5) cout << "\tchannel " << currentChannel << " is screened.." << endl;
			screenedChannel++;
			return false;
		}
		if (eventReader->isSaturated(det, currentChannel) == true) {
			if (verbosity > 5) cout << "\tchannel " << currentChannel << " has saturated.." << endl;
			saturatedChannel++;
			return false;
		}
	}
	return true;
}

// TODO: avoid wrong channel numbers (>128, <0)
/**
 *	returns the next channel number including a sign: + if pos - (int)pos <.5 else minus
 *	different approach:
 *	 ( 1+2*( (int)pos-(int)(pos+.5) ) ) * int (pos+.5)
 * @param position
 * @author Lukas Baeni
 * @return
 */
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
	UInt_t bins=512;
	vecvecResXChargeWeighted.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecResXHighest2Centroid.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecResXHighestHit.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecResXEtaCorrected.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecRelPos.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecRelPos2.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecEta.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecvecEtaCMNcorrected.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecVecLandau.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
	vecVecFidCutX.clear();
	vecVecFidCutY.clear();
	hEtaIntegrals.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));

	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		vecVecLandau.at(clusterSize).clear();
		vecvecEtaCMNcorrected.at(clusterSize).clear();
		vecvecRelPos.at(clusterSize).clear();
		vecvecRelPos2.at(clusterSize).clear();
		vecvecResXEtaCorrected.at(clusterSize).clear();
		vecvecResXChargeWeighted.at(clusterSize).clear();
		// TODO: take care of histogram names and bins!!
		stringstream histNameLandau,histNameLandau1Highest, histNameLandau2Highest, histNameEta, histNameResidualChargeWeighted, histNameResidualHighest2Centroid, histNameResidualEtaCorrected;
		// TODO: histogram naming!!
		histNameLandau << "hDiaTranspAnaPulseHeightOf" << clusterSize+1 << "Strips";
		histNameLandau1Highest<< "hDiaTranspAnaPulseHeightOfHighestIn"<<clusterSize+1<<"Strips";
		histNameLandau2Highest << "hDiaTranspAnaPulseHeightOf2HighestIn" << clusterSize+1 << "Strips";
		histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize+1 << "Strips";
		histNameResidualChargeWeighted << "hDiaTranspAnaResidualChargeWeightedIn" << clusterSize+1 << "StripsMinusPred";
		histNameResidualHighest2Centroid << "hDiaTranspAnaResidualHighest2CentroidIn" << clusterSize+1 << "StripsMinusPred";
		histNameResidualEtaCorrected << "hDiaTranspAnaResidualEtaCorrectedIn" << clusterSize+1 << "StripsMinusPred";
		TString nameResVsHitChargeWeighted = TString::Format("hDiaTransAnaResVsHitChargeWeightedIn%d",clusterSize+1);
		TString nameResVsHitEtaCor = TString::Format("hDiaTransAnaResVsHitEtaCorIn%d",clusterSize+1);
		TString nameResVsHitHeighest2Centroid = TString::Format("hDiaTransAnaResVsHitHigehst2CentroidIn%d",clusterSize+1);
		TString histNameResidualHighestHit = TString::Format("hDiaTranspAnaResidualHighestHitIn%dStripsMinusPred",clusterSize+1);

		hLandau.push_back(new TH1F(histNameLandau.str().c_str(),histNameLandau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		hLandau2Highest.push_back(new TH1F(histNameLandau2Highest.str().c_str(),histNameLandau2Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
		this->hLandau1Highest.push_back(new TH1F(histNameLandau1Highest.str().c_str(),histNameLandau1Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));

		hEta.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
		histNameEta.str("");
		histNameEta.clear();
		histNameEta << "hDiaTranspAnaEtaCMNcorrected2HighestIn" << clusterSize+1 << "Strips";
		hEtaCMNcorrected.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
		///TODO: PitchWidth Plot width
		hResidualChargeWeighted.push_back(new TH1F(histNameResidualChargeWeighted.str().c_str(),histNameResidualChargeWeighted.str().c_str(),bins,-2.5*50,2.5*50));
		hResidualHighest2Centroid.push_back(new TH1F(histNameResidualHighest2Centroid.str().c_str(),histNameResidualHighest2Centroid.str().c_str(),bins,-2.5*50,2.5*50));
		hResidualEtaCorrected.push_back(new TH1F(histNameResidualEtaCorrected.str().c_str(),histNameResidualEtaCorrected.str().c_str(),bins,-2.5*50,2.5*50));
		hResidualHighestHit.push_back(new TH1F(histNameResidualHighestHit,histNameResidualHighestHit,bins,-2.5*50,2.5*50));

		//		hResidualVsHitPositionChargeWeighted.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted));
//		hResidualVsHitPositionEtaCorrected.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted));
//		hResidualVsHitPositionHigehest2Centroid.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted,));

	}
	hLandauMean = new TH1F("hDiaTranspAnaPulseHeightMean","hDiaTranspAnaPulseHeightMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLandauMP = new TH1F("hDiaTranspAnaPulseHeightMP","hDiaTranspAnaPulseHeightMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLandau2HighestMean = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMean","hDiaTranspAnaPulseHeightOf2HighestMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hLandau2HighestMP = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMP","hDiaTranspAnaPulseHeightOf2HighestMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
	hPredictedPositionInStrip = new TH1F("hPredictedPositionInStrip","hPredictedPositionInStrip",2,-1.5,1.5);
}

void TTransparentAnalysis::fillHistograms() {
		vecVecFidCutX.push_back(eventReader->getFiducialValueX());
		vecVecFidCutY.push_back(eventReader->getFiducialValueY());
		vecPredictedChannel.push_back(positionInDetSystemChannelSpace);
		UInt_t maxSize = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);
	for (UInt_t clusterSize = 0; clusterSize < maxSize; clusterSize++) {
		Float_t charge = this->transparentClusters[clusterSize].getCharge();
		Float_t chargeOfTwo = this->transparentClusters[clusterSize].getCharge(2,false);
		vecVecLandau[clusterSize].push_back(charge);
		hLandau[clusterSize]->Fill(charge);
		hLandau2Highest[clusterSize]->Fill(chargeOfTwo);

		Float_t eta = this->transparentClusters[clusterSize].getEta();
		Float_t etaCMN = this->transparentClusters[clusterSize].getEta(true);
		if(eta>0&&eta<1)
			hEta[clusterSize]->Fill(eta);
		else if(verbosity>3)
			this->transparentClusters[clusterSize].Print();
		if(etaCMN>0&&etaCMN<1)
			hEtaCMNcorrected[clusterSize]->Fill(etaCMN);

		if (clusterSize>2){
			Int_t highestSignalClusterPos = this->transparentClusters[clusterSize].getHighestHitClusterPosition();
			if(highestSignalClusterPos<0){
				cout<<nEvent<<"errror5: highestSignalClusterPos: "<<highestSignalClusterPos<<endl;
				this->transparentClusters[clusterSize].Print();
			}
//			Float_t leftSig = transparentClusters[clusterSize].getSignal(highestSignalClusterPos-1);
//			Float_t rightSig= transparentClusters[clusterSize].getSignal(highestSignalClusterPos+1);

		}
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
		Float_t relPos =this->predPosition-(int)(this->predPosition+.5);
		Float_t residualCW =this->getResidual(this->transparentClusters[clusterSize],TCluster::chargeWeighted,hEtaIntegrals[clusterSize]);
		Float_t residualH2C = this->getResidual(this->transparentClusters[clusterSize],TCluster::highest2Centroid,hEtaIntegrals[clusterSize]);
		Float_t resXHighestHit= this->getResidual(this->transparentClusters[clusterSize],TCluster::maxValue,hEtaIntegrals[clusterSize]);

//		if(hResidualChargeWeightedVsEstimatedHitPosition==0)
//			hResidualChargeWeightedVsEstimatedHitPosition->Fill(residualCW,relPos,clusterSize);
//		if(hResidualHighest2CentroidVsEstimatedHitPosition)
//			hResidualHighest2CentroidVsEstimatedHitPosition>Fill(residualH2C,relPos,clusterSize);
		if (clusterSize+1 != transparentClusters[clusterSize].getClusterSize()) {
			cout << "wrong cluster size!" << endl;
			cout << "clusterSize+1 = " << clusterSize+1 << "\ttransparentClusters[clusterSize].getClusterSize() = " << transparentClusters[clusterSize].getClusterSize() << endl;
		}
		vecvecResXChargeWeighted[clusterSize].push_back(residualCW);
		vecvecResXHighest2Centroid[clusterSize].push_back(residualH2C);
		vecvecResXHighestHit[clusterSize].push_back(resXHighestHit);
		vecvecRelPos[clusterSize].push_back(relPos);
		vecvecRelPos2[clusterSize].push_back(relPos+.5);
		vecvecEta[clusterSize].push_back(eta);
		vecvecEtaCMNcorrected[clusterSize].push_back(etaCMN);

		hResidualChargeWeighted[clusterSize]->Fill(residualCW);
		hResidualHighest2Centroid[clusterSize]->Fill(residualH2C);
		hResidualHighestHit[clusterSize]->Fill(resXHighestHit);

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

TF1* TTransparentAnalysis::doDoubleGaussFit(TH1F *histo){
	if (!histo)return 0;
	if (histo->GetEntries()==0) return 0;
	Float_t mean = histo->GetMean();
	Float_t sigma = histo->GetRMS();
	Float_t max = histo->GetBinContent(histo->GetMaximumBin());
	Float_t pw = settings->getPitchWidth(subjectDetector);
	int nSigmas = 4;
	Float_t xmin = TMath::Min(-pw,mean-nSigmas*sigma);
	Float_t xmax = TMath::Max(pw,mean+nSigmas*sigma);
	TF1* histofitx = new TF1("fDoubleGaus","gaus(0)+gaus(3)",xmin,xmax);
	histofitx->SetParLimits(0,max/10,max);
	histofitx->SetParLimits(1,mean-2*sigma,mean+2*sigma);
	histofitx->SetParLimits(2,sigma/10,2*sigma);
	histofitx->SetParLimits(3,max/20,max/2);
	histofitx->SetParLimits(4,mean-2*sigma,mean+2*sigma);
	histofitx->SetParLimits(5,sigma/10,4*sigma);
	histofitx->SetParameters(.75*max,mean,sigma/5,.1*max,mean,sigma/4);
	histofitx->SetParNames("C_{0}","#mu_{0}","#sigma_{0}","C_{1}","#mu_{1}","#sigma_{1}");
	histofitx->SetLineColor(kBlue);
	histofitx->SetNpx(1000);
	histo->Fit(histofitx,"rq");
	return histofitx;
}

void TTransparentAnalysis::createEtaIntegrals() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		stringstream histName;
		histName << "hDiaTranspAnaEtaIntegral2HighestIn"<<clusterSize+1<<"Strips";
		if (hEtaIntegrals.at(clusterSize))
			delete hEtaIntegrals.at(clusterSize);
		hEtaIntegrals.at(clusterSize) = (TClustering::createEtaIntegral(hEta[clusterSize], histName.str()));
	}
}


void TTransparentAnalysis::createEfficiencyPlots(TH1F *hLandau){
	TString name = TString::Format("hEfficiency_%s",hLandau->GetName());
	TH1F* hEfficiency = new TH1F(name,name,settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector) );
	Float_t nentries = hLandau->GetEntries();
	Float_t integral = 0;
	for(Int_t bin = 1;bin <= hLandau->GetNbinsX();bin++){
		integral += hLandau->GetBinContent(bin);
		hEfficiency->SetBinContent(bin, (1.-integral/nentries)*100);
	}
	if(hEfficiency){
		hEfficiency->GetXaxis()->SetTitle("PH / adc counts");
		hEfficiency->GetYaxis()->SetTitle("efficientcy / %");
	}
//	TCutG *MP = new TCutG("gMP",1);
	histSaver->SaveHistogram(hEfficiency);
	if(hEfficiency) delete hEfficiency;
}


void TTransparentAnalysis::fitHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {

		createEfficiencyPlots(hLandau2Highest[clusterSize]);
		createEfficiencyPlots(hLandau[clusterSize]);
		createEfficiencyPlots(hLandau1Highest[clusterSize]);
		vector<Float_t> vecResChargeWeighted = vecvecResXChargeWeighted[clusterSize];

		stringstream name;
		name <<"hDiaTranspAnaResidualChargeWeightedIn" << clusterSize+1 << "StripsMinusPred";
		if(hResidualChargeWeighted[clusterSize])
			delete hResidualChargeWeighted[clusterSize];
		hResidualChargeWeighted[clusterSize] = histSaver->CreateDistributionHisto(name.str(), vecResChargeWeighted,8192,HistogrammSaver::maxWidth,-5000);
		Float_t plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
		hResidualChargeWeighted[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
		hResidualChargeWeighted[clusterSize]->GetXaxis()->SetTitle("Residual, ChargeWeighted / #mum");
		hResidualChargeWeighted[clusterSize]->GetYaxis()->SetTitle("number of entries #");

		name.str("");name.clear();
		name <<"hDiaTranspAnaResidualHighest2CentroidIn" << clusterSize+1 << "StripsMinusPred";
		if(hResidualHighest2Centroid[clusterSize] )
			delete hResidualHighest2Centroid[clusterSize] ;
		hResidualHighest2Centroid[clusterSize] = histSaver->CreateDistributionHisto(name.str(), vecvecResXHighest2Centroid[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
		plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
		hResidualHighest2Centroid[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
		hResidualHighest2Centroid[clusterSize]->GetXaxis()->SetTitle("Residual, highest 2 centroid / #mum");
		hResidualHighest2Centroid[clusterSize]->GetYaxis()->SetTitle("number of entries #");

		name.str("");name.clear();
		name <<"hDiaTranspAnaResidualEtaCorrectedIn" << clusterSize+1 << "StripsMinusPred";
		if(hResidualEtaCorrected[clusterSize])
			delete hResidualEtaCorrected[clusterSize];
		hResidualEtaCorrected[clusterSize] = histSaver->CreateDistributionHisto(name.str(), vecvecResXEtaCorrected[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
		plotWidth = settings->getPitchWidth(subjectDetector);
		hResidualEtaCorrected[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
		hResidualEtaCorrected[clusterSize]->GetXaxis()->SetTitle("Residual #eta corrected / #mum");
		hResidualEtaCorrected[clusterSize]->GetYaxis()->SetTitle("number of entries #");

		name.str("");name.clear();
		name <<"hDiaTranspAnaResidualHighestHitIn" << clusterSize+1 << "StripsMinusPred";
		if(hResidualHighestHit[clusterSize] )
			delete hResidualHighestHit[clusterSize] ;
		hResidualHighestHit[clusterSize] = histSaver->CreateDistributionHisto(name.str(), vecvecResXHighestHit[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
		plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
		hResidualHighestHit[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
		hResidualHighestHit[clusterSize]->GetXaxis()->SetTitle("Residual, highest Hit / #mum");
		hResidualHighestHit[clusterSize]->GetYaxis()->SetTitle("number of entries #");
		//,
		//,1024,HistogrammSaver::maxWidth,-6000);
		// fit histograms
		fitLandau.push_back(landauGauss->doLandauGaussFit(hLandau[clusterSize]));
		fitLandau2Highest.push_back(landauGauss->doLandauGaussFit(hLandau2Highest[clusterSize]));
		fitResidualChargeWeighted.push_back(doGaussFit(hResidualChargeWeighted[clusterSize]));
		fitResidualHighest2Centroid.push_back(doGaussFit(hResidualHighest2Centroid[clusterSize]));
		fitResidualEtaCorrected.push_back(doDoubleGaussFit(hResidualEtaCorrected[clusterSize]));
		
		// save fit parameters
		vecMPLandau.push_back(fitLandau[clusterSize]->GetParameter(1));
		vecMPLandau2Highest.push_back(fitLandau2Highest[clusterSize]->GetParameter(1));
		hLandauMP->SetBinContent(clusterSize+1,fitLandau[clusterSize]->GetParameter(1));
		hLandau2HighestMP->SetBinContent(clusterSize+1,fitLandau2Highest[clusterSize]->GetParameter(1));
		vecMeanLandau.push_back(hLandau[clusterSize]->GetMean());
		vecMeanLandau2Highest.push_back(hLandau2Highest[clusterSize]->GetMean());
		hLandauMean->SetBinContent(clusterSize+1,hLandau[clusterSize]->GetMean());
		hLandau2HighestMean->SetBinContent(clusterSize+1,hLandau2Highest[clusterSize]->GetMean());
		pair <Float_t,Float_t> tempPair,tempPair2;
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
			tempPair2.first = fitResidualEtaCorrected[clusterSize]->GetParameter(4);
			tempPair2.second = fitResidualEtaCorrected[clusterSize]->GetParameter(5);
		}
		else {
			tempPair.first = 0;
			tempPair.second = 0;
		}
		if (tempPair2.second>tempPair.second){
			vecResidualEtaCorrected.push_back(tempPair);
			vecResidualEtaCorrected_2ndGaus.push_back(tempPair2);
		}
		else{
			vecResidualEtaCorrected.push_back(tempPair2);
			vecResidualEtaCorrected_2ndGaus.push_back(tempPair);
		}
	}
	hLandauMean->Scale(1./hLandauMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLandauMP->Scale(1./hLandauMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLandau2HighestMean->Scale(1./hLandau2HighestMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
	hLandau2HighestMP->Scale(1./hLandau2HighestMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
}

void TTransparentAnalysis::analyseEtaDistribution(TH1F* hEtaDist){
	if(!hEtaDist)
		return;
	if(hEtaDist->GetEntries()==0)
		return;
	float threshold = 0.1;
	int n=0;
	int ntries = 0;
	int maxTries =30;
	while (n!=2&&ntries < maxTries){
		n = hEtaDist->ShowPeaks(3,"nobackground",threshold);
		if(n<2){
			threshold*=.9;
			//				cout<<ntries<<"-"<<n<<" ==> lowering threshold "<<threshold<<endl;
		}
		else if(n>2){
			threshold*=1.1;
			//				cout<<ntries<<"-"<<n<<" ==> higher threshold "<<threshold<<endl;
		}
		ntries++;
	}
	TList *functions = hEtaDist->GetListOfFunctions();
	TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	cout<<hEtaDist->GetName()<<" - "<<ntries<<endl;
	if (!pm){
		if(functions) functions->Print();
		return;
	}
	for (int i=0;i< pm->GetN();i++)
		cout<<"\t"<<i<<"\t"<<pm->GetX()[i]*100.<<": "<<pm->GetY()[i]<<"\n";
	if(pm->GetN()==2){
		Float_t x_0 = pm->GetX()[0];
		Float_t x_1 = pm->GetX()[1];
		Float_t y_0 = pm->GetY()[0];
		Float_t y_1 = pm->GetY()[1];
		if(x_0>x_1){
			x_1 = x_0;
			x_0 = pm->GetX()[1];
			y_0 = y_1;
			y_1 = pm->GetY()[0];
		}
		if(x_1>.5)
			x_1 = 1-x_1;
		x_0*=100.;
		x_1*=100.;
		cout<<"\t\t"<<x_0<<" - "<< x_1 <<"\t"<<(x_0-x_1)<<"\t"<<x_0/x_1<<"\t"<<(x_0-x_1)*100./x_0<<endl;
		cout<<"\t\t"<<y_0<<" - "<< y_1 <<"\t"<<(y_0-y_1)<<"\t"<<y_0/y_1<<"\t"<<(y_0-y_1)*100/y_0<<endl;
	}
	cout<<"\n"<<flush;
}

void TTransparentAnalysis::analyseEtaDistributions(){
	stringstream name;
	name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_etaOf2");
	if(vecRelatedEta2.size()!=vecDeltaEta.size()){
		cout<<"Something is wrong with vedDeltaEta Size"<<flush;
		char t;
		cin>>t;
	}
	TH2F* hDeltaEta = histSaver->CreateScatterHisto(name.str(),vecRelatedEta2,vecDeltaEta,400,400,-1,1,0,1);
	if(hDeltaEta){
		hDeltaEta->GetXaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
		hDeltaEta->GetYaxis()->SetTitle("#eta_{2 of 2}");
	}
	histSaver->SaveHistogram(hDeltaEta,false);

	if(hDeltaEta)delete hDeltaEta;

	name.str("");name.clear();
	name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_ResidualEtaCorrectedIn10");
	Float_t pw = settings->getDiamondPitchWidth();
	TH2F* hDeltaEtaVsResidual = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedResXEtaCorrected,512,400,-2*pw,2*pw,-1,1);
	if(hDeltaEtaVsResidual){
		hDeltaEtaVsResidual->GetXaxis()->SetTitle("Residual Eta Correct, 2 of 10 / #mum");
		hDeltaEtaVsResidual->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
	}
	histSaver->SaveHistogram(hDeltaEtaVsResidual,false);

	name.str("");name.clear();
	name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_etaOf10");
	TH2F* hDeltaEtaVsEta = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedEta10,512,512,0,1,0,1);
	if(hDeltaEtaVsEta){
		hDeltaEtaVsEta->GetXaxis()->SetTitle("#eta_{2 of 10}");
		hDeltaEtaVsEta->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
	}
	histSaver->SaveHistogram(hDeltaEtaVsEta,false);
	Float_t maxDelta = .199999;
	Int_t bin1 = hDeltaEtaVsEta->GetYaxis()->FindBin(-maxDelta);
	Int_t bin2 = hDeltaEtaVsEta->GetYaxis()->FindBin(+maxDelta);
	TString hName = TString::Format("hEtaIn10_DeltaEta_below_020");
	TH1F* hEtaBoundedEtaCorrected = (TH1F*)hDeltaEtaVsEta->ProjectionX(hName,bin1,bin2);
	hEtaBoundedEtaCorrected->SetTitle(TString::Format("#eta_{2 of 10}, |#Delta#eta| < 0.2"));
	histSaver->SaveHistogram(hEtaBoundedEtaCorrected);
	if(hEtaBoundedEtaCorrected)delete hEtaBoundedEtaCorrected;

//	name.str("");name.clear();
//	name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_ResidualEtaCorrectedIn10");
//	TH2F* hDeltaEtaVsResidual = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedResXEtaCorrected,512,400,-2*pw,2*pw,-1,1);
//	if(hDeltaEtaVsResidual){
//		hDeltaEtaVsResidual->GetXaxis()->SetTitle("Residual Eta Correct, 2 of 10 / #mum");
//		hDeltaEtaVsResidual->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
//	}
//	histSaver->SaveHistogram(hDeltaEtaVsResidual,false);

	maxDelta = .199999;
	bin1 = hDeltaEtaVsResidual->GetYaxis()->FindBin(-maxDelta);
	bin2 = hDeltaEtaVsResidual->GetYaxis()->FindBin(+maxDelta);
	hName = TString::Format("hResidualEtaCorrectedIn10_DeltaEta_below_020");
	TH1F* hEtaBoundedEtaCorrectedResidual = (TH1F*)hDeltaEtaVsResidual->ProjectionX(hName,bin1,bin2);
	hEtaBoundedEtaCorrectedResidual->SetTitle(TString::Format("Residual #eta corrected in 10 strips, |#Delta#eta| < 0.2"));
	histSaver->SaveHistogram(hEtaBoundedEtaCorrectedResidual);
	if(hDeltaEtaVsResidual)delete hDeltaEtaVsResidual;

	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		TH1F* hEtaDist = hEta[clusterSize];
		analyseEtaDistribution(hEtaDist);
		hEtaDist = hEtaCMNcorrected[clusterSize];
		analyseEtaDistribution(hEtaDist);
	}

	name.str("");name.clear();
	name<<"hEtaVsSignalLeftOfEta";
	TH2F* histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalLeftOfEta,this->vecEta);
	if(histo2d){
		histo2d->GetYaxis()->SetTitle("Signal left of #eta");
		histo2d->GetXaxis()->SetTitle("#eta");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}

	name.str("");name.clear();
	name<<"hEtaVsSignalRightOfEta";
	histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalRightOfEta,this->vecEta);
	if(histo2d){
		histo2d->GetYaxis()->SetTitle("Signal right of #eta");
		histo2d->GetXaxis()->SetTitle("#eta");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	vector<Float_t> vecRightFactor, vecLeftFactor,vecRightOverHighest,vecLeftOverHighest;
	for(UInt_t i=0;i<vecSignalLeftOfHighest.size()&&i<vecSignalRightOfHighest.size()&&i<vecClusterCharge.size();i++){
		vecRightFactor.push_back(vecSignalRightOfHighest.at(i)/vecClusterCharge.at(i));
		vecLeftFactor.push_back(vecSignalLeftOfHighest.at(i)/vecClusterCharge.at(i));
		vecLeftOverHighest.push_back(vecSignalLeftOfHighest.at(i)/vecHighestSignal.at(i));
		vecRightOverHighest.push_back(vecSignalRightOfHighest.at(i)/vecHighestSignal.at(i));
	}
	name.str("");name.clear();
	name<<"hSignalLeftOverHighestVsSignalRightOverHighest";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecRightOverHighest,vecLeftOverHighest);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("signal left of highest signal over highest signal");
		histo2d->GetYaxis()->SetTitle("signal right of highest signal over highest signal");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}

	name.str("");name.clear();
	name<<"hSignalRightVsHighestSignal";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalRightOfHighest,vecHighestSignal);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("highest signal");
		histo2d->GetYaxis()->SetTitle("signal right of highest signal ");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	name.str("");name.clear();
	name<<"hSignalLeftVsHighestSignal";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalLeftOfHighest,vecHighestSignal);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("highest signal");
		histo2d->GetYaxis()->SetTitle("signal left of highest signal ");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	name.str("");name.clear();
	name<<"hSignalLeftOfHighest";
	TH1F* histoLeft = histSaver->CreateDistributionHisto(name.str(),vecLeftFactor);
	histSaver->SaveHistogram(histoLeft);

	name.str("");name.clear();
	name<<"hSignalRightOfHighest";
	TH1F* histoRight = histSaver->CreateDistributionHisto(name.str(),vecRightFactor);
	histSaver->SaveHistogram(histoRight);

	name.str("");name.clear();
	name<<"cSignalNextToHighest";
	histoLeft->SetLineColor(kBlue);
	histoRight->SetLineColor(kRed);
	Float_t max = TMath::Max(histoLeft->GetMaximum(),histoRight->GetMaximum());
	histoLeft->SetMaximum(max);
	histoRight->SetMaximum(max);
	histSaver->SaveTwoHistos(name.str(),histoLeft,histoRight,1.,false);
	if(histoLeft) delete histoLeft;
	if(histoRight) delete histoRight;

	name.str("");name.clear();
	name<<"hSignalLeftOfEtaChannels";
	histoLeft = histSaver->CreateDistributionHisto(name.str(),this->vecSignalLeftOfEta);
	if(histoLeft){
		histoLeft->GetXaxis()->SetTitle("Signal left of #eta");
		histoLeft->GetYaxis()->SetTitle("number of entries #");
		histoLeft->SetLineColor(kBlue);
		histSaver->SaveHistogram(histoLeft);
	}
	name.str("");name.clear();
	name<<"hSignalRightOfEtaChannels";
	histoRight = histSaver->CreateDistributionHisto(name.str(),this->vecSignalRightOfEta);
	if(histoLeft){
		histoRight->GetXaxis()->SetTitle("Signal right of #eta");
		histoRight->GetYaxis()->SetTitle("number of entries #");
		histoRight->SetLineColor(kRed);
		histSaver->SaveHistogram(histoLeft);
	}
	name.str("");name.clear();
	name<<"cSignalOfSignalsAdjacentToEta";
	histSaver->SaveTwoHistos(name.str(),histoLeft,histoRight,1.,false);
	if(histoLeft) delete histoLeft;
	if(histoRight) delete histoRight;

}
void TTransparentAnalysis::saveHistograms() {

	stringstream name;
	name << "hLandauVsFidCutX_ClusterSize10";

	int i = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)-1;
	TH2F* htemp = histSaver->CreateScatterHisto(name.str(),vecVecFidCutX,vecVecLandau.at(i));
	if(htemp){
		htemp->GetXaxis()->SetTitle("pulse height, clusterSize 10");
		htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position /ch");
	}
	//else{
	else cout<<"couldn't ";
		cout<<"create scatterplot: "<<name<<"\t"<<vecVecFidCutX.size()<<" "<<vecVecLandau.at(i).size()<<" "<<i<<" "<<endl;
	if(!htemp){
		char t;
		cin>>t;
	}
	histSaver->SaveHistogram(htemp);
	if(htemp)delete htemp;

	name.str("");name.clear();
	name << "hLandauVsPredChannel_ClusterSize10";
	htemp = histSaver->CreateScatterHisto(name.str(),vecPredictedChannel,vecVecLandau.at(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)-1));
	if(htemp){
		htemp->GetXaxis()->SetTitle("pulse height, clusterSize 10");
		htemp->GetYaxis()->SetTitle("predicted channel position /ch");
	}
	histSaver->SaveHistogram(htemp);
	if(htemp) delete htemp;

	name.str("");name.clear();
	name << "hLandauVsFidCutY_ClusterSize10";
	htemp = histSaver->CreateScatterHisto(name.str(),vecVecFidCutY,vecVecLandau.at(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)-1));
	if(htemp){
		htemp->GetXaxis()->SetTitle("pulse height, clusterSize 10");
		htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position /ch");
	}
	histSaver->SaveHistogram(htemp);
	if (htemp) delete htemp;
	analyseEtaDistributions();
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		histSaver->SaveHistogram(hEta[clusterSize],0);
		histSaver->SaveHistogram(hEtaCMNcorrected[clusterSize],0);
//		if (clusterSize == 0) {
			histSaver->SaveHistogramLandau(hLandau[clusterSize]);
			histSaver->SaveHistogramLandau(hLandau2Highest[clusterSize]);
			histSaver->SaveHistogramLandau(hLandau1Highest[clusterSize]);
			histSaver->SaveHistogram(hResidualChargeWeighted[clusterSize]);
			histSaver->SaveHistogram(hResidualHighest2Centroid[clusterSize]);
			histSaver->SaveHistogram(hResidualHighestHit[clusterSize]);
//		}
//		else {
//			histSaver->SaveHistogramLandau(hLaundau[clusterSize],fitLandau[clusterSize]);
//			histSaver->SaveHistogramWithFit(hResidualChargeWeighted[clusterSize],fitResidualChargeWeighted[clusterSize]);
//			histSaver->SaveHistogramWithFit(hResidualHighest2Centroid[clusterSize],fitResidualHighest2Centroid[clusterSize]);
//		}
		histSaver->SaveHistogramWithFit(hResidualEtaCorrected[clusterSize],fitResidualEtaCorrected[clusterSize]);
		histSaver->SaveHistogram(hEtaIntegrals[clusterSize],0);
	}
	histSaver->SaveHistogram(hLandauMean);
	histSaver->SaveHistogram(hLandauMP);
	histSaver->SaveHistogram(hLandau2HighestMean);
	histSaver->SaveHistogram(hLandau2HighestMP);
	Float_t pw = settings->getPitchWidth(subjectDetector);
	for(UInt_t i=0;i<TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);i++){
		string name = (string)TString::Format("hRelPosVsResolutionEtaCorrectedIn%d",i+1);
		if(verbosity>6)cout<<"creating "<<name<<": "<<vecvecRelPos[i].size()<<"-"<<vecvecResXEtaCorrected[i].size()<<endl;
		TH2F* hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXEtaCorrected[i],512,512,-6000);
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		hist->GetYaxis()->SetTitle("Relative predicted Position ");
		hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		histSaver->SaveHistogram(hist);
		if (hist)delete hist;

		name = (string)TString::Format("hRelChPos2VsResChargeWeighted_In_%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecResXEtaCorrected[i],512,512,-6000);
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		hist->GetYaxis()->SetTitle("Relative predicted Position");
		hist->GetXaxis()->SetTitle("Residual, charge Weighted / #mum");
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		histSaver->SaveHistogram(hist);
		if (hist)delete hist;

		name = (string)TString::Format("hRelChPosVsResChargeWeighted_In_%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXChargeWeighted[i],512,512,-6000);
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		hist->GetYaxis()->SetTitle("Relative predicted Position");
		hist->GetXaxis()->SetTitle("Residual, charge Weighted / #mum");
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		histSaver->SaveHistogram(hist);
		if (hist)delete hist;

		name = (string)TString::Format("hRelChPosVsResHighest2Centroid_In_%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXHighest2Centroid[i],512,512,-6000);
		hist->GetYaxis()->SetTitle("Relative predicted Position");
		hist->GetXaxis()->SetTitle("Residual, Highest 2 Centorid / #mum");
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		histSaver->SaveHistogram(hist);
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		if (hist)delete hist;

		name = (string)TString::Format("hEtaVsResolutionEtaCorrectedIn%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecEta[i],vecvecResXEtaCorrected[i],512,512,-6000);
		hist->GetYaxis()->SetTitle("#eta");
		hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		histSaver->SaveHistogram(hist);
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		if (hist)delete hist;

		name = (string)TString::Format("hEtaCMNCorrectedVsResolutionEtaCorrectedIn%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecEtaCMNcorrected[i],vecvecResXEtaCorrected[i],512,512,-6000);
		hist->GetYaxis()->SetTitle("#eta_{CMN-corrected}");
		hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
		hist->GetXaxis()->SetRangeUser(-pw,pw);
		histSaver->SaveHistogram(hist);
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		if (hist)delete hist;

		name = (string)TString::Format("hRelChPosVsEta_In_%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecEta[i],512);
		hist->GetXaxis()->SetRangeUser(0,1);
		hist->GetYaxis()->SetTitle("Relative predicted Position");
		hist->GetXaxis()->SetTitle("#eta ");
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		histSaver->SaveHistogram(hist);
		if (hist)delete hist;

		name = (string)TString::Format("hRelChPosVsEtaCMN_In_%d",i+1);
		hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecEtaCMNcorrected[i],512);
		hist->GetXaxis()->SetRangeUser(0,1);
		hist->GetYaxis()->SetTitle("Relative predicted Position");
		hist->GetXaxis()->SetTitle("#eta_{CMN corrected}");
		if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
		histSaver->SaveHistogram(hist);
		if (hist)delete hist;
	}

	name.str("");
	name.clear();
	name <<"hPredictedChannelPositionVsChi2";
	TH2F* hPredictedPositionVsChi2 = histSaver->CreateScatterHisto(name.str(),vecChi2,vecPredictedPosition,2048,128,0,inf,0,20,0);
	if (hPredictedPositionVsChi2){
		hPredictedPositionVsChi2->GetXaxis()->SetTitle("Predicted Channel Position");
		hPredictedPositionVsChi2->GetYaxis()->SetTitle("Max. #chi^{2}_{X,Y}");
		histSaver->SaveHistogram(hPredictedPositionVsChi2,false);
		delete hPredictedPositionVsChi2;
	}
	name.str("");
	name.clear();
	name <<"hRelativePredictedChannelPositionVsChi2";
	hPredictedPositionVsChi2 = histSaver->CreateScatterHisto(name.str(),vecChi2,vecRelPredictedPosition,512,128,0,1,0,20,0);
	if (hPredictedPositionVsChi2){
		hPredictedPositionVsChi2->GetXaxis()->SetTitle("relative Predicted Channel Position");
		hPredictedPositionVsChi2->GetYaxis()->SetTitle("Max. #chi^{2}_{X,Y}");
		histSaver->SaveHistogram(hPredictedPositionVsChi2,false);
		delete hPredictedPositionVsChi2;
	}
	name.str("");
	name.clear();
	name<<"hPredictedChannelPosition";
	TH1F* hPredictedPosition = histSaver->CreateDistributionHisto(name.str(),vecPredictedPosition,2048,histSaver->maxWidth,0,inf);
	if (hPredictedPosition){
		hPredictedPosition->GetXaxis()->SetTitle("Predicted Channel Position");
		histSaver->SaveHistogram(hPredictedPosition);
		delete hPredictedPosition;
	}
	name.str("");
	name.clear();
	name<<"hRelativePredictedChannelPosition";
	hPredictedPosition = histSaver->CreateDistributionHisto(name.str(),vecRelPredictedPosition,512,histSaver->maxWidth,0,1);
	if (hPredictedPosition){
		hPredictedPosition->GetXaxis()->SetTitle("relative Predicted Channel Position");
		histSaver->SaveHistogram(hPredictedPosition);
		delete hPredictedPosition;
	}
}

void TTransparentAnalysis::deleteHistograms() {
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		if(hLandau[clusterSize]) delete hLandau[clusterSize];
		if(hLandau2Highest[clusterSize])delete hLandau2Highest[clusterSize];
		if(hLandau1Highest[clusterSize])delete hLandau1Highest[clusterSize];
		if ( hEta[clusterSize]) delete hEta[clusterSize];
		if ( hEtaCMNcorrected[clusterSize]) delete hEtaCMNcorrected[clusterSize];
		if (hResidualChargeWeighted[clusterSize]) delete hResidualChargeWeighted[clusterSize];
		if (hResidualHighest2Centroid[clusterSize]) delete hResidualHighest2Centroid[clusterSize];
		if (hResidualHighestHit[clusterSize]) delete hResidualHighestHit[clusterSize];
	}
	delete hLandauMean;
	delete hLandauMP;
	delete hLandau2HighestMean;
	delete hLandau2HighestMP;
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
	cout << "TTransparentAnalysis Cutflow" << endl;
	cout << "number of events\t " << setw(8) << nEvents << endl;
	cout << "no valid silicon track\t-" << setw(8) << noValidTrack << endl;
	cout << "not in fid cut region \t-" << setw(8) << noFidCutRegion << endl;
	cout << "used for alignment    \t-" << setw(8) << usedForAlignment << endl;
	cout << "too high Chi2 value   \t-" << setw(8) << highChi2 <<endl;
//	cout << "used for si alignment\t-" << setw(8) << usedForSiliconAlignment << endl;
	cout << "region not on plane\t-" << setw(8) << regionNotOnPlane << endl;
	cout << "screened channel\t-" << setw(8) << screenedChannel << endl;
	cout << "saturated channel\t-" << setw(8) << saturatedChannel << endl;
	cout << "\t\t\t---------" << endl;
	cout << "total analyzed events\t " << setw(8) << nAnalyzedEvents << endl;
}

void TTransparentAnalysis::printEvent() {
	cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << endl;
	if (eventReader->useForAnalysis() == 0) {
		cout << "this track is not used for the analysis.." << endl;
		return;
	}
	cout << "predicted pos in lab system:\t" << this->predPosition << "\tpredicted perp position:\t" << this->predPerpPosition << endl;
	cout << "predicted pos in det system:\t" << this->positionInDetSystemMetric << endl;
	cout << "clustered analysis position in lab system:\t" << eventReader->getStripXPosition(subjectPlane,this->predPerpPosition,clusterCalcMode) << endl;
	cout << "clustered analysis position in det system:\t" << eventReader->getMeasuredPositionMetricSpace(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
	if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystemMetric, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) {
		cout << "this track did not pass the check.." << endl;
		return;
	}
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		cout << "transparent cluster of size " << clusterSize+1 << ":" << endl;
		cout << "\tpulse height:\t" << this->transparentClusters[clusterSize].getCharge() << endl;
		cout << "\teta:\t" << this->transparentClusters[clusterSize].getEta() << endl;
		cout << "\tresidual:\t" << this->getResidual(this->transparentClusters[clusterSize],this->clusterCalcMode,hEtaIntegrals[clusterSize]) << endl;
		cout << "\tcluster pos in det system:\t" << this->transparentClusters[clusterSize].getPosition(this->clusterCalcMode) << endl;
		cout << "\tcluster pos in lab system:\t" << eventReader->getPositionOfCluster(subjectDetector, this->transparentClusters[clusterSize], this->predPerpPosition, this->clusterCalcMode) << endl;
	}
	return;
}

Float_t TTransparentAnalysis::getResidual(TCluster cluster, TCluster::calculationMode_t clusterCalculationMode, TH1F* hEtaInt) {
	return eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,clusterCalculationMode, hEtaInt)-this->predPosition;
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
	cout << "\n\tpredicted channel: " << positionInDetSystemMetric;
	cout << "\n\tpredicted position: " << predPosition;
	cout << "\n\tcharge weighted position (TCluster::chargeWeighted): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid position (TCluster::highest2Centroid): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::highest2Centroid);
	cout << "\n\tcharge weighted residual: " << getResidual(cluster,TCluster::chargeWeighted);
	cout << "\n\thighest 2 centroid residual: " << getResidual(cluster,TCluster::highest2Centroid);
	cout << "\n\teta: " << cluster.getEta();
	if (hEtaIntegrals.size() != 0) {
		cout << "\n\teta corrected position (TCluster::corEta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::corEta,hEtaIntegrals[cluster.getClusterSize()]);
		cout << "\n\teta corrected residual (TCluster::corEta): " << getResidual(cluster,TCluster::corEta,hEtaIntegrals[cluster.getClusterSize()-1]);
		cout << "\n\teta corrected position (TCluster::eta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::eta,hEtaIntegrals[cluster.getClusterSize()]);
//		cout << "\n\teta corrected residual (TCluster::eta): " << getResidual(cluster,TCluster::eta);
	}
	cout << "\n\t";
	cluster.Print();
}

/** returns difference between cluster position and calculated position for a given cluster
 * @param cluster
 * @param clusterCalculationMode
 * @param hEtaInt
 * @author Lukas Baeni
 * @return
 */
void TTransparentAnalysis::createEventVector(Int_t startEvent) {

	nAnalyzedEvents = 0;
	regionNotOnPlane = 0;
	saturatedChannel = 0;
	screenedChannel = 0;
	noValidTrack = 0;
	noFidCutRegion = 0;
	usedForAlignment = 0;
	highChi2 =0;
	vecTransparentClusters.clear();
	eventNumbers.clear();
	vecEvents.clear();
	cout<<"Creating  Event Vector "<<endl;
	for (nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
		//		if (verbosity > 4) cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << eventReader<<endl;
		eventReader->LoadEvent(nEvent);
		if (eventReader->isValidTrack() == 0) {
			//		if (eventReader->useForAnalysis() == 0) {
			if (verbosity > 6) printEvent();
			noValidTrack++;
			continue;
		}
		if (eventReader->isInFiducialCut() == 0) {
			noFidCutRegion++;
			continue;
		}
		if (eventReader->useForAlignment() == true) {
			usedForAlignment++;
			continue;
		}
		transparentClusters.clear();
		if (!this->predictPositions(true)){
			if (verbosity>4) cout<< nEvent << ": Chi2 to high: "<< positionPrediction->getChi2()<<endl;
			highChi2++;
			continue;
		}


		if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystemChannelSpace, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false)
			continue;
		for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
			transparentClusters.push_back(this->makeTransparentCluster(this->eventReader, this->settings,subjectDetector, this->positionInDetSystemChannelSpace, clusterSize));
		}
		nAnalyzedEvents++;
		this->fillHistograms();
		if (verbosity > 4) printEvent();

		// save clusters for eta corrected analysis
		vecTransparentClusters.push_back(transparentClusters);
		eventNumbers.push_back(nEvent);
		vecEvents.push_back(eventReader->getEvent());
	}
}

TCluster TTransparentAnalysis::makeTransparentCluster(TTracking *reader,TSettings* set, UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
	// get channel and direction for clustering
	if (reader==0){
		cerr<<" TTransparentAnalysis::makeTransparentCluster TTracking == 0 "<<endl;
		return TCluster();
	}
	if (set == 0 ){
		cerr<<" TTransparentAnalysis::makeTransparentCluster TSettings == 0 "<<endl;
		return TCluster();
	}

	UInt_t centerChannel;
	int direction;
	direction = getSignedChannelNumber(centerPosition);
//	cout << "centerPosition: " << centerPosition << "\tdirection: " << direction << endl;
	centerChannel = TMath::Abs(direction);
	if (direction < 0) direction = -1;
	else direction = 1;
	Float_t cmNoise = reader->getCMNoise();

	// make cluster
	TCluster transparentCluster = TCluster(reader->getEvent_number(), det, -99, -99, TPlaneProperties::getNChannels(det),cmNoise);
	int currentChannel = centerChannel;
	for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
		if( currentChannel < 0 || currentChannel >= TPlaneProperties::getNChannelsDiamond() )
			cout<<reader->getEvent_number()<<": Cannot create channel with: det"<<det<<", centerPoisition: "<<centerPosition<< ", direction: "<<direction<<", centerChannel: "<<centerChannel<<" "<<iChannel<<endl;
		direction *= -1;
		currentChannel += direction * iChannel;
		Int_t adcValue=reader->getAdcValue(det,currentChannel);
		Float_t pedMean = reader->getPedestalMean(det,currentChannel,false);
		Float_t pedMeanCMN = reader->getPedestalMean(det,currentChannel,true);
		Float_t pedSigma = reader->getPedestalSigma(det,currentChannel,false);
		Float_t pedSigmaCMN = reader->getPedestalSigma(det,currentChannel,true);
		bool isScreened = set->isDet_channel_screened(det,currentChannel);
		transparentCluster.addChannel(currentChannel,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened);
//		transparentCluster.addChannel(currentChannel, reader->getRawSignal(det,currentChannel), reader->getRawSignalInSigma(det,currentChannel), reader->getAdcValue(det,currentChannel), reader->isSaturated(det,currentChannel), settings->isDet_channel_screened(det,currentChannel));
	}
	return transparentCluster;
}

