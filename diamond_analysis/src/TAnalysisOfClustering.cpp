/*
 * TDeadChannels.cpp
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#include "../include/TAnalysisOfClustering.hh"

TAnalysisOfClustering::TAnalysisOfClustering(TSettings *settings) {
	cout<<"\n\n\n\n**********************************************************"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"*********TAnalysisOfClustering::TAnalysisOfClustering*****"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"**********************************************************\n\n\n"<<endl;
	if(settings==0)
		settings=new TSettings();
	setSettings(settings);
	UInt_t runNumber=settings->getRunNumber();
	sys = gSystem;
	htmlClus= new THTMLCluster(settings);

	settings->goToClusterTreeDir();
	eventReader=new TADCEventReader(settings->getClusterTreeFilePath(),runNumber);
	histSaver=new HistogrammSaver();


	settings->goToClusterAnalysisDir();
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	htmlClus->setFileGeneratingPath(sys->pwd());
	settings->goToClusterTreeDir();
	verbosity=0;
	initialiseHistos();
	cout<<"end initialise"<<endl;
	settings=0;
}

TAnalysisOfClustering::~TAnalysisOfClustering() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	htmlClus->createClusterSize(vecClusterSizes,vecClusterSeedSizes,vecNumberOfClusters);
	htmlClus->createPulseHeightPlots(this->vecPHMeans);
	htmlClus->createContent();
	htmlClus->generateHTMLFile();
	delete htmlClus;
	settings->goToOutputDir();
}

void TAnalysisOfClustering::setSettings(TSettings* settings){
	this->settings=settings;
}

void TAnalysisOfClustering::doAnalysis(int nEvents)
{
	cout<<"analyze clustering results..."<<endl;
	//	eventReader->checkADC();
	if(nEvents==0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
		analyseEvent();
	}
	saveHistos();
}

void TAnalysisOfClustering::analyseEvent(){
	checkForDeadChannels();
	checkForSaturatedChannels();
	//		getBiggestHit();//not working
	analyseForSeeds();
	analyseCluster();
	compareCentroid_ChargeWeightedMean();
	analyse2ndHighestHit();
	analyseClusterPosition();
	createPHDistribution();
	etaInvestigation();
}
void TAnalysisOfClustering::checkForDeadChannels()
{
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		int numberOfSeeds=0;
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			Float_t sigma=eventReader->getPedestalSigma(det,ch);
			if(sigma==0){
				//cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
				continue;
			};
			Float_t signalInSigma=eventReader->getSignalInSigma(det,ch);

			if(signalInSigma>settings->getClusterSeedFactor(det,ch)){
				hSeedMap[det]->Fill(ch);
				//cout<<"Found a Seed "<<det<<" "<<ch<<" "<<adcValueInSigma<<" "<<eventReader->getCurrent_event()<<endl;
				numberOfSeeds++;
			}
		}
		hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}

}
void TAnalysisOfClustering::analyseForSeeds(){
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		int nClusters = eventReader->getNClusters(det);
		if(nClusters==1)
			hSeedMap2[det]->Fill(eventReader->getCluster(det,0).getHighestSignalChannel());
	}
}

void TAnalysisOfClustering::checkForSaturatedChannels()
{
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			if(eventReader->getAdcValue(det,ch)>=TPlaneProperties::getMaxSignalHeight(det)){
				hSaturatedChannels[det]->Fill(ch);
			}
		}
}

void TAnalysisOfClustering::initialiseHistos()
{

	if(verbosity>3)cout<<"1"<<endl;
	{
		stringstream histoName;
		histoName<<"hDiamond_Delta_CWM_BiggestHit";
		histo_CWM_biggestHit=new TH2F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6,10,0,9);
		histoName.str("");
		histoName<<"hDiamond_Delta_highest2Centroid_BiggestHit";
		histo_H2C_biggestHit=new TH1F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6);
	}
	size_t nDet = TPlaneProperties::getNDetectors();
	vecvecSignalLeftLeft.resize(nDet);
	vecvecSignalRightRight.resize(nDet);
	vecvecLeftEtaSignal.resize(nDet);
	vecvecRightEtaSignal.resize(nDet);
	vecvecSignalLeftOfHighest.resize(nDet);
	vecvecSignalRightOfHighest.resize(nDet);
	vecvecSignalHighest.resize(nDet);
	vecvecEta.resize(nDet);

	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hClusterPositionRelativeToNextIntegerCWM_"<<TPlaneProperties::getStringForDetector(det);
		if (verbosity>2) cout<<histName.str()<<endl;
		hRelativeClusterPositionCWM[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,TPlaneProperties::getNChannels(det)-1,1024,-.5,.5);
		histName.str("");
		histName.clear();
		histName<<"hClusterPositionRelativeToNextIntegerCorEta_"<<TPlaneProperties::getStringForDetector(det);
		if (verbosity>2) cout<<histName.str()<<endl;
		hRelativeClusterPositionCorEta[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,TPlaneProperties::getNChannels(det)-1,512,-.5,.5);
		histName.str("");
		histName.clear();
		histName<<"hClusterPositionRelativeToNextIntegerEta_"<<TPlaneProperties::getStringForDetector(det);
		if (verbosity>2) cout<<histName.str()<<endl;
		hRelativeClusterPositionEta[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,TPlaneProperties::getNChannels(det)-1,512,-.5,.5);
		histName.str("");
		histName.clear();
		histName<<"hAbsoluteClusterPostion_"<<TPlaneProperties::getStringForDetector(det);;
		if (verbosity>2) cout<<histName.str()<<endl;
		hClusterPosition[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),4096,0,TPlaneProperties::getNChannels(det)-1);
		histName.str("");
		histName.clear();
		histName<<"hEtaDistribution_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistribution[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,0,1);
		hEtaDistribution[det]->GetXaxis()->SetTitle("#eta");
		hEtaDistribution[det]->GetYaxis()->SetTitle("number of entries");
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionCMN_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionCMN[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,0,1);
		hEtaDistributionCMN[det]->GetXaxis()->SetTitle("#eta_{CMN-corrected}");
		hEtaDistributionCMN[det]->GetYaxis()->SetTitle("number of entries");
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsLeftChannel_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionVsLeftChannel[det] = new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,1,256,0,255);
		hEtaDistributionVsLeftChannel[det]->GetXaxis()->SetTitle("#eta");
		hEtaDistributionVsLeftChannel[det]->GetYaxis()->SetTitle("left channel of #eta position");
		hEtaDistributionVsLeftChannel[det]->GetZaxis()->SetTitle("number of entries #");
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsClusterSize_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionVsClusterSize[det] = new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,1,10,-.5,9.5);
		hEtaDistributionVsClusterSize[det]->GetXaxis()->SetTitle("#eta");
		hEtaDistributionVsClusterSize[det]->GetYaxis()->SetTitle("ClusterSize");
		hEtaDistributionVsClusterSize[det]->GetZaxis()->SetTitle("number of entries #");
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsCharge_"<<TPlaneProperties::getStringForDetector(det);
		Int_t maxCharge = TPlaneProperties::isDiamondDetector(det)?4096:512;
		hEtaDistributionVsCharge[det] = new TH2F(histName.str().c_str(),histName.str().c_str(),512,0,1,512,0,maxCharge);
		hEtaDistributionVsCharge[det]->GetXaxis()->SetTitle("#eta");
		hEtaDistributionVsCharge[det]->GetYaxis()->SetTitle("Charge of two highest Channels /ADC counts");
		hEtaDistributionVsCharge[det]->GetYaxis()->SetTitle("number of entries");

		histName.str("");
		histName.clear();
		histName<<"hEtaDistribution5Percent_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistribution5Percent[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,0,1);
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsSignalLeft_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionVsSignalLeft[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),128,0,1,128,0,TPlaneProperties::getMaxSignalHeight(det));
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsSignalRight_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionVsSignalRight[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),128,0,1,128,0,TPlaneProperties::getMaxSignalHeight(det));
		histName.str("");
		histName.clear();
		histName<<"hEtaDistributionVsSignalSum_"<<TPlaneProperties::getStringForDetector(det);
		hEtaDistributionVsSignalSum[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),128,0,1,128,0,TPlaneProperties::getMaxSignalHeight(det)*2);
		histName.str("");
		histName.clear();
		histName<<"hSignalLeftVsSignalRight"<<TPlaneProperties::getStringForDetector(det);
		hSignalLeftVsSignalRight[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),128,0,TPlaneProperties::getMaxSignalHeight(det),128,0,TPlaneProperties::getMaxSignalHeight(det));
	}
	if(verbosity>3)cout<<"2"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"SaturatedChannels_"<<TPlaneProperties::getStringForDetector(det)<<"";
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		if(det==8)hSaturatedChannels[det]->GetXaxis()->SetRangeUser(0,128);
	}
	if(verbosity>3)cout<<"3"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPositionOfallSeeds_"<<TPlaneProperties::getStringForDetector(det);
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		if(det==8)hSeedMap[det]->GetXaxis()->SetRangeUser(0,128);
	}
	if(verbosity>3)cout<<"4"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPositionOfHighestSeed_"<<TPlaneProperties::getStringForDetector(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		hSeedMap2[det]->GetXaxis()->SetTitle("Position of Highest Seed of a Cluster");
		if(det==8)hSeedMap2[det]->GetXaxis()->SetRangeUser(0,128);
	}
	if(verbosity>3)cout<<"5"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hNumberOfSeeds_in_"<<TPlaneProperties::getStringForDetector(det);
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
		hNumberOfSeeds[det]->GetXaxis()->SetTitle("Number Of Seeds in Cluster");
		hNumberOfSeeds[det]->GetYaxis()->SetTitle("Entries #");
	}
	if(verbosity>3)cout<<"6"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHitChannelInSigma";
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	if(verbosity>3)cout<<"7"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHitNextToBiggestHit_ChannelInSigma";
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	if(verbosity>3)cout<<"8"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"Channel_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHit";
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for(int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hClusterSize_Seed"<<settings->getClusterSeedFactor(det,0)<<"-Hit"<<settings->getClusterHitFactor(det,0)<<"_"<<TPlaneProperties::getStringForDetector(det);
		hClusterSize[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
		hClusterSize[det]->GetXaxis()->SetTitle("Number of Seeds and Hits in Cluster");
		hClusterSize[det]->GetYaxis()->SetTitle("Entries #");
		histoName.str("");
		histoName.clear();
		histoName<<"hClusterSeedSize_Seed"<<settings->getClusterSeedFactor(det,0)<<"-Hit"<<settings->getClusterHitFactor(det,0)<<"_"<<TPlaneProperties::getStringForDetector(det);
		hClusterSeedSize[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
		hClusterSeedSize[det]->GetXaxis()->SetTitle("Number of Seeds in Cluster");
		hClusterSeedSize[det]->GetYaxis()->SetTitle("Entries #");
		histoName.str("");
		histoName.clear();
		histoName<<"NumberOfClusters_"<<TPlaneProperties::getStringForDetector(det);
		hNumberOfClusters[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
	}
	if(verbosity>3)cout<<"10"<<endl;
	for (int det = 0; det < 9; det++) {
		int nbins = 250;
		Float_t min = 0.;
		Float_t max = 250.;

		stringstream histoName;
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "BiggestHitChannelInSigma";
		histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigma";
		histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitMinusBiggestHitPosition";
		histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),2,-2.,2.);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigmaLeft";
		histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigmaRight";
		histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << TPlaneProperties::getStringForDetector(det) << "BiggestHitMap";
		histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0.,255.);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "LeftChipBiggestHitChannelInSigma";
		histo_pulseheight_left_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "RightChipBiggestHitChannelInSigma";
		histo_pulseheight_right_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "LeftChipSecondBiggestHitChannelInSigma";
		histo_pulseheight_left_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "RightChipSecondBiggestHitChannelInSigma";
		histo_pulseheight_right_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
	}
	if(verbosity>3)cout<<"11"<<endl;
	for(int det=0;det<9;det++){//analayse2ndHighestHit
		stringstream histName;
		histName<<"h2ndBiggestHitSignal_"<<TPlaneProperties::getStringForDetector(det);
		if(det<8)
			h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,200);
		else
			h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1024);
		h2ndBiggestHitSignal[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster");
		h2ndBiggestHitSignal[det]->GetYaxis()->SetTitle("Entries #");
		histName.str("");
		histName<<"h2ndBiggestHitOverCharge_"<<TPlaneProperties::getStringForDetector(det);
		h2ndBiggestHitOverCharge[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,0.5);
		h2ndBiggestHitOverCharge[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster over Sum of all signals of cluster");
		h2ndBiggestHitOverCharge[det]->GetYaxis()->SetTitle("Entries #");
		histName.str("");
		histName<<"h2ndBiggestHitPosition_"<<TPlaneProperties::getStringForDetector(det);
		h2ndBiggestHitPosition[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),3,-1.5,1.5);
		h2ndBiggestHitPosition[det]->GetXaxis()->SetTitle("position of snd biggest hit in respect to biggest Hit");
		h2ndBiggestHitPosition[det]->GetYaxis()->SetTitle("Entries #");
		histName.str("");
		histName<<"hLeftHitOverLeftAndRight_"<<TPlaneProperties::getStringForDetector(det);
		hLeftHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
		hLeftHitOverLeftAndRight[det]->GetXaxis()->SetTitle("Q_L/(Q_R +Q_L)");
		histName.str("");
		histName<<"hDeltaLeftRightHitOverLeftAndRight_"<<TPlaneProperties::getStringForDetector(det);
		hDeltaLeftRightHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,-1,1);
		hDeltaLeftRightHitOverLeftAndRight[det]->GetXaxis()->SetTitle("(Q_L-Q_R)/(Q_R +Q_L)");
		histName.str("");
		histName<<"hSignal2ndHighestOverSignalHighest_"<<TPlaneProperties::getStringForDetector(det);
		hSignal2ndHighestOverSignalHighestRatio[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
		hSignal2ndHighestOverSignalHighestRatio[det]->GetXaxis()->SetTitle("Q_{2ndHighest}/Q_{Highest}");
	}
	if(verbosity>3)cout<<"12"<<endl;
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		stringstream histName;
		histName<<"hPulseHeightDistribution_"<<TPlaneProperties::getStringForDetector(det);
		float max=0;
		if(det==TPlaneProperties::getDetDiamond())
			max = 4098;
		else max = 512;
		hPHDistribution[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),512,0,max,10,-.5,9.5);

	}

	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		stringstream name;
		name<<"hBiggestHitSizeVsClusterSize_"<<TPlaneProperties::getStringForDetector(det);
		hBiggestHitVsClusterSize[det] = new TH2F(name.str().c_str(),name.str().c_str(),1024,0,TPlaneProperties::getMaxSignalHeight(det)*2,8,0.5,8.5);
		hBiggestHitVsClusterSize[det]->GetXaxis()->SetTitle("Signal of Biggest Hit in Cluster [adc counts]");
		hBiggestHitVsClusterSize[det]->GetYaxis()->SetTitle("ClusterSize");
	}
}



void TAnalysisOfClustering::saveEtaIntegrals(){
	stringstream etaCorFileName;
	etaCorFileName<<"etaCorrection."<<settings->getRunNumber()<<".root";
	TFile* file = TFile::Open(etaCorFileName.str().c_str());
	if(file){
		cout<<"can read file: '"<<etaCorFileName.str()<<"'"<<endl;
		if(verbosity>4&&verbosity%2==1){
			cout<<"Press a key and enter to confirm.\t"<<flush;
			char t; cin>>t;
		}
		file->Close();
		return;
	}
	cout<<"cannot read file: '"<<etaCorFileName.str()<<"' ===> CREATE new Eta file"<<endl;
	if(verbosity>4&&verbosity%2==1){
		cout<<"Press a key and enter to confirm.]\t"<<flush;
		char t; cin>>t;
	}
	file = new TFile(etaCorFileName.str().c_str(),"RECREATE");
	file->cd();
	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hEtaIntegral_"<<det;
		if(!hEtaDistribution[det])
			continue;
		TH1F *histo= TClustering::createEtaIntegral(hEtaDistribution[det],histName.str());
		file->cd();
		histo->Write();
		hEtaDistribution[det]->Write();
	}
	file->Close();
}



void TAnalysisOfClustering::saveHistos(){
	if (verbosity>2) cout<<"plot histo "<<histo_CWM_biggestHit->GetName();
	histSaver->SaveHistogram(histo_CWM_biggestHit);
	histo_CWM_biggestHit->Delete();
	if (verbosity>2) cout<<"plot histo "<<histo_H2C_biggestHit->GetName();
	histSaver->SaveHistogram(histo_H2C_biggestHit);
	histo_H2C_biggestHit->Delete();
	for(int det=0;det<9;det++){//analyse 2nd biggest Hit
		if (verbosity>2) cout<<"plot histo "<<det<<"  h2ndBiggestHitSignal_"<<TPlaneProperties::getStringForDetector(det);
		histSaver->SaveHistogram(h2ndBiggestHitSignal[det]);
		delete h2ndBiggestHitSignal[det];
		if (verbosity>2) cout<<"plot histo "<<det<<"  h2ndBiggestHitOverCharge_"<<TPlaneProperties::getStringForDetector(det);
		histSaver->SaveHistogram(h2ndBiggestHitOverCharge[det]);
		delete h2ndBiggestHitOverCharge[det];
		if (verbosity>2) cout<<"plot histo "<<h2ndBiggestHitPosition[det]->GetName()<<endl;
		histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
		histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
		histSaver->SaveHistogram(hLeftHitOverLeftAndRight[det]);
		histSaver->SaveHistogram(hDeltaLeftRightHitOverLeftAndRight[det]);
		histSaver->SaveHistogram(hSignal2ndHighestOverSignalHighestRatio[det]);
	}

	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogram(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
		histSaver->SaveHistogram(hSeedMap[det]);
		hSeedMap[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
		histSaver->SaveHistogram(hSeedMap2[det]);
		hSeedMap2[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
		histSaver->SaveHistogram(hNumberOfSeeds[det]);
		hNumberOfSeeds[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(hPulsHeightBiggestHit[det]);
		hPulsHeightBiggestHit[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(hPulsHeightNextBiggestHit[det]);
		hPulsHeightNextBiggestHit[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(hChannelBiggestHit[det]);
		hChannelBiggestHit[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(this->hClusterSize[det]);
		if (verbosity>2) cout<<"save: "<<hClusterSeedSize[det]->GetName()<<endl;
		histSaver->SaveHistogram(this->hClusterSeedSize[det]);
		histSaver->SaveHistogram(this->hNumberOfClusters[det]);
		vecClusterSizes.push_back(hClusterSize[det]->GetMean());
		vecClusterSeedSizes.push_back(hClusterSeedSize[det]->GetMean());
		vecNumberOfClusters.push_back(hNumberOfClusters[det]->GetMean());
		delete hClusterSize[det];
		delete hClusterSeedSize[det];
		delete hNumberOfClusters[det];
	}
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		if (verbosity>2) cout<<"Print : "<<hClusterPosition[det]->GetTitle()<< " "<<hClusterPosition[det]->GetEntries()<<endl;
		histSaver->SaveHistogram(this->hClusterPosition[det]);
		histSaver->SaveHistogram(this->hRelativeClusterPositionCWM[det]);
		histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionCWM[det]->ProjectionY());
		//		histSaver->SaveHistogram(this->hRelativeClusterPositionCorEta[det]);
		histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionCorEta[det]->ProjectionY());
		histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionEta[det]->ProjectionY());
		delete hClusterPosition[det];
		delete hRelativeClusterPositionCWM[det];
		delete hRelativeClusterPositionEta[det];
		delete hRelativeClusterPositionCorEta[det];
	}
	saveEtaIntegrals();
	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hEtaIntegral_"<<TPlaneProperties::getStringForDetector(det);;
		TH1F *histo= TClustering::createEtaIntegral(hEtaDistribution[det],histName.str());
		histSaver->SaveHistogram(histo);
		if(histo) delete histo;
		histSaver->SaveHistogram(this->hEtaDistribution[det]);
		histSaver->SaveHistogram(this->hEtaDistributionCMN[det]);
		histSaver->SaveHistogram(this->hEtaDistributionVsCharge[det]);
		for(int area = 0; area < settings->getNDiaDetectorAreas() && TPlaneProperties::isDiamondDetector(det); area++){
			TString name = TString::Format("hEtaDistributionVsLeftChannel_%d_Area%d",det,area);
			TH2F *hEtaDistributionVsLeftChannelArea = (TH2F*)hEtaDistributionVsLeftChannel[det]->Clone();
			Float_t yMin = settings->getDiaDetectorArea(area).first;
			Float_t yMax = settings->getDiaDetectorArea(area).second;
			hEtaDistributionVsLeftChannelArea->GetYaxis()->SetRangeUser(yMin,yMax);
			histSaver->SaveHistogram(hEtaDistributionVsLeftChannelArea);
			if(hEtaDistributionVsLeftChannelArea)
				delete hEtaDistributionVsLeftChannelArea;
		}
		histSaver->SaveHistogram(this->hEtaDistributionVsLeftChannel[det]);

		for(int i=1;(i<6&&TPlaneProperties::isDiamondDetector(det))||i<3;i++){
			TString name =  TString::Format("hEtaDistribution_ClusterSize%d",i);
			Int_t bin = hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i);
			TH1F* hEtaOneClusterSize = (TH1F*)hEtaDistributionVsClusterSize[det]->ProjectionX(name,bin,bin);
			hEtaOneClusterSize->SetTitle(name);
			histSaver->SaveHistogram(hEtaOneClusterSize);
			delete hEtaOneClusterSize;
			if(i==1){
				name =  TString::Format("hEtaDistribution_ClusterSize_%d-%d",i,i+1);
				bin = hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i);
				Int_t bin2 =hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i+1);
				TH1F* hEtaTwoClusterSize = (TH1F*)hEtaDistributionVsClusterSize[det]->ProjectionX(name,bin,bin2);
				hEtaTwoClusterSize->SetTitle(name);
				histSaver->SaveHistogram(hEtaTwoClusterSize);
				delete hEtaTwoClusterSize;
			}
		}
		histSaver->SaveHistogram(this->hEtaDistributionVsClusterSize[det]);

		histSaver->SaveHistogram(this->hEtaDistribution5Percent[det]);
		histSaver->SaveHistogram(this->hEtaDistributionVsSignalLeft[det]);
		histSaver->SaveHistogram(this->hEtaDistributionVsSignalRight[det]);
		histSaver->SaveHistogram(this->hEtaDistributionVsSignalSum[det]);
		hSignalLeftVsSignalRight[det]->GetXaxis()->SetTitle("signalRight");
		hSignalLeftVsSignalRight[det]->GetYaxis()->SetTitle("signalLeft");
		histSaver->SaveHistogram(this->hSignalLeftVsSignalRight[det]);
		if(hEtaDistribution[det])delete hEtaDistribution[det];
	}
	savePHHistos();
	//    for (int det = 0; det < 9; det++) {
	//		cout << "saving histogram" << this->histo_pulseheight_sigma[det]->GetName() << ".." << endl;
	//        histSaver->SaveHistogram(this->histo_pulseheight_sigma[det]);
	//		cout << "saving histogram" << this->histo_pulseheight_sigma_second[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_pulseheight_sigma_second[det]);
	////		cout << "saving histogram" << this->histo_pulseheight_sigma125[det]->GetName() << ".." << endl;
	////		histSaver->SaveHistogram(this->histo_pulseheight_sigma125[det]);
	//		cout << "saving histogram" << this->histo_second_biggest_hit_direction[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_second_biggest_hit_direction[det]);
	//		cout << "saving histogram" << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_biggest_hit_map[det]);
	//		cout << "saving histogram" << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma[det]);
	//		cout << "saving histogram" << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma_second[det]);
	//		cout << "saving histogram" << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma[det]);
	//		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
	//		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma_second[det]);
	//        delete histo_pulseheight_sigma[det];
	//		delete histo_pulseheight_sigma_second[det];
	////		delete histo_pulseheight_sigma125[det];
	//		delete histo_second_biggest_hit_direction[det];
	//		delete histo_biggest_hit_map[det];
	//		delete histo_pulseheight_left_sigma[det];
	//		delete histo_pulseheight_left_sigma_second[det];
	//		delete histo_pulseheight_right_sigma[det];
	//		delete histo_pulseheight_right_sigma_second[det];
	//    }
	saveEtaInvestigationHistos();
}

void TAnalysisOfClustering::saveEtaDividedHistos(TH3F* histo3dLeft,TH3F* histo3dRight, TH2F* histo2Left, TH2F* histo2Right,string name_comparision, Float_t etaWidth){
	if(etaWidth<=0) {
		cout<<"etaWidth <= 0"<<endl;
		return;
	}
	int minZbin,maxZbin, minXbin,maxXbin;
	int minZ,minX,maxX,maxZ;
	Float_t eta1,eta2;
	stringstream name;
	name<<setfill('0');
	vector<TProfile *> vecLeftProfiles;
	vector<TProfile *> vecRightProfiles;

	for(Float_t eta=0;eta<=1;eta+=etaWidth){
		eta1 = eta<1?eta:-1000;
		eta2 = eta<1?eta+etaWidth:1000;
//		cout<<"eta in "<<eta1<<"-"<<eta2<<endl;
		/**** 2D  ***/

		name.str("");
		name.clear();
		minZbin = histo3dLeft->GetZaxis()->FindBin(eta1);
		maxZbin = histo3dLeft->GetZaxis()->FindBin(eta2);
//		cout<<"eta "<<eta<<": "<<minZbin<<"-"<<maxZbin<<endl;
		histo3dLeft->GetZaxis()->SetRange(minZbin, maxZbin);
		minZ = histo3dLeft->GetZaxis()->GetBinLowEdge(minZbin)*100;
		minZ=minZ<0?0:minZ;
		maxZ = histo3dLeft->GetZaxis()->GetBinLowEdge(maxZbin+1)*100;
		name<<histo3dLeft->GetName()<<"_eta_"<< setw(3)<<(minZ)<<"_"<<setw(3)<<(maxZ);
		histo3dLeft->GetZaxis()->SetBit(TAxis::kAxisRange);
		TH2F* histo2dLeft = (TH2F*)histo3dLeft->Project3D("yx");
		histo2dLeft->SetName(name.str().c_str());
		histo2dLeft->SetTitle(name.str().c_str());
		name.str("");
		name.clear();
		name<<histo2dLeft->GetName()<<"_pfx";
		TProfile * hLeft_pfx = histo2dLeft->ProfileX(name.str().c_str());
		hLeft_pfx->GetXaxis()->SetTitle("signal 'Left'");
		hLeft_pfx->GetYaxis()->SetTitle("avrg. signal 'LeftLeft'");
		hLeft_pfx->Draw("goff");
		if(hLeft_pfx->GetXaxis()->GetXmax()>2000)
			hLeft_pfx->GetXaxis()->SetRangeUser(0,2000);
		vecLeftProfiles.push_back(hLeft_pfx);

		name.str("");
		name.clear();
		minZbin = histo3dRight->GetZaxis()->FindBin(eta1);
		maxZbin = histo3dRight->GetZaxis()->FindBin(eta2);
//		cout<<eta<<": "<<minZbin<<"-"<<maxZbin<<endl;
		histo3dRight->GetZaxis()->SetRange(minZbin, maxZbin);
		minZ = histo3dRight->GetZaxis()->GetBinLowEdge(minZbin)*100;
		minZ=minZ<0?0:minZ;
		maxZ = histo3dRight->GetZaxis()->GetBinLowEdge(maxZbin+1)*100;
		name<<histo3dRight->GetName()<<"_eta_"<<setw(3)<<(minZ)<<"_"<<setw(3)<<(maxZ);
		histo3dRight->GetZaxis()->SetBit(TAxis::kAxisRange);

		TH2F* histo2dRight = (TH2F*)histo3dRight->Project3D("yx");
		histo2dRight->SetName(name.str().c_str());
		histo2dRight->SetTitle(name.str().c_str());
		name.str("");
		name.clear();
		name<<histo2dRight->GetName()<<"_pfx";
		TProfile * hRight_pfx = histo2dRight->ProfileX(name.str().c_str());
		hRight_pfx->GetXaxis()->SetTitle("signal 'Right'");
		hRight_pfx->GetYaxis()->SetTitle("avrg. signal 'RightRight'");
		hRight_pfx->Draw("goff");
		if(hRight_pfx->GetXaxis()->GetXmax()>2000)
			hRight_pfx->GetXaxis()->SetRangeUser(0,2000);
		vecRightProfiles.push_back(hRight_pfx);

		/**** 1D  ***/
		name.str("");
		name.clear();
		minXbin = histo2Left->GetXaxis()->FindBin(eta1);
		maxXbin = histo2Left->GetXaxis()->FindBin(eta2);
		minX = histo2Left->GetXaxis()->GetBinLowEdge(minXbin)*100;
		minX=minX<0?0:minX;
		maxX = histo2Left->GetXaxis()->GetBinLowEdge(maxXbin+1)*100;
		name<<histo2Left->GetName()<<"_eta_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
		TH1F *histo1dLeft = (TH1F*)histo2Left->ProjectionY(name.str().c_str(),minXbin,maxXbin);
		histo1dLeft->SetTitle(name.str().c_str());

		name.str("");
		name.clear();
		minXbin = histo2Right->GetXaxis()->FindBin(eta1);
		maxXbin = histo2Right->GetXaxis()->FindBin(eta2);
		minX = histo2Right->GetXaxis()->GetBinLowEdge(minXbin)*100;
		minX = minX<0?0:minX;
		maxX = histo2Right->GetXaxis()->GetBinLowEdge(maxXbin+1)*100;
		name<<histo2Right->GetName()<<"_eta_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
		TH1F *histo1dRight = (TH1F*)histo2Right->ProjectionY(name.str().c_str(),minXbin,maxXbin);
		histo1dRight->SetTitle(name.str().c_str());

//		cout<<"save "<<histo1dRight->GetName()<<endl;
		histSaver->SaveHistogram(histo1dRight);
//		cout<<"save "<<histo1dLeft->GetName()<<endl;
		histSaver->SaveHistogram(histo1dLeft);
//		cout<<"saving "<<histo2dLeft->GetName()<<endl;;
		histSaver->SaveHistogram(histo2dLeft,false);
//		cout<<"saving "<<histo2dRight->GetName()<<endl;
		histSaver->SaveHistogram(histo2dRight,false);
		name.str("");
		name.clear();
		name<<"cComparision_"<<name_comparision<<"_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
		histo1dRight->SetLineColor(kBlue);
		histo1dLeft->SetLineColor(kRed);
		histSaver->SaveTwoHistos(name.str(),histo1dRight,histo1dLeft,1);

		if (histo1dRight) delete histo1dRight;
		if(histo2dRight) delete histo2dRight;
		if (histo1dLeft) delete histo1dLeft;
		if(histo2dLeft) delete histo2dLeft;
	}
	for(UInt_t i=0; i<vecLeftProfiles.size()-1&&i<vecRightProfiles.size()-1;i++){
		TProfile* pLeft = vecLeftProfiles.at(i);
		Int_t j = vecRightProfiles.size()-2-i;
		TProfile* pRight = vecRightProfiles.at(j);
		Int_t nbins = pRight->GetNbinsX();
		Float_t xlow = pRight->GetBinLowEdge(1);
		Float_t xup = pRight->GetBinLowEdge(nbins+1);
		TString pname = TString::Format("pDiff_RightMinusLeft_%s_%d",name_comparision.c_str(),i);
		TProfile* pDiff = new TProfile(pname,pname,nbins,xlow,xup);
		pDiff->Add(pLeft,pRight,-1.,1.);
		pDiff->SetLineColor(kBlue);
		pLeft->SetLineColor(kRed);
		pLeft->SetMarkerColor(kRed);
		pRight->SetLineColor(kGreen);
		pRight->SetMarkerColor(kGreen);
		Float_t minY = pDiff->GetBinContent(pDiff->GetMinimumBin());
		TCanvas *c1 = new TCanvas(TString::Format("cProfile_%s_%d",name_comparision.c_str(),i));

//		cout<<"compare: "<<pLeft->GetName()<<" vs. "<<pRight->GetName()<<":\t"<<c1->GetName()<<endl;
		c1->cd();
		TString title = TString::Format("hProfile_%s_%d",name_comparision.c_str(),i);
		TH1F* pTitle = new TH1F(title,title,pLeft->GetNbinsX(),pLeft->GetXaxis()->GetXmin(),pLeft->GetXaxis()->GetXmax());
		Float_t maxX = 1000;
		Int_t maxBin = pLeft->FindBin(maxX);
		Double_t maxLeft =0;
		Double_t maxRight=0;
		for(Int_t i=0;i<maxBin;i++){
			maxLeft = TMath::Max(maxLeft,pLeft->GetBinContent(i));
			maxRight = TMath::Max(maxRight,pRight->GetBinContent(i));
		}
		pTitle->Draw();
		pTitle->GetXaxis()->SetRangeUser(0,maxX);
//		pRight->GetXaxis()->SetRangeUser(0,maxX);
//		pLeft->GetXaxis()->SetRangeUser(0,maxX);
		Float_t maxY = 1.1*TMath::Max(maxLeft,maxRight);
		pTitle->SetMaximum(maxY);
		pTitle->GetYaxis()->SetRangeUser(1.1*minY,maxY);

		TString rightName = pRight->GetTitle();
		TString leftName = pLeft->GetTitle();
		if (rightName.Contains("rightright",TString::kIgnoreCase)){
			pTitle->GetXaxis()->SetTitle("signal 'Left / Right'");
			pTitle->GetYaxis()->SetTitle("avrg. signal adjacent ('LeftLeft / RightRight')");
		}
		else if (rightName.Contains("highest",TString::kIgnoreCase)){
			pTitle->GetXaxis()->SetTitle("highest Signal");
			pTitle->GetYaxis()->SetTitle("avrg. signal adjacent (Left/Right)");
		}
		else{
			pTitle->GetXaxis()->SetTitle("signal 1");
			pTitle->GetYaxis()->SetTitle("avrg. signal adjacent (signal2)");

		}
		pTitle->Draw();
		pRight->Draw("same");
		pLeft->Draw("same");
		pDiff->Draw("same");

		TLegend *leg = c1->BuildLegend();
		leg->Clear();
		leg->AddEntry(pRight);
		leg->AddEntry(pLeft);
		leg->AddEntry(pDiff);
		leg->SetFillColor(kWhite);
		Float_t x1NDC = leg->GetX1NDC();
		Float_t x2NDC = leg->GetX2NDC();
		Float_t val = x1NDC;
		x1NDC = -(x2NDC-x1NDC) + x1NDC;
		x2NDC=val;
		leg->SetX1NDC(.15);//x1NDC);
		leg->SetX2NDC(.5);//x2NDC);

		leg->Draw();
		histSaver->SaveCanvas(c1);
		histSaver->SaveHistogram(pLeft);
		histSaver->SaveHistogram(pRight);
		histSaver->SaveHistogram(pDiff);
		delete c1;
		if (pDiff) delete pDiff;
	}

	for(UInt_t i=0; i<vecLeftProfiles.size()&&i<vecRightProfiles.size();i++){
		delete vecLeftProfiles.at(i);
		delete vecRightProfiles.at(i);
	}
	cout<<"saving "<<histo3dRight->GetName()<<endl;
	histSaver->SaveHistogramROOT(histo3dRight);
	cout<<"saving "<<histo3dLeft->GetName()<<endl;
	histSaver->SaveHistogramROOT(histo3dLeft);

//	cout<<"saving "<<histo2Left->GetName()<<endl;
	histSaver->SaveHistogram(histo2Left);
//	cout<<"saving "<<histo2Right->GetName()<<endl;
	histSaver->SaveHistogram(histo2Right);
}

void TAnalysisOfClustering::saveEtaInvestigationHistos(){
	stringstream name;
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
		name.str("");
		name.clear();
		name<<"hLeftLeft_"<<TPlaneProperties::getStringForDetector(det);
		TH2F* histo2Left = histSaver->CreateScatterHisto(name.str(),vecvecSignalLeftLeft.at(det),vecvecEta.at(det),512,0,1);
		histo2Left->GetXaxis()->SetTitle("#eta");
		histo2Left->GetYaxis()->SetTitle("Signal 'LeftLeft'");

		name.str("");
		name.clear();
		name<<"hRightRight_"<<TPlaneProperties::getStringForDetector(det);
		TH2F* histo2Right = histSaver->CreateScatterHisto(name.str(),vecvecSignalRightRight.at(det),vecvecEta.at(det),512,0,1);
		histo2Right->GetXaxis()->SetTitle("#eta");
		histo2Right->GetYaxis()->SetTitle("Signal 'RightRight'");

		name.str("");name.clear();
		name<<"hLeftLeft_vs_Left_3D_"<<TPlaneProperties::getStringForDetector(det);
		Float_t inf = std::numeric_limits<float>::infinity();
		TH3F* histo3dLeft = histSaver->Create3DHisto(name.str(),vecvecLeftEtaSignal.at(det),vecvecSignalLeftLeft.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
		histo3dLeft->GetXaxis()->SetTitle("Signal 'Left'");
		histo3dLeft->GetYaxis()->SetTitle("Signal 'LeftLeft'");
		histo3dLeft->GetZaxis()->SetTitle("#eta");
		/*****************/


		name.str("");
		name.clear();
		name<<"hRightRight_vs_Right_3D_"<<TPlaneProperties::getStringForDetector(det);
		TH3F* histo3dRight = histSaver->Create3DHisto(name.str(),vecvecRightEtaSignal.at(det),vecvecSignalRightRight.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
		histo3dRight->GetXaxis()->SetTitle("Signal 'Right'");
		histo3dRight->GetYaxis()->SetTitle("Signal 'RightRight'");
		histo3dRight->GetZaxis()->SetTitle("#eta");

		name.str("");
		name.str();
		name<<"LeftLeft-RightRight_"<<TPlaneProperties::getStringForDetector(det);
		this->saveEtaDividedHistos(histo3dLeft,histo3dRight,histo2Left,histo2Right,name.str());


		if (histo3dRight) delete histo3dRight;
		if (histo3dLeft) delete histo3dLeft;
		if (histo2Left) delete histo2Left;
		if (histo2Right) delete histo2Right;

		/**************************/

		name.str("");
		name.clear();
		name<<"hLeft_vs_Highest_3D_"<<TPlaneProperties::getStringForDetector(det);
		histo3dLeft = histSaver->Create3DHisto(name.str(),vecvecSignalHighest.at(det),vecvecSignalLeftOfHighest.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
		name.str("");
		name.clear();
		name<<"hRight_vs_Highest_3D_"<<TPlaneProperties::getStringForDetector(det);
		histo3dRight = histSaver->Create3DHisto(name.str(),vecvecSignalHighest.at(det),vecvecSignalRightOfHighest.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
		name.str("");
		name.clear();

		name<<"hLeftOfHighest_"<<TPlaneProperties::getStringForDetector(det);
		histo2Left = histSaver->CreateScatterHisto(name.str(),vecvecSignalLeftOfHighest.at(det),vecvecEta.at(det),512,0,1);
		histo2Left->GetXaxis()->SetTitle("#eta");
		histo2Left->GetYaxis()->SetTitle("Signal left of highest channel");

		name.str("");
		name.clear();
		name<<"hRightOfHighest_"<<TPlaneProperties::getStringForDetector(det);
		histo2Right = histSaver->CreateScatterHisto(name.str(),vecvecSignalRightOfHighest.at(det),vecvecEta.at(det),512,0,1);
		histo2Right->GetXaxis()->SetTitle("#eta");
		histo2Right->GetYaxis()->SetTitle("Signal right of highest channel");

		name.str("");
		name.clear();
		name<<"LeftHighest-RightHighest_"<<TPlaneProperties::getStringForDetector(det);
		this->saveEtaDividedHistos(histo3dLeft,histo3dRight,histo2Left,histo2Right,name.str());

		if (histo3dRight) delete histo3dRight;
		if (histo3dLeft) delete histo3dLeft;
		if (histo2Left) delete histo2Left;
		if (histo2Right) delete histo2Right;



	}
}
void TAnalysisOfClustering::compareCentroid_ChargeWeightedMean()
{
	bool check=true;
	for(int det=0;det<9;det++)
		check=eventReader->getNClusters(det)==1;
	if(check==true){
		TCluster cluster = eventReader->getCluster(8,0);
		Float_t xCWM=cluster.getChargeWeightedMean();
		Float_t xHit=(Float_t)cluster.getHighestSignalChannel();
		Float_t xH2C=(Float_t)cluster.getHighest2Centroid();
		Float_t delta=xCWM-xHit;
		this->histo_CWM_biggestHit->Fill(delta,cluster.size());
		delta = xH2C - xHit;
		this->histo_H2C_biggestHit->Fill(delta);
		//		if(eventReader->getNClusters(8)>=1){
		//			Float_t charge = eventReader->getCluster(8,0).getCharge();
		//			Float_t signal2ndHighestHit=eventReader->getCluster(8,0).getCharge(2)-eventReader->getCluster(8,0).getCharge(1);
		//			Float_t q =signal2ndHighestHit/charge;
		//		}


	}
}

void TAnalysisOfClustering::analyseClusterPosition()
{
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
			Float_t posCWM = eventReader->getEvent()->getPosition(det,cl,TCluster::chargeWeighted);
			Int_t chNo = (Int_t)(posCWM+0.5);
			Float_t relPos = posCWM - chNo;
			hClusterPosition[det]->Fill(posCWM);
			hRelativeClusterPositionCWM[det]->Fill(chNo+0.5,relPos);
			//			Float_t eta =
			TCluster cluster = eventReader->getCluster(det,cl);
			UInt_t highestClPos=cluster.getHighestHitClusterPosition();
			UInt_t nextHighestClPos=cluster.getHighestSignalNeighbourClusterPosition(highestClPos);
			if(nextHighestClPos==9999){
				//cout<<"\nnext highest=9999: "<<highestClPos<<" "<<eventReader->getCluster(det,cl).getClusterSize()<<endl;;
				//eventReader->getCluster(det,cl).Print(1);
				continue;
			}

			Float_t signalLeft,signalRight,adcRight,adcLeft,pedRight,pedLeft;
			UInt_t leftClPos,rightClPos;
			if(highestClPos>nextHighestClPos) {
				leftClPos=nextHighestClPos;
				rightClPos=highestClPos;
			}
			else{
				rightClPos=nextHighestClPos;
				leftClPos=highestClPos;
			}
			int clusterSize =cluster.getClusterSize();
			if(settings->isMaskedCluster(det,cluster,true))
				return;
			signalLeft= cluster.getSignal(leftClPos);
			signalRight=cluster.getSignal(rightClPos);
			adcLeft= cluster.getAdcValue(leftClPos);
			adcRight= cluster.getAdcValue(rightClPos);
			pedLeft= cluster.getPedestalMean(leftClPos);
			pedRight= cluster.getPedestalMean(rightClPos);
			Float_t a= 0.03;
			Float_t adcLeftReal=adcLeft/(1-a);
			Float_t adcRightReal=(-a*(1-a)*adcLeft+(1-a)*adcRight)/(1-a)/(1-a);
			Float_t signalAdcLeft = adcLeftReal-pedLeft;
			Float_t signalAdcRight= adcRightReal-pedRight;
			Float_t eta2=(signalAdcRight)/(signalAdcLeft+signalAdcRight);
			Int_t leftChannel=-1;
			Float_t eta= cluster.getEta(leftChannel);
			Float_t etaCmnCorrected = cluster.getEta(true);
			Float_t charge = cluster.getCharge(2,true);
			if(verbosity>3){
				cout<<"charge of 2: "<<charge<<"\t"<<flush;
				cluster.Print();
			}
			//			Float_t eta3= signalRightReal/(signalLeftReal+signalRightReal);
			//			cout<<nEvent<<" "<<eta<<" "<<eta1<<" "<<eta2<<" "<<eta3<<endl;
			if(hEtaDistribution[det])hEtaDistribution[det]->Fill(eta);
			if(hEtaDistributionCMN[det])hEtaDistributionCMN[det]->Fill(etaCmnCorrected);
			if(hEtaDistributionVsLeftChannel[det]) hEtaDistributionVsLeftChannel[det]->Fill(eta,leftChannel);
			if(hEtaDistributionVsClusterSize[det]) hEtaDistributionVsClusterSize[det]->Fill(eta,clusterSize);
			if(hEtaDistributionVsCharge[det]) hEtaDistributionVsCharge[det]->Fill(eta,charge);
			hEtaDistribution5Percent[det]->Fill(eta2);
			hSignalLeftVsSignalRight[det]->Fill(signalRight,signalLeft);

			//			cout<<nextHighestClPos<<"<"<<highestClPos<<"\t"<<signalLeft<<" < "<<signalRight<<"\t"<<eta<<endl;
			hEtaDistributionVsSignalLeft[det]->Fill(eta,signalLeft);
			hEtaDistributionVsSignalRight[det]->Fill(eta,signalRight);
			hEtaDistributionVsSignalSum[det]->Fill(eta,signalLeft+signalRight);
			TH1F *hEtaIntegral=eventReader->getEtaIntegral(det);
			Float_t posCorEta=  eventReader->getEvent()->getPosition(det,cl,TCluster::corEta,hEtaIntegral);
			chNo = (UInt_t)(posCorEta+0.5);
			relPos = posCorEta - chNo;
			if(verbosity) printf("%5d %3d %5.1f %5.1f\n",nEvent,chNo,posCWM,posCorEta);
			hRelativeClusterPositionCorEta[det]->Fill(chNo+0.5,relPos);
			Float_t posEta=  eventReader->getEvent()->getPosition(det,cl,TCluster::eta);
			chNo = (UInt_t)(posEta+0.5);
			relPos = posEta - chNo;
			hRelativeClusterPositionEta[det]->Fill(chNo+0.5,relPos);
		}
	}
}

void TAnalysisOfClustering::analyseCluster()
{
	for(int det=0;det<9;det++){
		hNumberOfClusters[det]->Fill(eventReader->getNClusters(det));
		for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
			UInt_t clSize = eventReader->getClusterSize(det,cl);
			hClusterSize[det]->Fill(clSize);
			hClusterSeedSize[det]->Fill(eventReader->getClusterSeedSize(det,cl));
			Float_t biggestSignal = eventReader->getCluster(det,cl).getHighestSignal();
			clSize=clSize>8?8:clSize;
			hBiggestHitVsClusterSize[det]->Fill(biggestSignal,clSize);
		}
	}

}


void TAnalysisOfClustering::analyse2ndHighestHit(){
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		//		if(nEvent==20){
		//			cout<<nEvent<<":"<<endl;
		//			eventReader->setVerbosity(10);
		//		}
		//		else eventReader->setVerbosity(0);
		Float_t nClusters=eventReader->getNClusters(det);

		for(UInt_t cl=0;cl<nClusters;cl++){
			TCluster cluster=eventReader->getCluster(det,cl);
			if(settings->isMaskedCluster(det,cluster,true))
				continue;
			if(cluster.size()==0){
				cout<<nEvent<<" "<<cl<<" "<<eventReader->getNClusters(det)<<" "<<cluster.size()<<endl;
				cluster.Print();
				continue;//why can this happen???
			}
			UInt_t highestChannel = cluster.getHighestSignalChannel();
			Float_t signalLeft = cluster.getSignalOfChannel(highestChannel-1);
			Float_t signalRight = cluster.getSignalOfChannel(highestChannel+1);
			if(signalLeft<0 && verbosity>2)
				cout<<"signalLeft is smaller than 0"<<endl;
			if(signalRight<0 && verbosity>2)
				cout<<"signalLeft is smaller than 0"<<endl;
			Float_t signalHighest = cluster.getHighestSignal();
			Float_t signal2ndHighest;
			Float_t deltaSignals = signalLeft-signalRight;
			if(deltaSignals<0)// right channel higher left channel
				signal2ndHighest=signalRight;
			else
				signal2ndHighest=signalLeft;
			Float_t sumSignals = signalLeft+signalRight;
			if(signalLeft==0&&signalRight==0)continue;
			if(sumSignals==0)continue;
			Float_t allCharge=cluster.getCharge(false);
			//			Float_t charge = cluster.getCharge(false);
			//			if (cluster.isHit())
			Float_t signalRatio=signal2ndHighest/signalHighest;
			if(signalRatio>1){
				cout<<"Ratio>1:"<<signal2ndHighest<<" "<<signalHighest<<endl;
				cluster.Print();
			}
			else
				hSignal2ndHighestOverSignalHighestRatio[det]->Fill(signalRatio);
			Float_t ratio;
			if(cluster.size()!=1){
				if(signalLeft>signalRight){
					ratio=signalLeft/allCharge;
					if(ratio>0.5||allCharge==0||ratio!=ratio){
						cout<<"\n2ndBiggestHitOverCharge>0.5: left "<<signalLeft<<" "<<allCharge<<endl;
						cluster.Print();
					}
					else{
						//					cout<<nEvent<<" "<<cl<<" "<<ratio<<endl;
						h2ndBiggestHitOverCharge[det]->Fill(ratio);
					}
					h2ndBiggestHitSignal[det]->Fill(signalLeft);
				}
				else{
					ratio=signalRight/allCharge;
					if(ratio>0.5||allCharge==0||ratio!=ratio){
						cout<<"\n2ndBiggestHitOverCharge>0.5: right"<<signalRight<<" "<<allCharge<<endl;
						cluster.Print();
					}
					else{
						//					cout<<nEvent<<" "<<cl<<" "<<ratio<<endl;
						h2ndBiggestHitOverCharge[det]->Fill(ratio);
					}

					h2ndBiggestHitSignal[det]->Fill(signalRight);
				}
			}
			if(signalLeft>signalRight){
				h2ndBiggestHitPosition[det]->Fill(-1);
			}
			else if(signalLeft<signalRight)
				h2ndBiggestHitPosition[det]->Fill(+1);
			else
				h2ndBiggestHitPosition[det]->Fill(0);
			ratio = (deltaSignals)/(sumSignals);
			if (ratio<1&&ratio>-1&&sumSignals!=0)
				hDeltaLeftRightHitOverLeftAndRight[det]->Fill(ratio);
			else {
				if(TMath::Abs(ratio)>1){
					//					cout<<"hDeltaLeftRightHitOverLeftAndRight "<<det<<" "<<cl<<" "<<deltaSignals<<" "<<sumSignals<<endl;
					//					cluster.Print();
				}
			}
			ratio = signalLeft/sumSignals;
			if(signalLeft>0&&sumSignals>0&&ratio<1&&ratio>0)hLeftHitOverLeftAndRight[det]->Fill((signalLeft)/(sumSignals));
		}
	}
}



void TAnalysisOfClustering::etaInvestigation(){
	for(int det=0;det<9;det++){
		for(UInt_t cl = 0; cl< eventReader->getNClusters(det);cl++){
			TCluster cluster = eventReader->getCluster(det,cl);
//			cluster.Print(1);
			Int_t leftEtaChannel =-1;
			Float_t eta = cluster.getEta(leftEtaChannel,false);
			Int_t leftEtaClusterPosition = cluster.getClusterPosition(leftEtaChannel);
			Float_t signalLeftLeft = cluster.getSignal(leftEtaClusterPosition-1);
			Float_t leftEtaSignal = cluster.getSignal(leftEtaClusterPosition);
			Float_t rightEtaSignal = cluster.getSignal(leftEtaClusterPosition+1);
			Float_t signalRightRight = cluster.getSignal(leftEtaClusterPosition+2);
			Float_t highestSignal = cluster.getHighestSignal();
			Int_t highestSignalChannel = cluster.getHighestSignalChannel();
			Int_t highestSignalClusterPos = cluster.getClusterPosition(highestSignalChannel);
			Float_t signalLeftOfHighest = cluster.getSignal(highestSignalClusterPos-1);
			Float_t signalRightOfHighest = cluster.getSignal(highestSignalClusterPos+1);
//			cout<< signalLeftOfEta<<" "<< leftEtaSignal<<" "<<rightEtaSignal<<" "<<signalRightOfEta<<"  ----  "<<signalLeftOfHighest<<" "<< highestSignal<<" "<<signalRightOfHighest<<"\n"<<endl;
			//		Float_t signalLeftOfEta = eventReader->
			if (signalLeftLeft < -100){
				int clusPos = leftEtaClusterPosition-1;
				cout<<"\n"<< nEvent<<", Problem signalLeftLeft: "<<signalLeftLeft<<" "<<leftEtaChannel<<": "<<cluster.getPedestalMean(clusPos)<<", "<<cluster.getAdcValue(clusPos)<<endl;
//				cluster.Print(1);
			}
			vecvecSignalLeftLeft.at(det).push_back(signalLeftLeft);
			vecvecSignalRightRight.at(det).push_back(signalRightRight);;
			vecvecLeftEtaSignal.at(det).push_back(leftEtaSignal);;
			vecvecRightEtaSignal.at(det).push_back(rightEtaSignal);;
			vecvecSignalLeftOfHighest.at(det).push_back(signalLeftOfHighest);
			vecvecSignalRightOfHighest.at(det).push_back(signalRightOfHighest);
			vecvecSignalHighest.at(det).push_back(highestSignal);
			vecvecEta.at(det).push_back(eta);
		}
	}
}


void TAnalysisOfClustering::savePHHistos()
{
	vector<Float_t> vecClusterSize,vecMVP,vecClusterSizeError,vecWidth;

	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		histSaver->SaveHistogram(hPHDistribution[det]);
		vecClusterSize.clear();
		vecMVP.clear();
		vecClusterSizeError.clear();
		vecWidth.clear();
		for(UInt_t nClusters=0;nClusters<10;nClusters++){
			stringstream histName;
			//    		string name = (string)hPHDistribution[det]->GetTitle();
			//    		name = name.substr(0,histName.str().size()-6);
			histName<<"hPulseHeightDistribution_";
			if(nClusters==0)
				histName<<"allClusterSizes";
			else
				histName<<"nClusters"<<nClusters;
			histName<<"_"<<TPlaneProperties::getStringForDetector(det);
			TObject *htemp2 = (TObject*)gROOT->FindObject(histName.str().c_str());
			if(htemp2!=0)delete htemp2;

			//CREATE HTEMP and ReBin it if necessary
			TH1F *htemp;
			htemp = (TH1F*)hPHDistribution[det]->ProjectionX(histName.str().c_str(),nClusters+1,nClusters+1);
			if(htemp==0) continue;

			//adjust binning if necessary
			UInt_t entries = htemp->GetEntries();
			UInt_t maximumEntries = htemp->GetMaximum();
			UInt_t nSteps =4;
			UInt_t nStep = 0;
			while((maximumEntries<50&&maximumEntries<entries*0.8)&&nStep<nSteps){
				htemp->Rebin(2);
				entries = htemp->GetEntries();
				maximumEntries = htemp->GetMaximum();
				nStep++;
			}
			TF1 *fit=0;
			htemp->SetTitle(htemp->GetName());
			htemp->GetXaxis()->SetTitle("Charge in ADC units");
			htemp->GetYaxis()->SetTitle("number of entries#");
			LandauGaussFit landauGauss;
			if(nClusters<4||det==TPlaneProperties::getDetDiamond())
				fit = landauGauss.doLandauGaussFit(htemp,nClusters==1&&det==TPlaneProperties::getDetDiamond());
			if(fit!=0){
				if (verbosity > 2) {
					cout<<"Width(scale): "<<fit->GetParameter(0)<<endl;
					cout<<"MostProb:     "<<fit->GetParameter(1)<<endl;
					cout<<"Area:         "<<fit->GetParameter(2)<<endl;
					cout<<"Width(sigma): "<<fit->GetParameter(3)<<endl;
				}
				if(nClusters==0)
					vecPHMeans.push_back(fit->GetParameter(1));
				vecClusterSize.push_back(nClusters);
				vecMVP.push_back(fit->GetParameter(1));
				vecClusterSizeError.push_back(0.5);
				vecWidth.push_back(fit->GetParameter(0));
				histSaver->SaveHistogramLandau(htemp);
			}
			else
				histSaver->SaveHistogram(htemp);
			delete htemp;
		}
		if(det==TPlaneProperties::getDetDiamond()){
			stringstream histTitle;
			histTitle<<"gChargeOfClusterVsClusterSize_"<<det;
			TGraphErrors graph = histSaver->CreateErrorGraph(histTitle.str(),vecClusterSize,vecMVP,vecClusterSizeError,vecWidth);
			graph.GetXaxis()->SetTitle("Cluster Size");
			graph.GetYaxis()->SetTitle("Charge of Cluster");
			histSaver->SaveGraph(&graph,histTitle.str());
		}
		delete hPHDistribution[det];
	}
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		histSaver->SaveHistogram(hBiggestHitVsClusterSize[det]);
		TProfile *profY = hBiggestHitVsClusterSize[det]->ProfileY();
		profY->GetXaxis()->SetTitle("ClusterSize");
		profY->GetYaxis()->SetTitle("mean of Biggest signal in Cl");
		string histoTitle = "mean of Biggest signal in Cluster vs. ClusterSize";
		histoTitle.append(TPlaneProperties::getStringForDetector(det).c_str());
		profY->SetTitle(histoTitle.c_str());
		histSaver->SaveHistogram(profY);
	}
}

void TAnalysisOfClustering::createPHDistribution(){
	bool isValid=true;
	for(UInt_t det =0;det<TPlaneProperties::getNDetectors();det++)
		isValid = (eventReader->getNClusters(det)==1)&&isValid;
	for(UInt_t det =0;det<TPlaneProperties::getNDetectors();det++){
		if(!isValid)
			continue;
		UInt_t nClusterSize = eventReader->getClusterSize(det,0);
		Float_t charge = eventReader->getCluster(det,0).getCharge(true);

		hPHDistribution[det]->Fill(charge,0);
		hPHDistribution[det]->Fill(charge,nClusterSize);
		//		cout<<"Fill PH histo with "<<charge<<" and Clustersize "<<nClusterSize<<endl;
	}
}




