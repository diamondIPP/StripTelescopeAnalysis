//
//  TAnalysisOfPedestal.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 30.11.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TAnalysisOfPedestal.hh"

TAnalysisOfPedestal::TAnalysisOfPedestal(TSettings* settings) {
	cout<<"**********************************************************"<<endl;
	cout<<"*********TAnalysisOfPedestal::TAnalysisOfPedestal*****"<<endl;
	cout<<"**********************************************************"<<endl;
	if(settings==0)
		exit(0);
	res=0;
	htmlPedestal= new THTMLPedestal(settings);
	sys = gSystem;
	UInt_t runNumber=settings->getRunNumber();
	sys->MakeDirectory(settings->getAbsoluteOuputPath().c_str());
	sys->cd(settings->getAbsoluteOuputPath().c_str());
	settings->goToPedestalTreeDir();

	eventReader=new TADCEventReader(settings->getPedestalTreeFilePath(),settings->getRunNumber());
	histSaver=new HistogrammSaver();
	histSaver->SetOptStat("ormen");
	histSaver->SetOptFit(111);
	settings->goToPedestalAnalysisDir();
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	htmlPedestal->setFileGeneratingPath(sys->pwd());
	settings->goToPedestalTreeDir();
	initialiseHistos();
	this->settings=settings;
	cout<<"end initialise"<<endl;

	pedestalMeanValue.resize(TPlaneProperties::getNDetectors());
	pedestalSigmaValue.resize(TPlaneProperties::getNDetectors());
	nPedestalHits.resize(TPlaneProperties::getNDetectors());
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		pedestalMeanValue.at(det).resize(TPlaneProperties::getNChannels(det),0);
		pedestalSigmaValue.at(det).resize(TPlaneProperties::getNChannels(det),0);
		nPedestalHits.at(det).resize(TPlaneProperties::getNChannels(det),0);
	}
	this->diaRawADCvalues.resize(TPlaneProperties::getNChannelsDiamond(),std::vector<UInt_t>());

	verbosity = 0;
}

TAnalysisOfPedestal::~TAnalysisOfPedestal() {
	htmlPedestal->createPageContent();
	//	htmlPedestal->createPedestalDistribution();
	htmlPedestal->generateHTMLFile();
	cout<<"Del htmlPed: "<<flush;
	delete htmlPedestal;
	cout<<"Del eventReader: "<<flush;
	delete eventReader;
	cout<<"Del histSaver: "<<flush;
	delete histSaver;
	settings->goToOutputDir();
}


void TAnalysisOfPedestal::doAnalysis(UInt_t nEvents)
{
	cout<<"TAnalysisOfPedestal::doAnalysis\nanalyze pedestal data..."<<endl;
	//	eventReader->checkADC();
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		eventReader->LoadEvent(nEvent);
		analyseEvent();

	}
	//	createPedestalMeanHistos();
	saveHistos();
}

void TAnalysisOfPedestal::analyseEvent(){
	analyseBiggestHit();
	analyseBiggestHit(true);
	checkForSaturatedChannels();
	updateMeanCalulation();
	//		UInt_t ch=30;
	//		cout<<nEvent<<"  "<<eventReader->getPedestalMean(8,ch,false)<<" "<<eventReader->getPedestalMean(8,ch,true)<<"\t";
	//		cout<<eventReader->getPedestalSigma(8,ch,false)<<" "<<eventReader->getPedestalSigma(8,ch,true)<<"\t"<<eventReader->getCMNoise()<<"\t\t";
	//		cout<<eventReader->getRawSignal(8,ch,false)<<" "<<eventReader->getRawSignal(8,ch,true)<<endl;
	/*cout<<nEvent;
	 for(unsigned int det=0;det< (eventReader->getCluster()->size());det++)
	 for(unsigned int cl=0;cl< eventReader->getCluster()->at(det).size();cl++)
	 cout<<" "<<eventReader->getCluster()->at(det).at(cl).getChargeWeightedMean()<<flush;
	 cout<<endl;//*/
	//		checkForDeadChannels();

	//		getBiggestHit();
	//		analyseForSeeds();
	//		analyseCluster();

}

void TAnalysisOfPedestal::checkForDeadChannels()
{
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		int numberOfSeeds=0;
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			if(settings->isDet_channel_screened(det,ch))
				continue;
			Float_t sigma=eventReader->getPedestalSigma(det,ch);
			Float_t signalInSigma = eventReader->getSignalInSigma(det,ch);
			if(sigma==0){
				//cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
				continue;
			};
			if(signalInSigma>settings->getClusterSeedFactor(det)){
				hSeedMap[det]->Fill(ch);
				numberOfSeeds++;
			}
		}
		hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}

}
void TAnalysisOfPedestal::analyseForSeeds(){
	//	for(int det=0;det<TPlaneProperties::getNDetectors();det++){
	//		int nClusters = eventReader->getNClusters(det);
	//		if(nClusters==1)
	//			hSeedMap2[det]->Fill(eventReader->getCluster(det,0).getHighestSignalChannel());
	//	}
}

void TAnalysisOfPedestal::checkForSaturatedChannels()
{

}
void TAnalysisOfPedestal::getBiggestHit(){
//
//	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
//		int biggestHit=0;
//		Float_t biggestHitInSigma=0;
//		int chBiggestHit;
//		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){//todo warum 70 bis 200
//			if(settings->isDet_channel_screened(det,ch))
//				continue;
//			Int_t adcValue=eventReader->getAdcValue(det,ch);
//			if (adcValue<TPlaneProperties::getMaxSignalHeight(det))
//				if(adcValue>biggestHit){
//					biggestHit=eventReader->getAdcValue(det,ch);
//					biggestHitInSigma=eventReader->getSignalInSigma(det,ch);
//					chBiggestHit=ch;
//				}
//		}
//		if(biggestHitInSigma>4){//removing all events where biggest hit is  <4 sigma, since probability of haveing one channel
//			//higher than 3 sigma is approx 50% (1-0.997^128)
//			hPulsHeightBiggestHit[det]->Fill(biggestHitInSigma);
//			hChannelBiggestHit[det]->Fill(chBiggestHit);
//			if(eventReader->getAdcValue(det,chBiggestHit-1)<eventReader->getAdcValue(det,chBiggestHit+1))
//				hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit+1));
//			else
//				hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit-1));
//		}
//	}

}

void TAnalysisOfPedestal::findPlotRangeForPHHisto(TH1F *histo, Float_t hitCut)
{

	Float_t max = 0;
	Int_t startBin = histo->FindBin(hitCut);
	for(Int_t binX=1;binX<startBin;binX++)
			if ( max < histo->GetBinContent(binX) )
				max = histo->GetBinContent(binX);
	Float_t oldMax = max;
	max=0;
	for(Int_t binX=startBin;binX<histo->GetNbinsX();binX++)
		if ( max < histo->GetBinContent(binX) )
			max = histo->GetBinContent(binX);
	max=max*1.1;
	if (oldMax>2*max){
		histo->SetFillStyle(3244);
		histo->SetFillColor(kGray);//Gray
	}
	else
		if(max<oldMax)
			max=oldMax;

	histo->GetYaxis()->SetRangeUser(0,max);

}

/**
 * create vector with biggest and biggest adjacent hit with PH in Sigma
 */
void TAnalysisOfPedestal::analyseBiggestHit(bool CMN_corrected) {
	//loop over all detectors
	for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++) {
		Float_t biggestSignal = -999;
		UInt_t biggestHitChannel =-999;
		//loop over all channels
		for (UInt_t ch = 0; ch < TPlaneProperties::getNChannels(det); ch++) {
			if(!settings->isDet_channel_screened(det,ch)){
				Float_t signal = eventReader->getSignal(det,ch, CMN_corrected);
				if (!eventReader->isSaturated(det,ch)){
					if (signal > biggestSignal) {
						biggestHitChannel = ch;
						biggestSignal = signal;
					}
				}
				else{
					if(!CMN_corrected){
						vecSaturatedChannels[det].push_back(ch);
						hSaturatedChannels[det]->Fill(ch);
					}
				}
			}
		}//end loop over channels
		Float_t leftSignal=eventReader->getSignal(det,biggestHitChannel-1, CMN_corrected);
		Float_t rightSignal=eventReader->getSignal(det,biggestHitChannel+1,CMN_corrected);
		Float_t biggestAdjacentSignal=0;
		Int_t biggestAdjacentHitChannel= -9999;
		if(leftSignal>rightSignal){
			biggestAdjacentSignal=leftSignal;
			biggestAdjacentHitChannel = biggestHitChannel-1;
		}
		else if (rightSignal>leftSignal){
			biggestAdjacentSignal = rightSignal;
			biggestAdjacentHitChannel = biggestHitChannel+1;
		}
		Int_t hitOrder = biggestAdjacentHitChannel!=-9999?biggestAdjacentHitChannel-biggestHitChannel:0;
//		if(det==6)cout<<nEvent<<" "<<setw(3)<<biggestHitChannel<<": "<<setw(8)<<biggestSignal<<" \t"<<setw(5)<<biggestAdjacentHitChannel<<": "<<setw(8)<<biggestAdjacentSignal<<"\t"<<setw(2)<<hitOrder<<endl;
		if(biggestAdjacentHitChannel!=-9999&&(biggestAdjacentHitChannel<0||biggestAdjacentHitChannel>=(Int_t)TPlaneProperties::getNChannels(det)))
			cout<<"something is wrong: biggestAdjacentHitChannel: "<<biggestAdjacentHitChannel<<" BiggestHitChannel:"<<biggestHitChannel<<" "<<leftSignal<<" "<<rightSignal<<endl;
		Float_t biggestSignalInSigma = eventReader->getSignalInSigma(det,biggestHitChannel,CMN_corrected);
		Float_t biggestAdjacentSignalInSigma =biggestAdjacentHitChannel>=0&&biggestAdjacentHitChannel< (Int_t)TPlaneProperties::getNChannels(det)? eventReader->getSignalInSigma(det,biggestAdjacentHitChannel,CMN_corrected):0;
		if(!CMN_corrected){
			vecBiggestHitChannel[det].push_back(biggestHitChannel);
			vecBiggestSignalInSigma[det].push_back(biggestSignalInSigma);
			if(biggestAdjacentHitChannel!=-9999){
				vecBiggestAdjacentSignal[det].push_back(biggestAdjacentSignalInSigma);
				vecBiggestAdjacentHitChannel[det].push_back(biggestAdjacentHitChannel);
				vecHitOrder[det].push_back(biggestAdjacentHitChannel-biggestHitChannel);
				hBiggestAdjacentSignalInSigma[det]->Fill(biggestAdjacentSignalInSigma);
			}
			else{
				vecHitOrder[det].push_back(0);
				vecBiggestAdjacentSignal[det].push_back(0);
				vecBiggestAdjacentHitChannel[det].push_back(-9999);
			}
			hHitOrderMap[det]->Fill(hitOrder);
			hBiggestHitChannelMap[det]->Fill(biggestHitChannel);
		}
		else if (TPlaneProperties::isDiamondDetector(det)){
			vecBiggestHitChannelCMN[det].push_back(biggestHitChannel);
			vecBiggestSignalInSigmaCMN[det].push_back(biggestSignalInSigma);
			if(biggestAdjacentHitChannel!=-9999){
				vecBiggestAdjacentSignalCMN[det].push_back(biggestAdjacentSignalInSigma);
				vecBiggestAdjacentHitChannelCMN[det].push_back(biggestAdjacentHitChannel);
//				vecHitOrder[det].push_back(biggestAdjacentHitChannel-biggestHitChannel);
				hBiggestAdjacentSignalInSigmaCMN[det]->Fill(biggestAdjacentSignalInSigma);
			}
			else{
//				vecHitOrder[det].push_back(0);
				vecBiggestAdjacentSignalCMN[det].push_back(0);
				vecBiggestAdjacentHitChannelCMN[det].push_back(-9999);
			}
//			hHitOrderMap[det]->Fill(hitOrder);
			hBiggestHitChannelMapCMN[det]->Fill(biggestHitChannel);
		}
	}//end loop over detectors
}


void TAnalysisOfPedestal::initialiseHistos()
{
	hCMNoiseDistribution= new TH1F("hCMNoiseDistribution","hCMNoiseDistribution",512,-20,20);
	hCMNoiseDistribution->GetXaxis()->SetTitle("Common Mode Noise [ADC]");
	hCMNoiseDistribution->GetYaxis()->SetTitle("number of entries [#]");
	for (UInt_t det =0;det<TPlaneProperties::getNDetectors();det++){
		stringstream histoName,histoTitle,xTitle,yTitle;
		histoName<<"hNoiseDistributionOfAllNonHitChannels_"<<TPlaneProperties::getStringForDetector(det);
		histoTitle<<"Noise Distribution  of all non hit channels in Plane"<<TPlaneProperties::getStringForDetector(det);
		xTitle<<"non hit Noise (Adc-Ped.) in ADC counts";
		yTitle<<"Number of Entries #";
		Float_t width = 8;
		UInt_t nBins =64;
		if (det==TPlaneProperties::getDetDiamond()){
			width = 32;
			nBins=128;
		}
		if(TPlaneProperties::isSiliconDetector(det)){
			hAllAdcNoise[det]= new TH1F(histoName.str().c_str(),histoTitle.str().c_str(),nBins,(-1)*width,width);
			hAllAdcNoise[det]->GetXaxis()->SetTitle(xTitle.str().c_str());
			hAllAdcNoise[det]->GetYaxis()->SetTitle(yTitle.str().c_str());
		}
		else{
			hDiaAllAdcNoise= new TH1F(histoName.str().c_str(),histoTitle.str().c_str(),nBins,(-1)*width,width);
			hDiaAllAdcNoise->GetXaxis()->SetTitle(xTitle.str().c_str());
			hDiaAllAdcNoise->GetYaxis()->SetTitle(yTitle.str().c_str());
			histoName<<"_CMNcorrected";
			histoTitle<<" - CMN corrected";
			hDiaAllAdcNoiseCMN = new TH1F(histoName.str().c_str(),histoTitle.str().c_str(),nBins,(-1)*width,width);
			xTitle<<" -  CMN corrected";
			hDiaAllAdcNoiseCMN->GetXaxis()->SetTitle(xTitle.str().c_str());
			hDiaAllAdcNoiseCMN->GetYaxis()->SetTitle(yTitle.str().c_str());
		}
	}

	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSaturatedChannels_"<<TPlaneProperties::getStringForDetector(det)<<"";
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSeedPosAllSeeds_"<<TPlaneProperties::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hMaxSeedPos_"<<TPlaneProperties::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hNumberOfSeeds_"<<TPlaneProperties::getStringForDetector(det);
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPulseHeight__BiggestHitChannelInSigma_"<<TPlaneProperties::getStringForDetector(det);
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){//todo why such a big histo?so big?
		stringstream histoName;
		histoName<<"hPulseHeight_BiggestHitNextToBiggestHit_ChannelInSigma"<<TPlaneProperties::getStringForDetector(det);
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hChannel_BiggestHit_"<<TPlaneProperties::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels);
	}

	for (UInt_t det = 0; det < 9; det++) {
		int nbins = 256;
		Float_t min = 0.;
		Float_t max = 64.;
		if(det==TPlaneProperties::getDetDiamond()){max=128;nbins=512;}

		stringstream histoName;


		histoName.str("");
		histoName << "hPulseHeight_SecondBiggestHitChannelInSigma_" << TPlaneProperties::getStringForDetector(det);
		hBiggestAdjacentSignalInSigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
		histoName << "_CMNcorrected";
		hBiggestAdjacentSignalInSigmaCMN[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hSecondBiggestHitMinusBiggestHitPosition_" << TPlaneProperties::getStringForDetector(det);
		hHitOrderMap[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),3,-1.5,1.5);

		histoName.str("");
		histoName << "hPulseHeightSecondBiggestHitChannelInSigmaLeft" << TPlaneProperties::getStringForDetector(det);
		histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hPulseHeightSecondBiggestHitChannelInSigmaRight" << TPlaneProperties::getStringForDetector(det) ;
		histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hBiggestHitMap"<< TPlaneProperties::getStringForDetector(det);
		hBiggestHitChannelMap[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),TPlaneProperties::getNChannels(det),0.,TPlaneProperties::getNChannels(det)-1);
		histoName << "_CMNcorrected";
		hBiggestHitChannelMapCMN[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),TPlaneProperties::getNChannels(det),0.,TPlaneProperties::getNChannels(det)-1);

		histoName.str("");
		histoName << "hPulseHeightLeftChipBiggestHitChannelInSigma" << TPlaneProperties::getStringForDetector(det);
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

}

void TAnalysisOfPedestal::savePHinSigmaHistos(){
	for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++) {
		TString histoName =TString::Format("hPulseHeight_BiggestHitChannelInSigma%s",TPlaneProperties::getStringForDetector(det).c_str()) ;
		double cut = settings->getClusterSeedFactor(det);
		cout << "saving histogram " << histoName<< ".. with CUT on " <<cut<< endl;
		Float_t max=0;
		Float_t mean=0;
		Float_t sigma=0;
		//Find max mean and sigma
		for(UInt_t i=0;i<vecBiggestSignalInSigma[det].size();i++){
			mean+=vecBiggestSignalInSigma[det][i];
			sigma+=vecBiggestSignalInSigma[det][i]*vecBiggestSignalInSigma[det][i];
			if (vecBiggestSignalInSigma[det][i]>max)
				max = vecBiggestSignalInSigma[det][i];
		}
		mean/=(Float_t)vecBiggestSignalInSigma[det].size();
		sigma=TMath::Sqrt(sigma/(Float_t)vecBiggestSignalInSigma[det].size()-mean*mean);
		cout<< "Mean: "<<mean<<" +/- "<<sigma<<"\tMaximum SNR: "<<max<<endl;
		//define xrange and nbins
		Float_t xRangeMax = TMath::Min(mean+3*sigma,max);
		Float_t xRangeMin = 0;
		UInt_t nbins = 512;
		cout<<"X Range: "<<xRangeMin<<","<<xRangeMax<<endl;
		hBiggestSignalInSigma[det] = new TH1F(histoName,histoName,nbins,xRangeMin,xRangeMax);
		histoName = histoName+"_CMNcorrected";
		hBiggestSignalInSigmaCMN[det] = new TH1F(histoName,histoName,nbins,xRangeMin,xRangeMax);
		//set Axis Titles
		this->hBiggestSignalInSigma[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->hBiggestSignalInSigma[det]->GetYaxis()->SetTitle("number of entries #");
		this->hBiggestSignalInSigmaCMN[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->hBiggestSignalInSigmaCMN[det]->GetYaxis()->SetTitle("number of entries #");
		histoName =TString::Format("hPulseHeight_BiggestHitChannelInSigma%s_2D",TPlaneProperties::getStringForDetector(det).c_str());
		Int_t nChannels = TPlaneProperties::getNChannels(det);
		TH2F* hBiggestSignalInSigma2D = new TH2F(histoName,histoName,1000,0,500,nChannels,0,nChannels-1);

		histoName =TString::Format("hPulseHeight_BiggestHitChannelInSigma%s_CMNcorrected_2D",TPlaneProperties::getStringForDetector(det).c_str()) ;
		TH2F* hBiggestSignalInSigma2DCMN = new TH2F(histoName,histoName,1000,0,500,nChannels,0,nChannels-1);
		//set Axis Titles
		hBiggestSignalInSigma2D->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		hBiggestSignalInSigma2D->GetYaxis()->SetTitle("Channel No of Biggest Hit");
		hBiggestSignalInSigma2D->GetZaxis()->SetTitle("number of entries #");

		hBiggestSignalInSigma2DCMN->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		hBiggestSignalInSigma2DCMN->GetYaxis()->SetTitle("Channel No of Biggest Hit");
		hBiggestSignalInSigma2DCMN->GetZaxis()->SetTitle("number of entries #");
		for(UInt_t i=0;i<vecBiggestSignalInSigma[det].size();i++){
			Float_t signal=vecBiggestSignalInSigma[det].at(i);
			Int_t channel = vecBiggestHitChannel[det].at(i);
			this->hBiggestSignalInSigma[det]->Fill(signal);
			hBiggestSignalInSigma2D->Fill(signal,channel);
		}
		for(UInt_t i=0;i<vecBiggestSignalInSigmaCMN[det].size();i++){
			Float_t signal=vecBiggestSignalInSigmaCMN[det].at(i);
			Int_t channel = vecBiggestHitChannelCMN[det].at(i);
			this->hBiggestSignalInSigmaCMN[det]->Fill(signal);
			hBiggestSignalInSigma2DCMN->Fill(signal,channel);
		}
		findPlotRangeForPHHisto(hBiggestSignalInSigma[det],settings->getClusterHitFactor(det));
		findPlotRangeForPHHisto(hBiggestSignalInSigmaCMN[det],settings->getClusterHitFactor(det));
		hBiggestSignalInSigma2D->GetXaxis()->SetRangeUser(hBiggestSignalInSigma[det]->GetXaxis()->GetXmin(),hBiggestSignalInSigma[det]->GetXaxis()->GetXmax());
		hBiggestSignalInSigma2DCMN->GetXaxis()->SetRangeUser(hBiggestSignalInSigma[det]->GetXaxis()->GetXmin(),hBiggestSignalInSigma[det]->GetXaxis()->GetXmax());
		hBiggestSignalInSigmaCMN[det]->GetXaxis()->SetRangeUser(hBiggestSignalInSigma[det]->GetXaxis()->GetXmin(),hBiggestSignalInSigma[det]->GetXaxis()->GetXmax());
		histSaver->SaveHistogram(hBiggestSignalInSigma2D);
		cout<<"BiggestHitSNR:\t"<<det<<"\t"<<hBiggestSignalInSigmaCMN[det]<<"\t"<<hBiggestSignalInSigmaCMN[det]->GetEntries()<<"\t"<<hBiggestSignalInSigmaCMN[det]->GetName()<<endl;
		histSaver->SaveHistogram(hBiggestSignalInSigmaCMN[det]);
		histSaver->SaveHistogram(hBiggestSignalInSigma2DCMN);

		TCanvas *c1 = new TCanvas(this->hBiggestSignalInSigma[det]->GetTitle(),this->hBiggestSignalInSigma[det]->GetTitle());
		c1->cd();
		this->hBiggestSignalInSigma[det]->Draw();
		double xCor[] = {cut,cut};
		double yCor[] = {0,this->hBiggestSignalInSigma[det]->GetMaximum()*2};
		TGraph* lineGraph = new TGraph(2,xCor,yCor);
		lineGraph->SetLineColor(kRed);
		lineGraph->SetLineWidth(2);
		lineGraph->Draw("Lsame");
		histSaver->SaveCanvas(c1);;
		//        histSaver->SaveHistogram(this->histo_pulseheight_sigma[det]);
		delete hBiggestSignalInSigma[det];
		delete lineGraph;
		delete c1;
		if (TPlaneProperties::isDiamondDetector(det)){
			for(Int_t i = 0;i < settings->getNDiaDetectorAreas();i++){
				std::pair<Int_t, Int_t> area = settings->getDiaDetectorArea(i);
				if(area.first<area.second){
					TString name = histoName;
					name.Append(TString::Format("_area%d_ch%d-ch%d",i,area.first,area.second));
					cout<<"Create "<<name<<endl;
					TH1F* hBiggestSignalArea = (TH1F*) hBiggestSignalInSigma2D->ProjectionX(name,area.first,area.second);
					findPlotRangeForPHHisto(hBiggestSignalArea,settings->getClusterSeedFactor(det));
					histSaver->SaveHistogram(hBiggestSignalArea);
					delete hBiggestSignalArea;
					name = histoName;
					name.Append(TString::Format("CMNcorrected_area%d_ch%d-ch%d",i,area.first,area.second));
					cout<<"Create "<<name<<endl;
					TH1F* hBiggestSignalAreaCMN = (TH1F*) hBiggestSignalInSigma2DCMN->ProjectionX(name,area.first,area.second);
					findPlotRangeForPHHisto(hBiggestSignalAreaCMN,settings->getClusterSeedFactor(det));
					histSaver->SaveHistogram(hBiggestSignalAreaCMN);
					delete hBiggestSignalAreaCMN;
				}
			}

		}
		delete hBiggestSignalInSigma2D;
	}
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){

		double cut = settings->getClusterHitFactor(det);
		//		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".. with CUT on " <<cut<< endl;
		TCanvas *c1 = new TCanvas(this->hBiggestAdjacentSignalInSigma[det]->GetTitle(),this->hBiggestAdjacentSignalInSigma[det]->GetTitle());
		c1->cd();
		this->hBiggestAdjacentSignalInSigma[det]->Draw();
		double xCor[] = {cut,cut};
		double yCor[] = {0,this->hBiggestAdjacentSignalInSigma[det]->GetMaximum()*2};
		this->hBiggestAdjacentSignalInSigma[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->hBiggestAdjacentSignalInSigma[det]->GetXaxis()->SetTitle("number of entries #");
		findPlotRangeForPHHisto(hBiggestAdjacentSignalInSigma[det],settings->getClusterHitFactor(det));
		TGraph* lineGraph = new TGraph(2,xCor,yCor);
		lineGraph->SetLineColor(kRed);
		lineGraph->SetLineWidth(2);
		lineGraph->Draw("Lsame");
		histSaver->SaveCanvas(c1);;
		delete hBiggestAdjacentSignalInSigma[det];
		delete lineGraph;
		delete c1;
		if(TPlaneProperties::isDiamondDetector(det)){
			double cut = settings->getClusterHitFactor(det);
			//		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".. with CUT on " <<cut<< endl;
			TCanvas *c2 = new TCanvas(this->hBiggestAdjacentSignalInSigmaCMN[det]->GetTitle(),this->hBiggestAdjacentSignalInSigmaCMN[det]->GetTitle());
			c2->cd();
			this->hBiggestAdjacentSignalInSigmaCMN[det]->Draw();
			double xCorCMN[] = {cut,cut};
			double yCorCMN[] = {0,this->hBiggestAdjacentSignalInSigmaCMN[det]->GetMaximum()*2};
			this->hBiggestAdjacentSignalInSigmaCMN[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
			this->hBiggestAdjacentSignalInSigmaCMN[det]->GetXaxis()->SetTitle("number of entries #");
			findPlotRangeForPHHisto(hBiggestAdjacentSignalInSigmaCMN[det],settings->getClusterHitFactor(det));
			TGraph* lineGraphCMN = new TGraph(2,xCorCMN,yCorCMN);
			lineGraphCMN->SetLineColor(kRed);
			lineGraphCMN->SetLineWidth(2);
			lineGraphCMN->Draw("Lsame");
			histSaver->SaveCanvas(c2);;
			delete hBiggestAdjacentSignalInSigmaCMN[det];
			delete lineGraphCMN;
			delete c2;
		}
	}
}


void TAnalysisOfPedestal::saveHistos(){
	createPedestalMeanHistos();
	savePHinSigmaHistos();
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSeedMap[det]);
		hSeedMap[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSeedMap2[det]);
		hSeedMap2[det]->Delete();
	}
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hNumberOfSeeds[det]);
		hNumberOfSeeds[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(hPulsHeightBiggestHit[det]);
		hPulsHeightBiggestHit[det]->Delete();
	}

	/**********************************************************************************
	 *
	 **********************************************************************************/
	TGraph *gAvrgPedCMN,*gAvrgPed,*gCMNoise,*gAvrgNoise,*gAvrgNoiseCMN;
	if(vecEventNo.size()==vecAvrgPed.size()&&vecEventNo.size()>0){
		if(verbosity>1)cout<<"Creating gAvrgPed with "<<vecEventNo.size()<<" Entries!"<<endl;
		gAvrgPed = new TGraph(vecEventNo.size(),&vecEventNo.at(0),&vecAvrgPed.at(0));
		gAvrgPed->SetTitle("Averg Pedestal of Diamond vs. eventNo");
		gAvrgPed->SetName("gAvrgPed");
		gAvrgPed->Draw("AL");
		gAvrgPed->GetXaxis()->SetTitle("EventNo.");
		gAvrgPed->GetYaxis()->SetTitle("Avrg. Ped Value [ADC counts]");
		histSaver->SaveGraph(gAvrgPed,"gAvrgPed","AL");
	}
	else{
		cout<<"Size for creatign gAvrgPed are wrong: "<<vecEventNo.size()<<"/"<<vecAvrgPed.size()<<endl;
	}
	if(vecEventNo.size()==vecAvrgPedCMN.size()){
		if(verbosity>1)cout<<"Creating gAvrgPedCMN with "<<vecEventNo.size()<<" Entries!"<<endl;
		gAvrgPedCMN = new TGraph(vecEventNo.size(),&vecEventNo.at(0),&vecAvrgPedCMN.at(0));
		gAvrgPedCMN->SetTitle("Averg CMN corrected Pedestal of Diamond vs. eventNo");
		gAvrgPedCMN->SetName("gAvrgPedCMN");
		gAvrgPedCMN->Draw("AL");
		gAvrgPedCMN->GetXaxis()->SetTitle("EventNo.");
		gAvrgPedCMN->GetYaxis()->SetTitle("Avrg. CMN corrected Ped Value [ADC counts]");
		histSaver->SaveGraph(gAvrgPedCMN,"gAvrgPedCMN","AL");
	}
	else{
		cout<<"Size for creatign gAvrgPedCMN are wrong: "<<vecEventNo.size()<<"/"<<vecAvrgPedCMN.size()<<endl;
	}
	if(vecEventNo.size()==vecAvrgSigma.size()){
		if(verbosity>1)cout<<"Creating gAvrgNoise with "<<vecEventNo.size()<<" Entries!"<<endl;
		gAvrgNoise = new TGraph(vecEventNo.size(),&vecEventNo.at(0),&vecAvrgSigma.at(0));
		gAvrgNoise->SetTitle("Avrg. Noise of Diamond vs. eventNo");
		gAvrgNoise->SetName("gAvrgNoise");
		gAvrgNoise->Draw("AL");
		gAvrgNoise->GetXaxis()->SetTitle("EventNo.");
		gAvrgNoise->GetYaxis()->SetTitle("Avrg Noise of diamond [ADC counts]");
		histSaver->SaveGraph(gAvrgNoise,"gAvrgNoise","AL");
	}
	else{
		cout<<"Size for creatign gAvrgNoise are wrong: "<<vecEventNo.size()<<"/"<<vecCMNoise.size()<<endl;
	}
	if(vecEventNo.size()==vecAvrgSigmaCMN.size()){
		if(verbosity>1)cout<<"Creating gAvrgNoiseCMN with "<<vecEventNo.size()<<" Entries!"<<endl;
		gAvrgNoiseCMN = new TGraph(vecEventNo.size(),&vecEventNo.at(0),&vecAvrgSigmaCMN.at(0));
		gAvrgNoiseCMN->SetTitle("Avrg. CMN correctedNoise of Diamond vs. eventNo");
		gAvrgNoiseCMN->SetName("gAvrgNoiseCMN");
		gAvrgNoiseCMN->Draw("AL");
		gAvrgNoiseCMN->GetXaxis()->SetTitle("EventNo.");
		gAvrgNoiseCMN->GetYaxis()->SetTitle("Avrg CMN corrected Noise of diamond [ADC counts]");
		histSaver->SaveGraph(gAvrgNoiseCMN,"gAvrgNoiseCMN","AL");
	}
	else{
		cout<<"Size for creatign gAvrgNoiseCMN are wrong: "<<vecEventNo.size()<<"/"<<vecCMNoise.size()<<endl;
	}
	if(vecEventNo.size()==vecAvrgPedCMN.size()){
		cout<<"Creating gAvrgPedCMN with "<<vecEventNo.size()<<" Entries!"<<endl;
		gCMNoise = new TGraph(vecEventNo.size(),&vecEventNo.at(0),&vecCMNoise.at(0));
		gCMNoise->SetTitle("Common Mode Noise vs. eventNo");
		gCMNoise->SetName("gCMNoise");
		gCMNoise->Draw("AL");
		gCMNoise->GetXaxis()->SetTitle("EventNo.");
		gCMNoise->GetYaxis()->SetTitle("Common Mode Noise[ADC counts]");
		histSaver->SaveGraph(gCMNoise,"gCMNoise","AL");
	}
	else{
		cout<<"Size for creatign gCMNoise are wrong: "<<vecEventNo.size()<<"/"<<vecCMNoise.size()<<endl;
	}
	//	if(gAvrgNoise!=0&&gAvrgPed!=0){
	//	  TMultiGraph* mg = new TMultiGraph('mgAvrgPedestal','Avrg Pedestal and Noise for diamond');
	//	}
	//	for(int det=0;det<9;det++){
	//		histSaver->SaveHistogramPNG(hPulsHeightNextBiggestHit[det]);
	//		hPulsHeightNextBiggestHit[det]->Delete();
	//	}
	//	for(int det=0;det<9;det++){
	//		histSaver->SaveHistogramPNG(hChannelBiggestHit[det]);
	//		hChannelBiggestHit[det]->Delete();
	//	}
	//	for(int det=0;det<9;det++){
	//		histSaver->SaveHistogramPNG(this->hClusterSize[det]);
	//		histSaver->SaveHistogramPNG(this->hNumberOfClusters[det]);
	//		delete hClusterSize[det];
	//		delete hNumberOfClusters[det];
	//	}


	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
		//		cout << "saving histogram" << this->histo_pulseheight_sigma125[det]->GetName() << ".." << endl;
		//		histSaver->SaveHistogramPNG(this->histo_pulseheight_sigma125[det]);
		//		cout << "saving histogram " << this->histo_second_biggest_hit_direction[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->hHitOrderMap[det]);
		//		cout << "saving histogram " << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->hBiggestHitChannelMap[det]);
		histSaver->SaveHistogram(this->hBiggestHitChannelMapCMN[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma_second[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma[det]);
		//		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma_second[det]);

		//		delete histo_pulseheight_sigma125[det];
		delete hHitOrderMap[det];
		delete hBiggestHitChannelMap[det];
		delete hBiggestHitChannelMapCMN[det];
		delete histo_pulseheight_left_sigma[det];
		delete histo_pulseheight_left_sigma_second[det];
		delete histo_pulseheight_right_sigma[det];
		delete histo_pulseheight_right_sigma_second[det];

		if(TPlaneProperties::isSiliconDetector(det)){
			TF1 histofitx("histofitx","gaus",hAllAdcNoise[det]->GetMean()-2*hAllAdcNoise[det]->GetRMS(),hAllAdcNoise[det]->GetMean()+2*hAllAdcNoise[det]->GetRMS());
			histofitx.SetLineColor(kBlue);
			hAllAdcNoise[det]->Fit(&histofitx,"rq");
			if(res!=0)res->SetNoise(det,histofitx.GetParameter(2));
			histSaver->SaveHistogram(hAllAdcNoise[det],true);

			delete hAllAdcNoise[det];
		}
	}
	cout<<"hFitX:"<<hDiaAllAdcNoise<<endl;
	cout<<hDiaAllAdcNoise->GetMean()<<endl;
	cout<<hDiaAllAdcNoise->GetRMS()<<endl;

	TF1 histofitx("histofitx","gaus",hDiaAllAdcNoise->GetMean()-2*hDiaAllAdcNoise->GetRMS(),hDiaAllAdcNoise->GetMean()+2*hDiaAllAdcNoise->GetRMS());
	histofitx.SetLineColor(kBlue);
	hDiaAllAdcNoise->Fit(&histofitx,"rq");
	if(res!=0)res->SetNoise(TPlaneProperties::getDetDiamond(),histofitx.GetParameter(2));
	histSaver->SaveHistogram(hDiaAllAdcNoise,true);
	TF1 histofitAllNoiseCMN("histofitx","gaus",hDiaAllAdcNoiseCMN->GetMean()-2*hDiaAllAdcNoiseCMN->GetRMS(),hDiaAllAdcNoiseCMN->GetMean()+2*hDiaAllAdcNoiseCMN->GetRMS());
	histofitAllNoiseCMN.SetLineColor(kBlue);
	hDiaAllAdcNoiseCMN->Fit(&histofitAllNoiseCMN,"rq");
	//    if(res!=0)res->SetNoise(TPlaneProperties::getDetDiamond(),histofitx.GetParameter(2));
	histSaver->SaveHistogram(hDiaAllAdcNoiseCMN,true);
	for(Int_t ntries=0;hCMNoiseDistribution->GetBinContent(hCMNoiseDistribution->GetMaximumBin())<.05*hCMNoiseDistribution->GetEntries()&&ntries<2;ntries++)
		hCMNoiseDistribution->Rebin(2);
	histSaver->SaveHistogram(hCMNoiseDistribution,true);
}

/**
 *
 */
void TAnalysisOfPedestal::createPedestalMeanHistos()
{
	for(UInt_t det = 0; det<TPlaneProperties::getNDetectors();det++){
		stringstream nameMean,titleMean,titleSigma,nameSigma,canvasTitle,graphTitle;
		nameMean<<"hMeanPedestal_Value_OfChannel_"<<TPlaneProperties::getStringForDetector(det);
		nameSigma<<"hMeanPedestal_Width_OfChannel_"<<TPlaneProperties::getStringForDetector(det);
		titleMean<<"mean of pedestalValue for each channel of "<<TPlaneProperties::getStringForDetector(det);
		titleSigma<<"mean of pedestalWidth for each channel of "<<TPlaneProperties::getStringForDetector(det);
		UInt_t nBins = pedestalMeanValue.at(det).size();
		TH1F *histoMean = new TH1F(nameMean.str().c_str(),titleMean.str().c_str(),nBins,-.5,nBins-.5);
		TH1F *histoSigma = new TH1F(nameSigma.str().c_str(),titleSigma.str().c_str(),nBins,-.5,nBins-.5);
		histoMean->GetXaxis()->SetTitle("channel No");
		histoMean->GetYaxis()->SetTitle("mean pedestal value");
		histoSigma->GetXaxis()->SetTitle("channel No");
		histoSigma->GetYaxis()->SetTitle("mean pedestal sigma");
		vector<Float_t> vecChNo,vecChError;
		for(UInt_t ch = 0; ch<TPlaneProperties::getNChannels(det);ch++){
			if(nPedestalHits.at(det).at(ch)!=0){
				this->pedestalMeanValue.at(det).at(ch)/=nPedestalHits.at(det).at(ch);
				this->pedestalSigmaValue.at(det).at(ch)/=nPedestalHits.at(det).at(ch);
			}
			else {
				cout<<"No channel non-hits in "<<det<<" "<<ch<<":"<<nPedestalHits.at(det).at(ch)<<endl;
				this->pedestalMeanValue.at(det).at(ch)=0;
				this->pedestalSigmaValue.at(det).at(ch)=0;
			}
			if(this->pedestalMeanValue.at(det).at(ch)!=this->pedestalMeanValue.at(det).at(ch))
				this->pedestalMeanValue.at(det).at(ch)=0;
			if(this->pedestalSigmaValue.at(det).at(ch)!=this->pedestalSigmaValue.at(det).at(ch))
				this->pedestalSigmaValue.at(det).at(ch)=0;
			histoMean->Fill(ch,pedestalMeanValue.at(det).at(ch));
			histoSigma->Fill(ch,pedestalSigmaValue.at(det).at(ch));
			vecChNo.push_back(ch+.1);
			vecChError.push_back(0);
		}
		Float_t max = histoMean->GetMaximum()*1.1;
		histoSigma->SetLineColor(kRed);
		histSaver->SaveHistogram(histoMean);
		histSaver->SaveHistogram(histoSigma);
		canvasTitle<<"cPedestalOfChannels_"<<TPlaneProperties::getStringForDetector(det);
		histSaver->SaveTwoHistos(canvasTitle.str(),histoMean,histoSigma,10);
		TGraphErrors *graph = new TGraphErrors(nBins,&vecChNo.at(0),&pedestalMeanValue.at(det).at(0),&vecChError.at(0),&pedestalSigmaValue.at(det).at(0));
		graph->Draw("APLgoff");
		graph->GetXaxis()->SetTitle("channel No.");
		graph->GetYaxis()->SetTitle("pedestalValue in ADC counts");
		graph->GetYaxis()->SetRangeUser(0,max);
		graph->GetXaxis()->SetRangeUser(0,nBins-1);
		graphTitle<<"gMeanPedestalValueOfChannelWithSigmaAsError_"<<TPlaneProperties::getStringForDetector(det);
		graph->SetTitle(graphTitle.str().c_str());
		histSaver->SaveGraph(graph,graphTitle.str(),"AP");
		delete histoMean;
	}


}

void TAnalysisOfPedestal::updateMeanCalulation(){
	Float_t sumPed =0;
	Float_t sumPedCMN=0;
	Float_t sumNoise=0;
	Float_t sumNoiseCMN=0;
	int nSumPed =0;
	int nSumPedCMN=0;
	int nSumNoiseCMN=0;
	int nSumNoise=0;
	Float_t cmNoise = eventReader->getCMNoise();
	vecCMNoise.push_back(cmNoise);
	vecEventNo.push_back(nEvent);
	hCMNoiseDistribution->Fill(cmNoise);
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			if(settings->isDet_channel_screened(det,ch))
				continue;
			Float_t snr = eventReader->getSignalInSigma(det,ch);
			Float_t pedestal = eventReader->getPedestalMean(det,ch,false);
			Float_t sigma = eventReader->getPedestalSigma(det,ch,false);
			UInt_t adc = eventReader->getAdcValue(det,ch);
			Float_t noise = eventReader->getRawSignal(det,ch,false);//adc-pedestal;
			Float_t pedestalCMN = eventReader->getPedestalMean(det,ch,true);
			////		  Float_t sigmaCMN = eventReader->getPedestalSigma(det,ch,true);
			//		  pedestalCMN-=cmNoise;
			Float_t noiseCMN =  eventReader->getRawSignal(det,ch,true);
			if(snr<settings->getClusterHitFactor(det)){
				pedestalMeanValue.at(det).at(ch) +=pedestal;
				pedestalSigmaValue.at(det).at(ch) +=sigma;
				nPedestalHits.at(det).at(ch)++;
				if(TPlaneProperties::isSiliconDetector(det))
					hAllAdcNoise[det]->Fill(noise);
				else if(TPlaneProperties::isDiamondDetector(det)){

					hDiaAllAdcNoise->Fill(noise);
					hDiaAllAdcNoiseCMN->Fill(noiseCMN);
					sumPed += pedestal;
					sumPedCMN+=pedestalCMN;
					nSumPed++;
					nSumPedCMN++;
					nSumNoise++;
					nSumNoiseCMN++;
					sumNoise+=noise;
					sumNoiseCMN+=noiseCMN;
				}
			}
			if(TPlaneProperties::getDetDiamond()==det){
				diaRawADCvalues.at(ch).push_back(adc);
			}

		}
	vecAvrgPed.push_back(sumPed/(float)nSumPed);
	vecAvrgPedCMN.push_back(sumPedCMN/(float)nSumPedCMN);
	vecAvrgSigma.push_back(sumNoise/(float)nSumNoise);
	vecAvrgSigmaCMN.push_back(sumNoiseCMN/(float)nSumNoiseCMN);
}


