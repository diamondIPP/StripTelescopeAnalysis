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
	//	stringstream  filepath;
	//	filepath.str("");
	//	filepath<<"pedestalData."<<runNumber<<".root";
	//	cout<<"currentPath: "<<sys->pwd()<<endl;
	//	cout<<filepath.str()<<endl;
	settings->goToPedestalTreeDir();

	eventReader=new TADCEventReader(settings->getPedestalTreeFilePath(),settings->getRunNumber());
	histSaver=new HistogrammSaver();
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
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
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
		checkForSaturatedChannels();
		//		getBiggestHit();
		//		analyseForSeeds();
		//		analyseCluster();
		analyseBiggestHit();
		updateMeanCalulation();
	}
	//	createPedestalMeanHistos();
	saveHistos();
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
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			if(settings->isDet_channel_screened(det,ch))
				continue;
			if(eventReader->getAdcValue(det,ch)>=TPlaneProperties::getMaxSignalHeight(det)){
				hSaturatedChannels[det]->Fill(ch);
			}
		}
}
void TAnalysisOfPedestal::getBiggestHit(){

	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		int biggestHit=0;
		Float_t biggestHitInSigma=0;
		int chBiggestHit;
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){//todo warum 70 bis 200
			if(settings->isDet_channel_screened(det,ch))
				continue;
			Int_t adcValue=eventReader->getAdcValue(det,ch);
			if (adcValue<TPlaneProperties::getMaxSignalHeight(det))
				if(adcValue>biggestHit){
					biggestHit=eventReader->getAdcValue(det,ch);
					biggestHitInSigma=eventReader->getSignalInSigma(det,ch);
					chBiggestHit=ch;
				}
		}
		if(biggestHitInSigma>4){//removing all events where biggest hit is  <4 sigma, since probability of haveing one channel
			//higher than 3 sigma is approx 50% (1-0.997^128)
			hPulsHeightBiggestHit[det]->Fill(biggestHitInSigma);
			hChannelBiggestHit[det]->Fill(chBiggestHit);
			if(eventReader->getAdcValue(det,chBiggestHit-1)<eventReader->getAdcValue(det,chBiggestHit+1))
				hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit+1));
			else
				hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit-1));
		}
	}

}

void TAnalysisOfPedestal::findPlotRangeForPHHisto(TH1F *histo, Float_t hitCut)
{

	// find a good Xaxis Range

	//first step find number of events which are smaller than hit factor
	UInt_t nBinsX=histo->GetXaxis()->GetNbins();
	UInt_t binx=0;
	UInt_t nEventsSmallerThanHitCut=0;
	for(binx=0;histo->GetBinCenter(binx)<hitCut&&binx<nBinsX;binx++)
		nEventsSmallerThanHitCut+= histo->GetBinContent(binx);

	//choose range of axis in such a way that the upper 20% are not included in the drawing
	UInt_t nMaxEvents = histo->GetEntries() - nEventsSmallerThanHitCut;
	nMaxEvents *= (1 - settings->getPHinSigmaPlotFactor());

	UInt_t nEventsfromBinX=0;
	for(binx= nBinsX;binx>0&&nEventsfromBinX<nMaxEvents;binx--)
		nEventsfromBinX+= histo->GetBinContent(binx);

	Float_t xHigh = histo->GetXaxis()->GetBinCenter(binx);
	histo->GetXaxis()->SetRangeUser(0,xHigh);

	Float_t max = 0;

	for(Int_t binX= histo->FindBin(hitCut);binX<histo->GetNbinsX();binX++)
		if(max<histo->GetBinContent(binX))max=histo->GetBinContent(binX);

	histo->GetYaxis()->SetRange(0,max*1.1);



}

/**
 * pedestalMeanValue,pedestalSigmaValue
 */
void TAnalysisOfPedestal::analyseBiggestHit() {
	for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++) {//todo change to 9
		Float_t biggestSignal =eventReader->getSignal(det,0);
		Float_t secondbiggestSignal =eventReader->getSignal(det,0);
		UInt_t biggest_hit_channel = 0;
		UInt_t second_biggest_hit_channel = 0;
		for (UInt_t ch = 0; ch < TPlaneProperties::getNChannels(det); ch++) {
			if(settings->isDet_channel_screened(det,ch))
				continue;
			Float_t signal = eventReader->getSignal(det,ch);
			//			cout << "event: " << nEvent << "\tdet: " << det << "\tchannel: " << i << "\tadc: " << adcValue << "\tPedMean: " << PedMean << "\tPedSigma: " << PedSigma << "\tsignal: " << signal << endl;
			if (/*i > 70 && i < 170 &&*/ !eventReader->isSaturated(det,ch)){
				if (signal > biggestSignal) {
					second_biggest_hit_channel = biggest_hit_channel;
					biggest_hit_channel = ch;
					secondbiggestSignal = biggestSignal;
					biggestSignal = signal;
				}
				else
					if (signal > secondbiggestSignal) {
						second_biggest_hit_channel = ch;
						secondbiggestSignal = signal;
					}
			}//end if fidcut region
		}//end loop over channels

		Float_t biggestSignalSigma = eventReader->getSignalInSigma(det,biggest_hit_channel);
		Float_t secondbiggestSignalSigma =  eventReader->getSignalInSigma(det,second_biggest_hit_channel);

		if (biggest_hit_channel > 0 && biggest_hit_channel < TPlaneProperties::getNChannels(det)) {
			// -- look for second biggest hit next to biggest hit
			Float_t leftHitSignal= eventReader->getSignal(det,biggest_hit_channel-1);
			Float_t rightHitSignal = eventReader->getSignal(det,biggest_hit_channel+1);
			if (leftHitSignal < rightHitSignal) {
				second_biggest_hit_channel = biggest_hit_channel + 1;
				histo_second_biggest_hit_direction[det]->Fill(1.);
			}
			else {
				second_biggest_hit_channel = biggest_hit_channel - 1;
				histo_second_biggest_hit_direction[det]->Fill(-1.);
			}
			secondbiggestSignal = eventReader->getSignal(det,second_biggest_hit_channel);
			secondbiggestSignalSigma = eventReader->getSignalInSigma(det,second_biggest_hit_channel);

			// -- start filling the histograms
			if(biggestSignalSigma<4)
				continue;
			histo_pulseheight_sigma[det]->Fill(biggestSignalSigma);
			histo_pulseheight_sigma_second[det]->Fill(secondbiggestSignalSigma);

			// -- left chip
			if (biggest_hit_channel < TPlaneProperties::getNChannels(det)/2) {
				histo_pulseheight_left_sigma[det]->Fill(biggestSignalSigma);
				histo_pulseheight_left_sigma_second[det]->Fill(secondbiggestSignalSigma);
			}
			// -- right chip
			else {
				histo_pulseheight_right_sigma[det]->Fill(biggestSignalSigma);
				histo_pulseheight_right_sigma_second[det]->Fill(secondbiggestSignalSigma);
			}

			histo_biggest_hit_map[det]->Fill(biggest_hit_channel);
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
		histoName << "hPulseHeight_BiggestHitChannelInSigma" << TPlaneProperties::getStringForDetector(det) ;
		histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);


		histoName.str("");
		histoName << "hPulseHeight_SecondBiggestHitChannelInSigma_" << TPlaneProperties::getStringForDetector(det);
		histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hSecondBiggestHitMinusBiggestHitPosition_" << TPlaneProperties::getStringForDetector(det);
		histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),3,-1.5,1.5);

		histoName.str("");
		histoName << "hPulseHeightSecondBiggestHitChannelInSigmaLeft" << TPlaneProperties::getStringForDetector(det);
		histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hPulseHeightSecondBiggestHitChannelInSigmaRight" << TPlaneProperties::getStringForDetector(det) ;
		histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

		histoName.str("");
		histoName << "hBiggestHitMap"<< TPlaneProperties::getStringForDetector(det);
		histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),TPlaneProperties::getNChannels(det),0.,TPlaneProperties::getNChannels(det)-1);

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
	for (int det = 0; det < 9; det++) {
		double cut = settings->getClusterSeedFactor(det);
		cout << "saving histogram " << this->histo_pulseheight_sigma[det]->GetName() << ".. with CUT on " <<cut<< endl;

		//set Axis Titiles
		this->histo_pulseheight_sigma[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->histo_pulseheight_sigma[det]->GetYaxis()->SetTitle("number of entries #");
		\
		findPlotRangeForPHHisto(histo_pulseheight_sigma[det],settings->getClusterHitFactor(det));

		TCanvas *c1 = new TCanvas(this->histo_pulseheight_sigma[det]->GetTitle(),this->histo_pulseheight_sigma[det]->GetTitle());
		c1->cd();
		this->histo_pulseheight_sigma[det]->Draw();
		double xCor[] = {cut,cut};
		double yCor[] = {0,this->histo_pulseheight_sigma[det]->GetMaximum()*2};
		TGraph* lineGraph = new TGraph(2,xCor,yCor);
		lineGraph->SetLineColor(kRed);
		lineGraph->SetLineWidth(2);
		lineGraph->Draw("Lsame");
		histSaver->SaveCanvas(c1);;
		//        histSaver->SaveHistogram(this->histo_pulseheight_sigma[det]);
		delete histo_pulseheight_sigma[det];
		delete lineGraph;
		delete c1;
	}
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
		double cut = settings->getClusterHitFactor(det);
		//		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".. with CUT on " <<cut<< endl;
		TCanvas *c1 = new TCanvas(this->histo_pulseheight_sigma_second[det]->GetTitle(),this->histo_pulseheight_sigma_second[det]->GetTitle());
		c1->cd();
		this->histo_pulseheight_sigma_second[det]->Draw();
		double xCor[] = {cut,cut};
		double yCor[] = {0,this->histo_pulseheight_sigma_second[det]->GetMaximum()*2};
		this->histo_pulseheight_sigma_second[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->histo_pulseheight_sigma_second[det]->GetXaxis()->SetTitle("number of entries #");
		findPlotRangeForPHHisto(histo_pulseheight_sigma_second[det],settings->getClusterHitFactor(det));
		TGraph* lineGraph = new TGraph(2,xCor,yCor);
		lineGraph->SetLineColor(kRed);
		lineGraph->SetLineWidth(2);
		lineGraph->Draw("Lsame");
		histSaver->SaveCanvas(c1);;
		delete histo_pulseheight_sigma_second[det];
		delete lineGraph;
		delete c1;
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
		histSaver->SaveHistogram(this->histo_second_biggest_hit_direction[det]);
		//		cout << "saving histogram " << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_biggest_hit_map[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma_second[det]);
		//		cout << "saving histogram " << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma[det]);
		//		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma_second[det]);

		//		delete histo_pulseheight_sigma125[det];
		delete histo_second_biggest_hit_direction[det];
		delete histo_biggest_hit_map[det];
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


