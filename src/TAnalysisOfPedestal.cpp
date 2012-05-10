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
	// TODO Auto-generated constructor stub
	if(settings==0)
		exit(0);
	sys = gSystem;
	stringstream  runString;
	UInt_t runNumber=settings->getRunNumber();
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());
	
	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<"pedestalData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str(),settings->getRunNumber());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("pedestalAnalysis");
	sys->cd("pedestalAnalysis");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
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

	htmlPedestal= new THTMLPedestal(settings);
	htmlPedestal->setPathName(plotsPath.str());
	htmlPedestal->setFileName("pedestal.html");
}

TAnalysisOfPedestal::~TAnalysisOfPedestal() {
	// TODO Auto-generated destructor stub


	htmlPedestal->setTitle("Pedestals");
	htmlPedestal->createTableOfCuts();
	htmlPedestal->createPedestalDistribution();
	htmlPedestal->generateHTMLFile();
	delete htmlPedestal;
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}


void TAnalysisOfPedestal::doAnalysis(UInt_t nEvents)
{
	cout<<"analyze pedestal data..."<<endl;
//	eventReader->checkADC();
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
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
	createPedestalMeanHistos();
	saveHistos();
}

void TAnalysisOfPedestal::checkForDeadChannels()
{
	for(int det=0;det<9;det++){
		int numberOfSeeds=0;
		deque< pair<int, Float_t> > seedQueue;
		UInt_t maxChannels=TPlaneProperties::getNChannels(det);
		for(UInt_t ch=0;ch<maxChannels;ch++){
			Float_t sigma=eventReader->getPedestalSigma(det,ch);
			Float_t signal = (Float_t)eventReader->getDet_ADC(det,ch)-eventReader->getPedestalMean(det,ch);
			if(sigma==0){
				//cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
				continue;
			};
			Float_t adcValueInSigma=signal/sigma;
			if(adcValueInSigma>10){
				hSeedMap[det]->Fill(ch);
				//cout<<"Found a Seed "<<det<<" "<<ch<<" "<<adcValueInSigma<<" "<<eventReader->getCurrent_event()<<endl;
				numberOfSeeds++;
				seedQueue.push_back(make_pair(ch,signal));
			}
		}
		hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}
	
}
void TAnalysisOfPedestal::analyseForSeeds(){
	for(int det=0;det<9;det++){
		int nClusters = eventReader->getNClusters(det);
		if(nClusters==1)
			hSeedMap2[det]->Fill(eventReader->getCluster(det,0).getHighestSignalChannel());
	}
}

void TAnalysisOfPedestal::checkForSaturatedChannels()
{
	for(int det=0;det<TPlaneProperties::getNDetectors();det++)
	for(int ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
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
		for(UInt_t ch=70;ch<200;ch++){//todo warum 70 bis 200
			UInt_t adcValue=eventReader->getAdcValue(det,ch);
			if (adcValue<TPlaneProperties::getMaxSignalHeight(det))
				if(adcValue>biggestHit){
					biggestHit=eventReader->getAdcValue(det,ch);
					biggestHitInSigma=eventReader->getSignalInSigma(det,ch);
					chBiggestHit=ch;
				}
		}
		hPulsHeightBiggestHit[det]->Fill(biggestHitInSigma);
		hChannelBiggestHit[det]->Fill(chBiggestHit);
		if(eventReader->getAdcValue(det,chBiggestHit-1)<eventReader->getAdcValue(det,chBiggestHit+1))
			hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit+1));
		else
			hPulsHeightNextBiggestHit[det]->Fill(eventReader->getSignalInSigma(det,chBiggestHit-1));
	}
	
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
		
		if (biggest_hit_channel > 0 && biggest_hit_channel < TPlaneProperties::getNChannels(det) && biggestSignalSigma > settings->getClusterHitFactor(det)) {
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
	for (UInt_t det =0;det<9;det++){
		stringstream histoName,histoTitle,xTitle,yTitle;
		histoName<<"hMeanSignalOfAllNonHitChannels_"<<TADCEventReader::getStringForDetector(det);
		histoTitle<<"signal (adc-pedestal) of all non hit channels in Plane"<<TADCEventReader::getStringForDetector(det);
		xTitle<<"non hit Noise in ADC counts";
		yTitle<<"Number of Entries #";
		Float_t width = 8;
		UInt_t nBins =512;
		if (det==TPlaneProperties::getDetDiamond())
			width = 20;
		hAllAdcNoise[det]= new TH1F(histoName.str().c_str(),histoTitle.str().c_str(),nBins,(-1)*width,width);
		hAllAdcNoise[det]->GetXaxis()->SetTitle(xTitle.str().c_str());
		hAllAdcNoise[det]->GetYaxis()->SetTitle(yTitle.str().c_str());
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSaturatedChannels_"<<TADCEventReader::getStringForDetector(det)<<"";
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSeedPosAllSeeds_"<<TADCEventReader::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hMaxSeedPos_"<<TADCEventReader::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hNumberOfSeeds_"<<TADCEventReader::getStringForDetector(det);
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPulseHeight__BiggestHitChannelInSigma_"<<TADCEventReader::getStringForDetector(det);
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){//todo why such a big histo?so big?
		stringstream histoName;
		histoName<<"hPulseHeight_BiggestHitNextToBiggestHit_ChannelInSigma"<<TADCEventReader::getStringForDetector(det);
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hChannel_BiggestHit_"<<TADCEventReader::getStringForDetector(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels);
	}
    
    for (UInt_t det = 0; det < 9; det++) {
		int nbins = 256;
		Float_t min = 0.;
		Float_t max = 64.;
		if(det==TPlaneProperties::getDetDiamond()){max=256;nbins=1024;}
		
        stringstream histoName;
        histoName << "hPulseHeight_BiggestHitChannelInSigma" << TADCEventReader::getStringForDetector(det) ;
        histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        

        histoName.str("");
        histoName << "hPulseHeight_SecondBiggestHitChannelInSigma_" << TADCEventReader::getStringForDetector(det);
        histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hSecondBiggestHitMinusBiggestHitPosition_" << TADCEventReader::getStringForDetector(det);
        histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),3,-1.5,1.5);
        
        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaLeft" << TADCEventReader::getStringForDetector(det);
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaRight" << TADCEventReader::getStringForDetector(det) ;
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hBiggestHitMap"<< TADCEventReader::getStringForDetector(det);
        histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),TPlaneProperties::getNChannels(det),0.,TPlaneProperties::getNChannels(det)-1);
        
        histoName.str("");
        histoName << "hPulseHeightLeftChipBiggestHitChannelInSigma" << TADCEventReader::getStringForDetector(det);
        histo_pulseheight_left_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForDetector(det) << "RightChipBiggestHitChannelInSigma";
        histo_pulseheight_right_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForDetector(det) << "LeftChipSecondBiggestHitChannelInSigma";
        histo_pulseheight_left_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForDetector(det) << "RightChipSecondBiggestHitChannelInSigma";
        histo_pulseheight_right_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
    }
	
}




void TAnalysisOfPedestal::saveHistos(){
	
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();
	}
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSeedMap[det]);
		hSeedMap[det]->Delete();
	}
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSeedMap2[det]);
		hSeedMap2[det]->Delete();
	}
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hNumberOfSeeds[det]);
		hNumberOfSeeds[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogram(hPulsHeightBiggestHit[det]);
		hPulsHeightBiggestHit[det]->Delete();
	}
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
    
    for (int det = 0; det < 9; det++) {
    	double cut = settings->getClusterSeedFactor(det);
		cout << "saving histogram " << this->histo_pulseheight_sigma[det]->GetName() << ".. with CUT on " <<cut<< endl;
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
    	TGraph* lineGraph = new TGraph(2,xCor,yCor);
    	lineGraph->SetLineColor(kRed);
    	lineGraph->SetLineWidth(2);
    	lineGraph->Draw("Lsame");
    	histSaver->SaveCanvas(c1);;
        delete histo_pulseheight_sigma_second[det];
        delete lineGraph;
        delete c1;
    }
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
		histSaver->SaveHistogram(hAllAdcNoise[det],true);
		delete hAllAdcNoise[det];
    }
}

/**
 *
 */
void TAnalysisOfPedestal::createPedestalMeanHistos()
{
	for(UInt_t det = 0; det<TPlaneProperties::getNDetectors();det++){
		stringstream nameMean,titleMean,titleSigma,nameSigma,canvasTitle,graphTitle;
		nameMean<<"hMeanPedestal_Value_OfChannel_"<<TADCEventReader::getStringForDetector(det);
		nameSigma<<"hMeanPedestal_Width_OfChannel_"<<TADCEventReader::getStringForDetector(det);
		titleMean<<"mean of pedestalValue for each channel of "<<TADCEventReader::getStringForDetector(det);
		titleSigma<<"mean of pedestalWidth for each channel of "<<TADCEventReader::getStringForDetector(det);
		UInt_t nBins = pedestalMeanValue.at(det).size();
		TH1F *histoMean = new TH1F(nameMean.str().c_str(),titleMean.str().c_str(),nBins,-.5,nBins-.5);
		TH1F *histoSigma = new TH1F(nameSigma.str().c_str(),titleSigma.str().c_str(),nBins,-.5,nBins-.5);
		histoMean->GetXaxis()->SetTitle("channel No");
		histoMean->GetYaxis()->SetTitle("mean pedestal value");
		histoSigma->GetXaxis()->SetTitle("channel No");
		histoSigma->GetYaxis()->SetTitle("mean pedestal sigma");
		vector<Float_t> vecChNo,vecChError;
		for(UInt_t ch = 0; ch<TPlaneProperties::getNChannels(det);ch++){
			this->pedestalMeanValue.at(det).at(ch)/=nPedestalHits.at(det).at(ch);
			this->pedestalSigmaValue.at(det).at(ch)/=nPedestalHits.at(det).at(ch);
			histoMean->Fill(ch,pedestalMeanValue.at(det).at(ch));
			histoSigma->Fill(ch,pedestalSigmaValue.at(det).at(ch));
			vecChNo.push_back(ch+.1);
			vecChError.push_back(0);
		}
		Float_t max = histoMean->GetMaximum()*1.1;
		histoSigma->SetLineColor(kRed);
		histSaver->SaveHistogram(histoMean);
		histSaver->SaveHistogram(histoSigma);
		canvasTitle<<"cPedestalOfChannels_"<<TADCEventReader::getStringForDetector(det);
		histSaver->SaveTwoHistos(canvasTitle.str(),histoMean,histoSigma,10);
		TGraphErrors *graph = new TGraphErrors(nBins,&vecChNo.at(0),&pedestalMeanValue.at(det).at(0),&vecChError.at(0),&pedestalSigmaValue.at(det).at(0));
		graph->Draw("APLgoff");
		graph->GetXaxis()->SetTitle("channel No.");
		graph->GetYaxis()->SetTitle("pedestalValue in ADC counts");
		graph->GetYaxis()->SetRangeUser(0,max);
		graph->GetXaxis()->SetRangeUser(0,nBins-1);
		graphTitle<<"gMeanPedestalValueOfChannelWithSigmaAsError_"<<TADCEventReader::getStringForDetector(det);
		graph->SetTitle(graphTitle.str().c_str());
		histSaver->SaveGraph(graph,graphTitle.str(),"AP");
		delete histoMean;
	}


}

void TAnalysisOfPedestal::updateMeanCalulation(){
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			Float_t snr = eventReader->getSignalInSigma(det,ch);
			Float_t pedestal = eventReader->getPedestalMean(det,ch);
			Float_t sigma = eventReader->getPedestalSigma(det,ch);
			UInt_t adc = eventReader->getAdcValue(det,ch);
			Float_t noise = adc-pedestal;
			if(snr<settings->getClusterHitFactor(det)){
//				cout<<"\n1 "<<det<<"/"<<ch<<" "<<pedestalMeanValue.size()<<"-"<<pedestalMeanValue.at(det).size()<<flush;
				pedestalMeanValue.at(det).at(ch) +=pedestal;
//				cout<<"\n2 "<<det<<"/"<<ch<<" "<<pedestalSigmaValue.size()<<"-"<<pedestalSigmaValue.at(det).size()<<flush;
				pedestalSigmaValue.at(det).at(ch) +=sigma;
//				cout<<"\n3 "<<det<<"/"<<ch<<" "<<nPedestalHits.size()<<"-"<<nPedestalHits.at(det).size()<<flush;
				nPedestalHits.at(det).at(ch)++;
//				cout<<noise<<endl;
				hAllAdcNoise[det]->Fill(noise);
//				cout<<"$"<<flush;
			}
			if(TPlaneProperties::getDetDiamond()==det){
//				cout<<"*"<<flush;
				diaRawADCvalues.at(ch).push_back(adc);
			}

		}

}


