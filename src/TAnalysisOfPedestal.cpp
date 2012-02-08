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
	eventReader=new TADCEventReader(filepath.str());
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
}

TAnalysisOfPedestal::~TAnalysisOfPedestal() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}


void TAnalysisOfPedestal::doAnalysis(int nEvents)
{
	cout<<"analyze pedestal data..."<<endl;
//	eventReader->checkADC();
	if(nEvents!=0) nEvents=eventReader->GetEntries();
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
//		checkForSaturatedChannels();
		//		getBiggestHit();
//		analyseForSeeds();
//		analyseCluster();
		analyseBiggestHit();
	}
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
	for(int det=0;det<9;det++){
		UInt_t maxChannel = TPlaneProperties::getNChannels(det);
		for(UInt_t ch=0;ch<maxChannel;ch++){
			if(eventReader->getAdcValue(det,ch)>=254){
				hSaturatedChannels[det]->Fill(ch);
			}
		}
	}
}
void TAnalysisOfPedestal::getBiggestHit(){
	
	for(int det=0;det<9;det++){
		int biggestHit=0;
		Float_t biggestHitInSigma=0;
		int chBiggestHit;
		for(int ch=70;ch<200;ch++){
			int adcValue=eventReader->getAdcValue(det,ch);
			if (adcValue<253)
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

void TAnalysisOfPedestal::analyseBiggestHit() {
    for (int det = 0; det < 8; det++) {//todo change to 9
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
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSaturatedChannels_"<<TADCEventReader::getStringForPlane(det)<<"";
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hSeedPosAllSeeds_"<<TADCEventReader::getStringForPlane(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hMaxSeedPos_"<<TADCEventReader::getStringForPlane(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels-1);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hNumberOfSeeds_"<<TADCEventReader::getStringForPlane(det);
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPulseHeight__BiggestHitChannelInSigma_"<<TADCEventReader::getStringForPlane(det);
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){//todo why such a big histo?so big?
		stringstream histoName;
		histoName<<"hPulseHeight_BiggestHitNextToBiggestHit_ChannelInSigma"<<TADCEventReader::getStringForPlane(det);
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hChannel_BiggestHit_"<<TADCEventReader::getStringForPlane(det);
		UInt_t nChannels=TPlaneProperties::getNChannels(det);
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels);
	}
    
    for (int det = 0; det < 9; det++) {
		int nbins = 256;
		Float_t min = 0.;
		Float_t max = 64.;
		
        stringstream histoName;
        histoName << "hPulseHeight_BiggestHitChannelInSigma" << TADCEventReader::getStringForPlane(det) ;
        histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        

        histoName.str("");
        histoName << "hPulseHeight_SecondBiggestHitChannelInSigma_" << TADCEventReader::getStringForPlane(det);
        histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hSecondBiggestHitMinusBiggestHitPosition_" << TADCEventReader::getStringForPlane(det);
        histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),3,-1.5,1.5);
        
        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaLeft" << TADCEventReader::getStringForPlane(det);
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaRight" << TADCEventReader::getStringForPlane(det) ;
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "hBiggestHitMap"<< TADCEventReader::getStringForPlane(det);
        histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),TPlaneProperties::getNChannels(det),0.,TPlaneProperties::getNChannels(det)-1);
        
        histoName.str("");
        histoName << "hPulseHeightLeftChipBiggestHitChannelInSigma" << TADCEventReader::getStringForPlane(det);
        histo_pulseheight_left_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "RightChipBiggestHitChannelInSigma";
        histo_pulseheight_right_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "LeftChipSecondBiggestHitChannelInSigma";
        histo_pulseheight_left_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "RightChipSecondBiggestHitChannelInSigma";
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
//	for(int det=0;det<9;det++){
//		histSaver->SaveHistogramPNG(hPulsHeightBiggestHit[det]);
//		hPulsHeightBiggestHit[det]->Delete();
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
    
    for (int det = 0; det < 9; det++) {
		cout << "saving histogram " << this->histo_pulseheight_sigma[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(this->histo_pulseheight_sigma[det]);
		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_sigma_second[det]);
		//		cout << "saving histogram" << this->histo_pulseheight_sigma125[det]->GetName() << ".." << endl;
		//		histSaver->SaveHistogramPNG(this->histo_pulseheight_sigma125[det]);
		cout << "saving histogram " << this->histo_second_biggest_hit_direction[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_second_biggest_hit_direction[det]);
		cout << "saving histogram " << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_biggest_hit_map[det]);
		cout << "saving histogram " << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma[det]);
		cout << "saving histogram " << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_left_sigma_second[det]);
		cout << "saving histogram " << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma[det]);
		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogram(this->histo_pulseheight_right_sigma_second[det]);
        delete histo_pulseheight_sigma[det];
		delete histo_pulseheight_sigma_second[det];
		//		delete histo_pulseheight_sigma125[det];
		delete histo_second_biggest_hit_direction[det];
		delete histo_biggest_hit_map[det];
		delete histo_pulseheight_left_sigma[det];
		delete histo_pulseheight_left_sigma_second[det];
		delete histo_pulseheight_right_sigma[det];
		delete histo_pulseheight_right_sigma_second[det];
    }
}





