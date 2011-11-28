/*
 * TDeadChannels.cpp
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#include "../include/TAnalysisOfClustering.hh"

TAnalysisOfClustering::TAnalysisOfClustering(int runNumber,int seedSigma,int hitSigma) {
	// TODO Auto-generated constructor stub
	sys = gSystem;
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());

	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<"clusterData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("deadChannels");
	sys->cd("deadChannels");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	initialiseHistos();
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
	cout<<"end initialise"<<endl;
}

TAnalysisOfClustering::~TAnalysisOfClustering() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}


void TAnalysisOfClustering::doAnalysis(int nEvents)
{
	cout<<"find dead channels..."<<endl;
	if(nEvents!=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		/*cout<<nEvent;
		for(unsigned int det=0;det< (eventReader->getCluster()->size());det++)
			for(unsigned int cl=0;cl< eventReader->getCluster()->at(det).size();cl++)
			cout<<" "<<eventReader->getCluster()->at(det).at(cl).getChargeWeightedMean()<<flush;
		cout<<endl;//*/
		checkForDeadChannels();
		checkForSaturatedChannels();
		getBiggestHit();
		analyseForSeeds();
		analyseCluster();
		analyseBiggestHit();
	}
	saveHistos();
}

void TAnalysisOfClustering::checkForDeadChannels()
{
	for(int det=0;det<8;det++){
		int numberOfSeeds=0;
		deque< pair<int, Float_t> > seedQueue;
	for(int ch=0;ch<N_DET_CHANNELS;ch++){
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
void TAnalysisOfClustering::analyseForSeeds(){
	for(int det=0;det<9;det++){
			int nClusters = eventReader->getCluster()->at(det).size();
			if(nClusters==1)
				hSeedMap2[det]->Fill(eventReader->getCluster()->at(det).at(0).getMaximumChannel());
	}
}

void TAnalysisOfClustering::checkForSaturatedChannels()
{
	for(int det=0;det<8;det++)
	for(int ch=0;ch<N_DET_CHANNELS;ch++){
		if(eventReader->getAdcValue(det,ch)>=254){
			hSaturatedChannels[det]->Fill(ch);
		}
	}
}
void TAnalysisOfClustering::getBiggestHit(){

	for(int det=0;det<8;det++){
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

void TAnalysisOfClustering::analyseBiggestHit() {
    for (int det = 0; det < 8; det++) {
		Float_t biggestSignal = (Float_t)eventReader->getAdcValue(det,0)-eventReader->getPedestalMean(det,0);
		Float_t secondbiggestSignal = (Float_t)eventReader->getAdcValue(det,0)-eventReader->getPedestalMean(det,0);
		int biggest_hit_channel = 0;
		int second_biggest_hit_channel = 0;
        for (int i = 0; i < 253; i++) {
			int adcValue = eventReader->getAdcValue(det,i);
			Float_t PedMean = eventReader->getPedestalMean(det,i);
			Float_t PedSigma = eventReader->getPedestalSigma(det,i);
			Float_t signal = (Float_t)adcValue-PedMean;
			if (i > 70 && i < 170 && adcValue < 255) {
				if (signal > biggestSignal) {
					second_biggest_hit_channel = biggest_hit_channel;
					biggest_hit_channel = i;
					secondbiggestSignal = biggestSignal;
					biggestSignal = signal;
				}
				else
					if (signal > secondbiggestSignal) {
						second_biggest_hit_channel = i;
						secondbiggestSignal = signal;
					}
//            if (det < 8 && i < 253 && i > 70 && i < 170) {
//                if (Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i]) && Det_ADC[det][i] < 255) {
//                    if (Det_ADC[det][i]-Det_PedMean[det][i] > Det_ADC[det][biggest_hit_channel]-Det_PedMean[det][biggest_hit_channel]/* && Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])*/) {
//                        second_biggest_hit_channel = biggest_hit_channel;
//                        biggest_hit_channel = i;
//                    }
//                    else {
//                        if (Det_ADC[det][i]-Det_PedMean[det][i] > Det_ADC[det][second_biggest_hit_channel]-Det_PedMean[det][second_biggest_hit_channel]) {
//                            second_biggest_hit_channel = i;
//                        }
//                    }
//                }
//                else {
//                    if (Det_ADC[det][i] < 255) {
//                        cout << endl << endl;
//                        cout << "current_event: " << current_event << " Det" << det << "\tchannel: " << i << " Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i]): " << Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i]) << endl;
//                        cout << "Det_ADC: " << (int)Det_ADC[det][i] << "\tDet_PedMean: " << Det_PedMean[det][i] << "\tDetPedWidth: " << Det_PedWidth[det][i] << "\tPH: " << (Det_ADC[det][i]-Det_PedMean[det][i]) / Det_PedWidth[det][i] << endl;
//                    }
//                }
//                if ((int)Det_Channels[det][i] == 125 && Det_ADC[det][i] < 255) {
//                    histo_pulseheight_sigma125[det]->Fill((Det_ADC[det][i]-Det_PedMean[det][i]) / Det_PedWidth[det][i]);
//                }
//            }
			}//end if fidcut region
        }//end loop over channels
		
		Float_t biggestSignalSigma = biggestSignal / eventReader->getPedestalSigma(det,biggest_hit_channel);
		Float_t secondbiggestSignalSigma = secondbiggestSignal / eventReader->getPedestalSigma(det,second_biggest_hit_channel);
		
		if (biggest_hit_channel > 0 && biggest_hit_channel < 255 && biggestSignalSigma > 10) {
			// -- look for second biggest hit next to biggest hit
			if ((Float_t)eventReader->getAdcValue(det,biggest_hit_channel-1)-eventReader->getPedestalMean(det,biggest_hit_channel-1) < (Float_t)eventReader->getAdcValue(det,biggest_hit_channel+1)-eventReader->getPedestalMean(det,biggest_hit_channel+1)) {
				second_biggest_hit_channel = biggest_hit_channel + 1;
				histo_second_biggest_hit_direction[det]->Fill(1.);
			}
			else {
				second_biggest_hit_channel = biggest_hit_channel - 1;
				histo_second_biggest_hit_direction[det]->Fill(-1.);
			}
			secondbiggestSignal = (Float_t)eventReader->getAdcValue(det,second_biggest_hit_channel)-eventReader->getPedestalMean(det,second_biggest_hit_channel);
			secondbiggestSignalSigma = secondbiggestSignal / eventReader->getPedestalSigma(det,second_biggest_hit_channel);
			
			// -- start filling the histograms
			histo_pulseheight_sigma[det]->Fill(biggestSignalSigma);
			histo_pulseheight_sigma_second[det]->Fill(secondbiggestSignalSigma);
			
			// -- left chip
			if (biggest_hit_channel < 128) {
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

void TAnalysisOfClustering::initialiseHistos()
{
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<"SaturatedChannels_"<<TADCEventReader::getStringForPlane(det)<<"";
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<"SeedPos_"<<TADCEventReader::getStringForPlane(det)<<"_allSeeds";
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"maxSeedPos_"<<TADCEventReader::getStringForPlane(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<TADCEventReader::getStringForPlane(det)<<"_NumberOfSeeds";
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHitChannelInSigma";
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHitNextToBiggestHit_ChannelInSigma";
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<"Channel_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHit";
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}

	for(int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"ClusterSize_"<<TADCEventReader::getStringForPlane(det);
		hClusterSize[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
		histoName.str("");
		histoName<<"NumberOfClusters_"<<TADCEventReader::getStringForPlane(det);
		hNumberOfClusters[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
	}
    
    for (int det = 0; det < 8; det++) {
        stringstream histoName;
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "BiggestHitChannelInSigma";
        histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigma";
        histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitMinusBiggestHitPosition";
        histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),2,-2.,2.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigmaLeft";
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigmaRight";
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << TADCEventReader::getStringForPlane(det) << "BiggestHitMap";
        histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0.,255.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "LeftChipBiggestHitChannelInSigma";
        histo_pulseheight_left_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "RightChipBiggestHitChannelInSigma";
        histo_pulseheight_right_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "LeftChipSecondBiggestHitChannelInSigma";
        histo_pulseheight_left_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "RightChipSecondBiggestHitChannelInSigma";
        histo_pulseheight_right_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),250,0.,50.);
    }

}




void TAnalysisOfClustering::saveHistos(){

	for (int det=0;det<8;det++){
		cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();
	}
	for (int det=0;det<8;det++){
		cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSeedMap[det]);
		hSeedMap[det]->Delete();
	}
	for (int det=0;det<9;det++){
			cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
			histSaver->SaveHistogramPNG(hSeedMap2[det]);
			hSeedMap2[det]->Delete();
		}
	for (int det=0;det<8;det++){
		cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hNumberOfSeeds[det]);
		hNumberOfSeeds[det]->Delete();
	}
	for(int det=0;det<8;det++){
			histSaver->SaveHistogramPNG(hPulsHeightBiggestHit[det]);
			hPulsHeightBiggestHit[det]->Delete();
	}
	for(int det=0;det<8;det++){
		histSaver->SaveHistogramPNG(hPulsHeightNextBiggestHit[det]);
		hPulsHeightNextBiggestHit[det]->Delete();
	}
	for(int det=0;det<8;det++){
		histSaver->SaveHistogramPNG(hChannelBiggestHit[det]);
		hChannelBiggestHit[det]->Delete();
	}
	for(int det=0;det<9;det++){
		histSaver->SaveHistogramPNG(this->hClusterSize[det]);
		histSaver->SaveHistogramPNG(this->hNumberOfClusters[det]);
		delete hClusterSize[det];
		delete hNumberOfClusters[det];
	}
    
    for (int det = 0; det < 8; det++) {
		cout << "saving histogram" << this->histo_pulseheight_sigma[det]->GetName() << ".." << endl;
        histSaver->SaveHistogramPNG(this->histo_pulseheight_sigma[det]);
		cout << "saving histogram" << this->histo_pulseheight_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_pulseheight_sigma_second[det]);
//		cout << "saving histogram" << this->histo_pulseheight_sigma125[det]->GetName() << ".." << endl;
//		histSaver->SaveHistogramPNG(this->histo_pulseheight_sigma125[det]);
		cout << "saving histogram" << this->histo_second_biggest_hit_direction[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_second_biggest_hit_direction[det]);
		cout << "saving histogram" << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_biggest_hit_map[det]);
		cout << "saving histogram" << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_pulseheight_left_sigma[det]);
		cout << "saving histogram" << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_pulseheight_left_sigma_second[det]);
		cout << "saving histogram" << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_pulseheight_right_sigma[det]);
		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
		histSaver->SaveHistogramPNG(this->histo_pulseheight_right_sigma_second[det]);
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

void TAnalysisOfClustering::analyseCluster()
{
	for(int det=0;det<9;det++){
		hNumberOfClusters[det]->Fill(eventReader->getCluster()->at(det).size());
		for(int cl=0;cl<eventReader->getCluster()->at(det).size();cl++){
			hClusterSize[det]->Fill(eventReader->getCluster()->at(det).at(cl).size());
		}
	}

}






