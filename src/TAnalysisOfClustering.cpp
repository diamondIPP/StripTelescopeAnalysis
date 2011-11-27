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
//		checkForDeadChannels();
//		checkForSaturatedChannels();
//		getBiggestHit();
//		analyseForSeeds();
		analyseCluster();
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






