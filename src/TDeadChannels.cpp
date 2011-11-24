/*
 * TDeadChannels.cpp
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#include "../include/TDeadChannels.hh"

TDeadChannels::TDeadChannels(int runNumber,int seedSigma,int hitSigma) {
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

TDeadChannels::~TDeadChannels() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}


void TDeadChannels::doAnalysis(int nEvents)
{
	cout<<"find dead channels..."<<endl;
	if(nEvents!=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		cout<<nEvent;
		for(unsigned int det=0;det< (eventReader->getCluster()->size());det++)
			for(unsigned int cl=0;cl< eventReader->getCluster()->at(det).size();cl++)
			cout<<" "<<eventReader->getCluster()->at(det).at(cl).getChargeWeightedMean()<<flush;//*/
		cout<<endl;
		checkForDeadChannels();
		checkForSaturatedChannels();
	}
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
	for (int det=0;det<8;det++){
		cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hNumberOfSeeds[det]);
		hNumberOfSeeds[det]->Delete();
	}
}

void TDeadChannels::checkForDeadChannels()
{
	for(int det=0;det<8;det++){
		int numberOfSeeds=0;
		deque< pair<int, Float_t> > seedQueue;
	for(int ch=0;ch<N_DET_CHANNELS;ch++){
		Float_t sigma=eventReader->getPedestalSigma(det,ch);
		Float_t signal = (Float_t)eventReader->getDet_ADC(det,ch)-eventReader->getPedestalMean(det,ch);
		if(sigma==0){
				cout<<nEvent<<" "<<det<<ch<<" sigma==0"<<endl;
				continue;
		};
		Float_t adcValueInSigma=signal/sigma;
		if(adcValueInSigma>10){
			hSeedMap[det]->Fill(ch);
//			cout<<"Found a Seed "<<det<<" "<<ch<<" "<<adcValueInSigma<<" "<<eventReader->getCurrent_event()<<endl;
			numberOfSeeds++;
			seedQueue.push_back(make_pair(ch,signal));
		}
	}
	hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}

}

void TDeadChannels::initialiseHistos()
{
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<TADCEventReader::getStringForPlane(det)<<"_SaturatedChannels";
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<TADCEventReader::getStringForPlane(det)<<"_SeedMap";
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	for (int det=0;det<8;det++){
		stringstream histoName;
		histoName<<TADCEventReader::getStringForPlane(det)<<"_NumberOfSeeds";
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
	}
}



void TDeadChannels::checkForSaturatedChannels()
{
	for(int det=0;det<8;det++)
	for(int ch=0;ch<N_DET_CHANNELS;ch++){
		if(eventReader->getDet_ADC(det,ch)>=255){
			hSaturatedChannels[det]->Fill(ch);
		}
	}
}






