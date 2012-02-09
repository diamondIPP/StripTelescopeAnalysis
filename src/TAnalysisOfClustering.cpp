/*
 * TDeadChannels.cpp
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#include "../include/TAnalysisOfClustering.hh"

TAnalysisOfClustering::TAnalysisOfClustering(int runNumber,int seedSigma,int hitSigma) {
	cout<<"\n\n\n\n**********************************************************"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"*********TAnalysisOfClustering::TAnalysisOfClustering*****"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"**********************************************************\n\n\n"<<endl;
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
	sys->MakeDirectory("clustering");
	sys->cd("clustering");
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
	cout<<"analyis clusterin results..."<<endl;
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
		checkForDeadChannels();
		checkForSaturatedChannels();
//		getBiggestHit();//not working
		analyseForSeeds();
		analyseCluster();
		compareCentroid_ChargeWeightedMean();
		analyse2ndHighestHit();
//		analyseBiggestHit(); // moved to TAnalysisOfPedestal.cpp
	}
	saveHistos();
}

void TAnalysisOfClustering::checkForDeadChannels()
{
	for(int det=0;det<9;det++){
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
			int nClusters = eventReader->getNClusters(det);
			if(nClusters==1)
				hSeedMap2[det]->Fill(eventReader->getCluster(det,0).getHighestSignalChannel());
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
	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		if(eventReader->getAdcValue(8,ch)>=1024){
			hSaturatedChannels[8]->Fill(ch);
		}
	}
}

void TAnalysisOfClustering::initialiseHistos()
{
	cout<<"1"<<endl;
	{
		stringstream histoName;
		histoName<<"hDiamond_Delta_CWM_BiggestHit";
		histo_CWM_biggestHit=new TH2F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6,10,0,9);
		histoName.str("");
		histoName<<"hDiamond_Delta_highest2Centroid_BiggestHit";
		histo_H2C_biggestHit=new TH1F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6);
	}
	cout<<"2"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"SaturatedChannels_"<<TADCEventReader::getStringForPlane(det)<<"";
		hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		if(det==8)hSaturatedChannels[det]->GetXaxis()->SetRangeUser(0,128);
	}
	cout<<"3"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPositionOfallSeeds_"<<TADCEventReader::getStringForPlane(det);
		hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		if(det==8)hSeedMap[det]->GetXaxis()->SetRangeUser(0,128);
	}
	cout<<"4"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hPositionOfHighestSeed_"<<TADCEventReader::getStringForPlane(det);
		hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
		hSeedMap2[det]->GetXaxis()->SetTitle("Position of Highest Seed of a Cluster");
		if(det==8)hSeedMap2[det]->GetXaxis()->SetRangeUser(0,128);
	}
	cout<<"5"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hNumberOfSeeds_in_"<<TADCEventReader::getStringForPlane(det);
		hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
		hNumberOfSeeds[det]->GetXaxis()->SetTitle("Number Of Seeds in Cluster");
		hNumberOfSeeds[det]->GetYaxis()->SetTitle("Entries #");
	}
	cout<<"6"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHitChannelInSigma";
		hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	cout<<"7"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"PulseHeight_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHitNextToBiggestHit_ChannelInSigma";
		hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
	}
	cout<<"8"<<endl;
	for (int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"Channel_"<<TADCEventReader::getStringForPlane(det)<<"_BiggestHit";
		hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,255);
	}
	cout<<"9"<<endl;
	for(int det=0;det<9;det++){
		stringstream histoName;
		histoName<<"hClusterSize_"<<TADCEventReader::getStringForPlane(det);
		hClusterSize[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
		hClusterSize[det]->GetXaxis()->SetTitle("Size Of Cluster in Channels");
		hClusterSize[det]->GetYaxis()->SetTitle("Entries #");
		histoName.str("");
		histoName<<"NumberOfClusters_"<<TADCEventReader::getStringForPlane(det);
		hNumberOfClusters[det]= new TH1F(histoName.str().c_str(),histoName.str().c_str(),10,-0.5,10.5);
	}
	cout<<"10"<<endl;
    for (int det = 0; det < 9; det++) {
		int nbins = 250;
		Float_t min = 0.;
		Float_t max = 250.;
		
        stringstream histoName;
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "BiggestHitChannelInSigma";
        histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigma";
        histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitMinusBiggestHitPosition";
        histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),2,-2.,2.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigmaLeft";
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "SecondBiggestHitChannelInSigmaRight";
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        
        histoName.str("");
        histoName << TADCEventReader::getStringForPlane(det) << "BiggestHitMap";
        histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0.,255.);
        
        histoName.str("");
        histoName << "PulseHeight" << TADCEventReader::getStringForPlane(det) << "LeftChipBiggestHitChannelInSigma";
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
    cout<<"11"<<endl;
    for(int det=0;det<9;det++){//analayse2ndHighestHit
    	stringstream histName;
    	histName<<"h2ndBiggestHitSignal_"<<TADCEventReader::getStringForPlane(det);
    	if(det<8)
    		h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,200);
    	else
    		h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1024);
    	h2ndBiggestHitSignal[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster");
    	h2ndBiggestHitSignal[det]->GetYaxis()->SetTitle("Entries #");
    	histName.str("");
    	histName<<"h2ndBiggestHitOverCharge_"<<TADCEventReader::getStringForPlane(det);
    	h2ndBiggestHitOverCharge[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,0.5);
    	h2ndBiggestHitOverCharge[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster over Sum of all signals of cluster");
    	h2ndBiggestHitOverCharge[det]->GetYaxis()->SetTitle("Entries #");
    	histName.str("");
    	histName<<"h2ndBiggestHitPosition_"<<TADCEventReader::getStringForPlane(det);
    	h2ndBiggestHitPosition[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),3,-1.5,1.5);
    	h2ndBiggestHitPosition[det]->GetXaxis()->SetTitle("position of snd biggest hit in respect to biggest Hit");
    	h2ndBiggestHitPosition[det]->GetYaxis()->SetTitle("Entries #");
    	histName.str("");
    	histName<<"hLeftHitOverLeftAndRight_"<<TADCEventReader::getStringForPlane(det);
    	hLeftHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
    	hLeftHitOverLeftAndRight[det]->GetXaxis()->SetTitle("Q_L/(Q_R +Q_L)");
    	histName.str("");
    	histName<<"hDeltaLeftRightHitOverLeftAndRight_"<<TADCEventReader::getStringForPlane(det);
    	hDeltaLeftRightHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,-1,1);
    	hDeltaLeftRightHitOverLeftAndRight[det]->GetXaxis()->SetTitle("(Q_L-Q_R)/(Q_R +Q_L)");
    	histName.str("");
    	histName<<"hSignal2ndHighestOverSignalHighest_"<<TADCEventReader::getStringForPlane(det);
    	hHighestTo2ndHighestSignalRatio[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
    	hHighestTo2ndHighestSignalRatio[det]->GetXaxis()->SetTitle("Q_{2ndHighest}/Q_{Highest}");
    }
    cout<<"12"<<endl;
}





void TAnalysisOfClustering::saveHistos(){
	cout<<"plot histo "<<histo_CWM_biggestHit->GetName();
	histSaver->SaveHistogram(histo_CWM_biggestHit);
	histo_CWM_biggestHit->Delete();
	cout<<"plot histo "<<histo_H2C_biggestHit->GetName();
	histSaver->SaveHistogram(histo_H2C_biggestHit);
	histo_H2C_biggestHit->Delete();
    for(int det=0;det<9;det++){//analyse 2nd biggest Hit
    	cout<<"plot histo "<<det<<"  h2ndBiggestHitSignal_"<<TADCEventReader::getStringForPlane(det);
    	histSaver->SaveHistogram(h2ndBiggestHitSignal[det]);
    	delete h2ndBiggestHitSignal[det];
    	cout<<"plot histo "<<det<<"  h2ndBiggestHitOverCharge_"<<TADCEventReader::getStringForPlane(det);
    	histSaver->SaveHistogram(h2ndBiggestHitOverCharge[det]);
    	delete h2ndBiggestHitOverCharge[det];
    	cout<<"plot histo "<<h2ndBiggestHitPosition[det]->GetName()<<endl;
    	histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
    	histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
    	histSaver->SaveHistogram(hLeftHitOverLeftAndRight[det]);
    	histSaver->SaveHistogram(hDeltaLeftRightHitOverLeftAndRight[det]);
    	histSaver->SaveHistogram(hHighestTo2ndHighestSignalRatio[det]);
    }

	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogram(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();
	}
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
		histSaver->SaveHistogram(hSeedMap[det]);
		hSeedMap[det]->Delete();
	}
	for (int det=0;det<9;det++){
			cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
			histSaver->SaveHistogram(hSeedMap2[det]);
			hSeedMap2[det]->Delete();
		}
	for (int det=0;det<9;det++){
		cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
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
		histSaver->SaveHistogram(this->hNumberOfClusters[det]);
		delete hClusterSize[det];
		delete hNumberOfClusters[det];
	}
    
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
}

void TAnalysisOfClustering::compareCentroid_ChargeWeightedMean()
{
	bool check=true;
	for(int det=0;det<9;det++)
		check=eventReader->getNClusters(det)==1;
	if(check==true){
		Float_t xCWM=eventReader->getCluster(8,0).getChargeWeightedMean();
		Float_t xHit=(Float_t)eventReader->getCluster(8,0).getHighestSignalChannel();
		Float_t xH2C=(Float_t)eventReader->getCluster(8,0).getHighest2Centroid();
		Float_t delta=xCWM-xHit;
		this->histo_CWM_biggestHit->Fill(delta,eventReader->getCluster(8,0).size());
		delta = xH2C - xHit;
		this->histo_H2C_biggestHit->Fill(delta);
//		if(eventReader->getNClusters(8)>=1){
//			Float_t charge = eventReader->getCluster(8,0).getCharge();
//			Float_t signal2ndHighestHit=eventReader->getCluster(8,0).getCharge(2)-eventReader->getCluster(8,0).getCharge(1);
//			Float_t q =signal2ndHighestHit/charge;
//		}


	}
}

void TAnalysisOfClustering::analyseCluster()
{
	for(int det=0;det<9;det++){
		hNumberOfClusters[det]->Fill(eventReader->getNClusters(det));
		for(int cl=0;cl<eventReader->getNClusters(det);cl++){
			hClusterSize[det]->Fill(eventReader->getClusterSize(det,cl));
		}
	}

}


void TAnalysisOfClustering::analyse2ndHighestHit(){
	for(int det=0;det<9;det++){
		for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
			TCluster cluster=eventReader->getCluster(det,cl);
			Float_t signalLeft = cluster.getSignalOfChannel(cluster.getHighestSignalChannel()-1);
			Float_t signalRight = cluster.getSignalOfChannel(cluster.getHighestSignalChannel()+1);
			if(signalLeft<0)
				cout<<"signalLeft is smaller than 0"<<endl;
			if(signalRight<0)
				cout<<"signalLeft is smaller than 0"<<endl;
			Float_t signalHighest = cluster.getHighestSignal();
			Float_t signal2ndHighest;
			Float_t deltaSignals = signalLeft-signalRight;
			if(deltaSignals<0)
				signal2ndHighest=signalRight;
			else
				signal2ndHighest=signalLeft;
			Float_t sumSignals = signalLeft+signalRight;
			Float_t allCharge=cluster.getCharge(true);
			Float_t signalRatio=signal2ndHighest/signalHighest;
			hHighestTo2ndHighestSignalRatio[det]->Fill(signalRatio);
			if(signalLeft>signalRight){
				h2ndBiggestHitSignal[det]->Fill(signalLeft);
				h2ndBiggestHitOverCharge[det]->Fill(signalLeft/allCharge);
			}
			else{
				h2ndBiggestHitSignal[det]->Fill(signalRight);
				h2ndBiggestHitOverCharge[det]->Fill(signalRight/allCharge);
			}
			if(signalLeft>signalRight){
				h2ndBiggestHitPosition[det]->Fill(-1);
			}
			else if(signalLeft<signalRight)
				h2ndBiggestHitPosition[det]->Fill(+1);
			else
				h2ndBiggestHitPosition[det]->Fill(0);
			if (TMath::Abs(deltaSignals/sumSignals)!=1)
				hDeltaLeftRightHitOverLeftAndRight[det]->Fill((deltaSignals)/(sumSignals));
			if(signalLeft!=0)hLeftHitOverLeftAndRight[det]->Fill((signalLeft)/(sumSignals));
		}
	}
}





