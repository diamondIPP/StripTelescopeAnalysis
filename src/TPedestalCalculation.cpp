/*
 * PedestalCalculation.cpp
 *
 *  Created on: 10.11.2011
 *      Author: bachmair
 */

#include "../include/TPedestalCalculation.hh"

TPedestalCalculation::TPedestalCalculation(int runNumber,int nEvents) {
	// TODO Auto-generated constructor stub
	slidingLength=1000;
	eventReader=NULL;
	pedestalTree=NULL;
	pedestalFile=NULL;
	this->runNumber=runNumber;
	sys = gSystem;
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());

	sys->cd(runString.str().c_str());

	stringstream rawfilepath;
	rawfilepath.str("");
	rawfilepath<<"rawData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<rawfilepath.str()<<endl;
	eventReader=new TADCEventReader(rawfilepath.str());
	cout<<eventReader->GetEntries()<<endl;
	createPedestalTree(nEvents);
}

TPedestalCalculation::~TPedestalCalculation() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	pedestalFile->cd();
	pedestalTree->Write();
	pedestalTree->Delete();
	pedestalFile->Close();
}


void TPedestalCalculation::calculatePedestals(){
	double meanSquared[9][256];
	cout<<"initialise arrays"<<endl;
	for(int det=0;det<9;det++)
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			meanValues[det][ch]=0;
			meanSquared[det][ch]=0;
		}

	cout<<"get mean and sigma"<<endl;
	/*
	 * calulate Pedestal mean and sigma with trick: ped=<x> and sigma=<x^2>-<x>^2
	 */
	for(int event=0;event<eventReader->GetEntries();event++){
		eventReader->GetEvent(event);
		for(int det=0;det <8;det++){
			for(int ch=0;ch<256;ch++){
				meanValues[det][ch]+=eventReader->getDet_ADC(det,ch);
				meanSquared[det][ch]+=eventReader->getDet_ADC(det,ch)*(int)eventReader->getDet_ADC(det,ch);
			}
		}
		for(int ch=0;ch<N_DIA_CHANNELS;ch++){
			meanValues[8][ch]+=eventReader->getDia_ADC(ch);
			meanSquared[8][ch]+=eventReader->getDia_ADC(ch)*eventReader->getDia_ADC(ch);
		}
	}
	for(int det=0;det<9;det++)
	for(int ch=0;ch<N_DET_CHANNELS;ch++){
		meanValues[det][ch]=meanValues[det][ch]/(double)eventReader->GetEntries();
		meanSquared[det][ch]=meanSquared[det][ch]/(double)eventReader->GetEntries();
		sigmaValues[det][ch]=meanSquared[det][ch]-meanValues[det][ch]*meanValues[det][ch];
		sigmaValues[det][ch]=TMath::Sqrt(sigmaValues[det][ch]);
	}
	cout<<"DONE"<<endl;

}

void TPedestalCalculation::calculateSlidingPedestals(int nEvents){
	cout<<"calculate Sliding Pedestals"<<endl;

	if(pedestalTree->GetEntries()>=nEvents){
		cout<<"NO Sliding PEdestal Calculation needed, calculations already done."<<endl;
		return;
	}
	TStopwatch watch;
	watch.Start(true);
	//clear adcValues
	for(int det=0;det <8;det++){
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			detAdcValues[det][ch].clear();
		}
	}
	for(int ch=0;ch<N_DIA_CHANNELS;ch++)
		diaAdcValues[ch].clear();

//	int nEvents=eventReader->GetEntries();
	TRawEventSaver::showStatusBar(0,nEvents,true);

	//initialise detAdcValues, diaAdcValues with values from rawTree
	for(nEvent=0;nEvent<slidingLength;nEvent++){
		eventReader->GetEvent(nEvent);
		for(int det=0;det <8;det++){
			for(int ch=0;ch<N_DET_CHANNELS;ch++){
				detAdcValues[det][ch].push_back(eventReader->getDet_ADC(det,ch));
			}
		}
		for(int ch=0;ch<N_DIA_CHANNELS;ch++)
			diaAdcValues[ch].push_back(eventReader->getDia_ADC(ch));
	}

	calculateFirstPedestals(detAdcValues,diaAdcValues);

	//save Sliding Pedestal Values for first slidingLength Events
	for(nEvent=0;nEvent<slidingLength;nEvent++){
		//Fill tree
		pedestalTree->Fill();
	}

	TRawEventSaver::showStatusBar(slidingLength,nEvents,true);

	//calculate sliding Pedestal Values for rest of Events and save them



	for(nEvent=slidingLength;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		//Add next Event to detAdcValues, diaAdcValues
		//Remove first Event from Queue
		eventReader->GetEvent(nEvent);
		for(int det=0;det <8;det++){
			for(int ch=0;ch<N_DET_CHANNELS;ch++){
				detAdcValues[det][ch].push_back(eventReader->getDet_ADC(det,ch));
				pair<float,float> values;
				values=	checkPedestalDet(det,ch);
				pedestalMean[det][ch]=values.first;
				pedestalSigma[det][ch]=values.second;

				detEventUsed[det][ch].pop_front();
				detAdcValues[det][ch].pop_front();
			}
		}
		for(int ch=0;ch<N_DIA_CHANNELS;ch++){
			diaAdcValues[ch].push_back(eventReader->getDia_ADC(ch));

			pair<float,float> values;
			values = checkPedestalDia(ch);
			pedestalMean[8][ch]=values.first;
			pedestalSigma[8][ch]=values.second;


			//if(ch==0&&(nEvent-slidingLength)%10000==0)cout<<nEvent<<". event, ch"<<ch<<"\t"<<values.first<<"+/-"<<values.second<<endl;

			diaEventUsed[ch].pop_front();
			diaAdcValues[ch].pop_front();
		}
		//calculateCurrentPedestals(detAdcValues,diaAdcValues);
		pedestalTree->Fill();
	}//end for
	watch.Stop();
	cout<<"\nStopWatch:"<<endl;
	watch.Print();
}

void TPedestalCalculation::calculateFirstPedestals(deque<UChar_t> DetAdcQueue[8][N_DET_CHANNELS], deque<UShort_t> DiaAdcQueue[N_DIA_CHANNELS]){
//	this->pedestalMean.clear();
//	this->pedestalMean.resize(9);
//	this->pedestalSigma.clear();
//	this->pedestalSigma.resize(9);
	for(int det=0;det <8;det++){
//		pedestalMean.at(det).resize(N_DET_CHANNELS);
//		pedestalSigma.at(det).resize(N_DET_CHANNELS);
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			pair<float,float> values;
			values=this->calculateFirstPedestalDet(det,ch,DetAdcQueue[det][ch],meanValues[det][ch],sigmaValues[det][ch]);
			pedestalMean[det][ch]=values.first;
			pedestalSigma[det][ch]=values.second;
//			pedestalMean.at(det).at(ch)=values.first;
//			pedestalSigma.at(det).at(ch)=values.second;
		}
	}
	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		pair<float,float> values;
		values=this->calculateFirstPedestalDia(ch,DiaAdcQueue[ch],meanValues[8][ch],sigmaValues[8][ch]);
		pedestalMean[8][ch]=values.first;
		pedestalSigma[8][ch]=values.second;
	}
}


pair <float,float> TPedestalCalculation::calculateFirstPedestalDet(int det,int ch,deque<UChar_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	detSUM[det][ch]=0;
	detSUM2[det][ch]=0;
	detEventsInSum[det][ch]=0;
	this->detEventUsed[det][ch].clear();
	for(nEvent=0;nEvent<adcQueue.size();nEvent++){
		if(adcQueue.at(nEvent)<meanChannel+sigmaChannel*maxSigma){
			detEventUsed[det][ch].push_back(true);
			detSUM[det][ch]+=adcQueue.at(nEvent);
			detSUM2[det][ch]+=adcQueue.at(nEvent)*adcQueue.at(nEvent);
			detEventsInSum[det][ch]++;
		}//end if
		else
			detEventUsed[det][ch].push_back(false);
	}//end for nEvent
	float mean=detSUM[det][ch]/(float)detEventsInSum[det][ch];
	float sigma=TMath::Sqrt(detSUM2[det][ch]/(float)detEventsInSum[det][ch]-mean*mean);
	pair<float,float> output = make_pair(mean,sigma);
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDet(det,ch,adcQueue,mean,sigma,iterations-1,maxSigma);
}

pair <float,float> TPedestalCalculation::calculateFirstPedestalDia(int ch,deque<UShort_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	diaSUM[ch]=0;
	diaSUM2[ch]=0;
	diaEventsInSum[ch]=0;
	this->diaEventUsed[ch].clear();
	for(int nEvent=0;nEvent<adcQueue.size();nEvent++){
		if(adcQueue.at(nEvent)<meanChannel+sigmaChannel*maxSigma){
			diaEventUsed[ch].push_back(true);
			diaSUM[ch]+=adcQueue.at(nEvent);
			diaSUM2[ch]+=(adcQueue.at(nEvent)*adcQueue.at(nEvent));
			diaEventsInSum[ch]++;
		}//end if
		else
			diaEventUsed[ch].push_back(false);
	}//end for nEvent
	float mean=(float)diaSUM[ch]/(float)diaEventsInSum[ch];
	float sigma=TMath::Sqrt( ((float)diaSUM2[ch]/(float)diaEventsInSum[ch])-mean*mean);
	pair<float,float> output = make_pair(mean,sigma);
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDia(ch,adcQueue,mean,sigma,iterations-1,maxSigma);
}

pair<float,float> TPedestalCalculation::checkPedestalDet(int det,int ch,int maxSigma){
	if(detEventUsed[det][ch].size()!=slidingLength)
		cout<<"detEventInUse has wrong length"<<detEventUsed[det][ch].size();
	if(detAdcValues[det][ch].size()!=slidingLength+1)
		cout<<"detAdcValues has wrong length..."<<detAdcValues[det][ch].size()<<endl;
	float mean =this->detSUM[det][ch]/(float)this->detEventsInSum[det][ch];
	float sigma=this->detSUM2[det][ch]/(float)this->detEventsInSum[det][ch]-mean*mean;
	if(this->detEventUsed[det][ch].front()){
		this->detSUM[det][ch]-=this->detAdcValues[det][ch].front();
		this->detSUM2[det][ch]-=this->detAdcValues[det][ch].front()*this->detAdcValues[det][ch].front();
		this->detEventsInSum[det][ch]--;
	}
	if(detAdcValues[det][ch].back()<mean+sigma*maxSigma){
		this->detSUM[det][ch]+=this->detAdcValues[det][ch].back();
		this->detSUM2[det][ch]+=this->detAdcValues[det][ch].back()*this->detAdcValues[det][ch].back();
		this->detEventsInSum[det][ch]++;
		this->detEventUsed[det][ch].push_back(true);
	}
	else
		this->detEventUsed[det][ch].push_back(false);


	mean =this->detSUM[det][ch]/(float)this->detEventsInSum[det][ch];
	sigma=TMath::Sqrt(this->detSUM2[det][ch]/(float)this->detEventsInSum[det][ch]-mean*mean);

	return make_pair(mean,sigma);
}


pair<float,float> TPedestalCalculation::checkPedestalDia(int ch,int maxSigma){
	float mean =this->diaSUM[ch]/(float)this->diaEventsInSum[ch];
	float sigma=this->diaSUM2[ch]/(float)this->diaEventsInSum[ch]-mean*mean;
	//cout<<mean<<" "<<sigma<<" "<<this->diaAdcValues[ch].front()<<" "<<this->diaAdcValues[ch].back()<<" "<<diaEventUsed[ch].front()<<" "<<(diaAdcValues[ch].back()<mean+sigma*maxSigma)<<endl;
	if(this->diaEventUsed[ch].front()){
		this->diaSUM[ch]-=this->diaAdcValues[ch].front();
		this->diaSUM2[ch]-=this->diaAdcValues[ch].front()*this->diaAdcValues[ch].front();
		this->diaEventsInSum[ch]--;
	}
	if(diaAdcValues[ch].back()<mean+sigma*maxSigma){
		this->diaSUM[ch]+=this->diaAdcValues[ch].back();
		this->diaSUM2[ch]+=this->diaAdcValues[ch].back()*this->diaAdcValues[ch].back();
		this->diaEventsInSum[ch]++;
		this->diaEventUsed[ch].push_back(true);
	}
	else
		this->diaEventUsed[ch].push_back(false);


	mean =this->diaSUM[ch]/(float)this->diaEventsInSum[ch];
	sigma=TMath::Sqrt(this->diaSUM2[ch]/(float)this->diaEventsInSum[ch]-mean*mean);
	return make_pair(mean,sigma);

}



bool TPedestalCalculation::createPedestalTree(int nEvents)
{
	stringstream pedestalfilepath;
	pedestalfilepath<<sys->pwd();
	pedestalfilepath<<"/pedestalData."<<runNumber<<".root";
	cout<<"Try to open \""<<pedestalfilepath.str()<<"\""<<endl;
	pedestalFile=TFile::Open(pedestalfilepath.str().c_str());
	if(pedestalFile==NULL){
		cout<<"pedestalfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		pedestalFile= new TFile(pedestalfilepath.str().c_str(),"CREATE");
		pedestalFile->cd();
	}
	else{
		createdNewFile=false;
		cout<<"File exists"<<endl;
	}
	pedestalFile->cd();
	stringstream treeDescription;
	treeDescription<<"Pedestal Data of run "<<runNumber;
	pedestalFile->GetObject("pedestalTree",pedestalTree);
	if(pedestalTree!=NULL){
		cout<<"File and Tree Exists... \t"<<flush;
		if(pedestalTree->GetEntries()>=nEvents){
			createdNewTree=false;
			setBranchAdresses();
			cout<<"tree has enough entries...."<<endl;
			return false;
		}
		else{
			pedestalTree->Delete();
			pedestalTree=NULL;
		}
	}
	if(pedestalTree==NULL){
		pedestalFile->Close();
		pedestalFile=new TFile(pedestalfilepath.str().c_str(),"RECREATE");
		this->pedestalTree=new TTree("pedestalTree",treeDescription.str().c_str());
		createdNewTree=true;
		cout<<"there exists no tree:\'pedestalTree\"\tcreate new one."<<pedestalTree<<endl;
	}
	setBranchAdresses();
	return true;
}

void TPedestalCalculation::setBranchAdresses(){
	pedestalTree->Branch("PedestalMean",&pedestalMean,"PedestalMean[9][256]/F");
	pedestalTree->Branch("PedestalSigma",&pedestalSigma,"PedestaSigma[9][256]/F");
	pedestalTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	pedestalTree->Branch("runNumber",&runNumber,"runNumber/i");
}
