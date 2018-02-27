/*
 * PedestalCalculation.cpp
 *
 *  Created on: 10.11.2011
 *      Author: bachmair
 */

#include "../include/TPedestalCalculation.hh"


TPedestalCalculation::TPedestalCalculation(TSettings *newSettings){
	if(newSettings==0)exit(0);
	this->settings=newSettings;
	verbosity = settings->getVerbosity();
	if(verbosity)cout<<"**********************************************************"<<endl;
	cout<<"*****TPedestalCalculation::TPedestalCalculation***********"<<endl;
	if(verbosity)cout<<"**********************************************************"<<endl;

	slidingLength=settings->getPedestalSildingLength();//1000;//settings->getSl
	eventReader=NULL;
	pedestalTree=NULL;
	pedestalFile=NULL;
	this->runNumber=settings->getRunNumber();
	sys = gSystem;
	settings->goToPedestalTreeDir();
	eventReader=new TADCEventReader(settings->getRawTreeFilePath(),settings);
	if(verbosity)cout<<"TPedestalCalculation::TPedestalCalculation -> Set HistoSaver: "<<settings <<endl;
	histSaver = new HistogrammSaver(settings);
	histSaver->SetPlotsPath(settings->getToPedestalAnalysisDir());
	histSaver->SetRunNumber(settings->getRunNumber());
	if(verbosity)cout<<eventReader->GetEntries()<<endl;
	MAXSDETSIGMA=settings->getSi_Pedestal_Hit_Factor();
	MAXDIASIGMA=settings->getDi_Pedestal_Hit_Factor();
	if(verbosity)cout<<"Pedestal Hit Factor Silicon: "<<MAXSDETSIGMA<<"\nPedestal Hit Factor Diamond: "<<MAXDIASIGMA<<endl;
	hCommonModeNoise = new TH1F("hCommonModeNoise","hCommonModeNoise",512,-32,32);
	doCMNCorrection= true;//settings->doCommonModeNoiseCorrection();
	if(verbosity)cout<<"DO Common Mode Noise Correction: ";
	if(doCMNCorrection){
		if(verbosity)cout<<"TRUE "<<endl;}
	else if(verbosity)cout<<"FALSE"<<endl;
	//char t; cin >>t;//test
	printChannel=1;
	//settings->doCommonModeNoiseCorrection();
	for (UInt_t i = 0; i < 8*2; i ++)
	    cmn_sil[i] = 0;
}

TPedestalCalculation::~TPedestalCalculation() {
	// TODO Auto-generated destructor stub
	histSaver->SaveHistogram(hCommonModeNoise,true);
	delete histSaver;
	if(createdNewTree){
		pedestalFile->cd();
		pedestalTree->AddFriend("rawTree",settings->getRawTreeFilePath().c_str());
		cout<<pedestalTree->GetListOfFriends()->GetEntries()<<endl;
		pedestalFile->cd();
		pedestalTree->Write();
		pedestalTree->Delete();
	}

	delete eventReader;
	pedestalFile->Close();
	settings->goToOutputDir();
	if(verbosity)cout<<"Closing TPedestalCalculation\n\n\n"<<endl;
}

/**
 * This Gives me the Starting Values for the pedestal Calculation
 * One gets meanValues[det][ch & sigmaValues[det][ch], it also fills the deques
 * with adc data of slidingLength Events
 */
void TPedestalCalculation::calculateStartingPedestal(int nEvents){
	initialiseDeques();
	histSaver->SetNumberOfEvents(nEvents);
	if(verbosity)cout<<"TPedestalCalculation::calculatePedestals:"<<nEvents<<endl;
	//	nEvents = eventReader->GetEntries();
	createPedestalTree(nEvents);
	if(pedestalTree->GetEntries()>=nEvents){
		if(verbosity)cout<<"NO Sliding Pedestal Calculation needed, calculations already done."<<endl;
		return;
	}
	double meanSquared[9][256];
	if(verbosity)cout<<"initialise arrays"<<endl;
	for(int det=0;det<9;det++)
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			meanValues[det][ch]=0;
			meanSquared[det][ch]=0;
		}

	if(verbosity)cout<<"get mean and sigma"<<endl;
	/*
	 * calulate Pedestal mean and sigma with trick: ped=<x> and sigma=<x^2>-<x>^2
	 */
	for(UInt_t event=0;event<slidingLength;event++){
		eventReader->LoadEvent(event);
		for(UInt_t det=0;det <TPlaneProperties::getNDetectors();det++){
			for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
				meanValues[det][ch]+=(int)eventReader->getAdcValue(det,ch);
				meanSquared[det][ch]+=(int)eventReader->getAdcValue(det,ch)*(int)eventReader->getAdcValue(det,ch);
			}
		}
		for(int det=0;det <8;det++){
			for(int ch=0;ch<N_DET_CHANNELS;ch++){
				detAdcValues[det][ch].push_back(eventReader->getDet_ADC(det,ch));
			}
		}
		for(int ch=0;ch<N_DIA_CHANNELS;ch++){
			diaAdcValues[ch].push_back(eventReader->getDia_ADC(ch));
			diaAdcValuesCMN[ch].push_back(eventReader->getDia_ADC(ch));
		}
	}
	for(int det=0;det<9;det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			meanValues[det][ch]=meanValues[det][ch]/(double)slidingLength;
			meanSquared[det][ch]=meanSquared[det][ch]/(double)slidingLength;
			sigmaValues[det][ch]=meanSquared[det][ch]-meanValues[det][ch]*meanValues[det][ch];
			sigmaValues[det][ch]=TMath::Sqrt(sigmaValues[det][ch]);
		}
	if(verbosity)cout<<"DONE"<<endl;

}

void TPedestalCalculation::calculateSlidingPedestals(UInt_t nEvents){
	cout<<"calculate Sliding Pedestals ";
	if(doCMNCorrection)cout<<"with CMN Correction";
	cout<<endl;
	calculateStartingPedestal(nEvents);
	if(pedestalTree->GetEntries()>=nEvents){
		if(verbosity)cout<<"NO Sliding PEdestal Calculation needed, calculations already done."<<endl;
		return;
	}
	TStopwatch watch;
	watch.Start(true);

	//	initialise detAdcValues, diaAdcValues with values from rawTree
	//	initialiseDeques();
	calculateFirstPedestals(detAdcValues,diaAdcValues,MAXSDETSIGMA);
	fillFirstEventsAndMakeDiaDeque();

	//char t; cin>>t;
	//calculate sliding Pedestal Values for rest of Events and save them

	for(nEvent=slidingLength;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		//Add next Event to detAdcValues, diaAdcValues
		//Remove first Event from Queue

		eventReader->LoadEvent(nEvent);
		//SILICON PLANES
		updateSiliconPedestals();
		doCmNoiseCalculation();
		//DIAMOND PLANE
		updateDiamondPedestals();
		//		printDiamond(30);
		//calculateCurrentPedestals(detAdcValues,diaAdcValues);
		pedestalTree->Fill();
	}//end for

	watch.Stop();
	if(verbosity)cout<<"\nTime needed for PedestalCalulation:"<<endl;
	if(verbosity)watch.Print();
}


/**
 *
 */
void TPedestalCalculation::calculateFirstPedestals(deque<Int_t> DetAdcQueue[8][N_DET_CHANNELS], deque<Int_t> DiaAdcQueue[N_DIA_CHANNELS], int maxSigma){
	if(verbosity)cout<<"calculate Pedestal for the first "<<slidingLength<<" Entries..."<<endl;
	for(int det=0;det <8;det++){
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			TRawEventSaver::showStatusBar(256*det+ch,256*8,10);
			pair<Float_t,Float_t> values;
			values=this->calculateFirstPedestalDet(det,ch,DetAdcQueue[det][ch],meanValues[det][ch],sigmaValues[det][ch],7,MAXSDETSIGMA);//7 iteration for first pedestal
			pedestalMean[det][ch]=RoundFloat(values.first);
			pedestalSigma[det][ch]=RoundFloat(values.second);
		}
	}

	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		pair<Float_t,Float_t> values;
		values=this->calculateFirstPedestalDia(ch,DiaAdcQueue[ch],meanValues[8][ch],sigmaValues[8][ch],7,MAXDIASIGMA);//7 iterations for first pedestal
		diaPedestalMeanStartValues[ch]=RoundFloat(values.first);
		diaPedestalSigmaStartValues[ch]=RoundFloat(values.second);
	}
	if(verbosity)cout<<"\tDONE!"<<endl;
}


pair <Float_t,Float_t> TPedestalCalculation::calculateFirstPedestalDet(int det,int ch,deque<Int_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	detSUM[det][ch]=0;
	detSUM2[det][ch]=0;
	detEventsInSum[det][ch]=0;

	this->detEventUsed[det][ch].clear();
	for(nEvent=0;nEvent<adcQueue.size();nEvent++){
		Float_t adc = (Float_t)adcQueue.at(nEvent);
		//mean - maxSigma*sigma<=adc<= mean + maxSigma*sigma <-- Than it is no hit/seed
		if(  ( adc >= getLowLimitPedestal(meanChannel,sigmaChannel,maxSigma) )
				&&( adc <= getHighLimitPedestal(meanChannel,sigmaChannel,maxSigma) ) ){
			detEventUsed[det][ch].push_back(true);
			detSUM[det][ch]+=adc;
			detSUM2[det][ch]+=adc*adc;
			detEventsInSum[det][ch]++;
		}//end if
		else
			detEventUsed[det][ch].push_back(false);
	}//end for nEvent
	if(detEventsInSum[det][ch]<0.1*adcQueue.size())
		cout<<"For the calculation of first pedestals in det "<<det<<" and ch "<<ch<<"there were only "<<detEventsInSum[det][ch]<<" Events used..."<<endl;
	Float_t mean  = detSUM[det][ch]/(Float_t)detEventsInSum[det][ch];
	Float_t sigma = TMath::Sqrt(detSUM2[det][ch]/(Float_t)detEventsInSum[det][ch]-mean*mean);
	pair<Float_t,Float_t> output = make_pair(mean,sigma);
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDet(det,ch,adcQueue,mean,sigma,iterations-1,maxSigma);
}

pair <Float_t,Float_t> TPedestalCalculation::calculateFirstPedestalDia(int ch,deque<Int_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	diaSUM[ch]=0;
	diaSUM2[ch]=0;
	diaEventsInSum[ch]=0;
	this->diaEventUsed[ch].clear();
	for(UInt_t nEvent=0;nEvent<adcQueue.size();nEvent++){
		Float_t adc = (Float_t) adcQueue.at(nEvent);
		if(   (adc >= getLowLimitPedestal(meanChannel,sigmaChannel,maxSigma)  )
				&& (adc <= getHighLimitPedestal(meanChannel,sigmaChannel,maxSigma) ) ){
			diaEventUsed[ch].push_back(true);
			diaSUM[ch]+=adc;
			diaSUM2[ch]+=adc*adc;
			diaEventsInSum[ch]++;
		}//end if
		else
			diaEventUsed[ch].push_back(false);
	}//end for nEvent
	if(diaEventsInSum[ch]<0.1*adcQueue.size())
		cout<<"For the calculation of first pedestals in diamond and ch "<<ch<<"there were only "<<diaEventsInSum[ch]<<" Events used..."<<endl;
	Float_t mean=diaSUM[ch]/(Float_t)diaEventsInSum[ch];
	Float_t sigma=TMath::Sqrt( diaSUM2[ch]/(Float_t)diaEventsInSum[ch]-mean*mean);

	diaPedestalMean[ch]=RoundFloat(mean);
	diaPedestalSigma[ch]=RoundFloat(sigma);
	pair<Float_t,Float_t> output = make_pair(mean,sigma);
	//	if(verbosity>4&& ch ==103)
	//		cout<<"diamond ch "<<ch<<", it "<<iterations<<", usedEvents " <<diaEventsInSum[ch]<<"\t"<<mean<<"+/-"<<sigma<<endl;
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDia(ch,adcQueue,mean,sigma,iterations-1,maxSigma);
}

pair<Float_t, Float_t> TPedestalCalculation::calculateFirstPedestalDiaCMN(int ch, deque<Float_t> adcQueue, float meanCMN, float sigmaCMN, int iterations, float maxSigma) {
	cout<<"calculateFirstPedestalDiaCMN "<<ch<<" "<<meanCMN<< " "<< sigmaCMN<< " "<<iterations<< " "<<maxSigma<<endl;
    diaSUMCmn[ch]=0;
	diaSUM2Cmn[ch]=0;
	diaEventsInSumCMN[ch]=0;
	//  if(ch==7)cout<<"calcFirstPedCMN:"<<ch<<" "<<meanCMN<<" "<<sigmaCMN<<" "<<diaEventsInSumCMN[ch]<<endl;
	this->diaEventUsedCMN[ch].clear();
	for(nEvent=0;nEvent<adcQueue.size();nEvent++){
		Float_t adc = adcQueue.at(nEvent);
		Float_t lowLimit = getLowLimitPedestal(meanCMN,sigmaCMN,maxSigma);
		Float_t highLimit = getHighLimitPedestal(meanCMN,sigmaCMN,maxSigma);
		if(ch ==0) cout<< nEvent<<" "<<ch<<" "<<adc<<" "<<lowLimit<<" "<<highLimit<<endl;
		if(   (adc >= lowLimit)	&& (adc <= highLimit) ){
			diaEventUsedCMN[ch].push_back(true);
			diaSUMCmn[ch]+=adc;
			diaSUM2Cmn[ch]+=adc*adc;
			diaEventsInSumCMN[ch]++;
		}//end if
		else
			diaEventUsedCMN[ch].push_back(false);
	}//end for nEvent
	//TODO!!! FIX!!! PROBLEM!!!!!!!!
	///WORK HERE!!!!!!!
	if(diaEventsInSumCMN[ch]==0)
		cout<<"events in sum=0: "<<nEvent<<" "<<ch<<" "<<diaEventsInSumCMN[ch]<<" "<<diaEventsInSum[ch]<<endl;
	meanCMN=diaSUMCmn[ch]/(Float_t)diaEventsInSumCMN[ch];
	//  if(ch==7)cout<<diaSUMCmn[ch]<<" "<<diaSUM2Cmn[ch]<<" "<<diaEventsInSumCMN[ch]<<" "<<meanCMN<<endl;
	sigmaCMN=TMath::Sqrt( (diaSUM2Cmn[ch]/(Float_t)diaEventsInSumCMN[ch])-meanCMN*meanCMN);

	diaPedestalMeanCMN[ch]=RoundFloat(meanCMN);
	diaPedestalSigmaCMN[ch]=RoundFloat(sigmaCMN);
	pair<Float_t,Float_t> output = make_pair(meanCMN,sigmaCMN);
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDiaCMN(ch,adcQueue,meanCMN,sigmaCMN,iterations-1,maxSigma);
}

pair<Float_t,Float_t> TPedestalCalculation::checkPedestalDet(int det,int ch,int maxSigma){
	if(detEventUsed[det][ch].size()!=slidingLength)
		cout<<"detEventInUse has wrong length "<<detEventUsed[det][ch].size();
	if(detAdcValues[det][ch].size()!=slidingLength+1)
		cout<<"detAdcValues has wrong length... "<<detAdcValues[det][ch].size()<<endl;

	float mean =this->detSUM[det][ch]/(float)this->detEventsInSum[det][ch];
	float sigma=TMath::Sqrt(this->detSUM2[det][ch]/(float)this->detEventsInSum[det][ch]-mean*mean);

	//	if(det==0&&ch==5&&nEvent<3490&&nEvent>3450)
	//		cout<<"\r"<<nEvent<<"\t"<<mean<<" +/- "<<sigma<<"\t"<<detSUM[det][ch]<<"\t"<<detSUM2[det][ch]<<"\t"<<(int)detAdcValues[det][ch].back()<<"\t"<<((detAdcValues[det][ch].back()<mean+sigma*maxSigma))<<flush;

	//the sum is calculated  from events 0-slidingLength-1
	if(this->detEventUsed[det][ch].front()){
		//		if(det==0&&ch==5&&nEvent<3490&&nEvent>3450)cout<<"det in use remove"<<(int)detAdcValues[det][ch].front()<<" "<<flush;
		this->detSUM[det][ch]-=(ULong_t)this->detAdcValues[det][ch].front();
		this->detSUM2[det][ch]-=(ULong_t)this->detAdcValues[det][ch].front()*(ULong_t)this->detAdcValues[det][ch].front();
		this->detEventsInSum[det][ch]--;
	}
	//now the sum is calculated from events 1-slidingLength-1
	//todo make it readable
	if((float)detAdcValues[det][ch].back()<=mean+max(sigma*maxSigma,(float)1.)&&(float)detAdcValues[det][ch].back()>=mean-max(sigma*maxSigma,(float)1.)){
		//		if(det==0&&ch==5&&nEvent<3490&&nEvent>3450)cout<<"new pedestalEvent "<<(int)detAdcValues[det][ch].back()<<" "<<flush;
		this->detSUM[det][ch]+=(ULong_t)this->detAdcValues[det][ch].back();
		this->detSUM2[det][ch]+=(ULong_t)this->detAdcValues[det][ch].back()*(ULong_t)this->detAdcValues[det][ch].back();
		this->detEventsInSum[det][ch]++;
		this->detEventUsed[det][ch].push_back(true);
	}
	else
		this->detEventUsed[det][ch].push_back(false);
	//now the sum is calculated for events 1-slidingLength

	mean =this->detSUM[det][ch]/(float)this->detEventsInSum[det][ch];
	sigma=TMath::Sqrt(this->detSUM2[det][ch]/(float)this->detEventsInSum[det][ch]-mean*mean);
	//	if(det==0&&ch==5&&nEvent<3490&&nEvent>3450)cout<<mean<<" "<<sigma<<" "<<detEventsInSum[det][ch]<<endl;
	return make_pair(mean,sigma);
}


pair<float,float> TPedestalCalculation::checkPedestalDia(int ch,int maxSigma){
	float mean =this->diaSUM[ch]/(float)this->diaEventsInSum[ch];//ok
	float meanCMN= this->diaSUMCmn[ch]/(float)diaEventsInSumCMN[ch];//ok
	float sigma=TMath::Sqrt(this->diaSUM2[ch]/(float)this->diaEventsInSum[ch]-mean*mean);//ok
	float sigmaCMN=TMath::Sqrt(this->diaSUM2Cmn[ch]/(float)diaEventsInSumCMN[ch]-meanCMN*meanCMN);//ok
	//cout<<mean<<" "<<sigma<<" "<<this->diaAdcValues[ch].front()<<" "<<this->diaAdcValues[ch].back()<<" "<<diaEventUsed[ch].front()<<" "<<(diaAdcValues[ch].back()<mean+sigma*maxSigma)<<endl;
	//NORMAL CALCULATION WAY //ok
	if(this->diaEventUsed[ch].front()){
		this->diaSUM[ch]-=this->diaAdcValues[ch].front();
		this->diaSUM2[ch]-=this->diaAdcValues[ch].front()*this->diaAdcValues[ch].front();
		this->diaEventsInSum[ch]--;
	}
	if(diaAdcValues[ch].back()<=mean+max(sigma*maxSigma,(float)1.)&&diaAdcValues[ch].back()>=mean-max(sigma*maxSigma,(float)1.)){
		this->diaSUM[ch]+=this->diaAdcValues[ch].back();
		this->diaSUM2[ch]+=this->diaAdcValues[ch].back()*this->diaAdcValues[ch].back();
		this->diaEventsInSum[ch]++;
		this->diaEventUsed[ch].push_back(true);
	}
	else
		this->diaEventUsed[ch].push_back(false);

	//COMMON MODE NOISE CALCULATION WAY
	if(this->diaEventUsedCMN[ch].front()){
		this->diaSUMCmn[ch]-=this->diaAdcValuesCMN[ch].front();
		this->diaSUM2Cmn[ch]-=this->diaAdcValuesCMN[ch].front()*this->diaAdcValuesCMN[ch].front();
		this->diaEventsInSumCMN[ch]--;
	}
	if(diaAdcValuesCMN[ch].back()<=meanCMN+max(sigmaCMN*maxSigma,(float)1.)&&diaAdcValuesCMN[ch].back()>=meanCMN-max(sigmaCMN*maxSigma,(float)1.)){
		this->diaSUMCmn[ch]+=this->diaAdcValuesCMN[ch].back();
		this->diaSUM2Cmn[ch]+=this->diaAdcValuesCMN[ch].back()*this->diaAdcValuesCMN[ch].back();
		this->diaEventsInSumCMN[ch]++;
		this->diaEventUsedCMN[ch].push_back(true);
	}
	else
		this->diaEventUsedCMN[ch].push_back(false);

	mean =this->diaSUM[ch]/(float)this->diaEventsInSum[ch];//ok
	meanCMN = diaSUMCmn[ch]/(float)diaEventsInSumCMN[ch];//ok
	sigma=TMath::Sqrt(this->diaSUM2[ch]/(float)this->diaEventsInSum[ch]-mean*mean);//ok
	sigmaCMN=TMath::Sqrt(this->diaSUM2Cmn[ch]/(float)this->diaEventsInSumCMN[ch]-meanCMN*meanCMN);//ok
	diaPedestalMeanCMN[ch]=RoundFloat(meanCMN);
	diaPedestalSigmaCMN[ch]=RoundFloat(sigmaCMN);
	diaPedestalMean[ch]=RoundFloat(mean);
	diaPedestalSigma[ch]=RoundFloat(sigma);
	//  if(diaPedestalSigma[ch]<diaPedestalSigmaCMN[ch])
	//    cout<<std::setw(5)<<nEvent<<" "<<std::setw(3)<<ch<<" "<<setw(6)<<diaPedestalSigma[ch]<<" "<<setw(6)<<diaPedestalSigmaCMN[ch]<<endl;
	//	if(ch==7) cout<<cmNoise<<" mean: "<<mean<<"/"<<meanCMN<<"\tsigma:"<<sigma<<"/"<<sigmaCMN<<"\t"<<diaEventsInSum[ch]<<"/"<<diaEventsInSumCMN[ch]<<endl;
	return make_pair(mean,sigma);

}



bool TPedestalCalculation::createPedestalTree(int nEvents)
{
	stringstream pedestalfilepath;
	pedestalfilepath<<sys->pwd();
	pedestalfilepath<<"/pedestalData."<<runNumber<<".root";
	if(verbosity)cout<<"Try to open \""<<pedestalfilepath.str()<<"\""<<endl;
	pedestalFile=TFile::Open(pedestalfilepath.str().c_str());
	if(pedestalFile==NULL){
		if(verbosity)cout<<"pedestalfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		pedestalFile= new TFile(pedestalfilepath.str().c_str(),"CREATE");
		pedestalFile->cd();
	}
	else{
		createdNewFile=false;
		if(verbosity)cout<<"File exists"<<endl;
	}
	pedestalFile->cd();
	stringstream treeDescription;
	treeDescription<<"Pedestal Data of run "<<runNumber;
	pedestalFile->GetObject("pedestalTree",pedestalTree);
	if(pedestalTree!=NULL){
		if(verbosity)cout<<"File and Tree Exists... \t"<<pedestalTree->GetEntries()<<" of "<< nEvents<<flush;
		if(pedestalTree->GetEntries()>=nEvents){
			createdNewTree=false;
			setBranchAdresses();
			if(verbosity)cout<<"tree has enough entries...."<<endl;
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
		pedestalFile->cd();
		this->pedestalTree=new TTree("pedestalTree",treeDescription.str().c_str());
		createdNewTree=true;
		if(verbosity)cout<<"there exists no tree:\'pedestalTree\"\tcreate new one."<<pedestalTree<<endl;
	}
	setBranchAdresses();
	return true;
}


void TPedestalCalculation::doCmNoiseCalculation()
{
	cmNoise=0;

	UInt_t nCmNoiseEvents=0;
	Float_t maxVal = TPlaneProperties::getMaxSignalHeightDiamond();
	UInt_t det = TPlaneProperties::getDetDiamond();
	for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
		if(nEvent>this->diaAdcValues[ch].size()&&nEvent<slidingLength){
			cerr<<"diaADCValues["<<ch<<"].size() = "<<diaAdcValues[ch].size()<<" < "<<nEvent<<"  --> BREAK"<<endl;
			exit(-1);
		}

		Float_t adc = (nEvent<slidingLength)?this->diaAdcValues[ch].at(nEvent):eventReader->getDia_ADC(ch);
        bool masked = settings->IsMasked(det,ch);
		Float_t mean =  (nEvent<slidingLength)?diaPedestalMeanStartValues[ch]:diaPedestalMeanCMN[ch];
        Float_t signal = adc-mean;
		Float_t sigma = (nEvent<slidingLength)?diaPedestalSigmaStartValues[ch]:diaPedestalSigmaCMN[ch];
		if(sigma<=0) {
		    if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<ch<<" sigma <= 0 : "<<sigma<<endl;
		    continue;
		}
		Float_t snr = TMath::Abs(signal/sigma);

        if(snr!=snr||adc!=adc||signal!=signal)
            continue;

        if(adc>=maxVal || adc<0||signal>maxVal){
            if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<ch<<" invalid adc/signal: "<<adc<<"/"<<signal<<endl;
            continue;
        }
        if (TMath::Abs(snr)>settings->getCMN_cut()){
            if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<ch <<" snr over cut: "<<snr<<endl;
            continue;
        }
        if (masked){
            if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<det<<"/"<<ch <<" Is Masked "<<settings->isDet_channel_screened(det,ch)<<endl;
            continue;
        }

		if(verbosity>10||(verbosity>4&&nEvent==0))cout<<" "<<ch<<"\t"<<adc<<" "<<mean<< " "<<sigma<<" "<<signal<<" "<<snr<<endl;
		cmNoise+=signal;
		nCmNoiseEvents++;
	}
    cmNoise = (nCMNoiseEvents != 0)? cmNoise/(Float_t)nCmNoiseEvents: 0;	// DA: division by zero!!!!
	if(verbosity>4)cout<<nEvent <<" cmNoise: "<<" "<<cmNoise<<" "<<nCmNoiseEvents<<" "<<eventReader->getCmnCreated(8)<<endl;
	hCommonModeNoise->Fill(cmNoise,true);
}

void TPedestalCalculation::fillFirstEventsAndMakeDiaDeque()
{
	for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
		//		diaAdcValues[ch].clear();
		diaAdcValuesCMN[ch].clear();
	}
	//	//save Sliding Pedestal Values for first slidingLength Events

	for(nEvent=0;nEvent<slidingLength;nEvent++){
		//Fill tree
		//		eventReader->LoadEvent(nEvent);
		doCmNoiseCalculation();
		cmnValues.push_back(cmNoise);
        cout<<cmNoise<<endl;
		for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
			Float_t adc = (nEvent<slidingLength)?this->diaAdcValues[ch].at(nEvent):eventReader->getDia_ADC(ch);;
			adc -=cmNoise;
			diaAdcValuesCMN[ch].push_back(adc);
			Float_t mean = RoundFloat(diaPedestalMeanStartValues[ch]);
			Float_t sigma= RoundFloat(diaPedestalSigmaStartValues[ch]);

			diaPedestalMean[ch]= RoundFloat(mean);
			diaPedestalSigma[ch]= RoundFloat(sigma);
			mean-=cmNoise;
			//if(ch==7)cout<<nEvent<<" deque "<<adc<<" "<<diaAdcValues[ch].size()<<endl;
			diaPedestalMeanCMN[ch]= RoundFloat(mean);;
			diaPedestalSigmaCMN[ch]=RoundFloat(sigma);

		}
	}
	if(verbosity)cout<<"update first Pedestal Calculation"<<endl;
	for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
		pair<Float_t, Float_t> values = calculateFirstPedestalDia(ch,diaAdcValues[ch],diaPedestalMeanStartValues[ch],diaPedestalSigmaStartValues[ch],7,MAXDIASIGMA); // DA: corrected typo: it was diaPedestalMeanStartValues in the position where it received the sigmas.
		values = calculateFirstPedestalDiaCMN(ch,diaAdcValuesCMN[ch],diaPedestalMeanStartValues[ch],diaPedestalSigmaStartValues[ch],7,MAXDIASIGMA);
		diaPedestalMeanCMN[ch] = values.first;
		diaPedestalSigmaCMN[ch] = values.second;
		if(ch==7&&verbosity>4){
			//      cout<<"PEDESTAL: ch: "<<ch<<" "<<values.first<<" "<<values.second<<endl;
			for(UInt_t i=0;i<diaAdcValues[ch].size()&&i<diaAdcValuesCMN[ch].size();i++){
				cout<<" "<<setw(3)<<i<<"  "<<diaAdcValues[ch].at(i)<<" "<<diaEventUsed[ch].at(i)<<" "<<diaAdcValuesCMN[ch].at(i)<<" "<<diaEventUsedCMN[ch].at(i)<<" ";
				cout<<std::setw(5)<<(diaAdcValues[ch].at(i)-diaAdcValuesCMN[ch].at(i))<<" "<<cmNoise<<" "<<diaEventsInSum[ch]<<" "<<diaEventsInSumCMN[ch]<<endl;
			}
			if(verbosity%2==1){
			cout<<"press a key and enter to continue..."<<endl;
			char t; cin>>t;
			}
		}
	}
	for(nEvent = 0; nEvent<slidingLength;nEvent++){
		cmNoise = cmnValues.at(nEvent);
		for (UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
			diaPedestalMean[ch]= RoundFloat(diaPedestalMean[ch]);
			diaPedestalSigma[ch]= RoundFloat(diaPedestalSigma[ch]);
			diaPedestalMeanCMN[ch] =  RoundFloat(diaPedestalMeanCMN[ch]);
			diaPedestalSigmaCMN[ch] = RoundFloat(diaPedestalSigmaCMN[ch]);
		}
		printDiamond(30);
		this->pedestalTree->Fill();
	}
}

void TPedestalCalculation::initialiseDeques()
{
	//clear adcValues
	for(int det=0;det <8;det++){
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			detAdcValues[det][ch].clear();
		}
	}
	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		diaAdcValues[ch].clear();
		diaAdcValuesCMN[ch].clear();
	}


}


void TPedestalCalculation::updateSiliconPedestals(){
	for(int det=0;det <8;det++){
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			detAdcValues[det][ch].push_back(eventReader->getDet_ADC(det,ch));
			pair<float,float> values;
			values= checkPedestalDet(det,ch,MAXSDETSIGMA);
			pedestalMean[det][ch]=RoundFloat(values.first);
			pedestalSigma[det][ch]=RoundFloat(values.second);

			if(detEventUsed[det][ch].size()>0)detEventUsed[det][ch].pop_front();
			if(detAdcValues[det][ch].size())  detAdcValues[det][ch].pop_front();
		}
	}
}

void TPedestalCalculation::updateDiamondPedestals(){
	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		Float_t adcValue = (Float_t)eventReader->getDia_ADC(ch);
		diaAdcValues[ch].push_back(eventReader->getDia_ADC(ch));
		adcValue-=cmNoise;
		diaAdcValuesCMN[ch].push_back(adcValue);//eventReader->getDia_ADC(ch));

		pair<float,float> values;
		values = checkPedestalDia(ch,MAXDIASIGMA);
		//        diaPedestalMean[ch]=doCMNCorrection?diaPedestalMeanCMN[ch]:diaPedestalMean[ch];
		//        diaPedestalSigma[ch]=doCMNCorrection?diaPedestalSigmaCMN[ch]:diaPedestalSigma[ch];
		//        if(ch==7&&nEvent%10==0)
		//          cout<<nEvent<<": "<<ch<<" "<<pedestalMean[8][ch]<<" "<<pedestalSigma[8][ch]<<" "<<cmNoise<<"\t"<<diaAdcValues[ch].size()<<" "<<diaEventUsed[ch].size()<<" "<<diaEventsInSum[ch]<<endl;

		//if(ch==0&&(nEvent-slidingLength)%10000==0)cout<<nEvent<<". event, ch"<<ch<<"\t"<<values.first<<"+/-"<<values.second<<endl;

		if(diaEventUsed[ch].size())diaEventUsed[ch].pop_front();
		if(diaEventUsedCMN[ch].size())diaEventUsedCMN[ch].pop_front();
		if(diaAdcValues[ch].size())diaAdcValues[ch].pop_front();
		if(diaAdcValuesCMN[ch].size())diaAdcValuesCMN[ch].pop_front();
	}
}

void TPedestalCalculation::setBranchAdresses(){
	pedestalTree->Branch("PedestalMean",&pedestalMean,"PedestalMean[8][256]/F");
	pedestalTree->Branch("PedestalSigma",&pedestalSigma,"PedestaSigma[8][256]/F");
	pedestalTree->Branch("cmn_sil",&cmn_sil,"cmn_sil[8]/F");

	pedestalTree->Branch("diaPedestalMean",&diaPedestalMean,"diaPedestalMean[128]/F");
	pedestalTree->Branch("diaPedestalSigma",&diaPedestalSigma,"diaPedestaSigma[128]/F");

	pedestalTree->Branch("diaPedestalMeanCMN",&diaPedestalMeanCMN,"diaPedestalMeanCMN[128]/F");
	pedestalTree->Branch("diaPedestalSigmaCMN",&diaPedestalSigmaCMN,"diaPedestaSigmaCMN[128]/F");
	pedestalTree->Branch("commonModeNoise",&cmNoise,"commonModeNoise/F");

	pedestalTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	pedestalTree->Branch("runNumber",&runNumber,"runNumber/i");
	pedestalTree->Branch("cmnCorrection",&doCMNCorrection,"cmnCorrection/O");
}



void TPedestalCalculation::printDiamond(UInt_t nChannel){
	if(!verbosity)
		return;
	if (nChannel<TPlaneProperties::getNChannelsDiamond()&&printChannel!=0&&nEvent%printChannel==0){
		cout<<nEvent<<"\t"<<diaEventsInSum[nChannel]<<" "<<diaEventsInSumCMN[nChannel]<<flush;
		cout<<" "<<diaPedestalMean[nChannel]<<" "<<diaPedestalMeanCMN[nChannel]<<"\t"<<diaPedestalSigma[nChannel]<<" "<<diaPedestalSigmaCMN[nChannel]<<endl;
	}
}

Float_t TPedestalCalculation::getLowLimitPedestal(Float_t pedMean,Float_t pedSigma, Float_t maxSigma) {
	return pedMean - TMath::Max(pedSigma*maxSigma,(Float_t)1.0);
}

void TPedestalCalculation::calculateCommonModeDet(int det) {
}

Float_t TPedestalCalculation::GetCommonModeNoise(int det, int ch) {
    if (TPlaneProperties::isDiamondDetector(det))
        return cmNoise;
    int  i = ch>=(TPlaneProperties::getNChannels(det)/2)?1:0;
    return cmn_sil[det*2+i];
}

Float_t TPedestalCalculation::getHighLimitPedestal(Float_t pedMean,Float_t pedSigma, Float_t maxSigma) {
	return pedMean + TMath::Max(pedSigma*maxSigma,(Float_t)1.0);
}
