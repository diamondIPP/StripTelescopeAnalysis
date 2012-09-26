/*
 * PedestalCalculation.cpp
 *
 *  Created on: 10.11.2011
 *      Author: bachmair
 */

#include "../include/TPedestalCalculation.hh"


TPedestalCalculation::TPedestalCalculation(TSettings *settings){

  cout<<"**********************************************************"<<endl;
  cout<<"*****TPedestalCalculation::TPedestalCalculation***********"<<endl;
  cout<<"**********************************************************"<<endl;
	if(settings==0)exit(0);
	this->settings=settings;
		slidingLength=settings->getPedestalSildingLength();//1000;//settings->getSl
		eventReader=NULL;
		pedestalTree=NULL;
		pedestalFile=NULL;
		this->runNumber=settings->getRunNumber();
		sys = gSystem;
		settings->goToPedestalTreeDir();
		eventReader=new TADCEventReader(settings->getRawTreeFilePath(),runNumber);
		histSaver = new HistogrammSaver();
		histSaver->SetPlotsPath(settings->getToPedestalAnalysisDir());
		histSaver->SetRunNumber(settings->getRunNumber());
		cout<<eventReader->GetEntries()<<endl;
		MAXSDETSIGMA=settings->getSi_Pedestal_Hit_Factor();
		MAXDIASIGMA=settings->getDi_Pedestal_Hit_Factor();
		cout<<"Pedestal Hit Factor Silicon: "<<MAXSDETSIGMA<<"\nPedestal Hit Factor Diamond: "<<MAXDIASIGMA<<endl;
		hCommonModeNoise = new TH1F("hCommonModeNoise","hCommonModeNoise",512,-32,32);
		doCMNCorrection= settings->doCommonModeNoiseCorrection();
		cout<<"DO Common Mode Noise Correction: ";
		if(doCMNCorrection)
			cout<<"TRUE "<<endl;
		else cout<<"FALSE"<<endl;
		//char t; cin >>t;//test
		printChannel=1;
		//settings->doCommonModeNoiseCorrection();
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
	cout<<"Closing TPedestalCalculation\n\n\n"<<endl;
}


void TPedestalCalculation::calculatePedestals(int nEvents){

  histSaver->SetNumberOfEvents(nEvents);
	cout<<"TPedestalCalculation::calculatePedestals:"<<nEvents<<endl;
//	nEvents = eventReader->GetEntries();
	createPedestalTree(nEvents);
	if(pedestalTree->GetEntries()>=nEvents){
		cout<<"NO Sliding Pedestal Calculation needed, calculations already done."<<endl;
		return;
	}
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
	for(UInt_t event=0;event<slidingLength;event++){
		eventReader->LoadEvent(event);
		for(UInt_t det=0;det <TPlaneProperties::getNDetectors();det++){
			for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
				meanValues[det][ch]+=(int)eventReader->getAdcValue(det,ch);
				meanSquared[det][ch]+=(int)eventReader->getAdcValue(det,ch)*(int)eventReader->getAdcValue(det,ch);
			}
		}
	}
	for(int det=0;det<9;det++)
		for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
			meanValues[det][ch]=meanValues[det][ch]/(double)slidingLength;
			meanSquared[det][ch]=meanSquared[det][ch]/(double)slidingLength;
			sigmaValues[det][ch]=meanSquared[det][ch]-meanValues[det][ch]*meanValues[det][ch];
			sigmaValues[det][ch]=TMath::Sqrt(sigmaValues[det][ch]);
		}
	cout<<"DONE"<<endl;

}

void TPedestalCalculation::calculateSlidingPedestals(UInt_t nEvents){
	cout<<"calculate Sliding Pedestals"<<endl;
	if(doCMNCorrection)cout<<"DO CMN Correction"<<endl;
	if(pedestalTree->GetEntries()>=nEvents){
		cout<<"NO Sliding PEdestal Calculation needed, calculations already done."<<endl;
		return;
	}
	TStopwatch watch;
	watch.Start(true);

	//initialise detAdcValues, diaAdcValues with values from rawTree
	initialiseDeques();

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
	cout<<"\nStopWatch:"<<endl;
	watch.Print();
}

void TPedestalCalculation::calculateFirstPedestals(deque<UChar_t> DetAdcQueue[8][N_DET_CHANNELS], deque<Int_t> DiaAdcQueue[N_DIA_CHANNELS], int maxSigma){
	cout<<"calculate Pedestal for the first "<<slidingLength<<" Entries..."<<endl;
	for(int det=0;det <8;det++){
		for(int ch=0;ch<N_DET_CHANNELS;ch++){
			TRawEventSaver::showStatusBar(256*det+ch,256*8,10);
			pair<float,float> values;
			values=this->calculateFirstPedestalDet(det,ch,DetAdcQueue[det][ch],meanValues[det][ch],sigmaValues[det][ch],7,MAXSDETSIGMA);//7 iteration for first pedestal
			pedestalMean[det][ch]=RoundFloat(values.first);
			pedestalSigma[det][ch]=RoundFloat(values.second);
		}
	}

	for(int ch=0;ch<N_DIA_CHANNELS;ch++){
		pair<float,float> values;
		values=this->calculateFirstPedestalDia(ch,DiaAdcQueue[ch],meanValues[8][ch],sigmaValues[8][ch],7,MAXDIASIGMA);//7 iterations for first pedestal
		diaPedestalMeanStartValues[ch]=RoundFloat(values.first);
		diaPedestalSigmaStartValues[ch]=RoundFloat(values.second);
	}
	cout<<"\tDONE!"<<endl;
}


pair <float,float> TPedestalCalculation::calculateFirstPedestalDet(int det,int ch,deque<UChar_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	detSUM[det][ch]=0;
	detSUM2[det][ch]=0;
	detEventsInSum[det][ch]=0;

	this->detEventUsed[det][ch].clear();
	for(nEvent=0;nEvent<adcQueue.size();nEvent++){
		if(  ( (float)adcQueue.at(nEvent) >= meanChannel-(max(sigmaChannel*maxSigma,(float)1.)) )
		   &&( (float)adcQueue.at(nEvent) <= meanChannel+(max(sigmaChannel*maxSigma,(float)1.)) ) ){
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

pair <float,float> TPedestalCalculation::calculateFirstPedestalDia(int ch,deque<Int_t> adcQueue, float meanChannel, float sigmaChannel,int iterations,float maxSigma){
	diaSUM[ch]=0;
	diaSUM2[ch]=0;
	diaEventsInSum[ch]=0;
	this->diaEventUsed[ch].clear();
	for(UInt_t nEvent=0;nEvent<adcQueue.size();nEvent++){

		if(   ((float)adcQueue.at(nEvent) >= (meanChannel-max(sigmaChannel*maxSigma,(float)1.)) )
		   && ((float)adcQueue.at(nEvent) <= (meanChannel+max(sigmaChannel*maxSigma,(float)1.))) ){
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

    diaPedestalMean[ch]=RoundFloat(mean);
    diaPedestalSigma[ch]=RoundFloat(sigma);
	pair<float,float> output = make_pair(mean,sigma);
	if(iterations==0)return output;
	else return this->calculateFirstPedestalDia(ch,adcQueue,mean,sigma,iterations-1,maxSigma);
}

pair<float, float> TPedestalCalculation::calculateFirstPedestalDiaCMN(int ch, deque<Float_t> adcQueue, float meanCMN, float sigmaCMN, int iterations, float maxSigma) {
  diaSUMCmn[ch]=0;
  diaSUM2Cmn[ch]=0;
  diaEventsInSumCMN[ch]=0;
//  if(ch==7)cout<<"calcFirstPedCMN:"<<ch<<" "<<meanCMN<<" "<<sigmaCMN<<" "<<diaEventsInSumCMN[ch]<<endl;
  this->diaEventUsedCMN[ch].clear();
  for(nEvent=0;nEvent<adcQueue.size();nEvent++){
    if(   ((float)adcQueue.at(nEvent) >= (meanCMN-max(sigmaCMN*maxSigma,(float)1.)) )
       && ((float)adcQueue.at(nEvent) <= (meanCMN+max(sigmaCMN*maxSigma,(float)1.))) ){
      diaEventUsedCMN[ch].push_back(true);
      diaSUMCmn[ch]+=adcQueue.at(nEvent);
      diaSUM2Cmn[ch]+=(Float_t)((Float_t)adcQueue.at(nEvent)*(Float_t)adcQueue.at(nEvent));
      diaEventsInSumCMN[ch]++;
    }//end if
    else
      diaEventUsedCMN[ch].push_back(false);
  }//end for nEvent
  //TODO!!! FIX!!! PROBLEM!!!!!!!!
  ///WORK HERE!!!!!!!
  if(diaEventsInSumCMN[ch]==0)
    cout<<"events in sum=0: "<<nEvent<<" "<<ch<<" "<<diaEventsInSumCMN[ch]<<" "<<diaEventsInSum[ch]<<endl;
  meanCMN=(float)diaSUMCmn[ch]/(float)diaEventsInSumCMN[ch];
//  if(ch==7)cout<<diaSUMCmn[ch]<<" "<<diaSUM2Cmn[ch]<<" "<<diaEventsInSumCMN[ch]<<" "<<meanCMN<<endl;
  sigmaCMN=TMath::Sqrt( ((float)diaSUM2Cmn[ch]/(float)diaEventsInSumCMN[ch])-meanCMN*meanCMN);

  diaPedestalMeanCMN[ch]=RoundFloat(meanCMN);
//  cout<<meanCMN<<" "<<diaPedestalMeanCMN[ch];
  diaPedestalSigmaCMN[ch]=RoundFloat(sigmaCMN);
  pair<float,float> output = make_pair(meanCMN,sigmaCMN);
  if(iterations==0)return output;
  else return this->calculateFirstPedestalDiaCMN(ch,adcQueue,meanCMN,sigmaCMN,iterations-1,maxSigma);
}

pair<float,float> TPedestalCalculation::checkPedestalDet(int det,int ch,int maxSigma){
	if(detEventUsed[det][ch].size()!=slidingLength)
		cout<<"detEventInUse has wrong length"<<detEventUsed[det][ch].size();
	if(detAdcValues[det][ch].size()!=slidingLength+1)
		cout<<"detAdcValues has wrong length..."<<detAdcValues[det][ch].size()<<endl;

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
		cout<<"File and Tree Exists... \t"<<pedestalTree->GetEntries()<<" of "<< nEvents<<flush;
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
		pedestalFile->cd();
		this->pedestalTree=new TTree("pedestalTree",treeDescription.str().c_str());
		createdNewTree=true;
		cout<<"there exists no tree:\'pedestalTree\"\tcreate new one."<<pedestalTree<<endl;
	}
	setBranchAdresses();
	return true;
}


void TPedestalCalculation::doCmNoiseCalculation()
{
  cmNoise=0;

  UInt_t nCmNoiseEvents=0;
  Float_t maxVal = TPlaneProperties::getMaxSignalHeightDiamond();
  for(int ch=0;ch<N_DIA_CHANNELS;ch++){
    Float_t adc = eventReader->getDia_ADC(ch);
    Float_t mean =  (nEvent<slidingLength)?diaPedestalMeanStartValues[ch]:diaPedestalMeanCMN[ch];
    Float_t sigma = (nEvent<slidingLength)?diaPedestalSigmaStartValues[ch]:diaPedestalSigmaCMN[ch];
    Float_t signal = adc-mean;
    Float_t snr = (sigma==0)?(-1.):TMath::Abs(signal/sigma);
    if(snr<0||snr>settings->getDi_Pedestal_Hit_Factor()||adc>=maxVal||adc<0||signal>maxVal)//settings->isDet_channel_screened(TPlaneProperties::getDetDiamond(),ch))
      continue;
    if(snr!=snr||adc!=adc||signal!=signal)
      continue;
   // cout<<" "<<ch<<"\t"<<adc<<" "<<mean<< " "<<sigma<<" "<<signal<<" "<<snr<<endl;
    cmNoise+=signal;
    nCmNoiseEvents++;
  }
  cmNoise = cmNoise/(Float_t)nCmNoiseEvents;
  //`cout<<nEvent <<" cmNoise: "<<" "<<cmNoise<<" "<<nCmNoiseEvents<<" "<<endl;
  hCommonModeNoise->Fill(cmNoise,true);
}

void TPedestalCalculation::fillFirstEventsAndMakeDiaDeque()
{
  for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
    diaAdcValues[ch].clear();
    diaAdcValuesCMN[ch].clear();
  }
  //save Sliding Pedestal Values for first slidingLength Events

  for(nEvent=0;nEvent<slidingLength;nEvent++){
    //Fill tree
    eventReader->LoadEvent(nEvent);
    doCmNoiseCalculation();
    for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
      Float_t adc = eventReader->getDia_ADC(ch);
      diaAdcValues[ch].push_back(eventReader->getDia_ADC(ch));
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
    printDiamond(30);
    pedestalTree->Fill();
  }
  cout<<"update first Pedestal Calculation"<<endl;
  for(UInt_t ch=0;ch<N_DIA_CHANNELS;ch++){
    pair<Float_t, Float_t> values = calculateFirstPedestalDia(ch,diaAdcValues[ch],diaPedestalMeanStartValues[ch],diaPedestalMeanStartValues[ch],7,MAXDIASIGMA);
    calculateFirstPedestalDiaCMN(ch,diaAdcValuesCMN[ch],values.first,values.second,7,3);
//    if(ch==7){
////      cout<<"PEDESTAL: ch: "<<ch<<" "<<values.first<<" "<<values.second<<endl;
//      for(UInt_t i;i<diaAdcValues[ch].size()&&i<diaAdcValuesCMN[ch].size();i++){
//        cout<<" "<<setw(3)<<i<<"  "<<diaAdcValues[ch].at(i)<<" "<<diaEventUsed[ch].at(i)<<" "<<diaAdcValuesCMN[ch].at(i)<<" "<<diaEventUsedCMN[ch].at(i)<<" ";
//        cout<<std::setw(5)<<(diaAdcValues[ch].at(i)-diaAdcValuesCMN[ch].at(i))<<" "<<cmNoise<<" "<<diaEventsInSum[ch]<<" "<<diaEventsInSumCMN[ch]<<endl;
//      }
//      cout<<"DDSKLAS"<<endl;
//      char t; cin>>t;
//    }

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

  for(nEvent=0;nEvent<slidingLength;nEvent++){
  eventReader->LoadEvent(nEvent);
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
  if (nChannel<TPlaneProperties::getNChannelsDiamond()&&printChannel!=0&&nEvent%printChannel==0){
    cout<<nEvent<<"\t"<<diaEventsInSum[nChannel]<<" "<<diaEventsInSumCMN[nChannel]<<flush;
    cout<<" "<<diaPedestalMean[nChannel]<<" "<<diaPedestalMeanCMN[nChannel]<<"\t"<<diaPedestalSigma[nChannel]<<" "<<diaPedestalSigmaCMN[nChannel]<<endl;
  }	
}



