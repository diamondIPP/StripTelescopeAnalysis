//
//  TAnalysisOfPedestal.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 30.11.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TAnalysisOfPedestal.hh"
using namespace std;
TAnalysisOfPedestal::TAnalysisOfPedestal(TSettings* newSettings) {
	cout<<"**********************************************************"<<endl;
	cout<<"*********TAnalysisOfPedestal::TAnalysisOfPedestal*****"<<endl;
	cout<<"**********************************************************"<<endl;
	if(newSettings==0)
		exit(0);

	this->settings=newSettings;
	res=0;
	htmlPedestal= new THTMLPedestal(settings);
	sys = gSystem;
	UInt_t runNumber=settings->getRunNumber();
	sys->MakeDirectory(settings->getAbsoluteOuputPath().c_str());
	sys->cd(settings->getAbsoluteOuputPath().c_str());
	settings->goToPedestalTreeDir();

	eventReader=new TADCEventReader(settings->getPedestalTreeFilePath(),settings);//->getRunNumber());
	histSaver=new HistogrammSaver(settings);
	histSaver->SetOptStat("ormen");
	histSaver->SetOptFit(111);
	settings->goToPedestalAnalysisDir();
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	htmlPedestal->setFileGeneratingPath(sys->pwd());
	settings->goToPedestalTreeDir();
	initialiseHistos();
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
	this->invalidBiggestHitChannels.resize(TPlaneProperties::getNDetectors());
	verbosity = 0;

	adcValues.clear();
	pedestalValues.clear();
	upperHitCutValues.clear();
	lowerHitCutValues.clear();
}

TAnalysisOfPedestal::~TAnalysisOfPedestal() {
	cout<<"Invalid SNR for Biggest Channel can be found in "<<endl;
	for(UInt_t det = 0; det < invalidBiggestHitChannels.size(); det++){
		cout<<"Detector: "<<det<<", endries: "<<invalidBiggestHitChannels[det].size()<<"\t{"<<flush;
		std::set<UInt_t>::iterator it;
		for (it=invalidBiggestHitChannels[det].begin(); it!=invalidBiggestHitChannels[det].end(); ++it)
			cout<<setw(3)<<*it<<",";
		cout<<"}"<<endl;

	}
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
		TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		eventReader->LoadEvent(nEvent);
		analyseEvent();

	}
//	createPedestalMeanHistos();
	saveHistos();
}

void TAnalysisOfPedestal::analyseEvent(){
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){

        TString name = "hADCProfiles_"+(TString)TPlaneProperties::getStringForDetector(det);
        if(det==TPlaneProperties::getDetDiamond()){
            cmn = eventReader->getCMNoise(det,0);
            checkCommonModeNoise();
        }
        biggestSignal = -999;
        biggestHitChannel =-1;
        biggestSignalCMN = -999;
        biggestHitChannelCMN =-1;
        for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
            cmn = eventReader->getCMNoise(det,ch);
            adc = eventReader->getAdcValue(det,ch);
            if (hHistoMap.count(name))
                if(hHistoMap[name])
                    ((TProfile2D*)hHistoMap[name])->Fill(nEvent,ch,adc);
            isSaturated = eventReader->isSaturated(det,ch);

            pedestal = eventReader->getPedestalMean(det,ch,false);
            sigma = eventReader->getPedestalSigma(det,ch,false);
            noise = eventReader->getRawSignal(det,ch,false);//adc-pedestal;
            signal = eventReader->getSignal(det,ch, false);
            snr = eventReader->getSignalInSigma(det,ch,false);

            pedestalCMN = eventReader->getPedestalMean(det,ch,true);
            sigmaCMN = eventReader->getPedestalSigma(det,ch,true);
            noiseCMN =  eventReader->getRawSignal(det,ch,true);
            signalCMN = eventReader->getSignal(det,ch, true);

            if(det==TPlaneProperties::getDetDiamond()&&verbosity>8)
                cout<<nEvent<<" "<<det<<" "<<ch<<" "<<adc<<" "<<pedestal<<" "<<pedestalCMN<<" "<<signal<<" "<<signalCMN<<endl;
            checkForDeadChannels(det,ch);
            updateMeanCalulation(det,ch);
            findBiggestSignalInDet(det,ch);
            hNoiseDistributionCMN[det]->Fill(ch,noiseCMN);
            hNoiseDistribution[det]->Fill(ch,noise);
            if(det == TPlaneProperties::getDetDiamond() && ch == settings->getNoisePlotChannel()){
                Float_t hitCut = settings->getClusterHitFactor(det,ch);
                Float_t seedCut = settings->getClusterSeedFactor(det,ch);
                adcValues.push_back(adc);
                pedestalValues.push_back(pedestal);
                upperHitCutValues.push_back(pedestal + sigma * hitCut);
                lowerHitCutValues.push_back(pedestal - sigma * hitCut);
                upperSeedCutValues.push_back(pedestal + sigma * seedCut);
                lowerSeedCutValues.push_back(pedestal - sigma * seedCut);
                pedestalValuesCMN.push_back(pedestalCMN);
                upperHitCutValuesCMN.push_back(pedestalCMN + sigmaCMN * hitCut);
                lowerHitCutValuesCMN.push_back(pedestalCMN - sigmaCMN * hitCut);
                upperSeedCutValuesCMN.push_back(pedestalCMN + sigmaCMN * seedCut);
                lowerSeedCutValuesCMN.push_back(pedestalCMN - sigmaCMN * seedCut);
                eventNumbers.push_back(nEvent);
            }
        }
        analyseBiggestHit(det,false);
        analyseBiggestHit(det,true);
    }
    for(int i = 0; i< settings->diamondPattern.getNPatterns();i++)
        findBiggestSignalInDia(i+1);
}

void TAnalysisOfPedestal::checkForDeadChannels(UInt_t det,UInt_t ch)
{
	if(ch==0)numberOfSeeds=0;
	if(settings->isDet_channel_screened(det,ch))
		return;
	if(sigma==0){
		//cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
		return;
	};
	if(snr>settings->getClusterSeedFactor(det,ch)){
		hSeedMap[det]->Fill(ch);
		numberOfSeeds++;
	}
	if(ch==TPlaneProperties::getNChannels(det)-1)
		hNumberOfSeeds[det]->Fill(numberOfSeeds);

}


/**
 *
 */
Float_t TAnalysisOfPedestal::findYPlotRangeForPHHisto(TH1F *histo, Float_t hitCut)
{
	histo->Draw();
	Float_t max = 0;
	Int_t startBin = histo->FindBin(hitCut);
	//looking in area ch 0-- startbin for higehestBin
	for(Int_t binX=1;binX<startBin;binX++)
		if ( max < histo->GetBinContent(binX) )
			max = histo->GetBinContent(binX);
	Float_t oldMax = max;
	max=0;
	//look for highest Bin in startBin: Nbins
	for(Int_t binX=startBin;binX<histo->GetNbinsX();binX++)
		if ( max < histo->GetBinContent(binX) )
			max = histo->GetBinContent(binX);

	max=max*1.1;
	//check which bin used for axis range
	if (oldMax>1.5*max){
		histo->SetFillStyle(3244);
		histo->SetFillColor(kGray);//Gray
	}
	else
		if(max<oldMax)
			max=oldMax;

	//set Range of y Axis
	histo->GetYaxis()->SetRangeUser(0,max);
	return max;
}


void TAnalysisOfPedestal::findBiggestSignalInDet(UInt_t det,UInt_t ch){
	if(!settings->isDet_channel_screened(det,ch)){
		if (!isSaturated){
			if (signal > biggestSignal) {
				biggestHitChannel = ch;
				biggestSignal = signal;
			}
			if(signalCMN > biggestSignalCMN){
				biggestHitChannelCMN = ch;
				biggestSignalCMN = signalCMN;
			}
		}
		else{
			vecSaturatedChannels[det].push_back(ch);
			hSaturatedChannels[det]->Fill(ch);
		}
	}
}


void TAnalysisOfPedestal::findBiggestSignalInDia(UInt_t area) {
    UInt_t det = TPlaneProperties::getDetDiamond();
    if (area > settings->diamondPattern.getNPatterns()||area == 0)
        return;
    pair<int,int> pattern = settings->diamondPattern.getPatternChannels(area);
    Float_t maxSignal = -1e9;
    Float_t maxSignalCMN = -1e9;
    Float_t adjacentSignal = 0;
    Float_t adjacentSignalCMN = 0;
    if (pattern.second-pattern.first<3)
        return;
    int ch = pattern.first;
    Float_t leftSignal = 0;
    Float_t leftSignalCMN = 0;
    Float_t signal = eventReader->getSignal(det,ch, false);
    Float_t signalCMN = eventReader->getSignal(det,ch, true);
    Float_t rightSignal = eventReader->getSignal(det,ch+1, false);
    Float_t rightSignalCMN = eventReader->getSignal(det,ch+1, true);
    Int_t channel = -1;
    Int_t channelCMN = -1;
    Int_t adjacentChannel = -1;
    Int_t adjacentChannelCMN = -1;
    ch++;
    Float_t max =0;
//    for(int i = 0; i< TPlaneProperties::getNChannels(det);i++)
//        max = TMath::Max(max,eventReader->getSignal(det,i,false));
//    cout<<nEvent<<": "<<max<<endl;
    for(;ch<=pattern.second;ch++){
//        if(area==1)
//            cout<<"\t"<<area<<"-"<<ch<<"\t"<<signal<<" "<<maxSignal<<endl;
        if(signal>maxSignal){
            maxSignal = signal;
            adjacentSignal = TMath::Max(leftSignal,rightSignal);
            if(leftSignal==rightSignal)adjacentChannel=-1;
            else adjacentChannel = leftSignal>rightSignal?ch-1:ch+1;
            channel  = ch;
        }
        if (signalCMN > maxSignalCMN){
            adjacentSignalCMN = TMath::Max(leftSignalCMN,rightSignalCMN);
            if(leftSignalCMN==rightSignalCMN)adjacentChannelCMN=-1;
            else adjacentChannelCMN = leftSignalCMN>rightSignalCMN?ch-1:ch+1;
            maxSignalCMN = signalCMN;
            channelCMN = ch;
        }
        signal = rightSignal;
        rightSignal = eventReader->getSignal(det,ch+1,false);
        signalCMN = rightSignalCMN;
        rightSignalCMN = eventReader->getSignal(det,ch+1,true);
    }
    Float_t snr = maxSignal / eventReader->getPedestalSigma(det,channel,false);
    Float_t snrCMN = maxSignalCMN / eventReader->getPedestalSigma(det,channelCMN,true);
//    if(area ==1 && maxSignal>0){
//    cout<<TString::Format("%d-%7d %6.1f/%5.1f \t %6.1f/%5.1f ",area,nEvent,maxSignal,adjacentSignal,maxSignalCMN,adjacentSignalCMN)<<"";
//    if(maxSignal>100) cout<<"\t*****";
//    cout<<endl;}
//    cout<<"\t"<<channel<<" "<<eventReader->getSignal(det,channel-1, false)<<"-"
//                            <<eventReader->getSignal(det,channel, false)<<"-"
//                            <<eventReader->getSignal(det,channel+1, false)<<" --> "<<snr<<endl;
//    cout<<"\t"<<channelCMN<<" "<<eventReader->getSignal(det,channelCMN-1, true)<<"-"
//                                <<eventReader->getSignal(det,channelCMN, true)<<"-"
//                                <<eventReader->getSignal(det,channelCMN+1, true)<<" --> "<<snrCMN<<endl;
    if (hBiggestSignalInSigmaDiaPattern.count(area)){
//        cout<<hBiggestSignalInSigmaDiaPattern[area]->GetName()<<endl;
        if(snr>3){
            hBiggestSignalInSigmaDiaPattern[area]->Fill(snr);
            if(adjacentChannel!=-1)
                hBiggestAdjacentSignalInSigmaDiaPattern[area]->Fill(adjacentSignal/eventReader->getPedestalSigma(det,adjacentChannel,false));
        }
    }
    else
        cout<<"hBiggestSignalInSigmaDiaPattern "<<area<<" doesn't exist"<<endl;
    if (hBiggestSignalInSigmaDiaPatternCMN.count(area)){
        if (snrCMN>3){
            hBiggestSignalInSigmaDiaPatternCMN[area]->Fill(snrCMN);
            if(adjacentChannelCMN != -1)
                hBiggestAdjacentSignalInSigmaDiaPatternCMN[area]->Fill(adjacentSignalCMN/eventReader->getPedestalSigma(det,adjacentChannelCMN,true));
        }
    }
    else
        cout<<"hBiggestSignalInSigmaDiaPatternCMN "<<area<<" doesn't exist"<<endl;
}


/**
 * create vector with biggest and biggest adjacent hit with PH in Sigma
 */
void TAnalysisOfPedestal::analyseBiggestHit(UInt_t det,bool CMN_corrected) {
	if (CMN_corrected&&TPlaneProperties::getDetDiamond()!=det)
		return;
	Int_t bigHit;
	Float_t bigSignal;
	if(!CMN_corrected){
		bigHit = biggestHitChannel;
		bigSignal=biggestSignal;
	}
	else{
		bigHit = biggestHitChannelCMN;
		bigSignal=biggestSignalCMN;
	}
	if(bigSignal<0)
		return;
	if(det==TPlaneProperties::getDetDiamond()&&verbosity>8)
		cout<<"Biggest Signal: "<<bigSignal<<"@"<<bigHit<<"-"<<CMN_corrected<<"\t"<<biggestSignal<<"@"<<biggestHitChannel<<"\t"<<biggestSignalCMN<<"@"<<biggestHitChannelCMN<<endl;

	Int_t leftCh = bigHit-1;
	Float_t leftSignal = 0;
	if(settings->isDet_channel_screened(det,bigHit)&&eventReader->isSaturated(det,bigHit))
		cout<<"Error1:"<<det<<"_"<<setw(5)<<nEvent<<"\t"<<settings->isDet_channel_screened(det,bigHit)<<eventReader->isSaturated(det,bigHit)<<endl;
	if(!settings->isDet_channel_screened(det,leftCh)&&!eventReader->isSaturated(det,leftCh))
		leftSignal = eventReader->getSignal(det,leftCh, CMN_corrected);
	Int_t rightCh = bigHit+1;
	Float_t rightSignal = 0;
	if(!settings->isDet_channel_screened(det,rightCh)&&!eventReader->isSaturated(det,rightCh))
		rightSignal=eventReader->getSignal(det,rightCh,CMN_corrected);

	Float_t biggestAdjacentSignal=0;
	Int_t biggestAdjacentHitChannel = -9999;

	if(leftSignal>rightSignal){
		biggestAdjacentSignal=leftSignal;
		biggestAdjacentHitChannel = bigHit-1;
	}
	else if (rightSignal>leftSignal){
		biggestAdjacentSignal = rightSignal;
		biggestAdjacentHitChannel = bigHit+1;
	}
	if(leftSignal>bigSignal||rightSignal>bigSignal)
		cout<<"Error2:"<<det<<"_"<<setw(5)<<nEvent<<"\t"<<leftSignal<<"-"<<bigSignal<<"-"<<rightSignal<<"\t"<<leftCh<<"_"<<bigHit<<"_"<<rightCh<<endl;

	Int_t hitOrder = biggestAdjacentHitChannel!=-9999?biggestAdjacentHitChannel-bigHit:0;

	if(biggestAdjacentHitChannel!=-9999&&(biggestAdjacentHitChannel<0||biggestAdjacentHitChannel>=(Int_t)TPlaneProperties::getNChannels(det)))
		cout<<"something is wrong: biggestAdjacentHitChannel: "<<biggestAdjacentHitChannel<<" BiggestHitChannel:"<<bigHit<<" "<<leftSignal<<" "<<rightSignal<<endl;

	Float_t biggestSignalInSigma = eventReader->getSignalInSigma(det,bigHit,CMN_corrected);
	Float_t biggestAdjacentSignalInSigma =biggestAdjacentHitChannel>=0&&biggestAdjacentHitChannel< (Int_t)TPlaneProperties::getNChannels(det)? eventReader->getSignalInSigma(det,biggestAdjacentHitChannel,CMN_corrected):0;
    Float_t S_L,S_R,eta,Q;
    Float_t SNR = biggestSignalInSigma;

    if (hitOrder == -1){
        S_L = biggestAdjacentSignal;
        S_R = bigSignal;
    }
    else if (hitOrder == 1){
        S_L = bigSignal;
        S_R = biggestAdjacentSignal;
    }
    if (hitOrder != 0){
        Q = S_R+S_L;
        eta = S_R/(Q);
        hLeftVsRightSignal[det]->Fill(S_L,S_R);
        hEtaVsCharge[det]->Fill(eta,Q);
        hEtaVsSNR[det]->Fill(eta,SNR);
    }

    if(!CMN_corrected){
        vecBiggestHitChannel[det].push_back(biggestHitChannel);
        vecBiggestSignalInSigma[det].push_back(biggestSignalInSigma);
        //		if(TPlaneProperties::isDiamondDetector(det)){
        //		    Int_t pattern = settings->diamondPattern.getPatternOfChannel(biggestHitChannel);
        //		    if (pattern!=-1)
        //		        if (hBiggestSignalInSigmaDiaPattern.count(pattern)>0)
        //		            hBiggestSignalInSigmaDiaPattern[pattern]->Fill(biggestSignalInSigma);
        //		}
        vecBiggestSignal[det].push_back(bigSignal);
        vecBiggestAdjacentSignal[det].push_back(biggestAdjacentSignal);
        if(biggestAdjacentHitChannel!=-9999){
            vecBiggestAdjacentSignalInSigma[det].push_back(biggestAdjacentSignalInSigma);
            vecBiggestAdjacentSignal[det].push_back(biggestAdjacentSignal);
            vecBiggestAdjacentHitChannel[det].push_back(biggestAdjacentHitChannel);
            vecHitOrder[det].push_back(biggestAdjacentHitChannel-bigHit);
            hBiggestAdjacentSignalInSigma[det]->Fill(biggestAdjacentSignalInSigma);
        }
        else{
            vecHitOrder[det].push_back(0);
            vecBiggestAdjacentSignalInSigma[det].push_back(0);
            vecBiggestAdjacentHitChannel[det].push_back(-9999);
        }
        if(vecBiggestAdjacentSignal[det].back()>vecBiggestSignal[det].back()&&false)
            cout<<setw(5)<<nEvent<<" "<<det<<" "<<CMN_corrected<<" "<<setw(3)<<vecBiggestHitChannel[det].back()<<":"<<vecBiggestSignal[det].back()<<":"<<vecBiggestSignalInSigma[det].back()
            <<"\t"<<setw(3)<<vecBiggestAdjacentHitChannel[det].back()<<":"<<vecBiggestAdjacentSignal[det].back()<<":"<<vecBiggestAdjacentSignalInSigma[det].back()<<endl;
        hHitOrderMap[det]->Fill(hitOrder);
        hBiggestHitChannelMap[det]->Fill(bigHit);
    }
    else if (TPlaneProperties::isDiamondDetector(det)){
        vecBiggestHitChannelCMN[det].push_back(bigHit);
        vecBiggestSignalInSigmaCMN[det].push_back(biggestSignalInSigma);
        if(biggestAdjacentHitChannel!=-9999){
            vecBiggestAdjacentSignalInSigmaCMN[det].push_back(biggestAdjacentSignalInSigma);
            vecBiggestAdjacentHitChannelCMN[det].push_back(biggestAdjacentHitChannel);
            //				vecHitOrder[det].push_back(biggestAdjacentHitChannel-bigHit);
            hBiggestAdjacentSignalInSigmaCMN[det]->Fill(biggestAdjacentSignalInSigma);
        }
        else{
            //				vecHitOrder[det].push_back(0);
            vecBiggestAdjacentSignalInSigmaCMN[det].push_back(0);
            vecBiggestAdjacentHitChannelCMN[det].push_back(-9999);
        }
        //			hHitOrderMap[det]->Fill(hitOrder);
        hBiggestHitChannelMapCMN[det]->Fill(bigHit);
    }
}


void TAnalysisOfPedestal::initialiseHistos()
{

    hRelCmnUncertainty = new TH1F("hRelCmnUncertainty","hRelCmnUncertainty",512,-200,+200);
    Int_t nChannels = TPlaneProperties::getNChannelsDiamond();
    hCmnNUsedChannels = new TH1F("hCmnNUsedChannels","hCmnNUsedChannels",nChannels,0,nChannels-1);
    hCmnUsedChannels = new TH1F("hCmnUsedChannels","hCmnUsedChannels",nChannels,0,nChannels-1);
    hNewComonModeNoise = new TH1F("hNewComonModeNoise","hNewComonModeNoise",512,-30,30);
    hCmnChannelWeight = new TH1F("hCmnChannelWeight","hCmnChannelWeight",512,-1000,+1000);

    hCmnChannelWeightVsChannel = new TH2F("hCmnChannelWeightVsChannel","hCmnChannelWeightVsChannel",nChannels,0,nChannels-1,512,-1000,+1000);
    hCmnFractionVsChannel = new TH2F("hCmnFractionVsChannel","hCmnFractionVsChannel",nChannels,0,nChannels-1,512,-100,+100);
    hCmnNewVsNUsedChannels = new TH2F("hCmnNewVsNUsedChannels","hCmnNewVsNUsedChannels",512,-30,30,nChannels,0,nChannels-1);
    hCmnVsNewCmn = new TH2F("hCmnVsNewCmn","hCmnVsNewCmn",512,-30,30,512,-30,30);
    Float_t xmax = settings->getNEvents();
    Int_t nbins = xmax/1e4;
    hNewCmnVsEventNo = new TH2F("hNewCmnVsEventNo","hNewCmnVsEventNo",nbins*10,0,xmax,512,-30,30);
    for(UInt_t i=0;i < settings->diamondPattern.getNIntervals();i++){
        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        TString name = TString::Format("hBiggestSignalInSigmaDiaPattern_%d_ch_%d_%d",i+1,channels.first,channels.second);
        TH1F* histo = new TH1F(name,name,196*2,4,200);
        histo->GetXaxis()->SetTitle("biggest hit in sigma");
        histo->GetYaxis()->SetTitle("number of entries #");
        hBiggestSignalInSigmaDiaPattern[i+1] = histo;

        name = TString::Format("hBiggestSignalInSigmaDiaPatternCMN_%d",i+1);
        histo = new TH1F(name,name,196*2,4,200);
        histo->GetXaxis()->SetTitle("biggest hit in sigma cmn");
        histo->GetYaxis()->SetTitle("number of entries #");
        hBiggestSignalInSigmaDiaPatternCMN[i+1] = histo;

        name = TString::Format("hBiggestAdjacentSignalInSigmaDiaPattern_%d",i+1);
        histo = new TH1F(name,name,196*2,4,200);
        histo->GetXaxis()->SetTitle("biggest hit in sigma");
        histo->GetYaxis()->SetTitle("number of entries #");
        hBiggestAdjacentSignalInSigmaDiaPattern[i+1] = histo;

        name = TString::Format("hBiggestAdjacentSignalInSigmaDiaPatternCMN_%d",i+1);
        histo = new TH1F(name,name,196*2,4,200);
        histo->GetXaxis()->SetTitle("biggest adjacent hit in sigma cmn");
        histo->GetYaxis()->SetTitle("number of entries #");
        hBiggestAdjacentSignalInSigmaDiaPatternCMN[i+1] = histo;

    }
    map<Int_t,TH1F*>::iterator it;
    for(it=hBiggestSignalInSigmaDiaPattern.begin();it!=hBiggestSignalInSigmaDiaPattern.end();it++)
        cout<<(*it).first<<": "<<(*it).second->GetName()<<endl;
    hCMNoiseDistribution= new TH1F("hCMNoiseDistribution","hCMNoiseDistribution",512,-20,20);
    hCMNoiseDistribution->GetXaxis()->SetTitle("Common Mode Noise [ADC]");
    hCMNoiseDistribution->GetYaxis()->SetTitle("number of entries [#]");
    Int_t nEvents = eventReader->GetEntries();
    Int_t nEventBins = (nEvents)/10000;
    hCMNoiseDistributionEventNo = new TH2F("hCMNoiseDistributionEventNo", "hCMNoiseDistributionEventNo", nEventBins, 0., nEvents, 512, -20, 20);
    hCMNoiseDistributionEventNo->GetXaxis()->SetTitle("Event number");
    hCMNoiseDistributionEventNo->GetYaxis()->SetTitle("Common Mode Noise [ADC]");
    hCMNoiseDistributionEventNo->GetZaxis()->SetTitle("number of entries [#]");
    for (UInt_t det =0;det<TPlaneProperties::getNDetectors();det++){
        TString name = "hADCProfiles_"+(TString)TPlaneProperties::getStringForDetector(det);
        Int_t ybins = TPlaneProperties::getNChannels(det);
        Float_t yend = ybins -1;
        Float_t ystart = 0;
        Float_t xstart =0;
        Float_t xend = settings->getNEvents();
        Int_t xbins = xend/1000;
        hHistoMap[name] = new TProfile2D(name,name,xbins,xstart,xend,ybins,ystart,yend);
        stringstream histoName,histoTitle,xTitle,yTitle;
        histoName<<"hNoiseDistributionOfAllNonHitChannels_"<<TPlaneProperties::getStringForDetector(det);
        histoTitle<<"Noise Distribution of all non hit channels in Plane"<<TPlaneProperties::getStringForDetector(det);
        xTitle<<"non hit Noise (Adc-Ped.) in ADC counts";
        yTitle<<"Number of Entries #";
        Float_t width = settings->getNoise_si_max     ();
        UInt_t nBins  = settings->getNoise_si_num_bins();
        if (det==TPlaneProperties::getDetDiamond()){
            width = settings->getNoise_di_max     ();
            nBins = settings->getNoise_di_num_bins();
        }
        if(TPlaneProperties::isSiliconDetector(det)){
            hAllAdcNoise[det]= new TH1F(histoName.str().c_str(),histoTitle.str().c_str(),nBins,(-1)*width,width);
            hAllAdcNoise[det]->GetXaxis()->SetTitle(xTitle.str().c_str());
            hAllAdcNoise[det]->GetYaxis()->SetTitle(yTitle.str().c_str());
        }
        else{
        	TString hName= histoName.str();
        	TString hTitle= histoTitle.str();
            hDiaAllAdcNoise= new TH1F(hName,hTitle,nBins,(-1)*width,width);
            hDiaAllAdcNoise->GetXaxis()->SetTitle(xTitle.str().c_str());
            hDiaAllAdcNoise->GetYaxis()->SetTitle(yTitle.str().c_str());

            hDiaAllAdcNoiseChannel = new TH2F(hName+"Channel",hTitle+" Channel",
            	TPlaneProperties::getNChannels(det),0,TPlaneProperties::getNChannels(det),
            	nBins,(-1)*width,width);
            hDiaAllAdcNoiseChannel->GetYaxis()->SetTitle(xTitle.str().c_str());
            hDiaAllAdcNoiseChannel->GetZaxis()->SetTitle(yTitle.str().c_str());
            hDiaAllAdcNoiseChannel->GetYaxis()->SetTitle("Channel no.");

			hDiaAllAdcNoiseEventNo = new TH2F(hName + "_EventNo", hTitle + " Event Number", nEventBins, 0., nEvents, nBins, (-1)*width, width);
			hDiaAllAdcNoiseEventNo->GetXaxis()->SetTitle("Event number");
			hDiaAllAdcNoiseEventNo->GetYaxis()->SetTitle(xTitle.str().c_str());
			hDiaAllAdcNoiseEventNo->GetZaxis()->SetTitle(yTitle.str().c_str());

            TString extName = "_CMNcorrected";
            TString extTitle = " CMN corrected";
            hDiaAllAdcNoiseCMN = new TH1F(hName+extName,
            		hTitle+"_CMNcorrected",nBins,(-1)*width,width);
            xTitle<<" -  CMN corrected";
            hDiaAllAdcNoiseCMN->GetXaxis()->SetTitle(xTitle.str().c_str());
            hDiaAllAdcNoiseCMN->GetYaxis()->SetTitle(yTitle.str().c_str());

            hDiaAllAdcNoiseCMNChannel = new TH2F(hName+"Channel"+extName,hTitle+" Channel"+extTitle,
                        	TPlaneProperties::getNChannels(det),0,TPlaneProperties::getNChannels(det),
                        	nBins,(-1)*width,width);
            hDiaAllAdcNoiseCMNChannel->GetYaxis()->SetTitle("Channel no.");
            hDiaAllAdcNoiseCMNChannel->GetYaxis()->SetTitle(xTitle.str().c_str());
            hDiaAllAdcNoiseCMNChannel->GetZaxis()->SetTitle(yTitle.str().c_str());

			hDiaAllAdcNoiseCMNEventNo = new TH2F(hName + "_EventNo" + extName, hTitle + " Event Number" + extTitle, nEventBins, 0., nEvents, nBins, (-1)*width, width);
			hDiaAllAdcNoiseCMNEventNo->GetXaxis()->SetTitle("Event number");
			hDiaAllAdcNoiseCMNEventNo->GetYaxis()->SetTitle(xTitle.str().c_str());
			hDiaAllAdcNoiseCMNEventNo->GetZaxis()->SetTitle(yTitle.str().c_str());
        }

        Int_t nbins = 512;
        Float_t xlow;
        Float_t xup;
        if (TPlaneProperties::isSiliconDetector(det)){
            xlow = -10;
            xup = 10;
        }
        else{
            xlow = -100;
            xup = 100;
        }


        name ="hNoiseDistributionPerChannelCMN_"+TPlaneProperties::getStringForDetector(det);
        TString title ="Noise Dist per Channel CMN "+TPlaneProperties::getStringForDetector(det)+"ch no;adc-ped-cmn";
        hNoiseDistributionCMN[det] = new TH2F(name,title,
                TPlaneProperties::getNChannels(det),0,TPlaneProperties::getNChannels(det),xbins,xlow,xup);

        name ="hNoiseDistributionPerChannel_"+TPlaneProperties::getStringForDetector(det);
        title ="Noise Dist per Channel"+TPlaneProperties::getStringForDetector(det)+"ch no;adc-ped";
        hNoiseDistribution[det] = new TH2F(name,title,
                TPlaneProperties::getNChannels(det),0,TPlaneProperties::getNChannels(det),xbins,xlow,xup);

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
        histoName<<"hPulseHeight__BiggestSignalChannelInSigma_"<<TPlaneProperties::getStringForDetector(det);
        hSNR_BiggestSignal[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
    }
    for (int det=0;det<9;det++){//todo why such a big histo?so big?
        stringstream histoName;
        histoName<<"hPulseHeightInSigma_BiggestSignalAdjacentToBiggestSignal"<<TPlaneProperties::getStringForDetector(det);
        hSNR_BiggestAdjacent[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
    }
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"hChannel_BiggestSignal_"<<TPlaneProperties::getStringForDetector(det);
        UInt_t nChannels=TPlaneProperties::getNChannels(det);
        hChannelBiggestSignal[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),nChannels,0,nChannels);
    }

    for (UInt_t det = 0; det < 9; det++) {
        int nbins = settings->getSnr_distribution_num_bins();
        Float_t min = 0.;
//        Float_t max = 64.;
//        if(det==TPlaneProperties::getDetDiamond()){max=128;nbins=512;}
		Float_t max;
		if (TPlaneProperties::isDiamondDetector(det)) max = settings->getSnr_distribution_di_max();
		else                                          max = settings->getSnr_distribution_si_max();

        stringstream histoName;

        histoName.str("");
        histoName << "hPulseHeight_BiggestAdjacentInSigma_" << TPlaneProperties::getStringForDetector(det);
        hBiggestAdjacentSignalInSigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);
        histoName << "_CMNcorrected";
        hBiggestAdjacentSignalInSigmaCMN[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << "hSecondBiggestHitMinusBiggestHitPosition_" << TPlaneProperties::getStringForDetector(det);
        hHitOrderMap[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),3,-1.5,1.5);

        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaLeft" << TPlaneProperties::getStringForDetector(det);
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << "hPulseHeightSecondBiggestHitChannelInSigmaRight" << TPlaneProperties::getStringForDetector(det) ;
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        TString hName = (TString)"hBiggestHitMap"+ (TString)TPlaneProperties::getStringForDetector(det);
        hBiggestHitChannelMap[det] = new TH1F(hName,hName,TPlaneProperties::getNChannels(det),-.5,TPlaneProperties::getNChannels(det)-.5);
        hName += (TString)"_CMNcorrected";
        hBiggestHitChannelMapCMN[det] = new TH1F(hName,hName,TPlaneProperties::getNChannels(det),-.5,TPlaneProperties::getNChannels(det)-.5);

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

        histoName.str("");
        //hLeftVsRightSignal
        max = TPlaneProperties::getMaxSignalHeight(det);
        if (max >1000)
            nbins = max/4;
        else
            nbins = max;
        histoName << "hLeftVsRightSignal_" << TPlaneProperties::getStringForDetector(det);
        hLeftVsRightSignal[det] = new TH2F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max,nbins,min,max);
        histoName << "_CMNcorrected";
        hLeftVsRightSignalCMN[det] = new TH2F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max,nbins,min,max);
        histoName.str("");
        histoName <<"hEtaVsCharge_" << TPlaneProperties::getStringForDetector(det);
        hEtaVsCharge[det] = new TH2F(histoName.str().c_str(),histoName.str().c_str(),nbins,0,1,nbins,min,max);
        hEtaVsCharge[det]->GetXaxis()->SetTitle("#eta = S_R/(S_L+S_R)");
        hEtaVsCharge[det]->GetYaxis()->SetTitle("Cluster Charge Q=S_L+S_R");
        histoName.str("");
        histoName <<"hEtaVsSNR_" << TPlaneProperties::getStringForDetector(det);
        max = 400;
        hEtaVsSNR[det] = new TH2F(histoName.str().c_str(),histoName.str().c_str(),nbins,0,1,nbins,min,max);
        hEtaVsSNR[det]->GetXaxis()->SetTitle("#eta = S_R/(S_L+S_R)");
        hEtaVsSNR[det]->GetYaxis()->SetTitle("SNR of biggest hit");
    }

}

void TAnalysisOfPedestal::savePHinSigmaHistos(){
	Float_t xMinBiggest,xMaxBiggest,xMinAdjacent,xMaxAdjacent;
	Float_t yMaxBiggest,yMaxAdjacent;
	//hBiggestSignalInSigma
	for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++) {
		TString histoName =TString::Format("hPulseHeight_BiggestSignalInSigma%s",TPlaneProperties::getStringForDetector(det).c_str()) ;
		double cut = settings->getClusterSeedFactor(det,0);
		cout << "saving histogram " << histoName<< ".. with CUT on " <<cut<< endl;
		Float_t max=0;
		Float_t mean=0;
		Float_t sigma=0;
		Int_t entries = 0;

		//Find max mean and sigma
		for(UInt_t i=0;i<vecBiggestSignalInSigma[det].size();i++){
		    if(vecBiggestSignalInSigma[det][i]<4)
		        continue;
		    entries ++;
			mean+=vecBiggestSignalInSigma[det][i];
			sigma+=vecBiggestSignalInSigma[det][i]*vecBiggestSignalInSigma[det][i];
			if (vecBiggestSignalInSigma[det][i]>max)
				max = vecBiggestSignalInSigma[det][i];
		}
		mean/=entries;//(Float_t)vecBiggestSignalInSigma[det].size();
		sigma=TMath::Sqrt(sigma/entries-mean*mean);
		cout<< "Mean: "<<mean<<" +/- "<<sigma<<"\tMaximum SNR: "<<max<<" "<<entries<<endl;
		//define xrange and nbins
//        Float_t xRangeMax = TMath::Min(mean+3*sigma,max);
		Float_t xRangeMax;
		if (TPlaneProperties::isDiamondDetector(det)) xRangeMax = settings->getSnr_distribution_di_max();
		else                                          xRangeMax = settings->getSnr_distribution_si_max();
        Float_t xRangeMin = 0;
        xMinBiggest = xRangeMin;
        xMaxBiggest = xRangeMax;
        UInt_t nbins = settings->getSnr_distribution_num_bins();
        cout<<"X Range: "<<xRangeMin<<","<<xRangeMax<<endl;
        hBiggestSignalInSigma[det] = new TH1F(histoName,histoName,nbins,xRangeMin,xRangeMax);
        histoName = histoName+"_CMNcorrected";
        hBiggestSignalInSigmaCMN[det] = new TH1F(histoName,histoName,nbins,xRangeMin,xRangeMax);
        //set Axis Titles
        this->hBiggestSignalInSigma[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
        this->hBiggestSignalInSigma[det]->GetYaxis()->SetTitle("number of entries #");
        this->hBiggestSignalInSigmaCMN[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
        this->hBiggestSignalInSigmaCMN[det]->GetYaxis()->SetTitle("number of entries #");
        histoName =TString::Format("hPulseHeight_BiggestSignalChannelInSigma%s_2D",TPlaneProperties::getStringForDetector(det).c_str());
        Int_t nChannels = TPlaneProperties::getNChannels(det);
        TH2F* hBiggestSignalInSigma2D = new TH2F(histoName,histoName,1000,0,500,nChannels,0,nChannels-1);

        histoName =TString::Format("hPulseHeight_BiggestSignalChannelInSigma%s_CMNcorrected_2D",TPlaneProperties::getStringForDetector(det).c_str()) ;
        TH2F* hBiggestSignalInSigma2DCMN = new TH2F(histoName,histoName,1000,0,500,nChannels,0,nChannels-1);
        //set Axis Titles
        hBiggestSignalInSigma2D->GetXaxis()->SetTitle("Biggest Signal PH in units of sigma");
        hBiggestSignalInSigma2D->GetYaxis()->SetTitle("Channel No of Biggest Hit");
        hBiggestSignalInSigma2D->GetZaxis()->SetTitle("number of entries #");

        hBiggestSignalInSigma2DCMN->GetXaxis()->SetTitle("Biggest Signal PH in units of sigma");
        hBiggestSignalInSigma2DCMN->GetYaxis()->SetTitle("Channel No of Biggest Hit");
        hBiggestSignalInSigma2DCMN->GetZaxis()->SetTitle("number of entries #");
        for(UInt_t i=0;i<vecBiggestSignalInSigma[det].size();i++){
            Float_t signal=vecBiggestSignalInSigma[det].at(i);
            Int_t channel = vecBiggestHitChannel[det].at(i);
            this->hBiggestSignalInSigma[det]->Fill(signal);
            hBiggestSignalInSigma2D->Fill(signal,channel);
        }
        for(UInt_t i=0;i<vecBiggestSignalInSigmaCMN[det].size();i++){
            Float_t signal=vecBiggestSignalInSigmaCMN[det].at(i);
            Int_t channel = vecBiggestHitChannelCMN[det].at(i);
            this->hBiggestSignalInSigmaCMN[det]->Fill(signal);
            hBiggestSignalInSigma2DCMN->Fill(signal,channel);
        }
        Float_t yMaxNorm = findYPlotRangeForPHHisto(hBiggestSignalInSigma[det],settings->getClusterSeedFactor(det,0));
        Float_t yMaxCMN = findYPlotRangeForPHHisto(hBiggestSignalInSigmaCMN[det],settings->getClusterSeedFactor(det,0));
        Float_t yMax = TMath::Max(yMaxNorm,yMaxCMN);
        //set X Axis Range to the same for all
        Float_t xmin = hBiggestSignalInSigma[det]->GetXaxis()->GetXmin();
        Float_t xmax = hBiggestSignalInSigma[det]->GetXaxis()->GetXmax();
        hBiggestSignalInSigma2DCMN->Draw("colz");
        hBiggestSignalInSigma2D->Draw("colz");
        hBiggestSignalInSigma2D->GetZaxis()->SetRangeUser(0,yMax);
        hBiggestSignalInSigma2DCMN->GetZaxis()->SetRangeUser(0,yMax);
        hBiggestSignalInSigma2D->GetXaxis()->SetRangeUser(xmin,xmax);
        hBiggestSignalInSigma2DCMN->GetXaxis()->SetRangeUser(xmin,xmax);

        hBiggestSignalInSigmaCMN[det]->GetXaxis()->SetRangeUser(xmin,xmax);

        cout<<"BiggestHitSNR:\t"<<det<<"\t"<<hBiggestSignalInSigmaCMN[det]<<"\t"<<hBiggestSignalInSigmaCMN[det]->GetEntries()<<"\t"<<hBiggestSignalInSigmaCMN[det]->GetName()<<endl;
        histSaver->SaveHistogram(hBiggestSignalInSigma2D);
        cout<<"save "<<hBiggestSignalInSigmaCMN[det]->GetTitle()<<endl;
        histSaver->SaveHistogram(hBiggestSignalInSigmaCMN[det]);
        cout<<"save "<<hBiggestSignalInSigma2DCMN->GetTitle()<<endl;
        histSaver->SaveHistogram(hBiggestSignalInSigma2DCMN);
        yMaxBiggest = hBiggestSignalInSigma[det]->GetYaxis()->GetXmax();
        TString canvasName = TString::Format("%s",this->hBiggestSignalInSigma[det]->GetTitle());
        TCanvas *c1 = new TCanvas(canvasName,canvasName);
        c1->cd();
        hBiggestSignalInSigma[det]->GetYaxis()->SetRangeUser(0,yMax);
        this->hBiggestSignalInSigma[det]->Draw();

        double xCor[] = {cut,cut};
        double yCor[] = {0,this->hBiggestSignalInSigma[det]->GetBinContent(hBiggestSignalInSigma[det]->GetMaximumBin())*2};
        TGraph* lineGraph = new TGraph(2,xCor,yCor);
        lineGraph->SetLineColor(kRed);
        lineGraph->SetLineWidth(2);
        lineGraph->SetTitle("Seed Cut");
        lineGraph->Draw("Lsame");
        cout<<"save "<<c1->GetTitle()<<endl;
        histSaver->SaveCanvas(c1);
        if(TPlaneProperties::isDiamondDetector(det)){
            TString canvasName1 = TString::Format("hPH_BiggestSignalInSigma_Compare_%s",TPlaneProperties::getStringForDetector(det).c_str());
            TCanvas *c2 = new TCanvas(canvasName1,canvasName1);
            c2->cd();
            cout<<"\n\n*****\n1: "<<hBiggestSignalInSigma[det]<<"\t"<<hBiggestSignalInSigma[det]->GetTitle()<<endl;
            this->hBiggestSignalInSigma[det]->Draw();
            this->hBiggestSignalInSigmaCMN[det]->SetLineColor(kBlue);
            hBiggestSignalInSigma[det]->GetYaxis()->SetRangeUser(0,yMax);
            cout<<"2: "<<hBiggestSignalInSigmaCMN[det]<<"\t"<<hBiggestSignalInSigmaCMN[det]->GetTitle()<<endl;
            this->hBiggestSignalInSigmaCMN[det]->Draw("same");
            cout<<"3: "<<lineGraph<<"\t"<<lineGraph->GetTitle()<<endl;
            lineGraph->Draw("Lsame");
            //x1,y1,x2,y2
            TLegend* leg = new TLegend(0.13,0.55,0.45,0.85);//c1->BuildLegend(0.13,0.55,0.45,0.85);
            leg->AddEntry(hBiggestSignalInSigma[det],"normal");
            leg->AddEntry(hBiggestSignalInSigmaCMN[det],"CMN corrected");
            leg->SetFillColor(kWhite);
            leg->Draw();
            cout<<"save Canvas: "<<c2->GetTitle()<<endl;
            histSaver->SaveCanvas(c2);
        }
        //        histSaver->SaveHistogram(this->histo_pulseheight_sigma[det]);
        delete hBiggestSignalInSigma[det];
        delete lineGraph;
        delete c1;
        if (TPlaneProperties::isDiamondDetector(det)){
            Int_t diamondAreas = settings->getNDiaDetectorAreas();
            for(Int_t i = 0;i < diamondAreas;i++){
                std::pair<Int_t, Int_t> area = settings->getDiaDetectorArea(i);
                if(area.first<area.second){
                    cout<<" save NDiaDetecotAreas "<<i<<"/"<<diamondAreas<<endl;
                    TString name = TString::Format("hPulseHeight_BiggestSignalInSigma%s_area_%d_ch%d-ch%d",TPlaneProperties::getStringForDetector(det).c_str(),i,area.first,area.second);
                    cout<<"Create "<<name<<endl;
                    Int_t binMin = hBiggestSignalInSigma2D->FindBin(area.first);
                    Int_t binMax = hBiggestSignalInSigma2D->FindBin(area.second);
                    TH1F* hBiggestPHinSigma_SubArea = (TH1F*) hBiggestSignalInSigma2D->ProjectionX(name,binMin,binMax);
                    if(!hBiggestPHinSigma_SubArea)continue;
                    hBiggestPHinSigma_SubArea->SetTitle(name);
                    hBiggestPHinSigma_SubArea->SetName(name);
                    //					Float_t ymax =
                    findYPlotRangeForPHHisto(hBiggestPHinSigma_SubArea,settings->getClusterSeedFactor(det,(area.second+area.first)/2));

                    name = TString::Format("hPulseHeight_BiggestSignalInSigma_area%d_ch%d-ch%d_CMNcorrected",i,area.first,area.second);
                    cout<<"Create "<<name<<endl;
                    binMin = hBiggestSignalInSigma2DCMN->FindBin(area.first);
                    binMax = hBiggestSignalInSigma2DCMN->FindBin(area.second);
                    TH1F* hBiggestPHInSigma_SubArea_CMNcorrected = (TH1F*) hBiggestSignalInSigma2DCMN->ProjectionX(name,binMin,binMax);
                    if(!hBiggestPHInSigma_SubArea_CMNcorrected) continue;
                    hBiggestPHInSigma_SubArea_CMNcorrected->SetTitle(name);
                    findYPlotRangeForPHHisto(hBiggestPHInSigma_SubArea_CMNcorrected,settings->getClusterSeedFactor(det,(area.second+area.first)/2));
                    histSaver->CopyAxisRangesToHisto(hBiggestPHInSigma_SubArea_CMNcorrected,hBiggestPHinSigma_SubArea);
                    histSaver->SaveHistogram(hBiggestPHinSigma_SubArea);
                    histSaver->SaveHistogram(hBiggestPHInSigma_SubArea_CMNcorrected);

                    //Draw Both histos
                    histSaver->CopyAxisRangesToHisto(hBiggestPHInSigma_SubArea_CMNcorrected,hBiggestPHinSigma_SubArea);
                    name = "c_PulseHeightInSigma_Compare";
                    name.Append(TString::Format("_area%d_ch%d-ch%d",i,area.first,area.second));
                    TCanvas *c1 = new TCanvas(name,name);
                    c1->cd();
                    hBiggestPHInSigma_SubArea_CMNcorrected->SetLineColor(kGreen);

                    hBiggestPHinSigma_SubArea->Draw();
                    hBiggestPHInSigma_SubArea_CMNcorrected->Draw("same");
                    hBiggestPHInSigma_SubArea_CMNcorrected->SetFillStyle(hBiggestPHinSigma_SubArea->GetFillStyle());
                    if(hBiggestPHinSigma_SubArea->GetFillColor()!=kWhite)
                        hBiggestPHInSigma_SubArea_CMNcorrected->SetFillColor(kGreen-8);
                    TLegend* leg1 = c1->BuildLegend(0.15,0.55,0.4,0.80);
                    leg1->SetFillColor(kWhite);
                    histSaver->SaveCanvas(c1);
                    Float_t xMinBiggest = 0;
                    Float_t xMinAdjacent = 0;
                    Float_t xMaxBiggest = hBiggestPHinSigma_SubArea->GetXaxis()->GetXmax();
                    Float_t xMaxAdjacent = xMaxBiggest;
                    TString histoName =TString::Format("hBiggestAndBiggestAdjacentSignal_in_SNR_%s_area%d_ch%d-%d",TPlaneProperties::getStringForDetector(det).c_str(),i,area.first,area.second) ;
                    TH2F *hBiggestAndBiggestAdjacentSignal_in_SNR = new TH2F(histoName,histoName,256,xMinBiggest,xMaxBiggest,256,xMinAdjacent,xMaxAdjacent);
                    histoName=histoName.Append("_CMNcorrected");
                    TH2F *hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected = new TH2F(histoName,histoName,256,xMinBiggest,xMaxBiggest,256,xMinAdjacent,xMaxAdjacent);
                    UInt_t bhc      = vecBiggestHitChannel[det].size();
                    UInt_t bhcCMN   = vecBiggestHitChannelCMN[det].size();
                    UInt_t bSNR     = vecBiggestSignalInSigma[det].size();
                    UInt_t bSNRCMN  = vecBiggestSignalInSigmaCMN[det].size();
                    UInt_t bahc     = vecBiggestAdjacentHitChannel[det].size();
                    UInt_t bahcCMN  = vecBiggestAdjacentHitChannelCMN[det].size();
                    UInt_t baSNR    = vecBiggestAdjacentSignalInSigma[det].size();
                    UInt_t baSNRCMN = vecBiggestAdjacentSignalInSigmaCMN[det].size();
                    if(bhc!=bhcCMN||bhcCMN!=bSNR||bSNR!=bSNRCMN||bSNRCMN!=bahc||bahc!=bahcCMN||bahcCMN!=baSNR||baSNR!=baSNRCMN)
                        cout<<"Vector Sizes do NOT AGREE: "<<bhc<<" "<<bhcCMN<<" "<<bSNR<<" "<<bSNRCMN<<" "<<bahc<<" "<<bahcCMN<<" "<<baSNR
                        <<" "<<baSNRCMN<<endl;
                    int len = 4;
                    UInt_t minArrayCMN[] = {bhcCMN,bSNRCMN,bahcCMN,baSNRCMN};
                    UInt_t minArray[] = {bhc,bSNR,bahc,baSNR};
                    UInt_t minSize = *min_element(minArray,minArray+4);
                    UInt_t minSizeCMN = *min_element(minArrayCMN,minArrayCMN+len);

                    if(verbosity){
                        cout<<"minSize: "<<minSize<<endl;
                        if(verbosity%2==1){char t;cin>>t;}
                    }
                    for(UInt_t i=0;i<minSize||i<minSizeCMN;i++){
                        if(i<minSize){
                            if(area.first<=vecBiggestHitChannel[det].at(i)&&vecBiggestHitChannel[det].at(i)<=area.second){
                                Float_t biggestSNR = vecBiggestSignalInSigma[det].at(i);
                                Float_t adjacentSNR = vecBiggestAdjacentSignalInSigma[det].at(i);
                                hBiggestAndBiggestAdjacentSignal_in_SNR->Fill(biggestSNR,adjacentSNR);
                            }
                        }
                        if(i<minSize&&i<minSizeCMN){
                            if(area.first<=vecBiggestHitChannelCMN[det].at(i)&&vecBiggestAdjacentHitChannelCMN[det].at(i)<=area.second){
                                Float_t biggestSNR = vecBiggestSignalInSigmaCMN[det].at(i);
                                Float_t adjacentSNR = vecBiggestAdjacentSignalInSigmaCMN[det].at(i);
                                Int_t biggestCh  = vecBiggestHitChannel[det].at(i);
                                Int_t adjacentCh = vecBiggestAdjacentHitChannel[det].at(i);
                                if(biggestSNR<adjacentSNR){
                                    if(TMath::Abs(adjacentSNR/biggestSNR) > 1.05){
                                        if(!settings->isDet_channel_screened(det,biggestCh)&&!settings->isDet_channel_screened(det,adjacentCh)){
                                            if(verbosity) cout<<"Error3:"<<i<<" "<<det<< " "<<setw(3)<<biggestCh<<" "<<biggestSNR<<"@"<<biggestCh<<" "<<" - "<<adjacentSNR<<"@"<<adjacentCh<<endl;
                                            invalidBiggestHitChannels.at(det).insert(biggestCh);
                                        }
                                    }
                                }

                                else
                                    hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->Fill(biggestSNR,adjacentSNR);
                            }
                        }
                    }
                    if(!hBiggestAndBiggestAdjacentSignal_in_SNR)continue;
                    if(!hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected)continue;
                    hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetXaxis()->SetTitle("Biggest Signal in SNR");
                    hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetYaxis()->SetTitle("Biggest adjacent Signal in SNR");
                    hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetZaxis()->SetTitle("number of entries #");
                    TH2F* hist = hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected;
                    if(!hist)continue;
                    Float_t zmax = 0;//hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetBinContent(hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetMaximumBin());
                    for(Int_t i = 0; i<hist->GetXaxis()->GetNbins();i++)
                        for(Int_t j=0;j<hist->GetYaxis()->GetNbins();j++){
                            Float_t xPos = hist->GetXaxis()->GetBinCenter(i);
                            Float_t binContent = hist->GetBinContent(j,j);
                            if(xPos>5&&binContent>zmax)
                                zmax=binContent;
                        }
                    hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetZaxis()->SetRangeUser(0,zmax*1.3);

                    hBiggestAndBiggestAdjacentSignal_in_SNR->GetXaxis()->SetTitle("Biggest Signal in SNR");
                    hBiggestAndBiggestAdjacentSignal_in_SNR->GetYaxis()->SetTitle("Biggest adjacent Signal in SNR");
                    hBiggestAndBiggestAdjacentSignal_in_SNR->GetZaxis()->SetTitle("number of entries #");
                    zmax = hBiggestAndBiggestAdjacentSignal_in_SNR->GetBinContent(hBiggestAndBiggestAdjacentSignal_in_SNR->GetMaximumBin());
                    hist = hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected;
                    zmax = 0;//hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetBinContent(hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected->GetMaximumBin());
                    for(Int_t i = 0; i<hist->GetXaxis()->GetNbins();i++)
                        for(Int_t j=0;j<hist->GetYaxis()->GetNbins();j++){
                            Float_t xPos = hist->GetXaxis()->GetBinCenter(i);
                            Float_t binContent = hist->GetBinContent(j,j);
                            if(xPos>5&&binContent>zmax)
                                zmax=binContent;
                        }
                    hBiggestAndBiggestAdjacentSignal_in_SNR->GetZaxis()->SetRangeUser(0,zmax*1.3);

                    histSaver->SaveHistogramLogZ(hBiggestAndBiggestAdjacentSignal_in_SNR_CMNcorrected);
                    histSaver->SaveHistogramLogZ(hBiggestAndBiggestAdjacentSignal_in_SNR);
                    if(hBiggestPHInSigma_SubArea_CMNcorrected)delete hBiggestPHInSigma_SubArea_CMNcorrected;
                    if(hBiggestPHinSigma_SubArea)delete hBiggestPHinSigma_SubArea;
                }
            }

        }
        delete hBiggestSignalInSigma2D;
    }

	//hBiggestAdjacentSignalInSigma
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){

		double cut = settings->getClusterHitFactor(det,0);
		//		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".. with CUT on " <<cut<< endl;
		TCanvas *c1 = new TCanvas(this->hBiggestAdjacentSignalInSigma[det]->GetTitle(),this->hBiggestAdjacentSignalInSigma[det]->GetTitle());
		c1->cd();
		this->hBiggestAdjacentSignalInSigma[det]->Draw();
		double xCor[] = {cut,cut};
		double yCor[] = {0,this->hBiggestAdjacentSignalInSigma[det]->GetMaximum()*2};
		this->hBiggestAdjacentSignalInSigma[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
		this->hBiggestAdjacentSignalInSigma[det]->GetYaxis()->SetTitle("number of entries #");
		findYPlotRangeForPHHisto(hBiggestAdjacentSignalInSigma[det],settings->getClusterHitFactor(det,0));
		xMinAdjacent = 0;
		xMaxAdjacent = hBiggestAdjacentSignalInSigma[det]->GetXaxis()->GetXmax();
		yMaxAdjacent = hBiggestAdjacentSignalInSigma[det]->GetYaxis()->GetXmax();
		TGraph* lineGraph = new TGraph(2,xCor,yCor);
		lineGraph->SetLineColor(kRed);
		lineGraph->SetLineWidth(2);
		lineGraph->Draw("Lsame");
		histSaver->SaveCanvas(c1);;
		if(TPlaneProperties::isDiamondDetector(det)){
			double cut = settings->getClusterHitFactor(det,0);
			//		cout << "saving histogram " << this->histo_pulseheight_sigma_second[det]->GetName() << ".. with CUT on " <<cut<< endl;
			TCanvas *c2 = new TCanvas(this->hBiggestAdjacentSignalInSigmaCMN[det]->GetTitle(),this->hBiggestAdjacentSignalInSigmaCMN[det]->GetTitle());
			c2->cd();
			this->hBiggestAdjacentSignalInSigmaCMN[det]->Draw();
			double xCorCMN[] = {cut,cut};
			double yCorCMN[] = {0,this->hBiggestAdjacentSignalInSigmaCMN[det]->GetMaximum()*2};
			this->hBiggestAdjacentSignalInSigmaCMN[det]->GetXaxis()->SetTitle("Biggest Hit PH in units of sigma");
			this->hBiggestAdjacentSignalInSigmaCMN[det]->GetXaxis()->SetTitle("number of entries #");
			hBiggestAdjacentSignalInSigmaCMN[det]->SetLineColor(kBlack);
			histSaver->CopyAxisRangesToHisto(hBiggestAdjacentSignalInSigmaCMN[det],hBiggestAdjacentSignalInSigma[det]);
			findYPlotRangeForPHHisto(hBiggestAdjacentSignalInSigmaCMN[det],cut);
			hBiggestAdjacentSignalInSigmaCMN[det]->Draw();
			TGraph* lineGraphCMN = new TGraph(2,xCorCMN,yCorCMN);
			lineGraphCMN->SetLineColor(kRed);
			lineGraphCMN->SetLineWidth(2);
			lineGraphCMN->Draw("Lsame");
			histSaver->SaveCanvas(c2);;
			delete hBiggestAdjacentSignalInSigmaCMN[det];
			delete lineGraphCMN;
			delete c2;
		}
		delete hBiggestAdjacentSignalInSigma[det];
		delete lineGraph;
		delete c1;
	}
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		TString histoName =TString::Format("hBiggestAndBiggestAdjacentSignal_in_SNR_%s",TPlaneProperties::getStringForDetector(det).c_str()) ;
		TH2F *hBiggestAndBiggestAdjacentSignal_in_SNR = new TH2F(histoName,histoName,256,xMinBiggest,xMaxBiggest,256,xMinAdjacent,xMaxAdjacent);
		hBiggestAndBiggestAdjacentSignal_in_SNR->GetXaxis()->SetTitle("Biggest Signal in SNR");
		hBiggestAndBiggestAdjacentSignal_in_SNR->GetYaxis()->SetTitle("Biggest Adjacent Signal in SNR");
		hBiggestAndBiggestAdjacentSignal_in_SNR->GetZaxis()->SetTitle("number of entries # ");
		if(TMath::Max(yMaxAdjacent,yMaxBiggest)>1)
			hBiggestAndBiggestAdjacentSignal_in_SNR->GetZaxis()->SetRangeUser(0,TMath::Max(yMaxAdjacent,yMaxBiggest));
		for(UInt_t i=0;i<vecBiggestAdjacentSignalInSigma[det].size()&&i<vecBiggestSignalInSigma[det].size();i++){
			Float_t biggestSNR = vecBiggestSignalInSigma[det].at(i);
			Float_t adjacentSNR= vecBiggestAdjacentSignalInSigma[det].at(i);
			hBiggestAndBiggestAdjacentSignal_in_SNR->Fill(biggestSNR,adjacentSNR);
		}
		histSaver->SaveHistogram(hBiggestAndBiggestAdjacentSignal_in_SNR);

		histoName = histoName.Append("_CMNcorrected");
		TH2F *hBiggestAndBiggestAdjacentSignal_in_SNR_CMN = new TH2F(histoName,histoName,256,xMinBiggest,xMaxBiggest,256,xMinAdjacent,xMaxAdjacent);
		if(TMath::Max(yMaxAdjacent,yMaxBiggest)>1)
			hBiggestAndBiggestAdjacentSignal_in_SNR_CMN->GetZaxis()->SetRangeUser(0,TMath::Max(yMaxAdjacent,yMaxBiggest));
		hBiggestAndBiggestAdjacentSignal_in_SNR_CMN->GetXaxis()->SetTitle("Biggest Signal in SNR");
		hBiggestAndBiggestAdjacentSignal_in_SNR_CMN->GetYaxis()->SetTitle("Biggest Adjacent Signal in SNR");
		hBiggestAndBiggestAdjacentSignal_in_SNR_CMN->GetZaxis()->SetTitle("number of entries # ");
		for(UInt_t i=0;i<vecBiggestSignalInSigmaCMN[det].size()&&i<vecBiggestSignalInSigmaCMN[det].size();i++){
			Float_t biggestSNR = vecBiggestSignalInSigmaCMN[det].at(i);
			Float_t adjacentSNR= vecBiggestAdjacentSignalInSigmaCMN[det].at(i);
			hBiggestAndBiggestAdjacentSignal_in_SNR_CMN->Fill(biggestSNR,adjacentSNR);
		}
		histSaver->SaveHistogram(hBiggestAndBiggestAdjacentSignal_in_SNR_CMN);
	}
}


void TAnalysisOfPedestal::saveHistos(){
    saveAdcVsEventProfiles();
    histSaver->SaveHistogram(hRelCmnUncertainty);
    delete hRelCmnUncertainty;
    histSaver->SaveHistogram(hCmnNUsedChannels);
    delete hCmnNUsedChannels;
    histSaver->SaveHistogram(hCmnChannelWeight);
    delete hCmnChannelWeight;
    hNewCmnVsEventNo->GetXaxis()->SetRangeUser(0,settings->getNEvents()/10.);
    histSaver->SaveHistogram(hNewCmnVsEventNo);
    TProfile *prof = hNewCmnVsEventNo->ProfileX();
    if (prof){
        Float_t mean =hNewCmnVsEventNo->GetMean(2);
        Float_t rms = hNewCmnVsEventNo->GetRMS(2);
        Float_t xmin = mean - 3 * rms;
        Float_t xmax = mean + 3 * rms;
        prof->GetYaxis()->SetRangeUser(xmin,xmax);
        cout<<"Set Range: "<<xmin<<"-"<<xmax<<endl;
        histSaver->SaveHistogram(prof,false,false,false);
        delete prof;
        prof=0;
    }
    hNewCmnVsEventNo->GetXaxis()->SetRangeUser(0,settings->getNEvents());
    hNewCmnVsEventNo->SetName("hNewCmnVsEventNo_all");
    histSaver->SaveHistogram(hNewCmnVsEventNo);

    prof = hNewCmnVsEventNo->ProfileX();
    if (prof){
        Float_t mean =hNewCmnVsEventNo->GetMean(2);
        Float_t rms = hNewCmnVsEventNo->GetRMS(2);
        prof->GetYaxis()->SetRangeUser(mean-3*rms, mean+3*rms);
        histSaver->SaveHistogram(prof,false,false,false);
        delete prof;
        prof=0;
    }
    TF1* fit= new TF1("fGauss_CMN","gaus",-30,30);
    histSaver->SaveHistogramWithFit(hNewComonModeNoise,fit);
    res->setFloatValue("Noise","CMN_Pos",fit->GetParameter(1));
    res->setFloatValue("Noise","CMN_Sigma",fit->GetParameter(2));
    delete hNewComonModeNoise;
    delete hNewCmnVsEventNo;
    histSaver->SaveHistogram(hCmnVsNewCmn);
    delete hCmnVsNewCmn;
    histSaver->SaveHistogram(hCmnNewVsNUsedChannels);
    delete hCmnNewVsNUsedChannels;
    histSaver->SaveHistogram(hCmnUsedChannels);
    delete hCmnUsedChannels;
    prof = hCmnChannelWeightVsChannel->ProfileX();
    if (prof){
        histSaver->SaveHistogram(prof,false,false,false);
        delete prof;
        prof=0;
    }
    histSaver->SaveHistogram(hCmnChannelWeightVsChannel);
    delete hCmnChannelWeightVsChannel;

    prof = hCmnFractionVsChannel->ProfileX();
    if (prof){
        histSaver->SaveHistogram(prof,false,false,false);
        delete prof;
        prof=0;
    }
    histSaver->SaveHistogram(hCmnFractionVsChannel);
    delete hCmnFractionVsChannel;

    map<Int_t,TH1F*>::iterator it;
    for(it=hBiggestSignalInSigmaDiaPattern.begin();it!=hBiggestSignalInSigmaDiaPattern.end();it++){
        SetYRangeForSignalInSigmaPlot((*it).second);
        histSaver->SaveHistogram((*it).second,false,false,true);
        delete (*it).second;
    }
    for(it=hBiggestSignalInSigmaDiaPatternCMN.begin();it!=hBiggestSignalInSigmaDiaPatternCMN.end();it++){
        SetYRangeForSignalInSigmaPlot((*it).second);
        histSaver->SaveHistogram((*it).second,false,false,true);
        delete (*it).second;
    }
    for(it=hBiggestAdjacentSignalInSigmaDiaPattern.begin();it!=hBiggestAdjacentSignalInSigmaDiaPattern.end();it++){
        SetYRangeForSignalInSigmaPlot((*it).second);
        histSaver->SaveHistogram((*it).second,false,false,true);
        delete (*it).second;
    }
    for(it=hBiggestAdjacentSignalInSigmaDiaPatternCMN.begin();it!=hBiggestAdjacentSignalInSigmaDiaPatternCMN.end();it++){
        SetYRangeForSignalInSigmaPlot((*it).second);
        histSaver->SaveHistogram((*it).second,false,false,true);
        delete (*it).second;
    }

	createPedestalMeanHistos();
	savePHinSigmaHistos();
	for (int det=0;det<9;det++){
		if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
		histSaver->SaveHistogramPNG(hSaturatedChannels[det]);
		hSaturatedChannels[det]->Delete();

        histSaver->SaveHistogram(hNoiseDistributionCMN[det]);
        delete hNoiseDistributionCMN[det];

        histSaver->SaveHistogram(hNoiseDistribution[det]);
        delete hNoiseDistribution[det];
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
        histSaver->SaveHistogram(hSNR_BiggestSignal[det]);
        hSNR_BiggestSignal[det]->Delete();
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
        histSaver->SaveHistogram(this->hHitOrderMap[det]);
        //		cout << "saving histogram " << this->histo_biggest_hit_map[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(hLeftVsRightSignal[det]);
        histSaver->SaveHistogram(hEtaVsCharge[det]);
        TH1D* hpx = hEtaVsCharge[det]->ProjectionX();
        histSaver->SaveHistogram(hpx);
        delete hpx;
        histSaver->SaveHistogram(hEtaVsSNR[det]);
        Int_t ybin = hEtaVsSNR[det]->GetYaxis()->FindBin(settings->getClusterSeedFactor(det,0));
        hpx =  hEtaVsSNR[det]->ProjectionX("",ybin);
        histSaver->SaveHistogram(hpx);
        delete hpx;

        histSaver->SaveHistogram(hLeftVsRightSignalCMN[det]);
        histSaver->SaveHistogram(this->hBiggestHitChannelMap[det]);
        histSaver->SaveHistogram(this->hBiggestHitChannelMapCMN[det]);
        if(TPlaneProperties::isDiamondDetector(det)){
            for (int pattern = 0; pattern < settings->getNDiaDetectorAreas();pattern++){
                pair<Int_t,Int_t> area = settings->getDiaDetectorArea(pattern);
                hBiggestHitChannelMap[det]->GetXaxis()->SetRangeUser(area.first-1,area.second+1);
                TString name = hBiggestHitChannelMap[det]->GetName();
                hBiggestHitChannelMap[det]->SetName(name+TString::Format("_pattern%d",pattern));
                histSaver->SaveHistogram(hBiggestHitChannelMap[det],false,false,true);
                hBiggestHitChannelMap[det]->SetName(name);
                //***
                hBiggestHitChannelMapCMN[det]->GetXaxis()->SetRangeUser(area.first-1,area.second+1);
                name = hBiggestHitChannelMapCMN[det]->GetName();
                hBiggestHitChannelMapCMN[det]->SetName(name+TString::Format("_pattern%d",pattern));
                histSaver->SaveHistogram(hBiggestHitChannelMapCMN[det],false,false,true);
                hBiggestHitChannelMapCMN[det]->SetName(name);
            }
        }
        //		cout << "saving histogram " << this->histo_pulseheight_left_sigma[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(this->histo_pulseheight_left_sigma[det]);
        //		cout << "saving histogram " << this->histo_pulseheight_left_sigma_second[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(this->histo_pulseheight_left_sigma_second[det]);
        //		cout << "saving histogram " << this->histo_pulseheight_right_sigma[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(this->histo_pulseheight_right_sigma[det]);
        //		cout << "saving histogram" << this->histo_pulseheight_right_sigma_second[det]->GetName() << ".." << endl;
        histSaver->SaveHistogram(this->histo_pulseheight_right_sigma_second[det]);

        //		delete histo_pulseheight_sigma125[det];
        delete hHitOrderMap[det];
        delete hLeftVsRightSignal[det];
        delete hLeftVsRightSignalCMN[det];
        delete hBiggestHitChannelMap[det];
        delete hBiggestHitChannelMapCMN[det];
        delete histo_pulseheight_left_sigma[det];
        delete histo_pulseheight_left_sigma_second[det];
        delete histo_pulseheight_right_sigma[det];
        delete histo_pulseheight_right_sigma_second[det];

        if(TPlaneProperties::isSiliconDetector(det)){
            TF1 histofitx("histofitx","gaus",hAllAdcNoise[det]->GetMean()-2*hAllAdcNoise[det]->GetRMS(),hAllAdcNoise[det]->GetMean()+2*hAllAdcNoise[det]->GetRMS());
            histofitx.SetLineColor(kBlue);
            hAllAdcNoise[det]->Fit(&histofitx,"rq");
            if(res!=0)res->setNoise(det,histofitx.GetParameter(2));
            histSaver->SaveHistogram(hAllAdcNoise[det],true);

			delete hAllAdcNoise[det];
		}
	}
	if(verbosity){
		cout<<"hFitX:"<<hDiaAllAdcNoise<<"\t"<<flush;
		cout<<"Mean: "<<hDiaAllAdcNoise->GetMean()<<"\t"<<flush;
		cout<<"'RMS': "<<hDiaAllAdcNoise->GetRMS()<<"\t"<<flush;
	}
	TF1 histofitx("histofitx","gaus",hDiaAllAdcNoise->GetMean()-2*hDiaAllAdcNoise->GetRMS(),hDiaAllAdcNoise->GetMean()+2*hDiaAllAdcNoise->GetRMS());
	histofitx.SetLineColor(kBlue);
	hDiaAllAdcNoise->Fit(&histofitx,"rq");

	Float_t diaNoise = histofitx.GetParameter(2);
	histSaver->SaveHistogram(hDiaAllAdcNoise,true);
	TF1 histofitAllNoiseCMN("histofitx","gaus",hDiaAllAdcNoiseCMN->GetMean()-2*hDiaAllAdcNoiseCMN->GetRMS(),hDiaAllAdcNoiseCMN->GetMean()+2*hDiaAllAdcNoiseCMN->GetRMS());
	histofitAllNoiseCMN.SetLineColor(kBlue);
	hDiaAllAdcNoiseCMN->Fit(&histofitAllNoiseCMN,"rq");
//	if(res!=0)res->SetNoise(TPlaneProperties::getDetDiamond(),histofitx.GetParameter(2));
	histSaver->SaveHistogram(hDiaAllAdcNoiseCMN,true);
    histSaver->SaveHistogram(hDiaAllAdcNoiseCMNChannel,true);
    histSaver->SaveHistogram(hDiaAllAdcNoiseChannel,true);
    histSaver->SaveHistogram(hDiaAllAdcNoiseEventNo, true);
    histSaver->SaveHistogram(hDiaAllAdcNoiseCMNEventNo, true);
    THStack *hstack = new THStack("hStackNoise","Noise per Area;Noise / ADC; no. of entries");
    for (UInt_t area = 0; area<settings->getNDiaDetectorAreas();area++){
    	TString ext = TString::Format("_area%d",area);
    	TString hname = TString::Format("hDiaAllAdcNoiseCMN_area%d",area);
    	TString htitle = TString::Format("Noise - CM corrected - Area %d",area);
    	Int_t binlow = settings->getDiaDetectorArea(area).first;
    	binlow = hDiaAllAdcNoiseCMNChannel->GetXaxis()->FindBin(binlow);
    	Int_t binhigh = settings->getDiaDetectorArea(area).second;
    	binhigh= hDiaAllAdcNoiseCMNChannel->GetXaxis()->FindBin(binhigh);
    	TH1D* hNoiseProjection = hDiaAllAdcNoiseCMNChannel->ProjectionY(hname,binlow,binhigh);
    	Float_t mean = hNoiseProjection->GetMean();
    	Float_t sigma = hNoiseProjection->GetRMS();
    	TF1* fGaus = new TF1((TString)"fGaus"+ext,"gaus",mean-2*sigma,mean+2*sigma);
    	hNoiseProjection->Fit(fGaus,"Q","",mean-sigma,mean+sigma);
    	histSaver->SaveHistogram(hNoiseProjection);
    	hstack->Add(hNoiseProjection);
    }
    if (settings->getNDiaDetectorAreas())
    	histSaver->SaveStack(hstack,"nostack",true,true,"Noise / ADC","no of entries");
    if (hstack)
    	delete hstack;

    Float_t noiseCMC = histofitAllNoiseCMN.GetParameter(2);
    for(Int_t ntries=0;hCMNoiseDistribution->GetBinContent(hCMNoiseDistribution->GetMaximumBin())<.05*hCMNoiseDistribution->GetEntries()&&ntries<2;ntries++)
        hCMNoiseDistribution->Rebin(2);
    TF1 histofitCMN("histofitx","gaus",hCMNoiseDistribution->GetMean()-2*hCMNoiseDistribution->GetRMS(),hCMNoiseDistribution->GetMean()+2*hCMNoiseDistribution->GetRMS());
    histofitCMN.SetLineColor(kBlue);
    hCMNoiseDistribution->Fit(&histofitCMN,"rq");

    Float_t cmn = histofitCMN.GetParameter(2);
    if(verbosity > 3){
        cout<<" diaNoise = "<<diaNoise<<endl;
        cout<<" CMN      = "<<cmn<<endl;
        cout<<" noiseCMC = "<<noiseCMC<<endl;
        char t;
        if(verbosity%2==1) cin>>t;
    }

	if(res!=0) res->setDiaNoiseCMcorrected(noiseCMC);
	if(res!=0) res->setCMN(cmn);
	if(res!=0) res->setNoise(TPlaneProperties::getDetDiamond(),diaNoise);
	histSaver->SaveHistogram(hCMNoiseDistribution,true);
	histSaver->SaveHistogram(hCMNoiseDistributionEventNo, true);
	stringstream name;
	name<< "gADC_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gADC = histSaver->CreateDipendencyGraph(name.str(),adcValues,eventNumbers,20000);
	gADC.SetMarkerStyle(7);

	name.str("");name.clear();
	name<< "gPed_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gPed = histSaver->CreateDipendencyGraph(name.str(),pedestalValues,eventNumbers,20000);
	gPed.SetLineColor(kBlue);
	gPed.SetMarkerColor(kBlue);
	gPed.SetLineWidth(2);


	name.str("");name.clear();
	name<< "gPedCMN_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gPedCMN = histSaver->CreateDipendencyGraph(name.str(),pedestalValuesCMN,eventNumbers,20000);
	gPedCMN.SetLineColor(kGreen);
	gPed.SetLineWidth(2);

	name.str("");name.clear();
	name<< "gUpperHitCut_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gUpperHitCut = histSaver->CreateDipendencyGraph(name.str(),upperHitCutValues,eventNumbers,20000);
	gUpperHitCut.SetLineColor(kRed);
	gUpperHitCut.SetMarkerColor(kRed);
	gUpperHitCut.SetLineWidth(2);

	name.str("");name.clear();
	name<< "gUpperSeedCut_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gUpperSeedCut = histSaver->CreateDipendencyGraph(name.str(),upperSeedCutValues,eventNumbers,20000);
	int n = gUpperSeedCut.GetN();
	double* y = gUpperSeedCut.GetY();
	Float_t max = MaxElement(n,y);
	gUpperSeedCut.SetLineColor(kMagenta);
	gUpperSeedCut.SetMarkerColor(kMagenta);
	gUpperSeedCut.SetLineWidth(2);


	name.str("");name.clear();
	name<< "gUpperHitCutCMN_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gUpperHitCutCMN = histSaver->CreateDipendencyGraph(name.str(),upperHitCutValuesCMN,eventNumbers,20000);
	gUpperHitCutCMN.SetLineColor(kOrange);

	name.str("");name.clear();
	name<< "gLowerHitCut_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gLowerHitCut = histSaver->CreateDipendencyGraph(name.str(),lowerHitCutValues,eventNumbers,20000);
	gLowerHitCut.SetLineColor(kRed);
	gLowerHitCut.SetMarkerColor(kRed);
	n = gLowerHitCut.GetN();
	y = gLowerHitCut.GetY();
	Float_t min = MinElement(n,y);
	gLowerHitCut.SetLineWidth(2);

	name.str("");name.clear();
	name << "gLowerHitCutCMN_ch"<<setw(3)<<setfill('0')<<settings->getNoisePlotChannel();
	TGraph gLowerHitCutCMN = histSaver->CreateDipendencyGraph(name.str(),lowerHitCutValuesCMN,eventNumbers,20000);
	gLowerHitCutCMN.SetLineColor(kOrange);

	name.str("");name.clear();
	name << "mgClusterCut_Ch_"<<settings->getNoisePlotChannel();
	TMultiGraph* mg = new TMultiGraph(name.str().c_str(),TString::Format("eventNumber vs adc - Ch %03d",settings->getNoisePlotChannel()));
	mg->Add(&gADC,"P");
	mg->Add(&gPed,"C");
	mg->Add(&gUpperSeedCut,"C");
	mg->Add(&gUpperHitCut,"C");
	mg->Add(&gLowerHitCut,"C");
	mg->Add(&gPedCMN,"C");
	mg->Add(&gUpperHitCutCMN,"C");
	mg->Add(&gLowerHitCutCMN,"C");

	name.str("");name.clear();
	name << "cClusterCut_Ch_"<<settings->getNoisePlotChannel();
	TCanvas *c1 = new TCanvas(name.str().c_str(),name.str().c_str(),1024,640);
	c1->cd();
	mg->Draw("AP");
	mg->GetXaxis()->SetTitle("event number");
	mg->GetXaxis()->SetRangeUser(0,10000);
	mg->GetYaxis()->SetTitle("Adc Value");
	mg->GetYaxis()->SetRangeUser(min-.1*(max-min),max+.3*(max-min));
	mg->Draw("A");
	TLegend* leg = c1->BuildLegend(0.5,0.75,0.88,0.99);
	leg->SetFillColor(kWhite);
	leg->Clear();
	leg->AddEntry(&gADC,"ADC","LP");
	leg->AddEntry(&gPed,"Pedestal","LP");
	Float_t nSigmasHit = settings->getClusterHitFactor(TPlaneProperties::getDetDiamond(),settings->getNoisePlotChannel());
	Float_t nSigmasSeed = settings->getClusterSeedFactor(TPlaneProperties::getDetDiamond(),settings->getNoisePlotChannel());
	leg->AddEntry(&gUpperSeedCut,  TString::Format("Seed Limit     - %1.0f #sigma",nSigmasSeed),"L");
	leg->AddEntry(&gUpperHitCut,   TString::Format("Hit Limit      -  %1.0f #sigma",nSigmasHit),"L");
	leg->AddEntry(&gLowerHitCutCMN,TString::Format("Hit Limit_{CMN} - %1.0f #sigma",nSigmasHit),"L");
	leg->AddEntry(&gPedCMN,"Pedestal_{CMN}","L");

	leg->Draw();
	histSaver->SaveCanvas(c1);
	delete c1;

	name.str("");name.clear();
	name << "mgClusterCut-CMN-Ch_"<<settings->getNoisePlotChannel();
	TMultiGraph* mgCMN = new TMultiGraph(name.str().c_str(),TString::Format("eventNumber vs adc - Ch %03d - CMN",settings->getNoisePlotChannel()));
	mgCMN->Add(&gADC);
	mgCMN->Add(&gPedCMN);
	mgCMN->Add(&gUpperHitCutCMN);
	mgCMN->Add(&gLowerHitCutCMN);
	name.str("");name.clear();
	name << "cClusterCut-CMN-Ch_"<<settings->getNoisePlotChannel();
	c1 = new TCanvas(name.str().c_str());
	c1->cd();
	mgCMN->Draw("APL");
	histSaver->SaveCanvas(c1);
	delete c1;

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

void TAnalysisOfPedestal::updateMeanCalulation(UInt_t det,UInt_t ch){
	if(det==0&&ch==0){
		sumPed =0;
		sumPedCMN=0;
		sumNoise=0;
		sumNoiseCMN=0;
		sumPed =0;
		nSumPedCMN=0;
		nSumNoiseCMN=0;
		nSumNoise=0;
        vecCMNoise.push_back(cmNoise);
        hCMNoiseDistribution->Fill(cmNoise);
        hCMNoiseDistributionEventNo->Fill(nEvent, cmNoise);
    }
    cmNoise = eventReader->getCMNoise(det,ch);

    if(settings->isDet_channel_screened(det,ch))
        return;
//    if(snr<settings->getClusterHitFactor(det,ch)){
    if(snr < settings->get_Pedestal_Hit_Factor(det)){
        pedestalMeanValue.at(det).at(ch) +=pedestal;
        pedestalSigmaValue.at(det).at(ch) +=sigma;
        nPedestalHits.at(det).at(ch)++;
        if(TPlaneProperties::isSiliconDetector(det))
            hAllAdcNoise[det]->Fill(noise);
        else if(TPlaneProperties::isDiamondDetector(det)){

            hDiaAllAdcNoise->Fill(noise);
            hDiaAllAdcNoiseCMN->Fill(noiseCMN);
            hDiaAllAdcNoiseChannel->Fill(ch,noise);
            hDiaAllAdcNoiseCMNChannel->Fill(ch,noiseCMN);
            hDiaAllAdcNoiseEventNo   ->Fill(nEvent, noise   );
            hDiaAllAdcNoiseCMNEventNo->Fill(nEvent, noiseCMN);
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
    if(det==TPlaneProperties::getNDetectors()-1 && ch <= TPlaneProperties::getNChannels(det) -1){
        vecAvrgPed.push_back(sumPed/(float)nSumPed);
        vecAvrgPedCMN.push_back(sumPedCMN/(float)nSumPedCMN);
        vecAvrgSigma.push_back(sumNoise/(float)nSumNoise);
        vecAvrgSigmaCMN.push_back(sumNoiseCMN/(float)nSumNoiseCMN);
        if(ch==0) vecEventNo.push_back(nEvent);
    }
}

void TAnalysisOfPedestal::saveAdcVsEventProfiles() {
    cout<<"TAnalysisOfPedestal::saveAdcVsEventProfiles"<<endl;
    for(UInt_t det =0; det<TPlaneProperties::getNDetectors();det++){
        cout<<"\tDetector: "<<det<<endl;
        TString name = (TString)"hADCProfiles_"+(TString)TPlaneProperties::getStringForDetector(det);
        if (!hHistoMap.count(name))
            continue;
        if(!hHistoMap[name])
            continue;
        TProfile2D* prof2d = (TProfile2D*) hHistoMap[name];
        prof2d->GetXaxis()->SetTitle("EventNo.");
        prof2d->GetYaxis()->SetTitle("Channel");
        prof2d->GetZaxis()->SetTitle("avrg adc. /1k events");
        histSaver->SaveHistogram(prof2d,false,false);
        if (TPlaneProperties::isSiliconDetector(det))
            continue;
        TF1* fit = new TF1("pol1fit","pol1",0,1e7);
        fit->SetLineWidth(1);
        fit->SetLineColor(kBlue);
        fit->SetLineStyle(2);
        TString name2 =name + (TString)"_ADCRatio";
        TH1F* hADCRatio = new TH1F(name2,name2, 512,0,2);
        hADCRatio->GetXaxis()->SetTitle("ratio adc/adc0");
        name2 =name + (TString)"_ADCRatio_Par0";
        TH1F* hADC_Fit_Par0 = new TH1F(name2,name2, 512,0,2);
        name2 =name + (TString)"_ADCRatio_Par1";
        TH1F* hADC_Fit_Par1 = new TH1F(name2,name2, 512,0,2);
        name2 =name + (TString)"_ADCRatio_ParChi2";
        TH1F* hADC_Fit_Chi2 = new TH1F(name2,name2, 80,0,40);
        vector<TH1F*>  vecADCPeakPositions;
        name2 =(TString)"hADCProfiles_"+(TString)TPlaneProperties::getStringForDetector(det);
        name2+=TString::Format("PeakPosition_SlidingWindow");
        TH1F* hPeakPos= new TH1F(name2,name2,prof2d->GetNbinsX(),0,prof2d->GetNbinsX()*1000);
        Int_t nChannels=0;
        for (UInt_t ch = 0; ch < TPlaneProperties::getNChannels(det); ch++){
            TRawEventSaver::showStatusBar(ch,TPlaneProperties::getNChannels(det),1,1,0);
            if(settings->IsMasked(det,ch))
                continue;
            TString name1 = name;
            name1 += TString::Format("_ch%03d",ch);
            TH1D* prof = prof2d->ProjectionX(name1,ch+1,ch+1);
            Float_t mean = 0;
            Int_t nAvrg = TMath::Min(7,prof->GetNbinsX());
            for(UInt_t i=1;i<=nAvrg;i++)mean+=prof->GetBinContent(i);
            TH1F* histo1 = doSlidingWindowAnalysis(prof,nAvrg);
            histo1->SetLineColor(kBlack);
            TH1F* histo3 = doSlidingWindowAnalysis(prof,10);
            histo3->SetLineColor(kRed);
            TH1F* histo4 = doSlidingWindowAnalysis(prof,15);
            mean/=nAvrg;
            prof->Scale(1./mean);
            histo4->SetLineColor(kGreen);
            histSaver->SaveHistogram(histo1,false,false,true,"PL");
            TCanvas* c1 = new TCanvas((TString)"c_"+histo1->GetName());
            int nPeaks = histo1->ShowPeaks(5,"",.1);
            histo4->Draw("PL");
            histo3->Draw("PLsame");
            histo1->Draw("PLsame");
            TLegend* leg = c1->BuildLegend();
            leg->SetFillColor(kWhite);
            histSaver->SaveCanvas(c1);
            TList *functions = histo1->GetListOfFunctions();
            TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
            if (pm){
                for (Int_t i =0; i<TMath::Min(nPeaks,4);i++){
                    Float_t position = pm->GetX()[i];
                    if(vecADCPeakPositions.size()<=i){
                        TString name2 =(TString)"hADCProfiles_"+(TString)TPlaneProperties::getStringForDetector(det);
                        name2+=TString::Format("PeakPosition_SlidingWindow_PeakNo%d",i+1);
                        vecADCPeakPositions.push_back( new TH1F(name2,name2,prof2d->GetNbinsX(),0,prof2d->GetNbinsX()*1000));
                    }
                    vecADCPeakPositions[i]->Fill(position);
                    hPeakPos->Fill(position);
                }
            }
            delete c1;
            if (fit)
                prof->Fit(fit,"Q");
            Float_t adc0 = prof->GetBinContent(1);
            for (UInt_t bin =1;bin <= prof->GetNbinsX();bin++){
                Float_t adc = prof->GetBinContent(bin);
                Float_t ratio = adc/adc0;
                hADCRatio->Fill(ratio);
            }
            if(fit){
                hADC_Fit_Par0->Fill(fit->GetParameter(0));
                hADC_Fit_Par1->Fill(fit->GetParameter(1));
                Float_t chi2 = fit->GetChisquare()/fit->GetNDF();
                //                cout<<det<<" "<<ch<<" "<<chi2<<endl;
                hADC_Fit_Chi2->Fill(chi2);
            }
            histSaver->SaveHistogram(prof,false,true,true);
            delete prof;
            nChannels++;
        } //for loop ch
        TCanvas* c1 = new TCanvas((TString)"c_ADCProfile_PeakPostions");
        for(UInt_t i =0; i< vecADCPeakPositions.size();i++){
            vecADCPeakPositions.at(i)->SetLineColor(i+1);
            if(i==0)
                vecADCPeakPositions.at(i)->Draw();
            else
                vecADCPeakPositions.at(i)->Draw("same");
        }
        histSaver->SaveHistogram(hPeakPos);
        Int_t nWindow = 1;
        Int_t nCandidates = 0;
        Int_t entries=0;
        for(UInt_t bin = 1; bin <= 2*nWindow+1 && bin <= hPeakPos->GetNbinsX();bin++)
            entries+=hPeakPos->GetBinContent(bin);
        bool counted = false;
        for(UInt_t bin=nWindow; bin<=hPeakPos->GetNbinsX()-nWindow;bin++){
            if(entries>.8*nChannels){
                cout<<"FOUND a Problematic candiate at bin"<<bin<<" --"<<hPeakPos->GetBinCenter(bin)<<endl;
                if(!counted)
                    nCandidates++;
                counted=true;
            }
            else
                counted = false;
            entries += hPeakPos->GetBinContent(bin+nWindow);
            entries -= hPeakPos->GetBinContent(bin-nWindow-1);
        }
        res->setIntValue("RunInfo","ProblematicAdcSteps",nCandidates);
        cout<<"\nThere are "<<nCandidates<<" candidates for adc jumps.";
        if (verbosity%2==1){
            cout<<"Press a key and enter"<<endl;
            char t; cin>>t;
        }

        delete hPeakPos;
        histSaver->SaveCanvas(c1);
        delete c1;
        vecADCPeakPositions.clear();
        hADCRatio->GetXaxis()->SetTitle("avrg adc / adc_{0}");
        hADCRatio->GetYaxis()->SetTitle("no. of entries");
        hADC_Fit_Par0->GetXaxis()->SetTitle("avrg adc fit par0");
        hADC_Fit_Par0->GetYaxis()->SetTitle("no of entries");
        hADC_Fit_Par1->GetXaxis()->SetTitle("avrg adc fit par1");
        hADC_Fit_Par1->GetYaxis()->SetTitle("no of entries");
        hADC_Fit_Chi2->GetXaxis()->SetTitle("avrg adc fit #chi^2");
        hADC_Fit_Chi2->GetYaxis()->SetTitle("no of entries");
        histSaver->SaveHistogram(hADCRatio,false,true,true);
        histSaver->SaveHistogram(hADC_Fit_Par0,false,true,true);
        histSaver->SaveHistogram(hADC_Fit_Par1,false,true,true);
        histSaver->SaveHistogram(hADC_Fit_Chi2);
        delete hADCRatio;
        delete hADC_Fit_Chi2;
        delete hADC_Fit_Par0;
        delete hADC_Fit_Par1;
        delete fit;
        Float_t adc,adc0;
        for(UInt_t binY = 1; binY <= prof2d->GetNbinsY();binY++){
            for(UInt_t binX=1; binX<=prof2d->GetNbinsX();binX++){
                adc=prof2d->GetBinContent(binX,binY);
                if (binX==1)
                    if (adc)
                        adc0 = adc;
                    else
                        adc0=1;
                Int_t bin = prof2d->GetBin(binX,binY);
                prof2d->SetBinContent(bin,adc/adc0);
                prof2d->SetBinEntries(bin,1);
            }
        }
        prof2d->SetName(prof2d->GetName()+(TString)"_rescaledChannels");
        histSaver->SaveHistogram(prof2d);
        cout<<"delete "<<name<<flush;
        delete hHistoMap[name];
        cout<<"."<<flush;
        hHistoMap.erase(name);
        cout<<" done."<<endl;
    }//for loop det
}

TH1F* TAnalysisOfPedestal::doSlidingWindowAnalysis(TH1D* histo,Int_t nAvrg,bool absolute) {
    if (!histo)
        return 0;
//    cout<<"do Sliding Window for '"<<histo->GetName()<<"'"<<endl;
    TString name  = histo->GetName()+TString::Format("_SlidingWindow_%dAvrg",nAvrg);
    TH1F* retHisto = (TH1F*)histo->Clone(name);
    retHisto->Reset();
    retHisto->SetTitle(name);
    Float_t meanPre=0;
    Float_t meanPast=0;
    Float_t max = -1e9;
    Float_t min = +1e9;
    for(UInt_t bin = 1; bin <= nAvrg&& bin+nAvrg <=histo->GetNbinsX(); bin++){
        meanPre+=histo->GetBinContent(bin);
        meanPast+=histo->GetBinContent(bin+nAvrg+1);
    }
    Float_t value = (meanPast-meanPre);
    if(absolute)    value = TMath::Abs(value);
    retHisto->SetBinContent(nAvrg+1,value);
    cout<< nAvrg+1<<"\t"<<retHisto->GetBinCenter(nAvrg+1)<< "\t"<<value;

    for (UInt_t bin = 2+nAvrg; bin <= histo->GetNbinsX()-nAvrg-1;bin++){
        meanPre -=  histo->GetBinContent(bin-nAvrg-1);
        meanPre +=  histo->GetBinContent(bin);
        meanPast -= histo->GetBinContent(bin);
        meanPast +=histo->GetBinContent(bin+nAvrg+1);
        value = (meanPast-meanPre);
        if(absolute)    value = TMath::Abs(value);
        retHisto->SetBinContent(bin-1,value);
    }
    return retHisto;
}

void TAnalysisOfPedestal::SetYRangeForSignalInSigmaPlot(TH1F* histo) {
    if (!histo) return;
    bool minFound=false;
    Int_t bin =1;
    Double_t content;
    Double_t max = 1e20;
    while(!minFound&&bin<histo->GetNbinsX()){
        content= histo->GetBinContent(bin);
        if(content>max)
            minFound = true;
        else max = content;
        bin++;
    }
    Int_t startbin = bin;
    max = content;
    for(;bin<histo->GetNbinsX();bin++){
        content = histo->GetBinContent(bin);
        if(content>max)max = content;
    }
    max *=1.1;
    histo->GetYaxis()->SetRangeUser(0,max);
}

void TAnalysisOfPedestal::checkCommonModeNoise(){
    UInt_t det = TPlaneProperties::getDetDiamond();
    Float_t maxVal =TPlaneProperties::getMaxSignalHeightDiamond();
    Float_t cmNoise = 0;
    UInt_t nCmNoiseEvents =0;
    vector<Float_t> channelWeight;
    for (UInt_t ch = 0; ch < TPlaneProperties::getNChannelsDiamond();ch++){
        channelWeight.push_back(0);

        Float_t adc = eventReader->getDia_ADC(ch);
        bool masked = settings->IsMasked(det,ch);
        Float_t mean = eventReader->getPedestalMean(det,ch,true);
        Float_t signal = adc - mean;
        Float_t sigma = eventReader->getPedestalSigma(det,ch,true);
        if(sigma<=0) {
            if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<ch<<" sigma < 0 : "<<sigma<<endl;
            continue;
        }
        Float_t snr = TMath::Abs(signal/sigma);

        if(snr!=snr||adc!=adc||signal!=signal)
            continue;

        if(adc>=maxVal || adc<0||signal>maxVal){
            if(verbosity>7)cout<<"cannot use "<<nEvent<<"/"<<ch<<" invalid adc/signal: "<<adc<<"/"<<signal<<endl;
            continue;
        }
        if (TMath::Abs(snr)>settings->getCMN_cut()){
            if(verbosity>7)cout<<"cannot use "<<nEvent<<"/"<<ch <<" snr over cut: "<<snr<<endl;
            continue;
        }

        if (masked){
            if(verbosity>7)cout<<"CMN: cannot use "<<nEvent<<"/"<<ch <<" Is Masked "<<endl;
            continue;
        }

        if(verbosity>10||(verbosity>4&&nEvent==0))cout<<" "<<ch<<"\t"<<adc<<" "<<mean<< " "<<sigma<<" "<<signal<<" "<<snr<<endl;
        cmNoise+=signal;
        channelWeight.back() = signal;
        nCmNoiseEvents++;
        hCmnUsedChannels->Fill(ch);
    }

    cmNoise = cmNoise/(Float_t)nCmNoiseEvents;
    Float_t relChange = (cmNoise-cmn)/cmn*100.;
    if (verbosity>6||(TMath::Abs(relChange)>100. && verbosity>1)||false)
        cout<<TString::Format("%7d - %3d %+6.2f %+6.2f - %6.1f",nEvent,nCmNoiseEvents,cmNoise,cmn,relChange)<<endl;
    if(verbosity>4)cout<<nEvent <<" cmNoise: "<<" "<<cmNoise<<" "<<nCmNoiseEvents<<" "<<eventReader->getCmnCreated(8)<<endl;
//    hCommonModeNoise->Fill(cmNoise,true);
    hRelCmnUncertainty->Fill(relChange);

    hCmnNUsedChannels->Fill(nCmNoiseEvents);
    hNewComonModeNoise->Fill(cmNoise);
    hCmnNewVsNUsedChannels->Fill(cmNoise,nCmNoiseEvents);
    for (UInt_t ch =0;ch< channelWeight.size();ch++)
        if (channelWeight.at(ch)!=0){
            Float_t weight = channelWeight.at(ch)/cmNoise*100.;
            hCmnChannelWeightVsChannel->Fill(ch,weight);
            hCmnFractionVsChannel->Fill(ch,channelWeight.at(ch));
            hCmnChannelWeight->Fill(weight);
            //cout<<nEvent<<"/"<<ch <<" weight: \t"<<channelWeight.at(ch)<< "\t"<<cmNoise<<"\t"<<weight<<endl;
        }
    hCmnVsNewCmn->Fill(cmn,cmNoise);
    hNewCmnVsEventNo->Fill(nEvent,cmNoise);
}
