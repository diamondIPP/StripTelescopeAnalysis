/*
 * TDeadChannels.cpp
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#include "../include/TAnalysisOfClustering.hh"

TAnalysisOfClustering::TAnalysisOfClustering(TSettings *newSettings) {
    cout<<"\n\n\n\n**********************************************************"<<endl;
    cout<<"**********************************************************"<<endl;
    cout<<"*********TAnalysisOfClustering::TAnalysisOfClustering*****"<<endl;
    cout<<"**********************************************************"<<endl;
    cout<<"**********************************************************\n\n\n"<<endl;
    if(newSettings==0)
        exit(-1);//todo
    //settings=new TSettings();
    setSettings(newSettings);
    res = 0;
    UInt_t runNumber=settings->getRunNumber();
    sys = gSystem;
    gROOT->ProcessLine("gErrorIgnoreLevel = 3000");
    htmlClus= new THTMLCluster(settings);

    settings->goToClusterTreeDir();
    eventReader=new TADCEventReader(settings->getClusterTreeFilePath(),settings);
    eventReader->SelectBranchesForClusteringAnalysis();
    histSaver=new HistogrammSaver(settings);


    settings->goToClusterAnalysisDir();
    stringstream plotsPath;
    plotsPath<<sys->pwd()<<"/";
    histSaver->SetPlotsPath(plotsPath.str().c_str());
    histSaver->SetRunNumber(runNumber);
    htmlClus->setFileGeneratingPath(sys->pwd());
    settings->goToClusterTreeDir();
    verbosity=settings->getVerbosity();
    initialiseHistos();
    cout<<"end initialise"<<endl;
    //	settings=0;
    vecVecClusters.resize(TPlaneProperties::getNDetectors());

    nMaxClusters = 40000;
    nInvalidReadout = 0;
}

TAnalysisOfClustering::~TAnalysisOfClustering() {
    // TODO Auto-generated destructor stub
    delete eventReader;
    delete histSaver;
    htmlClus->createClusterSize(vecClusterSizes,vecClusterSeedSizes,vecNumberOfClusters);
    htmlClus->createPulseHeightPlots(this->vecPHMeans);
    htmlClus->createContent();
    htmlClus->generateHTMLFile();
    delete htmlClus;
    settings->goToOutputDir();
}

void TAnalysisOfClustering::setSettings(TSettings* settings){
    this->settings=settings;
}

void TAnalysisOfClustering::doAnalysis(int nEvents)
{

    initPedestalAndNoiseHistos(nEvents);
    cout<<"analyze clustering results..."<<endl;
    //	eventReader->checkADC();
    if(nEvents==0) nEvents=eventReader->GetEntries();
    histSaver->SetNumberOfEvents(nEvents);
    for(nEvent=0;nEvent<nEvents;nEvent++){
        TRawEventSaver::showStatusBar(nEvent,nEvents,100);
        eventReader->LoadEvent(nEvent);
        analyseEvent();
    }
    cout<<"Save Histos!"<<endl;
    saveHistos();
}

void TAnalysisOfClustering::analyseEvent(){
//    if((nEvent != eventReader->getEvent_number()) || (nEvent != eventReader->getCurrent_event()))
//        cout<< "\nClustering analysis Event: " << int(nEvent) << ". Ev Reader Event Number: " << int(eventReader->getEvent_number()) << ". Ev Reader Current Event: " << int(eventReader->getCurrent_event()) << "\n" <<endl;
    if(eventReader->getEvent()->hasInvalidReadout()){
        //		cout<<nEvent<<": Invalid Readout..."<<endl;
        nInvalidReadout++;
        return;
    }
    checkForDeadChannels();
    checkForSaturatedChannels();
    //		getBiggestHit();//not working
    analyseForSeeds();
    analyseCluster();
    compareCentroid_ChargeWeightedMean();

    analyse2ndHighestHit();
    analyseClusterPosition();

    createPHDistribution();

    etaInvestigation();

    fillClusterVector();

    fillRelativeHitPosition();
    fillPedestalsAndNoiseHistos();
}

void TAnalysisOfClustering::fillRelativeHitPosition(){
    for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){

    }
}

void TAnalysisOfClustering::checkForDeadChannels()
{
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        int numberOfSeeds=0;
        for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
            Float_t sigma=eventReader->getPedestalSigma(det,ch);
            if(sigma==0){
                //cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
                continue;
            };
            Float_t signalInSigma=eventReader->getSignalInSigma(det,ch);

            if(signalInSigma>settings->getClusterSeedFactor(det,ch)){
                hSeedMap[det]->Fill(ch);
                //cout<<"Found a Seed "<<det<<" "<<ch<<" "<<adcValueInSigma<<" "<<eventReader->getCurrent_event()<<endl;
                numberOfSeeds++;
            }
        }
        hNumberOfSeeds[det]->Fill(numberOfSeeds);
    }

}
void TAnalysisOfClustering::analyseForSeeds(){
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        int nClusters = eventReader->getNClusters(det);
        if(nClusters==1)
            hSeedMap2[det]->Fill(eventReader->getCluster(det,0).getHighestSignalChannel());
    }
}

void TAnalysisOfClustering::checkForSaturatedChannels()
{
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
        for(UInt_t ch=0;ch<TPlaneProperties::getNChannels(det);ch++){
            if(eventReader->getAdcValue(det,ch)>=settings->getMaxSignalHeight(det)){ // DA
                hSaturatedChannels[det]->Fill(ch);
            }
        }
}


///**

void TAnalysisOfClustering::initPedestalAndNoiseHistos(UInt_t maxEvents) {
    cout<<"initPedestalHistos"<<endl;
    UInt_t start = 0;
    UInt_t nBins = (maxEvents-start)/10000;
    UInt_t diaDet = TPlaneProperties::getDetDiamond();
    for(UInt_t ch = 0; ch< TPlaneProperties::getNChannelsDiamond();ch++){
        if(settings->IsMasked(diaDet,ch))
            continue;
        TString name = TString::Format("hPedestalVsEventNo_det_%d_ch_%03d",diaDet,ch);
        TProfile* prof = new TProfile(name,name,nBins,start,maxEvents);
        prof->GetXaxis()->SetTitle("EventNo");
        TString title = TString::Format("pedestal_{ch %03d} /ADC",ch);
        if(settings->doCommonModeNoiseCorrection()) title.Append(" CM corrected");
        prof->GetYaxis()->SetTitle(title);
        hPedestalVsEvenNo[ch] = prof;

        //Noise of each channel
        name = TString::Format("hNoiseVsEventNo_det_%d_ch_%03d",diaDet,ch);
        prof = new TProfile(name,name,nBins,start,maxEvents);
        prof->GetXaxis()->SetTitle("EventNo");
        title = TString::Format("noise_{ch %03d} /ADC",ch);
        if(settings->doCommonModeNoiseCorrection()) title.Append(" CM corrected");
        prof->GetYaxis()->SetTitle(title);
        title =  TString::Format("noise of ch %03d vs. EvenNot",ch);
        if(settings->doCommonModeNoiseCorrection()) title.Append(" CM corrected");
        hNoiseVsEvenNo[ch] = prof;

        // raw ADC of each channel: hRawAdcVsEvenNo
        name = TString::Format("hADCVsEventNo_det_%d_ch_%03d",diaDet,ch);
        prof = new TProfile(name,name,nBins,start,maxEvents);
        prof->GetXaxis()->SetTitle("EventNo");
        title = TString::Format("raw ADC_{ch %03d} /ADC",ch);
        prof->GetYaxis()->SetTitle(title);
        title =  TString::Format("raw ADC of ch %03d vs. EvenNot",ch);
        hADCVsEvenNo[ch] = prof;
    }
    TString name = "hComonModeNoiseVsEventNo";
    hCmnVsEventNo = new TProfile(name,name,nBins,start,maxEvents);
    hCmnVsEventNo->GetXaxis()->SetTitle("EventNo");
    hCmnVsEventNo->GetYaxis()->SetTitle("common mode noise /ADC");
}

void TAnalysisOfClustering::fillPedestalsAndNoiseHistos() {
    UInt_t diaDet = TPlaneProperties::getDetDiamond();
    std::map<UInt_t,TProfile*>::iterator it;
    for(it=hPedestalVsEvenNo.begin(); it!=hPedestalVsEvenNo.end(); it++){
        UInt_t channel = (*it).first;
        Float_t pedestal = eventReader->getPedestalMean(diaDet,channel,settings->doCommonModeNoiseCorrection());
        (*it).second->Fill(nEvent,pedestal);
    }
    for(it=hNoiseVsEvenNo.begin(); it!=hNoiseVsEvenNo.end(); it++){
        UInt_t channel = (*it).first;
        Float_t noise = eventReader->getPedestalSigma(diaDet,channel,settings->doCommonModeNoiseCorrection());
        (*it).second->Fill(nEvent,noise);
    }
    //hRawAdcVsEvenNo

    for(it=hADCVsEvenNo.begin(); it!=hADCVsEvenNo.end(); it++){
        UInt_t channel = (*it).first;
        Float_t adc = eventReader->getAdcValue(diaDet,channel);
        (*it).second->Fill(nEvent,adc);
    }
    hCmnVsEventNo->Fill(nEvent,eventReader->getCMNoise());
}

void TAnalysisOfClustering::saveADCHistos() {
    saveVariableVsEventNoPlots(settings,histSaver,hADCVsEvenNo,"ADC",&vecRawADCSlope,&vecCh);
    return;
}

void TAnalysisOfClustering::saveNoiseHistos() {
    saveVariableVsEventNoPlots(settings,histSaver,hNoiseVsEvenNo,"Noise",&vecNoiseSlope,&vecCh);
    return;
}

void TAnalysisOfClustering::analysisSlopes(){
    TH1F* hPedestalSlopes = histSaver->CreateDistributionHisto("hPedestalSlopes",vecPedestalSlope,10);
    hPedestalSlopes->GetXaxis()->SetTitle("pedestal slope ADC/Event");
    hPedestalSlopes->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hPedestalSlopes);

    TH1F* hNoiseSlopes = histSaver->CreateDistributionHisto("hNoiseSlopes",vecNoiseSlope,10);
    hNoiseSlopes->GetXaxis()->SetTitle("noise slope ADC/Event");
    hNoiseSlopes->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hNoiseSlopes);

    TH1F* hADCSlopes = histSaver->CreateDistributionHisto("hADCSlopes",vecRawADCSlope,10);
    hADCSlopes->GetXaxis()->SetTitle("raw ADC slope ADC/Event");
    hADCSlopes->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hADCSlopes);

    TGraph gPedestalSlopeVsRawADCSlope = histSaver->CreateDipendencyGraph("gPedestalSlopeVsRawADCSlope",vecPedestalSlope,vecRawADCSlope);
    gPedestalSlopeVsRawADCSlope.GetXaxis()->SetTitle("slope of raw adc");
    gPedestalSlopeVsRawADCSlope.GetYaxis()->SetTitle("slope of pedestal");
    histSaver->SaveGraph(&gPedestalSlopeVsRawADCSlope,"gPedestalSlopeVsRawADCSlope");

    delete hNoiseSlopes;
    delete hPedestalSlopes;
    delete hADCSlopes;
}

void TAnalysisOfClustering::saveVariableVsEventNoPlots(TSettings* settings,HistogrammSaver*histSaver, std::map<UInt_t,TProfile*> mProf,
        TString nameOfVariable,vector <Float_t>* vec, vector<Float_t> *vecCh){
    TString name = (TString)"stack"+nameOfVariable+(TString)"VsEventNo";
    THStack * stack = new THStack(name,name);
    name = (TString)"h"+nameOfVariable+(TString)"SlopesVsChannel";
    TH1F* hVariableSlopesVsChannel = new TH1F(name,name,128,0,128);
    hVariableSlopesVsChannel->GetXaxis()->SetTitle("channel no");
    hVariableSlopesVsChannel->GetYaxis()->SetTitle("slope m = ADC/Event");
    UInt_t color = 0;
    std::map<UInt_t,TProfile*>::iterator it;
    TF1* pol1 = new TF1("pol1_fit","pol1",0,5e6);
    pol1->SetLineColor(kBlue);
    pol1->SetLineWidth(1);
    Double_t minStack = 1e9;
    Double_t maxStack = -1e9;
    name = (TString)"hRel"+nameOfVariable+(TString)"PedestalVsEventNo";
    it=mProf.begin();
    if(it==mProf.end()){
        cerr<<"Profile Map is empty..."<<endl;
        return;
    }
    TProfile* prof = (*it).second;
    if(!prof){
        cerr<<"Profile of iterator "<< (*it).first <<" is nullpointer"<<endl;
        return;
    }
    Float_t maxEvents =prof->GetXaxis()->GetXmax();
    Float_t start =prof->GetXaxis()->GetXmin();
    UInt_t nBins =prof->GetXaxis()->GetNbins();
    TH2D* hRelVariableVsEventNo = new TH2D(name,name,nBins,start,maxEvents,128,0,128);
    hRelVariableVsEventNo->GetXaxis()->SetTitle("Event No.");
    hRelVariableVsEventNo->GetYaxis()->SetTitle("channel No.");
    hRelVariableVsEventNo->GetZaxis()->SetTitle("relative change to pedestal_{0} in percent %% ");
    int i = 0;
    //Loop over all Histograms in Map
    for(it=mProf.begin(); it!=mProf.end(); it++){
        TProfile* prof = (*it).second;
        UInt_t channel = (*it).first;
        if(!prof) continue;
        TF1* fit = (TF1*)pol1->Clone(prof->GetName()+(TString)"_fit");

        if(channel%5==0){//add clone to stack if channel no % 5 ==0
            TProfile* prof2 =(TProfile*)prof->Clone();
            prof2->SetTitle(TString::Format("Channel %3d",channel));
            prof2->SetLineColor(color);
            prof2->SetMarkerColor(color);
            stack->Add(prof2);
            minStack = TMath::Min( minStack, prof->GetBinContent(prof->GetMinimumBin()));
            maxStack = TMath::Max( maxStack, prof->GetBinContent(prof->GetMaximumBin()));
            color++;
        }

        TProfile * prof3 = (TProfile*)prof->Clone();

        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        hVariableSlopesVsChannel->SetBinContent(hVariableSlopesVsChannel->FindBin(channel),fit->GetParameter(1));
        if ((int)vecCh->size()<=i){
            vecCh->push_back(channel);
            i++;
        }
        //        else if(vecCh.at(i)!=channel){
        //            cerr<<"something is wrong: channel: "<<channel<<" vecChannel @ "<<i<<": "<<vecCh.at(i)<<endl;
        //            return;
        //        }
        vec->push_back(fit->GetParameter(1));
        Float_t firstVarableValue = prof3->GetBinContent(1);
        if (false) cout<<channel<<" first"<<nameOfVariable<<"Value: "<<firstVarableValue<<endl;
        prof3->Scale(1./firstVarableValue);
        for(Int_t binx = 1;binx <= prof3->GetNbinsX();binx++){
            Float_t content = prof3->GetBinContent(binx);
            content = (content-1.)*100;
            Int_t bin = hRelVariableVsEventNo->GetBin(binx,channel+1);
            //            cout<<"Fill: "<<channel<<" "<<binx<<" "<<content<<endl;
            hRelVariableVsEventNo->SetBinContent(bin,content);
        }
        delete prof3;
        delete prof;
        (*it).second= 0;
        mProf.erase(it);
    }
    histSaver->SaveHistogram(hRelVariableVsEventNo,true);
    name = (TString)"g"+nameOfVariable+(TString)"SlopeVsChannel";
    TGraph graph = histSaver->CreateDipendencyGraph((string)name,*vec,*vecCh);
    graph.Draw("AP");
    graph.GetXaxis()->SetTitle("channel");
    graph.GetYaxis()->SetTitle(nameOfVariable+ (TString)" slope for channel");
    histSaver->SaveGraph(&graph,(string)name,"ABP");


    if(color!=0){
        cout<<"save stack "<<minStack<<"-"<<maxStack<<endl;
        stack->Draw("goff");
        stack->SetObjectStat(false);
        if(stack->GetXaxis()){
            stack->GetXaxis()->SetTitle("Event No");
            cout<<"Xaxis: "<<stack->GetXaxis()->GetTitle()<<endl;
        }
        if(stack->GetYaxis()){
            stack->GetYaxis()->SetTitle(nameOfVariable+ (TString)" /ADC");
            cout<<"Set range"<<endl;
            stack->GetYaxis()->SetRangeUser(minStack*.98,maxStack*1.05);
            cout<<"Yaxis: "<<stack->GetYaxis()->GetTitle()<<endl;
        }
        stack->SetMinimum(minStack*.98);
        stack->SetMaximum(maxStack*1.05);
        stack->SetObjectStat(false);
    }
    histSaver->SaveStack(stack,"nostack",true);
    cout<<"h"<<nameOfVariable<<"SlopesVsChannel:MIN: "<<hVariableSlopesVsChannel->GetBinContent(hVariableSlopesVsChannel->GetMinimumBin())<<endl;
    hVariableSlopesVsChannel->SetMinimum(hVariableSlopesVsChannel->GetBinContent(hVariableSlopesVsChannel->GetMinimumBin()));
    histSaver->SaveHistogram(hVariableSlopesVsChannel,false,false);
    delete hVariableSlopesVsChannel;
    delete stack;
}

void TAnalysisOfClustering::savePedestalHistos() {
    saveVariableVsEventNoPlots(settings,histSaver,hPedestalVsEvenNo,"Pedestal",&vecPedestalSlope,&vecCh);
    return;
}

///***

void TAnalysisOfClustering::initialiseHistos()
{
    if(verbosity>3)cout<<"1"<<endl;
    {
        stringstream histoName;
        histoName<<"hDiamond_Delta_CWM_BiggestHit";
        histo_CWM_biggestHit=new TH2F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6,10,0,9);
        histoName.str("");
        histoName<<"hDiamond_Delta_highest2Centroid_BiggestHit";
        histo_H2C_biggestHit=new TH1F(histoName.str().c_str(),histoName.str().c_str(),512,-0.6,0.6);
    }
    size_t nDet = TPlaneProperties::getNDetectors();
    vecvecSignalLeftLeft.resize(nDet);
    vecvecSignalRightRight.resize(nDet);
    vecvecLeftEtaSignal.resize(nDet);
    vecvecRightEtaSignal.resize(nDet);
    vecvecSignalLeftOfHighest.resize(nDet);
    vecvecSignalRightOfHighest.resize(nDet);
    vecvecSignalHighest.resize(nDet);
    vecvecEta.resize(nDet);

    TString hName = "hEtaDistributionDia";
    int nDiamonds = settings->getNDiamonds();
    hEtaDistributionDia = new TH2F(hName,hName,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),nDiamonds+2,-0.5,nDiamonds+1.5);

    for(UInt_t det=0;det<9;det++){
        TString name = "hClusterPositionRelativeToNextIntegerCWM_"+(TString)TPlaneProperties::getStringForDetector(det);
        if (verbosity>2) cout<<name<<endl;
        if(TPlaneProperties::isDiamondDetector(det)){
            pair<Int_t,Int_t> pattern = settings->diamondPattern.getTotalInterval();
            hRelativeClusterPositionCWM[det]=new TH2F(name,name,pattern.second-pattern.first,pattern.first,pattern.second,1024,-.5,.5);
        }
        else
            hRelativeClusterPositionCWM[det]=new TH2F(name,name,256,0,TPlaneProperties::getNChannels(det)-1,1024,-.5,.5);

        name = "hClusterPositionRelativeToNextIntegerCorEta_"+ (TString)TPlaneProperties::getStringForDetector(det);
        if (verbosity>2) cout<<name<<endl;
        if(TPlaneProperties::isDiamondDetector(det)){
            pair<Int_t,Int_t> pattern = settings->diamondPattern.getTotalInterval();
            hRelativeClusterPositionCorEta[det]=new TH2F(name,name,pattern.second-pattern.first,pattern.first,pattern.second,1024,-.5,.5);
        }
        else
            hRelativeClusterPositionCorEta[det]=new TH2F(name,name,256,0,TPlaneProperties::getNChannels(det)-1,1024,-.5,.5);

        name = "hClusterPositionRelativeToNextIntegerEta_"+(TString)TPlaneProperties::getStringForDetector(det);
        if (verbosity>2) cout<<name<<endl;
        if(TPlaneProperties::isDiamondDetector(det)){
            pair<Int_t,Int_t> pattern = settings->diamondPattern.getTotalInterval();
            hRelativeClusterPositionEta[det]=new TH2F(name,name,pattern.second-pattern.first,pattern.first,pattern.second,1024,-.5,.5);
        }
        else
            hRelativeClusterPositionEta[det]=new TH2F(name,name,256,0,TPlaneProperties::getNChannels(det)-1,1024,-.5,.5);

        name = "hAbsoluteClusterPostion_"+(TString)TPlaneProperties::getStringForDetector(det);;
        if (verbosity>2) cout<<name<<endl;
        hClusterPosition[det]=new TH1F(name,name,TPlaneProperties::getNChannels(det)*16,0,TPlaneProperties::getNChannels(det)-1);

        if(TPlaneProperties::isDiamondDetector(det)){
            name = "hAbsoluteClusterPostionActiveChannels_"+(TString)TPlaneProperties::getStringForDetector(det);
            if (verbosity>2) cout<<name<<endl;
            Float_t min = settings->diamondPattern.getTotalInterval().first;
            Float_t max = settings->diamondPattern.getTotalInterval().second;
            Int_t bins = (max - min)*20;
            hClusterPositionActiveChannels=new TH1F(name,name, bins, min,max);
        }

        name = "hEtaDistribution_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistribution[det]=new TH1F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS));
        hEtaDistribution[det]->GetXaxis()->SetTitle("#eta");
        hEtaDistribution[det]->GetYaxis()->SetTitle("number of entries");

        name = "hEtaDistributionCMN_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistributionCMN[det]=new TH1F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS));
        hEtaDistributionCMN[det]->GetXaxis()->SetTitle("#eta_{CMN-corrected}");
        hEtaDistributionCMN[det]->GetYaxis()->SetTitle("number of entries");

        name = "hEtaDistributionVsLeftChannel_"+ (TString) TPlaneProperties::getStringForDetector(det);
        if(TPlaneProperties::isDiamondDetector(det)){
            pair<Int_t,Int_t> pattern = settings->diamondPattern.getTotalInterval();
            hEtaDistributionVsLeftChannel[det] = new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),pattern.second-pattern.first,pattern.first,pattern.second);
        }
        else
            hEtaDistributionVsLeftChannel[det] = new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),256,0,256);
        hEtaDistributionVsLeftChannel[det]->GetXaxis()->SetTitle("#eta");
        hEtaDistributionVsLeftChannel[det]->GetYaxis()->SetTitle("left channel of #eta position");
        hEtaDistributionVsLeftChannel[det]->GetZaxis()->SetTitle("number of entries #");

        name = "hEtaDistributionVsClusterSize_"+ (TString) TPlaneProperties::getStringForDetector(det);

        hEtaDistributionVsClusterSize[det] = new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),10,-.5,9.5);
        hEtaDistributionVsClusterSize[det]->GetXaxis()->SetTitle("#eta");
        hEtaDistributionVsClusterSize[det]->GetYaxis()->SetTitle("ClusterSize");
        hEtaDistributionVsClusterSize[det]->GetZaxis()->SetTitle("number of entries #");

        name = "hEtaDistributionVsCharge_"+ (TString) TPlaneProperties::getStringForDetector(det);
        Int_t maxCharge = TPlaneProperties::isDiamondDetector(det)?settings->getPulse_height_di_max():settings->getPulse_height_si_max();
        hEtaDistributionVsCharge[det] = new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),512,0,maxCharge);
        hEtaDistributionVsCharge[det]->GetXaxis()->SetTitle("#eta");
        hEtaDistributionVsCharge[det]->GetYaxis()->SetTitle("Charge of two highest Channels /ADC counts");
        hEtaDistributionVsCharge[det]->GetYaxis()->SetTitle("number of entries");

        name = "hEtaDistribution5Percent_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistribution5Percent[det]=new TH1F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS));

        name = "hEtaDistributionVsSignalLeft_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistributionVsSignalLeft[det]=new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),128,0,settings->getMaxSignalHeight(det)); // DA

        name = "hEtaDistributionVsSignalRight_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistributionVsSignalRight[det]=new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),128,0,settings->getMaxSignalHeight(det)); // DA

        name = "hEtaDistributionVsSignalSum_"+ (TString) TPlaneProperties::getStringForDetector(det);
        hEtaDistributionVsSignalSum[det]=new TH2F(name,name,ETA_BINS+1,0-1.0/(2*ETA_BINS),1+1.0/(2*ETA_BINS),128,0,settings->getMaxSignalHeight(det)*2); // DA

        name = "hSignalLeftVsSignalRight"+ (TString) TPlaneProperties::getStringForDetector(det);
        hSignalLeftVsSignalRight[det]=new TH2F(name,name,128,0,settings->getMaxSignalHeight(det),128,0,settings->getMaxSignalHeight(det)); // DA

        name = "hDeltaLeftRightVsMaximum"+ (TString) TPlaneProperties::getStringForDetector(det);
        hDeltaLeftRightVsMaximum[det]=new TH2F(name,name,
                1024,-settings->getMaxSignalHeight(det),settings->getMaxSignalHeight(det),
                256,0,settings->getMaxSignalHeight(det)); // DA
//        hLeftVsMaximum[det]

        name = "hLeftVsMaximum"+ (TString) TPlaneProperties::getStringForDetector(det);
        hLeftVsMaximum[det]=new TH2F(name,name,
                1024,-settings->getMaxSignalHeight(det),settings->getMaxSignalHeight(det),
                256,0,settings->getMaxSignalHeight(det)); // DA

        name = "hRightVsMaximum"+ (TString) TPlaneProperties::getStringForDetector(det);
        hRightVsMaximum[det]=new TH2F(name,name,
                1024,-settings->getMaxSignalHeight(det),settings->getMaxSignalHeight(det),
                256,0,settings->getMaxSignalHeight(det)); // DA
    }
    //	for(int det = 0; det < 9; det++){
    //	    TString name = (TString)"hRelativeHitPosition"+(TString)TPlaneProperties::getStringForDetector(det);
    //	    hRelativeHitPosition[det] = new TH1F(name,name,512,-30,30);
    //	    hRelativeHitPosition[det]->GetXaxis()->SetTitle("relative Hit Position/#mum");
    //	    hRelativeHitPosition[det]->GetYaxis()->SetTitle("number of entries");
    //	}
    if(verbosity>3)cout<<"2"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"SaturatedChannels_"<<TPlaneProperties::getStringForDetector(det)<<"";
        hSaturatedChannels[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,256);
        if(det==8)hSaturatedChannels[det]->GetXaxis()->SetRangeUser(0,128);
    }
    if(verbosity>3)cout<<"3"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"hPositionOfallSeeds_"<<TPlaneProperties::getStringForDetector(det);
        hSeedMap[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,256);
        if(det==8)hSeedMap[det]->GetXaxis()->SetRangeUser(0,128);
    }
    if(verbosity>3)cout<<"4"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"hPositionOfHighestSeed_"<<TPlaneProperties::getStringForDetector(det);
        hSeedMap2[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,256);
        hSeedMap2[det]->GetXaxis()->SetTitle("Position of Highest Seed of a Cluster");
        if(det==8)hSeedMap2[det]->GetXaxis()->SetRangeUser(0,128);
    }
    if(verbosity>3)cout<<"5"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"hNumberOfSeeds_in_"<<TPlaneProperties::getStringForDetector(det);
        hNumberOfSeeds[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),31,0,30);
        hNumberOfSeeds[det]->GetXaxis()->SetTitle("Number Of Seeds in Cluster");
        hNumberOfSeeds[det]->GetYaxis()->SetTitle("Entries #");
    }
    if(verbosity>3)cout<<"6"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"PulseHeight_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHitChannelInSigma";
        hPulsHeightBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
    }
    if(verbosity>3)cout<<"7"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"PulseHeight_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHitNextToBiggestHit_ChannelInSigma";
        hPulsHeightNextBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),4000,0,400);
    }
    if(verbosity>3)cout<<"8"<<endl;
    for (int det=0;det<9;det++){
        stringstream histoName;
        histoName<<"Channel_"<<TPlaneProperties::getStringForDetector(det)<<"_BiggestHit";
        hChannelBiggestHit[det]=new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0,256);
    }
    for(int det=0;det<9;det++){
        Float_t seed = settings->getClusterSeedFactor(det,0);
        Float_t hit = settings->getClusterHitFactor(det,0);
        TString name = TString::Format("hClusterSize_Seed%d_Hit%d_",(int)seed,(int)hit)+(TString)TPlaneProperties::getStringForDetector(det);
        hClusterSize[det]= new TH1F(name,name,10,-0.5,10.5);
        hClusterSize[det]->GetXaxis()->SetTitle("Number of Seeds and Hits in Cluster");
        hClusterSize[det]->GetYaxis()->SetTitle("Entries #");
        name = TString::Format("hClusterSeedSize_Seed%d_Hit%d_",(int)seed,(int)hit)+(TString)TPlaneProperties::getStringForDetector(det);
        hClusterSeedSize[det]= new TH1F(name,name,10,-0.5,10.5);
        hClusterSeedSize[det]->GetXaxis()->SetTitle("Number of Seeds in Cluster");
        hClusterSeedSize[det]->GetYaxis()->SetTitle("Entries #");
        name = (TString)"NumberOfClusters_"+(TString)TPlaneProperties::getStringForDetector(det);
        hNumberOfClusters[det]= new TH1F(name,name,10,-0.5,10.5);
        name = (TString)"hClusterPos_"+(TString)TPlaneProperties::getStringForDetector(det);
        Int_t min = 0;
        Int_t max = TPlaneProperties::getNChannels(det);
        Int_t bins = max-min;
        hClusterPos[det] = new TH1F(name,name,bins,min,max);
        name = (TString)"hClusterPosCMN_"+(TString)TPlaneProperties::getStringForDetector(det);
        hClusterPosCMN[det] = new TH1F(name,name,bins,min,max);
    }
    if(verbosity>3)cout<<"10"<<endl;
    for (int det = 0; det < 9; det++) {
        int nbins = 250;
        Float_t min = 0.;
        Float_t max = 250.;

        stringstream histoName;
        histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "BiggestHitChannelInSigma";
        histo_pulseheight_sigma[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigma";
        histo_pulseheight_sigma_second[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitMinusBiggestHitPosition";
        histo_second_biggest_hit_direction[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),2,-2.,2.);

        histoName.str("");
        histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigmaLeft";
        histo_pulseheight_sigma_second_left[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "SecondBiggestHitChannelInSigmaRight";
        histo_pulseheight_sigma_second_right[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),nbins,min,max);

        histoName.str("");
        histoName << TPlaneProperties::getStringForDetector(det) << "BiggestHitMap";
        histo_biggest_hit_map[det] = new TH1F(histoName.str().c_str(),histoName.str().c_str(),256,0.,256.);

        histoName.str("");
        histoName << "PulseHeight" << TPlaneProperties::getStringForDetector(det) << "LeftChipBiggestHitChannelInSigma";
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
    }
    if(verbosity>3)cout<<"11"<<endl;
    for(int det=0;det<9;det++){//analayse2ndHighestHit
        stringstream histName;
        histName<<"h2ndBiggestHitSignal_"<<TPlaneProperties::getStringForDetector(det);
        if(det<8)
            h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,200);
        else
            h2ndBiggestHitSignal[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1024);
        h2ndBiggestHitSignal[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster");
        h2ndBiggestHitSignal[det]->GetYaxis()->SetTitle("Entries #");
        histName.str("");
        histName<<"h2ndBiggestHitOverCharge_"<<TPlaneProperties::getStringForDetector(det);
        h2ndBiggestHitOverCharge[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,0.5);
        h2ndBiggestHitOverCharge[det]->GetXaxis()->SetTitle("Signal of 2nd Biggest Hit of Cluster over Sum of all signals of cluster");
        h2ndBiggestHitOverCharge[det]->GetYaxis()->SetTitle("Entries #");
        histName.str("");
        histName<<"h2ndBiggestHitPosition_"<<TPlaneProperties::getStringForDetector(det);
        h2ndBiggestHitPosition[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),3,-1.5,1.5);
        h2ndBiggestHitPosition[det]->GetXaxis()->SetTitle("position of snd biggest hit in respect to biggest Hit");
        h2ndBiggestHitPosition[det]->GetYaxis()->SetTitle("Entries #");
        histName.str("");
        histName<<"hLeftHitOverLeftAndRight_"<<TPlaneProperties::getStringForDetector(det);
        hLeftHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
        hLeftHitOverLeftAndRight[det]->GetXaxis()->SetTitle("Q_L/(Q_R +Q_L)");
        histName.str("");
        histName<<"hDeltaLeftRightHitOverLeftAndRight_"<<TPlaneProperties::getStringForDetector(det);
        hDeltaLeftRightHitOverLeftAndRight[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,-1,1);
        hDeltaLeftRightHitOverLeftAndRight[det]->GetXaxis()->SetTitle("(Q_L-Q_R)/(Q_R +Q_L)");
        histName.str("");
        histName<<"hSignal2ndHighestOverSignalHighest_"<<TPlaneProperties::getStringForDetector(det);
        hSignal2ndHighestOverSignalHighestRatio[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),512,0,1);
        hSignal2ndHighestOverSignalHighestRatio[det]->GetXaxis()->SetTitle("Q_{2ndHighest}/Q_{Highest}");
    }
    if(verbosity>3)cout<<"12"<<endl;
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        stringstream histName;
        histName<<"hPulseHeightDistribution_"<<TPlaneProperties::getStringForDetector(det);
        float max=0;
        if(det==TPlaneProperties::getDetDiamond())
            max = settings->getPulse_height_di_max();
        else max = settings->getPulse_height_si_max();
        hPHDistribution[det]=new TH2F(histName.str().c_str(),histName.str().c_str(),512,0,max,10,-.5,9.5);
        histName.str("");
        histName<<"hPHDistIncreasingClustersize_"<<TPlaneProperties::getStringForDetector(det);
        hPHDistIncreasingClustersize[det] = new TH2F(histName.str().c_str(),histName.str().c_str(),
                                                    512,0,max,4,1,4);

    }

    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        stringstream name;
        name<<"hBiggestHitSizeVsClusterSize_"<<TPlaneProperties::getStringForDetector(det);
        hBiggestHitVsClusterSize[det] = new TH2F(name.str().c_str(),name.str().c_str(),1024,0,settings->getMaxSignalHeight(det)*2,8,0.5,8.5); // DA
        hBiggestHitVsClusterSize[det]->GetXaxis()->SetTitle("Signal of Biggest Hit in Cluster [adc counts]");
        hBiggestHitVsClusterSize[det]->GetYaxis()->SetTitle("ClusterSize");
    }
}



void TAnalysisOfClustering::saveEtaIntegrals(){
    stringstream etaCorFileName;
    etaCorFileName<<"etaCorrection."<<settings->getRunNumber()<<".root";
    TFile* file = TFile::Open(etaCorFileName.str().c_str());
    if(file){
        if(verbosity>4)cout<<"can read file: '"<<etaCorFileName.str()<<"'"<<endl;
        if(verbosity>4&&verbosity%2==1){
            cout<<"Press a key and enter to confirm.\t"<<flush;
            char t; cin>>t;
        }
        if(verbosity>4)cout<<"Close file"<<endl;
        file->Close();
        if(verbosity>4)cout<<"Return"<<endl;

        return;
    }
    if(verbosity>4)cout<<"cannot read file: '"<<etaCorFileName.str()<<"' ===> CREATE new Eta file"<<endl;
    if(verbosity>4&&verbosity%2==1){
        cout<<"Press a key and enter to confirm.]\t"<<flush;
        char t; cin>>t;
    }
    if(verbosity>4)cout <<"RECREATE file..."<<endl;
    file = new TFile(etaCorFileName.str().c_str(),"RECREATE");
    file->cd();
    for(UInt_t det=0;det<9;det++){
        stringstream histName;
        histName<<"hEtaIntegral_"<<det;
        if(!hEtaDistribution[det])
            continue;
        TH1F *histo= TClustering::createEtaIntegral(hEtaDistribution[det],histName.str());
        file->cd();
        histo->Write();
        hEtaDistribution[det]->Write();
    }
    hEtaDistributionDia->Write();
    file->Close();
}



void TAnalysisOfClustering::saveHistos(){
    analyseAsymmetricSample();
    //	char t; cin>>t;
    cout<<"1"<<flush;
    if (verbosity>2) cout<<"plot histo "<<histo_CWM_biggestHit->GetName();
    histSaver->SaveHistogram(histo_CWM_biggestHit);
    histo_CWM_biggestHit->Delete();
    if (verbosity>2) cout<<"plot histo "<<histo_H2C_biggestHit->GetName();
    histSaver->SaveHistogram(histo_H2C_biggestHit);
    histo_H2C_biggestHit->Delete();
    for(int det=0;det<9;det++){//analyse 2nd biggest Hit
        if (verbosity>2) cout<<"plot histo "<<det<<"  h2ndBiggestHitSignal_"<<TPlaneProperties::getStringForDetector(det);
        histSaver->SaveHistogram(h2ndBiggestHitSignal[det]);
        delete h2ndBiggestHitSignal[det];
        if (verbosity>2) cout<<"plot histo "<<det<<"  h2ndBiggestHitOverCharge_"<<TPlaneProperties::getStringForDetector(det);
        histSaver->SaveHistogram(h2ndBiggestHitOverCharge[det]);
        delete h2ndBiggestHitOverCharge[det];
        if (verbosity>2) cout<<"plot histo "<<h2ndBiggestHitPosition[det]->GetName()<<endl;
        histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
        histSaver->SaveHistogram(h2ndBiggestHitPosition[det]);
        histSaver->SaveHistogram(hLeftHitOverLeftAndRight[det]);
        histSaver->SaveHistogram(hDeltaLeftRightHitOverLeftAndRight[det]);
        histSaver->SaveHistogram(hSignal2ndHighestOverSignalHighestRatio[det]);
    }
    cout<<2<<flush;
    for (int det=0;det<9;det++){
        if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSaturatedChannels[det]->GetName()<<endl;
        histSaver->SaveHistogram(hSaturatedChannels[det]);
        hSaturatedChannels[det]->Delete();
    }
    for (int det=0;det<9;det++){
        if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap[det]->GetName()<<endl;
        histSaver->SaveHistogram(hSeedMap[det]);
        hSeedMap[det]->Delete();
    }
    cout<<3<<flush;
    for (int det=0;det<9;det++){
        if (verbosity>2) cout<<"plot histo"<<det<<" "<<hSeedMap2[det]->GetName()<<endl;
        histSaver->SaveHistogram(hSeedMap2[det]);
        hSeedMap2[det]->Delete();
    }
    for (int det=0;det<9;det++){
        if (verbosity>2) cout<<"plot histo"<<det<<" "<<hNumberOfSeeds[det]->GetName()<<endl;
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
        if (verbosity>2) cout<<"save: "<<hClusterSeedSize[det]->GetName()<<endl;
        histSaver->SaveHistogram(this->hClusterSeedSize[det]);
        histSaver->SaveHistogram(this->hNumberOfClusters[det]);
        vecClusterSizes.push_back(hClusterSize[det]->GetMean());
        vecClusterSeedSizes.push_back(hClusterSeedSize[det]->GetMean());
        vecNumberOfClusters.push_back(hNumberOfClusters[det]->GetMean());
        histSaver->SaveHistogram(hClusterPos[det]);
        histSaver->SaveHistogram(hClusterPosCMN[det]);
        delete hClusterPosCMN[det];
        delete hClusterPos[det];
        delete hClusterSize[det];
        delete hClusterSeedSize[det];
        delete hNumberOfClusters[det];
    }
    histSaver->SaveHistogram(hClusterPositionActiveChannels);
    delete hClusterPositionActiveChannels;
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        if (verbosity>2) cout<<"Print : "<<hClusterPosition[det]->GetTitle()<< " "<<hClusterPosition[det]->GetEntries()<<endl;
        histSaver->SaveHistogram(this->hClusterPosition[det]);
        histSaver->SaveHistogram(this->hRelativeClusterPositionCWM[det]);
        histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionCWM[det]->ProjectionY());
        //		histSaver->SaveHistogram(this->hRelativeClusterPositionCorEta[det]);
        histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionCorEta[det]->ProjectionY());
        histSaver->SaveHistogram((TH1F*)this->hRelativeClusterPositionEta[det]->ProjectionY());
        delete hClusterPosition[det];
        delete hRelativeClusterPositionCWM[det];
        delete hRelativeClusterPositionEta[det];
        delete hRelativeClusterPositionCorEta[det];
    }
    saveEtaIntegrals();
    if (verbosity) cout<<"create Eta Integrals"<<endl;
    histSaver->SaveHistogram(hEtaDistributionDia);
    saveEtaPerArea();
    for(UInt_t det=0;det<9;det++){
        stringstream histName;
        histName<<"hEtaIntegral_"<<TPlaneProperties::getStringForDetector(det);;
        TH1F *histo= TClustering::createEtaIntegral(hEtaDistribution[det],histName.str());
        if (verbosity) cout<<"Save: "<<histName.str()<<endl;
        histSaver->SaveHistogram(histo);
        //		if(histo) delete histo;
        if (verbosity) cout<<hEtaDistribution[det]->GetName()<<endl;
        histSaver->SaveHistogram(this->hEtaDistribution[det]);
        if (verbosity) cout<<hEtaDistributionCMN[det]->GetName()<<endl;
        histSaver->SaveHistogram(this->hEtaDistributionCMN[det]);
        if (verbosity) cout<<hEtaDistributionVsCharge[det]->GetName()<<endl;
        histSaver->SaveHistogram(this->hEtaDistributionVsCharge[det]);
        for(UInt_t area = 0; area < settings->getNDiaDetectorAreas() && TPlaneProperties::isDiamondDetector(det); area++){
            if (verbosity) cout<< "Save Eta Distributions plots for area "<<area<<endl;
            TString name = TString::Format("hEtaDistributionVsLeftChannel_%d_Area%d",det,area);
            TH2F *hEtaDistributionVsLeftChannelArea = (TH2F*)hEtaDistributionVsLeftChannel[det]->Clone();
            Float_t yMin = settings->getDiaDetectorArea(area).first;
            Float_t yMax = settings->getDiaDetectorArea(area).second;
            hEtaDistributionVsLeftChannelArea->GetYaxis()->SetRangeUser(yMin,yMax);
            histSaver->SaveHistogram(hEtaDistributionVsLeftChannelArea);
            if(hEtaDistributionVsLeftChannelArea)
                delete hEtaDistributionVsLeftChannelArea;
        }
        if (verbosity)  cout<<"Save "<< hEtaDistributionVsLeftChannel[det]->GetName()<<endl;
        histSaver->SaveHistogram(this->hEtaDistributionVsLeftChannel[det]);
        //		if(hEtaDistributionVsLeftChannel) delete hEtaDistributionVsLeftChannel;

        if (verbosity)  cout<<"Save ClusterSize Plotss"<<endl;
        for(int i=1;(i<6&&TPlaneProperties::isDiamondDetector(det))||i<3;i++){
            TString name =  TString::Format("hEtaDistribution_ClusterSize%d",i);
            Int_t bin = hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i);
            TH1F* hEtaOneClusterSize = (TH1F*)hEtaDistributionVsClusterSize[det]->ProjectionX(name,bin,bin);
            hEtaOneClusterSize->SetTitle(name);
            histSaver->SaveHistogram(hEtaOneClusterSize);
            delete hEtaOneClusterSize;
            if(i==1){
                name =  TString::Format("hEtaDistribution_ClusterSize_%d-%d",i,i+1);
                bin = hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i);
                Int_t bin2 =hEtaDistributionVsClusterSize[det]->GetYaxis()->FindBin(i+1);
                TH1F* hEtaTwoClusterSize = (TH1F*)hEtaDistributionVsClusterSize[det]->ProjectionX(name,bin,bin2);
                hEtaTwoClusterSize->SetTitle(name);
                histSaver->SaveHistogram(hEtaTwoClusterSize);
                delete hEtaTwoClusterSize;
            }
        }
        if (verbosity)  cout<<"Save "<< hEtaDistributionVsClusterSize[det]->GetName() << endl;
        histSaver->SaveHistogram(this->hEtaDistributionVsClusterSize[det]);

        histSaver->SaveHistogram(this->hEtaDistribution5Percent[det]);
        if (verbosity)  cout<<"Save "<< hEtaDistributionVsSignalLeft[det]->GetName() << endl;
        histSaver->SaveHistogram(this->hEtaDistributionVsSignalLeft[det]);
        histSaver->SaveHistogram(this->hEtaDistributionVsSignalRight[det]);
        if (verbosity)  cout<<"Save "<< hEtaDistributionVsSignalSum[det]->GetName() << endl;
        histSaver->SaveHistogram(this->hEtaDistributionVsSignalSum[det]);
        hSignalLeftVsSignalRight[det]->GetXaxis()->SetTitle("signalRight");
        hSignalLeftVsSignalRight[det]->GetYaxis()->SetTitle("signalLeft");
        if (verbosity)  cout<<"Save "<< hSignalLeftVsSignalRight[det]->GetName() << endl;
        histSaver->SaveHistogram(this->hSignalLeftVsSignalRight[det]);
        if(hEtaDistribution[det])delete hEtaDistribution[det];
        if(hDeltaLeftRightVsMaximum[det]){
            hDeltaLeftRightVsMaximum[det]->GetXaxis()->SetTitle("S_{Right} - S_{Left} / ADC");
            hDeltaLeftRightVsMaximum[det]->GetYaxis()->SetTitle("S_{Middle} / ADC");
            histSaver->SaveHistogram(this->hDeltaLeftRightVsMaximum[det]);
            TProfile* prof = hDeltaLeftRightVsMaximum[det]->ProfileY();
            if (TPlaneProperties::isDiamondDetector(det))
                prof->GetYaxis()->SetRangeUser(-500,500);
            else
                prof->GetYaxis()->SetRangeUser(-100,100);
            histSaver->SaveHistogram(prof,false,false,false);(hDeltaLeftRightVsMaximum[det],0,false);

            delete hDeltaLeftRightVsMaximum[det];
            hDeltaLeftRightVsMaximum[det]=0;
        }
        TProfile* profRight = 0;
        TProfile *profLeft = 0;
        if(hLeftVsMaximum[det]){
            hLeftVsMaximum[det]->GetXaxis()->SetTitle("S_{Left} / ADC");
            hLeftVsMaximum[det]->GetYaxis()->SetTitle("S_{Middle} / ADC");
            histSaver->SaveHistogram(this->hLeftVsMaximum[det]);
            profLeft = hLeftVsMaximum[det]->ProfileY();
            if (TPlaneProperties::isDiamondDetector(det))
                profLeft->GetYaxis()->SetRangeUser(-500,500);
            else
                profLeft->GetYaxis()->SetRangeUser(-100,100);
            histSaver->SaveHistogram(profLeft,false,false,false);
            delete hLeftVsMaximum[det];
            hLeftVsMaximum[det]=0;
        }
        if(hRightVsMaximum[det]){
            hRightVsMaximum[det]->GetXaxis()->SetTitle("S_{Right} / ADC");
            hRightVsMaximum[det]->GetYaxis()->SetTitle("S_{Middle} / ADC");
            histSaver->SaveHistogram(this->hRightVsMaximum[det]);
            profRight = hRightVsMaximum[det]->ProfileY();
            if (TPlaneProperties::isDiamondDetector(det))
                profRight->GetYaxis()->SetRangeUser(-500,500);
            else
                profRight->GetYaxis()->SetRangeUser(-100,100);
            histSaver->SaveHistogram(profRight,false,false,false);
            delete hRightVsMaximum[det];
            hRightVsMaximum[det]=0;
        }
        if (profLeft && profRight){
            histSaver->SaveTwoHistos("LeftRight_Vs_Maximum_profile",profLeft,profRight);
        }
        if (profLeft) delete profLeft;
        if (profRight) delete profRight;
    }
    if (verbosity)  cout<<"Save PH Histos"<<endl;
    savePHHistos();

    if (verbosity)  cout<<"Save Asymmetric Eta Sample analysis"<<endl;
//    analyseAsymmetricSample();
    savePedestalHistos();
    saveNoiseHistos();
    saveADCHistos();
    analysisSlopes();
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
    //	saveEtaInvestigationHistos();
}

void TAnalysisOfClustering::fillClusterHitPositions() {
    for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++){
    }

}

void TAnalysisOfClustering::saveEtaPerArea() {
    TH2F* histoNoCor = 0;
    int maxDia =(int)settings->getNDiamonds()+2;
    for ( Int_t dia=-1; hEtaDistributionDia->GetYaxis()->GetBinCenter(dia)<maxDia; dia++ ){
        TString histName = "hEtaDistribution_Dia_";
        if(dia==-1)
            histName.Append("All");
        else
            histName.Append(TString::Format("Area%d",(int)hEtaDistributionDia->GetYaxis()->GetBinCenter(dia)));

        TH1F* hist = 0;
        TH1F* histNoCor = 0;
        if(dia==-1){
            hist = (TH1F*)hEtaDistributionDia->ProjectionX(histName);
            if(histoNoCor) histNoCor =  (TH1F*)histoNoCor->ProjectionX();
        }
        else{
            hist = (TH1F*)hEtaDistributionDia->ProjectionX(histName,dia,dia);
            if (histoNoCor) histNoCor =  (TH1F*)histoNoCor->ProjectionX("_px",dia,dia);
        }
        TString title = "Eta Distribution Dia";
        if(dia==-1)
            title.Append(" All");
        else
            title.Append(TString::Format(" Area%d",(int)hEtaDistributionDia->GetYaxis()->GetBinCenter(dia)));
        int nbins = ETA_BINS;
        TString histName2 = histName;
        TString title2 = title;
        histName2.Append("inverted");
        title.Append(" inverted");
        TH1F* histInverted = new TH1F(histName2,title2,nbins+1,0-1.0/(2*nbins),1+1.0/(2*nbins));
        for(int bin = 1; bin <= nbins+1;bin++){
            histInverted->SetBinContent(bin,hist->GetBinContent(nbins+1+1-bin));
        }
        histInverted->SetLineColor(kBlue-7);
        if(hist){
            hist->SetTitle(title);
            string name = "c_";
            name.append(histName);
            if(hist->GetEntries()>0){
                TCanvas *c1 = new TCanvas(name.c_str());
                c1->cd();
                hist->Draw();
                if(histInverted) histInverted->Draw("same");
                hist->Draw("same");
                histSaver->SaveCanvas(c1);
                if(histNoCor&& histInverted) {
                    name.insert(1,"1");
                    name.append("_diff");
                    c1->Clear();
                    c1->SetName(name.c_str());
                    c1->SetTitle(name.c_str());
                    name = "stack_";
                    name.append(histName);
                    THStack *stack = new THStack(name.c_str(),hist->GetTitle());
                    stack->Add(histNoCor,"");//->Draw("AXIS");
                    stack->Add(histInverted,"");//->Draw("same");
                    stack->Add(hist,"");//->Draw("same");
                    stack->Draw("nostack");
                    histSaver->SaveCanvas(c1);
                }
            }
            else{
                TCanvas *c1 = new TCanvas(name.c_str());
                histSaver->SaveCanvas(c1);
            }
            histSaver->SaveHistogram(hist);
            if(hist) delete hist;

        }
    }
}

void TAnalysisOfClustering::saveProfileHistos(TProfile* pLeft, TProfile *pRight, Int_t etaLow, Int_t etaHigh,string name_comparision){
    cout << "Creating Profile Histos: " << etaLow << " - " << etaHigh << endl;
    Int_t nbins = pRight->GetNbinsX();
    Float_t xlow = pRight->GetBinLowEdge(1);
    Float_t xup = pRight->GetBinLowEdge(nbins+1);
    TString pname = TString::Format("pDiff_RightMinusLeft_%s_in_%03d_%03d",name_comparision.c_str(),etaLow,etaHigh);
    TProfile* pDiff = new TProfile(pname,pname,nbins,xlow,xup);
    for(int i = 0; i <= pLeft->GetNbinsX() && i <= pRight->GetNbinsX(); i++){
        Float_t value =  pRight->GetBinContent(i) - pLeft->GetBinContent(i);
        //		cout<< pDiff->GetBinCenter(i)<<": "<<value<<"\t\t"<<pRight->GetBinContent(i)<<" - "<<pLeft->GetBinContent(i)<<endl;
        pDiff->SetBinContent(i,value);
    }
    pDiff->SetLineColor(kBlue);
    pLeft->SetLineColor(kRed);
    pLeft->SetMarkerColor(kRed);
    pRight->SetLineColor(kGreen);
    pRight->SetMarkerColor(kGreen);
    Float_t minY = TMath::Max(-40.,pDiff->GetBinContent(pDiff->GetMinimumBin()));
    minY = minY>0?minY*.7:minY*1.1;
    TCanvas *c1 = new TCanvas(TString::Format("cProfile_%s_in_%03d_%03d",name_comparision.c_str(),etaLow,etaHigh));

    //		cout<<"compare: "<<pLeft->GetName()<<" vs. "<<pRight->GetName()<<":\t"<<c1->GetName()<<endl;
    c1->cd();
    TString title = TString::Format("hProfile_%s_eta_in_%03d_%03d",name_comparision.c_str(),etaLow,etaHigh);
    TH1F* pTitle = new TH1F(title,title,pLeft->GetNbinsX(),pLeft->GetXaxis()->GetXmin(),pLeft->GetXaxis()->GetXmax());
    Float_t maxX = 1000;
    Int_t maxBin = pLeft->FindBin(maxX);
    Double_t maxLeft =0;
    Double_t maxRight=0;
    for(Int_t i=0;i<maxBin;i++){
        maxLeft = TMath::Max(maxLeft,pLeft->GetBinContent(i));
        maxRight = TMath::Max(maxRight,pRight->GetBinContent(i));
    }
    pTitle->Draw();
    pTitle->GetXaxis()->SetRangeUser(0,maxX);
    //		pRight->GetXaxis()->SetRangeUser(0,maxX);
    //		pLeft->GetXaxis()->SetRangeUser(0,maxX);
    Float_t maxY = 1.1*TMath::Max(maxLeft,maxRight);
    pTitle->SetMaximum(maxY);
    pTitle->SetMinimum(minY);
    pTitle->GetYaxis()->SetRangeUser(minY,maxY);

    TString rightName = pRight->GetTitle();
    TString leftName = pLeft->GetTitle();
    if (rightName.Contains("rightright",TString::kIgnoreCase)){
        pTitle->GetXaxis()->SetTitle("signal 'Left / Right'");
        pTitle->GetYaxis()->SetTitle("avrg. signal adjacent ('LeftLeft / RightRight')");
    }
    else if (rightName.Contains("highest",TString::kIgnoreCase)){
        pTitle->GetXaxis()->SetTitle("highest Signal");
        pTitle->GetYaxis()->SetTitle("avrg. signal adjacent (Left/Right)");
    }
    else{
        pTitle->GetXaxis()->SetTitle("signal 1");
        pTitle->GetYaxis()->SetTitle("avrg. signal adjacent (signal2)");

    }
    pTitle->Draw();
    pRight->Draw("sameHIST");
    pLeft->Draw("sameHIST");
    //	pDiff->Draw("same");

    TLegend *leg = c1->BuildLegend();
    leg->Clear();
    leg->AddEntry(pRight);
    leg->AddEntry(pLeft);
    //	leg->AddEntry(pDiff);
    leg->SetFillColor(kWhite);
    Float_t x1NDC = leg->GetX1NDC();
    Float_t x2NDC = leg->GetX2NDC();
    Float_t val = x1NDC;
    x1NDC = -(x2NDC-x1NDC) + x1NDC;
    x2NDC=val;
    leg->SetX1NDC(.15);//x1NDC);
    leg->SetX2NDC(.5);//x2NDC);

    leg->Draw();
    histSaver->SaveCanvas(c1);
    histSaver->SaveHistogram(pLeft);
    histSaver->SaveHistogram(pRight);
    //	histSaver->SaveHistogram(pDiff);
    delete c1;
    if (pDiff) delete pDiff;
}

void TAnalysisOfClustering::saveEtaDividedHistos(TH3F* histo3dLeft,TH3F* histo3dRight, TH2F* histo2Left, TH2F* histo2Right,string name_comparision, Float_t etaWidth){
    if(etaWidth<=0) {
        cout<<"etaWidth <= 0"<<endl;
        return;
    }
    int minZbin,maxZbin, minXbin,maxXbin;
    int minZ,minX,maxX,maxZ;
    Float_t eta1,eta2;
    stringstream name;
    name<<setfill('0');
    vector<TProfile *> vecLeftProfiles;
    vector<TProfile *> vecRightProfiles;
    vector<Float_t> vecEtaLow;
    vector<Float_t> vecEtaHigh;
    cout<<"EtaWidth: "<<etaWidth<<endl;
    cout<<histo3dLeft->GetZaxis()->GetNbins()<<": "<<histo3dLeft->GetZaxis()->GetXmin()<<"-"<<histo3dLeft->GetZaxis()->GetXmax()<<endl;
    for(int i=0; i<=histo3dLeft->GetZaxis()->GetNbins()+1;i++)
        cout<<histo3dLeft->GetZaxis()->GetBinLowEdge(i)<<", "<<flush;
    cout<<endl;
    for(Float_t eta=0;eta<=1;eta+=etaWidth){
        eta1 = eta<1?eta:-1000;
        eta2 = eta<1?eta+etaWidth:1000;
        cout<<"eta in "<<eta1<<"-"<<eta2<<endl;
        /**** 2D  ***/

        name.str("");
        name.clear();
        minZbin = histo3dLeft->GetZaxis()->FindBin(eta1);
        maxZbin = histo3dLeft->GetZaxis()->FindBin(eta2);
        //		cout<<"eta "<<eta<<": "<<minZbin<<"-"<<maxZbin<<endl;
        histo3dLeft->GetZaxis()->SetRange(minZbin, maxZbin);
        minZ = histo3dLeft->GetZaxis()->GetBinLowEdge(minZbin)*100;
        minZ=minZ<0?0:minZ;
        maxZ = histo3dLeft->GetZaxis()->GetBinLowEdge(maxZbin)*100;
        maxZ=maxZ>100?100:maxZ;
        name<<histo3dLeft->GetName()<<"_eta_"<< setw(3)<<(minZ)<<"_"<<setw(3)<<(maxZ);
        histo3dLeft->GetZaxis()->SetBit(TAxis::kAxisRange);
        TH2F* histo2dLeft = (TH2F*)histo3dLeft->Project3D("yx");
        histo2dLeft->SetName(name.str().c_str());
        histo2dLeft->SetTitle(name.str().c_str());
        name.str("");
        name.clear();
        name<<histo2dLeft->GetName()<<"_pfx";
        TProfile * hLeft_pfx = histo2dLeft->ProfileX(name.str().c_str());
        cout << name.str() << endl;
        hLeft_pfx->GetXaxis()->SetTitle("signal 'Left'");
        hLeft_pfx->GetYaxis()->SetTitle("avrg. signal 'LeftLeft'");
        hLeft_pfx->Draw("goff");
        if(hLeft_pfx->GetXaxis()->GetXmax()>2000)
            hLeft_pfx->GetXaxis()->SetRangeUser(0,2000);
        vecLeftProfiles.push_back(hLeft_pfx);
        vecEtaLow.push_back(minZ);
        vecEtaHigh.push_back(maxZ);

        name.str("");
        name.clear();
        minZbin = histo3dRight->GetZaxis()->FindBin(eta1);
        maxZbin = histo3dRight->GetZaxis()->FindBin(eta2);
        //		cout<<eta<<": "<<minZbin<<"-"<<maxZbin<<endl;
        histo3dRight->GetZaxis()->SetRange(minZbin, maxZbin);
        minZ = histo3dRight->GetZaxis()->GetBinLowEdge(minZbin)*100;
        minZ=minZ<0?0:minZ;
        maxZ = histo3dRight->GetZaxis()->GetBinLowEdge(maxZbin)*100;
        name<<histo3dRight->GetName()<<"_eta_"<<setw(3)<<(minZ)<<"_"<<setw(3)<<(maxZ);
        histo3dRight->GetZaxis()->SetBit(TAxis::kAxisRange);

        TH2F* histo2dRight = (TH2F*)histo3dRight->Project3D("yx");
        histo2dRight->SetName(name.str().c_str());
        histo2dRight->SetTitle(name.str().c_str());
        name.str("");
        name.clear();
        name<<histo2dRight->GetName()<<"_pfx";
        TProfile * hRight_pfx = histo2dRight->ProfileX(name.str().c_str());
        cout << name.str() << endl;
        hRight_pfx->GetXaxis()->SetTitle("signal 'Right'");
        hRight_pfx->GetYaxis()->SetTitle("avrg. signal 'RightRight'");
        hRight_pfx->Draw("goff");
        if(hRight_pfx->GetXaxis()->GetXmax()>2000)
            hRight_pfx->GetXaxis()->SetRangeUser(0,2000);
        vecRightProfiles.push_back(hRight_pfx);

        /**** 1D  ***/
        name.str("");
        name.clear();
        minXbin = histo2Left->GetXaxis()->FindBin(eta1);
        maxXbin = histo2Left->GetXaxis()->FindBin(eta2);
        minX = histo2Left->GetXaxis()->GetBinLowEdge(minXbin)*100;
        minX=minX<0?0:minX;
        maxX = histo2Left->GetXaxis()->GetBinLowEdge(maxXbin+1)*100;
        name<<histo2Left->GetName()<<"_eta_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
        TH1F *histo1dLeft = (TH1F*)histo2Left->ProjectionY(name.str().c_str(),minXbin,maxXbin);
        histo1dLeft->SetTitle(name.str().c_str());

        name.str("");
        name.clear();
        minXbin = histo2Right->GetXaxis()->FindBin(eta1);
        maxXbin = histo2Right->GetXaxis()->FindBin(eta2);
        minX = histo2Right->GetXaxis()->GetBinLowEdge(minXbin)*100;
        minX = minX<0?0:minX;
        maxX = histo2Right->GetXaxis()->GetBinLowEdge(maxXbin+1)*100;
        name<<histo2Right->GetName()<<"_eta_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
        TH1F *histo1dRight = (TH1F*)histo2Right->ProjectionY(name.str().c_str(),minXbin,maxXbin);
        histo1dRight->SetTitle(name.str().c_str());

        //		cout<<"save "<<histo1dRight->GetName()<<endl;
        histSaver->SaveHistogram(histo1dRight);
        //		cout<<"save "<<histo1dLeft->GetName()<<endl;
        histSaver->SaveHistogram(histo1dLeft);
        //		cout<<"saving "<<histo2dLeft->GetName()<<endl;;
        histSaver->SaveHistogram(histo2dLeft,false);
        //		cout<<"saving "<<histo2dRight->GetName()<<endl;
        histSaver->SaveHistogram(histo2dRight,false);
        name.str("");
        name.clear();
        name<<"cComparision_"<<name_comparision<<"_"<<setw(3)<<minX<<"_"<<setw(3)<<maxX;
        histo1dRight->SetLineColor(kBlue);
        histo1dLeft->SetLineColor(kRed);
        histSaver->SaveTwoHistos(name.str(),histo1dRight,histo1dLeft,1);

        if (histo1dRight) delete histo1dRight;
        if(histo2dRight) delete histo2dRight;
        if (histo1dLeft) delete histo1dLeft;
        if(histo2dLeft) delete histo2dLeft;
    }

    TProfile *pLowLeft = 0;
    TProfile *pLowRight = 0;
    for(UInt_t i=0; i<vecLeftProfiles.size()&&i<vecRightProfiles.size();i++){
        TProfile* pLeft = vecLeftProfiles.at(i);
        Int_t j = vecRightProfiles.size()-1-i;
        TProfile* pRight = vecRightProfiles.at(j);
        Float_t etaLow = vecEtaLow.at(i);
        Float_t etaHigh = vecEtaHigh.at(i);
        cout<<"create profile:"<<"\n\t"<<pLeft->GetName()<<"\n\t"<<pRight->GetName()<<endl;
        if(pLowLeft == 0 && etaHigh <= etaWidth * 100){
            TString name = TString::Format("pLeft_%s_eta_%3.f_%3.f",name_comparision.c_str(),100-etaWidth*100,etaWidth*100);
            cout<<"Creating "<<name<<endl;
            pLowLeft = (TProfile*)pLeft->Clone(name);
        }
        if(pLowRight == 0 && etaHigh <= etaWidth * 100){
            TString name = TString::Format("pRight_%s_eta_%3.f_%3.f",name_comparision.c_str(),100-etaWidth*100,etaWidth*100);
            cout<<"Creating "<<name<<endl;
            pLowRight = (TProfile*)pRight->Clone(name);
        }
        if(pLowLeft && etaLow >= (1-etaWidth)*100){
            cout << "Adding " << pLeft->GetName() << " to " << pLowLeft->GetName() << endl;
            pLowLeft->Add(pLeft);
        }
        if(pLowRight&& etaLow >= (1-etaWidth)*100){
            cout << "Adding " << pRight->GetName() << " to " << pLowRight->GetName() << endl;
            pLowRight->Add(pRight);
        }
        saveProfileHistos(pLeft,pRight,etaLow,etaHigh,name_comparision);
    }
    if(pLowLeft&&pLowRight)
        saveProfileHistos(pLowLeft,pLowRight,(1-etaWidth)*100,etaWidth*100,name_comparision);

    for(UInt_t i=0; i<vecLeftProfiles.size()&&i<vecRightProfiles.size();i++){
        delete vecLeftProfiles.at(i);
        delete vecRightProfiles.at(i);
    }
    cout<<"saving "<<histo3dRight->GetName()<<endl;
    histSaver->SaveHistogramROOT(histo3dRight);
    cout<<"saving "<<histo3dLeft->GetName()<<endl;
    histSaver->SaveHistogramROOT(histo3dLeft);

    //	cout<<"saving "<<histo2Left->GetName()<<endl;
    histSaver->SaveHistogram(histo2Left);
    //	cout<<"saving "<<histo2Right->GetName()<<endl;
    histSaver->SaveHistogram(histo2Right);
}

void TAnalysisOfClustering::saveEtaInvestigationHistos(){
    stringstream name;
    for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
        name.str("");
        name.clear();
        name<<"hLeftLeft_"<<TPlaneProperties::getStringForDetector(det);
        TH2F* histo2Left = histSaver->CreateScatterHisto(name.str(),vecvecSignalLeftLeft.at(det),vecvecEta.at(det),512,0,1);
        histo2Left->GetXaxis()->SetTitle("#eta");
        histo2Left->GetYaxis()->SetTitle("Signal 'LeftLeft'");

        name.str("");
        name.clear();
        name<<"hRightRight_"<<TPlaneProperties::getStringForDetector(det);
        TH2F* histo2Right = histSaver->CreateScatterHisto(name.str(),vecvecSignalRightRight.at(det),vecvecEta.at(det),512,0,1);
        histo2Right->GetXaxis()->SetTitle("#eta");
        histo2Right->GetYaxis()->SetTitle("Signal 'RightRight'");

        name.str("");name.clear();
        name<<"hLeftLeft_vs_Left_3D_"<<TPlaneProperties::getStringForDetector(det);
        Float_t inf = std::numeric_limits<float>::infinity();
        TH3F* histo3dLeft = histSaver->Create3DHisto(name.str(),vecvecLeftEtaSignal.at(det),vecvecSignalLeftLeft.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
        histo3dLeft->GetXaxis()->SetTitle("Signal 'Left'");
        histo3dLeft->GetYaxis()->SetTitle("Signal 'LeftLeft'");
        histo3dLeft->GetZaxis()->SetTitle("#eta");
        /*****************/


        name.str("");
        name.clear();
        name<<"hRightRight_vs_Right_3D_"<<TPlaneProperties::getStringForDetector(det);
        TH3F* histo3dRight = histSaver->Create3DHisto(name.str(),vecvecRightEtaSignal.at(det),vecvecSignalRightRight.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
        histo3dRight->GetXaxis()->SetTitle("Signal 'Right'");
        histo3dRight->GetYaxis()->SetTitle("Signal 'RightRight'");
        histo3dRight->GetZaxis()->SetTitle("#eta");

        name.str("");
        name.str();
        name<<"LeftLeft-RightRight_"<<TPlaneProperties::getStringForDetector(det);
        this->saveEtaDividedHistos(histo3dLeft,histo3dRight,histo2Left,histo2Right,name.str());


        if (histo3dRight) delete histo3dRight;
        if (histo3dLeft) delete histo3dLeft;
        if (histo2Left) delete histo2Left;
        if (histo2Right) delete histo2Right;

        /**************************/

        name.str("");
        name.clear();
        name<<"hLeft_vs_Highest_3D_"<<TPlaneProperties::getStringForDetector(det);
        histo3dLeft = histSaver->Create3DHisto(name.str(),vecvecSignalHighest.at(det),vecvecSignalLeftOfHighest.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
        name.str("");
        name.clear();
        name<<"hRight_vs_Highest_3D_"<<TPlaneProperties::getStringForDetector(det);
        histo3dRight = histSaver->Create3DHisto(name.str(),vecvecSignalHighest.at(det),vecvecSignalRightOfHighest.at(det),vecvecEta.at(det),512,512,20,-100,inf,-100,inf,0,1);
        name.str("");
        name.clear();

        name<<"hLeftOfHighest_"<<TPlaneProperties::getStringForDetector(det);
        histo2Left = histSaver->CreateScatterHisto(name.str(),vecvecSignalLeftOfHighest.at(det),vecvecEta.at(det),512,0,1);
        histo2Left->GetXaxis()->SetTitle("#eta");
        histo2Left->GetYaxis()->SetTitle("Signal left of highest channel");

        name.str("");
        name.clear();
        name<<"hRightOfHighest_"<<TPlaneProperties::getStringForDetector(det);
        histo2Right = histSaver->CreateScatterHisto(name.str(),vecvecSignalRightOfHighest.at(det),vecvecEta.at(det),512,0,1);
        histo2Right->GetXaxis()->SetTitle("#eta");
        histo2Right->GetYaxis()->SetTitle("Signal right of highest channel");

        name.str("");
        name.clear();
        name<<"LeftHighest-RightHighest_"<<TPlaneProperties::getStringForDetector(det);
        this->saveEtaDividedHistos(histo3dLeft,histo3dRight,histo2Left,histo2Right,name.str());

        if (histo3dRight) delete histo3dRight;
        if (histo3dLeft) delete histo3dLeft;
        if (histo2Left) delete histo2Left;
        if (histo2Right) delete histo2Right;



    }
}
void TAnalysisOfClustering::compareCentroid_ChargeWeightedMean()
{
    bool check=true;
    for(int det=0;det<9;det++)
        check=eventReader->getNClusters(det)==1;
    if(check==true){
        TCluster cluster = eventReader->getCluster(8,0);
        Float_t xCWM=cluster.getChargeWeightedMean(settings->doCommonModeNoiseCorrection());
        Float_t xHit=(Float_t)cluster.getHighestSignalChannel();
        Float_t xH2C=(Float_t)cluster.getHighest2Centroid();
        Float_t delta=xCWM-xHit;
        this->histo_CWM_biggestHit->Fill(delta,cluster.size());
        delta = xH2C - xHit;
        this->histo_H2C_biggestHit->Fill(delta);
        //		if(eventReader->getNClusters(8)>=1){
        //			Float_t charge = eventReader->getCluster(8,0).getCharge();
        //			Float_t signal2ndHighestHit=eventReader->getCluster(8,0).getCharge(2)-eventReader->getCluster(8,0).getCharge(1);
        //			Float_t q =signal2ndHighestHit/charge;
        //		}


    }
}

void TAnalysisOfClustering::analyseClusterPosition()
{
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
            bool bCMN = TPlaneProperties::isDiamondDetector(det);
            Float_t posCWM = eventReader->getEvent()->getPosition(det,cl,bCMN,TCluster::chargeWeighted);
            Int_t chNo = (Int_t)(posCWM+0.5);
            Float_t relPos = posCWM - chNo;
            hClusterPosition[det]->Fill(posCWM);
            if(TPlaneProperties::isDiamondDetector(det))
                hClusterPositionActiveChannels->Fill(posCWM);
            hRelativeClusterPositionCWM[det]->Fill(chNo+0.5,relPos);
            //			Float_t eta =
            TCluster cluster = eventReader->getCluster(det,cl);

            UInt_t highestClPos=cluster.getHighestHitClusterPosition();
            UInt_t nextHighestClPos=cluster.getHighestSignalNeighbourClusterPosition(highestClPos);
            if(nextHighestClPos==9999){
                //cout<<"\nnext highest=9999: "<<highestClPos<<" "<<eventReader->getCluster(det,cl).getClusterSize()<<endl;;
                //eventReader->getCluster(det,cl).Print(1);
                continue;
            }

            Float_t signalLeft,signalRight,adcRight,adcLeft,pedRight,pedLeft;
            UInt_t leftClPos,rightClPos;
            if(highestClPos>nextHighestClPos) {
                leftClPos=nextHighestClPos;
                rightClPos=highestClPos;
            }
            else{
                rightClPos=nextHighestClPos;
                leftClPos=highestClPos;
            }
            int clusterSize =cluster.getClusterSize();
            if(settings->isMaskedCluster(det,cluster,true))
                return;
            signalLeft= cluster.getSignal(leftClPos);
            signalRight=cluster.getSignal(rightClPos);
            adcLeft= cluster.getAdcValue(leftClPos);
            adcRight= cluster.getAdcValue(rightClPos);
            pedLeft= cluster.getPedestalMean(leftClPos);
            pedRight= cluster.getPedestalMean(rightClPos);
            Float_t a= 0.03;
            Float_t adcLeftReal=adcLeft/(1-a);
            Float_t adcRightReal=(-a*(1-a)*adcLeft+(1-a)*adcRight)/(1-a)/(1-a);
            Float_t signalAdcLeft = adcLeftReal-pedLeft;
            Float_t signalAdcRight= adcRightReal-pedRight;
            Float_t eta2=(signalAdcRight)/(signalAdcLeft+signalAdcRight);
            Int_t leftChannel=-1;
            Float_t eta= cluster.getEta(leftChannel);
            Float_t etaCmnCorrected = cluster.getEta(true);
            Float_t charge = cluster.getCharge((UInt_t)2,true,true);
            Float_t S_M = cluster.getSignal(highestClPos);
            Float_t S_R = cluster.getSignal(highestClPos+1);
            Float_t S_L = cluster.getSignal(highestClPos-1);
            //					cluster.Print();
            //					cout<<"Cluster: "<<S_L<<"-"<<S_M<<"-"<<S_R<<endl;

            hDeltaLeftRightVsMaximum[det]->Fill(S_R-S_L,S_M);
            hLeftVsMaximum[det]->Fill(S_L,S_M);
            hRightVsMaximum[det]->Fill(S_R,S_M);
            if (cluster.isSaturatedCluster())
                return;
            if(verbosity>3){
                cout<<"charge of 2: "<<charge<<"\t"<<flush;
                cluster.Print();
            }
            //			Float_t eta3= signalRightReal/(signalLeftReal+signalRightReal);
            //			cout<<nEvent<<" "<<eta<<" "<<eta1<<" "<<eta2<<" "<<eta3<<endl;
            if(hEtaDistribution[det])
                hEtaDistribution[det]->Fill(eta);
            if (TPlaneProperties::isDiamondDetector(det)){
                Int_t nDia = settings->getDiaDetectorAreaOfChannel(leftChannel)+1;
                if(eta>0&&eta<1){
                    hEtaDistributionDia->Fill(eta,nDia);
                }
            }
            if(hEtaDistributionCMN[det])hEtaDistributionCMN[det]->Fill(etaCmnCorrected);
            if(hEtaDistributionVsLeftChannel[det]) hEtaDistributionVsLeftChannel[det]->Fill(eta,leftChannel);
            if(hEtaDistributionVsClusterSize[det]) hEtaDistributionVsClusterSize[det]->Fill(eta,clusterSize);
            if(hEtaDistributionVsCharge[det]) hEtaDistributionVsCharge[det]->Fill(eta,charge);
            hEtaDistribution5Percent[det]->Fill(eta2);
            hSignalLeftVsSignalRight[det]->Fill(signalRight,signalLeft);

            //			cout<<nextHighestClPos<<"<"<<highestClPos<<"\t"<<signalLeft<<" < "<<signalRight<<"\t"<<eta<<endl;
            hEtaDistributionVsSignalLeft[det]->Fill(eta,signalLeft);
            hEtaDistributionVsSignalRight[det]->Fill(eta,signalRight);
            hEtaDistributionVsSignalSum[det]->Fill(eta,signalLeft+signalRight);
            TH1F *hEtaIntegral=eventReader->getEtaIntegral(det);
            Float_t posCorEta=  eventReader->getEvent()->getPosition(det,cl,settings->doCommonModeNoiseCorrection(),TCluster::corEta,hEtaIntegral);
            chNo = (UInt_t)(posCorEta+0.5);
            relPos = posCorEta - chNo;
            if(verbosity>6) printf("%5d %3d %5.1f %5.1f\n",nEvent,chNo,posCWM,posCorEta);
            hRelativeClusterPositionCorEta[det]->Fill(chNo+0.5,relPos);
            Float_t posEta=  eventReader->getEvent()->getPosition(det,cl,TCluster::eta);
            chNo = (UInt_t)(posEta+0.5);
            relPos = posEta - chNo;
            hRelativeClusterPositionEta[det]->Fill(chNo+0.5,relPos);
        }
    }
}

void TAnalysisOfClustering::analyseCluster(){
    for(int det=0;det<9;det++){
        hNumberOfClusters[det]->Fill(eventReader->getNClusters(det));
        for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
            UInt_t clSize = eventReader->getClusterSize(det,cl);
            hClusterSize[det]->Fill(clSize);
            hClusterSeedSize[det]->Fill(eventReader->getClusterSeedSize(det,cl));
            Float_t biggestSignal = eventReader->getCluster(det,cl).getHighestSignal();
            clSize=clSize>8?8:clSize;
            hBiggestHitVsClusterSize[det]->Fill(biggestSignal,clSize);
            Float_t pos = eventReader->getCluster(det,cl).getPosition(false,TCluster::highest2Centroid);
            Float_t posCMN = eventReader->getCluster(det,cl).getPosition(true,TCluster::highest2Centroid);
            hClusterPos[det]->Fill(pos);
            hClusterPosCMN[det]->Fill(posCMN);
        }
    }

}


void TAnalysisOfClustering::analyse2ndHighestHit(){
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        //		if(nEvent==20){
        //			cout<<nEvent<<":"<<endl;
        //			eventReader->setVerbosity(10);
        //		}
        //		else eventReader->setVerbosity(0);
        Float_t nClusters=eventReader->getNClusters(det);

        for(UInt_t cl=0;cl<nClusters;cl++){
            TCluster cluster=eventReader->getCluster(det,cl);
            if(settings->isMaskedCluster(det,cluster,true))
                continue;
            if(cluster.size()==0){
                cout<<nEvent<<" "<<cl<<" "<<eventReader->getNClusters(det)<<" "<<cluster.size()<<endl;
                cluster.Print();
                continue;//why can this happen???
            }
            UInt_t highestChannel = cluster.getHighestSignalChannel();
            Float_t signalLeft = cluster.getSignalOfChannel(highestChannel-1);
            Float_t signalRight = cluster.getSignalOfChannel(highestChannel+1);
            if(signalLeft<0 && verbosity>2)
                cout<<"signalLeft is smaller than 0"<<endl;
            if(signalRight<0 && verbosity>2)
                cout<<"signalLeft is smaller than 0"<<endl;
            Float_t signalHighest = cluster.getHighestSignal();
            Float_t signal2ndHighest;
            Float_t deltaSignals = signalLeft-signalRight;
            if(deltaSignals<0)// right channel higher left channel
                signal2ndHighest=signalRight;
            else
                signal2ndHighest=signalLeft;
            Float_t sumSignals = signalLeft+signalRight;
            if(signalLeft==0&&signalRight==0)continue;
            if(sumSignals==0)continue;
            Float_t allCharge=cluster.getCharge(settings->doCommonModeNoiseCorrection());
            //			Float_t charge = cluster.getCharge(false);
            //			if (cluster.isHit())
            Float_t signalRatio=signal2ndHighest/signalHighest;
            if(signalRatio>1){
                cout<<"Ratio>1:"<<signal2ndHighest<<" "<<signalHighest<<endl;
                cluster.Print();
            }
            else
                hSignal2ndHighestOverSignalHighestRatio[det]->Fill(signalRatio);
            Float_t ratio;
            if(cluster.size()!=1){
                if(signalLeft>signalRight){
                    ratio=signalLeft/allCharge;
                    if(ratio>0.5||allCharge==0||ratio!=ratio){
//                        cout<<"\n2ndBiggestHitOverCharge>0.5: left "<<signalLeft<<" "<<allCharge<<endl;
//                        cluster.Print();
                    }
                    else{
                        //					cout<<nEvent<<" "<<cl<<" "<<ratio<<endl;
                        h2ndBiggestHitOverCharge[det]->Fill(ratio);
                    }
                    h2ndBiggestHitSignal[det]->Fill(signalLeft);
                }
                else{
                    ratio=signalRight/allCharge;
                    if(ratio>0.5||allCharge==0||ratio!=ratio){
//                        cout<<"\n2ndBiggestHitOverCharge>0.5: right"<<signalRight<<" "<<allCharge<<endl;
//                        cluster.Print();
                    }
                    else{
                        //					cout<<nEvent<<" "<<cl<<" "<<ratio<<endl;
                        h2ndBiggestHitOverCharge[det]->Fill(ratio);
                    }

                    h2ndBiggestHitSignal[det]->Fill(signalRight);
                }
            }
            if(signalLeft>signalRight){
                h2ndBiggestHitPosition[det]->Fill(-1);
            }
            else if(signalLeft<signalRight)
                h2ndBiggestHitPosition[det]->Fill(+1);
            else
                h2ndBiggestHitPosition[det]->Fill(0);
            ratio = (deltaSignals)/(sumSignals);
            if (ratio<1&&ratio>-1&&sumSignals!=0)
                hDeltaLeftRightHitOverLeftAndRight[det]->Fill(ratio);
            else {
                if(TMath::Abs(ratio)>1){
                    //					cout<<"hDeltaLeftRightHitOverLeftAndRight "<<det<<" "<<cl<<" "<<deltaSignals<<" "<<sumSignals<<endl;
                    //					cluster.Print();
                }
            }
            ratio = signalLeft/sumSignals;
            if(signalLeft>0&&sumSignals>0&&ratio<1&&ratio>0)hLeftHitOverLeftAndRight[det]->Fill((signalLeft)/(sumSignals));
        }
    }
}




void TAnalysisOfClustering::analyseAsymmetricSample(){
    cout<<"\nAnalyseAsymmetricSample:"<<vecVecClusters.size()<<endl;
    vector <Float_t> vecAlphas;
    vector <Float_t> vecNSteps;
    for(UInt_t det = 0 ; det<TPlaneProperties::getNDetectors();det ++){
        cout<<"***\n***\n Analyse: " << det << ": " << vecVecClusters.at(det).size() << endl;
        TAnalysisOfAsymmetricEta *asymmetricEta  = new TAnalysisOfAsymmetricEta(this->settings);
        if(asymmetricEta){
            asymmetricEta->setClusters(vecVecClusters.at(det));
            asymmetricEta->setDetector(det);
            vecNSteps.push_back(asymmetricEta->analyse());
            cout<<"Final Alpha: "<<asymmetricEta->getAlpha()*100 << " %"<<endl;
            vecAlphas.push_back(asymmetricEta->getAlpha());
            delete asymmetricEta;
        }
    }
    cout<<"****\n****\n FINAL RESULTS: "<<endl;
    std::ofstream ofs (settings->getCrossTalkFactorsFileName().c_str(), std::ofstream::out);
    for(UInt_t det = 0 ; det < vecAlphas.size();det ++){
        if (res){
            res->setSignalFeedOverCorrection(det, vecAlphas[det]);
            res->setFloatValue("Feed_Through_Correction",TPlaneProperties::getStringForDetector(det),vecAlphas[det]);
        }
        cout<<det<<": "<< TString::Format("%02.2f",vecAlphas.at(det)*100) << "%\t in "<<setw(2)<< vecNSteps.at(det) << " Steps" ;
        ofs<<det<<": "<<  TString::Format("%02.2f",vecAlphas.at(det)*100) << "%\t in "<<setw(2)<< vecNSteps.at(det) << " Steps";

        if(det<vecPHMeans.size()){
            cout<<TString::Format("\tMean: %4.1f, Noise: %2.2f",vecPHMeans.at(det),vecWidth.at(det))<<endl;
            ofs<<TString::Format("\tMean: %4.1f, Noise: %2.2f",vecPHMeans.at(det),vecWidth.at(det))<<endl;
        }
        else{
            cout<<endl;
            ofs<<endl;
        }
    }
    ofs.close();
}




void TAnalysisOfClustering::fillClusterVector(){
    //	if (!settings->isAsymmetricSample())
    //		return;
    for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
        //	UInt_t det = 8;
        if(eventReader->getNClusters(det)!=1)
            return;
        if(TPlaneProperties::isSiliconDetector(det))
            if (vecVecClusters.at(det).size()>nMaxClusters)
                continue;
        if(TPlaneProperties::isDiamondDetector(det))
            if (vecVecClusters.at(det).size()>3*nMaxClusters)
                continue;
        TCluster cluster = eventReader->getCluster(det,0);
        vecVecClusters.at(det).push_back(cluster);
    }
}


void TAnalysisOfClustering::etaInvestigation(){
    for(int det=0;det<9;det++){
        Float_t minSignal;
        if (TPlaneProperties::isSiliconDetector(det))
            minSignal = -10;
        else
            minSignal = -50;//todo
        for(UInt_t cl = 0; cl< eventReader->getNClusters(det);cl++){
            TCluster cluster = eventReader->getCluster(det,cl);
            //			cluster.Print(1);
            Int_t leftEtaChannel =-1;
            Float_t eta = cluster.getEta(leftEtaChannel,false);
            Int_t leftEtaClusterPosition = cluster.getClusterPosition(leftEtaChannel);
            Float_t signalLeftLeft = cluster.getSignal(leftEtaClusterPosition-1);
            Float_t leftEtaSignal = cluster.getSignal(leftEtaClusterPosition);
            Float_t rightEtaSignal = cluster.getSignal(leftEtaClusterPosition+1);
            Float_t signalRightRight = cluster.getSignal(leftEtaClusterPosition+2);
            Float_t highestSignal = cluster.getHighestSignal();
            Int_t highestSignalChannel = cluster.getHighestSignalChannel();
            Int_t highestSignalClusterPos = cluster.getClusterPosition(highestSignalChannel);
            Float_t signalLeftOfHighest = cluster.getSignal(highestSignalClusterPos-1);
            Float_t signalRightOfHighest = cluster.getSignal(highestSignalClusterPos+1);
            if(signalLeftLeft<minSignal||signalLeftOfHighest<minSignal||signalRightOfHighest<minSignal||signalRightRight<minSignal)
                continue;
            //			cout<< signalLeftOfEta<<" "<< leftEtaSignal<<" "<<rightEtaSignal<<" "<<signalRightOfEta<<"  ----  "<<signalLeftOfHighest<<" "<< highestSignal<<" "<<signalRightOfHighest<<"\n"<<endl;
            //		Float_t signalLeftOfEta = eventReader->
            if (signalLeftLeft < -100){
                int clusPos = leftEtaClusterPosition-1;
                cout<<"\n"<< nEvent<<", Problem signalLeftLeft: "<<signalLeftLeft<<" "<<leftEtaChannel<<": "<<cluster.getPedestalMean(clusPos)<<", "<<cluster.getAdcValue(clusPos)<<endl;
                cluster.Print(1);
            }
            vecvecSignalLeftLeft.at(det).push_back(signalLeftLeft);
            vecvecSignalRightRight.at(det).push_back(signalRightRight);;
            vecvecLeftEtaSignal.at(det).push_back(leftEtaSignal);;
            vecvecRightEtaSignal.at(det).push_back(rightEtaSignal);;
            vecvecSignalLeftOfHighest.at(det).push_back(signalLeftOfHighest);
            vecvecSignalRightOfHighest.at(det).push_back(signalRightOfHighest);
            vecvecSignalHighest.at(det).push_back(highestSignal);
            vecvecEta.at(det).push_back(eta);
        }
    }
}


void TAnalysisOfClustering::savePHHistos()
{
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        histSaver->SaveHistogram(hPHDistribution[det]);
        histSaver->SaveHistogram(hPHDistIncreasingClustersize[det]);
        vecClusterSize.clear();
        vecMPV.clear();
        vecClusterSizeError.clear();
        vecWidth.clear();
        for(UInt_t nClusters=0;nClusters<10;nClusters++){
            stringstream histName;
            //    		string name = (string)hPHDistribution[det]->GetTitle();
            //    		name = name.substr(0,histName.str().size()-6);
            TString hname = "hPulseHeightDistribution_";
            histName<<"hPulseHeightDistribution_";
            if(nClusters==0)
                hname += (TString) "allClusterSizes_";
            else
                hname += TString::Format("nClusters%d_",nClusters);

            hname+=(TString)TPlaneProperties::getStringForDetector(det);
            TObject *htemp2 = (TObject*)gROOT->FindObject(hname);
            if(htemp2!=0)delete htemp2;

            //CREATE HTEMP and ReBin it if necessary
            TH1F *htemp;
            //			if (nClusters == 0)
            //				htemp = (TH1F*)hPHDistribution[det]->ProjectionX(hname);//,nClusters+1,nClusters+1);
            //			else
            htemp = (TH1F*)hPHDistribution[det]->ProjectionX(hname,nClusters+1,nClusters+1);
            if(htemp==0) continue;

            //adjust binning if necessary
            UInt_t entries = htemp->GetEntries();
            UInt_t maximumEntries = htemp->GetMaximum();
            UInt_t nSteps =4;
            UInt_t nStep = 0;
            while((maximumEntries<50&&maximumEntries<entries*0.7)&&nStep<nSteps){
                htemp->Rebin(2);
                entries = htemp->GetEntries();
                maximumEntries = htemp->GetMaximum();
                nStep++;
            }
            TF1 *fit=0;
            htemp->SetTitle(htemp->GetName());
            htemp->GetXaxis()->SetTitle("Charge in ADC units");
            htemp->GetYaxis()->SetTitle("number of entries#");
            LandauGaussFit landauGauss;
            if(nClusters<4||det==TPlaneProperties::getDetDiamond())
                fit = landauGauss.doLandauGaussFit(htemp,nClusters==1&&det==TPlaneProperties::getDetDiamond());
            if(fit!=0){
                if (verbosity > 2) {
                    cout<<"Width(scale): "<<fit->GetParameter(0)<<endl;
                    cout<<"MostProb:     "<<fit->GetParameter(1)<<endl;
                    cout<<"Area:         "<<fit->GetParameter(2)<<endl;
                    cout<<"Width(sigma): "<<fit->GetParameter(3)<<endl;
                }
                if(nClusters==0)
                    vecPHMeans.push_back(fit->GetParameter(1));
                vecClusterSize.push_back(nClusters);
                vecMPV.push_back(fit->GetParameter(1));
                vecClusterSizeError.push_back(0.5);
                vecWidth.push_back(fit->GetParameter(0));
                histSaver->SaveHistogramLandau(htemp);
            }
            else
                histSaver->SaveHistogram(htemp);
            delete htemp;
        }
        if(det==TPlaneProperties::getDetDiamond()){
            stringstream histTitle;
            histTitle<<"gChargeOfClusterVsClusterSize_"<<det;
            TGraphErrors graph = histSaver->CreateErrorGraph(histTitle.str(),vecClusterSize,vecMPV,vecClusterSizeError,vecWidth);
            graph.GetXaxis()->SetTitle("Cluster Size");
            graph.GetYaxis()->SetTitle("Charge of Cluster");
            histSaver->SaveGraph(&graph,histTitle.str());
        }
        delete hPHDistribution[det];
        vector<TH1D*> vhtemps;
        TString hname = "hStackPHDist_ClusterSize";
        hname += TPlaneProperties::getStringForDetector(det);
        THStack* stack = new THStack(hname,hname);
        for (UInt_t cl = 1; cl < 5; cl++){
             hname = TString::Format("hPHDist_Cl%d_",cl);
            hname += TPlaneProperties::getStringForDetector(det);
            TH1D* htemp = hPHDistIncreasingClustersize[det]->ProjectionX(hname,cl,cl);
            vhtemps.push_back(htemp);
            histSaver->SaveHistogram(htemp);
            stack->Add(htemp);
        }
        histSaver->SaveStack(stack,"nostack",true,true,"Cluster Signal / ADC","number of entries");
        delete stack;
        for (UInt_t i = 0; i < vhtemps.size(); i++)
            delete vhtemps[i];
        delete hPHDistIncreasingClustersize[det];
    }
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        histSaver->SaveHistogram(hBiggestHitVsClusterSize[det]);
        TProfile *profY = hBiggestHitVsClusterSize[det]->ProfileY();
        profY->GetXaxis()->SetTitle("ClusterSize");
        profY->GetYaxis()->SetTitle("mean of Biggest signal in Cl");
        string histoTitle = "mean of Biggest signal in Cluster vs. ClusterSize";
        histoTitle.append(TPlaneProperties::getStringForDetector(det).c_str());
        profY->SetTitle(histoTitle.c_str());
        histSaver->SaveHistogram(profY);
    }
}

void TAnalysisOfClustering::createPHDistribution(){
    bool isValid=true;
    for(UInt_t det =0;det<TPlaneProperties::getNDetectors();det++)
        isValid = (eventReader->getNClusters(det)==1)&&isValid;
    for(UInt_t det =0;det<TPlaneProperties::getNDetectors();det++){
        if(!isValid)
            continue;
        UInt_t nClusterSize = eventReader->getClusterSize(det,0);
        Float_t charge = eventReader->getCluster(det,0).getCharge(TPlaneProperties::isDiamondDetector(det));

        hPHDistribution[det]->Fill(charge,0);
        hPHDistribution[det]->Fill(charge,nClusterSize);
        //		cout<<"Fill PH histo with "<<charge<<" and Clustersize "<<nClusterSize<<endl;
        for (UInt_t cl = 1; cl <= 5; cl++){
            charge = eventReader->getCluster(det,0).getCharge(cl,TPlaneProperties::isDiamondDetector(det),true);
            hPHDistIncreasingClustersize[det]->Fill(charge,cl);
        }
    }
}




