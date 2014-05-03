//
//  TTransparentAnalysis.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 05.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//


#include "../include/TTransparentAnalysis.hh"

TTransparentAnalysis::TTransparentAnalysis(TSettings* settings, TSettings::alignmentMode mode) {
    cout<<"**********************************************************"<<endl;
    cout<<"********TTransparentAnalysis::TTransparentAnalysis********"<<endl;
    cout<<"**********************************************************"<<endl;
    // TODO Auto-generated constructor stub
    if(settings==0){
        cerr<<"Settings invalid:"<<settings<<endl;
        exit(-1);
    }
    sys = gSystem;
    setSettings(settings);
    results=0;
    settings->goToAlignmentRootDir();
    eventReader = new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(mode),settings->getEtaDistributionPath(),settings);
    // TODO: load settings!!!
    histSaver = new HistogrammSaver(settings);
    //	settings->goToTransparentAnalysisDir();
    histSaver->SetPlotsPath(settings->getTransparentAnalysisDir(mode));
    histSaver->SetRunNumber(settings->getRunNumber());
    htmlTransAna = new THTMLTransparentAnalysis(settings);
    htmlTransAna->setFileGeneratingPath(settings->getTransparentAnalysisDir(mode));
    // TODO: move these setting to the proper place
    subjectDetector = 8;
    subjectPlane = subjectDetector/2;
    if (subjectDetector%2 == 0) {
        subjectDetectorCoordinate = TPlaneProperties::X_COR;
    }
    else {
        subjectDetectorCoordinate = TPlaneProperties::Y_COR;
    }
    for (int i = 0; i < 4; i++) {
        refPlanes.push_back(i);
    }
    clusterCalcMode = TCluster::highest2Centroid;
    verbosity = settings->getVerbosity();

    initHistograms();
    cout<<"end initialise"<<endl;
    positionPrediction = 0;
    inf  = std::numeric_limits<float>::infinity();
    alignMode = mode;
    gRandom->SetSeed(1);
    cmCorrected  = false;
    if(verbosity>5) settings->diamondPattern.Print();
    xDivisions = 3;
    yDivisions = 3;
}

TTransparentAnalysis::~TTransparentAnalysis() {
    // TODO Auto-generated destructor stub
    cout<<"\n\nClosing TTransparentAnalysis"<<endl;
    analyseNonHitEvents();
    fitHistograms();
    saveHistograms();
    deleteHistograms();
    deleteFits();

    vector<vector <Float_t> > meanPulseHeights;
    vector<vector <Float_t> > mpPulseHeights;
    vector<vector <pair <Float_t,Float_t> > > resolutions;

    mpPulseHeights.push_back(vecMPLandau);
    mpPulseHeights.push_back(vecMPLandau2Highest);
    meanPulseHeights.push_back(vecMeanLandau);
    //	this->settings->res
    meanPulseHeights.push_back(vecMeanLandau2Highest);
    if(results) this->results->setPH_NoutOfN(vecMeanLandau, alignMode);
    resolutions.push_back(vecResidualChargeWeighted);
    resolutions.push_back(vecResidualHighest2Centroid);
    resolutions.push_back(vecResidualEtaCorrected);
    resolutions.push_back(vecResidualEtaCorrected_2ndGaus);

    htmlTransAna->createPulseHeightPlots(meanPulseHeights, mpPulseHeights);
    htmlTransAna->createResolutionPlots(resolutions);
    htmlTransAna->createEtaPlots();
    htmlTransAna->createEtaIntegrals();
    htmlTransAna->generateHTMLFile();
    if (eventReader!=0) delete eventReader;
    if (histSaver!=0) delete histSaver;
    if (htmlTransAna!=0) delete htmlTransAna;
    settings->goToOutputDir();
}

void TTransparentAnalysis::analyze(UInt_t nEvents, UInt_t startEvent) {
    cout<<"\n\n******************************************\n";
    cout<<    "******Start Transparent Analysis...*******\n";
    cout<<"******************************************\n\n"<<endl;
    if(verbosity>10&&verbosity%2==1){
        cout<< "Press a Key and enter to continue..."<<flush;
        char t;
        cin >>t;
    }
    predXMin = predYMin = 1e9;
    predXMax = predYMax = -1e9;
    //	usedForSiliconAlignment = 0;
    if(verbosity>6)cout<<"Current Dir: "<<sys->pwd()<<endl;
    if (nEvents+startEvent > eventReader->GetEntries()) {
        cout << "only "<<eventReader->GetEntries()<<" in tree!\n";
        nEvents = eventReader->GetEntries()-startEvent;
    }
    UInt_t newstartEvent =0;
    if(settings->getAlignmentEvents(nEvents)>startEvent){
        cout<<"startEvent:      "<<startEvent<<endl;
        cout<<"alignmentEvents: "<<settings->getAlignmentEvents(nEvents)<<endl;
        newstartEvent = TMath::Max(settings->getAlignmentEvents(nEvents),startEvent);
        cout<<"newstartEvent: "<<newstartEvent<<endl;
        cout<<"nEvents: "<<nEvents<<endl;
        nEvents -= newstartEvent-startEvent;
        startEvent = newstartEvent;
        cout<<"\nnEvents = "<<nEvents<<endl;
        cout<<"startEvent= "<<startEvent<<endl;
    }
    initClusteredHistos(startEvent,nEvents+startEvent);
    initPedestalAndNoiseHistos(nEvents+startEvent);
    initPHvsEventNoAreaPlots(startEvent,nEvents+startEvent);
    this->nEvents = nEvents;
    createEventVector(startEvent);
    cout<<"X: "<<predXMin<<" - "<<predXMax<<endl;
    cout<<"Y: "<<predYMin<<" - "<<predYMax<<endl;

    usedForAlignment += newstartEvent;
    this->printCutFlow();
    createEtaIntegrals();
    calcEtaCorrectedResiduals();
}

void TTransparentAnalysis::calcEtaCorrectedResiduals() {
    cout << "\n\ncalculating eta corrected residuals..\n";
    if (eventNumbers.size() != vecTransparentClusters.size()) {
        cout << "now we are in deep trouble!! size of eventNumbers and transparentClusters do not match!" << endl;
        return;
    }
    if (vecTransparentClusters.size()==0 || eventNumbers.size()==0) {
        cout << "oh boy.. you didn't run the analysis!" << endl;
    }
    for(UInt_t clusterSize = 0; clusterSize <TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);clusterSize++){
        vecvecResXChargeWeighted[clusterSize].clear();
        vecvecRelPos[clusterSize].clear();
        vecvecRelPos2[clusterSize].clear();
        vecvecEta[clusterSize].clear();
        vecvecEtaCMNcorrected[clusterSize].clear();
        vecvecResXHighest2Centroid[clusterSize].clear();
        vecvecResXEtaCorrected[clusterSize].clear();
        vecvecResXHighestHit[clusterSize].clear();
    }
    //	vecPredictedPosition.clear();
    //	vecRelPredictedPosition.clear();
    //	vecPredictedChannel.clear();
    //	vecChi2.clear();
    UInt_t maxSize = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);
    for (UInt_t iEvent = 0; iEvent < eventNumbers.size(); iEvent++) {
        TRawEventSaver::showStatusBar(iEvent,eventNumbers.size(),100);
        nEvent = eventNumbers.at(iEvent);
        eventReader->LoadEvent(nEvent);
        if(!this->predictPositions(false))
            continue;
        Float_t etaClusSizeOf2 = -1;
        for (UInt_t clusterSize = 0; clusterSize < maxSize; clusterSize++) {
            vecTransparentClusters[iEvent].SetTransparentClusterSize(clusterSize+1);
            if (clusterSize == 2 && false) {
                cout << "using " << hEtaIntegrals[clusterSize]->GetName() << " to fill " << hResidualEtaCorrected[clusterSize]->GetName() << endl;
                printCluster(vecTransparentClusters.at(iEvent));
            }

            //			if (clusterSize == 1) printCluster(vecTransparentClusters.at(iEvent).at(clusterSize));
            Float_t metricPosInDetSystem = eventReader->getPositionInDetSystem(subjectDetector,predXPosition,predYPosition);
            Float_t channelPosInDetSystem = settings->convertMetricToChannelSpace(subjectDetector,metricPosInDetSystem);

            Float_t resXChargeWeighted = this->getResidual(this->vecTransparentClusters.at(iEvent),cmCorrected,TCluster::chargeWeighted,hEtaIntegrals[clusterSize]);
            Float_t resXEtaCorrected = this->getResidual(this->vecTransparentClusters.at(iEvent),cmCorrected,TCluster::corEta,hEtaIntegrals[clusterSize]);
            Float_t resXHighest2Centroid = this->getResidual(this->vecTransparentClusters.at(iEvent),cmCorrected,TCluster::highest2Centroid,hEtaIntegrals[clusterSize]);
            Float_t relChannelPos = channelPosInDetSystem - (int)(channelPosInDetSystem+.5);
            Float_t resXHighestHit= this->getResidual(this->vecTransparentClusters.at(iEvent),cmCorrected,TCluster::maxValue,hEtaIntegrals[clusterSize]);

            //			Float_t relHitPos = this->predPosition- (int)(predPosition+.5);
            hResidualEtaCorrected[clusterSize]->Fill(resXEtaCorrected);//getResidual(vecTransparentClusters.at(iEvent),TCluster::corEta,hEtaIntegrals[clusterSize]));
            Int_t leftChannel=-1;
            Float_t eta = vecTransparentClusters[iEvent].getEta(leftChannel);
            Float_t etaCMNcorrected = vecTransparentClusters[iEvent].getEta(true);
            if (clusterSize == 2){
                etaClusSizeOf2 = eta;
            }
            if(verbosity>4)
                cout<<nEvent<<": "<<clusterSize<<"clusterSize: "<<channelPosInDetSystem<<"-->"<<relChannelPos<<" <-> "<<resXChargeWeighted<<", "<<resXEtaCorrected<<", "<<resXHighest2Centroid<<endl;
            vecvecRelPos[clusterSize].push_back(relChannelPos);
            vecvecRelPos2[clusterSize].push_back(relChannelPos+.5);
            vecvecEta[clusterSize].push_back(eta);
            vecvecEtaCMNcorrected[clusterSize].push_back(etaCMNcorrected);
            //			if (resXChargeWeighted > -6000)
            vecvecResXChargeWeighted[clusterSize].push_back(resXChargeWeighted);
            //			if (resXHighest2Centroid > -6000)
            vecvecResXHighest2Centroid[clusterSize].push_back(resXHighest2Centroid);
            vecvecResXHighestHit[clusterSize].push_back(resXHighestHit);
            //			if (resXEtaCorrected > -6000)
            vecvecResXEtaCorrected[clusterSize].push_back(resXEtaCorrected);
            if(clusterSize==maxSize-1){
                Float_t deltaEta = eta-etaClusSizeOf2;
                vecDeltaEta.push_back(deltaEta);
                vecRelatedEta2.push_back(etaClusSizeOf2);
                vecRelatedEta10.push_back(eta);
                vecRelatedResXEtaCorrected.push_back(resXEtaCorrected);
                Float_t signalLeftOfEta = vecTransparentClusters[iEvent].getSignalOfChannel(leftChannel-1,cmCorrected);
                Float_t signalRightOfEta = vecTransparentClusters[iEvent].getSignalOfChannel(leftChannel+2);
                Int_t highestClusterPos =  vecTransparentClusters[iEvent].getHighestHitClusterPosition();
                Float_t leftOfHighestSignal =  vecTransparentClusters[iEvent].getSignal(highestClusterPos-1,cmCorrected);
                Float_t rightOfHighestSignal =  vecTransparentClusters[iEvent].getSignal(highestClusterPos+1,cmCorrected);
                Float_t charge = vecTransparentClusters[iEvent].getCharge((UInt_t)2,cmCorrected);
                Float_t highestSignal = vecTransparentClusters[iEvent].getHighestSignal(cmCorrected);
                this->vecSignalLeftOfEta.push_back(signalLeftOfEta);
                this->vecSignalRightOfEta.push_back(signalRightOfEta);
                this->vecSignalLeftOfHighest.push_back(leftOfHighestSignal);
                this->vecSignalRightOfHighest.push_back(rightOfHighestSignal);
                this->vecClusterCharge.push_back(charge);
                this->vecHighestSignal.push_back(highestSignal);
                this->vecEta.push_back(eta);

            }
        }
    }
}

bool TTransparentAnalysis::predictPositions(bool savePrediction) {
    if (positionPrediction) delete positionPrediction;
    positionPrediction = eventReader->predictPosition(subjectPlane,refPlanes,false);
    this->predXPosition = positionPrediction->getPositionX();
    this->predYPosition = positionPrediction->getPositionY();
    //		if (verbosity > 4) cout << "predicted x position:\t" << this->predXPosition << "\ty position:\t" << this->predYPosition << endl;
    if (subjectDetector%2 == 0) {
        this->predPerpPosition = this->predYPosition;
        this->predPosition = this->predXPosition;
    }
    else {
        this->predPerpPosition = this->predXPosition;
        this->predPosition = this->predYPosition;
    }
    // TODO: position in det system
    this->positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, this->predXPosition, this->predYPosition);
    this->positionInDetSystemMetricY = eventReader->getYPositionInDetSystem(subjectDetector,predXPosition,predYPosition);
    this->positionInDetSystemChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,positionInDetSystemMetric);
    if(verbosity>5){
        cout<<"\nEventNo: "<<nEvent<<":\t"<<predPosition<<"/"<<predPerpPosition<<"--->";
        cout<<positionInDetSystemMetric<<" mum --> "<<positionInDetSystemChannelSpace<<" ch"<<flush;
    }
    //		if (verbosity > 4) cout << "position in det system:\t" << this->positionInDetSystem << endl;
    //		if (verbosity > 4)
    //			cout << "clustered analysis strip position:\t" << eventReader->getMeasured(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
    Float_t chi2 = positionPrediction->getChi2();
    if(savePrediction){
        vecPredictedPosition.push_back(positionInDetSystemChannelSpace);
        vecRelPredictedPosition.push_back(positionInDetSystemChannelSpace-(int)(positionInDetSystemChannelSpace));
        vecChi2.push_back(chi2);}

    if(chi2>settings->getTransparentChi2())
        return false;
    return true;
}

bool TTransparentAnalysis::checkPredictedRegion(UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
    // get channel and direction for clustering
    UInt_t centerChannel;
    if(verbosity>3)
        cout<<"\ncheck Pred Region: "<<nEvent<< " "<<det<<" "<<centerPosition<<" "<<clusterSize<<endl;
    int direction;
    direction = getSignedChannelNumber(centerPosition);
    centerChannel = TMath::Abs(direction);
    if (direction < 0) direction = -1;
    else direction = 1;
    //	cout<<"\t"<<direction<<" x "<< centerChannel<<endl;

    // check predicted cluster channels
    UInt_t currentChannel = centerChannel;
    for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
        direction *= -1;
        currentChannel += direction * iChannel;
        if (currentChannel < 0) {
            if (verbosity > 5)
                cout << "\tchannel " << currentChannel << " is not on this detector.." << endl;
            regionNotOnPlane++;
            return false;
        }
        if (currentChannel > TPlaneProperties::getNChannels(det)-1) {
            if (verbosity > 5)
                cout << "\tchannel " << currentChannel << " is not on this detector.." << endl;
            regionNotOnPlane++;
            return false;
        }
        if (this->settings->getDet_channel_screen(det).isScreened(currentChannel) == true) {
            if (verbosity > 5) cout << "\tchannel " << currentChannel << " is screened.." << endl;
            screenedChannel++;
            return false;
        }
        if (eventReader->isSaturated(det, currentChannel) == true) {
            if (verbosity > 5) cout << "\tchannel " << currentChannel << " has saturated.." << endl;
            saturatedChannel++;
            return false;
        }
    }
    return true;
}

// TODO: avoid wrong channel numbers (>128, <0)
/**
 *	returns the next channel number including a sign: + if pos - (int)pos <.5 else minus
 *	different approach:
 *	 ( 1+2*( (int)pos-(int)(pos+.5) ) ) * int (pos+.5)
 * @param position
 * @author Lukas Baeni
 * @return
 */
int TTransparentAnalysis::getSignedChannelNumber(Float_t position) {
    if (position < 0) return -9999;
    UInt_t channel = 0;
    int direction;
    if (position-(int)position < 0.5) {
        channel = (UInt_t)position;
        direction = 1;
    }
    else {
        channel = (UInt_t)position+1;
        direction = -1;
    }
    return direction * channel;
}

void TTransparentAnalysis::setSettings(TSettings* settings) {
    this->settings=settings;
}

void TTransparentAnalysis::initHistograms() {
    cout<<"initHistos"<<endl;
    UInt_t bins=512;
    vecvecResXChargeWeighted.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecResXHighest2Centroid.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecResXHighestHit.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecResXEtaCorrected.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecRelPos.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecRelPos2.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecEta.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecvecEtaCMNcorrected.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecVecLandau.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecVecPh2Highest.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    vecVecFidCutX.clear();
    vecVecFidCutY.clear();
    vecPredX.clear();
    vecPredY.clear();
    hEtaIntegrals.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    hSelectedTracksAvrgSiliconHitPos = new TH2F("hSelectedTracksAvrgSiliconHitPos","hSelectedTracksAvrgSiliconHitPos",1024,0 ,256,1024,0,256);
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        vecVecLandau.at(clusterSize).clear();
        vecvecEtaCMNcorrected.at(clusterSize).clear();
        vecvecRelPos.at(clusterSize).clear();
        vecvecRelPos2.at(clusterSize).clear();
        vecVecPh2Highest.at(clusterSize).clear();
        vecvecResXEtaCorrected.at(clusterSize).clear();
        vecvecResXChargeWeighted.at(clusterSize).clear();
        // TODO: take care of histogram names and bins!!
        stringstream histNameLandau,histNameLandau1Highest, histNameLandau2Highest, histNameEta, histNameResidualChargeWeighted, histNameResidualHighest2Centroid, histNameResidualEtaCorrected;
        // TODO: histogram naming!!
        histNameLandau << "hDiaTranspAnaPulseHeightOf" << clusterSize+1 << "Strips";
        histNameLandau1Highest<< "hDiaTranspAnaPulseHeightOfHighestIn"<<clusterSize+1<<"Strips";
        histNameLandau2Highest << "hDiaTranspAnaPulseHeightOf2HighestIn" << clusterSize+1 << "Strips";
        histNameEta << "hDiaTranspAnaEta2HighestIn" << clusterSize+1 << "Strips";
        histNameResidualChargeWeighted << "hDiaTranspAnaResidualChargeWeightedIn" << clusterSize+1 << "StripsMinusPred";
        histNameResidualHighest2Centroid << "hDiaTranspAnaResidualHighest2CentroidIn" << clusterSize+1 << "StripsMinusPred";
        histNameResidualEtaCorrected << "hDiaTranspAnaResidualEtaCorrectedIn" << clusterSize+1 << "StripsMinusPred";
        TString nameResVsHitChargeWeighted = TString::Format("hDiaTransAnaResVsHitChargeWeightedIn%d",clusterSize+1);
        TString nameResVsHitEtaCor = TString::Format("hDiaTransAnaResVsHitEtaCorIn%d",clusterSize+1);
        TString nameResVsHitHeighest2Centroid = TString::Format("hDiaTransAnaResVsHitHigehst2CentroidIn%d",clusterSize+1);
        TString histNameResidualHighestHit = TString::Format("hDiaTranspAnaResidualHighestHitIn%dStripsMinusPred",clusterSize+1);

        hLandau.push_back(new TH1F(histNameLandau.str().c_str(),histNameLandau.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
        hLandau2Highest.push_back(new TH1F(histNameLandau2Highest.str().c_str(),histNameLandau2Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));
        this->hLandau1Highest.push_back(new TH1F(histNameLandau1Highest.str().c_str(),histNameLandau1Highest.str().c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector)));

        hEta.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
        histNameEta.str("");
        histNameEta.clear();
        histNameEta << "hDiaTranspAnaEtaCMNcorrected2HighestIn" << clusterSize+1 << "Strips";
        hEtaCMNcorrected.push_back(new TH1F(histNameEta.str().c_str(),histNameEta.str().c_str(),bins,0,1));
        ///TODO: PitchWidth Plot width
        hResidualChargeWeighted.push_back(new TH1F(histNameResidualChargeWeighted.str().c_str(),histNameResidualChargeWeighted.str().c_str(),bins,-2.5*50,2.5*50));
        hResidualHighest2Centroid.push_back(new TH1F(histNameResidualHighest2Centroid.str().c_str(),histNameResidualHighest2Centroid.str().c_str(),bins,-2.5*50,2.5*50));
        hResidualEtaCorrected.push_back(new TH1F(histNameResidualEtaCorrected.str().c_str(),histNameResidualEtaCorrected.str().c_str(),bins,-2.5*50,2.5*50));
        hResidualHighestHit.push_back(new TH1F(histNameResidualHighestHit,histNameResidualHighestHit,bins,-2.5*50,2.5*50));

        //		hResidualVsHitPositionChargeWeighted.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted));
        //		hResidualVsHitPositionEtaCorrected.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted));
        //		hResidualVsHitPositionHigehest2Centroid.push_back(new TH2F(nameResVsHitChargeWeighted,nameResVsHitChargeWeighted,));

    }
    hLandauMean = new TH1F("hDiaTranspAnaPulseHeightMean","hDiaTranspAnaPulseHeightMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
    hLandauMP = new TH1F("hDiaTranspAnaPulseHeightMP","hDiaTranspAnaPulseHeightMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
    hLandau2HighestMean = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMean","hDiaTranspAnaPulseHeightOf2HighestMean",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
    hLandau2HighestMP = new TH1F("hDiaTranspAnaPulseHeightOf2HighestMP","hDiaTranspAnaPulseHeightOf2HighestMP",TPlaneProperties::getMaxTransparentClusterSize(subjectDetector),0.5,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+0.5);
    hPredictedPositionInStrip = new TH1F("hPredictedPositionInStrip","hPredictedPositionInStrip",2,-1.5,1.5);
}

void TTransparentAnalysis::fillClusteredHistos(){
    Int_t nClusters = eventReader->getNClusters(subjectDetector);
    hNClusteres_Clustered->Fill(nClusters);
    if (nClusters != 1)
        return;
    TCluster clusteredCluster = eventReader->getCluster(subjectDetector,0);

    UInt_t clusterSize = clusteredCluster.size();
    Float_t clusteredPH = clusteredCluster.getCharge(cmCorrected,false);
    hLandauVsClusterSize_Clustered->Fill(clusteredPH,clusterSize);
    hClusterSize_Clustered->Fill(clusterSize);
    if (clusterSize <= 2)
        hLandauVsEventNo_Clustered->Fill(nEvent,clusteredPH);

}

void TTransparentAnalysis::fillHistograms() {
    hSelectedTracksAvrgSiliconHitPos->Fill(eventReader->getFiducialValueX(),eventReader->getFiducialValueY());
    vecVecFidCutX.push_back(eventReader->getFiducialValueX());
    vecVecFidCutY.push_back(eventReader->getFiducialValueY());
    vecPredX.push_back(predXPosition);
    vecPredY.push_back(predYPosition);
    vecPredictedChannel.push_back(positionInDetSystemChannelSpace);
    vecPredictedDetectorPositionY.push_back(positionInDetSystemMetricY);
    fillPedestalsAndNoiseHistos();
    UInt_t area = GetHitArea(settings,eventReader->getFiducialValueX(),eventReader->getFiducialValueY(),xDivisions,yDivisions);
    UInt_t maxSize = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);
    fillClusteredHistos();
    for (UInt_t clusterSize = 0; clusterSize < maxSize; clusterSize++) {
        transparentClusters.SetTransparentClusterSize(clusterSize+1);
        Float_t charge = this->transparentClusters.getCharge(cmCorrected);

        Float_t chargeOfTwo = this->transparentClusters.getCharge(2,cmCorrected);

        fillPHvsEventNoAreaPlots(area,clusterSize+1,charge,chargeOfTwo);
        vecVecLandau[clusterSize].push_back(charge);
        hLandau[clusterSize]->Fill(charge);
        hLandau2Highest[clusterSize]->Fill(chargeOfTwo);
        vecVecPh2Highest[clusterSize].push_back(chargeOfTwo);

        Float_t eta = this->transparentClusters.getEta();

        Float_t etaCMN = this->transparentClusters.getEta(true);
        if(eta>0&&eta<1)
            hEta[clusterSize]->Fill(eta);
        else if(verbosity>3)
            this->transparentClusters.Print();
        if(etaCMN>0&&etaCMN<1)
            hEtaCMNcorrected[clusterSize]->Fill(etaCMN);

        if (clusterSize>2){
            Int_t highestSignalClusterPos = this->transparentClusters.getHighestHitClusterPosition();
            if(highestSignalClusterPos<0){
                cout<<nEvent<<"errror5: highestSignalClusterPos: "<<highestSignalClusterPos<<endl;
                this->transparentClusters.Print();
            }
            //			Float_t leftSig = transparentClusters[clusterSize].getSignal(highestSignalClusterPos-1);
            //			Float_t rightSig= transparentClusters[clusterSize].getSignal(highestSignalClusterPos+1);

        }
        //		if (clusterSize == 1 /* && this->transparentClusters[clusterSize].getCharge() != this->transparentClusters[clusterSize].getCharge(2,false)*/) printCluster(this->transparentClusters[clusterSize]);
        //		if (clusterSize > 0 && this->transparentClusters[clusterSize].getEta() < 0) printCluster(this->transparentClusters[clusterSize]);
        // TODO: why is the eta distribution for 2 channel clusters more symmetric than for 3 and more channel clusters?
        //		if (clusterSize == 2 && this->transparentClusters[clusterSize-1].getEta() != this->transparentClusters[clusterSize].getEta()) {
        //			if (this->transparentClusters[clusterSize-1].getHighestSignalChannel()!=this->transparentClusters[clusterSize].getHighestSignalChannel()
        //				&&
        //				this->transparentClusters[clusterSize-1].getHighest2Centroid(cmCorrected)!=this->transparentClusters[clusterSize].getHighest2Centroid(cmCorrected)) {
        //			printCluster(this->transparentClusters[clusterSize-1]);
        //			printCluster(this->transparentClusters[clusterSize]);
        //			}
        //		}
        Float_t relPos =this->predPosition-(int)(this->predPosition+.5);
        Float_t residualCW =this->getResidual(this->transparentClusters,cmCorrected,TCluster::chargeWeighted,hEtaIntegrals[clusterSize]);
        Float_t residualH2C = this->getResidual(this->transparentClusters,cmCorrected,TCluster::highest2Centroid,hEtaIntegrals[clusterSize]);
        Float_t resXHighestHit= this->getResidual(this->transparentClusters,cmCorrected,TCluster::maxValue,hEtaIntegrals[clusterSize]);

        //		if(hResidualChargeWeightedVsEstimatedHitPosition==0)
        //			hResidualChargeWeightedVsEstimatedHitPosition->Fill(residualCW,relPos,clusterSize);
        //		if(hResidualHighest2CentroidVsEstimatedHitPosition)
        //			hResidualHighest2CentroidVsEstimatedHitPosition>Fill(residualH2C,relPos,clusterSize);
        if (clusterSize+1 != transparentClusters.getClusterSize()) {
            cout << "wrong cluster size!" << endl;
            cout << "clusterSize+1 = " << clusterSize+1 << "\ttransparentClusters[clusterSize].getClusterSize() = " << transparentClusters.getClusterSize() << endl;
        }
        vecvecResXChargeWeighted[clusterSize].push_back(residualCW);
        vecvecResXHighest2Centroid[clusterSize].push_back(residualH2C);
        vecvecResXHighestHit[clusterSize].push_back(resXHighestHit);
        vecvecRelPos[clusterSize].push_back(relPos);
        vecvecRelPos2[clusterSize].push_back(relPos+.5);
        vecvecEta[clusterSize].push_back(eta);
        vecvecEtaCMNcorrected[clusterSize].push_back(etaCMN);

        hResidualChargeWeighted[clusterSize]->Fill(residualCW);
        hResidualHighest2Centroid[clusterSize]->Fill(residualH2C);
        hResidualHighestHit[clusterSize]->Fill(resXHighestHit);
        if(predXPosition<predXMin)
            predXMin = predXPosition;
        if(predXPosition>predXMax)
            predXMax = predXPosition;
        if(predYPosition<predYMin)
            predYMin = predYPosition;
        if(predYPosition>predYMax)
            predYMax = predYPosition;

    }
    Float_t eventNo = eventReader->getEvent_number();
    vectorEventNo.push_back(eventNo);
    Float_t cmn = eventReader->getCMNoise();
    vectorCMN.push_back(cmn);
    //	hPredictedPositionInStrip->Fill();
}

TF1* TTransparentAnalysis::doGaussFit(TH1F *histo) {
    //	TH1* histo = (TH1*)htemp->Clone();
    if (histo->GetEntries()==0) return 0;
    TF1* histofitx = new TF1("histofitx","gaus",histo->GetMean()-2*histo->GetRMS(),histo->GetMean()+2*histo->GetRMS());
    histofitx->SetLineColor(kBlue);
    histo->Fit(histofitx,"rQ");
    return histofitx;
}

TF1* TTransparentAnalysis::doDoubleGaussFit(TH1F *histo){
    if (!histo)return 0;
    if (histo->GetEntries()==0) return 0;
    Float_t mean = histo->GetMean();
    Float_t sigma = histo->GetRMS();
    Float_t max = histo->GetBinContent(histo->GetMaximumBin());
    Float_t pw = settings->getPitchWidth(subjectDetector);
    int nSigmas = 4;
    Float_t xmin = TMath::Min(-pw,mean-nSigmas*sigma);
    Float_t xmax = TMath::Max(pw,mean+nSigmas*sigma);
    TF1* histofitx = new TF1("fDoubleGaus","gaus(0)+gaus(3)",xmin,xmax);
    histofitx->SetParLimits(0,max/10,max);
    histofitx->SetParLimits(1,mean-2*sigma,mean+2*sigma);
    histofitx->SetParLimits(2,sigma/10,2*sigma);
    histofitx->SetParLimits(3,max/20,max/2);
    histofitx->SetParLimits(4,mean-2*sigma,mean+2*sigma);
    histofitx->SetParLimits(5,sigma/10,4*sigma);
    histofitx->SetParameters(.75*max,mean,sigma/5,.1*max,mean,sigma/4);
    histofitx->SetParNames("C_{0}","#mu_{0}","#sigma_{0}","C_{1}","#mu_{1}","#sigma_{1}");
    histofitx->SetLineColor(kBlue);
    histofitx->SetNpx(1000);
    histo->Fit(histofitx,"rq");
    return histofitx;
}

void TTransparentAnalysis::createEtaIntegrals() {
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        stringstream histName;
        histName << "hDiaTranspAnaEtaIntegral2HighestIn"<<clusterSize+1<<"Strips";
        if (hEtaIntegrals.at(clusterSize))
            delete hEtaIntegrals.at(clusterSize);
        hEtaIntegrals.at(clusterSize) = (TClustering::createEtaIntegral(hEta[clusterSize], histName.str()));
    }
}


void TTransparentAnalysis::createEfficiencyPlots(TH1F *hLandau){
    if (!hLandau) cout<<"TTransparentAnalysis::createEfficiencyPlots: HLandau Histo Not Valid..."<<endl;
    TString name = TString::Format("hEfficiency_%s",hLandau->GetName());
    TH1F* hEfficiency = new TH1F(name,name,settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector) );
    Float_t nentries = hLandau->GetEntries();
    Float_t integral = 0;
    for(Int_t bin = 1;bin <= hLandau->GetNbinsX();bin++){
        integral += hLandau->GetBinContent(bin);
        hEfficiency->SetBinContent(bin, (1.-integral/nentries)*100);
    }
    if(hEfficiency){
        hEfficiency->GetXaxis()->SetTitle("PH / adc counts");
        hEfficiency->GetYaxis()->SetTitle("efficientcy / %");
    }
    //	TCutG *MP = new TCutG("gMP",1);
    histSaver->SaveHistogram(hEfficiency);
    if(hEfficiency) delete hEfficiency;
}


void TTransparentAnalysis::fitHistograms() {
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        if (verbosity) cout<<"Create EfficencyPlots "<<clusterSize<<endl;
        createEfficiencyPlots(hLandau2Highest[clusterSize]);
        createEfficiencyPlots(hLandau[clusterSize]);
        createEfficiencyPlots(hLandau1Highest[clusterSize]);
        vector<Float_t> vecResChargeWeighted = vecvecResXChargeWeighted[clusterSize];

        TString name;
        name = TString::Format("hDiaTranspAnaResidualChargeWeightedIn%02dStripsMinusPred",clusterSize+1);
        if(hResidualChargeWeighted[clusterSize])
            delete hResidualChargeWeighted[clusterSize];
        hResidualChargeWeighted[clusterSize] = histSaver->CreateDistributionHisto((string)name, vecResChargeWeighted,8192,HistogrammSaver::maxWidth,-5000);
        Float_t plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
        hResidualChargeWeighted[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
        hResidualChargeWeighted[clusterSize]->GetXaxis()->SetTitle("Residual, ChargeWeighted / #mum");
        hResidualChargeWeighted[clusterSize]->GetYaxis()->SetTitle("number of entries #");

        name = TString::Format("hDiaTranspAnaResidualHighest2CentroidIn%02dStripsMinusPred",clusterSize+1);
        if(hResidualHighest2Centroid[clusterSize] )
            delete hResidualHighest2Centroid[clusterSize] ;
        hResidualHighest2Centroid[clusterSize] = histSaver->CreateDistributionHisto((string)name, vecvecResXHighest2Centroid[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
        plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
        hResidualHighest2Centroid[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
        hResidualHighest2Centroid[clusterSize]->GetXaxis()->SetTitle("Residual, highest 2 centroid / #mum");
        hResidualHighest2Centroid[clusterSize]->GetYaxis()->SetTitle("number of entries #");

        name=TString::Format("hDiaTranspAnaResidualEtaCorrectedIn%02dStripsMinusPred",clusterSize+1);
        if(hResidualEtaCorrected[clusterSize])
            delete hResidualEtaCorrected[clusterSize];
        hResidualEtaCorrected[clusterSize] = histSaver->CreateDistributionHisto((string)name, vecvecResXEtaCorrected[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
        plotWidth = settings->getPitchWidth(subjectDetector);
        hResidualEtaCorrected[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
        hResidualEtaCorrected[clusterSize]->GetXaxis()->SetTitle("Residual #eta corrected / #mum");
        hResidualEtaCorrected[clusterSize]->GetYaxis()->SetTitle("number of entries #");

        name=TString::Format("hDiaTranspAnaResidualHighestHitIn%02dStripsMinusPred",clusterSize+1);
        if(hResidualHighestHit[clusterSize] )
            delete hResidualHighestHit[clusterSize] ;
        hResidualHighestHit[clusterSize] = histSaver->CreateDistributionHisto((string)name, vecvecResXHighestHit[clusterSize],8192,HistogrammSaver::maxWidth,-5000);
        plotWidth = 1.5 * settings->getPitchWidth(subjectDetector);
        hResidualHighestHit[clusterSize]->GetXaxis()->SetRangeUser(-plotWidth,plotWidth);
        hResidualHighestHit[clusterSize]->GetXaxis()->SetTitle("Residual, highest Hit / #mum");
        hResidualHighestHit[clusterSize]->GetYaxis()->SetTitle("number of entries #");
        //,
        //,1024,HistogrammSaver::maxWidth,-6000);
        // fit histograms
        name = hLandau[clusterSize]->GetName();
        name.Append("_fixedNoise");
        TH1F* hLandauFixedNoiseFit = (TH1F*) hLandau[clusterSize]->Clone(name);
        hLandauFixedNoiseFit->SetTitle(name);
        Float_t noise = cmCorrected?noiseWidthsCMN[clusterSize]:noiseWidths[clusterSize];
        cout<<"NOISE of "<<clusterSize<<": "<<noise<<endl;
        TF1* fitFixedNoise = landauGauss->doLandauGaussFitFixedNoise(hLandauFixedNoiseFit,noise,true);
        //		histSaver->SaveHistogram(hLandauFixedNoiseFit);
        hLandauFixedNoise.push_back(hLandauFixedNoiseFit);
        fitLandauFixedNoise.push_back(fitFixedNoise);
        cout<<"#"<<flush;

        name = hLandau2Highest[clusterSize]->GetName();
        name.Append("_fixedNoise");
        TH1F* hLandau2HighestFixedNoiseFit = (TH1F*) hLandau2Highest[clusterSize]->Clone(name);
        TF1* fit2HighestFixedNoise = landauGauss->doLandauGaussFitFixedNoise(hLandau2HighestFixedNoiseFit,noise,true);
        //		histSaver->SaveHistogram(hLandau2HighestFixedNoiseFit);
        hLandau2HighestFixedNoise.push_back(hLandau2HighestFixedNoiseFit);
        fitLandau2HighestFixedNoise.push_back(fit2HighestFixedNoise);
        cout<<"$"<<flush;



        TF1* fit = landauGauss->doLandauGaussFit(hLandau[clusterSize],true);
        if(fit==0){cout<<"PROBLEM with fit..."<<clusterSize<<endl;}
        fitLandau.push_back(fit);
        fitLandau2Highest.push_back(landauGauss->doLandauGaussFit(hLandau2Highest[clusterSize],true));
        if(clusterSize == TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)-1){
            Float_t mean;
            if(hLandau2Highest[clusterSize]) mean = hLandau2Highest[clusterSize]->GetMean();
            else mean = -1;
            Float_t mp;
            TF1* fit = fitLandau2Highest.back();
            if(fit) mp = fit->GetParameter(1);
            else mp = -1;
            Float_t width;
            if (fit) width = fit->GetParameter(0);
            else width = -1;
            Float_t gSigma;
            if (fit) gSigma = fit->GetParameter(3);
            else gSigma = -1;
            if(results){
                results->setPH_2outOf10(mean,mp,width,gSigma,alignMode);
            }
            else cout<<"setPH_2outOf10 DIDN'T WORK!!!"<<endl;
        }
        fitResidualChargeWeighted.push_back(doGaussFit(hResidualChargeWeighted[clusterSize]));
        fitResidualHighest2Centroid.push_back(doGaussFit(hResidualHighest2Centroid[clusterSize]));
        fitResidualEtaCorrected.push_back(doDoubleGaussFit(hResidualEtaCorrected[clusterSize]));
        if(clusterSize+1 >= TPlaneProperties::getMaxTransparentClusterSize(subjectDetector))
            saveResolutionPlot(hResidualEtaCorrected[clusterSize],clusterSize,"Normal");
        // save fit parameters
        vecMPLandau.push_back(fitLandau[clusterSize]->GetParameter(1));
        vecMPLandau2Highest.push_back(fitLandau2Highest[clusterSize]->GetParameter(1));
        hLandauMP->SetBinContent(clusterSize+1,fitLandau[clusterSize]->GetParameter(1));
        hLandau2HighestMP->SetBinContent(clusterSize+1,fitLandau2Highest[clusterSize]->GetParameter(1));
        vecMeanLandau.push_back(hLandau[clusterSize]->GetMean());
        vecMeanLandau2Highest.push_back(hLandau2Highest[clusterSize]->GetMean());
        hLandauMean->SetBinContent(clusterSize+1,hLandau[clusterSize]->GetMean());
        hLandau2HighestMean->SetBinContent(clusterSize+1,hLandau2Highest[clusterSize]->GetMean());
        pair <Float_t,Float_t> tempPair,tempPair2;
        if (fitResidualChargeWeighted[clusterSize]!=0) {
            tempPair.first = fitResidualChargeWeighted[clusterSize]->GetParameter(1);
            tempPair.second = fitResidualChargeWeighted[clusterSize]->GetParameter(2);
        }
        else {
            tempPair.first = 0;
            tempPair.second = 0;
        }
        vecResidualChargeWeighted.push_back(tempPair);
        if (fitResidualHighest2Centroid[clusterSize]!=0) {
            tempPair.first = fitResidualHighest2Centroid[clusterSize]->GetParameter(1);
            tempPair.second = fitResidualHighest2Centroid[clusterSize]->GetParameter(2);
        }
        else {
            tempPair.first = 0;
            tempPair.second = 0;
        }
        cout<<"%"<<flush;
        vecResidualHighest2Centroid.push_back(tempPair);
        if (fitResidualEtaCorrected[clusterSize]!=0) {
            tempPair.first = fitResidualEtaCorrected[clusterSize]->GetParameter(1);
            tempPair.second = fitResidualEtaCorrected[clusterSize]->GetParameter(2);
            tempPair2.first = fitResidualEtaCorrected[clusterSize]->GetParameter(4);
            tempPair2.second = fitResidualEtaCorrected[clusterSize]->GetParameter(5);
        }
        else {
            tempPair.first = 0;
            tempPair.second = 0;
        }
        if (tempPair2.second>tempPair.second){
            vecResidualEtaCorrected.push_back(tempPair);
            vecResidualEtaCorrected_2ndGaus.push_back(tempPair2);
        }
        else{
            vecResidualEtaCorrected.push_back(tempPair2);
            vecResidualEtaCorrected_2ndGaus.push_back(tempPair);
        }
    }
    hLandauMean->Scale(1./hLandauMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
    hLandauMP->Scale(1./hLandauMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
    hLandau2HighestMean->Scale(1./hLandau2HighestMean->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
    hLandau2HighestMP->Scale(1./hLandau2HighestMP->GetBinContent(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)));
    cout<<"#"<<endl;
}

void TTransparentAnalysis::analyseEtaDistribution(TH1F* hEtaDist){
    if(!hEtaDist)
        return;
    if(hEtaDist->GetEntries()==0)
        return;
    float threshold = 0.1;
    int n=0;
    int ntries = 0;
    int maxTries =30;
    while (n!=2&&ntries < maxTries){
        n = hEtaDist->ShowPeaks(3,"nobackground",threshold);
        if(n<2){
            threshold*=.9;
            //				cout<<ntries<<"-"<<n<<" ==> lowering threshold "<<threshold<<endl;
        }
        else if(n>2){
            threshold*=1.1;
            //				cout<<ntries<<"-"<<n<<" ==> higher threshold "<<threshold<<endl;
        }
        ntries++;
    }
    TList *functions = hEtaDist->GetListOfFunctions();
    TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
    if(verbosity) cout<<hEtaDist->GetName()<<" - "<<ntries<<endl;
    if (!pm){
        if(functions) functions->Print();
        return;
    }
    for (int i=0;i< pm->GetN();i++)
        if(verbosity) cout<<"\t"<<i<<"\t"<<pm->GetX()[i]*100.<<": "<<pm->GetY()[i]<<"\n";
    if(pm->GetN()==2){
        Float_t x_0 = pm->GetX()[0];
        Float_t x_1 = pm->GetX()[1];
        Float_t y_0 = pm->GetY()[0];
        Float_t y_1 = pm->GetY()[1];
        if(x_0>x_1){
            x_1 = x_0;
            x_0 = pm->GetX()[1];
            y_0 = y_1;
            y_1 = pm->GetY()[0];
        }
        if(x_1>.5)
            x_1 = 1-x_1;
        x_0*=100.;
        x_1*=100.;
        if(verbosity)cout<<"\t\t"<<x_0<<" - "<< x_1 <<"\t"<<(x_0-x_1)<<"\t"<<x_0/x_1<<"\t"<<(x_0-x_1)*100./x_0<<endl;
        if(verbosity)cout<<"\t\t"<<y_0<<" - "<< y_1 <<"\t"<<(y_0-y_1)<<"\t"<<y_0/y_1<<"\t"<<(y_0-y_1)*100/y_0<<endl;
    }
    if(verbosity)cout<<"\n"<<flush;
}

void TTransparentAnalysis::analyseEtaDistributions(){
    cout<<" TTransparentAnalysis::analyseEtaDistributions"<<endl;
    stringstream name;
    name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_etaOf2");
    if(vecRelatedEta2.size()!=vecDeltaEta.size()){
        cout<<"Something is wrong with vedDeltaEta Size"<<flush;
        char t;
        cin>>t;
    }
    TH2F* hDeltaEta = histSaver->CreateScatterHisto(name.str(),vecRelatedEta2,vecDeltaEta,400,400,-1,1,0,1);
    if(hDeltaEta){
        hDeltaEta->GetXaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
        hDeltaEta->GetYaxis()->SetTitle("#eta_{2 of 2}");
    }
    if(hDeltaEta)
        histSaver->SaveHistogram(hDeltaEta,false);

    if(hDeltaEta)delete hDeltaEta;

    name.str("");name.clear();
    name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_ResidualEtaCorrectedIn10");
    Float_t pw = settings->getDiamondPitchWidth();
    TH2F* hDeltaEtaVsResidual = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedResXEtaCorrected,512,400,-2*pw,2*pw,-1,1);
    if(hDeltaEtaVsResidual){
        hDeltaEtaVsResidual->GetXaxis()->SetTitle("Residual Eta Correct, 2 of 10 / #mum");
        hDeltaEtaVsResidual->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
    }
    if(hDeltaEtaVsResidual)
        histSaver->SaveHistogram(hDeltaEtaVsResidual,false);

    name.str("");name.clear();
    name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_etaOf10");
    TH2F* hDeltaEtaVsEta = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedEta10,512,512,0,1,0,1);
    if(hDeltaEtaVsEta){
        hDeltaEtaVsEta->GetXaxis()->SetTitle("#eta_{2 of 10}");
        hDeltaEtaVsEta->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
    }
    if(hDeltaEtaVsEta)
        histSaver->SaveHistogram(hDeltaEtaVsEta,false);
    Float_t maxDelta = .199999;
    Int_t bin1 = hDeltaEtaVsEta->GetYaxis()->FindBin(-maxDelta);
    Int_t bin2 = hDeltaEtaVsEta->GetYaxis()->FindBin(+maxDelta);
    TString hName = TString::Format("hEtaIn10_DeltaEta_below_020");
    TH1F* hEtaBoundedEtaCorrected = (TH1F*)hDeltaEtaVsEta->ProjectionX(hName,bin1,bin2);
    if(hEtaBoundedEtaCorrected)
        hEtaBoundedEtaCorrected->SetTitle(TString::Format("#eta_{2 of 10}, |#Delta#eta| < 0.2"));
    else
        cout<<"hEtaBoundedEtaCorrecte = 0"<<endl;
    histSaver->SaveHistogram(hEtaBoundedEtaCorrected);
    if(hEtaBoundedEtaCorrected)delete hEtaBoundedEtaCorrected;

    //	name.str("");name.clear();
    //	name<<TString::Format("hEtaOf10_minus_EtaOf2_vs_ResidualEtaCorrectedIn10");
    //	TH2F* hDeltaEtaVsResidual = histSaver->CreateScatterHisto(name.str(),vecDeltaEta,vecRelatedResXEtaCorrected,512,400,-2*pw,2*pw,-1,1);
    //	if(hDeltaEtaVsResidual){
    //		hDeltaEtaVsResidual->GetXaxis()->SetTitle("Residual Eta Correct, 2 of 10 / #mum");
    //		hDeltaEtaVsResidual->GetYaxis()->SetTitle("#Delta#eta = #eta_{2 of 10} - #eta_{2 of 2}");
    //	}
    //	histSaver->SaveHistogram(hDeltaEtaVsResidual,false);

    maxDelta = .199999;
    bin1 = hDeltaEtaVsResidual->GetYaxis()->FindBin(-maxDelta);
    bin2 = hDeltaEtaVsResidual->GetYaxis()->FindBin(+maxDelta);
    hName = TString::Format("hResidualEtaCorrectedIn10_DeltaEta_below_020");
    TH1F* hEtaBoundedEtaCorrectedResidual = (TH1F*)hDeltaEtaVsResidual->ProjectionX(hName,bin1,bin2);
    if(hEtaBoundedEtaCorrectedResidual)
        hEtaBoundedEtaCorrectedResidual->SetTitle(TString::Format("Residual #eta corrected in 10 strips, |#Delta#eta| < 0.2"));
    histSaver->SaveHistogram(hEtaBoundedEtaCorrectedResidual);
    if(hDeltaEtaVsResidual)delete hDeltaEtaVsResidual;
    if(hEtaBoundedEtaCorrectedResidual) delete hEtaBoundedEtaCorrectedResidual;

    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        TH1F* hEtaDist = hEta[clusterSize];
        analyseEtaDistribution(hEtaDist);
        hEtaDist = hEtaCMNcorrected[clusterSize];
        analyseEtaDistribution(hEtaDist);
    }

    name.str("");name.clear();
    name<<"hEtaVsSignalLeftOfEta";
    TH2F* histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalLeftOfEta,this->vecEta);
    if(histo2d){
        histo2d->GetYaxis()->SetTitle("Signal left of #eta");
        histo2d->GetXaxis()->SetTitle("#eta");
        histSaver->SaveHistogram(histo2d);
        delete histo2d;
    }

    name.str("");name.clear();
    name<<"hEtaVsSignalRightOfEta";
    histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalRightOfEta,this->vecEta);
    if(histo2d){
        histo2d->GetYaxis()->SetTitle("Signal right of #eta");
        histo2d->GetXaxis()->SetTitle("#eta");
        histSaver->SaveHistogram(histo2d);
        delete histo2d;
    }
    vector<Float_t> vecRightFactor, vecLeftFactor,vecRightOverHighest,vecLeftOverHighest;
    for(UInt_t i=0;i<vecSignalLeftOfHighest.size()&&i<vecSignalRightOfHighest.size()&&i<vecClusterCharge.size();i++){
        vecRightFactor.push_back(vecSignalRightOfHighest.at(i)/vecClusterCharge.at(i));
        vecLeftFactor.push_back(vecSignalLeftOfHighest.at(i)/vecClusterCharge.at(i));
        vecLeftOverHighest.push_back(vecSignalLeftOfHighest.at(i)/vecHighestSignal.at(i));
        vecRightOverHighest.push_back(vecSignalRightOfHighest.at(i)/vecHighestSignal.at(i));
    }
    name.str("");name.clear();
    name<<"hSignalLeftOverHighestVsSignalRightOverHighest";
    histo2d = histSaver->CreateScatterHisto(name.str(),vecRightOverHighest,vecLeftOverHighest);
    if(histo2d){
        histo2d->GetXaxis()->SetTitle("signal left of highest signal over highest signal");
        histo2d->GetYaxis()->SetTitle("signal right of highest signal over highest signal");
        histSaver->SaveHistogram(histo2d);
        delete histo2d;
    }

    name.str("");name.clear();
    name<<"hSignalRightVsHighestSignal";
    histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalRightOfHighest,vecHighestSignal);
    if(histo2d){
        histo2d->GetXaxis()->SetTitle("highest signal");
        histo2d->GetYaxis()->SetTitle("signal right of highest signal ");
        histSaver->SaveHistogram(histo2d);
        delete histo2d;
    }
    name.str("");name.clear();
    name<<"hSignalLeftVsHighestSignal";
    histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalLeftOfHighest,vecHighestSignal);
    if(histo2d){
        histo2d->GetXaxis()->SetTitle("highest signal");
        histo2d->GetYaxis()->SetTitle("signal left of highest signal ");
        histSaver->SaveHistogram(histo2d);
        delete histo2d;
    }
    name.str("");name.clear();
    name<<"hSignalLeftOfHighest";
    TH1F* histoLeft = histSaver->CreateDistributionHisto(name.str(),vecLeftFactor);
    histSaver->SaveHistogram(histoLeft);

    name.str("");name.clear();
    name<<"hSignalRightOfHighest";
    TH1F* histoRight = histSaver->CreateDistributionHisto(name.str(),vecRightFactor);
    histSaver->SaveHistogram(histoRight);

    name.str("");name.clear();
    name<<"cSignalNextToHighest";
    if(histoLeft)
        histoLeft->SetLineColor(kBlue);
    else
        cout<<"histoLeft = 0"<<endl;
    if (histoRight)
        histoRight->SetLineColor(kRed);
    else
        cout<<"histoRight = 0"<<endl;
    Float_t max = TMath::Max(histoLeft->GetMaximum(),histoRight->GetMaximum());
    if (histoLeft) histoLeft->SetMaximum(max);
    if (histoRight) histoRight->SetMaximum(max);
    histSaver->SaveTwoHistos(name.str(),histoLeft,histoRight,1.);
    if(histoLeft) delete histoLeft;
    if(histoRight) delete histoRight;

    name.str("");name.clear();
    name<<"hSignalLeftOfEtaChannels";
    histoLeft = histSaver->CreateDistributionHisto(name.str(),this->vecSignalLeftOfEta);
    if(histoLeft){
        histoLeft->GetXaxis()->SetTitle("Signal left of #eta");
        histoLeft->GetYaxis()->SetTitle("number of entries #");
        histoLeft->SetLineColor(kBlue);
        histSaver->SaveHistogram(histoLeft);
    }
    name.str("");name.clear();
    name<<"hSignalRightOfEtaChannels";
    histoRight = histSaver->CreateDistributionHisto(name.str(),this->vecSignalRightOfEta);
    if(histoLeft){
        histoRight->GetXaxis()->SetTitle("Signal right of #eta");
        histoRight->GetYaxis()->SetTitle("number of entries #");
        histoRight->SetLineColor(kRed);
        histSaver->SaveHistogram(histoLeft);
    }
    name.str("");name.clear();
    name<<"cSignalOfSignalsAdjacentToEta";
    if (histoLeft && histoRight)
        histSaver->SaveTwoHistos(name.str(),histoLeft,histoRight,1.);
    if(histoLeft) delete histoLeft;
    if(histoRight) delete histoRight;

}

void TTransparentAnalysis::SaveLandauVsEventNoPlots(UInt_t clusterSize){
    if (clusterSize==0)
        return;
    cout<<"SaveLandauVsEventNoPlots "<<clusterSize<<endl;
    TString name;
    TH2F* hLandau2OutOfXVsEventNo=0;
    TString section = "TimeDependence";
    if (alignMode == TSettings::transparentMode)
        section+="Trans";
    if(clusterSize-1 < vecVecPh2Highest.size()){
        name = (string)TString::Format("hLandauVsEventNo_2outOf%02d",clusterSize);
        hLandau2OutOfXVsEventNo = histSaver->CreateScatterHisto((string)name,vecVecPh2Highest.at(clusterSize-1),vectorEventNo,100,512,0,nEvents,0,3000);
        cout<<"Save "<<name<<" "<<hLandau2OutOfXVsEventNo;
        if(hLandau2OutOfXVsEventNo) cout<<" "<<hLandau2OutOfXVsEventNo->GetEntries();
        cout<<endl;
        if(vectorEventNo.size()!=vecVecPh2Highest.at(clusterSize-1).size())
            cerr<<"[TTransparentAnalysis::SaveLandauVsEventNoPlots]: Sizes of vectors are different for clusterSize "<<clusterSize<<endl;

        if (verbosity>3) cout<< name <<": "<<vectorEventNo.size()<<" "<<vecVecPh2Highest.at(clusterSize-1).size()<<endl;
        if(hLandau2OutOfXVsEventNo){
            hLandau2OutOfXVsEventNo->GetXaxis()->SetTitle("Event no.");
            hLandau2OutOfXVsEventNo->GetYaxis()->SetTitle("Pulse Height /ADC");
            histSaver->SaveHistogram(hLandau2OutOfXVsEventNo);
            TProfile* prof = histSaver->CreateAndSave1DProfileXWithFitAndInfluence(hLandau2OutOfXVsEventNo,"pol1");
            if (prof){
                TF1* fit = prof->GetFunction((TString)"fit_"+hLandau2OutOfXVsEventNo->GetName());

                Float_t value1 = 0;
                for (int i = 0; i < prof->GetNbinsX() && value1 ==0;i++)
                    value1 = prof->GetBinContent(i);

                Float_t value2 = 0;
                for (int i = prof->GetNbinsX(); i>0 && value2 ==0;i--)
                    value2 = prof->GetBinContent(i);
                TString key = TString::Format("DeltaLandauClusterSize%02d",clusterSize);
                results->setFloatValue(section,key,(value2-value1));
                key = TString::Format("LandauClusterBeginSize%02d",clusterSize);
                results->setFloatValue(section,key,(value1));
                key = TString::Format("LandauClusterEndSize%02d",clusterSize);
                results->setFloatValue(section,key,(value2));
                key = TString::Format("LandauClusterNEventsSize%02d",clusterSize);
                results->setIntValue(section,key,vectorEventNo.back() - vectorEventNo.front());
                if (fit){
                    key = TString::Format("LandauClusterFitOffsetSize%02d",clusterSize);
                    results->setFloatValue(section,key,fit->GetParameter(0));
                    key = TString::Format("LandauClusterFitSlopeSize%02d",clusterSize);
                    results->setFloatValue(section,key,fit->GetParameter(1));

                }
            }
            if (hLandau2OutOfXVsEventNo)
                delete hLandau2OutOfXVsEventNo;
        }
    }
}

void TTransparentAnalysis::saveLandausVsPositionPlots(UInt_t clusterSize){
    cout<<"saveLandausVsPositionPlots"<<endl;
    TString name;
    TH2F* htemp=0;
    if(clusterSize-1 < vecVecLandau.size()){
        name = TString::Format("hLandauVsFidCutX_ClusterSize%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecVecFidCutX,vecVecLandau[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecVecFidCutX.begin(),vecVecFidCutX.end()),
                *max_element(vecVecFidCutX.begin(),vecVecFidCutX.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position X/ch");
            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp)delete htemp;
        }
    }

    if(clusterSize-1 < vecVecPh2Highest.size() && clusterSize-1>=2){
        name = TString::Format("hLandauVsFidCutX_2OutOf%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecVecFidCutX,vecVecPh2Highest[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecVecFidCutX.begin(),vecVecFidCutX.end()),
                *max_element(vecVecFidCutX.begin(),vecVecFidCutX.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position X/ch");
            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp)delete htemp;
        }
    }


    if(clusterSize-1 < vecVecLandau.size()){
        name = TString::Format("hLandauVsPredChannel_ClusterSize%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredictedChannel,vecVecLandau[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecPredictedChannel.begin(),vecPredictedChannel.end()),
                *max_element(vecPredictedChannel.begin(),vecPredictedChannel.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted channel position (X)/ch");

            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp) delete htemp;
        }
    }


    if(clusterSize-1 < vecVecPh2Highest.size()&& clusterSize-1>=2){
        name = TString::Format("hLandauVsPredChannel_2OutOf%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredictedChannel,vecVecPh2Highest[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecPredictedChannel.begin(),vecPredictedChannel.end()),
                *max_element(vecPredictedChannel.begin(),vecPredictedChannel.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height,  2 out of %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted channel position (X)/ch");

            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp) delete htemp;
        }
    }

    if(clusterSize-1 < vecVecLandau.size()){
        name = TString::Format("hLandauVsPredDetPosY_ClusterSize%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredictedDetectorPositionY,vecVecLandau[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecPredictedChannel.begin(),vecPredictedChannel.end()),
                *max_element(vecPredictedChannel.begin(),vecPredictedChannel.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted det position (Y)/#mum");

            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp) delete htemp;
        }
    }

    if(clusterSize-1 < vecVecPh2Highest.size()){
        name = TString::Format("hLandauVsPredDetPosY_2OutOf%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredictedDetectorPositionY,vecVecPh2Highest[clusterSize-1],
                512,512,
                0,2800,
                *min_element(vecPredictedChannel.begin(),vecPredictedChannel.end()),
                *max_element(vecPredictedChannel.begin(),vecPredictedChannel.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted det position (Y)/#mum");

            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp) delete htemp;
        }
    }

    if(clusterSize-1<vecVecLandau.size() ) {
        name = TString::Format("hLandauVsFidCutY_ClusterSize%02d",clusterSize);
        Float_t miny = *min_element(vecVecFidCutY.begin(),vecVecFidCutY.end() );
        Float_t maxy = *max_element(vecVecFidCutY.begin(),vecVecFidCutY.end() );
        htemp = histSaver->CreateScatterHisto((string)name,
                vecVecFidCutY,vecVecLandau[clusterSize-1],
                512,512,
                0,2800,
                miny,maxy);
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position /ch");
            histSaver->SaveHistogram(htemp);
            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            if (htemp) delete htemp;
        }
    }

    if(clusterSize-1<vecVecPh2Highest.size() && clusterSize-1>=2) {
        name = TString::Format("hLandauVsFidCutY_2OutOf%02d",clusterSize);
        Float_t miny = *min_element(vecVecFidCutY.begin(),vecVecFidCutY.end() );
        Float_t maxy = *max_element(vecVecFidCutY.begin(),vecVecFidCutY.end() );
        htemp = histSaver->CreateScatterHisto((string)name,
                vecVecFidCutY,vecVecPh2Highest[clusterSize-1],
                512,512,
                0,2800,
                miny,maxy);
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("avrg. Silicon Hit position /ch");
            histSaver->SaveHistogram(htemp);
            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            if (htemp) delete htemp;
        }
    }

    if( clusterSize-1<vecVecLandau.size() ){
        name =TString::Format("hLandauVsPredX_ClusterSize%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredX,vecVecLandau[clusterSize-1],
                512,512,0,2800,
                *min_element(vecPredX.begin(),vecPredX.end()),
                *max_element(vecPredX.begin(),vecPredX.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted hit position X /#mum");

            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp)delete htemp;
        }
    }

    if( clusterSize-1<vecVecPh2Highest.size() && clusterSize-1>=2){
        name =TString::Format("hLandauVsPredX_2OutOf%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredX,vecVecPh2Highest[clusterSize-1],
                512,512,0,2800,
                *min_element(vecPredX.begin(),vecPredX.end()),
                *max_element(vecPredX.begin(),vecPredX.end()));
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, 2 out of %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted hit position X /#mum");
            histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
            histSaver->SaveHistogram(htemp);
            if(htemp)delete htemp;
        }
    }

    if(clusterSize-1<vecVecLandau.size()){
        name =TString::Format("hLandauVsPredY_ClusterSize%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredY,vecVecLandau[clusterSize-1],
                512,512,0,2800,
                *min_element(vecPredY.begin(),vecPredY.end()),
                *max_element(vecPredY.begin(),vecPredY.end())
        );
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, clusterSize %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted hit position Y /#mum");
        }
        histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
        histSaver->SaveHistogram(htemp);
        if(htemp)delete htemp;
    }

    if(clusterSize-1<vecVecPh2Highest.size() && clusterSize-1>=2){
        name =TString::Format("hLandauVsPredY_2OutOf%02d",clusterSize);
        htemp = histSaver->CreateScatterHisto((string)name,vecPredY,vecVecPh2Highest[clusterSize-1],
                512,512,0,2800,
                *min_element(vecPredY.begin(),vecPredY.end()),
                *max_element(vecPredY.begin(),vecPredY.end())
        );
        if(htemp){
            htemp->GetXaxis()->SetTitle(TString::Format("pulse height, 2 out of %02d",clusterSize));
            htemp->GetYaxis()->SetTitle("predicted hit position Y /#mum");
        }
        histSaver->Save1DProfileYWithFitAndInfluence(htemp,"pol1");
        histSaver->SaveHistogram(htemp);
        if(htemp)delete htemp;
    }

}

void TTransparentAnalysis::saveHistograms() {
    histSaver->SaveHistogram(hSelectedTracksAvrgSiliconHitPos);
    delete hSelectedTracksAvrgSiliconHitPos;
    string name;
    cout<<"&"<<endl;
    savePedestalHistos();
    saveNoiseHistos();
    analyseEtaDistributions();
    savePHvsEventNoAreaPlots();
    saveClusteredHistos();
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        this->saveLandausVsPositionPlots(clusterSize+1);
        this->SaveLandauVsEventNoPlots(clusterSize+1);
        histSaver->SaveHistogram(hEta[clusterSize],0);
        histSaver->SaveHistogram(hEtaCMNcorrected[clusterSize],0);
        //		if (clusterSize == 0) {
        histSaver->SaveHistogramLandau(hLandau[clusterSize]);
        histSaver->SaveHistogramLandau(hLandau2Highest[clusterSize]);
        histSaver->SaveHistogramLandau(hLandauFixedNoise[clusterSize]);
        histSaver->SaveHistogramLandau(hLandau2HighestFixedNoise[clusterSize]);
        histSaver->SaveHistogramLandau(hLandau1Highest[clusterSize]);
        histSaver->SaveHistogram(hResidualChargeWeighted[clusterSize]);
        histSaver->SaveHistogram(hResidualHighest2Centroid[clusterSize]);
        histSaver->SaveHistogram(hResidualHighestHit[clusterSize]);
        //			TCanvas *c1 = new TCanvas(TString::Format("cLandau_clusterSize%02d_both",clusterSize+1));
        //			c1->cd();
        //			hLandau[clusterSize]->Draw();
        //			fitLandau[clusterSize]->Draw("same");
        //			fitLandauFixedNoise[clusterSize]->SetLineColor(kRed);
        //			fitLandauFixedNoise[clusterSize]->Draw("same");
        //			histSaver->SaveCanvas(c1);
        ////		}
        //		else {
        //			histSaver->SaveHistogramLandau(hLaundau[clusterSize],fitLandau[clusterSize]);
        //			histSaver->SaveHistogramWithFit(hResidualChargeWeighted[clusterSize],fitResidualChargeWeighted[clusterSize]);
        //			histSaver->SaveHistogramWithFit(hResidualHighest2Centroid[clusterSize],fitResidualHighest2Centroid[clusterSize]);
        //		}
        histSaver->SaveHistogramWithFit(hResidualEtaCorrected[clusterSize],fitResidualEtaCorrected[clusterSize]);
        histSaver->SaveHistogram(hEtaIntegrals[clusterSize],0);
    }
    histSaver->SaveHistogram(hLandauMean);
    histSaver->SaveHistogram(hLandauMP);
    histSaver->SaveHistogram(hLandau2HighestMean);
    histSaver->SaveHistogram(hLandau2HighestMP);
    Float_t pw = settings->getPitchWidth(subjectDetector);
    for(UInt_t i=0;i<TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);i++){
        string name = (string)TString::Format("hRelPosVsResolutionEtaCorrectedIn%d",i+1);
        if(verbosity>6)cout<<"creating "<<name<<": "<<vecvecRelPos[i].size()<<"-"<<vecvecResXEtaCorrected[i].size()<<endl;
        TH2F* hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXEtaCorrected[i],512,512,-6000);
        hist->GetXaxis()->SetRangeUser(-pw,pw);
        hist->GetYaxis()->SetTitle("Relative predicted Position ");
        hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        histSaver->SaveHistogram(hist);
        if (hist)delete hist;

        name = (string)TString::Format("hRelChPos2VsResChargeWeighted_In_%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecResXEtaCorrected[i],512,512,-6000);
        hist->GetXaxis()->SetRangeUser(-pw,pw);
        hist->GetYaxis()->SetTitle("Relative predicted Position");
        hist->GetXaxis()->SetTitle("Residual, charge Weighted / #mum");
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        histSaver->SaveHistogram(hist);
        if (hist)delete hist;

        name = (string)TString::Format("hRelChPosVsResChargeWeighted_In_%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXChargeWeighted[i],512,512,-6000);
        hist->GetXaxis()->SetRangeUser(-pw,pw);
        hist->GetYaxis()->SetTitle("Relative predicted Position");
        hist->GetXaxis()->SetTitle("Residual, charge Weighted / #mum");
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        histSaver->SaveHistogram(hist);
        if (hist)delete hist;

        name = (string)TString::Format("hRelChPosVsResHighest2Centroid_In_%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecRelPos[i],vecvecResXHighest2Centroid[i],512,512,-6000);
        hist->GetYaxis()->SetTitle("Relative predicted Position");
        hist->GetXaxis()->SetTitle("Residual, Highest 2 Centorid / #mum");
        hist->GetXaxis()->SetRangeUser(-pw,pw);
        histSaver->SaveHistogram(hist);
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        if (hist)delete hist;

        name = (string)TString::Format("hEtaVsResolutionEtaCorrectedIn%d",i+1);
        Float_t inf = 1./0.;
        hist = histSaver->CreateScatterHisto(name,vecvecEta[i],vecvecResXEtaCorrected[i],512,512,-6000,inf,0,1);
        if(hist){
            hist->GetYaxis()->SetTitle("#eta");
            hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
            hist->GetXaxis()->SetRangeUser(-pw,pw);
            histSaver->SaveHistogram(hist);
            if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
            if ( i+1 >= TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)){
                Float_t minEta = settings->getMinimalAbsoluteEtaValue();
                TString hname = TString::Format("hResolutionEtaCorrectedIn%d_Eta_%02d_%02d",i+1,(int)(minEta*100),(int)((1-minEta)*100));
                Int_t minBin =  hist->GetYaxis()->FindBin(minEta);
                Int_t maxBin = hist->GetYaxis()->FindBin(1-minEta);
                if(verbosity) cout<<hname<<":"<<minEta<<" --> "<<minBin<<"-"<<maxBin<<" "<<flush;
                TString title = TString::Format("Resolution_{#eta-corrected} in %d channels, %.2f < #eta < %.2f",i+1,minEta,1-minEta);
                TH1F* hProj = (TH1F*)hist->ProjectionX(hname,minBin,maxBin);
                if (hProj) hProj->SetTitle(title);
                if (hProj) cout<<hProj->GetEntries()<<"/"<<	hist->GetEntries()<<endl;
                saveResolutionPlot(hProj,i,"SmallEtaRange");
                if (hProj) delete hProj;

            }
            if (hist)delete hist;
        }


        name = (string)TString::Format("hEtaCMNCorrectedVsResolutionEtaCorrectedIn%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecEtaCMNcorrected[i],vecvecResXEtaCorrected[i],512,512,-6000);
        hist->GetYaxis()->SetTitle("#eta_{CMN-corrected}");
        hist->GetXaxis()->SetTitle("Residual, Eta corrected / #mum");
        hist->GetXaxis()->SetRangeUser(-pw,pw);
        histSaver->SaveHistogram(hist);
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        if (hist)delete hist;

        name = (string)TString::Format("hRelChPosVsEta_In_%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecEta[i],512);
        hist->GetXaxis()->SetRangeUser(0,1);
        hist->GetYaxis()->SetTitle("Relative predicted Position");
        hist->GetXaxis()->SetTitle("#eta ");
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        histSaver->SaveHistogram(hist);
        if (hist)delete hist;

        name = (string)TString::Format("hRelChPosVsEtaCMN_In_%d",i+1);
        hist = histSaver->CreateScatterHisto(name,vecvecRelPos2[i],vecvecEtaCMNcorrected[i],512);
        hist->GetXaxis()->SetRangeUser(0,1);
        hist->GetYaxis()->SetTitle("Relative predicted Position");
        hist->GetXaxis()->SetTitle("#eta_{CMN corrected}");
        if(verbosity>6)cout<<hist<<" "<<hist->GetName()<<" --- > Entries:"<<hist->GetEntries()<<endl;
        histSaver->SaveHistogram(hist);
        if (hist)delete hist;
    }

    name = "hPredictedChannelPositionVsChi2";
    TH2F* hPredictedPositionVsChi2 = histSaver->CreateScatterHisto(name,vecChi2,vecPredictedPosition,2048,128,0,inf,0,20,0);
    if (hPredictedPositionVsChi2){
        hPredictedPositionVsChi2->GetXaxis()->SetTitle("Predicted Channel Position");
        hPredictedPositionVsChi2->GetYaxis()->SetTitle("Max. #chi^{2}_{X,Y}");
        histSaver->SaveHistogram(hPredictedPositionVsChi2,false);
        delete hPredictedPositionVsChi2;
    }
    name = "hRelativePredictedChannelPositionVsChi2";
    hPredictedPositionVsChi2 = histSaver->CreateScatterHisto(name,vecChi2,vecRelPredictedPosition,512,128,0,1,0,20,0);
    if (hPredictedPositionVsChi2){
        hPredictedPositionVsChi2->GetXaxis()->SetTitle("relative Predicted Channel Position");
        hPredictedPositionVsChi2->GetYaxis()->SetTitle("Max. #chi^{2}_{X,Y}");
        histSaver->SaveHistogram(hPredictedPositionVsChi2,false);
        delete hPredictedPositionVsChi2;
    }
    name = "hPredictedChannelPosition";
    TH1F* hPredictedPosition = histSaver->CreateDistributionHisto(name,vecPredictedPosition,2048,histSaver->maxWidth,0,inf);
    if (hPredictedPosition){
        hPredictedPosition->GetXaxis()->SetTitle("Predicted Channel Position");
        histSaver->SaveHistogram(hPredictedPosition);
        delete hPredictedPosition;
    }
    name = "hRelativePredictedChannelPosition";
    hPredictedPosition = histSaver->CreateDistributionHisto(name,vecRelPredictedPosition,512,histSaver->maxWidth,0,1);
    if (hPredictedPosition){
        hPredictedPosition->GetXaxis()->SetTitle("relative Predicted Channel Position");
        histSaver->SaveHistogram(hPredictedPosition);
        delete hPredictedPosition;
    }
}

void TTransparentAnalysis::deleteHistograms() {
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        if(hLandau[clusterSize]) delete hLandau[clusterSize];
        if(hLandau2Highest[clusterSize])delete hLandau2Highest[clusterSize];
        if(hLandau1Highest[clusterSize])delete hLandau1Highest[clusterSize];
        if ( hEta[clusterSize]) delete hEta[clusterSize];
        if ( hEtaCMNcorrected[clusterSize]) delete hEtaCMNcorrected[clusterSize];
        if (hResidualChargeWeighted[clusterSize]) delete hResidualChargeWeighted[clusterSize];
        if (hResidualHighest2Centroid[clusterSize]) delete hResidualHighest2Centroid[clusterSize];
        if (hResidualHighestHit[clusterSize]) delete hResidualHighestHit[clusterSize];
    }
    delete hLandauMean;
    delete hLandauMP;
    delete hLandau2HighestMean;
    delete hLandau2HighestMP;
}

void TTransparentAnalysis::deleteFits() {
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        delete fitLandau[clusterSize];
        delete fitLandau2Highest[clusterSize];
        delete fitResidualChargeWeighted[clusterSize];
        delete fitResidualHighest2Centroid[clusterSize];
    }
}

void TTransparentAnalysis::printCutFlow() {
    cout << "\n\n\n";
    cout << "TTransparentAnalysis Cutflow" << endl;
    cout << "number of events\t " << setw(8) << nEvents << endl;
    cout << "no valid silicon track\t-" << setw(8) << noValidTrack << endl;
    cout << "not in fid cut region \t-" << setw(8) << noFidCutRegion << endl;
    cout << "used for alignment    \t-" << setw(8) << usedForAlignment << endl;
    cout << "too high Chi2 value   \t-" << setw(8) << highChi2 <<endl;
    //	cout << "used for si alignment\t-" << setw(8) << usedForSiliconAlignment << endl;
    cout << "region not on plane\t-" << setw(8) << regionNotOnPlane << endl;
    cout << "screened channel\t-" << setw(8) << screenedChannel << endl;
    cout << "saturated channel\t-" << setw(8) << saturatedChannel << endl;
    cout << "\t\t\t---------" << endl;
    cout << "total analyzed events\t " << setw(8) << nAnalyzedEvents << endl;
    TString section = "CutFlowTransparentAnalysis";
    if(alignMode == TSettings::transparentMode)
        section+="Trans";
    results->setIntValue(section,"nEvents",nEvents);
    results->setIntValue(section,"noValidTrack",noValidTrack);
    results->setIntValue(section,"noFidCutRegion",noFidCutRegion);
    results->setIntValue(section,"usedForAlignment",usedForAlignment);
    results->setIntValue(section,"highChi2",highChi2);
    results->setIntValue(section,"regionNotOnPlane",regionNotOnPlane);
    results->setIntValue(section,"screenedChannel",screenedChannel);
    results->setIntValue(section,"saturatedChannel",saturatedChannel);
    results->setIntValue(section,"nAnalyzedEvents",nAnalyzedEvents);
}

void TTransparentAnalysis::printEvent() {
    cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << endl;
    if (eventReader->useForAnalysis() == 0) {
        cout << "this track is not used for the analysis.." << endl;
        return;
    }
    cout << "predicted pos in lab system:\t" << this->predPosition << "\tpredicted perp position:\t" << this->predPerpPosition << endl;
    cout << "predicted pos in det system:\t" << this->positionInDetSystemMetric << endl;
    cout << "clustered analysis position in lab system:\t" << eventReader->getStripXPosition(subjectPlane,this->predPerpPosition,clusterCalcMode) << endl;
    cout << "clustered analysis position in det system:\t" << eventReader->getMeasuredPositionMetricSpace(subjectDetectorCoordinate, subjectPlane, clusterCalcMode) << endl;
    if (this->checkPredictedRegion(subjectDetector, this->positionInDetSystemMetric, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)) == false) {
        cout << "this track did not pass the check.." << endl;
        return;
    }
    for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
        cout << "transparent cluster of size " << clusterSize+1 << ":" << endl;
        cout << "\tpulse height:\t" << this->transparentClusters.getCharge(cmCorrected) << endl;
        cout << "\teta:\t" << this->transparentClusters.getEta() << endl;
        cout << "\tresidual:\t" << this->getResidual(this->transparentClusters,cmCorrected,this->clusterCalcMode,hEtaIntegrals[clusterSize]) << endl;
        cout << "\tcluster pos in det system:\t" << this->transparentClusters.getPosition(this->clusterCalcMode) << endl;
        cout << "\tcluster pos in lab system:\t" << eventReader->getPositionOfCluster(subjectDetector, this->transparentClusters, this->predPerpPosition, this->clusterCalcMode) << endl;
    }
    return;
}

Float_t TTransparentAnalysis::getResidual(TCluster cluster,bool cmnCorrected, TCluster::calculationMode_t clusterCalculationMode, TH1F* hEtaInt) {
    if(clusterCalculationMode == TCluster::corEta && hEtaInt==0)
        cout<<"getResidual::EtaInt==0"<<endl;
    return eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,cmnCorrected,clusterCalculationMode, hEtaInt)-this->predPosition;
}

void TTransparentAnalysis::printCluster(TCluster cluster) {
    cout << "\n--- event " << nEvent;
    cout << "\n\tcluster size: " << cluster.getClusterSize();
    cout << "\n\tcharge: " << cluster.getCharge(cmCorrected);
    cout << "\n\tcharge of 2 highest centroid: " << cluster.getCharge((UInt_t)2,cmCorrected);
    cout << "\n\thighest channel: " << cluster.getHighestSignalChannel();
    cout << "\n\thighest 2 centroid: " << cluster.getHighest2Centroid(cmCorrected);
    cout << "\n\tcluster position of highest channel: " << cluster.getClusterPosition(cluster.getHighestSignalChannel());
    cout << "\n\thighest channel is seed? " << cluster.isSeed(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
    cout << "\n\thighest channel is hit? " << cluster.isHit(cluster.getClusterPosition(cluster.getHighestSignalChannel()));
    cout << "\n\tseed sigma: " << cluster.getSeedSigma();
    cout << "\n\thit sigma: " << cluster.getHitSigma();
    cout << "\n\tpredicted channel: " << positionInDetSystemMetric;
    cout << "\n\tpredicted position: " << predPosition;
    cout << "\n\tcharge weighted position (TCluster::chargeWeighted): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::chargeWeighted);
    cout << "\n\thighest 2 centroid position (TCluster::highest2Centroid): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,TCluster::highest2Centroid);
    cout << "\n\tcharge weighted residual: " << getResidual(cluster,cmCorrected,TCluster::chargeWeighted);
    cout << "\n\thighest 2 centroid residual: " << getResidual(cluster,cmCorrected,TCluster::highest2Centroid);
    cout << "\n\teta: " << cluster.getEta();
    if (hEtaIntegrals.size() != 0) {
        cout << "\n\teta corrected position (TCluster::corEta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,cmCorrected,TCluster::corEta,hEtaIntegrals[cluster.getClusterSize()]);
        cout << "\n\teta corrected residual (TCluster::corEta): " << getResidual(cluster,cmCorrected,TCluster::corEta,hEtaIntegrals[cluster.getClusterSize()-1]);
        cout << "\n\teta corrected position (TCluster::eta): " << eventReader->getPositionOfCluster(subjectDetector,cluster,this->predPerpPosition,cmCorrected,TCluster::eta,hEtaIntegrals[cluster.getClusterSize()]);
        //		cout << "\n\teta corrected residual (TCluster::eta): " << getResidual(cluster,TCluster::eta);
    }
    cout << "\n\t";
    cluster.Print();
}

/** returns difference between cluster position and calculated position for a given cluster
 * @param cluster
 * @param clusterCalculationMode
 * @param hEtaInt
 * @author Lukas Baeni
 * @return
 */
void TTransparentAnalysis::createEventVector(Int_t startEvent) {

    nAnalyzedEvents = 0;
    regionNotOnPlane = 0;
    saturatedChannel = 0;
    screenedChannel = 0;
    noValidTrack = 0;
    noFidCutRegion = 0;
    usedForAlignment = 0;
    highChi2 =0;
    vecTransparentClusters.clear();
    eventNumbers.clear();
    vecEvents.clear();
    settings->getSelectionFidCuts()->Print(1);
    settings->getSelectionFidCuts()->setRunDescription(settings->getRunDescription());
    cout<<"Rundes: "<<settings->getRunDescription()<<"\tIndex: "<<settings->getSelectionFidCuts()->getActiveIndex()<<endl;
    cout<<"Creating  Event Vector "<<endl;
    for (nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
        TRawEventSaver::showStatusBar(nEvent,nEvents+startEvent,100);
        //		if (verbosity > 4) cout << "-----------------------------\n" << "analyzing event " << nEvent << ".." << eventReader<<endl;
        if (settings->useForAlignment(nEvent,nEvents)){
            usedForAlignment++;
            continue;
        }
        eventReader->LoadEvent(nEvent);
        if (eventReader->isValidTrack() == 0) {
            //		if (eventReader->useForAnalysis() == 0) {
            if (verbosity > 6) printEvent();
            noValidTrack++;
            continue;
        }
        Float_t fiducialValueX = eventReader->getFiducialValueX();
        Float_t fiducialValueY = eventReader->getFiducialValueY();
        Int_t region = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
        cout<<"";
        if (!settings->getSelectionFidCuts()->IsInFiducialCut(fiducialValueX,fiducialValueY)) {
            noFidCutRegion++;
            continue;
        }
        transparentClusters.clear();
        if (!this->predictPositions(true)){
            if (verbosity>4) cout<< nEvent << ": Chi2 to high: "<< positionPrediction->getChi2()<<endl;
            highChi2++;
            continue;
        }
        //		cout<<"predRegion("<<nEvent<<");"<<endl;

        //		cout<<"add Event"<<nEvent<<endl;
        bool predRegion = this->checkPredictedRegion(subjectDetector, this->positionInDetSystemChannelSpace, TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
        if (predRegion == false)
            continue;

        //		cout<<"transparentClusters("<<nEvent<<");"<<endl;
        Int_t maxClusSize = TPlaneProperties::getMaxTransparentClusterSize(subjectDetector);
        transparentClusters = this->makeTransparentCluster(eventReader, settings,subjectDetector, positionInDetSystemChannelSpace, maxClusSize);

        Float_t pos = positionInDetSystemChannelSpace;
        Float_t channels = 15;
        Int_t direction = 2*(int)(gRandom->Uniform()>.5)-1;
        pos +=direction * channels;
        bool isMasked = settings->IsMasked(subjectDetector,pos);
        if(isMasked){
            //			cout<<"masked: "<<pos<<endl;
            pos -= 2*direction*channels;
            isMasked = settings->IsMasked(subjectDetector,pos);
        }
        //		cout<<"nonHit("<<nEvent<<");"<<endl;
        if(!isMasked){
            //			cout<<"not masked("<<nEvent<<");"<<endl;
            TCluster noHitCluster = makeTransparentCluster(eventReader,settings,subjectDetector,pos,TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
            //			cout<<"123"<<endl;
            Float_t charge =noHitCluster.getCharge(cmCorrected,true);
            //			cout<<"\n"<<nEvent<<" "<<noHitCluster.getCharge(true,true)<<"/"<<noHitCluster.getCharge(false,true);
            if(charge==0){
                //				cout<<"$#\n";
                noHitCluster.Print(11);
            }

            noHitClusters.push_back(noHitCluster);
            //			cout<<nEvent<<"Non Hit Cluster @ "<<pos <<" from "<<positionInDetSystemChannelSpace<<" "<<direction<<" "<<pos-positionInDetSystemChannelSpace<<endl;
        }
        else{
            //			cout<<nEvent<<"isMasked: "<< pos <<" from "<<positionInDetSystemChannelSpace<<" "<<direction<<" "<<pos-positionInDetSystemChannelSpace<<endl;
        }
        //		cout<<"analyse("<<nEvent<<");"<<endl;
        nAnalyzedEvents++;
        this->fillHistograms();
        if (verbosity > 4) printEvent();
        //		cout<<"push Back("<<nEvent<<");"<<endl;

        // save clusters for eta corrected analysis
        vecTransparentClusters.push_back(transparentClusters);
        eventNumbers.push_back(nEvent);
        vecEvents.push_back(eventReader->getEvent());

    }
}

TCluster TTransparentAnalysis::makeTransparentCluster(TTracking *reader,TSettings* set, UInt_t det, Float_t centerPosition, UInt_t clusterSize) {
    // get channel and direction for clustering
    //	cout<<"makeTransparentCluster"<<endl;
    if (reader==0){
        cerr<<" TTransparentAnalysis::makeTransparentCluster TTracking == 0 "<<endl;
        return TCluster();
    }
    if (set == 0 ){
        cerr<<" TTransparentAnalysis::makeTransparentCluster TSettings == 0 "<<endl;
        return TCluster();
    }
    //	cout<<"[TTransparentAnalysis::makeTransparentCluster]\t";
    UInt_t centerChannel;
    int direction;
    direction = getSignedChannelNumber(centerPosition);
    //	cout << "centerPosition: " << centerPosition << "\tdirection: " << direction << endl;
    centerChannel = TMath::Abs(direction);
    if (direction < 0) direction = -1;
    else direction = 1;
    Float_t cmNoise = reader->getCMNoise();

    // make cluster
    TCluster transparentCluster = TCluster(reader->getEvent_number(), det, -99, -99, TPlaneProperties::getNChannels(det),cmNoise);
    UInt_t currentChannel = centerChannel;
    for (UInt_t iChannel = 0; iChannel < clusterSize; iChannel++) {
        //		if( currentChannel < 0 || currentChannel >= TPlaneProperties::getNChannelsDiamond() )
        //			cout<<"\n"<<reader->getEvent_number()<<": Cannot create channel with: det"<<det<<", centerPoisition: "<<centerPosition<< ", direction: "<<direction<<", centerChannel: "<<centerChannel<<" "<<iChannel<<flush;
        direction *= -1;
        currentChannel += direction * iChannel;
        Int_t adcValue=reader->getAdcValue(det,currentChannel);
        Float_t pedMean = reader->getPedestalMean(det,currentChannel,false);
        Float_t pedMeanCMN = reader->getPedestalMean(det,currentChannel,true);
        Float_t pedSigma = reader->getPedestalSigma(det,currentChannel,false);
        Float_t pedSigmaCMN = reader->getPedestalSigma(det,currentChannel,true);
        bool isScreened = set->isDet_channel_screened(det,currentChannel);
        if (TPlaneProperties::IsValidChannel(det,currentChannel))
            transparentCluster.addChannel(currentChannel,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened);
        //		else
        //			cout<<"\t cannot add invalid channel"<<currentChannel<<" @ "<<reader->getEvent_number()<<endl;
        //		transparentCluster.addChannel(currentChannel, reader->getRawSignal(det,currentChannel), reader->getRawSignalInSigma(det,currentChannel), reader->getAdcValue(det,currentChannel), reader->isSaturated(det,currentChannel), settings->isDet_channel_screened(det,currentChannel));
    }
    //	cout<<"[done]"<<endl;
    transparentCluster.UpdateHighestSignalChannel();
    transparentCluster.SetTransparentCluster(centerPosition);
    transparentCluster.SetTransparentClusterSize(clusterSize);
    return transparentCluster;
}

void TTransparentAnalysis::clearEventVector() {
    while ( this->vecEvents.size() != 0)  {
        TEvent* event = vecEvents.back();
        if (event) delete event;
        vecEvents.pop_back();
    }

}

void TTransparentAnalysis::analyseNonHitEvents() {
    cout<<"[TTransparentAnalysis::analyseNonHitEvents] "<<noHitClusters.size()<<endl;
    vector <TH1F*> hNonHitNoiseDistributions;
    vector <TH1F*> hNonHitNoiseDistributionsCMN;
    vector <TH1F*> hNonHitNoiseDistributions2OutOfX;
    vector <TH1F*> hNonHitNoiseDistributions2OutOfXCMN;
    for (UInt_t i = 0; i < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); i++){
        TString name = TString::Format("hNonHitPulseHeightDitribution_ClusterSize%02d",i+1);
        TH1F *histo = new TH1F(name,name,1000,-499.5,499.5);
        histo->GetXaxis()->SetTitle("PH_{trans Clus - non Hit} / ADC");
        histo->GetYaxis()->SetTitle("number of entries #");
        hNonHitNoiseDistributions.push_back(histo);

        name = TString::Format("hNonHitPulseHeightDitributionCMN_ClusterSize%02d",i+1);
        histo = new TH1F(name,name,1000,-499.5,499.5);
        histo->GetXaxis()->SetTitle("PH_{trans Clus - non Hit} - cm corrected / ADC");
        histo->GetYaxis()->SetTitle("number of entries #");
        hNonHitNoiseDistributionsCMN.push_back(histo);

        name = TString::Format("hNonHitPulseHeightDitribution2OutOf%02d",i+1);
        histo = new TH1F(name,name,1000,-499.5,499.5);
        histo->GetXaxis()->SetTitle(TString::Format("PH_{trans Clus - non Hit - 2 out of %d }  / ADC",i+1));
        histo->GetYaxis()->SetTitle("number of entries #");
        hNonHitNoiseDistributions2OutOfX.push_back(histo);

        name = TString::Format("hNonHitPulseHeightDitribution2OutOf%02d",i+1);
        histo = new TH1F(name,name,1000,-499.5,499.5);
        histo->GetXaxis()->SetTitle(TString::Format("PH_{trans Clus - non Hit - 2 out of %d } - cm corrected / ADC",i+1));
        histo->GetYaxis()->SetTitle("number of entries #");
        hNonHitNoiseDistributions2OutOfXCMN.push_back(histo);
    }

    for (UInt_t i = 0; i< noHitClusters.size(); i++){
        for (UInt_t j = 0; j < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); j++){
            noHitClusters[i].SetTransparentClusterSize(j+1);
            Float_t chargeCMN = noHitClusters[i].getCharge(true,true);
            Float_t charge = noHitClusters[i].getCharge(false,true);
            hNonHitNoiseDistributions[j]->Fill(charge);
            hNonHitNoiseDistributionsCMN[j]->Fill(chargeCMN);
            hNonHitNoiseDistributions2OutOfX[j]->Fill(noHitClusters[i].getCharge(2,false,true));
            hNonHitNoiseDistributions2OutOfXCMN[j]->Fill(noHitClusters[i].getCharge(2,true,true));
        }
    }
    noiseWidths.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    noiseWidthsCMN.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));

    noiseWidths2OutOfX.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));
    noiseWidths2OutOfXCMN.resize(TPlaneProperties::getMaxTransparentClusterSize(subjectDetector));

    for (UInt_t j = 0; j < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); j++){
        for(UInt_t k = 0; k<4;k++){

            TH1F* histo;
            switch (k) {
                case 0: histo = hNonHitNoiseDistributions[j];break;
                case 1: histo = hNonHitNoiseDistributionsCMN[j];break;
                case 2: histo = hNonHitNoiseDistributions2OutOfX[j];break;
                case 3: histo = hNonHitNoiseDistributions2OutOfXCMN[j];break;
            }
            TString name;
            name = "fitGaus_";
            name.Append(histo->GetName());
            TF1* fit = new TF1(name,"gaus",-500,500);
            fit->SetLineColor(kBlue);
            histo->Fit(fit,"Q","",histo->GetMean()-2*histo->GetRMS(),histo->GetMean()+2*histo->GetRMS());
            histo->Draw("goff");
            Float_t xmin = fit->GetParameter(1)-4*fit->GetParameter(2);
            Float_t xmax = fit->GetParameter(1)+4*fit->GetParameter(2);
            histo->GetXaxis()->SetRangeUser(xmin,xmax);
            histSaver->SaveHistogram(histo);
            switch(k){
                case 0: noiseWidths[j] = fit->GetParameter(2);break;
                case 1: noiseWidthsCMN[j] = fit->GetParameter(2);break;
                case 2: noiseWidths2OutOfX[j] = fit->GetParameter(2);break;
                case 3: noiseWidths2OutOfXCMN[j] = fit->GetParameter(2);break;
            }
        }
        cout<<"ClusterSize "<<j<<": "<< noiseWidths[j]<<"/"<<noiseWidthsCMN[j]<<" "<<noiseWidths2OutOfX[j]<<"/"<<noiseWidths2OutOfXCMN[j]<<endl;
        TString name = TString::Format("cNonHitPulseHeightDitribution_ClusterSize%02d",j+1);
        histSaver->SaveTwoHistos((string)name,hNonHitNoiseDistributions[j],hNonHitNoiseDistributionsCMN[j]);
        name = TString::Format("cNonHitPulseHeightDitribution_2OutOf%02d",j+1);
        histSaver->SaveTwoHistos((string)name,hNonHitNoiseDistributions2OutOfX[j],hNonHitNoiseDistributions2OutOfXCMN[j]);
        delete hNonHitNoiseDistributions[j];
        delete hNonHitNoiseDistributionsCMN[j];
        delete hNonHitNoiseDistributions2OutOfX[j];
        delete hNonHitNoiseDistributions2OutOfXCMN[j];
    }
}

void TTransparentAnalysis::initPedestalAndNoiseHistos(UInt_t maxEvents) {
    cout<<"initPedestalHistos"<<flush;
    UInt_t start = settings->getAlignmentEvents(maxEvents);
    UInt_t nBins = (maxEvents-start)/1000;
    for(UInt_t ch = 0; ch< TPlaneProperties::getNChannelsDiamond();ch++){
        if(settings->IsMasked(subjectDetector,ch))
            continue;
        TString name = TString::Format("hPedestalVsEventNo_det_%d_ch_%03d",subjectDetector,ch);
        TProfile* prof = new TProfile(name,name,nBins,start,maxEvents);
        prof->GetXaxis()->SetTitle("EventNo");
        TString title = TString::Format("pedestal_{ch %03d} /ADC",ch);
        if(settings->doCommonModeNoiseCorrection()) title.Append(" CM corrected");
        prof->GetYaxis()->SetTitle(title);
        hPedestalVsEvenNo[ch] = prof;

        //Noise of each channel
        name = TString::Format("hNoiseVsEventNo_det_%d_ch_%03d",subjectDetector,ch);
        prof = new TProfile(name,name,nBins,start,maxEvents);
        prof->GetXaxis()->SetTitle("EventNo");
        title = TString::Format("noise_{ch %03d} /ADC",ch);
        if(settings->doCommonModeNoiseCorrection()) title.Append(" CM corrected");
        prof->GetYaxis()->SetTitle(title);
        hNoiseVsEvenNo[ch] = prof;
    }
    TString name = "hComonModeNoiseVsEventNo";
    hCmnVsEventNo = new TProfile(name,name,nBins,start,maxEvents);
    hCmnVsEventNo->GetXaxis()->SetTitle("EventNo");
    hCmnVsEventNo->GetYaxis()->SetTitle("common mode noise /ADC");
    cout<<"."<<endl;
}

void TTransparentAnalysis::fillPedestalsAndNoiseHistos() {
    std::map<UInt_t,TProfile*>::iterator it;
    for(it=hPedestalVsEvenNo.begin(); it!=hPedestalVsEvenNo.end(); it++){
        UInt_t channel = (*it).first;
        Float_t pedestal = eventReader->getPedestalMean(subjectDetector,channel,settings->doCommonModeNoiseCorrection());
        (*it).second->Fill(nEvent,pedestal);
    }
    for(it=hNoiseVsEvenNo.begin(); it!=hNoiseVsEvenNo.end(); it++){
        UInt_t channel = (*it).first;
        Float_t noise = eventReader->getPedestalSigma(subjectDetector,channel,settings->doCommonModeNoiseCorrection());
        (*it).second->Fill(nEvent,noise);
    }
    hCmnVsEventNo->Fill(nEvent,eventReader->getCMNoise());
}

void TTransparentAnalysis::saveNoiseHistos() {
    vector<Float_t> vecCh;
    vector<Float_t> vecSlope;
    TAnalysisOfClustering::saveVariableVsEventNoPlots(settings,histSaver,hNoiseVsEvenNo,"Noise",&vecSlope,&vecCh);
    return;
    THStack* stack = new THStack("hNoisesVsEventNo","Noises vs Event No");
    TH1F* hNoiseSlopesVsChannel = new TH1F("hNoiseSlopesVsChannel","slope of hNoiseVsEventNo for each ch",128,0,128);
    hNoiseSlopesVsChannel->GetXaxis()->SetTitle("channel no");
    hNoiseSlopesVsChannel->GetYaxis()->SetTitle("slope m = ADC/Event");
    UInt_t color = 0;
    std::map<UInt_t,TProfile*>::iterator it;
    TF1* pol1 = new TF1("pol1_fit","pol1",0,5e6);
    pol1->SetLineColor(kBlue);
    pol1->SetLineWidth(1);
    Double_t minStack = 1e9;
    Double_t maxStack = -1e9;
    for(it=hNoiseVsEvenNo.begin(); it!=hNoiseVsEvenNo.end(); it++){
        TProfile* prof = (*it).second;
        if(!prof) continue;
        TF1* fit = (TF1*)pol1->Clone(prof->GetName()+(TString)"_fit");
        if((*it).first%5==0){
            TProfile* prof2 =(TProfile*)prof->Clone();
            prof2->SetTitle(TString::Format("Channel %3d",(*it).first));
            prof2->SetLineColor(color);
            prof2->SetMarkerColor(color);
            stack->Add(prof2);
            minStack = TMath::Min( minStack, prof->GetBinContent(prof->GetMinimumBin()));
            maxStack = TMath::Max( maxStack, prof->GetBinContent(prof->GetMaximumBin()));
            color++;
        }
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        hNoiseSlopesVsChannel->SetBinContent((hNoiseSlopesVsChannel->FindBin((*it).first)),fit->GetParameter(1));
        vecCh.push_back((*it).first);
        vecSlope.push_back(fit->GetParameter(1));
        delete prof;
        (*it).second= 0;
        hNoiseVsEvenNo.erase(it);
    }
    TGraph graph = histSaver->CreateDipendencyGraph("gNoiseSlopeVsChannel",vecSlope,vecCh);
    graph.Draw("AP");
    graph.GetXaxis()->SetTitle("channel");
    graph.GetYaxis()->SetTitle("Noise slope for channel");
    histSaver->SaveGraph(&graph,"gNoiseSlopeVsChannel","ABP");

    TH1F* hSlopes = histSaver->CreateDistributionHisto("hNoiseSlopes",vecSlope,10);
    hSlopes->GetXaxis()->SetTitle("Noise slope ADC/Event");
    hSlopes->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hSlopes);
    delete hSlopes;
    if(hCmnVsEventNo) {
        cout<<"save "<<hCmnVsEventNo->GetName()<<endl;
        TF1* fit = (TF1*)pol1->Clone(hCmnVsEventNo->GetName()+(TString)"_fit");
        histSaver->Save1DProfileWithFitAndInfluence(hCmnVsEventNo,fit,true);
        delete hCmnVsEventNo;
        hCmnVsEventNo=0;
    }
    if(color!=0){
        cout<<"save stack "<<minStack<<"-"<<maxStack<<endl;
        stack->Draw("goff");
        stack->SetObjectStat(false);
        if(stack->GetXaxis()){
            stack->GetXaxis()->SetTitle("Event No");
            cout<<"Xaxis: "<<stack->GetXaxis()->GetTitle()<<endl;
        }
        if(stack->GetYaxis()){
            stack->GetYaxis()->SetTitle("Noise /ADC");
            cout<<"Set range"<<endl;
            stack->GetYaxis()->SetRangeUser(minStack*.98,maxStack*1.05);
            cout<<"Yaxis: "<<stack->GetYaxis()->GetTitle()<<endl;
        }
        stack->SetMinimum(minStack*.98);
        stack->SetMaximum(maxStack*1.05);
        stack->SetObjectStat(false);
    }
    histSaver->SaveStack(stack,"nostack",true);
    cout<<"hNoiseSlopesVsChannel:MIN: "<<hNoiseSlopesVsChannel->GetBinContent(hNoiseSlopesVsChannel->GetMinimumBin())<<endl;
    hNoiseSlopesVsChannel->SetMinimum(hNoiseSlopesVsChannel->GetBinContent(hNoiseSlopesVsChannel->GetMinimumBin()));
    histSaver->SaveHistogram(hNoiseSlopesVsChannel,false,false);
    delete hNoiseSlopesVsChannel;
    delete stack;
}


void TTransparentAnalysis::saveClusteredHistos(){
    if(hLandauVsEventNo_Clustered){
        TString name = "hLandau_Clustered";
        TH1F* hLandau_Clustered = (TH1F*)hLandauVsEventNo_Clustered->ProjectionY(name);
        hLandau_Clustered->SetTitle(name);
        TF1* fit = landauGauss->doLandauGaussFit(hLandau_Clustered,true);Float_t mean;
        if(hLandau_Clustered)
            mean = hLandau_Clustered->GetMean();
        else mean = -1;
        Float_t mp,width,gSigma;
        if(fit){
            cout<<"READ FIT Clustered"<<endl;
            mp = fit->GetParameter(1);
            width = fit->GetParameter(0);
            gSigma = fit->GetParameter(3);
        }
        else {
            cout<<"CANNOT READ FIT Clustered"<<endl;
            mp = -1;
            width = -1;
            gSigma = -1;
        }
        if(results){
            cout<<"SAVE PH clustered"<<mean<<" "<<mp<<" "<<width<<" "<<gSigma<<endl;
            results->setPH_clustered(mean,mp,width,gSigma,alignMode);
        }
        else
            cout<<" CANNOT SAVE PH clustered - results: "<<results<<endl;
        //char t;
        //cin >>t;
        histSaver->SaveHistogramLandau(hLandau_Clustered);
        delete hLandau_Clustered;
        histSaver->SaveHistogram(hLandauVsEventNo_Clustered);
        histSaver->CreateAndSave1DProfileXWithFitAndInfluence(hLandauVsEventNo_Clustered,"pol1");

        delete hLandauVsEventNo_Clustered;
        hLandauVsEventNo_Clustered =0;
    }
    if(hClusterSize_Clustered){
        histSaver->SaveHistogram(hClusterSize_Clustered);
        delete hClusterSize_Clustered;
        hClusterSize_Clustered =0;
    }
    if(hLandauVsClusterSize_Clustered){
        for (UInt_t i = 0; i < hLandauVsClusterSize_Clustered->GetNbinsX();i++){
            TString name = hLandauVsClusterSize_Clustered->GetName();
            name.Append(TString::Format("_clusterSize_%d",i));
            TH1F* hProjection = (TH1F*)hLandauVsClusterSize_Clustered->ProjectionX(name,i,i);
            histSaver->SaveHistogramLandau(hProjection);
            delete hProjection;
        }
        histSaver->SaveHistogram(hLandauVsClusterSize_Clustered);
        histSaver->CreateAndSave1DProfileXWithFitAndInfluence(hLandauVsClusterSize_Clustered,"");
        histSaver->Save1DProfileYWithFitAndInfluence(hLandauVsClusterSize_Clustered,"");
        delete hLandauVsClusterSize_Clustered;
        hLandauVsClusterSize_Clustered= 0;
    }
    if(hNClusteres_Clustered){
        histSaver->SaveHistogram(hNClusteres_Clustered);
        delete hNClusteres_Clustered;
        hNClusteres_Clustered = 0;
    }
}

void TTransparentAnalysis::savePedestalHistos() {

    vector<Float_t> vecCh;
    vector<Float_t> vecPed;
    TAnalysisOfClustering::saveVariableVsEventNoPlots(settings,histSaver,hPedestalVsEvenNo,"Pedestal",&vecPed,&vecCh);
    return;
    THStack* stack = new THStack("hPedestalsVsEventNo","pedestals vs Event No");
    TH1F* hPedestalSlopesVsChannel = new TH1F("hPedestalSlopesVsChannel","slope of hPedestalVsEventNo for each ch",128,0,128);
    hPedestalSlopesVsChannel->GetXaxis()->SetTitle("channel no");
    hPedestalSlopesVsChannel->GetYaxis()->SetTitle("slope m = ADC/Event");
    UInt_t color = 0;
    std::map<UInt_t,TProfile*>::iterator it;
    TF1* pol1 = new TF1("pol1_fit","pol1",0,5e6);
    pol1->SetLineColor(kBlue);
    pol1->SetLineWidth(1);
    Double_t minStack = 1e9;
    Double_t maxStack = -1e9;
    for(it=hPedestalVsEvenNo.begin(); it!=hPedestalVsEvenNo.end(); it++){
        TProfile* prof = (*it).second;
        if(!prof) continue;
        TF1* fit = (TF1*)pol1->Clone(prof->GetName()+(TString)"_fit");
        if((*it).first%5==0){
            TProfile* prof2 =(TProfile*)prof->Clone();
            prof2->SetTitle(TString::Format("Channel %3d",(*it).first));
            prof2->SetLineColor(color);
            prof2->SetMarkerColor(color);
            stack->Add(prof2);
            minStack = TMath::Min( minStack, prof->GetBinContent(prof->GetMinimumBin()));
            maxStack = TMath::Max( maxStack, prof->GetBinContent(prof->GetMaximumBin()));
            color++;
        }
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        hPedestalSlopesVsChannel->SetBinContent((hPedestalSlopesVsChannel->FindBin((*it).first)),fit->GetParameter(1));
        vecCh.push_back((*it).first);
        vecPed.push_back(fit->GetParameter(1));
        delete prof;
        (*it).second= 0;
        hPedestalVsEvenNo.erase(it);
    }
    TGraph graph = histSaver->CreateDipendencyGraph("gPedestalSlopeVsChannel",vecPed,vecCh);
    graph.Draw("AP");
    graph.GetXaxis()->SetTitle("channel");
    graph.GetYaxis()->SetTitle("pedestal slope for channel");
    histSaver->SaveGraph(&graph,"gPedestalSlopeVsChannel","ABP");

    TH1F* hSlopes = histSaver->CreateDistributionHisto("hPedestalSlopes",vecPed,10);
    hSlopes->GetXaxis()->SetTitle("pedestal slope ADC/Event");
    hSlopes->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hSlopes);
    delete hSlopes;
    if(color!=0){
        cout<<"save stack "<<minStack<<"-"<<maxStack<<endl;
        stack->Draw("goff");
        stack->SetObjectStat(false);
        if(stack->GetXaxis()){
            stack->GetXaxis()->SetTitle("Event No");
            cout<<"Xaxis: "<<stack->GetXaxis()->GetTitle()<<endl;
        }
        if(stack->GetYaxis()){
            stack->GetYaxis()->SetTitle("Pedestal /ADC");
            cout<<"Set range"<<endl;
            stack->GetYaxis()->SetRangeUser(minStack*.98,maxStack*1.05);
            cout<<"Yaxis: "<<stack->GetYaxis()->GetTitle()<<endl;
        }
        stack->SetMinimum(minStack*.98);
        stack->SetMaximum(maxStack*1.05);
        stack->SetObjectStat(false);
    }
    histSaver->SaveStack(stack,"nostack",true);
    cout<<"hPedestalSlopesVsChannel:MIN: "<<hPedestalSlopesVsChannel->GetBinContent(hPedestalSlopesVsChannel->GetMinimumBin())<<endl;
    hPedestalSlopesVsChannel->SetMinimum(hPedestalSlopesVsChannel->GetBinContent(hPedestalSlopesVsChannel->GetMinimumBin()));
    histSaver->SaveHistogram(hPedestalSlopesVsChannel,false,false);
    delete hPedestalSlopesVsChannel;
    delete stack;
}
std::pair<Float_t,Float_t >  TTransparentAnalysis::getFWCrossingPoint(TH1F* hRes,Float_t cM){
    Int_t firstBin = hRes->FindFirstBinAbove(cM);
    Float_t x1 = hRes->GetBinCenter(firstBin);
    Float_t y1 = hRes->GetBinContent(firstBin);
    Float_t x2 = hRes->GetBinCenter(firstBin-1);
    Float_t y2 = hRes->GetBinContent(firstBin-1);
    Float_t m = (y2-y1)/(x2-x1);
    Float_t b = y1-m*x1;
    Float_t start = (cM-b)/m;
    cout<<"START["<<cM<<"]: "<<x1<<"("<<y1<<") - "<<x2<<"("<<y2<<") => "<<start<<endl;
    Int_t lastBin = hRes->FindLastBinAbove(cM);
    x1 = hRes->GetBinCenter(lastBin);
    y1 = hRes->GetBinContent(lastBin);
    x2 = hRes->GetBinCenter(lastBin+1);
    y2 = hRes->GetBinContent(lastBin+1);
    m = (y2-y1)/(x2-x1);
    b = y1-m*x1;
    Float_t end = (cM-b)/m;
    cout<<"End[\"<<cM<<\"]: "<<x1<<"("<<y1<<") - "<<x2<<"("<<y2<<") => "<<end<<endl;
    //   char t; cin>>t;
    return std::make_pair(start,end);
}

void TTransparentAnalysis::saveResolutionPlot(TH1F* hRes, UInt_t clusterSize,TString additionalInformation) {
    if(!hRes)
        return;
    TString section = "Resolution";
    section+=additionalInformation;
    if(alignMode==TSettings::transparentMode)
        section+"Trans";
    //	Float_t mean = hRes->GetMean();
    //	Float_t sigma = hRes->GetRMS();
    TString hName;
    TFitResultPtr resPtr = hRes->Fit("gaus","NQS");
    Float_t mean = resPtr.Get()->GetParams()[1];//->Parameter(1);
    Float_t sigma = resPtr.Get()->GetParams()[2];//Parameter(2);
    //find fwhm
    Float_t max = hRes->GetBinContent(hRes->GetMaximumBin());
    Float_t start = hRes->GetBinLowEdge(hRes->FindFirstBinAbove(max/2));
    Float_t end =  hRes->GetBinLowEdge(hRes->FindLastBinAbove(max/2)+1);
    std::pair<Float_t,Float_t > fwhm = getFWCrossingPoint(hRes,max/2.);
    //std::make_pair(hRes->GetBinLowEdge(hRes->FindFirstBinAbove(max/2)),hRes->GetBinLowEdge(hRes->FindLastBinAbove(max/2)+1));
    std::pair<Float_t,Float_t > fwtm = getFWCrossingPoint(hRes,max/3.);
    //std::make_pair(hRes->GetBinLowEdge(hRes->FindFirstBinAbove(max/3)),hRes->GetBinLowEdge(hRes->FindLastBinAbove(max/3)+1));
    results->setFloatValue(section,"FWHM_width",fwhm.second-fwhm.first);
    results->setFloatValue(section,"FWHM_mean",(fwhm.second+fwhm.first)/2.);
    results->setFloatValue(section,"FWTM_width",fwtm.second-fwtm.first);
    results->setFloatValue(section,"FWTM_mean",(fwtm.second+fwtm.first)/2.);
    start = fwhm.first;
    end = fwhm.second;
    if (end<start){
        start = end;
        end = fwhm.first;
    }
    Float_t mean2 = (start+end)/2;
    Float_t sigma2 = end-mean2;
    if (sigma2 <0) sigma2*=-1;
    TString hTitle;
    hTitle = hRes->GetTitle();
    for(int i=0;i<5;i++){
        hName = hRes->GetName();
        switch (i){
            case 0: hName.Append("_SingleGausFit");hTitle.Append(" Single Gaus Fit 2x FWHM");break;
            case 1: hName.Append("_SingleGausFitFWHM");hTitle.Append(" Single Gaus Fit FWHM");break;
            case 2: hName.Append("_DoubleGausFit");hTitle.Append(" 2 x Gaus Fit");break;
            case 3: hName.Append("_FixedGausFit");hTitle.Append(" Single Gaus Fit -20#mum - 20 #mum");break;
            case 4: hName.Append("_SingleGausFitFWTM");hTitle.Append(" Single Gaus Fit Mean of 2/3");break;
        }
        TH1F* hClone = (TH1F*)hRes->Clone(hName);
        hClone->SetTitle(hTitle);
        Float_t gaus1=-1;
        Float_t gaus3=-1;
        Float_t mean1 = -9999;
        Float_t mean3 = -9999;
        TF1* fit;
        if(hClone) {
            switch(i){
                case 0:
                    resPtr=hClone->Fit("gaus","SQ","",mean2-2*sigma2,mean2+2*sigma2);
                    if (resPtr.Get()){//todo check wh neccessary
                        gaus1 = resPtr.Get()->GetParams()[2];
                        mean1 = resPtr.Get()->GetParams()[1];
                    }
                    break;
                case 1:
                    resPtr=hClone->Fit("gaus","SQ","",fwhm.first,fwhm.second);
                    if (resPtr.Get()){//todo check wh neccessary
                        gaus1 = resPtr.Get()->GetParams()[2];
                        mean1 = resPtr.Get()->GetParams()[1];
                    }
                    break;
                case 2:
                    fit = doDoubleGaussFit(hClone);
                    if (fit){//todo check wh neccessary
                        gaus1 = fit->GetParameter(2);
                        gaus3 = fit->GetParameter(5);
                        mean1 = fit->GetParameter(1);
                        mean3 = fit->GetParameter(4);
                        if(gaus1<gaus3){
                            Float_t gaus = gaus1;
                            gaus1 = gaus3;
                            gaus3 = gaus;
                            gaus = mean1;
                            mean1 = mean3;
                            mean3 = gaus;
                        }
                    }

                    break;
                case 3:
                    resPtr= hClone->Fit("gaus","SQ","",-20,20);
                    if (resPtr.Get()){//todo check wh neccessary
                        gaus1 = fit->GetParameter(2);
                        mean1 = fit->GetParameter(1);
                    }
                    break;
                case 4:
                    resPtr=hClone->Fit("gaus","SQ","",fwtm.first,fwtm.second);
                    if (resPtr.Get()){//todo check wh neccessary
                        gaus1 = resPtr.Get()->GetParams()[2];
                        mean1 = resPtr.Get()->GetParams()[1];
                    }
                    break;
            }
            if ( clusterSize == TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)-1 && results ){
                if( i == 0 ) {
                    results->setSingleGaussianResolution(gaus1,alignMode);
                    results->setFloatValue(section,"SingleGaus_Fit",gaus1);
                    results->setFloatValue(section,"SingleGaus_Mean",mean1);
                }
                else if (i == 1 ){
                    results->setSingleGaussianShortResolution(gaus1,alignMode);
                    results->setFloatValue(section,"FWHM_Fit",gaus1);
                    results->setFloatValue(section,"FWHM_FitMean",mean1);
                }
                else if (i == 2 ){
                    results->setDoubleGaussianResolution(gaus1,gaus3,alignMode);
                    results->setFloatValue(section,"DoubleGaus_Fit1",gaus1);
                    results->setFloatValue(section,"DoubleGaus_Fit2",gaus3);
                    results->setFloatValue(section,"DoubleGaus_Mean1",mean1);
                    results->setFloatValue(section,"DoubleGaus_Mean2",mean3);
                }
                else if (i == 3 ) {
                    results->setSingleGaussianFixedResolution(gaus1,alignMode);
                    results->setFloatValue(section,"FixedGaus_Fit",gaus1);
                    results->setFloatValue(section,"FixedGaus_Mean",mean1);
                }
                else if (i == 4 ) {
                    results->setSingleGaussianFWTMResolution(gaus1,alignMode);
                    results->setFloatValue(section,"FWTM_Fit",gaus1);
                    results->setFloatValue(section,"FWTM_FitMean",mean1);
                }
            }

            histSaver->SaveHistogram(hClone);
            delete hClone;
        }
    }
}

void TTransparentAnalysis::initPHvsEventNoAreaPlots(UInt_t nStart, UInt_t nEnd) {
    cout<<"initPHvsEventNoAreaPlots"<<flush;
    Int_t nentriesPerBin = 20000;
    UInt_t nBins = (nEnd-nStart)/nentriesPerBin;
    if((nEnd-nStart)%nentriesPerBin!=0)nBins++;
    if (nBins==0)nBins=1;

    TString name ="hPHVsEventNo_clustersize10";
    TString title = "PH_{clustersize = 10} vs eventNo";
    hPHVsEventNo = new TH2D(name,title,nBins,nStart,nEnd,settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector));
    hPHVsEventNo->GetXaxis()->SetTitle("event no.");
    hPHVsEventNo->GetXaxis()->SetTitle("pulse height");

    name ="hPHVsEventNo_2outOf10";
    title = "PH_{2 out of 10} vs eventNo";
    hPH2OutOf10VsEventNo= new TH2D(name,title,nBins,nStart,nEnd,settings->getPulse_height_num_bins(),0,settings->getPulse_height_max(subjectDetector));
    hPHVsEventNo->GetXaxis()->SetTitle("event no.");
    hPHVsEventNo->GetXaxis()->SetTitle("pulse height");
    for (UInt_t i = 0; i< TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); i++){
        name = TString::Format("hPHvsEventNoArea_clusterSize_%02d",i+1);
        title = TString::Format("ph vs eventNo, clustersize %d",i+1);
        UInt_t yBins = xDivisions*yDivisions;
        TProfile2D* prof = new TProfile2D(name,title,nBins,nStart,nEnd,yBins,0,yBins);
        prof->Draw();
        prof->GetXaxis()->SetTitle("event no");
        initDividedAreaAxis(prof->GetYaxis());
        prof->GetZaxis()->SetTitle(TString::Format("avrg. pulse height, cl size %d",i+1));
        vecPHVsEventNo_Areas.push_back(prof);
        name = TString::Format("hPHvsEventNo2HighestArea_clusterSize_%02d",i+1);
        title = TString::Format("ph vs eventNo 2Highest, clustersize %d",i+1);
        prof = new TProfile2D(name,title,nBins,nStart,nEnd,yBins,0,yBins);
        prof->GetXaxis()->SetTitle("event no");
        initDividedAreaAxis(prof->GetYaxis());
        prof->GetZaxis()->SetTitle(TString::Format("avrg. pulse height 2 out of %d",i+1));
        vecPH2HighestVsEventNo_Areas.push_back(prof);
    }
    cout<<"."<<endl;
}

void TTransparentAnalysis::initClusteredHistos(UInt_t startEvent,UInt_t maxEvents) {
    Int_t nBins = (maxEvents-startEvent)/1e4;
    Float_t min = 0;
    Float_t max = settings->getPulse_height_max(subjectDetector);
    Int_t bins = settings->getPulse_height_num_bins();
    cout<<"[TTransparentAnalysis::initClusteredHistos] "<<startEvent<<" - "<<maxEvents<<"-> "<<nBins<<"/\t" <<min<<"-"<<max<<": "<<bins<<endl;
    TString name ="hLandauVsEventNo_Clustered";
    cout<<"create \""<<name<<"\""<<endl;
    hLandauVsEventNo_Clustered = new TH2F(name,name,nBins,startEvent,maxEvents,bins,min,max);
    hLandauVsEventNo_Clustered->GetXaxis()->SetTitle("EventNo.");
    hLandauVsEventNo_Clustered->GetYaxis()->SetTitle("Pulse Height_{clustered} /ADC");
    hLandauVsEventNo_Clustered->GetZaxis()->SetTitle("number of entries #");

    name = "hClusterSize_Clustered";
    cout<<"create "<<name<<endl;
    hClusterSize_Clustered = new TH1F(name,name,10,-.5,9.5);
    hClusterSize_Clustered->GetXaxis()->SetTitle("cluster Size_{clustered}");
    hClusterSize_Clustered->GetYaxis()->SetTitle("no of entries #");

    name ="hLandauVsClusterSize_Clustered";
    cout<<"create "<<name<<endl;
    hLandauVsClusterSize_Clustered = new TH2F(name,name,bins,min,max,10,-.5,9.5);
    hLandauVsClusterSize_Clustered->GetXaxis()->SetTitle("Pulse Height_{clustered} /ADC");
    hLandauVsClusterSize_Clustered->GetYaxis()->SetTitle("cluster Size_{clustered}");
    hLandauVsClusterSize_Clustered->GetZaxis()->SetTitle("no of entries #");

    name = "hNClusteres_Clustered";
    cout<<"create "<<name<<endl;
    hNClusteres_Clustered = new TH1F(name,name,10,-.5,9.5);
    hNClusteres_Clustered->GetXaxis()->SetTitle("No of clusteres_{clustered}");
    hNClusteres_Clustered->GetYaxis()->SetTitle("number of entries");

}

void TTransparentAnalysis::savePHvsEventNoAreaPlots() {
    cout<<"[TTransparentAnalysis::savePHvsEventNoAreaPlots] "<<endl;
    histSaver->SaveHistogram(hPH2OutOf10VsEventNo);
    histSaver->Save1DProfileYWithFitAndInfluence(hPH2OutOf10VsEventNo,"pol1",true);
    delete hPH2OutOf10VsEventNo;
    histSaver->SaveHistogram(hPHVsEventNo);
    histSaver->Save1DProfileYWithFitAndInfluence(hPHVsEventNo,"pol1",true);
    delete hPHVsEventNo;
    for (UInt_t i = 0; i< vecPHVsEventNo_Areas.size(); i++){
        TProfile2D * prof2d = vecPHVsEventNo_Areas[i];
        if (!prof2d) continue;
        TAnalysisOfSelection::savePHvsEventNoAreaPlots(histSaver,settings,prof2d,xDivisions,yDivisions);
        //
        //        prof2d->Draw();
        //        TProfile* prof = histSaver->GetProfileX(prof2d);
        //        histSaver->Save1DProfileXWithFitAndInfluence(prof,0,false);
        //        histSaver->SaveProfile2DWithEntriesAsText(prof2d,true);//,false);
        //        TF1* pol1Fit = new TF1("pol1Fit","pol1",prof2d->GetXaxis()->GetXmin(),prof2d->GetXaxis()->GetXmax());
        //        vector<TProfile*> vecStack;
        //        for(int y = 0; y< prof2d->GetYaxis()->GetNbins();y++){
        //            TString name = prof2d->GetName()+(TString)"_"+GetNameOfArea(y%xDivisions,y/xDivisions);
        //            prof = histSaver->GetProfileX(prof2d,name,y+1,y+1);
        //            prof->Draw();
        //            prof->SetLineColor(y);
        //            prof->SetMarkerColor(y);
        //            TF1* fit = (TF1*) pol1Fit->Clone(prof->GetName()+(TString)"_fit");
        //            cout<<"Save: "<<prof->GetName()<<endl;
        //            histSaver->Save1DProfileXWithFitAndInfluence(prof,fit);
        //            vecStack.push_back(prof);
        //        }
        //        THStack* stack = new THStack("hPHvsEventNoAllAreas","Ph vs EventNo, All Areas");
        //        bool foundOneHisto = false;
        //        for(UInt_t i = 0; i< vecStack.size();i++){
        //            if(!vecStack[i]) continue;
        //            foundOneHisto = true;
        //            stack->Add(vecStack[i]);
        //        }
        //        if(foundOneHisto){
        //            stack->Draw();
        //            if(stack->GetXaxis())stack->GetXaxis()->SetTitle("Event No.");
        //            if(stack->GetYaxis())stack->GetYaxis()->SetTitle("avrg PH");
        //            cout<<"Save: "<<stack->GetName()<<endl;
        //            histSaver->SaveStack(stack);
        //        }
        //        delete stack;
        //        for(UInt_t i = 0; i< vecStack.size();i++){
        //            delete vecStack[i];
        //            vecStack[i]=0;
        //        }
    }
    //    char t; cin>>t;
}

TString TTransparentAnalysis::GetNameOfArea(Int_t x,Int_t y){

    vector<TString> xDirString;
    xDirString.push_back("Left");
    xDirString.push_back("Middle");
    xDirString.push_back("Right");
    vector<TString> yDirString;
    yDirString.push_back("low");
    yDirString.push_back("mid");
    yDirString.push_back("up");
    if(0<=x&&x<xDirString.size() && 0<=y&&y < yDirString.size())
        return yDirString[y]+(TString)"_"+xDirString[x];
    if (x==-1&&0<=y&&y < yDirString.size())
        return yDirString[y];
    if(y==-1&&0<=x&&x<xDirString.size())
        return xDirString[x];
    return "unkown area";
}
void TTransparentAnalysis::initDividedAreaAxis(TAxis* axis){
    if (!axis) return;
    Int_t bins = xDivisions*yDivisions;
    if (axis->GetNbins() != bins ){
        cerr<<"divided area axis has wrong number of bins: "<<axis->GetNbins()<<"/"<<bins<<endl;
        return;
    }
    for (UInt_t y = 0; y < yDivisions; y++)
        for (UInt_t x = 0; x < xDivisions; x++){
            axis->SetBinLabel(x+xDivisions*y+1, GetNameOfArea(x,y));
        }

}

void TTransparentAnalysis::fillPHvsEventNoAreaPlots(UInt_t area, UInt_t clusterSize, UInt_t charge, UInt_t chargeOfTwo) {
    if(clusterSize == TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)){
        hPH2OutOf10VsEventNo->Fill(nEvent,chargeOfTwo);
        hPHVsEventNo->Fill(nEvent,charge);
    }
    UInt_t i = clusterSize -1;
    if (i < vecPHVsEventNo_Areas.size())
        vecPHVsEventNo_Areas[i]->Fill(nEvent,area,charge);
    if (i < vecPH2HighestVsEventNo_Areas.size())
        vecPH2HighestVsEventNo_Areas[i]->Fill(nEvent,area,chargeOfTwo);
}



UInt_t TTransparentAnalysis::GetHitArea(TSettings* set,Float_t xVal,Float_t yVal,UInt_t xDivisions,UInt_t yDivisions){
    Int_t nFidCut = set->getSelectionFidCuts()->getFidCutRegion(xVal,yVal);
    TFiducialCut* fidcut = set->getSelectionFidCuts()->getFidCut(nFidCut);
    Float_t relX = (xVal-fidcut->GetXLow())/(fidcut->GetXHigh()-fidcut->GetXLow());
    Float_t relY = (yVal-fidcut->GetYLow())/(fidcut->GetYHigh()-fidcut->GetYLow());
    Int_t x = relX*xDivisions;
    Int_t y = relY*yDivisions;
    return x+xDivisions*y;
}
