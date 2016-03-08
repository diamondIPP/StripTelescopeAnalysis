/*
 * TAnalysisOfAnalysisDifferences.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: bachmair
 */

#include "../include/TAnalysisOfAnalysisDifferences.hh"

TAnalysisOfAnalysisDifferences::TAnalysisOfAnalysisDifferences(TSettings* settings,HistogrammSaver* histSaver,TString extension) {
    // TODO Auto-generated constructor stub
    cout<<"TAnalysisOfAnalysisDifferences:: "<<extension<<endl;
    this->settings = settings;

    if(settings==0 || histSaver==0)
        cerr<<"ERROR: invalid settings or histogram saver"<<endl;
    this->histSaver = new HistogrammSaver(settings);
    histSaver->SetNumberOfEvents(settings->getNEvents());
    oldPlotPath = histSaver->GetPlotsPath();
    this->histSaver->SetPlotsPath(oldPlotPath+(string)"/negativeCharges/");

    transparentMap = 0;
    clusteredMap = 0;
    predictedPositions = 0;
    stripHisto = 0;
    verbosity =settings->getVerbosity();
    negChargeCut = settings->getNegativeChargeCut();
    this->extension = "_"+extension;
    hPhantomLandau = 0;
}

TAnalysisOfAnalysisDifferences::~TAnalysisOfAnalysisDifferences() {
    delete histSaver;
    if (stripHisto)
        delete stripHisto;
    cout<<"Done with TAnalysisOfAnalysisDifferences"<<endl;
}

void TAnalysisOfAnalysisDifferences::setStripHistogram(TH1F* histo){
    this->stripHisto = (TH1F*)histo->Clone();
    if(stripHisto)
        stripHisto->SetLineColor(kBlue);
}
std::map<Int_t, TCluster>* TAnalysisOfAnalysisDifferences::getClusteredMap() const {
    return clusteredMap;
}

void TAnalysisOfAnalysisDifferences::setClusteredMap(
        std::map<Int_t, TCluster>* clusteredMap) {
    this->clusteredMap = clusteredMap;
}

std::map<Int_t, TCluster>* TAnalysisOfAnalysisDifferences::getTransparentMap() const {
    return transparentMap;
}

void TAnalysisOfAnalysisDifferences::setTransparentMap(
        std::map<Int_t, TCluster>* transparentMap) {
    this->transparentMap = transparentMap;
}


std::map<Int_t,std::pair<Float_t,Float_t> >* TAnalysisOfAnalysisDifferences::getPredictedPositions() const {
    return predictedPositions;
}

void TAnalysisOfAnalysisDifferences::setPredictedPositions(
        std::map<Int_t, pair<Float_t,Float_t> > *predictedPostionMap){
    this->predictedPositions = predictedPostionMap;
    cout<<"set predicted Positions"<<endl;
    cout<<predictedPostionMap->size()<<endl;
    cout<<predictedPositions->size()<<endl;
    cout<<"next"<<endl;
}

void TAnalysisOfAnalysisDifferences::Analysis() {
    char t;
    cout<<"\n\n[TAnalysisOfAnalysisDifferences]::Analysis"<<endl;
    if(!settings||!histSaver){
        return;
    }
    if(!transparentMap || !clusteredMap||!predictedPositions){
        cerr<<"ERROR: Maps not all set"<<endl;
        return;
    }
    cout<<"Clustered Analysis:  "<<clusteredMap->size()<<endl;
    cout<<"TransparentAnalysis: "<<transparentMap->size()<<endl;
    cout<<"predictedPostions: "<<predictedPositions->size()<<endl;
    InitHistograms();
    nSameEvents = 0;
    nOnlyClustered = 0;
    nOnlyTransparent = 0;
    LoopOverBothMaps();

//    cout<<"Save histos, Key Press:"<<flush; cin >>t;
    SaveHistograms();

    cout<<"\nSame Events:      "<<nSameEvents<<endl;
    cout<<"Only Clustered:   "<<nOnlyClustered<<endl;
    cout<<"Only Transparent: "<<nOnlyTransparent<<endl;
    if(verbosity%2==1){
        cout<<"Press a key"<<endl;
        cin>>t;
    }
}


void TAnalysisOfAnalysisDifferences::LoopOverBothMaps(){
    cout<<"set iterators"<<endl;
    itClustered = clusteredMap->begin();
    itTransparent = transparentMap->begin();
    itPredicted = predictedPositions->begin();
    bool endLoop = false;
    bool endClustered;
    bool endTransparent;
    endClustered =  itClustered== clusteredMap->end();
    endTransparent = itTransparent== transparentMap->end();
    cout<<"loop over events: "<<itPredicted->first<<endl;
//    char t; cout<<"Key Press:"<<flush; cin >>t;
    while (!endLoop){
        if (verbosity>2)
            cout<<"TAnalysisOfAnalysisDifferences::looping"<<flush;
        UpdatePredictedPosition();

        if (itClustered->first == itTransparent->first){
            AnalyseSameEvent();
            if( itClustered != clusteredMap->end()) itClustered++;
            if(itTransparent != transparentMap->end()) itTransparent++;

        }
        else if (itClustered->first < itTransparent->first){
            if(!endClustered){
                AnalyseOnlyClusteredEvent();
                itClustered++;
            }
            else{
                AnalyseOnlyTransparentEvent();
                itTransparent++;
            }
        }
        else if (itClustered->first > itTransparent->first){
            if (!endTransparent){
                AnalyseOnlyTransparentEvent();
                itTransparent++;
            }
            else{
                AnalyseOnlyClusteredEvent();
                itClustered++;
            }
        }
        endClustered =  itClustered== clusteredMap->end();
        endTransparent = itTransparent== transparentMap->end();
        endLoop = endClustered && endTransparent;
    }
    cout<<"done with looping"<<endl;
    if(verbosity%2==1){
        char t;
        cout<<"Press a key and enter"<<endl;
        cin>>t;
    }
}

void TAnalysisOfAnalysisDifferences::AnalyseSameEvent() {

    if (itClustered->second.getClusterSize()==0 ||
        itTransparent->second.getClusterSize()==0){
        cout<<"invalid cluster size"<<endl;
        return;
    }
    AnalyseTransparentEvent();

    if(itPredicted->first != itTransparent->first || itPredicted->first != itClustered->first)
        cout<<"predicted,transparent and clustered iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<" "<<itClustered->first<<endl;
    Float_t clusteredCharge = itClustered->second.getCharge(true);
    itTransparent->second.SetTransparentClusterSize(3);
    Float_t transparentCharge = itTransparent->second.getCharge(true);
    if (itClustered->second.isSaturatedCluster()||itTransparent->second.isSaturatedCluster())
        return;
    Float_t posCharge = itTransparent->second.getPositiveCharge(true);//100,true,true);
    Int_t eventNo = itClustered->first;
    mapHistos["hChargeDifference"]->Fill(transparentCharge-clusteredCharge);
    mapHistos["hClusteredPulseHeight"]->Fill(itClustered->second.getPositiveCharge());
    if ( posCharge - clusteredCharge < 0 && verbosity>0 && itClustered->second.getClusterSize()<4){
        cout<<"posCharge - clusteredCharge " <<eventNo<<"\t"<<posCharge-clusteredCharge<<"\n\tposCharge: "<<posCharge<<"\n\tclusCharge: "<<clusteredCharge<<endl;
        cout<<"\tclustered: ";
        itClustered->second.Print(2);
        cout<<"\ttransparent: ";
        itTransparent->second.Print(2);
        cout<<endl;
    }
    mapHistos["hPositiveTransparentCharge_Minus_ClusteredCharge"]->Fill(posCharge-clusteredCharge);
    //    cout<<"\t"<<transparentCharge<<"\t"<<posCharge<<" = "<<transparentCharge-posCharge<<" "<<clusteredCharge-posCharge<<endl;
    Float_t ratio = clusteredCharge/transparentCharge;
    //    cout<<itClustered->first<<" both\t"<<(ratio-1.)*100<<" %"<<endl;
    if( verbosity>3 && ( ratio-1 < -.1 || ratio-1 >.1)){
        cout<<"\tClustered\n\t";
        itClustered->second.Print(2);
        cout<<"\tTransparent\n\t";
        itTransparent->second.Print(2);
        cout<<"\n";
    }
    nSameEvents ++;
    mapHistos["hGoodCellsEventTypes"]->Fill(0);
}

void TAnalysisOfAnalysisDifferences::AnalyseTransparentEvent() {
    itTransparent->second.SetTransparentClusterSize(3);
    Float_t transparentCharge = itTransparent->second.getCharge(true);
    if (itTransparent->second.isSaturatedCluster())
        return;
    Int_t eventNo = itTransparent->first;
    Float_t posCharge = itTransparent->second.getPositiveCharge(true);//100,true,true);
    Float_t allCharge = itTransparent->second.getCharge(true,true);
    mapHistos["hTransparentPulseHeight"]->Fill(posCharge);
    mapHistos["hTransparentPulseHeightAllCharge"]->Fill(allCharge);
    if (1810<posCharge&&posCharge<1840){
        cout<<"PH: "<<posCharge;
        itTransparent->second.Print(1);
    }

    Int_t pos = 0;
    Float_t negCharge = 0;
    bool hasNegativeCharge = itTransparent->second.hasNegativeCharge(negCharge,pos,true);
    if (negCharge == 0)
        hasNegativeCharge = itTransparent->second.hasNegativeCharge(negCharge,pos,true,true);
    Float_t maxCharge = itTransparent->second.getHighestSignal(true);
    Float_t charge = itTransparent->second.getCharge(true);

    Int_t firstPos = itTransparent->second.getTransparentClusterPosition(0);
    Int_t secondPos = itTransparent->second.getTransparentClusterPosition(1);
    Int_t thirdPos = itTransparent->second.getTransparentClusterPosition(2);
    Float_t firstCharge = itTransparent->second.getSignal(firstPos);
    Float_t secondCharge = itTransparent->second.getSignal(secondPos);
    Float_t thirdCharge = itTransparent->second.getSignal(thirdPos);
//    cout<<"\nTransparent: "<<firstCharge<<" "<<secondCharge<<" "<<thirdCharge<<" "<<charge<<" "<<negCharge<<" "<<pos<<endl;
//    itTransparent->second.Print(1);

    Float_t lowThr = settings->getResponseWindow().first;
    Float_t highThr = settings->getResponseWindow().second;

    mapHistos["hNegativeChargeRatio"]->Fill( negCharge/charge,charge);
    mapHistos["hNegativeChargeRatioAbs"]->Fill( negCharge/charge,negCharge);
    mapHistos["hNegativeChargeRatioMax"]->Fill(negCharge/maxCharge,maxCharge);


    if(predictedPositions->count(eventNo)){
        Float_t xPredDet = predictedPositions->at(eventNo).first;
        Float_t yPredDet = predictedPositions->at(eventNo).second;
        pair<Float_t,Float_t> relPos =  settings->getRelativePositionInCell(xPredDet,yPredDet);
        if (lowThr < posCharge && posCharge < highThr){
            mapHistos["hChargeResponseWindowPosition"]->Fill(xPredDet,yPredDet);
            mapHistos["hChargeResponseWindowRelPosition"]->Fill(relPos.first,relPos.second);
        }
        mapHistos["hAllEvents_RelPosition"]->Fill(relPos.first,relPos.second);
        TProfile2D* prof = (TProfile2D*)mapHistos["hNegativeChargeProfileRelPosition"];
        cout<<"hNegativeChargeProfileRelPosition - FILL: "<<relPos.first<<","<<relPos.second<<": "<<negCharge<<endl;
        prof->Fill(relPos.first,relPos.second,negCharge);
        prof = (TProfile2D*)mapHistos["hNegativeChargeRatioOverlay"];
        prof->Fill(relPos.first,relPos.second,negCharge/maxCharge);
        if (firstCharge>secondCharge){
            prof = (TProfile2D*) mapHistos["hAdjacentChargeRatioOverlay"];
            prof->Fill(relPos.first,relPos.second,secondCharge/firstCharge);
            mapHistos["hAdjacentChargeRatio"]->Fill(secondCharge/firstCharge,maxCharge);
            if (!settings->IsOnTheEdgeOfCell(relPos)){
                prof = (TProfile2D*) mapHistos["hAdjacentChargeRatioOverlayNoBorder"];
                prof->Fill(relPos.first,relPos.second,secondCharge/firstCharge);
                mapHistos["hAdjacentChargeRatioNoBorder"]->Fill(secondCharge/firstCharge,maxCharge);
            }
//            hAdjacentChargeRatioOverlayNoBorder;
        }
    }
    if (hasNegativeCharge){
        Int_t pos_hit = itTransparent->second.getTransparentClusterPosition(0);
        Int_t pos_neg = itTransparent->second.getTransparentClusterPosition(pos-1);
        Int_t delta = pos_neg - pos_hit;

        mapHistos["hNegativeChargePosition"]->Fill(negCharge,pos);
        mapHistos["hNegativeChargePositionTransparent"]->Fill(negCharge,pos);
        mapHistos["hNegativeChargeChannelPositionTransparent"]->Fill(negCharge,delta);
        mapHistos["hNegativeCharge"]->Fill(negCharge);
        if(negCharge<negChargeCut){
            if(predictedPositions->count(eventNo)){
                Float_t xPredDet = predictedPositions->at(eventNo).first;
                Float_t yPredDet = predictedPositions->at(eventNo).second;
                mapHistos["hNegativeChargeAboveCut_Position"]->Fill(xPredDet,yPredDet);
                pair<Float_t,Float_t> relPos =  settings->getRelativePositionInCell(xPredDet,yPredDet);
                mapHistos["hNegativeChargeAboveCut_RelPosition"]->Fill(relPos.first,relPos.second);
            }
        }
    }
    else
        mapHistos["hNegativeChargePositionTransparent"]->Fill(negCharge,-1);
    if(hasNegativeCharge && negCharge < negChargeCut)
        mapHistos["hHasNegativeCharge_PulseHeight"]->Fill(posCharge);
    else{
        mapHistos["hNoNegativeCharge_PulseHeight"]->Fill(posCharge);

        if (predictedPositions->count(eventNo)){
            Float_t xPredDet = predictedPositions->at(eventNo).first;
            Float_t yPredDet = predictedPositions->at(eventNo).second;
            pair<Float_t,Float_t> relPos =  settings->getRelativePositionInCell(xPredDet,yPredDet);
            mapHistos["hNoNegativeCharge_RelPosition"]->Fill(relPos.first,relPos.second);
            cout<<"-->"<<((TH2F*)(mapHistos["hNoNegativeCharge_RelPosition"]))->GetEntries()<<endl;;
            if (posCharge < settings->getLowResponseThreshold()){
                mapHistos["hNoNegativeChargeLowResponsePosition"]->Fill(xPredDet,yPredDet);
            }
            if (lowThr < posCharge && posCharge < highThr)
                mapHistos["hNoNegativeChargeResponseWindowPosition"]->Fill(xPredDet,yPredDet);
        }
        else
            cout<<"Cannot find PredictedPosition of Event No: "<<eventNo<<endl;
    }
    mapHistos["hPositive_Minus_Negative_TransparentCharge"]->Fill(posCharge-transparentCharge);

}

void TAnalysisOfAnalysisDifferences::AnalyseOnlyTransparentEvent() {

    if (itTransparent->second.getClusterSize()==0){
        cout<<"invalid cluster size"<<endl;
        return;
    }
    AnalyseTransparentEvent();
    nOnlyTransparent++;
    if( verbosity>3)
        cout<<itTransparent->first<< " only in Transparent Analysis"<<endl;

    if(itPredicted->first != itTransparent->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<endl;
    mapHistos["hOnlyTranspClusterCharge"]->Fill(itTransparent->second.getPositiveCharge());

    mapHistos["hOnlyTranspClusterPosition"]->Fill(itPredicted->second.first,itPredicted->second.second);
    mapHistos["hGoodCellsEventTypes"]->Fill(1);
    Int_t pos = 0;
	Float_t charge = 0;
    bool hasNegativeCharge = itTransparent->second.hasNegativeCharge(charge,pos,true);
	if (hasNegativeCharge){
		mapHistos["hNegativeChargePositionTransparent"]->Fill(charge,pos);
	}
	else
		mapHistos["hNegativeChargePositionTransparent"]->Fill(charge,-1);

}

void TAnalysisOfAnalysisDifferences::set3DPhantomLandau(TH1F* hPhantomLandau) {
    this->hPhantomLandau = hPhantomLandau;
    if (!this->hPhantomLandau)
        return;
    this->hPhantomLandau->SetName("hPhantom");
    this->hPhantomLandau->SetTitle("Phantom");
}


void TAnalysisOfAnalysisDifferences::AnalyseOnlyClusteredEvent() {

    if (itClustered->second.getClusterSize()==0){
        cout<<"invalid cluster size"<<endl;
        return;
    }

    if( verbosity>3)
        cout<<itClustered->first<< " only in Clustered Analysis"<<endl;
    nOnlyClustered++;
    if(itPredicted->first != itClustered->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itClustered->first<<endl;

    mapHistos["hClusteredPulseHeight"]->Fill(itClustered->second.getPositiveCharge());
    mapHistos["hOnlyClusteredClusterCharge"]->Fill(itClustered->second.getPositiveCharge());
    mapHistos["hOnlyClusteredClusterPosition"]->Fill(itPredicted->second.first,itPredicted->second.second);
    mapHistos["hGoodCellsEventTypes"]->Fill(-1);
}

void TAnalysisOfAnalysisDifferences::InitHistograms() {
    this->InitSameHistos();
    this->InitTransparentHistos();
    this->InitClusteredHistos();
    TH1* histo;
    Int_t bins = stripHisto?stripHisto->GetNbinsX():256;
    Int_t xmin = stripHisto?stripHisto->GetXaxis()->GetXmin():1;
    Int_t xmax = stripHisto?stripHisto->GetXaxis()->GetXmax():2800;

    TString title;
    TString name = "hTransparentPulseHeight";
    TString hname = name+extension;
    histo = new TH1F(hname,hname,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("Charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    mapHistos[name] = histo;

    name = "hTransparentPulseHeightAllCharge";
    hname = name+extension;
    histo = (TH1F*) histo->Clone(hname);
    histo->GetXaxis()->SetTitle("Charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    mapHistos[name] = histo;

    name = "hClusteredPulseHeight";
    hname = name + extension;
    histo = (TH1F*) histo->Clone(hname);
    histo->SetLineColor(kBlack);
    mapHistos[name] = histo;

    name = "hGoodCellsEventTypes";
    hname=name+extension;
    histo = new TH1F(hname,hname,3,-1.5,1.5);
    histo->GetXaxis()->SetBinLabel(1,"Only Clustered");
    histo->GetXaxis()->SetBinLabel(2,"Same ");
    histo->GetXaxis()->SetBinLabel(3,"Only transparent");
    mapHistos[name] = histo;

    name ="hHasNegativeCharge_PulseHeight";
    hname=name+extension;
    histo = new TH1F(hname,hname,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    histo->SetLineColor(kGreen);
    mapHistos[name] = histo;

    name ="hNoNegativeCharge_PulseHeight";
    hname = name +extension;
    histo = new TH1F(hname,hname,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("charge / ADC");
    histo->GetYaxis()->SetTitle("entries a.u.");
    histo->SetLineColor(kBlack);
    mapHistos[name] = histo;

    name = "hNegativeChargeRatioOverlay";
    hname = name + extension;
    histo = settings->GetOverlayProfile(hname);
    histo->SetZTitle("avrg min. adjacent Signal");
    mapHistos[name] = histo;

    name = "hAdjacentChargeRatioOverlay";
    hname = name + extension;
    histo = settings->GetOverlayProfile(hname);
    histo->SetZTitle("avrg. adjacent Signal");
    mapHistos[name] = histo;

    name = "hAdjacentChargeRatioOverlayNoBorder";
    hname = name + extension;
    title = TString::Format("hAdjacentChargeRatioOverlayNo Border (%.0f #mum) ",settings->GetMinimumEdgeDistance());
    histo = settings->GetOverlayProfile(hname);
    histo->SetTitle(title);
    histo->SetZTitle("avrg. adjacent Signal");
    mapHistos[name] = histo;

    //hAdjacentChargeRatioOverlayNoBorder

    name = "hAdjacentChargeRatioNoBorder";
    hname = name +extension;
    title = TString::Format("Adjacent Charge Ratio No Border (%.0f #mum) ",settings->GetMinimumEdgeDistance());
    title+=extension;
    title+="; signal ratio: S_{Adjacent}/PH; Pulse Heigth / ADC;number of entries";
    histo = new TH2D(hname,title,1000,-.5,.5,bins/4,xmin,xmax);
    mapHistos[name] = histo;

    name = "hAdjacentChargeRatio";
    hname = name +extension;
    title = "Adjacent Charge Ratio "+extension;
    title+="; signal ratio: S_{Adjacent}/PH; Pulse Heigth / ADC;number of entries";
    histo = new TH2D(hname,title,1000,-.5,.5,bins/4,xmin,xmax);
    mapHistos[name] = histo;

    name = "hNegativeChargeRatio";
    hname = name +extension;
    title = "Negative Charge Ratio "+extension;
    title+="; signal ratio: S_{Min}/PH; Pulse Heigth / ADC;number of entries";
    histo = new TH2D(hname,title,1000,-.5,.5,bins/4,xmin,xmax);
    mapHistos[name] = histo;

    name = "hNegativeChargeRatioAbs";
    title = "Negative Charge Ratio "+extension;
    title+="; signal ratio: S_{Min}/PH; S_{Min} / ADC;number of entries";
    hname = name +extension;
    histo = new TH2D(hname,title,1000,-.5,.5,bins/4,-400,400);
    mapHistos[name] = histo;

    name = "hNegativeChargeRatioMax";
    hname = name+extension;
    title = "Negative Charge Ratio Max "+extension;
    title+="; signal ratio: S_{Min}/PH_{max}; Max. Pulse Height / ADC;number of entries";
    histo = new TH2D(hname,title,1000,-.5,.5,bins/4,xmin,xmax);
    mapHistos[name] = histo;
}


void TAnalysisOfAnalysisDifferences::SaveTransparentClusteredComparison() {
    TH1F* histo1 = (TH1F*)(((TH1F*)(mapHistos["hClusteredPulseHeight"]))->Clone());
    TH1F* histo2 = (TH1F*)(((TH1F*)(mapHistos["hTransparentPulseHeight"]))->Clone());
    histo1->SetTitle("Clustered");
    histo1->SetLineColor(kBlack);
    histo2->SetTitle("Transparent");
    histo2->SetLineColor(kRed);
    cout<<"histo1: "<<histo1->GetName()<<" | "<<histo1->GetTitle()<<": "<<histo1->GetEntries()<<endl;;
    cout<<"histo2: "<<histo2->GetName()<<" | "<<histo2->GetTitle()<<": "<<histo2->GetEntries()<<endl;;
    if (stripHisto)
        cout<<"strip"<<stripHisto->GetName()<<" | "<<stripHisto->GetTitle()<<": "<<stripHisto->GetEntries()<<endl;;
    TH1F* hPhantomScaled = 0;
    TH1F* hStripScaled = 0;
    if (stripHisto){
        stripHisto->SetLineColor(kBlue);
        stripHisto->SetTitle("Strip");
        hStripScaled = (TH1F*)stripHisto->Clone();
        Float_t scale = hStripScaled->GetBinContent(hStripScaled->GetMaximumBin());
        hStripScaled->Scale(1./scale);
    }
    if (hPhantomLandau){
        hPhantomScaled = (TH1F*)hPhantomLandau->Clone();
        Float_t scale = hPhantomScaled->GetBinContent(hPhantomScaled->GetMaximumBin());
        hPhantomScaled->Scale(1./scale);
    }
    histSaver->SaveTwoHistos("cComparisionPulseHeights"+extension,histo1,histo2);
    histSaver->SaveTwoHistosNormalized("cComparisionPulseHeightsNormalized"+extension,histo1,histo2);

    TString name = "cComparisionPulseHeightWithStrip"+extension;
    TString title = name+";Charge / ADC; number of entries";
    THStack *hs = new THStack(name,title);
    hs->Add(histo1);
    hs->Add(histo2);
    if(stripHisto){
        Float_t scale = histo2->GetBinContent(histo2->GetMaximumBin());
        scale += histo1->GetBinContent(histo1->GetMaximumBin());
        scale /= 2.;
        scale /= hStripScaled->GetBinContent(hStripScaled->GetMaximumBin());
        hStripScaled->Scale(scale);
        hs->Add(hStripScaled);
    }
    histSaver->SaveStack(hs,"nostack",true);

    name = "cComparisionPulseHeightWithStrip_scaled"+extension;
    title = name+";Charge / ADC; a.u.";
    THStack *hs1 = new THStack(name,title);
    TH1F* h1 = histo1;
    TH1F* h2 = histo2;
    h1->Scale(1./h1->GetBinContent(h1->GetMaximumBin()));
    h2->Scale(1./h2->GetBinContent(h2->GetMaximumBin()));
    hs1->Add(h1);
    hs1->Add(h2);
    if(stripHisto){
        Float_t scale = 1./hStripScaled->GetBinContent(hStripScaled->GetMaximumBin());
        hStripScaled->Scale(scale);
        hs1->Add(hStripScaled);
    }
    histSaver->SaveStack(hs1,"nostack",true);
    cout<<"Save6"<<endl;
    delete hs;
    delete histo1;
    delete histo2;
    cout<<"Save7"<<endl;
}


void TAnalysisOfAnalysisDifferences::SaveComparisonPlots(TString name, TH1* histo,bool includePhantom) {
    cout<<"TAnalysisOfAnalysisDifferences::SaveComparisonPlots"<<name<<endl;
    if (!histo) return;
    if (!hPhantomLandau) includePhantom=false;
    else{
        hPhantomLandau->SetTitle("3D Phantom");
        hPhantomLandau->SetLineStyle(4);
        hPhantomLandau->SetLineColor(kRed);
    }

    if (stripHisto){
        stripHisto->SetTitle("Strip");
        stripHisto->SetLineColor(kBlue);
        stripHisto->SetLineStyle(2);
    }
    TString hName = name+extension;
    hName.Insert(0,"c");
    TCanvas *c1 = new TCanvas(hName,hName);
    c1->SetTitle(hName);
    c1->SetName(hName);
    c1->Clear();
    histo->SetTitle("3D");
    histo->Draw();
    if (stripHisto) stripHisto->Draw("same");
    TLegend* leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->Draw();
    histSaver->SaveCanvas(c1);

    hName = "stack_"+name+extension;
    TString title = name +extension+";charge / ADC;number of entries";
    THStack *hs1 = new THStack(hName,title);
    hs1->Add(histo);
    if(stripHisto)hs1->Add(stripHisto);
    if (includePhantom) hs1->Add(hPhantomLandau);
    hs1->Draw("");
    histSaver->SaveStack(hs1,"nostack",true,false);

    hName = "stack_"+name+"_normalized"+extension;
    THStack *hs2 = new THStack(hName,title);
    hs2->Add(histo->DrawNormalized("goff"));
    if (includePhantom) hs2->Add(hPhantomLandau->DrawNormalized("goff"));
    if(stripHisto)hs2->Add(stripHisto->DrawNormalized("goff"));
    hs2->Draw("");

    hName = "stack_"+name+"_scaled"+extension;
    THStack *hs3 = new THStack(hName,title);
    TH1F* h_scaled = (TH1F*)histo->Clone();
    Float_t scale = h_scaled->GetBinContent(h_scaled->GetMaximumBin());
    h_scaled->Scale(1./scale);
    TH1F* hStrip_scaled =0;
    if (stripHisto){
        hStrip_scaled = (TH1F*)stripHisto->Clone();
        scale = hStrip_scaled->GetBinContent(hStrip_scaled->GetMaximumBin());
        hStrip_scaled->Scale(1./scale);
        hStrip_scaled->SetTitle("Strip");
    }
    TH1F* hPhantom_scaled =0;
    if (includePhantom){
        hPhantom_scaled = (TH1F*)hPhantomLandau->Clone();
        scale = hStrip_scaled->GetBinContent(hPhantom_scaled->GetMaximumBin());
        hPhantom_scaled->Scale(1./scale);
        hPhantom_scaled->SetTitle("Phantom");
    }
    hs3->Add(h_scaled);
    if (hStrip_scaled) hs3->Add(hStrip_scaled);
    if (hPhantom_scaled)hs3->Add(hPhantom_scaled);
    hs3->Draw("");
    histSaver->SaveStack(hs3,"nostack",true,false,"charge / ADC","rel. no of entries");
    cout<<"TAnalysisOfAnalysisDifferences::SaveComparisonPlots DONE"<<name<<endl;

}


void TAnalysisOfAnalysisDifferences::SaveHistograms() {
    SaveTransparentClusteredComparison();
    cout<<"TAnalysisOfAnalysisDifferences::SaveHistograms1"<<endl;
    TString name = "NoNegativeCharge_PulseHeight_Comparision";
    SaveComparisonPlots(name,mapHistos["hNoNegativeCharge_PulseHeight"]);
    name= "PositiveCharge_PulseHeight_Comparision";
    SaveComparisonPlots(name,mapHistos["hTransparentPulseHeight"]);
    name= "AllCharge_PulseHeight_Comparision";
    SaveComparisonPlots(name,mapHistos["hTransparentPulseHeightAllCharge"]);

    mapHistos["hHasNegativeCharge_PulseHeight"]->SetTitle("hasNegativeCharge");
    mapHistos["hNoNegativeCharge_PulseHeight"]->SetTitle("noNegativeCharge");
    name = "stackPulseHeights_NegativeChargeComparison"+extension;
    THStack *hs = new THStack(name,name);
    hs->Add((TH1F*)mapHistos["hHasNegativeCharge_PulseHeight"]->Clone());
    hs->Add((TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone());
    hs->Draw("");
    if(hs->GetXaxis())hs->GetXaxis()->SetTitle("charge / ADC");
    if(hs->GetYaxis())hs->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs,"",true,false);


    name = "stackPulseHeights_NegativeChargeComparison_nostack"+extension;
    hs->SetTitle(name);
    hs->SetName(name);
    histSaver->SaveStack(hs,"nostack",true,false);
    delete hs;

    name = "stackPulseHeights_NegativeChargeComparison_scaled"+extension;
    hs = new THStack(name,name);
    TH1F* hHasNegative_scaled = (TH1F*)mapHistos["hHasNegativeCharge_PulseHeight"]->Clone();
    Float_t factor = hHasNegative_scaled->GetBinContent(hHasNegative_scaled->GetMaximumBin());
    hHasNegative_scaled->Scale(1./factor);
    TH1F* hNoNegative_scaled = (TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone();
    factor = hNoNegative_scaled->GetBinContent(hNoNegative_scaled->GetMaximumBin());
    hNoNegative_scaled->Scale(1./factor);
    hs->Add(hHasNegative_scaled);
    hs->Add(hNoNegative_scaled);
    hs->Draw("");
    if(hs->GetXaxis())hs->GetXaxis()->SetTitle("charge / ADC");
    if(hs->GetYaxis())hs->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs,"nostack",true,false);
    delete hs;

    std::map<TString,TH1*>::iterator it = mapHistos.begin();
    cout<<"[TAnalysisOfAnalysisDifferences]Save Histos it map "<<mapHistos.size()<<endl;
    for(it;it!=mapHistos.end();it++){
        cout<<"Save: "<<it->first<<"\t"<<it->second->GetEntries()<<endl;
        if (it->second->GetEntries()==0)
            continue;
        TString className = it->second->ClassName();
        if (className.Contains("TH2") || className.Contains("TProfile2D"))
            histSaver->SaveHistogram((TH2F*)it->second);
        else
            histSaver->SaveHistogram(it->second);
        if (it->first.Contains("NegativeChargeProfileRelPosition") || it->first.Contains("RatioOverlay") )
            histSaver->SaveNegativeChargeOverlay((TProfile2D*)it->second);
        if (it->first == "hNegativeChargeRatio"  || it->first == "hNegativeChargeRatioMax"
                || it->first == "hAdjacentChargeRatio" || it->first == "hAdjacentChargeRatioNoBorder" ){
            histSaver->SaveProjectionX((TH2*)it->second);
            histSaver->SaveProjectionY((TH2*)it->second);
        }
        if (it->first.Contains("hNegativeChargePosition")){
            histSaver->SaveProjectionX((TH2*)it->second);
            histSaver->SaveBinnedProjectionX((TH2*)it->second);
        }
        if (it->first.Contains("RelPosition")){
            if (className.Contains("TProfile2D")){
                histSaver->SaveHistogram((TProfile2D*)it->second);
                histSaver->SaveOverlay((TProfile2D*)it->second);

            }
            else{
                histSaver->SaveOverlay((TH2*)it->second);
            }
        }
        if( it->first == "hNegativeChargeAboveCut_Position"||
                it->first == "hOnlyTranspClusterPosition" ||
                it->first == "hOnlyClusteredClusterPosition") {
            name = "c";
            name.Append(it->first);
            name.Append(extension);
            histSaver->SaveHistogramWithCellGrid((TH2*)it->second);
            if (it->second->GetEntries()){
                TCanvas *c1 = new TCanvas(name,name);
                c1->cd();
                it->second->Draw("colz");
                settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
                settings->DrawMetallisationGrid(c1,3);
                histSaver->SaveCanvas(c1);
                delete c1;
            }
        }
        if (it->first == "hNegativeChargeChannelPositionTransparent"){
            histSaver->SaveProjectionY((TH2*)it->second);
        }
        if (it->first == "hNoNegativeChargeLowResponsePosition" ||
            it->first == "hNoNegativeChargeResponseWindowPosition" ||
            it->first == "hChargeResponseWindowPosition"){
            name = "c";
            name.Append(it->first);
            name.Append(extension);
            histSaver->SaveHistogramWithCellGridAndMarkedCells((TH2*)it->second);
            TCanvas *c1 = histSaver->DrawHistogramWithCellGrid((TH2*)it->second);
            c1->cd();
//            it->second->Draw("colz");
            settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
            settings->DrawMetallisationGrid(c1,3);
            histSaver->AddMarkedCells(c1);
            histSaver->SaveCanvas(c1);
            delete c1;
        }
        delete it->second;
    }
}


void TAnalysisOfAnalysisDifferences::InitTransparentHistos() {
    TH1* histo;
    TString name = "hOnlyTranspClusterCharge";
    histo = new TH1F(name+extension,name+extension,256,0,2048);
    mapHistos[name] = histo;

    name = "hOnlyTranspClusterPosition";
    histo = histSaver->GetHistoBinedInCells(name+extension,4);
    mapHistos[name] = histo;

    name = "hNegativeChargePositionTransparent";
    histo = new TH2F(name+extension,name+extension,512,-512,0,6,-.5,5.5);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    histo->GetYaxis()->SetTitle("position of negative charge in transp. cluster");
    mapHistos[name] = histo;

    name = "hNegativeChargeChannelPositionTransparent";
    histo = new TH2F(name+extension,name+extension,512,-512,0,7,-3.5,3.5);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    histo->GetYaxis()->SetTitle("rel position of negative charge in transp. cluster");
    mapHistos[name] = histo;

}

void TAnalysisOfAnalysisDifferences::InitClusteredHistos() {
    TH1* histo;
    TString name = "hOnlyClusteredClusterCharge";
    TString hname = name +extension;
    histo = new TH1F(hname,hname,256,0,2048);
    mapHistos[name] = histo;

    name = "hOnlyClusteredClusterPosition";
    hname = name+extension;
    histo = histSaver->GetHistoBinedInCells(hname,4);
    mapHistos[name] = histo;

}

void TAnalysisOfAnalysisDifferences::InitSameHistos() {
    TH1* histo;
    TString name = "hChargeDifference";
    TString hname = name +extension;
    histo = new TH1F(hname,hname,1024,-512,512);
    histo->GetXaxis()->SetTitle("charge difference_{trans-clus} / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    mapHistos[name] = histo;

    name = "hNegativeChargePosition";
    hname =name+extension;
    histo = new TH2F(hname,hname,512,-512,0,6,-.5,5.5);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    histo->GetYaxis()->SetTitle("position of negative charge in transp. cluster");
    mapHistos[name] = histo;

    name = "hNegativeCharge";
    hname = name +extension;
    histo = new TH1F(hname,hname,256,-256,0);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    mapHistos[name] = histo;

    name = "hNegativeChargeAboveCut_Position";
    hname = name +extension;
    histo = histSaver->GetHistoBinedInCells(hname,8);
    mapHistos[name] = histo;

    name = "hAllEvents_RelPosition";
    hname = name +extension;
    histo = settings->GetOverlayHisto(hname);
    histo->SetMinimum(0);
    mapHistos[name] = histo;

    name = "hChargeResponseWindowRelPosition";
    hname = name +extension;
    histo = settings->GetOverlayHisto(hname);
    histo->SetTitle(TString::Format("ChargeResponseWindowPosition, Thrs %.0f - %0.f",
            settings->getResponseWindow().first,settings->getResponseWindow().first));
    mapHistos[name] = histo;

    name = "hNegativeChargeAboveCut_RelPosition";
    hname = name +extension;
    histo = settings->GetOverlayHisto(hname);
    mapHistos[name] = histo;

    name = "hNoNegativeCharge_RelPosition";
    hname = name +extension;
    histo = settings->GetOverlayHisto(hname);
    mapHistos[name] = histo;
// mapHistos["hNegativeChargePosition"]

    name = "hNegativeChargeProfileRelPosition";
    hname = name +extension;
    histo = settings->GetOverlayProfile(hname);
    mapHistos[name] = histo;

    name = "hPositive_Minus_Negative_TransparentCharge";
    hname  = name +extension;
    mapHistos[name] = new TH1F(hname,hname,1024,-64,512-64);

    name = "hPositiveTransparentCharge_Minus_ClusteredCharge";
    hname - name +extension;
    mapHistos[name] = new TH1F(hname,hname,512,-128,128);

    name = "hNoNegativeChargeLowResponsePosition";
    hname=name+extension;
    TString title = name;
    title +=TString::Format(" Thrs: %.0f",settings->getLowResponseThreshold())+extension;
    histo = histSaver->GetHistoBinedInCells(hname,8);
    histo->SetTitle(title);
    mapHistos[name] = histo;

    name = "hNoNegativeChargeResponseWindowPosition";
    hname=name+extension;
    title = name;
    title+=TString::Format(", Thrs: %.0f - %.0f",settings->getResponseWindow().first,
                                      settings->getResponseWindow().second)+extension;
    histo = histSaver->GetHistoBinedInCells(hname,8);
    histo->SetTitle(title);
    mapHistos[name] = histo;

    name = "hChargeResponseWindowPosition";
     title = name;
    hname=name+extension;
    title+=TString::Format(" Thrs %.0f - %.0f",settings->getResponseWindow().first,
                                      settings->getResponseWindow().second)+extension;
    histo = histSaver->GetHistoBinedInCells(hname,8);
    histo->SetTitle(title);
    mapHistos[name] = histo;
}

void TAnalysisOfAnalysisDifferences::UpdatePredictedPosition() {

    if (verbosity>2)cout<<"UPDATE Postion"<<endl;
    while (itClustered->first > itPredicted->first &&
            itTransparent->first > itPredicted->first &&
            itPredicted!= predictedPositions->end()){
        Int_t clustered = itClustered->first;;
        Int_t transparent = itTransparent->first;
        Int_t predicted = itPredicted->first;
        if(verbosity>8){
            cout<<"update Predicted Positions ";
            cout<<clustered<<" "<<transparent<<" "<<predicted<<" ";
            cout<< (clustered > predicted) <<" "<< (transparent>predicted) <<" "<< (itPredicted != predictedPositions->end())<<endl;;
        }
        itPredicted++;
    }
}
