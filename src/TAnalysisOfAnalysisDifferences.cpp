/*
 * TAnalysisOfAnalysisDifferences.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: bachmair
 */

#include "../include/TAnalysisOfAnalysisDifferences.hh"

TAnalysisOfAnalysisDifferences::TAnalysisOfAnalysisDifferences(TSettings* settings,HistogrammSaver* histSaver) {
    // TODO Auto-generated constructor stub
    this->settings = settings;
    this->histSaver = histSaver;
    if(settings==0 || histSaver==0)
        cerr<<"ERROR: invalid settings or histogram saver"<<endl;
    transparentMap = 0;
    clusteredMap = 0;
    predictedPositions = 0;
    stripHisto = 0;
    verbosity =settings->getVerbosity();
    negChargeCut = -50;
}

TAnalysisOfAnalysisDifferences::~TAnalysisOfAnalysisDifferences() {
//    delete stripHisto;
    // TODO Auto-generated destructor stub
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
    cout<<"nexnt"<<endl;
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
    cout<<"loop over events"<<endl;
    cout<<itPredicted->first<<endl;
    while (!endLoop){
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
    if(verbosity%2==1){
        char t;
        cout<<"Press a key and enter"<<endl;
        cin>>t;
    }
}

void TAnalysisOfAnalysisDifferences::AnalyseSameEvent() {
    if(itPredicted->first != itTransparent->first || itPredicted->first != itClustered->first)
        cout<<"predicted,transparent and clustered iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<" "<<itClustered->first<<endl;
    Float_t clusteredCharge = itClustered->second.getCharge(true);
    itTransparent->second.SetTransparentClusterSize(3);
    Float_t transparentCharge = itTransparent->second.getCharge(true);
    Float_t posCharge = itTransparent->second.getPositiveCharge(true);//100,true,true);
    mapHistos["hTransparentPulseHeight"]->Fill(posCharge);
    mapHistos["hClusteredPulseHeight"]->Fill(clusteredCharge);
    Int_t eventNo = itClustered->first;
    mapHistos["hChargeDifference"]->Fill(transparentCharge-clusteredCharge);
    Int_t pos = 0;
    Float_t charge = 0;
    bool hasNegativeCharge = itTransparent->second.hasNegativeCharge(charge,pos,true);
    if (hasNegativeCharge){
        mapHistos["hNegativeChargePosition"]->Fill(charge,pos);
        mapHistos["hNegativeCharge"]->Fill(charge);
        if(charge<negChargeCut){
            if(predictedPositions->count(eventNo)){
                Float_t xPredDet = predictedPositions->at(eventNo).first;
                Float_t yPredDet = predictedPositions->at(eventNo).second;
                mapHistos["hNegativeChargeAbove50_Position"]->Fill(xPredDet,yPredDet);
            }
        }
    }
    if(hasNegativeCharge && charge < negChargeCut)
        mapHistos["hHasNegativeCharge_PulseHeight"]->Fill(posCharge);
    else
        mapHistos["hNoNegativeCharge_PulseHeight"]->Fill(posCharge);
    mapHistos["hPositive_Minus_Negative_TransparentCharge"]->Fill(posCharge-transparentCharge);
    if ( posCharge - clusteredCharge < 0 ){
        cout<<eventNo<<"\t"<<posCharge-clusteredCharge<<"\n\tposCharge: "<<posCharge<<"\n\tclusCharge: "<<clusteredCharge<<endl;
        cout<<"\tclustered: ";
        itClustered->second.Print(2);
        cout<<"\ttransparent: ";
        itTransparent->second.Print(2);
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

void TAnalysisOfAnalysisDifferences::AnalyseOnlyTransparentEvent() {
    nOnlyTransparent++;
    if( verbosity>3)
        cout<<itTransparent->first<< " only in Transparent Analysis"<<endl;

    if(itPredicted->first != itTransparent->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<endl;
    mapHistos["hOnlyTranspClusterCharge"]->Fill(itTransparent->second.getPositiveCharge());

    mapHistos["hTransparentPulseHeight"]->Fill(itTransparent->second.getPositiveCharge());
    mapHistos["hOnlyTranspClusterPosition"]->Fill(itPredicted->second.first,itPredicted->second.second);
    mapHistos["hGoodCellsEventTypes"]->Fill(1);
}

void TAnalysisOfAnalysisDifferences::AnalyseOnlyClusteredEvent() {
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

    TString name = "hTransparentPulseHeight";
    histo = new TH1F(name,name,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("Charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    mapHistos[name] = histo;

    name = "hClusteredPulseHeight";
    histo = (TH1F*) histo->Clone(name);
    histo->SetLineColor(kBlack);
    mapHistos[name] = histo;

    name = "hGoodCellsEventTypes";
    histo = new TH1F(name,name,3,-1.5,1.5);
    histo->GetXaxis()->SetBinLabel(1,"Only Clustered");
    histo->GetXaxis()->SetBinLabel(2,"Same ");
    histo->GetXaxis()->SetBinLabel(3,"Only transparent");
    mapHistos[name] = histo;

    name ="hHasNegativeCharge_PulseHeight";
    histo = new TH1F(name,name,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    histo->SetLineColor(kGreen);
    mapHistos[name] = histo;

    name ="hNoNegativeCharge_PulseHeight";
    histo = new TH1F(name,name,bins,xmin,xmax);
    histo->GetXaxis()->SetTitle("charge / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    histo->SetLineColor(kBlack);
    mapHistos[name] = histo;
}

void TAnalysisOfAnalysisDifferences::SaveHistograms() {
    TH1F* histo1 = (TH1F*)mapHistos["hClusteredPulseHeight"]->Clone();
    TH1F* histo2 = (TH1F*)mapHistos["hTransparentPulseHeight"]->Clone();
    histSaver->SaveTwoHistos("cComparisionPulseHeights",histo1,histo2);
    histSaver->SaveTwoHistosNormalized("cComparisionPulseHeightsNormalized",histo1,histo2);
    TString name = "cComparisionPulseHeight";
    TCanvas *c1 = new TCanvas(name,name);
    histo1->Draw();
    histo2->Draw("same");
    if(stripHisto){
        Float_t scale = histo2->GetBinContent(histo2->GetMaximumBin());
        scale /= stripHisto->GetBinContent(stripHisto->GetMaximumBin());
        stripHisto->Scale(scale);
        stripHisto->SetLineColor(kBlue);
        stripHisto->Draw("same");
        stripHisto->SetTitle("Strip");
    }
    histo1->SetTitle("Clustered");
    histo2->SetTitle("Tranparent");
    TLegend* leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->Draw();
    histSaver->SaveCanvas(c1);
    cout<<"Save6"<<endl;


    name = "cNoNegativeCharge_PulseHeight_Comparision";
    c1->SetTitle(name);
    c1->SetName(name);
    c1->Clear();
    mapHistos["hNoNegativeCharge_PulseHeight"]->Draw();
//    if(histo2){
//        histo2->SetLineColor(kRed);
//        histo2->Draw("same");
//    }
    if(stripHisto)
        stripHisto->Draw("same");
    leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->Draw();
    histSaver->SaveCanvas(c1);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision";
    THStack *hs1 = new THStack(name,name);
    hs1->Add(mapHistos["hNoNegativeCharge_PulseHeight"]);
    if(stripHisto)hs1->Add(stripHisto);
    hs1->Draw("");
    if(hs1->GetXaxis())hs1->GetXaxis()->SetTitle("charge / ADC");
    if(hs1->GetYaxis())hs1->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs1,"nostack",true,false);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision_normalized";
    THStack *hs2 = new THStack(name,name);
    hs2->Add(mapHistos["hNoNegativeCharge_PulseHeight"]->DrawNormalized("goff"));
    if(stripHisto)hs2->Add(stripHisto->DrawNormalized("goff"));
    hs2->Draw("");
    if(hs2->GetXaxis())hs2->GetXaxis()->SetTitle("charge / ADC");
    if(hs2->GetYaxis())hs2->GetYaxis()->SetTitle("number of entries - normalized");
    histSaver->SaveStack(hs2,"nostack",true,false);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision_scaled";
    THStack *hs3 = new THStack(name,name);
    TH1F* histo3 = (TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone();
    Float_t scale = histo3->GetBinContent(histo3->GetMaximumBin());
    histo3->Scale(1./scale);
    histo3->SetTitle("Pulse Height - 3D");
    TH1F* histo4 =0;
    if (stripHisto){
        histo4 = (TH1F*)stripHisto->Clone();
        scale = histo4->GetBinContent(histo4->GetMaximumBin());
        histo4->Scale(1./scale);
        histo4->SetTitle("Pulse Height - Strip");
    }

    hs3->Add(histo3);
    if(histo4)
        hs3->Add(histo4);
    hs3->Draw("");
    if(hs3->GetXaxis())hs3->GetXaxis()->SetTitle("charge / ADC");
    if(hs3->GetYaxis())hs3->GetYaxis()->SetTitle("number of entries - rescaled");
    histSaver->SaveStack(hs3,"nostack",true,false,"charge / ADC","rel. no of entries");





    name = "stackPulseHeights_GoodCells";
    THStack *hs = new THStack(name,name);
    hs->Add((TH1F*)mapHistos["hHasNegativeCharge_PulseHeight"]->Clone());
    hs->Add((TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone());
    hs->Draw("");
    if(hs->GetXaxis())hs->GetXaxis()->SetTitle("charge / ADC");
    if(hs->GetYaxis())hs->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs,"",true,false);

    name = "stackPulseHeights_GoodCells";
    hs->SetTitle(name);
    hs->SetName(name);
    histSaver->SaveStack(hs,"",true,false);

    name = "stackPulseHeights_GoodCells_nostack";
    hs->SetTitle(name);
    hs->SetName(name);
    histSaver->SaveStack(hs,"nostack",true,false);

    delete hs;
    delete hs1;
    delete hs2;
    delete histo1;
    delete histo2;

    std::map<TString,TH1*>::iterator it = mapHistos.begin();
    for(it;it!=mapHistos.end();it++){
        TString className = it->second->ClassName();
        if( it->first == "hNegativeChargeAbove50_Position"||
                it->first == "hOnlyTranspClusterPosition" ||
                it->first == "hOnlyClusteredClusterPosition"){
            name = "c";
            name.Append(it->first);
            histSaver->SaveHistogramWithCellGrid((TH2*)it->second);
            TCanvas *c1 = new TCanvas(name,name);
            c1->cd();
            it->second->Draw("colz");
            settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
            settings->DrawMetallisationGrid(c1,3);
            histSaver->SaveCanvas(c1);
            delete c1;
        }
        if (className.Contains("TH2F"))
            histSaver->SaveHistogram((TH2F*)it->second);
        else
            histSaver->SaveHistogram(it->second);
        delete it->second;
    }
}


void TAnalysisOfAnalysisDifferences::InitTransparentHistos() {
    TH1* histo;
    TString name = "hOnlyTranspClusterCharge";
    histo = new TH1F(name,name,256,0,2048);
    mapHistos[name] = histo;

    name = "hOnlyTranspClusterPosition";
    histo = histSaver->GetHistoBinedInCells(name,4);
    mapHistos[name] = histo;

}

void TAnalysisOfAnalysisDifferences::InitClusteredHistos() {
    TH1* histo;
    TString name = "hOnlyClusteredClusterCharge";
    histo = new TH1F(name,name,256,0,2048);
    mapHistos[name] = histo;

    name = "hOnlyClusteredClusterPosition";
    histo = histSaver->GetHistoBinedInCells(name,4);
    mapHistos[name] = histo;

}

void TAnalysisOfAnalysisDifferences::InitSameHistos() {
    TH1* histo;
    TString name = "hChargeDifference";
    histo = new TH1F(name,name,1024,-512,512);
    histo->GetXaxis()->SetTitle("charge difference_{trans-clus} / ADC");
    histo->GetYaxis()->SetTitle("no of entries #");
    mapHistos[name] = histo;

    name = "hNegativeChargePosition";
    histo = new TH2F(name,name,512,-512,0,6,-.5,5.5);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    histo->GetYaxis()->SetTitle("position of negative charge in transp. cluster");
    mapHistos[name] = histo;

    name = "hNegativeCharge";
    histo = new TH1F(name,name,256,-256,0);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster / ADC");
    mapHistos[name] = histo;

    name = "hNegativeChargeAbove50_Position";
    histo = histSaver->GetHistoBinedInCells(name,4);
    mapHistos[name] = histo;

    name = "hPositive_Minus_Negative_TransparentCharge";
    mapHistos[name] = new TH1F(name,name,1024,-64,512-64);

    name = "hPositiveTransparentCharge_Minus_ClusteredCharge";
    mapHistos[name] = new TH1F(name,name,512,-128,128);

}

void TAnalysisOfAnalysisDifferences::UpdatePredictedPosition() {

    while (itClustered->first > itPredicted->first && itTransparent->first > itPredicted->first && itPredicted!= predictedPositions->end()){
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
