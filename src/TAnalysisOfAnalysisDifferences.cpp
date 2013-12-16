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
    verbosity =settings->getVerbosity();
}

TAnalysisOfAnalysisDifferences::~TAnalysisOfAnalysisDifferences() {
    // TODO Auto-generated destructor stub
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
    cout<<"Same Events:      "<<nSameEvents<<endl;
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
    Int_t eventNo = itClustered->first;
    mapHistos["hChargeDifference"]->Fill(transparentCharge-clusteredCharge);
    Int_t pos = 0;
    Float_t charge = 0;
    if (hasNegativeCharge(itTransparent,pos,charge)){
        //        cout<<itTransparent->first<<" "<<pos<<" "<<charge<<endl;
        mapHistos["hNegativeChargePosition"]->Fill(charge,pos);
        mapHistos["hNegativeCharge"]->Fill(charge);
        if(charge<-50){
            if(predictedPositions->count(eventNo)){
                Float_t xPredDet = predictedPositions->at(eventNo).first;
                Float_t yPredDet = predictedPositions->at(eventNo).second;
                mapHistos["hNegativeChargeAbove50_Position"]->Fill(xPredDet,yPredDet);
            }
            else{
                cerr<<"cannot find "<<eventNo<<endl;
            }
        }
        //        cout<<"FILL "<<charge<<" "<<pos<<endl;
    }
    Float_t posCharge = itTransparent->second.getPositiveCharge(true);//100,true,true);
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
}

void TAnalysisOfAnalysisDifferences::AnalyseOnlyTransparentEvent() {
    nOnlyTransparent++;
    if( verbosity>3)
        cout<<itTransparent->first<< " only in Transparent Analysis"<<endl;

    if(itPredicted->first != itTransparent->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<endl;
    mapHistos["hOnlyTranspClusterCharge"]->Fill(itTransparent->second.getPositiveCharge());
    mapHistos["hOnlyTranspClusterPosition"]->Fill(itPredicted->second.first,itPredicted->second.second);

}

void TAnalysisOfAnalysisDifferences::AnalyseOnlyClusteredEvent() {
    if( verbosity>3)
        cout<<itClustered->first<< " only in Clustered Analysis"<<endl;
    nOnlyClustered++;
    if(itPredicted->first != itClustered->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itClustered->first<<endl;
    mapHistos["hOnlyClusteredClusterCharge"]->Fill(itClustered->second.getPositiveCharge());
    mapHistos["hOnlyClusteredClusterPosition"]->Fill(itPredicted->second.first,itPredicted->second.second);

}

void TAnalysisOfAnalysisDifferences::InitHistograms() {
    this->InitSameHistos();
    this->InitTransparentHistos();
    this->InitClusteredHistos();

}

void TAnalysisOfAnalysisDifferences::SaveHistograms() {
    std::map<TString,TH1*>::iterator it = mapHistos.begin();
    for(it;it!=mapHistos.end();it++){
        TString className = it->second->ClassName();
        if( it->first == "hNegativeChargeAbove50_Position"||
                it->first == "hOnlyTranspClusterPosition" ||
                it->first == "hOnlyClusteredClusterPosition"){
            TString name = "c";
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

bool TAnalysisOfAnalysisDifferences::hasNegativeCharge(
        std::map<Int_t, TCluster>::iterator it, Int_t& pos, Float_t& charge) {
    TCluster* clus = &it->second;
    Float_t predPos  = clus->GetTransparentHitPosition();
    if (predPos<0)
        return false;
    Float_t oldCharge = 0;
    Float_t currentCharge;
    for(UInt_t i = 1; i<= clus->size();i++){
        currentCharge = clus->getCharge(i,true);
        charge = currentCharge - oldCharge;
        if (charge < 0){
            if( i+1<clus->size() && clus->getCharge(i+1,true)-currentCharge > charge){
                pos = i;
                //               cout << "found negative charge at "<< pos<<": "<<charge<<endl;
                return true;
            }else{
                pos = i+1;
                charge =  clus->getCharge(i+1,true)-currentCharge;
                //               cout << "found negative charge at "<< pos<<": "<<charge<<endl;
                return true;
            }

        }
        oldCharge = currentCharge;
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
    histo->GetXaxis()->SetTitle("PH Difference_{trans-clus} /ADC");
    histo->GetYaxis()->SetTitle("No of entries #");
    mapHistos[name] = histo;

    name = "hNegativeChargePosition";
    histo = new TH2F(name,name,512,-512,0,6,-.5,5.5);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster /adc");
    histo->GetYaxis()->SetTitle("position of negative charge in transp. cluster");
    mapHistos[name] = histo;

    name = "hNegativeCharge";
    histo = new TH1F(name,name,512,-512,0);
    histo->GetXaxis()->SetTitle("first neg. Charge in transparent Cluster /adc");
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
