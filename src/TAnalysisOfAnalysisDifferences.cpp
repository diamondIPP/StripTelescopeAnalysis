/*
 * TAnalysisOfAnalysisDifferences.cpp
 *
 *  Created on: Dec 13, 2013
 *      Author: bachmair
 */

#include "../include/TAnalysisOfAnalysisDifferences.hh"

TAnalysisOfAnalysisDifferences::TAnalysisOfAnalysisDifferences(TSettings* settings,HistogrammSaver* histSaver,TString extension) {
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
    negChargeCut = settings->getNegativeChargeCut();
    this->extension = "_"+extension;
    hPhantomLandau = 0;
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


    if(itPredicted->first != itTransparent->first || itPredicted->first != itClustered->first)
        cout<<"predicted,transparent and clustered iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<" "<<itClustered->first<<endl;
    Float_t clusteredCharge = itClustered->second.getCharge(true);
    itTransparent->second.SetTransparentClusterSize(3);
    Float_t transparentCharge = itTransparent->second.getCharge(true);
    if (itClustered->second.isSaturatedCluster()||itTransparent->second.isSaturatedCluster())
        continue;
    Float_t posCharge = itTransparent->second.getPositiveCharge(true);//100,true,true);
    mapHistos["hTransparentPulseHeight"]->Fill(posCharge);
    mapHistos["hClusteredPulseHeight"]->Fill(clusteredCharge);
    if (1810<posCharge&&posCharge<1840){
        cout<<"PH: "<<posCharge;
        itTransparent->second.Print(1);
    }

    Int_t eventNo = itClustered->first;
    mapHistos["hChargeDifference"]->Fill(transparentCharge-clusteredCharge);
    Int_t pos = 0;
    Float_t charge = 0;
    bool hasNegativeCharge = itTransparent->second.hasNegativeCharge(charge,pos,true);
//    cout<<eventNo<< " negCharge: "<<hasNegativeCharge<<" "<<charge<< " "<<pos<<endl;
    if (hasNegativeCharge){
        mapHistos["hNegativeChargePosition"]->Fill(charge,pos);
		mapHistos["hNegativeChargePositionTransparent"]->Fill(charge,pos);
        mapHistos["hNegativeCharge"]->Fill(charge);
        if(charge<negChargeCut){
            if(predictedPositions->count(eventNo)){
                Float_t xPredDet = predictedPositions->at(eventNo).first;
                Float_t yPredDet = predictedPositions->at(eventNo).second;
                mapHistos["hNegativeChargeAboveCut_Position"]->Fill(xPredDet,yPredDet);
            }
        }
    }
    else
    	mapHistos["hNegativeChargePositionTransparent"]->Fill(charge,-1);
    if(hasNegativeCharge && charge < negChargeCut)
        mapHistos["hHasNegativeCharge_PulseHeight"]->Fill(posCharge);
    else
        mapHistos["hNoNegativeCharge_PulseHeight"]->Fill(posCharge);
    mapHistos["hPositive_Minus_Negative_TransparentCharge"]->Fill(posCharge-transparentCharge);
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

void TAnalysisOfAnalysisDifferences::AnalyseOnlyTransparentEvent() {

    if (itTransparent->second.getClusterSize()==0){
        cout<<"invalid cluster size"<<endl;
        return;
    }

    nOnlyTransparent++;
    if( verbosity>3)
        cout<<itTransparent->first<< " only in Transparent Analysis"<<endl;

    if(itPredicted->first != itTransparent->first)
        cout<<"predicted and clusterer iterator do not agree"<<itPredicted->first<<" "<<itTransparent->first<<endl;
    mapHistos["hOnlyTranspClusterCharge"]->Fill(itTransparent->second.getPositiveCharge());

    mapHistos["hTransparentPulseHeight"]->Fill(itTransparent->second.getPositiveCharge());
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

    TString name = "hTransparentPulseHeight";
    TString hname = name+extension;
    histo = new TH1F(hname,hname,bins,xmin,xmax);
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
}

void TAnalysisOfAnalysisDifferences::SaveHistograms() {
    TH1F* histo1 = (TH1F*)mapHistos["hClusteredPulseHeight"]->Clone();
    TH1F* histo2 = (TH1F*)mapHistos["hTransparentPulseHeight"]->Clone();
    TH1F* hPhantomScaled = 0;
    if (hPhantomLandau){
        hPhantomScaled = (TH1F*)hPhantomLandau->Clone();
        Float_t scale = hPhantomScaled->GetBinContent(hPhantomScaled->GetMaximumBin());
        hPhantomScaled->Scale(1./scale);
    }
    histSaver->SaveTwoHistos("cComparisionPulseHeights"+extension,histo1,histo2);
    histSaver->SaveTwoHistosNormalized("cComparisionPulseHeightsNormalized"+extension,histo1,histo2);

    TString name = "cComparisionPulseHeightWithStrip"+extension;
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

    name = "cComparisionPulseHeightWithStrip_scaled"+extension;
    c1->SetName(name);
    TH1F* h1 = (TH1F*)histo1->Clone();
    TH1F* h2 = (TH1F*)histo2->Clone();
    h1->Scale(1./h1->GetBinContent(h1->GetMaximumBin()));
    h2->Scale(1./h2->GetBinContent(h2->GetMaximumBin()));
    h1->Draw();
    h2->Draw("same");
    TH1F* hS= 0;
    if(stripHisto){
        hS = (TH1F*)stripHisto->Clone();
        hS->Scale(1./hS->GetBinContent(hS->GetMaximumBin()));
        hS->SetLineColor(kBlue);
        hS->Draw("same");
        hS->SetTitle("Strip");
    }
    h1->SetTitle("Clustered");
    h2->SetTitle("Tranparent");
    if (leg) delete leg;
    leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->Draw();
    histSaver->SaveCanvas(c1);
    if (h2) delete h1;
    if (h1) delete h2;
    if (hS)
        delete hS;
    cout<<"Save6"<<endl;


    name = "cNoNegativeCharge_PulseHeight_Comparision"+extension;
    c1->SetTitle(name);
    c1->SetName(name);
    c1->Clear();
    mapHistos["hNoNegativeCharge_PulseHeight"]->SetTitle("3D");
    mapHistos["hNoNegativeCharge_PulseHeight"]->Draw();
//    if(histo2){
//        histo2->SetLineColor(kRed);
//        histo2->Draw("same");
//    }
    if(stripHisto){
        stripHisto->SetTitle("Strip");
        stripHisto->SetLineStyle(2);
        stripHisto->Draw("same");
    }
    leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->Draw();
    histSaver->SaveCanvas(c1);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision"+extension;
    THStack *hs1 = new THStack(name,name);
    hs1->Add(mapHistos["hNoNegativeCharge_PulseHeight"]);
    if(stripHisto)hs1->Add(stripHisto);
    hs1->Draw("");
    if(hs1->GetXaxis())hs1->GetXaxis()->SetTitle("charge / ADC");
    if(hs1->GetYaxis())hs1->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs1,"nostack",true,false);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision_normalized"+extension;
    THStack *hs2 = new THStack(name,name);
    hs2->Add(mapHistos["hNoNegativeCharge_PulseHeight"]->DrawNormalized("goff"));
    if(stripHisto)hs2->Add(stripHisto->DrawNormalized("goff"));
    hs2->Draw("");
    if(hs2->GetXaxis())hs2->GetXaxis()->SetTitle("charge / ADC");
    if(hs2->GetYaxis())hs2->GetYaxis()->SetTitle("number of entries - normalized");
//    histSaver->SaveStack(hs2,"nostack",true,false);

    name = "stack_NoNegativeCharge_PulseHeight_Comparision_scaled"+extension;
    THStack *hs3 = new THStack(name,name);
    TH1F* histo3 = (TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone();
    Float_t scale = histo3->GetBinContent(histo3->GetMaximumBin());
    histo3->Scale(1./scale);
    histo3->SetTitle("3D");
    TH1F* histo4 =0;
    if (stripHisto){
        histo4 = (TH1F*)stripHisto->Clone();
        scale = histo4->GetBinContent(histo4->GetMaximumBin());
        histo4->Scale(1./scale);
        histo4->SetTitle("Strip");
    }
    hs3->Add(histo3);
    if(histo4)
        hs3->Add(histo4);
    if (hPhantomScaled){
        hs3->Add(hPhantomScaled);
    }
    hs3->Draw("");
    if(hs3->GetXaxis())hs3->GetXaxis()->SetTitle("charge / ADC");
    if(hs3->GetYaxis())hs3->GetYaxis()->SetTitle("number of entries - rescaled");
    histSaver->SaveStack(hs3,"nostack",true,false,"charge / ADC","rel. no of entries");



    name = "stackPulseHeights_GoodCells"+extension;
    THStack *hs = new THStack(name,name);
    hs->Add((TH1F*)mapHistos["hHasNegativeCharge_PulseHeight"]->Clone());
    hs->Add((TH1F*)mapHistos["hNoNegativeCharge_PulseHeight"]->Clone());
    if (hPhantomScaled){
        hs3->Add(hPhantomScaled);
    }
    hs->Draw("");
    if(hs->GetXaxis())hs->GetXaxis()->SetTitle("charge / ADC");
    if(hs->GetYaxis())hs->GetYaxis()->SetTitle("number of entries");
    histSaver->SaveStack(hs,"",true,false);

    name = "stackPulseHeights_GoodCells"+extension;
    hs->SetTitle(name);
    hs->SetName(name);
    histSaver->SaveStack(hs,"",true,false);

    name = "stackPulseHeights_GoodCells_nostack"+extension;
    hs->SetTitle(name);
    hs->SetName(name);
    histSaver->SaveStack(hs,"nostack",true,false);

    delete hs;
    delete hs1;
    delete hs2;
    delete histo1;
    delete histo2;

    std::map<TString,TH1*>::iterator it = mapHistos.begin();
    cout<<"[TAnalysisOfAnalysisDifferences]Save Histos it map "<<mapHistos.size()<<endl;
    for(it;it!=mapHistos.end();it++){
        TString className = it->second->ClassName();
        if( it->first == "hNegativeChargeAboveCut_Position"||
                it->first == "hOnlyTranspClusterPosition" ||
                it->first == "hOnlyClusteredClusterPosition"){
            name = "c";
            name.Append(it->first);
            name.Append(extension);
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

    name = "hPositive_Minus_Negative_TransparentCharge";
    hname  = name +extension;
    mapHistos[name] = new TH1F(hname,hname,1024,-64,512-64);

    name = "hPositiveTransparentCharge_Minus_ClusteredCharge";
    hname - name +extension;
    mapHistos[name] = new TH1F(hname,hname,512,-128,128);

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
