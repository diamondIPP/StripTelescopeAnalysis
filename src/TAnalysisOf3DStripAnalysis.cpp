/*
 * TAnalysisOf3DStrip.cpp
 *
 *  Created on: Oct 7, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DStripAnalysis.hh>

TAnalysisOf3DStrip::TAnalysisOf3DStrip(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna) {
    this->histSaver = histSaver;
    this->settings = settings;
    verbosity = settings->getVerbosity();
    subjectPlane = TPlaneProperties::getDiamondPlane();
    subjectDetector = TPlaneProperties::getDetDiamond();
    useCMN = true;

    this->bTransAna = bTransAna;
    if (this->bTransAna)
        appendix ="_trans";
    else
        appendix ="";
    PulseHeightBins = 256;
    PulseHeightMin = 1;
    PulseHeightMax = 2800;
    PulseHeightMinMeanCharge = 1;
    PulseHeightMaxMeanCharge = 1500;
}

TAnalysisOf3DStrip::~TAnalysisOf3DStrip() {
    // TODO Auto-generated destructor stub
    delete hLandauStrip;
    delete hLandauStripFidCutXvsFidCutY;
}


/**
 * Analysis of events with hit in strip detector.
 * Checks if chi2 cut is fullfilled.
 * checks if right selection/area is hitted
 * checks if only one cluster exists
 * checks if is valid Cluster
 * checks that cluster is not saturated
 */
void TAnalysisOf3DStrip::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred,
        Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y) {
    if (verbosity)
        cout<<"TAnalysisOf3DStrip::addEvent"<<endl;
    predx = x_pred;
    predy = y_pred;
    fidx = x_fid;
    fidy = y_fid;
    this->chi2x = chi2x;
    this->chi2y = chi2y;
    if(!settings->do3dTransparentAnalysis()){
        diamondCluster = cluster;
    }
    else
    {
        diamondCluster = transparentCluster;

    }
    Float_t charge = diamondCluster->getPositiveCharge(false);
    fillPlots();
}

void TAnalysisOf3DStrip::initHistos() {
    TString name = "hLandauStrip";
    name.Append(appendix);
    hLandauStrip =  new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauStrip->GetXaxis()->SetTitle("charge / ADC");
    hLandauStrip->GetYaxis()->SetTitle("number of entries #");
    hLandauStrip->SetLineColor(kBlue);

    name = "hLandauStripFidCutXvsFidCutY";
    hLandauStripFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandauStripFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandauStripFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");
    hLandauStripFidCutXvsFidCutY->GetZaxis()->SetTitle("Charge ADC");
}

void TAnalysisOf3DStrip::saveHistos() {
    histSaver->SaveHistogram(hLandauStrip);
    histSaver->SaveHistogram(hLandauStripFidCutXvsFidCutY);
}

void TAnalysisOf3DStrip::fillPlots() {
    if (verbosity>5)
        cout<<"TAnalysisOf3DStrip::fillPlots"<<flush;
    if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())
        return;
    Int_t stripDetector = 1;
    if (verbosity>7)cout<<"1 "<<flush;
    if(settings->getSelectionFidCuts()->getFidCutRegion(fidx,fidy)!=stripDetector)
        return;
    if (verbosity>7)cout<<"2 "<<flush;
    if(eventReader->getNDiamondClusters()!=1)
        return;
    if (verbosity>7)cout<<"3 "<<flush;
    int areaStripDetector = 0;
    if (!settings->isClusterInDiaDetectorArea(diamondCluster,areaStripDetector) ){
        if (verbosity>7)cout<<"3.1"<<flush;
        return;
    }
    if (verbosity>7) cout<<"3.2"<<flush;
    if( !settings->diamondPattern.isValidCluster(diamondCluster)){
        if (verbosity>7)cout<<"4.1"<<flush;
        return;
    }
    if (verbosity>7)cout<<"4"<<flush;
    if (diamondCluster->isSaturatedCluster())
        return;
    if (verbosity>7)cout<<"4"<<flush;
    hLandauStrip->Fill(diamondCluster->getPositiveCharge());
    hLandauStripFidCutXvsFidCutY->Fill(fidx, fidy);
    if (verbosity>5) cout<<"done"<<flush;
}
