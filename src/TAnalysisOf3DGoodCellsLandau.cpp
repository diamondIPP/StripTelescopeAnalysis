/*
 * TAnalysisOf3DGoodCellsLandau.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DGoodCellsLandau.hh>

TAnalysisOf3DGoodCellsLandau::TAnalysisOf3DGoodCellsLandau(TSettings* settings,
        HistogrammSaver* histSaver, bool bTransAna) {
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

TAnalysisOf3DGoodCellsLandau::~TAnalysisOf3DGoodCellsLandau() {
    // TODO Auto-generated destructor stub
}

void TAnalysisOf3DGoodCellsLandau::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred,
        Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y) {
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
    FillGoodCellsLandaus(charge);
}

void TAnalysisOf3DGoodCellsLandau::initHistos() {
    mapPredictedPositionsGoodCells.clear();
//    mapClusteredAnalysisGoodCells.clear();
    mapTransparentAnalysisGoodCells.clear();
//    mapPredictedPositionsAllCells.clear();
//    mapClusteredAnalysisAllCells.clear();
//    mapTransparentAnalysisAllCells.clear();
    TString name = "hLandauGoodCells";
    name.Append(appendix);
    hLandauGoodCells = new TH1F(name,name,PulseHeightBins, PulseHeightMin,PulseHeightMax);
    hLandauGoodCells->GetXaxis()->SetTitle("PulseHeight [ADC]");
    hLandauGoodCells->GetYaxis()->SetTitle("number fo entries #");

    name = "hLandauGoodCellsWithoutEdges";
    name.Append(appendix);
    hLandauGoodCellsWithoutEdges = new TH1F(name,name,PulseHeightBins, PulseHeightMin,PulseHeightMax);
    hLandauGoodCellsWithoutEdges->GetXaxis()->SetTitle("PulseHeight [ADC]");
    hLandauGoodCellsWithoutEdges->GetYaxis()->SetTitle("number fo entries #");
    hLandauGoodCellsWithoutEdges->SetLineColor(kGreen);

    name = "hLandauGoodCellsWithoutColumns";
    name.Append(appendix);
    hLandauGoodCellsWithoutColumns= new TH1F(name,name,PulseHeightBins, PulseHeightMin,PulseHeightMax);
    hLandauGoodCellsWithoutColumns->GetXaxis()->SetTitle("PulseHeight [ADC]");
    hLandauGoodCellsWithoutColumns->GetYaxis()->SetTitle("number fo entries #");
    hLandauGoodCellsWithoutColumns->SetLineColor(kBlue);

    //hPulseHeightVsDetectorHitPostionXYGoodCells
    name = "hPulseHeightVsDetectorHitPostionXYGoodCells";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    Float_t factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXYGoodCells = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetXaxis()->SetTitle("#it{x} / #mum");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetYaxis()->SetTitle("#it{y} / #mum");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetZaxis()->SetTitle("charge /ADC");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
}

void TAnalysisOf3DGoodCellsLandau::saveHistos(TH1F* hLandauStrip,TH1F* hLandauPhantom) {
    cout<<"[ TAnalysisOf3DGoodCellsLandau::saveHistos]"<<endl;
    if (hLandauGoodCellsWithoutEdges){
        histSaver->SaveHistogramLandau(hLandauGoodCellsWithoutEdges);
    }
    if (hLandauGoodCellsWithoutColumns){
        histSaver->SaveHistogramLandau(hLandauGoodCellsWithoutColumns);
    }
    if (hLandauGoodCells)
        histSaver->SaveHistogramLandau(hLandauGoodCells);
    if (hLandauGoodCellsWithoutEdges && hLandauGoodCellsWithoutColumns && hLandauGoodCells){
        TCanvas *c1 = new TCanvas("cLandauGoodCellsAll");
        c1->cd();
        hLandauGoodCells->DrawNormalized();
        hLandauGoodCellsWithoutColumns->DrawNormalized("same");
        hLandauGoodCellsWithoutEdges->DrawNormalized("same");
        histSaver->SaveCanvas(c1);
        delete c1;
    }
    if(hLandauGoodCellsWithoutEdges) delete hLandauGoodCellsWithoutEdges;
    if(hLandauGoodCellsWithoutColumns) delete hLandauGoodCellsWithoutColumns;
    if(hLandauGoodCells) delete hLandauGoodCells;
    cout<<"[ TAnalysisOf3DGoodCellsLandau::saveHistos] DONE 1"<<endl;
}

void TAnalysisOf3DGoodCellsLandau::FillGoodCellsLandaus(Float_t charge) {
    if(settings->get3dMetallisationFidCuts()->getFidCutRegion(predx,predy)!=3)
        return;
    mapPredictedPositionsGoodCells[nEvent] = make_pair(predx,predy);
    if(validClusteredAnalysis)
        mapClusteredAnalysisGoodCells[nEvent] = *clusteredCluster;
    if (validTransparentAnalysis) //TODO
        mapTransparentAnalysisGoodCells[nEvent] = *transparentCluster;

    hPulseHeightVsDetectorHitPostionXYGoodCells->Fill(predx,predy,charge);
    hLandauGoodCells->Fill(charge);
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(predx,predy);
    bool isInEdgeRegion =  settings->IsOnTheEdgeOfCell(relPos.first,relPos.second);
    bool isInColumnRegion = settings->IsWithInTheColumnRadius(relPos.first,relPos.second);
    if (!isInEdgeRegion)
        hLandauGoodCellsWithoutEdges->Fill(charge);
    if(!isInColumnRegion)
        hLandauGoodCellsWithoutColumns->Fill(charge);
}
