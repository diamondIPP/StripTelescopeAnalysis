/*
 * TAnalysisOf3DResolutionStudies.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DResolutionStudies.hh>

TAnalysisOf3DResolutionStudies::TAnalysisOf3DResolutionStudies(TSettings* settings,
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
}

TAnalysisOf3DResolutionStudies::~TAnalysisOf3DResolutionStudies() {
    // TODO Auto-generated destructor stub
}

void TAnalysisOf3DResolutionStudies::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred,
        Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y) {
    predx = x_pred;
    predy = y_pred;
    fidx = x_fid;
    fidy = y_fid;
    this->chi2x = chi2x;
    this->chi2y = chi2y;
    if(!settings->do3dTransparentAnalysis()){
//        cout<<"cluster"<<flush;
        diamondCluster = cluster;
    }
    else
    {
//        cout<<"transparentcluster"<<flush;
        diamondCluster = transparentCluster;

    }
//    cout<<" "<<cluster<<endl;
    fillResolutionPlots();
}

void TAnalysisOf3DResolutionStudies::initHistos() {
    UInt_t nCells = 99;
    UInt_t nBins = 128;
    Float_t minX = - 1*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX = 1*settings->GetCellWidth(subjectDetector,2);
    for (UInt_t cell = 0; cell <nCells;cell++){
        TString name = TString::Format("hResolution_CellNo_%02d_maxValue",cell);
        TString title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        TH1F* histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_maxValue.push_back(histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_chargeWeighted.push_back(histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell);
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_highest2Centroid.push_back(histo);
    }
}

void TAnalysisOf3DResolutionStudies::saveHistos() {
    cout<<"TAnalysisOf3DResolutionStudies::saveHistos()"<<endl;
    CreateResolutionPlots(&vecHResolutionPerCell_chargeWeighted,"chargeWeighted");
    CreateResolutionPlots(&vecHResolutionPerCell_highest2Centroid,"highest2Centroid");
    CreateResolutionPlots(&vecHResolutionPerCell_maxValue,"maxValue");
    cout<<"TAnalysisOf3DResolutionStudies::saveHistos() DONE"<<flush;
}

void TAnalysisOf3DResolutionStudies::setEventReader(TTracking* eventReader) {
    this->eventReader=eventReader;
}

void TAnalysisOf3DResolutionStudies::fillResolutionPlots() {
    UInt_t cellNo = settings->getCellNo(predx,predy);
    Float_t cellWidth = settings->GetCellWidth(subjectDetector,2);
    diamondCluster->SetTransparentClusterSize(3);
    Float_t predPos = diamondCluster->GetTransparentHitPosition();
    Float_t pos = diamondCluster->getPosition(useCMN,TCluster::maxValue);
    Float_t delta = pos - predPos;
    if (cellNo< vecHResolutionPerCell_maxValue.size()){
        TH1F* histo  = vecHResolutionPerCell_maxValue.at(cellNo);
        if (histo);
        histo->Fill(delta*cellWidth);
    }
    /*********/
    pos = diamondCluster->getPosition(useCMN,TCluster::chargeWeighted);
    delta = pos - predPos;
    if (cellNo< vecHResolutionPerCell_chargeWeighted.size()){
        TH1F* histo  = vecHResolutionPerCell_chargeWeighted.at(cellNo);
        if (histo);
        histo->Fill(delta*cellWidth);
    }
    /*********/
    pos = diamondCluster->getPosition(useCMN,TCluster::highest2Centroid);
    delta = pos - predPos;
    if (cellNo< vecHResolutionPerCell_highest2Centroid.size()){
        TH1F* histo  = vecHResolutionPerCell_highest2Centroid.at(cellNo);
        if (histo);
        histo->Fill(delta*cellWidth);
    }
    if (cellNo<99&&false   ){

        cout<<"POS: "<<predPos<<" / "<<pos<<endl;
        diamondCluster->Print();
    }
}

void TAnalysisOf3DResolutionStudies::CreateResolutionPlots(vector<TH1F*>* vec, TString kind) {
    UInt_t nBins = 128;
    Float_t minX = -1*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX =settings->GetCellWidth(subjectDetector,2);
    TString name = "hResolutionGoodCells_"+kind;
    TH1F* hResolutionGoodCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionGoodCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionBadCells_"+kind;
    TH1F* hResolutionBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionBadCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllCells_"+kind;
    TH1F* hResolutionAllCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllButBadCells_"+kind;
    TH1F* hResolutionAllButBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllButBadCells->GetXaxis()->SetTitle("Residual / #mum");
    string plots_path = histSaver->GetPlotsPath();
    histSaver->SetPlotsPath(plots_path+(string)"/resolution/");
    for(UInt_t cell=0;cell< vec->size();cell++){
        TH1F* histo = vec->at(cell);
        if (!histo)
            continue;
        if (settings->IsGoodCell(3,cell))
            hResolutionGoodCells->Add(histo);
        if (settings->isBadCell(3,cell))
            hResolutionBadCells->Add(histo);
        else
            hResolutionAllButBadCells->Add(histo);
        hResolutionAllCells->Add(histo);
        histSaver->SaveHistogram(histo);
        vec->at(cell)= 0;
        delete histo;
    }
    histSaver->SetPlotsPath(plots_path);
    hResolutionGoodCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionGoodCells,false,false,false);
    hResolutionBadCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionBadCells);
    hResolutionAllCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionAllCells);
    hResolutionAllButBadCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionAllButBadCells);
    delete hResolutionGoodCells;
    delete hResolutionBadCells;
    delete hResolutionAllCells;
    delete hResolutionAllButBadCells;
}
