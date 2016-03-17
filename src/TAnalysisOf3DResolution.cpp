/*
 * TAnalysisOf3DResolution.cpp
 *
 *  Created on: Mar 7, 2016
 *      Author: bachmair
 */

#include <TAnalysisOf3DResolution.hh>

TAnalysisOf3DResolution::TAnalysisOf3DResolution(TSettings *settings,  HistogrammSaver *histSaver, TString appendix){
    // TODO Auto-generated constructor stub
    this->settings= settings;
    this->histSaver = histSaver;
    initialiseHistos();
    subjectDetector = TPlaneProperties::getDetDiamond();
    maxsnr = 50;
    useCMN = true;
    this->appendix =appendix;
}

TAnalysisOf3DResolution::~TAnalysisOf3DResolution() {
    // TODO Auto-generated destructor stub
    saveHistos();
    deleteHistos();
}

void TAnalysisOf3DResolution::Fill(TCluster* diamondCluster, Float_t xPredDet, Float_t yPredDet) {

    if (!settings->do3dTransparentAnalysis())
        return;
    UInt_t cellNo = settings->getCellNo(xPredDet,yPredDet);
    Float_t cellWidth = settings->GetCellWidth(subjectDetector,2);
    diamondCluster->SetTransparentClusterSize(5);
    Int_t highest_hit_pos = diamondCluster->getHighestHitClusterPosition();
    Int_t Second_highest_hit_pos = diamondCluster->getHighestSignalNeighbourClusterPosition(highest_hit_pos,useCMN,true);
    Float_t snr = diamondCluster->getSNR(Second_highest_hit_pos,useCMN);
    Float_t predPos = diamondCluster->GetTransparentHitPosition();
    if (snr > maxsnr) snr=maxsnr*109/110;
    if (hAdjacentChannels_SNR)
        hAdjacentChannels_SNR->Fill(diamondCluster->getSNR(highest_hit_pos-1,useCMN),
                diamondCluster->getSNR(highest_hit_pos+1,useCMN));
    if (hAdjacentSNR_vs_cellNo)
        hAdjacentSNR_vs_cellNo->Fill(cellNo,snr);
    if (hAdjacentChannels_Signal)
        hAdjacentChannels_Signal->Fill(diamondCluster->getSignal(highest_hit_pos-1,useCMN),
                diamondCluster->getSignal(highest_hit_pos+1,useCMN));

    Float_t pos_max = diamondCluster->getPosition(useCMN,TCluster::maxValue);
    Float_t delta_max = pos_max - predPos;
    Float_t pos_h2c = diamondCluster->getPosition(useCMN,TCluster::highest2Centroid);
    Float_t delta_h2C = pos_h2c - predPos;
    Float_t pos_weighted = diamondCluster->getPosition(useCMN,TCluster::chargeWeighted);
    Float_t delta_Weigthed = pos_weighted - predPos;
    Float_t relPredPos = fmod(predPos+.5,1)-.5;
    relPredPos*=cellWidth;
    pair<Float_t,Float_t> relPred = settings->getRelativePositionInCell(xPredDet,yPredDet);
    Float_t relPredPosY = relPred.second;
    relPred.first-=settings->GetCellWidth(subjectDetector,2)/2.;
    relPred.second-=settings->GetCellHeight()/2.;
    if (TMath::Abs(relPredPos-relPred.first) > 5)
        cout<<"\n "<<pos_max<<" "<<predPos<<" "<<delta_max<<" "<<delta_Weigthed<<" "<<delta_h2C<< " "<<relPredPos<<" "<<relPred.first<<"/"<<relPred.second<<endl;
    Int_t leftChannel;
    Float_t eta = diamondCluster->getEta(leftChannel,useCMN);
    //cout<<"\nETA: "<<eta<<": leftChannel:"<<leftChannel<<endl;
    //diamondCluster->Print(1);


    TH1* h;
    TString key = "h2C_vs_Eta";
    if (cellNo < cellHistos[key].size() && cellHistos[key].at(cellNo))
        cellHistos[key].at(cellNo)->Fill(delta_h2C*cellWidth,eta);

    if (!diamondCluster->getClusterSize() || snr < -100){
        cout<<"ERROR: POS: "<<predPos<<" / "<<pos_max<<"/"<<pos_h2c<<"/"<<pos_weighted
                <<" "<<Second_highest_hit_pos<<" - " << highest_hit_pos<<" "<<snr<<endl;
    }

    if (cellNo< vecHResolutionPerCell_maxValue.size()){
        TH1F* histo  = vecHResolutionPerCell_maxValue.at(cellNo);
        if (histo);
        histo->Fill(delta_max*cellWidth);
    }
    if (cellNo< vecHResolutionPerCell_maxValue_vs_SNR.size()){
        TH2F* histo  = (TH2F*)vecHResolutionPerCell_maxValue_vs_SNR.at(cellNo);
        if (histo)
            histo->Fill(delta_max*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_maxValue_vs_PredHit.size()){
        TH2F* histo  =  (TH2F*)vecHResolutionPerCell_maxValue_vs_PredHit.at(cellNo);
        if (histo)
            histo->Fill(delta_max*cellWidth,relPredPos);
    }
    if (cellNo< vecHResolutionPerCell_maxValue_vs_PredHitY.size()){
        TH2F* histo  =  (TH2F*)vecHResolutionPerCell_maxValue_vs_PredHitY.at(cellNo);
        if (histo)
            histo->Fill(delta_max*cellWidth,relPredPosY);
    }
    /*********/
    if (cellNo< vecHResolutionPerCell_chargeWeighted.size()){
        TH1F* histo  = vecHResolutionPerCell_chargeWeighted.at(cellNo);
        if (histo);
        histo->Fill(delta_Weigthed*cellWidth);
    }
    if (cellNo< vecHResolutionPerCell_chargeWeighted_vs_SNR.size()){
        TH2F* histo  = (TH2F*) vecHResolutionPerCell_chargeWeighted_vs_SNR.at(cellNo);
        if (histo)
            histo->Fill(delta_Weigthed*cellWidth,snr);
    }

    if (cellNo< vecHResolutionPerCell_chargeWeighted_vs_PredHit.size()){
        TH2F* histo  = (TH2F*) vecHResolutionPerCell_chargeWeighted_vs_PredHit.at(cellNo);
        //cout<<"FILL vecHResolutionPerCell_chargeWeighted_vs_PredHit:"<<cellNo<<"\t"<<relPredPos<<" --> "<<delta*cellWidth<<endl;
        if (histo)
            histo->Fill(delta_Weigthed*cellWidth,relPredPos);
    }
    else{
        cout<<"Cannot find "<<cellNo<< " in "<<vecHResolutionPerCell_chargeWeighted_vs_PredHit.size()<<endl;
    }

    if (cellNo< vecHResolutionPerCell_chargeWeighted_vs_PredHitY.size()){
        TH2F* histo  = (TH2F*) vecHResolutionPerCell_chargeWeighted_vs_PredHitY.at(cellNo);
        //cout<<"FILL vecHResolutionPerCell_chargeWeighted_vs_PredHit:"<<cellNo<<"\t"<<relPredPos<<" --> "<<delta*cellWidth<<endl;
        if (histo)
            histo->Fill(delta_Weigthed*cellWidth,relPredPosY);
    }
    /*********/
    if (cellNo< vecHResolutionPerCell_highest2Centroid.size()){
        TH1F* histo  = vecHResolutionPerCell_highest2Centroid.at(cellNo);
        if (histo);
        histo->Fill(delta_h2C*cellWidth);
    }
    if (cellNo< vecHResolutionPerCell_highest2Centroid_vs_SNR.size()){
        TH2F* histo  = (TH2F*) vecHResolutionPerCell_highest2Centroid_vs_SNR.at(cellNo);
        if (histo) histo->Fill(delta_h2C*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_highest2Centroid_vs_PredHit.size()){
        TH2F* histo  =  (TH2F*)vecHResolutionPerCell_highest2Centroid_vs_PredHit.at(cellNo);
        if (histo)  histo->Fill(delta_h2C*cellWidth,relPredPos);
    }
    if (cellNo< vecHResolutionPerCell_highest2Centroid_vs_PredHitY.size()){
        TH2F* histo  = (TH2F*)vecHResolutionPerCell_highest2Centroid_vs_PredHitY.at(cellNo);
        if (histo)  histo->Fill(delta_h2C*cellWidth,relPredPosY);
    }

    /*********/
    Float_t delta = snr>settings->GetResolutionSNR()?delta_h2C:delta_max;
    if (cellNo< vecHResolutionPerCell_h2C_WithCut.size()){
        TH1F* histo  = vecHResolutionPerCell_h2C_WithCut.at(cellNo);
        if (histo);
        histo->Fill(delta*cellWidth);
    }
    if (cellNo< vecHResolutionPerCell_h2C_WithCut_vs_SNR.size()){
        TH2F* histo  = (TH2F*)vecHResolutionPerCell_h2C_WithCut_vs_SNR.at(cellNo);
        if (histo) histo->Fill(delta*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_h2C_WithCut_vs_PredHit.size()){
        TH2F* histo  = (TH2F*)vecHResolutionPerCell_h2C_WithCut_vs_PredHit.at(cellNo);
        if (histo)  histo->Fill(delta*cellWidth,relPredPos);
    }

    if (cellNo< vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.size()){
        TH2F* histo  = (TH2F*)vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.at(cellNo);
        if (histo)  histo->Fill(delta*cellWidth,relPredPosY);
    }
}

void TAnalysisOf3DResolution::initialiseHistos() {

    TString key,  cellName;
    UInt_t nCells = settings->GetNCells3d();
    UInt_t nBins = 128;
    Float_t minX = - 1.5*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX = 1.5*settings->GetCellWidth(subjectDetector,2);
    TString name = "hAdjacentSNR_vs_cellNo"+appendix;
    TString title = "hAdjacentSNR_vs_cellNo"+appendix;

    title+=";Cell No;SNR adjacent Strip";
    hAdjacentSNR_vs_cellNo = new TH2F(name,title,nCells,0,nCells,160,-30,50);
    name = "hAdjacentChannels_Signal"+appendix;
    title = "hAdjacentChannels_Signal"+appendix;
    title+=";Signal left Strip / ADC No;Signal right Strip / ADC";
    hAdjacentChannels_Signal = new TH2F(name,title,300,-300,300,300,-300,300);
    TH1* histo;
    for (UInt_t cell = 0; cell <nCells;cell++){
        cellName = TString::Format("hResolution_CellNo_%02d_",cell);

        key = "h2C_vs_Eta";
        name = cellName+key+appendix;
        title = name+";residualt_{h2C}/#um;Eta = #frac{S_R}{S_L+S_R}";
        histo = new TH2F(name,name,nBins,minX,maxX,nBins,-1,1);
        if (cellHistos.find(key) == cellHistos.end() ){
            cellHistos[key] = vector<TH1*>();
//            cout<<"Add "<<key<<" to cellHistoMap"<<endl;
        }
        cellHistos[key].push_back((TH1F*)histo);
//        cout<<" * Add Histo: " <<name<<endl;

        name = TString::Format("hResolution_CellNo_%02d_maxValue",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_maxValue.push_back((TH1F*)histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_chargeWeighted.push_back((TH1F*)histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_highest2Centroid.push_back((TH1F*)histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_h2C_withCut",cell)+appendix;
        title = TString::Format("hResolution Cell %02d - h2C with SNR Cut: %2.1f",cell,settings->GetResolutionSNR());;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_h2C_WithCut.push_back((TH1F*)histo);

        /*******  RESOLUTION VS SNR *********/
        name = TString::Format("hResolution_CellNo_%02d_maxValue_vs_SNR",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        TH2F* histo2 = new TH2F(name,title,nBins,minX,maxX,110,-10,maxsnr);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("SNR 2nd hit");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_maxValue_vs_SNR.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted_vs_SNR",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,110,-10,maxsnr);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("SNR 2nd hit");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_chargeWeighted_vs_SNR.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid_vs_SNR",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,110,-10,maxsnr);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("SNR 2nd hit");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_highest2Centroid_vs_SNR.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_h2C_withCut_vs_SNR",cell)+appendix;
        title = TString::Format("hResolution Cell %02d - h2C with SNR Cut: %2.1f",cell,settings->GetResolutionSNR());;
        histo2 = new TH2F(name,title,nBins,minX,maxX,110,-10,maxsnr);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("SNR 2nd hit");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_h2C_WithCut_vs_SNR.push_back(histo2);

        /*******  RESOLUTION VS PredHit *********/

        name = TString::Format("hResolution_CellNo_%02d_maxValue_vs_PredHit",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_maxValue_vs_PredHit.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted_vs_PredHit",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_chargeWeighted_vs_PredHit.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid_vs_PredHit",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid_vs_PredHit",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_highest2Centroid_vs_PredHit.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_h2C_withCut_vs_PredHit",cell)+appendix;
        title = TString::Format("hResolution Cell %02d - h2C with SNR Cut: %2.1f",cell,settings->GetResolutionSNR());;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_h2C_WithCut_vs_PredHit.push_back(histo2);


        /*******  RESOLUTION VS PredHitY *********/

        name = TString::Format("hResolution_CellNo_%02d_maxValue_vs_PredHitY",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_maxValue_vs_PredHitY.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted_vs_PredHitY",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_chargeWeighted_vs_PredHitY.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid_vs_PredHitY",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid_vs_PredHitY",cell);;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_highest2Centroid_vs_PredHitY.push_back(histo2);

        name = TString::Format("hResolution_CellNo_%02d_h2C_withCut_vs_PredHitY",cell)+appendix;
        title = TString::Format("hResolution Cell %02d - h2C with SNR Cut: %2.1f",cell,settings->GetResolutionSNR());;
        histo2 = new TH2F(name,title,nBins,minX,maxX,160,-80,80);
        histo2->GetXaxis()->SetTitle("Residual / #mum");
        histo2->GetYaxis()->SetTitle("Pred Hit Pos / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.push_back(histo2);
    }
}

void TAnalysisOf3DResolution::saveHistos() {
    std::cout<<"[TAnalysisOf3DResolution::saveHistos()] "<<flush;
    if (!settings->do3dTransparentAnalysis())
        return;
    TString prefix = "hResolution";
    cout<<0<<flush;
    histSaver->CreateResolutionPlots(&vecHResolutionPerCell_chargeWeighted,"chargeWeighted",subjectDetector,appendix);
    histSaver->CreateResolutionPlots(&vecHResolutionPerCell_highest2Centroid,"highest2Centroid",subjectDetector,appendix);
    histSaver->CreateResolutionPlots(&vecHResolutionPerCell_maxValue,"maxValue",subjectDetector,appendix);
    histSaver->CreateResolutionPlots(&vecHResolutionPerCell_h2C_WithCut,"h2C_WithCut",subjectDetector,appendix);

    cout<<1<<flush;
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_SNR,"maxValue_SNR",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_SNR,"chargeWeighted_SNR",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_SNR,"highest2Centroid_SNR",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_SNR,"h2C_WithCut_SNR",prefix,appendix);

    cout<<2<<flush;
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_PredHit,"maxValue_PredHit",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_PredHit,"chargeWeighted_PredHit",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_PredHit,"highest2Centroid_PredHit",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_PredHit,"h2C_WithCut_PredHit",prefix,appendix);

    cout<<3<<flush;
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_PredHitY,"maxValue_PredHitY",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_PredHitY,"chargeWeighted_PredHitY",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_PredHitY,"highest2Centroid_PredHitY",prefix,appendix);
    histSaver->CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_PredHitY,"h2C_WithCut_PredHitY",prefix,appendix);

    for (TCellHistoMap::iterator it = cellHistos.begin();it != cellHistos.end();it++)
        histSaver->CreateTH2_CellPlots(&it->second,"h2C_WithCut_PredHitY",prefix,appendix);

    cout<<4<<flush;
    histSaver->SaveHistogram(hAdjacentSNR_vs_cellNo);
    histSaver->SaveHistogram(hAdjacentChannels_Signal);
    TH1D* pLeft = hAdjacentChannels_Signal->ProjectionX("hSNR_left");
    pLeft->SetTitle("Signal left / ADC");
    pLeft->SetLineColor(kRed);
    TH1D* pRight = hAdjacentChannels_Signal->ProjectionY("hSNR_right");
    pRight->SetTitle("Signal right / ADC");
    pRight->SetLineColor(kGreen);
    THStack *stack = new THStack("hAdjacentChannels_SignalAll","hAdjacentChannels_SignalAll");
    stack->Add(pLeft);
    stack->Add(pRight);
    histSaver->SaveStack(stack,"nostack",true,false,"Signal","number of entries");
    delete stack;
    delete pLeft;
    delete pRight;

    cout<<5<<flush;
    TH1D* pMax = hAdjacentSNR_vs_cellNo->ProjectionY("hAdjacentSNR");
    pMax->SetTitle("SNR max");
    histSaver->SaveHistogram(pMax);
    pLeft = hAdjacentChannels_SNR->ProjectionX("hSNR_left");
    pLeft->SetTitle("SNR left");
    pLeft->SetLineColor(kRed);
    pRight = hAdjacentChannels_SNR->ProjectionY("hSNR_right");
    pRight->SetTitle("SNR right");
    pRight->SetLineColor(kGreen);
    stack = new THStack("hAllSNRs","hAllSNRs");
    stack->Add(pLeft);
    stack->Add(pRight);
    stack->Add(pMax);
    histSaver->SaveStack(stack,"nostack",true,false,"SNR","number of entries");
    TCanvas * c5 = new TCanvas("ccAllSNRs");
    c5->cd();
    stack->Draw("nostack");
    histSaver->SaveCanvas(c5);
    histSaver->SaveHistogram(hAdjacentChannels_SNR);
    delete stack;
    delete pLeft;
    delete pRight;
    delete pMax;
    cout<<"  -> DONE "<<endl;
    //    LongAnalysis_SaveSNRPerCell();
}

void TAnalysisOf3DResolution::deleteHistos() {
    cout<<"[TAnalysisOf3DResolution::deleteHistos]"<<endl;
    delete hAdjacentSNR_vs_cellNo;
    delete hAdjacentChannels_SNR;
    delete hAdjacentChannels_Signal;
    cout<<"* vecHResolutionPerCell_maxValue"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_maxValue.size();i++)
        delete vecHResolutionPerCell_maxValue[i];
    cout<<"* vecHResolutionPerCell_chargeWeighted"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_chargeWeighted.size();i++)
        delete vecHResolutionPerCell_chargeWeighted[i];
    cout<<"* vecHResolutionPerCell_highest2Centroid"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_highest2Centroid.size();i++)
        delete vecHResolutionPerCell_highest2Centroid[i];
    cout<<"* vecHResolutionPerCell_h2C_WithCut"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_h2C_WithCut.size();i++)
        delete vecHResolutionPerCell_h2C_WithCut[i];
    cout<<"* vecHResolutionPerCell_maxValue_vs_SNR"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_maxValue_vs_SNR.size();i++)
        delete vecHResolutionPerCell_maxValue_vs_SNR[i];
    cout<<"* vecHResolutionPerCell_chargeWeighted_vs_SNR"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_chargeWeighted_vs_SNR.size();i++)
        delete vecHResolutionPerCell_chargeWeighted_vs_SNR[i];
    cout<<"* vecHResolutionPerCell_highest2Centroid_vs_SNR"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_highest2Centroid_vs_SNR.size();i++)
        delete vecHResolutionPerCell_highest2Centroid_vs_SNR[i];
    cout<<"* vecHResolutionPerCell_h2C_WithCut_vs_SNR"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_h2C_WithCut_vs_SNR.size();i++)
        delete vecHResolutionPerCell_h2C_WithCut_vs_SNR[i];
    cout<<"* vecHResolutionPerCell_maxValue_vs_PredHit"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_maxValue_vs_PredHit.size();i++)
        delete vecHResolutionPerCell_maxValue_vs_PredHit[i];
    cout<<"* vecHResolutionPerCell_maxValue_vs_PredHitY"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_maxValue_vs_PredHitY.size();i++)
        delete vecHResolutionPerCell_maxValue_vs_PredHitY[i];
    cout<<"* vecHResolutionPerCell_chargeWeighted_vs_PredHit"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_chargeWeighted_vs_PredHit.size();i++)
        delete vecHResolutionPerCell_chargeWeighted_vs_PredHit[i];
    cout<<"* vecHResolutionPerCell_chargeWeighted_vs_PredHitY"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_chargeWeighted_vs_PredHitY.size();i++)
        delete vecHResolutionPerCell_chargeWeighted_vs_PredHitY[i];
    cout<<"* vecHResolutionPerCell_highest2Centroid_vs_PredHit"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_highest2Centroid_vs_PredHit.size();i++)
        delete vecHResolutionPerCell_highest2Centroid_vs_PredHit[i];
    cout<<"* vecHResolutionPerCell_highest2Centroid_vs_PredHitY"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_highest2Centroid_vs_PredHitY.size();i++)
        delete vecHResolutionPerCell_highest2Centroid_vs_PredHitY[i];
    cout<<"* vecHResolutionPerCell_h2C_WithCut_vs_PredHit"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_h2C_WithCut_vs_PredHit.size();i++)
        delete vecHResolutionPerCell_h2C_WithCut_vs_PredHit[i];
    cout<<"* vecHResolutionPerCell_h2C_WithCut_vs_PredHitY"<<endl;
    for (UInt_t i = 0; i < vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.size();i++)
        delete vecHResolutionPerCell_h2C_WithCut_vs_PredHitY[i];
    cout<<"DONE"<<endl;
}
