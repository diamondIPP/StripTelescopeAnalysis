/*
 * TAnalysisOf3DLongAnalysis.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: bachmair
 */

#include "../include/TAnalysisOf3DLongAnalysis.hh"


TAnalysisOf3DLongAnalysis::TAnalysisOf3DLongAnalysis(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna) {
    // TODO Auto-generated constructor stub
    this->histSaver = histSaver;
    this->settings = settings;
    this->bTransAna = bTransAna;
    verbosity = settings->getVerbosity();
    subjectPlane = TPlaneProperties::getDiamondPlane();
    subjectDetector = TPlaneProperties::getDetDiamond();
    useCMN = true;
    PulseHeightBins = 256;
    PulseHeightMin = 1;
    PulseHeightMax = 2800;
    PulseHeightMinMeanCharge = 1;
    PulseHeightMaxMeanCharge = 1500;
    resolutionStudy = new TAnalysisOf3DResolutionStudies(settings,histSaver,bTransAna);
    goodCellsLandau = new TAnalysisOf3DGoodCellsLandau(settings,histSaver,bTransAna);

    this->bTransAna = bTransAna;
    if (this->bTransAna)
        appendix ="_trans";
    else
        appendix ="";
}

TAnalysisOf3DLongAnalysis::~TAnalysisOf3DLongAnalysis() {
    // TODO Auto-generated destructor stub
    delete resolutionStudy;
}

void TAnalysisOf3DLongAnalysis::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid,Float_t chi2x, Float_t chi2y){
    predx = x_pred;
    predy = y_pred;
    fidx = x_fid;
    fidy = y_fid;
    this->chi2x = chi2x;
    this->chi2y = chi2y;
    validClusteredAnalysis= true;
    validTransparentAnalysis = true;
    //check if clustered cluster exists
    checkClusteredAnalysis();
    checkTransparentAnalysis();
    goodCellsLandau->setValidClusterAnalysis(validClusteredAnalysis);
    goodCellsLandau->setValidTransparentAnalysis(validTransparentAnalysis);
    if (!validTransparentAnalysis && ! validClusteredAnalysis)
        return;
    if(!settings->do3dTransparentAnalysis()){
        diamondCluster = &clusteredCluster;
    }
    else
    {
        diamondCluster = transparentCluster;

    }
    this->setDiamondCluster(diamondCluster);
    resolutionStudy->setDiamondCluster(diamondCluster);
    goodCellsLandau->setDiamondCluster(diamondCluster);

    goodCellsLandau->setEventNo(nEvent);
    //    this->setDiamondCluster(diamondCluster);
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.3"<<endl;
    goodCellsLandau->setClusteredCluster(&clusteredCluster);
    resolutionStudy->setTransparentCluster(isTransparentCluster,transparentCluster);
    goodCellsLandau->setTransparentCluster(isTransparentCluster,transparentCluster);
    resolutionStudy->addEvent(cluster,predx,predy,fidx,fidy,chi2x,chi2y);//todo check if need to be moved

    //    if (settings->do3dTransparentAnalysis() && !validTransparentAnalysis)
    //        return;
    if(verbosity>5)cout<<"Cluster Formed."<<endl;
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.4"<<endl;
    pair<int,int> cell = settings->getCellAndQuarterNo(predx,predy);
    UInt_t cellNo = cell.first;
    UInt_t quarterNo = cell.second;

    Float_t charge;
    Int_t pos;
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.5"<<endl;
    if (diamondCluster->hasNegativeCharge(charge,pos,useCMN)){
        if(charge<-50)
            if (hNegativeChargePosition)
                hNegativeChargePosition->Fill(predx,predy);
            else
                cout<<"Invalid hhisto: hNegativeChargePosition"<<endl;
        if (charge < 0)
            if (hNegativeCharges)
                hNegativeCharges->Fill(charge);
            else
                cout<<"Invalid hhisto: hNegativeCharges"<<endl;
    }
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.6"<<endl;
    Int_t area3DwithColumns = 2;
    Int_t area3DwithoutColumns =1;
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.7"<<endl;
    if(!bTransAna){
        if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
            hInvalidCluster->Fill(predx,predy);
            return;
        }
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            hInvalidCluster->Fill(predx,predy);
            return;
        }
    }
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.8"<<endl;
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(predx,predy);

    charge = diamondCluster->getPositiveCharge(false);
    hPulseHeightVsDetectorHitPostionXY->Fill(predx,predy,charge);
    //  hDetXvsDetY3DEvents->Fill(predx,predy,1);
    //analyse Good Cells
    if(settings->IsGoodCell(3,cellNo)){
        goodCellsLandau->addEvent(diamondCluster,predx,predy,fidx,fidy,chi2x,chi2y);
    }
    ///*
    //if (nEvent>=10000) cout<<nEvent<<"Print 6.9"<<endl;
    //    if(settings->get3dMetallisationFidCuts()->getFidCutRegion(predx,predy)==3){
    //
    //        mapPredictedPositionsAllCells[nEvent] = make_pair(predx,predy);
    //        if(validClusteredAnalysis)
    //            mapClusteredAnalysisAllCells[nEvent] = clusteredCluster;
    //        if (validTransparentAnalysis)
    //            mapTransparentAnalysisAllCells[nEvent] = transparentCluster;
    //    }
    //    //*/
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.10"<<endl;
    //    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
    //        int badcell = settings->getDeadCell3D()[i];
    //        int relCellY = cellNo - badcell;
    //        //      cout<<i<<" "<<badcell<<" "<<cellNo<<endl;
    //        //      cout<<"sqrt(relCellY*relCellY): "<<sqrt(relCellY*relCellY)<<endl;
    //        if(relCellY ==0 || sqrt(relCellY*relCellY) <= 1){
    //            // @Iain: What is THAT?
    //            float relY = relPos.second+settings->GetCellHeight() + settings->GetCellHeight()*relCellY;
    //            if(verbosity>5) cout<<"DeadCellAnalysis: relCellY: "<<relCellY<<" relPos.second: "<<relPos.second<<" relY: "<<relY<<endl;
    //            hDeadCellCharge[i]->Fill(relY, charge);
    //            hDeadCellPositions[i]->Fill(predx,predy);
    //        }
    //    }
    //
    //    if(cellNo < hCellsLandau.size())
    //        hCellsLandau.at(cellNo)->Fill(charge);
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.11"<<endl;
    //    if(!settings->do3dTransparentAnalysis()){
    //        if(cellNo < hCellsClusteSize.size()){
    //            Int_t ClusterSize = diamondCluster->getClusterSize()-2;
    //            Int_t MaxHistoClusterSize = hCellsClusteSize.at(0)->GetNbinsX();
    //            if(ClusterSize>MaxHistoClusterSize)     //If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
    //                hCellsClusteSize.at(cellNo)->Fill(MaxHistoClusterSize);
    //            else
    //                hCellsClusteSize.at(cellNo)->Fill(ClusterSize);
    //        }
    //    }
    //
    //    if(nEvent >= 519000 &&nEvent<=520000 && settings->getRunNumber() == 17212)
    //        cout<<"QuarterCellLandaus"<<endl;
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.12"<<endl;
    //    if  (cellNo < hQuarterCellsLandau.size()){
    //        UInt_t size = hQuarterCellsLandau[cellNo].size();
    //        if(quarterNo < size){
    //            hQuarterCellsLandau[cellNo][quarterNo]->Fill(charge);
    //            if(verbosity>6)cout <<TString::Format("F%2d-%d,\t%3.1f",cellNo,quarterNo,charge)<<endl;
    //        }
    //        else{
    //            if(verbosity>6)
    //                cout << TString::Format("E%2d-%d(%d-%d),\t%3.1f",cellNo,quarterNo,
    //                        (int)hQuarterCellsLandau.size(),size,charge)<<endl;
    //        }
    //        if(!settings->do3dTransparentAnalysis()){
    //            if(quarterNo < hQuarterCellsLandau[cellNo].size()){
    //                Int_t ClusterSize = diamondCluster->size()-2;
    //                Int_t MaxHistoClusterSize = hQuarterCellsClusterSize[0].at(0)->GetNbinsX();
    //                if(ClusterSize>MaxHistoClusterSize)     //If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
    //                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(MaxHistoClusterSize);
    //                else
    //                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(ClusterSize);
    //            }
    //        }
    //    }
    //    else{
    //        if(verbosity>6) cout <<TString::Format("X%2d-%d,\t%3.1f",cellNo,quarterNo,charge)<<endl;
    //    }
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.13"<<endl;
    //    Int_t StartClusterSize;
    //    Int_t MaxOverlayClusterSize;
    //    if(settings->do3dTransparentAnalysis()){
    //        StartClusterSize = 1;
    //        MaxOverlayClusterSize = 3;
    //    }
    //    else{
    //        StartClusterSize = 3;
    //        MaxOverlayClusterSize = 3;
    //    }
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.14"<<endl;
    //    for(Int_t ClusterSize = StartClusterSize; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
    //        //cout<<"DiamondCLusterSize: "<<diamondCluster.getClusterSize()<<endl;
    //        //if (nEvent>=10000) cout<<nEvent<<"Print 6.14.0."<<ClusterSize<<endl;
    //        if(settings->do3dTransparentAnalysis()){
    //            diamondCluster->SetTransparentClusterSize(ClusterSize);
    //        }
    //        Float_t charge = diamondCluster->getPositiveCharge(false);
    //        //if (nEvent>=10000) cout<<nEvent<<"Print 6.14.1"<<endl;
    //        LongAnalysis_FillOverlayedHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize);
    //        //if (nEvent>=10000) cout<<nEvent<<"Print 6.14.2"<<endl;
    //        LongAnalysis_FillOverlayCentralColumnHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize, diamondCluster);
    //        //if (nEvent>=10000) cout<<nEvent<<"Print 6.14.3"<<endl;
    //        LongAnalysis_FillOverlayBiasColumnHistos(cellNo,relPos.first,relPos.second,diamondCluster->getPositiveCharge(false), ClusterSize, diamondCluster);
    //        //cout<<"After Fill Bias Column"<<endl;
    //
    //        //to check
    //    }
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.15"<<endl;
    //    LongAnalysis_FillEdgeFreeHistos(predx,predy,charge);
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.16"<<endl;
    //    LongAnalysis_FillRelativeAddedTransparentCharge();
    //    //if (nEvent>=10000) cout<<nEvent<<"Print 6.17"<<endl;
    //    if(settings->do3dTransparentAnalysis())
    //        LongAnalysis_FillResolutionPlots();
}


void TAnalysisOf3DLongAnalysis::initHistos() {
    resolutionStudy->initHistos();
    goodCellsLandau->initHistos();
    Float_t xmin = settings->get3dMetallisationFidCuts()->getXLow(0);
    Float_t xmax = settings->get3dMetallisationFidCuts()->getXHigh(0);
    Float_t ymin = settings->get3dMetallisationFidCuts()->getYLow(0);
    Float_t ymax = settings->get3dMetallisationFidCuts()->getYHigh(0);
    cout<<"RANGE: "<<xmin<<"/"<<xmax<<"\t"<<ymin<<"/"<<ymax<<endl;
    Float_t deltaX = xmax-xmin;
    Float_t deltaY = ymax - ymin;
    xmin = xmin - .2 * deltaX;
    xmin = xmax + .2 * deltaX;
    ymin = ymin - .2 * deltaY;
    ymax = ymax + .2 * deltaY;
    TString name = "hValidEventsDetSpace";
    name.Append(appendix);
    hValidEventsDetSpace = new TH2F(name,name,1024,xmin,xmax,1024,ymin,ymax);
    hValidEventsDetSpace->GetXaxis()->SetTitle("#it{X} / #mum");
    hValidEventsDetSpace->GetYaxis()->SetTitle("#it{Y} / #mum");
    hValidEventsDetSpace->GetZaxis()->SetTitle("number of entries #");

    hInvalidCellNo = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCellNo");
    hInvalidCellNo->SetTitle("hLongAnalysisInvalidCellNo");
    hInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCluster");
    hInvalidCellNo->SetTitle("hLongAnalysisInvalidCluster");

    hNegativeChargePosition = histSaver->GetHistoBinedInCells((TString)"hNegativeChargePositionAllCells"+appendix,4);
    name = "hNegativeCharges";
    name.Append(appendix);
    hNegativeCharges = new TH1F(name,name,8192,-4096,4096);

    name = "hPulseHeightVsDetectorHitPostionXY";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    UInt_t factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXY = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXY->GetXaxis()->SetTitle("#it{x} / #mum");
    hPulseHeightVsDetectorHitPostionXY->GetYaxis()->SetTitle("#it{y} / #mum");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetTitle("charge / ADC");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
}

void TAnalysisOf3DLongAnalysis::setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster) {
    resolutionStudy->setTransparentCluster(isTransparentCluster,transparentCluster);
    goodCellsLandau->setTransparentCluster(isTransparentCluster,transparentCluster);
    this->isTransparentCluster=isTransparentCluster;this->transparentCluster=transparentCluster;
}

void TAnalysisOf3DLongAnalysis::saveHistos(TH1F* hLandauStrip,TH1F* hLandauPhantom) {
    cout<<"TAnalysisOf3DLongAnalysis::saveHistos "<<hLandauStrip<<" "<<hLandauPhantom<<endl;
    resolutionStudy->saveHistos();
    goodCellsLandau->saveHistos(hLandauStrip,hLandauPhantom);
}

void TAnalysisOf3DLongAnalysis::checkClusteredAnalysis() {

    if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())
    {
        if(verbosity>8) cout<< "to high chi2: "<<chi2x<<" "<<chi2y<<endl;
        validClusteredAnalysis=false;
        return;
    }
    Int_t DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(predx,predy);
    if(DiamondPattern !=3){
        if(verbosity>8) cout<< "wrong diamond Pattern: "<<DiamondPattern<<endl;
        validClusteredAnalysis = false;
        return;
    }
    Int_t nClusters  = eventReader->getNDiamondClusters();
    if(nClusters!=1){
        if(verbosity>8) cout<< "wrong no clusters "<<nClusters<<endl;
        validClusteredAnalysis = false;
        return;
    }
    clusteredCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
    if(clusteredCluster.isSaturatedCluster()){
        validClusteredAnalysis = false;
        if(verbosity>8) cout<< "saturated Channel"<<endl;
        return;
    }
    Int_t clusterSize = clusteredCluster.size();
    if(clusterSize>3){
        validClusteredAnalysis = false;
        if(verbosity>8) cout<< "wrong clustersize >3: "<<clusterSize<<endl;
        return;
    }
    if(verbosity>8) cout<<"Valid Clustered Event"<<endl;
}

void TAnalysisOf3DLongAnalysis::checkTransparentAnalysis() {
    if (!settings->do3dTransparentAnalysis()){
        validTransparentAnalysis = false;
        return;
    }
    if(!isTransparentCluster){
        if(validClusteredAnalysis){
            //               if(verbosity>4)cout<<"\n"<<nEvent<<"\tvalid clustered analysis, but invalid transparentAnalysis: "<<predx<<"/"<<predy<<endl;
            if(verbosity>4)cout<<"\tpattern: "<<settings->get3dMetallisationFidCuts()->getFidCutRegion(predx,predy)<<endl;
            if(verbosity>4)settings->get3dMetallisationFidCuts()->Print(4);
            if(verbosity>4)cout<<"\tXdetChannelsSpaceInt"<<flush;
            cout<<settings->diamondPattern.convertMetricToIntChannel(predx);
            if(settings->diamondPattern.convertMetricToIntChannel(predx)<0){
                settings->diamondPattern.setVerbosity(8);
                if(verbosity>4)cout<<"\t"<<settings->diamondPattern.convertMetricToIntChannel(predx);
                settings->diamondPattern.setVerbosity(0);
            }
            if(verbosity>4)cout<<"\t"<<"isSaturated:"<<transparentCluster->isSaturatedCluster()<<endl;
        }

        validTransparentAnalysis = false;
        return;
    }
    Int_t DiamondPattern;
    DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(predx,predy);
    if (verbosity>5 && validClusteredAnalysis && DiamondPattern!=3){
        cout<<"\n"<<nEvent<<"\tvalid clustered analysis, but invalid transparentAnalysis: "<<predx<<"/"<<predy<<" "<<DiamondPattern<<endl;
        settings->get3dMetallisationFidCuts()->Print(1);
    }

    if(DiamondPattern !=3){
        validTransparentAnalysis = false;
    }
}
