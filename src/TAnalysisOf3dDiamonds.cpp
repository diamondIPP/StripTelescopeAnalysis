/*
 * TAnalysisOf3dDiamonds.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair,iain
 */

#include "../include/TAnalysisOf3dDiamonds.hh"

TAnalysisOf3dDiamonds::TAnalysisOf3dDiamonds(TSettings *newSettings) {
    cout<<"\n*********************************************************\n";
    cout<<"*****TAnalysisOf3dDiamonds::TAnalysisOf3dDiamonds********\n";
    cout<<"**********************************************************"<<endl;
    if(newSettings!=0)
        this->settings=newSettings;
    else exit(-1);
    verbosity = settings->getVerbosity();
    predictedPosition=0;
    UInt_t runNumber=settings->getRunNumber();

    html3D = new THTML3DAnalysis(settings);

    //htmlLandau=new THTMLLandaus(settings);

    settings->goTo3dDiamondTreeDir();
    cout<<"Selection Tree file path is: "<<settings->getSelectionTreeFilePath()<<endl;
    cout<<"Alignment file path is: "<<settings->getAlignmentFilePath()<<endl;
    cout<<"Eta distribution file path is: "<<settings->getEtaDistributionPath()<<endl;
    eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
    histSaver=new HistogrammSaver(settings);
    settings->goTo3dDiamondAnalysisDir();

    histSaver->SetPlotsPath(settings->get3dDiamondAnalysisPath());

    html3D->setFileGeneratingPath(settings->get3dDiamondAnalysisPath());
    cout<<"PATH: "<<settings->get3dDiamondAnalysisPath()<<endl;
    histSaver->SetRunNumber(runNumber);
    settings->goTo3dDiamondTreeDir();
    clusteredAnalysis = new TCellAnalysisClass(settings);
    cout<<"end initialise"<<endl;
    vecEdgePredX.resize(settings->get3dEdgeFidCuts()->getNFidCuts());
    vecEdgePredY.resize(settings->get3dEdgeFidCuts()->getNFidCuts());
    vecEdgePulseHeight.resize(settings->get3dEdgeFidCuts()->getNFidCuts());

    for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}//cout<<TPlaneProperties::getNSiliconPlanes()<<endl;}
    subjectPlane = TPlaneProperties::getDiamondPlane();
    subjectDetector = TPlaneProperties::getDetDiamond();

    PulseHeightBins = 256;
    PulseHeightMin = 1;
    PulseHeightMax = 2800;
    PulseHeightMinMeanCharge = 1;
    PulseHeightMaxMeanCharge = 1500;
    maxsnr = 40;

    maxClusterSize3d = 5;
    useCMN = true;
    if(settings->do3dTransparentAnalysis())
        appendix = "_trans";
    else
        appendix = "";

            

    PrintPositions();
    LongAnalysisSaveCellAndQuaterNumbering();
    hAdjacentSNR_vs_cellNo=0;
}

TAnalysisOf3dDiamonds::~TAnalysisOf3dDiamonds() {
    //htmlLandau->generateHTMLFile();

    html3D->createContent();
    html3D->generateHTMLFile();
    settings->GetMacro()->Clone()->Write();
    if(html3D!=0) delete html3D;

    if(eventReader!=0) delete eventReader;
    if(histSaver!=0)   delete histSaver;
    //if(htmlLandau!=0)  delete htmlLandau;
    settings->goToOutputDir();
}
void TAnalysisOf3dDiamonds::PrintPositions(){
    cout<<"Diamond Pattern: \n";
    settings->diamondPattern.Print();
    cout<<"\n";
    std::pair<Int_t,Int_t> bla = settings->diamondPattern.getInterval(2);
    std::cout.precision(7);
    cout<<std::fixed;
    for (int i = bla.first ; i <= bla.second;i++)
        cout<<setw(3)<<i<<":\t"<< std::fixed<<setw(10)<<settings->diamondPattern.convertChannelToMetric(i)<<endl;
    cout<<"\n";
    for(Int_t n = 0; n < settings->getNRows3d();n++){
        UInt_t cell = settings->get3DCellNo(n,(Int_t)2);
        settings->PrintCellPosition(cell,3);
    }
    cout<<endl;
    //cout<<"Press a key"<<endl;
    //char t;
    //cin>>t;
}
void TAnalysisOf3dDiamonds::doAnalysis(UInt_t nEvents) {
    FileNameEnd = "";
    cout<<"analyze selection data..."<<endl;

    //intialise vectors
    for(UInt_t i=0;i<settings->diamondPattern.getNIntervals();i++){
        vecPHDiamondHit.push_back(new vector<float>);
        vecXPredicted.push_back(new vector<float>);
        vecYPredicted.push_back(new vector<float>);
        ptrCanvas.push_back(new TCanvas); //To Create pointer to canvas for 3D Plot
        ptrCanvasEvents.push_back(new TCanvas);
        ptrCanvasMean.push_back(new TCanvas);
    }
    cout<<"Areas to be analysed:"<<endl;
    for(UInt_t i=0; i<settings->diamondPattern.getNIntervals(); i++){
        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        cout<<channels.first<<"-"<<channels.second<<endl;
    }

    initialiseHistos();

    if(nEvents<=0) nEvents=eventReader->GetEntries();
    cout<<"Number of Events: "<<eventReader->GetEntries()<<endl;
    settings->centralRegion3DnH->Print();
    histSaver->SetNumberOfEvents(nEvents);
    if(verbosity>5)settings->diamondPattern.Print();
    for(nEvent=0;nEvent<nEvents;nEvent++){
        TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
        eventReader->LoadEvent(nEvent);
        if(!eventValid()){
            if (verbosity > 7)	cout<<"don't use"<<endl;
            continue;
        }
        if (verbosity > 7)	cout<<"use Event"<<endl;
        // Analyse

        if(settings->do3dTransparentAnalysis()){
            isTransparentCluster = TransparentAnalysis();
        }
        else
            isTransparentCluster = false;

        if(!settings->do3dTransparentAnalysis()){
            diamondCluster = &clusteredCluster;
        }
        else
        {
            diamondCluster = &transparentCluster;
        }
        if (diamondCluster->isSaturatedCluster())
                continue;;
        //cout<<"Before Strip Analysis"<<endl;
        StripAnalysis();
        //cout<<"After Strip Analysis"<<endl;
        if(settings->do3dShortAnalysis() == 1){ShortAnalysis();}
        //cout<<"After Short Analysis"<<endl;
        if(settings->do3dLongAnalysis() == 1){LongAnalysis();}
        //cout<<"After Long Analysis"<<endl;

    }

    saveHistos();
    cout<< "ENTRIES: "<<clusteredAnalysis->getEntries()<<endl;
    createTreeTestHistos();
    clusteredAnalysis->cellAnalysisTree->SaveAs("analysis3d.root");
    TFile *file = new TFile("analysis3d-2.root","RECREATE");
    file->cd();
    TTree* tree = (TTree*)clusteredAnalysis->cellAnalysisTree->Clone("analysisTree");
    tree->Write();
    cout<<"tree: "<<tree->GetEntries()<<endl;
    cout<<gSystem->pwd()<<" "<<file->GetPath()<<" "<<file->GetName()<<endl;
    file->Close();
}

/**
 * checks if the event has one and only one valid silicon cluster in each
 * telescope plane
 * it predicts the position and than sets some variables as fiducialValue,
 * chi2, predictedPsosition and predictedDetectorPosition as global variables
 * since they are used quite often
 * @return wheater the event should be used for analysis or not
 */
bool TAnalysisOf3dDiamonds::eventValid(){
    if(!eventReader->isValidTrack()){
        if (verbosity > 7) cout<<nEvent<<" invalid Track"<<endl;
        return false;
    }
    fiducialValueX = eventReader->getFiducialValueX();
    fiducialValueY = eventReader->getFiducialValueY();
    if (!settings->isInRoughFiducialCut(fiducialValueX,fiducialValueY)){
        if (verbosity > 7)
            cout<<nEvent<<" not in rough fiducial cut: "<<fiducialValueX<<"/"<<fiducialValueY<<endl;
        return false;
    }

    if(predictedPosition) delete predictedPosition;
    predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
    if (!predictedPosition)
        return false;
    chi2x = predictedPosition->getChi2X();
    chi2y = predictedPosition->getChi2Y();

    if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())     //(chi2x>maxChi2||chi2y>maxChi2)
        return false;
    xPredicted = predictedPosition->getPositionX();	//Predicted positions in labframe
    yPredicted = predictedPosition->getPositionY();
    xPredDet = eventReader->getPositionInDetSystem( subjectDetector, xPredicted, yPredicted);
    yPredDet = eventReader->getYPositionInDetSystem(subjectDetector, xPredicted, yPredicted);
    //	cout<<nEvent<<"Valid Track"<<endl;
    hValidEventsFiducialSpace->Fill(fiducialValueX,fiducialValueY);
    hValidEventsDetSpace->Fill(xPredDet,yPredDet);
    return true;
}

void TAnalysisOf3dDiamonds::StripAnalysis() {

    if(!settings->do3dTransparentAnalysis()){
        if(!eventReader->getNDiamondClusters())
            return;
        clusteredCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
        diamondCluster = &clusteredCluster;

    }
    else{
        if(!isTransparentCluster)
            return;
        diamondCluster = & transparentCluster;
    }
    if (diamondCluster->isSaturatedCluster())
            return;

    //cout<<"Transparent Cluster, ";
    //cout<<"Chi Squared, ";
    /*
    if(transparentCluster.getFirstHitChannel()>23 && transparentCluster.getLastHitChannel()<40)
    	printf("Cluster Channels: %f - %f, Charge: %f \n", transparentCluster.getFirstHitChannel(),transparentCluster.getLastHitChannel(),diamondCluster->getPositiveCharge());
     */

    Int_t stripDetector = 1;
    if(settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY)!=stripDetector)
        return;
    //cout<<"Fiducial Cut, ";

    if(eventReader->getNDiamondClusters()!=1)
        return;
    //cout<<"Number Clusters, ";
    //	UInt_t diamond = TPlaneProperties::getDetDiamond();
    //TCluster diamondCluster = eventReader->getCluster(diamond,0);

    int areaStripDetector = 0;
    if (!settings->isClusterInDiaDetectorArea(diamondCluster,areaStripDetector) ){
        return;
    }

    //cout<<"Cluster completely in detector area"<<endl;
    if( !settings->diamondPattern.isValidCluster(diamondCluster)){
        return;
    }

    //cout<<"Is valid cluster"<<endl;
    if (diamondCluster->isSaturatedCluster())
        return;
    //cout<<"Is not saturated"<<endl;

    //cout<<"diamondCluster charge: "<<diamondCluster->getPositiveCharge()<<endl;

    //cout<<"Entry to Strip Histo."<<endl;

    Float_t charge = diamondCluster->getPositiveCharge();
    hLandauStripFidCutXvsFidCutY->Fill(fiducialValueX, fiducialValueY,charge);
    hLandauStripFiducialPosition->Fill(fiducialValueX, fiducialValueY);
    Float_t negativeCharge;
    Int_t clPos;
    UInt_t clsize;
    if (settings->do3dTransparentAnalysis()){
        clsize = diamondCluster->GetTransparentClusterSize();
        diamondCluster->SetTransparentClusterSize(3);
    }
    hLandauStrip->Fill(charge);
    if (!settings->do3dTransparentAnalysis())
        return;
    bool hasNegativeCharge = diamondCluster->hasNegativeCharge(negativeCharge,clPos,useCMN);
    diamondCluster->SetTransparentClusterSize(clsize);
    if (false||!hasNegativeCharge<0){
        cout<<"\nStrip: "<<hasNegativeCharge<<" "<<negativeCharge<<" "<<clPos<<" "<<useCMN<<" "<<nEvent;;
        Int_t ch_neg = diamondCluster->getChannel(clPos);
        Int_t ch_hit = diamondCluster->getTransparentClusterPosition(0);
        cout<<"Neg Position: "<<clPos<<endl;
        cout<<"Neg: "<<ch_neg<<"\t"<<ch_hit<<" = "<< ch_hit-ch_neg<<endl;
        diamondCluster->Print(1);
    }
    hLandauStripNegativeCharges->Fill(negativeCharge,charge);
    hLandauStripNegativeChargesFraction->Fill((int)(negativeCharge < settings->getNegativeChargeCut()));
    //                    " "<<setw(7)<<negativeCharge<<" "<<setw(5)<<settings->getNegativeChargeCut()<<endl;
//    if (negativeCharge < settings->getNegativeChargeCut())
//        hLandauStripNegativeChargesFraction->Fill(1);
//    else
//        hLandauStripNegativeChargesFraction->Fill(0);
    if (negativeCharge<0){
//        cout<<nEvent<<"\tFill: "<<hasNegativeCharge<<" "<<negativeCharge<< " " <<clPos<<" "<< charge<<endl;
        //    else
        //        hLandauStripNegativeCharges->Fill(0.0,charge);
        hLandauStripNegativeChargesClPos->Fill(negativeCharge,clPos);
        if (negativeCharge < settings->getNegativeChargeCut()){
            hLandauStripNegativeChargePosition->Fill(xPredDet,yPredDet);
        }
    }
    else{
//        hLandauStripNegativeCharges->Fill(1,charge);
    }
}

void TAnalysisOf3dDiamonds::ShortAnalysis() {

//    if(!settings->do3dTransparentAnalysis()){
        Float_t maxChi2 = settings->getChi2Cut3D();
        Int_t nClusters =eventReader->getNDiamondClusters();
        hNumberofClusters->Fill(nClusters);
        ClusterPlots(eventReader->getNDiamondClusters(),fiducialValueX,fiducialValueY);
        if (nClusters>0)
            hClusterEventsDetSpace->Fill(xPredDet,yPredDet);

        Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
        ShortAnalysis_FillEdgeAlignmentHistos();
        if(predictedDetector !=1 && predictedDetector !=2 && predictedDetector !=3){
            if (settings->do3dTransparentAnalysis())
                if (xPredDet>1650 && xPredDet<2200 && yPredDet>0 && yPredDet<1650){
                    cout<<"\nreject cluster: "<<predictedDetector<<" = "<<xPredDet<<"/"<<yPredDet<<"\t";
                    diamondCluster->Print();
                }
            return;
        }

        switch (eventReader->getNDiamondClusters()) {
            case 1:
                ShortAnalysis_Analyse1Cluster();
                hRelatviveNumberOfMultipleClusterEvents->Fill(predictedDetector,0);
                hRelatviveNumberOfMultipleClusterEventsSamePattern->Fill(predictedDetector,0);
                break;
            case 2:
                ShortAnalysis_Analyse2Cluster();
            default:
                hRelatviveNumberOfMultipleClusterEvents->Fill(predictedDetector,1);
        }
//    }
//    else
//        ShortAnalysis_Analyse1Cluster();
}

void TAnalysisOf3dDiamonds::ShortAnalysis_FillEdgeAlignmentHistos(){
    Int_t nClusters = eventReader->getNDiamondClusters();
    if (nClusters ==0)
        return;
    if ( nClusters >2)
        return;
    Float_t charge = 0;
    for (UInt_t i = 0; i<nClusters; i++)
        charge += eventReader->getCluster(TPlaneProperties::getDetDiamond(),i).getPositiveCharge(useCMN);
    hFidCutsVsMeanCharge->Fill(xPredDet,yPredDet,charge);
    ShortAnalysis_FillEdgeDistributions(charge);
}

void TAnalysisOf3dDiamonds::ShortAnalysis_Analyse1Cluster(UInt_t clusterNo){

    if(!settings->do3dTransparentAnalysis()){
        clusteredCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),clusterNo);
        diamondCluster = &clusteredCluster;
        if(diamondCluster->isSaturatedCluster())
            return;
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            //		cerr <<" Cluster is invalid: ";
            //		diamondCluster->Print(1);
            if (xPredDet>1650 && xPredDet < 2200)
            {
                cout<<"ShortAnalysis_Analyse1Cluster - invalid cluster: "<<xPredDet<<"/"<<yPredDet;
                diamondCluster->Print(1);
            }
            return;
        }
        HitandSeedCount(diamondCluster);
        Int_t clusterSize = diamondCluster->size()-2;
        vecClusterSize.push_back(clusterSize);
    }
    else{
        if(!isTransparentCluster)
            return;
        diamondCluster =& transparentCluster;
    }
    //Edge Finding
//    ShortAnalysis_FillEdgeDistributions(diamondCluster->getPositiveCharge(false));
    if (diamondCluster->isSaturatedCluster())
        return;
    Float_t clusterCharge = diamondCluster->getPositiveCharge(false);

    ShortAnalysis_FillMeanChargeVector(clusterCharge);
    hTotalAvrgChargeXY->Fill(xPredDet,yPredDet,clusterCharge);
    //Universal PHvsChannel Plot

    for(UInt_t i=0; i < settings->diamondPattern.getNIntervals();i++){

        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        //cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;
        Int_t HighestSignalChannel = diamondCluster->getHighestSignalChannel();

        if(HighestSignalChannel<=channels.second && HighestSignalChannel>=channels.first){

            TFiducialCut *cut = settings->getSelectionFidCuts()->getFidCut(i+1);
            if (!cut){
                cerr<<"Cannot get cut no "<<i+1<<" in "<<settings->getSelectionFidCuts()->getNFidCuts()<<endl;
                cout<<"[ShortAnalysis_Analyse1Cluster] ERROR:"<<flush;
                settings->getSelectionFidCuts()->Print(1);
                return;
            }
            if( cut->IsInFiducialCut(fiducialValueX,fiducialValueY))
                return;

            //hTransparentAnalysisValidClusterFidCutXvsFidCutY->Fill(fiducialValueX, fiducialValueY);
            Float_t charge = diamondCluster->getPositiveCharge(false);
            hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,charge);
            hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);

            hEventsvsChannel[i]->Fill(diamondCluster->getHighestSignalChannel());
            hPHvsChannel[i]->Fill(charge,diamondCluster->getHighestSignalChannel());
            hLandau[i]->Fill(charge);
            vecPHDiamondHit[i]->push_back(charge);
            hFidCutXvsFidCutY[i]->Fill(fiducialValueX,fiducialValueY);

            if(!settings->do3dTransparentAnalysis()){
                hHitandSeedCount[i]->Fill(HitCount,SeedCount);
            }
            hChi2XChi2Y[i]->Fill(chi2x, chi2y);
        }
    }		//End of for diamond patterns
}

void TAnalysisOf3dDiamonds::ShortAnalysis_Analyse2Cluster(){
    TCluster cluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
    TCluster cluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
    Float_t pos1 = cluster1.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t pos2 = cluster2.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t ph1 = cluster1.getPositiveCharge(useCMN);
    Float_t ph2 = cluster2.getPositiveCharge(useCMN);
    Int_t pattern1 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos1));
    Int_t pattern2 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos2));
    Int_t predictedArea = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
    Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
    if (pattern1 == pattern2)
        hRelatviveNumberOfMultipleClusterEventsSamePattern->Fill(predictedDetector,1);
    if(pattern1<0)pattern1=-2;
    if(pattern2<0)pattern2=-2;
    pattern1++;
    pattern2++;
    if (!settings->diamondPattern.isValidChannelPosition(pos1)){
        ShortAnalysis_Analyse1Cluster(1);//todo how can we use them for long analysis
        return;
    }
    if (!settings->diamondPattern.isValidChannelPosition(pos2)){
        ShortAnalysis_Analyse1Cluster(0);//todo how can we use them for long analysis
        return;
    }
    if (pos1>128||pos1<0||pos2>128||pos2<0){
        return;
    }
    ShortAnalysis_FillMeanChargeVector(ph1+ph2);

    Int_t delta1 = pattern1 - predictedDetector;
    Int_t delta2 = pattern2 - predictedDetector;
    hShortAnalysis2ClusterHitPattern_1stCluster->Fill(predictedDetector,delta1);
    hShortAnalysis2ClusterHitPattern_2ndCluster->Fill(predictedDetector,delta2);

    vecPH_Cluster1_ShortAna.push_back(ph1);
    vecPH_Cluster2_ShortAna.push_back(ph2);
    vecCh_Cluster1_ShortAna.push_back(pos1);
    vecCh_Cluster2_ShortAna.push_back(pos2);
    if (predictedDetector == pattern1||predictedDetector==pattern2){
        if(predictedDetector == 2|| true){
            Double_t relCharge = ph2/ph1;
            hRelativeChargeTwoClustersX->Fill(xPredDet,relCharge);
            hRelativeChargeTwoClustersY->Fill(yPredDet,relCharge);
            hRelativeChargeTwoClustersXY->Fill(xPredDet,yPredDet,relCharge);
            hShortAnalysis2TotalChargeXY->Fill(xPredDet,yPredDet,ph1+ph2);
        }
    }
    hTotalAvrgChargeXY->Fill(xPredDet,yPredDet,ph1+ph2);
}


void TAnalysisOf3dDiamonds::LongAnalysis_checkClusteredAnalysis(){
    Float_t maxChi2 = settings->getChi2Cut3D();

    if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())
    {
        if(verbosity>8) cout<< "to high chi2: "<<chi2x<<" "<<chi2y<<endl;
        validClusteredAnalysis=false;
        return;
    }
    Int_t DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
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
    if(nEvent >= 519000 &&nEvent<=520000 && settings->getRunNumber() == 17212)
        cout<<"getCluster "<<flush;
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

void TAnalysisOf3dDiamonds::LongAnalysis_checkTransparentAnalysis(){
    if (!settings->do3dTransparentAnalysis()){
        validTransparentAnalysis = false;
        return;
    }
    if(!isTransparentCluster){
        if(validClusteredAnalysis){
            if (verbosity>4){
                cout<<"[LongAnalysis_checkTransparentAnalysis]"<<flush;
                cout<<"\n"<<nEvent<<"\tvalid clustered analysis, but invalid transparentAnalysis: "<<xPredDet<<"/"<<yPredDet<<endl;
                cout<<"\tpattern: "<<settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet)<<endl;
                settings->get3dMetallisationFidCuts()->Print(4);
                cout<<"\tXdetChannelsSpaceInt"<<flush;
                cout<<settings->diamondPattern.convertMetricToIntChannel(xPredDet);
            }
            if(settings->diamondPattern.convertMetricToIntChannel(xPredDet)<0){
                settings->diamondPattern.setVerbosity(8);
                if(verbosity>4)cout<<"\t"<<settings->diamondPattern.convertMetricToIntChannel(xPredDet);
                settings->diamondPattern.setVerbosity(0);
            }
            if(verbosity>4)cout<<"\t"<<"isSaturated:"<<transparentCluster.isSaturatedCluster()<<endl;
        }

        validTransparentAnalysis = false;
        return;
    }
    Int_t DiamondPattern;
    DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
    if (verbosity>5 && validClusteredAnalysis && DiamondPattern!=3){
        cout<<"\n"<<nEvent<<"\tvalid clustered analysis, but invalid transparentAnalysis: "<<xPredDet<<"/"<<yPredDet<<" "<<DiamondPattern<<endl;
        settings->get3dMetallisationFidCuts()->Print(1);
    }

    if(DiamondPattern !=3){
        validTransparentAnalysis = false;
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis() {

    validClusteredAnalysis= true;
    validTransparentAnalysis = true;
    //check if clustered cluster exists
    LongAnalysis_checkClusteredAnalysis();
    LongAnalysis_checkTransparentAnalysis();
    if (!validTransparentAnalysis && ! validClusteredAnalysis)
        return;
    if(!settings->do3dTransparentAnalysis()){
        diamondCluster = &clusteredCluster;
    }
    else
    {
        diamondCluster = &transparentCluster;

    }
    if (diamondCluster->isSaturatedCluster())
            return;
    LongAnalysis_FillResolutionPlots();
    LongAnalysis_FillChargeSharingPlots();
//    if (settings->do3dTransparentAnalysis() && !validTransparentAnalysis)
//        return;
    if(verbosity>5)cout<<"Cluster Formed."<<endl;

    pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
    UInt_t cellNo = cell.first;
    UInt_t quarterNo = cell.second;

    Float_t charge = diamondCluster->getPositiveCharge(false);
    Float_t negCharge;
    Int_t pos;
    bool hasNegativeCharge = diamondCluster->hasNegativeCharge(negCharge,pos,useCMN);
    Float_t negativeChargeRatio = negCharge/charge;
    Float_t maxCharge = diamondCluster->getHighestSignal(useCMN);
    hNegativeChargeRatio->Fill(negativeChargeRatio,charge);
    hNegativeChargeRatioAbs->Fill(negativeChargeRatio,negCharge);
    hNegativeChargeRatioMax->Fill(negCharge/maxCharge,maxCharge);
    pair<Float_t,Float_t> relPos =  settings->getRelativePositionInCell(xPredDet,yPredDet);
    hNegativeChargeRatioOverlay->Fill(relPos.first,relPos.second,negCharge/maxCharge);
    Float_t dx = relPos.first;
    Float_t dy = relPos.second;
    if (dx > settings->GetCellWidth(subjectDetector,2)/2.)
        dx-=settings->GetCellWidth(subjectDetector,2);
    if (dy > settings->GetCellHeight()/2.)
        dy -= settings->GetCellHeight();
    dx /= settings->GetCellWidth(subjectDetector,2);
    dy /= settings->GetCellHeight();
    if (TMath::Abs(dx) + TMath::Abs(dy) < 0.4){
        hNegativeChargeFieldWireFraction->Fill(xPredDet,yPredDet,int(negCharge<settings->getNegativeChargeCut()));
        hNegativeChargeFieldWirePositions->Fill(xPredDet,yPredDet);
        hNegativeChargeFieldWirePositionsOverlay->Fill(relPos.first,relPos.second);
    }
    if (hasNegativeCharge){
        if(negCharge<settings->getNegativeChargeCut())
            hNegativeChargePosition->Fill(xPredDet,yPredDet);

        if (false&&!hasNegativeCharge<0){
            Int_t pos_neg = diamondCluster->getTransparentClusterPosition(pos-1);
            Int_t pos_hit = diamondCluster->getTransparentClusterPosition(0);
            Int_t delta = pos_neg - pos_neg;
            cout<<"\n"<<endl;
            diamondCluster->Print();
            cout<<"Delta: "<<pos_neg<<" "<<pos_hit<<": "<<delta<<endl;
        }
    }
    if(negCharge<settings->getNegativeChargeCut())
        hNegativeChargeFraction->Fill(1);
    else
        hNegativeChargeFraction->Fill(0);
    Int_t area3DwithColumns = 2;
    Int_t area3DwithoutColumns =1;
    if(!settings->do3dTransparentAnalysis()){
        if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
            hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
            return;
        }
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
            return;
        }
    }

    hPulseHeightVsDetectorHitPostionXY->Fill(xPredDet,yPredDet,charge);

    if (settings->do3dTransparentAnalysis()){
        UInt_t clusterSize = diamondCluster->GetTransparentClusterSize();
        for (UInt_t i = 1; i< diamondCluster->GetMaxTransparentClusterSize();i++){
            diamondCluster->SetTransparentClusterSize(i);
            Float_t transcharge = diamondCluster->getPositiveCharge(false);
            if (i-1<hPulseHeightVsDetectorHitPostionXY_trans.size())
                hPulseHeightVsDetectorHitPostionXY_trans[i-1]->Fill(xPredDet,yPredDet,transcharge);
            else
                cout<<"ERROR hPulseHeightVsDetectorHitPostionXY_trans["<<i-1<<"] doesn't exist"<<endl;
        }
        diamondCluster->SetTransparentClusterSize(clusterSize);
    }
    //	hDetXvsDetY3DEvents->Fill(xPredDet,yPredDet,1);
    //analyse Good Cells
    if(settings->IsGoodCell(3,cellNo)){
        this->LongAnalysis_FillGoodCellsLandaus(charge);
    }
    if(!settings->isBadCell(3,cellNo)){
        mapPredictedPositionsAllButBadCells[nEvent] = make_pair(xPredDet,yPredDet);
        if(validClusteredAnalysis)
            mapClusteredAnalysisAllButBadCells[nEvent] = clusteredCluster;
        if (validTransparentAnalysis)
            mapTransparentAnalysisAllButBadCells[nEvent] = transparentCluster;
    }
    ///*
     if(settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet)==3){
         if(verbosity > 7) cout<< "mapPredictedPositionsAllCells"<<endl;
    	 mapPredictedPositionsAllCells[nEvent] = make_pair(xPredDet,yPredDet);
    	 if(validClusteredAnalysis)
    		 mapClusteredAnalysisAllCells[nEvent] = clusteredCluster;
    	 if (validTransparentAnalysis)
    		 mapTransparentAnalysisAllCells[nEvent] = transparentCluster;
     }
     //*/

    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
        int badcell = settings->getDeadCell3D()[i];
        int relCellY = cellNo - badcell;
        //		cout<<i<<" "<<badcell<<" "<<cellNo<<endl;
        //		cout<<"sqrt(relCellY*relCellY): "<<sqrt(relCellY*relCellY)<<endl;
        if(relCellY ==0 || sqrt(relCellY*relCellY) <= 1){
            // @Iain: What is THAT?
            float relY = relPos.second+settings->GetCellHeight() + settings->GetCellHeight()*relCellY;
            if(verbosity>5) cout<<"DeadCellAnalysis: relCellY: "<<relCellY<<" relPos.second: "<<relPos.second<<" relY: "<<relY<<endl;
            hDeadCellCharge[i]->Fill(relY, charge);
            hDeadCellPositions[i]->Fill(xPredDet,yPredDet);
        }
    }

    if(cellNo < hCellsLandau.size())
        hCellsLandau.at(cellNo)->Fill(charge);

    if(!settings->do3dTransparentAnalysis()){
        if(cellNo < hCellsClusteSize.size()){
            Int_t ClusterSize = diamondCluster->getClusterSize()-2;
            Int_t MaxHistoClusterSize = hCellsClusteSize.at(0)->GetNbinsX();
            if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
                hCellsClusteSize.at(cellNo)->Fill(MaxHistoClusterSize);
            else
                hCellsClusteSize.at(cellNo)->Fill(ClusterSize);
        }
    }

    if(nEvent >= 519000 &&nEvent<=520000 && settings->getRunNumber() == 17212)
        cout<<"QuarterCellLandaus"<<endl;
    if  (cellNo < hQuarterCellsLandau.size()){
        UInt_t size = hQuarterCellsLandau[cellNo].size();
        if(quarterNo < size){
            hQuarterCellsLandau[cellNo][quarterNo]->Fill(charge);
            if(verbosity>6)cout <<TString::Format("F%2d-%d,\t%3.1f",cellNo,quarterNo,charge)<<endl;
        }
        else{
            if(verbosity>6)
                cout <<	TString::Format("E%2d-%d(%d-%d),\t%3.1f",cellNo,quarterNo,
                        (int)hQuarterCellsLandau.size(),size,charge)<<endl;
        }
        if(!settings->do3dTransparentAnalysis()){
            if(quarterNo < hQuarterCellsLandau[cellNo].size()){
                Int_t ClusterSize = diamondCluster->size()-2;
                Int_t MaxHistoClusterSize = hQuarterCellsClusterSize[0].at(0)->GetNbinsX();
                if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(MaxHistoClusterSize);
                else
                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(ClusterSize);
            }
        }
    }
    else{
        if(verbosity>6) cout <<TString::Format("Problem with QuarterCellsLandau: cellNo%2d - Quarter %d,\tCharge: %3.1f",cellNo,quarterNo,charge)<<endl;
    }

    Int_t StartClusterSize;
    Int_t MaxOverlayClusterSize;
    if(settings->do3dTransparentAnalysis()){
        StartClusterSize = 1;
        MaxOverlayClusterSize = 3;
    }
    else{
        StartClusterSize = 3;
        MaxOverlayClusterSize = 3;
    }

    for(Int_t ClusterSize = StartClusterSize; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        if(verbosity>6) cout<<"DiamondClusterSize: "<<diamondCluster->getClusterSize()<<endl;
        if(settings->do3dTransparentAnalysis()){
            diamondCluster->SetTransparentClusterSize(ClusterSize);
        }
        Float_t charge = diamondCluster->getPositiveCharge(false);
        LongAnalysis_FillOverlayedHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize);
        LongAnalysis_FillOverlayCentralColumnHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize, diamondCluster);
        LongAnalysis_FillOverlayBiasColumnHistos(cellNo,relPos.first,relPos.second,diamondCluster->getPositiveCharge(false), ClusterSize, diamondCluster);
        //cout<<"After Fill Bias Column"<<endl;

        //to check
    }
    if(verbosity>6) cout<<"LongAnalysis_FillEdgeFreeHistos"<<endl;
    LongAnalysis_FillEdgeFreeHistos(xPredDet,yPredDet,charge);
    if(verbosity>6) cout<<"LongAnalysis_FillRelativeAddedTransparentCharge"<<endl;
    LongAnalysis_FillRelativeAddedTransparentCharge();
};

bool TAnalysisOf3dDiamonds::TransparentAnalysis() {

    if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())
        return false;

    Int_t DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);

    bool bOutput = (verbosity >5 &&xPredDet >= 3660 && xPredDet <= 3710&& yPredDet >= 0 && yPredDet <= 1500);
    bOutput = bOutput;// || (xPredDet>1650 && xPredDet<2200 && yPredDet>0 && yPredDet<1650);
    TString s = TString::Format("[TAnalysisOf3dDiamonds::TransparentAnalysis] %6d, Hit in %5.1f,%5.1f ",nEvent, xPredDet,yPredDet);
    if(DiamondPattern !=1 && DiamondPattern !=2 && DiamondPattern !=3){
        //cout<<TString::Format("xPredDet: %.2f, yPredDet: %.2f", xPredDet, yPredDet)<<endl;
        hTransparentAnalysisInvalidCluster->Fill(xPredDet, yPredDet);
        if (bOutput) cout<<s<<" invalid Cluster: "<<DiamondPattern<<endl;
        return false;
    }

    Int_t XdetChannelSpaceInt =  settings->diamondPattern.convertMetricToIntChannel(xPredDet);
    Float_t XdetChannelSpace = settings->diamondPattern.convertMetricToChannel(xPredDet);

    if(XdetChannelSpaceInt<0){	//Returns -9998 when hit in channel 93.
        if (bOutput) cout<<s<<" invalid XdetChannelInSpaceInt: "<<XdetChannelSpaceInt<< "/" <<endl;
        return false;
    }

    hTransparentAnalysisValidCluster->Fill(xPredDet, yPredDet);
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->Fill(fiducialValueX, fiducialValueY);

    //Fill fiducial plot here.
    //Calculate Transparent charge Felix

    pair<int,int> channels = settings->diamondPattern.getPatternChannels(DiamondPattern);
    //cout<<TString::Format("3DwH Start Channel: %i, Hit Channel Int: %i, Hit Channel Float: %.2f, clusterSize: %i", channels.first, XdetChannelSpaceInt, XdetChannelSpace, clusterSize)<<endl;

    transparentCluster.clear();
    transparentCluster = TTransparentAnalysis::makeTransparentCluster(eventReader,settings,subjectDetector,XdetChannelSpace,maxClusterSize3d);

    if(transparentCluster.isSaturatedCluster()){
        if (bOutput) cout<<s<<"SaturatedCluster"<<endl;
        return false;
    }

    UInt_t  clusSize = transparentCluster.GetTransparentClusterSize();

//    if (bOutput) cout<<"Valid Cluster"<<endl;
    return true;


}

void TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos() {
    cout<<"[TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos()]"<<endl;
    //Universal histograms
    TString name = "hRelativeChargeTwoClustersX";
    name.Append(appendix);
    Float_t xmax = settings->get3dMetallisationFidCuts()->getXHigh();
    Float_t ymax = settings->get3dMetallisationFidCuts()->getYHigh();
    hRelativeChargeTwoClustersX = new TH2F(name,name,1024,0,xmax,100,0,20);
    hRelativeChargeTwoClustersX->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersX->GetYaxis()->SetTitle("rel. ph: ph_{clus1}/ph_{clus2}");

    name ="hFidCutsVsMeanCharge";
    name.Append(appendix);
    hFidCutsVsMeanCharge = new TProfile2D(name,name, 256,-.3*xmax,xmax*1.3,256,-.3*ymax,ymax*1.3);
    hFidCutsVsMeanCharge->GetXaxis()->SetTitle("X / #mum");
    hFidCutsVsMeanCharge->GetYaxis()->SetTitle("Y / #mum");
    hFidCutsVsMeanCharge->GetZaxis()->SetTitle("avrg total Charge of clusters - nClusters <= 2");

    name = "hRelativeChargeTwoClustersY";
    name.Append(appendix);
    hRelativeChargeTwoClustersY = new TH2F(name,name,1024,0,ymax,100,0,20);
    hRelativeChargeTwoClustersY->GetXaxis()->SetTitle("Y");
    hRelativeChargeTwoClustersY->GetYaxis()->SetTitle("rel. ph: ph_{clus1}/ph_{clus2}");

    name="hRelativeChargeTwoClustersXY";
    name.Append(appendix);
    hRelativeChargeTwoClustersXY = new TProfile2D(name,name,128,0,xmax,128,0,ymax);
    hRelativeChargeTwoClustersXY->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersXY->GetYaxis()->SetTitle("Y");
    hRelativeChargeTwoClustersXY->GetZaxis()->SetTitle("avrg. rel. ph: ph{clus1}/ph_{clus2}");

    name ="hShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    hShortAnalysis2TotalChargeXY = new TProfile2D(name,name,128,0,xmax,128,0,ymax);
    hRelativeChargeTwoClustersXY->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersXY->GetYaxis()->SetTitle("Y");
    hRelativeChargeTwoClustersXY->GetZaxis()->SetTitle("total charge: ph_{clus_1{}} + ph_{clus_{2}} /ADC");

    name = "hShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    hTotalAvrgChargeXY = new TProfile2D(name,name, xmax/10.,0,xmax,ymax/10.,0,ymax);
    hTotalAvrgChargeXY->Draw();
    hTotalAvrgChargeXY->GetXaxis()->SetTitle("#it{x} position / #mum");
    hTotalAvrgChargeXY->GetYaxis()->SetTitle("#it{y} position / #mum");
    hTotalAvrgChargeXY->GetZaxis()->SetTitle("avrg signal / ADC");
    hTotalAvrgChargeXY->GetZaxis()->SetTitleOffset(1.2);

    name = "hRelatviveNumberOfMultipleClusterEvents";
    name.Append(appendix);
    hRelatviveNumberOfMultipleClusterEvents = new TProfile(name,name,3,.5,3.5);
    hRelatviveNumberOfMultipleClusterEvents->GetXaxis()->SetTitle("predicted pattern");
    hRelatviveNumberOfMultipleClusterEvents->GetYaxis()->SetTitle("#_{multiple cluster events}/#_{total Events}");

    name = "hRelatviveNumberOfMultipleClusterEventsSamePattern";
    name.Append(appendix);
    hRelatviveNumberOfMultipleClusterEventsSamePattern = new TProfile(name,name,3,.5,3.5);
    hRelatviveNumberOfMultipleClusterEventsSamePattern->GetXaxis()->SetTitle("predicted pattern");
    hRelatviveNumberOfMultipleClusterEventsSamePattern->GetYaxis()->SetTitle("#_{multiple cluster events}/#_{total Events}");

    //hNumberofClusters
    name = "hNumberofClusters";
    hNumberofClusters = new TH1F(name,name,4,0,4);
    name.Append(appendix);
    hNumberofClusters->GetXaxis()->SetTitle("Number of Clusters");
    hNumberofClusters->GetYaxis()->SetTitle("Number of Entries #");

    //hEventsvsChannelCombined
    name = "hEventsvsChannelCombined";
    name.Append(appendix);
    hEventsvsChannelCombined = new TH1F(name,name,100,0,100);
    hEventsvsChannelCombined->GetXaxis()->SetTitle("Channel");
    hEventsvsChannelCombined->GetYaxis()->SetTitle("Number of Entries #");

    //hDoubleClusterPos
    name="hDoubleClusterPos";
    name.Append(appendix);
    hDoubleClusterPos = new TH1F(name,name,80,20,100);
    hDoubleClusterPos->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos->GetYaxis()->SetTitle("Number of Entries #");

    //hDoubleClusterPos0
    name ="hDoubleClusterPos0";
    name.Append(appendix);
    hDoubleClusterPos0 = new TH1F(name,name,80,20,100);
    hDoubleClusterPos0->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos0->GetYaxis()->SetTitle("Number of Entries #");
    hDoubleClusterPos0->SetFillColor(2);

    //hDoubleClusterPos1
    name = "hDoubleClusterPos1";
    name.Append(appendix);
    hDoubleClusterPos1 = new TH1F(name,name,80,20,100);
    hDoubleClusterPos1->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos1->GetYaxis()->SetTitle("Number of Entries #");
    hDoubleClusterPos1->SetFillColor(3);

    //hLandauCluster1
    name = "hLandauCluster1";
    name.Append(appendix);
    hLandauCluster1 = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauCluster1->GetXaxis()->SetTitle("Number of Clusters");
    hLandauCluster1->GetYaxis()->SetTitle("Number of Entries #");
    hLandauCluster1->SetFillColor(2);

    //hLandauCluster2
    name = "hLandauCluster2";
    name.Append(appendix);
    hLandauCluster2 = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauCluster2->GetXaxis()->SetTitle("Number of Clusters");
    hLandauCluster2->GetYaxis()->SetTitle("Number of Entries #");
    hLandauCluster2->SetFillColor(3);

    //hLandauDoubleCombined
    name = "hLandauDoubleCombined";
    name.Append(appendix);
    hLandauDoubleCombined = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauDoubleCombined->GetXaxis()->SetTitle("Number of Clusters");
    hLandauDoubleCombined->GetYaxis()->SetTitle("Number of Entries #");

    cout<<"loop over patterns"<<endl;
    for(UInt_t i=0; i<settings->diamondPattern.getNIntervals(); i++){
        if(verbosity) cout<<"Loop: "<<i<<endl;
        pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
        //hLandau
        name = TString::Format("hLandau_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second) +appendix;
        if(verbosity>1) cout<<"Create "<<name<<endl;
        hLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        if(hLandau.back()){
            hLandau.back()->GetXaxis()->SetTitle("PH of diamond cluster");
            hLandau.back()->GetYaxis()->SetTitle("number of entries #");
        }
        else
            cerr<<"hLandau:'"<<name<<"' wasn't created correctly"<<endl;

        //hPHvsChannel
        name = TString::Format("hEventsvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second) + appendix;
        hEventsvsChannel.push_back(new TH1F(name,name,100,0,100));
        if(hEventsvsChannel.back()){
            hEventsvsChannel.back()->GetXaxis()->SetTitle("HighestPH [ch]");
            hEventsvsChannel.back()->GetYaxis()->SetTitle("No. Events");
        }

        //hPHvsChannel
        name = TString::Format("hPHvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second) + appendix;
        hPHvsChannel.push_back(new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,100,0,100));
        if(hPHvsChannel.back()){
            hPHvsChannel.back()->Draw();
            hPHvsChannel.back()->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->GetXaxis()->SetLimits(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->SetAxisRange(PulseHeightMin, PulseHeightMax, "X");
            hPHvsChannel.back()->GetYaxis()->SetRangeUser(channels.first-1,channels.second+1);
            hPHvsChannel.back()->GetXaxis()->SetTitle("cluster charge /ADC");
            hPHvsChannel.back()->GetYaxis()->SetTitle("XPos /ch");
        }
        //hHitandSeedCount
        name = TString::Format("hHitandSeedCount_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hHitandSeedCount.push_back(new TH2F(name,name,10,-.5,9.5,10,-.5,9.5));
        hHitandSeedCount.back()->GetXaxis()->SetTitle("Hit Count");
        hHitandSeedCount.back()->GetYaxis()->SetTitle("Seed Count");

        //hChi2XChi2Y
        name = TString::Format("hChi2X_vs_chi2Y_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hChi2XChi2Y.push_back(new TH2F(name,name,60,0,30,60,0,30));
        hChi2XChi2Y.back()->GetXaxis()->SetTitle("#chi^{2}_{X}");
        hChi2XChi2Y.back()->GetYaxis()->SetTitle("#chi^{2}_{Y}");

        //hFidCutXvsFidCutY
        name = TString::Format("hFidCutXvsFidCutY_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hFidCutXvsFidCutY.push_back(new TH2F(name,name,160,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),120,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutY.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutY.at(i)->GetYaxis()->SetTitle("FidCutY");

        //hFidCutXvsFidCutYvsCharge		For TH2D

        name = TString::Format("hFidCutXvsFidCutYvsCharge_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hFidCutXvsFidCutYvsCharge.push_back(new TH2D(name,name,
                213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
                160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
        hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        name = "hFidCutXvsFidCutYvsEvents";
        name.Append(appendix);
        hFidCutXvsFidCutYvsEvents.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone(name));

        //hFidCutXvsFidCutYvsMeanCharge
        name = "hFidCutXvsFidCutYvsMeanCharge";
        name.Append(appendix);
        hFidCutXvsFidCutYvsMeanCharge.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone(name));

        //hXdetvsYdetvsCharge		For TH2D     -------> Error is around here!!!!!
        Int_t DiamondPattern = i+1;
        Float_t xLow = settings->get3dMetallisationFidCuts()->getXLow(DiamondPattern);
        Float_t xHigh = settings->get3dMetallisationFidCuts()->getXHigh(DiamondPattern);
        Float_t yLow = settings->get3dMetallisationFidCuts()->getYLow(DiamondPattern);
        Float_t yHigh = settings->get3dMetallisationFidCuts()->getYHigh(DiamondPattern);
        if(verbosity) cout<<"("<<xLow<<"-"<<xHigh<<", "<<yLow<<"-"<<yHigh<<")"<<endl;
        Float_t xDiv = (xHigh - xLow)/5;
        Float_t yDiv = (yHigh - yLow)/5;
        name = TString::Format("hXdetvsYdetvsCharge_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        //name = TString::Format("hXdetvsYdetvsCharge%02d_%02d%",channels.first,channels.second);
        name.Append(appendix);
        hXdetvsYdetvsCharge.push_back(new TH2D(name,name,xDiv,xLow,xHigh,yDiv,yLow,yHigh));
        hXdetvsYdetvsCharge.at(i)->GetXaxis()->SetTitle("X (um)");
        hXdetvsYdetvsCharge.at(i)->GetYaxis()->SetTitle("Y (um)");
        hXdetvsYdetvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        name ="hXdetvsYdetvsEvents";
        name.Append(appendix);
        hXdetvsYdetvsEvents.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone(name));

        //hFidCutXvsFidCutYvsMeanCharge
        name ="hXdetvsYdetvsMeanCharge";
        name.Append(appendix);
        hXdetvsYdetvsMeanCharge.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone(name));
    }

    //hFidCutXvsFidCutYvsMeanChargeAllDetectors
    hFidCutXvsFidCutYvsMeanChargeAllDetectors = (TH2D*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("hFidCutXvsFidCutYvsMeanChargeAllDetectors");

    for(int i=0;i<7;i++){
        //hFidCutXvsFidCutYClusters	For TH2D	{0,1,1_1Seed,2_FirstCluster,2_SecondCluster,3}
        name = TString::Format("hFidCutXvsFidCutYClusters%i",i);
        name.Append(appendix);
        hFidCutXvsFidCutYClusters.push_back(new TH2D(name,name,213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutYClusters.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutYClusters.at(i)->GetYaxis()->SetTitle("FidCutY");
        hFidCutXvsFidCutYClusters.at(i)->GetZaxis()->SetTitle("Events");
    }

    name ="hShortAnalysis2ClusterHitPattern_1stCluster" + appendix;
    int nFidCuts =settings->get3dMetallisationFidCuts()->getNFidCuts();
    hShortAnalysis2ClusterHitPattern_1stCluster =  new TH2F(name,name,
            nFidCuts,.5,nFidCuts+.5,
            5,-nFidCuts+.5,nFidCuts-.5);
    hShortAnalysis2ClusterHitPattern_1stCluster->GetXaxis()->SetTitle("predicteded hit pattern");
    hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()->SetTitle("cluster_{1}-hit-pattern - predicted-hit-pattern");

    name = "hShortAnalysis2ClusterHitPattern_2ndCluster" + appendix;
    hShortAnalysis2ClusterHitPattern_2ndCluster = (TH2F*)hShortAnalysis2ClusterHitPattern_1stCluster->Clone(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->SetTitle(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->GetYaxis()->SetTitle("cluster_{2}-hit-pattern - predicted-hit-pattern");

}

void TAnalysisOf3dDiamonds::initialise3DYAlignmentHistos() {
    //Fiducial Region with Edge Alignment Regions Highlighted

    //hDeadCellsProfile
    cout<<"settings->getDeadCell3D().size(): "<<settings->getDeadCell3D().size()<<endl;
    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
        TString name = TString::Format("hDeadCell_no_%d_MeanCharge",i);
        hDeadCellCharge.push_back(new TProfile(name,name,90,0,450));
        hDeadCellCharge.back()->GetXaxis()->SetTitle("pred hit position y in diamond / #mum");
        hDeadCellCharge.back()->GetYaxis()->SetTitle("Mean Charge [ADC]");
        name = TString::Format("hDeadCell_no_%d_HitPositions",i);
        hDeadCellPositions.push_back(histSaver->GetHistoBinedInCells(name,30));
        hDeadCellPositions.back()->GetXaxis()->SetTitle("pred hit position x in diamond / #mum");
        hDeadCellPositions.back()->GetYaxis()->SetTitle("pred hit position y in diamond / #mum");
    }
    cout<<"End initialise3DYAlignmentHistos()"<<endl;
};

void TAnalysisOf3dDiamonds::initialise3DOverviewHistos() {

    //hDetXvsDetY3DTotolCharge
    TString name = "hPulseHeightVsDetectorHitPostionXY";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    UInt_t factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXY = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXY->GetXaxis()->SetTitle("#it{x} / #mum");
    hPulseHeightVsDetectorHitPostionXY->GetYaxis()->SetTitle("#it{y} / #mum");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetTitle("charge / ADC");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    for (UInt_t i = 0; i< 6;i++){
        name = TString::Format("hPulseHeightVsDetectorHitPostionXY_clusterSize_%d",i+1);
        name+=appendix;
        hPulseHeightVsDetectorHitPostionXY_trans.push_back((TProfile2D*)hPulseHeightVsDetectorHitPostionXY->Clone(name));
    }

    //hPulseHeightVsDetectorHitPostionXYGoodCells
    name = "hPulseHeightVsDetectorHitPostionXYGoodCells";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXYGoodCells = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetXaxis()->SetTitle("#it{x} / #mum");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetYaxis()->SetTitle("#it{y} / #mum");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetZaxis()->SetTitle("charge /ADC");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    //hDetXvsDetY3DEvents
    //	hDetXvsDetY3DEvents = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Clone("hDetXvsDetY3DvsEvents");

    //hDetXvsDetY3DMeanCharge


    //hDetXvsDetY3DRebinnedMeanChargeRMS
    name = "h3DdetRebinnedRMS";
    name.Append(appendix);
    hDetXvsDetY3DRebinnedRMS = histSaver->GetHistoBinedInCells(name);
    //			new TH2D(hDetXvsDetY3DRebinnedRMSName.str().c_str(),hDetXvsDetY3DRebinnedRMSName.str().c_str(),
    //			settings->getNColumns3d(),getXMetalisationStart3d,getXMetalisationEnd3d,
    //			settings->getNRows3d(),getYMetalisationStart3d,getYMetalisationEnd3d);
    hDetXvsDetY3DRebinnedRMS->GetXaxis()->SetTitle("Xdet (um)");
    hDetXvsDetY3DRebinnedRMS->GetYaxis()->SetTitle("Ydet (um)");
    hDetXvsDetY3DRebinnedRMS->GetZaxis()->SetTitle("Charge ADC");

    //hBinnedMeanCharge
    name = "h3DdetCellMeanChargeBinned";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    hBinnedMeanCharge = new TH1F(name,name,9,400,1300);
    hBinnedMeanCharge->GetXaxis()->SetTitle("MeanCharge");
    hBinnedMeanCharge->GetYaxis()->SetTitle("Entries");

    //hDetXvsDetY3DOverview
    name = "hDetXvsDetY3DOverview";
    name.Append(appendix);
    hDetXvsDetY3DOverview = histSaver->GetHistoBinedInCells(name);
    hDetXvsDetY3DOverview->GetXaxis()->SetTitle("Xdet (#mum)");
    hDetXvsDetY3DOverview->GetYaxis()->SetTitle("Ydet (#mum)");
    //hDetXvsDetY3DOverview->GetZaxis()->SetTitle();



    //hCellsMeanClusteSize
    name = "hCellsMeanClusteSize";
    name.Append(appendix);
    hCellsMeanClusteSize = histSaver->GetHistoBinedInCells(name);

    //hQuarterCellsMeanClusteSize
    name = "hQuarterCellsMeanClusterSize";
    name.Append(appendix);
    hQuarterCellsMeanClusterSize = hCellsMeanClusteSize = histSaver->GetHistoBinedInQuarters(name);

    //RebinnedQuarterCellFails
    name = "3DdetNumberofQuarterCellFails";
    name.Append(appendix);
    RebinnedQuarterCellFails = histSaver->GetHistoBinedInCells(name);
    RebinnedQuarterCellFails->GetXaxis()->SetTitle("Xdet (um)");
    RebinnedQuarterCellFails->GetYaxis()->SetTitle("Ydet (um)");
    RebinnedQuarterCellFails->GetZaxis()->SetTitle("Quarter Fails");

    //hDetXvsDetY3DQuarterCellGrading
    for(int k=0; k<6; k++){
        name = TString::Format("hDetXvsDetY3DMeanChargeQuarterCellGrading_%d_Fail",k);
        name.Append(appendix);
        hDetXvsDetY3DMeanChargeQuarterCellGrading.push_back(histSaver->GetHistoBinedInQuarters(name));
        hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetXaxis()->SetTitle("Xdet (um)");
        hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetYaxis()->SetTitle("Ydet (um)");
        hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetZaxis()->SetTitle("Charge ADC");
    }

    //h3DdetQuarterCellFluctuation
    name = "h3DdetQuarterCellFluctuation";
    name.Append(appendix);
    h3DdetQuarterCellFluctuation = histSaver->GetHistoBinedInQuarters(name);
    h3DdetQuarterCellFluctuation->GetXaxis()->SetTitle("Xdet (um)");
    h3DdetQuarterCellFluctuation->GetYaxis()->SetTitle("Ydet (um)");
    h3DdetQuarterCellFluctuation->GetZaxis()->SetTitle("Fluctuation");

    //h3DdetQuarterCellFluctuation1
    name = "h3DdetQuarterCellFluctuation1";
    name.Append(appendix);
    h3DdetQuarterCellFluctuation1 = histSaver->GetHistoBinedInQuarters(name);
    if(h3DdetQuarterCellFluctuation1){
        h3DdetQuarterCellFluctuation1->GetXaxis()->SetTitle("Xdet (um)");
        h3DdetQuarterCellFluctuation1->GetYaxis()->SetTitle("Ydet (um)");
        h3DdetQuarterCellFluctuation1->GetZaxis()->SetTitle("Fluctuation");
    }
    else
        cerr<<name<<" ist not created correctly"<<endl;
    if(verbosity>3)cout<<"#(#"<<endl;

    //hDetXvsDetY3DMeanChargeQuarterCellGradingLandau
    for(int k=0;k<5;k++){
        //For Transparent analysis.
        //hQuarterCellGradedTransparentLandau
        name = TString::Format("hQuarterCellGradedTransparentLandau_Grade%d",k);
        name.Append(appendix);
        if(verbosity>1) cout<<"Create "<<name<<endl;
        hQuarterCellGradedTransparentLandau.push_back(new TH1F(name,name,256,0,2800));

        name = TString::Format("hDetXvsDetY3DMeanChargeQuarterCellGradingLandau_Grade%d",k);
        name.Append(appendix);
        if(verbosity>1)cout<<"Create "<<name<<endl;
        hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.push_back(new TH1F(name,name,256,0,2800));

        //hCellsDeltaXQuarterCellGrading
        name = TString::Format("hCellsDeltaXQuarterCellGrading_Grade%d",k);
        name.Append(appendix);
        if(verbosity>1)cout<<"Create "<<name<<endl;
        hCellsDeltaXQuarterCellGrading.push_back(new TH1F(name,name,100,-3,3));
    }

    //Define Landau function for Landau fit.
    Landau = new TF1("Landau","landau(0)",20,80);
    //
    hQuarterCellsLandau.resize(settings->GetNCells3d());
    hQuarterCellsClusterSize.resize(settings->GetNCells3d());

    for (UInt_t cell =0; cell < settings->GetNCells3d(); cell ++){
        TString name = TString::Format("hCellsLandau_no_%03d",cell);
        name.Append(appendix);
        if(verbosity>1)cout<<"Create "<<name<<endl;
        hCellsLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = "hCellsClusteSize";
        if(verbosity>1) cout<<"Create "<<name<<endl;
        hCellsClusteSize.push_back(new TH1F(name,name,6,0,6));

        //		hQuarterCellsLandau.push_back(vector<TH1F*>());
        hQuarterCellsLandau[cell].resize(settings->getNQuarters3d());
        for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){
            name = TString::Format("hQuaterCellsLandau_%d_%d",cell,quarter);
            name.Append(appendix);
            if(verbosity>1)cout<<"Create "<<name<<endl;
            hQuarterCellsLandau[cell][quarter] =
                    new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
            name = TString::Format("hQuarterCellsClusterSize_%d_%d",cell,quarter);
            name.Append(appendix);
            hQuarterCellsClusterSize[cell].push_back(new TH1F(name,name,6,0,6));
        }

    }

    //Transparent Analysis
    //hCellsTransparentHitPositionCellGraded
    for(int i=0;i<4;i++){
        name = TString::Format("hCellsTransparentHitPositionCellGraded%d",i);
        name.Append(appendix);
        //@iain: what are this hardcoded values? is it relative hit in cell?
        hCellsTransparentHitPositionCellGraded.push_back(new TH2D(name,name,30,0,150,30,0,150));
    }

    //hTransparentCharge3D
    name = "hTransparentCharge3D";
    name.Append(appendix);
    hTransparentCharge3D = new TH1F(name,name,256,0,2800);
    hTransparentCharge3D->GetXaxis()->SetTitle("Mean Charge [ADC]");
    hTransparentCharge3D->GetYaxis()->SetTitle("Entries");

    //hCellsLandauGraded    &&    hCellsLandauGradedNoColumn
    for(int i=0;i<12;i++){	//Group CellsLandaus within same ranges together. 0-100; 100-200; -> 1100-1200;
        name = TString::Format("hLandauCellsGraded_%d_to_%d",(i*100),(i+1)*100);
        name.Append(appendix);
        hCellsLandauGraded.push_back(new TH1F(name,name,256,0,2800));
        name = TString::Format("hLandauCellsGradedNoColumn_%d_to_%d",i*100,(i+1)*100);
        name.Append(appendix);
        hCellsLandauGradedNoColumn.push_back(new TH1F(name,name,256,0,2800));
    }

    //To create Strip detector Landau
    pair<int,int> channels =settings->diamondPattern.getPatternChannels(1);
    //hLandau
    name = TString::Format("hLandau_ch_%d_to_%d",channels.first, channels.second);
    name.Append(appendix);
    hLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
    hLandau.back()->GetXaxis()->SetTitle("PH of diamond cluster");
    hLandau.back()->GetYaxis()->SetTitle("number of entries #");

};

void TAnalysisOf3dDiamonds::initialiseEdgeFreeHistos(){
    Int_t factor = 2;
    TString name ="hPulseHeigthCentralRegion";
    name.Append(appendix);
    hPulseHeigthCentralRegion = histSaver->GetProfile2dBinedInCells(name,factor);
    name ="hPulseHeigthEdgeRegion";
    name.Append(appendix);
    hPulseHeigthEdgeRegion =histSaver->GetProfile2dBinedInCells(name,factor);

    hEventsEdgeRegion = new TH2F("hEventsEdgeRegion","hEventsEdgeRegion",150,0,150,150,0,150);
    hEventsEdgeRegion->GetXaxis()->SetTitle("rel x predicted / #mum");
    hEventsEdgeRegion->GetYaxis()->SetTitle("rel y predicted / #mum");

    hEventsCentralRegion = new TH2F("hEventsCentralRegion","hEventsCentralRegion",150,0,150,150,0,150);
    hEventsCentralRegion->GetXaxis()->SetTitle("rel x predicted / #mum");
    hEventsCentralRegion->GetYaxis()->SetTitle("rel y predicted / #mum");
    hEventsCentralRegion->GetXaxis()->SetLimits(0,150);
    hEventsCentralRegion->GetYaxis()->SetLimits(0,150);


}

void TAnalysisOf3dDiamonds::initialise3DCellOverlayHistos() {
    TH1F* hLandau = new TH1F("hLandau","hLandau",PulseHeightBins, PulseHeightMin,PulseHeightMax);
    hLandau->GetXaxis()->SetTitle("pulse height of cluster /ADC");
    hLandau->GetYaxis()->SetTitle("number of entries #");
    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2 += appendix;

        //Cell Overlay
        //Double_t xBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
        Double_t xBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
        Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

        //Double_t yBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
        Double_t yBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
        Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;
        //hCellsOverlayed
        TString name = "hCellsOverlayAvrgCharge";
        name.Append(appendix);
        /*hCellsOverlayAvrgCharge.push_back(new TProfile2D(name,name,
				30,0,settings->GetCellWidth(subjectDetector,2),
				30,0,settings->GetCellHeight()));*/
        hCellsOverlayAvrgCharge.push_back(new TProfile2D(name,name,
                nXbins,xBinEdges,
                nYbins,yBinEdges));
        hCellsOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("#it{x} position within a cell / #mum");
        hCellsOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("#it{y} position within a cell / #mum");
        hCellsOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster / ADC");
        //	hCellsOverlayAvrgCharge->SetContour(99);
        name = "hCellsOverlayAvrgCharge_noColumnHit";
        name.Append(appendix2);
        hCellsOverlayAvrgChargeNoColumnHit.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->SetTitle(name);

        if(verbosity>3)cout<<hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)<<" "<<hCellsOverlayAvrgCharge.at(ClusterSize)<<endl;
        if(verbosity>3)cout<<" "<<hCellsOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;

        //hCellsOverlayAvrgChargeMinusBadCells
        name = "hCellsOverlayAvrgChargeMinusBadCells";
        name.Append(appendix2);
        hCellsOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

        name = TString::Format( "hCellLandauChargeMinusBadCells_clustersize%0d",ClusterSize+1);
        name.Append(appendix2);
        hCellsLandauMinusBadCells.push_back((TH1F*)hLandau->Clone(name));
        hCellsLandauMinusBadCells.back()->SetTitle(name);


        //hCellsOverlayAvrgChargeGoodCells
        name = "hCellsOverlayAvrgChargeGoodCells";
        name.Append(appendix2);
        hCellsOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

        //hCellsOverlayAvrgChargeBadCells
        name = "hCellsOverlayAvrgChargeBadCells";
        name.Append(appendix2);
        hCellsOverlayAvrgChargeBadCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->SetTitle(name);

        //hCellsOverlayedColumnLandau
        name = "hCellOverlayWithColumnLandau";
        name.Append(appendix2);
        hCellOverlayWithColumnLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
        hCellOverlayWithColumnLandau.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
        hCellOverlayWithColumnLandau.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

        //hCellsOverlayedLandauNoColumn
        name = "hCellOverlayNoColumnLandau";
        name.Append(appendix);
        hCellOverlayNoColumnLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph

        name = "hCellOverlayColumnLandau";
        name.Append(appendix2);
        hCellOverlayColumnLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph


        /*for(int i=0;i<9;i++){
			hCellsOverlayedChargeBinAlignment.push_back(new TH2D("OverlayedCharge","OverlayedCharge",15,0,150,15,0,150));
			hCellsOverlayedChargeBinAlignment.at(i)->SetStats(kFALSE);
			hCellsOverlayedEventsBinAlignment.push_back(new TH2D("OverlayedEvents","OverlayedEvents",15,0,150,15,0,150));
			hCellsOverlayedEventsBinAlignment.at(i)->SetStats(kFALSE);
			stringstream hCellsOverlayedMeanChargeBinAlignmentName; hCellsOverlayedMeanChargeBinAlignmentName<<"hCellsOverlayedMeanChargeBinAlignment"<<FileNameEnd;
			hCellsOverlayedMeanChargeBinAlignment.push_back(new TH2D(hCellsOverlayedMeanChargeBinAlignmentName.str().c_str(),hCellsOverlayedMeanChargeBinAlignmentName.str().c_str(),15,0,150,15,0,150));
			hCellsOverlayedMeanChargeBinAlignment.at(i)->GetXaxis()->SetTitle("Xdet (um)");
			hCellsOverlayedMeanChargeBinAlignment.at(i)->GetYaxis()->SetTitle("Ydet (um)");
			hCellsOverlayedMeanChargeBinAlignment.at(i)->GetZaxis()->SetTitle("Charge ADC");
			hCellsOverlayedMeanChargeBinAlignment.at(i)->GetZaxis()->SetRangeUser(800,1200);
			hCellsOverlayedMeanChargeBinAlignment.at(i)->SetContour(99);
			hCellsOverlayedMeanChargeBinAlignment.at(i)->SetStats(kFALSE);
		}

		for(int i=0;i<9;i++){
			hCellsOverlayedChargeBinAlignment1.push_back(new TH2D("OverlayedCharge","OverlayedCharge",15,0,150,15,0,150));
			hCellsOverlayedChargeBinAlignment1.at(i)->SetStats(kFALSE);
			hCellsOverlayedEventsBinAlignment1.push_back(new TH2D("OverlayedEvents","OverlayedEvents",15,0,150,15,0,150));
			hCellsOverlayedEventsBinAlignment1.at(i)->SetStats(kFALSE);
			stringstream hCellsOverlayedChargeBinAlignment1Name; hCellsOverlayedChargeBinAlignment1Name<<"hCellsOverlayedMeanChargeBinAlignment1"<<FileNameEnd;
			hCellsOverlayedMeanChargeBinAlignment1.push_back(new TH2D(hCellsOverlayedChargeBinAlignment1Name.str().c_str(),hCellsOverlayedChargeBinAlignment1Name.str().c_str(),15,0,150,15,0,150));
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->GetXaxis()->SetTitle("Xdet (um)");
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->GetYaxis()->SetTitle("Ydet (um)");
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->GetZaxis()->SetTitle("Charge ADC");
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->GetZaxis()->SetRangeUser(800,1200);
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->SetContour(99);
			hCellsOverlayedMeanChargeBinAlignment1.at(i)->SetStats(kFALSE);
		}

		//hCellsColumnCheck55Name		//Cell histogram used to check whether hit is in column
		stringstream hCellsColumnCheck55Name; hCellsColumnCheck55Name<<"hCellsColumnCheck55"<<FileNameEnd;
		hCellsColumnCheck55 = new TH2D(hCellsColumnCheck55Name.str().c_str(),hCellsColumnCheck55Name.str().c_str(),30,0,150,30,0,150);

		//hCellsOverlayed55RMS
		stringstream hCellsOverlayed55RMSName; hCellsOverlayed55RMSName<<"hCellsOverlayed55RMS"<<FileNameEnd;
		hCellsOverlayed55RMS = new TH2D(hCellsOverlayed55RMSName.str().c_str(),hCellsOverlayed55RMSName.str().c_str(),30,0,150,30,0,150);

		//hCellsColumnCheck1010Name		//Cell histogram used to check whether hit is in column
		stringstream hCellsColumnCheck1010Name; hCellsColumnCheck1010Name<<"hCellsColumnCheck1010"<<FileNameEnd;
		hCellsColumnCheck1010 = new TH2D(hCellsColumnCheck1010Name.str().c_str(),hCellsColumnCheck1010Name.str().c_str(),15,0,150,15,0,150);        //30,0,150,30,0,150);

		//hCellsOverlayed1010RMS
		stringstream hCellsOverlayed1010RMSName; hCellsOverlayed1010RMSName<<"hCellsOverlayed1010RMS"<<FileNameEnd;
		hCellsOverlayed1010RMS = new TH2D(hCellsOverlayed1010RMSName.str().c_str(),hCellsOverlayed1010RMSName.str().c_str(),15,0,150,15,0,150);

		//hCellsOverlayed1010Significance
		stringstream hCellsOverlayed1010SignificanceName; hCellsOverlayed1010SignificanceName<<"hCellsOverlayed1010Significance"<<FileNameEnd;
		hCellsOverlayed1010Significance = new TH2D(hCellsOverlayed1010SignificanceName.str().c_str(),hCellsOverlayed1010SignificanceName.str().c_str(),15,0,150,15,0,150);

		//hCellsOverlayBinSpec55
		for(int i=0;i<30;i++){
			for(int j=0;j<30;j++){
				stringstream hCellsOverlayBinSpec55Name; hCellsOverlayBinSpec55Name<<"hCellsOverlayBinSpec55"<<(i*30+j)<<FileNameEnd;
				hCellsOverlayBinSpec55.push_back(new TH1F(hCellsOverlayBinSpec55Name.str().c_str(),hCellsOverlayBinSpec55Name.str().c_str(),256,0,2800));
			}
		}
		//hCellsOverlayBinSpec1010
		for(int i=0;i<15;i++){
			for(int j=0;j<15;j++){
				stringstream hCellsOverlayBinSpec1010Name; hCellsOverlayBinSpec1010Name<<"hCellsOverlayBinSpec1010"<<(i*15+j)<<FileNameEnd;
				hCellsOverlayBinSpec1010.push_back(new TH1F(hCellsOverlayBinSpec1010Name.str().c_str(),hCellsOverlayBinSpec1010Name.str().c_str(),256,0,2800));
			}
		}*/

    } //End of for ClusterSize


    if(verbosity>3)cout<<"End initialise3DCellOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DCellCentralColumnOverlayHistos() {

    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2 += appendix;
        Float_t CentralColumnXBins = settings->getCentralColumnOverlayXBins();
        Float_t CentralColumnXLow = settings->getCentralColumnOverlayXLow();
        Float_t CentralColumnXHigh = settings->getCentralColumnOverlayXHigh();
        Float_t CentralColumnYBins = settings->getCentralColumnOverlayYBins();
        Float_t CentralColumnYLow = settings->getCentralColumnOverlayYLow();
        Float_t CentralColumnYHigh = settings->getCentralColumnOverlayYHigh();

        //hCellsOverlayed
        TString name = "hCellsCentralColumnOverlayAvrgCharge";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgCharge.push_back(new TProfile2D(name,name,
                CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
                CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh));
        /*hCellsOverlayAvrgCharge = new TProfile2D(name,name,
				nXbins,xBinEdges,
				nYbins,yBinEdges);*/
        hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
        hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
        hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
        //	hCellsOverlayAvrgCharge->SetContour(99);

        //hCellsOverlayAvrgChargeMinusBadCells
        name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCells";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
        hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

        if(ClusterSize ==3){
            //hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells
            name = "hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells";
            name.Append(appendix2);
            /*hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells = (TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name);
			hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells->SetTitle(name);*/

            hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells = new TProfile2D(name,name,
                    CentralColumnXBins,0,CentralColumnXHigh-CentralColumnXLow,
                    CentralColumnYBins,0,CentralColumnYHigh-CentralColumnYLow);
            hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
            hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
            hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells->GetZaxis()->SetTitle("pulse height of cluster /ADC");
        }

        //hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut
        name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.push_back(new TH2F(name,name,
                CentralColumnXBins+2,CentralColumnXLow-2,CentralColumnXHigh+2,
                CentralColumnYBins+2,CentralColumnYLow-2,CentralColumnYHigh+2));
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetZaxis()->SetTitle("Events");

        //hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut
        name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->SetTitle(name);

        //hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis
        name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
        hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->SetTitle(name);

        name = (TString)"hLandauMinusBadCellsOffsetAnalysis" + appendix2;
        hLandauMinusBadCellsOffsetAnalysis.push_back(new TH1F(name,name,100,-50,50));
        hLandauMinusBadCellsOffsetAnalysis.back()->GetXaxis()->SetTitle("pulse height /ADC");
        hLandauMinusBadCellsOffsetAnalysis.back()->GetYaxis()->SetTitle("number of Entries #");

        //hCellsOverlayAvrgChargeGoodCells
        name = "hCellsCentralColumnOverlayAvrgChargeGoodCells";
        name.Append(appendix2);
        hCellsCentralColumnOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
        hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

        //hCellsOverlayedColumnLandau
        name = "hCellsCentralColumnOverlayLandau";
        name.Append(appendix2);
        hCellsCentralColumnOverlayLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
        hCellsCentralColumnOverlayLandau.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
        hCellsCentralColumnOverlayLandau.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

        //hCellsOverlayedLandauNoColumn
        name = "hCellsCentralColumnOverlayLandauMinusBadCells";
        name.Append(appendix2);
        hCellsCentralColumnOverlayLandauMinusBadCells.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
        hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
        hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

        //hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis
        name = "hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis";
        name.Append(appendix2);
        hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.push_back(new TH1F(name,name,256,-200,200));//todo iain ph
        hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
        hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

        name = "hCellsCentralColumnOverlayLandauGoodCells";
        name.Append(appendix2);
        hCellsCentralColumnOverlayLandauGoodCells.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
        hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
        hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

    } //End of for ClusterSize

    if(verbosity>3)cout<<"End initialise3DCellCentralColumnOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DCellBiasColumnOverlayHistos() {

    Float_t CentralColumnXBins = settings->getCentralColumnOverlayXBins();
    Float_t BiasColumnXLow = settings->getBiasColumnOverlayXLow();
    Float_t BiasColumnXHigh = settings->getBiasColumnOverlayXHigh();
    Float_t CentralColumnYBins = settings->getCentralColumnOverlayYBins();
    Float_t BiasColumnYLow = settings->getBiasColumnOverlayYLow();
    Float_t BiasColumnYHigh = settings->getBiasColumnOverlayYHigh();

    //hCellsOverlayed
    TString name = "hCellsBiasColumnOverlayAvrgCharge";
    name.Append(appendix);
    hCellsBiasColumnOverlayAvrgCharge = new TProfile2D(name,name,
            CentralColumnXBins,BiasColumnXLow,BiasColumnXHigh,
            CentralColumnYBins,BiasColumnYLow,BiasColumnYHigh);
    /*hCellsOverlayAvrgCharge = new TProfile2D(name,name,
				nXbins,xBinEdges,
				nYbins,yBinEdges);*/
    hCellsBiasColumnOverlayAvrgCharge->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
    hCellsBiasColumnOverlayAvrgCharge->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
    hCellsBiasColumnOverlayAvrgCharge->GetZaxis()->SetTitle("pulse height of cluster /ADC");
    //	hCellsOverlayAvrgCharge->SetContour(99);

    //hCellsOverlayAvrgChargeMinusBadCells
    name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCells";
    name.Append(appendix);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCells = (TProfile2D*)hCellsBiasColumnOverlayAvrgCharge->Clone(name);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCells->SetTitle(name);

    //hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells
    name = "hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells";
    name.Append(appendix);
    hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells = new TProfile2D(name,name,
            CentralColumnXBins,0,BiasColumnXHigh-BiasColumnXLow,
            CentralColumnYBins,0,BiasColumnYHigh-BiasColumnYLow);
    hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
    hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
    hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->GetZaxis()->SetTitle("pulse height of cluster /ADC");

    //hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut
    name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents";
    name.Append(appendix);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents = new TH2F(name,name,
            CentralColumnXBins+2,BiasColumnXLow-2,BiasColumnXHigh+2,
            CentralColumnYBins+2,BiasColumnYLow-2,BiasColumnYHigh+2);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->GetZaxis()->SetTitle("Events");

    //hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut
    name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut";
    name.Append(appendix);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut = (TProfile2D*)hCellsBiasColumnOverlayAvrgCharge->Clone(name);
    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->SetTitle(name);

    /*//hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis
		name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->SetTitle(name);

		//hCellsOverlayAvrgChargeGoodCells
		name = "hCellsCentralColumnOverlayAvrgChargeGoodCells";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
		hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

		//hCellsOverlayedColumnLandau
		name = "hCellsCentralColumnOverlayLandau";
		name.Append(appendix);
		hCellsCentralColumnOverlayLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
		hCellsCentralColumnOverlayLandau.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
		hCellsCentralColumnOverlayLandau.at(ClusterSize)->GetYaxis()->SetTitle("Entries");
     */
    //hCellsBiasColumnOverlayLandauMinusBadCells
    name = "hCellsBiasColumnOverlayLandauMinusBadCells";
    name.Append(appendix);
    hCellsBiasColumnOverlayLandauMinusBadCells = new TH1F(name,name,256,0,2800);//todo iain ph
    hCellsBiasColumnOverlayLandauMinusBadCells->GetXaxis()->SetTitle("charge / ADC");
    hCellsBiasColumnOverlayLandauMinusBadCells->GetYaxis()->SetTitle("Entries");

    /*//hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis
		name = "hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis";
		name.Append(appendix);
		hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.push_back(new TH1F(name,name,256,-200,200));//todo iain ph
		hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
		hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

		name = "hCellsCentralColumnOverlayLandauGoodCells";
		name.Append(appendix);
		hCellsCentralColumnOverlayLandauGoodCells.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
		hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
		hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->GetYaxis()->SetTitle("Entries");
     */
    //} //End of for ClusterSize

    cout<<"End initialise3DCellBiasColumnOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DCellOverlayIndividualBinHistos() {

    Int_t XBins = hCellsOverlayAvrgChargeMinusBadCells.at(0)->GetXaxis()->GetNbins();
    Int_t YBins = hCellsOverlayAvrgChargeMinusBadCells.at(0)->GetYaxis()->GetNbins();
    Int_t NBins = XBins*YBins;
    //Int_t hCellsOverlayAvrgChargeMinusBadCellsBins = XBins*YBins;
    Float_t BinWidth = 10; //hCellsOverlayAvrgChargeMinusBadCells->GetBin(1,1)->GetBinWidth();

    if (verbosity> 4) printf("XBins: %i, YBins: %i, BinWidth: %f \n", XBins, YBins, BinWidth);

    //hOverlayCellBinHitsBelowCut
    TString name = "hOverlayCellBinHitsBelowCut";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellBinHitsBelowCut = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellBinHitsBelowCut->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellBinHitsBelowCut->GetYaxis()->SetTitle("Entries");

    //hOverlayCellBinHits
    name = "hOverlayCellBinHits";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellBinHits = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellBinHits->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellBinHits->GetYaxis()->SetTitle("Entries");

    //hOverlayCellOffsetBinHitsBelowCut
    name = "hOverlayCellOffsetBinHitsBelowCut";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellOffsetBinHitsBelowCut = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellOffsetBinHitsBelowCut->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellOffsetBinHitsBelowCut->GetYaxis()->SetTitle("Entries");

    //hOverlayCellOffsetBinHits
    name = "hOverlayCellOffsetBinHits";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellOffsetBinHits = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellOffsetBinHits->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellOffsetBinHits->GetYaxis()->SetTitle("Entries");

    //hOverlayCellBinHitsBelowCut
    name = "hOverlayCellUnEvenBinningBinHitsBelowCut";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellUnEvenBinningBinHitsBelowCut = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellUnEvenBinningBinHitsBelowCut->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellUnEvenBinningBinHitsBelowCut->GetYaxis()->SetTitle("Entries");

    //hOverlayCellUnEvenBinningBinHits
    name = "hOverlayCellUnEvenBinningBinHits";
    //name = "hCellsCentralColumnOverlayLandau";
    name.Append(appendix);
    hOverlayCellUnEvenBinningBinHits = new TH1F(name,name,NBins+1,0,NBins);//todo iain ph
    hOverlayCellUnEvenBinningBinHits->GetXaxis()->SetTitle("OverlayCellBin");
    hOverlayCellUnEvenBinningBinHits->GetYaxis()->SetTitle("Entries");

    for(Int_t BinNumX =0; BinNumX<XBins; BinNumX++){
        for(Int_t BinNumY =0; BinNumY<YBins; BinNumY++){

            Int_t BinNum = BinNumX*YBins + BinNumY;

            Float_t XLow = BinNumX*BinWidth;
            Float_t XHigh = BinNumX*BinWidth + BinWidth;
            Float_t YLow = BinNumY*BinWidth;
            Float_t YHigh = BinNumY*BinWidth + BinWidth;

            if (verbosity>5) printf("XBin: %i, YBin: %i, BinNum: %i, XLow: %f, XHigh: %f, YLow: %f, YHigh: %f \n", BinNumX, BinNumY, BinNum, XLow, XHigh, YLow, YHigh);


            //hCellsOverlayed
            TString name = TString::Format("hCellsOverlayAvrgChargeMinusBadCells_Bin%i", BinNum);
            name.Append(appendix);
            VecOverlayCellBinHistos.push_back(new TProfile2D(name,name,
                    5,XLow,XHigh,
                    5,YLow,YHigh));
            VecOverlayCellBinHistos.at(BinNum)->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
            VecOverlayCellBinHistos.at(BinNum)->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
            VecOverlayCellBinHistos.at(BinNum)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
            //	hCellsOverlayAvrgCharge->SetContour(99);

            //hCellsOverlayedColumnLandau
            name = TString::Format("hCellsOverlayAvrgChargeMinusBadCellsLandau_Bin%i", BinNum);
            //name = "hCellsCentralColumnOverlayLandau";
            name.Append(appendix);
            VecOverlayCellBinLandaus.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
            VecOverlayCellBinLandaus.at(BinNum)->GetXaxis()->SetTitle("charge / ADC");
            VecOverlayCellBinLandaus.at(BinNum)->GetYaxis()->SetTitle("Entries");
        }
    }

    //} //End of for ClusterSize

    cout<<"End initialise3DCellOverlayIndividualBinHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DOffsetOverlayHistos() {

    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2 += appendix;

        /*With offset of 37.5, central column bin is at: 107.5 - 117.5.
    					 corner column bin is at: 32.5 - 42.5.*/

        //Double_t xBinEdges[] = {0,12.5,22.5,32.5,42.5,52.5,62.5,77.5,87.5,97.5,107.5,117.5,127.5,137.5,150};
        Double_t xBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
        //Double_t xBinEdges[] = {0,15,30,45,60,75,90,105,120,135,150};
        Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

        //Double_t yBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
        //Double_t yBinEdges[] = {0,12.5,22.5,32.5,42.5,52.5,62.5,77.5,87.5,97.5,107.5,117.5,127.5,137.5,150};
        Double_t yBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
        //Double_t yBinEdges[] = {0,15,30,45,60,75,90,105,120,135,150};
        Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;

        //hCellsOverlayed
        TString name = "hCellsOffsetOverlayAvrgCharge";
        name.Append(appendix2);
        /*hCellsCentralColumnOverlayAvrgCharge = new TProfile2D(name,name,
    		CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
    		CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh);*/
        hCellsOffsetOverlayAvrgCharge.push_back(new TProfile2D(name,name,
                nXbins,xBinEdges,
                nYbins,yBinEdges));
        hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
        hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
        hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
        //	hCellsOverlayAvrgCharge->SetContour(99);

        //hCellsOverlayAvrgChargeMinusBadCells
        name = "hCellsOffsetOverlayAvrgChargeMinusBadCells";
        name.Append(appendix);
        hCellsOffsetOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsOffsetOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

        //hCellsOverlayAvrgChargeGoodCells
        name = "hCellsOffsetOverlayAvrgChargeGoodCells";
        name.Append(appendix2);
        hCellsOffsetOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsOffsetOverlayAvrgCharge.at(0)->Clone(name));
        hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

        /*//hCellsOverlayedColumnLandau
	name = "hCellsCentralColumnOverlayLandau";
	name.Append(appendix2);
	hCellsCentralColumnOverlayLandau = new TH1F(name,name,256,0,2800);//todo iain ph
	hCellsCentralColumnOverlayLandau->GetXaxis()->SetTitle("charge / ADC");
	hCellsCentralColumnOverlayLandau->GetYaxis()->SetTitle("Entries");

    //hCellsOverlayedLandauNoColumn
    name = "hCellsCentralColumnOverlayLandauMinusBadCells";
    name.Append(appendix);
    hCellsCentralColumnOverlayLandauMinusBadCells = new TH1F(name,name,256,0,2800);//todo iain ph
    hCellsCentralColumnOverlayLandauMinusBadCells->GetXaxis()->SetTitle("charge / ADC");
    hCellsCentralColumnOverlayLandauMinusBadCells->GetYaxis()->SetTitle("Entries");

    name = "hCellsCentralColumnOverlayLandauGoodCells";
    name.Append(appendix);
    hCellsCentralColumnOverlayLandauGoodCells = new TH1F(name,name,256,0,2800);//todo iain ph
    hCellsCentralColumnOverlayLandauGoodCells->GetXaxis()->SetTitle("charge / ADC");
    hCellsCentralColumnOverlayLandauGoodCells->GetYaxis()->SetTitle("Entries");*/

    } // End of for ClusterSize

    Double_t xBinEdges[] = {0,5,15,25,35,45,55,65,70,80,90,100,110,120,130,140,150};
    Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

    Double_t yBinEdges[] = {0,5,15,25,35,45,55,65,70,80,90,100,110,120,130,140,150};
    Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;

    //hCellsOverlayed
    TString name = "hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells";
    name.Append(appendix);
    /*hCellsCentralColumnOverlayAvrgCharge = new TProfile2D(name,name,
	    		CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
	    		CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh);*/
    hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells = new TProfile2D(name,name,
            nXbins,xBinEdges,
            nYbins,yBinEdges);
    hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
    hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
    hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetZaxis()->SetTitle("pulse height of cluster /ADC");


    cout<<"End initialise3DCellCentralColumnOverlayHistos()"<<endl;
    if(verbosity>3)cout<<"End initialise3DCellCentralColumnOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DOffsetAlignmentOverlayHistos() {

    Float_t shift = 2;
    Int_t nShiftsX = 3;
    Int_t nShiftsY = 3;

    for (int i=0;i<=nShiftsX*2;i++){
        ShiftX.push_back(i*shift-nShiftsX*shift);
    }

    for (int i=0;i<=nShiftsY*2;i++){
        ShiftY.push_back(i*shift-nShiftsY*shift);
    }

    for(int i=0; i<ShiftX.size(); i++){
        Float_t OffsetX = settings->getOverlayOffsetX() + ShiftX.at(i);

        for(int j=0; j<ShiftY.size(); j++){

            Float_t OffsetY = settings->getOverlayOffsetY() + ShiftY.at(j);
            Int_t Shift = i*ShiftY.size() + j;

            TString appendix2 = TString::Format("Alignment%i_ShiftX%.1f_ShiftY%.1f",Shift,ShiftX[i],ShiftY[j]);
            appendix2 += appendix;

            Double_t xBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
            Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

            Double_t yBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
            Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;

            //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment
            TString name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment";
            name.Append(appendix2);
            /*hCellsCentralColumnOverlayAvrgCharge = new TProfile2D(name,name,
	    		CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
	    		CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh);*/
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.push_back(new TProfile2D(name,name,
                    nXbins,xBinEdges,
                    nYbins,yBinEdges));
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetXaxis()->SetTitle("rel. x position within a cell / #mum");
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetYaxis()->SetTitle("rel. y position within a cell / #mum");
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
            //	hCellsOverlayAvrgCharge->SetContour(99);

            name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut";
            name.Append(appendix);
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.push_back((TProfile2D*)hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(0)->Clone(name));
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Shift)->SetTitle(name);

            Int_t NBins = nXbins*nYbins;
            //hOverlayCellOffsetAlignmentBinHitsBelowCut			TString name = "hOverlayCellOffsetAlignmentBinHitsBelowCut";
            name = "hOverlayCellOffsetAlignmentBinHitsBelowCut";
            name.Append(appendix2);
            hOverlayCellOffsetAlignmentBinHitsBelowCut.push_back(new TH1F(name,name,NBins+1,0,NBins));//todo iain ph
            hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Shift)->GetXaxis()->SetTitle("OverlayCellBin");
            hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Shift)->GetYaxis()->SetTitle("Entries");

            //hOverlayCellOffsetBinHits
            name = "hOverlayCellOffsetAlignmentBinHits";
            //name = "hCellsCentralColumnOverlayLandau";
            name.Append(appendix);
            hOverlayCellOffsetAlignmentBinHits.push_back(new TH1F(name,name,NBins+1,0,NBins));//todo iain ph
            hOverlayCellOffsetAlignmentBinHits.at(Shift)->GetXaxis()->SetTitle("OverlayCellBin");
            hOverlayCellOffsetAlignmentBinHits.at(Shift)->GetYaxis()->SetTitle("Entries");
        }
    }
}

void TAnalysisOf3dDiamonds::initialiseTransparentAnalysisHistos() {
    hTransparentAnalysisInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hTransparentAnalysisInvalidCluster"+appendix);
    hTransparentAnalysisInvalidCluster->SetTitle("hTransparentAnalysisInvalidCluster");

    hTransparentAnalysisValidCluster = (TH2F*) hValidEventsDetSpace->Clone("hTransparentAnalysisValidCluster"+appendix);
    hTransparentAnalysisValidCluster->SetTitle("hTransparentAnalysisValidCluster");

    TString name = "hTransparentAnalysisValidClusterFidCutXvsFidCutY"+appendix;
    hTransparentAnalysisValidClusterFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetZaxis()->SetTitle("Charge ADC");

    for(int i=0; i<settings->diamondPattern.getNIntervals(); i++){
        pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
        //hLandau
        TString name = TString::Format("hLandauTransparent_pattern_%d_ch_%d_to_%d",i+1,channels.first,channels.second)+appendix;
        hLandauTransparent.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        hLandauTransparent.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
        hLandauTransparent.at(i)->GetYaxis()->SetTitle("number of entries #");
        hLandauTransparent.at(i)->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);

        //hLandauBadCellsRemoved
        name = TString::Format("hLandauTransparentBadCellsRemoved_pattern_%d_ch_%d_to_%d",i+1,channels.first,channels.second)+appendix;
        hLandauTransparentBadCellsRemoved.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        hLandauTransparentBadCellsRemoved.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
        hLandauTransparentBadCellsRemoved.at(i)->GetYaxis()->SetTitle("number of entries #");
        hLandauTransparentBadCellsRemoved.at(i)->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);

        //hXdetvsYdetvsCharge		For TH2D
        Int_t DiamondPattern = i+1;
        Float_t xLow = settings->get3dMetallisationFidCuts()->getXLow(DiamondPattern);
        Float_t xHigh = settings->get3dMetallisationFidCuts()->getXHigh(DiamondPattern);
        Float_t yLow = settings->get3dMetallisationFidCuts()->getYLow(DiamondPattern);
        Float_t yHigh = settings->get3dMetallisationFidCuts()->getYHigh(DiamondPattern);
        if(verbosity>3)cout<<"("<<xLow<<"-"<<xHigh<<", "<<yLow<<"-"<<yHigh<<")"<<endl;
        Float_t xDiv = (xHigh - xLow)/5;
        Float_t yDiv = (yHigh - yLow)/5;
        TString hXdetvsYdetvsChargeName = TString::Format("hXdetvsYdetvsCharge_%02d_%02d",channels.first,channels.second);
        hXdetvsYdetvsChargeName+= appendix + FileNameEnd;
        hXdetvsYdetvsCharge.push_back(new TH2D(hXdetvsYdetvsChargeName,hXdetvsYdetvsChargeName,xDiv,xLow,xHigh,yDiv,yLow,yHigh));
        hXdetvsYdetvsCharge.at(i)->GetXaxis()->SetTitle("X (um)");
        hXdetvsYdetvsCharge.at(i)->GetYaxis()->SetTitle("Y (um)");
        hXdetvsYdetvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        hXdetvsYdetvsEvents.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone("hXdetvsYdetvsEvents"+appendix));

        //hFidCutXvsFidCutYvsMeanCharge
        hXdetvsYdetvsMeanCharge.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone("hXdetvsYdetvsMeanCharge"+appendix));
        hXdetvsYdetvsMeanCharge.at(i)->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    }

    //hCellOverlay
    TString hCellOverlayvsChargeName ="hCellOverlayvsCharge"+appendix;
    hCellOverlayvsCharge = new TH2D(hCellOverlayvsChargeName,hCellOverlayvsChargeName,30,0,150,30,0,150);
    hCellOverlayvsCharge->GetXaxis()->SetTitle("X / #mum");
    hCellOverlayvsCharge->GetYaxis()->SetTitle("Y / #mum");
    hCellOverlayvsCharge->GetZaxis()->SetTitle("Charge ADC");

    //hCellOverlayvsEvents
    hCellOverlayvsEvents = (TH2D*)hCellOverlayvsCharge->Clone("hCellOverlayvsEvents"+appendix);

    //hCellOverlayvsMeanCharge
    hCellOverlayvsMeanCharge = (TH2D*)hCellOverlayvsCharge->Clone("hCellOverlayvsCharge"+appendix);

}

void TAnalysisOf3dDiamonds::SaveShortAnalysisHistos() {
    ShortAnalysis_SaveMeanChargeVector();
    ShortAnalysis_Save2ClusterPlots();
    vector<Float_t> xPred;
    vector<Float_t> yPred;
    vector<Float_t> charge;
    histSaver->SaveHistogram(hRelativeChargeTwoClustersX);
    histSaver->SaveHistogram(hRelativeChargeTwoClustersY);
    histSaver->SaveHistogram(hFidCutsVsMeanCharge);
    TString name = "cRelativeChargeTwoClustersXY";
    name.Append(appendix);
    TCanvas *c1 = new TCanvas(name,name);
    c1->cd();
    hRelativeChargeTwoClustersXY->Draw("colz");
    hRelativeChargeTwoClustersXY->GetZaxis()->SetRangeUser(0,20);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    histSaver->SaveCanvas(c1);
    delete c1;


    name = "cShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    c1 = new TCanvas(name,name);
    c1->cd();
    hShortAnalysis2TotalChargeXY->Draw("colz");
    hShortAnalysis2TotalChargeXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    histSaver->SaveCanvas(c1);
    delete c1;

    histSaver->SaveHistogram(hRelatviveNumberOfMultipleClusterEventsSamePattern);
    histSaver->SaveHistogram(hRelatviveNumberOfMultipleClusterEvents);
    if(hRelatviveNumberOfMultipleClusterEventsSamePattern) delete hRelatviveNumberOfMultipleClusterEventsSamePattern;
    if(hRelatviveNumberOfMultipleClusterEvents) delete hRelatviveNumberOfMultipleClusterEvents;

    name = "cTotalAvrgChargeXY";
    name.Append(appendix);
    if (settings->IsPaperMode())
        gStyle->SetCanvasDefW(gStyle->GetCanvasDefH()*2);
    c1 = new TCanvas(name,name);
    c1->cd();
    c1->SetRightMargin(.2);
    c1->SetObjectStat(false);
    hTotalAvrgChargeXY->Draw("colz");
//    hTotalAvrgChargeXY->GetZaxis()->SetTitleOffset(1.2);
    Float_t xmin = hTotalAvrgChargeXY->GetXaxis()->GetXmin();
    Float_t xmax = hTotalAvrgChargeXY->GetXaxis()->GetXmax();
    Float_t deltax = xmax - xmin;
    Float_t ymin = hTotalAvrgChargeXY->GetYaxis()->GetXmin();
    Float_t ymax = hTotalAvrgChargeXY->GetYaxis()->GetXmax();
    Float_t deltay = ymax-ymin;
    xmin = xmin -.05*deltax;
    xmax = xmax +.05*deltax;
    ymin = ymin -.05*deltay;
    ymax = ymax +.05*deltay;
    TH1F* histo = c1->DrawFrame(xmin,ymin,xmax,ymax, hTotalAvrgChargeXY->GetTitle());
    histo->GetXaxis()->SetTitle(hTotalAvrgChargeXY->GetXaxis()->GetTitle());
    histo->GetYaxis()->SetTitle(hTotalAvrgChargeXY->GetYaxis()->GetTitle());
    hTotalAvrgChargeXY->SetObjectStat(false);
    hTotalAvrgChargeXY->Draw("colz same");
    hTotalAvrgChargeXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    gPad->Update();
    c1->Update();
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    TCutG* centralRegion = settings->centralRegion3DnH->GetFiducialAreaCut();
    centralRegion->SetLineColor(kRed);
    centralRegion->Draw("same");
    histSaver->SaveCanvas(c1);

    delete c1;

    if (settings->IsPaperMode())
        gStyle->SetCanvasDefW(gStyle->GetCanvasDefH());

    for(UInt_t i = 0; i < vecEdgePredX.size(); i++){
        xPred.insert(  xPred.end(), vecEdgePredX[i].begin(), vecEdgePredX[i].end());
        yPred.insert(  yPred.end(), vecEdgePredY[i].begin(), vecEdgePredY[i].end());
        charge.insert(charge.end(), vecEdgePulseHeight[i].begin(), vecEdgePulseHeight[i].end());
    }
//    SaveEdgeFittingDistributions();
    //a.end(), b.begin(), b.end());
    histSaver->SaveHistogram(histSaver->CreateScatterHisto("hEdgeFittingPredictedPosition",yPred,xPred),false);
    TH3F* hEdgeFittingCharge = histSaver->Create3DHisto("hEdgeFittingCharge",xPred,yPred,charge);
    TH2F* hEdgeFittingAvrgCharge = (TH2F*) hEdgeFittingCharge->Project3DProfile("yx");
    hEdgeFittingAvrgCharge->SetName("hEdgeFittingAvrgCharge");
    hEdgeFittingAvrgCharge->SetTitle("hEdgeFittingAvrgCharge");
    histSaver->SaveHistogram(hEdgeFittingAvrgCharge);
    cout<<"vecEdgePredX: "<<vecEdgePredX.size()<<endl;
    if(verbosity%2==1){char t; cin>>t;}
    for(int i = 0; i < vecEdgePredX.size(); i++){
        cout<<"Edge no"<<i<<" "<<vecEdgePredX[i].size()<<" "<<vecEdgePredY[i].size()<<" "<<vecEdgePulseHeight[i].size()<<" "<<endl;
        name = "hEdgeFittingAvrgCharge_";
        name.Append(settings->getEdgePositionName(i));
        TH2F* hEdgeFittingAvrgCharge;
        if(settings->getEdgePositionType(i) == TPlaneProperties::X_COR)
            hEdgeFittingAvrgCharge = histSaver->CreateScatterHisto((string)name,vecEdgePulseHeight[i],vecEdgePredX[i],200);
        else
            hEdgeFittingAvrgCharge = histSaver->CreateScatterHisto((string)name,vecEdgePulseHeight[i],vecEdgePredY[i],200);

        hEdgeFittingAvrgCharge->GetYaxis()->SetTitle("Pulse Height /ADC");
        TString title = "predicted Position ";
        title.Append(TPlaneProperties::getCoordinateString(settings->getEdgePositionType(i)).c_str());
        title.Append(" / #mum");
        hEdgeFittingAvrgCharge->GetXaxis()->SetTitle(title);//"predicted Position X / #mum");
        histSaver->SaveHistogram(hEdgeFittingAvrgCharge);
        TH1F* hEdgeFittingAvrgCharge_pfx = (TH1F*)hEdgeFittingAvrgCharge->ProfileX();
        if(hEdgeFittingAvrgCharge_pfx){
            hEdgeFittingAvrgCharge_pfx->GetYaxis()->SetTitle("avrg. Charge / ADC");
            TCutG *cut = this->settings->getEdgePosition(i);
            name = "c";
            name.Append(hEdgeFittingAvrgCharge_pfx->GetName());
            name.Append(appendix);
            TCanvas *c1 = new TCanvas(name,name);
            hEdgeFittingAvrgCharge_pfx->Draw();
            if(cut)cut->Draw();
            histSaver->SaveCanvas(c1);
            cout<<"Saved "<<c1->GetName()<<endl;
            delete c1;
        }
        else
            cout<<" Cannot create ProfileX" << endl;
        TH1D* histo1st_py;
        TH1D* histo2nd_py;
        histSaver->SaveHistogram(hShortAnalysis2ClusterHitPattern_1stCluster);
        histSaver->SaveHistogram(hShortAnalysis2ClusterHitPattern_2ndCluster);
        for(int i = 0; i <=hShortAnalysis2ClusterHitPattern_1stCluster->GetNbinsX();i++){
            TString extension = TString::Format("_pattern%d",i);
            if (i==0)
                extension = "_all";
            name = hShortAnalysis2ClusterHitPattern_1stCluster->GetName();
            name.Append(extension);
            if(verbosity>3) cout<<name<<endl;
            if(i==0)
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name);
            else{
                //				int bin = hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name,i,i);
            }
            name = hShortAnalysis2ClusterHitPattern_2ndCluster->GetName();
            name.Append(extension);
            if(verbosity>3)cout<<name<<endl;

            if(i==0)
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name);
            else
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name,i,i);
            name = "h2ClusterAnalysis_ClusterPatterns" +appendix + extension;
            histSaver->SaveTwoHistos((string)name,histo1st_py,histo2nd_py);
            histSaver->SaveHistogram(histo1st_py);
            histSaver->SaveHistogram(histo2nd_py);
            delete histo1st_py;
            delete histo2nd_py;
        }
        name = "c2ClusterAnalysis_ClusterPatterns_px" +appendix;
        TH1D* histo1st_px = hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionX(name);
        name = "c2ClusterAnalysis_ClusterPatterns_py" +appendix;
        TH1D* histo2nd_px = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionX(name);
        name = "c2ClusterAnalysis_ClusterPatterns" +appendix;
        histSaver->SaveTwoHistos((string)name,histo1st_px,histo2nd_px);
        histSaver->SaveHistogram(histo1st_px);
        histSaver->SaveHistogram(histo2nd_px);
        if(histo1st_px) delete histo1st_px;
        if(histo2nd_px) delete histo2nd_px;
    }

    //	char t; cin>>t;
    //hNumberofClusters
    histSaver->SaveHistogram(hNumberofClusters);
    for(UInt_t i=0;i<settings->diamondPattern.getNIntervals();i++){
        hEventsvsChannelCombined->Add(hEventsvsChannel.at(i));
    }
    histSaver->SaveHistogram(hEventsvsChannelCombined);
    //	vector<TH1*> hLandauSorted;

    for(UInt_t i=0;i<settings->diamondPattern.getNIntervals();i++){

        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        name = "c_";
        name.Append(hLandau[i]->GetName());

        Float_t factor = hLandau[i]->GetBinContent(hLandau[i]->GetMaximumBin());
        factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
        name.Append("_normalized");
        histSaver->SaveTwoHistosNormalized((string)name,hLandau[i],hLandauStrip);

        Float_t max = hHitandSeedCount[i]->GetBinContent(hHitandSeedCount[i]->GetMaximumBin());
        hHitandSeedCount[i]->Scale(1./max);
        histSaver->SaveHistogram(hHitandSeedCount[i]);

        name = "c_"+(TString)hPHvsChannel[i]->GetName();
        TCanvas *c1 = new TCanvas(name,name);
        name = "h"+(TString)hPHvsChannel[i]->GetName();

        Int_t min = channels.first-1;
        max = channels.second+1;
        Int_t bins = max-min;
        TH2F* histo = new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,bins,min,max);
        c1->cd();
        histo->Draw();
        hPHvsChannel[i]->Draw("goff");
        hPHvsChannel[i]->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->GetXaxis()->SetLimits(PulseHeightMin,PulseHeightMax);
        //hPHvsChannel[i]->GetXaxis()->SetRange(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->Draw("colzsame");
        histSaver->SaveCanvas(c1);
        histSaver->SaveHistogram(hPHvsChannel[i]);
        //histSaver->SaveHistogram(hPHvsPredictedXPos.at(i));
        Float_t maxChi2 = 12;
        hChi2XChi2Y[i]->Draw("colz");
        hChi2XChi2Y[i]->GetXaxis()->SetRangeUser(0,maxChi2);
        hChi2XChi2Y[i]->GetYaxis()->SetRangeUser(0,maxChi2);
        histSaver->SaveHistogram(hChi2XChi2Y[i]);
        histSaver->SaveHistogram(hFidCutXvsFidCutY.at(i));
        //histSaver->SaveHistogram(hPHvsPredictedChannel.at(i));
        //histSaver->SaveHistogram(hFidCutXvsFidCutYvsCharge.at(i));

        //hFidCutXvsFidCutYvsMeanCharge
        ptrCanvasMean.at(i)->cd();
        *hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
        hFidCutXvsFidCutYvsMeanCharge.at(i)->SetEntries(hFidCutXvsFidCutYvsEvents.at(i)->Integral());
        hFidCutXvsFidCutYvsMeanCharge.at(i)->Draw("COLZ");
        hFidCutXvsFidCutYvsMeanCharge.at(i)->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
        TString hName  = TString::Format("cFidCutXvsFidCutYvsMeanCharge_%d_%d",channels.first,channels.second);
        hName.Append(appendix);
        ptrCanvasMean.at(i)->SetName(hName);
        histSaver->SaveCanvas(ptrCanvasMean[i]);

        /*//hXdetvsYdetvsEvents
		ptrCanvasXdetvsYdetMeanCharge.push_back(new TCanvas());
		ptrCanvasXdetvsYdetMeanCharge.at(i)->cd();
		hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
		hXdetvsYdetvsEvents.at(i)->SetEntries(hXdetvsYdetvsEvents.at(i)->Integral());
		hXdetvsYdetvsEvents.at(i)->Draw("COLZ");
		hName  = TString::Format("cXdetvsYdetMeanCharge_%d_%d",channels.first,channels.second);
		ptrCanvasXdetvsYdetMeanCharge.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvasXdetvsYdetMeanCharge[i]);
         */

    } //End of for loop
    TCanvas* cCombinedMeanCharge = new TCanvas();
    cCombinedMeanCharge->cd();

    // no Fiducial Cuts Drawn
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(0));
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(1));
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(2));
    name = "hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn";
    name.Append(appendix);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    // Selection Fiducial Cuts Drawn
    name = "hFidCutXvsFidCutYvsMeanChargeAllDetectors";
    name.Append(appendix);
    settings->getSelectionFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    // Edge F
    cCombinedMeanCharge->Clear();
    name = "hFidCutXvsFidCutYvsMeanChargeAllEdges";
    name.Append(appendix);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    settings->get3dEdgeFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge,true);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    for ( UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts(); i++){
        cCombinedMeanCharge->Clear();
        hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
        settings->get3dEdgeFidCuts()->getFidCut(i+1)->DrawFiducialCutToCanvas(cCombinedMeanCharge,true);
        TString name = "hFidCutXvsFidCutYvsMeanCharge_";
        name.Append(settings->getEdgePositionName(i));
        name.Append(appendix);
        cCombinedMeanCharge->SetName(name);
        histSaver->SaveCanvas(cCombinedMeanCharge);
    }

    for( UInt_t i=0; i < hFidCutXvsFidCutYClusters.size(); i++){
        if (hFidCutXvsFidCutYClusters[i])
            hFidCutXvsFidCutYClusters[i]->SetName(TString::Format("hFidCutXvsFidCutYClusters_%d",i));
        histSaver->SaveHistogram((TH2F*)hFidCutXvsFidCutYClusters[i]);
    }

}

void TAnalysisOf3dDiamonds::saveTransparentAnalysisHistos() {

    TCanvas* cTransparentAnalysisInvalidCluster = new TCanvas("cTransparentAnalysisInvalidCluster"+appendix, "cTransparentAnalysisInvalidCluster");
    cTransparentAnalysisInvalidCluster->cd();
    hTransparentAnalysisInvalidCluster->Draw("colz");
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisInvalidCluster);
    histSaver->SaveCanvas(cTransparentAnalysisInvalidCluster);

    TCanvas* cTransparentAnalysisValidCluster = new TCanvas("cTransparentAnalysisValidCluster"+appendix, "cTransparentAnalysisValidCluster");
    cTransparentAnalysisValidCluster->cd();
    hTransparentAnalysisValidCluster->Draw("colz");
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisValidCluster);
    histSaver->SaveCanvas(cTransparentAnalysisValidCluster);

    TFiducialCut *fidCut3dWithColumns = settings->get3dMetallisationFidCuts()->getFidCut(3);
    cout<<"3d FidCut: "<<endl;
    if (verbosity) fidCut3dWithColumns->Print(1);
    cout<<endl;
    Float_t xmin = fidCut3dWithColumns->GetXLow();
    Float_t xmax = fidCut3dWithColumns->GetXHigh();
    Float_t deltaX = TMath::Abs(.05*(xmax-xmin));
    Float_t ymin = fidCut3dWithColumns->GetYLow();
    Float_t ymax = fidCut3dWithColumns->GetYHigh();
    Float_t deltaY = TMath::Abs(.05*(ymax-ymin));

    hTransparentAnalysisValidCluster->GetXaxis()->SetRangeUser(xmin-deltaX,xmax+deltaX);
    hTransparentAnalysisValidCluster->GetYaxis()->SetRangeUser(ymin-deltaY,ymax+deltaY);
    hTransparentAnalysisValidCluster->Draw("colz");
    settings->DrawMetallisationGrid(cTransparentAnalysisValidCluster,3);
    cTransparentAnalysisValidCluster->Update();
    cTransparentAnalysisValidCluster->SetName("cTransparentAnalysisValidCluster_3DwH"+appendix);
    histSaver->SaveCanvas(cTransparentAnalysisValidCluster);
    delete cTransparentAnalysisValidCluster;

    TCanvas* cTransparentAnalysisValidClusterFidCutXvsFidCutY = new TCanvas("cTransparentAnalysisValidClusterFidCutXvsFidCutY"+appendix, "cTransparentAnalysisValidClusterFidCutXvsFidCutY");
    cTransparentAnalysisValidClusterFidCutXvsFidCutY->cd();
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->Draw();
    settings->getSelectionFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisValidClusterFidCutXvsFidCutY);
    histSaver->SaveCanvas(cTransparentAnalysisValidClusterFidCutXvsFidCutY);

}

void TAnalysisOf3dDiamonds::LongAnalysisSaveCellAndQuaterNumbering(){

    TString name = "h3DCellNumbering";
    name.Append(appendix);
    TH2D* hCellNumbering = new TH2D(name,name,
            settings->getNColumns3d(),settings->get3dMetallisationFidCuts()->getXLow(3),settings->get3dMetallisationFidCuts()->getXHigh(3),
            settings->getNRows3d(),settings->get3dMetallisationFidCuts()->getYLow(3),settings->get3dMetallisationFidCuts()->getYHigh(3));
    hCellNumbering->GetXaxis()->SetTitle("Xdet (#mum)");
    hCellNumbering->GetYaxis()->SetTitle("Ydet (#mum)");

    name = "hQuarterNumbering";
    name.Append(appendix);
    TH2D* hQuarterNumbering  = new TH2D(name,name,
            settings->getNColumns3d()*2,settings->get3dMetallisationFidCuts()->getXLow(3),settings->get3dMetallisationFidCuts()->getXHigh(3),
            settings->getNRows3d()*2,settings->get3dMetallisationFidCuts()->getYLow(3),settings->get3dMetallisationFidCuts()->getYHigh(3));
    hQuarterNumbering->GetXaxis()->SetTitle("Xdet (#mum)");
    hQuarterNumbering->GetYaxis()->SetTitle("Ydet (#mum)");
    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        for(UInt_t row=0;row<settings->getNRows3d();row++){
            Int_t cell =  settings->get3DCellNo((int)column,(int)row);
            hCellNumbering->SetBinContent(column+1,row+1,cell);
            for (UInt_t quarter = 0; quarter < settings->getNQuarters3d();quarter++){
                //				UInt_t quarterNo = settings->get3DQuarterNo(row,column,quarter);
                UInt_t quarterColumn = quarter/2;
                UInt_t quarterRow = quarter%2;
                Int_t binX = 2 * column + quarterColumn;
                Int_t binY = 2 * row + quarterRow;
                hQuarterNumbering->SetBinContent(binX+1,binY+1,cell+.1*quarter);
            }
        }
    }
    histSaver->SaveHistogramWithCellGrid(hCellNumbering,hCellNumbering);
    hCellNumbering->SetName(hCellNumbering->GetName()+(TString)"_markedCells");
    histSaver->SaveHistogramWithCellGridAndMarkedCells(hCellNumbering,hCellNumbering);
    histSaver->SaveHistogramWithCellGrid(hQuarterNumbering,hQuarterNumbering);
}

void TAnalysisOf3dDiamonds::SaveLongAnalysisHistos() {

    LongAnalysis_CompareTransparentAndClusteredAnalysis_Maps();
    //	LongAnalysis_SaveCellsOverlayMeanCharge();
    if(true){
        LongAnalysisSaveCellAndQuaterNumbering();
        LongAnalysis_SaveDeadCellProfile();
        LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors();
        LongAnalysis_SaveFailedQuarters();
        LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail();
    }
    LongAnalysis_SaveRawPulseHeightPlots();
    LongAnalysis_SaveGoodCellsLandaus();
    LongAnalysis_SaveEdgeFreeHistos();
    LongAnalysis_SaveGoodAndBadCellLandaus();
    LongAnalysis_SaveCellsOverlayMeanCharge();
    LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge();
    LongAnalysis_SaveCellsBiasColumnOverlayMeanCharge();
    LongAnalysis_SaveCellsOverlayBiasColumnAndCentralColumnStack();
    //LongAnalysis_Save3DCellOverlayIndividualBinHistos();
    LongAnalysis_Save3D3DOffsetOverlayBiasColumnAlignment();
    LongAnalysis_SaveCellsOverlayOffsetMeanCharge();
    LongAnalysis_SaveMeanChargePlots();
    LongAnalysis_CreateRelativeAddedTransparentChargeComparisonPlots();
    LongAnalysis_SaveRelativeAddedTransparentCharge();
    //LongAnalysis_SaveCellsClusterSize2DVsGrading();
    //LongAnalysis_SaveQuarterCellsClusterSize2DVsGrading();
    histSaver->SaveHistogram(hLongAnalysisInvalidCellNo);
    histSaver->SaveHistogram(hLongAnalysisInvalidCluster);
    if (!settings->do3dTransparentAnalysis())
        return;

    histSaver->SaveHistogram(hNegativeChargeFieldWireFraction);
    histSaver->SaveHistogramWithCellGrid(hNegativeChargeFieldWireFraction);
    TH1F* h = histSaver->SaveProjectionZ(hNegativeChargeFieldWireFraction,false,false,0,1.05,21);
    histSaver->SaveIntegral(h,true);
    histSaver->SaveIntegral(h,false);
    histSaver->SaveHistogram(h);
    delete h;
    TString name = hNegativeChargeFieldWireFraction->GetName();
    name+="_entries";
    TH2D* pxy = hNegativeChargeFieldWireFraction->ProjectionXY(name,"B");
    histSaver->SaveHistogram(pxy);
    histSaver->SaveHistogramWithCellGrid(pxy);
    delete pxy;
    histSaver->SaveHistogram(hNegativeChargeFieldWirePositions);
    histSaver->SaveOverlay(hNegativeChargeFieldWirePositionsOverlay);
    delete hNegativeChargeFieldWireFraction;
    delete hNegativeChargeFieldWirePositions;
    delete hNegativeChargeFieldWirePositionsOverlay;

    histSaver->SaveHistogram(hNegativeChargePosition);
    histSaver->SaveHistogram(hNegativeChargeRatio);
    TH1D* px = hNegativeChargeRatio->ProjectionX();
    histSaver->SaveHistogram(px);
    px->SetName(px->GetName()+(TString)"_logy");
    histSaver->SaveHistogram(px,false,false,true,"logy");
    delete px;
    histSaver->SaveHistogram(hNegativeChargeRatioAbs);
    histSaver->SaveHistogram(hNegativeChargeRatioMax);
    histSaver->SaveOverlay(hNegativeChargeRatioOverlay);
    histSaver->SaveNegativeChargeOverlay(hNegativeChargeRatioOverlay);
    delete hNegativeChargeRatioOverlay;
    histSaver->SaveHistogram(hNegativeChargeFraction);
    hLandauStripNegativeChargesFraction->SetTitle("Strip");
    hLandauStripNegativeChargesFraction->SetLineColor(kBlue);
    hNegativeChargeFraction->SetTitle("3D detector");
    histSaver->SaveTwoHistos("hNegativeChargeFractionComparision",hNegativeChargeFraction,hLandauStripNegativeChargesFraction);
    delete hNegativeChargeFraction;
     name = "hNegativeChargePosition"+appendix+"_px";
    TH1D* proj_px = hNegativeChargePosition->ProjectionX(name);
    histSaver->SaveHistogram(proj_px);
    name = "hNegativeChargePosition_"+appendix+"_py";
    TH1D* proj_py = hNegativeChargePosition->ProjectionY(name);
    histSaver->SaveHistogram(proj_py);
    delete proj_px;
    delete proj_py;
    TH2F* histo = (TH2F*)hNegativeChargePosition->Clone(hNegativeChargePosition->GetName() +(TString)"_grid");
    TCanvas * c1 = histSaver->DrawHistogramWithCellGrid(histo,histo);
    c1->SetName("chNegativeChargePosition"+appendix);
    histo->Draw("same colz");
    settings->DrawMetallisationGrid(c1,3);
//    histSaver->SaveCanvas(c1,c1->GetName());
    histo->SetTitle(histo->GetTitle()+(TString)" with grid");
    histSaver->SaveHistogramWithCellGrid(histo);
    histo->SetName(histo->GetName()+(TString)"_markedCells");
    histSaver->SaveHistogramWithCellGridAndMarkedCells(histo);
    delete c1;

    name = "chNegativeChargePositionGrid"+ appendix;
    c1 = new TCanvas(name,name);
    c1->cd();
    histo->Draw("colz");
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
    settings->DrawMetallisationGrid(c1,3);
    histSaver->SaveCanvas(c1);

    delete hNegativeChargePosition;
    hNegativeChargePosition = 0;
    delete histo;
    delete hNegativeChargeRatio;
    delete hNegativeChargeRatioAbs;
    histo =0;
    //histSaver->SaveHistogram(hNegativeChargePosition);
}

void TAnalysisOf3dDiamonds::LongAnalysis_InitChargeSharingPlots(){
    UInt_t nCells = settings->GetNCells3d();
    UInt_t nBinsX = 256;
    UInt_t nBinsY = 256;
    Float_t minX = 0;
    Float_t maxX = PulseHeightMax;
    Float_t minY = 0;
    Float_t maxY = PulseHeightMax/2.;
    for (UInt_t cell = 0; cell <nCells;cell++){
        TString name = TString::Format("hChargeSharing_CellNo_%02d",cell)+appendix;
        TString title = TString::Format("hChargeSharing_CellNo_%02d",cell);;
        TH2F* histo = new TH2F(name,title,nBinsX,minX,maxX,nBinsY,minY,maxY);
        histo->GetXaxis()->SetTitle("Charge highest signal / ADC");
        histo->GetYaxis()->SetTitle("Charge highest adjancent signal / ADC");
        vecHChargeSharing.push_back(histo);
    }
}


void TAnalysisOf3dDiamonds::LongAnalysis_FillChargeSharingPlots(){
    UInt_t cellNo = settings->getCellNo(xPredDet,yPredDet);
    Int_t highest_hit_pos = diamondCluster->getHighestHitClusterPosition();
    Int_t Second_highest_hit_pos = diamondCluster->getHighestSignalNeighbourClusterPosition(highest_hit_pos,useCMN,true);
    Float_t highestSignal = diamondCluster->getSignal(highest_hit_pos,useCMN);
    Float_t adjacentSignal = diamondCluster->getSignal(Second_highest_hit_pos,useCMN);
    vecHChargeSharing.at(cellNo)->Fill(highestSignal,adjacentSignal);
}
//LongAnalysis_SaveChargeSharingPlots
void TAnalysisOf3dDiamonds::LongAnalysis_SaveChargeSharingPlots(){
    LongAnalysis_CreateTH2_CellPlots(&vecHChargeSharing,"","hChargeSharing");
}

void TAnalysisOf3dDiamonds::LongAnalysis_InitResolutionPlots(){
    UInt_t nCells = settings->GetNCells3d();
    UInt_t nBins = 128;
    Float_t minX = - 1*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX = 1*settings->GetCellWidth(subjectDetector,2);
    TString name = "hAdjacentSNR_vs_cellNo"+appendix;
    TString title = "hAdjacentSNR_vs_cellNo"+appendix;
    title+=";Cell No;SNR adjacent Strip";
    hAdjacentSNR_vs_cellNo = new TH2F(name,title,nCells,0,nCells,160,-30,50);
    name = "hAdjacentChannels_Signal"+appendix;
    title = "hAdjacentChannels_Signal"+appendix;
    title+=";Signal left Strip / ADC No;Signal right Strip / ADC";
    hAdjacentChannels_Signal = new TH2F(name,title,300,-300,300,300,-300,300);
    for (UInt_t cell = 0; cell <nCells;cell++){
        name = TString::Format("hResolution_CellNo_%02d_maxValue",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_maxValue",cell);;
        TH1F* histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_maxValue.push_back(histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_chargeWeighted",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_chargeWeighted.push_back(histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell)+appendix;
        title = TString::Format("hResolution_CellNo_%02d_highest2Centroid",cell);;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_highest2Centroid.push_back(histo);
        /*******/
        name = TString::Format("hResolution_CellNo_%02d_h2C_withCut",cell)+appendix;
        title = TString::Format("hResolution Cell %02d - h2C with SNR Cut: %2.1f",cell,settings->GetResolutionSNR());;
        histo = new TH1F(name,title,nBins,minX,maxX);
        histo->GetXaxis()->SetTitle("Residual / #mum");
        vecHResolutionPerCell_h2C_WithCut.push_back(histo);

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

void TAnalysisOf3dDiamonds::LongAnalysis_FillResolutionPlots(){

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
    if (TMath::Abs(relPredPos-relPred.first) > 5)
        cout<<"\n "<<pos_max<<" "<<predPos<<" "<<delta_max<<" "<<delta_Weigthed<<" "<<delta_h2C<< " "<<relPredPos<<" "<<relPredPos.first<<"/"<<relPredPos.second<<endl;
    //diamondCluster->Print(1);

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
        TH2F* histo  = vecHResolutionPerCell_maxValue_vs_SNR.at(cellNo);
        if (histo)
            histo->Fill(delta_max*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_maxValue_vs_PredHit.size()){
        TH2F* histo  = vecHResolutionPerCell_maxValue_vs_PredHit.at(cellNo);
        if (histo)
            histo->Fill(delta_max*cellWidth,relPredPos);
    }
    if (cellNo< vecHResolutionPerCell_maxValue_vs_PredHitY.size()){
        TH2F* histo  = vecHResolutionPerCell_maxValue_vs_PredHitY.at(cellNo);
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
        TH2F* histo  = vecHResolutionPerCell_chargeWeighted_vs_SNR.at(cellNo);
        if (histo)
            histo->Fill(delta_Weigthed*cellWidth,snr);
    }

    if (cellNo< vecHResolutionPerCell_chargeWeighted_vs_PredHit.size()){
        TH2F* histo  = vecHResolutionPerCell_chargeWeighted_vs_PredHit.at(cellNo);
        //cout<<"FILL vecHResolutionPerCell_chargeWeighted_vs_PredHit:"<<cellNo<<"\t"<<relPredPos<<" --> "<<delta*cellWidth<<endl;
        if (histo)
            histo->Fill(delta_Weigthed*cellWidth,relPredPos);
    }
    else{
        cout<<"Cannot find "<<cellNo<< " in "<<vecHResolutionPerCell_chargeWeighted_vs_PredHit.size()<<endl;
    }

    if (cellNo< vecHResolutionPerCell_chargeWeighted_vs_PredHitY.size()){
        TH2F* histo  = vecHResolutionPerCell_chargeWeighted_vs_PredHitY.at(cellNo);
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
        TH2F* histo  = vecHResolutionPerCell_highest2Centroid_vs_SNR.at(cellNo);
        if (histo) histo->Fill(delta_h2C*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_highest2Centroid_vs_PredHit.size()){
        TH2F* histo  = vecHResolutionPerCell_highest2Centroid_vs_PredHit.at(cellNo);
        if (histo)  histo->Fill(delta_h2C*cellWidth,relPredPos);
    }
    if (cellNo< vecHResolutionPerCell_highest2Centroid_vs_PredHitY.size()){
        TH2F* histo  = vecHResolutionPerCell_highest2Centroid_vs_PredHitY.at(cellNo);
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
        TH2F* histo  = vecHResolutionPerCell_h2C_WithCut_vs_SNR.at(cellNo);
        if (histo) histo->Fill(delta*cellWidth,snr);
    }
    if (cellNo< vecHResolutionPerCell_h2C_WithCut_vs_PredHit.size()){
        TH2F* histo  = vecHResolutionPerCell_h2C_WithCut_vs_PredHit.at(cellNo);
        if (histo)  histo->Fill(delta*cellWidth,relPredPos);
    }

    if (cellNo< vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.size()){
        TH2F* histo  = vecHResolutionPerCell_h2C_WithCut_vs_PredHitY.at(cellNo);
        if (histo)  histo->Fill(delta*cellWidth,relPredPosY);
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateTH2_CellPlots(vector<TH2F*>*vec,TString kind, TString prefix){
    if (!vec) return;
    if (vec->size() == 0) return;
    TH2F* histo = vec->at(0);
    if (!histo) return;
    if (kind != "" && !kind.BeginsWith("_"))
        kind.Prepend("_");

    TString name = prefix+"GoodCells"+kind+appendix;
    TH2F* hGoodCells = (TH2F*)histo->Clone(name);
    hGoodCells->Reset();
    hGoodCells->SetTitle(prefix + " Good Cells "+ kind +appendix);

    name = prefix+"BadCells"+kind+appendix;
    TH2F* hBadCells =  (TH2F*)hGoodCells->Clone(name);
    hBadCells->Reset();
    hBadCells->SetTitle(prefix + " Bad Cells "+ kind +appendix);

    name = prefix+"AllButBadCells"+kind+appendix;
    TH2F* hAllButBadCells =  (TH2F*)hGoodCells->Clone(name);
    hAllButBadCells->Reset();
    hAllButBadCells->SetTitle(prefix + " AllButBad Cells "+ kind +appendix);

    name = prefix+"AllCells"+kind+appendix;
    TH2F* hAllCells =  (TH2F*)hGoodCells->Clone(name);
    hAllCells->Reset();
    hAllCells->SetTitle(prefix + " All Cells "+ kind +appendix);

    string plots_path = histSaver->GetPlotsPath();
    string new_plots_path = plots_path;
    new_plots_path+=(string)prefix;
    new_plots_path+="/";
    HistogrammSaver newHistSaver(settings);
    newHistSaver.SetPlotsPath(new_plots_path);
    for(UInt_t cell=0;cell< vec->size();cell++){
        TH2F* histo = vec->at(cell);
        if (!histo)
            continue;
        if (settings->IsGoodCell(3,cell))
            hGoodCells->Add(histo);
        if (settings->isBadCell(3,cell))
            hBadCells->Add(histo);
        else
            hAllButBadCells->Add(histo);
        hAllCells->Add(histo);
        newHistSaver.SaveHistogram(histo,true,false);
        vec->at(cell)= 0;
        delete histo;
    }
    hGoodCells->GetZaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hGoodCells);//,false,false,true);
    hBadCells->GetZaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hBadCells);//,false,false,true);
    hAllCells->GetZaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hAllCells);
    hAllButBadCells->GetZaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hAllButBadCells);//,false,false,true);
    delete hGoodCells;
    delete hBadCells;
    delete hAllCells;
    delete hAllButBadCells;
}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateResolutionPlots(vector<TH1F*>*vec,TString kind){
    UInt_t nBins = 128;
    Float_t minX = -1*settings->GetCellWidth(subjectDetector,2);
    Float_t maxX =settings->GetCellWidth(subjectDetector,2);
    TString name = "hResolutionGoodCells_"+kind+appendix;
    TH1F* hResolutionGoodCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionGoodCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionBadCells_"+kind+appendix;
    TH1F* hResolutionBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionBadCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllCells_"+kind+appendix;
    TH1F* hResolutionAllCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllCells->GetXaxis()->SetTitle("Residual / #mum");
    name = "hResolutionAllButBadCells_"+kind+appendix;
    TH1F* hResolutionAllButBadCells = new TH1F(name,name,nBins,minX,maxX);
    hResolutionAllButBadCells->GetXaxis()->SetTitle("Residual / #mum");
    string plots_path = histSaver->GetPlotsPath();
    HistogrammSaver newHistSaver(settings);
    newHistSaver.SetPlotsPath(plots_path+(string)"/resolution/");
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
        newHistSaver.SaveHistogram(histo);
        vec->at(cell)= 0;
        delete histo;
    }
    Float_t xmin = -100;
    Float_t xmax = +100;
    TF1* fitX = new TF1("fit","[0]*TMath::Sqrt(TMath::Pi()/2)*[1]*(TMath::Erf(([2]+[3]-x)/TMath::Sqrt(2)/[1])+TMath::Erf(([3]-[2]+x)/TMath::Sqrt(2)/[1]))",xmin,xmax);
    fitX->FixParameter(3,settings->GetCellWidth(subjectDetector,2)/2);//TODO
    fitX->SetParLimits(1,0,40);
    fitX->SetParNames("Integral","sigma of Gaus","position");
    fitX->SetParameter(2,0.);
    fitX->SetParameter(1,10);
    Int_t statOpt = gStyle->GetOptStat();
    gStyle->SetOptStat(1111);
    hResolutionGoodCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionGoodCells,false,false,true);
    hResolutionBadCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionBadCells,false,false,true);
    hResolutionAllCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionAllCells);
    hResolutionAllButBadCells->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hResolutionAllButBadCells,false,false,true);
    gStyle->SetOptStat(statOpt);
    delete hResolutionGoodCells;
    delete hResolutionBadCells;
    delete hResolutionAllCells;
    delete hResolutionAllButBadCells;
}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateResolutionPlots(){

    if (!settings->do3dTransparentAnalysis())
        return;
    LongAnalysis_CreateResolutionPlots(&vecHResolutionPerCell_chargeWeighted,"chargeWeighted");
    LongAnalysis_CreateResolutionPlots(&vecHResolutionPerCell_highest2Centroid,"highest2Centroid");
    LongAnalysis_CreateResolutionPlots(&vecHResolutionPerCell_maxValue,"maxValue");
    LongAnalysis_CreateResolutionPlots(&vecHResolutionPerCell_h2C_WithCut,"h2C_WithCut");

    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_SNR,"maxValue_SNR");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_SNR,"chargeWeighted_SNR");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_SNR,"highest2Centroid_SNR");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_SNR,"h2C_WithCut_SNR");

    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_PredHit,"maxValue_PredHit");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_PredHit,"chargeWeighted_PredHit");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_PredHit,"highest2Centroid_PredHit");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_PredHit,"h2C_WithCut_PredHit");


    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_maxValue_vs_PredHitY,"maxValue_PredHitY");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_chargeWeighted_vs_PredHitY,"chargeWeighted_PredHitY");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_highest2Centroid_vs_PredHitY,"highest2Centroid_PredHitY");
    LongAnalysis_CreateTH2_CellPlots(&vecHResolutionPerCell_h2C_WithCut_vs_PredHitY,"h2C_WithCut_PredHitY");


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
    LongAnalysis_SaveSNRPerCell();
    delete hAdjacentSNR_vs_cellNo;
    delete hAdjacentChannels_SNR;
    delete hAdjacentChannels_Signal;
    hAdjacentChannels_SNR=0;
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveSNRPerCell(){

    TH2D* h = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Clone();
    UInt_t xBins = h->GetXaxis()->GetNbins();
    UInt_t yBins = h->GetYaxis()->GetNbins();
    UInt_t xRebin = xBins/settings->getNColumns3d()/2;
    UInt_t yRebin = yBins/settings->getNRows3d()/2;
    TString name = "hPulseHeightVsDetectorHitPostionXY_perCell";
    h = (TH2D*)h->Rebin2D(xRebin*2,yRebin*2,name);
    TH2D* hNegativeSNRs = histSaver->GetHistoBinedInCells("hNegativeSNRs"+appendix,1);
    TH2D* hNegativeSNRsRelative = histSaver->GetHistoBinedInCells("hNegativeSNRsRelative"+appendix,1);
    if (h->GetNbinsX()!=hNegativeSNRs->GetNbinsX() ||
            h->GetNbinsY()!=hNegativeSNRs->GetNbinsY()){
        cout<<"ERROR; Bin Numbers do NOT Agree:"<<
                hNegativeSNRs->GetNbinsX()<<"/"<<hNegativeSNRs->GetNbinsY()<<"\t"<<
                h->GetNbinsX()<<"/"<<h->GetNbinsY()<<endl;
        char t;
        cin >>t;
    }
    vector<Float_t> nNegativeSNRs, rNegativeSNRs, fSignal;
    for (UInt_t i = 1; i <= hAdjacentSNR_vs_cellNo->GetNbinsX();i++){
        TH1D* p = hAdjacentSNR_vs_cellNo->ProjectionY("",i,i);
        Int_t maxbin = p->GetXaxis()->FindBin(0.);
        Float_t nNegatives = p->Integral(0,maxbin);
        Float_t entries = p->GetEntries();
        Int_t row = settings->getRowOfCell(i-1);
        Int_t column = settings->getColumnOfCell(i-1);
        hNegativeSNRs->SetBinContent(column+1,row+1,nNegatives);
        nNegativeSNRs.push_back(nNegatives);
        rNegativeSNRs.push_back(nNegatives/entries*100.);
        fSignal.push_back(h->GetBinContent(column+1,row+1));
        hNegativeSNRsRelative->SetBinContent(column+1,row+1,nNegatives/entries*100.);
        cout<<TString::Format("%3d -> %2d/%2d | %5.0f | %5.2f ",i,column,row,nNegatives,nNegatives/entries*100.)<<endl;
        delete p;
    }
    TGraph g = histSaver->CreateDipendencyGraph("gNegativeVsAvrgSignal",nNegativeSNRs,fSignal);
    g.SetMarkerStyle(5);
    histSaver->SaveGraph(&g,"gNegativeVsAvrgSignal");
    g.SetMarkerStyle(5);
    g = histSaver->CreateDipendencyGraph("gRelNegativeVsAvrgSignal",rNegativeSNRs,fSignal);
    g.SetMarkerStyle(5);
    histSaver->SaveGraph(&g,"gRelNegativeVsAvrgSignal");
    std::vector<Float_t>::iterator result = std::min_element(nNegativeSNRs.begin(), nNegativeSNRs.end());
    Float_t minX = *result;
    result = std::max_element(nNegativeSNRs.begin(), nNegativeSNRs.end());
    Float_t maxX = *result;
    result = std::min_element(rNegativeSNRs.begin(), rNegativeSNRs.end());
    Float_t minX2 = *result;
    result = std::max_element(rNegativeSNRs.begin(), rNegativeSNRs.end());
    Float_t maxX2 = *result;
    result = std::min_element(fSignal.begin(), fSignal.end());
    Float_t minY = *result;
    result = std::max_element(fSignal.begin(), fSignal.end());
    Float_t maxY = *result;
    cout<<"nNegativeSNR: "<<minX<<"-"<<maxX<<endl;
    cout<<"rNegativeSNR: "<<minX2<<"-"<<maxX2<<endl;
    cout<<"fSignal: "<<minY<<"-"<<maxY<<endl;
    Int_t xbins = 20;
    Int_t ybins = 20;
    Float_t xmin = minX - .1 * (maxX - minX);
    Float_t xmax = maxX + .1 * (maxX - minX);
    Float_t ymin = minY - .1 * (maxY - minY);
    Float_t ymax = maxY + .1 * (maxY - minY);
    name = "hNegativeVsAvrgSignal";
    TH2D* hNegativeVsAvrgSignal = new TH2D(name,name,xbins,xmin,xmax,ybins,ymin,ymax);
    xmin = minX2 - .1 * (maxX2 - minX2);
    xmax = maxX2 + .1 * (maxX2 - minX2);
    name = "hRelNegativeVsAvrgSignal";
    TH2D* hRelNegativeVsAvrgSignal = new TH2D(name,name,xbins,xmin,xmax,ybins,ymin,ymax);
    for (UInt_t i = 0; i < fSignal.size();i++){
        hNegativeVsAvrgSignal->Fill(nNegativeSNRs.at(i),fSignal.at(i));
        hRelNegativeVsAvrgSignal->Fill(rNegativeSNRs.at(i),fSignal.at(i));
    }
    hNegativeVsAvrgSignal->GetYaxis()->SetTitle("Avrg Signal of Cell");
    hNegativeVsAvrgSignal->GetXaxis()->SetTitle("No of adjacent negative SNRs in Cell");
    histSaver->SaveHistogram(hNegativeVsAvrgSignal);
    hRelNegativeVsAvrgSignal->GetYaxis()->SetTitle("Avrg Signal of Cell");
    hRelNegativeVsAvrgSignal->GetXaxis()->SetTitle("rel. no. of adjacent negative SNRs in Cell");
    histSaver->SaveHistogram(hRelNegativeVsAvrgSignal);

    histSaver->SaveHistogramWithCellGrid(hNegativeSNRs,hNegativeSNRs);
    histSaver->SaveHistogramWithCellGrid(hNegativeSNRsRelative,hNegativeSNRsRelative);
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveRawPulseHeightPlots(){
    TH1F* hStrip = (TH1F*)hLandauStrip->Clone();
    hStrip->SetLineStyle(2);
    histSaver->SaveHistogram(hLandau3DWithColumns);
    histSaver->SaveHistogram(hLandau3DPhantom);
    histSaver->SaveHistogram(hLandau3DPhantomCentral);
    histSaver->SaveHistogram(hLandau3DWithColumnsFidCutXvsFidCutY);
    histSaver->SaveHistogram(hLandau3DPhantomFidCutXvsFidCutY);

    hLandau3DWithColumns->SetTitle("3D");
    hLandau3DPhantom->SetTitle("3D Phantom");
    hLandau3DPhantomCentral->SetTitle("3D Phantom, central Region");
    TString name = "sAllPulseHeigthDistributions";
    name.Append(appendix);
    THStack sAllPulseHeigthDistributions(name,name);
    sAllPulseHeigthDistributions.Add(hLandauStrip);
    sAllPulseHeigthDistributions.Add(hLandau3DWithColumns);
    sAllPulseHeigthDistributions.Add(hLandau3DPhantom);
    sAllPulseHeigthDistributions.Add(hLandau3DPhantomCentral);
    histSaver->SaveStack(&sAllPulseHeigthDistributions,"nostack",true,false);


    Float_t max = hLandau3DWithColumns->GetBinContent(hLandau3DWithColumns->GetMaximumBin());
    hLandau3DWithColumns->Scale(1./max);

    max = hLandau3DPhantom->GetBinContent(hLandau3DPhantom->GetMaximumBin());
    hLandau3DPhantom->Scale(1./max);

    max = hLandau3DPhantomCentral->GetBinContent(hLandau3DPhantomCentral->GetMaximumBin());
    hLandau3DPhantomCentral->Scale(1./max);

    max = hStrip->GetBinContent(hStrip->GetMaximumBin());
    hStrip->Scale(1./max);
    hStrip->SetTitle("Strip");

    name = "sAllPulseHeigthDistributions_scaled";
    name.Append(appendix);
    THStack sAllPulseHeigthDistributions_normalized(name,name);
    sAllPulseHeigthDistributions_normalized.Add(hStrip);
    sAllPulseHeigthDistributions_normalized.Add(hLandau3DWithColumns);
    sAllPulseHeigthDistributions_normalized.Add(hLandau3DPhantom);
    sAllPulseHeigthDistributions_normalized.Add(hLandau3DPhantomCentral);
    sAllPulseHeigthDistributions_normalized.Draw("");
    gPad->Update();;
    if(sAllPulseHeigthDistributions_normalized.GetXaxis()) sAllPulseHeigthDistributions_normalized.GetXaxis()->SetTitle("charge / ADC");
    if(sAllPulseHeigthDistributions_normalized.GetYaxis()) sAllPulseHeigthDistributions_normalized.GetYaxis()->SetTitle("entries a.u.");
    histSaver->SaveStack(&sAllPulseHeigthDistributions_normalized,"nostack",true,false,"charge / ADC","entries a.u.");

    name = "ccAllPulseHeigthDistributions"+appendix;
    TCanvas *c1 = new TCanvas(name,name);
    c1->SetObjectStat(false);
    hLandau3DWithColumns->SetObjectStat(false);
    hLandau3DWithColumns->Draw();
    hLandau3DPhantom->Draw("same");
    hLandau3DPhantomCentral->Draw("same");
    hStrip->Draw("same");
    TLegend* leg = c1->BuildLegend();
    leg->SetFillColor(kWhite);
    histSaver->SaveCanvas(c1);

//    histSaver->SaveCanvas(c1);

    delete hStrip;
    delete hLandau3DWithColumns;
    delete hLandau3DPhantom;
    delete hLandau3DPhantomCentral;
    delete hLandau3DWithColumnsFidCutXvsFidCutY;
    delete hLandau3DPhantomFidCutXvsFidCutY;
}
vector<Float_t> TAnalysisOf3dDiamonds::LongAnalysis_GradeCellByQuarters(int quarterFailCriteriaTyp, vector<TH1F*> hQuarterLandaus){
    vector<Float_t> vecFluctuations;
    Float_t compareQuarterTo;
    //	if (true){
    vector<TH1F*> hQuarterCellsLandausSorted = hQuarterLandaus;
    sort(hQuarterCellsLandausSorted.begin(), hQuarterCellsLandausSorted.end(), TSettings::SorterForPulseHeightOfHisto);
    compareQuarterTo = hQuarterCellsLandausSorted.at(0)->GetMean();
    //	}
    for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){

        if(verbosity>3)cout<<"hQuarterCellsLandausSorted.at(quarter): "<<hQuarterCellsLandausSorted.at(quarter)->GetMean()<<"\t";
        Float_t QuarterMean = hQuarterLandaus[quarter]->GetMean();
        int entries = hQuarterLandaus[quarter]->GetEntries();
        Float_t fluctuation = (compareQuarterTo - QuarterMean)/compareQuarterTo;
        if(verbosity>3)cout<<TString::Format("\t%d: %.1f/%.1f with %3d (%2.1f%%)",
                quarter,QuarterMean,compareQuarterTo,entries,fluctuation)
        <<endl;
        vecFluctuations.push_back(fluctuation);
    }
    return vecFluctuations;

}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors(int quarterFailCriteriaTyp) {

    TString name = TString::Format("hQuarterFailCriteria_%d",quarterFailCriteriaTyp)+appendix;
    hLongAnalysisQuarterFluctuations = histSaver->GetHistoBinedInCells(name,2);
    cout<<"LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors"<<endl;
    cout<<" Criteria: "<<quarterFailCriteriaTyp<<endl;
    cout << " Mean of Landau - Strip:      " << (hLandauStrip->GetMean()) <<endl;
    vecQuarterCellsPassFail.clear();
    vecQuarterCellsPassFail.resize(settings->GetNCells3d());
    vector<vector<Float_t> > vecQuarterCellsFluctuation;
    vecQuarterCellsFluctuation.resize(settings->GetNCells3d());;
    Float_t DeadCellThreshold = 700;
    for (UInt_t cell = 0; cell < settings->GetNCells3d();cell++){
        Int_t Grading = 0;
        vecQuarterCellsFluctuation.at(cell) = LongAnalysis_GradeCellByQuarters(quarterFailCriteriaTyp, hQuarterCellsLandau[cell]);
        cout<<vecQuarterCellsFluctuation[cell].size()<<endl;

        if(LongAnalysis_IsDeadCell(hQuarterCellsLandau[cell], DeadCellThreshold)){
            vecDeadCells.push_back(cell);
            CellGrading.push_back(hQuarterCellsLandau[cell].size());
        }
        else{
            for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){
                Float_t Fluctuation = vecQuarterCellsFluctuation[cell][quarter];
                int row = settings->getRowOfCell(cell);
                int column = settings->getColumnOfCell(cell);
                cout<<"fill: "<<column<<"/"<<row<<" "<<quarter<<" "<<Fluctuation<<"\t\t"<<endl;
                Float_t FluctuationFactor = 0.1; //Needs to be added to settings file.
                bool isFailedQuarter = TMath::Abs(Fluctuation)>FluctuationFactor;
                cout<<TString::Format("%d - %d: %2.1f ==> %d",cell,quarter,100.*Fluctuation,isFailedQuarter);
                //vecQuarterCellsFluctuation[cell].push_back(Fluctuation);
                if(isFailedQuarter)
                    vecQuarterCellsPassFail[cell].push_back(1);
                else vecQuarterCellsPassFail[cell].push_back(0);
                Grading += vecQuarterCellsPassFail[cell].at(quarter);

                UInt_t quarterColumn = quarter/2;
                UInt_t quarterRow = quarter%2;
                Int_t binX = 2 * column + quarterColumn;
                Int_t binY = 2 * row + quarterRow;
                hLongAnalysisQuarterFluctuations->SetBinContent(binX+1,binY+1,Fluctuation);
            }
            cout<<"Cell: "<<cell<<"  Grading: "<<Grading<<endl;
            CellGrading.push_back(Grading);
        }
        if(settings->isBadCell(3,cell))
            HighlightGrading.push_back(2);
        else
            HighlightGrading.push_back(0);
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveEdgeFreeHistos() {
    hEventsCentralRegion->Draw("colz");
    hEventsCentralRegion->GetXaxis()->SetRangeUser(0,150);
    hEventsCentralRegion->GetYaxis()->SetRangeUser(0,150);
    histSaver->SaveHistogram(hEventsCentralRegion,true,false);
    histSaver->SaveHistogram(hEventsEdgeRegion,true,false);

    hPulseHeigthCentralRegion->Sumw2();
    hPulseHeigthEdgeRegion->Sumw2();

    histSaver->SaveHistogramWithCellGrid(hPulseHeigthCentralRegion);
    histSaver->SaveHistogramWithCellGrid(hPulseHeigthEdgeRegion);

    TProfile2D* hCompare;
    TString name = "hPulseHeigthCompareRegion";
    name.Append(appendix);
    hCompare = (TProfile2D*)hPulseHeigthCentralRegion->Clone(name);
    if(hCompare){
        hCompare->Draw("goffcolz");
        hCompare->SetTitle("comparing pulse height of edge and central region");
        hCompare->GetZaxis()->SetTitle("avrg ph_{central} / avrg. ph_{edge}");
        hCompare->Divide(hPulseHeigthEdgeRegion);
        histSaver->SaveHistogramWithCellGrid(hCompare);
        delete hCompare;
    }

    TProfile2D* hPulseHeigthCentralRegionCell = 0;
    TProfile2D* hPulseHeigthEdgeRegionCell = 0;
    if(hPulseHeigthCentralRegion){
        TString name = (TString)hPulseHeigthCentralRegion->GetName()+"_Cell"+appendix;
        hPulseHeigthCentralRegionCell = (TProfile2D*)hPulseHeigthCentralRegion->Rebin2D(2,2,name);
        hPulseHeigthCentralRegionCell->Sumw2();
        histSaver->SaveHistogramWithCellGrid(hPulseHeigthCentralRegionCell);
    }
    if(hPulseHeigthEdgeRegion){
        TString name =(TString)hPulseHeigthEdgeRegion->GetName()+"_Cell"+appendix;
        hPulseHeigthEdgeRegionCell = (TProfile2D*)hPulseHeigthEdgeRegion->Rebin2D(2,2,name);
        hPulseHeigthEdgeRegionCell->Sumw2();
        histSaver->SaveHistogramWithCellGrid(hPulseHeigthEdgeRegionCell);
    }
    if(hPulseHeigthEdgeRegionCell && hPulseHeigthCentralRegionCell){
        TString name = "hPulseHeigthCompareRegion_Cell"+appendix;
        hCompare = (TProfile2D*)hPulseHeigthCentralRegionCell->Clone(name);
        if(hCompare){
            hCompare->Draw("goffcolz");
            hCompare->SetTitle("comparing pulse height of edge and central region");
            hCompare->GetZaxis()->SetTitle("avrg ph_{central} / avrg. ph_{edge}");
            hCompare->Divide(hPulseHeigthEdgeRegionCell);
            histSaver->SaveHistogramWithCellGrid(hCompare);
            delete hCompare;
        }
    }
}

/**
 * Checks if all 4 Quarters are bellow a certain Threshold, this is an indicator for a dead cell, e.t. missing/broken readout column
 * @param nhQuarterCellsLandau
 * @param nThreshold
 * @return
 */
bool TAnalysisOf3dDiamonds::LongAnalysis_IsDeadCell(vector<TH1F*> nhQuarterCellsLandau, Float_t nThreshold){

    UInt_t NQuarterLow =0;
    for(UInt_t quarter=0;quarter<nhQuarterCellsLandau.size();quarter++){
        Float_t QuarterMean = nhQuarterCellsLandau.at(quarter)->GetMean();
        if(QuarterMean < nThreshold)
            NQuarterLow++;
        cout<<"QuarterMean: "<<QuarterMean<<" NQuarterLow: "<<NQuarterLow<<endl;
    }
    if(NQuarterLow == nhQuarterCellsLandau.size()){
        cout<<"Return True"<<endl;
        return true;
    }
    else{
        return false;
        cout<<"Return True"<<endl;
    }
    /*Float_t MeanSum = 0;
	for(UInt_t quarter=0;quarter<nhQuarterCellsLandau.size();quarter++){
		MeanSum += nhQuarterCellsLandau.at(quarter)->GetMean();
	}*/
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveFailedQuarters(){
    //	cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveFailedQuarters]"<<flush;
    vector < pair<Int_t,Int_t> > failedQuarters;
    for (UInt_t cell = 0; cell < vecQuarterCellsPassFail.size();cell++){
        for(UInt_t quarter = 0; quarter< vecQuarterCellsPassFail[cell].size(); quarter++){
            if ( vecQuarterCellsPassFail[cell][quarter])
                failedQuarters.push_back(make_pair((Int_t)cell,(Int_t)quarter));
        }
    }
    //	cout<<" with "<<failedQuarters.size()<<" failed quarters"<<endl;
    TH2F* histo = new TH2F();
    histo->SetName("hFailedQuarteters");
    histo->SetTitle("Failed Quarters");
    TCanvas *c1 = histSaver->DrawHistogramWithCellGrid(histo);
    histSaver->DrawFailedQuarters(failedQuarters,c1);
    histSaver->SaveCanvas(c1);
    delete c1;
    cout<<"Draw HistogramWithCellGrid"<<endl;
    c1 = histSaver->DrawHistogramWithCellGrid(hLongAnalysisQuarterFluctuations,hLongAnalysisQuarterFluctuations);
    cout<<"Draw DrawFailedQuarters"<<endl;
    histSaver->DrawFailedQuarters(failedQuarters,c1);
    if (histo){
        histo->Draw("sameTEXT");
    }
    histSaver->SaveCanvas(c1);
    delete c1;
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail() {

    Int_t yBins = hCellsLandau.at(0)->GetNbinsX();
    Float_t yMin = hCellsLandau.at(0)->GetBinLowEdge(1);
    Float_t yMax = hCellsLandau.at(0)->GetBinLowEdge(yBins) + hCellsLandau.at(0)->GetBinWidth(yBins);
    TString name = "hCellsLandau2DHighlightedQuarterFail"+appendix;
    name.Append(appendix);
    TH2D* hCellsLandau2DHighlightedQuarterFail = new TH2D(name,name,settings->GetNCells3d(),0,settings->GetNCells3d(),yBins,yMin,yMax);
    hCellsLandau2DHighlightedQuarterFail->GetXaxis()->SetTitle("Cell");
    hCellsLandau2DHighlightedQuarterFail->GetYaxis()->SetTitle("Charge ADC");
    name ="hCellsLandau2DHighlightedQuarterFailSorted"+appendix;
    name.Append(appendix);
    TH2D* hCellsLandau2DHighlightedQuarterFailSorted = (TH2D*)hCellsLandau2DHighlightedQuarterFail->Clone(name);
    vector<TH1F*> hCellLandausSorted = hCellsLandau;
    sort(hCellLandausSorted.begin(), hCellLandausSorted.end(), TSettings::SorterForPulseHeightOfHisto);

    name = "hHighlightedQuarterFail";
    name.Append(appendix);
    Int_t ncells = settings->GetNCells3d();
    TH2D* hHighlightedQuarterFail = new TH2D(name,name,ncells,0,ncells,1,yMin,yMax);
    name = "hHighlightedQuarterFailSorted";
    name.Append(appendix);
    TH2D* hHighlightedQuarterFailSorted = (TH2D*) hHighlightedQuarterFail->Clone(name);
    TH1D* histo = new TH1D("hMeanPulseHeightCells","hMeanPulseHeightCells",ncells,0,ncells);
    //hCellsLandau2DQuarterFail->GetXaxis()->SetTitle("Charge ADC");
    //hCellsLandau2DQuarterFail->GetYaxis()->SetTitle("Cell");
    for (UInt_t pos = 0 ; pos < ncells; pos ++){
        Int_t cellSorted = -1;
        string title = hCellLandausSorted[pos]->GetName();
        size_t found=title.find_last_of('_');
        if (found != string::npos){
            title = title.substr(found+1);
            cellSorted = atoi(title.c_str());
        }
        histo->Fill(pos,hCellsLandau.at(pos)->GetMean());

        TString binLabel = TString::Format("%3d",pos);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetBinLabel(pos,binLabel);

        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetLabelSize(0.02);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetBinLabel(pos,binLabel);

        binLabel = TString::Format("%3d",cellSorted);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetBinLabel(pos,binLabel);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetLabelSize(0.02);

        for(int yBin = 0; yBin<yBins; yBin++){
            Float_t ph = hCellsLandau[pos]->GetXaxis()->GetBinCenter(yBin+1);
            Float_t weight = hCellsLandau[pos]->GetBinContent(yBin);
            hCellsLandau2DHighlightedQuarterFail->Fill(pos,ph,weight);

            ph = hCellLandausSorted[pos]->GetXaxis()->GetBinCenter(yBin+1);
            weight =hCellLandausSorted[pos]->GetBinContent(yBin);
            hCellsLandau2DHighlightedQuarterFailSorted->Fill(pos,yBin,weight);
        }
        hHighlightedQuarterFail->SetBinContent(pos,1,HighlightGrading.at(pos));
        hHighlightedQuarterFailSorted->SetBinContent(pos,1,HighlightGrading.at(cellSorted));
    }
    histSaver->SaveHistogram(histo);
    delete histo;
    name = "cCellsLandau2DHighlightedQuarterFail";
    name.Append(appendix);
    TCanvas* cCellsLandau2DHighlightedQuarterFail = new TCanvas(name,name);
    cCellsLandau2DHighlightedQuarterFail->cd();
    //hCellsLandau2DHighlightedQuarterFail->SetEntries(hCellsLandau2DEntries);
    hCellsLandau2DHighlightedQuarterFail->SetStats(kFALSE);
    hHighlightedQuarterFail->SetStats(kFALSE);
    hCellsLandau2DHighlightedQuarterFail->Draw("COLZ");
    hHighlightedQuarterFailSorted->Draw("COLAHsame");
    hCellsLandau2DHighlightedQuarterFailSorted->Draw("sameCOLZ");
    histSaver->SaveCanvas(cCellsLandau2DHighlightedQuarterFail);
    delete cCellsLandau2DHighlightedQuarterFail;

    name = "cCellsLandau2DHighlightedQuarterFailSorted";
    name.Append(appendix);
    cCellsLandau2DHighlightedQuarterFail = new TCanvas(name,"cCellsLandau2DHighlightedQuarterFail - sorted"+appendix);
    cCellsLandau2DHighlightedQuarterFail->cd();
    //hCellsLandau2DHighlightedQuarterFail->SetEntries(hCellsLandau2DEntries);
    hCellsLandau2DHighlightedQuarterFailSorted->SetStats(kFALSE);
    hHighlightedQuarterFailSorted->SetStats(kFALSE);
    hCellsLandau2DHighlightedQuarterFailSorted->Draw("COLZ");
    hHighlightedQuarterFailSorted->Draw("COLAHsame");
    hCellsLandau2DHighlightedQuarterFailSorted->Draw("sameCOLZ");
    histSaver->SaveCanvas(cCellsLandau2DHighlightedQuarterFail);
    delete cCellsLandau2DHighlightedQuarterFail;

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsClusterSize2DVsGrading() {

    Int_t MaxClusterSize = hCellsClusteSize.at(0)->GetNbinsX() - 1; 		//Because Cluster Histo starts at 0.
    Int_t MaxGrading = settings->getNQuarters3d();
    Int_t MaxGradingBin = MaxGrading+1;

    //hCellsClusterSize2D
    stringstream hCellsClusterSize2DName; hCellsClusterSize2DName<<"hCellsClusterSize2D"<<FileNameEnd;
    TH2D* hCellsClusterSize2D = new TH2D(hCellsClusterSize2DName.str().c_str(),hCellsClusterSize2DName.str().c_str(),MaxGradingBin,0,MaxGradingBin,MaxClusterSize,0,MaxClusterSize);
    hCellsClusterSize2D->GetXaxis()->SetTitle("Grading");
    hCellsClusterSize2D->GetYaxis()->SetTitle("ClusterSize");

    for(int i=0;i<MaxClusterSize;i++){		//Set yAxis BinLabels
        stringstream jNumber;
        if(i==MaxClusterSize-1)
            jNumber<<MaxClusterSize<<"+";
        else
            jNumber<<(i+1);
        hCellsClusterSize2D->GetYaxis()->SetBinLabel(i+1,jNumber.str().c_str());
    }
    for(int i=0;i<MaxGradingBin;i++){		//Set xAxis BinLabels
        stringstream jNumber;
        jNumber<<(i);
        hCellsClusterSize2D->GetXaxis()->SetBinLabel(i+1,jNumber.str().c_str());
    }
    hCellsClusterSize2D->SetContour(99);

    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        for(UInt_t row=0;row<settings->getNRows3d();row++){

            Int_t cell = settings->get3DCellNo((int)column,row);
            Int_t Grading = CellGrading.at(cell);

            //histSaver->SaveHistogram(hCellsClusteSize.at(cell));

            Int_t xBins = hCellsClusteSize.at(0)->GetNbinsX();
            Int_t NumEvents = 0;

            for(int xBin =1; xBin<=xBins; xBin++){
                NumEvents = hCellsClusteSize.at(cell)->GetBinContent(xBin);
                Int_t ClusterSize = xBin-2;
                hCellsClusterSize2D->Fill(Grading,ClusterSize,NumEvents);
            }
        }
    }

    TCanvas* cCellsClusterSize2D = new TCanvas("cCellsClusterSize2D","cCellsClusterSize2D");
    cCellsClusterSize2D->cd();
    //hCellsClusterSize2D->SetStats(kFALSE);
    hCellsClusterSize2D->Draw("COLZ");
    /**h2DClusterSizeClone1 = *h2DClusterSize/(*h2DClusterSizeClone);
	gStyle->SetPaintTextFormat("3.2g");
	h2DClusterSizeClone1->Draw("sameTEXT");
	gStyle->SetPaintTextFormat("g");*/
    histSaver->SaveCanvas(cCellsClusterSize2D);

}

void TAnalysisOf3dDiamonds::LongAnalysis_FillEdgeFreeHistos(Float_t xPredDet,Float_t yPredDet, Float_t charge ){
    if(settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet)!=3)
        return;
    bool isInEdgeRegion = false;
    if (verbosity >8) cout<<"TAnalysisOf3dDiamonds::LongAnalysis_FillEdgeFreeHistos"<<xPredDet<<" "<<yPredDet<<" "<<charge<<endl;

    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);
    isInEdgeRegion =  settings->IsOnTheEdgeOfCell(relPos.first,relPos.second);
    if (isInEdgeRegion){
        if(hPulseHeigthEdgeRegion)
            hPulseHeigthEdgeRegion->Fill(xPredDet,yPredDet,charge);
        hEventsEdgeRegion->Fill(relPos.first,relPos.second);
    }
    else{
        if(hPulseHeigthCentralRegion)
            hPulseHeigthCentralRegion->Fill(xPredDet,yPredDet,charge);
        hEventsCentralRegion->Fill(relPos.first,relPos.second);
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayedHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize) {
    //	hCellsOverlayCharge->Fill(xRelPosDet,yRelPosDet,clusterCharge);
    //	hCellsOverlayEvents->Fill(xRelPosDet,yRelPosDet,1);

    Float_t pw = settings->getPitchWidth(subjectDetector,2);
    Float_t OffsetX = settings->getOverlayOffsetX();
    Float_t OffsetY = settings->getOverlayOffsetY();

    hCellsOverlayAvrgCharge.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
    hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);

    hCellOverlayWithColumnLandau.at(ClusterSize)->Fill(clusterCharge);

    if(!settings->isBadCell(3,cellNo)){
        hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
        hCellsLandauMinusBadCells.at(ClusterSize)->Fill(clusterCharge);
        hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);
        if(ClusterSize ==3)
            hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);
        //Int_t CellNumber = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge)
        if(ClusterSize ==3){
            LongAnalysis_Fill3DCellOverlayIndividualBinHistos(xRelPosDet, yRelPosDet, clusterCharge, ClusterSize);
            LongAnalysis_Fill3DOffsetOverlayBiasColumnAlignment(xRelPosDet, yRelPosDet, clusterCharge, ClusterSize);
        }

        /*Int_t XBin = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetXaxis()->FindBin(xRelPosDet);
    	Int_t YBin = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetYaxis()->FindBin(yRelPosDet);
    	Int_t YBins = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetYaxis()->GetNbins();
    	Int_t T2DBin = (XBin*YBins + YBin);
    	cout<<"T2DBin: "<<T2DBin<<endl;
    	cout<<"xbin + (nxbins+2)*ybin: "<<(XBin + (15+2)*YBin)<<endl;
    	cout<<"HELLO"<<endl;
    	hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetBinContent(XBin,YBin,clusterCharge);
    	cout<<"xRelPosDet: "<<xRelPosDet<<" yRelPosDet: "<<yRelPosDet<<" CellNumber: "<<CellNumber<<endl;*/
    }
    if(settings->IsGoodCell(3, cellNo)){
        hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
        hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);
    }
    if(settings->isBadCell(3,cellNo))
        hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);

    if(settings->IsWithInTheColumnRadius(xRelPosDet, yRelPosDet)){
        hCellOverlayColumnLandau.at(ClusterSize)->Fill(clusterCharge);
    }
    else{
        hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
        hCellOverlayNoColumnLandau.at(ClusterSize)->Fill(clusterCharge);
    }

}

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayCentralColumnHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize, TCluster *diamondCluster) {

    Float_t CentralColumnXLow = settings->getCentralColumnOverlayXLow();
    Float_t CentralColumnXHigh = settings->getCentralColumnOverlayXHigh();
    Float_t CentralColumnYLow = settings->getCentralColumnOverlayYLow();
    Float_t CentralColumnYHigh = settings->getCentralColumnOverlayYHigh();

    if(xRelPosDet>CentralColumnXLow && xRelPosDet<CentralColumnXHigh &&
            yRelPosDet>CentralColumnYLow && yRelPosDet<CentralColumnYHigh){

        hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
        hCellsCentralColumnOverlayLandau.at(ClusterSize)->Fill(clusterCharge);

        if(!settings->isBadCell(3,cellNo)){
            hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
            hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize)->Fill(clusterCharge);
            if(ClusterSize == 3){
                //To THStack column Histos, shift histo low edge back to (0,0).
                Float_t XShift = hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetXaxis()->GetBinLowEdge(1);
                Float_t YShift = hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetYaxis()->GetBinLowEdge(1);
                if(verbosity>5)cout<<"XShift: "<<XShift<<" YShift: "<<YShift<<endl;
                Float_t xRelPosDetOffsetShifted = xRelPosDet - XShift;
                Float_t yRelPosDetOffsetShifted = yRelPosDet - YShift;
                hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->Fill(xRelPosDetOffsetShifted,yRelPosDetOffsetShifted,clusterCharge);
            }
            //cout<<"xRelPosDet: "<<xRelPosDet<<" yRelPosDet: "<<yRelPosDet<<" clusterCharge: "<<clusterCharge<<endl;
            //cout<<"settings->getOverlayColumnPulseHeightCut()"<<settings->getOverlayColumnPulseHeightCut()<<endl;
            if(clusterCharge<settings->getOverlayColumnPulseHeightCut()){
                hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,1);
                hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
            }

        }
        if(settings->IsGoodCell(3, cellNo)){
            hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
            hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->Fill(clusterCharge);
        }

        //LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(cellNo,xRelPosDet,yRelPosDet, ClusterSize, diamondCluster);
    }

}

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayBiasColumnHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize, TCluster *diamondCluster) {

    Float_t BiasColumnXLow = settings->getBiasColumnOverlayXLow();
    Float_t BiasColumnXHigh = settings->getBiasColumnOverlayXHigh();
    Float_t BiasColumnYLow = settings->getBiasColumnOverlayYLow();
    Float_t BiasColumnYHigh = settings->getBiasColumnOverlayYHigh();

    Float_t pw = settings->getPitchWidth(subjectDetector,2);
    Float_t OffsetX = settings->getOverlayOffsetX();
    Float_t OffsetY = settings->getOverlayOffsetY();

    if(ClusterSize ==3){
        if(fmod(xRelPosDet+OffsetX,pw)>BiasColumnXLow && fmod(xRelPosDet+OffsetX,pw)<BiasColumnXHigh &&
                fmod(yRelPosDet+OffsetY,pw)>BiasColumnYLow && fmod(yRelPosDet+OffsetY,pw)<BiasColumnYHigh){

            //hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
            //hCellsCentralColumnOverlayLandau.at(ClusterSize)->Fill(clusterCharge);

            if(!settings->isBadCell(3,cellNo)){
                //printf("BiasColumnXLow: %f, BiasColumnXHigh %f, BiasColumnYLow %f, BiasColumnYHigh %f \n",
                //BiasColumnXLow, BiasColumnXHigh, BiasColumnYLow, BiasColumnYHigh);
                //printf("xRelPosDet: %f, yRelPosDet %f, OffsetxRelPosDet %f, OffsetyRelPosDet %f \n",
                //xRelPosDet, yRelPosDet, xRelPosDetOffset, yRelPosDetOffset);
                hCellsBiasColumnOverlayAvrgChargeMinusBadCells->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);
                //To THStack column Histos, shift histo low edge back to (0,0).
                Float_t XShift = hCellsBiasColumnOverlayAvrgChargeMinusBadCells->GetXaxis()->GetBinLowEdge(1);
                Float_t YShift = hCellsBiasColumnOverlayAvrgChargeMinusBadCells->GetYaxis()->GetBinLowEdge(1);
                if(verbosity>5)cout<<"XShift: "<<XShift<<" YShift: "<<YShift<<endl;
                Float_t xRelPosDetOffsetShifted = fmod(xRelPosDet+OffsetX,pw) - XShift;
                Float_t yRelPosDetOffsetShifted = fmod(yRelPosDet+OffsetY,pw) - YShift;
                hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->Fill(xRelPosDetOffsetShifted,yRelPosDetOffsetShifted,clusterCharge);
                hCellsBiasColumnOverlayLandauMinusBadCells->Fill(clusterCharge);
                //LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(cellNo,xRelPosDet,yRelPosDet,clusterCharge, ClusterSize, diamondCluster);
                //cout<<"xRelPosDet: "<<xRelPosDet<<" yRelPosDet: "<<yRelPosDet<<" clusterCharge: "<<clusterCharge<<endl;
                //cout<<"settings->getOverlayColumnPulseHeightCut()"<<settings->getOverlayColumnPulseHeightCut()<<endl;

                if(clusterCharge<settings->getOverlayColumnPulseHeightCut()){
                    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),1);
                    hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);
                }

            }
            /*if(settings->IsGoodCell(3, cellNo)){
			hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
			hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize)->Fill(clusterCharge);
		}*/

        }
    }

}


void TAnalysisOf3dDiamonds::LongAnalysis_Fill3DOffsetOverlayBiasColumnAlignment(Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize) {

    Float_t pw = settings->getPitchWidth(subjectDetector,2);

    for(int i=0; i<ShiftX.size(); i++){
        Float_t OffsetX = settings->getOverlayOffsetX() + ShiftX[i];
        for(int j=0; j<ShiftY.size(); j++){
            Float_t OffsetY = settings->getOverlayOffsetY() + ShiftY[j];
            Int_t Shift = i*ShiftY.size() + j;

            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->Fill(fmod(xRelPosDet+OffsetX,pw),fmod(yRelPosDet+OffsetY,pw),clusterCharge);

            Int_t XBin = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetXaxis()->FindBin(fmod(xRelPosDet+OffsetX,pw));
            Int_t YBin = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetYaxis()->FindBin(fmod(yRelPosDet+OffsetY,pw));
            Int_t YBins = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Shift)->GetYaxis()->GetNbins();
            Int_t BinNum = (XBin-1)*YBins + YBin-1;

            hOverlayCellOffsetAlignmentBinHits.at(Shift)->Fill(BinNum);
            if(clusterCharge<settings->getOverlayColumnPulseHeightCut())
                hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Shift)->Fill(BinNum);

        }
    }

}

void TAnalysisOf3dDiamonds::LongAnalysis_Fill3DCellOverlayIndividualBinHistos(Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize) {

    //hCellsOverlayAvrgChargeMinusBadCells
    if ( ClusterSize >= hCellsOverlayAvrgChargeMinusBadCells.size() ||
            hCellsOverlayAvrgChargeMinusBadCells.size() == 0){
        cerr<<"[ TAnalysisOf3dDiamonds::LongAnalysis_Fill3DCellOverlayIndividualBinHistos] "<<
                xRelPosDet<<"/"<<yRelPosDet<<" with charge: "<<clusterCharge<<" has Invalid ClusterSize "<<ClusterSize<<" but size is"<<hCellsOverlayAvrgChargeMinusBadCells.size()<<endl;
        return;
    }
    if (verbosity>6) cout<<"[ TAnalysisOf3dDiamonds::LongAnalysis_Fill3DCellOverlayIndividualBinHistos]  find "<<ClusterSize<<" in "<<hCellsOverlayAvrgChargeMinusBadCells.size()<<" sized hCellsOverlayAvrgChargeMinusBadCells"<<endl;
    Int_t XBin = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetXaxis()->FindBin(xRelPosDet);
    Int_t YBin = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetYaxis()->FindBin(yRelPosDet);
    Int_t YBins = hCellsOverlayAvrgChargeMinusBadCells.at(0)->GetYaxis()->GetNbins();
    Int_t BinNum = (XBin-1)*YBins + YBin-1;
    if (BinNum >= VecOverlayCellBinHistos.size()){
        cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_Fill3DCellOverlayIndividualBinHistos] invalid BinNum: "<<BinNum
                << "( Xbin: "<<XBin<<"/Ybin: "<<YBin<<" with BinNum = (XBin-1)*YBins + YBin-1) for a vector of size: "<<VecOverlayCellBinHistos.size()<<endl;
        cout<<"\t xRelPosDet: "<<xRelPosDet<<"\t yRelPosDet:"<<yRelPosDet<<endl;
        return;
    }

    VecOverlayCellBinHistos.at(BinNum)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
    VecOverlayCellBinLandaus.at(BinNum)->Fill(clusterCharge);
    hOverlayCellBinHits->Fill(BinNum);
    if(clusterCharge<settings->getOverlayColumnPulseHeightCut())
        hOverlayCellBinHitsBelowCut->Fill(BinNum);

    Float_t pw = settings->getPitchWidth(subjectDetector,2);
    Float_t OffsetX = settings->getOverlayOffsetX();
    Float_t OffsetY = settings->getOverlayOffsetY();

    //hCellsOffsetOverlayAvrgChargeMinusBadCells
    XBin = hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetXaxis()->FindBin(fmod(xRelPosDet+OffsetX,pw));
    YBin = hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetYaxis()->FindBin(fmod(yRelPosDet+OffsetY,pw));
    YBins = hCellsOffsetOverlayAvrgChargeMinusBadCells.at(0)->GetYaxis()->GetNbins();
    BinNum = (XBin-1)*YBins + YBin-1;

    hOverlayCellOffsetBinHits->Fill(BinNum);
    if(clusterCharge<settings->getOverlayColumnPulseHeightCut())
        hOverlayCellOffsetBinHitsBelowCut->Fill(BinNum);

    //hCellsOverlayUnEvenBinningAvrgChargeMinusBadCells
    XBin = hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetXaxis()->FindBin(fmod(xRelPosDet+OffsetX,pw));
    YBin = hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetYaxis()->FindBin(fmod(yRelPosDet+OffsetY,pw));
    YBins = hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetYaxis()->GetNbins();
    BinNum = (XBin-1)*YBins + YBin-1;

    hOverlayCellUnEvenBinningBinHits->Fill(BinNum);
    if(clusterCharge<settings->getOverlayColumnPulseHeightCut())
        hOverlayCellUnEvenBinningBinHitsBelowCut->Fill(BinNum);

}
void TAnalysisOf3dDiamonds::MakeGhostCluster(TCluster *diamondCluster, Int_t ClusterSize){
    Int_t ClusterStart = diamondCluster->getFirstHitChannel();
    Int_t ClusterEnd = diamondCluster->getLastHitChannel();
    if(verbosity>5)cout<<"\nClusterStart: "<<ClusterStart<<" ClusterEnd: "<<ClusterEnd<<endl;
    pair<int,int> channels = settings->diamondPattern.getPatternChannels(3);
    Float_t ghostHit;
    do{
        ghostHit = gRandom->Uniform(channels.first,channels.second);
    }
    while(ClusterStart<=ghostHit&&ghostHit<=ClusterEnd);
    ghostHit-=.5;
    if(verbosity>5)cout<<"new Ghost Hit at "<<ghostHit<<endl;
    GhostCluster = TTransparentAnalysis::makeTransparentCluster(eventReader,settings,subjectDetector,ghostHit,ClusterSize);
    if(verbosity>5) GhostCluster.Print(1);
}

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet, Float_t ClusterSize, TCluster *diamondCluster) {

//    diamondCluster->SetTransparentClusterSize(5);
    //    MakeGhostCluster(diamondCluster,5);
    Int_t cs = diamondCluster->GetTransparentClusterSize();
    diamondCluster->SetTransparentClusterSize(ClusterSize);
    GhostCluster.SetTransparentClusterSize(ClusterSize);
    Float_t GhostClusterCharge = GhostCluster.getPositiveCharge(useCMN);
    hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,GhostClusterCharge);
    hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->Fill(GhostClusterCharge);
    diamondCluster->SetTransparentClusterSize(cs);
}


void TAnalysisOf3dDiamonds::LongAnalysis_SaveQuarterCellsClusterSize2DVsGrading() {

    Int_t MaxClusterSize = hQuarterCellsClusterSize[0].at(0)->GetNbinsX() - 1; 		//Because Cluster Histo starts at 0.
    Int_t MaxGrading = settings->getNQuarters3d();
    Int_t MaxGradingBin = MaxGrading+1;

    stringstream hQuarterCellsClusterSize2DName; hQuarterCellsClusterSize2DName<<"hQuarterCellsClusterSize2D"<<FileNameEnd;
    TH2D* hQuarterCellsClusterSize2D = new TH2D(hQuarterCellsClusterSize2DName.str().c_str(),hQuarterCellsClusterSize2DName.str().c_str(),2*MaxGradingBin,0,2*MaxGradingBin,MaxClusterSize,0,MaxClusterSize);
    hQuarterCellsClusterSize2D->GetXaxis()->SetTitle("");
    hQuarterCellsClusterSize2D->GetYaxis()->SetTitle("ClusterSize");

    for(int i=0;i<MaxClusterSize;i++){		//Set yAxis BinLabels
        stringstream jNumber;
        if(i==MaxClusterSize-1)
            jNumber<<MaxClusterSize<<"+";
        else
            jNumber<<(i+1);
        hQuarterCellsClusterSize2D->GetYaxis()->SetBinLabel(i+1,jNumber.str().c_str());
    }

    for(int i=0;i<MaxGradingBin;i++){		//Set xAxis P/F BinLabels
        //stringstream iLetter; iLetter<<PassFail.at(0);
        stringstream Pass; Pass<<"P";
        hQuarterCellsClusterSize2D->GetXaxis()->SetBinLabel(2*i+1,Pass.str().c_str());
        stringstream Fail; Fail<<"F";
        hQuarterCellsClusterSize2D->GetXaxis()->SetBinLabel(2*i+2,Fail.str().c_str());
        cout<<"i: "<<i<<endl;
    }

    hQuarterCellsClusterSize2D->SetContour(99);
    //h2DClusterSize->SetTitleOffset(0.015,"X");


    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        for(UInt_t row=0;row<settings->getNRows3d();row++){
            Int_t cell = settings->get3DCellNo((int)column,row);
            for(int quarter=0;quarter<settings->getNQuarters3d();quarter++){
                Int_t Grading = 2*CellGrading.at(cell) + vecQuarterCellsPassFail[cell][quarter];;
                //histSaver->SaveHistogram(hQuarterCellsClusterSize[cell][quarter]);

                Int_t xBins = hQuarterCellsClusterSize[cell][0]->GetNbinsX();
                Int_t NumEvents = 0;

                for(int xBin =1; xBin<=xBins; xBin++){
                    NumEvents = hQuarterCellsClusterSize[cell][quarter]->GetBinContent(xBin);
                    Int_t ClusterSize = xBin-2;
                    cout<<"cell: "<<cell<<" Quarter: "<<quarter<<" Grading: "<<Grading<<" clusterSize: "<<ClusterSize<<endl;
                    NumEvents = hQuarterCellsClusterSize[cell][quarter]->GetBinContent(xBin);
                    hQuarterCellsClusterSize2D->Fill(Grading,ClusterSize,NumEvents);
                }
            }
        }
    }

    TCanvas* cQuarterCellsClusterSize2D = new TCanvas("cQuarterCellsClusterSize2D","cQuarterCellsClusterSize2D");
    cQuarterCellsClusterSize2D->cd();
    //hQuarterCellsClusterSize2D->SetStats(kFALSE);
    hQuarterCellsClusterSize2D->Draw("COLZ");
    cout<<"hQuarterCellsClusterSize2D"<<endl;
    /**h2DClusterSizeQuarterCellClone1 = *h2DClusterSizeQuarterCell/(*h2DClusterSizeQuarterCellClone);
			gStyle->SetPaintTextFormat("3.2g");
			h2DClusterSizeQuarterCellClone1->Draw("sameTEXT");
			gStyle->SetPaintTextFormat("g");*/
    //h2DClusterSizeXAxis->Draw("sameCOL");
    for(int i=0;i<MaxGradingBin+1;i++){			//Draw lines for Bin edges.
        TLine* BinEdge = new TLine(i*2,0,i*2,-.5);
        BinEdge->SetLineWidth(0.5);
        BinEdge->SetLineColor(kBlack);
        BinEdge->Draw("same");
        if(i<MaxGradingBin){					//set xAxis Grading bin labels
            stringstream Label; Label<<i;
            TText* Text = new TText(hQuarterCellsClusterSize2D->GetXaxis()->GetBinCenter(2*i+1)+.4,-.5,Label.str().c_str());
            Text->SetTextSize(0.04);
            Text->Draw("same");
        }
    }
    TText* XTitle = new TText((hQuarterCellsClusterSize2D->GetXaxis()->GetBinCenter(2*MaxGradingBin)-.75),-.80,"Grading");
    XTitle->SetTextSize(0.04);					//set xAxis title
    XTitle->Draw("same");
    cout<<"save cQuarterCellsClusterSize2D"<<endl;
    histSaver->SaveCanvas(cQuarterCellsClusterSize2D);

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveGoodAndBadCellLandaus() {
    cout<<"[TAnalysisOf3dDiamons::LongAnalysis_SaveGoodAndBadCellLandaus]"<<endl;
    TList *listGoodCells = new TList;
    TList *listBadCells = new TList;
    string plots_path = histSaver->GetPlotsPath();
    HistogrammSaver newHistSaver(settings);
    newHistSaver.SetPlotsPath(plots_path+(string)"/CellLandaus/");
    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        TString name = "hAllCellColumnLandau_Column_";
        name.Append(settings->getColumnChar(column));
        name+=appendix;
        TH1F* hColumnLandau = (TH1F*)hCellsLandau.at(0)->Clone(name);
        hColumnLandau->Reset();
        hColumnLandau->SetTitle(name);

        name = "hAllButBadCellColumnLandau_Column_";
        name.Append(settings->getColumnChar(column));
        name+=appendix;
        TH1F* hColumnLandauNotBad = (TH1F*)hCellsLandau.at(0)->Clone(name);
        hColumnLandauNotBad->Reset();
        hColumnLandauNotBad->SetTitle(name);

        name = "hGoodCellColumnLandau_Column_";
        name.Append(settings->getColumnChar(column));
        name+=appendix;
        TH1F* hColumnLandauGood = (TH1F*)hCellsLandau.at(0)->Clone(name);
        hColumnLandauGood->Reset();
        hColumnLandauGood->SetTitle(name);

        for(UInt_t row=0;row<settings->getNRows3d();row++){
            Int_t cell = settings->get3DCellNo((int)column,row);
            TH1F* h = hCellsLandau.at(cell);
            //hCellNumbering->SetBinContent(column+1,row+1,cell); //This should be a clone of the 2D Cell Mean Charge Plot, Wait till Felix has finished.
            newHistSaver.SaveHistogram(h);
            hColumnLandau->Add(h);
            if (!settings->isBadCell(3,cell))
                hColumnLandauNotBad->Add(h);
            if (settings->IsGoodCell(3,cell))
                hColumnLandau->Add(h);
            for(UInt_t i=0; i<settings->getBadCells3D().size(); i++)
                if(cell==settings->getBadCells3D().at(i)){
                    //                    Int_t Entries = hLandauBadCells->GetEntries();
                    //                    hLandauBadCells->Add(hCellsLandau.at(cell),1);	//Not working for some reason, ask Felix
                    //                    hLandauBadCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                    listBadCells->Add(h);
                    //                    cout<<"\tbad"<<endl;
                }
            for(UInt_t i=0; i<settings->getGoodCellRegions3d().size(); i++){
                for(UInt_t j=0; j<settings->getGoodCellRegions3d().at(i).size(); j++)
                    if(cell==settings->getGoodCellRegions3d().at(i).at(j)){
                        //                        Int_t Entries = hLandauGoodCells->GetEntries();
                        //                        hLandauGoodCells->Add(hCellsLandau.at(cell));	//Not working for some reason, ask Felix
                        //                        hLandauGoodCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                        listGoodCells->Add(h);
                        //                        cout<<"\tgood"<<endl;
                    }
            }//end good cells region
        }//end row
        newHistSaver.SaveHistogram(hColumnLandau);
        newHistSaver.SaveHistogram(hColumnLandauGood);
        newHistSaver.SaveHistogram(hColumnLandauNotBad);
        delete hColumnLandau;
    }//end columns
    cout<<"List Good Cells: "<<listGoodCells->GetEntries()<<endl;
    listGoodCells->Print();
    cout<<"\nList Bad Cells: "<<listBadCells->GetEntries()<<endl;
    listBadCells->Print();

    TString name = "hLandauBadCells";
    name.Append(appendix);
    TH1F* hLandauBadCells = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauBadCells->SetTitle("Landau of bad cells");
    hLandauBadCells->GetXaxis()->SetTitle("pulse height /adc");
    hLandauBadCells->GetXaxis()->SetTitle("number of entries #");
    hLandauBadCells->Reset();
    hLandauBadCells->Merge(listBadCells);
    histSaver->SaveHistogram(hLandauBadCells);

    name = "hLandauGoodCells";
    name.Append(appendix);
    TH1F* hLandauGoodCells = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauGoodCells->SetTitle("Landau of good cells");
    hLandauGoodCells->GetXaxis()->SetTitle("pulse height /adc");
    hLandauGoodCells->GetXaxis()->SetTitle("number of entries #");
    hLandauGoodCells->Reset();
    hLandauGoodCells->Merge(listGoodCells);
    //    histSaver->SaveHistogram(hLandauGoodCells);


    cout<<"hLandauGoodCells: "<<hLandauGoodCells->GetEntries()<<endl;
    cout<<"hLandauStrip:     "<<hLandauStrip->GetEntries()<<endl;

    Float_t factor = hLandauGoodCells->GetBinContent(hLandauGoodCells->GetMaximumBin());
    factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
    name = "cLandauGoodCells";
    name.Append(appendix);
    histSaver->SaveTwoHistos(name,hLandauGoodCells,hLandauStrip,factor,"right");

    name = "cLandauGoodCellsNormalized";
    name.Append(appendix);
    histSaver->SaveTwoHistosNormalized(name,hLandauGoodCells,hLandauStrip,1,"right");

    factor = hLandauBadCells->GetBinContent(hLandauBadCells->GetMaximumBin());
    factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
    name = "cLandauBadCells";
    name.Append(appendix);
    histSaver->SaveTwoHistos(name,hLandauBadCells,hLandauStrip,factor,"right");
    name = "cLandauBadCellsNormalized";
    name.Append(appendix);
    histSaver->SaveTwoHistosNormalized(name,hLandauBadCells,hLandauStrip,1,"right");
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveDeadCellProfile() {
    //hDeadCell
    //histSaver->SaveHistogram(hDeadCell);
    //histSaver->SaveHistogram(hDeadCellEvents);
    /*for(int i=0; i<settings->getDeadCell3D().size(); i++){
		cDeadCellMeanCharge.push_back(new TCanvas("cDeadCellMeanCharge","cDeadCellMeanCharge"));
		cDeadCellMeanCharge.at(i)->cd();
     *hDeadCellMeanCharge.at(i) = (*hDeadCell.at(i)/(*hDeadCellEvents.at(i)));
		hDeadCellMeanCharge.at(i)->SetEntries(hDeadCellEvents.at(i)->Integral());
		hDeadCellMeanCharge.at(i)->Draw();
		Float_t ymax1 = hDeadCellMeanCharge.at(i)->GetMaximum();
		TLine* CellEdge1 = new TLine(150,0,150,ymax1);
		TLine* CellEdge2 = new TLine(300,0,300,ymax1);
		CellEdge1->SetLineWidth(2);		CellEdge2->SetLineWidth(2);
		CellEdge1->SetLineColor(kRed);	CellEdge2->SetLineColor(kRed);
		CellEdge1->Draw("same");		CellEdge2->Draw("same");
		histSaver->SaveCanvas(cDeadCellMeanCharge.at(i));
	}*/

    //	vector<TCanvas*> cDeadCellMeanCharge;
    Float_t pw = settings->getPitchWidth(subjectDetector,2);
    Float_t maxSigma = pw/2;
    TF1* fitX = new TF1("fit",
            "[0]*TMath::Sqrt(TMath::Pi()/2)*[1]*(TMath::Erf(([5]+[2]+[3]-x)/TMath::Sqrt(2)/[1])+TMath::Erf(([3]-[5]-[2]+x)/TMath::Sqrt(2)/[1]))+[4]",
            0,300);
    fitX->SetParLimits(0,-500,0);		// Integral
    fitX->SetParameter(0,-100);			// Integral
    fitX->SetParLimits(1,0,maxSigma);	// Sigma Gaus
    fitX->SetParameter(1,2);			// Sigma Gaus
    fitX->SetParLimits(2,-pw,pw);		// rel. Pos wrt reference
    fitX->SetParameter(2,0);			// rel. Pos wrt reference
    fitX->FixParameter(3,pw/2);			// Pitch
    fitX->SetParLimits(4,100,3000);		// Offset
    fitX->SetParameter(4,1000);			// Offset
    fitX->FixParameter(5,1.5*pw);		// Reference position
    fitX->SetParNames("Integral","#Sigma_{Gaus}","Rel. position wrt reference","Pitch","Offset","Reference input");

    TCanvas *c1;
    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){

        if(!hDeadCellCharge[i]){
            cerr<<TString::Format("hDeadCellCharge[%d] invalid ", i)<<endl;
            continue;
        }
        TF1* fit = (TF1*)fitX->Clone(TString::Format("fDeadCell_%d",i));
        TString name = TString::Format("cDeadCellMeanCharge_%02d",i)+appendix;
        c1 = new TCanvas(name,name);
        c1->cd();
        hDeadCellCharge[i]->Fit(fit);
        hDeadCellCharge[i]->Draw("");
        TCutG* cellEdges = new TCutG(name,4);
        cellEdges->SetPoint(0,150,-1e9);
        cellEdges->SetPoint(1,150,1e9);
        cellEdges->SetPoint(2,300,1e9);
        cellEdges->SetPoint(3,300,-1e9);
        cellEdges->SetLineWidth(2);
        cellEdges->SetLineColor(kRed);
        cellEdges->Draw("same");
        histSaver->SaveCanvas(c1);
        delete c1;
        histSaver->SaveHistogramWithCellGrid(hDeadCellPositions[i]);
        delete hDeadCellCharge[i];
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayMeanCharge() {
    Float_t zmin = 100;
    Float_t zmax = 1200;
    //	return;
    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2 += appendix;

        cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayMeanCharge]"<<endl;
        cout<<hCellsOverlayAvrgCharge.at(ClusterSize)<<endl;
        if(hCellsOverlayAvrgCharge.at(ClusterSize)){
            cout<<hCellsOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOverlayAvrgCharge_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgCharge.at(ClusterSize)->Clone(name);
            if(verbosity)cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            delete histo;
            ////		hCellsOverlayAvrgCharge->SetName("hCellsOverlayAvrgCharge");
            //		cout<<"Set Name: "<<hCellsOverlayAvrgCharge<<endl;
            ////		hCellsOverlayAvrgCharge->SetTitle("Avrg PH - overlayed");
            //		histSaver->SaveHistogram(hCellsOverlayAvrgCharge);
            //		delete hCellsOverlayAvrgCharge ;
        }
        if(hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)){
            cout<<hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOverlayAvrgChargeGoodCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            delete histo;
        }
        if(hCellsOverlayAvrgChargeBadCells.at(ClusterSize)){
            cout<<hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOverlayAvrgChargeBadCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            delete histo;
        }
        if(hCellsLandauMinusBadCells.at(ClusterSize)){
            histSaver->SaveHistogram(hCellsLandauMinusBadCells[ClusterSize]);
            if(ClusterSize+1!=3)delete hCellsLandauMinusBadCells[ClusterSize];

        }
        if(hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)){
            cout<<hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOverlayAvrgChargeMinusBadCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            if(ClusterSize+1!=3)delete histo;
        }
        if(ClusterSize+1==3){
            TProfile2D* prof = hCellsOverlayAvrgChargeMinusBadCells[ClusterSize];
            DoMonteCarloOfAvrgChargePerBinInOverlay(prof,hCellsLandauMinusBadCells[ClusterSize]);
            delete hCellsLandauMinusBadCells[ClusterSize];
        }
        if(hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)){
            TString name = "hCellsOverlayAvrgChargeNoColumnHit_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->Clone(name);
            histo->SetTitle("Avrg PH - overlayed - no hit in columns");
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            delete histo;
            //		hCellsOverlayAvrgChargeNoColumnHit->SetName("hCellsOverlayAvrgChargeNoColumns");
            //		hCellsOverlayAvrgChargeNoColumnHit->SetTitle("Avrg PH - overlayed - no hit in columns");
            //		histSaver->SaveHistogram(hCellsOverlayAvrgChargeNoColumnHit);
            //		delete hCellsOverlayAvrgChargeNoColumnHit;
        }
        //		hCellsOverlayPulseHeight->Project3D("xy");

        if(hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)){
            cout<<hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOverlayAvrgChargeMinusBadCells_cl";
            name.Append(appendix);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(zmin,zmax);
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            delete histo;
            TH2D* hCellsOverlayAvrgChargeMinusBadCellsNentries = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->ProjectionXY("hCellsOverlayAvrgChargeMinusBadCellsNentries","b");
            hCellsOverlayAvrgChargeMinusBadCellsNentries->Draw("colz");
            hCellsOverlayAvrgChargeMinusBadCellsNentries->GetZaxis()->SetTitle("number of entries");
            histSaver->SaveHistogram(hCellsOverlayAvrgChargeMinusBadCellsNentries,false,false);
            delete hCellsOverlayAvrgChargeMinusBadCellsNentries;
        }
        if(hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)){
            TString name = "hCellsOverlayAvrgChargeNoColumnHit_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->Clone(name);
            histo->SetTitle("Avrg PH - overlayed - no hit in columns");
            cout<<"SAVE"<<endl;
            histSaver->SaveHistogram(histo);
            histSaver->SaveOverlay(histo);
            delete histo;
            TH2D* hCellsOverlayAvrgChargeNoColumnHitNentries = hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->ProjectionXY("hCellsOverlayAvrgChargeNoColumnHitNentries","b");
            hCellsOverlayAvrgChargeNoColumnHitNentries->Draw("colz");
            hCellsOverlayAvrgChargeNoColumnHitNentries->GetZaxis()->SetTitle("number of entries");
            histSaver->SaveHistogram(hCellsOverlayAvrgChargeNoColumnHitNentries,false,false);
            delete hCellsOverlayAvrgChargeNoColumnHitNentries;
            //		hCellsOverlayAvrgChargeNoColumnHit->SetName("hCellsOverlayAvrgChargeNoColumns");
            //		hCellsOverlayAvrgChargeNoColumnHit->SetTitle("Avrg PH - overlayed - no hit in columns");
            //		histSaver->SaveHistogram(hCellsOverlayAvrgChargeNoColumnHit);
            //		delete hCellsOverlayAvrgChargeNoColumnHit;
        }
        //		hCellsOverlayPulseHeight->Project3D("xy");

    } //End of for ClusterSize

    /*hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge() {
    //	return;

    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2+=appendix;

        cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge]"<<endl;
        cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)<<endl;
        if(hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgCharge_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(0,1200);
            histSaver->SaveHistogram(histo);
            delete histo;
            ////		hCellsOverlayAvrgCharge->SetName("hCellsOverlayAvrgCharge");
            //		cout<<"Set Name: "<<hCellsOverlayAvrgCharge<<endl;
            ////		hCellsOverlayAvrgCharge->SetTitle("Avrg PH - overlayed");
            //		histSaver->SaveHistogram(hCellsOverlayAvrgCharge);
            //		delete hCellsOverlayAvrgCharge ;

            /*	TProfile* pfx = (TProfile*)hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->ProfileX("Hello");
			pfx->GetYaxis()->SetTitle("Average PulseHeight [ADC]");
			pfx->GetXaxis()->SetTitle("Cell Pos X [um]");
			//histSaver->SaveHistogram(pfx);
			TCanvas* cCellsCentralColumnOverlayAvrgCharge_pfx = new TCanvas("cCellsCentralColumnOverlayAvrgCharge_pfx","cCellsCentralColumnOverlayAvrgCharge_pfx");
			cCellsCentralColumnOverlayAvrgCharge_pfx->cd();
			pfx->Draw();
			name = "cCellsCentralColumnOverlayAvrgCharge_pfx";
			name.Append(appendix);
			histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgCharge_pfx);
			delete pfx;
			TProfile*	TH2::ProfileX(const char* name = "_pfx", Int_t firstybin = 1, Int_t lastybin = -1, Option_t* option = "") constMENU
					TH1F* hEdgeFittingAvrgCharge_pfx = (TH1F*)hEdgeFittingAvrgCharge->ProfileX();*/

        }
        if(hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgChargeGoodCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeGoodCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(0,1200);
            histSaver->SaveHistogram(histo);
            delete histo;

            /*TProfile* pfx = (TProfile*)hCellsCentralColumnOverlayAvrgChargeGoodCells->ProfileX();
			pfx->GetYaxis()->SetTitle("Average PulseHeight [ADC]");
			pfx->GetXaxis()->SetTitle("Cell Pos X [um]");
			//histSaver->SaveHistogram(pfx);
			TCanvas* cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx = new TCanvas("cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx","cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx");
			cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx->cd();
			pfx->Draw();
			name = "cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx";
			name.Append(appendix);
			histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgChargeGoodCells_pfx);
			delete pfx;*/
        }
        if(hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolzTEXT");
            histo->GetZaxis()->SetRangeUser(0,1200);
            histSaver->SaveHistogram(histo);
            delete histo;

            /*TProfile* pfx = (TProfile*)hCellsCentralColumnOverlayAvrgChargeMinusBadCells->ProfileX();
			pfx->GetYaxis()->SetTitle("Average PulseHeight [ADC]");
			pfx->GetXaxis()->SetTitle("Cell Pos X [um]");
			//histSaver->SaveHistogram(pfx);
			TCanvas* cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx = new TCanvas("cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx","cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx");
			cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx->cd();
			pfx->Draw();
			name = "cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx";
			name.Append(appendix);
			histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgChargeMinusBadCells_pfx);
			delete pfx;*/
        }

        if(hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(0,1200);
            histSaver->SaveHistogram(histo);
            delete histo;

        }
        if (hLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)){
            cout<<hLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->GetName()<<endl;
            histSaver->SaveHistogram(hLandauMinusBadCellsOffsetAnalysis.at(ClusterSize));
            delete hLandauMinusBadCellsOffsetAnalysis.at(ClusterSize);

        }

        if(hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut_cl";
            name.Append(appendix2);
            TH2F* histo = (TH2F*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(0,10);
            //histSaver->SaveHistogram(histo);

            TCanvas* cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut = new TCanvas("cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut","cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut");
            cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut->cd();
            histo->Draw("colztext");

            Float_t xLow = settings->getCentralColumnOverlayXLow();
            Float_t xHigh = settings->getCentralColumnOverlayXHigh();
            Float_t yLow = settings->getCentralColumnOverlayYLow();
            Float_t yHigh = settings->getCentralColumnOverlayYHigh();

            TLine* BinLowerEdgeX = new TLine(xLow,yLow-2,xLow,yHigh+2);
            TLine* BinUpperEdgeX = new TLine(xHigh,yLow-2,xHigh,yHigh+2);
            TLine* BinLowerEdgeY = new TLine(xLow-2,yLow,xHigh+2,yLow);
            TLine* BinUpperEdgeY = new TLine(xLow-2,yHigh,xHigh+2,yHigh);

            BinLowerEdgeX->SetLineWidth(2);
            BinLowerEdgeX->SetLineColor(kRed);
            BinLowerEdgeX->Draw("same");
            BinUpperEdgeX->SetLineWidth(2);
            BinUpperEdgeX->SetLineColor(kRed);
            BinUpperEdgeX->Draw("same");
            BinLowerEdgeY->SetLineWidth(2);
            BinLowerEdgeY->SetLineColor(kRed);
            BinLowerEdgeY->Draw("same");
            BinUpperEdgeY->SetLineWidth(2);
            BinUpperEdgeY->SetLineColor(kRed);
            BinUpperEdgeY->Draw("same");

            histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut);

            delete histo;
        }

        if(hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)){
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->Clone(name);
            //TH2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->ProjectionXY("B");
            //histo->SetName(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            /*histo->GetXaxis()->SetRangeUser(settings->getCentralColumnOverlayXLow(),settings->getCentralColumnOverlayXHigh());
			histo->GetYaxis()->SetRangeUser(settings->getCentralColumnOverlayYLow(),settings->getCentralColumnOverlayYHigh());*/
            histSaver->SaveHistogram(histo);

            name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutScaled_cl";
            name.Append(appendix2);
            histo->SetName(name);
            histo->Scale(1/settings->getOverlayColumnPulseHeightCut());
            histo->GetZaxis()->SetRangeUser(0,1);
            histSaver->SaveHistogram(histo);

            /*TCanvas* cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut = new TCanvas("cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut","cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut");
			cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut->cd();
			histo->Draw("colztext");
			histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut);*/

            delete histo;
        }


        histSaver->SaveHistogram(hCellsCentralColumnOverlayLandau.at(ClusterSize));
        histSaver->SaveHistogram(hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize));
        histSaver->SaveHistogram(hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize));
        histSaver->SaveHistogram(hCellsCentralColumnOverlayLandauGoodCells.at(ClusterSize));

    } //End of for ClusterSize

    //		hCellsOverlayPulseHeight->Project3D("xy");

    /* hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsBiasColumnOverlayMeanCharge() {
    //	return;
    cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsBiasColumnOverlayMeanCharge]"<<endl;
    cout<<hCellsBiasColumnOverlayAvrgCharge<<endl;
    if(hCellsBiasColumnOverlayAvrgCharge){
        cout<<hCellsBiasColumnOverlayAvrgCharge->IsZombie()<<endl;
        cout<<hCellsBiasColumnOverlayAvrgCharge->GetEntries()<<endl;
        TString name = "hCellsBiasColumnOverlayAvrgCharge_cl";
        name.Append(appendix);
        TProfile2D* histo = (TProfile2D*)hCellsBiasColumnOverlayAvrgCharge->Clone(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolz");
        histo->GetZaxis()->SetRangeUser(0,1200);
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    if(hCellsBiasColumnOverlayAvrgChargeMinusBadCells){
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCells->IsZombie()<<endl;
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCells->GetEntries()<<endl;
        TString name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCells_cl";
        name.Append(appendix);
        TProfile2D* histo = (TProfile2D*)hCellsBiasColumnOverlayAvrgChargeMinusBadCells->Clone(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolzTEXT");
        histo->GetZaxis()->SetRangeUser(0,1200);
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    if(hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents){
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->IsZombie()<<endl;
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->GetEntries()<<endl;
        TString name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents_cl";
        name.Append(appendix);
        TH2F* histo = (TH2F*)hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents->Clone(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolz");
        histo->GetZaxis()->SetRangeUser(0,10);
        //histSaver->SaveHistogram(histo);

        TCanvas* cCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut = new TCanvas("cCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut","cCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut");
        cCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->cd();
        histo->Draw("colztext");

        Float_t xLow = settings->getBiasColumnOverlayXLow();
        Float_t xHigh = settings->getBiasColumnOverlayXHigh();
        Float_t yLow = settings->getBiasColumnOverlayYLow();
        Float_t yHigh = settings->getBiasColumnOverlayYHigh();

        TLine* BinLowerEdgeX = new TLine(xLow,yLow-2,xLow,yHigh+2);
        TLine* BinUpperEdgeX = new TLine(xHigh,yLow-2,xHigh,yHigh+2);
        TLine* BinLowerEdgeY = new TLine(xLow-2,yLow,xHigh+2,yLow);
        TLine* BinUpperEdgeY = new TLine(xLow-2,yHigh,xHigh+2,yHigh);

        BinLowerEdgeX->SetLineWidth(2);
        BinLowerEdgeX->SetLineColor(kRed);
        BinLowerEdgeX->Draw("same");
        BinUpperEdgeX->SetLineWidth(2);
        BinUpperEdgeX->SetLineColor(kRed);
        BinUpperEdgeX->Draw("same");
        BinLowerEdgeY->SetLineWidth(2);
        BinLowerEdgeY->SetLineColor(kRed);
        BinLowerEdgeY->Draw("same");
        BinUpperEdgeY->SetLineWidth(2);
        BinUpperEdgeY->SetLineColor(kRed);
        BinUpperEdgeY->Draw("same");

        histSaver->SaveCanvas(cCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut);

        delete histo;
    }

    if(hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut){
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->IsZombie()<<endl;
        cout<<hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->GetEntries()<<endl;
        TString name = "hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut_cl";
        name.Append(appendix);
        TProfile2D* histo = (TProfile2D*)hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut->Clone(name);
        //TH2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->ProjectionXY("B");
        //histo->SetName(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolz");
        histo->GetXaxis()->SetRangeUser(settings->getBiasColumnOverlayXLow(),settings->getBiasColumnOverlayXHigh());
        histo->GetYaxis()->SetRangeUser(settings->getBiasColumnOverlayYLow(),settings->getBiasColumnOverlayYHigh());
        histSaver->SaveHistogram(histo);

        /*name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutScaled_cl";
			name.Append(appendix);
			histo->SetName(name);
			histo->Scale(1/settings->getOverlayColumnPulseHeightCut());
			histo->GetZaxis()->SetRangeUser(0,1);
			histSaver->SaveHistogram(histo);

			TCanvas* cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut = new TCanvas("cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut","cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut");
			cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut->cd();
			histo->Draw("colztext");
			histSaver->SaveCanvas(cCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut);*/

        delete histo;
    }


    //histSaver->SaveHistogram(hCellsCentralColumnOverlayLandau.at(ClusterSize));
    histSaver->SaveHistogram(hCellsBiasColumnOverlayLandauMinusBadCells);


    //} //End of for ClusterSize

    //		hCellsOverlayPulseHeight->Project3D("xy");

    /* hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayBiasColumnAndCentralColumnStack() {

    //The Two histos i want to stack:

    /*hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(3)
	hCellsBiasColumnOverlayAvrgChargeMinusBadCells*/

    /*hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(3)
	hCellsBiasColumnOverlayAvrgChargeMinusBadCellsBelowCut;*/
    TString name = "sCellOverlayColumnBinOverlay";
    name.Append(appendix);
    THStack sCellOverlayColumnBinOverlay(name,name);

    TProfile2D* histo1;
    TProfile2D* histo2;

    /*if(hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells){
		//cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(3)->IsZombie()<<endl;
		//cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(3)->GetEntries()<<endl;
		histo1 = (TProfile2D*)hCellsCentralColumnOverlayShiftedAvrgChargeMinusBadCells->Clone();
		//cout<<"SAVE"<<endl;
		//histo->Draw("goffcolz");
		//histo->GetXaxis()->SetRangeUser(settings->getBiasColumnOverlayXLow(),settings->getBiasColumnOverlayXHigh());
		//histo->GetYaxis()->SetRangeUser(settings->getBiasColumnOverlayYLow(),settings->getBiasColumnOverlayYHigh());
		//histSaver->SaveHistogram(histo);
		histo1->GetZaxis()->SetRangeUser(300,1500);
		sCellOverlayColumnBinOverlay.Add(histo1);

		for(int i=0;i<5;i++){
			for(int j=0;j<5;j++){
				Int_t BinNum = i*5 + j +1;
				cout<<"Bin: "<<BinNum<<" Content: "<<histo1->GetBinContent(BinNum)<<endl;
			}
		}

		histo1->Draw("goffcolz");
		//delete histo;
	}*/
    if(hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells){
        //cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(3)->IsZombie()<<endl;
        //cout<<hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(3)->GetEntries()<<endl;
        histo2 = (TProfile2D*)hCellsBiasColumnOverlayShiftedAvrgChargeMinusBadCells->Clone();
        //cout<<"SAVE"<<endl;
        //histo->Draw("goffcolz");
        //histo->GetXaxis()->SetRangeUser(settings->getBiasColumnOverlayXLow(),settings->getBiasColumnOverlayXHigh());
        //histo->GetYaxis()->SetRangeUser(settings->getBiasColumnOverlayYLow(),settings->getBiasColumnOverlayYHigh());
        //histSaver->SaveHistogram(histo);
        histo2->GetZaxis()->SetRangeUser(300,1500);
        sCellOverlayColumnBinOverlay.Add(histo2);
        //histo2->Draw("goffcolz");
        //delete histo;
    }

    name = "cCellOverlayColumnBinOverlay";
    name.Append(appendix);
    TCanvas cCellOverlayColumnBinOverlay(name,name);
    cCellOverlayColumnBinOverlay.cd();

    sCellOverlayColumnBinOverlay.Draw("colzTEXT");
    //cCellOverlayColumnBinOverlay->Draw("nostack");
    histSaver->SaveCanvas(&cCellOverlayColumnBinOverlay);

    delete histo1;
    delete histo2;

    /*//hCellsCentralColumnOverlayLandauMinusBadCells.at(3)->SetFillColor(2);
	hCellsCentralColumnOverlayLandauMinusBadCells.at(3)->SetLineColor(2);

	hCellsCentralColumnOverlayLandauMinusBadCells.at(3)->Fill(2000);
	hCellsBiasColumnOverlayLandauMinusBadCells->Fill(2000);
	hCellsBiasColumnOverlayLandauMinusBadCells->Fill(2000);*/

    name = "sCellOverlayColumnBinOverlayLandau";
    name.Append(appendix);
    THStack sCellOverlayColumnBinOverlayLandau(name,name);
    sCellOverlayColumnBinOverlayLandau.Add(hCellsCentralColumnOverlayLandauMinusBadCells.at(3));
    sCellOverlayColumnBinOverlayLandau.Add(hCellsBiasColumnOverlayLandauMinusBadCells);

    name = "cCellOverlayColumnBinOverlayLandauCombined";
    name.Append(appendix);
    TCanvas cCellOverlayColumnBinOverlayLandauCombined(name,name);
    cCellOverlayColumnBinOverlayLandauCombined.cd();
    sCellOverlayColumnBinOverlayLandau.Draw();
    histSaver->SaveCanvas(&cCellOverlayColumnBinOverlayLandauCombined);

    hCellsCentralColumnOverlayLandauMinusBadCells.at(3)->SetLineColor(2);
    name = "cCellOverlayColumnBinOverlayLandau";
    name.Append(appendix);
    TCanvas cCellOverlayColumnBinOverlayLandau(name,name);
    cCellOverlayColumnBinOverlayLandau.cd();
    sCellOverlayColumnBinOverlayLandau.Draw("nostack");
    histSaver->SaveCanvas(&cCellOverlayColumnBinOverlayLandau);

}

void TAnalysisOf3dDiamonds::LongAnalysis_Save3DCellOverlayIndividualBinHistos() {
    //	return;
    cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_Save3DCellOverlayIndividualBinHistos]"<<endl;

    Int_t NBins = 225;
    for(Int_t BinNum=0;BinNum<NBins; BinNum++){

        if(BinNum == 126 || BinNum == 220 || BinNum ==20){

            cout<<VecOverlayCellBinHistos.at(BinNum)<<endl;

            if(VecOverlayCellBinHistos.at(BinNum)){
                cout<<VecOverlayCellBinHistos.at(BinNum)->IsZombie()<<endl;
                cout<<VecOverlayCellBinHistos.at(BinNum)->GetEntries()<<endl;
                TString name = TString::Format ("hOverlayCellsAvrgChargeBin%i_cl", BinNum);
                name.Append(appendix);
                TProfile2D* histo = (TProfile2D*)VecOverlayCellBinHistos.at(BinNum)->Clone(name);
                cout<<"SAVE"<<endl;
                histo->Draw("goffcolz");
                histo->GetZaxis()->SetRangeUser(0,1200);
                histSaver->SaveHistogram(histo);
                delete histo;
            }

            histSaver->SaveHistogram(VecOverlayCellBinLandaus.at(BinNum));
        }
        histSaver->SaveHistogram(hOverlayCellBinHitsBelowCut);
        histSaver->SaveHistogram(hOverlayCellBinHits);

        TString name = "hOverlayCellBinHitsBelowCutRelative";
        name.Append(appendix);
        TH1F* hOverlayCellBinHitsBelowCutRelative = (TH1F*)hOverlayCellBinHits->Clone(name);
        *hOverlayCellBinHitsBelowCutRelative = (*hOverlayCellBinHitsBelowCut)/(*hOverlayCellBinHits);

        vector<TCanvas*> ptrCanvas;
        vector<TH1F*> ptrHisto;
        name = "cOverlayCellBinHitsBelowCut";
        name.Append(appendix);
        TCanvas* cOverlayCellBinHitsBelowCut = new TCanvas(name,name);
        name = "cOverlayCellBinHitsBelowCutRelative";
        name.Append(appendix);
        TCanvas* cOverlayCellBinHitsBelowCutRelative = new TCanvas(name,name);
        ptrCanvas.push_back(cOverlayCellBinHitsBelowCut);
        ptrHisto.push_back(hOverlayCellBinHitsBelowCut);
        ptrCanvas.push_back(cOverlayCellBinHitsBelowCutRelative);
        ptrHisto.push_back(hOverlayCellBinHitsBelowCutRelative);

        for(int i=0; i<ptrCanvas.size(); i++){

            //Int_t MaxEntries = hOverlayCellBinHitsBelowCut->GetMaximum();  //Bin()->GetEntries();
            ptrCanvas.at(i)->cd();
            ptrHisto.at(i)->Draw();

            Int_t BinEntries = ptrHisto.at(i)->GetBinContent(1);
            TLine* Bias01LowerEdge = new TLine(0,0,0,BinEntries);
            TLine* Bias01UpperEdge = new TLine(1,0,1,BinEntries);
            BinEntries = ptrHisto.at(i)->GetBinContent(15);
            TLine* Bias02LowerEdge = new TLine(14,0,14,BinEntries);
            TLine* Bias02UpperEdge = new TLine(15,0,15,BinEntries);

            BinEntries = ptrHisto.at(i)->GetBinContent(113);
            TLine* ReadoutLowerEdge = new TLine(112,0,112,BinEntries);
            TLine* ReadoutUpperEdge = new TLine(113,0,113,BinEntries);

            BinEntries = ptrHisto.at(i)->GetBinContent(211);
            TLine* Bias03LowerEdge = new TLine(210,0,210,BinEntries);
            TLine* Bias03UpperEdge = new TLine(211,0,211,BinEntries);
            BinEntries = ptrHisto.at(i)->GetBinContent(225);
            TLine* Bias04LowerEdge = new TLine(224,0,224,BinEntries);
            TLine* Bias04UpperEdge = new TLine(225,0,225,BinEntries);

            Bias01LowerEdge->SetLineWidth(.2);
            Bias01LowerEdge->SetLineColor(kRed);
            Bias01LowerEdge->Draw("same");
            Bias01UpperEdge->SetLineWidth(.2);
            Bias01UpperEdge->SetLineColor(kRed);
            Bias01UpperEdge->Draw("same");

            Bias02LowerEdge->SetLineWidth(.2);
            Bias02LowerEdge->SetLineColor(kRed);
            Bias02LowerEdge->Draw("same");
            Bias02UpperEdge->SetLineWidth(.2);
            Bias02UpperEdge->SetLineColor(kRed);
            Bias02UpperEdge->Draw("same");

            ReadoutLowerEdge->SetLineWidth(.2);
            ReadoutLowerEdge->SetLineColor(kGreen);
            ReadoutLowerEdge->Draw("same");
            ReadoutUpperEdge->SetLineWidth(.2);
            ReadoutUpperEdge->SetLineColor(kGreen);
            ReadoutUpperEdge->Draw("same");

            Bias03LowerEdge->SetLineWidth(.2);
            Bias03LowerEdge->SetLineColor(kRed);
            Bias03LowerEdge->Draw("same");
            Bias03UpperEdge->SetLineWidth(.2);
            Bias03UpperEdge->SetLineColor(kRed);
            Bias03UpperEdge->Draw("same");

            Bias04LowerEdge->SetLineWidth(.2);
            Bias04LowerEdge->SetLineColor(kRed);
            Bias04LowerEdge->Draw("same");
            Bias04UpperEdge->SetLineWidth(.2);
            Bias04UpperEdge->SetLineColor(kRed);
            Bias04UpperEdge->Draw("same");

            histSaver->SaveCanvas(ptrCanvas.at(i));
        } //End of for ptrCanvas.size()
    }	//End of if BinNum ==

    histSaver->SaveHistogram(hOverlayCellUnEvenBinningBinHitsBelowCut);
    histSaver->SaveHistogram(hOverlayCellUnEvenBinningBinHits);

    TString name = "hOverlayCellUnEvenBinningBinHitsBelowCutRelative";
    name.Append(appendix);
    TH1F* hOverlayCellBinHitsBelowCutRelative = (TH1F*)hOverlayCellUnEvenBinningBinHits->Clone(name);
    *hOverlayCellBinHitsBelowCutRelative = (*hOverlayCellUnEvenBinningBinHitsBelowCut)/(*hOverlayCellUnEvenBinningBinHits);

    name = "cOverlayCellUnEvenBinningBinHitsBelowCutRelative";
    name.Append(appendix);
    TCanvas* cOverlayCellUnEvenBinningBinHitsBelowCutRelative = new TCanvas(name,name);
    cOverlayCellUnEvenBinningBinHitsBelowCutRelative->cd();
    hOverlayCellBinHitsBelowCutRelative->Draw();

    Float_t BinEntries = hOverlayCellBinHitsBelowCutRelative->GetBinContent(52);
    cout<<"Bin: 54 Content: "<<BinEntries<<endl;
    TLine* Bias01LowerEdge = new TLine(51,0,51,BinEntries);
    TLine* Bias01UpperEdge = new TLine(52,0,52,BinEntries);
    BinEntries = hOverlayCellBinHitsBelowCutRelative->GetBinContent(188);
    cout<<"Bin: 188 Content: "<<BinEntries<<endl;
    TLine* ReadoutLowerEdge = new TLine(187,0,187,BinEntries);
    TLine* ReadoutUpperEdge = new TLine(188,0,188,BinEntries);

    Bias01LowerEdge->SetLineWidth(.2);
    Bias01LowerEdge->SetLineColor(kRed);
    Bias01LowerEdge->Draw("same");
    Bias01UpperEdge->SetLineWidth(.2);
    Bias01UpperEdge->SetLineColor(kRed);
    Bias01UpperEdge->Draw("same");

    ReadoutLowerEdge->SetLineWidth(.2);
    ReadoutLowerEdge->SetLineColor(kGreen);
    ReadoutLowerEdge->Draw("same");
    ReadoutUpperEdge->SetLineWidth(.2);
    ReadoutUpperEdge->SetLineColor(kGreen);
    ReadoutUpperEdge->Draw("same");

    histSaver->SaveCanvas(cOverlayCellUnEvenBinningBinHitsBelowCutRelative);

    /*Int_t YBins = hCellsOverlayAvrgChargeMinusBadCells.at(0)->GetYaxis()->GetNbins();
				Int_t BinNum = BinX*YBins + BinX;
				VecOverlayCellBinHistos.at(BinNum)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
				VecOverlayCellBinLandaus.at(BinNum)->Fill(clusterCharge);*/

    //} //End of for ClusterSize

    //		hCellsOverlayPulseHeight->Project3D("xy");

    /* hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::LongAnalysis_Save3D3DOffsetOverlayBiasColumnAlignment() {

    cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_Save3D3DOffsetOverlayBiasColumnAlignment]"<<endl;
    vector<TH1F*> hOverlayCellOffsetAlignmentBinHitsBelowCutRelative;
    vector<TCanvas*> cOverlayCellOffsetAlignmentBinHitsBelowCutRelative;
    vector<Float_t> RelativeBinEntriesBelowCutVec;

    Int_t Alignments = ShiftX.size()*ShiftY.size();
    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut
    TString name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut";
    name.Append(appendix);
    TProfile2D* hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut = new TProfile2D(name,name,
            Alignments,0,Alignments,
            4,1,5);
    /*hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.push_back(new TProfile2D(name,name,
						nXbins,xBinEdges,
						nYbins,yBinEdges));*/
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->SetTitle("Alignment");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetYaxis()->SetTitle("Bias Column Bin");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetZaxis()->SetTitle("Relative Entries Below Cut");
    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->ChooseTimeFormat("Y");
    //	hCellsOverlayAvrgCharge->SetContour(99);

    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS
    name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS";
    name.Append(appendix);
    TH1F* hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS = new TH1F(name,name,
            Alignments,0,Alignments);
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS->GetXaxis()->SetTitle("Alignment");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS->GetYaxis()->SetTitle("RMS Between Bias Bins");

    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut
    name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut";
    name.Append(appendix);
    TH1F* hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut = new TH1F(name,name,
            Alignments,0,Alignments);
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut->GetXaxis()->SetTitle("Alignment");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut->GetYaxis()->SetTitle("Total Events Below Cut");

    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut
    name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut";
    name.Append(appendix);
    TH1F* hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut = new TH1F(name,name,
            Alignments,0,Alignments);
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut->GetXaxis()->SetTitle("Alignment");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut->GetYaxis()->SetTitle("Relative Events Below Cut");

    /*cout<<"hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.size(): "<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.size()<<endl;
	cout<<"hOverlayCellOffsetAlignmentBinHits.size(): "<<hOverlayCellOffsetAlignmentBinHits.size()<<endl;
	cout<<"hOverlayCellOffsetAlignmentBinHitsBelowCut.size(): "<<hOverlayCellOffsetAlignmentBinHitsBelowCut.size()<<endl;*/

    //for(int Alignment=0; Alignment<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.size(); Alignment++){

    Float_t TotalBiasEntries=0;
    Float_t TotalBiasEntriesBelowCut=0;
    Float_t RelativeBiasEntriesBelowCut=0;

    string plots_path = histSaver->GetPlotsPath();
    HistogrammSaver newHistSaver(settings);
    newHistSaver.SetPlotsPath(plots_path+(string)"/OverlayAlignment/");
    for(int i=0; i<ShiftX.size(); i++){
        Float_t OffsetX = settings->getOverlayOffsetX() + ShiftX[i];
        for(int j=0; j<ShiftY.size(); j++){
            Float_t OffsetY = settings->getOverlayOffsetY() + ShiftY[j];
            Int_t Alignment = i*ShiftY.size() + j;
            TotalBiasEntries=0;
            TotalBiasEntriesBelowCut=0;
            RelativeBiasEntriesBelowCut=0;

            //			cout<<"Alignment: "<<Alignment<<endl;
            //cout<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)<<endl;

            if(hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)){
                cout<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)->IsZombie()<<endl;
                cout<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)->GetEntries()<<endl;
                TString name = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)->GetName();
                //name.Append(appendix);
                TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignment.at(Alignment)->Clone(name);
                cout<<"SAVE"<<endl;
                histo->Draw("goffcolz");
                histo->GetZaxis()->SetRangeUser(700,1200);
                newHistSaver.SaveHistogram(histo);
                delete histo;
            }
            newHistSaver.SaveHistogram(hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Alignment));
            newHistSaver.SaveHistogram(hOverlayCellOffsetAlignmentBinHits.at(Alignment));

            TString name = TString::Format("%s_Relative", hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Alignment)->GetName());
            //name.Append(appendix);
            hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.push_back((TH1F*)hOverlayCellOffsetAlignmentBinHitsBelowCut.at(0)->Clone(name));
            *hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment) = (*hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Alignment))/(*hOverlayCellOffsetAlignmentBinHits.at(Alignment));

            name = TString::Format("c_%s", hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetName());
            //name.Append(appendix);
            cOverlayCellOffsetAlignmentBinHitsBelowCutRelative.push_back(new TCanvas(name,name));

            //Int_t MaxEntries = hOverlayCellBinHitsBelowCut->GetMaximum();  //Bin()->GetEntries();
            cOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->cd();
            hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->Draw();

            TString Binlabel = TString::Format("%.0f,%.0f", ShiftX[i], ShiftY[j]);
            /*if (Alignment % 2== 0){
				cout<<"Even: "<<Alignment<<endl;
				hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->GetBin(Alignment)->SetLabelOffset(0.015);
			}
			else{
				hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->GetBin(Alignment)->SetLabelOffset(0);
			}*/
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->SetBinLabel(Alignment,Binlabel);
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->SetLabelSize(0.035);

            vector<float> BiasColumnBinsVec;
            float BiasColumnBins[4] = {33,34,48,49};

            for(int i=0;i<4;i++){
                BiasColumnBinsVec.push_back(BiasColumnBins[i]);
            }

            RelativeBinEntriesBelowCutVec.empty(); RelativeBinEntriesBelowCutVec.clear();
            for(int BiasColumn =0; BiasColumn<BiasColumnBinsVec.size(); BiasColumn++){
                int BiasColumnBin = BiasColumnBinsVec.at(BiasColumn);

                Float_t BinEntries = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetBinContent(BiasColumnBin);

                RelativeBinEntriesBelowCutVec.push_back(BinEntries);
                TotalBiasEntriesBelowCut += hOverlayCellOffsetAlignmentBinHitsBelowCut.at(Alignment)->GetBinContent(BiasColumnBin);
                TotalBiasEntries += hOverlayCellOffsetAlignmentBinHits.at(Alignment)->GetBinContent(BiasColumnBin);

                hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->Fill(Alignment,BiasColumn+1,BinEntries);

                Double_t LowEdge = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetXaxis()->GetBinLowEdge(BiasColumnBin);
                Double_t UpperEdge = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetXaxis()->GetBinUpEdge(BiasColumnBin);
                TLine* BiasLowerEdge = new TLine(LowEdge,0,LowEdge,BinEntries);
                TLine* BiasUpperEdge = new TLine(UpperEdge,0,UpperEdge,BinEntries);

                //TLine* BiasLowerEdge = new TLine(LowEdge,0,LowEdge,.5);
                //TLine* BiasUpperEdge = new TLine(UpperEdge,0,UpperEdge,.5);

                BiasLowerEdge->SetLineWidth(.2);
                BiasLowerEdge->SetLineColor(kRed);
                BiasLowerEdge->Draw("same");
                BiasUpperEdge->SetLineWidth(.2);
                BiasUpperEdge->SetLineColor(kRed);
                BiasUpperEdge->Draw("same");
            }

            Float_t BinEntries = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetBinContent(161);
            Double_t LowEdge = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetXaxis()->GetBinLowEdge(161);
            Double_t UpperEdge = hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment)->GetXaxis()->GetBinUpEdge(161);

            TLine* ReadoutLowerEdge = new TLine(LowEdge,0,LowEdge,BinEntries);
            TLine* ReadoutUpperEdge = new TLine(UpperEdge,0,UpperEdge,BinEntries);

            ReadoutLowerEdge->SetLineWidth(.2);
            ReadoutLowerEdge->SetLineColor(kGreen);
            ReadoutLowerEdge->Draw("same");
            ReadoutUpperEdge->SetLineWidth(.2);
            ReadoutUpperEdge->SetLineColor(kGreen);
            ReadoutUpperEdge->Draw("same");

            RelativeBiasEntriesBelowCut = TotalBiasEntriesBelowCut/TotalBiasEntries;
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut->Fill(Alignment,TotalBiasEntriesBelowCut);
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut->Fill(Alignment,RelativeBiasEntriesBelowCut);
            Float_t RMS = LongAnalysis_CalculateRMS(&RelativeBinEntriesBelowCutVec);
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS->Fill(Alignment,RMS);

            newHistSaver.SaveCanvas(cOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment));
            LongAnalysis_Fill2DCellHitsBelowCutRelative(hOverlayCellOffsetAlignmentBinHitsBelowCutRelative.at(Alignment), Alignment);

        } //End of for Alignment for i
    } //End of for Alignment for j

    //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetXaxis()->CenterLabels();

    if(hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut){
        cout<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->IsZombie()<<endl;
        cout<<hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->GetEntries()<<endl;
        TString name = "hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut";
        name.Append(appendix);
        TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCut->Clone(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolz");
        //histo->GetZaxis()->SetRangeUser(0,1);
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    histSaver->SaveHistogram(hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentToTalBiasEntriesBelowCut);
    histSaver->SaveHistogram(hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeBiasEntriesBelowCut);
    histSaver->SaveHistogram(hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentBiasEntriesBelowCutRMS);

}

Float_t TAnalysisOf3dDiamonds::LongAnalysis_CalculateRMS(vector<Float_t>* nVector){

    Float_t SumSquares = 0;

    for(int i=0; i<nVector->size(); i++){
        SumSquares += (nVector->at(i))*(nVector->at(i));
    }

    SumSquares /= nVector->size();

    return sqrt(SumSquares);

}

void TAnalysisOf3dDiamonds::LongAnalysis_Fill2DCellHitsBelowCutRelative(TH1F* nHisto, Int_t Alignment){

    Int_t NXBins = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->GetXaxis()->GetNbins();
    Int_t NYBins = hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->GetYaxis()->GetNbins();
    Int_t NBins = NXBins*NYBins;
    if(verbosity>6)
        printf("NXBins: %i, NYBins: %i, NBins: %i \n", NXBins, NYBins, NBins);

    /*for(int Bin=0; Bin<NBins; Bin++){
		Float_t BinEntries = nHisto->GetBinContent(Bin);
		cout<<"Bin Entries: "<<BinEntries<<endl;
		hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->SetBinContent(Bin,BinEntries);
	}*/

    for(int BinX=0; BinX<NXBins; BinX++){
        for(int BinY=0; BinY<NYBins; BinY++){
            int Bin = BinX*NYBins+(BinY+1);
            //int Bin = (BinX+1)*(BinY+1)-1;
            //if(Bin == 33 || Bin == 34 || Bin == 48 || Bin == 49 || Bin == 161){
            float BinWidth = settings->getPitchWidth(subjectDetector,2)/NXBins;
            //cout<<"BinWidth: "<<BinWidth<<endl;
            float FillX = BinWidth*BinX + BinWidth/2;
            float FillY = BinWidth*BinY + BinWidth/2;
            Float_t BinEntries = nHisto->GetBinContent(Bin);
            //Float_t BinEntries = 1;
            //cout<<"Bin Entries: "<<BinEntries<<endl;
            hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->Fill(FillX,FillY,BinEntries);
            //}
            //hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->SetBinContent(Bin,BinEntries);
        }
    }

    TString name = TString::Format("c_%s",hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->GetName());
    TCanvas* c1 = new TCanvas(name,name);
    c1->cd();
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->Draw("TEXTcolz");
    hCellsOffsetOverlayAvrgChargeMinusBadCellsAlignmentRelativeEventsBelowCut.at(Alignment)->GetZaxis()->SetRangeUser(0.05,0.3);
    TCutG* BiasEdges = new TCutG("Bias",5);
    BiasEdges->SetPoint(0,20,20);
    BiasEdges->SetPoint(1,20,40);
    BiasEdges->SetPoint(2,40,40);
    BiasEdges->SetPoint(3,40,20);
    BiasEdges->SetPoint(4,20,20);
    BiasEdges->SetLineWidth(2);
    BiasEdges->SetLineColor(kRed);
    BiasEdges->Draw("same");

    TCutG* ReadoutEdges = new TCutG("Readout",5);
    ReadoutEdges->SetPoint(0,100,100);
    ReadoutEdges->SetPoint(1,100,110);
    ReadoutEdges->SetPoint(2,110,110);
    ReadoutEdges->SetPoint(3,110,100);
    ReadoutEdges->SetPoint(4,100,100);
    ReadoutEdges->SetLineWidth(2);
    ReadoutEdges->SetLineColor(kGreen);
    ReadoutEdges->Draw("same");

    histSaver->SaveCanvas(c1);
    delete c1;

}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayOffsetMeanCharge() {
    //	return;

    Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
    for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
        TString appendix2 = TString::Format("_ClusterSize%i", ClusterSize);
        appendix2+=appendix;

        //		cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayOffsetMeanCharge]"<<endl;
        cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)<<endl;
        if(hCellsOffsetOverlayAvrgCharge.at(ClusterSize)){
            cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOffsetOverlayAvrgCharge_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(700,1200);
            histSaver->SaveHistogram(histo);
            delete histo;
            ////		hCellsOverlayAvrgCharge->SetName("hCellsOverlayAvrgCharge");
            //		cout<<"Set Name: "<<hCellsOverlayAvrgCharge<<endl;
            ////		hCellsOverlayAvrgCharge->SetTitle("Avrg PH - overlayed");
            //		histSaver->SaveHistogram(hCellsOverlayAvrgCharge);
            //		delete hCellsOverlayAvrgCharge ;
        }
        if(hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)){
            cout<<hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOffsetOverlayAvrgChargeMinusBadCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(700,1200);
            histSaver->SaveHistogram(histo);
            delete histo;
        }
        if(hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)){
            cout<<hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->IsZombie()<<endl;
            cout<<hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->GetEntries()<<endl;
            TString name = "hCellsOffsetOverlayAvrgChargeGoodCells_cl";
            name.Append(appendix2);
            TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->Clone(name);
            cout<<"SAVE"<<endl;
            histo->Draw("goffcolz");
            histo->GetZaxis()->SetRangeUser(700,1200);
            histSaver->SaveHistogram(histo);
            delete histo;
        }

    } //End of for ClusterSize

    if(hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells){
        cout<<hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->IsZombie()<<endl;
        cout<<hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->GetEntries()<<endl;
        TString name = "hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells_cl";
        name.Append(appendix);
        TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayUnEvenBinningAvrgChargeMinusBadCells->Clone(name);
        cout<<"SAVE"<<endl;
        histo->Draw("goffcolz");
        histo->GetZaxis()->SetRangeUser(700,1200);
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    //		hCellsOverlayPulseHeight->Project3D("xy");

    /*hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::HitandSeedCount(TCluster* nCluster) {
    int Hit=0;int Seed=0;
    for (UInt_t i=0;i<nCluster->getClusterSize();i++){
        if(nCluster->isHit(i)) Hit++;
        if(nCluster->isSeed(i)) Seed++;
    }
    HitCount=(Hit-Seed);
    SeedCount=Seed;
}

void TAnalysisOf3dDiamonds::ClusterPlots(int nClusters, float nfiducialValueX, float nfiducialValueY) {
    if(nEvent >= 519000 &&nEvent<=520000 && settings->getRunNumber() == 17212)
        cout<<"ClusterPlots: "<<nClusters<<endl;
    if(nClusters==0){
        if(hFidCutXvsFidCutYClusters[0])
            hFidCutXvsFidCutYClusters.at(0)->Fill(nfiducialValueX,nfiducialValueY,1);
    }
    if(nClusters==1){
        if(hFidCutXvsFidCutYClusters[1])
            hFidCutXvsFidCutYClusters.at(1)->Fill(nfiducialValueX,nfiducialValueY,1);
    }
    if(nClusters==1){
        if(HitCount==0&&SeedCount==1){
            if(hFidCutXvsFidCutYClusters[2])
                hFidCutXvsFidCutYClusters.at(2)->Fill(nfiducialValueX,nfiducialValueY,1);
        }
    }
    if(nClusters==2){
        if(hFidCutXvsFidCutYClusters[3])
            hFidCutXvsFidCutYClusters.at(3)->Fill(nfiducialValueX,nfiducialValueY,1);
        TCluster diamondCluster0 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
        TCluster diamondCluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);

        hDoubleClusterPos->Fill(diamondCluster0.getHighestSignalChannel());
        hDoubleClusterPos->Fill(diamondCluster1.getHighestSignalChannel());
        if(diamondCluster0.getHighestSignalChannel()==85||diamondCluster1.getHighestSignalChannel()==85){
            hDoubleClusterPos0->Fill(diamondCluster0.getHighestSignalChannel());
            hDoubleClusterPos0->Fill(diamondCluster1.getHighestSignalChannel());
            hLandauCluster1->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
            hFidCutXvsFidCutYClusters.at(4)->Fill(nfiducialValueX,nfiducialValueY,1);
        }
        if(diamondCluster0.getHighestSignalChannel()==55||diamondCluster1.getHighestSignalChannel()==55){
            hDoubleClusterPos1->Fill(diamondCluster0.getHighestSignalChannel());
            hDoubleClusterPos1->Fill(diamondCluster1.getHighestSignalChannel());
            hLandauCluster2->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
            hFidCutXvsFidCutYClusters.at(5)->Fill(nfiducialValueX,nfiducialValueY,1);
        }
        if((!diamondCluster0.getHighestSignalChannel()==55&&!diamondCluster1.getHighestSignalChannel()==55)||(!diamondCluster0.getHighestSignalChannel()==85&&!diamondCluster1.getHighestSignalChannel()==85)){
            hLandauDoubleCombined->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
        }

    }
    if(nClusters==3){
        hFidCutXvsFidCutYClusters.at(6)->Fill(nfiducialValueX,nfiducialValueY,1);
    }
}


void TAnalysisOf3dDiamonds::createTreeTestHistos() {
    TH1F* histo = (TH1F*) clusteredAnalysis->getHistogram("hTest","pulseHeight","","");
    TH3F* histo2 = (TH3F*) clusteredAnalysis->getHistogram("hMeanPH","pulseHeight:nRow:nColumn","","");
    if(histo2){
        TH2F* hAvrgCharge = (TH2F*) histo2->Project3DProfile("yx");
        if(hAvrgCharge)hAvrgCharge->SetName("hAvrgChargeInCells");
        histSaver->SaveHistogram(hAvrgCharge);
        TCanvas *c1 = histSaver->DrawHistogramWithCellGrid(hAvrgCharge);
        histSaver->SaveCanvas(c1);


    }
    else{
        cout<<"PROBLEM: "<<histo2<<endl;
    }
    histSaver->SaveHistogram(histo);
}

/**
 * Cell Labeling
 * 			+-----------+-----------+
 * 			+           +           +
 * 			+     1     +     3     +       ^
 * 			+           +           +     Y |
 * 			+-----------+-----------+       |
 * 			+           +           +       |
 * 			+     0     +     2     +       |
 * 			+           +           +       |
 * 			+-----------+-----------+       +-------->
 * 			                                       X
 * @param xDet
 * @param yDet
 * @return
 */

void TAnalysisOf3dDiamonds::ShortAnalysis_FillEdgeDistributions(Float_t clusterCharge){
    for(UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts();i++ ){
        TFiducialCut* fidCut = settings->get3dEdgeFidCuts()->getFidCut(i+1);
        if(!fidCut)
            continue;
        if(fidCut->IsInFiducialCut(fiducialValueX,fiducialValueY)){
            vecEdgePredX[i].push_back(xPredDet);
            vecEdgePredY[i].push_back(yPredDet);
            vecEdgePulseHeight[i].push_back(clusterCharge);
            //				cout<<nEvent<<": "<<xPredDet<<"/"<<yPredDet<<" --> "<<i<<endl;
        }
    }

}

void TAnalysisOf3dDiamonds::ShortAnalysis_SaveEdgeDistributions() {
}

void TAnalysisOf3dDiamonds::ShortAnalysis_FillMeanChargeVector(
        Float_t clusterCharge) {
    vecPredDetX_ShortAna.push_back(xPredDet);
    vecPredDetY_ShortAna.push_back(yPredDet);
    vecPulseHeight_ShortAna.push_back(clusterCharge);

    Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
    if (predictedDetector == 2){
        hLandau3DPhantom->Fill(clusterCharge);
        hLandau3DPhantomFidCutXvsFidCutY->Fill(fiducialValueX,fiducialValueY);
        if(settings->centralRegion3DnH->IsInFiducialCut(xPredDet,yPredDet))
            hLandau3DPhantomCentral->Fill(clusterCharge);
    }
    else if (predictedDetector == 3){
        hLandau3DWithColumns->Fill(clusterCharge);
        hLandau3DWithColumnsFidCutXvsFidCutY->Fill(fiducialValueX,fiducialValueY);
    }

}

void TAnalysisOf3dDiamonds::ShortAnalysis_Save2ClusterPlots() {
    TH2F * hPH = histSaver->CreateScatterHisto("hPulseHeightComparision2Clusters",vecPH_Cluster2_ShortAna,vecPH_Cluster1_ShortAna,
            PulseHeightBins,PulseHeightBins,PulseHeightMin,PulseHeightMax,PulseHeightMin-1000,PulseHeightMax-1000);
    hPH->GetXaxis()->SetTitle("Pulse height cluster no. 1");
    hPH->GetYaxis()->SetTitle("Pulse height cluster no. 2");
    histSaver->SaveHistogram(hPH);
    delete hPH;

    TH2F* hCh = histSaver->CreateScatterHisto("hChannelComparision2Clusters",vecCh_Cluster2_ShortAna,vecCh_Cluster1_ShortAna,128,128,0,128,0,128);
    hCh->GetXaxis()->SetTitle("Channel no for cluster no. 1");
    hCh->GetYaxis()->SetTitle("Channel no for cluster no. 2");
    hCh->Draw("colz");
    //	Float_t xmax = hCh->GetZaxis()->GetXmax();
    //	Float_t xmin = hCh->GetZaxis()->GetXmin();
    //	hCh->GetZaxis()->SetRangeUser(xmin+.1*(xmax-xmin),xmax);
    histSaver->SaveHistogram(hCh,false);
    delete hCh;
}

void TAnalysisOf3dDiamonds::initialiseHistos() {
    // Initialise
    //Universal histograms

    TString name = "hValidEventsFiducialSpace";
    name.Append(appendix);
    hValidEventsFiducialSpace = new TH2F(name,name,1024,0,256,1024,0,256);
    hValidEventsFiducialSpace->GetXaxis()->SetTitle("avrg silicon position x /ch");
    hValidEventsFiducialSpace->GetYaxis()->SetTitle("avrg silicon position y /ch");
    hValidEventsFiducialSpace->GetZaxis()->SetTitle("number of entries #");
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
    name = "hValidEventsDetSpace";
    name.Append(appendix);
    hValidEventsDetSpace = new TH2F(name,name,1024,xmin,xmax,1024,ymin,ymax);
    hValidEventsDetSpace->GetXaxis()->SetTitle("#it{X} / #mum");
    hValidEventsDetSpace->GetYaxis()->SetTitle("#it{Y} / #mum");
    hValidEventsDetSpace->GetZaxis()->SetTitle("number of entries #");

    name = "hClusterEventsDetSpace";
    name.Append(appendix);
    hClusterEventsDetSpace = new TH2F(name,name,1024,xmin,xmax,1024,ymin,ymax);
    hClusterEventsDetSpace->GetXaxis()->SetTitle("#it{X} / #mum");
    hClusterEventsDetSpace->GetYaxis()->SetTitle("#it{Y} / #mum");
    hClusterEventsDetSpace->GetZaxis()->SetTitle("number of entries #");

    InitialiseStripAnalysisHistos();
    name = "hLandau3D";
    name.Append(appendix);
    hLandau3DWithColumns = new TH1F(name,"3D",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DWithColumns->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DWithColumns->GetYaxis()->SetTitle("number of entries #");
    hLandau3DWithColumns->SetLineColor(kBlack);

    name = "hLandau3DPhantom";
    name.Append(appendix);
    hLandau3DPhantom = new TH1F(name,"3D Phantom",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DPhantom->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DPhantom->GetYaxis()->SetTitle("number of entries #");
    hLandau3DPhantom->SetLineColor(kRed);
    hLandau3DPhantom->SetLineStyle(7);

    name = "hLandau3DPhantom_subset";
    name.Append(appendix);
    hLandau3DPhantomCentral = new TH1F(name,"3D Phantom, central Region",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DPhantomCentral->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DPhantomCentral->GetYaxis()->SetTitle("number of entries #");
    hLandau3DPhantomCentral->SetLineColor(kRed);
    hLandau3DPhantomCentral->SetLineWidth(2);

    name = "hLandau3DWithColumnsFidCutXvsFidCutY";
    hLandau3DWithColumnsFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandau3DWithColumnsFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandau3DWithColumnsFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");

    name = "hLandau3DPhantomFidCutXvsFidCutY";
    hLandau3DPhantomFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandau3DPhantomFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandau3DPhantomFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");

    if(settings->do3dShortAnalysis() == 1){
        initialiseShortAnalysisHistos();
    }
    if(settings->do3dLongAnalysis() == 1){
        initialiseLongAnalysisHistos();
    }
    if(settings->do3dTransparentAnalysis() == 1){
        initialiseTransparentAnalysisHistos();
        //initialiseTransparentAnalysisHistos();
    }
}

void TAnalysisOf3dDiamonds::saveHistos() {
//    if (settings->do3dTransparentAnalysis())

    histSaver->SaveHistogram(hValidEventsDetSpace);
    hValidEventsDetSpace->SetName("hValidEventsDetSpaceGrid"+appendix);
    hValidEventsDetSpace->Rebin2D(2,2);
    histSaver->SaveHistogramWithCellGrid(hValidEventsDetSpace);
    hValidEventsDetSpace->SetName("hValidEventsDetSpaceGridRebinned"+appendix);
    hValidEventsDetSpace->Rebin2D(2,2);
    histSaver->SaveHistogramWithCellGrid(hValidEventsDetSpace);

    histSaver->SaveHistogram(hClusterEventsDetSpace);
    hClusterEventsDetSpace->Rebin2D(2,2);
    hClusterEventsDetSpace->SetName("hClusterEventsDetSpaceRebin"+appendix);
    histSaver->SaveHistogram(hClusterEventsDetSpace);
    hClusterEventsDetSpace->SetName("hClusterEventsDetSpaceGrid"+appendix);
    histSaver->SaveHistogramWithCellGrid(hClusterEventsDetSpace);
    hClusterEventsDetSpace->SetName("hClusterEventsDetSpaceGridRebinned"+appendix);
    hClusterEventsDetSpace->Rebin2D(2,2);
    histSaver->SaveHistogramWithCellGrid(hClusterEventsDetSpace);

    histSaver->SaveHistogram(hValidEventsFiducialSpace);
    if(settings->do3dLongAnalysis() == 1){SaveLongAnalysisHistos();}

    LongAnalysis_CreateResolutionPlots();
    LongAnalysis_SaveChargeSharingPlots();
    SaveStripAnalysisHistos();
    // Save
    if(settings->do3dShortAnalysis() == 1){SaveShortAnalysisHistos();}
    if(settings->do3dTransparentAnalysis() == 1){saveTransparentAnalysisHistos();/*saveTransparentAnalysisHistos()*/}
}

void TAnalysisOf3dDiamonds::initialiseLongAnalysisHistos() {
    mapPredictedPositionsGoodCells.clear();
    mapClusteredAnalysisGoodCells.clear();
    mapTransparentAnalysisGoodCells.clear();
    mapPredictedPositionsAllCells.clear();
    mapClusteredAnalysisAllCells.clear();
    mapTransparentAnalysisAllCells.clear();
    mapPredictedPositionsAllButBadCells.clear();
    mapClusteredAnalysisAllButBadCells.clear();
    mapTransparentAnalysisAllButBadCells.clear();
    initialise3DYAlignmentHistos();
    initialise3DOverviewHistos();
    initialise3DCellOverlayHistos();
    initialise3DCellCentralColumnOverlayHistos();
    initialise3DCellBiasColumnOverlayHistos();
    initialise3DCellOverlayIndividualBinHistos();
    initialise3DOffsetOverlayHistos();
    initialise3DOffsetAlignmentOverlayHistos();
    initialiseEdgeFreeHistos();
    LongAnalysis_InitialiseRelativeAddedTransparentCharge();
    LongAnalysis_InitGoodCellsLandaus();

//    if (settings->do3dTransparentAnalysis())
    LongAnalysis_InitResolutionPlots();
    LongAnalysis_InitChargeSharingPlots();
    hLongAnalysisInvalidCellNo = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCellNo"+appendix);
    hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCellNo");
    hLongAnalysisInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCluster"+appendix);
    hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCluster");
    hNegativeChargePosition = histSaver->GetHistoBinedInCells("hNegativeChargePositionAllCells"+appendix,16);
    hNegativeChargePosition->SetTitle(TString::Format("hNegativeChargePositionAllCells, cut: %+3f",settings->getNegativeChargeCut()));

    TString name = "hNegativeChargeRatio"+appendix;
    TString title = "Negative Charge Ratio "+appendix;
    title+="; signal ratio: S_{Min}/PH; Pulse Heigth / ADC;number of entries";
    hNegativeChargeRatio = new TH2D(name,title,1000,.5,.5,PulseHeightBins/4,PulseHeightMin,PulseHeightMax);

    title = "Negative Charge Ratio "+appendix;
    title+="; signal ratio: S_{Min}/PH; S_{Min} / ADC;number of entries";
    name = "hNegativeChargeRatioAbs"+appendix;
    hNegativeChargeRatioAbs = new TH2D(name,title,1000,-.5,.5,PulseHeightBins/4,-400,400);

    name = "hNegativeChargeRatioMax"+appendix;
    title = "Negative Charge Ratio Max "+appendix;
    title+="; signal ratio: S_{Min}/PH_{max}; Max. Pulse Heigth / ADC;number of entries";
    hNegativeChargeRatioMax = new TH2D(name,title,1000,-.5,.5,PulseHeightBins/4,PulseHeightMin,PulseHeightMax);

    name = "hNegativeChargeRatioOverlay"+appendix;
    hNegativeChargeRatioOverlay = settings->GetOverlayProfile(name);
    hNegativeChargeRatioOverlay->SetZTitle("avrg min. adjacent Signal");

    name = "hAdjacentChannels_SNR"+appendix;
    title = "SNR of adjacent channels";
    hAdjacentChannels_SNR = new TH2F(name,title,160,-30,50,160,-30,50);
    hAdjacentChannels_SNR->GetXaxis()->SetTitle("SNR left strip");
    hAdjacentChannels_SNR->GetYaxis()->SetTitle("SNR right strip");
    hAdjacentChannels_SNR->GetZaxis()->SetTitle("number of entries");

    name = "hNegativeChargeFraction"+appendix;
    title = TString::Format("Negative Charge Fraction: Thr %d",(int)settings->getNegativeChargeCut());
    title+=";has Negative Charge below Thr;number of entries";
    hNegativeChargeFraction = new TH1F(name,title,2,0,1);

    name = "hNegativeChargeFieldWireFraction"+appendix;
    hNegativeChargeFieldWireFraction = histSaver->GetProfileBinnedAroundFieldWires(name,1);
    hNegativeChargeFieldWireFraction->GetZaxis()->SetTitle("Fraction of negative charges");

    name = "hNegativeChargeFieldWirePositions"+appendix;
    hNegativeChargeFieldWirePositions = histSaver->GetHistoBinnedAroundFieldWires(name,30);
    hNegativeChargeFieldWirePositions->GetZaxis()->SetTitle("Positions of entries in hNegativeChargeFieldWireFraction");

    name = "hNegativeChargeFieldWirePositionsOverlay"+appendix;
    hNegativeChargeFieldWirePositionsOverlay  = settings->GetOverlayHisto(name);
    hNegativeChargeRatioOverlay->SetZTitle("Positions of entries in hNegativeChargeFieldWireFraction");

}

void TAnalysisOf3dDiamonds::ShortAnalysis_SaveMeanChargeVector() {
    cout<<"SaveMeanChargeVector"<<endl;
    cout<<"vecPredDetX_ShortAna    "<<vecPredDetX_ShortAna.size()<<endl;
    cout<<"vecPredDetY_ShortAna    "<<vecPredDetY_ShortAna.size()<<endl;
    cout<<"vecPulseHeight_ShortAna "<<vecPulseHeight_ShortAna.size()<<endl;

    if (settings->do3dTransparentAnalysis())

        cout<<"Create Project3dProfile for hChargeDistribution3D ..."<<flush;

    TProfile2D* hMeanCharge = histSaver->CreateProfile2D("hChargeDistribution3D"+appendix,vecPredDetX_ShortAna,vecPredDetY_ShortAna,vecPulseHeight_ShortAna,1024,1024);
    if (!hMeanCharge){
        cerr<<" hChargDistribution3D: was not created:"<<endl;
        return;
    }
    else if (hMeanCharge->GetEntries() == 0){
        cerr<<" hChargDistribution3D: number of entries is 0"<<endl;
        for (int i = 0; i<vecPredDetX_ShortAna.size()&&i<100;i++)
            cout<<TString::Format("%3d %5.1f/%5.1f --> %6.1f",i,vecPredDetX_ShortAna.at(i),vecPredDetY_ShortAna.at(i),vecPulseHeight_ShortAna.at(i))<<endl;
        return;

    }
    hMeanCharge->GetXaxis()->SetTitle("#it{X} / #mum");
    hMeanCharge->GetYaxis()->SetTitle("#it{Y}/ #mum");
    hMeanCharge->GetYaxis()->SetTitleOffset(1.4);
    hMeanCharge->GetZaxis()->SetTitleOffset(1.3);
    cout<<"\t[done]"<<endl;
    //	if(!hMeanCharge)
    //		return;
    //	else
    //		cerr<<" hChargDistribution3D_pfyx: was not created:"<<endl;
    hMeanCharge->GetZaxis()->SetTitle("Avrg. pulse height /ADC");
    hMeanCharge->SetTitle("Avrg. pulse height in detector system");
    hMeanCharge->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    TString name = "hAvrgPulseHeigthDetSystem";
    name.Append(appendix);
    hMeanCharge->SetName(name);
    histSaver->SaveHistogram(hMeanCharge,false);
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer";
    name.Append(appendix);
    TCanvas *c1 = new TCanvas(name, name);
    c1->cd();
    c1->SetRightMargin(.15);
    hMeanCharge->Draw("colz");
    c1->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hMeanCharge->GetListOfFunctions()->FindObject("palette");
//    palette->SetY2NDC(0.7);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
    settings->DrawMetallisationGrid(c1,3);
    hMeanCharge->Draw("colz same");
    histSaver->SaveCanvas(c1);

    TFiducialCut *fidCut3dWithColumns = settings->get3dMetallisationFidCuts()->getFidCut(3);
    cout<<"3d FidCut: "<<endl;
    fidCut3dWithColumns->Print(1);
    cout<<endl;
    Float_t xmin = fidCut3dWithColumns->GetXLow();
    Float_t xmax = fidCut3dWithColumns->GetXHigh();
    Float_t deltaX = TMath::Abs(.05*(xmax-xmin));
    Float_t ymin = fidCut3dWithColumns->GetYLow();
    Float_t ymax = fidCut3dWithColumns->GetYHigh();
    Float_t deltaY = TMath::Abs(.05*(ymax-ymin));
    hMeanCharge->GetXaxis()->SetRangeUser(xmin-deltaX,xmax+deltaX);
    hMeanCharge->GetYaxis()->SetRangeUser(ymin-deltaY,ymax+deltaY);
    hMeanCharge->Draw("colz");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom";
    name.Append(appendix);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom_rebinned"+appendix;
    TProfile2D* hMeanCharge3D = histSaver->CreateProfile2D("hChargeDistribution3D_3D",
              vecPredDetX_ShortAna,vecPredDetY_ShortAna,vecPulseHeight_ShortAna,
              settings->getNColumns3d()*15,settings->getNRows3d()*15,
              xmin,xmax,ymin,ymax
      );
    hMeanCharge3D->Draw("colz");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom_rebinned";
    name.Append(appendix);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    pair<Float_t, Float_t> x = settings->getAllGoodCellsXpos();
    pair<Float_t, Float_t> y = settings->getAllGoodCellsYpos();
    cout<<"X: "<<x.first<<"-"<<x.second<<"   Y: "<<y.first<<"-"<<y.second<<endl;
    hMeanCharge->GetXaxis()->SetRangeUser(x.first,x.second);
    hMeanCharge->GetYaxis()->SetRangeUser(y.first,y.second);
    TH2D* hGridReferenceCellSpace = (TH2D*)histSaver->GetGridReferenceCellSpace()->Clone();
    hGridReferenceCellSpace->GetXaxis()->SetRangeUser(x.first,x.second);
    hGridReferenceCellSpace->GetYaxis()->SetRangeUser(y.first,y.second);
    hGridReferenceCellSpace->Draw("COL");
    c1->Update();
//    hMeanCharge->Rebin2D(4,4);
    hMeanCharge->Draw("sameCOLZ");
    settings->DrawMetallisationGrid(c1,3);
    hGridReferenceCellSpace->Draw("sameCOL");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_ZoomGoodCells"+appendix;
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    hGridReferenceCellSpace->Draw("COL");
    c1->Update();
    hMeanCharge->Rebin2D(2,2);
    hMeanCharge->Draw("sameCOLZ");
    settings->DrawMetallisationGrid(c1,3);
    hGridReferenceCellSpace->Draw("sameCOL");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_ZoomGoodCellsRebinned"+appendix;
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    hGridReferenceCellSpace->Draw("COL");
    c1->Update();
    hMeanCharge->Rebin2D(2,2);
    hMeanCharge->Draw("sameCOLZ");
    settings->DrawMetallisationGrid(c1,3);
    hGridReferenceCellSpace->Draw("sameCOL");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_ZoomGoodCellsRebinned2"+appendix;
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    if (hMeanCharge)
        delete hMeanCharge;
    if(hMeanCharge3D)
        delete hMeanCharge3D;
}

void TAnalysisOf3dDiamonds::InitialiseStripAnalysisHistos() {
    TString name = "hLandauStrip";
    name.Append(appendix);
    hLandauStrip =  new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauStrip->GetXaxis()->SetTitle("charge / ADC");
    hLandauStrip->GetYaxis()->SetTitle("number of entries #");
    hLandauStrip->SetLineColor(kBlue);

    Int_t xbins = 512;
    Float_t xlow = -PulseHeightMax/14;
    Float_t xup = .78*PulseHeightMax/5;
    name = "hLandauStripNegativeCharges"+appendix;
    hLandauStripNegativeCharges = new TH2F(name,name,xbins,xlow,xup,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauStripNegativeCharges->GetXaxis()->SetTitle("SmallesAdjacentNegativeCharge / ADC");
    hLandauStripNegativeCharges->GetYaxis()->SetTitle("charge / ADC");
    hLandauStripNegativeCharges->GetZaxis()->SetTitle("no of entries");

    name = "hLandauStripNegativeChargesFraction"+appendix;
    TString hName = TString::Format("Negative Charge Fraction: Thr %d",(int)settings->getNegativeChargeCutStrip());
    hName += ";has Negative charge below Thr";
    hName += ";number of entries";
    hLandauStripNegativeChargesFraction = new TH1F(name,hName,3,-1.5,1.5);

    Int_t ybins = 7;
    Float_t ylow = -2;
    Float_t yup = 5;
    name = "hLandauStripNegativeChargesClPos"+appendix;
    hLandauStripNegativeChargesClPos = new TH2F(name,name,xbins,xlow,xup,ybins,ylow,yup);//PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauStripNegativeChargesClPos->GetXaxis()->SetTitle("neg charge of cluster / ADC");
    hLandauStripNegativeChargesClPos->GetYaxis()->SetTitle("rel. Pos of negative in transp. cluster");
    hLandauStripNegativeChargesClPos->GetZaxis()->SetTitle("no of entries");

    name = "hLandauStripNegativeChargePosition";
    hName = name + TString::Format(" Thr: %d",int(settings->getNegativeChargeCutStrip()));
    hName+= ";pred pos x / #mum;pred pos y / #mum";
    name+=appendix;
    xlow = settings->get3dMetallisationFidCuts()->getXLow(1);
    xup = settings->get3dMetallisationFidCuts()->getXHigh(1);
    ylow = settings->get3dMetallisationFidCuts()->getYLow(1);
    yup = settings->get3dMetallisationFidCuts()->getYHigh(1);
    ybins = 128;
    xbins = 128;
    hLandauStripNegativeChargePosition  = new TH2F(name,hName,xbins,xlow,xup,ybins,ylow,yup);

    name = "hLandauStripFidCutXvsFidCutY"+appendix;
    hLandauStripFidCutXvsFidCutY = new TProfile2D(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandauStripFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandauStripFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");
    hLandauStripFidCutXvsFidCutY->GetZaxis()->SetTitle("Charge ADC");

    name = "hLandauStripFiducialPosition"+appendix;
    hLandauStripFiducialPosition= new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandauStripFiducialPosition->GetXaxis()->SetTitle("FidCutX");
    hLandauStripFiducialPosition->GetYaxis()->SetTitle("FidCutY");
    hLandauStripFiducialPosition->GetZaxis()->SetTitle("number of entries");
}

void TAnalysisOf3dDiamonds::SaveStripAnalysisHistos() {
    histSaver->SaveHistogram(hLandauStrip);
    histSaver->SaveHistogram(hLandauStripFidCutXvsFidCutY);
    histSaver->SaveHistogram(hLandauStripFiducialPosition);
    if (!settings->do3dTransparentAnalysis())
        return;
    histSaver->SaveHistogram(hLandauStripNegativeChargesClPos);
    histSaver->SaveHistogram(hLandauStripNegativeChargePosition);
    TCanvas *c1 = new TCanvas((TString)"c"+hLandauStripNegativeChargePosition->GetName());
    hLandauStripNegativeChargePosition->Draw("colz");
    settings->DrawMetallisationGrid(c1,1);
    histSaver->SaveCanvas(c1,"cLandauStripNegativeChargePosition");
    histSaver->SaveHistogram(hLandauStripNegativeCharges);
    histSaver->SaveHistogram(hLandauStripNegativeChargesFraction);
    TString name = "hLandauStripNegativeCharges"+appendix+"_px";
    TH1D* px = hLandauStripNegativeCharges->ProjectionX(name);
    Float_t mp = px->GetBinCenter(hLandauStripNegativeCharges->GetMaximumBin());
    Float_t max= hLandauStripNegativeCharges->GetBinContent(hLandauStripNegativeCharges->GetMaximumBin());
    Float_t xmin = hLandauStripNegativeCharges->FindFirstBinAbove(max/2);
    Float_t xmax = hLandauStripNegativeCharges->FindLastBinAbove(max/2);

    TF1* fit = new TF1("fit","gaus",xmin,xmax);
    px->Fit(fit,"RQ","+",xmin,xmax);
    histSaver->SaveHistogramWithFit(px,fit,xmin,xmax,false);
    px->SetName(name+"_logy");
    histSaver->SaveHistogram(px,false,false,true,"logy");
    delete px;
    delete hLandauStripNegativeCharges;
    delete hLandauStripNegativeChargePosition;
    delete hLandauStripNegativeChargesClPos;
    delete hLandauStripFidCutXvsFidCutY;
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveMeanChargePlots() {
    histSaver->SaveHistogram(hPulseHeightVsDetectorHitPostionXY);
    if (settings->do3dTransparentAnalysis()){
        for (UInt_t i = 0; i< hPulseHeightVsDetectorHitPostionXY_trans.size();i++){
            if( hPulseHeightVsDetectorHitPostionXY_trans[i]->GetEntries()){
                histSaver->SaveHistogram(hPulseHeightVsDetectorHitPostionXY_trans.at(i));
                histSaver->SaveHistogramWithCellGridAndMarkedCells(hPulseHeightVsDetectorHitPostionXY_trans.at(i));
            }
            delete hPulseHeightVsDetectorHitPostionXY_trans[i];
        }
    }
    histSaver->SaveHistogramWithCellGrid(hPulseHeightVsDetectorHitPostionXY);
    histSaver->SaveHistogramWithCellGrid(hPulseHeightVsDetectorHitPostionXYGoodCells);

    TString name = "hPulseHeightVsDetectorHitPostionXY_rebinned";
    name.Append(appendix);
    TProfile2D* profRebinned = (TProfile2D*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(2,2,name);
    profRebinned->SetContour(100);
    profRebinned->Draw();
    profRebinned->Draw("colz");
    histSaver->SaveHistogram(profRebinned);

    TCanvas *c1 =histSaver->DrawHistogramWithCellGrid(profRebinned);
    c1->SetName("cProfRebinned"+appendix);
    histSaver->DrawGoodCellsRegion(c1);
    histSaver->SaveCanvas(c1);
    histSaver->SaveHistogramWithCellGrid(profRebinned);
    c1->SetName("cProfRebinnedMarkedCells"+appendix);
    histSaver->AddMarkedCells(c1);

    c1->Update();
    histSaver->AddMarkedCells(c1);
    histSaver->SaveCanvas(c1);



    c1->cd();
    c1->SetRightMargin(.15);
    pair<Float_t, Float_t> x = settings->getAllGoodCellsXpos();
    pair<Float_t, Float_t> y = settings->getAllGoodCellsYpos();
    cout<<"X: "<<x.first<<"-"<<x.second<<"   Y: "<<y.first<<"-"<<y.second<<endl;
    TH2D* hGridReferenceDetSpace = (TH2D*)histSaver->GetGridReferenceDetSpace()->Clone("hGridReferenceDetSpace_cl");
    profRebinned->GetYaxis()->SetRangeUser(y.first,y.second);
    profRebinned->GetXaxis()->SetRangeUser(x.first,x.second);
    hGridReferenceDetSpace->GetYaxis()->SetRangeUser(y.first,y.second);
    hGridReferenceDetSpace->GetXaxis()->SetRangeUser(x.first,x.second);
    hGridReferenceDetSpace->SetTitle(profRebinned->GetTitle());        //Set title to require
    hGridReferenceDetSpace->Draw("COL");
    profRebinned->GetZaxis()->SetTitleOffset(1.3);
    profRebinned->GetZaxis()->SetLabelOffset(0);
    profRebinned->SetContour(100);
    profRebinned->Draw("sameCOLZ");
    hGridReferenceDetSpace->Draw("sameCOL");
    settings->DrawMetallisationGrid(c1, 3);
    hGridReferenceDetSpace->Draw("sameCOL");
    c1->Update();
    name = "cProfRebinned_ZoomGoodCells"+appendix;
    c1->SetName(name);
    histSaver->DrawGoodCellsRegion(c1);
    histSaver->SaveCanvas(c1);//*/
    delete profRebinned;
    delete c1;
    delete hGridReferenceDetSpace;

    UInt_t xBins = hPulseHeightVsDetectorHitPostionXY->GetXaxis()->GetNbins();
    UInt_t yBins = hPulseHeightVsDetectorHitPostionXY->GetYaxis()->GetNbins();
    UInt_t xRebin = xBins/settings->getNColumns3d()/2;
    UInt_t yRebin = yBins/settings->getNRows3d()/2;

    name = "hPulseHeightVsDetectorHitPostion_Quarters";
    name.Append(appendix);
    TH2F* hDetXvsDetY3DMeanChargeQuarters = (TH2F*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin,yRebin,name);
    histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeQuarters);
    histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeQuarters);
    delete hDetXvsDetY3DMeanChargeQuarters;
    //
    name = "hPulseHeightVsDetectorHitPostion_Cells";
    name.Append(appendix);
    TH2D* hDetXvsDetY3DMeanChargeCells = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin*2,yRebin*2,name);
    histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeCells,hDetXvsDetY3DMeanChargeCells);
    histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeCells);
    delete hDetXvsDetY3DMeanChargeCells;
    //
//    delete hPulseHeightVsDetectorHitPostionXY;

}

void TAnalysisOf3dDiamonds::LongAnalysis_FillRelativeAddedTransparentCharge() {
    if(!settings->do3dTransparentAnalysis())
        return;
    pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
    UInt_t diamondPattern = settings->get3dWithHolesDiamondPattern();
    transparentCluster.SetTransparentClusterSize(transparentCluster.size());
    Float_t charge5 = transparentCluster.getPositiveCharge(useCMN);
    Float_t charge = 0;
    Int_t realHitChannel = (Int_t)(transparentCluster.GetTransparentHitPosition()+.5);
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);
    bool isInEdgeRegion =  settings->IsOnTheEdgeOfCell(relPos.first,relPos.second);
    for(UInt_t i = 0; i < transparentCluster.size(); i++){
        UInt_t clusterSize = i+1;
        transparentCluster.SetTransparentClusterSize(clusterSize);
        Float_t addedCharge = - charge;
        charge = transparentCluster.getPositiveCharge(useCMN);
        Int_t relativeHitChannel = transparentCluster.GetHighestSignalChannelTransparentCluster();
        addedCharge += charge;
        Float_t relCharge = charge/charge5;
        Float_t relAddedCharge = addedCharge/charge5;

        hTransparentAnalysisTransparentRelativeHitChannel[i]->Fill(relativeHitChannel-realHitChannel);;
        hTransparentAnalysisTransparentRelativeHitChannelProfile[i]->Fill(xPredDet,yPredDet,relativeHitChannel-realHitChannel);;

        hTransparentAnalysisTransparentCharge[i]->Fill(charge);
        hTransparentAnalysisTransparentAddedCharge[i]->Fill(addedCharge);
        hTransparentAnalysisRelativeAddedCharge[i]->Fill(relAddedCharge);
        hTransparentAnalysisRelativeCharge[i]->Fill(relCharge);

        hTransparentAnalysisTransparentChargeProfile[i]->Fill(xPredDet,yPredDet,charge);
        hTransparentAnalysisTransparentAddedChargeProfile[i]->Fill(xPredDet,yPredDet,addedCharge);
        hTransparentAnalysisRelativeAddedChargeProfile[i]->Fill(xPredDet,yPredDet,relAddedCharge);
        hTransparentAnalysisRelativeChargeProfile[i]->Fill(xPredDet,yPredDet,relCharge);

        if (settings->isBadCell(diamondPattern, cell.first)){
            hTransparentAnalysisTransparentChargeBadCells[i]->Fill(charge);
            hTransparentAnalysisTransparentAddedChargeBadCells[i]->Fill(addedCharge);
            hTransparentAnalysisRelativeAddedChargeBadCells[i]->Fill(relAddedCharge);
            hTransparentAnalysisRelativeChargeBadCells[i]->Fill(relCharge);
            hTransparentAnalysisTransparentRelativeHitChannelBadCells[i]->Fill(relativeHitChannel-realHitChannel);
        }

        if (settings->IsGoodCell(diamondPattern, cell.first)){

            hTransparentAnalysisTransparentChargeGoodCells[i]->Fill(charge);
            hTransparentAnalysisTransparentAddedChargeGoodCells[i]->Fill(addedCharge);
            hTransparentAnalysisRelativeAddedChargeGoodCells[i]->Fill(relAddedCharge);
            hTransparentAnalysisRelativeChargeGoodCells[i]->Fill(relCharge);
            hTransparentAnalysisTransparentRelativeHitChannelGoodCells[i]->Fill(relativeHitChannel-realHitChannel);
        }
        if (!isInEdgeRegion){
            hTransparentAnalysisTransparentChargeWithoutEdge[i]->Fill(charge);
            hTransparentAnalysisTransparentAddedChargeWithoutEdge[i]->Fill(addedCharge);
            hTransparentAnalysisRelativeAddedChargeWithoutEdge[i]->Fill(relAddedCharge);

            hTransparentAnalysisTransparentChargeProfileWithoutEdge[i]->Fill(xPredDet,yPredDet,charge);
            hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge[i]->Fill(xPredDet,yPredDet,addedCharge);
            hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge[i]->Fill(xPredDet,yPredDet,relAddedCharge);

            if (settings->isBadCell(diamondPattern, cell.first)){
                hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]->Fill(charge);
                hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge[i]->Fill(addedCharge);
                hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge[i]->Fill(relAddedCharge);

                if(i==0 && charge < 500)
                    hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells->Fill(xPredDet,yPredDet);
            }

            if (settings->IsGoodCell(diamondPattern, cell.first)){
                hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]->Fill(charge);
                hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge[i]->Fill(addedCharge);
                hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge[i]->Fill(relAddedCharge);
                if(i==0 && charge < 500)
                    hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells->Fill(xPredDet,yPredDet);
            }
        }
        //			cout<<"\t"<<i<<"\t"<<transparentCluster.getPositiveCharge(useCMN)<<"\t"<<transparentCluster.getPositiveCharge(useCMN)/charge1<<endl;
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_InitialiseRelativeAddedTransparentCharge() {

    for(UInt_t i = 0; i<maxClusterSize3d;i++){
        UInt_t clusterSize = i+1;
        Float_t min = 0;
        Float_t max = 50;
        Int_t bins  =512;
        TString name;

        min = -1.2;
        max = 1.2;
        name = TString::Format("hTransparentAnalysisRelativeAddedCharge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedCharge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeGoodCells_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeGoodCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeBadCells_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeBadCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeProfile_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeProfile.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        name = TString::Format("hTransparentAnalysisRelativeCharge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeCharge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeChargeGoodCells_%02dover05",clusterSize);
        hTransparentAnalysisRelativeChargeGoodCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeChargeBadCells_%02dover05",clusterSize);
        hTransparentAnalysisRelativeChargeBadCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeChargeProfile_%02dover05",clusterSize);
        hTransparentAnalysisRelativeChargeProfile.push_back(histSaver->GetProfile2dBinedInCells(name,2));


        name = TString::Format("hTransparentAnalysisTransparentCharge_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentCharge.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeGoodCells_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeGoodCells.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeBadCells_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeBadCells.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeProfile_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeProfile.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        min = -50;
        max = 1500;
        name = TString::Format("hTransparentAnalysisTransparentAddedCharge_%02dvs01",clusterSize);
        hTransparentAnalysisTransparentAddedCharge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeGoodCells_%02dvs01",clusterSize);
        hTransparentAnalysisTransparentAddedChargeGoodCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeBadCells_%02dvs01",clusterSize);
        hTransparentAnalysisTransparentAddedChargeBadCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeProfile_%02dvs01",clusterSize);
        hTransparentAnalysisTransparentAddedChargeProfile.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        //without edges
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeWithoutEdge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge_%02dover05",clusterSize);
        hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        name = TString::Format("hTransparentAnalysisTransparentChargeWithoutEdge_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeWithoutEdge.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeBadCellsWithoutEdge_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        name = TString::Format("hTransparentAnalysisTransparentChargeProfileWithoutEdge_clusterSize%02d",clusterSize);
        hTransparentAnalysisTransparentChargeProfileWithoutEdge.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        min = -50;
        max = 1500;
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeWithoutEdge_%02dvs%02d",clusterSize,clusterSize-1);
        hTransparentAnalysisTransparentAddedChargeWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge_%02dvs%02d",clusterSize,clusterSize-1);
        hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge_%02dvs%02d",clusterSize,clusterSize-1);
        hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge_%02dvs%02d",clusterSize,clusterSize-1);
        hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge.push_back(histSaver->GetProfile2dBinedInCells(name,2));

        min = -3.5;
        max = 3.5;
        bins = 7;
        name = TString::Format("hTransparentAnalysisTransparentRelativeHitPosition_%02d",clusterSize);
        hTransparentAnalysisTransparentRelativeHitChannel.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentRelativeHitPositionGoodCells_%02d",clusterSize);
        hTransparentAnalysisTransparentRelativeHitChannelGoodCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentRelativeHitPositionBadCells_%02d",clusterSize);
        hTransparentAnalysisTransparentRelativeHitChannelBadCells.push_back(new TH1F(name,name,bins,min,max));
        name = TString::Format("hTransparentAnalysisTransparentRelativeHitPositionProfile_%02d",clusterSize);
        hTransparentAnalysisTransparentRelativeHitChannelProfile.push_back(histSaver->GetProfile2dBinedInCells(name,2));
    }

    TString name = TString::Format("hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells");
    hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells = histSaver->GetHistoBinedInCells(name,8);
    hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells->SetTitle("BadCell PH_{1ch}<500ADC");
    name = TString::Format("hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells");
    hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells = histSaver->GetHistoBinedInCells(name,8);
    hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells->SetTitle("GoodCell PH_{1ch}<500ADC");
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveRelativeAddedTransparentCharge() {
    for(UInt_t i = 0; i < hTransparentAnalysisRelativeAddedCharge.size(); i++){
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedCharge[i]);
        delete hTransparentAnalysisRelativeAddedCharge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedChargeBadCells[i]);
        delete hTransparentAnalysisRelativeAddedChargeBadCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedChargeGoodCells[i]);
        delete hTransparentAnalysisRelativeAddedChargeGoodCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisRelativeAddedChargeProfile[i]);
        delete hTransparentAnalysisRelativeAddedChargeProfile[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeCharge[i]);
        delete hTransparentAnalysisRelativeCharge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeChargeBadCells[i]);
        delete hTransparentAnalysisRelativeChargeBadCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeChargeGoodCells[i]);
        delete hTransparentAnalysisRelativeChargeGoodCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisRelativeChargeProfile[i]);
        delete hTransparentAnalysisRelativeChargeProfile[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentCharge[i]);
        delete hTransparentAnalysisTransparentCharge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentChargeBadCells[i]);
        delete hTransparentAnalysisTransparentChargeBadCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentChargeGoodCells[i]);
        delete hTransparentAnalysisTransparentChargeGoodCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisTransparentChargeProfile[i]);
        delete hTransparentAnalysisTransparentChargeProfile[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedCharge[i]);
        delete hTransparentAnalysisTransparentAddedCharge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedChargeBadCells[i]);
        delete hTransparentAnalysisTransparentAddedChargeBadCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedChargeGoodCells[i]);
        delete hTransparentAnalysisTransparentAddedChargeGoodCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisTransparentAddedChargeProfile[i]);
        delete hTransparentAnalysisTransparentAddedChargeProfile[i];

        // without edges
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedChargeWithoutEdge[i]);
        delete hTransparentAnalysisRelativeAddedChargeWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge[i]);
        delete hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge[i]);
        delete hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge[i]);
        delete hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentChargeWithoutEdge[i]);
        delete hTransparentAnalysisTransparentChargeWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]);
        delete hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]);
        delete hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisTransparentChargeProfileWithoutEdge[i]);
        delete hTransparentAnalysisTransparentChargeProfileWithoutEdge[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedChargeWithoutEdge[i]);
        delete hTransparentAnalysisTransparentAddedChargeWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge[i]);
        delete hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge[i]);
        delete hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge[i]);
        delete hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge[i];

        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentRelativeHitChannel[i]);
        delete hTransparentAnalysisTransparentRelativeHitChannel[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentRelativeHitChannelGoodCells[i]);
        delete hTransparentAnalysisTransparentRelativeHitChannelGoodCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogram(hTransparentAnalysisTransparentRelativeHitChannelBadCells[i]);
        delete hTransparentAnalysisTransparentRelativeHitChannelBadCells[i];
        if(settings->do3dTransparentAnalysis())
            histSaver->SaveHistogramWithCellGrid(hTransparentAnalysisTransparentRelativeHitChannelProfile[i]);
        delete hTransparentAnalysisTransparentRelativeHitChannelProfile[i];

    }

    if(settings->do3dTransparentAnalysis()){
        cout<<"hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells: "<<hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells<<endl;
        histSaver->SaveHistogram(hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells,false,false);
    }
    delete hTransparentAnalysisTransparentSmallChargeInFirstChannelBadCells;


    if(settings->do3dTransparentAnalysis()){
        cout<<"hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells: "<<hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells<<endl;
        histSaver->SaveHistogram(hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells,false,false);
    }
    delete hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells;

    hTransparentAnalysisRelativeAddedCharge.clear();
    hTransparentAnalysisRelativeAddedChargeBadCells.clear();
    hTransparentAnalysisRelativeAddedChargeGoodCells.clear();
    hTransparentAnalysisRelativeAddedChargeProfile.clear();

    hTransparentAnalysisRelativeCharge.clear();
    hTransparentAnalysisRelativeChargeBadCells.clear();
    hTransparentAnalysisRelativeChargeGoodCells.clear();
    hTransparentAnalysisRelativeChargeProfile.clear();

    hTransparentAnalysisTransparentCharge.clear();
    hTransparentAnalysisRelativeAddedChargeBadCells.clear();
    hTransparentAnalysisTransparentChargeGoodCells.clear();
    hTransparentAnalysisTransparentChargeProfile.clear();

    hTransparentAnalysisTransparentAddedCharge.clear();
    hTransparentAnalysisTransparentAddedChargeBadCells.clear();
    hTransparentAnalysisTransparentAddedChargeGoodCells.clear();
    hTransparentAnalysisTransparentAddedChargeProfile.clear();

    // without edge

    hTransparentAnalysisRelativeAddedChargeWithoutEdge.clear();
    hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge.clear();
    hTransparentAnalysisRelativeAddedChargeGoodCellsWithoutEdge.clear();
    hTransparentAnalysisRelativeAddedChargeProfileWithoutEdge.clear();

    hTransparentAnalysisTransparentChargeWithoutEdge.clear();
    hTransparentAnalysisRelativeAddedChargeBadCellsWithoutEdge.clear();
    hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge.clear();
    hTransparentAnalysisTransparentChargeProfileWithoutEdge.clear();

    hTransparentAnalysisTransparentAddedChargeWithoutEdge.clear();
    hTransparentAnalysisTransparentAddedChargeBadCellsWithoutEdge.clear();
    hTransparentAnalysisTransparentAddedChargeGoodCellsWithoutEdge.clear();
    hTransparentAnalysisTransparentAddedChargeProfileWithoutEdge.clear();
}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateRelativeAddedTransparentChargeComparisonPlots(){

    for(int i=1; i<hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge.size(); i++){
        TString Histo1Title = TString::Format("RelativeAddedTransparentCharge_Cluster Size %i",i-1);
        hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i-1]->SetTitle(Histo1Title);
        TString name1 = TString::Format("RelativeAddedTransparentCharge_Landau_GoodCellsWithoutEges_ClusterSize_%d",i-1);
        hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i-1]->SetName(name1);
        TString Histo2Title = TString::Format("RelativeAddedTransparentCharge_Cluster Size %i",i);
        hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]->SetTitle(Histo2Title);
        TString name2 = TString::Format("RelativeAddedTransparentCharge_Landau_GoodCellsWithoutEges_ClusterSize_%d",i);
        hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]->SetName(name2);

        TString name = TString::Format("hTransparentAnalysisTransparentChargeGoodCellsWithoutEdgeComparison%i_to_%i",i-1,i);
        histSaver->SaveTwoHistos((string)name,hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i-1],hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]);
    }

    for(int i=1; i<hTransparentAnalysisTransparentChargeBadCellsWithoutEdge.size(); i++){
        TString hString = TString::Format("Cluster Size %i/%i",i,i);
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i-1]->SetTitle(hString);
        hString = TString::Format("TransparentChargeBadCells_ClusterSize_%i_%i",i,i);
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i-1]->SetName(hString);
        hString = TString::Format("Cluster Size %i/%i",i,i+1);
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]->SetTitle(hString);
        hString = TString::Format("TransparentChargeBadCells_ClusterSize_%i_%i",i,i+1);
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]->SetName(hString);

        TString name = TString::Format("hTransparentAnalysisTransparentChargeBadCellsWithoutEdgeComparison%i_to_%i",i-1,i);
        histSaver->SaveTwoHistos((string)name,hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i-1],hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]);
    }

    if(settings->do3dShortAnalysis()){
        TString name = TString::Format("hTransparentAnalysisTransparentChargeWithoutEdgeBadCellsComparison_DiamondPattern2_to_ClusterSize1");
        hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]->SetLineColor(kRed);
        histSaver->SaveTwoHistos((string)name,hLandau[1],hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]);

        name = TString::Format("hTransparentAnalysisTransparentChargeWithoutEdgeBadCellsComparison_DiamondPattern2_to_ClusterSize1_normalized");
        histSaver->SaveTwoHistosNormalized((string)name,hLandau[1],hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]);
    }

}

void TAnalysisOf3dDiamonds::LongAnalysis_InitGoodCellsLandaus() {
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
}

void TAnalysisOf3dDiamonds::LongAnalysis_FillGoodCellsLandaus(Float_t charge) {
    if(settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet)!=3)
        return;
    mapPredictedPositionsGoodCells[nEvent] = make_pair(xPredDet,yPredDet);
    if(validClusteredAnalysis)
        mapClusteredAnalysisGoodCells[nEvent] = clusteredCluster;
    if (validTransparentAnalysis)
        mapTransparentAnalysisGoodCells[nEvent] = transparentCluster;

    hPulseHeightVsDetectorHitPostionXYGoodCells->Fill(xPredDet,yPredDet,charge);
    hLandauGoodCells->Fill(charge);
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);
    bool isInEdgeRegion =  settings->IsOnTheEdgeOfCell(relPos.first,relPos.second);
    bool isInColumnRegion = settings->IsWithInTheColumnRadius(relPos.first,relPos.second);
    if (!isInEdgeRegion)
        hLandauGoodCellsWithoutEdges->Fill(charge);
    if(!isInColumnRegion)
        hLandauGoodCellsWithoutColumns->Fill(charge);
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveGoodCellsLandaus() {
    cout<<"[ TAnalysisOf3dDiamonds::LongAnalysis_SaveGoodCellsLandaus]"<<endl;
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
        hLandauGoodCellsWithoutColumns->Draw("same");
        hLandauGoodCellsWithoutEdges->Draw("same");
        histSaver->SaveCanvas(c1);
        delete c1;
    }

    if(hLandauGoodCellsWithoutEdges) delete hLandauGoodCellsWithoutEdges;
    if(hLandauGoodCellsWithoutColumns) delete hLandauGoodCellsWithoutColumns;
    if(hLandauGoodCells) delete hLandauGoodCells;
}

void TAnalysisOf3dDiamonds::DoMonteCarloOfAvrgChargePerBinInOverlay(
        TProfile2D* profOverlay, TH1F* hLandauOfOverlay) {
    if(!profOverlay || !hLandauOfOverlay){
        cerr<<"[TAnalysisOf3dDiamonds::DoMonteCarloOfAvrgChargePerBinInOverlay] one of the histograms is not defined! ";
        cout<<profOverlay<<" "<<hLandauOfOverlay<<endl;
    }
    TH2D* hBinContents = histSaver->GetBinContentHisto(profOverlay);
    histSaver->SaveHistogram(hBinContents,false,false,(TString)"colzText");
    TString name = profOverlay->GetName()+(TString)"_Entries";
    Int_t max = hBinContents->GetBinContent(hBinContents->GetMaximumBin())+1;
    Int_t min = hBinContents->GetBinContent(hBinContents->GetMinimumBin())+1;
    max -=5;
    min -=5;
    Float_t cut = 500;
    Int_t bins = max-min;
    TH1F* hEntries = new TH1F(name,name,bins,min,max);
    for(Int_t binx = 1;binx < hBinContents->GetNbinsX();binx++)
        for(Int_t biny = 1;biny < hBinContents->GetNbinsY();biny++){
            Int_t bin = hBinContents->GetBin(binx,biny);
            hEntries->Fill(hBinContents->GetBinContent(bin));
        }
    histSaver->SaveHistogram(hEntries);
    name = "hMonteCarloAvrgChargePerBin";
    TH1D* hMonteCarloAvrgChargePerBin = new TH1D(name,name,PulseHeightBins,PulseHeightMinMeanCharge-100,PulseHeightMaxMeanCharge+100);

    name =TString::Format("hMonteCarloNBelowCut_%d",(int)cut);;
    TH1D* hNumberOfEntriesBelowCut = new TH1D(name,name,100,0,100);
    hNumberOfEntriesBelowCut->GetXaxis()->SetTitle(TString::Format("number of  entries below %.f",cut));
    hNumberOfEntriesBelowCut->GetYaxis()->SetTitle("number of  entries #");

    name =TString::Format("hMonteCarloRelativeBelowCut_%d",(int)cut);;
    TH1D* hRelativeNumberBelowCut = new TH1D(name,name,100,0,1);
    hRelativeNumberBelowCut->GetXaxis()->SetTitle(TString::Format("rel. number of  entries below %.f",cut));
    hRelativeNumberBelowCut->GetYaxis()->SetTitle("number of  entries #");
    Float_t mean =0;
    Int_t entries = 0;
    int nBelowCut = 0;
    for(UInt_t nMC = 0; nMC < 1e6; nMC++){
        entries = hEntries->GetRandom();
        mean = 0;
        nBelowCut = 0;
        for(Int_t entry =0; entry < entries; entry++){
            Float_t charge =  hLandauOfOverlay->GetRandom();
            mean += charge;
            if(charge<cut)
                nBelowCut ++;
        }
        mean /= (Float_t)entries;
        hNumberOfEntriesBelowCut->Fill(nBelowCut);
        hRelativeNumberBelowCut->Fill((Float_t)nBelowCut/(Float_t)entries);
        //        cout<<TString::Format("%6d --> %6.1f (%d)",nMC,mean,entries)<<endl;;
        hMonteCarloAvrgChargePerBin->Fill(mean);
    }
    hMonteCarloAvrgChargePerBin->GetXaxis()->SetTitle("mean charge per bin_{MonteCarlo}");
    hMonteCarloAvrgChargePerBin->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hMonteCarloAvrgChargePerBin,false,true,true);
    histSaver->SaveHistogram(hNumberOfEntriesBelowCut,false,true,true);
    histSaver->SaveHistogram(hRelativeNumberBelowCut,false,true,true);
    delete hMonteCarloAvrgChargePerBin;
    delete hBinContents;
    delete hEntries;

}


void TAnalysisOf3dDiamonds::LongAnalysis_CompareTransparentAndClusteredAnalysis_Maps(){
    if ( !settings->do3dTransparentAnalysis())
        return;
    TAnalysisOfAnalysisDifferences* differences = new TAnalysisOfAnalysisDifferences(settings,histSaver,"good");
    differences->setStripHistogram(this->hLandauStrip);
    differences->set3DPhantomLandau(this->hLandau3DPhantom);
    differences->setClusteredMap(&mapClusteredAnalysisGoodCells);
    differences->setTransparentMap(&mapTransparentAnalysisGoodCells);
    differences->setPredictedPositions(&mapPredictedPositionsGoodCells);
    differences->Analysis();
    delete differences;

    differences = new TAnalysisOfAnalysisDifferences(settings,histSaver,"all");
    differences->setStripHistogram(this->hLandauStrip);
    differences->set3DPhantomLandau(this->hLandau3DPhantom);
    differences->setClusteredMap(&mapClusteredAnalysisAllCells);
    differences->setTransparentMap(&mapTransparentAnalysisAllCells);
    differences->setPredictedPositions(&mapPredictedPositionsAllCells);
    differences->Analysis();
    delete differences;


    differences = new TAnalysisOfAnalysisDifferences(settings,histSaver,"allButBad");
    differences->setStripHistogram(this->hLandauStrip);
    differences->set3DPhantomLandau(this->hLandau3DPhantom);
    differences->setClusteredMap(&mapClusteredAnalysisAllButBadCells);
    differences->setTransparentMap(&mapTransparentAnalysisAllButBadCells);
    differences->setPredictedPositions(&mapPredictedPositionsAllButBadCells);
    differences->Analysis();
    delete differences;
}
