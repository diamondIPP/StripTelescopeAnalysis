/*
 * TAnalysisOf3dDiamonds.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair,iain
 */
//Plots to not be plotted by transparent Analysis:
/*c2ClusterAnalysis_ClusterPatterns_trans.png
cRelativeChargeTwoClustersXY_trans.png
cShortAnalysis2TotalChargeXY_trans.png
cTotalAvrgChargeXY_trans.png ???
c_pfx.png ???
h2ClusterAnalysis_ClusterPatterns_all.png
h2ClusterAnalysis_ClusterPatterns_pattern1.png
h2ClusterAnalysis_ClusterPatterns_pattern2.png
h2ClusterAnalysis_ClusterPatterns_pattern3.png*/

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

    //htmlLandau=new THTMLLandaus(settings);

    settings->goTo3dDiamondTreeDir();
    cout<<"Selection Tree file path is: "<<settings->getSelectionTreeFilePath()<<endl;
    cout<<"Alignment file path is: "<<settings->getAlignmentFilePath()<<endl;
    cout<<"Eta distribution file path is: "<<settings->getEtaDistributionPath()<<endl;
    eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
    histSaver=new HistogrammSaver(settings);
    settings->goTo3dDiamondAnalysisDir();

    histSaver->SetPlotsPath(settings->get3dDiamondAnalysisPath());
    histSaver->SetRunNumber(runNumber);
    //htmlLandau->setFileGeneratingPath(settings->getSelectionAnalysisPath());//todo Write html3dDiamond
    settings->goTo3dDiamondTreeDir();
    //initialiseHistos();
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

    maxClusterSize3d = 5;
}

TAnalysisOf3dDiamonds::~TAnalysisOf3dDiamonds() {
    //htmlLandau->generateHTMLFile();
    if(eventReader!=0) delete eventReader;
    if(histSaver!=0)   delete histSaver;
    //if(htmlLandau!=0)  delete htmlLandau;
    settings->goToOutputDir();
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
    histSaver->SetNumberOfEvents(nEvents);
    settings->diamondPattern.Print();
    for(nEvent=0;nEvent<nEvents;nEvent++){
        TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
        eventReader->LoadEvent(nEvent);
        if(!eventValid()){
            //			cout<<"don't use"<<endl;
            continue;
        }
        //		cout<<"use Event"<<endl;
        // Analyse

        //if(settings->do3dTransparentAnalysis() == 1){TransparentAnalysis();}
        if(settings->do3dTransparentAnalysis()){isTransparentCluster = TransparentAnalysis();}
        else
            isTransparentCluster = false;
        StripAnalysis();
        if(settings->do3dShortAnalysis() == 1){ShortAnalysis();}
        if(settings->do3dLongAnalysis() == 1){LongAnalysis();}

    }
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

    saveHistos();
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
        //		cout<<nEvent<<" invalid Track"<<endl;
        return false;
    }
    fiducialValueX = eventReader->getFiducialValueX();
    fiducialValueY = eventReader->getFiducialValueY();
    if (!settings->isInRoughFiducialCut(fiducialValueX,fiducialValueY)){
        //		cout<<nEvent<<" not in rough fiducial cut: "<<fiducialValueX<<"/"<<fiducialValueY<<endl;
        return false;
    }


    if(predictedPosition) delete predictedPosition;
    predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
    if (!predictedPosition)
        return false;
    chi2x = predictedPosition->getChi2X();
    chi2y = predictedPosition->getChi2Y();
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

    TCluster diamondCluster;
    if(!settings->do3dTransparentAnalysis()){
        if(!eventReader->getNDiamondClusters())
            return;
        diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);

    }
    else{
        if(!isTransparentCluster)
            return;
        diamondCluster = transparentCluster;
    }

    if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
        return;
    int stripDetector = 1;
    if(settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY)!=stripDetector)
        return;

    if(eventReader->getNDiamondClusters()!=1)
        return;
    //	UInt_t diamond = TPlaneProperties::getDetDiamond();
    //TCluster diamondCluster = eventReader->getCluster(diamond,0);
    int areaStripDetector = 0;
    if (!settings->isClusterInDiaDetectorArea(diamondCluster,areaStripDetector) ){
        return;
    }
    if( !settings->diamondPattern.isValidCluster(diamondCluster)){
        return;
    }
    if (diamondCluster.isSaturatedCluster())
        return;

    //cout<<"diamondCluster charge: "<<diamondCluster.getCharge()<<endl;


    hLandauStrip->Fill(diamondCluster.getCharge());
    hLandauStripFidCutXvsFidCutY->Fill(fiducialValueX, fiducialValueY);
}

void TAnalysisOf3dDiamonds::ShortAnalysis() {

    if(!settings->do3dTransparentAnalysis()){
        Float_t maxChi2 = settings->getChi2Cut3D();
        if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
            return;

        hNumberofClusters->Fill(eventReader->getNDiamondClusters());
        ClusterPlots(eventReader->getNDiamondClusters(),fiducialValueX,fiducialValueY);

        Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);

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
    }
    else
        ShortAnalysis_Analyse1Cluster();

    //	if(eventReader->getNDiamondClusters()==1){
    //		TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
    //		if(diamondCluster.isSaturatedCluster())
    //			return;
    //
    //		HitandSeedCount(&diamondCluster);
    //		Int_t clusterSize = diamondCluster.size()-2;
    //		vecClusterSize.push_back(clusterSize);
    //
    //		//hFidCutXvsFidCutYvsCharge.at(1)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
    //		//RemoveLumpyClusters(&diamondCluster);
    //
    //		//Edge Finding
    //
    //		Float_t clusterCharge = diamondCluster.getCharge(false);
    ////		this->fillEdgeDistributionmakes(clusterCharge);
    //
    //		//Universal PHvsChannel Plot
    //		for(UInt_t i=0; i < settings->diamondPattern.getNIntervals();i++){
    //
    //			pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
    //			//cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;
    //
    //			if((Int_t)diamondCluster.getHighestSignalChannel()<=channels.second&&(Int_t)diamondCluster.getHighestSignalChannel()>=channels.first){
    //
    //				hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
    //				hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);
    //				//if(!settings->getSelectionFidCuts()->getFidCut(i+1)->isInFiducialCut(fiducialValueX,fiducialValueY))
    //				//return;
    //
    //				if(RemoveEdgeHits(&diamondCluster,channels)==1)		//If cluster has hit in edge channel remove.
    //					return;
    //
    //				/*if(HitCount>0||SeedCount>1)
    //					return;
    //				 */
    //
    //				hEventsvsChannel[i]->Fill(diamondCluster.getHighestSignalChannel());
    //				hPHvsChannel[i]->Fill(diamondCluster.getCharge(false),diamondCluster.getHighestSignalChannel());
    //				//hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);
    //				//hPHvsPredictedXPos.at(i)->Fill(diamondCluster.getCharge(false),xPredicted);
    //				hLandau[i]->Fill(diamondCluster.getCharge(false));
    //				vecPHDiamondHit[i]->push_back(diamondCluster.getCharge(false));
    //				//vecXPredicted.at(i)->push_back(xPredicted);
    //				//vecYPredicted.at(i)->push_back(yPredicted);
    //				hChi2XChi2Y[i]->Fill(chi2x, chi2y);
    //				hFidCutXvsFidCutY[i]->Fill(fiducialValueX,fiducialValueY);
    //				//For hFidCutXvsFidCutYvsMeanCharge
    //				//hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);
    //
    //			}
    //		}		//End of for diamond patterns
    //	}		//End of if clusters = 1
}

void TAnalysisOf3dDiamonds::ShortAnalysis_Analyse1Cluster(UInt_t clusterNo){

    TCluster diamondCluster;
    if(!settings->do3dTransparentAnalysis()){
        diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),clusterNo);
        if(diamondCluster.isSaturatedCluster())
            return;
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            //		cerr <<" Cluster is invalid: ";
            //		diamondCluster.Print(1);
            return;
        }
        HitandSeedCount(&diamondCluster);
        Int_t clusterSize = diamondCluster.size()-2;
        vecClusterSize.push_back(clusterSize);
        //Edge Finding
        ShortAnalysis_FillEdgeDistributions(diamondCluster.getCharge(false));
    }
    else{
        if(!isTransparentCluster)
            return;
        diamondCluster = transparentCluster;
    }

    Float_t clusterCharge = diamondCluster.getCharge(false);
    //	cout<<"ClusterCharge: "<<clusterCharge<<endl;

    ShortAnalysis_FillMeanChargeVector(clusterCharge);
    hTotalAvrgChargeXY->Fill(xPredDet,yPredDet,clusterCharge);
    //Universal PHvsChannel Plot

    for(UInt_t i=0; i < settings->diamondPattern.getNIntervals();i++){

        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        //cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;
        Int_t HighestSignalChannel = diamondCluster.getHighestSignalChannel();

        if(HighestSignalChannel<=channels.second && HighestSignalChannel>=channels.first){

            if(!settings->do3dTransparentAnalysis()){
                /*if(!settings->getSelectionFidCuts()->getFidCut(i+1)->IsInFiducialCut(fiducialValueX,fiducialValueY))
								return;*/
                if(RemoveEdgeHits(&diamondCluster,channels)==1)		//If cluster has hit in edge channel remove.
                    return;
            }
            if(!settings->getSelectionFidCuts()->getFidCut(i+1)->IsInFiducialCut(fiducialValueX,fiducialValueY))
                return;

            //hTransparentAnalysisValidClusterFidCutXvsFidCutY->Fill(fiducialValueX, fiducialValueY);
            hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
            hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);

            /*if(HitCount>0||SeedCount>1)
						return;
             */

            hEventsvsChannel[i]->Fill(diamondCluster.getHighestSignalChannel());
            hPHvsChannel[i]->Fill(diamondCluster.getCharge(false),diamondCluster.getHighestSignalChannel());
            //hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);
            //hPHvsPredictedXPos.at(i)->Fill(diamondCluster.getCharge(false),xPredicted);
            hLandau[i]->Fill(diamondCluster.getCharge(false));
            vecPHDiamondHit[i]->push_back(diamondCluster.getCharge(false));
            //vecXPredicted.at(i)->push_back(xPredicted);
            //vecYPredicted.at(i)->push_back(yPredicted);
            hFidCutXvsFidCutY[i]->Fill(fiducialValueX,fiducialValueY);
            //For hFidCutXvsFidCutYvsMeanCharge
            //hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);

            if(!settings->do3dTransparentAnalysis()){
                hHitandSeedCount[i]->Fill(HitCount,SeedCount);
                hChi2XChi2Y[i]->Fill(chi2x, chi2y);
            }
        }
    }		//End of for diamond patterns
}

void TAnalysisOf3dDiamonds::ShortAnalysis_Analyse2Cluster(){
    TCluster cluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
    TCluster cluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
    Float_t pos1 = cluster1.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t pos2 = cluster2.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t ph1 = cluster1.getCharge();
    Float_t ph2 = cluster2.getCharge();
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

    //	cout<<TString::Format("%6d %2d/%2d %4.1f-->%2d %4.1f-->%2d\t\t%6.1f/%6.1f",
    //			nEvent,predictedArea,predictedDetector,pos1,pattern1+1,pos2,pattern2+1,ph1,ph2);
    //	if(pos2<85&&pos2>84)
    //		cluster2.Print(1);
    if (!settings->diamondPattern.isValidChannelPosition(pos1)){
        //		cout<<"\tAnalyze 1"<<endl;
        ShortAnalysis_Analyse1Cluster(1);//todo how can we use them for long analysis
        return;
    }
    if (!settings->diamondPattern.isValidChannelPosition(pos2)){
        //		cout<<"\tAnalyze 2"<<endl;
        ShortAnalysis_Analyse1Cluster(0);//todo how can we use them for long analysis
        return;
    }
    if (pos1>128||pos1<0||pos2>128||pos2<0){
        //		cout<<"\tInvalid"<<endl;
        //		cout<<"\n"<<nEvent<<" "<<pos1<<": "<<ph1<<"\t"<<pos2<<": "<<ph2<<endl;
        return;
    }

    Int_t delta1 = pattern1 - predictedDetector;
    Int_t delta2 = pattern2 - predictedDetector;
    hShortAnalysis2ClusterHitPattern_1stCluster->Fill(predictedDetector,delta1);
    hShortAnalysis2ClusterHitPattern_2ndCluster->Fill(predictedDetector,delta2);

    vecPH_Cluster1_ShortAna.push_back(ph1);
    vecPH_Cluster2_ShortAna.push_back(ph2);
    vecCh_Cluster1_ShortAna.push_back(pos1);
    vecCh_Cluster2_ShortAna.push_back(pos2);
    if (predictedDetector == pattern1||predictedDetector==pattern2){
        //		cout<<"\t1 "<<TString::Format("%1.3f",ph2/ph1);
        if(predictedDetector == 2|| true){
            Double_t relCharge = ph2/ph1;
            hRelativeChargeTwoClustersX->Fill(xPredDet,relCharge);
            hRelativeChargeTwoClustersY->Fill(yPredDet,relCharge);
            hRelativeChargeTwoClustersXY->Fill(xPredDet,yPredDet,relCharge);
            hShortAnalysis2TotalChargeXY->Fill(xPredDet,yPredDet,ph1+ph2);
            //			cout<<"Fill "<<xPredDet<<"/"<<yPredDet<<" "<<relCharge<<" ";
            //			cout<<" X";
        }
    }
    if (predictedDetector == pattern2){
        //		cout<<"\t2 "<<TString::Format("%1.3f",ph1/ph2);
    }
    hTotalAvrgChargeXY->Fill(xPredDet,yPredDet,ph1+ph2);
    //		cout<<endl;
}

void TAnalysisOf3dDiamonds::LongAnalysis() {

    TCluster diamondCluster;
    if(!settings->do3dTransparentAnalysis()){
        Float_t maxChi2 = settings->getChi2Cut3D();
        if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
            return;
        if(eventReader->getNDiamondClusters()!=1)
            return;
        diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
    }
    else{
        if(!isTransparentCluster)
            return;
        Int_t DiamondPattern;
        DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);

        if(DiamondPattern !=3)
            return;
        diamondCluster = transparentCluster;
    }

    //	cout<<"xPredDet, yPredDet"<<xPredDet<<", "<<yPredDet<<endl;
    pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
    //	cout<<cell.first<<" "<<cell.second<<endl;
    UInt_t cellNo = cell.first;
    //	cout<<"cellNo: "<<cellNo<<endl;
    UInt_t quarterNo = cell.second;
    //	if(cellNo>=0)cout<<"cell: "<<cellNo<<"--> row "<<row<<", column "<<column<<endl;

    if(!settings->do3dTransparentAnalysis()){
        Int_t area3DwithColumns = 2;
        if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
            hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
            return;
        }
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
            return;
        }
    }
    Int_t area3DwithColumns = 2;

    /*if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
		hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
		return;
	}*/
    /*if( !settings->diamondPattern.isValidCluster(diamondCluster)){
		hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
		return;
	}*/

    /*if(!settings->isValidCellNo(cellNo)){
		hLongAnalysisInvalidCellNo->Fill(xPredDet,yPredDet);
		return;
	}*/


    /*if(diamondCluster.isSaturatedCluster())
		return;*/

    //	cout<<cellNo<<" "<< flush;
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);

    Float_t charge = diamondCluster.getCharge(false);
    hPulseHeightVsDetectorHitPostionXY->Fill(xPredDet,yPredDet,charge);
    //	hDetXvsDetY3DEvents->Fill(xPredDet,yPredDet,1);
    if(settings->IsGoodCell(3,cellNo))
        hPulseHeightVsDetectorHitPostionXYGoodCells->Fill(xPredDet,yPredDet,charge);

    if(!settings->isBadCell(3,cellNo)){
        MakeGhostCluster(&diamondCluster,5);
        for(int i= 0; i< hLandauMinusBadCellsOffsetAnalysis.size();i++){
            GhostCluster.SetTransparentClusterSize(i+1);
            hLandauMinusBadCellsOffsetAnalysis.at(i)->Fill(GhostCluster.getCharge(false));
        }}

    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
        int badcell = settings->getDeadCell3D()[i];
        int relCellY = cellNo - badcell;
        //		cout<<i<<" "<<badcell<<" "<<cellNo<<endl;
        //		cout<<"sqrt(relCellY*relCellY): "<<sqrt(relCellY*relCellY)<<endl;
        if(relCellY ==0 || sqrt(relCellY*relCellY) <= 1){
            float relY = relPos.second+settings->GetCellHeight() + settings->GetCellHeight()*relCellY;
            if(verbosity) cout<<"DeadCellAnalysis: relCellY: "<<relCellY<<" relPos.second: "<<relPos.second<<" relY: "<<relY<<endl;
            hDeadCellCharge[i]->Fill(relY, charge);
            hDeadCellPositions[i]->Fill(xPredDet,yPredDet);
        }
    }
    if(cellNo < hCellsLandau.size())
        hCellsLandau.at(cellNo)->Fill(charge);
    //	cout<<cellNo<<"-"<<quarterNo<<endl;
    if(!settings->do3dTransparentAnalysis()){
        if(cellNo < hCellsClusteSize.size()){
            Int_t ClusterSize = diamondCluster.getClusterSize()-2;
            Int_t MaxHistoClusterSize = hCellsClusteSize.at(0)->GetNbinsX();
            if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
                hCellsClusteSize.at(cellNo)->Fill(MaxHistoClusterSize);
            else
                hCellsClusteSize.at(cellNo)->Fill(ClusterSize);
        }
    }
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
                Int_t ClusterSize = diamondCluster.size()-2;
                Int_t MaxHistoClusterSize = hQuarterCellsClusterSize[0].at(0)->GetNbinsX();
                if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(MaxHistoClusterSize);
                else
                    hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(ClusterSize);
            }
        }
    }
    else{
        if(verbosity>6) cout <<TString::Format("X%2d-%d,\t%3.1f",cellNo,quarterNo,charge)<<endl;
    }
    //cout<<"Transparent Cluster Channels: "<<diamondCluster.getFirstHitChannel()<<" - "<<diamondCluster.getLastHitChannel()<<endl;
    Int_t MaxOverlayClusterSize = 3;
    for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
    	diamondCluster.SetTransparentClusterSize(ClusterSize);
    	Float_t charge = diamondCluster.getCharge(false);
    	LongAnalysis_FillOverlayedHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize);
    	LongAnalysis_FillOverlayCentralColumnHistos(cellNo,relPos.first,relPos.second,charge, ClusterSize, diamondCluster);
    	LongAnalysis_FillOverlayOffsetHistos(cellNo,relPos.first,relPos.second,charge,ClusterSize);
        LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(cellNo,relPos.first,relPos.second,charge, ClusterSize, &diamondCluster);

    }
    LongAnalysis_FillEdgeFreeHistos(xPredDet,yPredDet,charge);
    LongAnalysis_FillRelativeAddedTransparentCharge();
};

bool TAnalysisOf3dDiamonds::TransparentAnalysis() {

    if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
        return false;

    Int_t DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
    //cout<<"DiamondPattern: "<<DiamondPattern<<endl;
    bool bOutput = (verbosity >5 &&xPredDet >= 3660 && xPredDet <= 3710&& yPredDet >= 0 && yPredDet <= 1500);
    if (bOutput)	cout<<TString::Format("%6d, Hit in %5.1f,%5.1f ",nEvent, xPredDet,yPredDet);
    if(DiamondPattern !=1 && DiamondPattern !=2 && DiamondPattern !=3){
        //cout<<TString::Format("xPredDet: %.2f, yPredDet: %.2f", xPredDet, yPredDet)<<endl;
        hTransparentAnalysisInvalidCluster->Fill(xPredDet, yPredDet);
        if (bOutput) cout<<" invalid Cluster: "<<DiamondPattern<<endl;
        return false;
    }

    /*if(DiamondPattern !=3)
		return;*/

    Int_t XdetChannelSpaceInt =  settings->diamondPattern.convertMetricToIntChannel(xPredDet);
    Float_t XdetChannelSpace = settings->diamondPattern.convertMetricToChannel(xPredDet);


    if(XdetChannelSpaceInt<0){	//Returns -9998 when hit in channel 93.
        if (bOutput) cout<<" invalid XdetChannelInSpaceInt"<<endl;
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
        if (bOutput) cout<<"SaturatedCluster"<<endl;
        return false;
    }

    UInt_t  clusSize = transparentCluster.GetTransparentClusterSize();
    //	pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
    //	transparentCluster.SetTransparentClusterSize(1);
    //	Float_t charge1 = transparentCluster.getCharge();
    //	if(cell.first == 36 || cell.first == 37 || cell.first == 47 || cell.first == 48){
    //		//		cout<<"\nnEvent: "<<nEvent<<" "<<cell.first<<"/"<<cell.second<<"\n";
    //		for(UInt_t i = 0; i < transparentCluster.size(); i++){
    //			transparentCluster.SetTransparentClusterSize(i+1);
    //			//			cout<<"\t"<<i<<"\t"<<transparentCluster.getCharge()<<"\t"<<transparentCluster.getCharge()/charge1<<endl;
    //		}
    //	}
    //	transparentCluster.SetTransparentClusterSize(transparentCluster.size());
    //	transparentCluster.SetTransparentClusterSize(transparentCluster.size()+1);
    //cout<<"transparentCluster.getCharge(): "<<transparentCluster.getCharge()<<endl;
    //	transparentCluster.Print();
    //cout<<"diamondCluster.getHighestSignalChannel(): "<<transparentCluster.getHighestSignalChannel()<<endl;
    if (bOutput) cout<<"Valid Cluster"<<endl;
    return true;

}

void TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    cout<<"[TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos()]"<<endl;
    //Universal histograms
    TString name = "hRelativeChargeTwoClustersX";
    name.Append(appendix);
    Float_t xmax = settings->get3dMetallisationFidCuts()->getXHigh();
    hRelativeChargeTwoClustersX = new TH2F(name,name,1024,0,xmax,100,0,20);
    hRelativeChargeTwoClustersX->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersX->GetYaxis()->SetTitle("rel. ph: ph_{clus1}/ph_{clus2}");

    Float_t ymax = settings->get3dMetallisationFidCuts()->getYHigh();
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
    hTotalAvrgChargeXY = new TProfile2D(name,name, 128,0,xmax,128,0,ymax);
    hTotalAvrgChargeXY->GetXaxis()->SetTitle("pred. 	det hit position X /#mum");
    hTotalAvrgChargeXY->GetYaxis()->SetTitle("pred. det hit position Y /#mum");
    hTotalAvrgChargeXY->GetZaxis()->SetTitle("total charge: ph_{total} = #sum_{i}{ph_{i}}/ADC");

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
        cout<<"Loop: "<<i<<endl;
        pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
        //hLandau
        name = TString::Format("hLandau_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        if(verbosity>1) cout<<"Create "<<name<<endl;
        hLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        if(hLandau.back()){
            hLandau.back()->GetXaxis()->SetTitle("PH of diamond cluster");
            hLandau.back()->GetYaxis()->SetTitle("number of entries #");
        }
        else
            cerr<<"hLandau:'"<<name<<"' wasn't created correctly"<<endl;

        //hPHvsChannel
        name = TString::Format("hEventsvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hEventsvsChannel.push_back(new TH1F(name,name,100,0,100));
        if(hEventsvsChannel.back()){
            hEventsvsChannel.back()->GetXaxis()->SetTitle("HighestPH [ch]");
            hEventsvsChannel.back()->GetYaxis()->SetTitle("No. Events");
        }

        //hPHvsChannel
        name = TString::Format("hPHvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hPHvsChannel.push_back(new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,100,0,100));
        if(hPHvsChannel.back()){
            hPHvsChannel.back()->Draw();
            hPHvsChannel.back()->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->GetXaxis()->SetLimits(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->SetAxisRange(PulseHeightMin, PulseHeightMax, "X");
            hPHvsChannel.back()->GetYaxis()->SetRangeUser(channels.first-1,channels.second+1);
            hPHvsChannel.back()->GetXaxis()->SetTitle("cluster charge /ADC");
            hPHvsChannel.back()->GetYaxis()->SetTitle("XPos /ch");
            //hPHvsChannel.back()->GetXaxis()->SetRange(PulseHeightMin,PulseHeightMax);
        }

        //hPHvsChannel.at(i)->SetMaximum(3000);
        //hPHvsChannel.at(i)->SetMinimum(0);

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
        cout<<"("<<xLow<<"-"<<xHigh<<", "<<yLow<<"-"<<yHigh<<")"<<endl;
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

    //hXEdgeCharge
    //xEdge
    //	stringstream hEdgeChargeName; hEdgeChargeName<<"hEdgeCharge"<<FileNameEnd;
    //	hEdgeCharge = new TH1F(hEdgeChargeName.str().c_str(),hEdgeChargeName.str().c_str(),250,3100,4100);
    //	hEdgeCharge->SetTitle(hEdgeChargeName.str().c_str());
    //	hEdgeCharge->GetXaxis()->SetTitle("X Diamond (um)");
    //	hEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");
    //
    //	hEdgeChargeEvents = (TH1F*)hEdgeCharge->Clone("hEdgeChargeEvents");
    //	hEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

    //	hEdgeMeanCharge = (TH1F*)hEdgeCharge->Clone("hEdgeMeanCharge");
    //	hEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

    //	//hYEdgeCharge
    //	stringstream hyEdgeChargeName; hyEdgeChargeName<<"hyEdgeCharge"<<FileNameEnd;
    //	hyEdgeCharge = new TH1F(hyEdgeChargeName.str().c_str(),hyEdgeChargeName.str().c_str(),2500,0000,6000);
    //	hyEdgeCharge->SetTitle(hyEdgeChargeName.str().c_str());
    //	hyEdgeCharge->GetXaxis()->SetTitle("Y Diamond (um)");
    //	hyEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

    //	hyEdgeChargeEvents = (TH1F*)hyEdgeCharge->Clone("hyEdgeChargeEvents");
    //	hyEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

    //	hyEdgeMeanCharge = (TH1F*)hyEdgeCharge->Clone("hyEdgeMeanCharge");
    //	hyEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

    name ="hShortAnalysis2ClusterHitPattern_1stCluster";
    name.Append(appendix);
    int nFidCuts =settings->get3dMetallisationFidCuts()->getNFidCuts();
    hShortAnalysis2ClusterHitPattern_1stCluster =  new TH2F(name,name,
            nFidCuts,.5,nFidCuts+.5,
            5,-nFidCuts+.5,nFidCuts-.5);
    hShortAnalysis2ClusterHitPattern_1stCluster->GetXaxis()->SetTitle("predicteded hit pattern");
    hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()->SetTitle("cluster_{1}-hit-pattern - predicted-hit-pattern");

    name = "hShortAnalysis2ClusterHitPattern_2ndCluster";
    name.Append(appendix);
    hShortAnalysis2ClusterHitPattern_2ndCluster = (TH2F*)hShortAnalysis2ClusterHitPattern_1stCluster->Clone(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->SetTitle(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->GetYaxis()->SetTitle("cluster_{2}-hit-pattern - predicted-hit-pattern");

}

void TAnalysisOf3dDiamonds::initialise3DGridReference() {

    stringstream hGridReferenceName; hGridReferenceName<<""<<FileNameEnd;
    TString nameDet = "hGridRefenrenceDetSpace";
    TString nameCell = "hGridRefenrenceCellSpace";
    Float_t xBins = settings->getNColumns3d();
    TFidCutRegions* mettalisationFidCuts = settings->get3dMetallisationFidCuts();
    TFiducialCut* fidCut = mettalisationFidCuts->getFidCut(3);
    Float_t xLow = fidCut->GetXLow();//getXMetalisationStart3d;
    Float_t xHigh = fidCut->GetXHigh();//getXMetalisationEnd3d;
    Float_t yBins = settings->getNRows3d();
    Float_t yLow = fidCut->GetYLow();
    Float_t yHigh = fidCut->GetYHigh();//getYMetalisationEnd3d;
    //	cout<<"nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh"<<endl;
    //	cout<<nameDet<<" "<<nameDet<<" "<<xBins<<" "<<xLow<<" "<<xHigh<<" "<<yBins<<" "<<yLow<<" "<<yHigh<<endl;
    hGridReferenceDetSpace = new TH2D(nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh);
    hGridReferenceCellSpace = new TH2D(nameCell,nameDet,xBins,0,xBins,yBins,0,yBins);

    for(int i=0;i<settings->getNRows3d();i++){
        hGridReferenceDetSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
        hGridReferenceCellSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
    }
    for(UInt_t j=0;j<settings->getNRows3d();j++){
        hGridReferenceDetSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
        hGridReferenceCellSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
    }
    hGridReferenceDetSpace->SetStats(kFALSE);
    hGridReferenceDetSpace->SetTickLength(0.0, "X");
    hGridReferenceDetSpace->SetTickLength(0.0, "Y");
    hGridReferenceCellSpace->SetStats(kFALSE);
    hGridReferenceCellSpace->SetTickLength(0.0, "X");
    hGridReferenceCellSpace->SetTickLength(0.0, "Y");

    cout<<"End initialise3DGridReference()"<<endl;

}

void TAnalysisOf3dDiamonds::initialise3DYAlignmentHistos() {

    //Fiducial Region with Edge Alignment Regions Highlighted

    //hDeadCellsProfile
    cout<<"settings->getDeadCell3D().size(): "<<settings->getDeadCell3D().size()<<endl;
    for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
        TString name = TString::Format("hDeadCell_no_%d_MeanCharge",i);
        hDeadCellCharge.push_back(new TProfile(name,name,90,0,450));
        hDeadCellCharge.back()->GetXaxis()->SetTitle("pred hit position y in diamond /#mum	");
        hDeadCellCharge.back()->GetYaxis()->SetTitle("Mean Charge [ADC]");
        name = TString::Format("hDeadCell_no_%d_HitPositions",i);
        hDeadCellPositions.push_back(histSaver->GetHistoBinedInCells(name,30));
        hDeadCellPositions.back()->GetXaxis()->SetTitle("pred hit position x in diamond /#mum	");
        hDeadCellPositions.back()->GetYaxis()->SetTitle("pred hit position y in diamond /#mum	");
    }
    cout<<"End initialise3DYAlignmentHistos()"<<endl;
};

void TAnalysisOf3dDiamonds::initialise3DOverviewHistos() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";

    //hDetXvsDetY3DTotolCharge
    TString name = "hPulseHeightVsDetectorHitPostionXY";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    UInt_t factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXY = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXY->GetXaxis()->SetTitle("predicted x position in detector /#mum");
    hPulseHeightVsDetectorHitPostionXY->GetYaxis()->SetTitle("predicted y position in detector /#mum");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetTitle("charge /ADC");
    hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    //hPulseHeightVsDetectorHitPostionXYGoodCells
    name = "hPulseHeightVsDetectorHitPostionXYGoodCells";
    name.Append(appendix);
    if(verbosity>1) cout<<"Create "<<name<<endl;
    factor = 10* (settings->getNQuarters3d()/2);
    hPulseHeightVsDetectorHitPostionXYGoodCells = histSaver->GetProfile2dBinedInCells(name,factor);
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetXaxis()->SetTitle("predicted x position in detector /#mum");
    hPulseHeightVsDetectorHitPostionXYGoodCells->GetYaxis()->SetTitle("predicted y position in detector /#mum");
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
    cout<<"#(#"<<endl;

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
    //	for(UInt_t column=0;column<settings->getNColumns3d();column++){
    //		for(UInt_t row=0;row<settings->getNRows3d();row++){
    //
    //			Int_t cell =Int_t cell = settings->get3DCellNo((int)column,row);
    //			TString name = TString::Format("hCellsLandau_no_%03d",cell);
    //			hCellsLandau.push_back(new TH1F(name,name,64,0,2800));
    //
    //			stringstream hCellsClusteSizeName; hCellsClusteSizeName<<"hCellsClusteSize"<<cell<<FileNameEnd;
    //			hCellsClusteSize.push_back(new TH1F(hCellsClusteSizeName.str().c_str(),hCellsClusteSizeName.str().c_str(),6,0,6));
    //
    //			hQuarterCellsLandau.push_back(vector<TH1F*>());
    //			hQuarterCellsClusterSize.push_back(vector<TH1F*>());
    //
    //			for(int quarter=0;quarter<settings->getNQuarters3d();quarter++){
    //				//cout<<"This should not repeat: "<<((i*settings->getNRows3d()*4)+j*4+k)<<endl;
    //				stringstream hQuaterCellsLandauName; hQuaterCellsLandauName<<"hQuaterCellsLandau"<<((column*settings->getNRows3d()*4)+row*4+quarter)<<FileNameEnd;
    //				hQuarterCellsLandau[cell].push_back(new TH1F(hQuaterCellsLandauName.str().c_str(),hQuaterCellsLandauName.str().c_str(),64,0,2800));
    //
    //				stringstream hQuarterCellsClusterSizeName; hQuarterCellsClusterSizeName<<"hQuarterCellsClusterSize"<<((column*settings->getNRows3d()*4)+row*4+quarter)<<FileNameEnd;
    //				hQuarterCellsClusterSize[cell].push_back(new TH1F(hQuarterCellsClusterSizeName.str().c_str(),hQuarterCellsClusterSizeName.str().c_str(),6,0,6));
    //			}
    //		}
    //	}


    //cout<<"The size of hQuarterCells is: "<<hQuaterCellsLandau.size()<<endl;

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

    //hCellsGoodandBad
    /*stringstream hCellsHarrisGoodName; hCellsHarrisGoodName<<"hCellsHarrisGood"<<FileNameEnd;
	hCellsHarrisGood = new TH1F(hCellsHarrisGoodName.str().c_str(),hCellsHarrisGoodName.str().c_str(),256,0,2800);
	hCellsHarrisGood->SetTitle(hCellsHarrisGoodName.str().c_str());
	hCellsHarrisGood->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hCellsHarrisGood->GetYaxis()->SetTitle("Entries");

	hCellsHarrisBad = (TH1F*)hCellsHarrisGood->Clone("hCellsHarrisBad");
	hCellsHarrisBad->SetTitle("hCellsHarrisBad");
     */
};

void TAnalysisOf3dDiamonds::initialiseEdgeFreeHistos(){
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    Int_t factor = 2;
    TString name ="hPulseHeigthCentralRegion";
    name.Append(appendix);
    hPulseHeigthCentralRegion = histSaver->GetProfile2dBinedInCells(name,factor);
    name ="hPulseHeigthEdgeRegion";
    name.Append(appendix);
    hPulseHeigthEdgeRegion =histSaver->GetProfile2dBinedInCells(name,factor);

    hEventsEdgeRegion = new TH2F("hEventsEdgeRegion","hEventsEdgeRegion",150,0,150,150,0,150);
    hEventsEdgeRegion->GetXaxis()->SetTitle("rel x predicted /#mum");
    hEventsEdgeRegion->GetYaxis()->SetTitle("rel y predicted /#mum");

    hEventsCentralRegion = new TH2F("hEventsCentralRegion","hEventsCentralRegion",150,0,150,150,0,150);
    hEventsCentralRegion->GetXaxis()->SetTitle("rel x predicted /#mum");
    hEventsCentralRegion->GetYaxis()->SetTitle("rel y predicted /#mum");
    hEventsCentralRegion->GetXaxis()->SetLimits(0,150);
    hEventsCentralRegion->GetYaxis()->SetLimits(0,150);


}

void TAnalysisOf3dDiamonds::initialise3DCellOverlayHistos() {
    TH1F* hLandau = new TH1F("hLandau","hLandau",PulseHeightBins, PulseHeightMin,PulseHeightMax);
    hLandau->GetXaxis()->SetTitle("pulse height of cluster /ADC");
    hLandau->GetYaxis()->SetTitle("number of entries #");
	Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
	for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		//Cell Overlay
		//Double_t xBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
		Double_t xBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
		Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

		//Double_t yBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
		Double_t yBinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
		Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;
		//Double_t BinEdges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
		//hCellsOverlayed
		TString name = "hCellsOverlayAvrgCharge";
		name.Append(appendix);
		/*hCellsOverlayAvrgCharge.push_back(new TProfile2D(name,name,
				30,0,settings->GetCellWidth(subjectDetector,2),
				30,0,settings->GetCellHeight()));*/
		hCellsOverlayAvrgCharge.push_back(new TProfile2D(name,name,
				nXbins,xBinEdges,
				nYbins,yBinEdges));
		hCellsOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell /#mum");
		hCellsOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell /#mum");
		hCellsOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
		//	hCellsOverlayAvrgCharge->SetContour(99);
		name = "hCellsOverlayAvrgCharge_noColumnHit";
		name.Append(appendix);
		hCellsOverlayAvrgChargeNoColumnHit.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->SetTitle(name);

		cout<<hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)<<" "<<hCellsOverlayAvrgCharge.at(ClusterSize)<<endl;
		cout<<" "<<hCellsOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;

		//hCellsOverlayAvrgChargeMinusBadCells
		name = "hCellsOverlayAvrgChargeMinusBadCells";
		name.Append(appendix);
		hCellsOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

		name = TString::Format( "hCellLandauChargeMinusBadCells_clustersize%0d",ClusterSize+1);
		name.Append(appendix);
		hCellsLandauMinusBadCells.push_back((TH1F*)hLandau->Clone(name));
		hCellsLandauMinusBadCells.back()->SetTitle(name);


		//hCellsOverlayAvrgChargeGoodCells
		name = "hCellsOverlayAvrgChargeGoodCells";
		name.Append(appendix);
		hCellsOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

		//hCellsOverlayAvrgChargeBadCells
		name = "hCellsOverlayAvrgChargeBadCells";
		name.Append(appendix);
		hCellsOverlayAvrgChargeBadCells.push_back((TProfile2D*)hCellsOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->SetTitle(name);

		//hCellsOverlayedColumnLandau
		name = "hCellOverlayWithColumnLandau";
		name.Append(appendix);
		hCellOverlayWithColumnLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
		hCellOverlayWithColumnLandau.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
		hCellOverlayWithColumnLandau.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

		//hCellsOverlayedLandauNoColumn
		name = "hCellOverlayNoColumnLandau";
		name.Append(appendix);
		hCellOverlayNoColumnLandau.push_back(new TH1F(name,name,256,0,2800));//todo iain ph

		name = "hCellOverlayColumnLandau";
		name.Append(appendix);
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


    cout<<"End initialise3DCellOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DCellCentralColumnOverlayHistos() {

	Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
	for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		Float_t CentralColumnXBins = settings->getCentralColumnOverlayXBins();
		Float_t CentralColumnXLow = settings->getCentralColumnOverlayXLow();
		Float_t CentralColumnXHigh = settings->getCentralColumnOverlayXHigh();
		Float_t CentralColumnYBins = settings->getCentralColumnOverlayYBins();
		Float_t CentralColumnYLow = settings->getCentralColumnOverlayYLow();
		Float_t CentralColumnYHigh = settings->getCentralColumnOverlayYHigh();

		//hCellsOverlayed
		TString name = "hCellsCentralColumnOverlayAvrgCharge";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgCharge.push_back(new TProfile2D(name,name,
				CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
				CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh));
		/*hCellsOverlayAvrgCharge = new TProfile2D(name,name,
				nXbins,xBinEdges,
				nYbins,yBinEdges);*/
		hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell /#mum");
		hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell /#mum");
		hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
		//	hCellsOverlayAvrgCharge->SetContour(99);

		//hCellsOverlayAvrgChargeMinusBadCells
		name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCells";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
		hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

		//hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut
		name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.push_back(new TH2F(name,name,
				CentralColumnXBins+2,CentralColumnXLow-2,CentralColumnXHigh+2,
				CentralColumnYBins+2,CentralColumnYLow-2,CentralColumnYHigh+2));
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell /#mum");
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell /#mum");
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCutEvents.at(ClusterSize)->GetZaxis()->SetTitle("Events");

		//hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut
		name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->SetTitle(name);

		//hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis
		name = "hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis";
		name.Append(appendix);
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.push_back((TProfile2D*)hCellsCentralColumnOverlayAvrgCharge.at(0)->Clone(name));
		hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->SetTitle(name);

		name = (TString)"hLandauMinusBadCellsOffsetAnalysis" + appendix;
		hLandauMinusBadCellsOffsetAnalysis.push_back(new TH1F(name,name,2* PulseHeightBins,-1.*PulseHeightMax,PulseHeightMax));
		hLandauMinusBadCellsOffsetAnalysis.back()->GetXaxis()->SetTitle("pulse height /ADC");
		hLandauMinusBadCellsOffsetAnalysis.back()->GetYaxis()->SetTitle("number of Entries #");

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

		//hCellsOverlayedLandauNoColumn
		name = "hCellsCentralColumnOverlayLandauMinusBadCells";
		name.Append(appendix);
		hCellsCentralColumnOverlayLandauMinusBadCells.push_back(new TH1F(name,name,256,0,2800));//todo iain ph
		hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize)->GetXaxis()->SetTitle("charge / ADC");
		hCellsCentralColumnOverlayLandauMinusBadCells.at(ClusterSize)->GetYaxis()->SetTitle("Entries");

		//hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis
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

	} //End of for ClusterSize

    cout<<"End initialise3DCellCentralColumnOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3DOffsetOverlayHistos() {

	Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
	for(Int_t ClusterSize = 0; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		/*With offset of 37.5, central column bin is at: 107.5 - 117.5.
    					 corner column bin is at: 32.5 - 42.5.*/

		Double_t xBinEdges[] = {0,12.5,22.5,32.5,42.5,52.5,62.5,77.5,87.5,97.5,107.5,117.5,127.5,137.5,150};
		//Double_t xBinEdges[] = {0,15,30,45,60,75,90,105,120,135,150};
		Int_t  nXbins = sizeof(xBinEdges)/sizeof(Double_t) - 1;

		//Double_t yBinEdges[] = {0,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,150};
		Double_t yBinEdges[] = {0,12.5,22.5,32.5,42.5,52.5,62.5,77.5,87.5,97.5,107.5,117.5,127.5,137.5,150};
		//Double_t yBinEdges[] = {0,15,30,45,60,75,90,105,120,135,150};
		Int_t nYbins = sizeof(yBinEdges)/sizeof(Double_t) - 1;

		//hCellsOverlayed
		TString name = "hCellsOffsetOverlayAvrgCharge";
		name.Append(appendix);
		/*hCellsCentralColumnOverlayAvrgCharge = new TProfile2D(name,name,
    		CentralColumnXBins,CentralColumnXLow,CentralColumnXHigh,
    		CentralColumnYBins,CentralColumnYLow,CentralColumnYHigh);*/
		hCellsOffsetOverlayAvrgCharge.push_back(new TProfile2D(name,name,
				nXbins,xBinEdges,
				nYbins,yBinEdges));
		hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetXaxis()->SetTitle("rel. x position within a cell /#mum");
		hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetYaxis()->SetTitle("rel. y position within a cell /#mum");
		hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetZaxis()->SetTitle("pulse height of cluster /ADC");
		//	hCellsOverlayAvrgCharge->SetContour(99);

		//hCellsOverlayAvrgChargeMinusBadCells
		name = "hCellsOffsetOverlayAvrgChargeMinusBadCells";
		name.Append(appendix);
		hCellsOffsetOverlayAvrgChargeMinusBadCells.push_back((TProfile2D*)hCellsOffsetOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->SetTitle(name);

		//hCellsOverlayAvrgChargeGoodCells
		name = "hCellsOffsetOverlayAvrgChargeGoodCells";
		name.Append(appendix);
		hCellsOffsetOverlayAvrgChargeGoodCells.push_back((TProfile2D*)hCellsOffsetOverlayAvrgCharge.at(0)->Clone(name));
		hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->SetTitle(name);

		/*//hCellsOverlayedColumnLandau
	name = "hCellsCentralColumnOverlayLandau";
	name.Append(appendix);
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

    cout<<"End initialise3DCellCentralColumnOverlayHistos()"<<endl;
}

void TAnalysisOf3dDiamonds::initialise3D2DLandauAndClustersizeHistos() {

    /*//hCellsLandau2D
	stringstream hCellsLandau2DName; hCellsLandau2DName<<"hCellsLandau2D"<<FileNameEnd;
	hCellsLandau2D = new TH2D(hCellsLandau2DName.str().c_str(),hCellsLandau2DName.str().c_str(),256,0,2800,99,0,99);
	hCellsLandau2D->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2D->GetYaxis()->SetTitle("Cell");
	hCellsLandau2DQuarterFail = new TH2D("hCellsLandau2DQuarterFail","hCellsLandau2DQuarterFail",1,0,2800,99,0,99);
	hCellsLandau2DQuarterFail->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2DQuarterFail->GetYaxis()->SetTitle("Cell");
	//hCellsLandau2D->SetCanExtend(TH1F::kAllAxes);


	//h2DClusterSizeClone
	h2DClusterSizeClone = (TH2D*)h2DClusterSize->Clone("h2DClusterSizeClone");

	//h2DClusterSizeClone1
	h2DClusterSizeClone1 = (TH2D*)h2DClusterSize->Clone("h2DClusterSizeClone1");

	//h2DClusterSizeQuarterCell
	stringstream h2DClusterSizeQuarterCellName; h2DClusterSizeQuarterCellName<<"h2DClusterSizeQuarterCell"<<FileNameEnd;
	h2DClusterSizeQuarterCell = new TH2D(h2DClusterSizeQuarterCellName.str().c_str(),h2DClusterSizeQuarterCellName.str().c_str(),20,0,20,5,0,5);
	h2DClusterSizeQuarterCell->GetXaxis()->SetTitle("");
	h2DClusterSizeQuarterCell->GetYaxis()->SetTitle("ClusterSize");
	vector<char> PassFail;
	PassFail.push_back('P');PassFail.push_back('P');PassFail.push_back('P');PassFail.push_back('P');
	PassFail.push_back('F');PassFail.push_back('P');PassFail.push_back('P');PassFail.push_back('P');
	PassFail.push_back('F');PassFail.push_back('F');PassFail.push_back('P');PassFail.push_back('P');
	PassFail.push_back('F');PassFail.push_back('F');PassFail.push_back('F');PassFail.push_back('P');
	PassFail.push_back('F');PassFail.push_back('F');PassFail.push_back('F');PassFail.push_back('F');
	for(int i=0;i<20;i++){
		stringstream iLetter; iLetter<<PassFail.at(i);
		h2DClusterSizeQuarterCell->GetXaxis()->SetBinLabel(i+1,iLetter.str().c_str());
	}
	for(int i=0;i<5;i++){
		stringstream jNumber;
		if(i==5)
			jNumber<<"5+";
		else{jNumber<<(i+1);}
		h2DClusterSizeQuarterCell->GetYaxis()->SetBinLabel(i+1,jNumber.str().c_str());
	}
	h2DClusterSizeQuarterCell->SetContour(99);
	//h2DClusterSize->SetTitleOffset(0.015,"X");

	//h2DClusterSizeQuarterCellClone
	h2DClusterSizeQuarterCellClone = (TH2D*)h2DClusterSizeQuarterCell->Clone("h2DClusterSizeQuarterCellClone");

	//h2DClusterSizeQuarterCellClone1
	h2DClusterSizeQuarterCellClone1 = (TH2D*)h2DClusterSizeQuarterCell->Clone("h2DClusterSizeQuarterCellClone1");

	//h2DMeanClusterSizeQuarterCell
	h2DMeanClusterSizeQuarterCell = (TH2D*)h2DClusterSizeQuarterCell->Clone("h2DMeanClusterSizeQuarterCell");

	//h2DMeanClusterSizeQuarterCellTotal
	h2DMeanClusterSizeQuarterCellTotal = (TH2D*)h2DClusterSizeQuarterCell->Clone("h2DMeanClusterSizeQuarterCellTotal");

	//h2DMeanClusterSizeQuarterCellEvents
	h2DMeanClusterSizeQuarterCellEvents = (TH2D*)h2DClusterSizeQuarterCell->Clone("h2DMeanClusterSizeQuarterCellEvents");

	h2DClusterSizeXAxis = new TH2D(h2DClusterSizeName.str().c_str(),h2DClusterSizeName.str().c_str(),5,0,20,5,0,5);
	for(int i=0;i<5;i++){
		stringstream jNumber; jNumber<<(i+1);
		h2DClusterSizeXAxis->GetXaxis()->SetBinLabel(i+1,jNumber.str().c_str());
	}
	h2DClusterSizeXAxis->GetXaxis()->SetTitle("Failed Quarters");
	h2DClusterSizeXAxis->SetLabelOffset(0.015,"X");
	//h2DClusterSizeXAxis->GetXaxis()->TitleOffset()

	cout<<"End initialise3D2DLandauClusterSizeHistos()"<<endl;
     */
}

void TAnalysisOf3dDiamonds::initialiseTransparentAnalysisHistos() {

    hTransparentAnalysisInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hTransparentAnalysisInvalidCluster");
    hTransparentAnalysisInvalidCluster->SetTitle("hTransparentAnalysisInvalidCluster");

    hTransparentAnalysisValidCluster = (TH2F*) hValidEventsDetSpace->Clone("hTransparentAnalysisValidCluster");
    hTransparentAnalysisValidCluster->SetTitle("hTransparentAnalysisValidCluster");

    TString name = "hTransparentAnalysisValidClusterFidCutXvsFidCutY";
    hTransparentAnalysisValidClusterFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->GetZaxis()->SetTitle("Charge ADC");

    /*stringstream hNonValidTransparentAnalysisName; hNonValidTransparentAnalysisName<<"hNonValidTransparentAnalysis"<<FileNameEnd;
	hNonValidTransparentAnalysis = new TH2D(hNonValidTransparentAnalysisName.str().c_str(),hNonValidTransparentAnalysisName.str().c_str(),840,0,4200,500,-500,2000);
	hNonValidTransparentAnalysis->GetXaxis()->SetTitle("xPredDet [um]");
	hNonValidTransparentAnalysis->GetYaxis()->SetTitle("yPredDet [um]");*/


    //Universal histograms
    /*
	//hNumberofClusters
	stringstream hNumberofClustersName; hNumberofClustersName<<"hNumberofClusters"<<FileNameEnd;
	hNumberofClusters = new TH1F(hNumberofClustersName.str().c_str(),hNumberofClustersName.str().c_str(),4,0,4);
	hNumberofClusters->SetTitle(hNumberofClustersName.str().c_str());
	hNumberofClusters->GetXaxis()->SetTitle("Number of Clusters");
	hNumberofClusters->GetYaxis()->SetTitle("Number of Entries #");

	//hEventsvsChannelCombined
	stringstream hEventsvsChannelCombinedName; hEventsvsChannelCombinedName<<"hEventsvsChannelCombined"<<FileNameEnd;
	hEventsvsChannelCombined = new TH1F(hEventsvsChannelCombinedName.str().c_str(),hEventsvsChannelCombinedName.str().c_str(),100,0,100);
	hEventsvsChannelCombined->SetTitle(hEventsvsChannelCombinedName.str().c_str());
	hEventsvsChannelCombined->GetXaxis()->SetTitle("Channel");
	hEventsvsChannelCombined->GetYaxis()->SetTitle("Number of Entries #");

	//hDoubleClusterPos
	stringstream hDoubleClusterPosName; hDoubleClusterPosName<<"hDoubleClusterPos"<<FileNameEnd;
	hDoubleClusterPos = new TH1F(hDoubleClusterPosName.str().c_str(),hDoubleClusterPosName.str().c_str(),80,20,100);
	hDoubleClusterPos->SetTitle(hDoubleClusterPosName.str().c_str());
	hDoubleClusterPos->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos->GetYaxis()->SetTitle("Number of Entries #");

	//hDoubleClusterPos0
	stringstream hDoubleClusterPos0Name; hDoubleClusterPos0Name<<"hDoubleClusterPos0"<<FileNameEnd;
	hDoubleClusterPos0 = new TH1F(hDoubleClusterPos0Name.str().c_str(),hDoubleClusterPos0Name.str().c_str(),80,20,100);
	hDoubleClusterPos0->SetTitle(hDoubleClusterPos0Name.str().c_str());
	hDoubleClusterPos0->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos0->GetYaxis()->SetTitle("Number of Entries #");
	hDoubleClusterPos0->SetFillColor(2);

	//hDoubleClusterPos1
	stringstream hDoubleClusterPos1Name; hDoubleClusterPos1Name<<"hDoubleClusterPos1"<<FileNameEnd;
	hDoubleClusterPos1 = new TH1F(hDoubleClusterPos1Name.str().c_str(),hDoubleClusterPos1Name.str().c_str(),80,20,100);
	hDoubleClusterPos1->SetTitle(hDoubleClusterPos1Name.str().c_str());
	hDoubleClusterPos1->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos1->GetYaxis()->SetTitle("Number of Entries #");
	hDoubleClusterPos1->SetFillColor(3);

	//hLandauCluster1
	stringstream hLandauCluster1Name; hLandauCluster1Name<<"hLandauCluster1"<<FileNameEnd;
	hLandauCluster1 = new TH1F(hLandauCluster1Name.str().c_str(),hLandauCluster1Name.str().c_str(),256,0,2800);
	hLandauCluster1->SetTitle(hLandauCluster1Name.str().c_str());
	hLandauCluster1->GetXaxis()->SetTitle("Number of Clusters");
	hLandauCluster1->GetYaxis()->SetTitle("Number of Entries #");
	hLandauCluster1->SetFillColor(2);

	//hLandauCluster2
	stringstream hLandauCluster2Name; hLandauCluster2Name<<"hLandauCluster2"<<FileNameEnd;
	hLandauCluster2 = new TH1F(hLandauCluster2Name.str().c_str(),hLandauCluster2Name.str().c_str(),256,0,2800);
	hLandauCluster2->SetTitle(hLandauCluster2Name.str().c_str());
	hLandauCluster2->GetXaxis()->SetTitle("Number of Clusters");
	hLandauCluster2->GetYaxis()->SetTitle("Number of Entries #");
	hLandauCluster2->SetFillColor(3);

	//hLandauDoubleCombined
	stringstream hLandauDoubleCombinedName; hLandauDoubleCombinedName<<"hLandauDoubleCombined"<<FileNameEnd;
	hLandauDoubleCombined = new TH1F(hLandauDoubleCombinedName.str().c_str(),hLandauDoubleCombinedName.str().c_str(),256,0,2800);
	hLandauDoubleCombined->SetTitle(hLandauDoubleCombinedName.str().c_str());
	hLandauDoubleCombined->GetXaxis()->SetTitle("Number of Clusters");
	hLandauDoubleCombined->GetYaxis()->SetTitle("Number of Entries #");
     */
    for(int i=0; i<settings->diamondPattern.getNIntervals(); i++){
        pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
        //hLandau
        TString name = TString::Format("hLandauTransparent_pattern_%d_ch_%d_to_%d",i+1,channels.first,channels.second);
        hLandauTransparent.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        hLandauTransparent.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
        hLandauTransparent.at(i)->GetYaxis()->SetTitle("number of entries #");
        hLandauTransparent.at(i)->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);

        //hLandauBadCellsRemoved
        name = TString::Format("hLandauTransparentBadCellsRemoved_pattern_%d_ch_%d_to_%d",i+1,channels.first,channels.second);
        hLandauTransparentBadCellsRemoved.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        hLandauTransparentBadCellsRemoved.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
        hLandauTransparentBadCellsRemoved.at(i)->GetYaxis()->SetTitle("number of entries #");
        hLandauTransparentBadCellsRemoved.at(i)->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);

        /*
		//hPHvsChannel
		stringstream hEventsvsChannelName; hEventsvsChannelName<<"hEventsvsChannel%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hEventsvsChannel.push_back(new TH1F(hEventsvsChannelName.str().c_str(),hEventsvsChannelName.str().c_str(),100,0,100));
		hEventsvsChannel.at(i)->GetXaxis()->SetTitle("HighestPH [ch]");
		hEventsvsChannel.at(i)->GetYaxis()->SetTitle("No. Events");

		//hPHvsChannel
		stringstream hPHvsChannelName; hPHvsChannelName<<"hPHvsChannel%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hPHvsChannel.push_back(new TH2F(hPHvsChannelName.str().c_str(),hPHvsChannelName.str().c_str(),150,0,2900,100,0,100));
		hPHvsChannel.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsChannel.at(i)->GetYaxis()->SetTitle("XPos(channel)");
		hPHvsChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsChannel.at(i)->GetXaxis()->SetRangeUser(1,2600);
		//hPHvsChannel.at(i)->SetMaximum(3000);
		//hPHvsChannel.at(i)->SetMinimum(0);



		//hFidCutXvsFidCutY
		stringstream hFidCutXvsFidCutYName; hFidCutXvsFidCutYName<<"hFidCutXvsFidCutY%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hFidCutXvsFidCutY.push_back(new TH2F(hFidCutXvsFidCutYName.str().c_str(),hFidCutXvsFidCutYName.str().c_str(),160,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),120,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
		hFidCutXvsFidCutY.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutY.at(i)->GetYaxis()->SetTitle("FidCutY");

		//hFidCutXvsFidCutYvsCharge		For TH2D
		stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hFidCutXvsFidCutYvsCharge.push_back(new TH2D(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
		hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

		//hFidCutXvsFidCutYvsEvents
		hFidCutXvsFidCutYvsEvents.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsEvents"));

		//hFidCutXvsFidCutYvsMeanCharge
		hFidCutXvsFidCutYvsMeanCharge.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsMeanCharge"));
		hFidCutXvsFidCutYvsMeanCharge.at(i)->GetZaxis()->SetRangeUser(0,1500);
         */


        //hXdetvsYdetvsCharge		For TH2D
        Int_t DiamondPattern = i+1;
        Float_t xLow = settings->get3dMetallisationFidCuts()->getXLow(DiamondPattern);
        Float_t xHigh = settings->get3dMetallisationFidCuts()->getXHigh(DiamondPattern);
        Float_t yLow = settings->get3dMetallisationFidCuts()->getYLow(DiamondPattern);
        Float_t yHigh = settings->get3dMetallisationFidCuts()->getYHigh(DiamondPattern);
        cout<<"("<<xLow<<"-"<<xHigh<<", "<<yLow<<"-"<<yHigh<<")"<<endl;
        Float_t xDiv = (xHigh - xLow)/5;
        Float_t yDiv = (yHigh - yLow)/5;
        stringstream hXdetvsYdetvsChargeName; hXdetvsYdetvsChargeName<<"hXdetvsYdetvsCharge%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
        hXdetvsYdetvsCharge.push_back(new TH2D(hXdetvsYdetvsChargeName.str().c_str(),hXdetvsYdetvsChargeName.str().c_str(),xDiv,xLow,xHigh,yDiv,yLow,yHigh));
        hXdetvsYdetvsCharge.at(i)->GetXaxis()->SetTitle("X (um)");
        hXdetvsYdetvsCharge.at(i)->GetYaxis()->SetTitle("Y (um)");
        hXdetvsYdetvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        hXdetvsYdetvsEvents.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone("hXdetvsYdetvsEvents"));

        //hFidCutXvsFidCutYvsMeanCharge
        hXdetvsYdetvsMeanCharge.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone("hXdetvsYdetvsMeanCharge"));
        hXdetvsYdetvsMeanCharge.at(i)->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    }

    //hCellOverlay
    stringstream hCellOverlayvsChargeName; hCellOverlayvsChargeName<<"hCellOverlayvsCharge"<<FileNameEnd;
    hCellOverlayvsCharge = new TH2D(hCellOverlayvsChargeName.str().c_str(),hCellOverlayvsChargeName.str().c_str(),30,0,150,30,0,150);
    hCellOverlayvsCharge->GetXaxis()->SetTitle("X (um)");
    hCellOverlayvsCharge->GetYaxis()->SetTitle("Y (um)");
    hCellOverlayvsCharge->GetZaxis()->SetTitle("Charge ADC");

    //hCellOverlayvsEvents
    hCellOverlayvsEvents = (TH2D*)hCellOverlayvsCharge->Clone("hCellOverlayvsEvents");

    //hCellOverlayvsMeanCharge
    hCellOverlayvsMeanCharge = (TH2D*)hCellOverlayvsCharge->Clone("hCellOverlayvsCharge");

    /*
	//hFidCutXvsFidCutYvsMeanChargeAllDetectors
	hFidCutXvsFidCutYvsMeanChargeAllDetectors = (TH2D*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("hFidCutXvsFidCutYvsMeanChargeAllDetectors");

	for(int i=0;i<7;i++){
		//hFidCutXvsFidCutYClusters	For TH2D	{0,1,1_1Seed,2_FirstCluster,2_SecondCluster,3}
		stringstream hFidCutXvsFidCutYClustersName; hFidCutXvsFidCutYClustersName<<"hFidCutXvsFidCutYClusters"<<i<<FileNameEnd;
		hFidCutXvsFidCutYClusters.push_back(new TH2D(hFidCutXvsFidCutYClustersName.str().c_str(),hFidCutXvsFidCutYClustersName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
		hFidCutXvsFidCutYClusters.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYClusters.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYClusters.at(i)->GetZaxis()->SetTitle("Events");
	}
     */
}

void TAnalysisOf3dDiamonds::SaveShortAnalysisHistos() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";

    ShortAnalysis_SaveMeanChargeVector();
    ShortAnalysis_Save2ClusterPlots();
    vector<Float_t> xPred;
    vector<Float_t> yPred;
    vector<Float_t> charge;
    histSaver->SaveHistogram(hRelativeChargeTwoClustersX);
    histSaver->SaveHistogram(hRelativeChargeTwoClustersY);
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
    c1 = new TCanvas(name,name);
    c1->cd();
    hTotalAvrgChargeXY->Draw("colz");
    hTotalAvrgChargeXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    histSaver->SaveCanvas(c1);
    delete c1;

    for(UInt_t i = 0; i < vecEdgePredX.size(); i++){
        xPred.insert(  xPred.end(), vecEdgePredX[i].begin(), vecEdgePredX[i].end());
        yPred.insert(  yPred.end(), vecEdgePredY[i].begin(), vecEdgePredY[i].end());
        charge.insert(charge.end(), vecEdgePulseHeight[i].begin(), vecEdgePulseHeight[i].end());
    }
    //a.end(), b.begin(), b.end());
    histSaver->SaveHistogram(histSaver->CreateScatterHisto("hEdgeFittingPredictedPosition",yPred,xPred),false);
    TH3F* hEdgeFittingCharge = histSaver->Create3DHisto("hEdgeFittingCharge",xPred,yPred,charge);
    TH2F* hEdgeFittingAvrgCharge = (TH2F*) hEdgeFittingCharge->Project3DProfile("yx");
    histSaver->SaveHistogram(hEdgeFittingAvrgCharge);

    for(int i = 0; i < vecEdgePredX.size(); i++){
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
        title.Append(" /#mum");
        hEdgeFittingAvrgCharge->GetXaxis()->SetTitle(title);//"predicted Position X /#mum");
        histSaver->SaveHistogram(hEdgeFittingAvrgCharge);
        TH1F* hEdgeFittingAvrgCharge_pfx = (TH1F*)hEdgeFittingAvrgCharge->ProfileX();
        if(hEdgeFittingAvrgCharge_pfx){
            hEdgeFittingAvrgCharge_pfx->GetYaxis()->SetTitle("avrg. Charge / ADC");
            TCutG *cut = this->settings->getEdgePosition(i);
            name = "c";
            name.Append(hEdgeFittingAvrgCharge_pfx->GetName());
            TCanvas *c1 = new TCanvas(name,name);
            hEdgeFittingAvrgCharge_pfx->Draw();
            if(cut)cut->Draw();
            histSaver->SaveCanvas(c1);
            delete c1;
        }
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
            cout<<name<<endl;
            if(i==0)
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name);
            else{
                //				int bin = hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name,i,i);
            }
            name = hShortAnalysis2ClusterHitPattern_2ndCluster->GetName();
            name.Append(extension);
            cout<<name<<endl;

            if(i==0)
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name);
            else
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name,i,i);
            name = "h2ClusterAnalysis_ClusterPatterns";
            name.Append(extension);
            histSaver->SaveTwoHistos((string)name,histo1st_py,histo2nd_py);
            histSaver->SaveHistogram(histo1st_py);
            histSaver->SaveHistogram(histo2nd_py);
            delete histo1st_py;
            delete histo2nd_py;
        }
        name = "c2ClusterAnalysis_ClusterPatterns";
        name.Append(appendix);
        TH1D* histo1st_px = hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionX();//name);
        TH1D* histo2nd_px = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionX();//name);
        histSaver->SaveTwoHistos((string)name,histo1st_px,histo2nd_px);
        histSaver->SaveHistogram(histo1st_px);
        histSaver->SaveHistogram(histo2nd_px);
        //		if(histo1st_px) delete histo1st_px;
        //		if(histo2nd_px) delete histo2nd_px;
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
        /*//hLandau
		stringstream hLandauName; hLandauName<<"hLandau%%"<<firstCh<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hLandau.push_back(HistogrammSaver::CreateDistributionHisto(hLandauName.str().c_str(),*vecPHDiamondHit.at(i),256));
		hLandau.at(i)->SetTitle(hLandauName.str().c_str());
		hLandau.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
		hLandau.at(i)->GetYaxis()->SetTitle("number of entries #");
         */
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

    TCanvas* cTransparentAnalysisInvalidCluster = new TCanvas("cTransparentAnalysisInvalidCluster", "cTransparentAnalysisInvalidCluster");
    cTransparentAnalysisInvalidCluster->cd();
    hTransparentAnalysisInvalidCluster->Draw();
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisInvalidCluster);
    histSaver->SaveCanvas(cTransparentAnalysisInvalidCluster);

    TCanvas* cTransparentAnalysisValidCluster = new TCanvas("cTransparentAnalysisValidCluster", "cTransparentAnalysisValidCluster");
    cTransparentAnalysisValidCluster->cd();
    hTransparentAnalysisValidCluster->Draw();
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisValidCluster);
    histSaver->SaveCanvas(cTransparentAnalysisValidCluster);

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

    hTransparentAnalysisValidCluster->GetXaxis()->SetRangeUser(xmin-deltaX,xmax+deltaX);
    hTransparentAnalysisValidCluster->GetYaxis()->SetRangeUser(ymin-deltaY,ymax+deltaY);
    hTransparentAnalysisValidCluster->Draw("colz");
    settings->DrawMetallisationGrid(cTransparentAnalysisValidCluster,3);
    cTransparentAnalysisValidCluster->Update();
    cTransparentAnalysisValidCluster->SetName("cTransparentAnalysisValidCluster_3DwH");
    histSaver->SaveCanvas(cTransparentAnalysisValidCluster);
    delete cTransparentAnalysisValidCluster;

    TCanvas* cTransparentAnalysisValidClusterFidCutXvsFidCutY = new TCanvas("cTransparentAnalysisValidClusterFidCutXvsFidCutY", "cTransparentAnalysisValidClusterFidCutXvsFidCutY");
    cTransparentAnalysisValidClusterFidCutXvsFidCutY->cd();
    hTransparentAnalysisValidClusterFidCutXvsFidCutY->Draw();
    settings->getSelectionFidCuts()->DrawFiducialCutsToCanvas(cTransparentAnalysisValidClusterFidCutXvsFidCutY);
    histSaver->SaveCanvas(cTransparentAnalysisValidClusterFidCutXvsFidCutY);

}

void TAnalysisOf3dDiamonds::LongAnalysisSaveCellAndQuaterNumbering(){

    TString name = "h3DCellNumbering";
    TH2D* hCellNumbering = new TH2D(name,name,
            settings->getNColumns3d(),settings->get3dMetallisationFidCuts()->getXLow(3),settings->get3dMetallisationFidCuts()->getXHigh(3),
            settings->getNRows3d(),settings->get3dMetallisationFidCuts()->getYLow(3),settings->get3dMetallisationFidCuts()->getYHigh(3));
    hCellNumbering->GetXaxis()->SetTitle("Xdet (#mum)");
    hCellNumbering->GetYaxis()->SetTitle("Ydet (#mum)");

    name = "hQuarterNumbering";
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
    histSaver->SaveHistogramWithCellGrid(hQuarterNumbering,hQuarterNumbering);
}

void TAnalysisOf3dDiamonds::SaveLongAnalysisHistos() {
    //	LongAnalysis_SaveCellsOverlayMeanCharge();
    if(!settings->do3dTransparentAnalysis()){
        LongAnalysisSaveCellAndQuaterNumbering();
        LongAnalysis_SaveDeadCellProfile();
        LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors();
        LongAnalysis_SaveFailedQuarters();
        LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail();
    }
    LongAnalysis_SaveEdgeFreeHistos();
    LongAnalysis_SaveGoodAndBadCellLandaus();
    LongAnalysis_SaveCellsOverlayMeanCharge();
    LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge();
    LongAnalysis_SaveCellsOverlayOffsetMeanCharge();
    LongAnalysis_SaveMeanChargePlots();
    LongAnalysis_CreateRelativeAddedTransparentChargeComparisonPlots();
    LongAnalysis_SaveRelativeAddedTransparentCharge();
    //LongAnalysis_SaveCellsClusterSize2DVsGrading();
    //LongAnalysis_SaveQuarterCellsClusterSize2DVsGrading();

    histSaver->SaveHistogram(hLongAnalysisInvalidCellNo);
    histSaver->SaveHistogram(hLongAnalysisInvalidCluster);
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

        cout<<"hQuarterCellsLandausSorted.at(quarter): "<<hQuarterCellsLandausSorted.at(quarter)->GetMean()<<"\t";
        Float_t QuarterMean = hQuarterLandaus[quarter]->GetMean();
        int entries = hQuarterLandaus[quarter]->GetEntries();
        Float_t fluctuation = (compareQuarterTo - QuarterMean)/compareQuarterTo;
        cout<<TString::Format("\t%d: %.1f/%.1f with %3d (%2.1f%%)",
                quarter,QuarterMean,compareQuarterTo,entries,fluctuation)
        <<endl;
        vecFluctuations.push_back(fluctuation);
    }
    return vecFluctuations;

}

void TAnalysisOf3dDiamonds::LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors(int quarterFailCriteriaTyp) {

    TString name = TString::Format("hQuarterFailCriteria_%d",quarterFailCriteriaTyp);
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
    }
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveEdgeFreeHistos() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix = "_trans";
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
        TString name = (TString)hPulseHeigthCentralRegion->GetName()+"_Cell";
        hPulseHeigthCentralRegionCell = (TProfile2D*)hPulseHeigthCentralRegion->Rebin2D(2,2,name);
        hPulseHeigthCentralRegionCell->Sumw2();
        histSaver->SaveHistogramWithCellGrid(hPulseHeigthCentralRegionCell);
    }
    if(hPulseHeigthEdgeRegion){
        TString name =(TString)hPulseHeigthEdgeRegion->GetName()+"_Cell";
        hPulseHeigthEdgeRegionCell = (TProfile2D*)hPulseHeigthEdgeRegion->Rebin2D(2,2,name);
        hPulseHeigthEdgeRegionCell->Sumw2();
        histSaver->SaveHistogramWithCellGrid(hPulseHeigthEdgeRegionCell);
    }
    if(hPulseHeigthEdgeRegionCell && hPulseHeigthCentralRegionCell){
        TString name = "hPulseHeigthCompareRegion_Cell";
        name.Append(appendix);
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
    c1 = histSaver->DrawHistogramWithCellGrid(hLongAnalysisQuarterFluctuations,hLongAnalysisQuarterFluctuations);
    histSaver->DrawFailedQuarters(failedQuarters,c1);
    histo->Draw("sameTEXT");
    histSaver->SaveCanvas(c1);
    delete c1;
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail() {

    Int_t yBins = hCellsLandau.at(0)->GetNbinsX();
    Float_t yMin = hCellsLandau.at(0)->GetBinLowEdge(1);
    Float_t yMax = hCellsLandau.at(0)->GetBinLowEdge(yBins) + hCellsLandau.at(0)->GetBinWidth(yBins);
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix+="_trans";
    TString name = "hCellsLandau2DHighlightedQuarterFail";
    name.Append(appendix);
    TH2D* hCellsLandau2DHighlightedQuarterFail = new TH2D(name,name,settings->GetNCells3d(),0,settings->GetNCells3d(),yBins,yMin,yMax);
    hCellsLandau2DHighlightedQuarterFail->GetXaxis()->SetTitle("Cell");
    hCellsLandau2DHighlightedQuarterFail->GetYaxis()->SetTitle("Charge ADC");
    name ="hCellsLandau2DHighlightedQuarterFailSorted";
    name.Append(appendix);
    TH2D* hCellsLandau2DHighlightedQuarterFailSorted = (TH2D*)hCellsLandau2DHighlightedQuarterFail->Clone(name);
    vector<TH1F*> hCellLandausSorted = hCellsLandau;
    sort(hCellLandausSorted.begin(), hCellLandausSorted.end(), TSettings::SorterForPulseHeightOfHisto);

    name = "hHighlightedQuarterFail";
    name.Append(appendix);
    TH2D* hHighlightedQuarterFail = new TH2D(name,name,settings->GetNCells3d(),0,settings->GetNCells3d(),1,yMin,yMax);
    name = "hHighlightedQuarterFailSorted";
    name.Append(appendix);
    TH2D* hHighlightedQuarterFailSorted = (TH2D*) hHighlightedQuarterFail->Clone(name);
    //hCellsLandau2DQuarterFail->GetXaxis()->SetTitle("Charge ADC");
    //hCellsLandau2DQuarterFail->GetYaxis()->SetTitle("Cell");
    for (UInt_t pos = 0 ; pos < settings->GetNCells3d(); pos ++){
        Int_t cellSorted = -1;
        string title = hCellLandausSorted[pos]->GetName();
        size_t found=title.find_last_of('_');
        if (found != string::npos){
            title = title.substr(found+1);
            cellSorted = atoi(title.c_str());
        }
        TString binLabel = TString::Format("%3d",cellSorted);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetBinLabel(pos,binLabel);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetLabelSize(0.02);
        binLabel = TString::Format("%3d",pos);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetBinLabel(pos,binLabel);
        hCellsLandau2DHighlightedQuarterFailSorted->GetXaxis()->SetLabelSize(0.02);

        for(int yBin =0; yBin<yBins; yBin++){
            hCellsLandau2DHighlightedQuarterFail->SetBinContent(pos,yBin,hCellsLandau[pos]->GetBinContent(yBin));
            hCellsLandau2DHighlightedQuarterFailSorted->SetBinContent(pos,yBin,hCellLandausSorted[pos]->GetBinContent(yBin));
        }
        hHighlightedQuarterFail->SetBinContent(pos,1,CellGrading.at(pos));
        hHighlightedQuarterFailSorted->SetBinContent(pos,1,CellGrading.at(cellSorted));
    }

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
    hCellsOverlayAvrgCharge.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
    hCellOverlayWithColumnLandau.at(ClusterSize)->Fill(clusterCharge);

    if(!settings->isBadCell(3,cellNo)){
        hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
        hCellsLandauMinusBadCells.at(ClusterSize)->Fill(clusterCharge);
    }
    if(settings->IsGoodCell(3, cellNo))
        hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,clusterCharge);
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
        Float_t clusterCharge, Float_t ClusterSize, TCluster diamondCluster) {

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

	}

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
void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayCentralColumnHistosOffsetAnalysis(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize, TCluster *diamondCluster) {

    diamondCluster->SetTransparentClusterSize(5);
//    MakeGhostCluster(diamondCluster,5);

	diamondCluster->SetTransparentClusterSize(ClusterSize);
	GhostCluster.SetTransparentClusterSize(ClusterSize);
	Float_t GhostClusterCharge = GhostCluster.getCharge(false);
	hCellsCentralColumnOverlayAvrgChargeMinusBadCellsOffsetAnalysis.at(ClusterSize)->Fill(xRelPosDet,yRelPosDet,GhostClusterCharge);
	hCellsCentralColumnOverlayLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->Fill(GhostClusterCharge);
    hLandauMinusBadCellsOffsetAnalysis.at(ClusterSize)->Fill(GhostClusterCharge);
}

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayOffsetHistos(Int_t cellNo,Float_t xRelPosDet,Float_t yRelPosDet,
        Float_t clusterCharge, Float_t ClusterSize) {

	int row = settings->getRowOfCell(cellNo);
	int column = settings->getColumnOfCell(cellNo);
	//cout<<"Row of cell 98: "<<settings->getRowOfCell(98)<<"  Column of cell 98: "<<settings->getColumnOfCell(98)<<endl;
	//if cell not in top row or right column then continue
	if(row != 10 && column !=8){

		Float_t pw = settings->getPitchWidth(subjectDetector,2);
		//cout<<"pw: "<<pw<<endl;
		Float_t xRelPosDetOffset = xRelPosDet + settings->getOverlayOffsetX();
		if(xRelPosDetOffset>pw)
			xRelPosDetOffset = xRelPosDet + settings->getOverlayOffsetX() - pw;
		Float_t yRelPosDetOffset = yRelPosDet + settings->getOverlayOffsetY();
		if(yRelPosDetOffset>pw)
			yRelPosDetOffset = yRelPosDet + settings->getOverlayOffsetY() - pw;

		hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->Fill(xRelPosDetOffset,yRelPosDetOffset,clusterCharge);
		//hCellsOffsetOverlayLandau->Fill(clusterCharge);

		if(!settings->isBadCell(3,cellNo)){
			hCellsOffsetOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Fill(xRelPosDetOffset,yRelPosDetOffset,clusterCharge);
			//hCellsOffsetOverlayLandauMinusBadCells->Fill(clusterCharge);
		}
		if(settings->IsGoodCell(3, cellNo)){
			hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->Fill(xRelPosDetOffset,yRelPosDetOffset,clusterCharge);
			//hCellsOffsetOverlayLandauGoodCells->Fill(clusterCharge);
		}
	}
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
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    TList *listGoodCells = new TList;
    TList *listBadCells = new TList;

    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        for(UInt_t row=0;row<settings->getNRows3d();row++){
            Int_t cell = settings->get3DCellNo((int)column,row);
            //hCellNumbering->SetBinContent(column+1,row+1,cell); //This should be a clone of the 2D Cell Mean Charge Plot, Wait till Felix has finished.
            histSaver->SaveHistogram(hCellsLandau.at(cell));
            for(UInt_t i=0; i<settings->getBadCells3D().size(); i++)
                if(cell==settings->getBadCells3D().at(i)){
//                    Int_t Entries = hLandauBadCells->GetEntries();
//                    hLandauBadCells->Add(hCellsLandau.at(cell),1);	//Not working for some reason, ask Felix
//                    hLandauBadCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                    listBadCells->Add(hCellsLandau.at(cell));
//                    cout<<"\tbad"<<endl;
                }
            for(UInt_t i=0; i<settings->getGoodCellRegions3d().size(); i++){
                for(UInt_t j=0; j<settings->getGoodCellRegions3d().at(i).size(); j++)
                    if(cell==settings->getGoodCellRegions3d().at(i).at(j)){
//                        Int_t Entries = hLandauGoodCells->GetEntries();
//                        hLandauGoodCells->Add(hCellsLandau.at(cell));	//Not working for some reason, ask Felix
//                        hLandauGoodCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                        listGoodCells->Add(hCellsLandau.at(cell));
//                        cout<<"\tgood"<<endl;
                    }
            }
        }
    }
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
    histSaver->SaveHistogram(hLandauGoodCells);


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
        TString name = TString::Format("cDeadCellMeanCharge_%d",i);
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
    //	return;
	Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
	for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayMeanCharge]"<<endl;
		cout<<hCellsOverlayAvrgCharge.at(ClusterSize)<<endl;	if(hCellsOverlayAvrgCharge.at(ClusterSize)){
			cout<<hCellsOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
			cout<<hCellsOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
			TString name = "hCellsOverlayAvrgCharge_cl";
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgCharge.at(ClusterSize)->Clone(name);
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
		if(hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)){
			cout<<hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->IsZombie()<<endl;
			cout<<hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->GetEntries()<<endl;
			TString name = "hCellsOverlayAvrgChargeGoodCells_cl";
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeGoodCells.at(ClusterSize)->Clone(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->GetZaxis()->SetRangeUser(700,1200);
			histSaver->SaveHistogram(histo);
			delete histo;
		}
		if(hCellsOverlayAvrgChargeBadCells.at(ClusterSize)){
			cout<<hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->IsZombie()<<endl;
			cout<<hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->GetEntries()<<endl;
			TString name = "hCellsOverlayAvrgChargeBadCells_cl";
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeBadCells.at(ClusterSize)->Clone(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->GetZaxis()->SetRangeUser(300,1200);
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
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->GetZaxis()->SetRangeUser(700,1200);
			histSaver->SaveHistogram(histo);
			if(ClusterSize+1!=3)delete histo;
		}
		if(ClusterSize+1==3){
		    TProfile2D* prof = hCellsOverlayAvrgChargeMinusBadCells[ClusterSize];
		    DoMonteCarloOfAvrgChargePerBinInOverlay(prof,hCellsLandauMinusBadCells[ClusterSize]);
		    delete hCellsLandauMinusBadCells[ClusterSize];
		}
		if(hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)){
			TString name = "hCellsOverlayAvrgChargeNoColumnHit_cl";
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->Clone(name);
			histo->SetTitle("Avrg PH - overlayed - no hit in columns");
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->GetZaxis()->SetRangeUser(700,1200);
			histSaver->SaveHistogram(histo);
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
			histSaver->SaveHistogram(histo);
			delete histo;
			TH2D* hCellsOverlayAvrgChargeMinusBadCellsNentries = hCellsOverlayAvrgChargeMinusBadCells.at(ClusterSize)->ProjectionXY("hCellsOverlayAvrgChargeMinusBadCellsNentries","b");
			hCellsOverlayAvrgChargeMinusBadCellsNentries->Draw("colz");
			hCellsOverlayAvrgChargeMinusBadCellsNentries->GetZaxis()->SetTitle("number of entries");
			histSaver->SaveHistogram(hCellsOverlayAvrgChargeMinusBadCellsNentries,false,false);
			delete hCellsOverlayAvrgChargeMinusBadCellsNentries;
		}
		if(hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)){
			TString name = "hCellsOverlayAvrgChargeNoColumnHit_cl";
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOverlayAvrgChargeNoColumnHit.at(ClusterSize)->Clone(name);
			histo->SetTitle("Avrg PH - overlayed - no hit in columns");
			cout<<"SAVE"<<endl;
			histSaver->SaveHistogram(histo);
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
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsCentralColumnOverlayMeanCharge]"<<endl;
		cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)<<endl;
		if(hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)){
			cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
			cout<<hCellsCentralColumnOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
			TString name = "hCellsCentralColumnOverlayAvrgCharge_cl";
			name.Append(appendix);
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
			name.Append(appendix);
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
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCells.at(ClusterSize)->Clone(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
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
			name.Append(appendix);
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
			name.Append(appendix);
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
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->Clone(name);
			//TH2D* histo = (TProfile2D*)hCellsCentralColumnOverlayAvrgChargeMinusBadCellsBelowCut.at(ClusterSize)->ProjectionXY("B");
			//histo->SetName(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->Scale(1/settings->getOverlayColumnPulseHeightCut());
			/*histo->GetXaxis()->SetRangeUser(settings->getCentralColumnOverlayXLow(),settings->getCentralColumnOverlayXHigh());
			histo->GetYaxis()->SetRangeUser(settings->getCentralColumnOverlayYLow(),settings->getCentralColumnOverlayYHigh());*/
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

void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayOffsetMeanCharge() {
	//	return;

	Int_t MaxOverlayClusterSize = settings->getMaxOverlayClusterSize();
	for(Int_t ClusterSize = 1; ClusterSize<=MaxOverlayClusterSize; ClusterSize++){
		TString appendix = TString::Format("_ClusterSize%i", ClusterSize);
		if (settings->do3dTransparentAnalysis())
			appendix = TString::Format("_ClusterSize%i_trans", ClusterSize);

		cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayOffsetMeanCharge]"<<endl;
		cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)<<endl;
		if(hCellsOffsetOverlayAvrgCharge.at(ClusterSize)){
			cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->IsZombie()<<endl;
			cout<<hCellsOffsetOverlayAvrgCharge.at(ClusterSize)->GetEntries()<<endl;
			TString name = "hCellsOffsetOverlayAvrgCharge_cl";
			name.Append(appendix);
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
			name.Append(appendix);
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
			name.Append(appendix);
			TProfile2D* histo = (TProfile2D*)hCellsOffsetOverlayAvrgChargeGoodCells.at(ClusterSize)->Clone(name);
			cout<<"SAVE"<<endl;
			histo->Draw("goffcolz");
			histo->GetZaxis()->SetRangeUser(700,1200);
			histSaver->SaveHistogram(histo);
			delete histo;
		}

	} //End of for ClusterSize
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
            hLandauCluster1->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
            hFidCutXvsFidCutYClusters.at(4)->Fill(nfiducialValueX,nfiducialValueY,1);
        }
        if(diamondCluster0.getHighestSignalChannel()==55||diamondCluster1.getHighestSignalChannel()==55){
            hDoubleClusterPos1->Fill(diamondCluster0.getHighestSignalChannel());
            hDoubleClusterPos1->Fill(diamondCluster1.getHighestSignalChannel());
            hLandauCluster2->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
            hFidCutXvsFidCutYClusters.at(5)->Fill(nfiducialValueX,nfiducialValueY,1);
        }
        if((!diamondCluster0.getHighestSignalChannel()==55&&!diamondCluster1.getHighestSignalChannel()==55)||(!diamondCluster0.getHighestSignalChannel()==85&&!diamondCluster1.getHighestSignalChannel()==85)){
            hLandauDoubleCombined->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
        }

    }
    if(nClusters==3){
        hFidCutXvsFidCutYClusters.at(6)->Fill(nfiducialValueX,nfiducialValueY,1);
    }
}

int TAnalysisOf3dDiamonds::RemoveEdgeHits(TCluster* nCluster, pair<int,int> nDetector) {
    //Removes edge detector hits
    int Remove =0;
    for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
        if(nCluster->isHit(clPos)){
            if(nCluster->getChannel(clPos)>=nDetector.second||nCluster->getChannel(clPos)<=nDetector.first)
                Remove = 1;
        }
    }
    return Remove;
}

void TAnalysisOf3dDiamonds::createTreeTestHistos() {
    TH1F* histo = (TH1F*) clusteredAnalysis->getHistogram("hTest","pulseHeight","","");
    TH3F* histo2 = (TH3F*) clusteredAnalysis->getHistogram("hMeanPH","pulseHeight:nRow:nColumn","","");
    if(histo2){
        TH2F* hAvrgCharge = (TH2F*) histo2->Project3DProfile("yx");
        if(hAvrgCharge)hAvrgCharge->SetName("hAvrgChargeInCells");
        histSaver->SaveHistogram(hAvrgCharge);
        TCanvas *c1 = new TCanvas("cAvrgChargeInCells","cAvrgChargeInCells");
        c1->cd();
        hGridReferenceCellSpace->Draw("COL");
        hAvrgCharge->Draw("sameCOLZAH");
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
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    //Universal histograms

    TString name = "hValidEventsFiducialSpace";
    name.Append(appendix);
    hValidEventsFiducialSpace = new TH2F(name,name,1024,0,256,1024,0,256);
    hValidEventsFiducialSpace->GetXaxis()->SetTitle("avrg silicon position x /ch");
    hValidEventsFiducialSpace->GetYaxis()->SetTitle("avrg silicon position y /ch");
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
    hValidEventsDetSpace->GetXaxis()->SetTitle("predicted Hit position x in det space /#mum");
    hValidEventsDetSpace->GetYaxis()->SetTitle("predicted Hit position y in det space /#mum");
    InitialiseStripAnalysisHistos();

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
    histSaver->SaveHistogram(hValidEventsDetSpace);
    histSaver->SaveHistogram(hValidEventsFiducialSpace);
    SaveStripAnalysisHistos();
    // Save
    if(settings->do3dShortAnalysis() == 1){SaveShortAnalysisHistos();}
    if(settings->do3dLongAnalysis() == 1){SaveLongAnalysisHistos();}
    if(settings->do3dTransparentAnalysis() == 1){saveTransparentAnalysisHistos();/*saveTransparentAnalysisHistos()*/}
}

void TAnalysisOf3dDiamonds::initialiseLongAnalysisHistos() {
    initialise3DGridReference();
    initialise3DYAlignmentHistos();
    initialise3DOverviewHistos();
    initialise3D2DLandauAndClustersizeHistos();
    initialise3DCellOverlayHistos();
    initialise3DCellCentralColumnOverlayHistos();
    initialise3DOffsetOverlayHistos();
    initialiseEdgeFreeHistos();
    LongAnalysis_InitialiseRelativeAddedTransparentCharge();
    hLongAnalysisInvalidCellNo = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCellNo");
    hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCellNo");
    hLongAnalysisInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCluster");
    hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCluster");
}

void TAnalysisOf3dDiamonds::ShortAnalysis_SaveMeanChargeVector() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())

        cout<<"Create Project3dProfile for hChargeDistribution3D ..."<<flush;
    TProfile2D* hMeanCharge = histSaver->CreateProfile2D("hChargeDistribution3D",vecPredDetX_ShortAna,vecPredDetY_ShortAna,vecPulseHeight_ShortAna,1024,1024);
    if (!hMeanCharge){
        cerr<<" hChargDistribution3D: was not created:"<<endl;
        return;
    }
    hMeanCharge->GetXaxis()->SetTitle("Predicted Hit Position X in Detector /#mum");
    hMeanCharge->GetYaxis()->SetTitle("Predicted Hit Position Y in Detector /#mum");

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
    hMeanCharge->Draw("colz");
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
    settings->DrawMetallisationGrid(c1,3);
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
    delete c1;
    if (hMeanCharge)
        delete hMeanCharge;
}

void TAnalysisOf3dDiamonds::InitialiseStripAnalysisHistos() {
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    TString name = "hLandauStrip";
    name.Append(appendix);
    hLandauStrip =  new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauStrip->GetXaxis()->SetTitle("charge /ADC");
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

void TAnalysisOf3dDiamonds::SaveStripAnalysisHistos() {
    histSaver->SaveHistogram(hLandauStrip);
    histSaver->SaveHistogram(hLandauStripFidCutXvsFidCutY);
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveMeanChargePlots() {

    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";

    histSaver->SaveHistogram(hPulseHeightVsDetectorHitPostionXY);
    histSaver->SaveHistogramWithCellGrid(hPulseHeightVsDetectorHitPostionXY);
    histSaver->SaveHistogramWithCellGrid(hPulseHeightVsDetectorHitPostionXYGoodCells);

    UInt_t xBins = hPulseHeightVsDetectorHitPostionXY->GetXaxis()->GetNbins();
    UInt_t yBins = hPulseHeightVsDetectorHitPostionXY->GetYaxis()->GetNbins();
    UInt_t xRebin = xBins/settings->getNColumns3d()/2;
    UInt_t yRebin = yBins/settings->getNRows3d()/2;

    TString name = "hPulseHeightVsDetectorHitPostion_Quarters";
    name.Append(appendix);
    TH2F* hDetXvsDetY3DMeanChargeQuarters = (TH2F*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin,yRebin,name);
    histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeQuarters);
    histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeQuarters);
    delete hDetXvsDetY3DMeanChargeQuarters;
    //
    name = "hPulseHeightVsDetectorHitPostion_Cells";
    name.Append(appendix);
    TH2D* hDetXvsDetY3DMeanChargeCells = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin*2,yRebin*2,name);
    histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeCells);
    histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeCells);
    delete hDetXvsDetY3DMeanChargeCells;
    //
    delete hPulseHeightVsDetectorHitPostionXY;

}

void TAnalysisOf3dDiamonds::LongAnalysis_FillRelativeAddedTransparentCharge() {
    if(!settings->do3dTransparentAnalysis())
        return;
    pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
    UInt_t diamondPattern = settings->get3dWithHolesDiamondPattern();
    transparentCluster.SetTransparentClusterSize(transparentCluster.size());
    Float_t charge5 = transparentCluster.getCharge();
    Float_t charge = 0;
    Int_t realHitChannel = (Int_t)(transparentCluster.GetTransparentHitPosition()+.5);
    pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);
    bool isInEdgeRegion =  settings->IsOnTheEdgeOfCell(relPos.first,relPos.second);
    for(UInt_t i = 0; i < transparentCluster.size(); i++){
        UInt_t clusterSize = i+1;
        transparentCluster.SetTransparentClusterSize(clusterSize);
        Float_t addedCharge = - charge;
        charge = transparentCluster.getCharge();
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
        //			cout<<"\t"<<i<<"\t"<<transparentCluster.getCharge()<<"\t"<<transparentCluster.getCharge()/charge1<<endl;
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
		TString Histo1Title = TString::Format("Cluster Size %i",i-1);
		hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i-1]->SetName(Histo1Title);
		TString Histo2Title = TString::Format("Cluster Size %i",i);
		hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]->SetName(Histo2Title);

		TString name = TString::Format("hTransparentAnalysisTransparentChargeGoodCellsWithoutEdgeComparison%i_to_%i",i-1,i);
		histSaver->SaveTwoHistos((string)name,hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i-1],hTransparentAnalysisTransparentChargeGoodCellsWithoutEdge[i]);
	}

	for(int i=1; i<hTransparentAnalysisTransparentChargeBadCellsWithoutEdge.size(); i++){
		TString Histo1Title = TString::Format("Cluster Size %i",i);
		hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i-1]->SetName(Histo1Title);
		TString Histo2Title = TString::Format("Cluster Size %i",i+1);
		hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[i]->SetName(Histo2Title);

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
    Int_t bins = max-min;
    TH1F* hEntries = new TH1F(name,name,bins,min,max);
    for(Int_t binx = 1;binx < hBinContents->GetNbinsX();binx++)
    for(Int_t biny = 1;biny < hBinContents->GetNbinsY();biny++){
        Int_t bin = hBinContents->GetBin(binx,biny);
        hEntries->Fill(hBinContents->GetBinContent(bin));
    }
    histSaver->SaveHistogram(hEntries);
    name = "hMonteCarloAvrgChargePerBin";
    TH1D* hMonteCarloAvrgChargePerBin = new TH1D(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    for(UInt_t nMC = 0; nMC < 1e6; nMC++){
        Int_t entries = hEntries->GetRandom();
        Float_t mean = 0;
        for(Int_t entry =0; entry < entries; entry++){
            mean += hLandauOfOverlay->GetRandom();
        }
        mean /= (Float_t)entries;
        cout<<TString::Format("%6d --> %6.1f (%d)",nMC,mean,entries)<<endl;;
        hMonteCarloAvrgChargePerBin->Fill(mean);
    }
    hMonteCarloAvrgChargePerBin->GetXaxis()->SetTitle("mean charge per bin_{MonteCarlo}");
    hMonteCarloAvrgChargePerBin->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveHistogram(hMonteCarloAvrgChargePerBin);
    delete hMonteCarloAvrgChargePerBin;
    delete hBinContents;
    delete hEntries;

}
