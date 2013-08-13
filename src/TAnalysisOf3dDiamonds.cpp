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
	PulseHeightMin = 0;
	PulseHeightMax = 2800;
	PulseHeightMaxMeanCharge = 1500;

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
		StripAnalysis();
		if(settings->do3dShortAnalysis() == 1){ShortAnalysis();}
		if(settings->do3dLongAnalysis() == 1){LongAnalysis();}
		if(settings->do3dTransparentAnalysis() == 1){TransparentAnalysis();}
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
	if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
		return;
	int stripDetector = 1;
	if(settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY)!=stripDetector)
		return;

	if(eventReader->getNDiamondClusters()!=1)
		return;
	UInt_t diamond = TPlaneProperties::getDetDiamond();
	TCluster diamondCluster = eventReader->getCluster(diamond,0);
	int areaStripDetector = 0;
	if (!settings->isClusterInDiaDetectorArea(diamondCluster,areaStripDetector) ){
		return;
	}
	if( !settings->diamondPattern.isValidCluster(diamondCluster)){
		return;
	}
	if (diamondCluster.isSaturatedCluster())
		return;
	hLandauStrip->Fill(diamondCluster.getCharge());
}


void TAnalysisOf3dDiamonds::ShortAnalysis() {
	Float_t maxChi2 = settings->getChi2Cut3D();
	if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
		return;

	hNumberofClusters->Fill(eventReader->getNDiamondClusters());
	ClusterPlots(eventReader->getNDiamondClusters(),fiducialValueX,fiducialValueY);
	switch (eventReader->getNDiamondClusters()) {
	case 1:
		ShortAnalysis_Analyse1Cluster();
		break;
	case 2:
		ShortAnalysis_Analyse2Cluster();
		break;
	default:
		break;
	}
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
	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),clusterNo);
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

	//hFidCutXvsFidCutYvsCharge.at(1)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
	//RemoveLumpyClusters(&diamondCluster);

	//Edge Finding

	Float_t clusterCharge = diamondCluster.getCharge(false);
	ShortAnalysis_FillEdgeDistributions(clusterCharge);
	ShortAnalysis_FillMeanChargeVector(clusterCharge);
	//Universal PHvsChannel Plot
	for(UInt_t i=0; i < settings->diamondPattern.getNIntervals();i++){

		pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
		//cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;

		if((Int_t)diamondCluster.getHighestSignalChannel()<=channels.second&&(Int_t)diamondCluster.getHighestSignalChannel()>=channels.first){

			hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
			hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);
			if(!settings->getSelectionFidCuts()->getFidCut(i+1)->IsInFiducialCut(fiducialValueX,fiducialValueY))
				return;

			if(RemoveEdgeHits(&diamondCluster,channels)==1)		//If cluster has hit in edge channel remove.
				return;

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
			hHitandSeedCount[i]->Fill(HitCount,SeedCount);
			hChi2XChi2Y[i]->Fill(chi2x, chi2y);
			hFidCutXvsFidCutY[i]->Fill(fiducialValueX,fiducialValueY);
			//For hFidCutXvsFidCutYvsMeanCharge
			//hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);

		}
	}		//End of for diamond patterns
}


void TAnalysisOf3dDiamonds::ShortAnalysis_Analyse2Cluster(){
	TCluster cluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
	TCluster cluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
	Float_t pos1 = cluster1.getPosition();
	Float_t pos2 = cluster2.getPosition();
	Float_t ph1 = cluster1.getCharge();
	Float_t ph2 = cluster2.getCharge();
	if (!settings->diamondPattern.isValidChannelPosition(pos1)){
		ShortAnalysis_Analyse1Cluster(1);
		return;
	}
	if (!settings->diamondPattern.isValidChannelPosition(pos2)){
		ShortAnalysis_Analyse1Cluster(0);
		return;
	}
	if (pos1>128||pos1<0||pos2>128||pos2<0){
		cout<<"\n"<<nEvent<<" "<<pos1<<": "<<ph1<<"\t"<<pos2<<": "<<ph2<<endl;
		return;
	}
	vecPH_Cluster1_ShortAna.push_back(ph1);
	vecPH_Cluster2_ShortAna.push_back(ph2);
	vecCh_Cluster1_ShortAna.push_back(pos1);
	vecCh_Cluster2_ShortAna.push_back(pos2);
	Int_t pattern1 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos1));
	Int_t pattern2 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos2));
	Int_t predictedArea = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
	//	cout<<nEvent<<" "<<predictedArea<<" "<<pos1<<"-->"<<pattern1+1<<" "<<pos2<<"-->"<<pattern2+1<<endl;

}



void TAnalysisOf3dDiamonds::LongAnalysis() {

	Float_t maxChi2 = settings->getChi2Cut3D();
	if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
		return;

	if(eventReader->getNDiamondClusters()!=1)
		return;

	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());


	//	cout<<"xPredDet, yPredDet"<<xPredDet<<", "<<yPredDet<<endl;
	pair<int,int> cell = settings->getCellAndQuarterNo(xPredDet,yPredDet);
	//	cout<<cell.first<<" "<<cell.second<<endl;
	UInt_t cellNo = cell.first;
	//	cout<<"cellNo: "<<cellNo<<endl;
	UInt_t quarterNo = cell.second;
	//	if(cellNo>=0)cout<<"cell: "<<cellNo<<"--> row "<<row<<", column "<<column<<endl;

	Int_t area3DwithColumns = 2;
	if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
		hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
		return;
	}
	if( !settings->diamondPattern.isValidCluster(diamondCluster)){
		hLongAnalysisInvalidCluster->Fill(xPredDet,yPredDet);
		return;
	}
	if(!settings->isValidCellNo(cellNo)){
		hLongAnalysisInvalidCellNo->Fill(xPredDet,yPredDet);
		return;
	}

	if(diamondCluster.isSaturatedCluster())
		return;

	//	cout<<cellNo<<" "<< flush;
	pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);

	Float_t charge = diamondCluster.getCharge(false);
	hPulseHeightVsDetectorHitPostionXY->Fill(xPredDet,yPredDet,charge);
	//	hDetXvsDetY3DEvents->Fill(xPredDet,yPredDet,1);

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
	if(cellNo < hCellsClusteSize.size()){
		Int_t ClusterSize = diamondCluster.getClusterSize()-2;
		Int_t MaxHistoClusterSize = hCellsClusteSize.at(0)->GetNbinsX();
		if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
			hCellsClusteSize.at(cellNo)->Fill(MaxHistoClusterSize);
		else
			hCellsClusteSize.at(cellNo)->Fill(ClusterSize);
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
		if(quarterNo < hQuarterCellsLandau[cellNo].size()){
			Int_t ClusterSize = diamondCluster.size()-2;
			Int_t MaxHistoClusterSize = hQuarterCellsClusterSize[0].at(0)->GetNbinsX();
			if(ClusterSize>MaxHistoClusterSize)		//If ClusterSize is > hCellsClusterSize2D MaxClusterSize then goes into MaxClusterBin.
				hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(MaxHistoClusterSize);
			else
				hQuarterCellsClusterSize[cellNo][quarterNo]->Fill(ClusterSize);
		}
	}
	else{
		if(verbosity>6) cout <<TString::Format("X%2d-%d,\t%3.1f",cellNo,quarterNo,charge)<<endl;
	}
	LongAnalysis_FillOverlayedHistos(relPos.first,relPos.second,charge);
};

void TAnalysisOf3dDiamonds::TransparentAnalysis() {


	//Float_t maxChi2 = settings->getChi2Cut3D();
	if(chi2x>5||chi2y>20)     //(chi2x>maxChi2||chi2y>maxChi2)
		return;

	//	Float_t xPredDet = eventReader->getPositionInDetSystem(subjectDetector, xPredicted, yPredicted);
	//cout<<"YOffset: "<<settings->get3DYOffset()<<endl;

	//	vecXPredicted.push_back(xPredicted);
	//	vecYPredicted.push_back(yPredicted);
	//	vecXPreditedDetector.push_back(Xdet);
	//	vecYPreditedDetector.push_back(Ydet);

	//if(!settings->get3dMetallisationFidCuts()->isInFiducialCut(xPredDet,yPredDet))		// return if not in any of the metallisation regions
	//return;

	Int_t DiamondPattern;
	DiamondPattern = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPredDet,yPredDet);
	cout<<"DiamondPattern: "<<DiamondPattern<<endl;

	if (DiamondPattern!=3)
		return;
	vecChi2X.push_back(chi2x);
	vecChi2Y.push_back(chi2y);
	pair<int,int> channels = settings->diamondPattern.getPatternChannels(DiamondPattern);

	//Int_t XdetChannelSpaceInt = settings->convertMetricToChannelSpace(subjectDetector,Xdet)+0.5;
	//Float_t XdetChannelSpaceFloat = settings->convertMetricToChannelSpace(subjectDetector,Xdet);
	//cout<<"Diamond Pattern: "<<DiamondPattern<<endl;
	//cout<<"Hit coordinates: "<<"("<<Xdet<<", "<<Ydet<<")"<<endl;
	//cout<<"Channel Hit: "<<XdetChannelSpaceInt<<"   "<<XdetChannelSpaceFloat<<endl;
	Int_t XdetChannelSpaceInt =  settings->diamondPattern.convertMetricToIntChannel(xPredDet);
	Float_t XdetChannelSpace = settings->diamondPattern.convertMetricToChannel(xPredDet);
	UInt_t clusterSize =5;
	if(DiamondPattern == 1)
		clusterSize = 5;
	else{clusterSize = 3;}
	if(!TPlaneProperties::IsValidChannel(subjectDetector,XdetChannelSpaceInt))
		return;
	TCluster transparentCluster = TTransparentAnalysis::makeTransparentCluster(eventReader,settings,subjectDetector,XdetChannelSpace,clusterSize);
	transparentCluster.Print();
	Float_t TransparentCharge = getTransparentCharge(DiamondPattern, XdetChannelSpaceInt);
	cout<<TString::Format("%7d: Event in diamond pattern %d:\t %05.1f/%05.1f --> %3d --> %5.1f/%5.1f  ADC counts",
			nEvent,DiamondPattern,xPredDet,yPredDet,XdetChannelSpaceInt-channels.first,TransparentCharge,transparentCluster.getCharge())<<endl;

	//** WHY?
	if(TransparentCharge == -9999)
		return;
	if (DiamondPattern == 3)
		cout << TString::Format("%7d: %f - %f = %f\t%d\t%d",
				nEvent, TransparentCharge, transparentCluster.getCharge(true),TransparentCharge-transparentCluster.getCharge(true),
				DiamondPattern,settings->isMaskedCluster(subjectDetector,transparentCluster) ) << endl;

	hXdetvsYdetvsCharge[DiamondPattern-1]->Fill(xPredDet,yPredDet,TransparentCharge);
	hXdetvsYdetvsEvents[DiamondPattern-1]->Fill(xPredDet,yPredDet,1);
	hLandauTransparent[DiamondPattern-1]->Fill(TransparentCharge);

	if(DiamondPattern ==2 || DiamondPattern ==3){
		if(settings->isBadCell(DiamondPattern,xPredDet,yPredDet) == 1){
			//cout<<settings->getCellNo(Xdet,Ydet).first<<endl;
			return;
		}
		hLandauTransparentBadCellsRemoved.at(DiamondPattern-1)->Fill(TransparentCharge);
	}

	//For overlay of 3DwH
	if(DiamondPattern == 3){
		pair<Float_t,Float_t> relPos = settings->getRelativePositionInCell(xPredDet,yPredDet);

		hCellOverlayvsCharge->Fill(relPos.first,relPos.second,TransparentCharge);
		hCellOverlayvsEvents->Fill(relPos.first,relPos.second,1);
	}

	/*
	//For Transparent Analysis.
	Int_t HitCell;
	float hEntries0 = 0;
	float TransparentCharge = 0;
	float TransparentChargeAddition = 0;
	int SaturatedEvent = 0;
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){

				for(int k=1;k<3;k++){	//First memeber of array is the number of channels to readout.
					if(eventReader->isSaturated(subjectDetector,CellToChannel(HitCell,Xcell)[k]))
						SaturatedEvent = 1;
					//					cout<<"Hit channel: "<<CellToChannel(HitCell,Xcell)[k]<<endl;
					//TransparentCharge = eventReader->getAdcValue(subjectDetector,HitChannel);
					TransparentChargeAddition = TransparentCharge;
					TransparentCharge = eventReader->getSignal(subjectDetector,CellToChannel(HitCell,Xcell)[k]) + TransparentChargeAddition;
				}	//End of for HitChannels
				if(SaturatedEvent == 0){	// Only plot events with no saturated channels.
					//					cout<<"Transparent charge is: "<<TransparentCharge<<endl;
					hTransparentCharge3D->Fill(TransparentCharge);
					hCellTransparentLandau.at(i*11+j)->Fill(TransparentCharge);
					if(TransparentCharge<700)
						hCellsTransparentHitPosition.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
					//cout<<eventReader->getSignal(subjectDetector,HitChannel)<<endl;
				}
			}
		}
	}

		//Universal PHvsChannel Plot
		for(int i=0; i<settings->diamondPattern.getNIntervals();i++){      //settings->diamondPattern.getNIntervals(); i++){

			pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
			//cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;

			if(diamondCluster.getHighestSignalChannel()<=channels.second&&diamondCluster.getHighestSignalChannel()>=channels.first){

				hEventsvsChannel.at(i)->Fill(diamondCluster.getHighestSignalChannel());
				hPHvsChannel.at(i)->Fill(diamondCluster.getCharge(false),diamondCluster.getHighestSignalChannel());
				//hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);
				//hPHvsPredictedXPos.at(i)->Fill(diamondCluster.getCharge(false),xPredicted);
				hLandau.at(i)->Fill(diamondCluster.getCharge(false));
				vecPHDiamondHit.at(i)->push_back(diamondCluster.getCharge(false));
				//vecXPredicted.at(i)->push_back(xPredicted);
				//vecYPredicted.at(i)->push_back(yPredicted);
				hHitandSeedCount.at(i)->Fill(HitCount,SeedCount);
				//hChi2XChi2Y.at(i)->Fill(chi2x, chi2y);
				hFidCutXvsFidCutY.at(i)->Fill(fiducialValueX,fiducialValueY);
				//For hFidCutXvsFidCutYvsMeanCharge
				//hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);

			}
		}		//End of for diamond patterns
	//eventReader->getAdcValue(subjectDetector,CellToChannel());

	//Int_t getAdcValue(UInt_t det,UInt_t ch);

	//CellToChannel();
	 *
	 */
}

void TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos() {
	cout<<"[TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos()]"<<endl;
	//Universal histograms

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

	cout<<"loop over patterns"<<endl;
	for(UInt_t i=0; i<settings->diamondPattern.getNIntervals(); i++){
		pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
		//hLandau
		TString name = TString::Format("hLandau_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
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
		hEventsvsChannel.push_back(new TH1F(name,name,100,0,100));
		if(hEventsvsChannel.back()){
			hEventsvsChannel.back()->GetXaxis()->SetTitle("HighestPH [ch]");
			hEventsvsChannel.back()->GetYaxis()->SetTitle("No. Events");
		}
		//hPHvsChannel
		name = TString::Format("hPHvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
		hPHvsChannel.push_back(new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,100,0,100));
		if(hPHvsChannel.back()){
			hPHvsChannel.back()->GetXaxis()->SetTitle("Charge in ADC counts");
			hPHvsChannel.back()->GetYaxis()->SetTitle("XPos(channel)");
			hPHvsChannel.back()->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
		}
		//hPHvsChannel.at(i)->SetMaximum(3000);
		//hPHvsChannel.at(i)->SetMinimum(0);

		//hHitandSeedCount
		name = TString::Format("hHitandSeedCount_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
		hHitandSeedCount.push_back(new TH2F(name,name,10,-.5,9.5,10,-.5,9.5));
		hHitandSeedCount.back()->GetXaxis()->SetTitle("Hit Count");
		hHitandSeedCount.back()->GetYaxis()->SetTitle("Seed Count");

		//hChi2XChi2Y
		name = TString::Format("hChi2X_vs_chi2Y_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
		hChi2XChi2Y.push_back(new TH2F(name,name,60,0,30,60,0,30));
		hChi2XChi2Y.back()->GetXaxis()->SetTitle("#chi^{2}_{X}");
		hChi2XChi2Y.back()->GetYaxis()->SetTitle("#chi^{2}_{Y}");

		//hFidCutXvsFidCutY
		name = TString::Format("hFidCutXvsFidCutY_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
		hFidCutXvsFidCutY.push_back(new TH2F(name,name,160,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),120,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
		hFidCutXvsFidCutY.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutY.at(i)->GetYaxis()->SetTitle("FidCutY");

		//hFidCutXvsFidCutYvsCharge		For TH2D

		name = TString::Format("hFidCutXvsFidCutYvsCharge_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
		hFidCutXvsFidCutYvsCharge.push_back(new TH2D(name,name,
				213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
				160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
		hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

		//hFidCutXvsFidCutYvsEvents
		hFidCutXvsFidCutYvsEvents.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsEvents"));

		//hFidCutXvsFidCutYvsMeanCharge
		hFidCutXvsFidCutYvsMeanCharge.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsMeanCharge"));

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
	}

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
	//hDetXvsDetY3DTotolCharge
	TString name = "hPulseHeightVsDetectorHitPostionXY";
	if(verbosity>1) cout<<"Create "<<name<<endl;
	UInt_t factor = 10* (settings->getNQuarters3d()/2);
	hPulseHeightVsDetectorHitPostionXY = histSaver->GetProfile2dBinedInCells(name,factor);
	hPulseHeightVsDetectorHitPostionXY->GetXaxis()->SetTitle("predicted x position in detector /#mum");
	hPulseHeightVsDetectorHitPostionXY->GetYaxis()->SetTitle("predicted y position in detector /#mum");
	hPulseHeightVsDetectorHitPostionXY->GetZaxis()->SetTitle("charge /ADC");

	//hDetXvsDetY3DEvents
	//	hDetXvsDetY3DEvents = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Clone("hDetXvsDetY3DvsEvents");

	//hDetXvsDetY3DMeanCharge


	//hDetXvsDetY3DRebinnedMeanChargeRMS
	name = "h3DdetRebinnedRMS";
	hDetXvsDetY3DRebinnedRMS = histSaver->GetHistoBinedInCells(name);
	//			new TH2D(hDetXvsDetY3DRebinnedRMSName.str().c_str(),hDetXvsDetY3DRebinnedRMSName.str().c_str(),
	//			settings->getNColumns3d(),getXMetalisationStart3d,getXMetalisationEnd3d,
	//			settings->getNRows3d(),getYMetalisationStart3d,getYMetalisationEnd3d);
	hDetXvsDetY3DRebinnedRMS->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinnedRMS->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinnedRMS->GetZaxis()->SetTitle("Charge ADC");

	//hBinnedMeanCharge
	name = "h3DdetCellMeanChargeBinned";
	if(verbosity>1) cout<<"Create "<<name<<endl;
	hBinnedMeanCharge = new TH1F(name,name,9,400,1300);
	hBinnedMeanCharge->GetXaxis()->SetTitle("MeanCharge");
	hBinnedMeanCharge->GetYaxis()->SetTitle("Entries");

	//hDetXvsDetY3DOverview
	name = "hDetXvsDetY3DOverview";
	hDetXvsDetY3DOverview = histSaver->GetHistoBinedInCells(name);
	hDetXvsDetY3DOverview->GetXaxis()->SetTitle("Xdet (#mum)");
	hDetXvsDetY3DOverview->GetYaxis()->SetTitle("Ydet (#mum)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();



	//hCellsMeanClusteSize
	name = "hCellsMeanClusteSize";
	hCellsMeanClusteSize = histSaver->GetHistoBinedInCells(name);

	//hQuarterCellsMeanClusteSize
	name = "hQuarterCellsMeanClusterSize";
	hQuarterCellsMeanClusterSize = hCellsMeanClusteSize = histSaver->GetHistoBinedInQuarters(name);

	//RebinnedQuarterCellFails
	name = "3DdetNumberofQuarterCellFails";
	RebinnedQuarterCellFails = histSaver->GetHistoBinedInCells(name);
	RebinnedQuarterCellFails->GetXaxis()->SetTitle("Xdet (um)");
	RebinnedQuarterCellFails->GetYaxis()->SetTitle("Ydet (um)");
	RebinnedQuarterCellFails->GetZaxis()->SetTitle("Quarter Fails");

	//hDetXvsDetY3DQuarterCellGrading
	for(int k=0; k<6; k++){
		name = TString::Format("hDetXvsDetY3DMeanChargeQuarterCellGrading_%d_Fail",k);
		hDetXvsDetY3DMeanChargeQuarterCellGrading.push_back(histSaver->GetHistoBinedInQuarters(name));
		hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetXaxis()->SetTitle("Xdet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetYaxis()->SetTitle("Ydet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.back()->GetZaxis()->SetTitle("Charge ADC");
	}

	//h3DdetQuarterCellFluctuation
	name = "h3DdetQuarterCellFluctuation";
	h3DdetQuarterCellFluctuation = histSaver->GetHistoBinedInQuarters(name);
	h3DdetQuarterCellFluctuation->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetQuarterCellFluctuation->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetQuarterCellFluctuation->GetZaxis()->SetTitle("Fluctuation");

	//h3DdetQuarterCellFluctuation1
	name = "h3DdetQuarterCellFluctuation1";
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
		if(verbosity>1) cout<<"Create "<<name<<endl;
		hQuarterCellGradedTransparentLandau.push_back(new TH1F(name,name,256,0,2800));

		name = TString::Format("hDetXvsDetY3DMeanChargeQuarterCellGradingLandau_Grade%d",k);
		if(verbosity>1)cout<<"Create "<<name<<endl;
		hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.push_back(new TH1F(name,name,256,0,2800));

		//hCellsDeltaXQuarterCellGrading
		name = TString::Format("hCellsDeltaXQuarterCellGrading_Grade%d",k);
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
		if(verbosity>1)cout<<"Create "<<name<<endl;
		hCellsLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
		name = "hCellsClusteSize";
		if(verbosity>1) cout<<"Create "<<name<<endl;
		hCellsClusteSize.push_back(new TH1F(name,name,6,0,6));

		//		hQuarterCellsLandau.push_back(vector<TH1F*>());
		hQuarterCellsLandau[cell].resize(settings->getNQuarters3d());
		for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){
			name = TString::Format("hQuaterCellsLandau_%d_%d",cell,quarter);
			if(verbosity>1)cout<<"Create "<<name<<endl;
			hQuarterCellsLandau[cell][quarter] =
					new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
			name = TString::Format("hQuarterCellsClusterSize_%d_%d",cell,quarter);
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
		//@iain: what are this hardcoded values? is it relative hit in cell?
		hCellsTransparentHitPositionCellGraded.push_back(new TH2D(name,name,30,0,150,30,0,150));
	}

	//hTransparentCharge3D
	name = "hTransparentCharge3D";
	hTransparentCharge3D = new TH1F(name,name,256,0,2800);
	hTransparentCharge3D->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hTransparentCharge3D->GetYaxis()->SetTitle("Entries");

	//hCellsLandauGraded    &&    hCellsLandauGradedNoColumn
	for(int i=0;i<12;i++){	//Group CellsLandaus within same ranges together. 0-100; 100-200; -> 1100-1200;
		name = TString::Format("hLandauCellsGraded_%d_to_%d",(i*100),(i+1)*100);
		hCellsLandauGraded.push_back(new TH1F(name,name,256,0,2800));
		name = TString::Format("hLandauCellsGradedNoColumn_%d_to_%d",i*100,(i+1)*100);
		hCellsLandauGradedNoColumn.push_back(new TH1F(name,name,256,0,2800));
	}

	//To create Strip detector Landau
	pair<int,int> channels =settings->diamondPattern.getPatternChannels(1);
	//hLandau
	name = TString::Format("hLandau_ch_%d_to_%d",channels.first, channels.second);
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

void TAnalysisOf3dDiamonds::initialise3DCellOverlayHistos() {

	//Cell Overlay
	//hCellsOverlayed
	TString name = "hCellsOverlayAvrgCharge";
	hCellsOverlayAvrgCharge = new TProfile2D(name,name,
			30,0,settings->GetCellWidth(subjectDetector,2),
			30,0,settings->GetCellHeight());
	hCellsOverlayAvrgCharge->GetXaxis()->SetTitle("rel. x position within a cell /#mum");
	hCellsOverlayAvrgCharge->GetYaxis()->SetTitle("rel. y position within a cell /#mum");
	hCellsOverlayAvrgCharge->GetZaxis()->SetTitle("pulse height of cluster /ADC");
	//	hCellsOverlayAvrgCharge->SetContour(99);
	name = "hCellsOverlayAvrgCharge_noColumnHit";
	hCellsOverlayAvrgChargeNoColumnHit = (TProfile2D*)hCellsOverlayAvrgCharge->Clone(name);
	hCellsOverlayAvrgChargeNoColumnHit->SetTitle(name);

	cout<<hCellsOverlayAvrgChargeNoColumnHit<<" "<<hCellsOverlayAvrgCharge<<endl;
	cout<<" "<<hCellsOverlayAvrgCharge->IsZombie()<<endl;

	//hCellsOverlayedColumnLandau
	name = "hCellOverlayWithColumnLandau";
	hCellOverlayWithColumnLandau = new TH1F(name,name,256,0,2800);//todo iain ph
	hCellOverlayWithColumnLandau->GetXaxis()->SetTitle("charge / ADC");
	hCellOverlayWithColumnLandau->GetYaxis()->SetTitle("Entries");

	//hCellsOverlayedLandauNoColumn
	name = "hCellOverlayNoColumnLandau";
	hCellOverlayNoColumnLandau = new TH1F(name,name,256,0,2800);//todo iain ph

	name = "hCellOverlayColumnLandau";
	hCellOverlayColumnLandau = new TH1F(name,name,256,0,2800);//todo iain ph


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
	}
	 */

	cout<<"End initialise3DCellOverlayHistos()"<<endl;
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
	ShortAnalysis_SaveMeanChargeVector();
	ShortAnalysis_Save2ClusterPlots();
	vector<Float_t> xPred;
	vector<Float_t> yPred;
	vector<Float_t> charge;

	for(int i = 0; i < vecEdgePredX.size(); i++){
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
		TString name = "hEdgeFittingAvrgCharge_";
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
		TString name = "c_";
		name.Append(hLandau[i]->GetName());

		Float_t factor = hLandau[i]->GetBinContent(hLandau[i]->GetMaximumBin());
		factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
		name.Append("_normalized");
		histSaver->SaveTwoHistosNormalized((string)name,hLandau[i],hLandauStrip);

		Float_t max = hHitandSeedCount[i]->GetBinContent(hHitandSeedCount[i]->GetMaximumBin());
		hHitandSeedCount[i]->Scale(1./max);
		histSaver->SaveHistogram(hHitandSeedCount[i]);
		histSaver->SaveHistogram(hPHvsChannel.at(i));
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
		TString hName  = TString::Format("cFidCutXvsFidCutYvsMeanCharge_%d_%d",channels.first,channels.second);
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
	TString name = "hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn";
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
	cCombinedMeanCharge->SetName(name);
	histSaver->SaveCanvas(cCombinedMeanCharge);

	// Selection Fiducial Cuts Drawn
	name = "hFidCutXvsFidCutYvsMeanChargeAllDetectors";
	settings->getSelectionFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge);
	cCombinedMeanCharge->SetName(name);
	histSaver->SaveCanvas(cCombinedMeanCharge);

	// Edge F
	cCombinedMeanCharge->Clear();
	name = "hFidCutXvsFidCutYvsMeanChargeAllEdges";
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
	settings->get3dEdgeFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge,true);
	cCombinedMeanCharge->SetName(name);
	histSaver->SaveCanvas(cCombinedMeanCharge);

	for ( UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts(); i++){
		cCombinedMeanCharge->Clear();
		hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
		settings->get3dEdgeFidCuts()->getFidCut(i+1)->DrawFiducialCutToCanvas(cCombinedMeanCharge,true);
		TString name = "hFidCutXvsFidCutYvsMeanCharge_";
		name.Append(settings->getEdgePositionName(i));
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

	/*//hNumberofClusters
	histSaver->SaveHistogram(hNumberofClusters);
	for(int i=0;i<settings->diamondPattern.getNIntervals();i++){
		hEventsvsChannelCombined->Add(hEventsvsChannel.at(i));
	}
	histSaver->SaveHistogram(hEventsvsChannelCombined);
	 */

	for(int i=0;i<settings->diamondPattern.getNIntervals();i++){
		/*//hLandau
		stringstream hLandauName; hLandauName<<"hLandau%%"<<firstCh<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hLandau.push_back(HistogrammSaver::CreateDistributionHisto(hLandauName.str().c_str(),*vecPHDiamondHit.at(i),256));
		hLandau.at(i)->SetTitle(hLandauName.str().c_str());
		hLandau.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
		hLandau.at(i)->GetYaxis()->SetTitle("number of entries #");
		 */
		pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
		TString name = "c_";
		name.Append(hLandauTransparent[i]->GetName());
		Float_t factor = hLandauTransparent[i]->GetBinContent(hLandau[i]->GetMaximumBin());
		factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
		histSaver->SaveTwoHistos((string)name,hLandauTransparent[i],hLandauStrip,factor);
		name = "c_";
		name.Append(hLandauTransparentBadCellsRemoved[i]->GetName());
		histSaver->SaveTwoHistos((string)name,hLandauTransparentBadCellsRemoved[i],hLandauStrip);

		//histSaver->SaveHistogram(hPHvsChannel.at(i));
		//histSaver->SaveHistogram(hPHvsPredictedXPos.at(i));
		//histSaver->SaveHistogram(hChi2XChi2Y.at(i));
		//histSaver->SaveHistogram(hFidCutXvsFidCutY.at(i));
		//histSaver->SaveHistogram(hPHvsPredictedChannel.at(i));
		//histSaver->SaveHistogram(hFidCutXvsFidCutYvsCharge.at(i));

		/*//hFidCutXvsFidCutYvsMeanCharge
		ptrCanvasMean.at(i)->cd();
		 *hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
		hFidCutXvsFidCutYvsMeanCharge.at(i)->SetEntries(hFidCutXvsFidCutYvsEvents.at(i)->Integral());
		hFidCutXvsFidCutYvsMeanCharge.at(i)->Draw("COLZ");
		TString hName  = TString::Format("cFidCutXvsFidCutYvsMeanCharge_%d_%d",channels.first,channels.second);
		ptrCanvasMean.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvasMean[i]);
		 */

		//hXdetvsYdetvsEvents
		Int_t DiamondPattern = i+1;
		ptrCanvasXdetvsYdetMeanCharge.push_back(new TCanvas());
		ptrCanvasXdetvsYdetMeanCharge.at(i)->cd();
		*hXdetvsYdetvsMeanCharge.at(i) = (*hXdetvsYdetvsCharge.at(i)/(*hXdetvsYdetvsEvents.at(i)));
		hXdetvsYdetvsMeanCharge.at(i)->SetEntries(hXdetvsYdetvsEvents.at(i)->Integral());
		hXdetvsYdetvsMeanCharge.at(i)->Draw("COLZ");
		settings->DrawMetallisationGrid(ptrCanvasXdetvsYdetMeanCharge.at(i), DiamondPattern);
		TString hName  = TString::Format("cXdetvsYdetMeanCharge_%d_%d",channels.first,channels.second);
		ptrCanvasXdetvsYdetMeanCharge.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvasXdetvsYdetMeanCharge[i]);

	} //End of for loop

	cCellOverlayvsMeanCharge = new TCanvas();
	cCellOverlayvsMeanCharge->cd();
	*hCellOverlayvsMeanCharge = (*hCellOverlayvsCharge/(*hCellOverlayvsEvents));
	hCellOverlayvsMeanCharge->SetEntries(hCellOverlayvsEvents->Integral());
	hCellOverlayvsMeanCharge->Draw("COLZ");
	cCellOverlayvsMeanCharge->SetName("cCellOverlayvsMeanChargeTransparent");
	histSaver->SaveCanvas(cCellOverlayvsMeanCharge);


	/*//For all Diamond hFidCutXvsFidCutYvsMeanCharge
	hCombinedMeanCharge = new TCanvas();
	hCombinedMeanCharge->cd();
	//FidCudBoundMetric->Draw("same"); //To draw the fiducial cut regions
	hCombinedMeanCharge->SetName("hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn");
	histSaver->SaveCanvas(hCombinedMeanCharge);
	settings->get3dFidCuts()->drawFiducialCutsToCanvas(hCombinedMeanCharge);
	hCombinedMeanCharge->SetName("hFidCutXvsFidCutYvsMeanChargeAllDetectors");
	histSaver->SaveCanvas(hCombinedMeanCharge);

	for(int i=0;i<7;i++){
		stringstream cClustersName; cClustersName<<"cFidCutXvsFidCutYClusters"<<i;
		cClusters.push_back(new TCanvas(cClustersName.str().c_str(),cClustersName.str().c_str()));
		cClusters.at(i)->cd();
		hFidCutXvsFidCutYClusters.at(i)->Draw("COLZ");
		histSaver->SaveCanvas(cClusters.at(i));
	}
	 */

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
	LongAnalysisSaveCellAndQuaterNumbering();
	LongAnalysis_SaveGoodAndBadCellLandaus();
	LongAnalysis_SaveDeadCellProfile();
	LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors();
	LongAnalysis_SaveFailedQuarters();
	LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail();
	LongAnalysis_SaveCellsOverlayMeanCharge();
	LongAnalysis_SaveMeanChargePlots();
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
		//Float_t Fluctuation = TMath::Abs((hLandauGoodCellsMean - QuarterMean)/MeanOfLandauGoodCells);
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
	TH2D* histo = histSaver->GetHistoBinedInCells(name,2);
	cout<<"LongAnalysis_CreateQuarterCellsPassFailAndCellGradingVectors"<<endl;
	cout<<" Criteria: "<<quarterFailCriteriaTyp<<endl;
	cout<<" Mean of Landau - Good Cells: "<<MeanOfLandauGoodCells<<endl;
	cout<<" Mean of Landau - Strip:      "<<hLandauStrip->GetMean()<<endl;

	vecQuarterCellsPassFail.clear();
	vecQuarterCellsPassFail.resize(settings->GetNCells3d());
	vector<vector<Float_t> > vecQuarterCellsFluctuation;
	vecQuarterCellsFluctuation.resize(settings->GetNCells3d());;
	Float_t DeadCellThreshold = 700;
	for (UInt_t cell = 0; cell < settings->GetNCells3d();cell++){
		Int_t Grading = 0;
		vecQuarterCellsFluctuation.at(cell) = LongAnalysis_GradeCellByQuarters(quarterFailCriteriaTyp, hQuarterCellsLandau[cell]);
		cout<<vecQuarterCellsFluctuation[cell].size()<<endl;
		if(LongAnalysis_IsDeadCell(hQuarterCellsLandau[cell], DeadCellThreshold))
			CellGrading.push_back(hQuarterCellsLandau[cell].size());
		else{
			for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){
				Float_t Fluctuation = vecQuarterCellsFluctuation[cell][quarter];
				Float_t FluctuationFactor = 0.1; //Needs to be added to settings file.
				bool isFailedQuarter = TMath::Abs(Fluctuation)>FluctuationFactor;
				cout<<TString::Format("%d - %d: %2.1f ==> %d",cell,quarter,100.*Fluctuation,isFailedQuarter);
				//vecQuarterCellsFluctuation[cell].push_back(Fluctuation);
				if(isFailedQuarter)
					vecQuarterCellsPassFail[cell].push_back(1);
				else vecQuarterCellsPassFail[cell].push_back(0);
				Grading += vecQuarterCellsPassFail[cell].at(quarter);
			}
			cout<<"Cell: "<<cell<<"  Grading: "<<Grading<<endl;
			CellGrading.push_back(Grading);
		}
	}
}

bool TAnalysisOf3dDiamonds::LongAnalysis_IsDeadCell(vector<TH1F*> nhQuarterCellsLandau, Float_t nThreshold){

	Int_t NQuarterLow =0;
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
}


void TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsLandau2DHighlightedQuarterFail() {

	Int_t yBins = hCellsLandau.at(0)->GetNbinsX();
	Float_t yMin = hCellsLandau.at(0)->GetBinLowEdge(1);
	Float_t yMax = hCellsLandau.at(0)->GetBinLowEdge(yBins) + hCellsLandau.at(0)->GetBinWidth(yBins);

	TString name = "hCellsLandau2DHighlightedQuarterFail";
	TH2D* hCellsLandau2DHighlightedQuarterFail = new TH2D(name,name,settings->GetNCells3d(),0,settings->GetNCells3d(),yBins,yMin,yMax);
	hCellsLandau2DHighlightedQuarterFail->GetXaxis()->SetTitle("Cell");
	hCellsLandau2DHighlightedQuarterFail->GetYaxis()->SetTitle("Charge ADC");
	TH2D* hCellsLandau2DHighlightedQuarterFailSorted = (TH2D*)hCellsLandau2DHighlightedQuarterFail->Clone("hCellsLandau2DHighlightedQuarterFailSorted");
	vector<TH1F*> hCellLandausSorted = hCellsLandau;
	sort(hCellLandausSorted.begin(), hCellLandausSorted.end(), TSettings::SorterForPulseHeightOfHisto);

	TH2D* hHighlightedQuarterFail = new TH2D("hHighlightedQuarterFail","",settings->GetNCells3d(),0,settings->GetNCells3d(),1,yMin,yMax);
	TH2D* hHighlightedQuarterFailSorted = (TH2D*) hHighlightedQuarterFail->Clone("hHighlightedQuarterFailSorted");
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

	TCanvas* cCellsLandau2DHighlightedQuarterFail = new TCanvas("cCellsLandau2DHighlightedQuarterFail","cCellsLandau2DHighlightedQuarterFail");
	cCellsLandau2DHighlightedQuarterFail->cd();
	//hCellsLandau2DHighlightedQuarterFail->SetEntries(hCellsLandau2DEntries);
	hCellsLandau2DHighlightedQuarterFail->SetStats(kFALSE);
	hHighlightedQuarterFail->SetStats(kFALSE);
	hCellsLandau2DHighlightedQuarterFail->Draw("COLZ");
	hHighlightedQuarterFailSorted->Draw("COLAHsame");
	hCellsLandau2DHighlightedQuarterFailSorted->Draw("sameCOLZ");
	histSaver->SaveCanvas(cCellsLandau2DHighlightedQuarterFail);
	delete cCellsLandau2DHighlightedQuarterFail;

	cCellsLandau2DHighlightedQuarterFail = new TCanvas("cCellsLandau2DHighlightedQuarterFailSorted","cCellsLandau2DHighlightedQuarterFail - sorted");
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

void TAnalysisOf3dDiamonds::LongAnalysis_FillOverlayedHistos(Float_t xRelPosDet,Float_t yRelPosDet,
		Float_t clusterCharge) {
	//	hCellsOverlayCharge->Fill(xRelPosDet,yRelPosDet,clusterCharge);
	//	hCellsOverlayEvents->Fill(xRelPosDet,yRelPosDet,1);
	hCellsOverlayAvrgCharge->Fill(xRelPosDet,yRelPosDet,clusterCharge);

	hCellOverlayWithColumnLandau->Fill(clusterCharge);
	if(settings->IsWithInTheColumnRadius(xRelPosDet, yRelPosDet)){
		hCellOverlayColumnLandau->Fill(clusterCharge);
	}
	else{
		hCellsOverlayAvrgChargeNoColumnHit->Fill(xRelPosDet,yRelPosDet,clusterCharge);
		hCellOverlayNoColumnLandau->Fill(clusterCharge);
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

	TH1F* hLandauBadCells = (TH1F*)hCellsLandau.at(0)->Clone("hLandauBadCells");
	hLandauBadCells->SetTitle("Landau of bad cells");
	TH1F* hLandauGoodCells = (TH1F*)hCellsLandau.at(0)->Clone("hLandauGoodCells");
	hLandauGoodCells->SetTitle("Landau of good cells");

	for(UInt_t column=0;column<settings->getNColumns3d();column++){
		for(UInt_t row=0;row<settings->getNRows3d();row++){
			Int_t cell = settings->get3DCellNo((int)column,row);
			//hCellNumbering->SetBinContent(column+1,row+1,cell); //This should be a clone of the 2D Cell Mean Charge Plot, Wait till Felix has finished.
			//histSaver->SaveHistogram(hCellsLandau.at(cell));

			for(UInt_t i=0; i<settings->getBadCells3D().size(); i++)
				if(cell==settings->getBadCells3D().at(i)){
					Int_t Entries = hLandauBadCells->GetEntries();
					hLandauBadCells->Add(hCellsLandau.at(cell));	//Not working for some reason, ask Felix
					hLandauBadCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
				}
			for(UInt_t i=0; i<settings->getGoodCells3D().size(); i++)
				if(cell==settings->getGoodCells3D().at(i)){
					Int_t Entries = hLandauGoodCells->GetEntries();
					hLandauGoodCells->Add(hCellsLandau.at(cell));	//Not working for some reason, ask Felix
					hLandauGoodCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
				}
		}
	}
	Float_t factor = hLandauGoodCells->GetBinContent(hLandauGoodCells->GetMaximumBin());
	factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
	histSaver->SaveTwoHistos("hLandauGoodCells",hLandauGoodCells,hLandauStrip,factor);
	histSaver->SaveTwoHistosNormalized("hLandauGoodCellsNormalized",hLandauGoodCells,hLandauStrip);

	factor = hLandauBadCells->GetBinContent(hLandauBadCells->GetMaximumBin());
	factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
	histSaver->SaveTwoHistos("cLandauBadCells",hLandauBadCells,hLandauStrip,factor);
	//Histogram(hLandauBadCells);
	histSaver->SaveTwoHistosNormalized("cLandauBadCellsNormalized",hLandauBadCells,hLandauStrip);

	MeanOfLandauGoodCells = hLandauGoodCells->GetMean();
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
	TCanvas *c1;
	for(UInt_t i=0; i<settings->getDeadCell3D().size(); i++){
		TString name = TString::Format("cDeadCellMeanCharge_%d",i);
		c1 = new TCanvas(name,name);
		c1->cd();
		if(hDeadCellCharge[i])
			hDeadCellCharge[i]->Draw("");
		else
			cerr<<TString::Format("hDeadCellCharge[%d] invalid ", i)<<endl;
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
	cout<<"[TAnalysisOf3dDiamonds::LongAnalysis_SaveCellsOverlayMeanCharge]"<<endl;
	cout<<hCellsOverlayAvrgCharge<<endl;
	if(hCellsOverlayAvrgCharge){
		cout<<hCellsOverlayAvrgCharge->IsZombie()<<endl;
		cout<<hCellsOverlayAvrgCharge->GetEntries()<<endl;
		TString name = "hCellsOverlayAvrgCharge_cl";
		TH2F* histo = (TH2F*)hCellsOverlayAvrgCharge->Clone(name);
		cout<<"SAVE"<<endl;
		histSaver->SaveHistogram(histo);
		delete histo;
		////		hCellsOverlayAvrgCharge->SetName("hCellsOverlayAvrgCharge");
		//		cout<<"Set Name: "<<hCellsOverlayAvrgCharge<<endl;
		////		hCellsOverlayAvrgCharge->SetTitle("Avrg PH - overlayed");
		//		histSaver->SaveHistogram(hCellsOverlayAvrgCharge);
		//		delete hCellsOverlayAvrgCharge ;
	}
	if(hCellsOverlayAvrgChargeNoColumnHit){
		TString name = "hCellsOverlayAvrgChargeNoColumnHit_cl";
		TH2F* histo = (TH2F*)hCellsOverlayAvrgChargeNoColumnHit->Clone(name);
		histo->SetTitle("Avrg PH - overlayed - no hit in columns");
		cout<<"SAVE"<<endl;
		histSaver->SaveHistogram(histo);
		delete histo;
		//		hCellsOverlayAvrgChargeNoColumnHit->SetName("hCellsOverlayAvrgChargeNoColumns");
		//		hCellsOverlayAvrgChargeNoColumnHit->SetTitle("Avrg PH - overlayed - no hit in columns");
		//		histSaver->SaveHistogram(hCellsOverlayAvrgChargeNoColumnHit);
		//		delete hCellsOverlayAvrgChargeNoColumnHit;
	}
	//		hCellsOverlayPulseHeight->Project3D("xy");

	/*hCellsOverlayEvents->Draw("sameTEXT");
		cCellsOverlayMeanCharge->SetName("cCellsOverlayMeanChargeWithEntries");
		histSaver->SaveCanvas(cCellsOverlayMeanCharge);*/

}

void TAnalysisOf3dDiamonds::SaveLongAnalysisHistos2() {




	/*//hCellNumbering
	cCellNumbering = new TCanvas("cCellNumbering","cCellNumbering");
	cCellNumbering->cd();
	hGridReferenceDetSpace->SetTitle("hCellNumbering");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hCellNumbering->Draw("sameCOLZAH");
	hCellNumbering->Draw("sameTEXTAH");
	settings->DrawMetallisationGrid(cCellNumbering, 3);
	histSaver->SaveCanvas(cCellNumbering);
	 */

	//For h3DdetMeanCharge
	//	TString name = "c3DdetMeanCharge";
	//	cDetXvsDetY3DMeanCharge = new TCanvas(name,name);
	//	cDetXvsDetY3DMeanCharge->cd();

	/*
	//hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanCharge");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hGridReference->Draw("COL");
	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanCharge, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanCharge);
//*/
	//hCellsLandau2DHighlightedQuarterFail
	/*	Int_t Ncells = settings->getNColumns3d()*settings->getNRows3d();
	xBins = hCellsLandau.at(0)->GetNbinsX();
	Float_t xMin = hCellsLandau.at(0)->GetBinLowEdge(1);
	Float_t xMax = hCellsLandau.at(0)->GetBinLowEdge(xBins) + hCellsLandau.at(0)->GetBinWidth(xBins);*/
	//RebinnedMeanCharge
	//	cDetXvsDetY3DMeanChargeRebinned = new TCanvas("c3DdetMeanChargeRebinned","c3DdetMeanChargeRebinned");
	//	cDetXvsDetY3DMeanChargeRebinned->cd();
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	//	 *hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DChargeRebinned/(*hDetXvsDetY3DEventsRebinned));
	//	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DEventsRebinned->Integral());
	//	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	//	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanChargeRebinned");		//Set title to require
	//	hGridReferenceDetSpace->Draw("COL");
	//	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	//	hDetXvsDetY3DEventsRebinned->Draw("sameTEXTAH");
	//	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	//	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeRebinned, 3);
	//	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeRebinned);
	/*c3DdetEventsRebinned = new TCanvas("c3DdetEventsRebinned","c3DdetEventsRebinned");
	//c3DdetEventsRebinned->cd();
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");			//("COLZ");
	//settings->DrawMetallisationGrid();
	histSaver->SaveCanvas(c3DdetEventsRebinned);
//	 */

	//cout<<"Ncells: "<<Ncells<<" xBins: "<<xBins<<" xMin: "<<xMin<<" xMax: "<<xMax<<endl;


	/*

	for(int column=0;column<settings->getNColumns3d();column++){
		for(int row=0;row<settings->getNRows3d();row++){

			Int_t cell = column*settings->getNRows3d() + row;

			stringstream binLabel; binLabel<<cell;
//			hCellsLandau2DHighlightedQuarterFail->GetYaxis()->SetBinLabel((cell),binLabel.str().c_str());
//			hCellsLandau2DHighlightedQuarterFail->GetYaxis()->SetLabelSize(0.02);
//			for(int xBin =0; xBin<xBins; xBin++)
//				hCellsLandau2DHighlightedQuarterFail->SetBinContent(xBin,cell,hCellsLandau.at(cell)->GetBinContent(xBin));
//			hHighlightedQuarterFail->SetBinContent(1,cell,hCellGrading.at(cell));
		}
	}
	TCanvas* cCellsLandau2DHighlightedQuarterFail = new TCanvas("cCellsLandau2DHighlightedQuarterFail","cCellsLandau2DHighlightedQuarterFail");
	cCellsLandau2DHighlightedQuarterFail->SetCanvasSize(1500,3000);
	cCellsLandau2DHighlightedQuarterFail->cd();
	//hCellsLandau2DHighlightedQuarterFail->SetEntries(hCellsLandau2DEntries);
	 */
	//	hCellsLandau2DHighlightedQuarterFail->SetStats(kFALSE);
	//	hHighlightedQuarterFail->SetStats(kFALSE);
	//	hCellsLandau2DHighlightedQuarterFail->Draw("COLZ");
	//	//hHighlightedQuarterFail->Draw("COLAH");
	//	//hCellsLandau2DHighlightedQuarterFail->Draw("sameCOLZ");
	//	histSaver->SaveCanvas(cCellsLandau2DHighlightedQuarterFail);
	//
	//	int QuaterCellEntrySum = 0;
	//	int QuaterCellEntrySumAddition = 0;
	//	for(int column=0;column<settings->getNColumns3d();column++){
	//		for(int row=0;row<settings->getNRows3d();row++){
	//			for(int quarter=0;quarter<settings->getNQuarters3d();quarter++){	//For each quater cell
	//				Float_t xBin = 2*column + (quarter/2)+1;
	//				Float_t yBin = 2*row +quarter+1-2*(quarter/2);
	//				Int_t entry = column*settings->getNRows3d()*4+row*4+quarter;
	//				if (entry<hQuaterCellsLandau.size()){
	//					Float_t mean = hQuaterCellsLandau.at(entry)->GetMean();
	//					Float_t rms = hQuaterCellsLandau.at(entry)->GetRMS();
	//					if(hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetNbinsX()>xBin && hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetNbinsY()>yBin)
	//						hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetBinContent(xBin,yBin,mean);
	//					else
	//						cout<<" hDetXvsDetY3DMeanChargeRebinnedQuarterCell, bin does not exist: " << xBin<<"/"<<yBin<<endl;
	//					if(hDetXvsDetY3DRebinnedQuarterCellRMS->GetNbinsX()>xBin && hDetXvsDetY3DRebinnedQuarterCellRMS->GetNbinsY()>yBin)
	//						hDetXvsDetY3DRebinnedQuarterCellRMS->SetBinContent(xBin, yBin, rms);
	//					else
	//						cout<<" hDetXvsDetY3DRebinnedQuarterCellRMS, bin does not exist: " << xBin<<"/"<<yBin<<endl;
	//					QuaterCellEntrySumAddition = QuaterCellEntrySum;
	//					QuaterCellEntrySum = QuaterCellEntrySumAddition + hQuaterCellsLandau.at(entry)->GetEntries();
	//					//histSaver->SaveHistogram(hQuaterCellsLandau.at(entry));
	//					if(entry<hQuarterCellsClusterSize.size()){
	//						cout<<"save: hQuarterCellsClusterSize["<<entry<<"] / "<<hQuarterCellsClusterSize.size()<<endl;
	//						cout<<"hQuarterCellsClusterSize: "<<flush<<" "<<hQuarterCellsClusterSize[entry]<<endl;
	//						//histSaver->SaveHistogram(hQuarterCellsClusterSize.at(entry));
	//					}
	//				}
	//				else{
	//					cout<<"size of hQuaterCellsLandau is to small"<<entry<<"/"<<hQuaterCellsLandau.size()<<endl;
	//				}
	//			}
	//		}
	//	}
	//
	//	cout<<"done"<<endl;
	//	cout<<"cDetXvsDetY3DMeanChargeRebinnedQuarterCell"<<endl;
	//	TCanvas* cDetXvsDetY3DMeanChargeRebinnedQuarterCell = new TCanvas("c3DdetQuarterCellMeanCharge","c3DdetQuarterCellMeanCharge");
	//	cDetXvsDetY3DMeanChargeRebinnedQuarterCell->cd();
	//	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetEntries(QuaterCellEntrySum);
	//	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetStats(kFALSE);
	//	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellMeanCharge");		//Set title to require
	//	hGridReferenceDetSpace->Draw("COL");
	//	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameCOLZAH");
	//	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameTEXTAH");
	//	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeRebinnedQuarterCell, 3);
	//	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeRebinnedQuarterCell);

	//Cell Overlay

	/*//cout<<"CellArraySize is: "<<hCellsCharge.size()<<endl;
	//for(int i=0;i<hCellsCharge.size();i++){
	for(int column=0;column<settings->getNColumns3d();column++){		//run over each cell
		for(int row=0;row<settings->getNRows3d();row++){
			//histSaver->SaveHistogram(hCellsLandau.at(i*settings->getNRows3d()+j));	//Save each cells Landau;
			LandauGaussFit landauGauss;
			//TF1* fit = landauGauss.doLandauGaussFit(hCellsDeltaX.at(i*settings->getNRows3d()+j));
			//histSaver->SaveHistogram(hCellsDeltaX.at(i*settings->getNRows3d()+j));		//Save each cells X res in Ch.

			int Continue = 1;
			for(int badCell=0;badCell<settings->getBadCells3D().size();badCell++){
				//if((i*settings->getNRows3d()+j)==CellsHarris10Bad[k]){
				if((column*settings->getNRows3d()+row)==settings->getBadCells3D().at(badCell)){
					Continue = 0;
					//					int Entries = hCellsHarris10Bad->GetEntries();;
					//					hCellsHarris10Bad->Add(hCellsLandau.at(column*settings->getNRows3d()+row));
					//					hCellsHarris10Bad->SetEntries((Entries+hCellsLandau.at(column*settings->getNRows3d()+row)->GetEntries()));

					RebinnedQuarterCellFails->SetBinContent(column+1,row+1,4);

					for(int l=0;l<4;l++){		//To fill highlighted grading, plot 5.
						cout<<"h2DMeanClusterSizeQuarterCell"<<column<<" "<<row<<" "<<badCell<<" "<<l<<endl;
						//h2DMeanClusterSizeQuarterCell
						Float_t xBin = 4*4+SortArrayPointer[badCell];
						//						Float_t yBin = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+badCell)->GetMean();
						//						Float_t zBin = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+badCell)->GetMean();
						//						h2DMeanClusterSizeQuarterCell->Fill(xBin,yBin,zBin);
						//h2DMeanClusterSizeQuarterCellTotal->Fill(4*4+SortArrayPointer[k],hQuaterCellsLandau.at(i*settings->getNRows3d()*4+j*4+k)->GetMean(),1);
						//h2DMeanClusterSizeQuarterCellEvents->Fill(4*4+SortArrayPointer[k],1);

						//NumberQuarterFails++;
						if(badCell<2)			//To make it more obvious which quarters have failed.
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*column+1),(2*row+1+badCell),500);
						if(badCell>1)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*column+2),(2*row+1+badCell-2),500);
					}
					//For Transparent Analysis
					hQuarterCellGradedTransparentLandau.at(4)->Add(hCellTransparentLandau.at(column*settings->getNRows3d()+row));

					if(settings->getBadCells3D().at(badCell)>10&&settings->getBadCells3D().at(badCell)<88){		//Only non edge cells used
						//hCellsDeltaXQuarterCellGrading
						hCellsDeltaXQuarterCellGrading.at(4)->Add(hCellsDeltaX.at(column*settings->getNRows3d()+row));
					}
					int Events =0;
					int ContentSum =0;
					for(int ll=2;ll<7;ll++){
						if(ll==6){
							int EventsSum = 0;
							for(int m=6;m<20;m++){
								EventsSum = Events;
								Events = hCellsClusteSize.at(column*settings->getNRows3d()+row)->GetBinContent(m+1)+EventsSum;
							}
						}
						else{
							Events = hCellsClusteSize.at(column*settings->getNRows3d()+row)->GetBinContent(ll);
						}
						//cout<<"Cell: "<<i*settings->getNRows3d()*4+j*4+SortArrayPointer[l]<<"Quarter Fails: "<<NumberQuarterFails<<"ClusterSize: "<<ll-1<<"Events: "<<Events<<endl;
						if(Events>0){
							//printArray(QuarterMeanCharge, 4, "-");
							//printArray(QuarterMeanChargeUnSorted, 4, "-");
							//printArray(SortArrayPointer, 4, "-");

						}
						ContentSum = h2DClusterSize->GetBinContent(4+1,ll-1); //,Events); //4 quarters failed
						h2DClusterSize->SetBinContent(4+1,ll-1,Events+ContentSum);
						for(int m=0;m<5;m++){
							ContentSum = h2DClusterSizeClone->GetBinContent(4+1,m+1); //,Events);
							h2DClusterSizeClone->SetBinContent(4+1,m+1,Events+ContentSum);
						}
					}
					for(int k=0;k<4;k++){	//Run over 4 quarters
						int Events =0;
						int ContentSum =0;
						for(int l=2;l<7;l++){
							if(l==6){
								int EventsSum = 0;
								for(int m=6;m<20;m++){
									EventsSum = Events;
									//									Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+k)->GetBinContent(m+1)+EventsSum;
								}
							}
							else{
								//								Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+k)->GetBinContent(l);
							}
							ContentSum = h2DClusterSizeQuarterCell->GetBinContent(4*4+k+1,l-1);	//To fill 4 failed quarters area of plot
							h2DClusterSizeQuarterCell->SetBinContent(4*4+k+1,l-1,Events+ContentSum);
							for(int m=0;m<5;m++){
								ContentSum = h2DClusterSizeQuarterCellClone->GetBinContent(4*4+k+1,m+1); //,Events);
								h2DClusterSizeQuarterCellClone->SetBinContent(4*4+k+1,m+1,Events+ContentSum);
								//cout<<"Here"<<endl;
							}
						}
					}	//End of running over quarters.
				}
			}
			if(Continue == 1){		//If a bad cell, the code won't continue.

				//hCellsHarrisGoodCells
				for(int k=0;k<18;k++){
					if((column*settings->getNRows3d()+row)==settings->getGoodCells3D().at(k)){
						//						int Entries = hCellsHarris18Good->GetEntries();
						//						hCellsHarris18Good->Add(hCellsLandau.at(column*settings->getNRows3d()+row));
						//						hCellsHarris18Good->SetEntries((Entries+hCellsLandau.at(column*settings->getNRows3d()+row)->GetEntries()));
					}
				}

				int NumberQuarterFails = 0;
				float QuarterMeanCharge[4];
				float QuarterMeanChargeUnSorted[4];
				for(int k=0;k<4;k++){
					SortArrayPointer[k] = k;
					//					QuarterMeanCharge[k] = hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+k)->GetMean();
					//					QuarterMeanChargeUnSorted[k] = hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+k)->GetMean();
				}
				SortArrayBtoS(QuarterMeanCharge, 4);
				float FlucFail = 0.1;		//If fluctuation greater than this -> Quarter Fails
				float OneFail = (QuarterMeanCharge[3]-QuarterMeanCharge[0])/QuarterMeanCharge[3];
				float TwoFail = (QuarterMeanCharge[3]-QuarterMeanCharge[1])/QuarterMeanCharge[3];
				float ThreeFail = (QuarterMeanCharge[3]-QuarterMeanCharge[2])/QuarterMeanCharge[3];
				if(ThreeFail>FlucFail)
					NumberQuarterFails = 3;
				else{
					if(TwoFail>FlucFail)
						NumberQuarterFails = 2;
					else{
						if(OneFail>FlucFail)
							NumberQuarterFails = 1;
					}
				}
				RebinnedQuarterCellFails->SetBinContent(column+1,row+1,NumberQuarterFails);
				if(NumberQuarterFails ==0){
					//hCellsTransparentHitPositionCellGraded.at(NumberQuarterFails)->Add(hCellsTransparentHitPosition.at(i*settings->getNRows3d()+j));
					int Events =0;
					int ContentSum =0;
					for(int ll=2;ll<7;ll++){
						if(ll==6){
							int EventsSum = 0;
							for(int m=6;m<20;m++){
								EventsSum = Events;
								Events = hCellsClusteSize.at(i*settings->getNRows3d()+j)->GetBinContent(m+1)+EventsSum;
							}
						}
						else{
							Events = hCellsClusteSize.at(i*settings->getNRows3d()+j)->GetBinContent(ll);
						}
						//cout<<"Cell: "<<i*settings->getNRows3d()*4+j*4+SortArrayPointer[l]<<"Quarter Fails: "<<NumberQuarterFails<<"ClusterSize: "<<ll-1<<"Events: "<<Events<<endl;
						if(Events>0){
							//printArray(QuarterMeanCharge, 4, "-");
							//printArray(QuarterMeanChargeUnSorted, 4, "-");
							//printArray(SortArrayPointer, 4, "-");

						}
						ContentSum = h2DClusterSize->GetBinContent(NumberQuarterFails+1,ll-1,Events);
						h2DClusterSize->SetBinContent(NumberQuarterFails+1,ll-1,Events+ContentSum);
					}

					for(int k=0;k<4;k++){
						int Events =0;
						int ContentSum =0;
						for(int l=2;l<7;l++){
							if(l==6){
								int EventsSum = 0;
								for(int m=6;m<20;m++){
									EventsSum = Events;
									Events = hQuarterCellsClusterSize.at(i*settings->getNRows3d()*4+j*4+k)->GetBinContent(m+1)+EventsSum;
								}
							}
							else{
								Events = hQuarterCellsClusterSize.at(i*settings->getNRows3d()*4+j*4+k)->GetBinContent(l);
							}
							ContentSum = h2DClusterSizeQuarterCell->GetBinContent(k+1,l-1,Events);
							h2DClusterSizeQuarterCell->SetBinContent(k+1,l-1,Events+ContentSum);
						}
					}	//End of running over quarters.
				}
	 *//*
				for(int k=0;k<4;k++){		//To fill fluctuation plot.
					//NumberQuarterFails++;
					float Fluctuation = 0;
					float Fluctuation1 = 0;
					if(QuarterMeanCharge[3]-QuarterMeanCharge[k] == 0)
						Fluctuation = 0;
					else{
						Fluctuation = (QuarterMeanCharge[3]-QuarterMeanCharge[k])/QuarterMeanCharge[3];
					}
					Fluctuation1 = (hCellsHarris18Good->GetMean()-QuarterMeanCharge[k])/hCellsHarris18Good->GetMean();

					//cout<<QuarterMeanCharge[3]<<" "<<QuarterMeanCharge[k]<<" "<<Fluctuation<<endl;
					if(SortArrayPointer[k]<2){			//To make it more obvious which quarters have failed.
						h3DdetQuarterCellFluctuation->SetBinContent((2*column+1),(2*row+1+SortArrayPointer[k]),Fluctuation);
						h3DdetQuarterCellFluctuation1->SetBinContent((2*column+1),(2*row+1+SortArrayPointer[k]),Fluctuation1);
					}
					if(SortArrayPointer[k]>1){
						h3DdetQuarterCellFluctuation->SetBinContent((2*column+2),(2*row+1+SortArrayPointer[k]-2),Fluctuation);
						h3DdetQuarterCellFluctuation1->SetBinContent((2*column+2),(2*row+1+SortArrayPointer[k]-2),Fluctuation1);
					}
				}
				for(int k=0;k<NumberQuarterFails;k++){		//To fill highlighted grading, plot 5.
						//NumberQuarterFails++;
					if(SortArrayPointer[k]<2)			//To make it more obvious which quarters have failed.
						hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*column+1),(2*row+1+SortArrayPointer[k]),500);
					if(SortArrayPointer[k]>1)
						hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*column+2),(2*row+1+SortArrayPointer[k]-2),500);
				}
					else{
						if(k<2)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+1),(2*j+1+k),1000);
						if(k>1)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+2),(2*j+1+k-2),1000);
					}
				}

				for(int k=0;k<4;k++){	//For each Grading
					if(k==NumberQuarterFails){
						//For Transparent Analysis
						hCellsTransparentHitPositionCellGraded.at(NumberQuarterFails)->Add(hCellsTransparentHitPosition.at(column*settings->getNRows3d()+row));
						hQuarterCellGradedTransparentLandau.at(k)->Add(hCellTransparentLandau.at(column*settings->getNRows3d()+row));
						hCellsDeltaXQuarterCellGrading.at(k)->Add(hCellsDeltaX.at(column*settings->getNRows3d()+row));

						int Events =0;
						int ContentSum =0;
						for(int ll=2;ll<7;ll++){
							if(ll==6){
								int EventsSum = 0;
								for(int m=6;m<20;m++){
									EventsSum = Events;
									Events = hCellsClusteSize.at(column*settings->getNRows3d()+row)->GetBinContent(m+1)+EventsSum;
								}
							}
							else{
								Events = hCellsClusteSize.at(column*settings->getNRows3d()+row)->GetBinContent(ll);
							}
							//cout<<"Cell: "<<i*settings->getNRows3d()*4+j*4+SortArrayPointer[l]<<"Quarter Fails: "<<NumberQuarterFails<<"ClusterSize: "<<ll-1<<"Events: "<<Events<<endl;
							if(Events>0){
								//printArray(QuarterMeanCharge, 4, "-");
								//printArray(QuarterMeanChargeUnSorted, 4, "-");
								//printArray(SortArrayPointer, 4, "-");

							}
							ContentSum = h2DClusterSize->GetBinContent(NumberQuarterFails+1,ll-1); //,Events);
							h2DClusterSize->SetBinContent(NumberQuarterFails+1,ll-1,Events+ContentSum);
							for(int m=0;m<5;m++){
								ContentSum = h2DClusterSizeClone->GetBinContent(NumberQuarterFails+1,m+1); //,Events);
								h2DClusterSizeClone->SetBinContent(NumberQuarterFails+1,m+1,Events+ContentSum);
							}
						}

						for(int l=0;l<4;l++){	//To fill each quarter cell of the graded cell.

							h2DMeanClusterSizeQuarterCell->Fill(NumberQuarterFails*4+SortArrayPointer[l],hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+l)->GetMean(),hQuarterCellsClusterSize.at(column*11*4+row*4+k)->GetMean());
							//h2DMeanClusterSizeQuarterCellTotal->Fill(NumberQuarterFails*4+SortArrayPointer[l],hQuaterCellsLandau.at(i*settings->getNRows3d()*4+j*4+l)->GetMean());
							//h2DMeanClusterSizeQuarterCellEvents->Fill(NumberQuarterFails*4+SortArrayPointer[l],1);

							hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.at(k)->Add(hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+l));
							if(l<2){
								hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->SetBinContent((2*column+1),(2*row+1+l),hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+l)->GetMean());
							}
							if(l>1){
								hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->SetBinContent((2*column+2),(2*row+1+l-2),hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+l)->GetMean());
							}

							int Events =0;
							int ContentSum =0;
							for(int ll=2;ll<7;ll++){
								if(ll==6){
									int EventsSum = 0;
									//cout<<"Mean Charge Array: "<<endl;
									//printArray(QuarterMeanChargeUnSorted, 4, "-");
									//cout<<"Pointer to Array: "<<endl;
									//printArray(SortArrayPointer, 4, "-");
									for(int m=6;m<20;m++){
										EventsSum = Events;
										Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+SortArrayPointer[l])->GetBinContent(m+1)+EventsSum;
									}
								}
								else{
									Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+SortArrayPointer[l])->GetBinContent(ll);
								}
								//cout<<"Cell: "<<i*settings->getNRows3d()*4+j*4+SortArrayPointer[l]<<"Quarter Fails: "<<NumberQuarterFails<<"ClusterSize: "<<ll-1<<"Events: "<<Events<<endl;
								if(Events>0){
									//printArray(QuarterMeanCharge, 4, "-");
									//printArray(QuarterMeanChargeUnSorted, 4, "-");
									//printArray(SortArrayPointer, 4, "-");

								}
								ContentSum = h2DClusterSizeQuarterCell->GetBinContent(NumberQuarterFails*4+SortArrayPointer[l]+1,ll-1);
								h2DClusterSizeQuarterCell->SetBinContent(NumberQuarterFails*4+SortArrayPointer[l]+1,ll-1,Events+ContentSum);
								for(int m=0;m<5;m++){
									ContentSum = h2DClusterSizeQuarterCellClone->GetBinContent(NumberQuarterFails*4+SortArrayPointer[l]+1,m+1); //,Events);
									h2DClusterSizeQuarterCellClone->SetBinContent(NumberQuarterFails*4+SortArrayPointer[l]+1,m+1,Events+ContentSum);
									//cout<<"Here"<<endl;
								}
							}
						}	//End of running over number of quarters
					}
				}

				//if(i!=BadCells[k]){
				if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1)>1000)&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1)<1100)){
					//cout<<"Cell with average charge > 1000: "<<(i*settings->getNRows3d()+j)<<endl;
					//if(NumberQuarterFails<3){
					hCellsOverlayedCharge->Add(hCellsCharge.at(column*settings->getNRows3d()+row));
					hCellsOverlayedEvents->Add(hCellsEvents.at(column*settings->getNRows3d()+row));
					hCellsOverlayedEventsNoColumns->Add(hCellsEventsNoColumn.at(column*settings->getNRows3d()+row));
					hCellsOverlayedLandauNoColumn->Add(hCellsLandauNoColumn.at(column*settings->getNRows3d()+row));
					//cout<<"Cell: "<<(i*settings->getNRows3d()+j)<<" is overlayed."<<endl;

					hDetXvsDetY3DOverview->SetBinContent(column+1,row+1,1);
					//cout<<"Bin content is: "<<hDetXvsDetY3DOverview->GetBinContent(1,1)<<endl;
					//histSaver->SaveHistogram(hCellsCharge.at(i));
					//gStyle->SetOptStat("irm");
					hCellsEvents.at(column*settings->getNRows3d()+row)->SetEntries(hCellsEvents.at(column*settings->getNRows3d()+row)->Integral()); //Set the number of entries to the integral of hEvents
					//histSaver->SaveHistogram(hCellsEvents.at(i*settings->getNRows3d()+j));
					//cout<<"The integral of Events"<<(i*settings->getNRows3d()+j)<<" is: "<<hCellsEvents.at(i*settings->getNRows3d()+j)->Integral()<<endl;
				}
				for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
					if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1)>(k*100))&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1)<((k+1)*100))){
						hCellsLandauGraded.at(k)->Add(hCellsLandau.at(column*settings->getNRows3d()+row));
						hCellsLandauGradedNoColumn.at(k)->Add(hCellsLandauNoColumn.at(column*settings->getNRows3d()+row));
					}
				}
			}	//End of if Continue == 1.
		//}
		}	//End of for j
	}	//End of for i
	//To draw RMS of each bin
	for(int i=0;i<30;i++){
		for(int j=0;j<30;j++){
			hCellsOverlayed55RMS->SetBinContent(i+1,j+1,hCellsOverlayBinSpec55.at(i*30+j)->GetRMS());
		}
	}
	cCellsOverlayed = new TCanvas("cCellsOverlayedMeanCharge","cCellsOverlayedMeanCharge");
	cCellsOverlayed->cd();
	  *hCellsOverlayedMeanCharge = (*hCellsOverlayedCharge/(*hCellsOverlayedEvents));
	hCellsOverlayedMeanCharge->SetEntries(hCellsOverlayedEvents->Integral());

	hCellsOverlayedMeanCharge->GetZaxis()->SetRangeUser(800,1200);
	hCellsOverlayedMeanCharge->SetContour(99);
	hCellsOverlayedMeanCharge->Draw("COLZ");
	histSaver->SaveCanvas(cCellsOverlayed);
	//hCellsOverlayed55RMS->Draw("sameTEXT");

	hCellsOverlayedEvents->Draw("sameTEXT");
	 cCellsOverlayed->SetName("cCellsOverlayedMeanChargeWithEntries");
	 histSaver->SaveCanvas(cCellsOverlayed);
	histSaver->SaveHistogram(hBinnedMeanCharge);	//*//*

	//To rebin and replot Overlayed plot.
	int minus0X[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
	int minus5X[] = {30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
	int plus5X[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1};

	int minus0Y[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
	int minus5Y[] = {30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
	int plus5Y[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1};

	int* ptrAxisX[3];
	int* ptrAxisY[3];

	ptrAxisX[0] = minus0X; ptrAxisX[1] = minus5X; ptrAxisX[2] = plus5X;
	ptrAxisY[0] = minus0Y; ptrAxisY[1] = minus5Y; ptrAxisY[2] = plus5Y;

	int temp;
	for(int i=0;i<3;i++){		//increment ptrAxisX
		for(int ii=0;ii<3;ii++){		//increment ptrAxisY
			for(int j=0;j<15;j++){		//Cells in x
				for(int k=0;k<15;k++){		//Cells in y
					for(int l=0;l<4;l++){		//Four subcells.
						if(l<2){
							temp = hCellsOverlayedChargeBinAlignment1.at(i*3+ii)->GetBinContent(j+1,k+1);
							hCellsOverlayedChargeBinAlignment1.at(i*3+ii)->SetBinContent(j+1,k+1,(temp+hCellsOverlayedCharge->GetBinContent(ptrAxisX[i][j*2],ptrAxisY[ii][k*2+l])));
							//cout<<ptrAxisX[i][j*2]<<" "<<ptrAxisY[ii][k*2+l]<<endl;
							temp = hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->GetBinContent(j+1,k+1);
							hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->SetBinContent(j+1,k+1,(temp+hCellsOverlayedEvents->GetBinContent(ptrAxisX[i][j*2],ptrAxisY[ii][k*2+l])));
						}
						if(l>1){
							temp = hCellsOverlayedChargeBinAlignment1.at(i*3+ii)->GetBinContent(j+1,k+1);
							hCellsOverlayedChargeBinAlignment1.at(i*3+ii)->SetBinContent(j+1,k+1,(temp+hCellsOverlayedCharge->GetBinContent(ptrAxisX[i][j*2+1],ptrAxisY[ii][k*2+l-2])));
							//cout<<ptrAxisX[i][i*2+1]<<" "<<ptrAxisY[ii][j*2+l-2]<<endl;
							temp = hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->GetBinContent(j+1,k+1);
							hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->SetBinContent(j+1,k+1,(temp+hCellsOverlayedEvents->GetBinContent(ptrAxisX[i][j*2+1],ptrAxisY[ii][k*2+l-2])));
						}
					}	//subcells
				}	//y cells
			}	//x cells
			float Significance = 0;
			if((i*3+ii)==0){
				for(int j=0;j<15;j++){
					for(int k=0;k<15;k++){
						float OverlayedMean = hCellsOverlayedLandauNoColumn->GetMean();
						float BinMean = hCellsOverlayBinSpec1010.at(j*15+k)->GetMean();
						float BinRMS = hCellsOverlayBinSpec1010.at(j*15+k)->GetRMS();
						Significance = (OverlayedMean-BinMean)/BinRMS;
						hCellsOverlayed1010RMS->SetBinContent(j+1,k+1,hCellsOverlayBinSpec1010.at(j*15+k)->GetRMS());
						hCellsOverlayed1010Significance->SetBinContent(j+1,k+1,Significance);
						//histSaver->SaveHistogram(hCellsOverlayBinSpec1010.at(j*15+k));
					}
				}
			}	//end of if first plot
			TString hName = TString::Format("hCellsOverlayedMeanChargeBinAlignment_%03d",i*3+ii);
			cCellsOverlayedBinAlignment.push_back(new TCanvas(hName,hName));
			cCellsOverlayedBinAlignment.at(i*3+ii)->cd();
	 *hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii) = (*hCellsOverlayedChargeBinAlignment1.at(i*3+ii)/(*hCellsOverlayedEventsBinAlignment1.at(i*3+ii)));
			hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii)->SetEntries(hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->Integral());
			hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii)->Draw("COLZ");
			if((i*3+ii)==0)
				hCellsOverlayed1010Significance->Draw("sameTEXT");
			histSaver->SaveCanvas(cCellsOverlayedBinAlignment.at(i*3+ii));
			//hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->Draw("sameTEXT");
		}	//increment y
	}	//increment x

	//hCellsOverlayedEventsNoColumns
	for(int i=0;i<15;i++){
		for(int j=0;j<15;j++){
			hCellsOverlayedEntriesNoColumns->Fill(hCellsOverlayedEventsNoColumns->GetBinContent(i+1,j+1));
		}
	}
	LandauGaussFit landauGauss;
	histSaver->SaveHistogram(hCellsOverlayedEntriesNoColumns);
	TF1* fit00 = landauGauss.doLandauGaussFit(hCellsOverlayedLandauNoColumn);
	histSaver->SaveHistogram(hCellsOverlayedLandauNoColumn);

	//hCellsHarrisandAlexGoodandBadCells
	//TF1* fit01 = landauGauss.doLandauGaussFit(hLandau.at(0));
	pair<int,int> channels = settings->diamondPattern.getPatternChannels(1);

	//TF1* fit0 = landauGauss.doLandauGaussFit(hCellsHarris18Good);
	histSaver->SaveHistogram(hCellsHarris18Good);
	TF1* fit1 = landauGauss.doLandauGaussFit(hCellsHarris10Bad);
	histSaver->SaveHistogram(hCellsHarris10Bad);
	TF1* fit4 = landauGauss.doLandauGaussFit(hCellsOverlayedColumnLandau);
	histSaver->SaveHistogram(hCellsOverlayedColumnLandau);

	//TF1* fit5 = landauGauss.doLandauGaussFit(hTransparentCharge3D);
	histSaver->SaveHistogram(hTransparentCharge3D);

	//To normalise to histograms onto the same figure
	TLegend* Legend = new 	TLegend(0.5,0.6,0.79,0.79);  //const char* header = "", Option_t* option = "brNDC")
	cHarrisGoodandStripNormailsed = new TCanvas("cHarrisGoodandStripNormailsed","cHarrisGoodandStripNormailsed");
	cHarrisGoodandStripNormailsed->cd();
	Legend->AddEntry(hLandau.at(0),"Strip Detector at 500V");   //,"1");  //, Option_t* option = "lpf")
	hLandau.at(0)->GetYaxis()->SetTitle("Number of Entries Normalised");
	hLandau.at(0)->SetStats(kFALSE);
	hLandau.at(0)->SetLineColor(2);
	hLandau.at(0)->DrawNormalized();
	Legend->AddEntry(hCellsHarris18Good,"3D Detector at 25V");   //"1");
	hCellsHarris18Good->SetStats(kFALSE);
	hCellsHarris18Good->SetLineColor(4);
	hCellsHarris18Good->DrawNormalized("same");
	Legend->Draw();
	histSaver->SaveCanvas(cHarrisGoodandStripNormailsed);

	//To normalise to histograms onto the same figure
	TLegend* Legend1 = new 	TLegend(0.5,0.6,0.79,0.79);  //const char* header = "", Option_t* option = "brNDC")
	cGoodGradedandStripNormailsed = new TCanvas("cGoodGradedandStripNormailsed","cGoodGradedandStripNormailsed");
	cGoodGradedandStripNormailsed->cd();
	Legend1->AddEntry(hLandau.at(0),"Strip Detector at 500V");   //,"1");  //, Option_t* option = "lpf")
	hLandau.at(0)->GetYaxis()->SetTitle("Number of Entries Normalised");
	hLandau.at(0)->SetStats(kFALSE);
	hLandau.at(0)->SetLineColor(2);
	hLandau.at(0)->DrawNormalized();
	Legend1->AddEntry(hCellsLandauGraded.at(0),"3D Detector at 25V");   //"1");
	hCellsLandauGraded.at(0)->SetStats(kFALSE);
	hCellsLandauGraded.at(0)->SetLineColor(4);
	hCellsLandauGraded.at(0)->DrawNormalized("same");
	Legend1->Draw();
	histSaver->SaveCanvas(cGoodGradedandStripNormailsed);

	//To show mean charge around Strip detector

	//hDetXvsDetY3DMeanChargeQuarterCellGrading
	//vector<TF1*> GradedFits;
	for(int k=0;k<6;k++){
		TString hName = TString::Format("hDetXvsDetY3DMeanChargeQuarterCellGrading%d",k);
		cDetXvsDetY3DMeanChargeQuarterCellGrading.push_back(new TCanvas(hName,hName));
		cDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->cd();
		hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanChargeQuarterCellGrading");		//Set title to require
		hGridReferenceDetSpace->Draw("COL");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->Draw("sameCOLZAH");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->Draw("sameTEXTAH");
		//gStyle->SetPaintTextFormat(3.2g);
		settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeQuarterCellGrading.at(k), 3);
		histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeQuarterCellGrading.at(k));
		if(k!=5){
			//For Transparent Analysis
			TF1* fit2 = landauGauss.doLandauGaussFit(hQuarterCellGradedTransparentLandau.at(k));
			histSaver->SaveHistogram(hQuarterCellGradedTransparentLandau.at(k));

			TF1* fit1 = landauGauss.doLandauGaussFit(hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.at(k));
			histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.at(k));

			//hCellsDeltaXQuarterCellGrading
			//TF1* fit2 = landauGauss.doLandauGaussFit(hCellsDeltaXQuarterCellGrading.at(k));
			histSaver->SaveHistogram(hCellsDeltaXQuarterCellGrading.at(k));
		}
	}

	//Transparent Analysis
	//hCellsTransparentHitPositionCellGraded
	for(int i=0;i<4;i++){
		TString hName = TString::Format("cCellsTransparentHitPositionCellGraded%d",i);
		cCellsTransparentHitPositionCellGraded.push_back(new TCanvas(hName, hName));
		cCellsTransparentHitPositionCellGraded.at(i)->cd();
		hCellsTransparentHitPositionCellGraded.at(i)->SetStats(kFALSE);
		hCellsTransparentHitPositionCellGraded.at(i)->Draw("COLZ");
		hCellsTransparentHitPositionCellGraded.at(i)->Draw("sameTEXT");
		histSaver->SaveCanvas(cCellsTransparentHitPositionCellGraded.at(i));
	}
	//hDetXvsDetY3DMeanChargeHighlightedQuarters
	cDetXvsDetY3DMeanChargeHighlightedQuarters = new TCanvas("cDetXvsDetY3DMeanChargeHighlightedQuarters","cDetXvsDetY3DMeanChargeHighlightedQuarters");
	cDetXvsDetY3DMeanChargeHighlightedQuarters->cd();
	hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanChargeHighlightedQuarters");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hDetXvsDetY3DMeanCharge->Draw("sameTEXTAH");
	//hDetXvsDetY3DMeanCharge->Draw("sameCOLZ");
	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeHighlightedQuarters, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeHighlightedQuarters);

	//h3DdetQuarterCellFluctuation
	c3DdetQuarterCellFluctuation1 = new TCanvas("c3DdetQuarterCellFluctuation","c3DdetQuarterCellFluctuation");
	c3DdetQuarterCellFluctuation1->cd();
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellFluctuation");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	gStyle->SetPaintTextFormat("3.2g");
	h3DdetQuarterCellFluctuation->Draw("sameTEXT");
	settings->DrawMetallisationGrid(c3DdetQuarterCellFluctuation1, 3);
	histSaver->SaveCanvas(c3DdetQuarterCellFluctuation1);


	//h3DdetQuarterCellFluctuation1
	c3DdetQuarterCellFluctuation = new TCanvas("c3DdetQuarterCellFluctuation1","c3DdetQuarterCellFluctuation1");
	c3DdetQuarterCellFluctuation->cd();
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellFluctuation1");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	gStyle->SetPaintTextFormat("3.2g");
	h3DdetQuarterCellFluctuation1->Draw("sameTEXT");
	settings->DrawMetallisationGrid(c3DdetQuarterCellFluctuation, 3);
	histSaver->SaveCanvas(c3DdetQuarterCellFluctuation1);

	gStyle->SetPaintTextFormat("g");



	//hCellsOverlayedEventsNoColumns
	cCellsEventsNoColumn = new TCanvas("cCellsEventsNoColumn","cCellsEventsNoColumn");
	cCellsEventsNoColumn->cd();
	hCellsOverlayedEventsNoColumns->Draw("COLZ");
	hCellsOverlayedEventsNoColumns->Draw("sameTEXT");
	histSaver->SaveCanvas(cCellsEventsNoColumn);

	//hDetXvsDetY3DRebinnedRMS
	cDetXvsDetY3DRebinnedRMS = new TCanvas("cDetXvsDetY3DRebinnedRMS","cDetXvsDetY3DRebinnedRMS");
	cDetXvsDetY3DRebinnedRMS->cd();
	hGridReferenceDetSpace->SetTitle("h3DdetRebinnedRMS");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DRebinnedRMS->Draw("sameCOLZAH");
	hDetXvsDetY3DRebinnedRMS->Draw("sameTEXTAH");
	settings->DrawMetallisationGrid(cDetXvsDetY3DRebinnedRMS, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DRebinnedRMS);

	//hCellsLandauGraded			//Save each of the histograms
	for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
		//histSaver->SaveHistogram(hCellsLandauGraded.at(k));
		//histSaver->SaveHistogram(hCellsLandauGradedNoColumn.at(k));
	}

	//hCellsLandau2D
	TH2D* hDetXvsDetY3DMeanChargeRebinnedForCellOrdering = (TH2D*)hDetXvsDetY3DMeanChargeRebinned->Clone("hDetXvsDetY3DMeanChargeRebinnedForCellOrdering");
	int hCellsLandau2DEntries =0;
	int hCellsLandau2DEntriesPlus =0;
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){
			hCellsLandau2DEntriesPlus = hCellsLandau2DEntries;
			hCellsLandau2DEntries = hCellsLandau2DEntriesPlus + hCellsLandau.at(i*settings->getNRows3d()+j)->GetEntries();

			//hCellsMeanClusteSize
			hCellsMeanClusteSize->SetBinContent(i+1,j+1,hCellsClusteSize.at(i*settings->getNRows3d()+j)->GetMean());
			//hQuarterCellsMeanClusteSize
			for(int l=0;l<4;l++){	//To fill each quarter cell of the graded cell.
				if(l<2){
					hQuarterCellsMeanClusterSize->SetBinContent((2*i+1),(2*j+1+l),hQuarterCellsClusterSize.at(i*settings->getNRows3d()*4+j*4+l)->GetMean());
				}
				if(l>1){
					hQuarterCellsMeanClusterSize->SetBinContent((2*i+2),(2*j+1+l-2),hQuarterCellsClusterSize.at(i*settings->getNRows3d()*4+j*4+l)->GetMean());
				}
			}

			//Fit Landau to each cells landau.
			//hCellsLandau.at(i*settings->getNRows3d()+j)->Fit("landau");
			/*Int_t minX, minY, minZ;
			Int_t maxX, maxY, maxZ;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMaximumBin(maxX, maxY, maxZ);
			cout<<"The Maximum bin is: "<<maxX<<", "<<maxY<<endl;
			cout<<"The Maximum ADC is: "<<hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMaximum()<<endl;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimumBin(minX, minY, minZ);
			cout<<"The Minimum bin is: "<<minX<<", "<<minY<<endl;;
			cout<<"The Minimum ADC is: "<<hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimum()<<endl;
	 *//*
			Int_t minX, minY, minZ;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimumBin(minX, minY, minZ);	//Returns cell with Minimum average.
			//cout<<"hCellsLandau2D bin in y: "<<(i*settings->getNRows3d()+j+1)<<".   Filled with Cell: "<<((minX-1)*settings->getNRows3d()+(minY-1))<<endl;
			stringstream binLabel; binLabel<<((minX-1)*settings->getNRows3d()+(minY-1));
			hCellsLandau2D->GetYaxis()->SetBinLabel(i*settings->getNRows3d()+j+1,binLabel.str().c_str());
			hCellsLandau2D->GetYaxis()->SetLabelSize(0.02);
			hCellsLandau2DQuarterFail->GetYaxis()->SetBinLabel(i*settings->getNRows3d()+j+1,binLabel.str().c_str());
			hCellsLandau2DQuarterFail->GetYaxis()->SetLabelSize(0.02);
			hCellsLandau2DQuarterFail->SetBinContent(1,(i*settings->getNRows3d()+j+1),RebinnedQuarterCellFails->GetBinContent(minX,minY));
			for(int k=0;k<256;k++){
				hCellsLandau2D->SetBinContent((k+1),(i*settings->getNRows3d()+j+1),hCellsLandau.at((minX-1)*settings->getNRows3d()+(minY-1))->GetBinContent((k+1)));		//Set Cell in Y, set all pulseheights in x.
			}	//End of cell Landau loop
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->SetBinContent(minX, minY, 99999);	//The minimum bin now becomes maximum.
		}
	}	//End of looping over all cells


	//RebinnedQuarterCellFails
	cRebinnedQuarterCellFails = new TCanvas("c3DdetNumberofQuarterCellFails","c3DdetNumberofQuarterCellFails");
	cRebinnedQuarterCellFails->cd();
	RebinnedQuarterCellFails->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetNumberofQuarterCellFails");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	RebinnedQuarterCellFails->Draw("sameCOLZAH");
	RebinnedQuarterCellFails->Draw("sameTEXTAH");
	settings->DrawMetallisationGrid(cRebinnedQuarterCellFails, 3);
	histSaver->SaveCanvas(cRebinnedQuarterCellFails);

	//hDetXvsDetY3DOverview
	hOverview = new TCanvas("c3DdetOverview","c3DdetOverview");
	hOverview->cd();
	hDetXvsDetY3DOverview->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetOverview");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DOverview->Draw("sameCOLZAH");
	histSaver->SaveCanvas(hOverview);


	}
	  */


	/////////////////////////  Felix is working on ///////////////////////////
	/*
	//For h3DdetMeanCharge
	cDetXvsDetY3DMeanCharge = new TCanvas("c3DdetMeanCharge","c3DdetMeanCharge");
	cDetXvsDetY3DMeanCharge->cd();
	 *hDetXvsDetY3DMeanCharge = (*hDetXvsDetY3DCharge/(*hDetXvsDetY3DEvents));
	hDetXvsDetY3DMeanCharge->SetEntries(hDetXvsDetY3DEvents->Integral());
	//hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanCharge");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hGridReference->Draw("COL");
	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanCharge, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanCharge);

	//RebinnedMeanCharge
	cDetXvsDetY3DMeanChargeRebinned = new TCanvas("c3DdetMeanChargeRebinned","c3DdetMeanChargeRebinned");
	cDetXvsDetY3DMeanChargeRebinned->cd();
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	 *hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DChargeRebinned/(*hDetXvsDetY3DEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DMeanChargeRebinned");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	hDetXvsDetY3DEventsRebinned->Draw("sameTEXTAH");
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeRebinned, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeRebinned);
	/*c3DdetEventsRebinned = new TCanvas("c3DdetEventsRebinned","c3DdetEventsRebinned");
	//c3DdetEventsRebinned->cd();
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");			//("COLZ");
	//settings->DrawMetallisationGrid();
	histSaver->SaveCanvas(c3DdetEventsRebinned);
	 *//*

	//RebinnedMeanChargeQuarterCell
	for(int column=0;column<settings->getNColumns3d();column++){
		for(int row=0;row<settings->getNRows3d();row++){

			Int_t cell = column*settings->getNRows3d() + row;
			hCellNumbering->SetBinContent(column+1,row+1,cell);
			hCellMeanCharge->Fill(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1));	//Fill 1D hist with each cells mean charge


			for(int quarter=0;quarter<settings->getNQuarters3d();quarter++){	//For each quater cell

				Float_t xBin = 2*column + (quarter/2)+1;
				Float_t yBin = 2*row +quarter+1-2*(quarter/2);

				Float_t mean = hQuarterCellsLandau[cell][quarter]->GetMean();
				Float_t rms = hQuarterCellsLandau[cell][quarter]->GetRMS();

				hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetBinContent(xBin,yBin,mean);
				hDetXvsDetY3DRebinnedQuarterCellRMS->SetBinContent(xBin, yBin, rms);

			}
		}
	}


	/*int QuaterCellEntrySum = 0;
	int QuaterCellEntrySumAddition = 0;
	for(int column=0;column<settings->getNColumns3d();column++){
		for(int row=0;row<settings->getNRows3d();row++){
			for(int quarter=0;quarter<settings->getNQuarters3d();quarter++){	//For each quater cell

				Float_t xBin = 2*column + (quarter/2)+1;
				Float_t yBin = 2*row +quarter+1-2*(quarter/2);
				Int_t entry = column*settings->getNRows3d()*4+row*4+quarter;
				if (entry<hQuaterCellsLandau.size()){
					Float_t mean = hQuaterCellsLandau.at(entry)->GetMean();
					Float_t rms = hQuaterCellsLandau.at(entry)->GetRMS();
					if(hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetNbinsX()>xBin && hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetNbinsY()>yBin)
						hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetBinContent(xBin,yBin,mean);
					else
						cout<<" hDetXvsDetY3DMeanChargeRebinnedQuarterCell, bin does not exist: " << xBin<<"/"<<yBin<<endl;
					if(hDetXvsDetY3DRebinnedQuarterCellRMS->GetNbinsX()>xBin && hDetXvsDetY3DRebinnedQuarterCellRMS->GetNbinsY()>yBin)
						hDetXvsDetY3DRebinnedQuarterCellRMS->SetBinContent(xBin, yBin, rms);
					else
						cout<<" hDetXvsDetY3DRebinnedQuarterCellRMS, bin does not exist: " << xBin<<"/"<<yBin<<endl;
					QuaterCellEntrySumAddition = QuaterCellEntrySum;
					QuaterCellEntrySum = QuaterCellEntrySumAddition + hQuaterCellsLandau.at(entry)->GetEntries();
					//histSaver->SaveHistogram(hQuaterCellsLandau.at(entry));
					if(entry<hQuarterCellsClusterSize.size()){
						cout<<"save: hQuarterCellsClusterSize["<<entry<<"] / "<<hQuarterCellsClusterSize.size()<<endl;
						cout<<"hQuarterCellsClusterSize: "<<flush<<" "<<hQuarterCellsClusterSize[entry]<<endl;
						//histSaver->SaveHistogram(hQuarterCellsClusterSize.at(entry));
					}
				}
				else{
					cout<<"size of hQuaterCellsLandau is to small"<<entry<<"/"<<hQuaterCellsLandau.size()<<endl;
				}
			}
		}
	}

		cout<<"done"<<endl;
	cout<<"cDetXvsDetY3DMeanChargeRebinnedQuarterCell"<<endl;
	cDetXvsDetY3DMeanChargeRebinnedQuarterCell = new TCanvas("c3DdetQuarterCellMeanCharge","c3DdetQuarterCellMeanCharge");
	cDetXvsDetY3DMeanChargeRebinnedQuarterCell->cd();
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetEntries(QuaterCellEntrySum);
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellMeanCharge");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameCOLZAH");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameTEXTAH");
	settings->DrawMetallisationGrid(cDetXvsDetY3DMeanChargeRebinnedQuarterCell, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeRebinnedQuarterCell);

	cout<<"cDetXvsDetY3DRebinnedRMSQuarterCell"<<endl;
	cDetXvsDetY3DRebinnedRMSQuarterCell = new TCanvas("c3DdetQuarterCellRMS","c3DdetQuarterCellRMS");
	cDetXvsDetY3DRebinnedQuarterCellRMSCanvas->cd();
	hDetXvsDetY3DRebinnedRMSQuarterCell->SetEntries(QuaterCellEntrySum);
	hGridReferenceDetSpace->SetTitle("hDetXvsDetY3DRebinnedRMSQuarterCell");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DRebinnedRMSQuarterCell->Draw("sameCOLZAH");
	hDetXvsDetY3DRebinnedRMSQuarterCell->Draw("sameTEXTAH");
	settings->DrawMetallisationGrid(cDetXvsDetY3DRebinnedRMSQuarterCell, 3);
	histSaver->SaveCanvas(cDetXvsDetY3DRebinnedRMSQuarterCell);



			}
		}
	}*/
}

float* TAnalysisOf3dDiamonds::SortArrayBtoS(float* nArray, int nSize) {

	bool bDone = false; // this flag will be used to check whether we have to continue the algorithm
	int length = nSize; //This actually isn't needed but I need something to hold the value of size so I can print it properly.

	//printArray(nArray, nSize, "-"); // print the initial array
	while (!bDone){
		bDone = true; // assume that the array is currently sorted
		for (int i = 0; i != length - 1; ++i){ // for every element in the array     ( notice: this can be a bit optimized, see http://en.wikipedia.org/wiki/Bubblesort#Optimizing_bubble_sort )
			if ( nArray[i] > nArray[i + 1] ){ // compare the current element with the following one
				// They are in the wrong order, swap them
				float Temp = nArray[i];
				float Temp2 = SortArrayPointer.at(i);
				nArray[i] = nArray[i+1];
				SortArrayPointer.at(i) = SortArrayPointer.at(i+1);
				nArray[i+1] = Temp;
				SortArrayPointer.at(i+1) = Temp2;
				bDone = false; // since we performed a swap, the array needs to be checked to see if it is sorted
				// this is done in the next iteration of the while
				//printArray(nArray, nSize, "-"); // print current array to show the steps done by the algorithm
			}
		}
		length -= 1; //After each full iteration of the while loop, we know that the last element in the array is the highest number
		//Because we know this, we can lower the amount of iterations in the array we need to do and simply subtract 1.
	}
	return nArray;
}

void TAnalysisOf3dDiamonds::printArray(float* nArray, int nSize, const std::string& space){
	for (int i = 0; i != nSize; ++i){
		if (i != nSize - 1)
			std::cout << nArray[i] << space;
		else
			std::cout << nArray[i];
	}
	std::cout << std::endl;
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
		hFidCutXvsFidCutYClusters.at(0)->Fill(nfiducialValueX,nfiducialValueY,1);
	}
	if(nClusters==1){
		hFidCutXvsFidCutYClusters.at(1)->Fill(nfiducialValueX,nfiducialValueY,1);
	}
	if(nClusters==1){
		if(HitCount==0&&SeedCount==1){
			hFidCutXvsFidCutYClusters.at(2)->Fill(nfiducialValueX,nfiducialValueY,1);
		}
	}
	if(nClusters==2){
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


void TAnalysisOf3dDiamonds::RemoveLumpyClusters(TCluster* nCluster) {
	//Removes clusters with inbetween hits
	Int_t Left=0; Int_t Span=0; Int_t Right=0;
	for (UInt_t clPos=0;clPos <= nCluster->getClusterSize();clPos++){
		if(nCluster->isSeed(clPos)) Left=1;
		if(Left==1&&nCluster->isHit(clPos)&&!nCluster->isSeed(clPos)) Span=1;
		if(Left==1&&Span==1&&nCluster->isSeed(clPos)){
			Right=1;
			cout<<"Removed: "<<endl;
			nCluster->Print(1);
		}
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

int TAnalysisOf3dDiamonds::RemoveClustersWithHitOutside3D(TCluster* nCluster) {
	int Remove =0;
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)){
			if(nCluster->getChannel(clPos)>93||nCluster->getChannel(clPos)<85)
				Remove =1;
		}
	}
	return Remove;
}

int TAnalysisOf3dDiamonds::RemoveEdgeClusters(TCluster* nCluster, int nDetector) {
	pair<int,int> channels = settings->diamondPattern.getPatternChannels(nDetector);
	//cout<<channels.first<<"  "<<channels.second<<endl;
	int Remove =0;
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)){
			if(nCluster->getChannel(clPos)>=channels.second||nCluster->getChannel(clPos)<=channels.first)
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

Float_t TAnalysisOf3dDiamonds::getTransparentCharge(Int_t nDiamondPattern, Int_t nChannelHit) {

	Float_t TransparentCharge = -9999;
	pair<int,int> channels = settings->diamondPattern.getPatternChannels(nDiamondPattern);
	Int_t ChannelsToRead;
	if(nDiamondPattern==1)
		ChannelsToRead = 5;
	else{ChannelsToRead = 3;}

	vector<Float_t> ChannelsCharge;
	float* ptrChannelsChargeSorted;

	//cout<<XdetChannelSpaceInt<<endl;
	Int_t StripStart = -(ChannelsToRead/2);
	//cout<<"StripStart: "<<StripStart<<endl;
	Int_t ChannelStart = nChannelHit + StripStart;
	//cout<<"ChannelStart: "<<ChannelStart<<endl;
	Int_t DetectorLeftMostChannel = channels.first;
	Int_t DetectorRightMostChannel = channels.second;

	SortArrayPointer.clear();
	for(int i=0;i<ChannelsToRead;i++){
		Int_t Channel = ChannelStart + i;
		SortArrayPointer.push_back(Channel);
		ChannelsCharge.push_back(eventReader->getSignal(subjectDetector,Channel));
		//cout<<"Channel: "<<Channel<<"   Charge: "<<ChannelsCharge.at(i)<<endl;
	}

	ptrChannelsChargeSorted = SortArrayBtoS(VectorToArray(ChannelsCharge), ChannelsCharge.size());

	//cout<<"DiamondPattern"<<nDiamondPattern<<endl;
	for(int i=0; i<ChannelsCharge.size();i++){
		//cout<<SortArrayPointer.at(i)<<"   "<<ptrChannelsChargeSorted[i]<<endl;
	}

	Float_t HighestCharge = ptrChannelsChargeSorted[ChannelsCharge.size()-1];
	Int_t HighestChargeChannel = SortArrayPointer.at(ChannelsCharge.size()-1);
	//cout<<HighestChargeChannel<<"   "<<HighestCharge<<endl;
	Float_t SecondHighestCharge = ptrChannelsChargeSorted[ChannelsCharge.size()-2];
	Int_t SecondHighestChargeChannel = SortArrayPointer.at(ChannelsCharge.size()-2);
	//cout<<SecondHighestChargeChannel<<"   "<<SecondHighestCharge<<endl;

	Int_t Saturated =0;
	if(eventReader->isSaturated(subjectDetector,HighestChargeChannel)||eventReader->isSaturated(subjectDetector,SecondHighestChargeChannel))
		Saturated = 1;
	if(Saturated == 1||HighestChargeChannel == DetectorLeftMostChannel||HighestChargeChannel == DetectorRightMostChannel)    //||SecondHighestChargeChannel == DetectorLeftMostChannel||SecondHighestChargeChannel == DetectorRightMostChannel)
		TransparentCharge = -9999;
	else{
		TransparentCharge = HighestCharge + SecondHighestCharge;
	}

	return TransparentCharge;

	/*
		Int_t Channel = ChannelStart + i;
		//cout<<"Channel: "<<Channel<<endl;
		TransparentChargeAddition = TransparentCharge;

		if(Channel<DetectorLeftMostChannel||Channel>DetectorRightMostChannel)		//Channel is outside detector
			TransparentCharge = 0 + TransparentChargeAddition;
		TransparentCharge = eventReader->getSignal(subjectDetector,Channel) + TransparentChargeAddition;

	}




	Float_t TransparentCharge = 0;
	Float_t TransparentChargeAddition = 0;

	for(int i=0;i<ChannelsToRead;i++){
		Int_t Channel = ChannelStart + i;
		//cout<<"Channel: "<<Channel<<endl;
		TransparentChargeAddition = TransparentCharge;

		if(Channel<DetectorLeftMostChannel||Channel>DetectorRightMostChannel)		//Channel is outside detector
			TransparentCharge = 0 + TransparentChargeAddition;
		TransparentCharge = eventReader->getSignal(subjectDetector,Channel) + TransparentChargeAddition;
	}
	 */
}

float* TAnalysisOf3dDiamonds::VectorToArray(vector<float> nvector){
	float* ret=new float[nvector.size()];
	for (int i=0;i<nvector.size();i++){
		ret[i]=nvector.at(i);
	}
	return ret;
}

/*void TAnalysisOf3dDiamonds::fillEdgeDistributions(Float_t clusterCharge){
	for(UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts();i++ ){
		TFiducialCut* fidCut = settings->get3dEdgeFidCuts()->getFidCut(i+1);
		if(!fidCut)
			continue;
		if(fidCut->isInFiducialCut(fiducialValueX,fiducialValueY)){
			vecEdgePredX[i].push_back(xPredDet);
			vecEdgePredY[i].push_back(yPredDet);
			vecEdgePulseHeight[i].push_back(clusterCharge);
			//				cout<<nEvent<<": "<<xPredDet<<"/"<<yPredDet<<" --> "<<i<<endl;
		}
	}

}*/



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
	Float_t xmax = hCh->GetZaxis()->GetXmax();
	Float_t xmin = hCh->GetZaxis()->GetXmin();
//	hCh->GetZaxis()->SetRangeUser(xmin+.1*(xmax-xmin),xmax);
	histSaver->SaveHistogram(hCh,false);
	delete hCh;
}

void TAnalysisOf3dDiamonds::initialiseHistos() {
	// Initialise
	hValidEventsFiducialSpace = new TH2F("hValidEventsFiducialSpace","hValidEventsFiducialSpace",1024,0,256,1024,0,256);
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
	hValidEventsDetSpace = new TH2F("hValidEventsDetSpace","hValidEventsDetSpace",1024,xmin,xmax,1024,ymin,ymax);
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
	if(settings->do3dLongAnalysis() == 1){SaveLongAnalysisHistos();SaveLongAnalysisHistos2();}
	if(settings->do3dTransparentAnalysis() == 1){saveTransparentAnalysisHistos();/*saveTransparentAnalysisHistos()*/}
}

void TAnalysisOf3dDiamonds::initialiseLongAnalysisHistos() {
	initialise3DGridReference();
	initialise3DYAlignmentHistos();
	initialise3DOverviewHistos();
	initialise3D2DLandauAndClustersizeHistos();
	initialise3DCellOverlayHistos();
	hLongAnalysisInvalidCellNo = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCellNo");
	hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCellNo");
	hLongAnalysisInvalidCluster = (TH2F*) hValidEventsDetSpace->Clone("hLongAnalysisInvalidCluster");
	hLongAnalysisInvalidCellNo->SetTitle("hLongAnalysisInvalidCluster");
}

void TAnalysisOf3dDiamonds::ShortAnalysis_SaveMeanChargeVector() {
	TH3F* hChargeDistribution = histSaver->Create3DHisto("hChargeDistribution3D",vecPredDetX_ShortAna,vecPredDetY_ShortAna,vecPulseHeight_ShortAna,1024,1024);
	hChargeDistribution->GetXaxis()->SetTitle("Predicted Hit Position X in Detector /#mum");
	hChargeDistribution->GetYaxis()->SetTitle("Predicted Hit Position Y in Detector /#mum");
	hChargeDistribution->GetZaxis()->SetTitle("Pulse Height of Cluster  / ADC");
	if (!hChargeDistribution)
		return;
	else
		cerr<<" hChargDistribution3D: was not created:"<<endl;
	cout<<"Create Project3dProfile for "<<hChargeDistribution->GetName()<<"..."<<flush;
	TH2F* hMeanCharge = (TH2F*)hChargeDistribution->Project3DProfile("yx");
	cout<<"\t[done]"<<endl;
	//	if(!hMeanCharge)
	//		return;
	//	else
	//		cerr<<" hChargDistribution3D_pfyx: was not created:"<<endl;
	hMeanCharge->GetXaxis()->SetTitle("Predicted Hit Position X in Detector /#mum");
	hMeanCharge->GetYaxis()->SetTitle("Predicted Hit Position Y in Detector /#mum");
	hMeanCharge->GetZaxis()->SetTitle("Avrg. pulse height /ADC");
	hMeanCharge->SetTitle("Avrg. pulse height in detector system");
	hMeanCharge->SetName("hAvrgPulseHeigthDetSystem");
	histSaver->SaveHistogram(hMeanCharge,false);
	TCanvas *c1 = new TCanvas("cAvrgPulseHeigthDetSystem_MetalizationLayer");
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
	c1->SetName("cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom");
	histSaver->SaveCanvas(c1);
	delete c1;
	if (hChargeDistribution)
		delete hChargeDistribution;
}

void TAnalysisOf3dDiamonds::InitialiseStripAnalysisHistos() {
	TString name = "hLandauStrip";
	hLandauStrip =  new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
	hLandauStrip->GetXaxis()->SetTitle("charge /ADC");
	hLandauStrip->GetYaxis()->SetTitle("number of entries #");
	hLandauStrip->SetLineColor(kBlue);
}

void TAnalysisOf3dDiamonds::SaveStripAnalysisHistos() {
	histSaver->SaveHistogram(hLandauStrip);
}

void TAnalysisOf3dDiamonds::LongAnalysis_SaveMeanChargePlots() {

	histSaver->SaveHistogram(hPulseHeightVsDetectorHitPostionXY);
	histSaver->SaveHistogramWithCellGrid(hPulseHeightVsDetectorHitPostionXY);

	UInt_t xBins = hPulseHeightVsDetectorHitPostionXY->GetXaxis()->GetNbins();
	UInt_t yBins = hPulseHeightVsDetectorHitPostionXY->GetYaxis()->GetNbins();
	UInt_t xRebin = xBins/settings->getNColumns3d()/2;
	UInt_t yRebin = yBins/settings->getNRows3d()/2;

	TString name = "hPulseHeightVsDetectorHitPostion_Quarters";
	TH2F* hDetXvsDetY3DMeanChargeQuarters = (TH2F*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin,yRebin,name);
	histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeQuarters);
	histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeQuarters);
	delete hDetXvsDetY3DMeanChargeQuarters;
	//
	name = "hPulseHeightVsDetectorHitPostion_Cells";
	TH2D* hDetXvsDetY3DMeanChargeCells = (TH2D*)hPulseHeightVsDetectorHitPostionXY->Rebin2D(xRebin*2,yRebin*2,name);
	histSaver->SaveHistogramWithCellGrid(hDetXvsDetY3DMeanChargeCells);
	histSaver->SaveHistogram(hDetXvsDetY3DMeanChargeCells);
	delete hDetXvsDetY3DMeanChargeCells;
	//
	delete hPulseHeightVsDetectorHitPostionXY;

}
