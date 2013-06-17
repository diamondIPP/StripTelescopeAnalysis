/*
 * TAnalysisOf3dDiamonds.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair,iain
 */

#include "../include/TAnalysisOf3dDiamonds.hh"

TAnalysisOf3dDiamonds::TAnalysisOf3dDiamonds(TSettings *newSettings) {
	if(newSettings!=0)
		this->settings=newSettings;
	else exit(0);

	UInt_t runNumber=settings->getRunNumber();

	//htmlLandau=new THTMLLandaus(settings);

	settings->goTo3dDiamondTreeDir();
	eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
	histSaver=new HistogrammSaver();
	settings->goTo3dDiamondAnalysisDir();

	histSaver->SetPlotsPath(settings->get3dDiamondAnalysisPath());
	histSaver->SetRunNumber(runNumber);
	//htmlLandau->setFileGeneratingPath(settings->getSelectionAnalysisPath());//todo Write html3dDiamond
	settings->goTo3dDiamondTreeDir();
	//initialiseHistos();

	cout<<"end initialise"<<endl;
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
	//Print out detectors to be analysed
	int Strip[] = {26,38};
	int threeDnH[] = {56,62};
	int threeDwH[] = {86,92};

	//int Strip[] = {32,32};
	//int threeDnH[] = {59,59};
	//int threeDwH[] = {88,88};

	Detector.push_back(Strip); //Strip
	Detector.push_back(threeDnH); //3D no holes
	Detector.push_back(threeDwH); //3D with holes
	//Fiducial Cut Regions       (xHigh, xLow, yHigh, yLow,     xChHigh, xChLow, yChHigh, yChLow)
	int StripFid[] =             {102,   92,   110,   68,       104,     93,     96,      68};
	int threeDnHFid[] =          {102,   92,   110,   68,       130,     112,    105,     98};
	int threeDwHFid[] =          {102,   92,   110,   68,       161,     141,    107,     70};
	FidCut.push_back(StripFid); //Strip
	FidCut.push_back(threeDnHFid); //3D no holes
	FidCut.push_back(threeDwHFid); //3D with holes

	FiducialMetric=0;			// Which Fiducial cut to apply
	FiducialChannel=1;

	/*for(int i=0;i<FidCut.size();i++){
		for(int j=0; j<8;j++){
			cout<<FidCut.at(i)[j]<<endl;
		}
	}
		*/

	//For YAlignment Fiducial cuts
	//Fiducial Cut Regions    (xHigh, xLow, yHigh, yLow,     xChHigh, xChLow, yChHigh, yChLow)
	int xEdge[] =             {102,   92,   110,   68,       170,     162,    100,     83};
	int yEdge[] =             {102,   92,   110,   68,       164,     140,    112,     103};
	FidCutYAlignment.push_back(xEdge);
	FidCutYAlignment.push_back(yEdge);
	for(int i=0;i<FidCutYAlignment.size();i++){
		for(int j=0; j<8;j++){
			cout<<FidCutYAlignment.at(i)[j]<<endl;
		}
	}

	//intialise vectors
	for(int i=0;i<Detector.size();i++){
		vecPHDiamondHit.push_back(new vector<float>);
		vecXPredicted.push_back(new vector<float>);
		vecYPredicted.push_back(new vector<float>);
		ptrCanvas.push_back(new TCanvas); //To Create pointer to canvas for 3D Plot
		ptrCanvasEvents.push_back(new TCanvas);
		ptrCanvasMean.push_back(new TCanvas);
	}
	cout<<"Areas to be analysed:"<<endl;
	for(int i=0; i<Detector.size(); i++){
		cout<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<endl;
	}
	//Moved initialiseHistos to here
	initialiseYAlignmentHistos();
	/*initialiseHistos();
	 	 */
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		eventReader->LoadEvent(nEvent);
		/*analyseEvent();
		 *
		 */
		YAlignment();
	}
	saveYAlignmentHistos();
	/*saveHistos();
	 	 */
}

void TAnalysisOf3dDiamonds::initialiseHistos() {
	//Universal histograms

	//hNumberofClusters
	stringstream hNumberofClustersName; hNumberofClustersName<<"hNumberofClusters%%"<<FileNameEnd;
	hNumberofClusters = new TH1F(hNumberofClustersName.str().c_str(),hNumberofClustersName.str().c_str(),4,0,4);
	hNumberofClusters->SetTitle(hNumberofClustersName.str().c_str());
	hNumberofClusters->GetXaxis()->SetTitle("Number of Clusters");
	hNumberofClusters->GetYaxis()->SetTitle("Number of Entries #");

	//hDoubleClusterPos
	stringstream hDoubleClusterPosName; hDoubleClusterPosName<<"hDoubleClusterPos%%"<<FileNameEnd;
	hDoubleClusterPos = new TH1F(hDoubleClusterPosName.str().c_str(),hDoubleClusterPosName.str().c_str(),80,20,100);
	hDoubleClusterPos->SetTitle(hDoubleClusterPosName.str().c_str());
	hDoubleClusterPos->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos->GetYaxis()->SetTitle("Number of Entries #");

	//hDoubleClusterPos0
	stringstream hDoubleClusterPos0Name; hDoubleClusterPos0Name<<"hDoubleClusterPos0%%"<<FileNameEnd;
	hDoubleClusterPos0 = new TH1F(hDoubleClusterPos0Name.str().c_str(),hDoubleClusterPos0Name.str().c_str(),80,20,100);
	hDoubleClusterPos0->SetTitle(hDoubleClusterPos0Name.str().c_str());
	hDoubleClusterPos0->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos0->GetYaxis()->SetTitle("Number of Entries #");
	hDoubleClusterPos0->SetFillColor(2);

	//hDoubleClusterPos1
	stringstream hDoubleClusterPos1Name; hDoubleClusterPos1Name<<"hDoubleClusterPos1%%"<<FileNameEnd;
	hDoubleClusterPos1 = new TH1F(hDoubleClusterPos1Name.str().c_str(),hDoubleClusterPos1Name.str().c_str(),80,20,100);
	hDoubleClusterPos1->SetTitle(hDoubleClusterPos1Name.str().c_str());
	hDoubleClusterPos1->GetXaxis()->SetTitle("HighestPH Channel Hit");
	hDoubleClusterPos1->GetYaxis()->SetTitle("Number of Entries #");
	hDoubleClusterPos1->SetFillColor(3);

	//hLandauCluster1
	stringstream hLandauCluster1Name; hLandauCluster1Name<<"hLandauCluster1%%"<<FileNameEnd;
	hLandauCluster1 = new TH1F(hLandauCluster1Name.str().c_str(),hLandauCluster1Name.str().c_str(),256,0,2800);
	hLandauCluster1->SetTitle(hLandauCluster1Name.str().c_str());
	hLandauCluster1->GetXaxis()->SetTitle("Number of Clusters");
	hLandauCluster1->GetYaxis()->SetTitle("Number of Entries #");
	hLandauCluster1->SetFillColor(2);

	//hLandauCluster2
	stringstream hLandauCluster2Name; hLandauCluster2Name<<"hLandauCluster2%%"<<FileNameEnd;
	hLandauCluster2 = new TH1F(hLandauCluster2Name.str().c_str(),hLandauCluster2Name.str().c_str(),256,0,2800);
	hLandauCluster2->SetTitle(hLandauCluster2Name.str().c_str());
	hLandauCluster2->GetXaxis()->SetTitle("Number of Clusters");
	hLandauCluster2->GetYaxis()->SetTitle("Number of Entries #");
	hLandauCluster2->SetFillColor(3);

	//hLandauDoubleCombined
	stringstream hLandauDoubleCombinedName; hLandauDoubleCombinedName<<"hLandauDoubleCombined%%"<<FileNameEnd;
	hLandauDoubleCombined = new TH1F(hLandauDoubleCombinedName.str().c_str(),hLandauDoubleCombinedName.str().c_str(),256,0,2800);
	hLandauDoubleCombined->SetTitle(hLandauDoubleCombinedName.str().c_str());
	hLandauDoubleCombined->GetXaxis()->SetTitle("Number of Clusters");
	hLandauDoubleCombined->GetYaxis()->SetTitle("Number of Entries #");

	for(int i=0; i<Detector.size(); i++){

		//hLandau
		stringstream hLandauName; hLandauName<<"hLandau%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hLandau.push_back(new TH1F(hLandauName.str().c_str(),hLandauName.str().c_str(),256,0,2800));
		hLandau.at(i)->SetTitle(hLandauName.str().c_str());
		hLandau.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
		hLandau.at(i)->GetYaxis()->SetTitle("number of entries #");

		//hPHvsChannel
		stringstream hEventsvsChannelName; hEventsvsChannelName<<"hEventsvsChannel%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hEventsvsChannel.push_back(new TH1F(hEventsvsChannelName.str().c_str(),hEventsvsChannelName.str().c_str(),90,10,100));
		hEventsvsChannel.at(i)->GetXaxis()->SetTitle("HighestPH [ch]");
		hEventsvsChannel.at(i)->GetYaxis()->SetTitle("No. Events");

		//hPHvsChannel
		stringstream hPHvsChannelName; hPHvsChannelName<<"hPHvsChannel%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hPHvsChannel.push_back(new TH2F(hPHvsChannelName.str().c_str(),hPHvsChannelName.str().c_str(),150,0,2900,80,20,100));
		hPHvsChannel.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsChannel.at(i)->GetYaxis()->SetTitle("XPos(channel)");
		hPHvsChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsChannel.at(i)->GetXaxis()->SetRangeUser(1,2600);
		//hPHvsChannel.at(i)->SetMaximum(3000);
		//hPHvsChannel.at(i)->SetMinimum(0);

		//hPHvsPredictedChannel
		stringstream hPHvsPredictedChannelName; hPHvsPredictedChannelName<<"hPHvsPredictedChannel%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hPHvsPredictedChannel.push_back(new TH2F(hPHvsPredictedChannelName.str().c_str(),hPHvsPredictedChannelName.str().c_str(),150,0,2900,320,20,100));
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsPredictedChannel.at(i)->GetYaxis()->SetTitle("XPos(Predicted channel)");
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetRangeUser(1,2600);
		//hPHvsPredictedChannel.at(i)->SetMaximum(3000);
		//hPHvsPredictedChannel.at(i)->SetMinimum(0);

		//hPHvsPredictedXPos
		stringstream hPHvsPredictedXPosName; hPHvsPredictedXPosName<<"hPHvsPredictedXPos%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hPHvsPredictedXPos.push_back(new TH2F(hPHvsPredictedXPosName.str().c_str(),hPHvsPredictedXPosName.str().c_str(),150,0,2900,50,5800,10200));
		hPHvsPredictedXPos.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsPredictedXPos.at(i)->GetYaxis()->SetTitle("XPos(um)");
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsPredictedXPos.at(i)->GetXaxis()->SetRange(1,2600);
		//hPHvsPredictedXPos.at(i)->SetMaximum(3000);
		//hPHvsPredictedXPos.at(i)->SetMinimum(0);

		//hHitandSeedCount
		stringstream hHitandSeedCountName; hHitandSeedCountName<<"hHitandSeedCount%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hHitandSeedCount.push_back(new TH2F(hHitandSeedCountName.str().c_str(),hHitandSeedCountName.str().c_str(),10,0,10,10,0,10));
		hHitandSeedCount.at(i)->GetXaxis()->SetTitle("Hit Count");
		hHitandSeedCount.at(i)->GetYaxis()->SetTitle("Seed Count");

		//hChi2XChi2Y
		stringstream hChi2XChi2YName; hChi2XChi2YName<<"hChi2XChi2Y%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hChi2XChi2Y.push_back(new TH2F(hChi2XChi2YName.str().c_str(),hChi2XChi2YName.str().c_str(),60,0,60,60,0,60));
		hChi2XChi2Y.at(i)->GetXaxis()->SetTitle("Chi2X");
		hChi2XChi2Y.at(i)->GetYaxis()->SetTitle("Chi2Y");

		//hFidCutXvsFidCutY
		stringstream hFidCutXvsFidCutYName; hFidCutXvsFidCutYName<<"hFidCutXvsFidCutY%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hFidCutXvsFidCutY.push_back(new TH2F(hFidCutXvsFidCutYName.str().c_str(),hFidCutXvsFidCutYName.str().c_str(),160,90,170,120,60,120));
		hFidCutXvsFidCutY.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutY.at(i)->GetYaxis()->SetTitle("FidCutY");

		//hFidCutXvsFidCutYvsCharge		For TH2D
		stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hFidCutXvsFidCutYvsCharge.push_back(new TH2D(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),213,90,170,160,60,120));
		hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

		//hFidCutXvsFidCutYvsEvents
		hFidCutXvsFidCutYvsEvents.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsEvents"));

		//hFidCutXvsFidCutYvsMeanCharge
		hFidCutXvsFidCutYvsMeanCharge.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone("hFidCutXvsFidCutYvsMeanCharge"));

		/*//hFidCutXvsFidCutYvsCharge		For TH3F
		stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%";
		hFidCutXvsFidCutYvsCharge.push_back(new TH3F(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),160,90,170,120,60,120,30,0,3000));
		hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");
		*/

	}

	//To calculate Efficiency
	//hFidCutXvsFidCutYvsPredictedEvents	For TH2D
	stringstream hFidCutXvsFidCutYvsPredictedEventsName; hFidCutXvsFidCutYvsPredictedEventsName<<"hFidCutXvsFidCutYvsPredictedEvents"<<FileNameEnd;
	hFidCutXvsFidCutYvsPredictedEvents = new TH2D(hFidCutXvsFidCutYvsPredictedEventsName.str().c_str(),hFidCutXvsFidCutYvsPredictedEventsName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutYvsPredictedEvents->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsPredictedEvents->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsPredictedEvents->GetZaxis()->SetTitle("Predicted Events");

	//hFidCutXvsFidCutYvsSeenEventsName	For TH2D
	stringstream hFidCutXvsFidCutYvsSeenEventsName; hFidCutXvsFidCutYvsSeenEventsName<<"hFidCutXvsFidCutYvsSeenEvents"<<FileNameEnd;
	hFidCutXvsFidCutYvsSeenEvents = new TH2D(hFidCutXvsFidCutYvsSeenEventsName.str().c_str(),hFidCutXvsFidCutYvsSeenEventsName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutYvsSeenEvents->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsSeenEvents->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsSeenEvents->GetZaxis()->SetTitle("Seen Events");

	//hEfficiency
	hEfficiency = (TH2D*)hFidCutXvsFidCutYvsSeenEvents->Clone("hEfficiency");

	//hFidCutXvsFidCutY0Clusters	For TH2D
	stringstream hFidCutXvsFidCutY0ClustersName; hFidCutXvsFidCutY0ClustersName<<"hFidCutXvsFidCutY0Clusters"<<FileNameEnd;
	hFidCutXvsFidCutY0Clusters = new TH2D(hFidCutXvsFidCutY0ClustersName.str().c_str(),hFidCutXvsFidCutY0ClustersName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY0Clusters->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY0Clusters->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY0Clusters->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY1_1Clusters	For TH2D
	stringstream hFidCutXvsFidCutY1_1ClustersName; hFidCutXvsFidCutY1_1ClustersName<<"hFidCutXvsFidCutY1_1Clusters"<<FileNameEnd;
	hFidCutXvsFidCutY1_1Clusters = new TH2D(hFidCutXvsFidCutY1_1ClustersName.str().c_str(),hFidCutXvsFidCutY1_1ClustersName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY1_1Clusters->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY1_1Clusters->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY1_1Clusters->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY1Clusters	For TH2D
	stringstream hFidCutXvsFidCutY1ClustersName; hFidCutXvsFidCutY1ClustersName<<"hFidCutXvsFidCutY1Clusters"<<FileNameEnd;
	hFidCutXvsFidCutY1Clusters = new TH2D(hFidCutXvsFidCutY1ClustersName.str().c_str(),hFidCutXvsFidCutY1ClustersName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY1Clusters->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY1Clusters->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY1Clusters->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY2Clusters	For TH2D
	stringstream hFidCutXvsFidCutY2ClustersName; hFidCutXvsFidCutY2ClustersName<<"hFidCutXvsFidCutY2Clusters"<<FileNameEnd;
	hFidCutXvsFidCutY2Clusters = new TH2D(hFidCutXvsFidCutY2ClustersName.str().c_str(),hFidCutXvsFidCutY2ClustersName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY2Clusters->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY2Clusters->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY2Clusters->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY2Clusters0	For TH2D
	stringstream hFidCutXvsFidCutY2Clusters0Name; hFidCutXvsFidCutY2Clusters0Name<<"hFidCutXvsFidCutY2Clusters0"<<FileNameEnd;
	hFidCutXvsFidCutY2Clusters0 = new TH2D(hFidCutXvsFidCutY2Clusters0Name.str().c_str(),hFidCutXvsFidCutY2Clusters0Name.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY2Clusters0->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY2Clusters0->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY2Clusters0->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY2Clusters1	For TH2D
	stringstream hFidCutXvsFidCutY2Clusters1Name; hFidCutXvsFidCutY2Clusters1Name<<"hFidCutXvsFidCutY2Clusters1"<<FileNameEnd;
	hFidCutXvsFidCutY2Clusters1 = new TH2D(hFidCutXvsFidCutY2Clusters1Name.str().c_str(),hFidCutXvsFidCutY2Clusters1Name.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY2Clusters1->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY2Clusters1->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY2Clusters1->GetZaxis()->SetTitle("Events");

	//hFidCutXvsFidCutY3Clusters	For TH2D
	stringstream hFidCutXvsFidCutY3ClustersName; hFidCutXvsFidCutY3ClustersName<<"hFidCutXvsFidCutY3Clusters"<<FileNameEnd;
	hFidCutXvsFidCutY3Clusters = new TH2D(hFidCutXvsFidCutY3ClustersName.str().c_str(),hFidCutXvsFidCutY3ClustersName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutY3Clusters->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutY3Clusters->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutY3Clusters->GetZaxis()->SetTitle("Seen Events");

	//hFidCutXvsFidCutYvsMeanChargeAllDetectors
	hFidCutXvsFidCutYvsMeanChargeAllDetectors = (TH2D*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("hFidCutXvsFidCutYvsMeanChargeAllDetectors");

	//hXPosvsYPosvsMeanCharge		For TH2D
	stringstream hXPosvsYPosvsMeanChargeName; hXPosvsYPosvsMeanChargeName<<"hXPosvsYPosvsMeanCharge";
	hXPosvsYPosvsMeanCharge = new TH2D(hXPosvsYPosvsMeanChargeName.str().c_str(),hXPosvsYPosvsMeanChargeName.str().c_str(),213,5800,9600,160,3400,5900);
	hXPosvsYPosvsMeanCharge->GetXaxis()->SetTitle("XPos");
	hXPosvsYPosvsMeanCharge->GetYaxis()->SetTitle("YPos");
	hXPosvsYPosvsMeanCharge->GetZaxis()->SetTitle("Charge ADC");
	hXPosvsYPosvsCharge = (TH2D*)hXPosvsYPosvsMeanCharge->Clone("hXPosvsYPosvsCharge");
	hXPosvsYPosvsEvents = (TH2D*)hXPosvsYPosvsMeanCharge->Clone("hXPosvsYPosvsEvents");
	//hFidCutXvsFidCutYvsCharge
	/*stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%";   //<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%";
	//hFidCutXvsFidCutYvsCharge.push_back(new TH3F(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),2800,90,170,2800,60,120,2800,0,3000));
	hFidCutXvsFidCutYvsCharge = new TH3F(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),160,90,170,120,60,120,30,0,3000);
	hFidCutXvsFidCutYvsCharge->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsCharge->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsCharge->GetZaxis()->SetTitle("Charge ADC");
	*/
}

void TAnalysisOf3dDiamonds::saveHistos() {

	//hNumberofClusters
	histSaver->SaveHistogram(hNumberofClusters);
	histSaver->SaveHistogram(hDoubleClusterPos);
	histSaver->SaveHistogram(hDoubleClusterPos0);
	histSaver->SaveHistogram(hDoubleClusterPos1);
	histSaver->SaveHistogram(hLandauCluster1);
	histSaver->SaveHistogram(hLandauCluster2);
	histSaver->SaveHistogram(hLandauDoubleCombined);

	DoubleCluster = new TCanvas();
	DoubleCluster->cd();
	hDoubleClusterPos0->Draw();
	hDoubleClusterPos1->Draw("samehist");
	stringstream str111; str111<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hDoubleClusterCombinedPos%%"<<FileNameEnd<<".png";
	DoubleCluster->SaveAs(str111.str().c_str());


	for(int i=0;i<Detector.size();i++){
		/*//hLandau
		stringstream hLandauName; hLandauName<<"hLandau%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hLandau.push_back(HistogrammSaver::CreateDistributionHisto(hLandauName.str().c_str(),*vecPHDiamondHit.at(i),256));
		hLandau.at(i)->SetTitle(hLandauName.str().c_str());
		hLandau.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
		hLandau.at(i)->GetYaxis()->SetTitle("number of entries #");
		*/
		histSaver->SaveHistogram(hLandau.at(i));

		//hPredictedPositionDiamondHit
		stringstream hPredictedPositionDiamondHitName; hPredictedPositionDiamondHitName<<"hPredictedPositionDiamondHit%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd;
		hPredictedPositionDiamondHit.push_back(HistogrammSaver::CreateScatterHisto(hPredictedPositionDiamondHitName.str().c_str(),*vecXPredicted.at(i),*vecYPredicted.at(i),512));
		hPredictedPositionDiamondHit.at(i)->SetTitle(hPredictedPositionDiamondHitName.str().c_str());
		hPredictedPositionDiamondHit.at(i)->GetXaxis()->SetTitle("predicted x Position");
		hPredictedPositionDiamondHit.at(i)->GetYaxis()->SetTitle("predicted y Position");
		histSaver->SaveHistogram(hPredictedPositionDiamondHit.at(i));

		histSaver->SaveHistogram(hHitandSeedCount.at(i));
		histSaver->SaveHistogram(hEventsvsChannel.at(i));
		histSaver->SaveHistogram(hPHvsChannel.at(i));
		histSaver->SaveHistogram(hPHvsPredictedXPos.at(i));
		histSaver->SaveHistogram(hChi2XChi2Y.at(i));
		histSaver->SaveHistogram(hFidCutXvsFidCutY.at(i));
		histSaver->SaveHistogram(hPHvsPredictedChannel.at(i));
		//histSaver->SaveHistogram(hFidCutXvsFidCutYvsCharge.at(i));

		//hFidCutXvsFidCutYvsCharge
		ptrCanvas.at(i)->cd();
		hFidCutXvsFidCutYvsCharge.at(i)->Draw("COLZ");
		stringstream str; str<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsCharge%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd<<".png";
		ptrCanvas.at(i)->SaveAs(str.str().c_str());

		//hFidCutXvsFidCutYvsEvents
		ptrCanvasEvents.at(i)->cd();
		hFidCutXvsFidCutYvsEvents.at(i)->Draw("COLZ");
		stringstream str1; str1<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsChargeClone%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd<<".png";
		ptrCanvasEvents.at(i)->SaveAs(str1.str().c_str());

		//hFidCutXvsFidCutYvsMeanCharge
		ptrCanvasMean.at(i)->cd();
		*hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
		hFidCutXvsFidCutYvsMeanCharge.at(i)->Draw("COLZ");
		stringstream str2; str2<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanCharge%%"<<Detector.at(i)[0]<<"-"<<Detector.at(i)[1]<<"%%"<<FileNameEnd<<".png";
		ptrCanvasMean.at(i)->SaveAs(str2.str().c_str());

	} //End of for loop
	//For all Diamond hFidCutXvsFidCutYvsMeanCharge
	hCombinedMeanCharge = new TCanvas();
	hCombinedMeanCharge->cd();
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(0));
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(1));
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(2));
	hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
	/*DrawFidCutRegions(); //Draw Fiducial Cut Regions
	 	 */
	//FidCudBoundMetric->Draw("same"); //To draw the fiducial cut regions
	stringstream str1; str1<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn"<<"%%"<<FileNameEnd<<".png";
	hCombinedMeanCharge->SaveAs(str1.str().c_str());

	DrawFidCutRegions(); //Draw Fiducial Cut Regions
	stringstream str11; str11<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanChargeAllDetectors"<<"%%"<<FileNameEnd<<".png";
	hCombinedMeanCharge->SaveAs(str11.str().c_str());

	//For efficiency
	//Predicted Events
	hPredictedEvents = new TCanvas();
	hPredictedEvents->cd();
	hFidCutXvsFidCutYvsPredictedEvents->Draw("COLZ");
	stringstream str; str<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsPredictedEvents"<<"%%"<<FileNameEnd<<".png";
	hPredictedEvents->SaveAs(str.str().c_str());
	//Seen Events
	hSeenEvents = new TCanvas();
	hSeenEvents->cd();
	hFidCutXvsFidCutYvsSeenEvents->Draw("COLZ");
	stringstream str3; str3<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsSeenEvents"<<"%%"<<FileNameEnd<<".png";
	hSeenEvents->SaveAs(str3.str().c_str());
	//Efficiency
	TCanvas* hEfficiencyMap = new TCanvas();
	*hEfficiency = (*hFidCutXvsFidCutYvsSeenEvents/(*hFidCutXvsFidCutYvsPredictedEvents));
	hEfficiency->Draw("COLZ");
	stringstream str4; str4<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hEfficiencyMap"<<"%%"<<FileNameEnd<<".png";
	hEfficiencyMap->SaveAs(str4.str().c_str());

	//0Clusters
	h0Clusters = new TCanvas();
	h0Clusters->cd();
	hFidCutXvsFidCutY0Clusters->Draw("COLZ");
	stringstream str7; str7<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY0Clusters"<<"%%"<<FileNameEnd<<".png";
	h0Clusters->SaveAs(str7.str().c_str());

	//1_1Clusters
	h1_1Clusters = new TCanvas();
	h1_1Clusters->cd();
	hFidCutXvsFidCutY1_1Clusters->Draw("COLZ");
	stringstream str8; str8<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY1_1Clusters"<<"%%"<<FileNameEnd<<".png";
	h1_1Clusters->SaveAs(str8.str().c_str());

	//1Clusters
	h1Clusters = new TCanvas();
	h1Clusters->cd();
	hFidCutXvsFidCutY1Clusters->Draw("COLZ");
	stringstream str9; str9<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY1Clusters"<<"%%"<<FileNameEnd<<".png";
	h1Clusters->SaveAs(str9.str().c_str());

	//2Clusters
	h2Clusters = new TCanvas();
	h2Clusters->cd();
	hFidCutXvsFidCutY2Clusters->Draw("COLZ");
	stringstream str5; str5<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY2Clusters"<<"%%"<<FileNameEnd<<".png";
	h2Clusters->SaveAs(str5.str().c_str());

	//2Clusters0
	h2Clusters0 = new TCanvas();
	h2Clusters0->cd();
	hFidCutXvsFidCutY2Clusters0->Draw("COLZ");
	stringstream str50; str50<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY2Clusters0"<<"%%"<<FileNameEnd<<".png";
	h2Clusters0->SaveAs(str50.str().c_str());

	//2Clusters1
	h2Clusters1 = new TCanvas();
	h2Clusters1->cd();
	hFidCutXvsFidCutY2Clusters1->Draw("COLZ");
	stringstream str51; str51<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY2Clusters1"<<"%%"<<FileNameEnd<<".png";
	h2Clusters1->SaveAs(str51.str().c_str());

	//3Clusters
	h3Clusters = new TCanvas();
	h3Clusters->cd();
	hFidCutXvsFidCutY3Clusters->Draw("COLZ");
	stringstream str6; str6<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutY3Clusters"<<"%%"<<FileNameEnd<<".png";
	h3Clusters->SaveAs(str6.str().c_str());

	//For hXPosvsYPosvsMeanCharge
	TCanvas* XPosvsYPosvsMean = new TCanvas();
	XPosvsYPosvsMean->cd();
	*hXPosvsYPosvsMeanCharge = (*hXPosvsYPosvsCharge/(*hXPosvsYPosvsEvents));
	hXPosvsYPosvsMeanCharge->Draw("COLZ");
	stringstream str2; str2<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"XPosvsYPosvsMeanCharge"<<"%%"<<FileNameEnd<<".png";
	XPosvsYPosvsMean->SaveAs(str2.str().c_str());
}

void TAnalysisOf3dDiamonds::analyseEvent() {

	if(!eventReader->isValidTrack()) return;
	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	if(!eventReader->isInFiducialCut())	//This is a larger fiducial cut around silicon
		return;

	/*if(eventReader->getNDiamondClusters()==0||eventReader->getNDiamondClusters()==2||eventReader->getNDiamondClusters()==3)
		return;
		*/

	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster())
		return;

	TPositionPrediction *predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	if(predictedPosition->getChi2X()>100||predictedPosition->getChi2Y()>100)
		return;

	Float_t xPos = predictedPosition->getPositionX();
	Float_t yPos = predictedPosition->getPositionY();

	//Predict Channel hit
	Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
	Float_t positionInDetSystemChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,positionInDetSystemMetric);
	//Float_t CorrectedpositionInDetSystemChannelSpace = positionInDetSystemChannelSpace+0.5;
	//cout<<"The pridicted hit channel is: "<<positionInDetSystemChannelSpace<<endl;

	Float_t Chi2X = predictedPosition->getChi2X();
	Float_t Chi2Y = predictedPosition->getChi2Y();
	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();

	HitandSeedCount(&diamondCluster, 0);
	hNumberofClusters->Fill(eventReader->getNDiamondClusters());

	if(eventReader->getNDiamondClusters()==1){
		if(HitCount==0&&SeedCount==1){
			//cout<<"1 and 1Seed 0Hits"<<endl;
			hFidCutXvsFidCutY1_1Clusters->Fill(fiducialValueX,fiducialValueY,1);
		}
	}
	if(eventReader->getNDiamondClusters()==1){
			//cout<<"1"<<endl;
			hFidCutXvsFidCutY1Clusters->Fill(fiducialValueX,fiducialValueY,1);
	}
	if(eventReader->getNDiamondClusters()==0){
		//cout<<"0"<<endl;
		hFidCutXvsFidCutY0Clusters->Fill(fiducialValueX,fiducialValueY,1);
	}
	if(eventReader->getNDiamondClusters()==2){
		cout<<"2"<<endl;
		hFidCutXvsFidCutY2Clusters->Fill(fiducialValueX,fiducialValueY,1);
		TCluster diamondCluster0 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
		TCluster diamondCluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
		//TCluster diamondCluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),2);
		cout<<"Cl 0 Highest Channel: "<<diamondCluster0.getHighestSignalChannel()<<endl;
		cout<<"Cl 1 Highest Channel: "<<diamondCluster1.getHighestSignalChannel()<<endl;
		//cout<<"Cl 2 Highest Channel: "<<diamondCluster2.getHighestSignalChannel()<<endl;

		hDoubleClusterPos->Fill(diamondCluster0.getHighestSignalChannel());
		hDoubleClusterPos->Fill(diamondCluster1.getHighestSignalChannel());
		if(diamondCluster0.getHighestSignalChannel()==85||diamondCluster1.getHighestSignalChannel()==85){
			hDoubleClusterPos0->Fill(diamondCluster0.getHighestSignalChannel());
			hDoubleClusterPos0->Fill(diamondCluster1.getHighestSignalChannel());
			hLandauCluster1->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
			hFidCutXvsFidCutY2Clusters0->Fill(fiducialValueX,fiducialValueY,1);
		}
		if(diamondCluster0.getHighestSignalChannel()==55||diamondCluster1.getHighestSignalChannel()==55){
			hDoubleClusterPos1->Fill(diamondCluster0.getHighestSignalChannel());
			hDoubleClusterPos1->Fill(diamondCluster1.getHighestSignalChannel());
			hLandauCluster2->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
			hFidCutXvsFidCutY2Clusters1->Fill(fiducialValueX,fiducialValueY,1);
		}
		//hLandauCluster1->Fill(diamondCluster0.getCharge(false));
		//hLandauCluster2->Fill(diamondCluster1.getCharge(false));
		if((!diamondCluster0.getHighestSignalChannel()==55&&!diamondCluster1.getHighestSignalChannel()==55)||(!diamondCluster0.getHighestSignalChannel()==85&&!diamondCluster1.getHighestSignalChannel()==85)){
			hLandauDoubleCombined->Fill((diamondCluster0.getCharge(false)+diamondCluster1.getCharge(false)));
		}

	}
	if(eventReader->getNDiamondClusters()==3){
		//cout<<"3"<<endl;
		hFidCutXvsFidCutY3Clusters->Fill(fiducialValueX,fiducialValueY,1);
	}

	/*if(!eventReader->getNDiamondClusters()==1)
		return;
		*/

	//RemoveLumpyClusters(&diamondCluster);

	//Universal PHvsChannel Plot
	for(int i=0; i<Detector.size(); i++){
		//RemoveEdgeHits(&diamondCluster, Detector.at(i));
		//hFidCutXvsFidCutYvsCharge->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
		if(eventReader->getNDiamondClusters()==1){

			if(diamondCluster.getHighestSignalChannel()<=Detector.at(i)[1]&&diamondCluster.getHighestSignalChannel()>=Detector.at(i)[0]){

				if(FiducialCut(i)==1)		//Cut on fiducial cuts
					return;

				/*HitandSeedCount(&diamondCluster, i);
				if(HitCount>0||SeedCount>1)
					return;
					*/
				hEventsvsChannel.at(i)->Fill(diamondCluster.getHighestSignalChannel());
				hPHvsChannel.at(i)->Fill(diamondCluster.getCharge(false),diamondCluster.getHighestSignalChannel());

				hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);

				hPHvsPredictedXPos.at(i)->Fill(diamondCluster.getCharge(false),xPos);
				hLandau.at(i)->Fill(diamondCluster.getCharge(false));
				vecPHDiamondHit.at(i)->push_back(diamondCluster.getCharge(false));
				vecXPredicted.at(i)->push_back(xPos);
				vecYPredicted.at(i)->push_back(yPos);
				HitandSeedCount(&diamondCluster, i);
				hHitandSeedCount.at(i)->Fill(HitCount,SeedCount);
				hChi2XChi2Y.at(i)->Fill(Chi2X, Chi2Y);
				hFidCutXvsFidCutY.at(i)->Fill(fiducialValueX,fiducialValueY);
				//For hFidCutXvsFidCutYvsMeanCharge
				hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
				hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);
				hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);

				//For XPosvsYPosvsMeanCharge
				//cout<<"XPos = "<<xPos<<"   YPos = "<<yPos<<endl;
				hXPosvsYPosvsCharge->Fill(xPos,yPos,diamondCluster.getCharge(false));
				//hXPosvsYPosvsMeanCharge->Fill(xPos,yPos,diamondCluster.getCharge(false));
				hXPosvsYPosvsEvents->Fill(xPos,yPos,1);
			}        			//if(diamondCluster.getHighestSignalChannel()<=Detector.at(i)[1]&&diamondCluster.getHighestSignalChannel()>=Detector.at(i)[0]){

		}     		//if(eventReader->getNDiamondClusters()==1){

		// Part of loop for Transparent analysis
		// Calculate Detector efficiency
		if(eventReader->getNDiamondClusters()==1||eventReader->getNDiamondClusters()==0){
			if(positionInDetSystemChannelSpace<=(Detector.at(i)[1]+0.5)&&positionInDetSystemChannelSpace>=(Detector.at(i)[0]-0.5)){

				if(FiducialCut(i)==1)		//Cut on fiducial cuts
					return;

				hFidCutXvsFidCutYvsPredictedEvents->Fill(fiducialValueX,fiducialValueY,1);

				/*hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);
				hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
				*/
				//hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);

			}
		}     		//if(eventReader->getNDiamondClusters()==1||eventReader->getNDiamondClusters()==0){

	}

}
void TAnalysisOf3dDiamonds::initialiseYAlignmentHistos() {

	//hGridReference
	vector<char> GridReference;
	GridReference.push_back('A');GridReference.push_back('B');GridReference.push_back('C');
	GridReference.push_back('D');GridReference.push_back('E');GridReference.push_back('F');
	GridReference.push_back('G');GridReference.push_back('H');GridReference.push_back('I');
	int GridReferenceY[] = {1,2,3,4,5,6,7,8,9,10,11};
	stringstream hGridReferenceName; hGridReferenceName<<""<<FileNameEnd;
	hGridReference = new TH2D(hGridReferenceName.str().c_str(),hGridReferenceName.str().c_str(),9,2365,3715,11,0,1650);
	for(int i=0;i<9;i++){
		stringstream iLetter; iLetter<<GridReference.at(i);
		hGridReference->GetXaxis()->SetBinLabel(i+1,iLetter.str().c_str());
	}
	for(int j=0;j<11;j++){
		stringstream jNumber; jNumber<<GridReferenceY[j];
		hGridReference->GetYaxis()->SetBinLabel(j+1,jNumber.str().c_str());
	}
	hGridReference->SetStats(kFALSE);
	hGridReference->SetTickLength(0.0, "X");
	hGridReference->SetTickLength(0.0, "Y");

	//hFidCutXvsFidCutYvsChargeYAlignment
	stringstream hFidCutXvsFidCutYvsChargeYAlignmentName; hFidCutXvsFidCutYvsChargeYAlignmentName<<"hFidCutXvsFidCutYvsChargeYAlignment"<<FileNameEnd;
	hFidCutXvsFidCutYvsChargeYAlignment = new TH2D(hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutYvsChargeYAlignment->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsChargeYAlignment->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsChargeYAlignment->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignment
	hFidCutXvsFidCutYvsEventsYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsEventsYAlignment");

	//hFidCutXvsFidCutYvsMeanChargeYAlignment
	hFidCutXvsFidCutYvsMeanChargeYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsMeanChargeYAlignment");

	//hDetXvsDetY3D
	stringstream hDetXvsDetY3DName; hDetXvsDetY3DName<<"hFidCutXvsFidCutYvsChargeYAlignment"<<FileNameEnd;
	hDetXvsDetY3D = new TH2D(hDetXvsDetY3DName.str().c_str(),hDetXvsDetY3DName.str().c_str(),270,2365,3715,330,0,1650);
	hDetXvsDetY3D->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3D->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3D->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignment
	hDetXvsDetY3DvsEvents = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DvsEvents");

	//hFidCutXvsFidCutYvsMeanChargeYAlignment
	hDetXvsDetY3DMeanCharge = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DMeanCharge");

	//hCellsMeanClusteSize
	stringstream hCellsMeanClusteSizeName; hCellsMeanClusteSizeName<<"hCellsMeanClusteSize"<<FileNameEnd;
	hCellsMeanClusteSize = new TH2D(hCellsMeanClusteSizeName.str().c_str(),hCellsMeanClusteSizeName.str().c_str(),9,2365,3715,11,0,1650);

	//hDetXvsDetY3DRebinned
	stringstream hDetXvsDetY3DRebinnedName; hDetXvsDetY3DRebinnedName<<"hFidCutXvsFidCutYvsChargeRebinnedYAlignment"<<FileNameEnd;
	hDetXvsDetY3DRebinned = new TH2D(hDetXvsDetY3DRebinnedName.str().c_str(),hDetXvsDetY3DRebinnedName.str().c_str(),9,2365,3715,11,0,1650);
	hDetXvsDetY3DRebinned->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinned->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinned->GetZaxis()->SetTitle("Charge ADC");

	//RebinnedQuarterCellFails
	stringstream RebinnedQuarterCellFailsName; RebinnedQuarterCellFailsName<<"3DdetNumberofQuarterCellFails"<<FileNameEnd;
	RebinnedQuarterCellFails = new TH2D(RebinnedQuarterCellFailsName.str().c_str(),RebinnedQuarterCellFailsName.str().c_str(),9,2365,3715,11,0,1650);
	RebinnedQuarterCellFails->GetXaxis()->SetTitle("Xdet (um)");
	RebinnedQuarterCellFails->GetYaxis()->SetTitle("Ydet (um)");
	RebinnedQuarterCellFails->GetZaxis()->SetTitle("Quarter Fails");

	//hDetXvsDetY3DRebinnedRMS
	stringstream hDetXvsDetY3DRebinnedRMSName; hDetXvsDetY3DRebinnedRMSName<<"h3DdetRebinnedRMS"<<FileNameEnd;
	hDetXvsDetY3DRebinnedRMS = new TH2D(hDetXvsDetY3DRebinnedRMSName.str().c_str(),hDetXvsDetY3DRebinnedRMSName.str().c_str(),9,2365,3715,11,0,1650);
	hDetXvsDetY3DRebinnedRMS->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinnedRMS->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinnedRMS->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DMeanChargeRebinnedQuarterCell
	stringstream hDetXvsDetY3DMeanChargeRebinnedQuarterCellName; hDetXvsDetY3DMeanChargeRebinnedQuarterCellName<<"h3DdetQuarterCellMeanCharge"<<FileNameEnd;
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell = new TH2D(hDetXvsDetY3DMeanChargeRebinnedQuarterCellName.str().c_str(),hDetXvsDetY3DMeanChargeRebinnedQuarterCellName.str().c_str(),18,2365,3715,22,0,1650);
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DRebinnedQuarterCellRMS
	stringstream hDetXvsDetY3DRebinnedQuarterCellRMSName; hDetXvsDetY3DRebinnedQuarterCellRMSName<<"h3DdetQuarterCellRMS"<<FileNameEnd;
	hDetXvsDetY3DRebinnedQuarterCellRMS = new TH2D(hDetXvsDetY3DRebinnedQuarterCellRMSName.str().c_str(),hDetXvsDetY3DRebinnedQuarterCellRMSName.str().c_str(),18,2365,3715,22,0,1650);
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DMeanChargeQuarterCellGrading
	for(int k=0; k<6; k++){
		stringstream hDetXvsDetY3DMeanChargeQuarterCellGradingName; hDetXvsDetY3DMeanChargeQuarterCellGradingName<<"hDetXvsDetY3DMeanChargeQuarterCellGrading"<<k<<"%%Fail"<<FileNameEnd;
		hDetXvsDetY3DMeanChargeQuarterCellGrading.push_back(new TH2D(hDetXvsDetY3DMeanChargeQuarterCellGradingName.str().c_str(),hDetXvsDetY3DMeanChargeQuarterCellGradingName.str().c_str(),18,2365,3715,22,0,1650));
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetXaxis()->SetTitle("Xdet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetYaxis()->SetTitle("Ydet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetZaxis()->SetTitle("Charge ADC");
	}

	//hDetXvsDetY3DMeanChargeQuarterCellGradingLandau
	for(int k=0;k<5;k++){

		//For Transparent analysis.
		//hQuarterCellGradedTransparentLandau
		stringstream hQuarterCellGradedTransparentLandauName; hQuarterCellGradedTransparentLandauName<<"hQuarterCellGradedTransparentLandau"<<k<<FileNameEnd;
		hQuarterCellGradedTransparentLandau.push_back(new TH1F(hQuarterCellGradedTransparentLandauName.str().c_str(),hQuarterCellGradedTransparentLandauName.str().c_str(),256,0,2800));

		stringstream hDetXvsDetY3DMeanChargeQuarterCellGradingLandauName; hDetXvsDetY3DMeanChargeQuarterCellGradingLandauName<<"hDetXvsDetY3DMeanChargeQuarterCellGradingLandau"<<k<<FileNameEnd;
		hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.push_back(new TH1F(hDetXvsDetY3DMeanChargeQuarterCellGradingLandauName.str().c_str(),hDetXvsDetY3DMeanChargeQuarterCellGradingLandauName.str().c_str(),256,0,2800));

		//hCellsDeltaXQuarterCellGrading
		stringstream hCellsDeltaXQuarterCellGradingName; hCellsDeltaXQuarterCellGradingName<<"hCellsDeltaXQuarterCellGrading"<<k<<FileNameEnd;
		hCellsDeltaXQuarterCellGrading.push_back(new TH1F(hCellsDeltaXQuarterCellGradingName.str().c_str(),hCellsDeltaXQuarterCellGradingName.str().c_str(),100,-3,3));
	}

	//hFidCutXvsFidCutYvsEventsYAlignmentRebinned
	hDetXvsDetY3DvsEventsRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DvsEventsRebinned");

	//hFidCutXvsFidCutYvsMeanChargeYAlignmentRebinned
	hDetXvsDetY3DMeanChargeRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DMeanChargeRebinned");
	//hDetXvsDetY3DMeanChargeRebinned->SetBins(9,2365,3715,11,0,1650,12,0,1200);

	//hDetXvsDetY3DOverview
	stringstream hDetXvsDetY3DOverviewName; hDetXvsDetY3DOverviewName<<"hDetXvsDetY3DOverview"<<FileNameEnd;
	hDetXvsDetY3DOverview = new TH2D(hDetXvsDetY3DOverviewName.str().c_str(),hDetXvsDetY3DOverviewName.str().c_str(),9,2365,3715,11,0,1650);
	hDetXvsDetY3DOverview->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DOverview->GetYaxis()->SetTitle("Ydet (um)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();

	//hCellNumbering
	stringstream hCellNumberingName; hCellNumberingName<<"h3DdetCellNumbering"<<FileNameEnd;
	hCellNumbering = new TH2D(hCellNumberingName.str().c_str(),hCellNumberingName.str().c_str(),9,2365,3715,11,0,1650);
	hCellNumbering->GetXaxis()->SetTitle("Xdet (um)");
	hCellNumbering->GetYaxis()->SetTitle("Ydet (um)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();

	//hBinnedMeanCharge
	stringstream hBinnedMeanChargeName; hBinnedMeanChargeName<<"h3DdetCellMeanChargeBinned"<<FileNameEnd;
	hBinnedMeanCharge = new TH1F(hBinnedMeanChargeName.str().c_str(),hBinnedMeanChargeName.str().c_str(),9,400,1300);
	hBinnedMeanCharge->SetTitle(hBinnedMeanChargeName.str().c_str());
	hBinnedMeanCharge->GetXaxis()->SetTitle("MeanCharge");
	hBinnedMeanCharge->GetYaxis()->SetTitle("Entries");

	//h3DdetDeltaXChannel
	stringstream h3DdetDeltaXChannelName; h3DdetDeltaXChannelName<<"h3DdetDeltaXChannel"<<FileNameEnd;
	h3DdetDeltaXChannel = new TH2D(h3DdetDeltaXChannelName.str().c_str(),h3DdetDeltaXChannelName.str().c_str(),270,2365,3715,330,0,1650);
	h3DdetDeltaXChannel->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetDeltaXChannel->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetDeltaXChannel->GetZaxis()->SetTitle("Delta X (ch)");

	//h3DdetDeltaXChannelAbove1000
	stringstream h3DdetDeltaXChannelAbove1000Name; h3DdetDeltaXChannelAbove1000Name<<"h3DdetDeltaXChannelAbove1000"<<FileNameEnd;
	h3DdetDeltaXChannelAbove1000 = new TH2D(h3DdetDeltaXChannelAbove1000Name.str().c_str(),h3DdetDeltaXChannelAbove1000Name.str().c_str(),270,2365,3715,330,0,1650);
	h3DdetDeltaXChannelAbove1000->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetDeltaXChannelAbove1000->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetDeltaXChannelAbove1000->GetZaxis()->SetTitle("Delta X (ch)");

	//Cell Overlay
	//hCellsOverlayed
	hCellsOverlayedCharge = new TH2D("OverlayedCharge","OverlayedCharge",30,0,150,30,0,150);        //30,0,150,30,0,150);
	hCellsOverlayedEvents = new TH2D("OverlayedEvents","OverlayedEvents",30,0,150,30,0,150);        //30,0,150,30,0,150);
	hCellsOverlayedEventsNoColumns = new TH2D("hCellsOverlayedEventsNoColumns","hCellsOverlayedEventsNoColumns",15,0,150,15,0,150);        //30,0,150,30,0,150);
	stringstream hCellsOverlayedMeanChargeName; hCellsOverlayedMeanChargeName<<"hCellsOverlayedMeanCharge"<<FileNameEnd;
	hCellsOverlayedMeanCharge = new TH2D(hCellsOverlayedMeanChargeName.str().c_str(),hCellsOverlayedMeanChargeName.str().c_str(),30,0,150,30,0,150);        //30,0,150,30,0,150);
	hCellsOverlayedMeanCharge->GetXaxis()->SetTitle("Xdet (um)");
	hCellsOverlayedMeanCharge->GetYaxis()->SetTitle("Ydet (um)");
	hCellsOverlayedMeanCharge->GetZaxis()->SetTitle("Charge ADC");
	hCellsOverlayedMeanCharge->GetZaxis()->SetRangeUser(800,1200);
	hCellsOverlayedMeanCharge->SetContour(99);

	for(int i=0;i<9;i++){
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

	//Define Landau function for Landau fit.
	Landau = new TF1("Landau","landau(0)",20,80);
	//
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			/*
			float xLow = 2365 + i*150;
			float yLow = j*150;
			float xHigh = xLow+150;
			float yHigh = yLow+150;
				*/
			//For Transparent analysis.
			//hCellTransparentLandau
			stringstream hCellTransparentLandauName; hCellTransparentLandauName<<"hCellTransparentLandau"<<((i*11)+j)<<FileNameEnd;
			hCellTransparentLandau.push_back(new TH1F(hCellTransparentLandauName.str().c_str(),hCellTransparentLandauName.str().c_str(),256,0,2800));
			//hCellsTransparentHitPosition
			stringstream hCellsTransparentHitPositionName; hCellsTransparentHitPositionName<<"hCellsTransparentHitPosition"<<((i*11)+j)<<FileNameEnd;
			hCellsTransparentHitPosition.push_back(new TH2D(hCellsTransparentHitPositionName.str().c_str(),hCellsTransparentHitPositionName.str().c_str(),30,0,150,30,0,150));

			stringstream hCellsChargeName; hCellsChargeName<<"TotalCharge"<<((i*11)+j)<<FileNameEnd;
			hCellsCharge.push_back(new TH2D(hCellsChargeName.str().c_str(),hCellsChargeName.str().c_str(),30,0,150,30,0,150));        //30,0,150,30,0,150);
			stringstream hCellsEventsName; hCellsEventsName<<"Events"<<((i*11)+j)<<FileNameEnd;
			hCellsEvents.push_back(new TH2D(hCellsEventsName.str().c_str(),hCellsEventsName.str().c_str(),30,0,150,30,0,150));			//30,0,150,30,0,150);

			//hCellsEventsCheck
			stringstream hCellsEventsCheckName; hCellsEventsCheckName<<"hCellsEventsCheck"<<((i*11)+j)<<FileNameEnd;
			hCellsEventsCheck.push_back(new TH2D(hCellsEventsCheckName.str().c_str(),hCellsEventsCheckName.str().c_str(),30,0,150,30,0,150));

			stringstream hCellsLandauName; hCellsLandauName<<"hLandauCell"<<((i*11)+j)<<FileNameEnd;
			hCellsLandau.push_back(new TH1F(hCellsLandauName.str().c_str(),hCellsLandauName.str().c_str(),256,0,2800));

			stringstream hCellsClusteSizeName; hCellsClusteSizeName<<"hCellsClusteSize"<<((i*11)+j)<<FileNameEnd;
			hCellsClusteSize.push_back(new TH1F(hCellsClusteSizeName.str().c_str(),hCellsClusteSizeName.str().c_str(),20,0,20));

			//hCellsDeltaX
			stringstream hCellsDeltaXName; hCellsDeltaXName<<"hCellsDeltaX"<<((i*11)+j)<<FileNameEnd;
			hCellsDeltaX.push_back(new TH1F(hCellsDeltaXName.str().c_str(),hCellsDeltaXName.str().c_str(),100,-3,3));

			stringstream hCellsLandauNoColumnName; hCellsLandauNoColumnName<<"hCellsLandauNoColumn"<<((i*11)+j)<<FileNameEnd;
			hCellsLandauNoColumn.push_back(new TH1F(hCellsLandauNoColumnName.str().c_str(),hCellsLandauNoColumnName.str().c_str(),256,0,2800));

			stringstream hCellsEventsNoColumnName; hCellsEventsNoColumnName<<"hCellsEventsNoColumn"<<((i*11)+j)<<FileNameEnd;
			hCellsEventsNoColumn.push_back(new TH2D(hCellsEventsNoColumnName.str().c_str(),hCellsEventsNoColumnName.str().c_str(),15,0,150,15,0,150));		//30,0,150,30,0,150);

			for(int k=0;k<4;k++){
				//cout<<"This should not repeat: "<<((i*11*4)+j*4+k)<<endl;
				stringstream hQuaterCellsLandauName; hQuaterCellsLandauName<<"hQuaterCellsLandau"<<((i*11*4)+j*4+k)<<FileNameEnd;
				hQuaterCellsLandau.push_back(new TH1F(hQuaterCellsLandauName.str().c_str(),hQuaterCellsLandauName.str().c_str(),256,0,2800));
			}
		}
	}
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

	//cout<<"The size of hQuarterCells is: "<<hQuaterCellsLandau.size()<<endl;

	//Transparent Analysis
	//hCellsTransparentHitPositionCellGraded0
	stringstream hCellsTransparentHitPositionCellGraded0Name; hCellsTransparentHitPositionCellGraded0Name<<"hCellsTransparentHitPositionCellGraded0"<<FileNameEnd;
	hCellsTransparentHitPositionCellGraded0 = new TH2D(hCellsTransparentHitPositionCellGraded0Name.str().c_str(),hCellsTransparentHitPositionCellGraded0Name.str().c_str(),30,0,150,30,0,150);

	//hCellsOverlayedLandauNoColumn
	stringstream hCellsOverlayedEntriesNoColumnsName; hCellsOverlayedEntriesNoColumnsName<<"hCellsOverlayedEntriesNoColumns"<<FileNameEnd;
	hCellsOverlayedEntriesNoColumns = new TH1F(hCellsOverlayedEntriesNoColumnsName.str().c_str(),hCellsOverlayedEntriesNoColumnsName.str().c_str(),100,0,100);

	stringstream hCellsOverlayedLandauNoColumnName; hCellsOverlayedLandauNoColumnName<<"hCellsOverlayedLandauNoColumn"<<FileNameEnd;
	hCellsOverlayedLandauNoColumn = new TH1F(hCellsOverlayedLandauNoColumnName.str().c_str(),hCellsOverlayedLandauNoColumnName.str().c_str(),256,0,2800);

	//hCellsLandauGraded    &&    hCellsLandauGradedNoColumn
	for(int i=0;i<12;i++){	//Group CellsLandaus within same ranges together. 0-100; 100-200; -> 1100-1200;
		stringstream hCellsLandauGradedName; hCellsLandauGradedName<<"hLandauCellsGraded"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGraded.push_back(new TH1F(hCellsLandauGradedName.str().c_str(),hCellsLandauGradedName.str().c_str(),256,0,2800));
		stringstream hCellsLandauGradedNoColumnName; hCellsLandauGradedNoColumnName<<"hLandauCellsGradedNoColumn%%"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGradedNoColumn.push_back(new TH1F(hCellsLandauGradedNoColumnName.str().c_str(),hCellsLandauGradedNoColumnName.str().c_str(),256,0,2800));
	}

	//hCellsGoodandBad
	stringstream hCellsHarris18GoodName; hCellsHarris18GoodName<<"hCellsHarris18Good"<<FileNameEnd;
	hCellsHarris18Good = new TH1F(hCellsHarris18GoodName.str().c_str(),hCellsHarris18GoodName.str().c_str(),256,0,2800);
	hCellsHarris18Good->SetTitle(hCellsHarris18GoodName.str().c_str());
	hCellsHarris18Good->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hCellsHarris18Good->GetYaxis()->SetTitle("Entries");

	hCellsHarris10Bad = (TH1F*)hCellsHarris18Good->Clone("hCellsHarris10Bad");

	//hCellsOverlayedColumnLandau
	stringstream hCellsOverlayedColumnLandauName; hCellsOverlayedColumnLandauName<<"hCellsOverlayedColumnLandau"<<FileNameEnd;
	hCellsOverlayedColumnLandau = new TH1F(hCellsOverlayedColumnLandauName.str().c_str(),hCellsOverlayedColumnLandauName.str().c_str(),256,0,2800);
	hCellsOverlayedColumnLandau->SetTitle(hCellsOverlayedColumnLandauName.str().c_str());
	hCellsOverlayedColumnLandau->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hCellsOverlayedColumnLandau->GetYaxis()->SetTitle("Entries");

	//hTransparentCharge3D
	stringstream hTransparentCharge3DName; hTransparentCharge3DName<<"hTransparentCharge3D"<<FileNameEnd;
	hTransparentCharge3D = new TH1F(hTransparentCharge3DName.str().c_str(),hTransparentCharge3DName.str().c_str(),256,0,2800);
	hTransparentCharge3D->SetTitle(hTransparentCharge3DName.str().c_str());
	hTransparentCharge3D->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hTransparentCharge3D->GetYaxis()->SetTitle("Entries");

	//hCellsLandau2D
	stringstream hCellsLandau2DName; hCellsLandau2DName<<"hCellsLandau2D"<<FileNameEnd;
	hCellsLandau2D = new TH2D(hCellsLandau2DName.str().c_str(),hCellsLandau2DName.str().c_str(),256,0,2800,99,0,99);
	hCellsLandau2D->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2D->GetYaxis()->SetTitle("Cell");
	hCellsLandau2DQuarterFail = new TH2D("hCellsLandau2DQuarterFail","hCellsLandau2DQuarterFail",1,0,2800,99,0,99);
	hCellsLandau2DQuarterFail->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2DQuarterFail->GetYaxis()->SetTitle("Cell");
	//hCellsLandau2D->SetCanExtend(TH1F::kAllAxes);

	//hEdgeCharge
	stringstream hEdgeChargeName; hEdgeChargeName<<"hEdgeCharge"<<FileNameEnd;
	hEdgeCharge = new TH1F(hEdgeChargeName.str().c_str(),hEdgeChargeName.str().c_str(),250,3100,4100);
	hEdgeCharge->SetTitle(hEdgeChargeName.str().c_str());
	hEdgeCharge->GetXaxis()->SetTitle("X Diamond (um)");
	hEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hEdgeChargeEvents = (TH1F*)hEdgeCharge->Clone("hEdgeChargeEvents");
	hEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

	hEdgeMeanCharge = (TH1F*)hEdgeCharge->Clone("hEdgeMeanCharge");
	hEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	//hEdgeCharge
	stringstream hyEdgeChargeName; hyEdgeChargeName<<"hyEdgeCharge"<<FileNameEnd;
	hyEdgeCharge = new TH1F(hyEdgeChargeName.str().c_str(),hyEdgeChargeName.str().c_str(),250,5000,6000);
	hyEdgeCharge->SetTitle(hyEdgeChargeName.str().c_str());
	hyEdgeCharge->GetXaxis()->SetTitle("Y Diamond (um)");
	hyEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hyEdgeChargeEvents = (TH1F*)hyEdgeCharge->Clone("hyEdgeChargeEvents");
	hyEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

	hyEdgeMeanCharge = (TH1F*)hyEdgeCharge->Clone("hyEdgeMeanCharge");
	hyEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	/////////////////
	// Alignment using dead cells profile
	/////////////////

	//hDeadCellsProfile
	//hEdgeCharge
	stringstream hDeadCellName; hDeadCellName<<"hDeadCell"<<FileNameEnd;
	hDeadCell = new TH1F(hDeadCellName.str().c_str(),hDeadCellName.str().c_str(),90,600,1050);
	hDeadCell->SetTitle(hDeadCellName.str().c_str());
	hDeadCell->GetXaxis()->SetTitle("Y Diamond (um)");
	hDeadCell->GetYaxis()->SetTitle("Mean Charge [ADC]");

	hDeadCellEvents = (TH1F*)hDeadCell->Clone("hDeadCellEvents");
	hDeadCellEvents->GetYaxis()->SetTitle("Entries");

	hDeadCellMeanCharge = (TH1F*)hDeadCell->Clone("hDeadCellMeanCharge");
	hDeadCellMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	/*
	stringstream hDeadCellsProfileChargeName; hDeadCellsProfileChargeName<<"hEdgeCharge%%"<<FileNameEnd;
	hDeadCellsProfileCharge = new TH1F(hDeadCellsProfileChargeName.str().c_str(),hDeadCellsProfileChargeName.str().c_str(),90,0,450);
	hDeadCellsProfileCharge->SetTitle(hDeadCellsProfileChargeName.str().c_str());
	hDeadCellsProfileCharge->GetXaxis()->SetTitle("X Diamond (um)");
	hDeadCellsProfileCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hDeadCellsProfileEvents = (TH1F*)hDeadCellsProfileCharge->Clone("hEdgeChargeEvents");
	hDeadCellsProfileEvents->GetYaxis()->SetTitle("Entries");

	hDeadCellsProfileMeanCharge = (TH1F*)hDeadCellsProfileCharge->Clone("hEdgeMeanCharge");
	hDeadCellsProfileMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	int DeadCellsArray[] = {13, 27, 53, 97};
	DeadCellsArrayPointer = DeadCellsArray;
	for(int i=0;i<4;i++){
		stringstream hDeadCellsChargeName; hDeadCellsChargeName<<"TotalCharge"<<(i)<<FileNameEnd;
		hDeadCellsCharge.push_back(new TH1F(hDeadCellsChargeName.str().c_str(),hDeadCellsChargeName.str().c_str(),90,0,450));
		stringstream hDeadCellsEventsName; hDeadCellsEventsName<<"Events"<<(i)<<FileNameEnd;
		hDeadCellsEvents.push_back(new TH1F(hDeadCellsEventsName.str().c_str(),hDeadCellsEventsName.str().c_str(),90,0,450));
	}
		*/
}
void TAnalysisOf3dDiamonds::saveYAlignmentHistos() {

	//For all Diamond hFidCutXvsFidCutYvsMeanCharge
	hCombinedMeanChargeYAlignment = new TCanvas();
	hCombinedMeanChargeYAlignment->cd();
	*hFidCutXvsFidCutYvsMeanChargeYAlignment = (*hFidCutXvsFidCutYvsChargeYAlignment/(*hFidCutXvsFidCutYvsEventsYAlignment));
	hFidCutXvsFidCutYvsMeanChargeYAlignment->SetEntries(hFidCutXvsFidCutYvsEventsYAlignment->Integral());
	hFidCutXvsFidCutYvsMeanChargeYAlignment->Draw("COLAH");
	/*DrawFidCutRegions(); //Draw Fiducial Cut Regions
		 */
	stringstream str1; str1<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanChargeYAlignmentNoFidDrawn"<<FileNameEnd<<".png";
	hCombinedMeanChargeYAlignment->SaveAs(str1.str().c_str());

	DrawYAlignmentFidCutRegions(); //Draw Fiducial Cut Regions
	stringstream str11; str11<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanCharge"<<FileNameEnd<<".png";
	hCombinedMeanChargeYAlignment->SaveAs(str11.str().c_str());

	//For h3DdetMeanCharge
	h3DdetMeanCharge = new TCanvas();
	h3DdetMeanCharge->cd();
	*hDetXvsDetY3DMeanCharge = (*hDetXvsDetY3D/(*hDetXvsDetY3DvsEvents));
	hDetXvsDetY3DMeanCharge->SetEntries(hDetXvsDetY3DvsEvents->Integral());
	hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetMeanCharge");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hGridReference->Draw("COL");
	DrawMetallisationGrid(h3DdetMeanCharge);
	stringstream str12; str12<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetMeanCharge"<<FileNameEnd<<".png";
	h3DdetMeanCharge->SaveAs(str12.str().c_str());

	//RebinnedMeanCharge
	h3DdetMeanChargeRebinned = new TCanvas();
	h3DdetMeanChargeRebinned->cd();
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	*hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DRebinned/(*hDetXvsDetY3DvsEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DvsEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetMeanChargeRebinned");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXTAH");
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	DrawMetallisationGrid(h3DdetMeanChargeRebinned);
	stringstream str121; str121<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetMeanChargeRebinned"<<FileNameEnd<<".png";
	h3DdetMeanChargeRebinned->SaveAs(str121.str().c_str());
	/*h3DdetEventsRebinned = new TCanvas();
	//h3DdetEventsRebinned->cd();
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");			//("COLZ");
	//DrawMetallisationGrid();
	stringstream str122; str122<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetEventsRebinned"<<"%%"<<FileNameEnd<<".png";
	h3DdetMeanChargeRebinned->SaveAs(str122.str().c_str());
		*/

	//h3DdetDeltaXChannel
	h3DdetDeltaXChannelCanvas = new TCanvas();
	h3DdetDeltaXChannelCanvas->cd();
	h3DdetDeltaXChannel->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetDeltaXChannelCanvas");		//Set title to require
	hGridReference->Draw("COL");
	h3DdetDeltaXChannel->Draw("sameCOLZAH");
	//h3DdetDeltaXChannelCanvas->Draw("sameTEXTAH");
	DrawMetallisationGrid(h3DdetDeltaXChannelCanvas);
	stringstream str12000; str12000<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetDeltaXChannel"<<FileNameEnd<<".png";
	h3DdetDeltaXChannelCanvas->SaveAs(str12000.str().c_str());

	//h3DdetDeltaXChannelAbove1000
	h3DdetDeltaXChannelAbove1000Canvas = new TCanvas();
	h3DdetDeltaXChannelAbove1000Canvas->cd();
	h3DdetDeltaXChannelAbove1000->SetStats(kFALSE);
	h3DdetDeltaXChannelAbove1000->Draw("COLZ");
	DrawMetallisationGrid(h3DdetDeltaXChannelAbove1000Canvas);
	stringstream str120000; str120000<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetDeltaXChannelAbove1000"<<FileNameEnd<<".png";
	h3DdetDeltaXChannelAbove1000Canvas->SaveAs(str120000.str().c_str());

	//RebinnedMeanChargeQuarterCell
	int QuaterCellEntrySum = 0;
	int QuaterCellEntrySumAddition = 0;
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			for(int k=0;k<4;k++){	//For each quater cell
				if(k<2){
					hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetBinContent((2*i+1),(2*j+1+k),hQuaterCellsLandau.at(i*11*4+j*4+k)->GetMean());
					hDetXvsDetY3DRebinnedQuarterCellRMS->SetBinContent((2*i+1),(2*j+1+k),hQuaterCellsLandau.at(i*11*4+j*4+k)->GetRMS());
				}
				if(k>1){
					hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetBinContent((2*i+2),(2*j+1+k-2),hQuaterCellsLandau.at(i*11*4+j*4+k)->GetMean());
					hDetXvsDetY3DRebinnedQuarterCellRMS->SetBinContent((2*i+2),(2*j+1+k-2),hQuaterCellsLandau.at(i*11*4+j*4+k)->GetRMS());
				}
				QuaterCellEntrySumAddition = QuaterCellEntrySum;
				QuaterCellEntrySum = QuaterCellEntrySumAddition + hQuaterCellsLandau.at(i*11*4+j*4+k)->GetEntries();
				//histSaver->SaveHistogram(hQuaterCellsLandau.at(i*11*4+j*4+k));
			}
		}
	}
	hDetXvsDetY3DMeanChargeRebinnedQuarterCellCanvas = new TCanvas();
	hDetXvsDetY3DMeanChargeRebinnedQuarterCellCanvas->cd();
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetEntries(QuaterCellEntrySum);
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetQuarterCellMeanCharge");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameCOLZAH");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->Draw("sameTEXTAH");
	DrawMetallisationGrid(hDetXvsDetY3DMeanChargeRebinnedQuarterCellCanvas);
	stringstream str128; str128<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetQuarterCellMeanCharge"<<FileNameEnd<<".png";
	hDetXvsDetY3DMeanChargeRebinnedQuarterCellCanvas->SaveAs(str128.str().c_str());

	hDetXvsDetY3DRebinnedQuarterCellRMSCanvas = new TCanvas();
	hDetXvsDetY3DRebinnedQuarterCellRMSCanvas->cd();
	hDetXvsDetY3DRebinnedQuarterCellRMS->SetEntries(QuaterCellEntrySum);
	hGridReference->SetTitle("h3DdetQuarterCellRMS");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DRebinnedQuarterCellRMS->Draw("sameCOLZAH");
	hDetXvsDetY3DRebinnedQuarterCellRMS->Draw("sameTEXTAH");
	DrawMetallisationGrid(hDetXvsDetY3DRebinnedQuarterCellRMSCanvas);
	stringstream str129; str129<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetQuarterCellRMS"<<FileNameEnd<<".png";
	hDetXvsDetY3DRebinnedQuarterCellRMSCanvas->SaveAs(str129.str().c_str());

	//Cell Overlay
	int CellsHarris18Good[] = {39,40,41,   50,51,52,   61,62,63,   72,73,74,   83,84,85,   94,95,96};
	int CellsHarris10Bad[] = {1,6,7,8,13,27,45,53,77,90,97};

	//cout<<"CellArraySize is: "<<hCellsCharge.size()<<endl;
	//for(int i=0;i<hCellsCharge.size();i++){
	for(int i=0;i<9;i++){		//run over each cell
		for(int j=0;j<11;j++){
			//histSaver->SaveHistogram(hCellsLandau.at(i*11+j));	//Save each cells Landau;
			LandauGaussFit landauGauss;
			//TF1* fit = landauGauss.doLandauGaussFit(hCellsDeltaX.at(i*11+j));
			//histSaver->SaveHistogram(hCellsDeltaX.at(i*11+j));		//Save each cells X res in Ch.
			hCellNumbering->SetBinContent((i+1),(j+1),(i*11+j));
			hDetXvsDetY3DRebinnedRMS->SetBinContent((i+1),(j+1),(hCellsLandau.at(i*11+j)->GetRMS()));
			hBinnedMeanCharge->Fill(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1));	//Fill 1D hist with each cells mean charge

			int Continue = 1;
			for(int k=0;k<11;k++){
				if((i*11+j)==CellsHarris10Bad[k]){
					Continue = 0;
					int Entries = hCellsHarris10Bad->GetEntries();;
					hCellsHarris10Bad->Add(hCellsLandau.at(i*11+j));
					hCellsHarris10Bad->SetEntries((Entries+hCellsLandau.at(i*11+j)->GetEntries()));

					RebinnedQuarterCellFails->SetBinContent(i+1,j+1,4);

					for(int k=0;k<4;k++){		//To fill highlighted grading, plot 5.
						//NumberQuarterFails++;
						if(k<2)			//To make it more obvious which quarters have failed.
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+1),(2*j+1+k),500);
						if(k>1)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+2),(2*j+1+k-2),500);
					}
					//For Transparent Analysis
					hQuarterCellGradedTransparentLandau.at(4)->Add(hCellTransparentLandau.at(i*11+j));

					if(CellsHarris10Bad[k]>10&&CellsHarris10Bad[k]<88){		//Only non edge cells used
						//hCellsDeltaXQuarterCellGrading
						hCellsDeltaXQuarterCellGrading.at(4)->Add(hCellsDeltaX.at(i*11+j));
					}
				}
			}
			if(Continue == 1){		//If a bad cell, the code won't continue.

				//hCellsHarrisGoodCells
				for(int k=0;k<18;k++){
					if((i*11+j)==CellsHarris18Good[k]){
						int Entries = hCellsHarris18Good->GetEntries();
						hCellsHarris18Good->Add(hCellsLandau.at(i*11+j));
						hCellsHarris18Good->SetEntries((Entries+hCellsLandau.at(i*11+j)->GetEntries()));
					}
				}

				int NumberQuarterFails = 0;
				float QuarterMeanCharge[4];
				for(int k=0;k<4;k++){
					SortArrayPointer[k] = k;
					QuarterMeanCharge[k] = hQuaterCellsLandau.at(i*11*4+j*4+k)->GetMean();
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
				RebinnedQuarterCellFails->SetBinContent(i+1,j+1,NumberQuarterFails);
				if(NumberQuarterFails ==0)
					hCellsTransparentHitPositionCellGraded0->Add(hCellsTransparentHitPosition.at(i*11+j));
				for(int k=0;k<NumberQuarterFails;k++){		//To fill highlighted grading, plot 5.
						//NumberQuarterFails++;
					if(SortArrayPointer[k]<2)			//To make it more obvious which quarters have failed.
						hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+1),(2*j+1+SortArrayPointer[k]),500);
					if(SortArrayPointer[k]>1)
						hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+2),(2*j+1+SortArrayPointer[k]-2),500);
				}
					/*else{
						if(k<2)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+1),(2*j+1+k),1000);
						if(k>1)
							hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->SetBinContent((2*i+2),(2*j+1+k-2),1000);
					}
				}
					*/
				for(int k=0;k<4;k++){	//For each Grading
					if(k==NumberQuarterFails){
						//For Transparent Analysis
						hQuarterCellGradedTransparentLandau.at(k)->Add(hCellTransparentLandau.at(i*11+j));
						hCellsDeltaXQuarterCellGrading.at(k)->Add(hCellsDeltaX.at(i*11+j));
						for(int l=0;l<4;l++){	//To fill each quarter cell of the graded cell.
							hDetXvsDetY3DMeanChargeQuarterCellGradingLandau.at(k)->Add(hQuaterCellsLandau.at(i*11*4+j*4+l));
							if(l<2){
								hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->SetBinContent((2*i+1),(2*j+1+l),hQuaterCellsLandau.at(i*11*4+j*4+l)->GetMean());
							}
							if(l>1){
								hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->SetBinContent((2*i+2),(2*j+1+l-2),hQuaterCellsLandau.at(i*11*4+j*4+l)->GetMean());
							}
						}
					}
				}

				//if(i!=BadCells[k]){
				if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)>1000)&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)<1100)){
					cout<<"Cell with average charge > 1000: "<<(i*11+j)<<endl;
					//if(NumberQuarterFails<3){
					hCellsOverlayedCharge->Add(hCellsCharge.at(i*11+j));
					hCellsOverlayedEvents->Add(hCellsEvents.at(i*11+j));
					hCellsOverlayedEventsNoColumns->Add(hCellsEventsNoColumn.at(i*11+j));
					hCellsOverlayedLandauNoColumn->Add(hCellsLandauNoColumn.at(i*11+j));
					//cout<<"Cell: "<<(i*11+j)<<" is overlayed."<<endl;

					hDetXvsDetY3DOverview->SetBinContent(i+1,j+1,1);
					//cout<<"Bin content is: "<<hDetXvsDetY3DOverview->GetBinContent(1,1)<<endl;
					//histSaver->SaveHistogram(hCellsCharge.at(i));
					//gStyle->SetOptStat("irm");
					hCellsEvents.at(i*11+j)->SetEntries(hCellsEvents.at(i*11+j)->Integral()); //Set the number of entries to the integral of hEvents
					//histSaver->SaveHistogram(hCellsEvents.at(i*11+j));
					//cout<<"The integral of Events"<<(i*11+j)<<" is: "<<hCellsEvents.at(i*11+j)->Integral()<<endl;
				}
				for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
					if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)>(k*100))&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)<((k+1)*100))){
						hCellsLandauGraded.at(k)->Add(hCellsLandau.at(i*11+j));
						hCellsLandauGradedNoColumn.at(k)->Add(hCellsLandauNoColumn.at(i*11+j));
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
	hCellsOverlayedCanvas = new TCanvas();
	hCellsOverlayedCanvas->cd();
	*hCellsOverlayedMeanCharge = (*hCellsOverlayedCharge/(*hCellsOverlayedEvents));
	hCellsOverlayedMeanCharge->SetEntries(hCellsOverlayedEvents->Integral());

	hCellsOverlayedMeanCharge->GetZaxis()->SetRangeUser(800,1200);
	hCellsOverlayedMeanCharge->SetContour(99);
	hCellsOverlayedMeanCharge->Draw("COLZ");
	//hCellsOverlayed55RMS->Draw("sameTEXT");
	stringstream str13; str13<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedMeanCharge"<<FileNameEnd<<".png";
	hCellsOverlayedCanvas->SaveAs(str13.str().c_str());

	/*hCellsOverlayedEvents->Draw("sameTEXT");
	stringstream str133; str133<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedMeanChargeWithEntries"<<FileNameEnd<<".png";
	hCellsOverlayedCanvas->SaveAs(str133.str().c_str());
	histSaver->SaveHistogram(hBinnedMeanCharge);
		*/

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
			if((i*3+ii)==0){
				for(int j=0;j<15;j++){
					for(int k=0;k<15;k++){
						hCellsOverlayed1010RMS->SetBinContent(j+1,k+1,hCellsOverlayBinSpec1010.at(j*15+k)->GetRMS());
						//histSaver->SaveHistogram(hCellsOverlayBinSpec1010.at(j*15+k));
					}
				}
			}	//end of if first plot
			hCellsOverlayedBinAlignmentCanvas1.push_back(new TCanvas());
			hCellsOverlayedBinAlignmentCanvas1.at(i*3+ii)->cd();
			*hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii) = (*hCellsOverlayedChargeBinAlignment1.at(i*3+ii)/(*hCellsOverlayedEventsBinAlignment1.at(i*3+ii)));
			hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii)->SetEntries(hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->Integral());
			hCellsOverlayedMeanChargeBinAlignment1.at(i*3+ii)->Draw("COLZ");
			if((i*3+ii)==0)
				hCellsOverlayed1010RMS->Draw("sameTEXT");
			//hCellsOverlayedEventsBinAlignment1.at(i*3+ii)->Draw("sameTEXT");
			stringstream str13001; str13001<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedMeanCharge"<<"BinAlignment"<<i*3+ii<<FileNameEnd<<".png";
			hCellsOverlayedBinAlignmentCanvas1.at(i*3+ii)->SaveAs(str13001.str().c_str());
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
	TF1* fit0 = landauGauss.doLandauGaussFit(hCellsHarris18Good);
	histSaver->SaveHistogram(hCellsHarris18Good);
	TF1* fit1 = landauGauss.doLandauGaussFit(hCellsHarris10Bad);
	histSaver->SaveHistogram(hCellsHarris10Bad);
	TF1* fit4 = landauGauss.doLandauGaussFit(hCellsOverlayedColumnLandau);
	histSaver->SaveHistogram(hCellsOverlayedColumnLandau);

	//TF1* fit5 = landauGauss.doLandauGaussFit(hTransparentCharge3D);
	histSaver->SaveHistogram(hTransparentCharge3D);

	//hDetXvsDetY3DMeanChargeQuarterCellGrading
	//vector<TF1*> GradedFits;
	for(int k=0;k<6;k++){
		hDetXvsDetY3DMeanChargeQuarterCellGradingCanvas.push_back(new TCanvas());
		hDetXvsDetY3DMeanChargeQuarterCellGradingCanvas.at(k)->cd();
		hGridReference->SetTitle("hDetXvsDetY3DMeanChargeQuarterCellGrading");		//Set title to require
		hGridReference->Draw("COL");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->Draw("sameCOLZAH");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->Draw("sameTEXTAH");
		//gStyle->SetPaintTextFormat(3.2g);
		DrawMetallisationGrid(hDetXvsDetY3DMeanChargeQuarterCellGradingCanvas.at(k));
		stringstream str1401; str1401<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hDetXvsDetY3DMeanChargeQuarterCellGrading"<<k<<FileNameEnd<<".png";
		hDetXvsDetY3DMeanChargeQuarterCellGradingCanvas.at(k)->SaveAs(str1401.str().c_str());
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
	//hCellsTransparentHitPositionCellGraded0
	hCellsTransparentHitPositionCellGraded0Canvas = new TCanvas();
	hCellsTransparentHitPositionCellGraded0Canvas->cd();
	hCellsTransparentHitPositionCellGraded0->SetStats(kFALSE);
	hCellsTransparentHitPositionCellGraded0->Draw("COLZ");
	hCellsTransparentHitPositionCellGraded0->Draw("sameTEXT");
	stringstream str14020; str14020<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsTransparentHitPositionCellGraded0"<<FileNameEnd<<".png";
	hCellsTransparentHitPositionCellGraded0Canvas->SaveAs(str14020.str().c_str());

	//hDetXvsDetY3DMeanChargeHighlightedQuarters
	hDetXvsDetY3DMeanChargeHighlightedQuartersCanvas = new TCanvas();
	hDetXvsDetY3DMeanChargeHighlightedQuartersCanvas->cd();
	hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReference->SetTitle("hDetXvsDetY3DMeanChargeHighlightedQuarters");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hDetXvsDetY3DMeanCharge->Draw("sameTEXTAH");
	//hDetXvsDetY3DMeanCharge->Draw("sameCOLZ");
	DrawMetallisationGrid(hDetXvsDetY3DMeanChargeHighlightedQuartersCanvas);
	stringstream str1402; str1402<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hDetXvsDetY3DMeanChargeHighlightedQuarters"<<FileNameEnd<<".png";
	hDetXvsDetY3DMeanChargeHighlightedQuartersCanvas->SaveAs(str1402.str().c_str());

	//hCellNumbering
	hCellNumberingCanvas = new TCanvas();
	hCellNumberingCanvas->cd();
	hGridReference->SetTitle("hCellNumbering");		//Set title to require
	hGridReference->Draw("COL");
	hCellNumbering->Draw("sameCOLZAH");
	hCellNumbering->Draw("sameTEXTAH");
	DrawMetallisationGrid(hCellNumberingCanvas);
	stringstream str1400; str1400<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellNumbering"<<FileNameEnd<<".png";
	hCellNumberingCanvas->SaveAs(str1400.str().c_str());

	//hCellsOverlayedEventsNoColumns
	hCellsEventsNoColumnCanvas = new TCanvas();
	hCellsEventsNoColumnCanvas->cd();
	hCellsOverlayedEventsNoColumns->Draw("COLZ");
	hCellsOverlayedEventsNoColumns->Draw("sameTEXT");
	stringstream str14; str14<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedEventsNoColumns"<<FileNameEnd<<".png";
	hCellsEventsNoColumnCanvas->SaveAs(str14.str().c_str());

	//hDetXvsDetY3DRebinnedRMS
	hDetXvsDetY3DRebinnedRMSCanvas = new TCanvas();
	hDetXvsDetY3DRebinnedRMSCanvas->cd();
	hGridReference->SetTitle("h3DdetRebinnedRMS");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DRebinnedRMS->Draw("sameCOLZAH");
	hDetXvsDetY3DRebinnedRMS->Draw("sameTEXTAH");
	DrawMetallisationGrid(hDetXvsDetY3DRebinnedRMSCanvas);
	stringstream str145; str145<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetRebinnedRMS"<<FileNameEnd<<".png";
	hDetXvsDetY3DRebinnedRMSCanvas->SaveAs(str145.str().c_str());

	//hCellsLandauGraded			//Save each of the histograms
	for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
		//histSaver->SaveHistogram(hCellsLandauGraded.at(k));
		//histSaver->SaveHistogram(hCellsLandauGradedNoColumn.at(k));
	}

	//hCellsLandau2D
	TH2D* hDetXvsDetY3DMeanChargeRebinnedForCellOrdering = (TH2D*)hDetXvsDetY3DMeanChargeRebinned->Clone("hDetXvsDetY3DMeanChargeRebinnedForCellOrdering");
	int hCellsLandau2DEntries =0;
	int hCellsLandau2DEntriesPlus =0;
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			hCellsLandau2DEntriesPlus = hCellsLandau2DEntries;
			hCellsLandau2DEntries = hCellsLandau2DEntriesPlus + hCellsLandau.at(i*11+j)->GetEntries();

			//hCellsMeanClusteSize
			hCellsMeanClusteSize->SetBinContent(i+1,j+1,hCellsClusteSize.at(i*11+j)->GetMean());
			//Fit Landau to each cells landau.
			//hCellsLandau.at(i*11+j)->Fit("landau");
			/*Int_t minX, minY, minZ;
			Int_t maxX, maxY, maxZ;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMaximumBin(maxX, maxY, maxZ);
			cout<<"The Maximum bin is: "<<maxX<<", "<<maxY<<endl;
			cout<<"The Maximum ADC is: "<<hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMaximum()<<endl;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimumBin(minX, minY, minZ);
			cout<<"The Minimum bin is: "<<minX<<", "<<minY<<endl;;
			cout<<"The Minimum ADC is: "<<hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimum()<<endl;
				*/
			Int_t minX, minY, minZ;
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->GetMinimumBin(minX, minY, minZ);	//Returns cell with Minimum average.
			//cout<<"hCellsLandau2D bin in y: "<<(i*11+j+1)<<".   Filled with Cell: "<<((minX-1)*11+(minY-1))<<endl;
			stringstream binLabel; binLabel<<((minX-1)*11+(minY-1));
			hCellsLandau2D->GetYaxis()->SetBinLabel(i*11+j+1,binLabel.str().c_str());
			hCellsLandau2D->GetYaxis()->SetLabelSize(0.02);
			hCellsLandau2DQuarterFail->GetYaxis()->SetBinLabel(i*11+j+1,binLabel.str().c_str());
			hCellsLandau2DQuarterFail->GetYaxis()->SetLabelSize(0.02);
			hCellsLandau2DQuarterFail->SetBinContent(1,(i*11+j+1),RebinnedQuarterCellFails->GetBinContent(minX,minY));
			for(int k=0;k<256;k++){
				hCellsLandau2D->SetBinContent((k+1),(i*11+j+1),hCellsLandau.at((minX-1)*11+(minY-1))->GetBinContent((k+1)));		//Set Cell in Y, set all pulseheights in x.
			}	//End of cell Landau loop
			hDetXvsDetY3DMeanChargeRebinnedForCellOrdering->SetBinContent(minX, minY, 99999);	//The minimum bin now becomes maximum.
		}
	}	//End of looping over all cells
	hCellsLandau2DCanvas = new TCanvas();
	hCellsLandau2DCanvas->SetCanvasSize(1500,3000);
	hCellsLandau2DCanvas->cd();
	//hCellsLandau2D->LabelsDeflate();
	hCellsLandau2D->SetEntries(hCellsLandau2DEntries);
	//hCellsLandau2DQuarterFail->Draw("COLZ");
	hCellsLandau2D->SetStats(kFALSE);
	hCellsLandau2D->Draw("COLZ");
	//hCellsLandau2D->Draw("sameTEXT");
	stringstream str144; str144<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsLandau2D"<<FileNameEnd<<".png";
	hCellsLandau2DCanvas->SaveAs(str144.str().c_str());

	hCellsLandau2DCanvas1 = new TCanvas();
	hCellsLandau2DCanvas1->SetCanvasSize(1500,3000);
	hCellsLandau2DCanvas1->cd();
	hCellsLandau2DQuarterFail->SetStats(kFALSE);
	hCellsLandau2DQuarterFail->Draw("COL");
	hCellsLandau2D->Draw("sameCOLZ");
	stringstream str1440; str1440<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsLandau2DHighlightedQuarters"<<FileNameEnd<<".png";
	hCellsLandau2DCanvas1->SaveAs(str1440.str().c_str());

	//For h3DdetMeanCharge with Mean ClusterSize.
	h3DdetMeanChargeWithMeanClusterSize = new TCanvas();
	h3DdetMeanChargeWithMeanClusterSize->cd();
	*hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DRebinned/(*hDetXvsDetY3DvsEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DvsEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetMeanChargeWithMeanClusterSize");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	hCellsMeanClusteSize->Draw("sameTEXTAH");
	DrawMetallisationGrid(h3DdetMeanChargeWithMeanClusterSize);
	stringstream str1200; str1200<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetMeanChargeWithMeanClusterSize"<<FileNameEnd<<".png";
	h3DdetMeanChargeWithMeanClusterSize->SaveAs(str1200.str().c_str());

	//RebinnedQuarterCellFails
	RebinnedQuarterCellFailsCanvas = new TCanvas();
	RebinnedQuarterCellFailsCanvas->cd();
	RebinnedQuarterCellFails->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetNumberofQuarterCellFails");		//Set title to require
	hGridReference->Draw("COL");
	RebinnedQuarterCellFails->Draw("sameCOLZAH");
	RebinnedQuarterCellFails->Draw("sameTEXTAH");
	DrawMetallisationGrid(RebinnedQuarterCellFailsCanvas);
	stringstream str146; str146<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetNumberofQuarterCellFails"<<FileNameEnd<<".png";
	RebinnedQuarterCellFailsCanvas->SaveAs(str146.str().c_str());

	//hDetXvsDetY3DOverview
	hOverview = new TCanvas();
	hOverview->cd();
	hDetXvsDetY3DOverview->SetStats(kFALSE);
	hGridReference->SetTitle("h3DdetOverview");		//Set title to require
	hGridReference->Draw("COL");
	hDetXvsDetY3DOverview->Draw("sameCOLZAH");
	stringstream str131; str131<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetOverview"<<FileNameEnd<<".png";
	hOverview->SaveAs(str131.str().c_str());

	//hDeadCell
	//histSaver->SaveHistogram(hDeadCell);
	//histSaver->SaveHistogram(hDeadCellEvents);

	hDeadCellMeanChargeCanvas = new TCanvas();
	hDeadCellMeanChargeCanvas->cd();
	*hDeadCellMeanCharge = (*hDeadCell/(*hDeadCellEvents));
	hDeadCellMeanCharge->SetEntries(hDeadCellEvents->Integral());
	hDeadCellMeanCharge->Draw();
	Float_t ymax1 = hDeadCellMeanCharge->GetMaximum();
	TLine* CellEdge1 = new TLine(750,0,750,ymax1);
	TLine* CellEdge2 = new TLine(900,0,900,ymax1);
	CellEdge1->SetLineWidth(2);		CellEdge2->SetLineWidth(2);
	CellEdge1->SetLineColor(kRed);	CellEdge2->SetLineColor(kRed);
	CellEdge1->Draw("same");		CellEdge2->Draw("same");

	stringstream hDeadCellMeanChargeCanvasName; hDeadCellMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hDeadCellMeanCharge"<<FileNameEnd<<".png";
	hDeadCellMeanChargeCanvas->SaveAs(hDeadCellMeanChargeCanvasName.str().c_str());

	//hxEdgeCharge
	//histSaver->SaveHistogram(hEdgeCharge);
	//histSaver->SaveHistogram(hEdgeChargeEvents);

	hEdgeMeanChargeCanvas = new TCanvas();
	hEdgeMeanChargeCanvas->cd();
	*hEdgeMeanCharge = (*hEdgeCharge/(*hEdgeChargeEvents));
	hEdgeMeanCharge->SetEntries(hEdgeChargeEvents->Integral());
	hEdgeMeanCharge->Draw();
	Float_t ymax = hEdgeMeanCharge->GetMaximum();
	TLine* DetectorEdge = new TLine(3715,0,3715,ymax);
	DetectorEdge->SetLineWidth(2);
	DetectorEdge->SetLineColor(kRed);
	DetectorEdge->Draw("same");
	stringstream hEdgeMeanChargeCanvasName; hEdgeMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hEdgeMeanCharge"<<FileNameEnd<<".png";
	hEdgeMeanChargeCanvas->SaveAs(hEdgeMeanChargeCanvasName.str().c_str());
	//histSaver->SaveHistogram(hEdgeMeanCharge);

	//hyEdgeCharge
	//histSaver->SaveHistogram(hyEdgeCharge);
	//histSaver->SaveHistogram(hyEdgeChargeEvents);

	hyEdgeMeanChargeCanvas = new TCanvas();
	hyEdgeMeanChargeCanvas->cd();
	*hyEdgeMeanCharge = (*hyEdgeCharge/(*hyEdgeChargeEvents));
	hyEdgeMeanCharge->SetEntries(hyEdgeChargeEvents->Integral());
	hyEdgeMeanCharge->Draw();
	//Float_t ymax = hyEdgeMeanCharge->GetMaximum();
	//TLine* DetectorEdge = new TLine(3715,0,3715,ymax);
	//DetectorEdge->SetLineWidth(2);
	//DetectorEdge->SetLineColor(kRed);
	//DetectorEdge->Draw("same");
	stringstream hyEdgeMeanChargeCanvasName; hyEdgeMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hyEdgeMeanCharge"<<FileNameEnd<<".png";
	hyEdgeMeanChargeCanvas->SaveAs(hyEdgeMeanChargeCanvasName.str().c_str());

}
void TAnalysisOf3dDiamonds::YAlignment() {

	if(!eventReader->isValidTrack())
		return;
	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	TPositionPrediction *predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	if(predictedPosition->getChi2X()>100||predictedPosition->getChi2Y()>100)
		return;

	Float_t xPos = predictedPosition->getPositionX();	//Predicted positions in labframe
	Float_t yPos = predictedPosition->getPositionY();
	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	Float_t Xdet = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos);
	Float_t Ydet = GetYPositionInDetSystem()-3890;	//3890 is the calculated offset from comparing x and y edge charge profiles.

	if(!eventReader->isInFiducialCut())	//This is a larger fiducial cut around silicon
		return;

	//For Transparent Analysis.
	float hEntries0 = 0;
	float TransparentCharge = 0;
	float TransparentChargeAddition = 0;
	int SaturatedEvent = 0;
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			hEntries0 = hCellsEventsCheck.at(i*11+j)->Integral();
			float xminus = 2365+i*150+5;		//2365 is the start of the 3D detector in x
			float yminus = j*150;
			hCellsEventsCheck.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
			int Xcell = (Xdet-xminus);
			Int_t HitCell;
			//Int_t* HitChannels;
			if(hEntries0 != hCellsEventsCheck.at(i*11+j)->Integral()){	//If entries in this histogram range have increased fill corresponding hCellLandau
				HitCell = (i*11+j);
				cout<<"Hit cell: "<<HitCell<<endl;
				//HitChannels = CellToChannel(HitCell);
				cout<<CellToChannel(HitCell,Xcell)[0]<<endl;
				for(int k=1;k<CellToChannel(HitCell,Xcell)[0]+1;k++){	//First memeber of array is the number of channels to readout.
					if(eventReader->isSaturated(subjectDetector,CellToChannel(HitCell,Xcell)[k]))
						SaturatedEvent = 1;
					cout<<"Hit channel: "<<CellToChannel(HitCell,Xcell)[k]<<endl;
					//TransparentCharge = eventReader->getAdcValue(subjectDetector,HitChannel);
					TransparentChargeAddition = TransparentCharge;
					TransparentCharge = eventReader->getSignal(subjectDetector,CellToChannel(HitCell,Xcell)[k]) + TransparentChargeAddition;
				}	//End of for HitChannels
				if(SaturatedEvent == 0){	// Only plot events with no saturated channels.
					cout<<"Transparent charge is: "<<TransparentCharge<<endl;
					hTransparentCharge3D->Fill(TransparentCharge);
					hCellTransparentLandau.at(i*11+j)->Fill(TransparentCharge);
					if(TransparentCharge<700)
						hCellsTransparentHitPosition.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
					//cout<<eventReader->getSignal(subjectDetector,HitChannel)<<endl;
				}
			}
		}
	}

	//eventReader->getAdcValue(subjectDetector,CellToChannel());

	//Int_t getAdcValue(UInt_t det,UInt_t ch);

	//CellToChannel();

	if(eventReader->getNDiamondClusters()==0||eventReader->getNDiamondClusters()==2||eventReader->getNDiamondClusters()==3)	//Only single cluster events
		return;


	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster())
		return;

	//HitandSeedCount(&diamondCluster, 0);

	//if(HitCount==1&&SeedCount==1){
	//if(SeedCount<3){
	if(RemoveClustersWithHitOutside3D(&diamondCluster) == 0){	//Only clusters with all hit channels within 3D detector continue

		//ClusterShape(&diamondCluster);

		/*if(YAlignmentFiducialCut())
			return;
				*/

		hFidCutXvsFidCutYvsChargeYAlignment->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
		hFidCutXvsFidCutYvsEventsYAlignment->Fill(fiducialValueX,fiducialValueY,1);

		for(int i=0;i<2;i++){		//For both x and y edge charge profiles
			if(i==0){
				if(!YAlignmentFiducialCut(0)){
					Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
					hEdgeChargeEvents->Fill(positionInDetSystemMetric);
					hEdgeCharge->Fill(positionInDetSystemMetric, diamondCluster.getCharge(false));
				}
			}
			if(i==1){
				if(!YAlignmentFiducialCut(1)){
					hyEdgeChargeEvents->Fill(GetYPositionInDetSystem());
					hyEdgeCharge->Fill(GetYPositionInDetSystem(), diamondCluster.getCharge(false));
				}
			}
		}	//End of for

		hDetXvsDetY3D->Fill(Xdet,Ydet,diamondCluster.getCharge(false));
		hDetXvsDetY3DvsEvents->Fill(Xdet,Ydet,1);
		hDetXvsDetY3DRebinned->Fill(Xdet,Ydet,diamondCluster.getCharge(false));
		hDetXvsDetY3DvsEventsRebinned->Fill(Xdet,Ydet,1);

		//h3DdetDeltaXChannel
		//if(RemoveEdge3DClusters(&diamondCluster)==0){
			Float_t XdetChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,Xdet);
			//h3DdetDeltaXChannel->Fill(Xdet,Ydet,sqrt((diamondCluster.getHighestSignalChannel()-XdetChannelSpace)*(diamondCluster.getHighestSignalChannel()-XdetChannelSpace)));
			if(sqrt((diamondCluster.getHighestSignalChannel()-XdetChannelSpace)*(diamondCluster.getHighestSignalChannel()-XdetChannelSpace))<1000){
				h3DdetDeltaXChannel->Fill(Xdet,Ydet,sqrt((diamondCluster.getHighestSignalChannel()-XdetChannelSpace)*(diamondCluster.getHighestSignalChannel()-XdetChannelSpace)));
				//cout<<"Delta X is: "<<sqrt((diamondCluster.getHighestSignalChannel()-XdetChannelSpace)*(diamondCluster.getHighestSignalChannel()-XdetChannelSpace))<<endl;
				//cout<<"Xdet is: "<<Xdet<<" Ydet is: "<<Ydet<<endl;
				//cout<<"HighestChannel: "<<diamondCluster.getHighestSignalChannel()<<" PredictedChannel: "<<XdetChannelSpace<<endl;
				//ClusterChannels(&diamondCluster);
			}
		//}
		//if(sqrt((diamondCluster.getHighestSignalChannel()-XdetChannelSpace)*(diamondCluster.getHighestSignalChannel()-XdetChannelSpace))>1000)
			//h3DdetDeltaXChannelAbove1000->Fill(Xdet, Ydet, 1);

		if(Xdet>2665&&Xdet<2815){	//Dead cell chosen to be: (2,5) in cell space
			hDeadCell->Fill(Ydet, diamondCluster.getCharge(false));
			hDeadCellEvents->Fill(Ydet, 1);
		}
		//	To fill multiple cell histograms
		//for(int i=0;i<hCellsCharge.size();i++){
		int CellsAboveThousand[]={0,12,14,15,16,17,18,20,21,22,23,25,26,28,31,32,33,35,36,38,39,40,41,42,43,48,49,50,51,52,54,55,56,58,61,62,66,69,70,72,73,76,81,83,84,85,86,87,92,93,94,95,96,98};	//Cells predetermined to be above 1000 ADC.
		float hEntries = 0;
		for(int i=0;i<9;i++){
			for(int j=0;j<11;j++){
				hEntries = hCellsEvents.at(i*11+j)->Integral();
				float xminus = 2365+i*150+5;		//2365 is the start of the 3D detector in x
				float yminus = j*150;
				hCellsCharge.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),diamondCluster.getCharge(false));
				hCellsEvents.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
				if(hEntries != hCellsEvents.at(i*11+j)->Integral()){	//If entries in this histogram range have increased fill corresponding hCellLandau
					hCellsLandau.at(i*11+j)->Fill(diamondCluster.getCharge(false));
					hCellsClusteSize.at(i*11+j)->Fill((diamondCluster.getClusterSize()-2));		//Cluster seems to be 2 smaller than ClusterSize for some reason?
					hCellsDeltaX.at(i*11+j)->Fill((diamondCluster.getHighestSignalChannel()-XdetChannelSpace));
					hCellsColumnCheck55->Fill((Xdet-xminus),(Ydet-yminus),1);	//Fill Column Check histo 5*5
					hCellsColumnCheck1010->Fill((Xdet-xminus),(Ydet-yminus),1);	//Fill Column Check histo 10*10
					//To fill quarter cell histograms
					if((Xdet-xminus)>0&&(Xdet-xminus)<75  &&  (Ydet-yminus)>0&&(Ydet-yminus)<75){	//bottom left
						hQuaterCellsLandau.at(i*11*4+j*4)->Fill(diamondCluster.getCharge(false));
					}
					if((Xdet-xminus)>0&&(Xdet-xminus)<75  &&  (Ydet-yminus)>75&&(Ydet-yminus)<150){	//top left
						hQuaterCellsLandau.at(i*11*4+j*4+1)->Fill(diamondCluster.getCharge(false));
					}
					if((Xdet-xminus)>75&&(Xdet-xminus)<150  &&  (Ydet-yminus)>0&&(Ydet-yminus)<75){	//bottom right
						hQuaterCellsLandau.at(i*11*4+j*4+2)->Fill(diamondCluster.getCharge(false));
					}
					if((Xdet-xminus)>75&&(Xdet-xminus)<150  &&  (Ydet-yminus)>75&&(Ydet-yminus)<150){	//top right
						hQuaterCellsLandau.at(i*11*4+j*4+3)->Fill(diamondCluster.getCharge(false));
					}
					//
					for(int k=0;k<30;k++){		//looping over number of cell bins in x		//5*5um
						for(int l=0;l<30;l++){		//looping over number of cell bins in y
							int ColumnEntriesX[] = {1,1,2,2,1,1,2,2,16,16,17,17,29,29,30,30,29,29,30,30};	//Column cell x coordinates
							int ColumnEntriesY[] = {1,2,1,2,29,30,29,30,15,16,15,16,1,2,1,2,29,30,29,30};	//Column cell y coordinates
							if(hCellsColumnCheck55->GetBinContent((k+1),(l+1))==1){	//True when in hit bin
								/*int inColumn =0;
								for(int m=0;m<20;m++){	//Loop over number of column cells
									if((k+1)==ColumnEntriesX[m]&&(l+1)==ColumnEntriesY[m]) inColumn =1;
								}
								if(inColumn==0){
									hCellsLandauNoColumn.at(i*11+j)->Fill(diamondCluster.getCharge(false));
									hCellsEventsNoColumn.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
								}
									*/
								//Fill 5*5 overlay bin spec.
								for(int m=0;m<54;m++){
									if(CellsAboveThousand[m]==(i*11+j))
										hCellsOverlayBinSpec55.at(k*30+l)->Fill(diamondCluster.getCharge(false));
								}
							}	//End of in hit bin
							hCellsColumnCheck55->SetBinContent((k+1),(l+1),0);	//Empty column check histo
						}	//End of y cells loop
					}	//End of x cells loop

					for(int k=0;k<15;k++){		//looping over number of cell bins in x		//10*10um
						for(int l=0;l<15;l++){		//looping over number of cell bins in y
							int ColumnEntriesX[] = {1,8,9};	//Column cell x coordinates
							int ColumnEntriesY[] = {1,8,8};	//Column cell y coordinates
							if(hCellsColumnCheck1010->GetBinContent((k+1),(l+1))==1){	//True when in hit bin
								int inColumn =0;
								for(int m=0;m<3;m++){	//Loop over number of column cells
									if((k+1)==ColumnEntriesX[m]&&(l+1)==ColumnEntriesY[m]) inColumn =1;
								}
								if(inColumn==1)
									hCellsOverlayedColumnLandau->Fill(diamondCluster.getCharge(false));	//Fill with column Spectrum.
								if(inColumn==0){
									hCellsLandauNoColumn.at(i*11+j)->Fill(diamondCluster.getCharge(false));
									hCellsEventsNoColumn.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
								}
								//Fill 10*10 overlay bin spec.
								for(int m=0;m<54;m++){
									if(CellsAboveThousand[m]==(i*11+j))
										hCellsOverlayBinSpec1010.at(k*15+l)->Fill(diamondCluster.getCharge(false));
								}
							}	//End of in hit bin
							hCellsColumnCheck1010->SetBinContent((k+1),(l+1),0);	//Empty column check histo
						}	//End of y cells loop
					}	//End of x cells loop
				}
			}
		}
	}	//End of if SeedCount





}
int TAnalysisOf3dDiamonds::YAlignmentFiducialCut(int nRegion){

	//cout<<"IN FID CUT"<<endl;
	//FidYAlignment.push_back(162);FidYAlignment.push_back(170);FidYAlignment.push_back(83);FidYAlignment.push_back(100);	//Fiducial Cut region end of 3DWH detector
	//cout<<nRegion<<endl;
	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	int Throw=0;

	if(fiducialValueX<FidCutYAlignment.at(nRegion)[5]||fiducialValueX>FidCutYAlignment.at(nRegion)[4])
		Throw =1;
	if(fiducialValueY<FidCutYAlignment.at(nRegion)[7]||fiducialValueY>FidCutYAlignment.at(nRegion)[6])
		Throw =1;

	//if(Throw==1)
		//cout<<"Thrown: "<<fiducialValueX<<",    "<<fiducialValueY<<endl;
	return Throw;
}
void TAnalysisOf3dDiamonds::DrawYAlignmentFidCutRegions() {

	hCombinedMeanChargeYAlignment->cd();
	for(int i=0;i<FidCutYAlignment.size();i++){
		FidCutChannelYAlignmentTBox.push_back(new TBox(FidCutYAlignment.at(i)[5], FidCutYAlignment.at(i)[7], FidCutYAlignment.at(i)[4], FidCutYAlignment.at(i)[6]));

		FidCutChannelYAlignmentTBox.at(i)->SetFillStyle(0);
		FidCutChannelYAlignmentTBox.at(i)->SetLineWidth(2);
		FidCutChannelYAlignmentTBox.at(i)->SetLineColor(kRed);
		FidCutChannelYAlignmentTBox.at(i)->Draw("same");
	}

}

Int_t* TAnalysisOf3dDiamonds::CellToChannel(int nCell, float nXcell) {
	Int_t* ChannelsToRead;
	Int_t Channel = -9999;
	Int_t nCellsPerRow = 11;
	Int_t nColumn = nCell / nCellsPerRow;
	Int_t startChannel = 85;
	Channel = startChannel + nColumn;
//	if(nCell<11 && nCell>=0)
//		Channel = 85;
//	if(nCell<22 && nCell>=11)
//		Channel = 86;
//	if(nCell<33 && nCell>=22)
//		Channel = 87;
//	if(nCell<44 && nCell>=33)
//		Channel = 88;
//	if(nCell<55 && nCell>=44)
//		Channel = 89;
//	if(nCell<66 && nCell>=55)
//		Channel = 90;
//	if(nCell<77 && nCell>=66)
//		Channel = 91;
//	if(nCell<88 && nCell>=77)
//		Channel = 92;
//	if(nCell<99 && nCell>=88)
//		Channel = 93;

//	int BadCells[] = {1,6,7,8,13,27,45,53,77,90,97};
	vector<int> BadCells = settings->getBadCells3D();
	int BadCell = 0;
	for(int i=0; i<BadCells.size(); i++){
		if(nCell ==  BadCells[i]){
			BadCell = 1;
			if(nCell<11 && nCell>=0){	//If cell on left edge
				Int_t A[] = {2,Channel,(Channel+1)};
				ChannelsToRead = A;
			}
			if(nCell<99 && nCell>=88){	//If cell on right edge
				Int_t A[] = {2,(Channel-1),Channel};
				ChannelsToRead = A;
			}
			if(nCell<88 && nCell>=11){	//If cell in the central region
				Int_t A[] = {3,(Channel-1),Channel,(Channel+1)};
				ChannelsToRead = A;
			}
		}
	}
	if(BadCell == 0){		//To account for charge sharing in edge 50um
		if(nXcell<5){	//If Xpos in cell left charge share
			if(nCell<11 && nCell>=0){	//If on leftmost column of detector
				Int_t A[] = {1,Channel};
				ChannelsToRead = A;
			}
			else{
				Int_t A[] = {2,(Channel-1),Channel};
				ChannelsToRead = A;
			}
		}
		if(nXcell>145){	//If Xpos in cell right charge share
			if(nCell<99 && nCell>=88){	//If on rightmost column of detector
				Int_t A[] = {1,Channel};
				ChannelsToRead = A;
			}
			else{
				Int_t A[] = {2,Channel,(Channel+1)};
				ChannelsToRead = A;
			}
		}
		if(nXcell>=5 && nXcell<=145){	//If Xpos in cell centre
			Int_t A[] = {1,Channel};
			ChannelsToRead = A;
		}
	}

	/*if(BadCell == 0){
		Int_t A[] = {1,Channel};
		ChannelsToRead = A;
	}
		*/
	return ChannelsToRead;
}

void TAnalysisOf3dDiamonds::DrawMetallisationGrid(TCanvas* nCanvas) {

	nCanvas->cd();
	//vector<TBox*> Grid;
	TBox* Grid;
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			float xLow = 2365 + i*150;
			float yLow = j*150;
			float xHigh = xLow+150;
			float yHigh = yLow+150;

			Grid = new TBox(xLow,yLow,xHigh,yHigh);
			Grid->SetFillStyle(0);
			Grid->SetLineWidth(1);
			Grid->SetLineColor(kBlack);
			Grid->Draw("same");
		}
	}
}

Float_t TAnalysisOf3dDiamonds::GetYPositionInDetSystem() {

	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	TPositionPrediction *predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	Float_t xPos = predictedPosition->getPositionX();	//Predicted positions in labframe
	Float_t yPos = predictedPosition->getPositionY();

	TDetectorAlignment* detectorAlignment = eventReader->getAlignment();
	Double_t xOffset = detectorAlignment->GetXOffset(subjectPlane);
	Double_t PhiOffset = detectorAlignment->GetPhiXOffset(subjectPlane);

	Float_t yDet = yPos*TMath::Cos(PhiOffset)+(xPos-xOffset)*TMath::Sin(PhiOffset);

	return yDet;

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
				float Temp2 = SortArrayPointer[i];
				nArray[i] = nArray[i+1];
				SortArrayPointer[i] = SortArrayPointer[i+1];
				nArray[i+1] = Temp;
				SortArrayPointer[i+1] = Temp2;
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

/*
void TAnalysisOf3dDiamonds::CreateDeadCellsProfile() {

	for(int i=0;i<4;i++){
		for(int j=0;j<3;j++){
			hCellsCharge.at(DeadCellsArrayPointer.at(i-1+j))->GetBinContent()
		}




		stringstream hDeadCellsChargeName; hDeadCellsChargeName<<"TotalCharge"<<(i)<<FileNameEnd;
		hDeadCellsCharge.push_back(new TH1F(hDeadCellsChargeName.str().c_str(),hDeadCellsChargeName.str().c_str(),90,0,450));
		stringstream hDeadCellsEventsName; hDeadCellsEventsName<<"Events"<<(i)<<FileNameEnd;
		hDeadCellsEvents.push_back(new TH1F(hDeadCellsEventsName.str().c_str(),hDeadCellsEventsName.str().c_str(),90,0,450));
	}


}
*/


void TAnalysisOf3dDiamonds::HitandSeedCount(TCluster* nCluster, int ni) {
	int Hit=0;int Seed=0;
	for (UInt_t i=0;i<nCluster->getClusterSize();i++){
		if(nCluster->isHit(i)) Hit++;
		if(nCluster->isSeed(i)) Seed++;
	}
	HitCount=(Hit-Seed);
	SeedCount=Seed;
}

void TAnalysisOf3dDiamonds::ClusterShape(TCluster* nCluster) {
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)&&!nCluster->isSeed(clPos)) cout<<" HIT ";
		if(nCluster->isSeed(clPos)) cout<<" SEED ";
	}
	cout<<endl;
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
			ClusterShape(nCluster);
		}
	}
}

void TAnalysisOf3dDiamonds::RemoveEdgeHits(TCluster* nCluster, int* nEdgeChannel) {
	//Removes edge detector hits
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)){
			if(nCluster->getChannel(clPos)>=nEdgeChannel[0]||nCluster->getChannel(clPos)<=nEdgeChannel[1]) return;
		}
	}
}

void TAnalysisOf3dDiamonds::ClusterChannels(TCluster* nCluster) {
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)){
			cout<<nCluster->getChannel(clPos)<<", ";
		}
		cout<<endl;
	}
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

int TAnalysisOf3dDiamonds::RemoveEdge3DClusters(TCluster* nCluster) {
	int Remove =0;
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)){
			if(nCluster->getChannel(clPos)>=93||nCluster->getChannel(clPos)<=85)
				Remove =1;
		}
	}
	return Remove;
}

int TAnalysisOf3dDiamonds::FiducialCut(int Detector){

	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	int Throw=0;

	//cout<<"X Less than: "<<FidCut.at(Detector)[5]<<" Greater than: "<<FidCut.at(Detector)[4]<<endl;
	//cout<<"Y Less than: "<<FidCut.at(Detector)[7]<<" Greater than: "<<FidCut.at(Detector)[6]<<endl;
	if(FiducialChannel){
		if(fiducialValueX<FidCut.at(Detector)[5]||fiducialValueX>FidCut.at(Detector)[4])
			Throw =1;
		if(fiducialValueY<FidCut.at(Detector)[7]||fiducialValueY>FidCut.at(Detector)[6])
			Throw =1;
	}
	if(Throw==1)
		cout<<"Detector: "<<Detector<<"Thrown: "<<fiducialValueX<<",    "<<fiducialValueY<<endl;
	return Throw;
}

void TAnalysisOf3dDiamonds::DrawFidCutRegions() {

	hCombinedMeanCharge->cd();

	for(int i=0;i<FidCut.size();i++){
		FidCutChannelTBox.push_back(new TBox(FidCut.at(i)[5], FidCut.at(i)[7], FidCut.at(i)[4], FidCut.at(i)[6]));
		FidCutChannelTBox.at(i)->SetFillStyle(0);
		FidCutChannelTBox.at(i)->SetLineWidth(2);
		FidCutChannelTBox.at(i)->SetLineColor(kRed);
		FidCutChannelTBox.at(i)->Draw("same");
	}

/*
	FidCudBoundMetric = (TH2F*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("FidCudBoundMetric");
	FidCudBoundMetric->SetFillStyle(3444);
	//FidCudBoundMetric->SetBins(80,90,170,60,60,120);
	FidCudBoundChannelCanvas = new TCanvas();
	cout<<"FitCut.size() = "<<FidCut.size()<<endl;
	cout<<"Y boundaries: "<<FidCut.at(0)[6]<<", "<<FidCut.at(0)[7]<<endl;
	cout<<"X boundaries: "<<FidCut.at(0)[4]<<", "<<FidCut.at(0)[5]<<endl;

	float FillXCh, FillYCh;
	for(int i=0;i<FidCut.size();i++){
		cout<<"1i is: "<<i<<endl;
		for(int j=0;j<(FidCut.at(i)[6]-FidCut.at(i)[7]);j++){ //Y boundaries
			FillYCh = FidCut.at(i)[7]+j+0.5; //Fill value of y
			cout<<"2i is: "<<i<<endl;

			for(int k=0;k<(FidCut.at(i)[4]-FidCut.at(i)[5]);k++){ //X boundaries
				cout<<"3i is: "<<i<<endl;

				FillXCh = FidCut.at(i)[5]+k; //Fill value of x
				cout<<"FillYCh = "<<FillYCh<<"FillXCh = "<<FillXCh<<endl;
				FidCudBoundMetric->Fill(FillXCh, FillYCh, 1);
			}
		}
	}
	// Highlight the maximum

	TBox b(FidCut.at(0)[5], FidCut.at(0)[7], FidCut.at(0)[4], FidCut.at(0)[6]);
	b.SetFillStyle(0);
	b.SetLineWidth(2);
	b.SetLineColor(kRed);
	b.Draw();
	FidCudBoundChannelCanvas->cd();
	FidCudBoundMetric->SetTitle("");
	FidCudBoundMetric->GetXaxis()->SetTitle("");
	FidCudBoundMetric->GetYaxis()->SetTitle("");
	FidCudBoundMetric->SetStats(kFALSE);
	FidCudBoundMetric->Draw();
	b.Draw();
	stringstream str; str<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"FidCudBoundMetric.png";
	FidCudBoundChannelCanvas->SaveAs(str.str().c_str());
	*/
}

void TAnalysisOf3dDiamonds::PredictChannelHit(TCluster* nCluster) {
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)&&!nCluster->isSeed(clPos)) cout<<" HIT ";
		if(nCluster->isSeed(clPos)) cout<<" SEED ";
	}
	cout<<endl;
}

