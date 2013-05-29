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

	FileNameEnd = "YAlignment17";
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

	//hFidCutXvsFidCutYvsChargeYAlignment
	stringstream hFidCutXvsFidCutYvsChargeYAlignmentName; hFidCutXvsFidCutYvsChargeYAlignmentName<<"hFidCutXvsFidCutYvsChargeYAlignment%%"<<FileNameEnd;
	hFidCutXvsFidCutYvsChargeYAlignment = new TH2D(hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),213,90,170,160,60,120);
	hFidCutXvsFidCutYvsChargeYAlignment->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsChargeYAlignment->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsChargeYAlignment->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignment
	hFidCutXvsFidCutYvsEventsYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsEventsYAlignment");

	//hFidCutXvsFidCutYvsMeanChargeYAlignment
	hFidCutXvsFidCutYvsMeanChargeYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsMeanChargeYAlignment");

	//hDetXvsDetY3D
	stringstream hDetXvsDetY3DName; hDetXvsDetY3DName<<"hFidCutXvsFidCutYvsChargeYAlignment%%"<<FileNameEnd;
	hDetXvsDetY3D = new TH2D(hDetXvsDetY3DName.str().c_str(),hDetXvsDetY3DName.str().c_str(),270,2365,3715,330,0,1650);
	hDetXvsDetY3D->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3D->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3D->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignment
	hDetXvsDetY3DvsEvents = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DvsEvents");

	//hFidCutXvsFidCutYvsMeanChargeYAlignment
	hDetXvsDetY3DMeanCharge = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DMeanCharge");

	//hDetXvsDetY3DRebinned
	stringstream hDetXvsDetY3DRebinnedName; hDetXvsDetY3DRebinnedName<<"hFidCutXvsFidCutYvsChargeRebinnedYAlignment%%"<<FileNameEnd;
	hDetXvsDetY3DRebinned = new TH2D(hDetXvsDetY3DRebinnedName.str().c_str(),hDetXvsDetY3DRebinnedName.str().c_str(),9,2365,3715,11,0,1650);
	hDetXvsDetY3DRebinned->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinned->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinned->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignmentRebinned
	hDetXvsDetY3DvsEventsRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DvsEventsRebinned");

	//hFidCutXvsFidCutYvsMeanChargeYAlignmentRebinned
	hDetXvsDetY3DMeanChargeRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DMeanChargeRebinned");
	//hDetXvsDetY3DMeanChargeRebinned->SetBins(9,2365,3715,11,0,1650,12,0,1200);

	//hDetXvsDetY3DOverview
	stringstream hDetXvsDetY3DOverviewName; hDetXvsDetY3DOverviewName<<"hDetXvsDetY3DOverview%%"<<FileNameEnd;
	hDetXvsDetY3DOverview = new TH2D(hDetXvsDetY3DOverviewName.str().c_str(),hDetXvsDetY3DOverviewName.str().c_str(),9,2365,3715,11,0,1650);
	hDetXvsDetY3DOverview->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DOverview->GetYaxis()->SetTitle("Ydet (um)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();

	//hBinnedMeanCharge
	stringstream hBinnedMeanChargeName; hBinnedMeanChargeName<<"hBinnedMeanCharge%%"<<FileNameEnd;
	hBinnedMeanCharge = new TH1F(hBinnedMeanChargeName.str().c_str(),hBinnedMeanChargeName.str().c_str(),9,400,1300);
	hBinnedMeanCharge->SetTitle(hBinnedMeanChargeName.str().c_str());
	hBinnedMeanCharge->GetXaxis()->SetTitle("MeanCharge");
	hBinnedMeanCharge->GetYaxis()->SetTitle("Entries");

	//Cell Overlay
	//hCellsOverlayed
	hCellsOverlayedCharge = new TH2D("OverlayedCharge","OverlayedCharge",30,0,150,30,0,150);
	hCellsOverlayedEvents = new TH2D("OverlayedEvents","OverlayedEvents",30,0,150,30,0,150);
	hCellsOverlayedEventsNoColumns = new TH2D("hCellsOverlayedEventsNoColumns","hCellsOverlayedEventsNoColumns",30,0,150,30,0,150);
	stringstream hCellsOverlayedMeanChargeName; hCellsOverlayedMeanChargeName<<"hCellsOverlayedMeanCharge%%"<<FileNameEnd;
	hCellsOverlayedMeanCharge = new TH2D(hCellsOverlayedMeanChargeName.str().c_str(),hCellsOverlayedMeanChargeName.str().c_str(),30,0,150,30,0,150);
	hCellsOverlayedMeanCharge->GetXaxis()->SetTitle("Xdet (um)");
	hCellsOverlayedMeanCharge->GetYaxis()->SetTitle("Ydet (um)");
	hCellsOverlayedMeanCharge->GetZaxis()->SetTitle("Charge ADC");

	//hCellsColumnCheckName		//Cell histogram used to check whether hit is in column
	stringstream hCellsColumnCheckName; hCellsColumnCheckName<<"hCellsColumnCheck%%"<<FileNameEnd;
	hCellsColumnCheck = new TH2D(hCellsColumnCheckName.str().c_str(),hCellsColumnCheckName.str().c_str(),30,0,150,30,0,150);

	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			/*
			float xLow = 2365 + i*150;
			float yLow = j*150;
			float xHigh = xLow+150;
			float yHigh = yLow+150;
				*/
			stringstream hCellsChargeName; hCellsChargeName<<"TotalCharge"<<((i*11)+j)<<FileNameEnd;
			hCellsCharge.push_back(new TH2D(hCellsChargeName.str().c_str(),hCellsChargeName.str().c_str(),30,0,150,30,0,150));
			stringstream hCellsEventsName; hCellsEventsName<<"Events"<<((i*11)+j)<<FileNameEnd;
			hCellsEvents.push_back(new TH2D(hCellsEventsName.str().c_str(),hCellsEventsName.str().c_str(),30,0,150,30,0,150));

			stringstream hCellsLandauName; hCellsLandauName<<"hLandauCell%%"<<((i*11)+j)<<FileNameEnd;
			hCellsLandau.push_back(new TH1F(hCellsLandauName.str().c_str(),hCellsLandauName.str().c_str(),256,0,2800));
			stringstream hCellsLandauNoColumnName; hCellsLandauNoColumnName<<"hCellsLandauNoColumn%%"<<((i*11)+j)<<FileNameEnd;
			hCellsLandauNoColumn.push_back(new TH1F(hCellsLandauNoColumnName.str().c_str(),hCellsLandauNoColumnName.str().c_str(),256,0,2800));

			stringstream hCellsEventsNoColumnName; hCellsEventsNoColumnName<<"hCellsEventsNoColumn%%"<<((i*11)+j)<<FileNameEnd;
			hCellsEventsNoColumn.push_back(new TH2D(hCellsEventsNoColumnName.str().c_str(),hCellsEventsNoColumnName.str().c_str(),30,0,150,30,0,150));
		}
	}
	//hCellsLandauGraded    &&    hCellsLandauGradedNoColumn
	for(int i=0;i<12;i++){	//Group CellsLandaus within same ranges together. 0-100; 100-200; -> 1100-1200;
		stringstream hCellsLandauGradedName; hCellsLandauGradedName<<"hLandauCellsGraded%%"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGraded.push_back(new TH1F(hCellsLandauGradedName.str().c_str(),hCellsLandauGradedName.str().c_str(),256,0,2800));
		stringstream hCellsLandauGradedNoColumnName; hCellsLandauGradedNoColumnName<<"hLandauCellsGradedNoColumn%%"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGradedNoColumn.push_back(new TH1F(hCellsLandauGradedNoColumnName.str().c_str(),hCellsLandauGradedNoColumnName.str().c_str(),256,0,2800));
	}

	//hEdgeCharge
	stringstream hEdgeChargeName; hEdgeChargeName<<"hEdgeCharge%%"<<FileNameEnd;
	hEdgeCharge = new TH1F(hEdgeChargeName.str().c_str(),hEdgeChargeName.str().c_str(),250,3100,4100);
	hEdgeCharge->SetTitle(hEdgeChargeName.str().c_str());
	hEdgeCharge->GetXaxis()->SetTitle("X Diamond (um)");
	hEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hEdgeChargeEvents = (TH1F*)hEdgeCharge->Clone("hEdgeChargeEvents");
	hEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

	hEdgeMeanCharge = (TH1F*)hEdgeCharge->Clone("hEdgeMeanCharge");
	hEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	//hEdgeCharge
	stringstream hyEdgeChargeName; hyEdgeChargeName<<"hyEdgeCharge%%"<<FileNameEnd;
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
	stringstream hDeadCellName; hDeadCellName<<"hDeadCell%%"<<FileNameEnd;
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
	hFidCutXvsFidCutYvsMeanChargeYAlignment->Draw("COLZ");
	/*DrawFidCutRegions(); //Draw Fiducial Cut Regions
		 */
	stringstream str1; str1<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanChargeYAlignmentNoFidDrawn"<<"%%"<<FileNameEnd<<".png";
	hCombinedMeanChargeYAlignment->SaveAs(str1.str().c_str());

	DrawYAlignmentFidCutRegions(); //Draw Fiducial Cut Regions
	stringstream str11; str11<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hFidCutXvsFidCutYvsMeanCharge"<<"%%"<<FileNameEnd<<".png";
	hCombinedMeanChargeYAlignment->SaveAs(str11.str().c_str());

	//For h3DdetMeanCharge
	h3DdetMeanCharge = new TCanvas();
	h3DdetMeanCharge->cd();
	*hDetXvsDetY3DMeanCharge = (*hDetXvsDetY3D/(*hDetXvsDetY3DvsEvents));
	hDetXvsDetY3DMeanCharge->SetEntries(hDetXvsDetY3DvsEvents->Integral());
	hDetXvsDetY3DMeanCharge->Draw("COLZ");
	DrawMetallisationGrid();
	stringstream str12; str12<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetMeanCharge"<<"%%"<<FileNameEnd<<".png";
	h3DdetMeanCharge->SaveAs(str12.str().c_str());

	//RebinnedMeanCharge
	h3DdetMeanChargeRebinned = new TCanvas();
	h3DdetMeanChargeRebinned->cd();
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	*hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DRebinned/(*hDetXvsDetY3DvsEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DvsEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->Draw("COLZ");
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	DrawMetallisationGrid();
	stringstream str121; str121<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetMeanChargeRebinned"<<"%%"<<FileNameEnd<<".png";
	h3DdetMeanChargeRebinned->SaveAs(str121.str().c_str());
	/*h3DdetEventsRebinned = new TCanvas();
	//h3DdetEventsRebinned->cd();
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");			//("COLZ");
	//DrawMetallisationGrid();
	stringstream str122; str122<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"h3DdetEventsRebinned"<<"%%"<<FileNameEnd<<".png";
	h3DdetMeanChargeRebinned->SaveAs(str122.str().c_str());
		*/

	//Cell Overlay
	//int BadCells[] = {0,1,2,3,4,5,6,7,8,9,10,13,19,24,27,34,35,44,45,48,53,57,64,68,75,77,78,88,90,97};
	cout<<"CellArraySize is: "<<hCellsCharge.size()<<endl;
	//for(int i=0;i<hCellsCharge.size();i++){
	for(int i=0;i<9;i++){
		for(int j=0;j<11;j++){
			histSaver->SaveHistogram(hCellsLandau.at(i*11+j));	//Save each cells Landau;
			//cout<<"The content of bin: "<<(i*11+j)<<" is: " <<(hDetXvsDetY3DvsEventsRebinned->GetBinContent(i+1,j+1))<<endl;
			hBinnedMeanCharge->Fill(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1));
			//for(int k=0;k<30;k++){
			//if(i!=BadCells[k]){
			if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)>1000)&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)<1100)){
				hCellsOverlayedCharge->Add(hCellsCharge.at(i*11+j));
				hCellsOverlayedEvents->Add(hCellsEvents.at(i*11+j));
				hCellsOverlayedEventsNoColumns->Add(hCellsEventsNoColumn.at(i*11+j));
				cout<<"Cell: "<<(i*11+j)<<" is overlayed."<<endl;

				hDetXvsDetY3DOverview->SetBinContent(i+1,j+1,1);
				//cout<<"Bin content is: "<<hDetXvsDetY3DOverview->GetBinContent(1,1)<<endl;
				//histSaver->SaveHistogram(hCellsCharge.at(i));
				//gStyle->SetOptStat("irm");
				hCellsEvents.at(i*11+j)->SetEntries(hCellsEvents.at(i*11+j)->Integral()); //Set the number of entries to the integral of hEvents
				histSaver->SaveHistogram(hCellsEvents.at(i*11+j));
				cout<<"The integral of Events"<<(i*11+j)<<" is: "<<hCellsEvents.at(i*11+j)->Integral()<<endl;
			}
			for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
				if((hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)>(k*100))&&(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(i+1,j+1)<((k+1)*100))){
					hCellsLandauGraded.at(k)->Add(hCellsLandau.at(i*11+j));
					hCellsLandauGradedNoColumn.at(k)->Add(hCellsLandauNoColumn.at(i*11+j));
				}
			}
		//}
		}
	}
	hCellsOverlayedCanvas = new TCanvas();
	hCellsOverlayedCanvas->cd();
	*hCellsOverlayedMeanCharge = (*hCellsOverlayedCharge/(*hCellsOverlayedEvents));
	hCellsOverlayedMeanCharge->SetEntries(hCellsOverlayedEvents->Integral());
	hCellsOverlayedMeanCharge->Draw("COLZ");
	stringstream str13; str13<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedMeanCharge"<<"%%"<<FileNameEnd<<".png";
	hCellsOverlayedCanvas->SaveAs(str13.str().c_str());

	hCellsOverlayedEvents->Draw("sameTEXT");
	stringstream str133; str133<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedMeanChargeWithEntries"<<"%%"<<FileNameEnd<<".png";
	hCellsOverlayedCanvas->SaveAs(str133.str().c_str());
	histSaver->SaveHistogram(hBinnedMeanCharge);

	//hCellsOverlayedEventsNoColumns
	hCellsEventsNoColumnCanvas = new TCanvas();
	hCellsEventsNoColumnCanvas->cd();
	hCellsOverlayedEventsNoColumns->Draw("COLZ");
	hCellsOverlayedEventsNoColumns->Draw("sameTEXT");
	stringstream str14; str14<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hCellsOverlayedEventsNoColumns"<<"%%"<<FileNameEnd<<".png";
	hCellsEventsNoColumnCanvas->SaveAs(str14.str().c_str());

	//hCellsLandauGraded			//Save each of the histograms
	for(int k=0;k<12;k++){	// to fill each of the sub ranges in 100's of ADC
		histSaver->SaveHistogram(hCellsLandauGraded.at(k));
		histSaver->SaveHistogram(hCellsLandauGradedNoColumn.at(k));
	}

	//hDetXvsDetY3DOverview
	hOverview = new TCanvas();
	hOverview->cd();
	hDetXvsDetY3DOverview->Draw("COLZ");
	stringstream str131; str131<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hOverview"<<"%%"<<FileNameEnd<<".png";
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

	stringstream hDeadCellMeanChargeCanvasName; hDeadCellMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hDeadCellMeanCharge"<<"%%"<<FileNameEnd<<".png";
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
	stringstream hEdgeMeanChargeCanvasName; hEdgeMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hEdgeMeanCharge"<<"%%"<<FileNameEnd<<".png";
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
	stringstream hyEdgeMeanChargeCanvasName; hyEdgeMeanChargeCanvasName<<"/Users/iainhaughton/3D_diamond_analysis/output/17107/3dDiamondAnalysis/"<<"hyEdgeMeanCharge"<<"%%"<<FileNameEnd<<".png";
	hyEdgeMeanChargeCanvas->SaveAs(hyEdgeMeanChargeCanvasName.str().c_str());

}
void TAnalysisOf3dDiamonds::YAlignment() {

	if(!eventReader->isValidTrack())
		return;
	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	if(!eventReader->isInFiducialCut())	//This is a larger fiducial cut around silicon
		return;

	if(eventReader->getNDiamondClusters()==0||eventReader->getNDiamondClusters()==2||eventReader->getNDiamondClusters()==3)	//Only single cluster events
		return;


	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster())
		return;

	HitandSeedCount(&diamondCluster, 0);

	//if(HitCount==0&&SeedCount==1){

		ClusterShape(&diamondCluster);

		TPositionPrediction *predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
		if(predictedPosition->getChi2X()>100||predictedPosition->getChi2Y()>100)
			return;

		Float_t xPos = predictedPosition->getPositionX();	//Predicted positions in labframe
		Float_t yPos = predictedPosition->getPositionY();
		float fiducialValueX= eventReader->getFiducialValueX();
		float fiducialValueY= eventReader->getFiducialValueY();

		/*if(YAlignmentFiducialCut())
			return;
				*/

		hFidCutXvsFidCutYvsChargeYAlignment->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
		hFidCutXvsFidCutYvsEventsYAlignment->Fill(fiducialValueX,fiducialValueY,1);

		for(int i=0;i<2;i++){
			if(i==0){
				if(!YAlignmentFiducialCut(0)){
					//cout<<"HERE"<<endl;
					Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
					hEdgeChargeEvents->Fill(positionInDetSystemMetric);
					hEdgeCharge->Fill(positionInDetSystemMetric, diamondCluster.getCharge(false));
				}
			}
			if(i==1){
				if(!YAlignmentFiducialCut(1)){
					//cout<<"HERE"<<endl;
					/*TDetectorAlignment* detectorAlignment = eventReader->getAlignment();
					Double_t PhiOffset = detectorAlignment->GetPhiXOffset(subjectPlane);
					cout<<"The PhiOffset is: "<<PhiOffset<<endl;
					Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
					Float_t positionYInDetSystemMetric = yPos - positionInDetSystemMetric*TMath::Sin(PhiOffset);
					 */
					hyEdgeChargeEvents->Fill(GetYPositionInDetSystem());
					hyEdgeCharge->Fill(GetYPositionInDetSystem(), diamondCluster.getCharge(false));
				}
			}
		}	//End of for

		Float_t Xdet = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos);
		Float_t Ydet = GetYPositionInDetSystem()-3890;
		hDetXvsDetY3D->Fill(Xdet,Ydet,diamondCluster.getCharge(false));
		hDetXvsDetY3DvsEvents->Fill(Xdet,Ydet,1);
		hDetXvsDetY3DRebinned->Fill(Xdet,Ydet,diamondCluster.getCharge(false));
		hDetXvsDetY3DvsEventsRebinned->Fill(Xdet,Ydet,1);

		if(Xdet>2665&&Xdet<2815){	//Dead cell chosen to be: (2,5) in cell space
			hDeadCell->Fill(Ydet, diamondCluster.getCharge(false));
			hDeadCellEvents->Fill(Ydet, 1);
		}

		//	To fill multiple cell histograms
		//for(int i=0;i<hCellsCharge.size();i++){
		float hEntries = 0;
		for(int i=0;i<9;i++){
			for(int j=0;j<11;j++){
				hEntries = hCellsEvents.at(i*11+j)->Integral();
				float xminus = 2365+i*150;		//2365 is the start of the 3D detector in x
				float yminus = j*150;
				hCellsCharge.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),diamondCluster.getCharge(false));
				hCellsEvents.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
				if(hEntries != hCellsEvents.at(i*11+j)->Integral()){	//If entries in this histogram range have increased fill corresponding hCellLandau
					hCellsLandau.at(i*11+j)->Fill(diamondCluster.getCharge(false));
					hCellsColumnCheck->Fill((Xdet-xminus),(Ydet-yminus),1);	//Fill Column Check histo
					for(int k=0;k<30;k++){		//looping over number of cell bins in x
						for(int l=0;l<30;l++){		//looping over number of cell bins in y
							int ColumnEntriesX[] = {1,1,2,2,1,1,2,2,16,16,17,17,29,29,30,30,29,29,30,30};	//Column cell x coordinates
							int ColumnEntriesY[] = {1,2,1,2,29,30,29,30,15,16,15,16,1,2,1,2,29,30,29,30};	//Column cell y coordinates
							if(hCellsColumnCheck->GetBinContent((k+1),(l+1))==1){	//True when in hit bin
								int inColumn =0;
								for(int m=0;m<20;m++){	//Loop over number of column cells
									if((k+1)==ColumnEntriesX[m]&&(l+1)==ColumnEntriesY[m]) inColumn =1;
								}
								if(inColumn==0){
									hCellsLandauNoColumn.at(i*11+j)->Fill(diamondCluster.getCharge(false));
									hCellsEventsNoColumn.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
								}
							}	//End of in hit bin
							hCellsColumnCheck->SetBinContent((k+1),(l+1),0);	//Empty column check histo
						}	//End of y cells loop
					}	//End of x cells loop
				}
			}
		}
	//}	//End of if SeedCount

/*
	//Predict Channel hit
	Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
	//cout<<"PrintOut position in Diamond plane: "<<positionInDetSystemMetric<<endl;
	//Float_t positionInDetSystemChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,positionInDetSystemMetric);
	hEdgeChargeEvents->Fill(positionInDetSystemMetric);
	hEdgeCharge->Fill(positionInDetSystemMetric, diamondCluster.getCharge(false));
		*/

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
void TAnalysisOf3dDiamonds::DrawMetallisationGrid() {

	h3DdetMeanCharge->cd();
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
			Grid->SetLineColor(kRed);
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

