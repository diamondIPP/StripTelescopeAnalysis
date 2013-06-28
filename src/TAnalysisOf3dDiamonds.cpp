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
	else exit(-1);
	verbosity = settings->getVerbosity();
	predictedPosition=0;
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
	clusteredAnalysis = new TCellAnalysisClass(settings);
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
//	int Strip[] = {26,38};
//	int threeDnH[] = {56,62};
//	int threeDwH[] = {86,92};

	//int Strip[] = {32,32};
	//int threeDnH[] = {59,59};
	//int threeDwH[] = {88,88};

//	Detector.push_back(Strip); //Strip
//	Detector.push_back(threeDnH); //3D no holes
//	Detector.push_back(threeDwH); //3D with holes
//	//Fiducial Cut Regions       (xHigh, xLow, yHigh, yLow,     xChHigh, xChLow, yChHigh, yChLow)
//	int StripFid[] =             {102,   92,   110,   68,       104,     93,     96,      68};
//	int threeDnHFid[] =          {102,   92,   110,   68,       130,     112,    105,     98};
//	int threeDwHFid[] =          {102,   92,   110,   68,       161,     141,    107,     70};

	/*TFiducialCut* fidCutStrip = new TFiducialCut(0,92,106,75,95);		//These are the fiducial regions used to plot events within the detector.
	TFiducialCut* fidCut3DNoHoles = new TFiducialCut(1,108,135,73,104);
	TFiducialCut* fidCut3DWiHoles = new TFiducialCut(2,137,166,73,104);
	FidCut.push_back(fidCutStrip);
	FidCut.push_back(fidCut3DNoHoles);
	FidCut.push_back(fidCut3DWiHoles);

	fidCuts.addFiducialCut(fidCutStrip);
	fidCuts.addFiducialCut(fidCut3DNoHoles);
	fidCuts.addFiducialCut(fidCut3DWiHoles);

	cout<<settings->get3dFidCuts()->getFidCut(1)->GetXLow()<<"    "<<settings->get3dFidCuts()->getFidCut(1)->GetXHigh()<<endl;
	cout<<settings->get3dFidCuts()->getFidCut(2)->GetXLow()<<"    "<<settings->get3dFidCuts()->getFidCut(2)->GetXHigh()<<endl;
	cout<<settings->get3dFidCuts()->getFidCut(3)->GetXLow()<<"    "<<settings->get3dFidCuts()->getFidCut(3)->GetXHigh()<<endl;
*/

	FiducialMetric=0;			// Which Fiducial cut to apply
	FiducialChannel=1;


	//For YAlignment Fiducial cuts
	//Fiducial Cut Regions    (xHigh, xLow, yHigh, yLow,     xChHigh, xChLow, yChHigh, yChLow)
	int xEdge[] =             {102,   92,   110,   68,       170,     162,    100,     83};
	int yEdge[] =             {102,   92,   110,   68,       164,     140,    112,     103};
	TFiducialCut* fidCutYAlignmentX = new TFiducialCut(0,162,170,83,100);
	TFiducialCut* fidCutYAlignmentY = new TFiducialCut(0,140,164,103,112);

	FidCutYAlignment.push_back(xEdge);
	FidCutYAlignment.push_back(yEdge);
	fidCutYAlignmentX->Print();
	fidCutYAlignmentY->Print();
//	for(int i=0;i<FidCutYAlignment.size();i++){
//		for(int j=0; j<8;j++){
//			cout<<FidCutYAlignment.at(i)[j]<<endl;
//		}
//	}

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
	//Moved initialiseHistos to here
	/*initialise3DGridReference();
	initialise3DYAlignmentHistos();
	initialise3DOverviewHistos();
	initialise3D2DLandauAndClustersizeHistos();
	initialise3DCellOverlayHistos();
		*/
	initialiseHistos();
	 	 //*/
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		eventReader->LoadEvent(nEvent);
		analyseEvent();
		// *
		// */
		//YAlignment();
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

	//saveYAlignmentHistos();
	saveHistos();
	 	 //*/
}

void TAnalysisOf3dDiamonds::initialiseHistos() {
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

	for(int i=0; i<settings->diamondPattern.getNIntervals(); i++){
		pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
		//hLandau
		stringstream hLandauName; hLandauName<<"hLandau%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hLandau.push_back(new TH1F(hLandauName.str().c_str(),hLandauName.str().c_str(),256,0,2800));
		hLandau.at(i)->SetTitle(hLandauName.str().c_str());
		hLandau.at(i)->GetXaxis()->SetTitle("PH of diamond cluster");
		hLandau.at(i)->GetYaxis()->SetTitle("number of entries #");

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

		/*
		//hPHvsPredictedChannel
		stringstream hPHvsPredictedChannelName; hPHvsPredictedChannelName<<"hPHvsPredictedChannel%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hPHvsPredictedChannel.push_back(new TH2F(hPHvsPredictedChannelName.str().c_str(),hPHvsPredictedChannelName.str().c_str(),150,0,2900,320,20,100));
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsPredictedChannel.at(i)->GetYaxis()->SetTitle("XPos(Predicted channel)");
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetRangeUser(1,2600);
		//hPHvsPredictedChannel.at(i)->SetMaximum(3000);
		//hPHvsPredictedChannel.at(i)->SetMinimum(0);
		 	 */
		/*
		//hPHvsPredictedXPos
		stringstream hPHvsPredictedXPosName; hPHvsPredictedXPosName<<"hPHvsPredictedXPos%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hPHvsPredictedXPos.push_back(new TH2F(hPHvsPredictedXPosName.str().c_str(),hPHvsPredictedXPosName.str().c_str(),150,0,2900,50,5800,10200));
		hPHvsPredictedXPos.at(i)->GetXaxis()->SetTitle("Charge in ADC counts");
		hPHvsPredictedXPos.at(i)->GetYaxis()->SetTitle("XPos(um)");
		hPHvsPredictedChannel.at(i)->GetXaxis()->SetLimits(0.,2900.);
		hPHvsPredictedXPos.at(i)->GetXaxis()->SetRange(1,2600);
		//hPHvsPredictedXPos.at(i)->SetMaximum(3000);
		//hPHvsPredictedXPos.at(i)->SetMinimum(0);
		 	 */

		//hHitandSeedCount
		stringstream hHitandSeedCountName; hHitandSeedCountName<<"hHitandSeedCount%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hHitandSeedCount.push_back(new TH2F(hHitandSeedCountName.str().c_str(),hHitandSeedCountName.str().c_str(),10,0,10,10,0,10));
		hHitandSeedCount.at(i)->GetXaxis()->SetTitle("Hit Count");
		hHitandSeedCount.at(i)->GetYaxis()->SetTitle("Seed Count");

		//hChi2XChi2Y
		stringstream hChi2XChi2YName; hChi2XChi2YName<<"hChi2XChi2Y%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hChi2XChi2Y.push_back(new TH2F(hChi2XChi2YName.str().c_str(),hChi2XChi2YName.str().c_str(),60,0,60,60,0,60));
		hChi2XChi2Y.at(i)->GetXaxis()->SetTitle("Chi2X");
		hChi2XChi2Y.at(i)->GetYaxis()->SetTitle("Chi2Y");

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

		/*//hFidCutXvsFidCutYvsCharge		For TH3F
		stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%"<<firstCh<<"-"<<channels.second<<"%%";
		hFidCutXvsFidCutYvsCharge.push_back(new TH3F(hFidCutXvsFidCutYvsChargeName.str().c_str(),hFidCutXvsFidCutYvsChargeName.str().c_str(),160,90,170,120,60,120,30,0,3000));
		hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
		hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
		hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");
		*/

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

	/*
	//To calculate Efficiency
	//hFidCutXvsFidCutYvsPredictedEvents	For TH2D
	stringstream hFidCutXvsFidCutYvsPredictedEventsName; hFidCutXvsFidCutYvsPredictedEventsName<<"hFidCutXvsFidCutYvsPredictedEvents"<<FileNameEnd;
	hFidCutXvsFidCutYvsPredictedEvents = new TH2D(hFidCutXvsFidCutYvsPredictedEventsName.str().c_str(),hFidCutXvsFidCutYvsPredictedEventsName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
	hFidCutXvsFidCutYvsPredictedEvents->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsPredictedEvents->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsPredictedEvents->GetZaxis()->SetTitle("Predicted Events");

	//hFidCutXvsFidCutYvsSeenEventsName	For TH2D
	stringstream hFidCutXvsFidCutYvsSeenEventsName; hFidCutXvsFidCutYvsSeenEventsName<<"hFidCutXvsFidCutYvsSeenEvents"<<FileNameEnd;
	hFidCutXvsFidCutYvsSeenEvents = new TH2D(hFidCutXvsFidCutYvsSeenEventsName.str().c_str(),hFidCutXvsFidCutYvsSeenEventsName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
	hFidCutXvsFidCutYvsSeenEvents->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsSeenEvents->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsSeenEvents->GetZaxis()->SetTitle("Seen Events");

	//hEfficiency
	hEfficiency = (TH2D*)hFidCutXvsFidCutYvsSeenEvents->Clone("hEfficiency");
		*/

	//hXPosvsYPosvsMeanCharge		For TH2D
	/*
	stringstream hXPosvsYPosvsMeanChargeName; hXPosvsYPosvsMeanChargeName<<"hXPosvsYPosvsMeanCharge";
	hXPosvsYPosvsMeanCharge = new TH2D(hXPosvsYPosvsMeanChargeName.str().c_str(),hXPosvsYPosvsMeanChargeName.str().c_str(),213,5800,9600,160,3400,5900);
	hXPosvsYPosvsMeanCharge->GetXaxis()->SetTitle("XPos");
	hXPosvsYPosvsMeanCharge->GetYaxis()->SetTitle("YPos");
	hXPosvsYPosvsMeanCharge->GetZaxis()->SetTitle("Charge ADC");
	hXPosvsYPosvsCharge = (TH2D*)hXPosvsYPosvsMeanCharge->Clone("hXPosvsYPosvsCharge");
	hXPosvsYPosvsEvents = (TH2D*)hXPosvsYPosvsMeanCharge->Clone("hXPosvsYPosvsEvents");
			*/
	//hFidCutXvsFidCutYvsCharge
	/*stringstream hFidCutXvsFidCutYvsChargeName; hFidCutXvsFidCutYvsChargeName<<"hFidCutXvsFidCutYvsCharge%%";   //<<firstCh<<"-"<<channels.second<<"%%";
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
	for(int i=0;i<settings->diamondPattern.getNIntervals();i++){
		hEventsvsChannelCombined->Add(hEventsvsChannel.at(i));
	}
	histSaver->SaveHistogram(hEventsvsChannelCombined);
	//histSaver->SaveHistogram(hDoubleClusterPos);
	//histSaver->SaveHistogram(hDoubleClusterPos0);
	//histSaver->SaveHistogram(hDoubleClusterPos1);
	//histSaver->SaveHistogram(hLandauCluster1);
	//histSaver->SaveHistogram(hLandauCluster2);
	//histSaver->SaveHistogram(hLandauDoubleCombined);
	/*
	cDoubleCluster = new TCanvas("cDoubleClusterCombinedPos","cDoubleClusterCombinedPos");
	cDoubleCluster->cd();
	hDoubleClusterPos0->Draw();
	hDoubleClusterPos1->Draw("samehist");
	histSaver->SaveCanvas(cDoubleCluster);
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
		histSaver->SaveHistogram(hLandau.at(i));

		/*//hPredictedPositionDiamondHit
		stringstream hPredictedPositionDiamondHitName; hPredictedPositionDiamondHitName<<"hPredictedPositionDiamondHit%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
		hPredictedPositionDiamondHit.push_back(HistogrammSaver::CreateScatterHisto(hPredictedPositionDiamondHitName.str().c_str(),*vecXPredicted.at(i),*vecYPredicted.at(i),512));
		hPredictedPositionDiamondHit.at(i)->SetTitle(hPredictedPositionDiamondHitName.str().c_str());
		hPredictedPositionDiamondHit.at(i)->GetXaxis()->SetTitle("predicted x Position");
		hPredictedPositionDiamondHit.at(i)->GetYaxis()->SetTitle("predicted y Position");
		histSaver->SaveHistogram(hPredictedPositionDiamondHit.at(i));
			*/
		histSaver->SaveHistogram(hHitandSeedCount.at(i));
		histSaver->SaveHistogram(hPHvsChannel.at(i));
		//histSaver->SaveHistogram(hPHvsPredictedXPos.at(i));
		histSaver->SaveHistogram(hChi2XChi2Y.at(i));
		histSaver->SaveHistogram(hFidCutXvsFidCutY.at(i));
		//histSaver->SaveHistogram(hPHvsPredictedChannel.at(i));
		//histSaver->SaveHistogram(hFidCutXvsFidCutYvsCharge.at(i));

		/*//hFidCutXvsFidCutYvsCharge
		ptrCanvas.at(i)->cd();
		hFidCutXvsFidCutYvsCharge.at(i)->Draw("COLZ");
		TString hName = TString::Format("cFidCutXvsFidCutYvsCharge_%d_%d",channels.first,channels.second);
		ptrCanvas.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvas.at(i));

		//hFidCutXvsFidCutYvsEvents
		ptrCanvasEvents.at(i)->cd();
		hFidCutXvsFidCutYvsEvents.at(i)->Draw("COLZ");
		hName  = TString::Format("cFidCutXvsFidCutYvsChargeClone_%d_%d",channels.first,channels.second);
		ptrCanvasEvents.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvasEvents.at(i));
			*/

		//hFidCutXvsFidCutYvsMeanCharge
		ptrCanvasMean.at(i)->cd();
		*hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
		hFidCutXvsFidCutYvsMeanCharge.at(i)->SetEntries(hFidCutXvsFidCutYvsEvents.at(i)->Integral());
		hFidCutXvsFidCutYvsMeanCharge.at(i)->Draw("COLZ");
		TString hName  = TString::Format("cFidCutXvsFidCutYvsMeanCharge_%d_%d",channels.first,channels.second);
		ptrCanvasMean.at(i)->SetName(hName);
		histSaver->SaveCanvas(ptrCanvasMean[i]);
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
	hCombinedMeanCharge->SetName("hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn");
	histSaver->SaveCanvas(hCombinedMeanCharge);
	settings->get3dFidCuts()->drawFiducialCutsToCanvas(hCombinedMeanCharge);
	hCombinedMeanCharge->SetName("hFidCutXvsFidCutYvsMeanChargeAllDetectors");
	histSaver->SaveCanvas(hCombinedMeanCharge);

	/*
	//For efficiency
	//Predicted Events
	cPredictedEvents = new TCanvas("cFidCutXvsFidCutYvsPredictedEvents","cFidCutXvsFidCutYvsPredictedEvents");
	cPredictedEvents->cd();
	hFidCutXvsFidCutYvsPredictedEvents->Draw("COLZ");
	histSaver->SaveCanvas(cPredictedEvents);

	//Seen Events
	cSeenEvents = new TCanvas("cFidCutXvsFidCutYvsSeenEvents","cFidCutXvsFidCutYvsSeenEvents");
	cSeenEvents->cd();
	hFidCutXvsFidCutYvsSeenEvents->Draw("COLZ");
	histSaver->SaveCanvas(cSeenEvents);

	//Efficiency
	TCanvas* cEfficiencyMap = new TCanvas("cEfficiencyMap","cEfficiencyMap");
	*hEfficiency = (*hFidCutXvsFidCutYvsSeenEvents/(*hFidCutXvsFidCutYvsPredictedEvents));
	hEfficiency->Draw("COLZ");
	histSaver->SaveCanvas(cEfficiencyMap);
	 	 */

	for(int i=0;i<7;i++){
		stringstream cClustersName; cClustersName<<"cFidCutXvsFidCutYClusters"<<i;
		cClusters.push_back(new TCanvas(cClustersName.str().c_str(),cClustersName.str().c_str()));
		cClusters.at(i)->cd();
		hFidCutXvsFidCutYClusters.at(i)->Draw("COLZ");
		histSaver->SaveCanvas(cClusters.at(i));
	}

	/*
	//For hXPosvsYPosvsMeanCharge
	TCanvas* cXPosvsYPosvsMean = new TCanvas("cXPosvsYPosvsMeanCharge","cXPosvsYPosvsMeanCharge");
	cXPosvsYPosvsMean->cd();
	*hXPosvsYPosvsMeanCharge = (*hXPosvsYPosvsCharge/(*hXPosvsYPosvsEvents));
	hXPosvsYPosvsMeanCharge->Draw("COLZ");
	histSaver->SaveCanvas(cXPosvsYPosvsMean);
		*/
}

void TAnalysisOf3dDiamonds::analyseEvent() {

	if(!eventReader->isValidTrack()) return;
	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	if(!eventReader->isInFiducialCut())	//This is a larger fiducial cut around silicon
		return;

	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster())
		return;
	if(predictedPosition)
		delete predictedPosition;
	predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	Float_t chi2x = predictedPosition->getChi2X();
	Float_t chi2y = predictedPosition->getChi2Y();

	if(chi2x>20||chi2y>20)
		return;
//	cout<< nEvent<<" "<<chi2x<<" "<<chi2y<<endl;

	Float_t xPos = predictedPosition->getPositionX();
	Float_t yPos = predictedPosition->getPositionY();
	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	//Predict Channel hit
	Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
	Float_t positionInDetSystemChannelSpace = settings->convertMetricToChannelSpace(subjectDetector,positionInDetSystemMetric);

	HitandSeedCount(&diamondCluster);
	Int_t clusterSize = diamondCluster.size()-2;
	vecClusterSize.push_back(clusterSize);

	hNumberofClusters->Fill(clusterSize);
	ClusterPlots(eventReader->getNDiamondClusters(),fiducialValueX,fiducialValueY);

	//hFidCutXvsFidCutYvsCharge.at(1)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));

	if(eventReader->getNDiamondClusters()==1){

		//RemoveLumpyClusters(&diamondCluster);

		//Universal PHvsChannel Plot
		for(int i=0; i<settings->diamondPattern.getNIntervals();i++){      //settings->diamondPattern.getNIntervals(); i++){

			hFidCutXvsFidCutYvsCharge.at(i)->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
			hFidCutXvsFidCutYvsEvents.at(i)->Fill(fiducialValueX,fiducialValueY,1);

			pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
			if(diamondCluster.getHighestSignalChannel()<=channels.second&&diamondCluster.getHighestSignalChannel()>=channels.first){

				if(FiducialCut(i+1)==1)	//Cut on fiducial cuts
					return;
				if(RemoveEdgeHits(&diamondCluster,channels)==1)		//If cluster has hit in edge channel remove.
					return;

				//cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;

				/*HitandSeedCount(&diamondCluster);
				if(HitCount>0||SeedCount>1)
					return;
				 */

				hEventsvsChannel.at(i)->Fill(diamondCluster.getHighestSignalChannel());
				hPHvsChannel.at(i)->Fill(diamondCluster.getCharge(false),diamondCluster.getHighestSignalChannel());
				//hPHvsPredictedChannel.at(i)->Fill(diamondCluster.getCharge(false),positionInDetSystemChannelSpace);
				//hPHvsPredictedXPos.at(i)->Fill(diamondCluster.getCharge(false),xPos);
				hLandau.at(i)->Fill(diamondCluster.getCharge(false));
				vecPHDiamondHit.at(i)->push_back(diamondCluster.getCharge(false));
				vecXPredicted.at(i)->push_back(xPos);
				vecYPredicted.at(i)->push_back(yPos);
				hHitandSeedCount.at(i)->Fill(HitCount,SeedCount);
				hChi2XChi2Y.at(i)->Fill(chi2x, chi2y);
				hFidCutXvsFidCutY.at(i)->Fill(fiducialValueX,fiducialValueY);
				//For hFidCutXvsFidCutYvsMeanCharge
				//hFidCutXvsFidCutYvsSeenEvents->Fill(fiducialValueX,fiducialValueY,1);

			}
		}
	}
}

void TAnalysisOf3dDiamonds::initialise3DGridReference() {

	stringstream hGridReferenceName; hGridReferenceName<<""<<FileNameEnd;
	TString nameDet = "hGridRefenrenceDetSpace";
	TString nameCell = "hGridRefenrenceCellSpace";
	Float_t xBins = settings->getNColumns3d();
	Float_t xLow = settings->getXMetalisationStart3d();
	Float_t xHigh =settings->getXMetalisationEnd3d();
	Float_t yBins = settings->getNRows3d();
	Float_t yLow = 0;
	Float_t yHigh = settings->getYMetalisationEnd3d();
	cout<<"nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh"<<endl;
	cout<<nameDet<<" "<<nameDet<<" "<<xBins<<" "<<xLow<<" "<<xHigh<<" "<<yBins<<" "<<yLow<<" "<<yHigh<<endl;
	hGridReferenceDetSpace = new TH2D(nameDet,nameDet,xBins,xLow,xHigh,yBins,yLow,yHigh);
	hGridReferenceCellSpace = new TH2D(nameCell,nameDet,xBins,0,xBins,yBins,0,yBins);

	for(int i=0;i<settings->getNRows3d();i++){
		hGridReferenceDetSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
		hGridReferenceCellSpace->GetXaxis()->SetBinLabel(i+1,TString::Format("%c",(char)('A'+i)));//iLetter.str().c_str());
	}
	for(int j=0;j<settings->getNRows3d();j++){
		hGridReferenceDetSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
		hGridReferenceCellSpace->GetYaxis()->SetBinLabel(j+1,TString::Format("%d",j+1));
	}
	hGridReferenceDetSpace->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTickLength(0.0, "X");
	hGridReferenceDetSpace->SetTickLength(0.0, "Y");
	hGridReferenceCellSpace->SetStats(kFALSE);
	hGridReferenceCellSpace->SetTickLength(0.0, "X");
	hGridReferenceCellSpace->SetTickLength(0.0, "Y");

}

void TAnalysisOf3dDiamonds::initialise3DYAlignmentHistos() {

	//Fiducial Region with Edge Alignment Regions Highlighted
	//hFidCutXvsFidCutYvsChargeYAlignment
	stringstream hFidCutXvsFidCutYvsChargeYAlignmentName; hFidCutXvsFidCutYvsChargeYAlignmentName<<"hFidCutXvsFidCutYvsChargeYAlignment"<<FileNameEnd;
	hFidCutXvsFidCutYvsChargeYAlignment = new TH2D(hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
	hFidCutXvsFidCutYvsChargeYAlignment->GetXaxis()->SetTitle("FidCutX");
	hFidCutXvsFidCutYvsChargeYAlignment->GetYaxis()->SetTitle("FidCutY");
	hFidCutXvsFidCutYvsChargeYAlignment->GetZaxis()->SetTitle("Charge ADC");

	//hFidCutXvsFidCutYvsEventsYAlignment
	hFidCutXvsFidCutYvsEventsYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsEventsYAlignment");

	//hFidCutXvsFidCutYvsMeanChargeYAlignment
	hFidCutXvsFidCutYvsMeanChargeYAlignment = (TH2D*)hFidCutXvsFidCutYvsChargeYAlignment->Clone("hFidCutXvsFidCutYvsMeanChargeYAlignment");

	//hXEdgeCharge
	//xEdge
	stringstream hEdgeChargeName; hEdgeChargeName<<"hEdgeCharge"<<FileNameEnd;
	hEdgeCharge = new TH1F(hEdgeChargeName.str().c_str(),hEdgeChargeName.str().c_str(),250,3100,4100);
	hEdgeCharge->SetTitle(hEdgeChargeName.str().c_str());
	hEdgeCharge->GetXaxis()->SetTitle("X Diamond (um)");
	hEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hEdgeChargeEvents = (TH1F*)hEdgeCharge->Clone("hEdgeChargeEvents");
	hEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

	hEdgeMeanCharge = (TH1F*)hEdgeCharge->Clone("hEdgeMeanCharge");
	hEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	//hYEdgeCharge
	stringstream hyEdgeChargeName; hyEdgeChargeName<<"hyEdgeCharge"<<FileNameEnd;
	hyEdgeCharge = new TH1F(hyEdgeChargeName.str().c_str(),hyEdgeChargeName.str().c_str(),250,5000,6000);
	hyEdgeCharge->SetTitle(hyEdgeChargeName.str().c_str());
	hyEdgeCharge->GetXaxis()->SetTitle("Y Diamond (um)");
	hyEdgeCharge->GetYaxis()->SetTitle("Total Charge [ADC]");

	hyEdgeChargeEvents = (TH1F*)hyEdgeCharge->Clone("hyEdgeChargeEvents");
	hyEdgeChargeEvents->GetYaxis()->SetTitle("Entries");

	hyEdgeMeanCharge = (TH1F*)hyEdgeCharge->Clone("hyEdgeMeanCharge");
	hyEdgeMeanCharge->GetYaxis()->SetTitle("Mean Charge [ADC]");

	//hDeadCellsProfile
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

void TAnalysisOf3dDiamonds::initialise3DOverviewHistos() {

	//hDetXvsDetY3DTotolCharge
	stringstream hDetXvsDetY3DName; hDetXvsDetY3DName<<"hFidCutXvsFidCutYvsChargeYAlignment"<<FileNameEnd;
	hDetXvsDetY3D = new TH2D(hDetXvsDetY3DName.str().c_str(),hDetXvsDetY3DName.str().c_str(),270,settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),330,0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3D->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3D->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3D->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DEvents
	hDetXvsDetY3DvsEvents = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DvsEvents");

	//hDetXvsDetY3DMeanCharge
	hDetXvsDetY3DMeanCharge = (TH2D*)hDetXvsDetY3D->Clone("hDetXvsDetY3DMeanCharge");

	//hDetXvsDetY3DRebinnedChargeTotal
	stringstream hDetXvsDetY3DRebinnedName; hDetXvsDetY3DRebinnedName<<"hFidCutXvsFidCutYvsChargeRebinnedYAlignment"<<FileNameEnd;
	hDetXvsDetY3DRebinned = new TH2D(hDetXvsDetY3DRebinnedName.str().c_str(),hDetXvsDetY3DRebinnedName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3DRebinned->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinned->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinned->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DRebinnedEvents
	hDetXvsDetY3DvsEventsRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DvsEventsRebinned");

	//hDetXvsDetY3DRebinnedMeanCharge
	hDetXvsDetY3DMeanChargeRebinned = (TH2D*)hDetXvsDetY3DRebinned->Clone("hDetXvsDetY3DMeanChargeRebinned");
	//hDetXvsDetY3DMeanChargeRebinned->SetBins(9,2365,3715,11,0,1650,12,0,1200);

	//hDetXvsDetY3DRebinnedMeanChargeRMS
	stringstream hDetXvsDetY3DRebinnedRMSName; hDetXvsDetY3DRebinnedRMSName<<"h3DdetRebinnedRMS"<<FileNameEnd;
	hDetXvsDetY3DRebinnedRMS = new TH2D(hDetXvsDetY3DRebinnedRMSName.str().c_str(),hDetXvsDetY3DRebinnedRMSName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3DRebinnedRMS->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinnedRMS->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinnedRMS->GetZaxis()->SetTitle("Charge ADC");

	//hBinnedMeanCharge
	stringstream hBinnedMeanChargeName; hBinnedMeanChargeName<<"h3DdetCellMeanChargeBinned"<<FileNameEnd;
	hBinnedMeanCharge = new TH1F(hBinnedMeanChargeName.str().c_str(),hBinnedMeanChargeName.str().c_str(),9,400,1300);
	hBinnedMeanCharge->SetTitle(hBinnedMeanChargeName.str().c_str());
	hBinnedMeanCharge->GetXaxis()->SetTitle("MeanCharge");
	hBinnedMeanCharge->GetYaxis()->SetTitle("Entries");

	//hDetXvsDetY3DOverview
	stringstream hDetXvsDetY3DOverviewName; hDetXvsDetY3DOverviewName<<"hDetXvsDetY3DOverview"<<FileNameEnd;
	hDetXvsDetY3DOverview = new TH2D(hDetXvsDetY3DOverviewName.str().c_str(),hDetXvsDetY3DOverviewName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3DOverview->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DOverview->GetYaxis()->SetTitle("Ydet (um)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();

	//hCellNumbering
	stringstream hCellNumberingName; hCellNumberingName<<"h3DdetCellNumbering"<<FileNameEnd;
	hCellNumbering = new TH2D(hCellNumberingName.str().c_str(),hCellNumberingName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hCellNumbering->GetXaxis()->SetTitle("Xdet (um)");
	hCellNumbering->GetYaxis()->SetTitle("Ydet (um)");
	//hDetXvsDetY3DOverview->GetZaxis()->SetTitle();

	//hCellsMeanClusteSize
	stringstream hCellsMeanClusteSizeName; hCellsMeanClusteSizeName<<"hCellsMeanClusteSize"<<FileNameEnd;
	hCellsMeanClusteSize = new TH2D(hCellsMeanClusteSizeName.str().c_str(),hCellsMeanClusteSizeName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());

	///// For Quarter Cell /////

	//hDetXvsDetY3DQuarterCellMeanCharge
	stringstream hDetXvsDetY3DMeanChargeRebinnedQuarterCellName; hDetXvsDetY3DMeanChargeRebinnedQuarterCellName<<"h3DdetQuarterCellMeanCharge"<<FileNameEnd;
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell = new TH2D(hDetXvsDetY3DMeanChargeRebinnedQuarterCellName.str().c_str(),hDetXvsDetY3DMeanChargeRebinnedQuarterCellName.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DMeanChargeRebinnedQuarterCell->GetZaxis()->SetTitle("Charge ADC");

	//hDetXvsDetY3DQuarterCellRMS
	stringstream hDetXvsDetY3DRebinnedQuarterCellRMSName; hDetXvsDetY3DRebinnedQuarterCellRMSName<<"h3DdetQuarterCellRMS"<<FileNameEnd;
	hDetXvsDetY3DRebinnedQuarterCellRMS = new TH2D(hDetXvsDetY3DRebinnedQuarterCellRMSName.str().c_str(),hDetXvsDetY3DRebinnedQuarterCellRMSName.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetXaxis()->SetTitle("Xdet (um)");
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetYaxis()->SetTitle("Ydet (um)");
	hDetXvsDetY3DRebinnedQuarterCellRMS->GetZaxis()->SetTitle("Charge ADC");

	//hQuarterCellsMeanClusteSize
	stringstream hQuarterCellsMeanClusterSizeName; hQuarterCellsMeanClusterSizeName<<"hQuarterCellsMeanClusterSize"<<FileNameEnd;
	hQuarterCellsMeanClusterSize = new TH2D(hQuarterCellsMeanClusterSizeName.str().c_str(),hQuarterCellsMeanClusterSizeName.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d());

	//RebinnedQuarterCellFails
	stringstream RebinnedQuarterCellFailsName; RebinnedQuarterCellFailsName<<"3DdetNumberofQuarterCellFails"<<FileNameEnd;
	RebinnedQuarterCellFails = new TH2D(RebinnedQuarterCellFailsName.str().c_str(),RebinnedQuarterCellFailsName.str().c_str(),settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	RebinnedQuarterCellFails->GetXaxis()->SetTitle("Xdet (um)");
	RebinnedQuarterCellFails->GetYaxis()->SetTitle("Ydet (um)");
	RebinnedQuarterCellFails->GetZaxis()->SetTitle("Quarter Fails");

	//hDetXvsDetY3DQuarterCellGrading
	for(int k=0; k<6; k++){
		stringstream hDetXvsDetY3DMeanChargeQuarterCellGradingName; hDetXvsDetY3DMeanChargeQuarterCellGradingName<<"hDetXvsDetY3DMeanChargeQuarterCellGrading"<<k<<"%%Fail"<<FileNameEnd;
		hDetXvsDetY3DMeanChargeQuarterCellGrading.push_back(new TH2D(hDetXvsDetY3DMeanChargeQuarterCellGradingName.str().c_str(),hDetXvsDetY3DMeanChargeQuarterCellGradingName.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d()));
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetXaxis()->SetTitle("Xdet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetYaxis()->SetTitle("Ydet (um)");
		hDetXvsDetY3DMeanChargeQuarterCellGrading.at(k)->GetZaxis()->SetTitle("Charge ADC");
	}

	//h3DdetQuarterCellFluctuation
	stringstream h3DdetQuarterCellFluctuationName; h3DdetQuarterCellFluctuationName<<"h3DdetQuarterCellFluctuation"<<FileNameEnd;
	h3DdetQuarterCellFluctuation = new TH2D(h3DdetQuarterCellFluctuationName.str().c_str(),h3DdetQuarterCellFluctuationName.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	h3DdetQuarterCellFluctuation->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetQuarterCellFluctuation->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetQuarterCellFluctuation->GetZaxis()->SetTitle("Fluctuation");

	//h3DdetQuarterCellFluctuation1
	stringstream h3DdetQuarterCellFluctuation1Name; h3DdetQuarterCellFluctuation1Name<<"h3DdetQuarterCellFluctuation1"<<FileNameEnd;
	h3DdetQuarterCellFluctuation1 = new TH2D(h3DdetQuarterCellFluctuation1Name.str().c_str(),h3DdetQuarterCellFluctuation1Name.str().c_str(),2*settings->getNColumns3d(),settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),2*settings->getNRows3d(),0,settings->getYMetalisationEnd3d());
	h3DdetQuarterCellFluctuation1->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetQuarterCellFluctuation1->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetQuarterCellFluctuation1->GetZaxis()->SetTitle("Fluctuation");

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

	//Define Landau function for Landau fit.
	Landau = new TF1("Landau","landau(0)",20,80);
	//
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){
			/*
			float xLow = 2365 + i*150;
			float yLow = j*150;
			float xHigh = xLow+150;
			float yHigh = yLow+150;
				*/
			//For Transparent analysis.
			//hCellTransparentLandau
			stringstream hCellTransparentLandauName; hCellTransparentLandauName<<"hCellTransparentLandau"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellTransparentLandau.push_back(new TH1F(hCellTransparentLandauName.str().c_str(),hCellTransparentLandauName.str().c_str(),256,0,2800));
			//hCellsTransparentHitPosition
			stringstream hCellsTransparentHitPositionName; hCellsTransparentHitPositionName<<"hCellsTransparentHitPosition"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsTransparentHitPosition.push_back(new TH2D(hCellsTransparentHitPositionName.str().c_str(),hCellsTransparentHitPositionName.str().c_str(),30,0,150,30,0,150));

			stringstream hCellsChargeName; hCellsChargeName<<"TotalCharge"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsCharge.push_back(new TH2D(hCellsChargeName.str().c_str(),hCellsChargeName.str().c_str(),30,0,150,30,0,150));        //30,0,150,30,0,150);
			stringstream hCellsEventsName; hCellsEventsName<<"Events"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsEvents.push_back(new TH2D(hCellsEventsName.str().c_str(),hCellsEventsName.str().c_str(),30,0,150,30,0,150));			//30,0,150,30,0,150);

			//hCellsEventsCheck
			stringstream hCellsEventsCheckName; hCellsEventsCheckName<<"hCellsEventsCheck"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsEventsCheck.push_back(new TH2D(hCellsEventsCheckName.str().c_str(),hCellsEventsCheckName.str().c_str(),30,0,150,30,0,150));

			stringstream hCellsLandauName; hCellsLandauName<<"hLandauCell"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsLandau.push_back(new TH1F(hCellsLandauName.str().c_str(),hCellsLandauName.str().c_str(),256,0,2800));

			stringstream hCellsClusteSizeName; hCellsClusteSizeName<<"hCellsClusteSize"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsClusteSize.push_back(new TH1F(hCellsClusteSizeName.str().c_str(),hCellsClusteSizeName.str().c_str(),20,0,20));

			//hCellsDeltaX
			stringstream hCellsDeltaXName; hCellsDeltaXName<<"hCellsDeltaX"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsDeltaX.push_back(new TH1F(hCellsDeltaXName.str().c_str(),hCellsDeltaXName.str().c_str(),100,-3,3));

			stringstream hCellsLandauNoColumnName; hCellsLandauNoColumnName<<"hCellsLandauNoColumn"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsLandauNoColumn.push_back(new TH1F(hCellsLandauNoColumnName.str().c_str(),hCellsLandauNoColumnName.str().c_str(),256,0,2800));

			stringstream hCellsEventsNoColumnName; hCellsEventsNoColumnName<<"hCellsEventsNoColumn"<<((i*settings->getNRows3d())+j)<<FileNameEnd;
			hCellsEventsNoColumn.push_back(new TH2D(hCellsEventsNoColumnName.str().c_str(),hCellsEventsNoColumnName.str().c_str(),15,0,150,15,0,150));		//30,0,150,30,0,150);

			for(int k=0;k<4;k++){
				//cout<<"This should not repeat: "<<((i*settings->getNRows3d()*4)+j*4+k)<<endl;
				stringstream hQuaterCellsLandauName; hQuaterCellsLandauName<<"hQuaterCellsLandau"<<((i*settings->getNRows3d()*4)+j*4+k)<<FileNameEnd;
				hQuaterCellsLandau.push_back(new TH1F(hQuaterCellsLandauName.str().c_str(),hQuaterCellsLandauName.str().c_str(),64,0,2800));

				stringstream hQuarterCellsClusterSizeName; hQuarterCellsClusterSizeName<<"hQuarterCellsClusterSize"<<((i*settings->getNRows3d()*4)+j*4+k)<<FileNameEnd;
				hQuarterCellsClusterSize.push_back(new TH1F(hQuarterCellsClusterSizeName.str().c_str(),hQuarterCellsClusterSizeName.str().c_str(),20,0,20));
			}
		}
	}

	//cout<<"The size of hQuarterCells is: "<<hQuaterCellsLandau.size()<<endl;

	//Transparent Analysis
	//hCellsTransparentHitPositionCellGraded
	for(int i=0;i<4;i++){
		stringstream hCellsTransparentHitPositionCellGradedName; hCellsTransparentHitPositionCellGradedName<<"hCellsTransparentHitPositionCellGraded"<<i<<FileNameEnd;
		hCellsTransparentHitPositionCellGraded.push_back(new TH2D(hCellsTransparentHitPositionCellGradedName.str().c_str(),hCellsTransparentHitPositionCellGradedName.str().c_str(),30,0,150,30,0,150));
	}

	//hTransparentCharge3D
	stringstream hTransparentCharge3DName; hTransparentCharge3DName<<"hTransparentCharge3D"<<FileNameEnd;
	hTransparentCharge3D = new TH1F(hTransparentCharge3DName.str().c_str(),hTransparentCharge3DName.str().c_str(),256,0,2800);
	hTransparentCharge3D->SetTitle(hTransparentCharge3DName.str().c_str());
	hTransparentCharge3D->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hTransparentCharge3D->GetYaxis()->SetTitle("Entries");

	//hCellsLandauGraded    &&    hCellsLandauGradedNoColumn
	for(int i=0;i<12;i++){	//Group CellsLandaus within same ranges together. 0-100; 100-200; -> 1100-1200;
		stringstream hCellsLandauGradedName; hCellsLandauGradedName<<"hLandauCellsGraded"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGraded.push_back(new TH1F(hCellsLandauGradedName.str().c_str(),hCellsLandauGradedName.str().c_str(),256,0,2800));
		stringstream hCellsLandauGradedNoColumnName; hCellsLandauGradedNoColumnName<<"hLandauCellsGradedNoColumn%%"<<(i*100)<<" - "<<((i+1)*100)<<FileNameEnd;
		hCellsLandauGradedNoColumn.push_back(new TH1F(hCellsLandauGradedNoColumnName.str().c_str(),hCellsLandauGradedNoColumnName.str().c_str(),256,0,2800));
	}

	//To create Strip detector Landau
	pair<int,int> channels =settings->diamondPattern.getPatternChannels(1);
	//hLandau
	stringstream hLandauName; hLandauName<<"hLandau%%"<<channels.first<<"-"<<channels.second<<"%%"<<FileNameEnd;
	hLandau.push_back(new TH1F(hLandauName.str().c_str(),hLandauName.str().c_str(),256,0,2800));
	hLandau.at(0)->SetTitle(hLandauName.str().c_str());
	hLandau.at(0)->GetXaxis()->SetTitle("PH of diamond cluster");
	hLandau.at(0)->GetYaxis()->SetTitle("number of entries #");

	//hStripFidCutXFidCutYvsCharge
	stringstream hFidCutXvsFidCutYvsChargeYAlignmentName; hFidCutXvsFidCutYvsChargeYAlignmentName<<"hStripFidCutXFidCutYvsCharge"<<FileNameEnd;
	hStripFidCutXFidCutYvsCharge = new TH2D(hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),hFidCutXvsFidCutYvsChargeYAlignmentName.str().c_str(),213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
	hStripFidCutXFidCutYvsCharge->GetXaxis()->SetTitle("FidCutX");
	hStripFidCutXFidCutYvsCharge->GetYaxis()->SetTitle("FidCutY");
	hStripFidCutXFidCutYvsCharge->GetZaxis()->SetTitle("Charge ADC");

	//hStripFidCutXFidCutYvsChargeEvents
	hStripFidCutXFidCutYvsEvents = (TH2D*)hStripFidCutXFidCutYvsCharge->Clone("hStripFidCutXFidCutYvsChargeEvents");

	//hStripFidCutXFidCutYvsMeanCharge
	hStripFidCutXFidCutYvsMeanCharge = (TH2D*)hStripFidCutXFidCutYvsCharge->Clone("hStripFidCutXFidCutYvsMeanCharge");

	//hCellsGoodandBad
	stringstream hCellsHarris18GoodName; hCellsHarris18GoodName<<"hCellsHarris18Good"<<FileNameEnd;
	hCellsHarris18Good = new TH1F(hCellsHarris18GoodName.str().c_str(),hCellsHarris18GoodName.str().c_str(),256,0,2800);
	hCellsHarris18Good->SetTitle(hCellsHarris18GoodName.str().c_str());
	hCellsHarris18Good->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hCellsHarris18Good->GetYaxis()->SetTitle("Entries");

	hCellsHarris10Bad = (TH1F*)hCellsHarris18Good->Clone("hCellsHarris10Bad");

	//h3DdetDeltaXChannel
	stringstream h3DdetDeltaXChannelName; h3DdetDeltaXChannelName<<"h3DdetDeltaXChannel"<<FileNameEnd;
	h3DdetDeltaXChannel = new TH2D(h3DdetDeltaXChannelName.str().c_str(),h3DdetDeltaXChannelName.str().c_str(),270,settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),330,0,settings->getYMetalisationEnd3d());
	h3DdetDeltaXChannel->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetDeltaXChannel->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetDeltaXChannel->GetZaxis()->SetTitle("Delta X (ch)");

	//h3DdetDeltaXChannelAbove1000
	stringstream h3DdetDeltaXChannelAbove1000Name; h3DdetDeltaXChannelAbove1000Name<<"h3DdetDeltaXChannelAbove1000"<<FileNameEnd;
	h3DdetDeltaXChannelAbove1000 = new TH2D(h3DdetDeltaXChannelAbove1000Name.str().c_str(),h3DdetDeltaXChannelAbove1000Name.str().c_str(),270,settings->getXMetalisationStart3d(),settings->getXMetalisationEnd3d(),330,0,settings->getYMetalisationEnd3d());
	h3DdetDeltaXChannelAbove1000->GetXaxis()->SetTitle("Xdet (um)");
	h3DdetDeltaXChannelAbove1000->GetYaxis()->SetTitle("Ydet (um)");
	h3DdetDeltaXChannelAbove1000->GetZaxis()->SetTitle("Delta X (ch)");

}

void TAnalysisOf3dDiamonds::initialise3DCellOverlayHistos() {

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

	//hCellsOverlayedColumnLandau
	stringstream hCellsOverlayedColumnLandauName; hCellsOverlayedColumnLandauName<<"hCellsOverlayedColumnLandau"<<FileNameEnd;
	hCellsOverlayedColumnLandau = new TH1F(hCellsOverlayedColumnLandauName.str().c_str(),hCellsOverlayedColumnLandauName.str().c_str(),256,0,2800);
	hCellsOverlayedColumnLandau->SetTitle(hCellsOverlayedColumnLandauName.str().c_str());
	hCellsOverlayedColumnLandau->GetXaxis()->SetTitle("Mean Charge [ADC]");
	hCellsOverlayedColumnLandau->GetYaxis()->SetTitle("Entries");

	//hCellsOverlayedLandauNoColumn
	stringstream hCellsOverlayedEntriesNoColumnsName; hCellsOverlayedEntriesNoColumnsName<<"hCellsOverlayedEntriesNoColumns"<<FileNameEnd;
	hCellsOverlayedEntriesNoColumns = new TH1F(hCellsOverlayedEntriesNoColumnsName.str().c_str(),hCellsOverlayedEntriesNoColumnsName.str().c_str(),100,0,100);

	stringstream hCellsOverlayedLandauNoColumnName; hCellsOverlayedLandauNoColumnName<<"hCellsOverlayedLandauNoColumn"<<FileNameEnd;
	hCellsOverlayedLandauNoColumn = new TH1F(hCellsOverlayedLandauNoColumnName.str().c_str(),hCellsOverlayedLandauNoColumnName.str().c_str(),256,0,2800);

}

void TAnalysisOf3dDiamonds::initialise3D2DLandauAndClustersizeHistos() {

	//hCellsLandau2D
	stringstream hCellsLandau2DName; hCellsLandau2DName<<"hCellsLandau2D"<<FileNameEnd;
	hCellsLandau2D = new TH2D(hCellsLandau2DName.str().c_str(),hCellsLandau2DName.str().c_str(),256,0,2800,99,0,99);
	hCellsLandau2D->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2D->GetYaxis()->SetTitle("Cell");
	hCellsLandau2DQuarterFail = new TH2D("hCellsLandau2DQuarterFail","hCellsLandau2DQuarterFail",1,0,2800,99,0,99);
	hCellsLandau2DQuarterFail->GetXaxis()->SetTitle("Charge ADC");
	hCellsLandau2DQuarterFail->GetYaxis()->SetTitle("Cell");
	//hCellsLandau2D->SetCanExtend(TH1F::kAllAxes);

	//h2DClusterSize
	stringstream h2DClusterSizeName; h2DClusterSizeName<<"h2DClusterSize"<<FileNameEnd;
	h2DClusterSize = new TH2D(h2DClusterSizeName.str().c_str(),h2DClusterSizeName.str().c_str(),5,0,5,5,0,5);
	h2DClusterSize->GetXaxis()->SetTitle("Failed Quarters");
	h2DClusterSize->GetYaxis()->SetTitle("ClusterSize");
	for(int i=0;i<5;i++){
		stringstream jNumber;
		if(i==5)
			jNumber<<"5+";
		else
			jNumber<<(i+1);
		h2DClusterSize->GetYaxis()->SetBinLabel(i+1,jNumber.str().c_str());
	}
	for(int i=0;i<5;i++){
		stringstream jNumber;
		jNumber<<(i);
		h2DClusterSize->GetXaxis()->SetBinLabel(i+1,jNumber.str().c_str());
	}
	h2DClusterSize->SetContour(99);

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

}

void TAnalysisOf3dDiamonds::saveYAlignmentHistos() {

	//For all Diamond hFidCutXvsFidCutYvsMeanCharge
	cCombinedMeanChargeYAlignment = new TCanvas("cFidCutXvsFidCutYvsMeanChargeYAlignmentNoFidDrawn","cFidCutXvsFidCutYvsMeanChargeYAlignmentNoFidDrawn");
	cCombinedMeanChargeYAlignment->cd();
	*hFidCutXvsFidCutYvsMeanChargeYAlignment = (*hFidCutXvsFidCutYvsChargeYAlignment/(*hFidCutXvsFidCutYvsEventsYAlignment));
	hFidCutXvsFidCutYvsMeanChargeYAlignment->SetEntries(hFidCutXvsFidCutYvsEventsYAlignment->Integral());
	hFidCutXvsFidCutYvsMeanChargeYAlignment->Draw("COLAH");
	histSaver->SaveCanvas(cCombinedMeanChargeYAlignment);

	/*DrawFidCutRegions(); //Draw Fiducial Cut Regions
		 */

	DrawYAlignmentFidCutRegions(); //Draw Fiducial Cut Regions
	cCombinedMeanChargeYAlignment->SetName("hFidCutXvsFidCutYvsMeanCharge");
	histSaver->SaveCanvas(cCombinedMeanChargeYAlignment);

	//For h3DdetMeanCharge
	c3DdetMeanCharge = new TCanvas("c3DdetMeanCharge","c3DdetMeanCharge");
	c3DdetMeanCharge->cd();
	*hDetXvsDetY3DMeanCharge = (*hDetXvsDetY3D/(*hDetXvsDetY3DvsEvents));
	hDetXvsDetY3DMeanCharge->SetEntries(hDetXvsDetY3DvsEvents->Integral());
	//hDetXvsDetY3DMeanCharge->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("c3DdetMeanCharge");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanCharge->Draw("sameCOLZAH");
	//hGridReference->Draw("COL");
	DrawMetallisationGrid(c3DdetMeanCharge);
	histSaver->SaveCanvas(c3DdetMeanCharge);

	//RebinnedMeanCharge
	c3DdetMeanChargeRebinned = new TCanvas("c3DdetMeanChargeRebinned","c3DdetMeanChargeRebinned");
	c3DdetMeanChargeRebinned->cd();
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	*hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DRebinned/(*hDetXvsDetY3DvsEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DvsEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetMeanChargeRebinned");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXTAH");
	//hDetXvsDetY3DvsEventsRebinned->Draw("TEXT");
	DrawMetallisationGrid(c3DdetMeanChargeRebinned);
	histSaver->SaveCanvas(c3DdetMeanChargeRebinned);
	/*c3DdetEventsRebinned = new TCanvas("c3DdetEventsRebinned","c3DdetEventsRebinned");
	//c3DdetEventsRebinned->cd();
	hDetXvsDetY3DvsEventsRebinned->Draw("sameTEXT");			//("COLZ");
	//DrawMetallisationGrid();
	histSaver->SaveCanvas(c3DdetEventsRebinned);
		*/

	//h3DdetDeltaXChannel
	c3DdetDeltaXChannel = new TCanvas("c3DdetDeltaXChannel","c3DdetDeltaXChannel");
	c3DdetDeltaXChannel->cd();
	h3DdetDeltaXChannel->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetDeltaXChannelCanvas");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	h3DdetDeltaXChannel->Draw("sameCOLZAH");
	//h3DdetDeltaXChannelCanvas->Draw("sameTEXTAH");
	DrawMetallisationGrid(c3DdetDeltaXChannel);
	histSaver->SaveCanvas(c3DdetDeltaXChannel);

	//h3DdetDeltaXChannelAbove1000
	c3DdetDeltaXChannelAbove1000 = new TCanvas("c3DdetDeltaXChannelAbove1000","c3DdetDeltaXChannelAbove1000");
	c3DdetDeltaXChannelAbove1000->cd();
	h3DdetDeltaXChannelAbove1000->SetStats(kFALSE);
	h3DdetDeltaXChannelAbove1000->Draw("COLZ");
	DrawMetallisationGrid(c3DdetDeltaXChannelAbove1000);
	histSaver->SaveCanvas(c3DdetDeltaXChannelAbove1000);

	//RebinnedMeanChargeQuarterCell
	int QuaterCellEntrySum = 0;
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
	DrawMetallisationGrid(cDetXvsDetY3DMeanChargeRebinnedQuarterCell);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeRebinnedQuarterCell);

	cout<<"cDetXvsDetY3DRebinnedQuarterCellRMSCanvas"<<endl;
	cDetXvsDetY3DRebinnedQuarterCellRMSCanvas = new TCanvas("c3DdetQuarterCellRMS","c3DdetQuarterCellRMS");
	cDetXvsDetY3DRebinnedQuarterCellRMSCanvas->cd();
	hDetXvsDetY3DRebinnedQuarterCellRMS->SetEntries(QuaterCellEntrySum);
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellRMS");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DRebinnedQuarterCellRMS->Draw("sameCOLZAH");
	hDetXvsDetY3DRebinnedQuarterCellRMS->Draw("sameTEXTAH");
	DrawMetallisationGrid(cDetXvsDetY3DRebinnedQuarterCellRMSCanvas);
	histSaver->SaveCanvas(cDetXvsDetY3DRebinnedQuarterCellRMSCanvas);

	//Cell Overlay

	//cout<<"CellArraySize is: "<<hCellsCharge.size()<<endl;
	//for(int i=0;i<hCellsCharge.size();i++){
	for(int column=0;column<settings->getNColumns3d();column++){		//run over each cell
		for(int row=0;row<settings->getNRows3d();row++){
			//histSaver->SaveHistogram(hCellsLandau.at(i*settings->getNRows3d()+j));	//Save each cells Landau;
			LandauGaussFit landauGauss;
			//TF1* fit = landauGauss.doLandauGaussFit(hCellsDeltaX.at(i*settings->getNRows3d()+j));
			//histSaver->SaveHistogram(hCellsDeltaX.at(i*settings->getNRows3d()+j));		//Save each cells X res in Ch.
			hCellNumbering->SetBinContent((column+1),(row+1),(column*settings->getNRows3d()+row));
			hDetXvsDetY3DRebinnedRMS->SetBinContent((column+1),(row+1),(hCellsLandau.at(column*settings->getNRows3d()+row)->GetRMS()));
			hBinnedMeanCharge->Fill(hDetXvsDetY3DMeanChargeRebinned->GetBinContent(column+1,row+1));	//Fill 1D hist with each cells mean charge

			int Continue = 1;
			for(int badCell=0;badCell<settings->getBadCells3D().size();badCell++){
				//if((i*settings->getNRows3d()+j)==CellsHarris10Bad[k]){
				if((column*settings->getNRows3d()+row)==settings->getBadCells3D().at(badCell)){
					Continue = 0;
					int Entries = hCellsHarris10Bad->GetEntries();;
					hCellsHarris10Bad->Add(hCellsLandau.at(column*settings->getNRows3d()+row));
					hCellsHarris10Bad->SetEntries((Entries+hCellsLandau.at(column*settings->getNRows3d()+row)->GetEntries()));

					RebinnedQuarterCellFails->SetBinContent(column+1,row+1,4);

					for(int l=0;l<4;l++){		//To fill highlighted grading, plot 5.
						cout<<"h2DMeanClusterSizeQuarterCell"<<column<<" "<<row<<" "<<badCell<<" "<<l<<endl;
						//h2DMeanClusterSizeQuarterCell
						Float_t xBin = 4*4+SortArrayPointer[badCell];
						Float_t yBin = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+badCell)->GetMean();
						Float_t zBin = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+badCell)->GetMean();
						h2DMeanClusterSizeQuarterCell->Fill(xBin,yBin,zBin);
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
									Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+k)->GetBinContent(m+1)+EventsSum;
								}
							}
							else{
								Events = hQuarterCellsClusterSize.at(column*settings->getNRows3d()*4+row*4+k)->GetBinContent(l);
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
						int Entries = hCellsHarris18Good->GetEntries();
						hCellsHarris18Good->Add(hCellsLandau.at(column*settings->getNRows3d()+row));
						hCellsHarris18Good->SetEntries((Entries+hCellsLandau.at(column*settings->getNRows3d()+row)->GetEntries()));
					}
				}

				int NumberQuarterFails = 0;
				float QuarterMeanCharge[4];
				float QuarterMeanChargeUnSorted[4];
				for(int k=0;k<4;k++){
					SortArrayPointer[k] = k;
					QuarterMeanCharge[k] = hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+k)->GetMean();
					QuarterMeanChargeUnSorted[k] = hQuaterCellsLandau.at(column*settings->getNRows3d()*4+row*4+k)->GetMean();
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
				/*if(NumberQuarterFails ==0){
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
					*/
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

	/*hCellsOverlayedEvents->Draw("sameTEXT");
	 cCellsOverlayed->SetName("cCellsOverlayedMeanChargeWithEntries");
	 histSaver->SaveCanvas(cCellsOverlayed);
	histSaver->SaveHistogram(hBinnedMeanCharge);	//*/

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
	histSaver->SaveHistogram(hLandau.at(0));

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
	cStripFidCutXFidCutYvsMeanCharge = new TCanvas("cStripFidCutXFidCutYvsMeanCharge","cStripFidCutXFidCutYvsMeanCharge");
	cStripFidCutXFidCutYvsMeanCharge->cd();
	*hStripFidCutXFidCutYvsMeanCharge = *hStripFidCutXFidCutYvsCharge/(*hStripFidCutXFidCutYvsEvents);
	hStripFidCutXFidCutYvsMeanCharge->Draw("COLZ");
	histSaver->SaveCanvas(cStripFidCutXFidCutYvsMeanCharge);

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
		DrawMetallisationGrid(cDetXvsDetY3DMeanChargeQuarterCellGrading.at(k));
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
	DrawMetallisationGrid(cDetXvsDetY3DMeanChargeHighlightedQuarters);
	histSaver->SaveCanvas(cDetXvsDetY3DMeanChargeHighlightedQuarters);

	//h3DdetQuarterCellFluctuation
	c3DdetQuarterCellFluctuation1 = new TCanvas("c3DdetQuarterCellFluctuation","c3DdetQuarterCellFluctuation");
	c3DdetQuarterCellFluctuation1->cd();
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellFluctuation");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	gStyle->SetPaintTextFormat("3.2g");
	h3DdetQuarterCellFluctuation->Draw("sameTEXT");
	DrawMetallisationGrid(c3DdetQuarterCellFluctuation1);
	histSaver->SaveCanvas(c3DdetQuarterCellFluctuation1);


	//h3DdetQuarterCellFluctuation1
	c3DdetQuarterCellFluctuation = new TCanvas("c3DdetQuarterCellFluctuation1","c3DdetQuarterCellFluctuation1");
	c3DdetQuarterCellFluctuation->cd();
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellFluctuation1");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	gStyle->SetPaintTextFormat("3.2g");
	h3DdetQuarterCellFluctuation1->Draw("sameTEXT");
	DrawMetallisationGrid(c3DdetQuarterCellFluctuation);
	histSaver->SaveCanvas(c3DdetQuarterCellFluctuation1);

	gStyle->SetPaintTextFormat("g");

	//hCellNumbering
	cCellNumbering = new TCanvas("cCellNumbering","cCellNumbering");
	cCellNumbering->cd();
	hGridReferenceDetSpace->SetTitle("hCellNumbering");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hCellNumbering->Draw("sameCOLZAH");
	hCellNumbering->Draw("sameTEXTAH");
	DrawMetallisationGrid(cCellNumbering);
	histSaver->SaveCanvas(cCellNumbering);

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
	DrawMetallisationGrid(cDetXvsDetY3DRebinnedRMS);
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
				*/
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
	cCellsLandau2D = new TCanvas("cCellsLandau2D","cCellsLandau2D");
	cCellsLandau2D->SetCanvasSize(1500,3000);
	cCellsLandau2D->cd();
	//hCellsLandau2D->LabelsDeflate();
	hCellsLandau2D->SetEntries(hCellsLandau2DEntries);
	//hCellsLandau2DQuarterFail->Draw("COLZ");
	hCellsLandau2D->SetStats(kFALSE);
	hCellsLandau2D->Draw("COLZ");
	histSaver->SaveCanvas(cCellsLandau2D);
	//hCellsLandau2D->Draw("sameTEXT");

	cCellsLandau2DHighlightedQuarters = new TCanvas("cCellsLandau2DHighlightedQuarters","cCellsLandau2DHighlightedQuarters");
	cCellsLandau2DHighlightedQuarters->SetCanvasSize(1500,3000);
	cCellsLandau2DHighlightedQuarters->cd();
	hCellsLandau2DQuarterFail->SetStats(kFALSE);
	hCellsLandau2DQuarterFail->Draw("COL");
	hCellsLandau2D->Draw("sameCOLZ");
	histSaver->SaveCanvas(cCellsLandau2DHighlightedQuarters);

	//h2DClusterSize
	c2DClusterSizeCanvas = new TCanvas("c2DClusterSizeCanvas","c2DClusterSizeCanvas");
	c2DClusterSizeCanvas->cd();
	h2DClusterSize->SetStats(kFALSE);
	h2DClusterSize->Draw("COLZ");
	*h2DClusterSizeClone1 = *h2DClusterSize/(*h2DClusterSizeClone);
	gStyle->SetPaintTextFormat("3.2g");
	h2DClusterSizeClone1->Draw("sameTEXT");
	gStyle->SetPaintTextFormat("g");
	histSaver->SaveCanvas(c2DClusterSizeCanvas);

	//h2DClusterSizeQuarterCell
	c2DClusterSizeQuarterCell = new TCanvas("c2DClusterSizeQuarterCell","c2DClusterSizeQuarterCell");
	//h2DClusterSizeQuarterCellCanvas->SetCanvasSize(1500,3000);
	c2DClusterSizeQuarterCell->cd();
	//h2DClusterSize->SetEntries(hCellsLandau2DEntries);
	//hCellsLandau2DQuarterFail->Draw("COLZ");
	h2DClusterSizeQuarterCell->SetStats(kFALSE);
	h2DClusterSizeQuarterCell->Draw("COLZ");
	cout<<"h2DClusterSizeQuarterCell"<<endl;
	*h2DClusterSizeQuarterCellClone1 = *h2DClusterSizeQuarterCell/(*h2DClusterSizeQuarterCellClone);
	gStyle->SetPaintTextFormat("3.2g");
	h2DClusterSizeQuarterCellClone1->Draw("sameTEXT");
	gStyle->SetPaintTextFormat("g");
	//h2DClusterSizeXAxis->Draw("sameCOL");
	for(int i=0;i<6;i++){	//Draw lines for Bin edges.
		TLine* BinEdge = new TLine(i*4,0,i*4,-.5);
		BinEdge->SetLineWidth(0.5);
		BinEdge->SetLineColor(kBlack);
		BinEdge->Draw("same");
		if(i<5){
			stringstream Label; Label<<i;
			TText* Text = new TText(h2DClusterSizeXAxis->GetXaxis()->GetBinCenter(i+1),-.5,Label.str().c_str());
			Text->SetTextSize(0.04);
			Text->Draw("same");
		}
	}
	TText* XTitle = new TText((h2DClusterSizeXAxis->GetXaxis()->GetBinCenter(5)-1.75),-.80,"Failed Quarters");
	XTitle->SetTextSize(0.04);
	XTitle->Draw("same");
	cout<<"save c2DClusterSizeQuarterCell"<<endl;
	histSaver->SaveCanvas(c2DClusterSizeQuarterCell);
	//hCellsLandau2D->Draw("sameTEXT");

	//h2DMeanClusterSizeQuarterCell
	c2DMeanClusterSizeQuarterCell = new TCanvas("c2DMeanClusterSizeQuarterCell","c2DMeanClusterSizeQuarterCell");
	c2DMeanClusterSizeQuarterCell->cd();
	h2DMeanClusterSizeQuarterCell->SetStats(kFALSE);
	//*h2DMeanClusterSizeQuarterCell = *h2DMeanClusterSizeQuarterCellTotal/(*h2DMeanClusterSizeQuarterCellEvents);
	h2DMeanClusterSizeQuarterCell->Draw("COLZ");
	gStyle->SetPaintTextFormat("3.2g");
	h2DMeanClusterSizeQuarterCell->Draw("sameTEXT");
	gStyle->SetPaintTextFormat("g");
	for(int i=0;i<6;i++){	//Draw lines for Bin edges.
		TLine* BinEdge = new TLine(i*4,0,i*4,-.5);
		BinEdge->SetLineWidth(0.5);
		BinEdge->SetLineColor(kBlack);
		BinEdge->Draw("same");
		if(i<5){
			stringstream Label; Label<<i;
			TText* Text = new TText(h2DClusterSizeXAxis->GetXaxis()->GetBinCenter(i+1),-.5,Label.str().c_str());
			Text->SetTextSize(0.04);
			Text->Draw("same");
		}
	}
	//TText* XTitle = new TText((h2DClusterSizeXAxis->GetXaxis()->GetBinCenter(5)-1.75),-.80,"Failed Quarters");
	//XTitle->SetTextSize(0.04);
	XTitle->Draw("same");
	histSaver->SaveCanvas(c2DClusterSizeQuarterCell);
	//hCellsLandau2D->Draw("sameTEXT");

	//For c3DdetMeanCharge with Mean ClusterSize.
	c3DdetMeanChargeWithMeanClusterSize = new TCanvas("c3DdetMeanChargeWithMeanClusterSize","c3DdetMeanChargeWithMeanClusterSize");
	c3DdetMeanChargeWithMeanClusterSize->cd();
	*hDetXvsDetY3DMeanChargeRebinned = (*hDetXvsDetY3DRebinned/(*hDetXvsDetY3DvsEventsRebinned));
	hDetXvsDetY3DMeanChargeRebinned->SetEntries(hDetXvsDetY3DvsEventsRebinned->Integral());
	hDetXvsDetY3DMeanChargeRebinned->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetMeanChargeWithMeanClusterSize");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeRebinned->Draw("sameCOLZAH");
	hCellsMeanClusteSize->Draw("sameTEXTAH");
	DrawMetallisationGrid(c3DdetMeanChargeWithMeanClusterSize);
	histSaver->SaveCanvas(c3DdetMeanChargeWithMeanClusterSize);

	//Quarter Cell Mean ClusterSize, with highlighted Fails.
	c3DdetQuarterCellClusterSize = new TCanvas("cQuarterCellsMeanClusteSize","cQuarterCellsMeanClusteSize");
	c3DdetQuarterCellClusterSize->cd();
	hQuarterCellsMeanClusterSize->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetQuarterCellMeanClusterSize");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DMeanChargeQuarterCellGrading.at(5)->Draw("sameCOLAH");
	gStyle->SetPaintTextFormat("3.2g");
	hQuarterCellsMeanClusterSize->Draw("sameTEXT");
	DrawMetallisationGrid(c3DdetQuarterCellClusterSize);
	histSaver->SaveCanvas(c3DdetQuarterCellClusterSize);
	gStyle->SetPaintTextFormat("g");

	//RebinnedQuarterCellFails
	cRebinnedQuarterCellFails = new TCanvas("c3DdetNumberofQuarterCellFails","c3DdetNumberofQuarterCellFails");
	cRebinnedQuarterCellFails->cd();
	RebinnedQuarterCellFails->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetNumberofQuarterCellFails");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	RebinnedQuarterCellFails->Draw("sameCOLZAH");
	RebinnedQuarterCellFails->Draw("sameTEXTAH");
	DrawMetallisationGrid(cRebinnedQuarterCellFails);
	histSaver->SaveCanvas(cRebinnedQuarterCellFails);

	//hDetXvsDetY3DOverview
	hOverview = new TCanvas("c3DdetOverview","c3DdetOverview");
	hOverview->cd();
	hDetXvsDetY3DOverview->SetStats(kFALSE);
	hGridReferenceDetSpace->SetTitle("h3DdetOverview");		//Set title to require
	hGridReferenceDetSpace->Draw("COL");
	hDetXvsDetY3DOverview->Draw("sameCOLZAH");
	histSaver->SaveCanvas(hOverview);

	//hDeadCell
	//histSaver->SaveHistogram(hDeadCell);
	//histSaver->SaveHistogram(hDeadCellEvents);

	cDeadCellMeanCharge = new TCanvas("cDeadCellMeanCharge","cDeadCellMeanCharge");
	cDeadCellMeanCharge->cd();
	*hDeadCellMeanCharge = (*hDeadCell/(*hDeadCellEvents));
	hDeadCellMeanCharge->SetEntries(hDeadCellEvents->Integral());
	hDeadCellMeanCharge->Draw();
	Float_t ymax1 = hDeadCellMeanCharge->GetMaximum();
	TLine* CellEdge1 = new TLine(750,0,750,ymax1);
	TLine* CellEdge2 = new TLine(900,0,900,ymax1);
	CellEdge1->SetLineWidth(2);		CellEdge2->SetLineWidth(2);
	CellEdge1->SetLineColor(kRed);	CellEdge2->SetLineColor(kRed);
	CellEdge1->Draw("same");		CellEdge2->Draw("same");
	histSaver->SaveCanvas(cDeadCellMeanCharge);

	//hxEdgeCharge
	//histSaver->SaveHistogram(hEdgeCharge);
	//histSaver->SaveHistogram(hEdgeChargeEvents);

	cEdgeMeanCharge = new TCanvas("cEdgeMeanCharge","cEdgeMeanCharge");
	cEdgeMeanCharge->cd();
	*hEdgeMeanCharge = (*hEdgeCharge/(*hEdgeChargeEvents));
	hEdgeMeanCharge->SetEntries(hEdgeChargeEvents->Integral());
	hEdgeMeanCharge->Draw();
	Float_t ymax = hEdgeMeanCharge->GetMaximum();
	TLine* DetectorEdge = new TLine(3715,0,3715,ymax);
	DetectorEdge->SetLineWidth(2);
	DetectorEdge->SetLineColor(kRed);
	DetectorEdge->Draw("same");
	histSaver->SaveCanvas(cEdgeMeanCharge);
	//histSaver->SaveHistogram(hEdgeMeanCharge);

	//hyEdgeCharge
	//histSaver->SaveHistogram(hyEdgeCharge);
	//histSaver->SaveHistogram(hyEdgeChargeEvents);

	cyEdgeMeanCharge = new TCanvas("cyEdgeMeanCharge","cyEdgeMeanCharge");
	cyEdgeMeanCharge->cd();
	*hyEdgeMeanCharge = (*hyEdgeCharge/(*hyEdgeChargeEvents));
	hyEdgeMeanCharge->SetEntries(hyEdgeChargeEvents->Integral());
	hyEdgeMeanCharge->Draw();
	//Float_t ymax = hyEdgeMeanCharge->GetMaximum();
	//TLine* DetectorEdge = new TLine(3715,0,3715,ymax);
	//DetectorEdge->SetLineWidth(2);
	//DetectorEdge->SetLineColor(kRed);
	//DetectorEdge->Draw("same");
	histSaver->SaveCanvas(cyEdgeMeanCharge);

}
void TAnalysisOf3dDiamonds::YAlignment() {

	if(!eventReader->isValidTrack())
		return;
	vector<UInt_t> vecSilPlanes;

	for(UInt_t pl=0;pl<TPlaneProperties::getNSiliconPlanes();pl++){vecSilPlanes.push_back(pl);}
	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();

	if(predictedPosition)
		delete predictedPosition;
	predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	Float_t chi2x = predictedPosition->getChi2X();
	Float_t chi2y = predictedPosition->getChi2Y();
	Float_t maxChi2 = settings->getChi2Cut3D();
	if(chi2x>100||chi2y>100)     //(chi2x>maxChi2||chi2y>maxChi2)
		return;

	Float_t xPos = predictedPosition->getPositionX();	//Predicted positions in labframe
	Float_t yPos = predictedPosition->getPositionY();
	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	Float_t Xdet = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos);
	//cout<<"YOffset: "<<settings->get3DYOffset()<<endl;
	Float_t Ydet = GetYPositionInDetSystem()-settings->get3DYOffset();	//3890 is the calculated offset from comparing x and y edge charge profiles.
	if(!eventReader->isInFiducialCut())	//This is a larger fiducial cut around silicon
		return;
	vecChi2X.push_back(chi2x);
	vecChi2Y.push_back(chi2y);
//	vecXPredicted.push_back(xPos);
//	vecYPredicted.push_back(yPos);
//	vecXPreditedDetector.push_back(Xdet);
//	vecYPreditedDetector.push_back(Ydet);

	//For Transparent Analysis.
	float hEntries0 = 0;
	float TransparentCharge = 0;
	float TransparentChargeAddition = 0;
	int SaturatedEvent = 0;
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){
			hEntries0 = hCellsEventsCheck.at(i*11+j)->Integral();
			float xminus = 2365+i*150; //+5;		//2365 is the start of the 3D detector in x
			float yminus = j*150;
			hCellsEventsCheck.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
			int Xcell = (Xdet-xminus);
			Int_t HitCell;
			//Int_t* HitChannels;
			if(hEntries0 != hCellsEventsCheck.at(i*11+j)->Integral()){	//If entries in this histogram range have increased fill corresponding hCellLandau
				HitCell = (i*11+j);
				//				cout<<"Hit cell: "<<HitCell<<endl;
				//HitChannels = CellToChannel(HitCell);
				//				cout<<CellToChannel(HitCell,Xcell)[0]<<endl;
				for(int k=1;k<CellToChannel(HitCell,Xcell)[0]+1;k++){	//First memeber of array is the number of channels to readout.
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

	//eventReader->getAdcValue(subjectDetector,CellToChannel());

	//Int_t getAdcValue(UInt_t det,UInt_t ch);

	//CellToChannel();

	if(eventReader->getNDiamondClusters()!=1)
		return;

	TCluster diamondCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond());
	if(diamondCluster.isSaturatedCluster())
		return;

	//HitandSeedCount(&diamondCluster, 0);

	//To create Strip detector Spectrum
	if(RemoveEdgeClusters(&diamondCluster, 1)==0){
		//if(fiducialValueX<101&&fiducialValueX>94&&fiducialValueY<102&&fiducialValueY>72){
		if(fiducialValueX<104&&fiducialValueX>93&&fiducialValueY<96&&fiducialValueY>70){
			hLandau.at(0)->Fill(diamondCluster.getCharge(false));
			hStripFidCutXFidCutYvsCharge->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
			hStripFidCutXFidCutYvsEvents->Fill(fiducialValueX,fiducialValueY,1);
			//cout<<"Made it into loop"<<endl;
		}
	}

	//if(HitCount==1&&SeedCount==1){
	//if(SeedCount<3){
	pair<int,int> cell = getCellNo(Xdet,Ydet);
	Int_t cellNo = cell.first;
	Int_t quarterNo = cell.second;
	Int_t row = cellNo% settings->getNColumns3d();
	Int_t column = cellNo / settings->getNColumns3d();
//	if(cellNo>=0)cout<<"cell: "<<cellNo<<"--> row "<<row<<", column "<<column<<endl;
	Float_t cellWidth = 150;
	Float_t cellHight = 150;
	Float_t startOf3dDetectorX = 2365;
	Float_t xminus = startOf3dDetectorX+column*cellWidth; //+5;		//2365 is the start of the 3D detector in x
	Float_t yminus = row*cellHight;

	Float_t relPosX =Xdet - xminus;
	Float_t relPosY = Ydet - yminus;
	Int_t area3DwithColumns = 2;
	if (!settings->isClusterInDiaDetectorArea(diamondCluster,area3DwithColumns)){
		return;
	}
//	if(cellNo<0)
//		return;

	//clusteredAnalysis->addEvent(xPos,yPos,cellNo,quarterNo,relPosX,relPosY,diamondCluster);
	//ClusterShape(&diamondCluster);

	/*if(YAlignmentFiducialCut())
			return;
	 */
	hFidCutXvsFidCutYvsChargeYAlignment->Fill(fiducialValueX,fiducialValueY,diamondCluster.getCharge(false));
	hFidCutXvsFidCutYvsEventsYAlignment->Fill(fiducialValueX,fiducialValueY,1);

	//xEdge
	if(!YAlignmentFiducialCut(settings->getxEdgeFicucialRegion())){
		Float_t positionInDetSystemMetric = eventReader->getPositionInDetSystem(subjectDetector, xPos, yPos); //This takes x and y predicted in lab frame and gives corresponding position in detector frame.
		hEdgeChargeEvents->Fill(positionInDetSystemMetric);
		hEdgeCharge->Fill(positionInDetSystemMetric, diamondCluster.getCharge(false));
	}
	//yEdge
	if(!YAlignmentFiducialCut(settings->getyEdgeFicucialRegion())){
		hyEdgeChargeEvents->Fill(GetYPositionInDetSystem());
		hyEdgeCharge->Fill(GetYPositionInDetSystem(), diamondCluster.getCharge(false));
	}

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
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){
			hEntries = hCellsEvents.at(i*11+j)->Integral();
			float xminus = 2365+i*150; //+5;		//2365 is the start of the 3D detector in x
			float yminus = j*150;
			hCellsCharge.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),diamondCluster.getCharge(false));
			hCellsEvents.at((i*11+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
			if(hEntries != hCellsEvents.at(i*11+j)->Integral()){	//If entries in this histogram range have increased fill corresponding hCellLandau
//				cout<< i<<" "<<j<<" "<<i*11+j<<" ?= "<< cellNo<<" "<<xPos<<" "<<yPos<<endl;
				if(cellNo != i*11+j)
					cerr<<"something is wrong "<<i<<" "<<j<<endl;
				hCellsLandau.at(i*11+j)->Fill(diamondCluster.getCharge(false));
				hCellsClusteSize.at(i*11+j)->Fill((diamondCluster.getClusterSize()-2));		//Cluster seems to be 2 smaller than ClusterSize for some reason?
				hCellsDeltaX.at(i*11+j)->Fill((diamondCluster.getHighestSignalChannel()-XdetChannelSpace));
				hCellsColumnCheck55->Fill((Xdet-xminus),(Ydet-yminus),1);	//Fill Column Check histo 5*5
				hCellsColumnCheck1010->Fill((Xdet-xminus),(Ydet-yminus),1);	//Fill Column Check histo 10*10
				//To fill quarter cell histograms
				if((Xdet-xminus)>0&&(Xdet-xminus)<75  &&  (Ydet-yminus)>0&&(Ydet-yminus)<75){	//bottom left
					hQuaterCellsLandau.at(i*11*4+j*4)->Fill(diamondCluster.getCharge(false));
					hQuarterCellsClusterSize.at(i*11*4+j*4)->Fill((diamondCluster.getClusterSize()-2));
				}
				if((Xdet-xminus)>0&&(Xdet-xminus)<75  &&  (Ydet-yminus)>75&&(Ydet-yminus)<150){	//top left
					hQuaterCellsLandau.at(i*11*4+j*4+1)->Fill(diamondCluster.getCharge(false));
					hQuarterCellsClusterSize.at(i*11*4+j*4+1)->Fill((diamondCluster.getClusterSize()-2));
				}
				if((Xdet-xminus)>75&&(Xdet-xminus)<150  &&  (Ydet-yminus)>0&&(Ydet-yminus)<75){	//bottom right
					hQuaterCellsLandau.at(i*11*4+j*4+2)->Fill(diamondCluster.getCharge(false));
					hQuarterCellsClusterSize.at(i*11*4+j*4+2)->Fill((diamondCluster.getClusterSize()-2));
				}
				if((Xdet-xminus)>75&&(Xdet-xminus)<150  &&  (Ydet-yminus)>75&&(Ydet-yminus)<150){	//top right
					hQuaterCellsLandau.at(i*11*4+j*4+3)->Fill(diamondCluster.getCharge(false));
					hQuarterCellsClusterSize.at(i*11*4+j*4+3)->Fill((diamondCluster.getClusterSize()-2));
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
									hCellsLandauNoColumn.at(i*settings->getNRows3d()+j)->Fill(diamondCluster.getCharge(false));
									hCellsEventsNoColumn.at((i*settings->getNRows3d()+j))->Fill((Xdet-xminus),(Ydet-yminus),1);
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

}

int TAnalysisOf3dDiamonds::YAlignmentFiducialCut(vector<int> nFiducialRegion){

	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	int Throw=0;
	if(fiducialValueX<nFiducialRegion.at(1)||fiducialValueX>nFiducialRegion.at(0))
		Throw =1;
	if(fiducialValueY<nFiducialRegion.at(3)||fiducialValueY>nFiducialRegion.at(2))
		Throw =1;
	return Throw;
}


void TAnalysisOf3dDiamonds::DrawYAlignmentFidCutRegions() {

	cCombinedMeanChargeYAlignment->cd();
	//Draw xEdge Region
	FidCutChannelYAlignmentTBox.push_back(new TBox(settings->getxEdgeFicucialRegion().at(1), settings->getxEdgeFicucialRegion().at(3), settings->getxEdgeFicucialRegion().at(0), settings->getxEdgeFicucialRegion().at(2)));
	//Draw yEdge Region
	FidCutChannelYAlignmentTBox.push_back(new TBox(settings->getyEdgeFicucialRegion().at(1), settings->getyEdgeFicucialRegion().at(3), settings->getyEdgeFicucialRegion().at(0), settings->getyEdgeFicucialRegion().at(2)));
	for(int i=0;i<2;i++){
		FidCutChannelYAlignmentTBox.at(i)->SetFillStyle(0);
		FidCutChannelYAlignmentTBox.at(i)->SetLineWidth(2);
		FidCutChannelYAlignmentTBox.at(i)->SetLineColor(kRed);
		FidCutChannelYAlignmentTBox.at(i)->Draw("same");
	}
}

Int_t* TAnalysisOf3dDiamonds::CellToChannel(int nCell, float nXcell) {
	Int_t* ChannelsToRead;
	Int_t Channel = -9999;
	Int_t nCellsPerColumn = settings->getNRows3d();
	Int_t nColumn = nCell / nCellsPerColumn;
	int pattern3d = settings->get3dWithHolesDiamondPattern();
	pair<int,int> channels = settings->diamondPattern.getPatternChannels(pattern3d);
//	cout<<"\tPattern3D "<<pattern3d<<": channels "<<channels.first<<"-"<<channels.second<<endl;
	Int_t startChannel = channels.first;//85;
//	cout<<"First 3DwH channel: "<<settings->diamondPattern.getPatternChannels(settings->get3dWithHolesDiamondPattern()).first<<endl;
	/// This part of code is not working properly.
//	cout<<"nColumn hit is: "<<nColumn<<endl;
	Channel = startChannel + nColumn;

//	int BadCells[] = {1,6,7,8,13,27,45,53,77,90,97};
	vector<int> BadCells = settings->getBadCells3D();
	int BadCell = 0;
	for(int i=0; i<BadCells.size(); i++){
		if(nCell ==  BadCells[i]){
			BadCell = 1;
			if(nCell<settings->getNRows3d() && nCell>=0){	//If cell on left edge
				Int_t A[] = {2,Channel,(Channel+1)};
				ChannelsToRead = A;
			}
			if(nCell<9*settings->getNRows3d() && nCell>=8*settings->getNRows3d()){	//If cell on right edge
				Int_t A[] = {2,(Channel-1),Channel};
				ChannelsToRead = A;
			}
			if(nCell<8*settings->getNRows3d() && nCell>=settings->getNRows3d()){	//If cell in the central region
				Int_t A[] = {3,(Channel-1),Channel,(Channel+1)};
				ChannelsToRead = A;
			}
		}
	}
	if(BadCell == 0){		//To account for charge sharing in edge 50um
		if(nXcell<5){	//If Xpos in cell left charge share
			if(nCell<settings->getNRows3d() && nCell>=0){	//If on leftmost column of detector
				Int_t A[] = {1,Channel};
				ChannelsToRead = A;
			}
			else{
				Int_t A[] = {2,(Channel-1),Channel};
				ChannelsToRead = A;
			}
		}
		if(nXcell>145){	//If Xpos in cell right charge share
			if(nCell<9*settings->getNRows3d() && nCell>=8*settings->getNRows3d()){	//If on rightmost column of detector
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
	if (!nCanvas)
		return;
	nCanvas->cd();
	//vector<TBox*> Grid;
	TCutG* gridPoint;
	for(int i=0;i<settings->getNColumns3d();i++){
		for(int j=0;j<settings->getNRows3d();j++){
			float xLow = settings->getXMetalisationStart3d() + i*150;
			float yLow = j*150;
			float xHigh = xLow+150;
			float yHigh = yLow+150;
			TString name = nCanvas->GetName();
			name.Append(TString::Format("_CellGrid%d_%d",i,j));
			TCutG * gridPoint = new TCutG(name,5);
			gridPoint->SetPoint(0,xLow,yLow);
			gridPoint->SetPoint(1,xLow,yHigh);
			gridPoint->SetPoint(2,xHigh,yHigh);
			gridPoint->SetPoint(3,xHigh,yLow);
			gridPoint->SetPoint(4,xLow,yLow);
			gridPoint->SetFillStyle(0);
			gridPoint->SetLineWidth(1);
			gridPoint->SetLineColor(kBlack);
			gridPoint->Draw("same");
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

int TAnalysisOf3dDiamonds::FiducialCut(int Detector){

	float fiducialValueX= eventReader->getFiducialValueX();
	float fiducialValueY= eventReader->getFiducialValueY();
	int Throw=0;
	//cout<<"FidX: "<<fiducialValueX<<endl;
	//cout<<"FidY: "<<fiducialValueY<<endl;
	//cout<<"X Less than: "<<FidCut.at(Detector)->GetXHigh()<<" Greater than: "<<FidCut.at(Detector)->GetXLow()<<endl;
	//cout<<"Y Less than: "<<FidCut.at(Detector)->GetYHigh()<<" Greater than: "<<FidCut.at(Detector)->GetYLow()<<endl;
	if(FiducialChannel){
		if(fiducialValueX<settings->get3dFidCuts()->getFidCut(Detector)->GetXLow()||fiducialValueX>settings->get3dFidCuts()->getFidCut(Detector)->GetXHigh())
			Throw =1;
		if(fiducialValueY<settings->get3dFidCuts()->getFidCut(Detector)->GetYLow()||fiducialValueY>settings->get3dFidCuts()->getFidCut(Detector)->GetYHigh())
			Throw =1;
	}
	//if(Throw==1)
		//cout<<"Detector: "<<Detector<<"Thrown: "<<fiducialValueX<<",    "<<fiducialValueY<<endl;
	return Throw;
}

void TAnalysisOf3dDiamonds::DrawFidCutRegions() {

	hCombinedMeanCharge->cd();
	for(int i=0;i<settings->get3dFidCuts()->getNFidCuts();i++){
		FidCutChannelTBox.push_back(new TBox(settings->get3dFidCuts()->getFidCut(i+1)->GetXLow(), settings->get3dFidCuts()->getFidCut(i+1)->GetYLow(), settings->get3dFidCuts()->getFidCut(i+1)->GetXHigh(), settings->get3dFidCuts()->getFidCut(i+1)->GetYHigh()));
		FidCutChannelTBox.at(i)->SetFillStyle(0);
		FidCutChannelTBox.at(i)->SetLineWidth(2);
		FidCutChannelTBox.at(i)->SetLineColor(kRed);
		FidCutChannelTBox.at(i)->Draw("same");
	}

/*
	FidCudBoundMetric = (TH2F*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("FidCudBoundMetric");
	FidCudBoundMetric->SetFillStyle(3444);
	//FidCudBoundMetric->SetBins(80,90,170,60,60,120);
	cFidCudBoundChannel = new TCanvas("cFidCudBoundMetric","cFidCudBoundMetric");
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
	cFidCudBoundChannel->cd();
	FidCudBoundMetric->SetTitle("");
	FidCudBoundMetric->GetXaxis()->SetTitle("");
	FidCudBoundMetric->GetYaxis()->SetTitle("");
	FidCudBoundMetric->SetStats(kFALSE);
	FidCudBoundMetric->Draw();
	b.Draw();
	histSaver->SaveCanvas(cFidCudBoundChannel);
	//*/
}

void TAnalysisOf3dDiamonds::PredictChannelHit(TCluster* nCluster) {
	for (UInt_t clPos=0;clPos < nCluster->getClusterSize();clPos++){
		if(nCluster->isHit(clPos)&&!nCluster->isSeed(clPos)) cout<<" HIT ";
		if(nCluster->isSeed(clPos)) cout<<" SEED ";
	}
	cout<<endl;
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
 * 			+     0     +     1     +       ^
 * 			+           +           +     Y |
 * 			+-----------+-----------+       |
 * 			+           +           +       |
 * 			+     2     +     3     +       |
 * 			+           +           +       |
 * 			+-----------+-----------+       +-------->
 * 			                                       X
 * @param xDet
 * @param yDet
 * @return
 */

pair<int,int> TAnalysisOf3dDiamonds::getCellNo(Float_t xDet, Float_t yDet) {
	// i column
	// j row

	Float_t startOf3dDetectorX = 2365;
	Float_t startOf3dDetectorY = 0;
	Float_t cellWidth =150;
	Float_t cellHight = 150;
	Int_t column = (xDet-startOf3dDetectorX)/cellWidth;
	Int_t row = (yDet - startOf3dDetectorY)/cellHight;

	Float_t xminus = startOf3dDetectorX+column*cellWidth; //+5;		//2365 is the start of the 3D detector in x
	Float_t yminus = row*cellHight;
	Float_t deltaX = xDet - xminus;
	Float_t deltaY = yDet - yminus;
	Int_t cell = -1;
	if(column >= 0 && column < settings->getNColumns3d() && row >= 0 && row < settings->getNRows3d())
		cell = column * settings->getNRows3d()+row;
	Int_t quarter = -1;
	if (deltaY>=0&&deltaX>=0&&deltaX<=cellWidth&&deltaY<=cellHight){
		int quarterX = deltaX/(cellWidth/2);
		int quarterY = deltaY/(cellHight/2);
		quarter = quarterX +quarterY*2;
		if(verbosity>4)cout << "\t"<<deltaX <<"/"<<deltaY << " "<<quarterX<< "/"<<quarterY<<" "<<cellWidth<<"/"<<cellHight<<endl;
	}
	if(verbosity>4 && column >= 0 && row >= 0){
		cout<<" TAnalysisOf3dDiamonds::getCellNo " << xDet <<"/"<<yDet<<endl;
		cout << "\tcolumn: " << column << ", row: " << row << endl;
		cout << "\tdeltaX: " << deltaX << ", deltaY: " << deltaY <<endl;
		cout<<"\t cell: "<< cell << ", quarter: " << quarter <<endl;
	}
//	i*11+j
	return make_pair(cell,quarter);
}
