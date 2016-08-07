/*
 * TAnalysisOfSelection.cpp
 *
 *  Created on: May 18, 2012
 *      Author: bachmair
 */

#include "../include/TAnalysisOfSelection.hh"

TAnalysisOfSelection::TAnalysisOfSelection(TSettings *newSettings) {
	if(newSettings!=0)
		this->settings=newSettings;
	else exit(0);
    cout<<"\n\n\n**********************************************************"<<endl;
    cout<<"*****************Analysis of Selection *******************"<<endl;
    cout<<"**********************************************************"<<endl;
	verbosity = settings->getVerbosity();
	UInt_t runNumber=settings->getRunNumber();
	htmlLandau=new THTMLLandaus(settings);
	htmlSelection = new THTMLSelectionAnalysis(settings);
	settings->goToSelectionTreeDir();
	eventReader=new TADCEventReader(settings->getSelectionTreeFilePath(),settings);
	histSaver=new HistogrammSaver(settings);
	settings->goToSelectionAnalysisDir();
	//	htmlPedestal->setSubdirPath("selectionAnalysis");
	histSaver->SetPlotsPath(settings->getSelectionAnalysisPath());
	histSaver->SetRunNumber(runNumber);
	htmlLandau->setFileGeneratingPath(settings->getSelectionAnalysisPath());
	htmlSelection->setFileGeneratingPath(settings->getSelectionAnalysisPath());
	settings->goToSelectionTreeDir();
	initialiseHistos();

	if (verbosity) cout<<"end initialise"<<endl;
    xDivisions = 3;
    yDivisions = 3;

    if (verbosity)
        settings->diamondPattern.showPatterns();
}

TAnalysisOfSelection::~TAnalysisOfSelection() {
	htmlLandau->generateHTMLFile();
	htmlSelection->addSelectionPlots();
	htmlSelection->addAreaPlots();
	htmlSelection->addFiducialCutPlots();

	htmlSelection->generateHTMLFile();

	if(eventReader!=0) delete eventReader;
	if(histSaver!=0)   delete histSaver;
	if(htmlLandau!=0)  delete htmlLandau;
	if(htmlSelection!=0) delete htmlSelection;
	settings->goToOutputDir();
}

void TAnalysisOfSelection::doAnalysis(UInt_t nEvents)
{
    initPHvsEventNoAreaPlots(0,nEvents);
	cout<<"analyze selection data..."<<endl;
	if(nEvents<=0) nEvents=eventReader->GetEntries();
	histSaver->SetNumberOfEvents(nEvents);
	for(nEvent=0;nEvent<nEvents;nEvent++){
	    TRawEventSaver::showStatusBar(nEvent,nEvents,1000);
		bool skip = false;
		for(Int_t i=0; i< settings->getSkipEvents().size(); i++) {
			if ((settings->getSkipEvents().at(i).first < nEvent) && (nEvent < settings->getSkipEvents().at(i).second)){
				skip = true;
				break;
			}
		}
		if(!skip) {
			eventReader->LoadEvent(nEvent);
			analyseEvent();
		}
	}
	saveHistos();

}

void TAnalysisOfSelection::initialiseHistos()
{
	histoLandauDistribution = new TH2F("hLandauDiamond_OneCluster","hLandauDiamond_OneCluster",512,0,4096,8,0.5,8.5);
	histoLandauDistribution->GetXaxis()->SetTitle("Charge in ADC counts");
	histoLandauDistribution->GetYaxis()->SetTitle("ClusterSize");

	histoLandauDistribution2D = new TH2F("histoLandauDistribution2D_Clustersize_1_2","histoLandauDistribution2D_Clustersize_1_2",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2D->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2D->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2D->GetZaxis()->SetTitle("number of entries");

	histoLandauDistribution2DNoBorderSeed = new TH2F("histoLandauDist2DNoBorderSeed","histoLandauDist2DNoBorderSeed",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2DNoBorderSeed->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2DNoBorderSeed->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2DNoBorderSeed->GetZaxis()->SetTitle("number of entries");

	histoLandauDistribution2DNoBorderHit = new TH2F("histoLandauDist2D_Clustersize_1_2_noBorderHit","histoLandauDist2D_Clustersize_1_2_noBorderHit",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2DNoBorderHit->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2DNoBorderHit->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2DNoBorderHit->GetZaxis()->SetTitle("number of entries");

	histoLandauDistribution2D_unmasked = new TH2F("histoLandauDistribution2D_Clustersize_1_2_unmasked","histoLandauDistribution2D_Clustersize_1_2_unmasked",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2D_unmasked->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2D_unmasked->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2D_unmasked->GetZaxis()->SetTitle("number of entries");

	histoLandauDistribution2DNoBorderSeed_unmasked = new TH2F("hLandauDist2D_Clustersize_1_2NoBorderSeed-unmasked","hLandauDist2D_Clustersize_1_2NoBorderSeed",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2DNoBorderSeed_unmasked->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2DNoBorderSeed_unmasked->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2DNoBorderSeed_unmasked->GetZaxis()->SetTitle("number of entries");

	histoLandauDistribution2DNoBorderHit_unmasked = new TH2F("hLandauDist2D_Clustersize_1_2NoBorderHit-unmasked","hLandauDist2D_Clustersize_1_2NoBorderHit",512,0,4096,TPlaneProperties::getNChannelsDiamond(),0,TPlaneProperties::getNChannelsDiamond()-1);
	histoLandauDistribution2DNoBorderHit_unmasked->GetXaxis()->SetTitle("Charge of Cluster in ADC counts");
	histoLandauDistribution2DNoBorderHit_unmasked->GetYaxis()->SetTitle("channel of highest Signal");
	histoLandauDistribution2DNoBorderHit_unmasked->GetZaxis()->SetTitle("number of entries");

	hValidSiliconAndDiamondHit = new TH2F("hValidSiliconAndDiamondHit","hValidSiliconAndDiamondHit",256*4,0,255,256*4,0,255);
	hValidSiliconAndDiamondHit->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndDiamondHit->GetYaxis()->SetTitle("FidCutValue in Y");

	hValidSiliconAndOneDiamondHit = new TH2F("hValidSiliconAndOneDiamondHit","hValidSiliconAndOneDiamondHit",256*4,0,255,256*4,0,255);
	hValidSiliconAndOneDiamondHit->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndOneDiamondHit->GetYaxis()->SetTitle("FidCutValue in Y");

	hValidSiliconAndOneDiamondHitNotMasked  = new TH2F("hValidSiliconAndOneDiamondHitNotMasked","hValidSiliconAndOneDiamondHitNotMasked	",256*4,0,255,256*4,0,255);
	hValidSiliconAndOneDiamondHitNotMasked->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndOneDiamondHitNotMasked->GetYaxis()->SetTitle("FidCutValue in Y");

	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels= new TH2F("hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels","hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels",256*4,0,255,256*4,0,255);
	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels->GetYaxis()->SetTitle("FidCutValue in Y");

	hValidSiliconAndOneDiamondHitInOneArea  = new TH2F("hValidSiliconAndOneDiamondHitInOneArea","hValidSiliconAndOneDiamondHitInOneArea	",256*4,0,255,256*4,0,255);
	hValidSiliconAndOneDiamondHitInOneArea->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndOneDiamondHitInOneArea->GetYaxis()->SetTitle("FidCutValue in Y");

	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut = new TH2F("hValidSiliconAndOneDiamondHitInSameAreaAndFidCut","hValidSiliconAndOneDiamondHitInSameAreaAndFidCut",256*4,0,255,256*4,0,255);
	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->GetXaxis()->SetTitle("FidCutValue in X");
	hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->GetYaxis()->SetTitle("FidCutValue in Y");

	hOneClusterHitChannel = new TH1F("hOneClusterHitChannel","hOneClusterHitChannel",128,0,128);
	hOneClusterHitChannel->GetXaxis()->SetTitle("Diamond Cluster Channel no");
	hOneClusterHitChannel->GetYaxis()->SetTitle("Number of entries");

	hOneClusterHitChannelVsArea = new TH2F("hOneClusterHitChannelVsArea","hOneClusterHitChannelVsArea",128,0,128,
	        settings->getNDiaDetectorAreas(),0,settings->getNDiaDetectorAreas());
	hOneClusterHitChannelVsArea->GetXaxis()->SetTitle("Diamond Cluster Channel no");
	hOneClusterHitChannelVsArea->GetYaxis()->SetTitle("Area no");

	hOneClusterHitChannelVsFiducialArea = new TH2F("hOneClusterHitChannelVsFiducialArea","hOneClusterHitChannelVsArea",128,0,128,
	        settings->getSelectionFidCuts()->getNFidCuts()+2,-1, settings->getSelectionFidCuts()->getNFidCuts()+1);
	hOneClusterHitChannelVsFiducialArea->GetXaxis()->SetTitle("Diamond Cluster Channel no");
	hOneClusterHitChannelVsFiducialArea->GetYaxis()->SetTitle("Fiducial Area");

	hOneClusterHitChannelAreaVsFiducialArea = new TH2F("hOneClusterHitChannelAreaVsFiducialArea","hOneClusterHitChannelVsArea",
	        settings->getNDiaDetectorAreas()+1,-1,settings->getNDiaDetectorAreas(),
            settings->getSelectionFidCuts()->getNFidCuts()+2,-1, settings->getSelectionFidCuts()->getNFidCuts()+1);
	hOneClusterHitChannelAreaVsFiducialArea->GetXaxis()->SetTitle("areas");
	hOneClusterHitChannelAreaVsFiducialArea->GetYaxis()->SetTitle("Fiducial Area");

	hFidCut= new TH2F("hFidCut","hFidCut",256*4,0,255,256*4,0,255);
	hFidCut->GetXaxis()->SetTitle("FidCutValue in X");
	hFidCut->GetYaxis()->SetTitle("FidCutValue in Y");

	hFidCutOneDiamondCluster= new TH2F("hFidCut_oneDiamondCluster","hFidCut_oneDiamondCluster",256*4,0,255,256*4,0,255);
	hFidCutOneDiamondCluster->GetXaxis()->SetTitle("FidCutValue in X");
	hFidCutOneDiamondCluster->GetYaxis()->SetTitle("FidCutValue in Y");

	hClusterPosition=new TH1F("hClusterPositionDia","Events which have a valid Silicon Track",128,0,127);
	hClusterPosition->GetXaxis()->SetTitle("highes Cluster Channel Position");
	hClusterPosition->GetYaxis()->SetTitle("number of Events #");

	hClusterSizeVsChannelPos = new TH2F("hClusterSizeVsChannelPos","hClusterSizeVsChannelPos",10,-.5,9.5,128,0,127);
	hClusterSizeVsChannelPos->GetXaxis()->SetTitle("cluster size");
	hClusterSizeVsChannelPos->GetYaxis()->SetTitle("channel of highest signal in cluster");

	h3dDiamond = new TH1F("h3dDiamond","Sum of Charge for all 18 3d-channels",4096,0,4095);
	h3dDiamond_hit = new TH1F("h3dDiamond_hit","Sum of Charge for all 18 3d-channels with a Hit",4096,0,4095);
	hNoDiamond = new TH1F("hNoDiamond","Sum of Charge for all 18 no-channels",4096,0,4095);
	hNoDiamond_hit = new TH1F("hNoDiamond_hit","Sum of Charge for all 18 no-channels with a Hit",4096,0,4095);


	TString name = TString::Format("hEtaVsLeftChannelNo");
	hEtaVsLeftChannelNo = new TH2F(name,name,256,0,1,128,0,127);
	hEtaVsLeftChannelNo->GetXaxis()->SetTitle("#eta");
	hEtaVsLeftChannelNo->GetYaxis()->SetTitle("left channel of #eta");
	name = TString::Format("hEtaCMNcorrectedVsLeftChannelNo");
	hEtaCMNcorrectedVsLeftChannelNo = new TH2F(name,name,256,0,1,128,0,127);
	hEtaCMNcorrectedVsLeftChannelNo->GetXaxis()->SetTitle("#eta_{CMNcorrected}");
	hEtaCMNcorrectedVsLeftChannelNo->GetYaxis()->SetTitle("left channel of #eta ");
	for (UInt_t i=1;i<=settings->getSelectionFidCuts()->getNFidCuts();i++){
		Float_t xMin = settings->getSelectionFidCuts()->getMinFiducialX(i);
		Float_t xMax = settings->getSelectionFidCuts()->getMaxFiducialX(i);
		UInt_t xBins = (Int_t)(2*(xMax-xMin));
		name = TString::Format("hChargeVsFidX_HitInFidCutNo%d",i);
		TH2F* histo = new TH2F(name,name,4096,0,4095,xBins,xMax,xMin);
		histo->GetXaxis()->SetTitle("charge of diamond Hit");
		histo->GetYaxis()->SetTitle("FiducialCut X/ch");
		hChargeVsFidX.push_back(histo);
		if(verbosity){
			cout<<"added Histogram: "<<name
					<< " X:"<<histo->GetNbinsX()<<"_"<<histo->GetXaxis()->GetXmin()<<":"<<histo->GetXaxis()->GetXmax()
					<<" -Y:"<<histo->GetNbinsY()<<"_"<<histo->GetYaxis()->GetXmin()<<":"<<histo->GetYaxis()->GetXmax()<<endl;
		}
		Float_t yMin = settings->getSelectionFidCuts()->getMinFiducialY(i);
		Float_t yMax =  settings->getSelectionFidCuts()->getMaxFiducialY(i);
		xBins = (Int_t)(2*(yMax-yMin));
		name = TString::Format("hChargeVsFidY_HitInFidCutNo%d",i);
		histo = new TH2F(name,name,4096,0,4095,xBins,yMax,yMin);
		histo->GetXaxis()->SetTitle("charge of diamond Hit");
		histo->GetYaxis()->SetTitle("FiducialCut Y/ch");
		hChargeVsFidY.push_back(histo);
		if(verbosity){
			cout<<"added Histogram: "<<name
					<< " X:"<<histo->GetNbinsX()<<"_"<<histo->GetXaxis()->GetXmin()<<":"<<histo->GetXaxis()->GetXmax()
					<<" -Y:"<<histo->GetNbinsY()<<"_"<<histo->GetYaxis()->GetXmin()<<":"<<histo->GetYaxis()->GetXmax()<<endl;
		}
		if(verbosity>8&&verbosity%2==1){
			cout<<"Press a key and enter to continue.\t"<<flush;
			char t;
			cin>>t;
		}
	}
	name = TString::Format("hChargeVsFidCut");
	Float_t xmin = settings->getSelectionFidCuts()->getMinFiducialX();
	Float_t xmax = settings->getSelectionFidCuts()->getMaxFiducialX();
	Float_t ymin = settings->getSelectionFidCuts()->getMinFiducialY();
	Float_t ymax = settings->getSelectionFidCuts()->getMaxFiducialY();
	Float_t deltaX = xmax - xmin;
	xmin = xmin - .1 * deltaX;
	xmax = xmax + .1 * deltaX;
	Float_t deltaY = ymax - ymin;
	ymin = ymin - .1 * deltaY;
	ymax = ymax + .1 * deltaY;
	Int_t xBins = (Int_t)(3*(xmax-xmin));
	Int_t yBins = (Int_t)(3*(ymax-ymin));
	hChargeVsFidCut = new TH3F(name,name,xBins,xmin,xmax,yBins,ymin,ymax,4096,0,4096);
	hChargeVsFidCut->GetXaxis()->SetTitle("FiducialValue X/ch");
	hChargeVsFidCut->GetYaxis()->SetTitle("FiducialValue Y/ch");
	hChargeVsFidCut->GetZaxis()->SetTitle("charge");
	if(verbosity)cout<<"added Histogram: "<<name<<endl;
	Float_t chBegin = settings->getMinDiamondChannel();
	Float_t chEnd = settings->getMaxDiamondChannel();
	Float_t deltaCh = chEnd - chBegin;
	chBegin -= 0.1 * deltaCh;
	chEnd   += 0.1 * deltaCh;
	name = TString::Format("hFidCutXvsDiamondClusterChannelPos");
	yBins=3*(chEnd-chBegin);
	hFidCutXvsChannelPos = new TH2F(name,name,xBins,xmin,xmax,yBins,chBegin,chEnd);
	hFidCutXvsChannelPos->GetXaxis()->SetTitle("FiducialValue X / ch");
	hFidCutXvsChannelPos->GetYaxis()->SetTitle("DiamondCluster Channel Pos / ch");
	hFidCutXvsChannelPos->GetYaxis()->SetTitle("# number of entries");

	hTwoClustersArea = new TH2F("hTwoClustersArea","hTwoClustersArea",6,-1.5,4.5,6,-1.5,4.5);
	hTwoClustersArea->GetXaxis()->SetTitle("Area of Cluster No. 1");
	hTwoClustersArea->GetYaxis()->SetTitle("Area of CLuster No. 2");

	hNDiaClusters = new TH1F("hNDiaClusters","hNDiaClusters",10,-.5,9.5);
	hNDiaClusters->GetXaxis()->SetTitle("no of Dia Clusters");
	hNDiaClusters->GetYaxis()->SetTitle("no of entries #");
}

void TAnalysisOfSelection::saveDiamondAreaHistos(){
	settings->diamondPattern.showPatterns();
	cout << "Save: "<<histoLandauDistribution2DNoBorderSeed_unmasked->GetName() << endl;
	histSaver->SaveHistogram(histoLandauDistribution2DNoBorderSeed_unmasked);
	for(Int_t area=0;area<settings->getNDiaDetectorAreas();area++){

		Int_t chLow = settings->getDiaDetectorArea(area).first;
		Int_t chHigh =  settings->getDiaDetectorArea(area).second;

		cout <<" Area: " << area << ", ch: : "<< chLow <<" - "<< chHigh <<endl;

		TString name = TString::Format("hEtaVsLeftChannelNoArea%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<"save: "<<name<<endl;
		TH2F* hEtaVsLeftChannelNoArea = (TH2F*)hEtaVsLeftChannelNo->Clone(name);
		hEtaVsLeftChannelNoArea->SetTitle(name);
		hEtaVsLeftChannelNoArea->Draw();
		hEtaVsLeftChannelNoArea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);
		histSaver->SaveHistogram(hEtaVsLeftChannelNoArea);
		delete hEtaVsLeftChannelNoArea;

		Int_t minBin = hClusterSizeVsChannelPos->GetYaxis()->FindBin(chLow);
		Int_t maxBin = hClusterSizeVsChannelPos->GetYaxis()->FindBin(chHigh);
		name = TString::Format("hEta_Area_%d_ch_%d_%d",area,chLow,chHigh);
		TH1F* hEtaArea = (TH1F*)hEtaVsLeftChannelNo->ProjectionX(name,minBin,maxBin);
		if(!hEtaArea){
			cout<<"Projection does not work: "<<name;
			if(verbosity%2==1){
				cout<<".\tPress a key and enter to continue..."<<flush;
				char t;cin>>t;
			}
			else
				cout<<endl;
		}
		else
			hEtaVsLeftChannelNo->SetTitle(name);



		name = TString::Format("hEtaCMNcorrectedVsLeftChannelNoArea%d_ch_%d_%d",area,chLow,chHigh);
		TH2F* hEtaCMNcorrectedVsLeftChannelNoArea = (TH2F*)hEtaCMNcorrectedVsLeftChannelNo->Clone(name);
		if(hEtaCMNcorrectedVsLeftChannelNoArea){
			hEtaCMNcorrectedVsLeftChannelNoArea->SetTitle(name);
			hEtaCMNcorrectedVsLeftChannelNoArea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);
		}
		minBin = hEtaCMNcorrectedVsLeftChannelNoArea->GetYaxis()->FindBin(chLow);
		maxBin = hEtaCMNcorrectedVsLeftChannelNoArea->GetYaxis()->FindBin(chHigh);
		name = TString::Format("hEtaCMNcorrected_Area_%d_ch_%d_%d",area,chLow,chHigh);
		TH1F* hEtaCMNcorrectedArea = (TH1F*)hEtaCMNcorrectedVsLeftChannelNoArea->ProjectionX(name,minBin,maxBin);
		if(!hEtaCMNcorrectedArea){
			cout<<"Projection does not work: "<<name;
			if(verbosity%2==1){
				cout<<".\t Press a key and enter to continue..."<<flush;
				char t;cin>>t;
			}
			else
				cout<<endl;
		}
		histSaver->SaveHistogram(hEtaCMNcorrectedArea);
		histSaver->SaveHistogram(hEtaCMNcorrectedVsLeftChannelNoArea);

		histSaver->SaveHistogram(hEtaArea);
		if(hEtaArea)delete hEtaArea;

		//		delete hEtaCMNcorrectedVsLeftChannelNoArea;
		//2d area channel vs cluster size
		name = TString::Format("hClusterSizeVsChannelPos_Area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH2F* hClusterSizeVsChannelPosArea = (TH2F*)hClusterSizeVsChannelPos->Clone(name);
		if(hClusterSizeVsChannelPosArea){
			hClusterSizeVsChannelPosArea->SetTitle(name);
			hClusterSizeVsChannelPosArea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);
			histSaver->SaveHistogram(hClusterSizeVsChannelPosArea,false);
		}
		if(hClusterSizeVsChannelPosArea) delete hClusterSizeVsChannelPosArea;

		Float_t yMax = hClusterSizeVsChannelPos->GetYaxis()->GetXmax();
		Float_t yMin = hClusterSizeVsChannelPos->GetYaxis()->GetXmin();
		int binMax = hClusterSizeVsChannelPos->GetYaxis()->FindBin(yMax);
		int binMin = hClusterSizeVsChannelPos->GetYaxis()->FindBin(yMin);
		name = TString::Format("hClusterSize_Area_%d_ch_%d_%d",area,chLow,chHigh);
		TH1F* hClusterSizeVsChannelPosAreaProjection = (TH1F*)hClusterSizeVsChannelPos->ProjectionX(name,binMin,binMax);
		if(hClusterSizeVsChannelPosAreaProjection){
			hClusterSizeVsChannelPosAreaProjection->SetTitle(name);
			hClusterSizeVsChannelPosAreaProjection->GetXaxis()->SetTitle("cluster size");
			hClusterSizeVsChannelPosAreaProjection->GetYaxis()->SetTitle("number of entries");
			histSaver->SaveHistogram(hClusterSizeVsChannelPosAreaProjection);
		}
		if(hClusterSizeVsChannelPosAreaProjection)delete hClusterSizeVsChannelPosAreaProjection;

		// 2d area clusterSize 1or 2 normal
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH2F* histoLandauDistribution2Darea = (TH2F*)histoLandauDistribution2D->Clone(name);
		if(histoLandauDistribution2Darea){
			histoLandauDistribution2Darea->SetTitle(name);
			histoLandauDistribution2Darea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);
			name = name + (TString)"_pfx";
			TProfile *prof = histoLandauDistribution2Darea->ProfileY(name);
			prof->GetXaxis()->SetRangeUser(chLow,chHigh);
			TF1* fit = new TF1("fit","pol1",chLow,chHigh);
			histSaver->Save1DProfileWithFitAndInfluence(prof, fit,true);
			if (prof) delete prof;
		}
		histSaver->SaveHistogram(histoLandauDistribution2Darea);
		if(histoLandauDistribution2Darea)delete histoLandauDistribution2Darea;

		// 2d area clusterSize 1or 2 noBorderSeeds
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_noBorderSeed_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH2F* histoLandauDistribution2DNoBorderSeedarea = (TH2F*)histoLandauDistribution2DNoBorderSeed->Clone(name);
		if(histoLandauDistribution2DNoBorderSeedarea){
		histoLandauDistribution2DNoBorderSeedarea->SetTitle(name);
		histoLandauDistribution2DNoBorderSeedarea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);

        name = name + (TString)"_pfx";
        TProfile *prof = histoLandauDistribution2DNoBorderSeedarea->ProfileY(name);
        prof->GetXaxis()->SetRangeUser(chLow,chHigh);
        TF1* fit = new TF1("fit","pol1",chLow,chHigh);
        histSaver->Save1DProfileWithFitAndInfluence(prof, fit,true);
        if (prof) delete prof;
		}
		histSaver->SaveHistogram(histoLandauDistribution2DNoBorderSeedarea);
		if(histoLandauDistribution2DNoBorderSeedarea)
			delete histoLandauDistribution2DNoBorderSeedarea;


		// 2d area clusterSize 1or 2 noBorderHits
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_noBorderHit_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;


		TH2F* histoLandauDistribution2DNoBorderHitarea = (TH2F*)histoLandauDistribution2DNoBorderHit->Clone(name);
		histoLandauDistribution2DNoBorderHitarea->SetTitle(name);
		histoLandauDistribution2DNoBorderHitarea->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);


        name = name + (TString)"_pfx";
        TProfile *prof = histoLandauDistribution2DNoBorderHitarea->ProfileY(name);
        prof->GetXaxis()->SetRangeUser(chLow,chHigh);
        TF1* fit = new TF1("fit","pol1",chLow,chHigh);
        histSaver->Save1DProfileWithFitAndInfluence(prof, fit,true);
        if (prof) delete prof;

		histSaver->SaveHistogram(histoLandauDistribution2DNoBorderHitarea);
		Int_t firstBin = histoLandauDistribution2DNoBorderHitarea->GetYaxis()->GetFirst();
		Int_t lastBin = histoLandauDistribution2DNoBorderHitarea->GetYaxis()->GetLast();
		vector<TH1F*> vecHistos;
		Double_t max =0;
		int k = 0;
		THStack* stack = new THStack("hStack",TString::Format("Charge of Cluster per Channel - Area %d",area));
		for(Int_t bin = firstBin;bin<lastBin;bin++){
			Int_t ch = histoLandauDistribution2DNoBorderHitarea->GetYaxis()->GetBinCenter(bin);
			if(!settings->isInDiaDetectorArea(ch,area))
				continue;
			if(settings->isDet_channel_screened(TPlaneProperties::getDetDiamond(),ch))
				continue;
			name = TString::Format("hChargeOfCluster_ch%d_area%d",ch,area);
			TH1F * histo = (TH1F*) histoLandauDistribution2DNoBorderHitarea->ProjectionX(name,bin,bin);
			if(verbosity>4)cout<<name<<" - "<<histo<<endl;
			if(!histo)
				continue;
			name = TString::Format("ch %d, area %d",ch,area);
			histo->SetTitle(name);
			histo->Draw();
			histo->Rebin();
			histo->Rebin();
			histo->GetXaxis()->SetTitle("Charge of Cluster");
			histo->GetYaxis()->SetTitle("number of entries #");
			histo->SetLineColor(settings->GetColor(k));
			k++;
			max = TMath::Max(max,histo->GetBinContent(histo->GetMaximumBin()));
			vecHistos.push_back(histo);
			stack->Add(histo);
		}
		name = TString::Format("cChargePerChannel_area%d",area);
		max*=1.1;
		TCanvas *c1 = new TCanvas(name,name,1024, 768);
		c1->cd();
		//		Float_t xmin = histoLandauDistribution2DNoBorderHit->GetXaxis()->GetXmin();
		//		Float_t xmax = histoLandauDistribution2DNoBorderHit->GetXaxis()->GetXmax();
		//		Int_t bins = histoLandauDistribution2DNoBorderHit->GetXaxis()->GetNbins();
		stack->Draw("nostack");
		TLegend *leg = c1->BuildLegend();
		leg->SetFillColor(kWhite);
		leg->Draw();
		histSaver->SaveCanvas(c1);


		//		delete histoLandauDistribution2DNoBorderSeedarea;

		// 2d area clusterSize 1 or 2 unmasked
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_unmasked_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		UInt_t binLow=0;
		UInt_t binHigh=-1;
		TH2F* histoLandauDistribution2DareaUnmasked = (TH2F*)histoLandauDistribution2D_unmasked->Clone(name);
		if(histoLandauDistribution2DareaUnmasked){
			histoLandauDistribution2DareaUnmasked->SetTitle(name);
			histoLandauDistribution2DareaUnmasked->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);


	        name = name + (TString)"_pfx";
	        TProfile *prof = histoLandauDistribution2DareaUnmasked->ProfileY(name);
	        prof->GetXaxis()->SetRangeUser(chLow,chHigh);
	        TF1* fit = new TF1("fit","pol1",chLow,chHigh);
	        histSaver->Save1DProfileWithFitAndInfluence(prof, fit,true);
	        if (prof) delete prof;

            binLow = histoLandauDistribution2DareaUnmasked->GetYaxis()->FindBin(chLow);
			binHigh = histoLandauDistribution2DareaUnmasked->GetYaxis()->FindBin(chHigh);
		}
		histSaver->SaveHistogram(histoLandauDistribution2DareaUnmasked);

		if(histoLandauDistribution2DareaUnmasked)delete histoLandauDistribution2DareaUnmasked;

		// 2d area clusterSize 1or 2 unmasked noBorderSeed
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_noBorderSeed_unmasked_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH2F* histoLandauDistribution2DareaUnmaskedNoBorderSeed = (TH2F*)histoLandauDistribution2DNoBorderSeed_unmasked->Clone(name);
		if(histoLandauDistribution2DareaUnmaskedNoBorderSeed){
			histoLandauDistribution2DareaUnmaskedNoBorderSeed->SetTitle(name);
			histoLandauDistribution2DareaUnmaskedNoBorderSeed->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);


            name = name + (TString)"_pfx";
            TProfile *prof = histoLandauDistribution2DareaUnmaskedNoBorderSeed->ProfileY(name);
            prof->GetXaxis()->SetRangeUser(chLow,chHigh);
            TF1* fit = new TF1("fit","pol1",chLow,chHigh);
            histSaver->Save1DProfileWithFitAndInfluence(prof, fit,true);
            if (prof) delete prof;

		}
		histSaver->SaveHistogram(histoLandauDistribution2DareaUnmaskedNoBorderSeed);
		if(histoLandauDistribution2DareaUnmaskedNoBorderSeed) delete histoLandauDistribution2DareaUnmaskedNoBorderSeed;

		// 2d area clusterSize 1or 2 unmasked noBorderSeed
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_noBorderHit_unmasked_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH2F* histoLandauDistribution2DareaUnmaskedNoBorderHit = (TH2F*)histoLandauDistribution2DNoBorderHit_unmasked->Clone(name);
		if(histoLandauDistribution2DareaUnmaskedNoBorderHit){
			histoLandauDistribution2DareaUnmaskedNoBorderHit->SetTitle(name);
			histoLandauDistribution2DareaUnmaskedNoBorderHit->GetYaxis()->SetRangeUser(chLow-1,chHigh+1);
            name = name + (TString)"_pfx";
            TProfile *prof = histoLandauDistribution2DareaUnmaskedNoBorderHit->ProfileY(name);
            prof->GetXaxis()->SetRangeUser(chLow,chHigh);
            TF1* fit = new TF1("fit","pol1",chLow,chHigh);
            histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
		}
		histSaver->SaveHistogram(histoLandauDistribution2DareaUnmaskedNoBorderHit);
		if(histoLandauDistribution2DareaUnmaskedNoBorderHit)delete histoLandauDistribution2DareaUnmaskedNoBorderHit;

		LandauGaussFit landaugaus;
		/******* PROJECTIONS ******/
		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		TH1F* hProjection = (TH1F*)histoLandauDistribution2D->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;

		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_NoBorderSeed_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		hProjection = (TH1F*)histoLandauDistribution2DNoBorderSeed->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;

		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_NoBorderHit_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		hProjection = (TH1F*)histoLandauDistribution2DNoBorderHit->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;

		//unmasked projections
		name = TString::Format("hChargeOfCluster_ClusterSizeUnmasked_1_2_area_%d_ch_%d_%d",area,chLow,chHigh);
		if(verbosity)cout<<name<<endl;
		hProjection = (TH1F*)histoLandauDistribution2D_unmasked->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;

		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_NoBorderSeedUnmasked_area_%d_ch_%d_%d",area,chLow,chHigh);
		hProjection = (TH1F*)histoLandauDistribution2DNoBorderSeed_unmasked->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;

		name = TString::Format("hChargeOfCluster_ClusterSize_1_2_NoBorderHitUnmasked_area_%d_ch_%d_%d",area,chLow,chHigh);
		hProjection = (TH1F*)histoLandauDistribution2DNoBorderHit_unmasked->ProjectionX(name,binLow,binHigh);
		if(hProjection){
			hProjection->SetTitle(name);
			hProjection->GetXaxis()->SetTitle(TString::Format("ChargeOfCluster in area %d",area));
			hProjection->GetYaxis()->SetTitle("number of entries");
			landaugaus.doLandauGaussFit(hProjection);
			histSaver->SaveHistogramLandau(hProjection);
			delete hProjection;
		}
		hProjection=0;
	}
	if(hEtaVsLeftChannelNo){
		if(verbosity)cout<<"save: "<<hEtaVsLeftChannelNo->GetName()<<endl;;
		histSaver->SaveHistogram(hEtaVsLeftChannelNo);
		delete hEtaVsLeftChannelNo;
	}
	if(hEtaCMNcorrectedVsLeftChannelNo){
		histSaver->SaveHistogram(hEtaCMNcorrectedVsLeftChannelNo);
		delete hEtaCMNcorrectedVsLeftChannelNo;
	}
}
void TAnalysisOfSelection::saveFidCutHistos(){
	int bins = hChargeVsFidCut->GetNbinsX() * hChargeVsFidCut->GetNbinsY() * hChargeVsFidCut->GetNbinsZ();
	if(verbosity)cout<<"create hChargeVsFidCutProfile (this take some time since is doing a Project3D Profile for "<<bins<<endl;
	TH2F* hChargeVsFidCutProfile = (TH2F*)hChargeVsFidCut->Project3DProfile("yx");
	if(verbosity)cout<<"draw hChargeVsFidCutProfile"<<endl;
	if(	hChargeVsFidCutProfile){
		hChargeVsFidCutProfile->Draw();
		hChargeVsFidCutProfile->SetName("hMeanChargeVsFiducialCutPosition");
		hChargeVsFidCutProfile->SetTitle("MeanChargeVsFiducialCutPosition");
		hChargeVsFidCutProfile->GetXaxis()->SetTitle("fiducialCut X/ch");
		hChargeVsFidCutProfile->GetYaxis()->SetTitle("fiducialCut Y/ch");
		hChargeVsFidCutProfile->GetZaxis()->SetTitle("mean charge");
		TCanvas *c3 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hChargeVsFidCutProfile);
		c3->SetName(hChargeVsFidCutProfile->GetName());
		if(verbosity)cout<<"save hChargeVsFidCutProfile"<<endl;
		histSaver->SaveCanvas(c3);
	}
	//	Histogram(hChargeVsFidCutProfile);
	//	cout<<"save hChargeVsFidCut"<<endl;
	//	histSaver->SaveHistogramROOT(hChargeVsFidCut);Profile;
	histSaver->SaveHistogram(hFidCutXvsChannelPos);
	if(hFidCutXvsChannelPos)delete hFidCutXvsChannelPos;
	if(verbosity)cout<<"save hChargeVsFidCut x/Y"<<endl;
	Float_t xmin,xmax,xCorrect,ymin,ymax,yCorrect;
	for (UInt_t i=0;i<settings->getSelectionFidCuts()->getNFidCuts();i++){
		if(i<hChargeVsFidX.size()){
			if(verbosity)cout<<"save: "<<hChargeVsFidX[i]->GetTitle()<<endl;
			histSaver->OptimizeXYRange(hChargeVsFidX[i]);
			histSaver->SaveHistogram(hChargeVsFidX[i]);
		}
		if(i<hChargeVsFidY.size() ){
			if(verbosity)cout<<"save: "<<hChargeVsFidY[i]->GetTitle()<<endl;
			histSaver->OptimizeXYRange(hChargeVsFidY[i]);
			histSaver->SaveHistogram(hChargeVsFidY[i]);
		}
		TString name  = TString::Format("hMeanChargeFiducialCutNo%d",i+1);
		if(verbosity)cout<<name<<" "<<hChargeVsFidCutProfile<<flush;
		TH2F* hMeanChargeArea =  (TH2F*)hChargeVsFidCutProfile->Clone(name);
		hMeanChargeArea->SetTitle(name);
		if(verbosity)cout<<"."<<flush;
		xmin = settings->getSelectionFidCuts()->getMinFiducialX(i+1);
		xmax = settings->getSelectionFidCuts()->getMaxFiducialX(i+1);
		xCorrect = 0.1*(xmax-xmin);
		if(verbosity)cout<<"."<<flush;
		ymin = settings->getSelectionFidCuts()->getMinFiducialY(i+1);
		ymax = settings->getSelectionFidCuts()->getMaxFiducialY(i+1);
		yCorrect = 0.1*(ymax-ymin);
		if(verbosity)cout<<"."<<flush;
		hMeanChargeArea->GetXaxis()->SetRangeUser(xmin-xCorrect,xmax+xCorrect);
		hMeanChargeArea->GetYaxis()->SetRangeUser(ymin-yCorrect,ymax+yCorrect);
		if(verbosity)cout<<"."<<flush;
		histSaver->SaveHistogram(hMeanChargeArea);
		if(verbosity)cout<<"."<<flush;
		if(hMeanChargeArea)delete hMeanChargeArea;
		if(verbosity)cout<<"#"<<endl;
	}
//	histSaver->SaveHistogramROOT(hChargeVsFidCut);
	if(verbosity)cout<<"create hEntriesOfMeanChargeHisto "<<flush;
	TH2F* hEntriesOfMeanChargeHisto = (TH2F*)hChargeVsFidCut->Project3D("yx");
	if(hEntriesOfMeanChargeHisto){
		if(verbosity)cout<<"."<<flush;
		hEntriesOfMeanChargeHisto->SetName("hEntriesOfMeanChargeHisto");
		hEntriesOfMeanChargeHisto->SetTitle("hEntriesOfMeanChargeHisto");
		if(verbosity)cout<<"."<<flush;
		hEntriesOfMeanChargeHisto->GetXaxis()->SetTitle("fidCut X/ch");
		hEntriesOfMeanChargeHisto->GetYaxis()->SetTitle("fidCut Y/ch");
		hEntriesOfMeanChargeHisto->GetZaxis()->SetTitle("number of entries #");
		if(verbosity)cout<<"."<<flush;
		histSaver->SaveHistogram(hEntriesOfMeanChargeHisto);
		if(verbosity)cout<<"#"<<flush;
	}

	if(verbosity)cout<<endl;

	if(hChargeVsFidCut) delete hChargeVsFidCut;
	if(hChargeVsFidCutProfile) delete hChargeVsFidCutProfile;
}
void TAnalysisOfSelection::saveHistos()
{
	TH1F *histo=0;
	cout<<"\n\nSAVE HISTOGRAMS!!!!!"<<endl;
	savePHvsEventNoAreaPlots(histSaver,settings,hPHVsEventNo_Areas,xDivisions,yDivisions);
	saveDiamondAreaHistos();
	saveFidCutHistos();
	//	cout<<"AREAS: "<<settings->getNDiaDetectorAreas()<<endl;
	//	char t;cin >>t;
	histSaver->SaveHistogram(hTwoClustersArea);
	if(hTwoClustersArea)delete hTwoClustersArea;
	LandauGaussFit landauGauss;
	histSaver->OptimizeXYRange(histoLandauDistribution2D_unmasked);
	histSaver->OptimizeXYRange(histoLandauDistribution2D);
	if(verbosity)cout<<"unmasked: "<<histoLandauDistribution2D_unmasked->GetEntries()<<"\nmasked: "<<histoLandauDistribution2D->GetEntries()<<endl;
	if(verbosity)cout<< "save "<<histoLandauDistribution2D->GetName()<<endl;
	histSaver->SaveHistogram(histoLandauDistribution2D);
	if(verbosity)cout<< "save "<<histoLandauDistribution2D_unmasked->GetName()<<endl;
	histSaver->SaveHistogram(histoLandauDistribution2D_unmasked);

	stringstream name;
	name<<"hEtaVsSignalLeftOfEta";
	TH2F* histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalLeftOfEta,this->vecEta);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("Signal left of #eta");
		histo2d->GetYaxis()->SetTitle("#eta");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}

	name.str("");name.clear();
	name<<"hEtaVsSignalRightOfEta";
	histo2d = histSaver->CreateScatterHisto(name.str(),this->vecSignalRightOfEta,this->vecEta);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("Signal right of #eta");
		histo2d->GetYaxis()->SetTitle("#eta");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	name.str("");name.clear();
	name<<"hTest";
	histo = histSaver->CreateDistributionHisto(name.str(),vecTest,512);
	histSaver->SaveHistogram(histo);
	delete histo;

	vector<Float_t> vecRightFactor, vecLeftFactor,vecRightOverHighest,vecLeftOverHighest;
	for(UInt_t i=0;i<vecSignalLeftOfHighest.size()&&i<vecSignalRightOfHighest.size()&&i<vecClusterCharge.size();i++){
		vecRightFactor.push_back(vecSignalRightOfHighest.at(i)/vecClusterCharge.at(i));
		vecLeftFactor.push_back(vecSignalLeftOfHighest.at(i)/vecClusterCharge.at(i));
		vecLeftOverHighest.push_back(vecSignalLeftOfHighest.at(i)/vecHighestSignal.at(i));
		vecRightOverHighest.push_back(vecSignalRightOfHighest.at(i)/vecHighestSignal.at(i));
	}
	name.str("");name.clear();
	name<<"hSignalLeftOverHighestVsSignalRightOverHighest";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecRightOverHighest,vecLeftOverHighest);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("signal left of highest signal over highest signal");
		histo2d->GetYaxis()->SetTitle("signal right of highest signal over highest signal");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}

	name.str("");name.clear();
	name<<"hSignalRightVsHighestSignal";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalRightOfHighest,vecHighestSignal);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("highest signal");
		histo2d->GetYaxis()->SetTitle("signal right of highest signal ");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	name.str("");name.clear();
	name<<"hSignalLeftVsHighestSignal";
	histo2d = histSaver->CreateScatterHisto(name.str(),vecSignalLeftOfHighest,vecHighestSignal);
	if(histo2d){
		histo2d->GetXaxis()->SetTitle("highest signal");
		histo2d->GetYaxis()->SetTitle("signal left of highest signal ");
		histSaver->SaveHistogram(histo2d);
		delete histo2d;
	}
	name.str("");name.clear();
	name<<"hSignalLeftOfHighest";
	TH1F* histoLeft = histSaver->CreateDistributionHisto(name.str(),vecLeftFactor);
	histSaver->SaveHistogram(histoLeft);
	name.str("");name.clear();
	name<<"hSignalRightOfHighest";
	TH1F* histoRight = histSaver->CreateDistributionHisto(name.str(),vecRightFactor);
	histSaver->SaveHistogram(histoRight);

	name.str("");name.clear();
	name<<"cSignalNextToHighest";
	histoLeft->SetLineColor(kBlue);
	histoRight->SetLineColor(kRed);
	Float_t max = TMath::Max(histoLeft->GetMaximum(),histoRight->GetMaximum());
	histoLeft->SetMaximum(max);
	histoRight->SetMaximum(max);
	histSaver->SaveTwoHistos(name.str(),histoLeft,histoRight,1);
	if(histoLeft) delete histoLeft;
	if(histoRight) delete histoRight;

	histSaver->SaveHistogram(hClusterSizeVsChannelPos,false);

	vector <Float_t> vecMP;
	vector <Float_t> vecClusSize;
	vector <Float_t> vecWidth;
	vector <Float_t> vecXError;
	vector <Float_t> vecHistoMax;
	vector <Float_t> vecHistoMean;
	vector <Float_t> vecHistoMeanGaus;
	vector <Float_t> vecHistoMeanLandau;
	TH1F* histoClusSize=0;
	if(histoLandauDistribution){
		histoClusSize = (TH1F*)histoLandauDistribution->ProjectionY("ClusterSizeDiamond",0,4096);
		histo = (TH1F*)histoLandauDistribution->ProjectionX("hPulseHeightDiamondAll",0,8);
	}

	if(histo==0){

		cout<<"Couldn't creaye histo hPulseHeightDiamondAll"<<flush;
		if(verbosity%2==1){cout<<"\tPress a key to continue..."<<endl;char t;cin>>t;}
	}

	Float_t histoMean,histoMax,histoRMS,histoMeanGausFit;
	TF1* gausFit=0;


	Float_t width=0;
	Float_t MP=0;
	Float_t area=0;
	Float_t gWidth=0;
	Float_t xmax,xmin;

	if (histo){
		if(verbosity)cout<< "preparing "<<histo->GetName()<<endl;

		histo->GetYaxis()->SetTitle("number of Entries #");
		Double_t xmin,xmax;
		TF1* fit=0;
		histoMean = histo->GetMean();
		histoMax = histo->GetBinCenter(histo->GetMaximumBin());
		histoRMS = histo->GetRMS();
		xmin=histoMax-histoRMS;
		xmax=histoMax+histoRMS;
		gausFit = new TF1("gausFit","gaus",xmin,xmax);
		//	cout<<"gausFit: "<<gausFit<<endl;
		histo->Fit(gausFit,"0Q+","goff",xmin,xmax);
		fit = landauGauss.doLandauGaussFit(histo);
		//	cout <<"gausFit:"<<gausFit->GetTitle()<<" is a:"<< gausFit->ClassName()<<" "<<gausFit->GetNpar()<<endl;
		histoMeanGausFit = gausFit->GetParameter(1);
		vecWidth.push_back(fit->GetParameter(0));
		vecHistoMax.push_back(histoMax);
		vecHistoMean.push_back(histoMean);
		vecHistoMeanGaus.push_back(histoMeanGausFit);
		vecHistoMeanLandau.push_back(fit->GetParameter(1));

		if(verbosity)cout<< "saving "<<histo->GetName()<<endl;
		histSaver->SaveHistogramLandau(histo);
		width=fit->GetParameter(0);
		MP = fit->GetParameter(1);
		area = fit->GetParameter(2);
		gWidth = fit->GetParameter(3);
		//	vecMP.push_back(MP);
		//	vecWidth.push_back(width);
		//vecClusSize.push_back(0);
		//vecXError.push_back(0);
	}

	name.str("");
	name.clear();
	name<< "hPulseHeigthDiamond_1_2_ClusterSize";
	TH1F* histo12 =0;

	if (histoLandauDistribution){
		Int_t binLow = histoLandauDistribution->GetYaxis()->FindBin(1);
		Int_t binHigh = histoLandauDistribution->GetYaxis()->FindBin(2);
		histo12 = (TH1F*)histoLandauDistribution->ProjectionX(name.str().c_str(),binLow,binHigh);
		landauGauss.doLandauGaussFit(histo12);
	}
	else{
		cerr<<"histoLandauDistribution is not valid, Press a key to continue..."<<flush;
		if(verbosity%2==1){char t; cin>>t;}
	}
	if(histo12==0){
		cout<<"'hPulseHeigthDiamond_1_2_ClusterSize' == 0"<<flush;
		if(verbosity%2==1){char t; cin>>t;}
	}
	else if(histo12->GetEntries()==0){
		cout<<"'hPulseHeigthDiamond_1_2_ClusterSize' has 0 entries"<<flush;
		if(verbosity%2==1){char t; cin>>t;}
	}
	else
		cout<<"'hPulseHeigthDiamond_1_2_ClusterSize' has "<< histo12->GetEntries()<<" entries"<<endl;
	if(histo==0){
		cout<<"histo does not exist!";
		if(verbosity%2==1){
			cout<<"\t Press a key and enter to confirm."<<flush;
			char t;cin>>t;
		}
		else
			cout<<endl;
	}
	if(histo12!=0 && histo!=0){
		if(verbosity) cout<<"histo: "<<histo <<" or histo12: "<<histo12<<" are valid..."<<endl;
		histo12->SetTitle(name.str().c_str());
		histo12->GetYaxis()->SetTitle("number of Entries #");
		if(histo12->GetEntries()>0){
			TF1* fitCS12=0;
			int nTries=0;
			while(histo->GetMaximum()<histo->GetEntries()*0.1&&nTries<5)
				histo->Rebin(),nTries++;
			histoMean = histo->GetMean();
			histoMax = histo->GetBinCenter(histo->GetMaximumBin());
			histoRMS = histo->GetRMS();
			xmin=histoMax-histoRMS, xmax=histoMax+histoRMS;
			if(gausFit!=0)delete gausFit;
			if (verbosity) cout<< "Fitting LandauGaus in "<<histo->GetName()<<endl;
			gausFit = new TF1("gausFit","gaus",xmin,xmax);
			histo->Fit(gausFit,"0+Q","goff",xmin,xmax);
			histoMeanGausFit = gausFit->GetParameter(1);
			if(fitCS12!=0)delete fitCS12;
			fitCS12 = landauGauss.doLandauGaussFit(histo);
		}
		else{
			cout<<"1_2 Cluster plot is empty....."<<endl;
		}
		if(histo12){
			cout<<"Save HISTOGRAM: "<<histo12->GetName()<<endl;
			histSaver->SaveHistogramLandau(histo12);
			delete histo12;
		}
	}
	else
		cerr<<"histo: "<<histo <<" or histo12: "<<histo12<<" are not initialized correctly"<<endl;

	for(UInt_t clusSize=1;clusSize<8;clusSize++){

		stringstream name;
		name<< "hPulseHeigthDiamond_"<<clusSize<<"_ClusterSize";
		if(verbosity>2)cout<<"save "<<name.str()<<endl;
		TH1F* histo = (TH1F*)histoLandauDistribution->ProjectionX(name.str().c_str(),clusSize,clusSize);
		if(histo==0) {
			cout<<"TAnalysisOfSelection:: saverHistos ==> oooh Boy, something went terribly wrong, Lukas you better fix it! NOW!"<<endl;
			continue;
		}
		histo->SetTitle(name.str().c_str());
		histo->GetYaxis()->SetTitle("number of Entries #");
		TF1* fitCS=0;
		if(clusSize<5&& histo->GetEntries()>0){
			int nTries=0;
			while(histo->GetMaximum()<histo->GetEntries()*0.1&&nTries<5)
				histo->Rebin(),nTries++;
			histoMean = histo->GetMean();
			histoMax = histo->GetBinCenter(histo->GetMaximumBin());
			histoRMS = histo->GetRMS();
			xmin=histoMax-histoRMS, xmax=histoMax+histoRMS;
			if(gausFit!=0)delete gausFit;
			gausFit = new TF1("gausFit","gaus",xmin,xmax);
			histo->Fit(gausFit,"0+Q","sames+",xmin,xmax);
			histoMeanGausFit = gausFit->GetParameter(1);
			if(fitCS!=0)delete fitCS;
			fitCS = landauGauss.doLandauGaussFit(histo);
			vecMP.push_back(fitCS->GetParameter(1));
			vecClusSize.push_back(clusSize);
			vecXError.push_back(.5);
			vecWidth.push_back(fitCS->GetParameter(0));
			vecHistoMax.push_back(histoMax);
			vecHistoMean.push_back(histoMean);
			vecHistoMeanGaus.push_back(histoMeanGausFit);
			vecHistoMeanLandau.push_back(fitCS->GetParameter(1));
			histSaver->SaveHistogramLandau(histo);
		}
		else
			histSaver->SaveHistogram(histo);

		if(histo) delete histo;
	}


	if(vecMP.size()==0)
		cout<<" size of vecMP is == 0"<<endl;
	else{
		if(verbosity)cout<<"Create ErrorGraph."<<vecClusSize.size()<<"."<<vecMP.size()<<"."<<vecXError.size()<<"."<<vecWidth.size()<<"."<<flush;
		TGraphErrors* graph = new TGraphErrors(vecMP.size(),&vecClusSize.at(0),&vecMP.at(0),&vecXError.at(0),&vecWidth.at(0));
		if(verbosity)cout<<"."<<flush;
		name.str("");name.clear();
		if(verbosity)cout<<"-"<<flush;
		name<<"MPV of Landau for one ClusterSizes";
		if(verbosity)cout<<"."<<flush;
		graph->SetTitle(name.str().c_str());
		if(verbosity)cout<<"+"<<flush;
		name.clear();name.str("");name.clear();
		name<<"hMPV_Landau_diff_ClusterSizes";
		if(verbosity)cout<<"!"<<flush;
		graph->SetName(name.str().c_str());
		if(verbosity)cout<<"#"<<flush;
		graph->Draw("APLE1 goff");
		if(verbosity)cout<<"$"<<flush;
		graph->GetXaxis()->SetTitle("Cluster Size");
		if(verbosity)cout<<"%"<<flush;
		graph->GetYaxis()->SetTitle("Most Probable Value of Landau");
		if(verbosity)cout<<"^"<<flush;
		graph->SetMarkerColor(kGreen);
		if(verbosity)cout<<"*"<<flush;
		graph->SetMarkerStyle(22);
		if(verbosity)cout<<"&"<<flush;
		graph->SetFillColor(kWhite);
		if(verbosity)cout<<"?"<<flush;
		graph->SetLineWidth(2);
		if(verbosity)cout<<"done"<<endl;
		if(verbosity)cout<<"Create Canvas"<<endl;
		TCanvas *c1= new TCanvas("cMVP_Landau_vs_ClusterSize","cMVP_Landau_vs_ClusterSize",800,600);
		c1->cd();
		Float_t xVal[] = {0,5};
		Float_t exVal[] = {0.5,0.5};
		Float_t yVal[] = {MP,MP};
		Float_t eyVal[]= {width,width};
		if(verbosity)cout<<"Create ErrorGraph MEAN"<<flush;
		TGraphErrors *gMVP = new TGraphErrors(2,xVal,yVal,exVal,eyVal);
		if(verbosity)cout<<"."<<flush;
		gMVP->SetName("gMPV_ALL");
		if(verbosity)cout<<"."<<flush;
		gMVP->SetTitle("MVP of all Clusters");
		if(verbosity)cout<<"."<<flush;
		gMVP->SetFillColor(kRed);
		if(verbosity)cout<<"."<<flush;
		gMVP->SetFillStyle(3002);
		if(verbosity)cout<<"."<<flush;
		gMVP->SetLineColor(kBlue);
		if(verbosity)cout<<"DONE"<<endl;
		if(verbosity)cout<<"Create MultiGraph"<<endl;
		TMultiGraph *mg = new TMultiGraph("mgMVP_ClusterSize","MVP of Landau vs. ClusterSize");
		mg->Add(gMVP,"3L");
		mg->Add(graph,"PLE1");
		if(verbosity)cout<<"Draw Canvas"<<endl;
		mg->Draw("a");
		mg->GetXaxis()->SetTitle("Cluster Size of Diamond");
		mg->GetYaxis()->SetTitle("MPV of Landau  ");
		mg->GetXaxis()->SetRangeUser(0.5,4.5);
		TLegend *leg = c1->BuildLegend(0.15,0.55,0.6,0.8);
		leg->SetFillColor(kWhite);
		if(verbosity)cout<<"Save Canvas"<<endl;
		histSaver->SaveCanvas(c1);

		//	TLine *lMVP = new TLine(graph->GetXaxis()->GetXmin(),MP,graph->GetXaxis()->GetXmax(),MP);
		//	TLine *lMVPplus = new TLine(graph->GetXaxis()->GetXmin(),MP+width,graph->GetXaxis()->GetXmax(),MP+width);
		//	TLine *lMVPminus = new TLine(graph->GetXaxis()->GetXmin(),MP-width,graph->GetXaxis()->GetXmax(),MP-width);
		if (graph)
			histSaver->SaveGraph(graph,name.str(),"APLE1");
		htmlLandau->addLandauDiamond(width,MP,area,gWidth);
		htmlLandau->addLandauDiamondTable(vecHistoMean,vecHistoMax,vecHistoMeanGaus,vecHistoMeanLandau);

		if(histoClusSize){
			histoClusSize->SetTitle("ClusterSize Diamond");
			histoClusSize->GetXaxis()->SetTitle("ClusterSize");
			histoClusSize->GetYaxis()->SetTitle("Number of Entries #");
			histSaver->SaveHistogram(histoClusSize);
			delete histoClusSize;
		}
		htmlLandau->addSection("ClusterSize Diamond",htmlLandau->putImageOfPath("ClusterSizeDiamond","png",50));
		//	if(fit!=0)delete fit;
		if (histo)delete histo;

		histSaver->SaveHistogram(histoLandauDistribution);
		if(histoLandauDistribution)
			delete histoLandauDistribution;
		if (histoLandauDistribution2D)
			delete histoLandauDistribution2D;
		if (histoLandauDistribution2D_unmasked)
			delete histoLandauDistribution2D_unmasked;
		if (mg)
			delete mg;
		if (c1)
			delete c1;
	}
	TCanvas* c1;
	if (hValidSiliconAndDiamondHit){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndDiamondHit,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndDiamondHit->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		if(hValidSiliconAndDiamondHit)delete hValidSiliconAndDiamondHit;
	}

	if (hValidSiliconAndOneDiamondHit){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndOneDiamondHit,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndOneDiamondHit->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		if(hValidSiliconAndOneDiamondHit)delete hValidSiliconAndOneDiamondHit;
	}

	if (hValidSiliconAndOneDiamondHitNotMasked){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndOneDiamondHitNotMasked,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndOneDiamondHitNotMasked->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		if(hValidSiliconAndOneDiamondHitNotMasked)delete hValidSiliconAndOneDiamondHitNotMasked;
	}
	if (hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		if(hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels)delete hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels;
	}
	if (hValidSiliconAndOneDiamondHitInOneArea){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndOneDiamondHitInOneArea,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndOneDiamondHitInOneArea->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		delete hValidSiliconAndOneDiamondHitInOneArea;
	}
	if (hValidSiliconAndOneDiamondHitInSameAreaAndFidCut){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hValidSiliconAndOneDiamondHitInSameAreaAndFidCut,true);
		c1->SetName(TString::Format("c%s",hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		delete hValidSiliconAndOneDiamondHitInSameAreaAndFidCut;
	}
	if (hOneClusterHitChannelVsArea){
        histSaver->SaveHistogram(hOneClusterHitChannelVsArea);
        delete hOneClusterHitChannelVsArea;
	}
	if (hOneClusterHitChannel){
        histSaver->SaveHistogram(hOneClusterHitChannel);
        delete hOneClusterHitChannel;
	}
	if (hOneClusterHitChannelVsFiducialArea){
        histSaver->SaveHistogram(hOneClusterHitChannelVsFiducialArea);
        delete hOneClusterHitChannelVsFiducialArea;
	}
	if (hOneClusterHitChannelAreaVsFiducialArea){
	    histSaver->SaveHistogram(hOneClusterHitChannelAreaVsFiducialArea);
	    delete hOneClusterHitChannelAreaVsFiducialArea;
	}


	if(hFidCut){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hFidCut,true);
		c1->SetName(TString::Format("c%s",hFidCut->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		delete hFidCut;
	}

	if(hFidCutOneDiamondCluster){
		c1 = settings->getSelectionFidCuts()->getAllFiducialCutsCanvas(hFidCutOneDiamondCluster,true);
		c1->SetName(TString::Format("c%s",hFidCutOneDiamondCluster->GetName()));
		histSaver->SaveCanvas(c1);
		delete c1;
		delete hFidCutOneDiamondCluster;
	}

	if(hClusterPosition){
		histSaver->SaveHistogram(hClusterPosition,0,1);
		delete hClusterPosition;
	}
	if(verbosity)cout<<"Save histos "<<endl;
	if(h3dDiamond&&verbosity)cout<<h3dDiamond->GetEntries()<<endl;
	if(h3dDiamond)histSaver->SaveHistogram(h3dDiamond,0,1);
	if(hNoDiamond&&verbosity)cout<<hNoDiamond->GetEntries()<<endl;
	if(hNoDiamond)histSaver->SaveHistogram(hNoDiamond,0,1);
	if(h3dDiamond_hit&&verbosity)cout<<h3dDiamond_hit->GetEntries()<<endl;
	histSaver->SaveHistogram(h3dDiamond_hit,0,1);
	if (hNoDiamond_hit&&verbosity)cout<<hNoDiamond_hit->GetEntries()<<endl;
	histSaver->SaveHistogram(hNoDiamond_hit,0,1);
	if(h3dDiamond_hit&&h3dDiamond_hit){
		if(verbosity)cout<<"Save canvas"<<endl;
		TCanvas *c2= new TCanvas("c2","c2",1024,800);
		c2->cd();
		h3dDiamond_hit->Draw();
		hNoDiamond_hit->SetLineColor(kBlue);
		hNoDiamond_hit->Draw("same");
		histSaver->SaveCanvas(c2);
		delete c2;
		delete h3dDiamond;
		delete hNoDiamond;
	}
	if(hNDiaClusters){
		histSaver->SaveHistogram(hNDiaClusters);
		delete hNDiaClusters;
	}
}

void TAnalysisOfSelection::analyseEvent()
{

	if(!eventReader->isValidTrack()) //just Tracks with Valid Silicon Track
		return;

	Float_t fiducialValueX= eventReader->getFiducialValueX();
	Float_t fiducialValueY = eventReader->getFiducialValueY();
	Int_t nDiaClusters = eventReader->getNClusters(TPlaneProperties::getDetDiamond());
	if(nDiaClusters<=0) // at least one Diamond Cluster
		return;
	TCluster cluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
	Float_t charge = cluster.getCharge(false);
	UInt_t clustSize = cluster.size();
	bool isMasked = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),cluster,false);
	bool isMaskedAdjacentChannels = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),cluster,true);
	Float_t pos = cluster.getPosition(settings->doCommonModeNoiseCorrection(),TCluster::maxValue,0);
	Int_t fidRegionIndex = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY)-1;
	Int_t area = settings->getDiaDetectorAreaOfChannel(pos,0);
	bool isInOneArea = !(area==-1);
//	if(verbosity && pos <= 20) settings->diamondPattern.showPatterns();
//	if(verbosity && pos <= 13 && pos>0) cout<<"\n*****"<<(int)pos<< " "<< fidRegionIndex << " / " << area <<" "<<nDiaClusters<<" "<<eventReader->isInOneFiducialArea()<<"\t"<<flush;

	hValidSiliconAndDiamondHit ->Fill(fiducialValueX,fiducialValueY);
	if (nDiaClusters==1){
		hValidSiliconAndOneDiamondHit ->Fill(fiducialValueX,fiducialValueY);
		if(!isMasked)
			hValidSiliconAndOneDiamondHitNotMasked ->Fill(fiducialValueX,fiducialValueY);
		if(!isMaskedAdjacentChannels)
			hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels ->Fill(fiducialValueX,fiducialValueY);
		if(isInOneArea)
			hValidSiliconAndOneDiamondHitInOneArea->Fill(fiducialValueX,fiducialValueY);
		if (area==fidRegionIndex)
			hValidSiliconAndOneDiamondHitInSameAreaAndFidCut->Fill(fiducialValueX,fiducialValueY);
		hOneClusterHitChannel->Fill(pos);
		hOneClusterHitChannelVsArea->Fill(pos,area);
	}
	Int_t fiducial_area = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
	hOneClusterHitChannelVsFiducialArea->Fill(pos,fiducial_area);
	hOneClusterHitChannelAreaVsFiducialArea->Fill(area,fiducial_area);
	if(!eventReader->isInOneFiducialArea())
		return;
	Int_t hitArea = TTransparentAnalysis::GetHitArea(settings,fiducialValueX,fiducialValueY,xDivisions,yDivisions);

	Float_t charge1 = cluster.getCharge(settings->doCommonModeNoiseCorrection());
	Float_t charge2 = cluster.getCharge(2,settings->doCommonModeNoiseCorrection());
	        fillPHvsEventNoAreaPlots(hitArea,charge1,charge2);
//	if(fidRegionIndex>0) cout<<"FidRegion: "<<fidRegionIndex<<" "<<endl;
	hFidCut->Fill(fiducialValueX,fiducialValueY);
	hNDiaClusters -> Fill(nDiaClusters);
	if (nDiaClusters > 1){
		TCluster cluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
		Float_t pos2 = cluster2.getPosition(settings->doCommonModeNoiseCorrection(),TCluster::maxValue,0);
		Float_t area2 = settings->getDiaDetectorAreaOfChannel(pos2);
		hTwoClustersArea->Fill(area,area2);
//		if(verbosity && pos <= 20) {cluster.Print();cluster2.Print();}
	}
	else if(nDiaClusters>1)
		return;
	hFidCutOneDiamondCluster->Fill(fiducialValueX,fiducialValueY);


	if(clustSize>8) clustSize=8;
	if(cluster.isSaturatedCluster())
		return;
	histoLandauDistribution->Fill(charge,clustSize);

	hClusterPosition->Fill(pos);
	//	if(isMaskedAdjacentChannels)
	//		return;
	//
	hClusterSizeVsChannelPos->Fill(clustSize,pos);
	if(area!=fidRegionIndex){
		if(verbosity>2)cout<<"\r"<<nEvent<<" Problem: "<<fiducialValueX<<"/"<<fiducialValueY<<"-->"<<fidRegionIndex<<" \t"<<pos<<"-->"<<area<<"_\n"<<flush;
		return;
	}
	if(clustSize<=2&&nDiaClusters==1&&area==fidRegionIndex){
		if(verbosity && pos <= 20) cout<<"FILL AT "<<pos<< " "<< fidRegionIndex << " / " << area <<" nClus"<< nDiaClusters<<" size"<<clustSize<<" "<<flush;
		//
		//		else{
		//			if(verbosity>5)cout<<nEvent<<" Good"<<endl;
		//		}
		vecTest.push_back(pos);
		Int_t leftChannel =-1;

		Float_t eta = cluster.getEta(leftChannel,false);
		Float_t etaCMNCor = cluster.getEta(true);
		Float_t signalLeftOfEta = cluster.getSignalOfChannel(leftChannel-1);
		Float_t signalRightOfEta = cluster.getSignalOfChannel(leftChannel+2);
		Int_t highestClusterPos = cluster.getHighestHitClusterPosition();
		Float_t highestSignal = cluster.getHighestSignal();
		Float_t leftOfHighestSignal = cluster.getSignal(highestClusterPos-1);
		Float_t rightOfHighestSignal = cluster.getSignal(highestClusterPos+1);
		this->vecSignalLeftOfEta.push_back(signalLeftOfEta);
		this->vecSignalRightOfEta.push_back(signalRightOfEta);
		this->vecSignalLeftOfHighest.push_back(leftOfHighestSignal);
		this->vecSignalRightOfHighest.push_back(rightOfHighestSignal);
		this->vecClusterCharge.push_back(charge);
		this->vecHighestSignal.push_back(highestSignal);
		this->vecEta.push_back(eta);
//		Float_t relPos = pos - (Int_t) (pos-.5);
//		if(clusSize==2) hEtaVsRelPos->Fill(eta,relPos);
		if(clustSize==2) hEtaVsLeftChannelNo ->Fill(eta,leftChannel);
		if(clustSize==2) hEtaCMNcorrectedVsLeftChannelNo ->Fill(etaCMNCor,leftChannel);
		if(fidRegionIndex<(Int_t)hChargeVsFidX.size()&&fidRegionIndex>=0){
			hChargeVsFidX[fidRegionIndex]->Fill(charge,fiducialValueX);
			hChargeVsFidY[fidRegionIndex]->Fill(charge,fiducialValueY);
		}
		else{
			cout<<"fidRegion not valid: "<<fidRegionIndex<<" "<<fiducialValueX<<"/"<<fiducialValueY<<endl;
			settings->getSelectionFidCuts()->Print();
		}
		hChargeVsFidCut->Fill(fiducialValueX,fiducialValueY,charge);
		hFidCutXvsChannelPos->Fill(fiducialValueX,pos);
		histoLandauDistribution2D_unmasked->Fill(charge,pos);
		bool isBorderSeedCluster = settings->hasBorderSeed(TPlaneProperties::getDetDiamond(),cluster);
		bool isBorderHitCluster = settings->hasBorderHit(TPlaneProperties::getDetDiamond(),cluster);
		if (!isBorderSeedCluster)
			histoLandauDistribution2DNoBorderSeed_unmasked->Fill(charge,pos);
		if (!isBorderHitCluster)
			histoLandauDistribution2DNoBorderHit_unmasked->Fill(charge,pos);
		bool isMaskedCluster = settings->isMaskedCluster(TPlaneProperties::getDetDiamond(),cluster,false);
		if(verbosity && pos <= 20) cout<<"masked:"<< isMaskedCluster;
		if(!isMaskedCluster){
			histoLandauDistribution2D->Fill(charge,pos);
			if (!isBorderSeedCluster)
				histoLandauDistribution2DNoBorderSeed->Fill(charge,pos);
			if (!isBorderHitCluster)
				histoLandauDistribution2DNoBorderHit->Fill(charge,pos);
		}
		if(verbosity && pos <= 20) cout<<endl;
	}
}


void TAnalysisOfSelection::initDividedAreaAxis(TAxis* axis){
    if (!axis) return;
    Int_t bins = xDivisions*yDivisions;
    if (axis->GetNbins() != bins ){
        cerr<<"divided area axis has wrong number of bins: "<<axis->GetNbins()<<"/"<<bins<<endl;
        return;
    }
    for (UInt_t y = 0; y < yDivisions; y++)
    for (UInt_t x = 0; x < xDivisions; x++){
        axis->SetBinLabel(x+xDivisions*y+1, TTransparentAnalysis::GetNameOfArea(x,y));
    }
}


void TAnalysisOfSelection::initPHvsEventNoAreaPlots(UInt_t nStart, UInt_t nEnd) {
    if(verbosity) cout<<"initPHvsEventNoAreaPlots"<<endl;
    Int_t nentriesPerBin = 20000;
//    Int_t subjectDetector = TPlaneProperties::getDetDiamond();
    TString name = TString::Format("hPHvsEventNoArea");
    TString title = TString::Format("ph vs eventNo");
    UInt_t nBins = (nEnd-nStart)/nentriesPerBin;
    if((nEnd-nStart)%nentriesPerBin!=0)nBins++;
    if (nBins==0)nBins=1;
    UInt_t yBins = xDivisions*yDivisions;
    TProfile2D* prof = new TProfile2D(name,title,nBins,nStart,nEnd,yBins,0,yBins);
    prof->Draw();
    prof->GetXaxis()->SetTitle("event no");
    initDividedAreaAxis(prof->GetYaxis());
    prof->GetZaxis()->SetTitle(TString::Format("avrg. pulse height"));
    hPHVsEventNo_Areas = prof;
    name = TString::Format("hPHvsEventNo2HighestArea");
    title = TString::Format("ph vs eventNo 2Highest");
    prof = new TProfile2D(name,title,nBins,nStart,nEnd,yBins,0,yBins);
    prof->GetXaxis()->SetTitle("event no");
    initDividedAreaAxis(prof->GetYaxis());
    prof->GetZaxis()->SetTitle(TString::Format("avrg. pulse height 2 higehest"));
    hPH2HighestVsEventNo_Areas = prof;
}

void TAnalysisOfSelection::fillPHvsEventNoAreaPlots(UInt_t area, UInt_t charge, UInt_t chargeOfTwo) {
    hPHVsEventNo_Areas->Fill(nEvent,area,charge);
    hPH2HighestVsEventNo_Areas->Fill(nEvent,area,chargeOfTwo);
}



void TAnalysisOfSelection::savePHvsEventNoAreaPlots(HistogrammSaver* histSaver,TSettings* settings,TProfile2D * prof2d,UInt_t xDivisions, UInt_t yDivisions) {
    cout<<"[TAnalysisOfSelection::savePHvsEventNoAreaPlots] "<<endl;
    if (!prof2d)return;
    if(xDivisions == 0 || yDivisions ==0 )
        return
    prof2d->Draw();
    TProfile* prof = histSaver->GetProfileX(prof2d);
    TF1* pol1Fit = new TF1("pol1Fit","pol1",prof2d->GetXaxis()->GetXmin(),prof2d->GetXaxis()->GetXmax());
    pol1Fit->SetLineWidth(1);
    pol1Fit->SetLineStyle(2);
    histSaver->Save1DProfileWithFitAndInfluence(prof,(TF1*)pol1Fit->Clone(),true);
    histSaver->SaveProfile2DWithEntriesAsText(prof2d,true);//,false);
    vector<TProfile*> vecStack;
    Int_t bins = prof2d->GetXaxis()->GetNbins();
    Float_t min = prof2d->GetXaxis()->GetXmin();
    Float_t max = prof2d->GetXaxis()->GetXmax();
    TProfile* hPHVsEventNo_AreasX [xDivisions];
    TProfile* hPHVsEventNo_AreasY [yDivisions];

    //create
    for (UInt_t x = 0; x < xDivisions; x++){
        TString name = prof2d->GetName()+(TString)"_X_"+TTransparentAnalysis::GetNameOfArea(x,-1);
        TString title = (TString)"avrg ph vs event no, X "+TTransparentAnalysis::GetNameOfArea(x,-1);
        hPHVsEventNo_AreasX[x] = new TProfile(name,title,bins,min,max);
    }
    for(UInt_t y = 0; y< yDivisions; y++){
        TString name = prof2d->GetName()+(TString)"_Y_"+TTransparentAnalysis::GetNameOfArea(-1,y);
        TString title = (TString)"avrg ph vs event no, Y "+TTransparentAnalysis::GetNameOfArea(-1,y);
        hPHVsEventNo_AreasY[y] = new TProfile(name,title,bins,min,max);
    }

    //fill
    for(int y = 0; y< prof2d->GetYaxis()->GetNbins();y++){
        TString name = prof2d->GetName()+(TString)"_"+TTransparentAnalysis::GetNameOfArea(y%xDivisions,y/xDivisions);
        prof = histSaver->GetProfileX(prof2d,name,y+1,y+1);
        prof->Draw();
        prof->SetLineColor(settings->GetColor(y));
        prof->SetMarkerColor(settings->GetColor(y));
        prof->Sumw2();
        TF1* fit = (TF1*) pol1Fit->Clone(prof->GetName()+(TString)"_fit");
        fit->SetLineColor(settings->GetColor(y));
        prof->SetTitle((TString)"Ph vs EventNo, " + TTransparentAnalysis::GetNameOfArea(y%xDivisions,y/xDivisions));
        cout<<"Save: "<<prof->GetName()<<" "<<y<<" "<<prof->GetEntries()<< endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit);
        TString title = (TString)"Ph vs EventNo, ";
        title.Append(TString::Format(", slope: %6.1f adc/1M, X",fit->GetParameter(1)*1.e6));
        title.Append(TTransparentAnalysis::GetNameOfArea(y%xDivisions,y/xDivisions));
        prof->SetTitle(title);
        cout<<"save: "<<prof->GetName()<<endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        vecStack.push_back(prof);
        TProfile* profX = hPHVsEventNo_AreasX[y%xDivisions];
        TProfile* profY = hPHVsEventNo_AreasY[y/xDivisions];
        cout<<profX->GetNbinsX()<<"/"<<prof->GetNbinsX()<<" "<<profX->GetNbinsY()<<"/"<<prof->GetNbinsY()<<" "<<profX->GetNbinsZ()<<"/"<<prof->GetNbinsZ()<<" "<<endl;
        profX->Add(prof);
        profY->Add(prof);
        cout<<"hPHVsEventNo_AreasX: "<<y%xDivisions<<" "<<profX->GetEntries()<<endl;
        cout<<"hPHVsEventNo_AreasY: "<<y/xDivisions<<" "<<profY->GetEntries()<<endl;
    }

    TString name = (TString)"hStack"+prof2d->GetName() + (TString)"AllAreasX";
    THStack* stack = new THStack(name,"Ph vs EventNo, All Areas X");
    for (UInt_t x = 0; x < xDivisions; x++){
        TF1* fit= (TF1*)pol1Fit->Clone();
        fit->SetLineColor(settings->GetColor(x));
        TProfile* prof = hPHVsEventNo_AreasX[x];
        cout<<"save: "<<prof->GetName()<<endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        TString title = (TString)"Ph vs EventNo, ";
        title.Append(TString::Format(", slope: %6.1f adc/1M, X",fit->GetParameter(1)*1.e6));
        title.Append(TTransparentAnalysis::GetNameOfArea(x,-1));
        prof->SetTitle(title);
        cout<<"save: "<<prof->GetName()<<endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        prof->SetLineColor(settings->GetColor(x));
        prof->SetMarkerColor(settings->GetColor(x));
        if(stack)stack->Add((TProfile*)prof->Clone());
    }
    histSaver->SaveStack(stack,"nostack",true);
    delete stack;
    stack = 0;

    name = (TString)"hStack"+prof2d->GetName() + (TString)"AllAreasY";
    stack = new THStack(name,"Ph vs EventNo, All Areas Y");
    for(UInt_t y = 0; y< yDivisions; y++){
        TF1* fit= (TF1*)pol1Fit->Clone();
        fit->SetLineColor(settings->GetColor(y));
        TProfile* prof = hPHVsEventNo_AreasY[y];
        cout<<"save: "<<prof->GetName()<<endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        TString title = (TString)"Ph vs EventNo, ";
        title.Append(TString::Format(", slope: %6.1f adc/1M, y",fit->GetParameter(1)*1.e6));
        title.Append(TTransparentAnalysis::GetNameOfArea(-1,y));
        prof->SetTitle(title);
        cout<<"save: "<<prof->GetName()<<endl;
        histSaver->Save1DProfileWithFitAndInfluence(prof,fit,true);
        prof->SetLineColor(settings->GetColor(y));
        prof->SetMarkerColor(settings->GetColor(y));
        if(stack)stack->Add((TProfile*)prof->Clone());
    }
    histSaver->SaveStack(stack,"nostack",true);
    delete stack;
    stack = 0;

    name = (TString)"hStack"+prof2d->GetName() + (TString)"AllAreas";
    stack = new THStack(name,"Ph vs EventNo, All Areas");
    bool foundOneHisto = false;
    for(UInt_t i = 0; i< vecStack.size();i++){
        if(!vecStack[i]) continue;
        foundOneHisto = true;
        vecStack[i]->SetLineColor(settings->GetColor(i));
        if(stack)stack->Add(vecStack[i]);
    }
    if(stack &&foundOneHisto){
        stack->Draw();
        if(stack->GetXaxis())stack->GetXaxis()->SetTitle("Event No.");
        cout<<"Save: "<<stack->GetName()<<endl;
        Float_t min = stack->GetMinimum("nostack");
        Float_t max = stack->GetMaximum("nostack");
        if(stack->GetYaxis()){
            stack->GetYaxis()->SetTitle("avrg PH");
            stack->GetYaxis()->SetRangeUser(min-.1*(max-min),max+.1*(max-min));
        }
        histSaver->SaveStack(stack,"nostack",true);
    }
    if(stack)delete stack;
    for(UInt_t i = 0; i< vecStack.size();i++){
        delete vecStack[i];
        vecStack[i]=0;
    }
//    char t; cin>>t;
}
