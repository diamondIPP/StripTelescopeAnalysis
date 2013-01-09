/*
 * TAnalysisOfAlignment.cpp
 *
 *  Created on: Mar 2, 2012
 *      Author: bachmair
 */

#include "../include/TAnalysisOfAlignment.hh"

TAnalysisOfAlignment::TAnalysisOfAlignment(TSettings *settings) {
	cout<<"\n\n\n\n**********************************************************"<<endl;
		cout<<"**********************************************************"<<endl;
		cout<<"*********TAnalysisOfClustering::TAnalysisOfAlignment*****"<<endl;
		cout<<"**********************************************************"<<endl;
		cout<<"**********************************************************\n\n\n"<<endl;
		setSettings(settings);
		UInt_t runNumber=settings->getRunNumber();
		sys = gSystem;
		settings->goToAlignmentRootDir();
		eventReader=new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
		histSaver=new HistogrammSaver();
		settings->goToAlignmentAnalysisDir();
		htmlAlignment=new THTMLAlignment(settings);
		htmlAlignment->setAlignment(eventReader->getAlignment());
		htmlAlignment->setFileGeneratingPath(settings->getAlignmentAnalysisFilePath());
		htmlAlignment->createContent();
		htmlAlignment->generateHTMLFile();
		settings->goToAlignmentAnalysisDir();
		sys->cd("anaAlignmnet");
		histSaver->SetPlotsPath(settings->getAlignmentAnalysisFilePath());
		histSaver->SetRunNumber(runNumber);
		settings->goToAlignmentRootDir();
		initialiseHistos();
		cout<<"end initialise"<<endl;
		settings=0;
		verbosity=0;
}

TAnalysisOfAlignment::~TAnalysisOfAlignment() {
	delete eventReader;
	delete histSaver;
	settings->goToOutputDir();
}

void TAnalysisOfAlignment::doAnalysis(UInt_t nEvents){
	if (nEvents==0||eventReader->GetEntries()<nEvents)
		nEvents=eventReader->GetEntries();
	cout<<"Do Analysis After Alignment...."<<endl;
	DoEtaCorrection(1);

//	TH1F *histo = new TH1F("hPredictedStripPositionDet2","hPredictedStripPositionDet2",128,-0.501,0.501);
//	for(nEvent=0;nEvent<nEvents;nEvent++){
//		eventReader->LoadEvent(nEvent);
//		if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//			continue;
//		TCluster xClus = eventReader->getCluster(2,0);
//		vector<UInt_t>refPlanes;
//		refPlanes.push_back(0);
//		refPlanes.push_back(2);
//		refPlanes.push_back(3);
//		TPositionPrediction* pred = eventReader->predictPosition(1,refPlanes,false);
//		TH1F *hEtaIntegral=eventReader->getEtaIntegral(2);
//		if(nEvent==0)histSaver->SaveHistogram(eventReader->getEtaIntegral(2));
//		Float_t predictedStripPosition = eventReader->getPositionInDetSystem(2,pred->getPositionX(),pred->getPositionY());
//		UInt_t stripMiddle=(UInt_t) (predictedStripPosition+0.5);
//		Float_t delta = predictedStripPosition-stripMiddle;
//		histo->Fill(delta);
//	}
//	histSaver->SaveHistogram(histo);
//	Int_t minBin = histo->GetMinimumBin();
//	UInt_t nMinEntries =histo->GetBinContent(minBin);
//	TH1F * histo2=new TH1F("hPredictedStripPositionDet2_flattend","hPredictedStripPositionDet2_flattend",128,-0.501,0.501);
//	vector<UInt_t> vecEventNo;
//	vecEventNo.clear();
//	for(nEvent=0;nEvent<nEvents;nEvent++){
//			eventReader->LoadEvent(nEvent);
//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//				continue;
//			TCluster xClus = eventReader->getCluster(2,0);
//			vector<UInt_t>refPlanes;
//			refPlanes.push_back(0);
//			refPlanes.push_back(2);
//			refPlanes.push_back(3);
//			TPositionPrediction* pred = eventReader->predictPosition(1,refPlanes,false);
//			TH1F *hEtaIntegral=eventReader->getEtaIntegral(2);
//			if(nEvent==0)histSaver->SaveHistogram(eventReader->getEtaIntegral(2));
//			Float_t predictedStripPosition = eventReader->getPositionInDetSystem(2,pred->getPositionX(),pred->getPositionY());
//			UInt_t stripMiddle=(UInt_t) (predictedStripPosition+0.5);
//			Float_t delta = predictedStripPosition-stripMiddle;
//			Int_t bin =histo2->FindBin(delta);
//			if(histo2->GetBinContent(bin)<nMinEntries){
//				vecEventNo.push_back(nEvent);
//				histo2->Fill(delta);
//			}
//	}
//	histSaver->SaveHistogram(histo2);
//
//	TH1F *hEta=new TH1F("hEtaDistributionDet2_flattend","hEtaDistribution on det 2 from flatted StripPositionDistribution",128,0,1);
//	for(UInt_t i=0;i<vecEventNo.size();i++){
//		nEvent= vecEventNo.at(i);
//		eventReader->LoadEvent(nEvent);
//
//		if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//			continue;
//		Float_t eta = eventReader->getCluster(2,0).getEta();
//		hEta->Fill(eta);
//
//	}
//	histSaver->SaveHistogram(hEta);
}
void TAnalysisOfAlignment::setSettings(TSettings *settings)
{
	this->settings=settings;
}

void TAnalysisOfAlignment::initialiseHistos()
{
  this->hChi2DistributionX=new TH1F("hChi2X","#chi^{2}_{X}-Distribution valid SiliconTrack",512,0,5);
  this->hChi2DistributionX->GetXaxis()->SetTitle("chi^{2}_{X}");
  this->hChi2DistributionX->GetYaxis()->SetTitle("Number of Entries #");
  this->hChi2DistributionY=new TH1F("hChi2Y","#chi^{2}_{Y}-Distribution valid SiliconTrack",512,0,5);
  this->hChi2DistributionY->GetXaxis()->SetTitle("chi^{2}_{Y}");
  this->hChi2DistributionY->GetYaxis()->SetTitle("Number of Entries #");
  this->hChi2DistributionXY=new TH1F("hChi2XY","#chi^{2}_{X}+#chi^{2}_{Y}-Distribution valid SiliconTrack",512,0,5);
  this->hChi2DistributionXY->GetXaxis()->SetTitle("#chi^{2}_{X}+#chi^{2}_{Y}");
  this->hChi2DistributionXY->GetYaxis()->SetTitle("Number of Entries #");
  this->hChi2DistributionXY2D=new TH2F("hChi2XY2D","#chi^{2}_{X}+#chi^{2}_{Y}-Distribution valid SiliconTrack",512,0,5,512,0,5);
  this->hChi2DistributionXY2D->GetXaxis()->SetTitle("#chi^{2}_{X}");
  this->hChi2DistributionXY2D->GetYaxis()->SetTitle("#chi^{2}_{Y}");
  this->hChi2DistributionXY->GetZaxis()->SetTitle("number of entries");
  this->hAngleDistribution=new TH2F("hAngularDistribution","hAngularDistribution",256,-5,5,256,-5,5);
  this->hAngleDistribution->GetXaxis()->SetTitle("#Phi_{X}");
  this->hAngleDistribution->GetYaxis()->SetTitle("#Phi_{Y}");
}


void TAnalysisOfAlignment::DoEtaCorrection(UInt_t correctionStep){
cout<<"****************************************"<<endl;
	cout<<"Do Eta correction for all silicon planes "<<correctionStep<< endl;
	cout<<"****************************************"<<endl;
	vector<TH1F*> histoStripDistribution,histoStripDistributionFlattned;
	vector<TH1F *>vecHEta;
//	stringstream fileName;
//	fileName<<	"etaIntegral_Step"<<correctionStep<<"."<<settings->getRunNumber()<<".root";
	TFile* correctedEtaFile = new TFile(settings->getEtaDistributionPath(correctionStep).c_str(),"RECREATE");
	cout<<"create Histos..."<<endl;
	for(UInt_t det=0; det<TPlaneProperties::getNSiliconDetectors();det++){
		stringstream histoTitle;
		histoTitle<<"hPredictedStripPosition"<<"_step"<<correctionStep<<"_"<<TPlaneProperties::getStringForDetector(det);
		histoStripDistribution.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,-0.501,0.501));
		histoTitle<<"_flattend";
		histoStripDistributionFlattned.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,-0.501,0.501));
		histoTitle.str("");
		histoTitle.clear();
		histoTitle<<"hCorrectedEtaDistribution"<<"_step"<<correctionStep<<"_"<<TPlaneProperties::getStringForDetector(det);
		vecHEta.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,0,1));
	}

	cout<<"fill first strip hit histo"<<eventReader->GetEntries()<<endl;

  for( nEvent=0;nEvent<eventReader->GetEntries();nEvent++){
		TRawEventSaver::showStatusBar(nEvent,eventReader->GetEntries());
		eventReader->LoadEvent(nEvent);
		if(!eventReader->useForAnalysis()&&!eventReader->useForAlignment())
			continue;
		for(UInt_t subjectPlane =0; subjectPlane<TPlaneProperties::getNSiliconPlanes();subjectPlane++){
//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//				continue;
			vector<UInt_t>vecRefPlanes;
			for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
				if(subjectPlane!=refPlane)vecRefPlanes.push_back(refPlane);
			TPositionPrediction* pred = eventReader->predictPosition(subjectPlane,vecRefPlanes,false);
			Float_t predictedStripPositionX = eventReader->getPositionInDetSystem(subjectPlane*2,pred->getPositionX(),pred->getPositionY());
			Float_t predictedStripPositionY = eventReader->getPositionInDetSystem(subjectPlane*2+1,pred->getPositionX(),pred->getPositionY());
			UInt_t stripMiddleX=(UInt_t) (predictedStripPositionX+0.5);
			Float_t deltaX = predictedStripPositionX-stripMiddleX;
			UInt_t stripMiddleY=(UInt_t) (predictedStripPositionY+0.5);
			Float_t deltaY = predictedStripPositionY-stripMiddleY;
//			cout<<nEvent<<": "<<subjectPlane<<"Fill "<<deltaX<<" "<<deltaY<<endl;
			histoStripDistribution.at(subjectPlane*2)->Fill(deltaX);
			histoStripDistribution.at(subjectPlane*2+1)->Fill(deltaY);
			chi2Distribution();
		}
	}
  saveHistos();
	vector<UInt_t> vecMinEntries;
	cout<<"Minimal Entries in a bin of historgram:"<<endl;
	for(UInt_t det=0;det<TPlaneProperties::getNSiliconDetectors();det++){
		TH1F* histo=histoStripDistribution.at(det);
		Int_t minBin =      histo->GetMinimumBin();
		Int_t nMinEntries = histo->GetBinContent(minBin);
		cout<<endl<< det<< ": "<<minBin<<" "<<nMinEntries<<"\t";
		vecMinEntries.push_back(nMinEntries);
	}
	cout<<"\n\n"<<endl;

	vector<UInt_t>vecEventNo[9];
	cout<<"create flattened strip hit histo "<<eventReader->GetEntries()<<endl;
	for( nEvent=0;nEvent<eventReader->GetEntries();nEvent++){
			TRawEventSaver::showStatusBar(nEvent,eventReader->GetEntries());
			eventReader->LoadEvent(nEvent);
			if(!eventReader->useForAnalysis()&&!eventReader->useForAlignment())
				continue;

		for(UInt_t subjectPlane=0; subjectPlane<TPlaneProperties::getNSiliconPlanes();subjectPlane++){
//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//				continue;
			vector<UInt_t>refPlanes;
			for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
				if(subjectPlane!=refPlane)refPlanes.push_back(refPlane);
			TPositionPrediction* pred = eventReader->predictPosition(subjectPlane,refPlanes,false);

			Float_t predictedStripPositionX = eventReader->getPositionInDetSystem(subjectPlane*2,pred->getPositionX(),pred->getPositionY());
			Float_t predictedStripPositionY = eventReader->getPositionInDetSystem(subjectPlane*2+1,pred->getPositionX(),pred->getPositionY());
			UInt_t stripMiddleX=(UInt_t) (predictedStripPositionX+0.5);
			Float_t deltaX = predictedStripPositionX-stripMiddleX;
//			histoStripDistribution.at(subjectPlane*2)->Fill(deltaX);
			UInt_t stripMiddleY=(UInt_t) (predictedStripPositionY+0.5);
			Float_t deltaY = predictedStripPositionY-stripMiddleY;
//			histoStripDistribution.at(subjectPlane*2)->Fill(deltaY);

			Int_t binX=histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaX);
			Int_t binY=histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaY);
			if(histoStripDistributionFlattned.at(subjectPlane*2)->GetBinContent(binX)<vecMinEntries.at(subjectPlane*2)){
				vecEventNo[subjectPlane*2].push_back(nEvent);
				histoStripDistributionFlattned.at(subjectPlane*2)->Fill(deltaX);
			}
			if(histoStripDistributionFlattned.at(subjectPlane*2+1)->GetBinContent(binY)<vecMinEntries.at(subjectPlane*2+1)){
				vecEventNo[subjectPlane*2+1].push_back(nEvent);
				histoStripDistributionFlattned.at(subjectPlane*2+1)->Fill(deltaY);
			}
		}
	}
	for(UInt_t det =0; det<TPlaneProperties::getNSiliconDetectors();det++){
		cout<<"save histogram: "<<det<<"  "<<histoStripDistributionFlattned.at(det)->GetTitle()<<"  "<<histoStripDistribution.at(det)->GetTitle()<<endl;
		histSaver->SaveHistogram(histoStripDistributionFlattned.at(det));
		correctedEtaFile->Add(histoStripDistributionFlattned.at(det));
		histSaver->SaveHistogram(histoStripDistribution.at(det));
		correctedEtaFile->Add(histoStripDistribution.at(det));
	}

	cout<<"\n\ncreate eta correction histo"<<endl;
	for(UInt_t det =0; det<TPlaneProperties::getNSiliconDetectors();det++){
		for(UInt_t i=0;i<vecEventNo[det].size();i++){
			nEvent= vecEventNo[det].at(i);
			eventReader->LoadEvent(nEvent);

			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
				continue;
			Float_t eta = eventReader->getCluster(det,0).getEta();
			vecHEta.at(det)->Fill(eta);

		}

			stringstream histName;
			histName<<"hEtaIntegral"<<"_step"<<correctionStep<<"_"<<TPlaneProperties::getStringForDetector(det);;
//			UInt_t nBins = vecHEta.at(det)->GetNbinsX();
			TH1F *histo= TClustering::createEtaIntegral(vecHEta.at(det),histName.str());
			cout<<"save "<<vecHEta.at(det)->GetTitle()<<" "<<vecHEta.at(det)->GetEntries()<<endl;
			histSaver->SaveHistogram(vecHEta.at(det));
			correctedEtaFile->Add(vecHEta.at(det));
			cout<<"save "<<histo->GetTitle()<<" "<<histo->GetEntries()<<endl;
			histSaver->SaveHistogram(histo);
			correctedEtaFile->Add((TH1F*)histo->Clone());
	}
	correctedEtaFile->Write();
	cout<<"Closing "<<correctedEtaFile->GetName()<<endl;
	for(UInt_t i=0;i<histoStripDistribution.size();i++)
	  cout<<histoStripDistribution.at(i)<<" "<<flush;
	cout<<endl;
	for(UInt_t i =0;i<histoStripDistribution.size();i++)
	    cout<<histoStripDistributionFlattned.at(i)<<" "<<flush;
	cout<<endl;
	for(UInt_t i =0;i<histoStripDistribution.size();i++)
	  cout<<vecHEta.at(i)<<" "<<flush;

}

void TAnalysisOfAlignment::chi2Distribution()
{
  vector<UInt_t> vecRefPlanes;
  for(UInt_t i=0;i<TPlaneProperties::getNSiliconPlanes();i++)vecRefPlanes.push_back(i);
  TPositionPrediction *prediction = eventReader->predictPosition(TPlaneProperties::getDiamondPlane(),vecRefPlanes,false);

  Float_t phiX  = prediction->getPhiX();
  Float_t phiY  = prediction->getPhiY();
  Float_t chi2X = prediction->getChi2X();
  Float_t chi2Y = prediction->getChi2Y();
//  cout<<nEvent<<"\t"<<chi2X<<" "<<chi2Y<<endl;
  this->hAngleDistribution->Fill(phiX,phiY);
  this->hChi2DistributionX->Fill(chi2X);
  this->hChi2DistributionY->Fill(chi2Y);
  this->hChi2DistributionXY->Fill(chi2X+chi2Y);
  this->hChi2DistributionXY2D->Fill(chi2X,chi2Y);
}

void TAnalysisOfAlignment::saveHistos()
{
  cout<<"Save Histos"<<endl;
  histSaver->SaveHistogram(hAngleDistribution);
  if(hAngleDistribution!=0)delete hAngleDistribution;
  histSaver->SaveHistogram(hChi2DistributionX);
  if(hChi2DistributionX!=0)delete hChi2DistributionX;
  histSaver->SaveHistogram(hChi2DistributionY);
  if(hChi2DistributionY!=0)delete hChi2DistributionY;
  histSaver->SaveHistogram(hChi2DistributionXY);
  if(hChi2DistributionXY!=0)delete hChi2DistributionXY;
  histSaver->SaveHistogram(hChi2DistributionXY2D);
  if(hChi2DistributionXY2D!=0)delete hChi2DistributionXY2D;

  cout<<"DONE"<<endl;
}




