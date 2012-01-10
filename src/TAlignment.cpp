/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(TSettings* settings) {
	cout<<"\n\n\n**********************************************************"<<endl;
	cout<<"*************TAlignment::TAlignment***********************"<<endl;
	cout<<"**********************************************************"<<endl;

	// TODO Auto-generated constructor stub
	sys = gSystem;
	setSettings(settings);
	runNumber=settings->getRunNumber();
	cout<<runNumber<<endl;
	stringstream  runString;
	runString.str("");
	runString<<settings->getRunNumber();

	sys->MakeDirectory(runString.str().c_str());
	gErrorIgnoreLevel=1001;
	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<sys->pwd();
	filepath<<"/selectionData."<<runNumber<<".root";
//	cout<<"currentPath: "<<sys->pwd()<<endl;
//	cout<<filepath.str()<<endl;
	cout<<"OPEN eventReader with file \""<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str());
	//eventReader->checkADC();
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("alignment");
	sys->cd("alignment");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	initialiseHistos();
	cout<<"end initialise"<<endl;
	alignmentPercentage=0.1;
	Float_t stripSize=1.;// 50./10000.;//mu m
	detectorD0Z = 0.725/stripSize; // by definition in cm
	detectorD1Z = 1.625/stripSize; // by definition in cm
	detectorD2Z = 18.725/stripSize; // by definition in cm
	detectorD3Z = 19.625/stripSize; // by definition in cm
	detectorDiaZ = 10.2/stripSize; // by definition in cm
	this->runNumber=runNumber;
	verbosity=2;
	nAlignSteps=1;
	res_keep_factor = 1000;//todo anpassen
	align=NULL;
	myTrack=NULL;

}

TAlignment::~TAlignment() {
	// TODO Auto-generated destructor stub
	cout<<"TAlignment deconstructor"<<endl;
	saveHistos();
	delete histSaver;
	delete eventReader;
	sys->cd("..");
}

void TAlignment::setSettings(TSettings* settings){
	this->settings=settings;
}

void TAlignment::createVectors(){
	createVectors(eventReader->GetEntries());
}
void TAlignment::createVectors(UInt_t nEvents){
	int noHitDet=0;
	int falseClusterSizeDet=0;
	int noHitDia=0;
	int falseClusterSizeDia=0;
	int nCandidates=0;
	int nScreened=0;
	cout<<"ANALYSE VECTORS...."<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
		//check if every plane has exactly one cluster
		if(!eventReader->isValidTrack()){
			noHitDet++;
			continue;
		}
		if (eventReader->isDetMasked()){
			nScreened++;
			continue;
		}
		if(eventReader->getNDiamondClusters()!=1){
			falseClusterSizeDia++;
			continue;
		}
		nCandidates++;
		this->addEventToTracks();
	}
	cout<<"\n\nDetAnalysed "<<nEvents<<": "<<nEvents-noHitDet-falseClusterSizeDet<<" Candidates while "<< noHitDet<<" Events have not exactly one Cluster and "<<falseClusterSizeDet<<" Events have wrong cluster size"<<endl;
	cout<<"\n\nDiaAnalysed "<<nEvents-noHitDet-falseClusterSizeDet<<": "<<nCandidates<<" Candidates while "<< noHitDia<<" Events have not exactly one Cluster and "<<falseClusterSizeDia<<" Events have wrong cluster size"<<endl;
	cout<<"EVENTS SCREENED:"<<nScreened<<endl;
	cout<<tracks.size()<<" "<<tracks_masked.size()<<" "<<tracks_fidcut.size()<<" "<<tracks_masked_fidcut.size()<<endl;

//	for(UInt_t trackNo=0;trackNo<tracks.size();trackNo++){
//		Float_t xPos=tracks.at(trackNo).GetDetectorHitPosition(0);
//		Float_t yPos=tracks.at(trackNo).GetDetectorHitPosition(1);
//		cout<<trackNo<<" "<<tracks.at(trackNo).GetEventNumber()<<flush;
//		for(int no=1;no<4;no++){
//			Float_t deltaX = tracks.at(trackNo).GetDetectorHitPosition(no*2)-xPos;
//			cout<<" "<<deltaX;
//			hXPositionDifference[no-1]->Fill(deltaX);
//			hXXPositionDifference[no-1]->Fill(deltaX,tracks.at(trackNo).GetDetectorHitPosition(no*2));
//			hXYPositionDifference[no-1]->Fill(deltaX,tracks.at(trackNo).GetDetectorHitPosition(no*2+1));
//		}
//		cout<<"\ty:  ";
//		for(int no=1;no<4;no++){
//			Float_t deltaY = tracks.at(trackNo).GetDetectorHitPosition(no*2+1)-yPos;
//			cout<<" "<<deltaY;
//			hYPositionDifference[no-1]->Fill(deltaY);
//			hYXPositionDifference[no-1]->Fill(deltaY,tracks.at(trackNo).GetDetectorHitPosition(no*2));
//			hYYPositionDifference[no-1]->Fill(deltaY,tracks.at(trackNo).GetDetectorHitPosition(no*2+1));
//		}
//		cout<<endl;
//	}
}

void TAlignment::initialiseHistos(){
	for(int no = 0;no <4;no++){
		stringstream histName;
		histName<<"ScatterPlot_VectorCreation_8Hits_Plane_"<<no;
		this->hScatterPosition[no] = new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,255,256,0,255);
	}
	for(int no=0;no<3;no++){
		//*******************XPOS******************************************
		stringstream histName;
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no+1)*2);
		this->hXPositionDifference[no]=new TH1F(histName.str().c_str(),histName.str().c_str(),2048,-256,256);

		histName.str("");
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no)*2+2)<<"_vsX";
		this->hXXPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		histName.str("");
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no)*2+2)<<"_vsY";
		this->hXYPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		//*********************YPOS***********************************************

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3);
		this->hYPositionDifference[no]=new TH1F(histName.str().c_str(),histName.str().c_str(),2048,-256,256);

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3)<<"_vsX";
		this->hYXPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3)<<"_vsY";
		this->hYYPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

	}
}

void TAlignment::saveHistos(){

	for(int no=0;no<4;no++){
		this->histSaver->SaveHistogramPNG(this->hScatterPosition[no]);
		delete this->hScatterPosition[no];
	}
	for(int no=0;no<3;no++){
		//** XPOS***
		float mean = hXPositionDifference[no]->GetMean();
		float sigma= hXPositionDifference[no]->GetRMS();
		hXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		this->histSaver->SaveHistogramPNG(hXPositionDifference[no]);
		delete hXPositionDifference[no];
		hXXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXXPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		hXYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		this->histSaver->SaveHistogramPNG(hXXPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hXYPositionDifference[no]);
		delete this->hXXPositionDifference[no];
		delete this->hXYPositionDifference[no];

		//** YPOS***

		mean = hYPositionDifference[no]->GetMean();
		sigma= hYPositionDifference[no]->GetRMS();
		hYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hYXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hYYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		hYYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		this->histSaver->SaveHistogramPNG(hYPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hYXPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hYYPositionDifference[no]);
		delete hYPositionDifference[no];
		delete this->hYXPositionDifference[no];
		delete this->hYYPositionDifference[no];
	}
}

void TAlignment::addEventToTracks()
{
	TDetectorPlane D0;
	TDetectorPlane D1;
	TDetectorPlane D2;
	TDetectorPlane D3;
	TDetectorPlane Dia;
	D0.SetZ(detectorD0Z);
	D1.SetZ(detectorD1Z);
	D2.SetZ(detectorD2Z);
	D3.SetZ(detectorD3Z);
	Dia.SetZ(detectorDiaZ);
//	for(int det=0;det<eventReader->getCluster()->size();det++)
//		cout<<eventReader->getCluster()->at(det).size();
//	cout<<endl;
	D0.SetX(eventReader->getCluster()->at(0).at(0).getPosition());
	D1.SetX(eventReader->getCluster()->at(2).at(0).getPosition());
	D2.SetX(eventReader->getCluster()->at(4).at(0).getPosition());
	D3.SetX(eventReader->getCluster()->at(4).at(0).getPosition());

	D0.SetY(eventReader->getCluster()->at(1).at(0).getPosition());
	D1.SetY(eventReader->getCluster()->at(3).at(0).getPosition());
	D2.SetY(eventReader->getCluster()->at(5).at(0).getPosition());
	D3.SetY(eventReader->getCluster()->at(7).at(0).getPosition());

	Dia.SetX(eventReader->getCluster()->at(8).at(0).getPosition());

	TDiamondTrack newDiamondTrack(nEvent,D0,D1,D2,D3,Dia);

	Float_t xPos=eventReader->getCluster()->at(0).at(0).getPosition();
	Float_t yPos=eventReader->getCluster()->at(1).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(0,xPos);
	newDiamondTrack.SetDetectorHitPosition(1,yPos);
	hScatterPosition[0]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(2).at(0).getPosition();
	yPos=eventReader->getCluster()->at(3).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(2,xPos);
	newDiamondTrack.SetDetectorHitPosition(3,yPos);
	hScatterPosition[1]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(4).at(0).getPosition();
	yPos=eventReader->getCluster()->at(5).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(4,xPos);
	newDiamondTrack.SetDetectorHitPosition(5,yPos);
	hScatterPosition[2]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(6).at(0).getPosition();
	yPos=eventReader->getCluster()->at(7).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(6,xPos);
	newDiamondTrack.SetDetectorHitPosition(7,yPos);
	hScatterPosition[3]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(8).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(8,xPos);

	tracks.push_back(newDiamondTrack);
	tracks_masked.push_back((rand.Rndm()>this->alignmentPercentage));
	tracks_fidcut.push_back(newDiamondTrack);
	tracks_masked_fidcut.push_back((rand.Rndm()>this->alignmentPercentage));

}


int TAlignment::Align()
{
	if(tracks.size()==0)
		createVectors();

	if (verbosity){
		cout<<"\n\n\nTAlignment::Align:Starting \""<<histSaver->GetPlotsPath()<<"\""<<endl;
		cout<<"\t\t"<<tracks.size()<<" "<<tracks_masked.size()<< " ";
		cout<<		  tracks_fidcut.size()<<" "<<tracks_masked_fidcut.size()<<endl;
		cout << "\t\t "<<eventReader<<" ."<<endl;
	}
	if(tracks.size()==0) {
		cout<<"TAlignment::Align: No tracks found; need to CallClustering::ClusterRun first..."<<endl;
		return -1;
		//ClusterRun(plots); // doesn't use alternative clustering
	}

	if (tracks.size() == 0) {
		cout << "TAlignment::Align:No tracks available. Alignment not possible. (tracks.size() = " << tracks.size() << ")" << endl;
		return 0;
	}

	if(align==NULL){
		align = new TDetectorAlignment();//histSaver->GetPlotsPath(), tracks, tracks_masked);
		cout<<"TAlignment::Align::Detectoralignment did not exist, so created new DetectorAlignment"<<endl;
		align->SetZOffset(0,detectorD0Z);
		align->SetZOffset(1,detectorD1Z);
		align->SetZOffset(2,detectorD2Z);
		align->SetZOffset(3,detectorD3Z);
		align->SetZOffset(4,detectorDiaZ);
	}
	if(myTrack==NULL){
		myTrack = new TTrack(*align);
		cout<<"TAlignment::Align::created new TTrack";
	}

	AlignDetector(TPlane::XY_COR,1,0,0,true);
	AlignDetector(TPlane::XY_COR,2,0,0,true);
	AlignDetector(TPlane::XY_COR,3,0,0,true);
	myTrack->setDetectorAlignment(*align);
	cout<<endl;
	for(UInt_t nSteps=0;nSteps<5;nSteps++){
		AlignDetector(TPlane::XY_COR,3,0,2);
		myTrack->setDetectorAlignment(*align);
		AlignDetector(TPlane::XY_COR,1,0,2);
		myTrack->setDetectorAlignment(*align);
		AlignDetector(TPlane::XY_COR,2,0,3);
		myTrack->setDetectorAlignment(*align);
		cout<<endl;
	}
	AlignDetector(TPlane::XY_COR,3,0,2,true);
	AlignDetector(TPlane::XY_COR,1,0,2,true);
	AlignDetector(TPlane::XY_COR,2,0,3,true);

	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
//	for (int alignStep = 0; alignStep < nAlignSteps; alignStep++) {
//		//TDetectorAlignment* align = new TDetectorAlignment(plots_path, tracks, tracks_mask);
//		align->LoadTracks(tracks,tracks_masked);
//		if (verbosity) cout<<"TAlignment::Align:start with alignmentStep no. "<<alignStep+1 << " of " <<nAlignSteps<<endl;
//		doDetAlignmentStep();
//		cout << "AlignmentClass::Align:Intrinsic silicon resolution " << align->GetSiResolution() << " strips or " << align->GetSiResolution() * 50 << "um" << endl;
//		doDiaAlignmentStep();
//	}
	return 1;
}


void TAlignment::AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::XY_COR,subjectPlane,refPlane1,refPlane2);
}



void TAlignment::AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::X_COR,subjectPlane,refPlane1,refPlane2);
}



void TAlignment::AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::Y_COR,subjectPlane,refPlane1,refPlane2);
}

void TAlignment::AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot)
{
	cout<<"TAlignment::AlignDetector\n\t align "<<TPlane::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with Plane "<<refPlane1<<" and "<<refPlane2<<endl;
	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
		return;
	}

	TResidual res = getResidual(cor,subjectPlane,refPlane1,refPlane2);
	if(refPlane1==refPlane2){
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR)
			align->AddToXOffset(subjectPlane,res.resXMean);
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR)
			align->AddToYOffset(subjectPlane,res.resYMean);
		return;
	}
	else{
		Float_t variableDif = (res.nUsedTracks * res.sumV2y - res.sumVy * res.sumVy);
		Float_t y_offset = (res.sumRy * res.sumV2y - res.sumVRy * res.sumVy) / variableDif;
		//phiy_offset = (observedX.size() * sumvr - sumr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
		Float_t phiy_offset = TMath::ATan((res.nUsedTracks* res.sumVRy - res.sumRy * res.sumVy) / (res.nUsedTracks * res.sumV2y - res.sumVy * res.sumVy));
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR){
			align->AddToYOffset(subjectPlane,y_offset);
			align->AddToPhiYOffset(subjectPlane,phiy_offset);
			cout<<"/tY:"<<y_offset<<" "<<phiy_offset<<endl;
		}

		variableDif= (res.nUsedTracks * res.sumV2x - res.sumVx * res.sumVx);
		Float_t x_offset = (res.sumRx * res.sumV2x - res.sumVRx * res.sumVx) / variableDif;
		//phix_offset = -(observedX.size() * sumvr - sumr * sumv) / (observedX.size() * sumv2 - sumv * sumv);
		Float_t phix_offset = TMath::ATan(-(res.nUsedTracks * res.sumVRx - res.sumRx * res.sumVx) / (res.nUsedTracks * res.sumV2x - res.sumVx * res.sumVx));
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR){
			align->AddToXOffset(subjectPlane,x_offset);
			align->AddToPhiXOffset(subjectPlane,phix_offset);
			cout<<"/tX:"<<x_offset<<" "<<phix_offset<<endl;
		}
	}
}

/**
 * @brief:
 */
TResidual TAlignment::getResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool plot){
	if(verbosity)cout<<"TAlignment::getResidual"<<subjectPlane<<" "<<refPlane1<<" "<<refPlane2<<" "<<plot<<" with "<<alignmentPercentage<<endl;
	vector<Float_t> vecXPred;
	vector<Float_t> vecYPred;
	vector<Float_t> vecXObs;
	vector<Float_t> vecYObs;
	vector<Float_t> vecDeltaX;
	vector<Float_t> vecDeltaY;

	vector<UInt_t> vecRefPlanes;

	vecRefPlanes.push_back(refPlane1);
	if(refPlane1!=refPlane2)
		vecRefPlanes.push_back(refPlane2);

	for(Int_t nEvent=0; nEvent<eventReader->GetEntries(); nEvent++)
	{
		TRawEventSaver::showStatusBar(nEvent,eventReader->GetEntries());
		if (rand.Rndm()>this->alignmentPercentage) continue; //skip randomly to use only a part alignmentPercentage for alignment //todo improve method!
		eventReader->LoadEvent(nEvent);
		if (eventReader->getEvent()->isMasked()) continue; // skip tracks not selected for determining alignment constants
		if (!eventReader->getEvent()->isValidSiliconEvent()) continue;

		myTrack->setEvent(eventReader->getEvent());
		//if(verbosity>3)	cout<<myTrack->getEvent()->getEventNumber()<<" "<<myTrack->getPosition(cor,subjectPlane)<<" "<<eventReader->getEvent()->getPlane(subjectPlane).getXPosition(0)<<endl;;
		Float_t xPositionObserved  = myTrack->getPosition(TPlane::X_COR,subjectPlane);
		Float_t yPositionObserved  = myTrack->getPosition(TPlane::Y_COR,subjectPlane);
		vecXObs.push_back(xPositionObserved);
		vecYObs.push_back(yPositionObserved);
		TPositionPrediction predictedPostion = myTrack->predictPosition(subjectPlane,vecRefPlanes);
		if(verbosity>3)	predictedPostion.Print();
		if(verbosity>3)	cout<<xPositionObserved<<" / "<<yPositionObserved<<endl;
		Float_t deltaX = xPositionObserved-predictedPostion.getPositionX();//X_OBS-X_Pred
		Float_t deltaY = yPositionObserved-predictedPostion.getPositionY();//Y_OBS-Y_Pred
		vecDeltaX.push_back(deltaX);
		vecDeltaY.push_back(deltaY);
		vecXPred.push_back(predictedPostion.getPositionX());
		vecYPred.push_back(predictedPostion.getPositionY());
		if(verbosity>3)	cout<< deltaX<<" "<<deltaY<<endl;
	}

	if(verbosity>2)cout<<vecDeltaX.size()<<" "<<vecDeltaY.size()<<" "<< vecXPred.size()<<" "<<vecYPred.size()<<endl;
	//first estimate residuals widths
	TResidual residuum = calculateResidual(cor,vecXPred,vecDeltaX,vecYPred,vecDeltaY);
	if(plot){
		stringstream histName;
		histName<<"hScatterPlot_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2<<"_alignment";
		histSaver->SaveHistogramPNG((TH2F*)histSaver->CreateScatterHisto(histName.str(),vecXObs,vecYObs).Clone());
		histName.str("");
		histName<<"hDistributionPlot_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2<<"DeltaY_alignment";
		//histSaver->SaveHistogramPNG((TH2F*)histSaver->CreateDistributionHisto(histName.str(),vecYObs,vecDeltaX).Clone());
		histName.str("");
		histName<<"hDistributionPlot_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2<<"DeltaY_alignment";
	//	histSaver->SaveHistogramPNG((TH2F*)histSaver->CreateDistributionHisto(histName.str(),vecYObs,vecDeltaX).Clone());
	}
	return residuum;
}


/**
 * calculateResidual if there is no residuals calculatet to cut on the res_keep_factor;
 * this funcition opens the other calculateResidual function but uses a TResidual res
 * which is using all items...
 */
TResidual TAlignment::calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY){
	TResidual res;
	res.resXMean=0;
	res.resXSigma=10000;
	res.resYMean=0;
	res.resYSigma=10000;
	if(verbosity>2)cout<<"\t"<<deltaX.size()<<" "<<deltaY.size()<<" "<< xPred.size()<<" "<<yPred.size()<<endl;
	return calculateResidual(cor,xPred,deltaX,yPred,deltaY,res);
}

/**
 * calculating residual:
 * uses the input TResidual res to use only events with
 * resxtest= TMath::Abs(deltaX.at(i)-res.resXMean)/res.resXSigm < res_keep_factor
 * resytest= TMath::Abs(deltaY.at(i)-res.resYMean)/res.resYSigma < res_keep_factor
 *
 */
TResidual TAlignment::calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t>yPred,vector<Float_t> deltaY,TResidual res)
{
	Float_t resxmean = 0;
	Float_t resXsigma = 0;
	Float_t resymean = 0;
	Float_t resYsigma = 0;
	Float_t resxtest;
	Float_t resytest;
	UInt_t nUsedTracks=0;
	Float_t sumRx=0;
	Float_t sumRy=0;
	Float_t sumVx=0;
	Float_t sumVy=0;
	Float_t sumV2x=0;
	Float_t sumV2y=0;
	Float_t sumVRx=0;
	Float_t sumVRy=0;
	if(verbosity>2)cout<<"\tcalculate Residuum "<<res_keep_factor<<endl;
	if(verbosity>2)cout<<"\t"<<deltaX.size()<<" "<<deltaY.size()<<" "<< xPred.size()<<" "<<yPred.size()<<endl;
	for(UInt_t i=0;i<deltaX.size();i++){
		resxtest= TMath::Abs(deltaX.at(i)-res.resXMean)/res.resXSigma;
		resytest= TMath::Abs(deltaY.at(i)-res.resYMean)/res.resYSigma;
		bool used=false;
		//only add if restest is smaller than res_keep_factor
		if((cor==TPlane::X_COR)&&resxtest<res_keep_factor){
			resxmean += deltaX.at(i);
			resXsigma +=deltaX.at(i)*deltaX.at(i);
			used=true;
		}//end if
		else if((cor==TPlane::X_COR)&&resytest<res_keep_factor){
			resymean += deltaY.at(i);
			resYsigma += deltaY.at(i)*deltaY.at(i);
			used=true;
		}//end else if
		else if((cor==TPlane::XY_COR)&&resxtest<res_keep_factor&&resytest<res_keep_factor){
			resymean += deltaY.at(i);
			resYsigma += deltaY.at(i)*deltaY.at(i);
			resxmean += deltaX.at(i);
			resXsigma +=deltaX.at(i)*deltaX.at(i);
			used=true;
		}//end else if

		//if track can be used (restest<res_keep_factor), add values
		if(used){
			sumRx+=deltaX.at(i);
			sumRy+=deltaY.at(i);
			sumVx+=yPred.at(i);
			sumVy+=xPred.at(i);
			sumV2x+=yPred.at(i)*yPred.at(i);
			sumV2y+=xPred.at(i)*xPred.at(i);
			sumVRx+=yPred.at(i)*deltaX.at(i);
			sumVRy+=xPred.at(i)*deltaY.at(i);
			nUsedTracks++;
		}//end if
	}//end for loop
	resxmean = resxmean /(Double_t)nUsedTracks;
	resymean = resymean /(Double_t)nUsedTracks;
	resXsigma = TMath::Sqrt(resXsigma / (Double_t)nUsedTracks - resxmean*resxmean);
	resYsigma = TMath::Sqrt(resYsigma / (Double_t)nUsedTracks- resymean*resymean);

	if(verbosity>0)	cout<<"\tused "<<nUsedTracks<<" Tracks"<<endl;
	if(verbosity>0)	cout<<"\tX: "<<std::setprecision(4)<<resxmean<<"+/-"<<resXsigma<<endl;
	if(verbosity>0)	cout<<"\tY: "<<resymean<<"+/-"<<resYsigma<<endl;
	//set values
	TResidual residual;
	residual.resXMean=resxmean;
	residual.resXSigma=resXsigma;
	residual.resYMean=resymean;
	residual.resYSigma=resYsigma;
	residual.sumRx=sumRx;
	residual.sumRy=sumRy;
	residual.sumV2x=sumV2x;
	residual.sumV2y=sumV2y;
	residual.sumVRx=sumVRx;
	residual.sumVRy=sumVRy;
	residual.sumVx =sumVx;
	residual.sumVy = sumVy;
	residual.nUsedTracks=nUsedTracks;

	return residual;
}


void TAlignment::PrintEvents(UInt_t maxEvent){
	cout<<"\n\n\n\n\n\nPrintEVENTS"<<maxEvent<<"\n\n\n"<<flush;
	for(UInt_t event=0;event<maxEvent;event++){
		eventReader->LoadEvent(event);
		cout<<event<<":\t"<<flush;
		eventReader->getEvent()->Print(1);
	}
	cout<<"\n\n\n\n\n\n"<<flush;
	sleep(10);
}
