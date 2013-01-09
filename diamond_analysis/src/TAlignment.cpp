/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(TSettings* inputSettings) {
	cout << "\n\n\n**********************************************************" << endl;
	cout << "*************TAlignment::TAlignment***********************" << endl;
	cout << "**********************************************************" << endl;

	sys = gSystem;
	setSettings(inputSettings);
	runNumber = settings->getRunNumber();
	cout << runNumber << endl;
	stringstream runString;
	settings->goToSelectionTreeDir();
	htmlAlign = new THTMLAlignment(settings);
	eventReader = new TADCEventReader(settings->getSelectionTreeFilePath(), settings->getRunNumber(),settings->getVerbosity()<15?0:settings->getVerbosity()-15);
	eventReader->setEtaDistributionPath(settings->getEtaDistributionPath());
	cout<<"Eta dist path: "<<eventReader->getEtaDistributionPath()<<endl;
	histSaver = new HistogrammSaver();
	settings->goToAlignmentDir();
	histSaver->SetPlotsPath(settings->getAlignmentDir());
	histSaver->SetRunNumber(runNumber);
	histSaver->SetNumberOfEvents(eventReader->GetEntries());
	htmlAlign->setFileGeneratingPath(settings->getAlignmentDir());
	settings->goToAlignmentRootDir();
	cout << "end initialise" << endl;
	alignmentPercentage = settings->getAlignment_training_track_fraction();
	Float_t stripSize = 1.;    // 50./10000.;//mu m
	detectorD0Z = 0.725 / stripSize;    // by definition in cm
	detectorD1Z = 1.625 / stripSize;    // by definition in cm
	detectorD2Z = 18.725 / stripSize;    // by definition in cm
	detectorD3Z = 19.625 / stripSize;    // by definition in cm
	detectorDiaZ = 10.2 / stripSize;    // by definition in cm
	verbosity = settings->getVerbosity();
	cout<<"Verbosity is: "<<verbosity<<" "<<settings->getVerbosity()<<endl;
	res_keep_factor = settings->getRes_keep_factor();
	cout << "Res Keep factor is set to " << res_keep_factor << endl;
	align = NULL;
	myTrack = NULL;
	nAlignmentStep = -1;
	nAlignSteps = settings->GetSiliconAlignmentSteps();
	nDiaAlignmentStep = -1;
	nDiaAlignSteps = settings->GetDiamondAlignmentSteps();
	
	diaCalcMode = TCluster::maxValue;    //todo
	silCalcMode = TCluster::corEta;    //todo

	bPlotAll = settings->doAllAlignmentPlots()||verbosity>6;
	results=0;

}

TAlignment::~TAlignment() {
	if(verbosity)cout << "TAlignment deconstructor" << endl;
	if (results!=0)results->setAlignment(this->align);
	htmlAlign->setAlignment(align);
	htmlAlign->createContent();
	htmlAlign->generateHTMLFile();
	this->saveAlignment();
	if(htmlAlign!=0)delete htmlAlign;
	if (myTrack) delete myTrack;
	if (histSaver) delete histSaver;
	//	if(eventReader)delete eventReader;
	settings->goToOutputDir();
	if(verbosity)cout<<"Closed TAlginment"<<endl;
}

void TAlignment::setSettings(TSettings* settings) {
	if(settings==0)
	{
		cerr<<"Input settings == 0. BREAK!"<<endl;
		exit(-1);
	}
	this->settings = settings;
}

/**
 * initialise the variable align and set the Z Offsets
 */
void TAlignment::initialiseDetectorAlignment() {
	if (align == NULL) {
		align = new TDetectorAlignment();
		align->setVerbosity(verbosity-4>0?verbosity-4:0);
		htmlAlign->setAlignment(align);
		cout << "TAlignment::Align::Detectoralignment did not exist, so created new DetectorAlignment" << endl;
		align->SetZOffset(0, detectorD0Z);
		align->SetZOffset(1, detectorD1Z);
		align->SetZOffset(2, detectorD2Z);
		align->SetZOffset(3, detectorD3Z);
		align->SetZOffset(4, detectorDiaZ);
	}
}

void TAlignment::loadDetectorAlignment(){
	if (align == NULL) {
		settings->goToAlignmentRootDir();
		TFile *alignmentFile = new TFile(settings->getAlignmentFilePath().c_str(), "READ");
		cout << "TAlignment:loadDetectorAlignment(): \""<<settings->getAlignmentFilePath()<<"\""<<endl;
		if (alignmentFile==0){
			cout<<"couldn't open alignment File.Creating new Alignment"<<endl;
			initialiseDetectorAlignment();
		}
		else {
			alignmentFile->cd();
			align = (TDetectorAlignment*)alignmentFile->Get("alignment");
			if (align!=0){
				cout<<"Found DetectorAlignment!"<<endl;
				if(verbosity)align->Print();
				alignmentFile->Close();
				align = (TDetectorAlignment*)align->Clone();
			}
			else{
				cout<<"Couldn't find detectorAlignment in alignmentFile: "<<endl;
				alignmentFile->Close();
				initialiseDetectorAlignment();
			}
		}
		if(alignmentFile) delete alignmentFile;
	}
	//Break if it couldn't create/open an Detecotr Alignment <-Shouldn't happen
	if(!align){
		cerr<<"Problem with creating/loading DetectorAlignment. EXIT!"<<endl;
		exit(-1);
	}
	//be sure that verbosity is set correctly
	align->setVerbosity(verbosity-4>0?verbosity-4:0);
	if(myTrack!=0){
		myTrack->setDetectorAlignment(align);
	}

}
/**
 * Fills a vector of TEvents which are candidates to be used in the Alignment
 * @param nEvents
 * @param startEvent
 */
void TAlignment::createEventVectors(UInt_t nEvents, UInt_t startEvent) {
	initialiseDetectorAlignment();
	if (nEvents == 0) nEvents = eventReader->GetEntries() - startEvent;
	if (nEvents + startEvent > eventReader->GetEntries()) nEvents = eventReader->GetEntries() - startEvent;
	int noHitDet = 0;
	//	int falseClusterSizeDet=0;
	//int noHitDia=0;
	int falseClusterSizeDia = 0;
	int nCandidates = 0;
	int nScreened = 0;
	int nNotInFidCut=0;
	cout << "CREATING VECTOR OF VALID EVENTS..." << endl;

	for (nEvent = startEvent; nEvent < nEvents + startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent - startEvent, nEvents, 1000);
		eventReader->LoadEvent(nEvent);
		if (!eventReader->isValidTrack()) {
			noHitDet++;
			continue;
		}
		if (eventReader->isDetMasked()) {
			nScreened++;
			continue;
		}
		if (eventReader->getNDiamondClusters() != 1) {
			falseClusterSizeDia++;
			continue;
		}
		if(!eventReader->isInFiducialCut()){
			nNotInFidCut++;
			//    	float fiducialValueX=0;
			//    	float fiducialValueY=0;
			//    	for(UInt_t plane=0;plane<4;plane++){
			//    		fiducialValueX+=eventReader->getCluster(plane,TPlaneProperties::X_COR,0).getPosition();
			//    		fiducialValueY+=eventReader->getCluster(plane,TPlaneProperties::Y_COR,0).getPosition();
			//    	}
			//    	fiducialValueX/=4.;
			//    	fiducialValueY/=4.;
			//    	cout<<nEvent<<" "<<eventReader->isInFiducialCut()<<" "<<fiducialValueX<<"/"<<fiducialValueY<<endl;
			continue;
		}
		if(nEvent==startEvent&&verbosity>4)
			cout<<"\nEvent\tvalid\tnClus\tmasked\tFidCut\tAlign"<<endl;
		if(verbosity>20)
			cout<<nEvent<<"\t"<<eventReader->isValidTrack()<<"\t"<<eventReader->getNDiamondClusters()
			<<"\t"<<eventReader->isDetMasked()<<"\t"<<eventReader->isInFiducialCut()<<"\t"<<eventReader->useForAlignment()<<endl;
		if (eventReader->useForAlignment()) {
			bool bBreak = false;
			for(UInt_t det=0;det<TPlaneProperties::getNDetectors()&&!bBreak;det++){
				Float_t clusPos = eventReader->getCluster(det,0).getPosition();
				if(clusPos<0||clusPos>=TPlaneProperties::getNChannels(det)||eventReader->getNClusters(det)!=1){
					bBreak = true;
					cout<<"Do not take event clusPos is not valid...."<<endl;
				}
			}
			if(bBreak)
				continue;
			nCandidates++;

			this->events.push_back(*eventReader->getEvent());
		}
		if (eventReader->useForAnalysis()) {
			cout << "\nFound first Event for Analysis ->BREAK" << endl;
			break;
		}
	}
	if(verbosity){
		cout<<"\nCreated Vector of Events with one and only one Hit in Silicon  + one diamond Cluster + in Fiducial Cut Area"<<endl;
		cout<<"first Event: "<<startEvent<<"\t last Event: "<<nEvent<<"\t total Events: "<<endl;//todo
		cout<<"Cut Flow:\n";
		cout<<"\tTotal Events looked at:     "<<setw(7)<<nEvent-startEvent<<"\n";
		cout<<"\tPlane with no Silicon Hit:   "<<setw(7)<<noHitDet<<"\n";
		cout<<"\tSil Track not in Fid Cut:    "<<setw(7)<<nNotInFidCut<<"\n";
		cout<<"\tNo of Diamond Clust. != 1:   "<<setw(7)<<falseClusterSizeDia<<"\n";
		cout<<"\t                             "<<"-------\n";
		cout<<"\t                             "<<setw(7)<<nCandidates<<" = "<<(Float_t)nCandidates/(Float_t)(nEvent-startEvent)*100.<<"%\n"<<endl;

	}
	align->addEventIntervall(startEvent, nEvent);
	align->setNUsedEvents((UInt_t) events.size());
	align->setAlignmentTrainingTrackFraction(settings->getAlignment_training_track_fraction());
	align->setRunNumber(settings->getRunNumber());
	if(events.size()<=0){
		cerr<<"The Event vector has the size 0. No Alignment can be performed. EXIT."<<endl;
		exit(-1);
	}
}

/**
 * @brief main function of TAlignment, creates eventVector and makes alignemnt
 *
 * @param nEvents	number of Events form which the alignment events are taken
 * @param startEvent first event number
 * @return 1 if everything worked else return 0
 *
 * @todo work on return values
 */
int TAlignment::Align(UInt_t nEvents, UInt_t startEvent,enumDetectorsToAlign detToAlign) {
	if (verbosity>2) {
		cout << "\n\n\nTAlignment::Align:Starting \"" << histSaver->GetPlotsPath() << "\"" << endl;
		cout << "\t\t" << events.size() << "\t";
		cout << "\t\t " << eventReader << " ." << endl;
	}

	if(detToAlign!=diaAlignment&&settings->resetAlignment())
		initialiseDetectorAlignment();
	else
		loadDetectorAlignment();
	if (events.size() == 0) createEventVectors(nEvents, startEvent);

	//create an TTRack object and set the eta distributions.
	if (myTrack == NULL) {
		if(verbosity>2)cout << "TAlignment::Align::create new TTrack" << endl;
		myTrack = new TTrack(align,settings);
		if(verbosity>2)cout << "TAlignment::Align::created new TTrack" << endl;
		for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++){
			TH1F* etaInt = eventReader->getEtaIntegral(det);
			if(etaInt==0){char t;cout<<"eta Int ==0"<<det<<endl;cin >>t;}
			myTrack->setEtaIntegral(det,etaInt );
		}

		myTrack->setVerbosity(verbosity-4>0?verbosity-4:0);
	}
	if(!myTrack){
		cerr<<"could not create my Track ----> EXIT"<<endl;
		exit(-1);
	}
	myTrack->setDetectorAlignment(align);

	/*
	 * distinguish here between diamond and silicon alignment
	 */
	if(detToAlign==silAlignment||detToAlign==bothAlignment){
		AlignSiliconPlanes();
	}
	if(detToAlign==diaAlignment||detToAlign==bothAlignment) {
		if(settings->resetAlignment())
			align->ResetAlignment(TPlaneProperties::getDiamondPlane());
		AlignDiamondPlane();
	}
	align->PrintResults((UInt_t) 1);
	cout<<"Done with alignment of ";
	if (detToAlign==diaAlignment)cout<<"diamond"<<endl;
	else if(detToAlign==silAlignment)cout<<"silicon"<<endl;
	else if(detToAlign==bothAlignment)cout<<"silicon and diamond"<<endl;
	return (1);
}




int TAlignment::AlignSilicon(UInt_t nEvents,UInt_t startEvent){
	return Align(nEvents,startEvent,silAlignment);
}

int TAlignment::AlignDiamond(UInt_t nEvents,UInt_t startEvent){
	return Align(nEvents,startEvent,diaAlignment);
}
/**
 *
 */
void TAlignment::AlignSiliconPlanes() {
	cout<<"Alignment of Silicon Planes. max. Alignment Steps: "<<nAlignSteps<<endl;
	nAlignmentStep = -1;
	if(!align->isPreAligned()){
		if (verbosity)cout << "\nCheck Detector Alignment:" << endl;
		resPlane1 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 1, 0, 3, true);
		resPlane2 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 2, 0, 3, true);
		resPlane3 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 3, 1, 2, true);
	}
	else
		if(verbosity)
			cout<<"Detectors are already pre alignend do not create new Prealigment plots."<<endl;


	doPreAlignment();

	if (verbosity) cout << endl;
	if (verbosity)cout << "Start with Alignment Steps" << endl;
	bool bAlignmentGood = false;
	for (nAlignmentStep = 0; nAlignmentStep < nAlignSteps&&bAlignmentGood!=true; nAlignmentStep++) {
		bAlignmentGood = siliconAlignmentStep(nAlignmentStep-1>=nAlignSteps||bPlotAll);
	}
	int neededAlignmentSteps = nAlignmentStep;
	nAlignmentStep=nAlignSteps;
	if(neededAlignmentSteps<nAlignSteps-1){
		if(verbosity) cout<<"Creating Post Alignment Plots, since the alignment was good after "<< neededAlignmentSteps <<" Steps."<<endl;
		siliconAlignmentStep(nAlignmentStep-1>=nAlignSteps||bPlotAll,false);
	}
	cout<<"Alignment step: "<<neededAlignmentSteps<<endl;
	cout<<"\n\n\n*******************";
	cout<<"\nAlignment of Silicon Planes is done after "<<neededAlignmentSteps<<" steps. Now get final Silicon Alignment Results..."<<endl;
	cout<<"*******************\n\n"<<endl;
	nAlignmentStep = nAlignSteps;
	getFinalSiliconAlignmentResuluts();
	myTrack->setVerbosity(verbosity-4>0?verbosity-4:0);

}

/** "Rough Alignment" of Plane 3 in order to have an aligned system when starting
 *  with alignment of the other planes
 *  It alignes plane 3 up to 4 times when the Offsets are not small enough
 *
 */
void TAlignment::doPreAlignment(){
	if (verbosity)cout << "\nStart Pre-Alignment:" << endl;
	nAlignmentStep=0;
	Float_t xOff3,yOff3,xPhiOff3,yPhiOff3;
	bool xAlignment3,yAlignment3;
	int nTry =0;
	int nTries=4;
	do{
		if(verbosity) cout<< "\n\nPreAlignment Step "<<nTry<<endl;
		resPlane3 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 3, 0, 0, true);
		alignDetector(TPlaneProperties::XY_COR, 3, 0, 0, true,resPlane3);
		xOff3 = TMath::Abs(align->GetLastXOffset(3));
		yOff3 = TMath::Abs(align->GetLastYOffset(3));
		xPhiOff3  = TMath::Abs(align->GetLastPhiXOffset(3));
		xPhiOff3 = TMath::Abs(align->GetLastPhiYOffset(3));
		xAlignment3 = xOff3<settings->getAlignmentPrecisionOffset() && xPhiOff3<settings->getAlignmentPrecisionAngle();
		yAlignment3 = yOff3<settings->getAlignmentPrecisionOffset() && yPhiOff3<settings->getAlignmentPrecisionAngle();
		nTry++;
	}
	while ((!xAlignment3||!yAlignment3)&&nTry<nTries);
}

bool TAlignment::siliconAlignmentStep(bool bPlot, bool bUpdateAlignment) {
	TPlaneProperties::enumCoordinate coordinates = TPlaneProperties::XY_COR;//((nAlignmentStep < nAlignSteps) ? TPlaneProperties::XY_COR : TPlaneProperties::Y_COR);

	cout << "\n\n\nALIGNMENT STEP:\t" << nAlignmentStep + 1 << " of " << nAlignSteps << "\n\n" << endl;
	bool isAlignmentDone = true;

	//Plane 1 with plane 0 and 3
	UInt_t subjectPlane1 = 1;
	UInt_t refPlane1_1 = 0, refPlane1_2 = 3;
	cout << "\n\nAlign Plane " << subjectPlane1 << " with Plane " << refPlane1_1 << " and " << refPlane1_2 << endl;
	resPlane1 = CheckDetectorAlignment(coordinates, subjectPlane1, refPlane1_1, refPlane1_2, false);
	resPlane1 = CheckDetectorAlignment(coordinates, subjectPlane1, refPlane1_1, refPlane1_2, bUpdateAlignment==false && bPlot,resPlane1);
	if(bUpdateAlignment)
		alignDetector(TPlaneProperties::XY_COR, subjectPlane1, refPlane1_1, refPlane1_2, bPlot || bPlotAll,resPlane1);
//	resPlane1 = CheckDetectorAlignment(coordinates, subjectPlane1, refPlane1_1, refPlane1_2, false,resPlane1);
//	alignDetector(TPlaneProperties::Y_COR, subjectPlane1, refPlane1_1, refPlane1_2, bPlot || plotAll,resPlane1);


//	//Plane 2 with plane 0 and 3
	UInt_t subjectPlane2 = 2;
	UInt_t refPlane2_1 = 0, refPlane2_2 = 3;
	cout << "\n\n\nAlign Plane " << subjectPlane2 << " with Plane " << refPlane2_1 << " and " << refPlane2_2 << endl;
	resPlane2 = CheckDetectorAlignment(coordinates, subjectPlane2, refPlane2_1, refPlane2_2, false);
	resPlane2 = CheckDetectorAlignment(coordinates, subjectPlane2, refPlane2_1, refPlane2_2, bUpdateAlignment==false && bPlot,resPlane2);
	if (bUpdateAlignment)
		alignDetector(TPlaneProperties::XY_COR, subjectPlane2, refPlane2_1, refPlane2_2, bPlot || bPlotAll, resPlane2);
//	resPlane2 = CheckDetectorAlignment(coordinates, subjectPlane2, refPlane2_1, refPlane2_2, false,resPlane2);
//	alignDetector(TPlaneProperties::Y_COR, subjectPlane2, refPlane2_1, refPlane2_2, bPlot || plotAll, resPlane2);

////
////	//Plane 3 with 0 and 2
//	UInt_t subjectPlane3 = 3;vector<UInt_t> vecRefPlanes;
//	UInt_t refPlane3_1 = 1, refPlane3_2 = 2;
//	vecRefPlanes.push_back(refPlane3_1);
//	vecRefPlanes.push_back(refPlane3_2);
//	cout << "\n\nAlign Plane " << subjectPlane3 << " with Plane " << refPlane3_1 << " and " << refPlane3_2 << endl;
//	resPlane3 = CheckDetectorAlignment(coordinates, subjectPlane3, vecRefPlanes, false);
//	resPlane3 = CheckDetectorAlignment(coordinates, subjectPlane3, vecRefPlanes, false,resPlane3);
//	alignDetector(TPlaneProperties::XY_COR, subjectPlane3, vecRefPlanes, bPlot || plotAll, resPlane3);
////	resPlane3 = CheckDetectorAlignment(coordinates, subjectPlane3, vecRefPlanes, false,resPlane3);
////	alignDetector(TPlaneProperties::Y_COR, subjectPlane3, vecRefPlanes, bPlot || plotAll, resPlane3);


	Float_t xOff1 = TMath::Abs(align->GetLastXOffset(1));
	Float_t yOff1 = TMath::Abs(align->GetLastYOffset(1));
	bool test1 = (xOff1 > settings->getAlignmentPrecisionOffset() || yOff1 >  settings->getAlignmentPrecisionOffset()||xOff1 == -9999||yOff1 == -9999);
	Float_t xOff2 = TMath::Abs(align->GetLastXOffset(2));
	Float_t yOff2 = TMath::Abs(align->GetLastYOffset(2));
	bool test2 =  (xOff2 > settings->getAlignmentPrecisionOffset() || yOff2 > settings->getAlignmentPrecisionOffset()||xOff2 == -9999||yOff2 == -9999);
	isAlignmentDone = !test1&&!test2;

	if(verbosity>3)	align->PrintResults(0);
	if(verbosity) cout<< xOff1<< "/" << yOff1 << "\t" << xOff2 <<"/"<< yOff2<<endl;
	cout<<test1<<" "<<test2<<" "<<flush;
	if (isAlignmentDone)
		cout << "AlignmentStep was successfully done" << endl;
	else
		cout << "AlignmentStep was NOT the finally alignment step" << endl;
	cout << endl;
	return (isAlignmentDone);
}

void TAlignment::AlignDiamondPlane() {
	cout << "\n\n\n*******************************************************" << endl;
	cout << "*******************************************************" << endl;
	cout << "***************** Align Diamond ***********************" << endl;
	cout << "*******************************************************" << endl;
	cout << "*******************************************************\n" << endl;

	//create ReferencePlane Vector: using all Silicon Planes for Alignment 0,1,2,3
	UInt_t diaPlane = 4;
	vector<UInt_t> vecRefPlanes;
	for (UInt_t i = 0; i < 4; i++)
		if (i != diaPlane) vecRefPlanes.push_back(i);
	nDiaAlignmentStep = -1;
	//checking Residual
	TResidual resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false, true);
	bool diaAlignmentDone = false;

	for (nDiaAlignmentStep = 0; (nDiaAlignmentStep < nDiaAlignSteps) && (!diaAlignmentDone||nDiaAlignmentStep<2)	 ; nDiaAlignmentStep++) {
		cout << "\n\n " << nDiaAlignmentStep << " of " << nDiaAlignSteps << " Steps..." << endl;
		//do Alignment using the resDia Residual
		alignStripDetector(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false||bPlotAll, resDia);
		//creating new Residual resDia
		resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes);

		//check if this was the alignment has changed less than certain values read from the settings file
		diaAlignmentDone = (TMath::Abs(align->GetLastXOffset(TPlaneProperties::getDiamondPlane())) < settings->getAlignmentPrecisionOffset());
		diaAlignmentDone =diaAlignmentDone &&
				(TMath::Abs(align->GetLastPhiXOffset(TPlaneProperties::getDiamondPlane())) < settings->getAlignmentPrecisionAngle());
		cout<<" "<<( TMath::Abs(align->GetLastXOffset(TPlaneProperties::getDiamondPlane())) < settings->getAlignmentPrecisionOffset())<<flush;
		cout<<" "<<( TMath::Abs(align->GetLastPhiXOffset(TPlaneProperties::getDiamondPlane())) < settings->getAlignmentPrecisionAngle())<<endl;
		if(diaAlignmentDone)
			cout<<"dia Alignment is Done! after " <<nDiaAlignmentStep <<" Steps\t"<<flush;
		else if(verbosity)
			cout<<"dia Alignment is not done, another step!\t"<<flush;
		if(verbosity||diaAlignmentDone)
			cout<<align->GetLastXOffset(TPlaneProperties::getDiamondPlane())<<" "<<
				align->GetLastPhiXOffset(TPlaneProperties::getDiamondPlane())*1000<<" mrad"<<endl;
	}
	nDiaAlignmentStep = nDiaAlignSteps;
	resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, true, true, resDia);
	CheckStripDetectorAlignmentChi2(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, true, true, settings->getAlignment_chi2());
	align->setDiamondDate();
}


/**
 *
 * @param cor
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 * @param bPlot
 * @param resOld
 * @return
 */
TResidual TAlignment::alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld) {
	vector<UInt_t> vecRefPlanes;
	vecRefPlanes.push_back(refPlane1);
	if (refPlane1 != refPlane2) vecRefPlanes.push_back(refPlane2);
	return (alignDetector(cor, subjectPlane, vecRefPlanes, bPlot, resOld));
}



TResidual TAlignment::alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {

	if (true) {    //some outputs
		cout << "TAlignment::AlignDetector::\taligning Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << vecRefPlanes.size() << " Plane(s): ";
		for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
			cout << vecRefPlanes.at(i) << " ";
		cout << "\t";
		if (bPlot) cout << " plots are created\t";
		if (resOld.isTestResidual())
			cout << "resOld is a testResidual" << endl;
		else {
			if(verbosity>2)
			cout <<"\nfollowing Residuals are used as Input residuals:"<< endl;
			else
				cout<<", with a real residual as an input"<<endl;
			if(verbosity>2)resOld.Print(2);
		}

		if(verbosity>1)
		for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
			if (vecRefPlanes.at(i) == subjectPlane) cerr << "Plane " << subjectPlane << " is used as a reference Plane and as a subject Plane" << endl;
	}

	//get Residual
	TResidual res = this->getResidual(cor, subjectPlane, vecRefPlanes, bPlot, resOld,getClusterCalcMode(subjectPlane));
	Float_t x_offset = res.getXOffset();
	Float_t phix_offset = res.getPhiXOffset();
	Float_t y_offset = res.getYOffset();
	Float_t phiy_offset = res.getPhiYOffset();

	if(verbosity)printf("Plane: %d\tCorrection Values: X: %2.6f,  PhiX: %2.6f,   Y: %2.6f,  PhiY: %2.6f\n",subjectPlane, x_offset, phix_offset, y_offset, phiy_offset);
	if(verbosity&&!res.isTestResidual())res.Print();

	//save corrections to alignment
	if( res.getUsedTracks()/(Float_t)(vecXPred.size()) < 0.1){
		cout<<"Something is wrong used less than 10% of the tracks for alignment: "<<res.getUsedTracks()<<" Tracks out of "<<vecXPred.size()<<endl;
		char t;
		cin >> t;
	}
	if(verbosity)cout<<"Set Alignment of Detector "<<subjectPlane<<endl;

	if (TPlaneProperties::X_COR == cor || TPlaneProperties::XY_COR == cor) {
		if ((x_offset != N_INVALID)&& (phiy_offset != N_INVALID) ) {
			align->AddToXOffset(subjectPlane, x_offset);
			align->AddToPhiXOffset(subjectPlane, phix_offset);
		}
	}
	if (TPlaneProperties::Y_COR == cor || TPlaneProperties::XY_COR == cor) {
		if ((y_offset != N_INVALID) && (phix_offset != N_INVALID)){
			align->AddToYOffset(subjectPlane, y_offset);
			align->AddToPhiYOffset(subjectPlane, phiy_offset);
		}
	}
	if(verbosity>3){
		cout<<"Press a Key to confirm and press enter!"<<endl;
		char t; cin>>t;
	}
	this->myTrack->setDetectorAlignment(align);
	return (res);
}

TResidual TAlignment::alignStripDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {
	if (verbosity) {
		cout << "TAlignment::AlignStripDetector::\taligning Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << vecRefPlanes.size() << " Planes: ";

		for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
			cout << vecRefPlanes.at(i) << " ";
		cout << "\t";
		if (bPlot) cout << " plots are created\t";
		if (resOld.isTestResidual())
			cout << "resOld is a testResidual" << endl;
		else {
			cout << endl;
			resOld.Print(2);
		}
		for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
			if (vecRefPlanes.at(i) == subjectPlane) cerr << "Plane " << subjectPlane << " is used as a reference Plane and as a subject Plane" << endl;
	}
	//get residual
	TResidual res = this->getStripResidual(cor, subjectPlane, vecRefPlanes, false, bPlot, resOld);
	Float_t x_offset = res.getXOffset();
	Float_t phix_offset = res.getPhiXOffset();
	Float_t y_offset = res.getYOffset();
	Float_t phiy_offset = res.getPhiYOffset();

	if (verbosity||true) cout << "Correction Values: X:" << x_offset << ",  PhiX: " << phix_offset << ",   Y: " << y_offset << ",  PhiY: " << phiy_offset << "\n" << endl;

	if (subjectPlane == 4 && nDiaAlignmentStep == 0) {
		align->AddToXOffset(subjectPlane, x_offset);
		return (res);
	}
	if (TPlaneProperties::X_COR == cor || TPlaneProperties::XY_COR == cor) {
		align->AddToXOffset(subjectPlane, x_offset);
		align->AddToPhiXOffset(subjectPlane, phix_offset);
	}
	if (TPlaneProperties::Y_COR == cor || TPlaneProperties::XY_COR == cor) {
		align->AddToYOffset(subjectPlane, y_offset);
		align->AddToPhiYOffset(subjectPlane, phiy_offset);
	}

	return res;
}

TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld, TCluster::calculationMode_t mode,resCalcMode calcMode,Float_t maxChi2) {
	vector<UInt_t> vecRefPlanes;
	vecRefPlanes.push_back(refPlane1);
	if (refPlane1 != refPlane2) vecRefPlanes.push_back(refPlane2);
	return getResidual(cor, subjectPlane, vecRefPlanes, bPlot, resOld, mode,calcMode,maxChi2);
}

TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld, TCluster::calculationMode_t mode,resCalcMode calcMode,Float_t maxChi2) {
	for(UInt_t refPlane1=0;refPlane1<vecRefPlanes.size()-1;refPlane1++){
		for(UInt_t refPlane2=refPlane1+1;refPlane2<vecRefPlanes.size();refPlane2++){
			if(vecRefPlanes.at(refPlane1)==vecRefPlanes.at(refPlane2)){
				cout<<" Removing "<<refPlane1<<" "<<refPlane2<<":"<< vecRefPlanes.at(refPlane1)<<" from"<<vecRefPlanes.size()<<"-";
				vecRefPlanes.erase(vecRefPlanes.begin()+refPlane1);
				cout<<vecRefPlanes.size()<<endl;
			}
		}
	}
	stringstream refPlaneString;
	for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
		if (i == 0)
			refPlaneString << vecRefPlanes.at(i);
		else if (i + 1 < vecRefPlanes.size())
			refPlaneString << "_" << vecRefPlanes.at(i);
		else
			refPlaneString << "_and_" << vecRefPlanes.at(i);
	if (verbosity>0)
		cout << "\tTAlignment::getResidual of Plane " << subjectPlane<<" " << TPlaneProperties::getCoordinateString(cor)
	         << " with " << refPlaneString.str() << ", plotting: " << bPlot << "  with " << alignmentPercentage << "\t"
	         << resOld.isTestResidual()<< endl;
	if(verbosity>0&&resOld.isTestResidual())resOld.Print(1);

	clearMeasuredVectors();
	Float_t xDelta,xRes,yRes, yDelta,xLabMeasMetric, yLabMeasMetric, xDetMeasuredMetric, yDetMeasuredMetric,  xLabPredictedMetric, yLabPredictedMetric, resxtest, resytest, xPhi,yPhi,chi2x,chi2y;
	TPositionPrediction* predictedPostionMetric = 0;
	UInt_t nUsedEvents=0;
	for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent, events.size(),100,false,false);
		myTrack->setEvent(&events.at(nEvent));
		if (verbosity > 5) cout << "Sil Alignment - Event No.:"<< nEvent << endl;
		xLabMeasMetric = myTrack->getPositionInLabFrame(TPlaneProperties::X_COR, subjectPlane, mode,myTrack->getEtaIntegral(subjectPlane*2));
		yLabMeasMetric = myTrack->getPositionInLabFrame(TPlaneProperties::Y_COR, subjectPlane, mode,myTrack->getEtaIntegral(subjectPlane*2+1));
		if (verbosity > 5) cout<<"Predict Position: "<<endl;
		predictedPostionMetric = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, verbosity>8);
		if(verbosity>6)	predictedPostionMetric->Print(1);
		xLabPredictedMetric = predictedPostionMetric->getPositionX();
		yLabPredictedMetric = predictedPostionMetric->getPositionY();
		xRes = predictedPostionMetric->getSigmaX();
		yRes = predictedPostionMetric->getSigmaY();
		xPhi = predictedPostionMetric->getPhiX()*360./TMath::TwoPi();
		yPhi = predictedPostionMetric->getPhiY()*360./TMath::TwoPi();
		xDelta = xLabMeasMetric - xLabPredictedMetric;    //X_OBS-X_Pred
		yDelta = yLabMeasMetric - yLabPredictedMetric;    //Y_OBS-Y_Pred
		resxtest = TMath::Abs(xDelta - resOld.getXMean()) / TMath::Max((Double_t)resOld.getXSigma(),align->getResolution(subjectPlane,cor));
		resytest = TMath::Abs(yDelta - resOld.getYMean()) / resOld.getYSigma();
		xDetMeasuredMetric = myTrack->getXMeasuredClusterPositionMetricSpace(subjectPlane, mode,myTrack->getEtaIntegral(subjectPlane*2));
		yDetMeasuredMetric = myTrack->getYMeasuredClusterPositionMetricSpace(subjectPlane, mode,myTrack->getEtaIntegral(subjectPlane*2+1));
		chi2x = predictedPostionMetric->getChi2X();
		chi2y = predictedPostionMetric->getChi2Y();
		bool useEvent=false;
		if(calcMode==normalCalcMode){
			if(cor==TPlaneProperties::X_COR)		useEvent = resxtest < res_keep_factor;
			else if(cor==TPlaneProperties::Y_COR) 	useEvent = resytest < res_keep_factor;
			else if(cor==TPlaneProperties::XY_COR)	useEvent = resxtest < res_keep_factor && resytest < res_keep_factor;
		}
		else if(calcMode==chi2CalcMode){
			if(cor==TPlaneProperties::X_COR)		useEvent = chi2x < maxChi2;
			else if(cor==TPlaneProperties::Y_COR) 	useEvent = chi2y < maxChi2;
			else if(cor==TPlaneProperties::XY_COR)	useEvent = chi2x < maxChi2 && chi2y < maxChi2;
		}
		if (useEvent) {
			if(verbosity>4)cout<<nEvent<<" add to vector"<<endl;
			vecXObs.push_back(xLabMeasMetric);
			vecYObs.push_back(yLabMeasMetric);
			vecXDelta.push_back(xDelta);
			vecYDelta.push_back(yDelta);
			vecXPred.push_back(xLabPredictedMetric);
			vecYPred.push_back(yLabPredictedMetric);
			vecXMeasured.push_back(xDetMeasuredMetric);
			vecYMeasured.push_back(yDetMeasuredMetric);
			vecXPhi.push_back(xPhi);
			vecYPhi.push_back(yPhi);
			vecXResPrediction.push_back(xRes);
			vecYResPrediction.push_back(yRes);
			vecXChi2.push_back(chi2x);
			vecYChi2.push_back(chi2y);
			int det = subjectPlane*2+cor==TPlaneProperties::X_COR?0:1;
			vecClusterSize.push_back(myTrack->getClusterSize(det,0));
			nUsedEvents++;
		}
		if (verbosity > 4) cout << "Measured: " <<xDetMeasuredMetric << " / " << yDetMeasuredMetric << endl;
		if (verbosity > 4) cout << "Observed: " << xLabMeasMetric << " / " << yLabMeasMetric << endl;
		if (verbosity > 4) cout << "Predicted: " << predictedPostionMetric->getPositionX() << " / " << predictedPostionMetric->getPositionY() << endl;
		if (verbosity > 4) cout << "Predicted: " << xLabPredictedMetric << " / " << yLabPredictedMetric << endl;
		if (verbosity > 4) cout << "Delta:    " << xDelta << " / " << yDelta << endl;
		if (verbosity > 4) cout << "ResTest:  " << resxtest << " / " << resytest << "\n\n" << endl;
		predictedPostionMetric->Delete();
	}
	cout<<"using "<<vecXDelta.size() <<" Events"<<endl;
	if(nUsedEvents==0){
		cerr<< "cannot calculate Residual/Alignment for 0 Events. BREAK!"<<endl;
		align->PrintResults(1);
		cerr<< "cannot calculate Residual/Alignment for 0 Events. BREAK!"<<endl;
		exit(-1);
	}
	if (verbosity > 2) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size()
												  <<" "<<vecXMeasured.size()<<" "<<vecYMeasured.size()<< endl;
	//first estimate residuals widths
	TResidual res;
	res.setVerbosity(verbosity>4?verbosity-4:0);
	res.setResKeepFactor(res_keep_factor);
	res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta);
	this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot,calcMode==chi2CalcMode&&bPlot,calcMode==chi2CalcMode);
	clearMeasuredVectors();
	return res;

}
/**
 * @brief creates element TResidual to adjust the alignment
 *
 * creates a vector of pedicted X positions, predicted Y positions, delta X and delta Y
 * and use the function calculateResidual to get the residual with this vectors
 * @param	cor coordindate for which the residual is calculated
 * @param	subjectPlane plane for which the residual should be calculated
 * @param	refPlane1 first reference plane
 * @param	refPlane2 second reference plane
 * @param	bPlot	variable to create plots or not
 *
 * @return
 */
TResidual TAlignment::getStripResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld,
									   TCluster::calculationMode_t mode,resCalcMode calcMode,Float_t maxChi2)
{
	stringstream refPlaneString;
	for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
		if (i == 0)
			refPlaneString << vecRefPlanes.at(i);
		else if (i + 1 < vecRefPlanes.size())
			refPlaneString << "_" << vecRefPlanes.at(i);
		else
			refPlaneString << "_and_" << vecRefPlanes.at(i);
	if (verbosity)
		cout << "TAlignment::getStripResidual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << refPlaneString.str() << ", plotting: " << bPlot<< "using "<<events.size()<<" Events\t" << resOld.isTestResidual() << endl;
	clearMeasuredVectors();

	Float_t deltaXMetric, deltaY;
	Float_t xMeasuredMetric;
	Float_t xPositionObservedMetric,yPositionObservedMetric;
	Float_t xPredSigma,yPredSigma;
	TPositionPrediction* predictedPostion = 0;
	Float_t resxtest, resytest;
	Float_t xPhi,yPhi;
	Float_t xPredictedMetric,yPredictedMetric;
	Int_t clusterSize;
	Float_t chi2x,chi2y;
	for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent, events.size());
		myTrack->setEvent(&events.at(nEvent));
		if (verbosity > 3) cout << "Event no.: " << nEvent << endl;
		predictedPostion = myTrack->predictPosition(subjectPlane, vecRefPlanes, silCalcMode, verbosity > 7);
		xMeasuredMetric = myTrack->getXMeasuredClusterPositionMetricSpace(subjectPlane,diaCalcMode);
		xPositionObservedMetric = myTrack->getStripXPosition(subjectPlane, predictedPostion->getPositionY(), diaCalcMode);
		yPositionObservedMetric = myTrack->getPositionInLabFrame(TPlaneProperties::Y_COR, subjectPlane, diaCalcMode);
		xPredictedMetric=predictedPostion->getPositionX();
		yPredictedMetric=predictedPostion->getPositionY();
		deltaXMetric = xPositionObservedMetric - xPredictedMetric; //X_OBS-X_Pred
		deltaY = yPositionObservedMetric - yPredictedMetric; //Y_OBS-Y_Pred
		resxtest = TMath::Abs(TMath::Abs(deltaXMetric - resOld.getXMean()) / resOld.getXSigma());
		resytest = TMath::Abs(TMath::Abs(deltaY - resOld.getYMean()) / resOld.getYSigma());
		xPhi = predictedPostion->getPhiX()*360./TMath::TwoPi();
		yPhi = predictedPostion->getPhiY()*360./TMath::TwoPi();
		xPredSigma = predictedPostion->getSigmaX();
		yPredSigma = predictedPostion->getSigmaY();
		clusterSize = myTrack->getClusterSize(subjectPlane*2+cor==TPlaneProperties::X_COR?0:1,0);
		chi2x = predictedPostion->getChi2X();
		chi2y = predictedPostion->getChi2Y();
		int oldVerb=verbosity;

		if(verbosity>3)	predictedPostion->Print();
		if (verbosity > 3) events.at(nEvent).getPlane(subjectPlane).Print();
		if (verbosity > 3) cout << "MeasuredChannel: " << myTrack->getXMeasuredClusterPositionChannelSpace(subjectPlane) << "/" << myTrack->getYMeasuredClusterPositionChannelSpace(subjectPlane) << endl;
		if (verbosity > 3) cout << "MeasuredMetric: " << myTrack->getXMeasuredClusterPositionMetricSpace(subjectPlane) << "/" << myTrack->getYMeasuredClusterPositionMetricSpace(subjectPlane) << endl;
		if (verbosity > 3) cout << "Observed: " << xPositionObservedMetric << " / " << yPositionObservedMetric << endl;
		if (verbosity > 3) cout << "Predicted: " << xPredictedMetric << "/" << yPredictedMetric << endl;
		if (verbosity > 3) cout << "Delta:    " << deltaXMetric << " / " << yPositionObservedMetric << endl;
		if (verbosity > 3) cout << "ResTest:  " << resxtest << " / " << resytest << "\n\n" << endl;

		bool takeTrack = false;
		if (resxtest < res_keep_factor && resytest < res_keep_factor && calcMode == normalCalcMode) { // take Track/Event if residualtest is small enough and in normal calc mode
			takeTrack = true;
		}
		else if (calcMode == chi2CalcMode && chi2x <maxChi2 && chi2y < maxChi2)
			takeTrack = true;
		if(takeTrack ){
			vecXObs.push_back(xPositionObservedMetric);
			vecYObs.push_back(yPositionObservedMetric);
			vecXMeasured.push_back(xMeasuredMetric);
			vecXDelta.push_back(deltaXMetric);
			vecYDelta.push_back(deltaY);
			vecXPred.push_back(xPredictedMetric);//predictedPostion->getPositionX());
			vecYPred.push_back(yPredictedMetric);//predictedPostion->getPositionY());
			vecXPhi.push_back(xPhi);
			vecYPhi.push_back(yPhi);
			vecXChi2.push_back(chi2x);
			vecYChi2.push_back(predictedPostion->getChi2Y());
			vecXResPrediction.push_back(xPredSigma);
			vecYResPrediction.push_back(yPredSigma);
			vecClusterSize.push_back(clusterSize);
		}
		if (verbosity > 3) cout << deltaXMetric << " " << deltaY << endl;
		predictedPostion->Delete();
	}

	if (verbosity > 2) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << endl;

	//first estimate residuals widths
	TResidual res;
	res.setResKeepFactor(res_keep_factor);
	res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta);
	this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot, bAlign);
	clearMeasuredVectors();
	return res;
}
//
//TResidual TAlignment::getStripResidualChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2, TCluster::calculationMode_t mode) {
//
//	TResidual resOld;
//	stringstream refPlaneString;
//	for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
//		if (i == 0)
//			refPlaneString << vecRefPlanes.at(i);
//		else if (i + 1 < vecRefPlanes.size())
//			refPlaneString << "_" << vecRefPlanes.at(i);
//		else
//			refPlaneString << "_and_" << vecRefPlanes.at(i);
//	if (verbosity) cout << "TAlignment::getStripResidual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << refPlaneString.str() << ", plotting: " << bPlot << "  with " << alignmentPercentage << "\t" << resOld.isTestResidual() << endl;
//	clearMeasuredVectors();
//
//	Float_t deltaX, deltaY;
//	Float_t xPositionObserved, yPositionObserved;
//	Float_t xPredSigma,yPredSigma;
//	Float_t xMeasured, yMeasured;
//	TPositionPrediction* predictedPosition = 0;
//	Float_t resxtest, resytest;
//	Float_t chi2x, chi2y;
//	Float_t xPhi,yPhi;
//	Float_t clusterSize;
//	for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
//		TRawEventSaver::showStatusBar(nEvent, events.size());
//		myTrack->setEvent(&events.at(nEvent));
//		predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, nEvent < 1&&verbosity>2);
//		chi2x = predictedPosition->getChi2X();
//		chi2y = predictedPosition->getChi2Y();
//		xPositionObserved = myTrack->getStripXPosition(subjectPlane, predictedPosition->getPositionY(), TCluster::maxValue);
//		yPositionObserved = myTrack->getPositionInLabFrame(TPlaneProperties::Y_COR, subjectPlane, TCluster::maxValue);
//		deltaX = xPositionObserved - predictedPosition->getPositionX();    //X_OBS-X_Pred
//		deltaY = yPositionObserved - predictedPosition->getPositionY();    //Y_OBS-Y_Pred
//		xMeasured = myTrack->getMeasured(TPlaneProperties::X_COR, subjectPlane, TCluster::maxValue);
//		yMeasured = myTrack->getMeasured(TPlaneProperties::Y_COR, subjectPlane, TCluster::maxValue);
//		resxtest = TMath::Abs(deltaX - resOld.getXMean()) / resOld.getXSigma();
//		resytest = TMath::Abs(deltaY - resOld.getYMean()) / resOld.getYSigma();
//		xPhi=predictedPosition->getPhiX()*360./TMath::TwoPi();
//		yPhi=predictedPosition->getPhiY()*360./TMath::TwoPi();
//		xPredSigma = predictedPosition->getSigmaX();
//		yPredSigma = predictedPosition->getSigmaY();
//		clusterSize = eventReader->getClusterSize(subjectPlane*2,0);
//		if (verbosity > 3) cout << "Event no.: " << nEvent << endl;
//		if (verbosity > 3) events.at(nEvent).getPlane(subjectPlane).Print();
//		if (verbosity > 3) cout << "Measured: " << myTrack->getXMeasured(subjectPlane) << "/" << myTrack->getYMeasured(subjectPlane) << endl;
//		if (verbosity > 3) cout << "Observed: " << xPositionObserved << " / " << yPositionObserved << endl;
//		if (verbosity > 3) cout << "Predicted: " << predictedPosition->getPositionX() << " / " << predictedPosition->getPositionY() << endl;
//		if (verbosity > 3) cout << "Delta:    " << deltaX << " / " << yPositionObserved << endl;
//		if (verbosity > 3) cout << "Chi2:     " << chi2x << " / " << chi2y << endl;
//		if (verbosity > 3) cout << "ResTest:  " << resxtest << " / " << resytest << endl;
//
//		if (chi2x < settings->getAlignment_chi2() && chi2y < settings->getAlignment_chi2()) {
//			if (verbosity > 3) cout << "Take event for Residual" << endl;
//			vecXObs.push_back(xPositionObserved);
//			vecYObs.push_back(yPositionObserved);
//			vecXDelta.push_back(deltaX);
//			vecYDelta.push_back(deltaY);
//			vecXPred.push_back(predictedPosition->getPositionX());
//			vecYPred.push_back(predictedPosition->getPositionY());
//			vecXMeasured.push_back(xMeasured);
//			vecYMeasured.push_back(yMeasured);
//			vecXPhi.push_back(xPhi);
//			vecYPhi.push_back(yPhi);
//			vecXResPrediction.push_back(xPredSigma);
//			vecYResPrediction.push_back(yPredSigma);
//			vecXChi2.push_back(predictedPosition->getChi2X());
//			vecYChi2.push_back(predictedPosition->getChi2Y());
//			vecClusterSize.push_back(clusterSize);
//		} else if (verbosity > 3) cout << "through event away..." << endl;
//		if (verbosity > 3) cout << "\n\n" << endl;
//		predictedPosition->Delete();
//	}
//	cout<<"\tused "<<vecXDelta.size()<<" Events."<<endl;
//	if (bAlign) {
//		align->setNDiamondAlignmentEvents((UInt_t) vecXDelta.size());
//		align->setDiaChi2(settings->getAlignment_chi2());
//	}
//
//	if (verbosity) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << " " << vecXObs.size() << " " << vecYObs.size() << endl;
//	if (verbosity) cout << "used " << vecXDelta.size() << " Events of " << events.size() << " with a Chi2 < " << settings->getAlignment_chi2() << endl;
//	//first estimate residuals widths
//	TResidual res;
//	res.setResKeepFactor(res_keep_factor);
//	res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta);
//	this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot, bAlign);
//
//	clearMeasuredVectors();
//	return res;
//}
/**
 * calculateResidual if there is no residuals calculatet to cut on the res_keep_factor;
 * this funcition opens the other calculateResidual function but uses a TResidual res
 * which is using all items...
 */
TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld) {
	vector<UInt_t> vecRefPlanes;
	vecRefPlanes.push_back(refPlane1);
	vecRefPlanes.push_back(refPlane2);
	return CheckDetectorAlignment(cor, subjectPlane, vecRefPlanes, bPlot, resOld);
}

TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {
	if (verbosity) {
		cout << "TAlignment::checkDetectorAlignment\n\t    check " << TPlaneProperties::getCoordinateString(cor) << " coordinate of Plane " << subjectPlane << " with Planes: ";
		for(UInt_t i=0;i<vecRefPlanes.size();i++)
			cout<<vecRefPlanes.at(i)<<" ";
		cout<<endl;
	}
	TResidual res = getResidual(cor, subjectPlane, vecRefPlanes, bPlot, resOld,getClusterCalcMode(subjectPlane));
	if (verbosity) res.Print();
	return res;
}

TResidual TAlignment::CheckStripDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld) {
	TResidual res = getStripResidual(cor, subjectPlane, vecRefPlanes, false, false, resOld);
	if (verbosity) cout << endl;
	res.SetTestResidual(false);
	res = getStripResidual(cor, subjectPlane, vecRefPlanes, bAlign, bPlot, res);
	if (verbosity) res.Print();
	return res;
}

TResidual TAlignment::CheckStripDetectorAlignmentChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2) {
	//	getStripResidual(cor,subjectPlane,vecRefPlanes,false,false,maxChi2);
	if (verbosity) cout << endl;
	TResidual res = getStripResidual(cor, subjectPlane, vecRefPlanes, bAlign, bPlot,TResidual(),diaCalcMode,chi2CalcMode, maxChi2);
	if (verbosity) res.Print();
	return res;
}
/**
 *
 */
void TAlignment::saveAlignment() {
	stringstream fileName;
	settings->goToAlignmentRootDir();
	TFile *alignmentFile = new TFile(settings->getAlignmentFilePath().c_str(), "RECREATE");
	cout << "TAlignment:saveAlignment(): path: \"" << settings->getAlignmentFilePath().c_str() << "\"" << endl;
	alignmentFile->cd();
	align->SetName("alignment");
	stringstream title;
	title << "alignment of Run " << settings->getRunNumber();
	align->SetTitle(title.str().c_str());
	align->Write();
	alignmentFile->Write();
	delete alignmentFile;
}

/**
 *
 */
void TAlignment::getChi2Distribution(Float_t maxChi2) {
	vecXChi2.clear();
	vecYChi2.clear();
	vecXPhi.clear();
	vecYPhi.clear();
	//	if(verbosity)
	cout << "\nTAlignment::getChi2Distribution final" << endl;
	//	Float_t xPositionObserved,yPositionObserved,deltaX,deltaY,resxtest,resytest;
	TPositionPrediction* predictedPosition = 0;
	vector<UInt_t> vecRefPlanes;

	for (UInt_t i = 0; i < 4; i++)
		vecRefPlanes.push_back(i);
	//	UInt_t oldVerbosity=myTrack->getVerbosity();
	vector<Float_t> vecSumDeltaX;
	vector<Float_t> vecSumDeltaY;
	Float_t xPhi, yPhi;
	Float_t chi2X, chi2Y;
	Float_t xPredSigma,yPredSigma;
	for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent, events.size());
		myTrack->setEvent(&events.at(nEvent));
		Float_t sumDeltaX = 0;
		Float_t sumDeltaY = 0;
		UInt_t subjectPlane = 0;
		predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
		chi2X = predictedPosition->getChi2X();
		chi2Y = predictedPosition->getChi2Y();
		xPhi = predictedPosition->getPhiX()*360./TMath::TwoPi();
		yPhi = predictedPosition->getPhiY()*360./TMath::TwoPi();
		xPredSigma = predictedPosition->getSigmaX();
		yPredSigma = predictedPosition->getSigmaY();
		for (subjectPlane = 0; subjectPlane < 4; subjectPlane++) {
			if (subjectPlane != 0) {
				predictedPosition->Delete();
				predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
			}
			Float_t deltaX = myTrack->getXPositionMetric(subjectPlane);
			Float_t deltaY = myTrack->getYPositionMetric(subjectPlane);

			deltaX -= predictedPosition->getPositionX();
			deltaY -= predictedPosition->getPositionY();
			sumDeltaX += TMath::Abs(deltaX);
			sumDeltaY += TMath::Abs(deltaY);
		}    //for loop over subjectPlane

		//		if (predictedPosition->getChi2X() < maxChi2 && predictedPosition->getChi2Y() < maxChi2) {
		vecXChi2.push_back(chi2X);
		vecYChi2.push_back(chi2Y);
		vecSumDeltaX.push_back(sumDeltaX);
		vecSumDeltaY.push_back(sumDeltaY);
		vecXPhi.push_back(xPhi);
		vecYPhi.push_back(yPhi);
		vecXResPrediction.push_back(xPredSigma);
		vecYResPrediction.push_back(yPredSigma);
		//		}    //end if chi2x<maxCi && chi2y<maxChi
		predictedPosition->Delete();
	}    //end for loop over nEvent

	//	myTrack->setVerbosity(oldVerbosity);
	stringstream histName;
	//Chi2X Distribution
	histName.str("");
	if (nAlignmentStep == -1)
		histName << "hPreAlignment";
	else if (nAlignmentStep >= nAlignSteps - 1)
		histName << "hPostAlignment";
	else
		histName << "h_" << nAlignmentStep << "_Step";
	histName << "_Chi2X_Distribution";
	histSaver->SetRange(0,20);
	TH1F *histoChi2X = histSaver->CreateDistributionHisto(histName.str(), vecXChi2, 1024, HistogrammSaver::manual,0,20);
	histoChi2X->GetXaxis()->SetTitle("#chi^{2} of X plane");
	histoChi2X->GetYaxis()->SetTitle("number of entries");
	histoChi2X->GetXaxis()->SetRangeUser(0,20);
	histSaver->SaveHistogram(histoChi2X);

	//Chi2Y Distribution
	histName.str("");
	if (nAlignmentStep == -1)
		histName << "hPreAlignment";
	else if (nAlignmentStep >= nAlignSteps - 1)
		histName << "hPostAlignment";
	else
		histName << "h_" << nAlignmentStep << "_Step";
	histName << "_Chi2Y_Distribution";
	histSaver->SetRange(0,20);
	TH1F *histoChi2Y = histSaver->CreateDistributionHisto(histName.str(), vecYChi2, 1024, HistogrammSaver::manual,0,20);
	histoChi2Y->GetXaxis()->SetTitle("#chi^{2} of Y plane");
	histoChi2Y->GetYaxis()->SetTitle("number of entries");
	histoChi2Y->GetXaxis()->SetRangeUser(0,20);
	histSaver->SaveHistogram(histoChi2Y);

	histName.str("");
	if (nAlignmentStep == -1)
		histName << "hPreAlignment";
	else if (nAlignmentStep >= nAlignSteps - 1)
		histName << "hPostAlignment";
	else
		histName << "h_" << nAlignmentStep << "_Step";
	histName << "_Chi2X_vs_SumDeltaX";
	cout << "CREATE: " << histName.str() << endl;
	TH2F *histo = histSaver->CreateScatterHisto(histName.str(),vecSumDeltaX, vecXChi2 );
	histo->GetYaxis()->SetTitle("Sum of |Delta X| / #mum ");
	histo->GetXaxis()->SetTitle("#chi^{2} of X fit");
	histSaver->SaveHistogram(histo);

	histName.str("");
	if (nAlignmentStep == -1)
		histName << "hPreAlignment";
	else if (nAlignmentStep >= nAlignSteps - 1)
		histName << "hPostAlignment";
	else
		histName << "h_" << nAlignmentStep << "_Step";
	histName << "_Chi2Y_vs_SumDeltaY";
	cout << "CREATE PLOT: \"" << histName.str() << "\"" << endl;
	TH2F *histo1 = histSaver->CreateScatterHisto(histName.str(), vecSumDeltaY,vecYChi2);
	histo1->GetYaxis()->SetTitle("Sum of |Delta Y| / #mum");
	histo1->GetXaxis()->SetTitle("#chi^{2} of Y fit");
	histName << "_graph";
	TGraph graph1 = histSaver->CreateDipendencyGraph(histName.str(), vecSumDeltaY, vecYChi2);
	graph1.Draw("AP");
	graph1.GetYaxis()->SetTitle("Sum of Delta Y / #mum");
	graph1.GetXaxis()->SetTitle("Chi2 in Y");

	histSaver->SaveHistogram(histo1);
	histSaver->SaveGraph(&graph1, histName.str(), "AP");


	histName.str("");
	if (nAlignmentStep == -1)
		histName << "hPreAlignment";
	else if (nAlignmentStep >= nAlignSteps - 1)
		histName << "hPostAlignment";
	else
		histName << "h_" << nAlignmentStep << "_Step";
	histName << "_PhiXY";
	TH2F *hPhiXY = histSaver->CreateScatterHisto(histName.str(), vecXPhi,vecYPhi,1024);
	hPhiXY->GetXaxis()->SetTitle("Phi X / degree");
	hPhiXY->GetYaxis()->SetTitle("Phi Y / degree");
	histSaver->SaveHistogram(hPhiXY);
}

void TAlignment::getFinalSiliconAlignmentResuluts() {
	Float_t maxChi2 = settings->getAlignment_chi2();
	cout << "*****\n*****get Final Silicon Resolution with a  maximum chi2 of " << maxChi2 <<"\n*****"<< endl;
	//set Resolutions
	//	setDetectorResolution(maxChi2);
	setSiliconDetectorResolution(maxChi2);
	getChi2Distribution(15);
	TResidual res = CheckDetectorAlignment(TPlaneProperties::XY_COR,3,0,0,true);
	CheckDetectorAlignment(TPlaneProperties::XY_COR,3,0,0,true,res);
}

void TAlignment::setSiliconDetectorResolution(Float_t maxChi2) {
	//get  something like a aprroximate Sigma with calculating the residual
	for (UInt_t plane = 0; plane < 4; plane++) {
		vector<UInt_t> vecRefPlanes;
		for (UInt_t i = 0; i < 4; i++)
			if (i != plane) vecRefPlanes.push_back(i);
		TResidual res = CheckDetectorAlignment(TPlaneProperties::XY_COR, plane, vecRefPlanes, false);
		res = CheckDetectorAlignment(TPlaneProperties::XY_COR, plane, vecRefPlanes, false, res);
		if(verbosity) cout<<"Results of silicon plane "<<plane<<" for resolution: "<<endl;
		if(verbosity) res.Print(1);
		align->setXResolution(res.getXSigma()/2, plane);
		align->setYResolution(res.getYSigma()/2, plane);
		align->setXMean(res.getDeltaXMean(),plane);
		align->setYMean(res.getDeltaYMean(),plane);
	}

	for (UInt_t plane = 0; plane < 4; plane++) {
		vector<UInt_t> vecRefPlanes;
		for (UInt_t i = 0; i < 4; i++)
			if (i != plane) vecRefPlanes.push_back(i);
		TResidual res = getResidual(TPlaneProperties::XY_COR, plane, vecRefPlanes,true,TResidual(),getClusterCalcMode(plane),chi2CalcMode,maxChi2);
	}
}

void TAlignment::CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, string refPlaneString, bool bPlot, bool bUpdateResolution, bool bChi2) {
	if (!bPlot && !bUpdateResolution) return;
	if (bPlot)
		if(verbosity>3)cout<<"Save Histograms: "<<  vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << " " << vecXObs.size() << " " << vecYObs.size() << endl;
	// define preName
	stringstream preName,postName;
	if (subjectPlane == 4) {
		preName << "hDiamond_";
		if (nDiaAlignmentStep == -1)
			preName << "PreAlignment";
		else if (nDiaAlignmentStep == nDiaAlignSteps)
			preName << "PostAlignment";
		else
			preName << nDiaAlignmentStep << "Step";
	} else {
		preName << "hSilicon_";
		if (nAlignmentStep == -1)
			preName << "PreAlignment";
		else if (nAlignmentStep == nAlignSteps)
			preName << "PostAlignment";
		else
			preName << nAlignmentStep << "_Step";
	}
	if(bChi2){
		postName<<"with_Chi2_cut_on_"<<settings->getAlignment_chi2();
	}
	else
		postName.str("");


	stringstream histName;
	if(verbosity){cout << "\nCreatePlots with " << preName.str() << " " << (subjectPlane!=4?nAlignmentStep:nDiaAlignmentStep) <<" Step" << flush;
	if (bUpdateResolution)
		cout << "\twith Alignment Resolution Update\n" << endl;
	else
		cout << endl;
	}
	Float_t xPredictionSigma=0;
	if (bUpdateResolution&&(cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    																  //SigmaOfPredictionX
			histName.str("");
			histName.clear();
			histName << preName.str() << "_SigmaOfPredictionX"<< "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
			TH1F* histo = histSaver->CreateDistributionHisto(histName.str(), vecXResPrediction, 512, HistogrammSaver::threeSigma);
			if (!histo)
				cerr<<"Could not CreateDistributionHisto: "<<histName.str()<<endl;
			else{
				histo->Draw("goff");
				histo->GetXaxis()->SetTitle("Delta X/ #mum");
				histo->GetYaxis()->SetTitle("Number of entries");
				xPredictionSigma = 	histSaver->GetMean(vecXResPrediction);
				cout<<"X-Prediction of sigma (*100): "<<xPredictionSigma*100<<"\thisto:"<<histo->GetMean()*100<<endl;

				if (bPlot) histSaver->SaveHistogram(histo);
				delete histo;
			}
	}

	if (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) {    																  //DistributionPlot DeltaX
		histName.str("");
		histName.clear();
		histName << preName.str() << "_Distribution_DeltaX"<< "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		TH1F* histo=0;
		if(TPlaneProperties::isDiamondPlane(subjectPlane)&&(nDiaAlignmentStep == nDiaAlignSteps)){
			Float_t pitchWidth = settings->getDiamondPitchWidth();
			histo = histSaver->CreateDistributionHisto(histName.str(), vecXDelta, 512, HistogrammSaver::manual,-1.1*pitchWidth,1.1*pitchWidth);
		}
		else
			histo = histSaver->CreateDistributionHisto(histName.str(), vecXDelta, 512, HistogrammSaver::threeSigma);
		if (!histo)
			cerr<<"Could not CreateDistributionHisto: "<<histName.str()<<endl;
		else{
			histo->Draw("goff");
			Float_t sigma = histo->GetRMS();
			Float_t fitWidth = sigma * 1.5;
			Float_t mean = histo->GetMean();
			cout << "Alignment for plane" << subjectPlane << endl;
			TF1* fitX=0;
			if(TPlaneProperties::isDiamondPlane(subjectPlane)){
				fitWidth = 3*sigma;
				fitX = new TF1("fit","[0]*TMath::Sqrt(TMath::Pi()/2)*[1]*(TMath::Erf(([2]+[3]-x)/TMath::Sqrt(2)/[1])+TMath::Erf(([3]-[2]+x)/TMath::Sqrt(2)/[1]))",mean-fitWidth,mean+fitWidth);
				fitX->FixParameter(3,settings->getDiamondPitchWidth()/2);//TODO
				fitX->SetParLimits(1,0,sigma);
				fitX->SetParNames("Integral","sigma of Gaus","position");
				fitX->SetParameter(2,0);
				fitX->SetParameter(1,0.1);
			}
			else{
				fitX = new TF1("fitGausX", "gaus", mean - fitWidth , mean + fitWidth);
				fitX->SetParameter(1,mean);
				fitX->SetParameter(2,sigma);
			}
			histo->Fit(fitX, "Q", "",mean-fitWidth, mean+fitWidth);
			Float_t xRes=0;

			if(TPlaneProperties::isDiamondPlane(subjectPlane)){
				xRes = TMath::Max(fitX->GetParameter(1),fitX->GetParError(2));
				mean = fitX->GetParameter(2);
			}
			else{
				mean = fitX->GetParameter(1);
				xRes = fitX->GetParameter(2);
			}

			if(xRes>0&&bUpdateResolution&&histo->GetEntries()>0){
				if(TPlaneProperties::isDiamondPlane(subjectPlane)){
					cout << "set Resolution via Gaus fit for diamond: " << xRes*100 << " with " << vecXDelta.size() << " Events" << endl;
					align->setXResolution(xRes,subjectPlane);
					align->setXMean(mean,subjectPlane);
				}
				else{
					cout << "\n\nset Resolution via Gaus-Fit: " << xRes*100 << " with " << vecXDelta.size() << " Events" << endl;
					cout << "xRes: "<<xRes*100<<endl;
					cout << "xPre: "<<xPredictionSigma*100<<endl;
					Float_t xres2 =xRes;
					if(xRes>xPredictionSigma){
						xres2 = xRes*xRes-xPredictionSigma*xPredictionSigma;
						xres2 = TMath::Sqrt(xres2);
					}
					else{
						xres2 = xres2/TMath::Sqrt2();
						cout<<" .... xRes < xPredictionSigma .....Update xRes to "<<xres2<<endl;
					}
					cout<< "xDet: "<<xres2*100 <<" = Sqrt("<<xRes*100<<"^2 + "<<xPredictionSigma*100<<"^2)"<<endl;
					align->setXResolution(xres2, subjectPlane);
					align->setXMean(mean,subjectPlane);
				}
				if(xRes<fitWidth/3&&TPlaneProperties::isDiamondPlane(subjectPlane))
					xRes = fitWidth/3;
				histo->GetXaxis()->SetRangeUser(mean-4 * xRes, mean + 4 * xRes);
			}
			histo->GetXaxis()->SetTitle("Delta X / #mum");
			histo->GetYaxis()->SetTitle("Number of entries");
			if (bPlot) histSaver->SaveHistogram(histo);
			delete fitX;
			delete histo;
		}
	}

	if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {   														  //ScatterPlot DeltaX vs Ypred
		histName.str("");
		histName.clear();
		histName << preName.str()<< "_ScatterPlot_YPred_vs_DeltaX"<< "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(),vecXDelta, vecYPred, 256);
		//    histo.Draw("goff");
		if(histo){
			histo->GetXaxis()->SetTitle("Y Predicted / #mum");
			histo->GetYaxis()->SetTitle("Delta X / #mum");
			histSaver->SaveHistogram(histo);
			delete histo;
		}
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecYPred);
		if(&graph==0) cerr<<"Could not create DipendencyGraph vecXDelta,vecYPred"<<endl;
		else {
			graph.Draw("APL");
			graph.GetXaxis()->SetTitle("predicted Y position / #mum");
			graph.GetYaxis()->SetTitle("delta X / #mum");
			histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		}
	}

	if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) && subjectPlane == TPlaneProperties::getDiamondPlane()) {    //ScatterPlot DeltaX vs Xpred
		histName.str("");
		histName.clear();
		histName << preName.str() << "_ScatterPlot_XPred_vs_DeltaX" << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXDelta,vecXPred, 256);
		if(!histo)
			cerr<<"Could not CreateScatterHisto: vecXDelta,vecXPred"<<endl;
		else{
			histo->GetXaxis()->SetTitle("X Predicted / #mum");
			histo->GetYaxis()->SetTitle("Delta X / #mum");
			histSaver->SaveHistogram(histo);
			delete histo;
		}
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXPred);
		if(&graph==0)
			cerr<<"Could not CreateDipendencyHisto: vecXDelta,vecXPred"<<endl;
		else{
			graph.Draw("APL");
			graph.GetXaxis()->SetTitle("predicted X position / #mum");
			graph.GetYaxis()->SetTitle("delta X / #mum");
			histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		}
	}

	if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) && subjectPlane == TPlaneProperties::getDiamondPlane()) {    //ScatterPlot DeltaX vs XMeas
		histName.str("");
		histName.clear();
		histName << preName.str()<< "_ScatterPlot_XMeasured_vs_DeltaX" << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		TH2F* histo = histSaver->CreateScatterHisto(histName.str(),vecXDelta, vecXMeasured,256);
		if(!histo)
			cerr<<"Could not CreateScatterHisto: "<<histName.str()<<endl;
		else{
			histo->GetXaxis()->SetTitle("X Measured / #mum");
			histo->GetYaxis()->SetTitle("Delta X / #mum");
			histSaver->SaveHistogram(histo);
			delete histo;
		}
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXMeasured);
		if(&graph==0)
			cerr<<"Could not CreateDipendencyGraph: "<<histName.str()<<endl;
		else{
			graph.Draw("APL");
			graph.GetXaxis()->SetTitle("measured X  / #mum");
			graph.GetYaxis()->SetTitle("delta X / #mum");
			histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		}
	}
	Float_t yPredictionSigma=0;
	if (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR) {    																  //DistributionPlot DeltaX
		histName.str("");
		histName.clear();
		histName << preName.str() << "_SigmaOfPredictionY"<< "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		TH1F* histo = histSaver->CreateDistributionHisto(histName.str(), vecYResPrediction, 512, HistogrammSaver::threeSigma);
		if (!histo)
			cerr<<"Could not CreateDistributionHisto: "<<histName.str()<<endl;
		else{
			histo->Draw("goff");
			histo->GetXaxis()->SetTitle("Delta X / #mum");
			histo->GetYaxis()->SetTitle("Number of entries");
			yPredictionSigma = 	histSaver->GetMean(vecYResPrediction);
			cout<<"Y-Prediction of sigma (*100): "<<yPredictionSigma*100<<"\thisto:"<<histo->GetMean()*100<<endl;
			if (bPlot) histSaver->SaveHistogram(histo);
			delete histo;
		}
	}

	if (subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {   											  //DistributionPlot DeltaY
		histName.str("");
		histName.clear();
		histName << preName.str()<<"_Distribution_DeltaY"<< "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3)cout<<"Save: "<<histName.str()<<flush;

		TH1F* histo = (TH1F*) histSaver->CreateDistributionHisto(histName.str(), vecYDelta, 512, HistogrammSaver::threeSigma);
		if(!histo)
			cerr<<"Could not CreateDistributionHisto: "<<histName.str()<<endl;
		else{
			histo->Draw("goff");
			Float_t sigma = histo->GetRMS();
			Float_t fitWidth = sigma *1.5;
			Float_t mean = histo->GetMean();
			TF1* fitGausY = new TF1("fitGausY", "gaus", mean-fitWidth, mean+fitWidth);
			fitGausY->SetParameter(1,mean);
			fitGausY->SetParameter(2,sigma);
			histo->Fit(fitGausY, "Q", "", mean-fitWidth, mean+fitWidth);
			TF1* fitGausY2 = (TF1*) fitGausY->Clone();
			fitGausY2->SetRange(mean-fitWidth*2,mean+fitWidth*2);
			fitGausY2->SetLineWidth(1);
			fitGausY2->SetLineStyle(3);
			histo->Fit(fitGausY2, "Q+", "", mean-fitWidth, mean+fitWidth);
			fitGausY2->SetParameter(0,fitGausY->GetParameter(0));
			fitGausY2->SetParameter(1,fitGausY->GetParameter(1));
			fitGausY2->SetParameter(2,fitGausY->GetParameter(2));
			Float_t yRes = fitGausY->GetParameter(2);
			mean = fitGausY->GetParameter(1);
			if (bUpdateResolution&&histo->GetEntries()>0 && yRes > 0) {
				cout << "\n\nset Y-Resolution via Gaus-Fit: " << yRes*100 << " with " << vecYDelta.size() << " Events" << endl;
				cout << "yRes: "<<yRes*100<<endl;
				cout << "yPre: "<<yPredictionSigma*100<<endl;
				Float_t yres2 =yRes;
				if(yRes>yPredictionSigma){
					yres2 = yRes*yRes-yPredictionSigma*yPredictionSigma;
					yres2 = TMath::Sqrt(yres2);
				}
				else{
					yres2 = yres2/TMath::Sqrt2();
					cout<<" .... yRes < yPredictionSigma .....Update yRes to "<<yres2<<endl;
				}
				cout<< "yDet: "<<yres2*100<<" = sqrt( "<<yRes*100<<"^2 + "<<yPredictionSigma*100<<"^2)"<<endl;
				align->setYResolution(yres2, subjectPlane);
				align->setYMean(mean,subjectPlane);
				histo->GetXaxis()->SetRangeUser(mean-4 * yRes,mean +4 * yRes);
			}
			histo->GetXaxis()->SetTitle("Delta Y / #mum");
			histo->GetYaxis()->SetTitle("Number of entries");
			if (bPlot) histSaver->SaveHistogram(histo);
			if(verbosity>3)cout<<" DONE"<<endl;
			delete fitGausY;
			delete histo;
		}
	}

	if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    									  //ScatterPlot DeltaY vs Xpred
		histName.str("");
		histName.clear();
		histName << preName.str() << "_ScatterPlot_XPred_vs_DeltaY" << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(),vecYDelta,  vecXPred, 256);
		if(!histo)
			cerr<<"Could not create "<<histName.str()<<endl;
		else{
			histo->GetXaxis()->SetTitle("X Predicted / #mum");
			histo->GetYaxis()->SetTitle("Delta Y / #mum");
			histSaver->SaveHistogram((TH2F*) histo->Clone());
			delete histo;
		}
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecYDelta, vecXPred);
		if(&graph==0)
			cerr<<"Could not CreateDipendencyGraph "<<histName.str()<<endl;
		else{
			graph.Draw("APL");
			graph.GetXaxis()->SetTitle("X Predicted / #mum");
			graph.GetYaxis()->SetTitle("Delta Y / #mum");
			histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		}
		if(verbosity>3)cout<<" DONE"<<endl;
	}

	if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR)) {   									  								  //ScatterHisto XObs vs YObs
		histName.str("");
		histName.clear();
		histName << preName.str() << "_ScatterPlot_XObs_vs_YObs" << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecYObs, vecXObs,512);
		if(!histo)
			cerr<<"Could not create CreateScatterHisto: "<<histName.str()<<endl;
		else{
			histo->GetXaxis()->SetTitle("XObs / #mum");
			histo->GetYaxis()->SetTitle("YObs / #mum");
			histSaver->SaveHistogram(histo);    //,histName.str());
			delete histo;
		}
		if(verbosity>3)cout<<"DONE"<<endl;
	}

	if (bPlot && nAlignmentStep > -1 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    								  //ScatterHisto DeltaX vs Chi2X
		histName.str("");
		histName.clear();
		histName << preName.str();
		histName << "_ScatterPlot_DeltaX_vs_Chi2X";
		histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXDelta, vecXChi2, 256);
		//    histo->Draw("goff");
		histo->GetXaxis()->SetTitle("Delta X / #mum");
		histo->GetYaxis()->SetTitle("Chi2 X");
		histSaver->SaveHistogram((TH2F*) histo->Clone());
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXChi2);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("#chi^{2} per NDF");
		graph.GetYaxis()->SetTitle("Delta X / #mum");
		histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		delete histo;
		if(verbosity>3)cout<<"DONE"<<endl;
	}

	if (bPlot && nAlignmentStep > -1 && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {   			  //ScatterHisto DeltaY vs Chi2Y
		histName.str("");
		histName << preName.str();
		histName << "_ScatterPlot_DeltaX_vs_Chi2X";
		histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecYDelta, vecYChi2, 256);
		//    histo->Draw("goff");
		histo->GetYaxis()->SetTitle("Sum of Delta Y / #mum");
		histo->GetXaxis()->SetTitle("Chi2 Y");

		histSaver->SaveHistogram((TH2F*) histo->Clone());
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecYDelta, vecYChi2);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("#chi^{2} per NDF");
		graph.GetYaxis()->SetTitle("Delta Y / #mum");
		histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		delete histo;
		if(verbosity>3)cout<<"DONE"<<endl;
	}

	if (bPlot && subjectPlane == 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    								  //predX vs deltaX <-> Diamond
		histName.str("");
		histName << preName.str();
		histName << "_ScatterPlot_XPred_vs_DeltaX";
		histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXDelta, vecXPred, 512);
		//    histo->Draw("goff");
		histo->GetXaxis()->SetTitle("X Predicted / #mum");
		histo->GetYaxis()->SetTitle("Delta X / #mum");

		histSaver->SaveHistogram((TH2F*) histo->Clone());
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXPred);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("X Predicted / #mum");
		graph.GetYaxis()->SetTitle("Delta X / #mum");

		histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		delete histo;
		if(verbosity>3)cout<<"DONE"<<endl;
	}

	if (bPlot ){																																  //hAngularDistribution
			histName.str("");histName.clear();
			histName<<preName.str()<<"_AngularDistribution_for_"<<subjectPlane<<"_with_"<<refPlaneString<<postName.str();
			if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
			TH2F *histo = histSaver->CreateScatterHisto(histName.str(),vecXPhi,vecYPhi,512);
			histo->GetXaxis()->SetTitle("PhiX / degree");
			histo->GetYaxis()->SetTitle("PhiY / degree");
			histSaver->SaveHistogram((TH2F*) histo->Clone());
			histName << "_graph";
			TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecYPhi, vecXPhi);
			graph.Draw("APL");
			graph.GetXaxis()->SetTitle("xPhi / degree");
			graph.GetYaxis()->SetTitle("yPhi / degree");
			histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
			delete histo;
			if(verbosity>3)cout<<"DONE"<<endl;
		}

	if (bPlot && subjectPlane == 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    //DeltaX vs ClusterSize
		histName.str("");
		histName << preName.str() << "_ScatterPlot_ClusterSize_vs_DeltaX_-_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName.str();
		if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
		TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXDelta, vecClusterSize, 512);
		histo->Draw("goff");
		histo->GetXaxis()->SetTitle("Cluster Size");
		histo->GetYaxis()->SetTitle("Delta X / #mum");
		histSaver->SaveHistogram((TH2F*) histo->Clone());
		histName << "_graph";
		TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecClusterSize);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Cluster Size");
		graph.GetYaxis()->SetTitle("Delta X / #mum");
		histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
		delete histo;
		if(verbosity>3)cout<<"DONE"<<endl;
	}

}

void TAlignment::DoEtaCorrectionSilicon(UInt_t correctionStep) {
	cout << "****************************************" << endl;
	cout << "Do Eta correction for all silicon planes " << correctionStep << endl;
	cout << "****************************************" << endl;
	vector<TH1F*> histoStripDistribution, histoStripDistributionFlattned;
	vector<TH1F *> vecHEta;
	stringstream fileName;
	fileName << "etaIntegral_Step" << correctionStep << "." << settings->getRunNumber() << ".root";
	TFile* correctedEtaFile = new TFile(fileName.str().c_str(), "RECREATE");
	cout << "create Histos..." << endl;
	for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
		stringstream histoTitle;
		histoTitle << "hPredictedStripPosition" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
		histoStripDistribution.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, -0.501, 0.501));
		histoTitle << "_flattend";
		histoStripDistributionFlattned.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, -0.501, 0.501));
		histoTitle.str("");
		histoTitle.clear();
		histoTitle << "hCorrectedEtaDistribution" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
		vecHEta.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, 0, 1));
	}

	cout << "fill first strip hit histo" << events.size() << endl;

	for (nEvent = 0; nEvent < this->events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent, events.size());
		myTrack->setEvent(&events.at(nEvent));

		for (UInt_t subjectPlane = 0; subjectPlane < TPlaneProperties::getNSiliconPlanes(); subjectPlane++) {
			//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
			//				continue;
			vector<UInt_t> refPlanes;
			for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
				if (subjectPlane != refPlane) refPlanes.push_back(refPlane);
			TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, TCluster::corEta, false);
			Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane * 2, pred->getPositionX(), pred->getPositionY());
			Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane * 2 + 1, pred->getPositionX(), pred->getPositionY());
			UInt_t stripMiddleX = (UInt_t) (predictedStripPositionX + 0.5);
			Float_t deltaX = predictedStripPositionX - stripMiddleX;
			UInt_t stripMiddleY = (UInt_t) (predictedStripPositionY + 0.5);
			Float_t deltaY = predictedStripPositionY - stripMiddleY;
			//			cout<<nEvent<<": "<<subjectPlane<<"Fill "<<deltaX<<" "<<deltaY<<endl;
			histoStripDistribution.at(subjectPlane * 2)->Fill(deltaX);
			histoStripDistribution.at(subjectPlane * 2 + 1)->Fill(deltaY);
		}
	}
	vector<UInt_t> vecMinEntries;
	cout << "Minimal Entries in a bin of historgram:" << endl;
	for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
		TH1F* histo = histoStripDistribution.at(det);
		Int_t minBin = histo->GetMinimumBin();
		Int_t nMinEntries = histo->GetBinContent(minBin);
		cout << endl << det << ": " << minBin << " " << nMinEntries << "\t";
		vecMinEntries.push_back(nMinEntries);
	}
	cout << "\n\n" << endl;

	vector<UInt_t> vecEventNo[9];
	cout << "create flattened strip hit histo " << events.size() << endl;
	for (nEvent = 0; nEvent < this->events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent, events.size());
		myTrack->setEvent(&events.at(nEvent));

		for (UInt_t subjectPlane = 0; subjectPlane < TPlaneProperties::getNSiliconPlanes(); subjectPlane++) {
			vector<UInt_t> refPlanes;
			for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
				if (subjectPlane != refPlane) refPlanes.push_back(refPlane);
			TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, TCluster::corEta, false);

			Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane * 2, pred->getPositionX(), pred->getPositionY());
			Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane * 2 + 1, pred->getPositionX(), pred->getPositionY());
			UInt_t stripMiddleX = (UInt_t) (predictedStripPositionX + 0.5);
			Float_t deltaX = predictedStripPositionX - stripMiddleX;
			UInt_t stripMiddleY = (UInt_t) (predictedStripPositionY + 0.5);
			Float_t deltaY = predictedStripPositionY - stripMiddleY;

			Int_t binX = histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaX);
			Int_t binY = histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaY);
			if (histoStripDistributionFlattned.at(subjectPlane * 2)->GetBinContent(binX) < vecMinEntries.at(subjectPlane * 2)) {
				vecEventNo[subjectPlane * 2].push_back(nEvent);
				histoStripDistributionFlattned.at(subjectPlane * 2)->Fill(deltaX);
			}
			if (histoStripDistributionFlattned.at(subjectPlane * 2 + 1)->GetBinContent(binY) < vecMinEntries.at(subjectPlane * 2 + 1)) {
				vecEventNo[subjectPlane * 2 + 1].push_back(nEvent);
				histoStripDistributionFlattned.at(subjectPlane * 2 + 1)->Fill(deltaY);
			}
		}
	}
	for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
		cout << "save histogram: " << det << "  " << histoStripDistributionFlattned.at(det)->GetTitle() << "  " << histoStripDistribution.at(det)->GetTitle() << endl;
		histSaver->SaveHistogram(histoStripDistributionFlattned.at(det));
		correctedEtaFile->Add(histoStripDistributionFlattned.at(det));
		histSaver->SaveHistogram(histoStripDistribution.at(det));
		correctedEtaFile->Add(histoStripDistribution.at(det));
	}

	cout << "\n\ncreate eta correction histo" << endl;
	for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
		for (UInt_t i = 0; i < vecEventNo[det].size(); i++) {
			nEvent = vecEventNo[det].at(i);
			eventReader->LoadEvent(nEvent);

			if (!eventReader->useForAlignment() && !eventReader->useForAnalysis()) continue;
			Float_t eta = eventReader->getCluster(det, 0).getEta();
			vecHEta.at(det)->Fill(eta);

		}

		stringstream histName;
		histName << "hEtaIntegral" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
		;
		UInt_t nBins = vecHEta.at(det)->GetNbinsX();
		TH1F *histo = new TH1F(histName.str().c_str(), histName.str().c_str(), nBins, 0, 1);
		Int_t entries = vecHEta.at(det)->GetEntries();
		entries -= vecHEta.at(det)->GetBinContent(0);
		entries -= vecHEta.at(det)->GetBinContent(nBins + 1);
		Int_t sum = 0;
		for (UInt_t bin = 1; bin < nBins + 1; bin++) {
			Int_t binContent = vecHEta.at(det)->GetBinContent(bin);
			sum += binContent;
			Float_t pos = vecHEta.at(det)->GetBinCenter(bin);
			histo->Fill(pos, (Float_t) sum / (Float_t) entries);
		}
		cout << "save " << vecHEta.at(det)->GetTitle() << " " << vecHEta.at(det)->GetEntries() << endl;
		histSaver->SaveHistogram(vecHEta.at(det));
		correctedEtaFile->Add(vecHEta.at(det));
		cout << "save " << histo->GetTitle() << " " << histo->GetEntries() << endl;
		histSaver->SaveHistogram(histo);
		correctedEtaFile->Add(histo->Clone());
	}
	correctedEtaFile->Write();
	correctedEtaFile->Close();
}

/**
 * @brief: Print Informations for all Events smaller than maxEvent, starting with startEvent
 * @param maxEvent  maximum Event No to print event
 * @param:  startEvent  first Event where to start with printing Event Information
 *
 */
void TAlignment::PrintEvents(UInt_t maxEvent, UInt_t startEvent) {
	if (maxEvent == 0) maxEvent = eventReader->GetEntries();
	if (maxEvent > eventReader->GetEntries()) maxEvent = eventReader->GetEntries();
	cout << "\n\n\n\n\n\nPrintEVENTS" << maxEvent << "\n\n\n" << flush;
	for (UInt_t event = startEvent; event < maxEvent; event++) {
		eventReader->LoadEvent(event);
		cout << event << ":\t" << flush;
		eventReader->getEvent()->Print(1);
	}
	cout << "\n\n\n\n\n\n" << flush;

}

void TAlignment::clearMeasuredVectors() {
	vecXPred.clear();
	vecYPred.clear();
	vecXObs.clear();
	vecYObs.clear();
	vecXDelta.clear();
	vecYDelta.clear();
	vecXChi2.clear();
	vecYChi2.clear();
	vecXPhi.clear();
	vecYPhi.clear();
	vecXMeasured.clear();
	vecClusterSize.clear();
	vecYMeasured.clear();
	vecXResPrediction.clear();
	vecYResPrediction.clear();
}

