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
	histSaver->SetNumberOfEvents(eventReader->GetEntries());//todo anpassen
	sys->cd("..");
	cout<<"end initialise"<<endl;
	alignmentPercentage=0.1;
	Float_t stripSize=1.;// 50./10000.;//mu m
	detectorD0Z = 0.725/stripSize; // by definition in cm
	detectorD1Z = 1.625/stripSize; // by definition in cm
	detectorD2Z = 18.725/stripSize; // by definition in cm
	detectorD3Z = 19.625/stripSize; // by definition in cm
	detectorDiaZ = 10.2/stripSize; // by definition in cm
	this->runNumber=runNumber;
	verbosity=1;
	nAlignSteps=1;
	res_keep_factor = 2;//todo anpassen
	align=NULL;
	myTrack=NULL;
	alignmentSteps=5;
	nAlignmentStep=-1;

}

TAlignment::~TAlignment() {
	// TODO Auto-generated destructor stub
	cout<<"TAlignment deconstructor"<<endl;
	delete histSaver;
	delete eventReader;
	sys->cd("..");
}

void TAlignment::setSettings(TSettings* settings){
	this->settings=settings;
}


/**
 *
 * @param nEvents
 * @param startEvent
 * @todo	chekc of fiduccial cut and other stuff...
 */
void TAlignment::createEventVectors(UInt_t nEvents, UInt_t startEvent){
	if (nEvents == 0) nEvents = eventReader->GetEntries()-startEvent;
	if (nEvents+startEvent>eventReader->GetEntries())nEvents=eventReader->GetEntries()-startEvent;
	int noHitDet=0;
	int falseClusterSizeDet=0;
	int noHitDia=0;
	int falseClusterSizeDia=0;
	int nCandidates=0;
	int nScreened=0;
	cout << "CREATING VECTOR OF VALID EVENTS..." << endl;
	
	for (nEvent = startEvent; nEvent < nEvents+startEvent; nEvent++) {
		TRawEventSaver::showStatusBar(nEvent-startEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
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
		this->events.push_back(*eventReader->getEvent());
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
int TAlignment::Align(UInt_t nEvents,UInt_t startEvent)
{
	if(tracks.size()==0)
		createEventVectors(nEvents,startEvent);

	if (verbosity){
		cout<<"\n\n\nTAlignment::Align:Starting \""<<histSaver->GetPlotsPath()<<"\""<<endl;
		cout<<"\t\t"<<events.size()<<"\t";
		cout << "\t\t "<<eventReader<<" ."<<endl;
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
		cout<<"TAlignment::Align::create new TTrack"<<endl;
		myTrack = new TTrack(align);
		cout<<"TAlignment::Align::created new TTrack"<<endl;
	}
	cout<<"~"<<endl;
	cout<<"~"<<endl;
	cout<<"~"<<endl;

	nAlignmentStep=-1;

	CheckDetectorAlignment(TPlane::XY_COR,1,0,3,true);
	CheckDetectorAlignment(TPlane::XY_COR,2,0,3,true);
	CheckDetectorAlignment(TPlane::XY_COR,3,1,2,true);
	AlignDetector(TPlane::XY_COR,1,0,0,false);
	AlignDetector(TPlane::XY_COR,2,0,0,false);
	AlignDetector(TPlane::XY_COR,3,0,0,false);

	TResidual res103=CheckDetectorAlignment(TPlane::XY_COR,1,0,3,false);
	TResidual res203=CheckDetectorAlignment(TPlane::XY_COR,2,0,3,false);
	TResidual res312=CheckDetectorAlignment(TPlane::XY_COR,3,1,2,false);
	cout<<endl;
	for(nAlignmentStep=0;nAlignmentStep<alignmentSteps;nAlignmentStep++){
		cout<<"\n\n\nALIGNMENT STEP:\t"<<nAlignmentStep+1<<" of "<<alignmentSteps<<"\n\n"<<endl;
		printf("\n1 with 0 and 3:\n");
		res103=AlignDetector(TPlane::XY_COR,1,0,3,false,res103);
		printf("\n2 with 0 and 3:\n");
		res203=AlignDetector(TPlane::XY_COR,2,0,3,false,res203);
		printf("\n3 with 1 and 2:\n");
		res312=AlignDetector(TPlane::XY_COR,3,1,2,false,res312);
		cout<<endl;
	}
	this->CheckDetectorAlignment(TPlane::XY_COR,1,0,3,true,res103);
	this->CheckDetectorAlignment(TPlane::XY_COR,2,0,3,true,res203);
	this->CheckDetectorAlignment(TPlane::XY_COR,3,1,2,true,res312);

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

	cout<<"**********************************************"<<endl;
	cout<<"**********************************************"<<endl;
	cout<<"******** TAlignment::Align:RESULTS ***********"<<endl;
	cout<<"**********************************************"<<endl;
	cout<<"**********************************************"<<endl;
	align->PrintResults(1);
	printf("\n\n\n");
	for(int k=0;k<(int)vecRes103.size();k++){
		printf("\n\n%d of %d Residuals",k+1,(int)vecRes103.size());
		vecRes103.at(k).Print(1);
	}
	return 1;
}

/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::XY_COR,subjectPlane,refPlane1,refPlane2);
}


/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::X_COR,subjectPlane,refPlane1,refPlane2);
}


/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlane::Y_COR,subjectPlane,refPlane1,refPlane2);
}

/**
 *
 * @param cor
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 * @param bPlot
 */
TResidual TAlignment::AlignDetector(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot,TResidual resOld)
{
	//cout<<"\n\nTAlignment::AlignDetector\n\t align "<<TPlane::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with Plane "<<refPlane1<<" and "<<refPlane2<<endl;
	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
		return resOld;
	}

	TResidual res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot);
	//res.Print(0);
	if(refPlane1==refPlane2){//todo ersetzen durch "echtes alignment?"
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR)
			align->AddToXOffset(subjectPlane,res.getXMean());
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR)
			align->AddToYOffset(subjectPlane,res.getYMean());
		return res;
	}
	else{
		Float_t x_offset = res.getXOffset();
		Float_t phix_offset = res.getPhiXOffset();
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR){
			align->AddToXOffset(subjectPlane,x_offset);
			align->AddToPhiXOffset(subjectPlane,phix_offset);
			if(verbosity)
				printf("\tUpdate  Plane %d, X: %2.6f, PhiX: %2.6f\n",subjectPlane,x_offset,phix_offset);
		}
		if(subjectPlane==1)vecRes103.push_back(res);
		res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot);
		if(subjectPlane==1)vecRes103.push_back(res);
		Float_t y_offset = res.getYOffset();
		Float_t phiy_offset = res.getPhiYOffset();
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR){
			align->AddToYOffset(subjectPlane,y_offset);
			align->AddToPhiYOffset(subjectPlane,phiy_offset);
			if(verbosity)
				printf("\tUpdate Plane %d, Y: %2.6f, PhiY: %2.6f\n",subjectPlane,y_offset,phiy_offset);
		}
		//
		//res.Print();

	}
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
TResidual TAlignment::getResidual(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot,TResidual resOld){
	if(verbosity)cout<<"TAlignment::getResidual of Plane "<<subjectPlane<<" with "<<refPlane1<<" and "<<refPlane2<<", plotting"<<bPlot<<" with "<<alignmentPercentage<<endl;
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

	for(UInt_t nEvent=0; nEvent<events.size(); nEvent++)
	{
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));
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
	TResidual res = calculateResidual(cor,vecXPred,vecDeltaX,vecYPred,vecDeltaY,resOld);
	if(bPlot){//todo
		stringstream histName;
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR){
			//DistributionPlot DeltaX
			histName.str("");
			if(nAlignmentStep==-1)histName<<"hPreAlignment";
			else if(nAlignmentStep==alignmentSteps)histName<<"hPostAlignment";
			else histName<<"h_"<<nAlignmentStep<<"_Step";
			histName<<"_DistributionPlot_DeltaX";
			histName<<"_-_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2;
			histSaver->SaveHistogramPNG((TH1F*)histSaver->CreateDistributionHisto(histName.str(),vecDeltaX,4096).Clone());
		}
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR){
			//DistributionPlot DeltaY
			histName.str("");
			if(nAlignmentStep==-1)histName<<"hPreAlignment";
						else if(nAlignmentStep==alignmentSteps)histName<<"hPostAlignment";
						else histName<<"h_"<<nAlignmentStep<<"_Step";
			histName<<"_DistributionPlot_DeltaY";
			histName<<"_-_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2;
			TH1F* histo = (TH1F*)histSaver->CreateDistributionHisto(histName.str(),vecDeltaY,4096).Clone();
			histSaver->SaveHistogramPNG(histo);
			histSaver->SaveHistogramROOT(histo);
			delete histo;
		}
		if(cor==TPlane::XY_COR||cor==TPlane::X_COR){
			//DistributionPlot DeltaX vs Ypred
			histName.str("");
			if(nAlignmentStep==-1)histName<<"hPreAlignment";
						else if(nAlignmentStep==alignmentSteps)histName<<"hPostAlignment";
						else histName<<"h_"<<nAlignmentStep<<"_Step";
			histName<<"_ScatterPlot_YPred_vs_DeltaX";
			histName<<"_-_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2;
			TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecYPred,vecDeltaX,4096);
			histSaver->SaveHistogramPNG(&histo);
			histSaver->SaveHistogramROOT(&histo);
			histName<<"_graph";
			TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecYPred);
			histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
		}
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR){
			//DistributionPlot DeltaY vs Xpred
			histName.str("");
			if(nAlignmentStep==-1)histName<<"hPreAlignment";
			else if(nAlignmentStep==alignmentSteps)histName<<"hPostAlignment";
			else histName<<"h_"<<nAlignmentStep<<"_Step";
			histName<<"_ScatterPlot_XPred_vs_DeltaY";
			histName<<"_-_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2;
			histSaver->SaveHistogramPNG((TH2F*)histSaver->CreateScatterHisto(histName.str(),vecXPred,vecDeltaY,4096).Clone());
			histName<<"_graph";
			TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecYPred);
			histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
		}
		if(cor==TPlane::XY_COR||cor==TPlane::Y_COR){
			//DistributionPlot DeltaY vs Xpred
			histName.str("");
			if(nAlignmentStep==-1)histName<<"hPreAlignment";
			else if(nAlignmentStep==alignmentSteps)histName<<"hPostAlignment";
			else histName<<"h_"<<nAlignmentStep<<"_Step";
			histName<<"_ScatterPlot_YPred_vs_YObs";
			histName<<"_-_Plane_"<<subjectPlane<<"_with"<<refPlane1<<"_and_"<<refPlane2;
			TH2F histo =histSaver->CreateScatterHisto(histName.str(),vecYObs,vecYPred);
			histo.GetXaxis()->SetTitle("YObs");
			histo.GetYaxis()->SetTitle("YPred");
			histSaver->SaveHistogram((TH2F*)histo.Clone());//,histName.str());
		}
	}
	return res;
}


/**
 * calculateResidual if there is no residuals calculatet to cut on the res_keep_factor;
 * this funcition opens the other calculateResidual function but uses a TResidual res
 * which is using all items...
 */
TResidual TAlignment::calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t> yPred,vector<Float_t> deltaY){
	TResidual res;
	res.SetTestResidual();
	if(verbosity>2)cout<<"\t"<<deltaX.size()<<" "<<deltaY.size()<<" "<< xPred.size()<<" "<<yPred.size()<<endl;
	return calculateResidual(cor,xPred,deltaX,yPred,deltaY,res);
}

TResidual TAlignment::CheckDetectorAlignment(TPlane::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot,TResidual resOld)
{
	cout<<"\n\nTAlignment::checkDetectorAlignment\n\t check "<<TPlane::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with Plane "<<refPlane1<<" and "<<refPlane2<<endl;
	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
		return TResidual(true);
	}

	TResidual res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot,resOld);
	res.Print();
	return res;
}

/**
 * calculating residual:
 * uses the input TResidual res to use only events with
 * resxtest= TMath::Abs(deltaX.at(i)-res.resXMean)/res.resXSigm < res_keep_factor
 * resytest= TMath::Abs(deltaY.at(i)-res.resYMean)/res.resYSigma < res_keep_factor
 *
 * @param 	cor		which coordinate (X,Y,X&Y) should be used
 * @param	xPred	vector of predicted x Positions
 * @param	deltaX	vector of difference between predicted and measured x Positions
 * @param	yPred	vector of predicted y Position
 * @param	deltaY	vector of difference between predicted and measured y Positions
 * @param	res		current residuals for use if res_keep_factor
 *
 * @return calculated Residual
 */
TResidual TAlignment::calculateResidual(TPlane::enumCoordinate cor,vector<Float_t>xPred,vector<Float_t> deltaX,vector<Float_t>yPred,vector<Float_t> deltaY,TResidual res)
{
	//cout<<"\nTAlignment::calculateResidual"<<endl;
	TResidual residual;
	Float_t resxtest;
	Float_t resytest;
	if(verbosity>2)cout<<"\tcalculate Residual "<<res_keep_factor<<endl;
	if(verbosity>2)cout<<"\t"<<deltaX.size()<<" "<<deltaY.size()<<" "<< xPred.size()<<" "<<yPred.size()<<endl;
	for(UInt_t i=0;i<deltaX.size();i++){
		resxtest= TMath::Abs(deltaX.at(i)-res.getXMean())/res.getXSigma();
		resytest= TMath::Abs(deltaY.at(i)-res.getYMean())/res.getYSigma();
		//only add if restest is smaller than res_keep_factor
		if((cor==TPlane::X_COR)&&resxtest<res_keep_factor){
			residual.addDataPoint(deltaX.at(i),xPred.at(i),deltaY.at(i),yPred.at(i));
		}//end if
		else if((cor==TPlane::X_COR)&&resytest<res_keep_factor){
			residual.addDataPoint(deltaX.at(i),xPred.at(i),deltaY.at(i),yPred.at(i));
		}//end else if
		else if((cor==TPlane::XY_COR)&&resxtest<res_keep_factor&&resytest<res_keep_factor){
			residual.addDataPoint(deltaX.at(i),xPred.at(i),deltaY.at(i),yPred.at(i));
		}//end else if
	}//end for loop

	if(verbosity>0)	cout<<"\n\tused "<<residual.getUsedTracks()<<" Tracks"<<endl;
	if(verbosity>0)	cout<<"\tX: "<<std::setprecision(4)<<residual.getXMean()<<"+/-"<<residual.getXSigma()<<endl;
	if(verbosity>0)	cout<<"\tY: "<<residual.getYMean()<<"+/-"<<residual.getYSigma()<<"\n"<<endl;
	//set values
	return residual;
}


/**
 * @brief: Print Informations for all Events smaller than maxEvent, starting with startEvent
 * @param	maxEvent	maximum Event No to print event
 * @param: 	startEvent	first Event where to start with printing Event Information
 *
 */
void TAlignment::PrintEvents(UInt_t maxEvent,UInt_t startEvent){
	if(maxEvent==0)maxEvent=eventReader->GetEntries();
	if(maxEvent>eventReader->GetEntries())maxEvent=eventReader->GetEntries();
	cout<<"\n\n\n\n\n\nPrintEVENTS"<<maxEvent<<"\n\n\n"<<flush;
	for(UInt_t event=startEvent;event<maxEvent;event++){
		eventReader->LoadEvent(event);
		cout<<event<<":\t"<<flush;
		eventReader->getEvent()->Print(1);
	}
	cout<<"\n\n\n\n\n\n"<<flush;

}
