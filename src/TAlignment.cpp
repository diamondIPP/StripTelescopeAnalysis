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
	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<sys->pwd();
	filepath<<"/selectionData."<<runNumber<<".root";
	cout<<"OPEN eventReader with file \""<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str(),settings->getRunNumber());
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
	alignmentPercentage=settings->getAlignment_training_track_fraction();
	Float_t stripSize=1.;// 50./10000.;//mu m
	detectorD0Z = 0.725/stripSize; // by definition in cm
	detectorD1Z = 1.625/stripSize; // by definition in cm
	detectorD2Z = 18.725/stripSize; // by definition in cm
	detectorD3Z = 19.625/stripSize; // by definition in cm
	detectorDiaZ = 10.2/stripSize; // by definition in cm
	verbosity=1;
	res_keep_factor = settings->getRes_keep_factor();
	cout<<"Res Keep factor is set to "<<res_keep_factor<<endl;
	align=NULL;
	myTrack=NULL;
	nAlignmentStep=-1;
	nAlignSteps=5;
	nDiaAlignmentStep=-1;
	nDiaAlignSteps=5;

	diaCalcMode = TCluster::corEta;
	silCalcMode = TCluster::maxValue;

}

TAlignment::~TAlignment() {
	// TODO Auto-generated destructor stub
	cout<<"TAlignment deconstructor"<<endl;

	if(myTrack)delete myTrack;
	if(histSaver)delete histSaver;
	//todo this is not correct but otherwise it stops working
//	if(eventReader)delete eventReader;
	sys->cd("..");
}

void TAlignment::setSettings(TSettings* settings){
	this->settings=settings;
}

/**
 * initialise the variable align and set the Z Offsets
 */
void TAlignment::initialiseDetectorAlignment(){
	if(align==NULL){
		align = new TDetectorAlignment();
		cout<<"TAlignment::Align::Detectoralignment did not exist, so created new DetectorAlignment"<<endl;
		align->SetZOffset(0,detectorD0Z);
		align->SetZOffset(1,detectorD1Z);
		align->SetZOffset(2,detectorD2Z);
		align->SetZOffset(3,detectorD3Z);
		align->SetZOffset(4,detectorDiaZ);
	}
}
/**
 *
 * @param nEvents
 * @param startEvent
 * @todo	chekc of fiduccial cut and other stuff...
 */
void TAlignment::createEventVectors(UInt_t nEvents, UInt_t startEvent){
	initialiseDetectorAlignment();
	if (nEvents == 0) nEvents = eventReader->GetEntries()-startEvent;
	if (nEvents+startEvent>eventReader->GetEntries())nEvents=eventReader->GetEntries()-startEvent;
	int noHitDet=0;
//	int falseClusterSizeDet=0;
	//int noHitDia=0;
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
		if(eventReader->useForAlignment()){
			nCandidates++;
			this->events.push_back(*eventReader->getEvent());
		}
	}
	align->addEventIntervall(startEvent,startEvent+nEvents);
	align->setNUsedEvents(events.size());
	align->setAlignmentTrainingTrackFraction(settings->getAlignment_training_track_fraction());
	align->setRunNumber(settings->getRunNumber());

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
/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorXY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlaneProperties::XY_COR,subjectPlane,refPlane1,refPlane2);
}


/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorX(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlaneProperties::X_COR,subjectPlane,refPlane1,refPlane2);
}


/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::AlignDetectorY(UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2)
{
	AlignDetector(TPlaneProperties::Y_COR,subjectPlane,refPlane1,refPlane2);
}


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

int TAlignment::Align(UInt_t nEvents,UInt_t startEvent)
{
	if (verbosity){
		cout<<"\n\n\nTAlignment::Align:Starting \""<<histSaver->GetPlotsPath()<<"\""<<endl;
		cout<<"\t\t"<<events.size()<<"\t";
		cout << "\t\t "<<eventReader<<" ."<<endl;
	}

	initialiseDetectorAlignment();
	if(events.size()==0)
		createEventVectors(nEvents,startEvent);
	if(myTrack==NULL){
		cout<<"TAlignment::Align::create new TTrack"<<endl;
		myTrack = new TTrack(align);
		cout<<"TAlignment::Align::created new TTrack"<<endl;
		for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
			myTrack->setEtaIntergral(det,eventReader->getEtaIntegral(det));
	}
	AlignSiliconPlanes();
	AlignDiamondPlane();
	align->PrintResults(1);

	this->saveAlignment();

	return 1;
}

TResidual TAlignment::AlignStripDetector(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane,vector<UInt_t> vecRefPlanes,bool bPlot,TResidual resOld){
	if(verbosity){
		cout<<"TAlignment::AlignStripDetector::\taligning Plane "<<subjectPlane<<TPlaneProperties::getCoordinateString(cor)<<" with "<<vecRefPlanes.size()<< " Planes: ";

		for(UInt_t i=0;i<vecRefPlanes.size();i++)
			cout<<vecRefPlanes.at(i)<<" ";
		cout<<"\t";
		if(bPlot)cout<<" plots are created\t";
		if(resOld.isTestResidual())
			cout<<"resOld is a testResidual"<<endl;
		else{
			cout<<endl;
			resOld.Print(2);
		}
		for(UInt_t i=0;i<vecRefPlanes.size();i++)
			if(vecRefPlanes.at(i)==subjectPlane)cerr<<"Plane "<<subjectPlane<<" is used as a reference Plane and as a subject Plane"<<endl;
	}
	//get residual
	TResidual res = this->getStripResidual(cor,subjectPlane,vecRefPlanes,false,bPlot,resOld);
	Float_t x_offset = res.getXOffset();
	Float_t phix_offset = res.getPhiYOffset();
	Float_t y_offset = res.getYOffset();
	Float_t phiy_offset = res.getPhiXOffset();

	if(verbosity)printf("Correction Values: X: %2.6f,  PhiX: %2.6f,   Y: %2.6f,  PhiY: %2.6f\n",x_offset,phix_offset,y_offset,phiy_offset);
	//save corrections to alignment
	if(vecRefPlanes.size()==1){
		if(TPlaneProperties::X_COR==cor||TPlaneProperties::XY_COR==cor)
			align->AddToXOffset(subjectPlane,x_offset);
		if(TPlaneProperties::Y_COR==cor||TPlaneProperties::XY_COR==cor)
			align->AddToYOffset(subjectPlane,y_offset);
	}
	else{
		if(subjectPlane==4&&nDiaAlignmentStep==0){
			align->AddToXOffset(subjectPlane,x_offset);
			return res;
		}
		if(TPlaneProperties::X_COR==cor||TPlaneProperties::XY_COR==cor){
			align->AddToXOffset(subjectPlane,x_offset);
			align->AddToPhiXOffset(subjectPlane,phix_offset);
		}
		if(TPlaneProperties::Y_COR==cor||TPlaneProperties::XY_COR==cor){
			align->AddToYOffset(subjectPlane,y_offset);
			align->AddToPhiYOffset(subjectPlane,phiy_offset);
		}
	}
	return res;
}

/**
 *
 * @param cor
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 * @param bPlot
 * @todo upschreiben dass die andere funktion aufgreufen wird...
 */
TResidual TAlignment::AlignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot,TResidual resOld)
{
	//cout<<"\n\nTAlignment::AlignDetector\n\t align "<<TPlaneProperties::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with Plane "<<refPlane1<<" and "<<refPlane2<<endl;
	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
		return resOld;
	}

	TResidual res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot,resOld);
	if(refPlane1==refPlane2){//todo ersetzen durch "echtes alignment?"
		if(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)
			align->AddToXOffset(subjectPlane,res.getXMean());
		if(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)
			align->AddToYOffset(subjectPlane,res.getYMean());
		return res;
	}
	else{
		Float_t x_offset = res.getXOffset();
		Float_t phix_offset = res.getPhiXOffset();
		if(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR){
			align->AddToXOffset(subjectPlane,x_offset);
			align->AddToPhiXOffset(subjectPlane,phix_offset);
			if(verbosity)
				printf("\tUpdate  Plane %d, X: %2.6f, PhiX: %2.6f\n",subjectPlane,x_offset,phix_offset);
		}
//		if(resOld.isTestResidual())
//			res=getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot,resOld);
//		else{
//			resOld = checkDetectorAlignment(TPlaneProperties::XY_COR,subjectPlane,refPlane1)
//			res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot,res);
//		}
		Float_t y_offset = res.getYOffset();
		Float_t phiy_offset = res.getPhiYOffset();
		if(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR){
			align->AddToYOffset(subjectPlane,y_offset);
			align->AddToPhiYOffset(subjectPlane,phiy_offset);
			if(verbosity)
				printf("\tUpdate Plane %d, Y: %2.6f, PhiY: %2.6f\n",subjectPlane,y_offset,phiy_offset);
		}
		//

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
TResidual TAlignment::AlignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld)
{
	if(verbosity)cout<<"TAlignment::AlignDetector::\taligning Plane "<<subjectPlane<<TPlaneProperties::getCoordinateString(cor)<<" with "<<vecRefPlanes.size()<< " Planes: ";
	for(UInt_t i=0;i<vecRefPlanes.size();i++)
		if(verbosity)cout<<vecRefPlanes.at(i)<<" ";
	if(verbosity)cout<<"\t";
	if(bPlot)if(verbosity)cout<<" plots are created\t";
	if(resOld.isTestResidual())
		if(verbosity)cout<<"resOld is a testResidual"<<endl;
	else{
		if(verbosity)cout<<endl;
		if(verbosity)resOld.Print(2);
	}
	for(UInt_t i=0;i<vecRefPlanes.size();i++)
		if(vecRefPlanes.at(i)==subjectPlane)cerr<<"Plane "<<subjectPlane<<" is used as a reference Plane and as a subject Plane"<<endl;

	//get residual
	TResidual res = this->getResidual(cor,subjectPlane,vecRefPlanes,bPlot,resOld);
	Float_t x_offset = res.getXOffset();
	Float_t phix_offset = res.getPhiXOffset();
	Float_t y_offset = res.getYOffset();
	Float_t phiy_offset = res.getPhiYOffset();

	printf("Correction Values: X: %2.6f,  PhiX: %2.6f,   Y: %2.6f,  PhiY: %2.6f\n",x_offset,phix_offset,y_offset,phiy_offset);
	//save corrections to alignment
	if(vecRefPlanes.size()==1){
		if(TPlaneProperties::X_COR==cor||TPlaneProperties::XY_COR==cor)
			align->AddToXOffset(subjectPlane,x_offset);
		if(TPlaneProperties::Y_COR==cor||TPlaneProperties::XY_COR==cor)
			align->AddToYOffset(subjectPlane,y_offset);
	}
	else{
		if(subjectPlane==4&&nDiaAlignmentStep==0){
			align->AddToXOffset(subjectPlane,x_offset);
			return res;
		}
		if(TPlaneProperties::X_COR==cor||TPlaneProperties::XY_COR==cor){
			align->AddToXOffset(subjectPlane,x_offset);
			align->AddToPhiYOffset(subjectPlane,phiy_offset);
		}
		if(TPlaneProperties::Y_COR==cor||TPlaneProperties::XY_COR==cor){
			align->AddToYOffset(subjectPlane,y_offset);
			align->AddToPhiXOffset(subjectPlane,phix_offset);
		}
	}
	return res;
}

bool TAlignment::SiliconAlignmentStep(bool bPlot){


	cout<<"\n\n\nALIGNMENT STEP:\t"<<nAlignmentStep+1<<" of "<<nAlignSteps<<"\n\n"<<endl;
	bool isAlignmentDone=true;
	if (nAlignmentStep<2){
		for(UInt_t plane=0;plane<TPlaneProperties::getNSiliconPlanes();plane++){
			vector<UInt_t> vecRefPlanes;
			for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
				if(refPlane!=plane)vecRefPlanes.push_back(refPlane);
			TResidual res=CheckDetectorAlignment(TPlaneProperties::XY_COR,plane,vecRefPlanes,false);
			res=CheckDetectorAlignment(TPlaneProperties::XY_COR,plane,vecRefPlanes,false,res);
			align->setXResolution(res.getXSigma(),plane);
			align->setYResolution(res.getYSigma(),plane);
		}
		cout<<endl;
	}
//	for(UInt_t subjectPlane=1;subjectPlane<TPlaneProperties::getNSiliconPlanes()-1;subjectPlane++){
//		vector<UInt_t>vecRefPlanes;
//		cout<<"\nAlign Plane "<<subjectPlane<<" with Plane "<<flush;
//		for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
//			if(refPlane!=subjectPlane){
//				vecRefPlanes.push_back(refPlane);
//				cout<<" "<<refPlane<<flush;
//			}
//		cout<<endl;
//		TResidual res=CheckDetectorAlignment(TPlaneProperties::XY_COR,subjectPlane,vecRefPlanes,false);
//		res = AlignDetector(TPlaneProperties::XY_COR,subjectPlane,vecRefPlanes,false,res);
//		res.getXSigma();
//		align->setXResolution(res.getXSigma(),subjectPlane);
//		align->setYResolution(res.getYSigma(),subjectPlane);
//		Float_t xOff= align->GetLastXOffset(subjectPlane);
//			Float_t yOff= align->GetLastYOffset(subjectPlane);
//			cout<<"addXoffset: "<<xOff<<"/"<<align->GetLastXOffset(subjectPlane)<<"\taddYoffset:"<<yOff<<"/"<<align->GetLastYOffset(subjectPlane)<<endl;
//			if(xOff>0.01||yOff>0.01)
//				isAlignmentDone= false;
//	}


	TResidual res103=CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,false);
	if(verbosity)printf("\n1 with 0 and 3:%d,%1.2f+/-%1.2f\n",res103.isTestResidual(),res103.getXMean(),res103.getXSigma());
	AlignDetector(TPlaneProperties::XY_COR,1,0,3,bPlot,res103);
	Float_t xOff= align->GetLastXOffset(1);
	Float_t yOff= align->GetLastYOffset(1);
	if(xOff>0.01||yOff>0.01)
		isAlignmentDone= false;

	cout<<"\nAlign Plane 2 with Plane 0 and 3"<<endl;
	TResidual res203=CheckDetectorAlignment(TPlaneProperties::XY_COR,2,0,3,false);
	if(verbosity)printf("\n2 with 0 and 3:%d,%1.2f+/-%1.2f\n",res203.isTestResidual(),res203.getXMean(),res203.getXSigma());
	AlignDetector(TPlaneProperties::XY_COR,2,0,3,bPlot,res203);
	xOff= align->GetLastXOffset(2);
	yOff= align->GetLastYOffset(2);
	if(xOff>0.01||yOff>0.01)
		isAlignmentDone= false;


	if(nAlignmentStep<nAlignSteps-2){
	cout<<"\nAlign Plane 3 with Plane 1 and 2"<<endl;
	TResidual res312=CheckDetectorAlignment(TPlaneProperties::XY_COR,3,1,2,false);
	if(verbosity)printf("\n3 with 1 and 2:%d,%1.2f+/-%1.2f\n",res312.isTestResidual(),res312.getXMean(),res312.getXSigma());
	res312=AlignDetector(TPlaneProperties::XY_COR,3,1,2,bPlot,res312);
	xOff= align->GetLastXOffset(3);
	yOff= align->GetLastYOffset(3);
	if(xOff>0.01||yOff>0.01)
		isAlignmentDone= false;
	}

//	align->PrintResults(0);
		if(isAlignmentDone)
			cout<<"AlignmentStep was successfully done"<<endl;
		else
			cout<<"AlignmentStep was NOT the finally alignment step"<<endl;
	cout<<endl;
	return isAlignmentDone;
}

void TAlignment::AlignSiliconPlanes(){
	verbosity=0;
	nAlignmentStep=-1;
	//getChi2Distribution(10000);
	CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,true);
	CheckDetectorAlignment(TPlaneProperties::XY_COR,2,0,3,true);
	CheckDetectorAlignment(TPlaneProperties::XY_COR,3,1,2,true);
	AlignDetector(TPlaneProperties::XY_COR,1,0,0,false);
	AlignDetector(TPlaneProperties::XY_COR,2,0,0,false);
	AlignDetector(TPlaneProperties::XY_COR,3,0,0,false);

//	TResidual res103;//=CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,false);
//	TResidual res203;//=CheckDetectorAlignment(TPlaneProperties::XY_COR,2,0,3,false);
//	TResidual res312;//=CheckDetectorAlignment(TPlaneProperties::XY_COR,3,1,2,false);
//	res103=CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,false);//test
//	AlignDetector(TPlaneProperties::XY_COR,1,0,3,false,res103);//test

	if(verbosity)cout<<endl;

	bool bAlignmentGood=true;
	for(nAlignmentStep=0;nAlignmentStep<nAlignSteps;nAlignmentStep++){//||(nAlignmentStep<10&&!bAlignmentGood);nAlignmentStep++){
		bAlignmentGood=true;
		bAlignmentGood=SiliconAlignmentStep(false);
	}
//	cout<<"\n\n***************************"<<endl;
//	cout<<"DO Eta Corrections Step 1"<<endl;
//	cout<<"***************************"<<endl;
//	DoEtaCorrectionSilicon(1);
//	cout<<"DO new alignment steps"<<endl;
//	bAlignmentGood=true;
//	UInt_t nAlignStep=nAlignmentStep;
//	for(nAlignmentStep;nAlignmentStep<nAlignStep+2;nAligmentStep++){
//		bAlignmentGood=SiliconAlignmentStep(nAlignmentStep==nAlignStep+1);
//	}
//	cout<<"DO Eta Corrections Step 2"<<endl;
//	DoEtaCorrectionSilicon(2);
//	align->PrintResults();
//
////	res103=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,false,res);
////	res203=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,2,0,3,false);
//////	res312=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,3,1,2,false);
////////
////////
////	res103=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,1,0,3,true,res103);
////	res203=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,2,0,3,true,res203);
////	res312=this->CheckDetectorAlignment(TPlaneProperties::XY_COR,3,1,2,true,res312);
////	verbosity=0;
	getFinalSiliconAlignmentResuluts();
}

void TAlignment::AlignDiamondPlane(){
	verbosity=1;
	cout<<"\n\n\n*******************************************************"<<endl;
	cout<<"*******************************************************"<<endl;
	cout<<"***************** Align Diamond ***********************"<<endl;
	cout<<"*******************************************************"<<endl;
	cout<<"*******************************************************\n"<<endl;

	//create ReferencePlane Vector:
	UInt_t diaPlane=4;
	vector<UInt_t>vecRefPlanes;
	for(UInt_t i=0;i<4;i++)if(i!=diaPlane)vecRefPlanes.push_back(i);
	//align->setVerbosity(5);
//	this->verbosity=5;
	nDiaAlignmentStep=-1;
	TResidual resDia=CheckStripDetectorAlignment(TPlaneProperties::X_COR,diaPlane,vecRefPlanes,false,true);
	//myTrack->setVerbosity(1);
	for(nDiaAlignmentStep=0;nDiaAlignmentStep<nDiaAlignSteps;nDiaAlignmentStep++){
			cout<<"\n\n "<<nDiaAlignmentStep<<" of "<<nDiaAlignSteps<<" Steps..."<<endl;
//		if(nDiaAlignmentStep==0)verbosity=4;
//		else verbosity=0;
		AlignStripDetector(TPlaneProperties::X_COR,diaPlane,vecRefPlanes,false,resDia);
		resDia=CheckStripDetectorAlignment(TPlaneProperties::X_COR,diaPlane,vecRefPlanes);
	}
	verbosity=0;
	resDia=CheckStripDetectorAlignment(TPlaneProperties::X_COR,diaPlane,vecRefPlanes,true,true,resDia);
//	verbosity=100;
	CheckStripDetectorAlignmentChi2(TPlaneProperties::X_COR,diaPlane,vecRefPlanes,true,true,settings->getAlignment_chi2());
	align->setDiamondDate();
}

TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2,bool bPlot,TResidual resOld,TCluster::calculationMode_t mode){
	vector<UInt_t>vecRefPlanes;
	vecRefPlanes.push_back(refPlane1);
	if(refPlane1!=refPlane2)
		vecRefPlanes.push_back(refPlane2);
	return getResidual(cor,subjectPlane,vecRefPlanes,bPlot,resOld,mode);
}



TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld,TCluster::calculationMode_t mode)
{
	stringstream  refPlaneString;
	for(UInt_t i=0;i<vecRefPlanes.size();i++)
		if(i==0)
			refPlaneString<<vecRefPlanes.at(i);
		else if(i+1<vecRefPlanes.size())
			refPlaneString<<"_"<<vecRefPlanes.at(i);
		else
			refPlaneString<<"_and_"<<vecRefPlanes.at(i);
	if(verbosity)cout<<"TAlignment::getResidual of Plane "<<subjectPlane<<TPlaneProperties::getCoordinateString(cor)<<" with "<<refPlaneString.str()<<", plotting: "<<bPlot<<"  with "<<alignmentPercentage<<"\t"<<resOld.isTestResidual()<<endl;
	vecXPred.clear();
	vecYPred.clear();
	vecXObs.clear();
	vecYObs.clear();
	vecDeltaX.clear();
	vecDeltaY.clear();
	vecXChi2.clear();
	vecYChi2.clear();

	Float_t deltaX, deltaY;
	Float_t xPositionObserved;
	Float_t yPositionObserved;
	TPositionPrediction* predictedPostion=0;
	Float_t resxtest,resytest;
	for(UInt_t nEvent=0; nEvent<events.size(); nEvent++)
	{
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));
		xPositionObserved  = myTrack->getPosition(TPlaneProperties::X_COR,subjectPlane,mode);
		yPositionObserved  = myTrack->getPosition(TPlaneProperties::Y_COR,subjectPlane,mode);
		predictedPostion  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,false);
		deltaX = xPositionObserved-predictedPostion->getPositionX();//X_OBS-X_Pred
		deltaY = yPositionObserved-predictedPostion->getPositionY();//Y_OBS-Y_Pred
		resxtest= TMath::Abs(deltaX-resOld.getXMean())/resOld.getXSigma();
		resytest= TMath::Abs(deltaY-resOld.getYMean())/resOld.getYSigma();
		if(verbosity>3)cout<<nEvent<<endl;
		//if(verbosity>3)	predictedPostion->Print();
		if(verbosity>3)	cout<<"Measured: "<<myTrack->getXMeasured(subjectPlane,mode)<<"/"<<myTrack->getYMeasured(subjectPlane,mode)<<endl;
		if(verbosity>3)	cout<<"Observed: "<<xPositionObserved<<" / "<<yPositionObserved<<endl;
		if(verbosity>3)	cout<<"Predicted: "<<predictedPostion->getPositionX()<<" / "<<predictedPostion->getPositionY()<<endl;
		if(verbosity>3)	cout<<"Delta:    "<<deltaX<<" / "<<yPositionObserved<<endl;
		if(verbosity>3)	cout<<"ResTest:  "<<resxtest<<" / "<<resytest<<"\n\n"<<endl;

		if(resxtest<res_keep_factor&&resytest<res_keep_factor){
			vecXObs.push_back(xPositionObserved);
			vecYObs.push_back(yPositionObserved);
			vecDeltaX.push_back(deltaX);
			vecDeltaY.push_back(deltaY);
			vecXPred.push_back(predictedPostion->getPositionX());
			vecYPred.push_back(predictedPostion->getPositionY());

			vecXChi2.push_back(predictedPostion->getChi2X());
			vecYChi2.push_back(predictedPostion->getChi2Y());
		}
		if(verbosity>3)	cout<< deltaX<<" "<<deltaY<<endl;
		predictedPostion->Delete();
	}

	if(verbosity>2)cout<<vecDeltaX.size()<<" "<<vecDeltaY.size()<<" "<< vecXPred.size()<<" "<<vecYPred.size()<<endl;
	//first estimate residuals widths
	TResidual res = calculateResidual(cor,&vecXPred,&vecDeltaX,&vecYPred,&vecDeltaY,resOld);
	this->CreatePlots(cor, subjectPlane,refPlaneString.str(),bPlot);

	vecXObs.clear();
	vecYObs.clear();
	vecDeltaX.clear();
	vecDeltaY.clear();
	vecXPred.clear();
	vecYPred.clear();
	return res;

}


TResidual TAlignment::getStripResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld,TCluster::calculationMode_t mode){

	stringstream  refPlaneString;
		for(UInt_t i=0;i<vecRefPlanes.size();i++)
			if(i==0)
				refPlaneString<<vecRefPlanes.at(i);
			else if(i+1<vecRefPlanes.size())
				refPlaneString<<"_"<<vecRefPlanes.at(i);
			else
				refPlaneString<<"_and_"<<vecRefPlanes.at(i);
		if(verbosity)cout<<"TAlignment::getStripResidual of Plane "<<subjectPlane<<TPlaneProperties::getCoordinateString(cor)<<" with "<<refPlaneString.str()<<", plotting: "<<bPlot<<"  with "<<alignmentPercentage<<"\t"<<resOld.isTestResidual()<<endl;
		vecXPred.clear();
		vecYPred.clear();
		vecXObs.clear();
		vecYObs.clear();
		vecDeltaX.clear();
		vecDeltaY.clear();
		vecXChi2.clear();
		vecYChi2.clear();

		Float_t deltaX, deltaY;
		Float_t xPositionObserved;
		Float_t yPositionObserved;
		TPositionPrediction* predictedPostion=0;
		Float_t resxtest,resytest;
		for(UInt_t nEvent=0; nEvent<events.size(); nEvent++)
		{
			TRawEventSaver::showStatusBar(nEvent,events.size());
			myTrack->setEvent(&events.at(nEvent));
			predictedPostion  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,verbosity>1);
			xPositionObserved  = myTrack->getStripXPosition(subjectPlane,predictedPostion->getPositionY(),TCluster::maxValue);
			yPositionObserved  = myTrack->getPosition(TPlaneProperties::Y_COR,subjectPlane,TCluster::maxValue);
			deltaX = xPositionObserved-predictedPostion->getPositionX();//X_OBS-X_Pred
			deltaY = yPositionObserved-predictedPostion->getPositionY();//Y_OBS-Y_Pred
			resxtest= TMath::Abs(deltaX-resOld.getXMean())/resOld.getXSigma();
			resytest= TMath::Abs(deltaY-resOld.getYMean())/resOld.getYSigma();
			if(verbosity>3)cout<<"Event no.: "<<nEvent<<endl;
			//if(verbosity>3)	predictedPostion->Print();
			if(verbosity>3) events.at(nEvent).getPlane(subjectPlane).Print();
			if(verbosity>3)	cout<<"Measured: "<<myTrack->getXMeasured(subjectPlane)<<"/"<<myTrack->getYMeasured(subjectPlane)<<endl;
			if(verbosity>3)	cout<<"Observed: "<<xPositionObserved<<" / "<<yPositionObserved<<endl;
			if(verbosity>3)	cout<<"Predicted: "<<predictedPostion->getPositionX()<<" / "<<predictedPostion->getPositionY()<<endl;
			if(verbosity>3)	cout<<"Delta:    "<<deltaX<<" / "<<yPositionObserved<<endl;
			if(verbosity>3)	cout<<"ResTest:  "<<resxtest<<" / "<<resytest<<"\n\n"<<endl;

			if(resxtest<res_keep_factor&&resytest<res_keep_factor){
				vecXObs.push_back(xPositionObserved);
				vecYObs.push_back(yPositionObserved);
				vecDeltaX.push_back(deltaX);
				vecDeltaY.push_back(deltaY);
				vecXPred.push_back(predictedPostion->getPositionX());
				vecYPred.push_back(predictedPostion->getPositionY());

				vecXChi2.push_back(predictedPostion->getChi2X());
				vecYChi2.push_back(predictedPostion->getChi2Y());
			}
			if(verbosity>3)	cout<< deltaX<<" "<<deltaY<<endl;
			predictedPostion->Delete();
		}

		if(verbosity>2)cout<<vecDeltaX.size()<<" "<<vecDeltaY.size()<<" "<< vecXPred.size()<<" "<<vecYPred.size()<<endl;

		//first estimate residuals widths
		TResidual res = calculateResidual(cor,&vecXPred,&vecDeltaX,&vecYPred,&vecDeltaY,resOld);
		this->CreatePlots(cor, subjectPlane,refPlaneString.str(),bPlot,bAlign);

		vecXObs.clear();
		vecYObs.clear();
		vecDeltaX.clear();
		vecDeltaY.clear();
		vecXPred.clear();
		vecYPred.clear();
		return res;

}


TResidual TAlignment::getStripResidualChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bAlign,bool bPlot,Float_t maxChi2,TCluster::calculationMode_t mode){

	TResidual resOld;
	stringstream  refPlaneString;
		for(UInt_t i=0;i<vecRefPlanes.size();i++)
			if(i==0)
				refPlaneString<<vecRefPlanes.at(i);
			else if(i+1<vecRefPlanes.size())
				refPlaneString<<"_"<<vecRefPlanes.at(i);
			else
				refPlaneString<<"_and_"<<vecRefPlanes.at(i);
		if(verbosity)cout<<"TAlignment::getStripResidual of Plane "<<subjectPlane<<TPlaneProperties::getCoordinateString(cor)<<" with "<<refPlaneString.str()<<", plotting: "<<bPlot<<"  with "<<alignmentPercentage<<"\t"<<resOld.isTestResidual()<<endl;
		vecXPred.clear();
		vecYPred.clear();
		vecXObs.clear();
		vecYObs.clear();
		vecDeltaX.clear();
		vecDeltaY.clear();
		vecXChi2.clear();
		vecYChi2.clear();
		vecMeasuredX.clear();
		vecMeasuredY.clear();

		Float_t deltaX, deltaY;
		Float_t xPositionObserved;
		Float_t yPositionObserved;
		Float_t xMeasured,yMeasured;
		TPositionPrediction* predictedPosition=0;
		Float_t resxtest,resytest;
		Float_t chi2x,chi2y;
		for(UInt_t nEvent=0; nEvent<events.size(); nEvent++)
		{
			TRawEventSaver::showStatusBar(nEvent,events.size());
			myTrack->setEvent(&events.at(nEvent));
			predictedPosition  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,nEvent<1);
			chi2x=predictedPosition->getChi2X();
			chi2y=predictedPosition->getChi2Y();
			xPositionObserved  = myTrack->getStripXPosition(subjectPlane,predictedPosition->getPositionY(),TCluster::maxValue);
			yPositionObserved  = myTrack->getPosition(TPlaneProperties::Y_COR,subjectPlane,TCluster::maxValue);
			deltaX = xPositionObserved-predictedPosition->getPositionX();//X_OBS-X_Pred
			deltaY = yPositionObserved-predictedPosition->getPositionY();//Y_OBS-Y_Pred
			xMeasured = myTrack->getMeasured(TPlaneProperties::X_COR,subjectPlane,TCluster::maxValue);
			yMeasured = myTrack->getMeasured(TPlaneProperties::Y_COR,subjectPlane,TCluster::maxValue);
			resxtest= TMath::Abs(deltaX-resOld.getXMean())/resOld.getXSigma();
			resytest= TMath::Abs(deltaY-resOld.getYMean())/resOld.getYSigma();
			if(verbosity>3)cout<<"Event no.: "<<nEvent<<endl;
			//if(verbosity>3)	predictedPostion->Print();
			if(verbosity>3) events.at(nEvent).getPlane(subjectPlane).Print();
			if(verbosity>3)	cout<<"Measured: "<<myTrack->getXMeasured(subjectPlane)<<"/"<<myTrack->getYMeasured(subjectPlane)<<endl;
			if(verbosity>3)	cout<<"Observed: "<<xPositionObserved<<" / "<<yPositionObserved<<endl;
			if(verbosity>3)	cout<<"Predicted: "<<predictedPosition->getPositionX()<<" / "<<predictedPosition->getPositionY()<<endl;
			if(verbosity>3)	cout<<"Delta:    "<<deltaX<<" / "<<yPositionObserved<<endl;
			if(verbosity>3) cout<<"Chi2:     "<<chi2x<<" / "<<chi2y<<endl;
			if(verbosity>3)	cout<<"ResTest:  "<<resxtest<<" / "<<resytest<<endl;


			if(chi2x<settings->getAlignment_chi2()&&chi2y<settings->getAlignment_chi2()){
				if(verbosity>3)cout<<"Take event for Residual"<<endl;
				vecXObs.push_back(xPositionObserved);
				vecYObs.push_back(yPositionObserved);
				vecDeltaX.push_back(deltaX);
				vecDeltaY.push_back(deltaY);
				vecXPred.push_back(predictedPosition->getPositionX());
				vecYPred.push_back(predictedPosition->getPositionY());
				vecMeasuredX.push_back(xMeasured);
				vecMeasuredY.push_back(yMeasured);
				vecXChi2.push_back(predictedPosition->getChi2X());
				vecYChi2.push_back(predictedPosition->getChi2Y());
			}
			else
				if(verbosity>3)cout<<"through event away..."<<endl;
			if(verbosity>3)cout<<"\n\n"<<endl;
//			if(verbosity>3)	cout<< deltaX<<" "<<deltaY<<endl;
			predictedPosition->Delete();
		}
		if(bAlign){
			align->setNDiamondAlignmentEvents((UInt_t)vecDeltaX.size());
			align->setDiaChi2(settings->getAlignment_chi2());
		}

		if(verbosity)cout<<vecDeltaX.size()<<" "<<vecDeltaY.size()<<" "<< vecXPred.size()<<" "<<vecYPred.size()<<" "<<vecXObs.size()<<" "<<vecYObs.size()<<endl;
		if(verbosity)cout<<"used "<<vecDeltaX.size()<<" Events of "<<events.size()<<" with a Chi2 < "<<settings->getAlignment_chi2()<<endl;
		//first estimate residuals widths
		TResidual res = calculateResidual(cor,&vecXPred,&vecDeltaX,&vecYPred,&vecDeltaY);
		this->CreatePlots(cor, subjectPlane,refPlaneString.str(),bPlot,bAlign);

		vecXObs.clear();
		vecYObs.clear();
		vecDeltaX.clear();
		vecDeltaY.clear();
		vecXPred.clear();
		vecYPred.clear();
		return res;
}
/**
 * calculateResidual if there is no residuals calculatet to cut on the res_keep_factor;
 * this funcition opens the other calculateResidual function but uses a TResidual res
 * which is using all items...
 */
TResidual TAlignment::calculateResidual(TPlaneProperties::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>* yPred,vector<Float_t>* deltaY){
	TResidual res;
	res.SetTestResidual();
	if(verbosity>2)cout<<"\t"<<deltaX->size()<<" "<<deltaY->size()<<" "<< xPred->size()<<" "<<yPred->size()<<endl;
	return calculateResidual(cor,xPred,deltaX,yPred,deltaY,res);
}

TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot,TResidual resOld)
{
//	if(verbosity)cout<<"\n\nTAlignment::checkDetectorAlignment\n\t check "<<TPlaneProperties::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with Plane "<<refPlane1<<" and "<<refPlane2<<endl;
//	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
//		return TResidual(true);
//	}
//	int verb=verbosity;
//	verbosity=0;
//	TResidual res = getResidual(cor,subjectPlane,refPlane1,refPlane2,false,resOld);
//	if(verbosity)cout<<endl;
//	res.SetTestResidual(false);
//	res = getResidual(cor,subjectPlane,refPlane1,refPlane2,bPlot,res);
//	verbosity=verb;
//	if(verbosity)res.Print();
//	return res;
	vector<UInt_t>vecRefPlanes;
	vecRefPlanes.push_back(refPlane1);
	vecRefPlanes.push_back(refPlane2);
	return CheckDetectorAlignment(cor,subjectPlane,vecRefPlanes,bPlot,resOld);
}


TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld)
{
	if(verbosity)cout<<"\n\nTAlignment::checkDetectorAlignment\n\t check "<<TPlaneProperties::getCoordinateString(cor)<<" coordinate of Plane "<<subjectPlane<<" with "<<vecRefPlanes.size()<<" Planes"<<endl;
//	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
//		return TResidual(true);
//	}
	int verb=verbosity;
	verbosity=0;
	TResidual res = getResidual(cor,subjectPlane,vecRefPlanes,false,resOld);
	if(verbosity)cout<<endl;
	res.SetTestResidual(false);
	res = getResidual(cor,subjectPlane,vecRefPlanes,bPlot,res);
	verbosity=verb;
	if(verbosity)res.Print();
	return res;
}

TResidual TAlignment::CheckStripDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld){
	int verb=verbosity;
	verbosity=0;
	TResidual res = getStripResidual(cor,subjectPlane,vecRefPlanes,false,false,resOld);
	if(verbosity)cout<<endl;
	res.SetTestResidual(false);
	res = getStripResidual(cor,subjectPlane,vecRefPlanes,bAlign,bPlot,res);
	verbosity=verb;
	if(verbosity)res.Print();
	return res;
}

TResidual TAlignment::CheckStripDetectorAlignmentChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2){
	int verb=verbosity;
	verbosity=1;

//	getStripResidual(cor,subjectPlane,vecRefPlanes,false,false,maxChi2);
	if(verbosity)cout<<endl;
	TResidual res = getStripResidualChi2(cor,subjectPlane,vecRefPlanes,bAlign,bPlot,maxChi2);
	if(verbosity)res.Print();
	verbosity=verb;
	return res;
}
/**
 *
 */
void TAlignment::saveAlignment()
{
	stringstream fileName;
	fileName<<"alignment."<<settings->getRunNumber()<<".root";
	TFile *alignmentFile=new TFile(fileName.str().c_str(),"RECREATE");//todo anpassen runnumber
	cout<<"TAlignment:saveAlignment(): path: \""<<sys->pwd()<<"\", file Name:\""<<fileName.str()<<"\""<<endl;
	alignmentFile->cd();
	align->SetName("alignment");
	stringstream title;
	title<<"alignment of Run "<<settings->getRunNumber();
	align->SetTitle(title.str().c_str());
	align->Write();
	alignmentFile->Write();
	alignmentFile->Close();
}


/**
 *
 */
void TAlignment::getChi2Distribution(Float_t maxChi2)
{
	vecXChi2.clear();
	vecYChi2.clear();
//	if(verbosity)
		cout<<"TAlignment::getChi2Distribution"<<endl;
//	Float_t xPositionObserved,yPositionObserved,deltaX,deltaY,resxtest,resytest;
	TPositionPrediction* predictedPosition=0;
	vector<UInt_t> vecRefPlanes;

	for(UInt_t i=0;i<4;i++)
		vecRefPlanes.push_back(i);
//	UInt_t oldVerbosity=myTrack->getVerbosity();
//	myTrack->setVerbosity(4);
	vector<Float_t> vecSumDeltaX;
	vector<Float_t> vecSumDeltaY;
	for(UInt_t nEvent=0; nEvent<events.size(); nEvent++) {
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));
		Float_t sumDeltaX=0;
		Float_t sumDeltaY=0;
		Float_t chi2X=0;
		Float_t chi2Y=0;
		UInt_t subjectPlane=0;
		predictedPosition  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,false);
		chi2X=predictedPosition->getChi2X();
		chi2Y= predictedPosition->getChi2Y();
		if(predictedPosition->getChi2X()<maxChi2&&predictedPosition->getChi2Y()<maxChi2){
			for(subjectPlane=0;subjectPlane<4;subjectPlane++){
				if(subjectPlane!=0){
					predictedPosition->Delete();
					predictedPosition  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,false);
				}
				Float_t deltaX=myTrack->getXPosition(subjectPlane);;
				Float_t deltaY=myTrack->getYPosition(subjectPlane);;

				deltaX-=predictedPosition->getPositionX();
				deltaY-=predictedPosition->getPositionY();
				sumDeltaX+=TMath::Abs(deltaX);
				sumDeltaY+=TMath::Abs(deltaY);
			}//for loop over subjectPlane
			vecXChi2.push_back(chi2X);
			vecYChi2.push_back(chi2Y);
			vecSumDeltaX.push_back(sumDeltaX);
			vecSumDeltaY.push_back(sumDeltaY);
		}//end if chi2x<maxCi && chi2y<maxChi
		predictedPosition->Delete();
	}//end for loop over nEvent

//	myTrack->setVerbosity(oldVerbosity);
	stringstream histName;
	//Chi2X Distribution
	histName.str("");
	if(nAlignmentStep==-1)histName<<"hPreAlignment";
	else if(nAlignmentStep>=nAlignSteps-1)histName<<"hPostAlignment";
	else histName<<"h_"<<nAlignmentStep<<"_Step";
	histName<<"_Chi2X_Distribution";
	TH1F histoChi2X = histSaver->CreateDistributionHisto(histName.str(),vecXChi2,4096,HistogrammSaver::positiveSigma);
	histoChi2X.GetXaxis()->SetTitle("Chi^2/NDF of X plane");
	histoChi2X.GetYaxis()->SetTitle("number of entries");

	histSaver->SaveHistogram(&histoChi2X);

	//Chi2Y Distribution
	histName.str("");
	if(nAlignmentStep==-1)histName<<"hPreAlignment";
	else if(nAlignmentStep>=nAlignSteps-1)histName<<"hPostAlignment";
	else histName<<"h_"<<nAlignmentStep<<"_Step";
	histName<<"_Chi2Y_Distribution";
	TH1F histoChi2Y=histSaver->CreateDistributionHisto(histName.str(),vecYChi2,4096,HistogrammSaver::positiveSigma);
	histoChi2Y.GetXaxis()->SetTitle("Chi^2/NDF of Y plane");
	histoChi2Y.GetYaxis()->SetTitle("number of entries");
	histSaver->SaveHistogram(&histoChi2Y);

	histName.str("");
	if(nAlignmentStep==-1)histName<<"hPreAlignment";
	else if(nAlignmentStep>=nAlignSteps-1)histName<<"hPostAlignment";
	else histName<<"h_"<<nAlignmentStep<<"_Step";
	histName<<"_Chi2X_vs_SumDeltaX";
	cout<<"CREATE: "<<histName.str()<<endl;
	TH2F histo= histSaver->CreateScatterHisto(histName.str(),vecSumDeltaX,vecXChi2);
	histo.GetYaxis()->SetTitle("Sum of Delta X");
	histo.GetXaxis()->SetTitle("Chi2 in X");
	histSaver->SaveHistogram(&histo);


	histName.str("");
	if(nAlignmentStep==-1)histName<<"hPreAlignment";
	else if(nAlignmentStep>=nAlignSteps-1)histName<<"hPostAlignment";
	else histName<<"h_"<<nAlignmentStep<<"_Step";
	histName<<"_Chi2Y_vs_SumDeltaY";
	cout<<"CREATE PLOT: \""<<histName.str()<<"\""<<endl;
	TH2F histo1= histSaver->CreateScatterHisto(histName.str(),vecSumDeltaY,vecYChi2);
	histName<<"_graph";
	TGraph graph1 = histSaver->CreateDipendencyGraph(histName.str(),vecSumDeltaY,vecYChi2);
	graph1.Draw("AP");
	graph1.GetYaxis()->SetTitle("Sum of Delta Y");
	graph1.GetXaxis()->SetTitle("Chi2 in Y");

	histo1.GetYaxis()->SetTitle("Sum of Delta Y");
	histo1.GetXaxis()->SetTitle("Chi2 in Y");
	histSaver->SaveHistogram(&histo1);
	histSaver->SaveGraph(&graph1,histName.str(),"AP");
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
 * @TODO POINTER UEBERGEBEN
 */
TResidual TAlignment::calculateResidual(TPlaneProperties::enumCoordinate cor,vector<Float_t>*xPred,vector<Float_t>* deltaX,vector<Float_t>*yPred,vector<Float_t>* deltaY,TResidual resOld)
{
	if(verbosity)cout<<"\n\tTAlignment::calculateResidual with "<<resOld.getXMean()<<"+/-"<<resOld.getXSigma()<<"\t"<<resOld.getYMean()<<"+/-"<<resOld.getYSigma();
	if(verbosity)if(resOld.isTestResidual())cout<<"\tno residual Correction, inputRes is a TestResidual";
	if(verbosity)cout<<endl;
	TResidual residual(false);
	Float_t resxtest;
	Float_t resytest;
	if(verbosity>2)cout<<"\tcalculate Residual "<<res_keep_factor<<endl;
	if(verbosity>2)cout<<"\t"<<deltaX->size()<<" "<<deltaY->size()<<" "<< xPred->size()<<" "<<yPred->size()<<endl;
	for(UInt_t i=0;i<deltaX->size();i++){
		resxtest= TMath::Abs(deltaX->at(i)-resOld.getXMean())/resOld.getXSigma();
		resytest= TMath::Abs(deltaY->at(i)-resOld.getYMean())/resOld.getYSigma();

		//only add if restest is smaller than res_keep_factor
		if((cor==TPlaneProperties::X_COR)&&resxtest<res_keep_factor){
			residual.addDataPoint(deltaX->at(i),xPred->at(i),deltaY->at(i),yPred->at(i));
		}//end if
		else if((cor==TPlaneProperties::X_COR)&&resytest<res_keep_factor){
			residual.addDataPoint(deltaX->at(i),xPred->at(i),deltaY->at(i),yPred->at(i));
		}//end else if
		else if((cor==TPlaneProperties::XY_COR)&&resxtest<res_keep_factor&&resytest<res_keep_factor){
			residual.addDataPoint(deltaX->at(i),xPred->at(i),deltaY->at(i),yPred->at(i));
		}//end else if

	}//end for loop
	if(verbosity)cout<<"\n";
	if(!resOld.isTestResidual()&&verbosity)printf("\tresidual with x:%1.2f+/-%1.2f and y:%1.2f+/-%1.2f\n",resOld.getXMean(),resOld.getXSigma(),resOld.getYMean(),resOld.getYSigma());
	if(verbosity>0)	cout<<"\tused "<<residual.getUsedTracks()<<" Tracks"<<endl;
	if(verbosity>0)	cout<<"\tX: "<<std::setprecision(4)<<residual.getXMean()<<"+/-"<<residual.getXSigma()<<endl;
	if(verbosity>0)	cout<<"\tY: "<<residual.getYMean()<<"+/-"<<residual.getYSigma()<<"\n"<<endl;
	//set values
	residual.SetTestResidual(false);
	return residual;
}


/**
 * @brief: Print Informations for all Events smaller than maxEvent, starting with startEvent
 * @param	maxEvent	maximum Event No to print event
 * @param: 	startEvent	first Event where to start with printing Event Information
 *
 */
void TAlignment::getFinalSiliconAlignmentResuluts()
{
	Float_t maxChi2=settings->getAlignment_chi2();
	cout<<"get Final Silicon Resolution with a  maximum chi2 of "<<maxChi2<<endl;
	//set Resolutions
//	setDetectorResolution(maxChi2);
	setSiliconDetectorResolution(maxChi2);

//	vector<UInt_t>vecRefPlanes;
//	for(UInt_t plane=0;plane<4;plane++){
//		vecRefPlanes.clear();
//		for(UInt_t i=0;i<4;i++)
//			if(i!=plane)vecRefPlanes.push_back(i);
//			TResidual res0=calculateResidualWithChi2(TPlaneProperties::XY_COR,plane,vecRefPlanes,maxChi2,false,true);
//	}

	getChi2Distribution(15);
}

void TAlignment::setSiliconDetectorResolution(Float_t maxChi2)
{
	//get  something like a aprroximate Sigma with calculating the residual
	for(UInt_t plane=0;plane<4;plane++){
		vector<UInt_t>vecRefPlanes;
		for(UInt_t i=0;i<4;i++)
			if(i!=plane)vecRefPlanes.push_back(i);
		TResidual res=CheckDetectorAlignment(TPlaneProperties::XY_COR,plane,vecRefPlanes,false);
		res=CheckDetectorAlignment(TPlaneProperties::XY_COR,plane,vecRefPlanes,false,res);
		align->setXResolution(res.getXSigma(),plane);
		align->setYResolution(res.getYSigma(),plane);
	}

	for(UInt_t plane=0;plane<4;plane++){
		vector<UInt_t>vecRefPlanes;
		for(UInt_t i=0;i<4;i++)
			if(i!=plane)vecRefPlanes.push_back(i);
		TResidual res=calculateResidualWithChi2(TPlaneProperties::XY_COR,plane,vecRefPlanes,maxChi2,true,true);
	}
}



void TAlignment::CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,string refPlaneString, bool bPlot, bool bUpdateAlignment,bool bChi2)
{
	if(bPlot)
		if(verbosity)cout<<vecDeltaX.size()<<" "<<vecDeltaY.size()<<" "<< vecXPred.size()<<" "<<vecYPred.size()<<" "<<vecXObs.size()<<" "<<vecYObs.size()<<endl;
	if(!bPlot&&!bUpdateAlignment)return;
	stringstream preName;
	if(subjectPlane==4){
		preName<<"hDiamond_";
		if(nDiaAlignmentStep==-1)preName<<"PreAlignment";
				else if(nDiaAlignmentStep==nDiaAlignSteps)preName<<"PostAlignment";
				else preName<<nDiaAlignmentStep<<"Step";
	}
	else{
		preName<<"hSilicon_";
		if(nAlignmentStep==-1)preName<<"PreAlignment";
		else if(nAlignmentStep==nAlignSteps)preName<<"PostAlignment";
		else preName<<nAlignmentStep<<"Step";
	}
	stringstream histName;
	cout<<"\nCreatePlots with "<<preName.str()<<flush;
	if(bUpdateAlignment)cout<<"\twith Alignment Resolution Update\n"<<endl;
	else cout<<endl;

	if(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR){//DistributionPlot DeltaX
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_DistributionPlot_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH1F* histo=(TH1F*)histSaver->CreateDistributionHisto(histName.str(),vecDeltaX,512,HistogrammSaver::threeSigma).Clone();
		TF1* fitGausX=new TF1("fitGaus","gaus",-1,1);
		if(bUpdateAlignment){
			cout<<"Alignment for plane"<<subjectPlane<<endl;
			histo->Draw("goff");
			histo->Fit(fitGausX,"Q","",0.5,0.5);
			Float_t xRes=fitGausX->GetParameter(2);
			cout<<"set Resolution via Gaus-Fit: "<<xRes<<" with "<<vecDeltaX.size()<<" Events"<<endl;
			align->setXResolution(xRes,subjectPlane);
			histo->GetXaxis()->SetRangeUser(-5*xRes,+5*xRes);
		}
		histo->GetXaxis()->SetTitle("Delta X in Channels");
		histo->GetYaxis()->SetTitle("Number of entries");
		if(bPlot)histSaver->SaveHistogram(histo);
		delete fitGausX;
		delete histo;
	}

	if(bPlot&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)){//DistributionPlot DeltaX vs Ypred
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_YPred_vs_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecYPred,vecDeltaX,256);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("Y Predicted");
		histo.GetYaxis()->SetTitle("Delta X");
		histSaver->SaveHistogram(&histo);
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecYPred);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("predicted Y position in ChannelNo.");
		graph.GetYaxis()->SetTitle("delta X in Channel No.");
		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

	if(bPlot&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)&&subjectPlane==TPlaneProperties::getDiamondPlane()){//DistributionPlot DeltaX vs Xpred
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_XPred_vs_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecXPred,vecDeltaX,256);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("X Predicted");
		histo.GetYaxis()->SetTitle("Delta X");
		histSaver->SaveHistogram(&histo);
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecXPred);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("predicted X position in ChannelNo.");
		graph.GetYaxis()->SetTitle("delta X in Channel No.");
		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

	if(bPlot&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)&&subjectPlane==TPlaneProperties::getDiamondPlane()){//DistributionPlot DeltaX vs Xpred
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_XMeasured_vs_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecMeasuredX,vecDeltaX,256);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("X Measured");
		histo.GetYaxis()->SetTitle("Delta X");
		histSaver->SaveHistogram(&histo);
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecMeasuredX);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("measured X position in ChannelNo.");
		graph.GetYaxis()->SetTitle("delta X in Channel No.");
		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

	if(subjectPlane<4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)){//DistributionPlot DeltaY
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_DistributionPlot_DeltaY";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH1F* histo = (TH1F*)histSaver->CreateDistributionHisto(histName.str(),vecDeltaY,512,HistogrammSaver::threeSigma).Clone();
		TF1* fitGausY=new TF1("fitGaus","gaus",-1,1);
		if(bUpdateAlignment){
			histo->Draw("goff");
			histo->Fit(fitGausY,"Q","",0.5,0.5);
			Float_t yRes=fitGausY->GetParameter(2);
			align->setYResolution(yRes,subjectPlane);
			histo->GetXaxis()->SetRangeUser(-5*yRes,+5*yRes);
		}

		histo->GetXaxis()->SetTitle("Delta Y in Channels");
		histo->GetYaxis()->SetTitle("Number of entries");
		if(bPlot)histSaver->SaveHistogram(histo);
		delete fitGausY;
		delete histo;
	}

	if(bPlot&&subjectPlane<4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)){//ScatterPlot DeltaY vs Xpred
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_XPred_vs_DeltaY";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecXPred,vecDeltaY,256);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("X Predicted");
		histo.GetYaxis()->SetTitle("Delta Y");
		histSaver->SaveHistogram((TH2F*)histo.Clone());
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaY,vecXPred);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Predicted X position in channel no");
		graph.GetYaxis()->SetTitle("Delta Y in channel no");

		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

	if(bPlot&&subjectPlane<4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)){//ScatterHisto XObs vs YObs
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_YPred_vs_YObs";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo =histSaver->CreateScatterHisto(histName.str(),vecXObs,vecYObs);
		histo.GetXaxis()->SetTitle("XObs");
		histo.GetYaxis()->SetTitle("YObs");
		histSaver->SaveHistogram((TH2F*)histo.Clone());//,histName.str());
	}

	if(bPlot&&nAlignmentStep>-1&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)){//ScatterHisto DeltaX vs Chi2X
		histName.str("");
		histName.clear();
		histName<<preName.str();
		histName<<"_ScatterPlot_DeltaX_vs_Chi2X";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateScatterHisto(histName.str(),vecDeltaX,vecXChi2,256);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("Delta X");
		histo.GetYaxis()->SetTitle("Chi2 X");
		histSaver->SaveHistogram((TH2F*)histo.Clone());
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecXChi2);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Chi^2 per NDF");
		graph.GetYaxis()->SetTitle("Delta X in channel no");
		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}
	if(bPlot&&nAlignmentStep>-1&&subjectPlane<4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::Y_COR)){//ScatterHisto DeltaY vs Chi2Y
		histName.str("");
		histName<<preName.str();
		histName<<"_ScatterPlot_DeltaX_vs_Chi2X";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateScatterHisto(histName.str(),vecDeltaY,vecYChi2,256);
		histo.Draw("goff");
		histo.GetYaxis()->SetTitle("Delta Y");
		histo.GetXaxis()->SetTitle("Chi2 Y");
		histSaver->SaveHistogram((TH2F*)histo.Clone());
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaY,vecYChi2);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Chi^2 per NDF");
		graph.GetYaxis()->SetTitle("Delta Y in channel no");
		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}
	if(bPlot&&subjectPlane==4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)){//predX vs deltaX
		histName.str("");
		histName<<preName.str();
		histName<<"_ScatterPlot_XPred_vs_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecXPred,vecDeltaX,512);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("X Predicted");
		histo.GetYaxis()->SetTitle("Delta X");
		histSaver->SaveHistogram((TH2F*)histo.Clone());
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecXPred);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Predicted X position in channel no");
		graph.GetYaxis()->SetTitle("Delta X in channel no");

		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

	//TODO!!!!!!
	if(bPlot&&subjectPlane==4&&(cor==TPlaneProperties::XY_COR||cor==TPlaneProperties::X_COR)){//DeltaX vs ClusterSize
		histName.str("");
		histName<<preName.str();
		histName<<"_ScatterPlot_ClusterSize_	vs_DeltaX";
		histName<<"_-_Plane_"<<subjectPlane<<"_with_"<<refPlaneString;;//<<"_with"<<refPlane1<<"_and_"<<refPlane2;
		TH2F histo=histSaver->CreateDipendencyHisto(histName.str(),vecClusterSize,vecDeltaX,512);
		histo.Draw("goff");
		histo.GetXaxis()->SetTitle("Cluster Size");
		histo.GetYaxis()->SetTitle("Delta X");
		histSaver->SaveHistogram((TH2F*)histo.Clone());
		histName<<"_graph";
		TGraph graph =histSaver->CreateDipendencyGraph(histName.str(),vecDeltaX,vecClusterSize);
		graph.Draw("APL");
		graph.GetXaxis()->SetTitle("Cluster Size");
		graph.GetYaxis()->SetTitle("Delta X in channel no");

		histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
	}

}


void TAlignment::DoEtaCorrectionSilicon(UInt_t correctionStep)
{
	cout<<"****************************************"<<endl;
	cout<<"Do Eta correction for all silicon planes "<<correctionStep<< endl;
	cout<<"****************************************"<<endl;
	vector<TH1F*> histoStripDistribution,histoStripDistributionFlattned;
	vector<TH1F *>vecHEta;
	stringstream fileName;
	fileName<<	"etaIntegral_Step"<<correctionStep<<"."<<settings->getRunNumber()<<".root";
	TFile* correctedEtaFile = new TFile(fileName.str().c_str(),"RECREATE");
	cout<<"create Histos..."<<endl;
	for(UInt_t det=0; det<TPlaneProperties::getNSiliconDetectors();det++){
		stringstream histoTitle;
		histoTitle<<"hPredictedStripPosition"<<"_step"<<correctionStep<<"_"<<TADCEventReader::getStringForDetector(det);
		histoStripDistribution.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,-0.501,0.501));
		histoTitle<<"_flattend";
		histoStripDistributionFlattned.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,-0.501,0.501));
		histoTitle.str("");
		histoTitle.clear();
		histoTitle<<"hCorrectedEtaDistribution"<<"_step"<<correctionStep<<"_"<<TADCEventReader::getStringForDetector(det);
		vecHEta.push_back(new TH1F(histoTitle.str().c_str(),histoTitle.str().c_str(),128,0,1));
	}

	cout<<"fill first strip hit histo"<<events.size()<<endl;

	for( nEvent=0;nEvent<this->events.size();nEvent++){
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));

		for(UInt_t subjectPlane =0; subjectPlane<TPlaneProperties::getNSiliconPlanes();subjectPlane++){
//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//				continue;
			vector<UInt_t>refPlanes;
			for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
				if(subjectPlane!=refPlane)refPlanes.push_back(refPlane);
			TPositionPrediction* pred = myTrack->predictPosition(subjectPlane,refPlanes,TCluster::corEta,false);
			Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane*2,pred->getPositionX(),pred->getPositionY());
			Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane*2+1,pred->getPositionX(),pred->getPositionY());
			UInt_t stripMiddleX=(UInt_t) (predictedStripPositionX+0.5);
			Float_t deltaX = predictedStripPositionX-stripMiddleX;
			UInt_t stripMiddleY=(UInt_t) (predictedStripPositionY+0.5);
			Float_t deltaY = predictedStripPositionY-stripMiddleY;
//			cout<<nEvent<<": "<<subjectPlane<<"Fill "<<deltaX<<" "<<deltaY<<endl;
			histoStripDistribution.at(subjectPlane*2)->Fill(deltaX);
			histoStripDistribution.at(subjectPlane*2+1)->Fill(deltaY);
		}
	}
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
	cout<<"create flattened strip hit histo "<<events.size()<<endl;
	for( nEvent=0;nEvent<this->events.size();nEvent++){
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));

		for(UInt_t subjectPlane=0; subjectPlane<TPlaneProperties::getNSiliconPlanes();subjectPlane++){
//			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
//				continue;
			vector<UInt_t>refPlanes;
			for(UInt_t refPlane=0;refPlane<TPlaneProperties::getNSiliconPlanes();refPlane++)
				if(subjectPlane!=refPlane)refPlanes.push_back(refPlane);
			TPositionPrediction* pred = myTrack->predictPosition(subjectPlane,refPlanes,TCluster::corEta,false);

			Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane*2,pred->getPositionX(),pred->getPositionY());
			Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane*2+1,pred->getPositionX(),pred->getPositionY());
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
			histName<<"hEtaIntegral"<<"_step"<<correctionStep<<"_"<<TADCEventReader::getStringForDetector(det);;
			UInt_t nBins = vecHEta.at(det)->GetNbinsX();
			TH1F *histo=new TH1F(histName.str().c_str(),histName.str().c_str(),nBins,0,1);
			Int_t entries = vecHEta.at(det)->GetEntries();
			entries -= vecHEta.at(det)->GetBinContent(0);
			entries -=  vecHEta.at(det)->GetBinContent(nBins+1);
			Int_t sum =0;
			for(UInt_t bin=1;bin<nBins+1;bin++){
				Int_t binContent = vecHEta.at(det)->GetBinContent(bin);
				sum +=binContent;
				Float_t pos =  vecHEta.at(det)->GetBinCenter(bin);
				histo->Fill(pos, (Float_t)sum/(Float_t)entries);
			}
			cout<<"save "<<vecHEta.at(det)->GetTitle()<<" "<<vecHEta.at(det)->GetEntries()<<endl;
			histSaver->SaveHistogram(vecHEta.at(det));
			correctedEtaFile->Add(vecHEta.at(det));
			cout<<"save "<<histo->GetTitle()<<" "<<histo->GetEntries()<<endl;
			histSaver->SaveHistogram(histo);
			correctedEtaFile->Add(histo->Clone());
	}
	correctedEtaFile->Write();
	correctedEtaFile->Close();
}

TResidual TAlignment::calculateResidualWithChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,Float_t maxChi2, bool bAlign,bool bPlot)
{
	if(verbosity)cout<<"\n\nTAlignment::calculateResidualWithChi2"<<endl;
	vecXObs.clear();
	vecYObs.clear();
	vecDeltaX.clear();
	vecDeltaY.clear();
	vecXPred.clear();
	vecYPred.clear();
	vecXChi2.clear();
	vecYChi2.clear();
	stringstream  refPlaneString;
	for(UInt_t i=0;i<vecRefPlanes.size();i++){
		if(i==0)
			refPlaneString<<vecRefPlanes.at(i);
		else if(i+1<vecRefPlanes.size())
			refPlaneString<<"_"<<vecRefPlanes.at(i);
		else
			refPlaneString<<"_and_"<<vecRefPlanes.at(i);
	}
	refPlaneString<<"with_Chi2_cut_on_"<<maxChi2;

	//	Float_t xPositionObserved,yPositionObserved,deltaX,deltaY,resxtest,resytest;
	TPositionPrediction* predictedPosition=0;

	//	UInt_t oldVerbosity=myTrack->getVerbosity();
	//	myTrack->setVerbosity(4);
	Float_t chi2x,chi2y,xPositionObserved,yPositionObserved,deltaX,deltaY;;//_
	for(UInt_t nEvent=0; nEvent<events.size(); nEvent++)
	{
		TRawEventSaver::showStatusBar(nEvent,events.size());
		myTrack->setEvent(&events.at(nEvent));
		predictedPosition  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,false);
		chi2x=predictedPosition->getChi2X();
		chi2y=predictedPosition->getChi2Y();
		xPositionObserved  = myTrack->getPosition(TPlaneProperties::X_COR,subjectPlane);
		yPositionObserved  = myTrack->getPosition(TPlaneProperties::Y_COR,subjectPlane);
		predictedPosition  = myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta);
		if(verbosity>3)	predictedPosition->Print();
		if(verbosity>3)	cout<<xPositionObserved<<" / "<<yPositionObserved<<endl;
		deltaX = xPositionObserved-predictedPosition->getPositionX();//X_OBS-X_Pred
		deltaY = yPositionObserved-predictedPosition->getPositionY();//Y_OBS-Y_Pred

		if(chi2x<maxChi2&&chi2y<maxChi2){
			vecXObs.push_back(xPositionObserved);
			vecYObs.push_back(yPositionObserved);
			vecDeltaX.push_back(deltaX);
			vecDeltaY.push_back(deltaY);
			vecXPred.push_back(predictedPosition->getPositionX());
			vecYPred.push_back(predictedPosition->getPositionY());
			vecXChi2.push_back(predictedPosition->getChi2X());
			vecYChi2.push_back(predictedPosition->getChi2Y());
		}

		predictedPosition->Delete();
	}
	cout<< "use "<<vecDeltaX.size()<<" of "<<events.size()<< " Events ( "<<(Float_t)vecDeltaX.size()/(Float_t)events.size()*100.<<" %)"<<endl;

	TResidual res = this->calculateResidual(cor,&vecXPred,&vecDeltaX,&vecYPred,&vecDeltaY);
	res.SetTestResidual(false);
	this->CreatePlots(cor,subjectPlane,refPlaneString.str(),bPlot,bAlign);


	vecXObs.clear();
	vecYObs.clear();
	vecDeltaX.clear();
	vecDeltaY.clear();
	vecXPred.clear();
	vecYPred.clear();
	vecXChi2.clear();
	vecYChi2.clear();
	if(verbosity)res.Print(1);
	return res;
}




