/*
 * TSelectionClass.cpp
 *
 *  Created on: 02.12.2011
 *      Author: bachmair
 */

#include "../include/TSelectionClass.hh"

TSelectionClass::TSelectionClass(TSettings* settings) {
	// TODO Auto-generated constructor stub
	cout<<"\n\n\n**********************************************************"<<endl;
	cout<<"************TSelectionClass::TSelectionClass**************"<<endl;
	cout<<"**********************************************************"<<endl;
	this->settings=settings;
	cout<<settings->getRunNumber()<<endl;

	// TODO Auto-generated constructor stub
	sys = gSystem;
	runString.str("");
	runString<<settings->getRunNumber();
	sys->MakeDirectory(runString.str().c_str());

	createdNewTree=false;
	createdNewFile=false;
	selectionTree=NULL;
	selectionFile=NULL;

	sys->cd(runString.str().c_str());
	rawfilepath<<"rawData."<<settings->getRunNumber()<<".root";
	pedestalfilepath.str("");
	pedestalfilepath<<"pedestalData."<<settings->getRunNumber()<<".root";
	clusterfilepath<<sys->pwd()<<"/clusterData."<<settings->getRunNumber()<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<"\""<<clusterfilepath.str()<<"\""<<endl;
	cout<<"OPEN TADCEventReader"<<flush;
	eventReader=new TADCEventReader(clusterfilepath.str(),settings->getRunNumber());
	cout<<" DONE"<<endl;
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("selections");
	sys->cd("selections");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	sys->cd("..");
	cout<<"HISTSAVER:"<<sys->pwd()<<endl;

			verbosity=1;

	createdTree=false;
	cout<<"Fiducial Cut:\n\n\tAccept following Range in Silicon Planes: "<<endl;
	cout<<"\t\tX: "<<settings->getSi_avg_fidcut_xlow()<<"/"<<settings->getSi_avg_fidcut_xhigh()<<endl;
	cout<<"\t\tY: "<<settings->getSi_avg_fidcut_ylow()<<"/"<<settings->getSi_avg_fidcut_yhigh()<<endl;
	cout<<" for Alignment use "<<settings->getAlignment_training_track_fraction()*100<<" % of the events." <<endl;
	nUseForAlignment=0;
	nUseForAnalysis=0;
	nUseForSiliconAlignment=0;
	nValidButMoreThanOneDiaCluster=0;
	nValidSiliconNoDiamondHit=0;
	nNoValidSiliconTrack=0;
	nValidSiliconAndDiamondCluster=0;
	nValidSiliconTrack=0;
	nValidDiamondTrack=0;
	initialiseHistos();
}

TSelectionClass::~TSelectionClass() {
	// TODO Auto-generated destructor stub
	saveHistos();
	cout<<"\n\nClosing TSelectionClass"<<endl;
	selectionFile->cd();
	if(selectionTree!=NULL&&this->createdTree){
		cout<<"CLOSING TREE"<<endl;
		cout<<eventReader->getTree()->GetName()<<" "<<clusterfilepath.str().c_str()<<endl;
		selectionTree->AddFriend(eventReader->getTree()->GetName(),clusterfilepath.str().c_str());
		cout<<"pedestalTree"<<" "<<pedestalfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("pedestalTree",pedestalfilepath.str().c_str());
		cout<<"rawTree"<<" "<<rawfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("rawTree",rawfilepath.str().c_str());
		cout<<"save selectionTree: "<<selectionTree->GetListOfFriends()->GetEntries()<<endl;
		selectionTree->Write();
		createCutFlowDiagramm();
	}

	selectionFile->Close();
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TSelectionClass::MakeSelection()
{
	MakeSelection(eventReader->GetEntries());
}



void TSelectionClass::MakeSelection(UInt_t nEvents)
{
	if(nEvents==0)
		this->nEvents=eventReader->GetEntries();
	else if(nEvents>eventReader->GetEntries()){
		cerr<<"nEvents is bigger than entries in eventReader tree: \""<<eventReader->getTree()->GetName()<<"\""<<endl;
	}
	else
		this->nEvents=nEvents;
	createdTree=createSelectionTree(nEvents);
	if(!createdTree) return;
	this->setBranchAdressess();
	nUseForAlignment=0;
	nUseForAnalysis=0;
	nUseForSiliconAlignment=0;
	nValidButMoreThanOneDiaCluster=0;
	nValidSiliconNoDiamondHit=0;
	nNoValidSiliconTrack=0;
	nValidSiliconTrack=0;
	nValidSiliconAndDiamondCluster=0;
	nValidDiamondTrack=0;
	cout<<"start selection with "<<nEvents<<" Events, training fraction: "<<settings->getAlignment_training_track_fraction()*100.<<"%"<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100,verbosity>=20);
		eventReader->LoadEvent(nEvent);
		if(verbosity>10)cout<<"Loaded Event "<<nEvent<<flush;
		resetVariables();
		if(verbosity>10)cout<<"."<<flush;
		setVariables();
		if(verbosity>10)cout<<"."<<flush;
		selectionTree->Fill();
		if(verbosity>10)cout<<"DONE"<<endl;
	}
}
bool TSelectionClass::createSelectionTree(int nEvents)
{

	cout<<"TSelectionClass::checkTree"<<endl;
	bool createdNewFile=false;
	bool createdNewTree=false;
	stringstream selectionfilepath;
	sys->cd(	runString.str().c_str());
	selectionfilepath<<sys->pwd();
	selectionfilepath<<"/selectionData."<<settings->getRunNumber()<<".root";
	cout<<"Try to open \""<<selectionfilepath.str()<<"\""<<endl;
	selectionFile=new TFile(selectionfilepath.str().c_str(),"READ");
	if(selectionFile->IsZombie()){
		cout<<"selectionfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		selectionFile= new TFile(selectionfilepath.str().c_str(),"CREATE");
		cout<<"DONE"<<flush;
		selectionFile->cd();
	}
	else{
		createdNewFile=false;
		cout<<"File exists"<<endl;
	}
	selectionFile->cd();
	cout<<"get Tree"<<endl;
	stringstream treeDescription;
	treeDescription<<"Selection Data of run "<<settings->getRunNumber();
	cout<<"get Tree2"<<endl;
	selectionFile->GetObject("selectionTree",selectionTree);
	cout<<"check Selection Tree:"<<selectionTree<<endl;
	if(selectionTree!=NULL){
		cout<<"File and Tree Exists... \t"<<flush;
		if(selectionTree->GetEntries()>=nEvents){
			createdNewTree=false;
			selectionTree->GetEvent(0);
				return false;
		}
		else{
			cout<<"selectionTree.events !- nEvents"<<flush;
			selectionTree->Delete();
			selectionTree=NULL;
		}
	}

	if(selectionTree==NULL){
		this->nEvents=nEvents;
		cout<<"selectionTree does not exists, close file"<<endl;
		delete selectionFile;
		cout<<"."<<endl;
		cout<<selectionfilepath.str().c_str()<<endl;
		selectionFile=new TFile(selectionfilepath.str().c_str(),"RECREATE");
		cout<<"."<<endl;
		this->selectionTree=new TTree("selectionTree",treeDescription.str().c_str());
		cout<<"."<<endl;
		createdNewTree=true;
		cout<<"\n\n\n***************************************************************\n";
		cout<<"***************************************************************\n";
		cout<<"there exists no tree:\'selectionTree\"\tcreate new one."<<selectionTree<<"\n";
		cout<<"***************************************************************\n";
		cout<<"***************************************************************\n"<<endl;
	}

	return createdNewTree;
}



void TSelectionClass::resetVariables(){

	isDetMasked = false;//one of the Silicon Planes contains a Cluster with a masked channel
	isDiaMasked.clear();
	nDiamondHits=0;
	hasValidSiliconTrack=true;; //One and only one cluster in each silicon plane;
	useForAnalysis=false;
	useForAlignment=false;
	useForSiliconAlignment=false;
}


/**
 * In this function all
 */
void TSelectionClass::setVariables(){
	if(verbosity>1)cout<<"setVariables\n"<<endl;

	for(UInt_t det=0;det<TPlaneProperties::getNSiliconDetectors();det++){
		bool oneAndOnlyOne = (eventReader->getNClusters(det)==1);
		if(verbosity>10)cout<<"DET "<<det<<": "<<oneAndOnlyOne<<" "<<flush;
		hasValidSiliconTrack=hasValidSiliconTrack&&oneAndOnlyOne;
		if(verbosity>10)cout<<"."<<checkDetMasked(det)<<flush;
		isDetMasked+=checkDetMasked(det);
		if(verbosity>10)cout<<"."<<isDetMasked<<endl;
	}
	for(UInt_t cl=0;cl<eventReader->getNClusters(TPlaneProperties::getDetDiamond());cl++){
		isDiaMasked.push_back(checkDetMasked(8,cl));
		if(verbosity>10)cout<<isDiaMasked[cl]<<endl;
	}

	nDiamondHits=eventReader->getNClusters(TPlaneProperties::getDetDiamond());
	isInFiducialCut=true;
	Float_t fiducialValueX=0;
	Float_t fiducialValueY=0;
	if(hasValidSiliconTrack){
		for(UInt_t plane=0;plane<4;plane++){
			fiducialValueX+=eventReader->getCluster(plane,TPlaneProperties::X_COR,0).getPosition();
			fiducialValueY+=eventReader->getCluster(plane,TPlaneProperties::Y_COR,0).getPosition();
		}
		fiducialValueX/=4.;
		fiducialValueY/=4.;
		isInFiducialCut=isInFiducialCut&&fiducialValueX>settings->getSi_avg_fidcut_xlow();
		isInFiducialCut=isInFiducialCut&&fiducialValueX<settings->getSi_avg_fidcut_xhigh();
		isInFiducialCut=isInFiducialCut&&fiducialValueY>settings->getSi_avg_fidcut_ylow();
		isInFiducialCut=isInFiducialCut&&fiducialValueY<settings->getSi_avg_fidcut_yhigh();
		if(verbosity>10)cout<<"fidCut:"<<fiducialValueX<<"/"<<fiducialValueY<<": Fidcut:"<<isInFiducialCut<<endl;
		hFiducialCutSilicon->Fill(fiducialValueX,fiducialValueY);
	}
	else
		isInFiducialCut=false;
	bool isValidSiliconTrack = isInFiducialCut&&hasValidSiliconTrack&&!isDetMasked;
	bool isValidDiamondEvent = nDiamondHits==1&&!checkDetMasked(TPlaneProperties::getDetDiamond());
	bool validMoreThanOneClusterDiamondevent = nDiamondHits>=1&&!checkDetMasked(TPlaneProperties::getDetDiamond());
	useForSiliconAlignment= isValidSiliconTrack;
	useForAlignment=useForSiliconAlignment&&isValidDiamondEvent;
	float relativeEventNumber = (float)nEvent/(float)nEvents;

//	useForAlignment=isInFiducialCut&&hasValidSiliconTrack&&nDiamondHits==1&&!checkDetMasked(TPlaneProperties::getDetDiamond());
	useForAnalysis=useForAlignment;

//	if(hasValidSiliconTrack&&isInFiducialCut){
//		if(verbosity)cout<<nEvent<<":\t"<<flush;
//		if(verbosity)printf("%5.1f %5.1f\t",fiducialValueX,fiducialValueY);
//		//if(verbosity)cout<<<<" "<<<<" "<<isInFiducialCut<<"\t"<<flush;
//		if(verbosity)cout<<isDetMasked<<" "<<eventReader->getNClusters(8)<<" "<<checkDetMasked(8)<<" "<<nDiamondHits<<endl;
//	}

	double fraction =settings->getAlignment_training_track_fraction();
	useForSiliconAlignment = isValidSiliconTrack&& !isValidDiamondEvent;// one and only one hit in silicon but not exactly one hit in diamond
	useForAlignment = isValidDiamondEvent && isValidSiliconTrack//one and only one hit in all detectors (also diamond)
					  && relativeEventNumber<fraction;			//and reltative event Number smaller than fraction
	useForAnalysis=isValidDiamondEvent && isValidSiliconTrack&& relativeEventNumber>fraction;
	if(useForAnalysis){
		nUseForAnalysis++;
		nValidSiliconAndDiamondCluster++;
	}
	if(useForAlignment){
		nUseForAlignment++;
		nValidSiliconAndDiamondCluster++;
	}
	if(useForSiliconAlignment)
		nUseForSiliconAlignment++;
	if(useForSiliconAlignment&&validMoreThanOneClusterDiamondevent){
		nValidButMoreThanOneDiaCluster++;
		nValidSiliconAndDiamondCluster++;
	}
	if(isValidSiliconTrack&&!validMoreThanOneClusterDiamondevent&&!isValidDiamondEvent)
		nValidSiliconNoDiamondHit++;
	if(isValidSiliconTrack&&isValidDiamondEvent)
		nValidDiamondTrack++;

	if(!isValidSiliconTrack)
		nNoValidSiliconTrack++;
	else
		nValidSiliconTrack++;
	//else cout<<nEvent<<"\tuseNOTforAlignemnt..."<<endl;
//		UInt_t nDiamondHits; //number of  in diamond plane;
//		bool isInFiducialCut; //if hasValidSiliconTrack avarage of x and y of all planes is in fidcut region
}



bool TSelectionClass::checkDetMasked(UInt_t det){
	bool isMasked=false;

	for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
		isMasked=isMasked||checkDetMasked(det,cl);
	}
	return isMasked;
}

bool TSelectionClass::checkDetMasked(UInt_t det,UInt_t cl){
	if(verbosity>20)cout<<"check if det Masked"<<endl;
	bool isMasked=false;

	if(cl<eventReader->getNClusters(det)){
		if(verbosity>20)cout<<"getCLuster"<<flush;
		TCluster cluster = eventReader->getCluster(det,cl);
		if(verbosity>20)cout<<"."<<flush;
		UInt_t min = cluster.getSmallestChannelNumber();
		if(verbosity>20)cout<<"."<<min<<flush;
		UInt_t max = cluster.getHighestChannelNumber();
		if(verbosity>20)cout<<":"<<max<<flush;
		for(UInt_t ch = max; ch<=max;ch++){
			if(verbosity>20)cout<<"ch"<<ch<<" "<<flush;
			isMasked=isMasked||settings->getDet_channel_screen(det).isScreened(ch);
		}
	}
	else
		return true;
	return isMasked;

	return false;
}

void TSelectionClass::setBranchAdressess(){
	selectionTree->Branch("nDiamondHits",&nDiamondHits,"nDiamondHits/i");
	selectionTree->Branch("isInFiducialCut",&isInFiducialCut,"isInFiducialCut/O");
	selectionTree->Branch("isDetMasked",&isDetMasked,"isDetMasked/O");
	selectionTree->Branch("hasValidSiliconTrack",&hasValidSiliconTrack,"hasValidSiliconTrack/O");
	selectionTree->Branch("isDiaMasked",&this->isDiaMasked,"isDiaMasked");
	selectionTree->Branch("useForSiliconAlignment",&this->useForSiliconAlignment,"useForSiliconAlignment/O");
	selectionTree->Branch("useForAlignment",&this->useForAlignment,"useForAlignment/O");
	selectionTree->Branch("useForAnalysis",&this->useForAnalysis,"useForAnalysis/O");
}

void TSelectionClass::initialiseHistos()
{
	std::string name = "hFidCutSilicon_OneAndOnlyOneCluster";
	hFiducialCutSilicon = new TH2F(name.c_str(),name.c_str(),512,0,256,512,0,256);
}



void TSelectionClass::createCutFlowDiagramm()
{
	char output[4000];
	int n=0;
	n+=sprintf(&output[n],"Finished with Selection with alignment training fraction of %f%%\n",settings->getAlignment_training_track_fraction()*100.);
	n+=sprintf(&output[n],"Selection Result: \n\tfor Silicon Alignment: %4.1f %%  %6d\n",((float)nUseForSiliconAlignment*100./(Float_t)nEvents),nUseForSiliconAlignment);
	n+=sprintf(&output[n],"\tfor Diamond Alignment: %4.1f %%  %6d\n",(float)nUseForAlignment*100./(Float_t)nEvents,nUseForAlignment);
	n+=sprintf(&output[n],"\tfor Diamond  Analysis: %4.1f %%  %6d\n",(float)nUseForAnalysis*100./(Float_t)nEvents,nUseForAnalysis);
	n+=sprintf(&output[n],"\nCUT-FLOW:\n");
	n+=sprintf(&output[n],"AllEvents: %6d ------>%6d (%4.1f%%) no valid Silicon Track\n",nEvents,nNoValidSiliconTrack,(float)nNoValidSiliconTrack*100./(float)nEvents);
	n+=sprintf(&output[n],"                    |\n");
	n+=sprintf(&output[n],"                    L--->%6d (%4.1f%%) valid Silicon Track\n",nValidSiliconTrack,(float)nValidSiliconTrack*100./(float)nEvents);
	n+=sprintf(&output[n],"                              |\n");
	n+=sprintf(&output[n],"                              L--->%6d (%4.1f%%) no Diamond Hit (absolute %4.1f%%)\n",nValidSiliconNoDiamondHit,(float)nValidSiliconNoDiamondHit*100./(float)nValidSiliconTrack,(float)nValidSiliconNoDiamondHit*100./(float)nEvents);
	n+=sprintf(&output[n],"                              |\n");
	n+=sprintf(&output[n],"                              L--->%6d (%4.1f%%) at least one Diamond Hit\n",nValidSiliconAndDiamondCluster,(float)nValidSiliconAndDiamondCluster*100./(float)nValidSiliconTrack);
	n+=sprintf(&output[n],"                                        |\n");
	n+=sprintf(&output[n],"                                        L--->%6d (%4.1f%%) more than one Diamond Hit\n",nValidButMoreThanOneDiaCluster,(float)nValidButMoreThanOneDiaCluster*100./(float)nValidSiliconAndDiamondCluster);
	n+=sprintf(&output[n],"                                        |\n");
	n+=sprintf(&output[n],"                                        L--->%6d (%4.1f%%) exactly one Diamond Hit\n",nValidDiamondTrack,(float)nValidDiamondTrack*100./(float)nValidSiliconAndDiamondCluster);
	n+=sprintf(&output[n],"                                                  |\n");
	n+=sprintf(&output[n],"                                                  L--->%6d (%4.1f%%) Alignment (absolute: %4.1f%%)\n",nUseForAlignment,(float)nUseForAlignment*100./(float)nValidDiamondTrack,(float)nUseForAlignment*100./(float)nEvents);
	n+=sprintf(&output[n],"                                                  |\n");
	n+=sprintf(&output[n],"                                                  L--->%6d (%4.1f%%) Analysis (absolute: %4.1f%%)\n",nUseForAnalysis,(float)nUseForAnalysis*100./(float)nValidDiamondTrack,(float)nUseForAnalysis*100./(float)nEvents);
	cout<<output<<endl;
	histSaver->SaveStringToFile("cutFlow.txt",output);
}

void TSelectionClass::saveHistos()
{
	cout<<"save Histo: "<<hFiducialCutSilicon->GetTitle()<<endl;
	histSaver->SaveHistogram(hFiducialCutSilicon);
	delete hFiducialCutSilicon;
}


