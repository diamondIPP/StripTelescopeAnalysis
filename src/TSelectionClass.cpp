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
	if(settings==0)exit(-1);
	this->settings=settings;
	this->results=0;
	cout<<settings->getRunNumber()<<endl;

	// TODO Auto-generated constructor stub
	sys = gSystem;
	cout<<"goToClusterTree"<<endl;
  settings->goToClusterTreeDir();
  fiducialCuts=0;
	createdNewTree=false;
	createdNewFile=false;
	selectionTree=NULL;
	selectionFile=NULL;
	htmlSelection = new THTMLSelection(settings);
  cout<<"OPEN TADCEventReader"<<flush;
  cout<<"\ngoToSelectionTreeDir"<<endl;
  settings->goToSelectionTreeDir();
  cout<<"open Tree:"<<endl;
	eventReader=new TADCEventReader(settings->getClusterTreeFilePath(),settings->getRunNumber());
	cout<<" DONE"<<endl;

	histSaver=new HistogrammSaver();
	cout<<"goToSelectionDir"<<endl;
  settings->goToSelectionDir();
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
  htmlSelection->setFileGeneratingPath(sys->pwd());
  cout<<"goToSelectionTREEDir"<<endl;
  settings->goToSelectionTreeDir();
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
	nSiliconTrackNotFiducialCut=0;
	nValidDiamondTrack=0;
	initialiseHistos();
	htmlSelection->createFiducialCuts();
}

TSelectionClass::~TSelectionClass() {
	// TODO Auto-generated destructor stub
	cout<<"\n\nClosing TSelectionClass"<<endl;
	selectionFile->cd();
	if(selectionTree!=NULL&&this->createdTree){
		saveHistos();
		cout<<"CLOSING TREE"<<endl;
		settings->goToAlignmentRootDir();
		cout<<"\t"<<eventReader->getTree()->GetName()<<" "<<settings->getClusterTreeFilePath()<<endl;
		selectionTree->AddFriend("clusterTree",settings->getClusterTreeFilePath().c_str());

		cout<<"\t"<<"pedestalTree"<<" "<<pedestalfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("pedestalTree",settings->getPedestalTreeFilePath().c_str());

		cout<<"\t"<<"rawTree"<<" "<<rawfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("rawTree",settings->getRawTreeFilePath().c_str());

		cout<<"\n\n\t"<<"save selectionTree: "<<selectionTree->GetListOfFriends()->GetEntries()<<endl;
		selectionFile->cd();
		cout<<"\t"<<"WRITE TREE: "<<flush;
		int retVal = selectionTree->Write();
		cout<<retVal<<endl;
		htmlSelection->generateHTMLFile();
	}
	selectionFile->Close();
	delete eventReader;
	delete histSaver;
	delete htmlSelection;
	cout<<"goToOutputDir"<<endl;
	settings->goToOutputDir();
}

void TSelectionClass::MakeSelection()
{
	MakeSelection(eventReader->GetEntries());
}



void TSelectionClass::MakeSelection(UInt_t nEvents)
{cout<<"Make Selection"<<endl;
	if(nEvents==0)
		this->nEvents=eventReader->GetEntries();
	else if(nEvents>eventReader->GetEntries()){
		cerr<<"nEvents is bigger than entries in eventReader tree: \""<<eventReader->getTree()->GetName()<<"\""<<endl;
	}
	else
		this->nEvents=nEvents;
	cout<<"goToSelectionTreeDir"<<endl;
	settings->goToSelectionTreeDir();
	histSaver->SetNumberOfEvents(this->nEvents);
	createdTree=createSelectionTree(nEvents);
	if(!createdTree) return;
	this->setBranchAdressess();
	  createFiducialCut();
	hFiducialCutSilicon->Reset();
	hFiducialCutSiliconDiamondHit->Reset();
	hAnalysisFraction->Reset();
	nUseForAlignment=0;
	nUseForAnalysis=0;
	nUseForSiliconAlignment=0;
	nValidButMoreThanOneDiaCluster=0;
	nValidSiliconNoDiamondHit=0;
	nNoValidSiliconTrack=0;
	nValidSiliconTrack=0;
	nValidSiliconAndDiamondCluster=0;
	nValidDiamondTrack=0;
	nSiliconTrackNotFiducialCut=0;
	nToBigDiamondCluster=0;
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
	createCutFlowDiagramm();
}
bool TSelectionClass::createSelectionTree(int nEvents)
{

	cout<<"TSelectionClass::checkTree"<<endl;
	bool createdNewFile=false;
	bool createdNewTree=false;
	cout<<"goToSelection Tree:"<<endl;
	settings->goToSelectionTreeDir();
	selectionFile=new TFile(settings->getSelectionTreeFilePath().c_str(),"READ");
	if(selectionFile->IsZombie()){
		cout<<"selectionfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		selectionFile= new TFile(settings->getSelectionTreeFilePath().c_str(),"CREATE");
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
	cout<<sys->pwd()<<endl;
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
		selectionFile=new TFile(settings->getSelectionTreeFilePath().c_str(),"RECREATE");
		selectionFile->cd();
		cout<<"."<<endl;
		this->selectionTree=new TTree("selectionTree",treeDescription.str().c_str());
		cout<<"."<<endl;
		createdNewTree=true;
		cout<<"\n***************************************************************\n";
		cout<<"there exists no tree:\'selectionTree\"\tcreate new one."<<selectionTree<<"\n";
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
	bool isSilSaturated=false;
	if(hasValidSiliconTrack){
	  for(UInt_t det=0;det<TPlaneProperties::getNSiliconDetectors();det++)
	    isSilSaturated = isSilSaturated ||isSaturated(det);
	}
//	if(isSilSaturated&&hasValidSiliconTrack){
//	  hasValidSiliconTrack==false;
//	}
	if(hasValidSiliconTrack){
		for(UInt_t plane=0;plane<4;plane++){
			fiducialValueX+=eventReader->getCluster(plane,TPlaneProperties::X_COR,0).getPosition();
			fiducialValueY+=eventReader->getCluster(plane,TPlaneProperties::Y_COR,0).getPosition();
		}
		fiducialValueX/=4.;
		fiducialValueY/=4.;
		isInFiducialCut = fiducialCuts->isInFiducialCut(fiducialValueX,fiducialValueY);
//		isInFiducialCut=isInFiducialCut&&fiducialValueX>settings->getSi_avg_fidcut_xlow();
//		isInFiducialCut=isInFiducialCut&&fiducialValueX<settings->getSi_avg_fidcut_xhigh();
//		isInFiducialCut=isInFiducialCut&&fiducialValueY>settings->getSi_avg_fidcut_ylow();
//		isInFiducialCut=isInFiducialCut&&fiducialValueY<settings->getSi_avg_fidcut_yhigh();
		if(verbosity>10)cout<<"fidCut:"<<fiducialValueX<<"/"<<fiducialValueY<<": Fidcut:"<<isInFiducialCut<<endl;
		hFiducialCutSilicon->Fill(fiducialValueX,fiducialValueY);
	}
	else
		isInFiducialCut=false;
	bool isSiliconTrackNotFiducialCut= !isInFiducialCut&&hasValidSiliconTrack&&!isDetMasked;
	isValidSiliconTrack = isInFiducialCut&&hasValidSiliconTrack&&!isDetMasked;
	isDiaSaturated=this->isSaturated(TPlaneProperties::getDetDiamond());
	bool isValidDiamondEvent = nDiamondHits==1&&!checkDetMasked(TPlaneProperties::getDetDiamond())&&!isDiaSaturated;
	bool isDiaOneHitButToBigCluster = false;
	bool validMoreThanOneClusterDiamondevent = nDiamondHits>=1&&!checkDetMasked(TPlaneProperties::getDetDiamond());
	useForSiliconAlignment= isValidSiliconTrack;

	if(nDiamondHits==1){
	  nDiaClusterSize =eventReader->getClusterSize(TPlaneProperties::getDetDiamond(),0);
	  if(nDiaClusterSize>=3){
	    isValidDiamondEvent=false;
	    isDiaOneHitButToBigCluster=true;
	  }
	}
	else nDiaClusterSize = -1;

	useForAlignment=useForSiliconAlignment&&isValidDiamondEvent;
	float relativeEventNumber = (float)nEvent/(float)nEvents;
	isDiaOneHitButToBigCluster = isDiaOneHitButToBigCluster&&useForSiliconAlignment;
	useForAnalysis=useForAlignment;


	double fraction =settings->getAlignment_training_track_fraction();
	useForSiliconAlignment = isValidSiliconTrack&& !isValidDiamondEvent;// one and only one hit in silicon but not exactly one hit in diamond
	useForAlignment = isValidDiamondEvent && isValidSiliconTrack//one and only one hit in all detectors (also diamond)
					  && relativeEventNumber<fraction;			//and reltative event Number smaller than fraction
	useForAnalysis=isValidDiamondEvent && isValidSiliconTrack&& relativeEventNumber>fraction;
	if((isValidSiliconTrack||isSiliconTrackNotFiducialCut)&&nDiamondHits>0)
		hFiducialCutSiliconDiamondHit->Fill(fiducialValueX,fiducialValueY);
	if(isSiliconTrackNotFiducialCut)
		nSiliconTrackNotFiducialCut++;
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
	if(isDiaOneHitButToBigCluster)
	  nToBigDiamondCluster++;

	if(!isValidSiliconTrack)
		nNoValidSiliconTrack++;
	else
		nValidSiliconTrack++;
	if(useForAnalysis||useForAlignment)
	  hAnalysisFraction->Fill(nEvent);
	//else cout<<nEvent<<"\tuseNOTforAlignemnt..."<<endl;
//		UInt_t nDiamondHits; //number of  in diamond plane;
//		bool isInFiducialCut; //if hasValidSiliconTrack avarage of x and y of all planes is in fidcut region
//	)
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
	selectionTree->Branch("diaClusterSize",&this->nDiaClusterSize,"diaClusterSize/I");
	selectionTree->Branch("isDiaSaturated",&this->isDiaSaturated,"isDiaSaturated/O");
}

void TSelectionClass::initialiseHistos()
{
	std::string name = "hFidCutSilicon_OneAndOnlyOneCluster";
	hFiducialCutSilicon = new TH2F(name.c_str(),name.c_str(),512,0,256,512,0,256);
	hFiducialCutSilicon->GetYaxis()->SetTitle("yCoordinate in Channels");
	hFiducialCutSilicon->GetXaxis()->SetTitle("xCoordinate in Channels");

	std::string name2 = "hFidCutSilicon_OneAndOnlyOneCluster_DiamondCluster";
	hFiducialCutSiliconDiamondHit = new TH2F(name2.c_str(),name2.c_str(),512,0,256,512,0,256);

	hFiducialCutSiliconDiamondHit->GetYaxis()->SetTitle("yCoordinate in Channels");
	hFiducialCutSiliconDiamondHit->GetXaxis()->SetTitle("xCoordinate in Channels");
	int nEvents= eventReader->GetEntries();

	int i=nEvents/1000;
	i++;
	nEvents = (i)*1000;
	hAnalysisFraction = new TH1F("hAnalysisFraction","hAnalysisFraction",i,0,nEvents);
	hAnalysisFraction->SetTitle("Fraction of Events for Analysis");
	hAnalysisFraction->GetXaxis()->SetTitle("event no");
	hAnalysisFraction->GetYaxis()->SetTitle("fraction of daimond + silicon hit events (%)");
}



void TSelectionClass::createCutFlowDiagramm()
{
  if(results!=0){
    if(!results->IsZombie()){
      results->setAllEvents(nEvents);
      results->setNoSiliconHit(nNoValidSiliconTrack-nSiliconTrackNotFiducialCut);
      results->setOneAndOnlyOneSiliconNotFiducialCut(nSiliconTrackNotFiducialCut);
      results->setValidSiliconTrack(nValidSiliconTrack);
      results->setNoDiamondHit(nValidSiliconNoDiamondHit);
      results->setMoreThanOneDiamondHit(nValidButMoreThanOneDiaCluster);
      results->setExactlyOneDiamondHit(nValidDiamondTrack);
      results->setUseForAlignment(nUseForAlignment);
      results->setUseForAnalysis(nUseForAnalysis);
    }
  }
	char output[4000];
	int n=0;
	n+=sprintf(&output[n],"Finished with Selection with alignment training fraction of %f%%\n",settings->getAlignment_training_track_fraction()*100.);
	n+=sprintf(&output[n],"Selection Result: \n\tfor Silicon Alignment: %4.1f %%  %6d\n",((float)nUseForSiliconAlignment*100./(Float_t)nEvents),nUseForSiliconAlignment);
	n+=sprintf(&output[n],"\tfor Diamond Alignment: %4.1f %%  %6d\n",(float)nUseForAlignment*100./(Float_t)nEvents,nUseForAlignment);
	n+=sprintf(&output[n],"\tfor Diamond  Analysis: %4.1f %%  %6d\n",(float)nUseForAnalysis*100./(Float_t)nEvents,nUseForAnalysis);
	n+=sprintf(&output[n],"\nCUT-FLOW:\n");
	n+=sprintf(&output[n],"AllEvents: %6d ------>%6d (%4.1f%%) no only one and only one Silicon Hit\n",nEvents,(nNoValidSiliconTrack-nSiliconTrackNotFiducialCut),(float)(nNoValidSiliconTrack-nSiliconTrackNotFiducialCut)*100./(float)nEvents);
	n+=sprintf(&output[n],"                    |\n");
	n+=sprintf(&output[n],"                    L--->%6d (%4.1f%%) one and only one silicon hit, not in Fiducial Cut\n",nSiliconTrackNotFiducialCut,(float)nSiliconTrackNotFiducialCut*100./(float)nEvents);
	n+=sprintf(&output[n],"                    |\n");
	n+=sprintf(&output[n],"                    L--->%6d (%4.1f%%) valid Silicon Track\n",nValidSiliconTrack,(float)nValidSiliconTrack*100./(float)nEvents);
	n+=sprintf(&output[n],"                              |\n");
	n+=sprintf(&output[n],"                              L--->%6d (%4.1f%%) no Diamond Hit (absolute %4.1f%%)\n",nValidSiliconNoDiamondHit,(float)nValidSiliconNoDiamondHit*100./(float)nValidSiliconTrack,(float)nValidSiliconNoDiamondHit*100./(float)nEvents);
	n+=sprintf(&output[n],"                              |\n");
	n+=sprintf(&output[n],"                              L--->%6d (%4.1f%%) at least one Diamond Hit\n",nValidSiliconAndDiamondCluster,(float)nValidSiliconAndDiamondCluster*100./(float)nValidSiliconTrack);
	n+=sprintf(&output[n],"                                        |\n");
	n+=sprintf(&output[n],"                                        L--->%6d (%4.1f%%) more than one Diamond Hit\n",nValidButMoreThanOneDiaCluster,(float)nValidButMoreThanOneDiaCluster*100./(float)nValidSiliconAndDiamondCluster);
	n+=sprintf(&output[n],"                                        |\n");
	n+=sprintf(&output[n],"                                        L--->%6d (%4.1f%%) toBigClusters (absolute: %4.1f%%)\n",nToBigDiamondCluster,(float)nToBigDiamondCluster*100./(float)nValidSiliconAndDiamondCluster,(float)nToBigDiamondCluster*100./(float)nEvents);
	n+=sprintf(&output[n],"                                        |\n");
	n+=sprintf(&output[n],"                                        L--->%6d (%4.1f%%) exactly one Diamond Hit\n",nValidDiamondTrack,(float)nValidDiamondTrack*100./(float)nValidSiliconAndDiamondCluster);

	n+=sprintf(&output[n],"                                                  |\n");
	n+=sprintf(&output[n],"                                                  L--->%6d (%4.1f%%) Alignment (absolute: %4.1f%%)\n",nUseForAlignment,(float)nUseForAlignment*100./(float)nValidDiamondTrack,(float)nUseForAlignment*100./(float)nEvents);
	n+=sprintf(&output[n],"                                                  |\n");
	n+=sprintf(&output[n],"                                                  L--->%6d (%4.1f%%) Analysis (absolute: %4.1f%%)\n",nUseForAnalysis,(float)nUseForAnalysis*100./(float)nValidDiamondTrack,(float)nUseForAnalysis*100./(float)nEvents);
	cout<<output<<endl;
	histSaver->SaveStringToFile("cutFlow.txt",output);
	string cutFlow;
	cutFlow = output;
	htmlSelection->createCutFlowGraph(cutFlow);



	Double_t values [] = {(nNoValidSiliconTrack-nSiliconTrackNotFiducialCut),nSiliconTrackNotFiducialCut,nValidSiliconNoDiamondHit+nValidButMoreThanOneDiaCluster,nUseForAlignment,nUseForAnalysis};
	Int_t colors[] = {kRed,kRed+1,kBlue,kYellow,kGreen};
	Int_t nvals = sizeof(values)/sizeof(values[0]);
	TCanvas *cpieMain = new TCanvas("cMainCutFlow","Main Cut Flow",700,700);
	cpieMain->cd();
	TPie *pie4 = new TPie("pieCutFlow","cutFlow Silicon",nvals,values,colors);
	pie4->SetEntryLabel(0,"noSilTrack");
	pie4->SetEntryLabel(1,"notInFidCut");
	pie4->SetEntryLabel(2,"notExactlyOneDiamondCluster");
	pie4->SetEntryLabel(3,"useForAlignment");
	pie4->SetEntryLabel(4,"useForAnalysis");
	pie4->SetRadius(.3);
	pie4->SetLabelsOffset(.02);
	pie4->SetTextSize(pie4->GetTextSize()*0.4);
	pie4->SetLabelFormat("#splitline{%val (%perc)}{%txt}");
	pie4->SetValueFormat("%d");
	pie4->Draw("nol ");
	TLegend* legend = pie4->MakeLegend(0.7,0.7,0.99,0.95,"Silicon CutFlow");
	legend->SetFillColor(kWhite);
	legend ->Draw();
	histSaver->SaveCanvas(cpieMain);

	Double_t  valuesSilTracks[] = {nValidSiliconNoDiamondHit,nValidButMoreThanOneDiaCluster,nUseForAlignment,nUseForAnalysis};
	Int_t colorsSilTracks[] = {kBlue+1,kBlue+2,kYellow,kGreen};
	Int_t nvalsSilTracks = sizeof(valuesSilTracks)/sizeof(valuesSilTracks[0]);
	TCanvas *cPieValidSilicontTrack = new TCanvas("cPieValidSiliconTrack","CutFlow Valid Silicon Track",700,700);
	cPieValidSilicontTrack->cd();
	TPie *pieValidSilicontTrack = new TPie("pieValidSilicontTrack","CutFlow Valid Silicon Track",nvalsSilTracks,valuesSilTracks,colorsSilTracks);
	pieValidSilicontTrack->SetEntryLabel(0,"noDiamondHit");
	pieValidSilicontTrack->SetEntryLabel(1,"moreThanOneDiamondHit");
	pieValidSilicontTrack->SetEntryLabel(2,"useForAlignment");
	pieValidSilicontTrack->SetEntryLabel(3,"useForAnalysis");
	pieValidSilicontTrack->SetRadius(.25);
	pieValidSilicontTrack->SetLabelsOffset(.03);
	pieValidSilicontTrack->SetTextSize(pieValidSilicontTrack->GetTextSize()*0.35);
	pieValidSilicontTrack->SetLabelFormat("%val (%perc) - %txt ");
	pieValidSilicontTrack->SetValueFormat("%.0f");
	pieValidSilicontTrack->SetX(0.4);
	pieValidSilicontTrack->Draw("nol");
	TLegend* legendSilTrack = pieValidSilicontTrack->MakeLegend(0.76,0.65,0.99,0.97,"Diamond CutFlow");
	legendSilTrack->SetFillColor(kWhite);

	legendSilTrack ->Draw();
	histSaver->SaveCanvas(cPieValidSilicontTrack);
}

bool TSelectionClass::isSaturated(UInt_t det,UInt_t cl)
{
  if(eventReader->getNClusters(det)<=cl)
    return true;
  TCluster cluster = eventReader->getCluster(det,cl);
//  if(cluster.hasSaturatedChannels() )cout<<nEvent<<" hasSaturatedChannel"<<flush;
  for(UInt_t clPos=0;clPos<cluster.size();clPos++)
    if(cluster.getAdcValue(clPos)>=TPlaneProperties::getMaxSignalHeight(det)){
//      cout<<"\t"<<this->nEvent<<" confirmed"<<endl;
      return true;
    }
  return false;
}

void TSelectionClass::saveHistos()
{
	cout<<"save Histo: "<<hFiducialCutSilicon->GetTitle()<<endl;
	std::string name = "c";
	name.append(hFiducialCutSilicon->GetName());
	TCanvas *c1= fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSilicon);
	    /*new TCanvas(name.c_str(),hFiducialCutSilicon->GetTitle());
	c1->cd();
	hFiducialCutSilicon->Draw("colz");
	double xLow = settings->getSi_avg_fidcut_xlow();
	double xHigh = settings->getSi_avg_fidcut_xhigh();
	double yLow = settings->getSi_avg_fidcut_ylow();
	double yHigh = settings->getSi_avg_fidcut_yhigh();
	TLine* lXlower = new TLine(xLow,yLow,xLow,yHigh);
	TLine* lXhigher = new TLine(xHigh,yLow,xHigh,yHigh);
	TLine* lYlower = new TLine(xLow,yLow,xHigh,yLow);
	TLine* lYhigher = new TLine(xLow,yHigh,xHigh,yHigh);
	lXlower->Draw();
	lXhigher->Draw();
	lYlower->Draw();
	lYhigher->Draw();*/
	histSaver->SaveCanvas(c1);
	std::string name2 = "c";
	name2.append(hFiducialCutSiliconDiamondHit->GetName());
	TCanvas *c2= fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSiliconDiamondHit);
	 //   new TCanvas(name2.c_str(),hFiducialCutSiliconDiamondHit->GetTitle());
	/*c2->cd();
	hFiducialCutSiliconDiamondHit->Draw("colz");
	lXlower->Draw();
	lXhigher->Draw();
	lYlower->Draw();
	lYhigher->Draw();*/
	histSaver->SaveCanvas(c2);

//	histSaver->SaveHistogram(hFiducialCutSilicon);
	delete hFiducialCutSilicon;
	delete hFiducialCutSiliconDiamondHit;
	hAnalysisFraction->Scale(.1);
	hAnalysisFraction->SetStats(false);
	hAnalysisFraction->GetYaxis()->SetRangeUser(0,100);
	histSaver->SaveHistogram(hAnalysisFraction);
	delete hAnalysisFraction;
}

void TSelectionClass::createFiducialCut(){

  std::vector<std::pair<Float_t,Float_t> > xInt,yInt;
  xInt.push_back( make_pair(settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh()));
  yInt.push_back( make_pair(settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
  fiducialCuts = new TFidCutRegions(xInt,yInt,1);
  cout<<"Create AutoFidCut"<<endl;
  UInt_t nEvents = settings->getAutoFidCutEvents();
  if(nEvents>eventReader->GetEntries())nEvents=eventReader->GetEntries();
  cout<<" "<<nEvents<<endl;
  for(nEvent=0;nEvent<nEvents;nEvent++){
      TRawEventSaver::showStatusBar(nEvent,nEvents,100,verbosity>=20);
      eventReader->LoadEvent(nEvent);
      if(verbosity>10)cout<<"Loaded Event "<<nEvent<<flush;
      resetVariables();
      if(verbosity>10)cout<<"."<<flush;
      setVariables();
  }
//  findFiducialCut(hFiducialCutSiliconDiamondHit);
  delete fiducialCuts;
  fiducialCuts = new TFidCutRegions(hFiducialCutSiliconDiamondHit,settings->getNDiamonds(),settings->getAutoFidCutPercentage());

  histSaver->SaveCanvas(fiducialCuts->getFiducialCutCanvas(TPlaneProperties::X_COR));
  histSaver->SaveCanvas(fiducialCuts->getFiducialCutCanvas(TPlaneProperties::Y_COR));
  histSaver->SaveCanvas(fiducialCuts->getFiducialCutCanvas(TPlaneProperties::XY_COR));
  //  fiducialCuts->getFiducialCutCanvas();
}

//void TSelectionClass::findFiducialCut(TH2F* hFidCut){
//
//  TH1D* hProjX = hFiducialCutSiliconDiamondHit->ProjectionX("hFiducialCutOneDiamondHitProjX");
//  hProjX->GetXaxis()->SetTitle("Mean Silicon Value in X[strips]");
//  hProjX->GetYaxis()->SetTitle("Number Of Entries #");
//  TH1D* hProjY = hFiducialCutSiliconDiamondHit->ProjectionY("hFiducialCutOneDiamondHitProjY");
//  hProjY->GetXaxis()->SetTitle("Mean Silicon Value in Y[strips]");
//  hProjY->GetYaxis()->SetTitle("Number Of Entries #");
//
//  hProjX->SetName("hProjX");
//  hProjY->SetName("hProjY");
//  if(hProjX==0||hProjY==0)
//    return;
//  if(hProjX->IsZombie()||hProjY->IsZombie())
//    return;
//
////
////  vector<pair <Float_t, Float_t> >xIntervals = findFiducialCutIntervall( hProjX);
////  vector<pair <Float_t, Float_t> >yIntervals = findFiducialCutIntervall( hProjY);
////  chooseFidCut(xIntervals,yIntervals);
////  cout<<"HERE"<<endl;
////  TLegend *leg = new TLegend();
//
//  histSaver->SaveHistogram(hProjX);
//  cout<<"%"<<endl;
//  histSaver->SaveHistogram(hProjY);
////  cout<<"SAVED"<<endl;
//}
//
//
//void TSelectionClass::chooseFidCut(std::vector<std::pair <Float_t, Float_t> > xInt,std::vector<std::pair <Float_t, Float_t> >yInt){
//  cout<<"Found "<<xInt.size()<<" Intervals in x direction."<<endl;
//  cout<<"Found "<<yInt.size()<<" Intervals in y direction."<<endl;
//  cout<<"The RunDescription is \""<<settings->getRunDescription()<<"\""<<endl;
//  cout<<"Finding Areas for "<<settings->getNDiamonds()<<" Diamonds"<<endl;
//  if(fiducialCuts!=0) delete fiducialCuts;
//  fiducialCuts= new TFidCutRegions(xInt,yInt,settings->getNDiamonds());
//  fiducialCuts->setRunDescription(settings->getRunDescription());
//  TFiducialCut *fidCut = fiducialCuts->getFidCut(settings->getRunDescription());
//  fidCut->Print();
//  cout<<"old Fiducial Cut is "<<settings->getSi_avg_fidcut_xlow()<<":"<<settings->getSi_avg_fidcut_xhigh()<<endl;
//  settings->setFidCut(fidCut);
//  cout<<"new Fiducial Cut is "<<settings->getSi_avg_fidcut_xlow()<<":"<<settings->getSi_avg_fidcut_xhigh()<<endl;
//}
//
//void TSelectionClass::DrawFiduciaCuts(TH1D* hProj,vector< pair<Float_t,Float_t> > intervals){
//  if(hProj==0)
//    return;
//  if(hProj->IsZombie())
//    return;
//  stringstream canvasName;
//  canvasName<<"c"<<hProj->GetName();
//  TCanvas *c1 =new TCanvas(canvasName.str().c_str(),canvasName.str().c_str(),800,600);
//  c1->cd();
//  hProj->Draw();
//  vector<TBox* > boxes;
//  for(UInt_t i=0;i<intervals.size();i++){
//    TPaveText *box = new TPaveText(intervals.at(i).first,0,intervals.at(i).second,hProj->GetMaximum());
//    box->SetFillColor(kRed+i);
//    box->SetFillStyle(3013);
//    box->AddText("");
//    box->AddText("");
//    box->AddText(Form("Area\n\n%i",i));
//    boxes.push_back(box);
//    box->Draw("same");
//  }
//  histSaver->SaveCanvas(c1);
//
//}
