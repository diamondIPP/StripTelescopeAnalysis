/*
 * TSelectionClass.cpp
 *
 *  Created on: 02.12.2011
 *      Author: bachmair
 */

#include "../include/TSelectionClass.hh"

TSelectionClass::TSelectionClass(TSettings* newSettings) {
    // TODO Auto-generated constructor stub
    cout<<"\n\n\n**********************************************************"<<endl;
    cout<<"************TSelectionClass::TSelectionClass**************"<<endl;
    cout<<"**********************************************************"<<endl;
    if(newSettings==0)exit(-1);
    this->settings=newSettings;
    settings->diamondPattern.getNPatterns();
    verbosity=settings->getVerbosity();
    this->results=0;
    if(verbosity)cout<<settings->getRunNumber()<<endl;

    // TODO Auto-generated constructor stub
    sys = gSystem;
    if(verbosity)cout<<"goToClusterTree"<<endl;
    settings->goToClusterTreeDir();
    fiducialCuts=0;
    createdNewTree=false;
    createdNewFile=false;
    selectionTree=NULL;
    selectionFile=NULL;
    settings->goToSelectionTreeDir();
    htmlSelection = new THTMLSelection(settings);
    if(verbosity)cout<<"OPEN TADCEventReader"<<flush;
    if(verbosity)cout<<"\ngoToSelectionTreeDir"<<endl;
    settings->goToSelectionTreeDir();
    if(verbosity)cout<<"open Tree:"<<endl;
    eventReader=new TADCEventReader(settings->getClusterTreeFilePath(),settings);
    //settings->getRunNumber());
    if(verbosity)cout<<" DONE"<<endl;

    histSaver=new HistogrammSaver(settings);
    if(verbosity)cout<<"goToSelectionDir"<<endl;
    settings->goToSelectionDir();
    stringstream plotsPath;
    plotsPath<<sys->pwd()<<"/";
    histSaver->SetPlotsPath(plotsPath.str().c_str());
    histSaver->SetRunNumber(settings->getRunNumber());
    htmlSelection->setFileGeneratingPath(sys->pwd());
    if(verbosity)cout<<"goToSelectionTREEDir"<<endl;
    settings->goToSelectionTreeDir();
    if(verbosity)cout<<"HISTSAVER:"<<sys->pwd()<<endl;

    createdTree=false;
    cout<<"\nFiducial Cut:\n\tAccept following Range in Silicon Planes: "<<endl;
    cout<<"\t\tX: "<<settings->getSi_avg_fidcut_xlow()<<"/"<<settings->getSi_avg_fidcut_xhigh()<<endl;
    cout<<"\t\tY: "<<settings->getSi_avg_fidcut_ylow()<<"/"<<settings->getSi_avg_fidcut_yhigh()<<endl;
    cout<<"for Alignment use "<<settings->getAlignment_training_track_fraction()*100<<" % of the events." <<endl;
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
        if(verbosity)cout<<"CLOSING TREE"<<endl;
        settings->goToAlignmentRootDir();
        if(verbosity)cout<<"\t"<<eventReader->getTree()->GetName()<<" "<<settings->getClusterTreeFilePath()<<endl;
        selectionTree->AddFriend("clusterTree",settings->getClusterTreeFilePath().c_str());

        if(verbosity)cout<<"\t"<<"pedestalTree"<<" "<<pedestalfilepath.str().c_str()<<endl;
        selectionTree->AddFriend("pedestalTree",settings->getPedestalTreeFilePath().c_str());

        if(verbosity)cout<<"\t"<<"rawTree"<<" "<<rawfilepath.str().c_str()<<endl;
        selectionTree->AddFriend("rawTree",settings->getRawTreeFilePath().c_str());

        if(verbosity)cout<<"\n\n\t"<<"save selectionTree: "<<selectionTree->GetListOfFriends()->GetEntries()<<endl;
        selectionFile->cd();
        if(verbosity)cout<<"\t"<<"WRITE TREE: "<<flush;
        int retVal = selectionTree->Write();
        if(verbosity)cout<<retVal<<endl;
        htmlSelection->generateHTMLFile();
    }
    selectionFile->Close();
    delete eventReader;
    delete histSaver;
    delete htmlSelection;
    if(verbosity)cout<<"goToOutputDir"<<endl;
    settings->goToOutputDir();
}

void TSelectionClass::MakeSelection()
{
    MakeSelection(eventReader->GetEntries());
}



void TSelectionClass::MakeSelection(UInt_t nEvents){
    if(nEvents==0)
        this->nEvents=eventReader->GetEntries();
    else if(nEvents>eventReader->GetEntries()){
        cerr<<"nEvents is bigger than entries in eventReader tree: \""<<eventReader->getTree()->GetName()<<"\""<<endl;
    }
    else
        this->nEvents=nEvents;

    if(verbosity)cout<<"Make Selection"<<endl;
    if(verbosity)cout<<"goToSelectionTreeDir"<<endl;
    settings->goToSelectionTreeDir();
    histSaver->SetNumberOfEvents(this->nEvents);
    createdTree=createSelectionTree(nEvents);
    if(!createdTree) return;
    this->setBranchAdressess();
    createFiducialCut();
    hFiducialCutSilicon->Reset();
    hFiducialCutSiliconDiamondHit->Reset();
    hFiducialCutSiliconOneAndOnlyOneDiamondHit->Reset();
    hAnalysisFraction->Reset();
    hSelectedEvents->Reset();
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
    cout<<"start selection in  "<<nEvents<<" Events."<<endl;
    if(settings->getTrainingMethod()==TSettings::enumFraction)
        cout<<"Use Fraction Training, with  fraction: "<<settings->getAlignment_training_track_fraction()*100.<<"%"<<endl;
    else
        cout<<"Use the first "<<settings->getAlignmentTrainingTrackNumber()<<" Events for Alignment!"<<endl;
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
    if(verbosity)cout<<"TSelectionClass::checkTree"<<endl;
    bool createdNewFile=false;
    bool createdNewTree=false;
    if(verbosity)cout<<"\tgoToSelection Tree:"<<endl;
    settings->goToSelectionTreeDir();
    selectionFile=new TFile(settings->getSelectionTreeFilePath().c_str(),"READ");
    if(selectionFile->IsZombie() || settings->RerunSelection()){
        if(verbosity && selectionFile->IsZombie())
            cout<<"\tselectionfile does not exist, create new one..."<<endl;
        else if (settings->RerunSelection())
            cout<<"\tRerun Selection activated"<<endl;
        createdNewFile =true;
        selectionFile= new TFile(settings->getSelectionTreeFilePath().c_str(),"CREATE");
        if(verbosity)cout<<"DONE"<<flush;
        selectionFile->cd();
    }
    else{
        createdNewFile=false;
        if(verbosity)cout<<"\tFile exists"<<endl;
    }
    selectionFile->cd();
    if(verbosity)cout<<"\tget Tree"<<endl;
    stringstream treeDescription;
    treeDescription<<"Selection Data of run "<<settings->getRunNumber();
    if(verbosity)cout<<"\tget Tree2"<<endl;
    selectionFile->GetObject("selectionTree",selectionTree);
    if(verbosity)cout<<"\tcheck Selection Tree:"<<selectionTree<<endl;
    if(verbosity)cout<<sys->pwd()<<endl;
    if(selectionTree!=NULL){
        if(verbosity)cout<<"\tFile and Tree Exists... \t"<<selectionTree->GetEntries()<<" Events\t"<<flush;
        if(selectionTree->GetEntries()>=nEvents){
            createdNewTree=false;
            selectionTree->GetEvent(0);
            return false;
        }
        else{
            if(verbosity)cout<<"\tselectionTree.events !- nEvents"<<flush;
            selectionTree->Delete();
            selectionTree=NULL;
        }
    }

    if(selectionTree==NULL){
        this->nEvents=nEvents;
        if(verbosity)cout<<"selectionTree does not exists, close file"<<endl;
        delete selectionFile;
        if(verbosity)cout<<"."<<endl;
        selectionFile=new TFile(settings->getSelectionTreeFilePath().c_str(),"RECREATE");
        selectionFile->cd();
        if(verbosity)cout<<"."<<endl;
        this->selectionTree=new TTree("selectionTree",treeDescription.str().c_str());
        if(verbosity)cout<<"."<<endl;
        createdNewTree=true;
        cout<<"\n***************************************************************\n";
        cout<<"there exists no tree:\'selectionTree\"\tcreate new one."<<selectionTree<<"\n";
        cout<<"***************************************************************\n"<<endl;
    }

    return createdNewTree;
}



void TSelectionClass::resetVariables(){

    isDetMasked = false;//one of the Silicon Planes contains a Cluster with a masked channel
    nDiamondClusters=0;
    oneAndOnlyOneSiliconCluster=true;; //One and only one cluster in each silicon plane;
    useForAnalysis=false;
    useForAlignment=false;
    useForSiliconAlignment=false;
    atLeastOneValidDiamondCluster=false;
    oneAndOnlyOneDiamondCluster = false;
    hasBigDiamondCluster=false;
    fiducialRegion =-1;
}


/**
 * Checks if there is One and Only One Cluster in each Silicon Detector
 * This includes a check if a channel is Masked or Saturated
 */
bool TSelectionClass::isOneAndOnlyOneClusterSiliconEvent(){
    bool oneAndOnlyOneClusterInAllSilicon = true;
    for(UInt_t det=0;det<TPlaneProperties::getNSiliconDetectors()&&oneAndOnlyOneClusterInAllSilicon==true;det++){
        bool oneAndOnlyOne = (eventReader->getNClusters(det)==1);
        if(verbosity>10)cout<<"DET "<<det<<": "<<oneAndOnlyOne<<" "<<checkDetMasked(det)<<" "<<isSaturated(det)<<flush;
        oneAndOnlyOne = oneAndOnlyOne && !checkDetMasked(det) && !isSaturated(det);
        oneAndOnlyOneClusterInAllSilicon=oneAndOnlyOneClusterInAllSilicon&&oneAndOnlyOne;
    }
    return oneAndOnlyOneClusterInAllSilicon;
}


void TSelectionClass::checkSiliconTrackInFiducialCut(){
    //if not one and only one cluster one cannot check the Silicon Track if is fullfills the fiducia cut
    if (!oneAndOnlyOneSiliconCluster){
        IsInFiducialCut=false;
        fiducialValueY=N_INVALID;
        fiducialValueX=N_INVALID;
        return;
    }
    IsInFiducialCut=true;
    //Calculate Fiducial Values
    fiducialValueX=0;
    fiducialValueY=0;
    for(UInt_t plane=0;plane<4;plane++){
        fiducialValueX+=eventReader->getCluster(plane,TPlaneProperties::X_COR,0).getPosition(settings->doCommonModeNoiseCorrection());
        fiducialValueY+=eventReader->getCluster(plane,TPlaneProperties::Y_COR,0).getPosition(settings->doCommonModeNoiseCorrection());
    }
    fiducialValueX/=4.;
    fiducialValueY/=4.;
    //check the fiducial Values
    IsInFiducialCut = fiducialCuts->IsInFiducialCut(fiducialValueX,fiducialValueY);
    if(IsInFiducialCut){
        fiducialRegion = fiducialCuts->getFiducialCutIndex(fiducialValueX,fiducialValueY);
        if(verbosity>6){
            cout<< setw(6) << nEvent <<" fidCut: "
                    << std::right<<setw(6)<<std::setprecision(2)<<fiducialValueX<<" / "
                    << std::left <<setw(6)<<std::setprecision(2)<<fiducialValueY
                    <<" - "<<std::right<<IsInFiducialCut<< " -> "<<fiducialRegion<<"   "<<flush;
            fiducialCuts->getFidCut(fiducialRegion)->Print();
        }
    }
    //	if(verbosity>4)cout<<"fidCut:"<<fiducialValueX<<"/"<<fiducialValueY<<": Fidcut:"<<IsInFiducialCut<<endl;
}


void TSelectionClass::checkSiliconTrack(){

    oneAndOnlyOneSiliconCluster = isOneAndOnlyOneClusterSiliconEvent();
    checkSiliconTrackInFiducialCut();
    if(verbosity>3)
        cout<<"\t"<<fiducialValueX<<"/"<<fiducialValueY<<"\tSilicon: oneAndOnlyOne:"<<oneAndOnlyOneSiliconCluster<<"\tinFidCut:"<<IsInFiducialCut<<endl;
    //if IsInFiducialCut it has one and Only One CLuster in each silicon detector
    isValidSiliconTrack = IsInFiducialCut;
    isSiliconTrackNotFiducialCut = !IsInFiducialCut&&oneAndOnlyOneSiliconCluster;
}

void TSelectionClass::checkDiamondTrack(){
    //Do not look at tracks where is not at least one and only one cluster in each sil det
    if(!oneAndOnlyOneSiliconCluster)
        return;
    UInt_t diaDet =TPlaneProperties::getDetDiamond();

    nDiamondClusters=eventReader->getNClusters(diaDet);

    isDiaSaturated=this->isSaturated(diaDet);
    if(verbosity>4&&nDiamondClusters>0&&!IsInFiducialCut)
        cout<<"\nThis event has diamond hit which is not in fid Cut"<<endl;

    atLeastOneValidDiamondCluster = nDiamondClusters >0 && !checkDetMasked(diaDet) && !isDiaSaturated;
    oneAndOnlyOneDiamondCluster = atLeastOneValidDiamondCluster && nDiamondClusters==1;
    hasBigDiamondCluster =false;
    for(UInt_t cl=0;cl<nDiamondClusters;cl++){
        nDiaClusterSize = eventReader->getClusterSize(diaDet,cl);
        if(nDiaClusterSize>=3){
            hasBigDiamondCluster=true;
            oneAndOnlyOneDiamondCluster=false;
            atLeastOneValidDiamondCluster=false;
        }
    }
    if(verbosity>3){
        cout<<"\tDiamond: nClusters:"<<nDiamondClusters<<"\tSaturated:"<<isDiaSaturated<<"\tmasked:"<<checkDetMasked(TPlaneProperties::getDetDiamond());
        cout<<"\tatLeastOne:"<<atLeastOneValidDiamondCluster<<"\texactlyOne"<<oneAndOnlyOneDiamondCluster;
        cout<<"\thasBigCluster: "<<hasBigDiamondCluster<<endl;
    }

}
/**
 * In this function all
 */
void TSelectionClass::setVariables(){
    if(verbosity>3)cout<<"\nsetVariables: "<<nEvent<<endl;

    checkSiliconTrack();
    checkDiamondTrack();


    useForSiliconAlignment = isValidSiliconTrack&& !oneAndOnlyOneDiamondCluster&&IsInFiducialCut;//isValidDiamondEvent;// one and only one hit in silicon but not exactly one hit in diamond
    //	useForAlignment = oneAndOnlyOneDiamondCluster&&settings->useForAlignment(nEvent,nEvents)&&IsInFiducialCut;//one and only one hit in all detectors (also diamond)
    //	useForAnalysis=oneAndOnlyOneDiamondCluster&&!useForAlignment&&IsInFiducialCut;;
    useForAlignment = atLeastOneValidDiamondCluster&&settings->useForAlignment(nEvent,nEvents)&&IsInFiducialCut;//one and only one hit in all detectors (also diamond)
    useForAnalysis=atLeastOneValidDiamondCluster&&!useForAlignment&&IsInFiducialCut;;
    validMoreThanOneClusterDiamondevent = atLeastOneValidDiamondCluster && !oneAndOnlyOneDiamondCluster&&IsInFiducialCut;
    doEventCounting();
    fillHitOccupancyPlots();
}

void TSelectionClass::fillHitOccupancyPlots(){
    if(!oneAndOnlyOneSiliconCluster)
        return;
    hFiducialCutSilicon->Fill(fiducialValueX,fiducialValueY);
    //	if((isValidSiliconTrack||isSiliconTrackNotFiducialCut)&&nDiamondClusters>0)
    if (settings->isInRoughFiducialCut(fiducialValueX,fiducialValueY))
        hFiducialCutSiliconRoughCut->Fill(fiducialValueX,fiducialValueY);
    FillHitOccupancyPlotsSamePattern();
    if(!atLeastOneValidDiamondCluster)
        return;
    hFiducialCutSiliconDiamondHit->Fill(fiducialValueX,fiducialValueY);
    if(!oneAndOnlyOneDiamondCluster)
        return;
    hFiducialCutSiliconOneAndOnlyOneDiamondHit->Fill(fiducialValueX,fiducialValueY);
    if(!IsInFiducialCut)
        return;
    hAnalysisFraction->Fill(nEvent);
    hSelectedEvents->Fill(fiducialValueX,fiducialValueY);
    //	if (verbosity>4)
    //		  cout<<nEvent<<" selected Event for ana or alignment @ "<<setprecision(4) <<std::setw(5)<<fiducialValueX<<"/"<<setprecision (4) <<std::setw(5)<<fiducialValueY<<":\t"<<std::setw(3)<<fiducialCuts->getFidCutRegion(fiducialValueX,fiducialValueY)<<endl;;
    //	}
    //	else
    //		cout<<nEvent<<" not Use for ana or alignment @ "<<setprecision(4) <<std::setw(5)<<fiducialValueX<<"/"<<setprecision (4) <<std::setw(5)<<fiducialValueY<<":\t"<<std::setw(3)<<fiducialCuts->getFidCutRegion(fiducialValueX,fiducialValueY)<<endl;;
}
void TSelectionClass::doEventCounting(){
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
    if(isValidSiliconTrack&&!atLeastOneValidDiamondCluster)
        nValidSiliconNoDiamondHit++;
    if(isValidSiliconTrack&&oneAndOnlyOneDiamondCluster)
        nValidDiamondTrack++;
    if(hasBigDiamondCluster)
        nToBigDiamondCluster++;
    if(!isValidSiliconTrack)
        nNoValidSiliconTrack++;
    else
        nValidSiliconTrack++;
    Int_t fidCutNo = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
    while (fidCutNo >=0 && fidCutNo >= vSiliconEventsInFidCutRegions.size() ){
        vSiliconEventsInFidCutRegions.push_back(0);
        cout<<fidCutNo<<" "<<vSiliconEventsInFidCutRegions.size()<<endl;
    }
    if (fidCutNo>=0){
//        cout<<" Fill "<<fidCutNo<<" from: "<<vSiliconEventsInFidCutRegions[fidCutNo];
        vSiliconEventsInFidCutRegions[fidCutNo]++;
//        cout<<" --> "<<vSiliconEventsInFidCutRegions[fidCutNo]<<endl;
    }
}

/**
 * checks if a detector has a masked cluster
 * opens checkDetMasked(det,cl)
 * gives verbosity output if verb>5 for diamond
 * and verb>7 for silicon
 */
bool TSelectionClass::checkDetMasked(UInt_t det){
    bool isMasked=false;
    if((verbosity>5 && TPlaneProperties::isDiamondDetector(det))||verbosity>7)
        cout<<"checkDetMasked("<<det<<"):\t";
    for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
        TCluster cluster = eventReader->getCluster(det,cl);
        bool clusterMasked = settings->isMaskedCluster(det,cluster,settings->checkAdjacentChannelsMasked());
        isMasked=isMasked||clusterMasked;
        if((verbosity>5 && TPlaneProperties::isDiamondDetector(det))||verbosity>7)
            cout<<cl<<":"<<clusterMasked<<" ";
    }
    if((verbosity>5 && TPlaneProperties::isDiamondDetector(det))||verbosity>7)
        cout<<endl;
    return isMasked;
}


/**
 * @todo: probably move to  Tevent?
 */
bool TSelectionClass::checkDetMasked(UInt_t det,UInt_t cl){
    if(verbosity>20)cout<<"masked of "<<det<<"-"<<cl<<":\t"<<flush;
    bool isMasked=false;

    if(cl<eventReader->getNClusters(det)){
        //		if(verbosity>20)cout<<"getCLuster"<<flush;
        TCluster cluster = eventReader->getCluster(det,cl);
        if(verbosity>20)cout<<"."<<flush;
        UInt_t min = cluster.getSmallestChannelNumber();
        if(verbosity>20)cout<<"."<<min<<flush;
        UInt_t max = cluster.getHighestChannelNumber();
        if(verbosity>20)cout<<":"<<max<<"\t"<<flush;
        for(UInt_t ch = min; ch<=max;ch++){
            bool channelScreened = settings->getDet_channel_screen(det).isScreened(ch);
            isMasked=isMasked|| channelScreened;
            if(verbosity>20)cout<<"ch"<<ch<<":"<<channelScreened<<"=>"<<isMasked<<" "<<flush;
        }
        if(verbosity>20)cout<<endl;
    }
    else
        return true;
    return isMasked;

    return false;
}

void TSelectionClass::setBranchAdressess(){
    selectionTree->Branch("nDiamondHits",&nDiamondClusters,"nDiamondHits/i");
    selectionTree->Branch("IsInFiducialCut",&IsInFiducialCut,"IsInFiducialCut/O");
    selectionTree->Branch("isDetMasked",&isDetMasked,"isDetMasked/O");
    selectionTree->Branch("hasValidSiliconTrack",&oneAndOnlyOneSiliconCluster,"hasValidSiliconTrack/O");
    selectionTree->Branch("useForSiliconAlignment",&this->useForSiliconAlignment,"useForSiliconAlignment/O");
    selectionTree->Branch("useForAlignment",&this->useForAlignment,"useForAlignment/O");
    selectionTree->Branch("useForAnalysis",&this->useForAnalysis,"useForAnalysis/O");
    selectionTree->Branch("diaClusterSize",&this->nDiaClusterSize,"diaClusterSize/I");
    selectionTree->Branch("isDiaSaturated",&this->isDiaSaturated,"isDiaSaturated/O");
    selectionTree->Branch("fiducialRegion",&this->fiducialRegion,"fiducialRegion/I");
    selectionTree->Branch("fiducialValueX",&this->fiducialValueX,"fiducialValueX/F");
    selectionTree->Branch("fiducialValueY",&this->fiducialValueY,"fiducialValueY/F");
}

void TSelectionClass::initialiseHistos()
{
    TString name = "hFidCutSilicon_OneAndOnlyOneCluster";
    hFiducialCutSilicon = new TH2F(name,name,512,0,256,512,0,256);
    hFiducialCutSilicon->GetYaxis()->SetTitle("yCoordinate in Channels");
    hFiducialCutSilicon->GetXaxis()->SetTitle("xCoordinate in Channels");
    hFiducialCutSilicon->GetXaxis()->SetTitle("no. of entries");

    name = "hFidCutSilicon_OneAndOnlyOneCluster_RoughCut";
    hFiducialCutSiliconRoughCut = (TH2F*)hFiducialCutSilicon->Clone(name);
    hFiducialCutSiliconRoughCut->SetTitle(name);

    name = "hFidCutSilicon_OneAndOnlyOneCluster_DiamondCluster";
    hFiducialCutSiliconDiamondHit = (TH2F*)hFiducialCutSilicon->Clone(name);
    hFiducialCutSiliconDiamondHit->SetTitle(name);

    name = "hFiducialCutSiliconOneAndOnlyOneDiamondHit";
    hFiducialCutSiliconOneAndOnlyOneDiamondHit = (TH2F*)hFiducialCutSilicon->Clone(name);
    hFiducialCutSiliconOneAndOnlyOneDiamondHit->SetTitle(name);

    name = "hSelectedEvents";
    TString title  = "ValidSiliconTrackInFidCutAndOneAndOnlyOneDiamondHit";
    hSelectedEvents = (TH2F*)hFiducialCutSilicon->Clone(name);
    hSelectedEvents->SetTitle(title);

    int nEvents= eventReader->GetEntries();
    int i=nEvents/1000;
    i++;
    nEvents = (i)*1000;
    hAnalysisFraction = new TH1F("hAnalysisFraction","hAnalysisFraction",i,0,nEvents);
    hAnalysisFraction->SetTitle("Fraction of Events for Analysis");
    hAnalysisFraction->GetXaxis()->SetTitle("event no");
    hAnalysisFraction->GetYaxis()->SetTitle("fraction of daimond + silicon hit events (%)");

    UInt_t area =  settings->getSelectionFidCuts()->getNFidCuts();
    area = TMath::Max(area,settings->diamondPattern.getNPatterns());
    for(UInt_t i =0;i <area; i++){
        name = TString::Format("hFiducialCutSiliconDiamondHitSamePattern_%d",i+1);
        TH2F* histo= (TH2F*)hFiducialCutSilicon->Clone(name);
        histo->SetTitle(name);
        mapFiducialCutSiliconDiamondHitSamePattern[i+1] =histo;
    }

    name = "hDiamondPatternFiducialPattern";
    hDiamondPatternFiducialPattern = new TH1F(name,name,settings->diamondPattern.getNPatterns()+1,0.5,settings->diamondPattern.getNPatterns()+1.5);
    hDiamondPatternFiducialPattern->GetXaxis()->SetTitle("pattern");
    hDiamondPatternFiducialPattern->GetYaxis()->SetTitle("number of correct patterns");

    name = "hDiamondPatternFiducialPatternNoMapping";
    hDiamondPatternFiducialPatternNoMapping = new TH1F(name,name,settings->diamondPattern.getNPatterns()+1,0.5,settings->diamondPattern.getNPatterns()+1.5);
    hDiamondPatternFiducialPatternNoMapping->GetXaxis()->SetTitle("pattern selectionCut");
    hDiamondPatternFiducialPatternNoMapping->GetYaxis()->SetTitle("number of no mapping found");

    name = "pDiamondPatternFiducialPatternProfile";
    pDiamondPatternFiducialPatternProfile = new TProfile(name,name,settings->diamondPattern.getNPatterns()+1,05,settings->diamondPattern.getNPatterns()+1.5);
    pDiamondPatternFiducialPatternProfile->GetXaxis()->SetTitle("pattern selectionCut");
    pDiamondPatternFiducialPatternProfile->GetYaxis()->SetTitle("rel. number of mappings found");



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
    TString name = hFiducialCutSilicon->GetName();
    name.Insert(0,"c");
    TCanvas *c1= fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSilicon);
    TH1D* hProjection = hFiducialCutSilicon->ProjectionX("hFidCutSilicon_OneAndOnlyOneCluster_px");
    cout<<"\nSaving: "<<hProjection->GetName()<<" "<<hProjection->GetEntries()<<endl;
    histSaver->SaveHistogram(hProjection);
    delete hProjection;
    hProjection = hFiducialCutSilicon->ProjectionY("hFidCutSilicon_OneAndOnlyOneCluster_py");
    cout<<"Saving: "<<hProjection->GetName()<<" "<<hProjection->GetEntries()<<endl;
    histSaver->SaveHistogram(hProjection);
    delete hProjection;
    c1->SetName(name);
    histSaver->SaveCanvas(c1);
    delete c1;
    c1 = 0;
    delete hFiducialCutSilicon;



    name = hFiducialCutSiliconRoughCut->GetName();
    name.Insert(0,"c");
    c1 = fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSiliconRoughCut,true);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);
    hProjection = hFiducialCutSiliconRoughCut->ProjectionX("hFidCutSilicon_OneAndOnlyOneCluster_RoughCut_px");
    hProjection->SetTitle("One and onlu one silicon in Rough FidCut - Projection X");
    cout<<"\nSaving: "<<hProjection->GetName()<<" "<<hProjection->GetEntries()<<endl;
    histSaver->SaveHistogram(hProjection);
    delete hProjection;
    hProjection = hFiducialCutSiliconRoughCut->ProjectionY("hFidCutSilicon_OneAndOnlyOneCluster_RoughCut_py");
    hProjection->SetTitle("One and onlu one silicon in Rough FidCut - Projection Y");
    cout<<"Saving: "<<hProjection->GetName()<<" "<<hProjection->GetEntries()<<endl;
    histSaver->SaveHistogram(hProjection);
    delete hProjection;
    delete c1;
    c1=0;
    delete hFiducialCutSiliconRoughCut;

    name = hFiducialCutSiliconDiamondHit->GetName();
    name.Insert(0,"c");
    c1 = fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSiliconDiamondHit,true);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);
    delete c1;
    c1=0;
    delete hFiducialCutSiliconDiamondHit;

    name = "c";
    name.Append(hFiducialCutSiliconOneAndOnlyOneDiamondHit->GetName());
    c1= fiducialCuts->getAllFiducialCutsCanvas(hFiducialCutSiliconOneAndOnlyOneDiamondHit,true);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);
    delete c1;
    c1=0;
    delete hFiducialCutSiliconOneAndOnlyOneDiamondHit;

    name = "c";
    name.Append(hSelectedEvents->GetName());
    c1 = fiducialCuts->getAllFiducialCutsCanvas(hSelectedEvents,true);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);
    delete c1;
    c1=0;
    delete hSelectedEvents;

    map<Int_t,TH2F*>::iterator it;
    for (it = mapFiducialCutSiliconDiamondHitSamePattern.begin(); it!=mapFiducialCutSiliconDiamondHitSamePattern.end(); it++){
        TH2F* histo = (*it).second;
        name = histo->GetName();
        name.Replace(0,1,"c");
        c1 = fiducialCuts->getAllFiducialCutsCanvas(histo,true);
        c1->SetName(name);
        histSaver->SaveCanvas(c1);
        delete c1;
        c1=0;
        delete histo;
    }
    hAnalysisFraction->Scale(.1);
    hAnalysisFraction->SetStats(false);
    //	hAnalysisFraction->GetYaxis()->SetRangeUser(0,100);
    histSaver->SaveHistogram(hAnalysisFraction);
    delete hAnalysisFraction;


    histSaver->SaveHistogram(hDiamondPatternFiducialPattern);

    histSaver->SaveHistogram(hDiamondPatternFiducialPatternNoMapping);
    name = "stackPatternMapping";
    hDiamondPatternFiducialPattern->SetLineColor(kGreen);
    hDiamondPatternFiducialPatternNoMapping->SetLineColor(kRed);
    THStack *stack = new THStack(name,name);
    stack->Add(hDiamondPatternFiducialPattern);
    stack->Add(hDiamondPatternFiducialPatternNoMapping);
    stack->Draw();
    if(stack->GetXaxis()) stack->GetXaxis()->SetTitle("pattern no.");
    if(stack->GetYaxis()) stack->GetYaxis()->SetTitle("number of entries #");
    histSaver->SaveStack(stack,"hist",true);
    if(stack) delete stack;
    histSaver->SaveHistogram(pDiamondPatternFiducialPatternProfile);
    if(hDiamondPatternFiducialPatternNoMapping) delete hDiamondPatternFiducialPatternNoMapping;
    if(hDiamondPatternFiducialPattern) delete hDiamondPatternFiducialPattern;
    if (pDiamondPatternFiducialPatternProfile) delete pDiamondPatternFiducialPatternProfile;

    name = "hSiliconEventsInCutRegions";
    Int_t bins = vSiliconEventsInFidCutRegions.size()+1;
    TH1F *htemp = new TH1F(name,name,bins,0,bins);
    for (UInt_t i = 0; i < vSiliconEventsInFidCutRegions.size(); i++){
        cout<<"Set Region: "<<i<<": "<<vSiliconEventsInFidCutRegions.at(i)<<endl;
        htemp->Fill(i,vSiliconEventsInFidCutRegions[i]);
        htemp->SetEntries(htemp->GetEntries()+vSiliconEventsInFidCutRegions[i]);
    }
    htemp->GetXaxis()->SetTitle("SelectionCutRegion");
    htemp->GetYaxis()->SetTitle("no of entries");
    histSaver->SaveHistogram(htemp,false,false,true,"histTEXT");
    name.Append("_Normalized");
    TH1F* htemp2 = (TH1F*)htemp->DrawNormalized(name,100);
    htemp2->SetName(name);
    htemp2->SetTitle(name);
    htemp2->GetYaxis()->SetTitle("rel. no of entries %");
    histSaver->SaveHistogram(htemp2,false,false,true,"histTEXT");
    delete htemp;
    delete htemp2;
}

void TSelectionClass::FillHitOccupancyPlotsSamePattern() {
    UInt_t diamondDet = TPlaneProperties::getDetDiamond();
    Int_t selectionPattern = settings->getSelectionFidCuts()->getFidCutRegion(fiducialValueX,fiducialValueY);
    if(verbosity>4) cout<<"\n";
    if(verbosity>4) settings->diamondPattern.showPatterns();
    bool mappingFound = false;
    for(UInt_t cl = 0; cl< nDiamondClusters;cl++){
        if(verbosity>4)cout<<"Cluster: "<<cl+1<<"/"<<nDiamondClusters<<endl;
        TCluster cluster = eventReader->getCluster(diamondDet,cl);
        if(verbosity>4)cluster.Print();
        Int_t diamondPattern = settings->diamondPattern.getClusterPattern(&cluster);
        if(verbosity>4)cout<<"selectionPattern: "<<selectionPattern<<"\tdiamondPattern: "<<diamondPattern<<endl;
        if(selectionPattern == diamondPattern && diamondPattern!=-1){
            mappingFound = true;
            mapFiducialCutSiliconDiamondHitSamePattern[diamondPattern]->Fill(fiducialValueX,fiducialValueY);
            hDiamondPatternFiducialPattern->Fill(selectionPattern);
            pDiamondPatternFiducialPatternProfile->Fill(selectionPattern,1);
        }
    }
    if(!mappingFound&& selectionPattern>0){
        pDiamondPatternFiducialPatternProfile->Fill(selectionPattern,0);
        hDiamondPatternFiducialPatternNoMapping->Fill(selectionPattern);
    }
}

void TSelectionClass::createFiducialCut(){

    //	std::vector<std::pair<Float_t,Float_t> > xInt,yInt;
    //	xInt.push_back( make_pair(settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh()));
    //	yInt.push_back( make_pair(settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
    //	fiducialCuts = new TFidCutRegions(xInt,yInt,1);
    fiducialCuts = settings->getSelectionFidCuts();
    cout<<"Create AutoFidCut with "<<endl;
    fiducialCuts->Print(1);
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

    if(settings->getUseAutoFidCut()==true){
        delete fiducialCuts;
        fiducialCuts = new TFidCutRegions(hFiducialCutSiliconDiamondHit,settings->getNDiamonds(),settings->getAutoFidCutPercentage());
        fiducialCuts->setRunDescription(settings->getRunDescription());
    }
    else{
        fiducialCuts->setHistogramm(hFiducialCutSiliconDiamondHit);
    }
    histSaver->SaveCanvas(fiducialCuts->getFiducialCutCanvas(TPlaneProperties::X_COR));
    histSaver->SaveCanvas(fiducialCuts->getFiducialCutCanvas(TPlaneProperties::Y_COR));
    TCanvas *c1 = fiducialCuts->getFiducialCutCanvas(TPlaneProperties::XY_COR);
    c1->SetTitle(Form("Fiducial Cut of Run %i with \"%s\"",settings->getRunNumber(),settings->getRunDescription().c_str()));
    c1->SetName("cFidCutCanvasXY");
    histSaver->SaveCanvas(c1);
    if(verbosity)
        fiducialCuts->Print(1);
    if (verbosity>3&&verbosity%2==1){
        cout<<"Press a key and enter to continue..."<<flush;
        char t;
        cin >>t;
    }
}
