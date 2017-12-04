/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(TSettings* inputSettings,TSettings::alignmentMode mode) {
    if (inputSettings) verbosity = inputSettings->getVerbosity();
    cout << "\n\n\n**********************************************************" << endl;
    cout << "*************TAlignment::TAlignment***********************" << endl;
    cout << "**********************************************************" << endl;
    this->mode = mode;
    sys = gSystem;
    setSettings(inputSettings);
    runNumber = settings->getRunNumber();
    if(verbosity) cout << runNumber << endl;
    stringstream runString;
    settings->goToSelectionTreeDir();
    htmlAlign = new THTMLAlignment(settings);
    if (mode == TSettings::transparentMode){
        settings->goToAlignmentRootDir();
        eventReader = new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);
    }
    else{
//TTracking* trackReader = new TTracking(settings->getSelectionTreeFilePath(),settings->getAlignmentFilePath(),settings->getEtaDistributionPath(),settings);

        eventReader = (TTracking*) new TADCEventReader(settings->getSelectionTreeFilePath(),settings);
//	settings->getRunNumber(),settings->getVerbosity()<15?0:settings->getVerbosity()-15);
        eventReader->setEtaDistributionPath(settings->getEtaDistributionPath());
//		cout<<"Eta dist path: "<<eventReader->getEtaDistributionPath()<<endl;
    }
    histSaver = new HistogrammSaver(inputSettings);
    settings->goToAlignmentDir(mode);
    histSaver->SetPlotsPath(settings->getAlignmentDir(mode));
    histSaver->SetRunNumber(runNumber);
    histSaver->SetNumberOfEvents(eventReader->GetEntries());
    htmlAlign->setFileGeneratingPath(settings->getAlignmentDir(mode));
    settings->goToAlignmentRootDir(mode);
    if(verbosity) cout << "end initialise" << endl;
    alignmentPercentage = settings->getAlignment_training_track_fraction();
    Float_t stripSize = 1.;    // 50./10000.;//mu m
    detectorD0Z = 0.725 / stripSize;    // by definition in cm
    detectorD1Z = 1.625 / stripSize;    // by definition in cm
    detectorD2Z = 18.725 / stripSize;    // by definition in cm
    detectorD3Z = 19.625 / stripSize;    // by definition in cm
    detectorDiaZ = 10.2 / stripSize;    // by definition in cm
//	cout<<"Verbosity is: "<<verbosity<<" "<<settings->getVerbosity()<<endl;
    res_keep_factor = settings->getRes_keep_factor();
    if(verbosity) cout << "Res Keep factor is set to " << res_keep_factor << endl;
    align = NULL;
    myTrack = NULL;
    nAlignmentStep = -1;
    nAlignSteps = settings->GetSiliconAlignmentSteps();
    nDiaAlignmentStep = -1;
    nDiaAlignSteps = settings->GetDiamondAlignmentSteps();

    if(mode == TSettings::transparentMode)
        diaCalcMode = TCluster::corEta;
    else
        diaCalcMode = TCluster::maxValue;    //todo
    silCalcMode = TCluster::corEta;    //todo

    bPlotAll = settings->doAllAlignmentPlots()||verbosity>6;
    results=0;
    gausFitValuesX.resize(4);
    gausFitValuesY.resize(4);
    cout<<"Selection Fidcuts: "<<endl;
    settings->getSelectionFidCuts()->Print();
    cout<<"\nAlignment FIdcuts: "<<settings->getAlignmentFidCuts()<<endl;
    silicon_chi2_cut = false;

}

TAlignment::~TAlignment() {
    while (events.size()) {
//		events.back().Delete();
        events.pop_back();
    }
    while (telescopeAlignmentEvent.size())
        telescopeAlignmentEvent.pop_back();
    if(verbosity)cout << "TAlignment deconstructor" << endl;
    if (results!=0)results->setAlignment(this->align);
    htmlAlign->setAlignment(align);
    htmlAlign->createContent();
    htmlAlign->generateHTMLFile();
    this->saveAlignment(this->mode);
    if(htmlAlign!=0)delete htmlAlign;
    if (myTrack) delete myTrack;
    if (histSaver) delete histSaver;
//	if(eventReader)delete eventReader;
    settings->goToOutputDir();
    if(verbosity)cout<<"Closed TAlginment"<<endl;
}

void TAlignment::setSettings(TSettings* inputsettings) {
    if(inputsettings==0)
    {
        cerr<<"Input settings == 0. BREAK!"<<endl;
        exit(-1);
    }
    this->settings = inputsettings;
}

/**
 * initialise the variable align and set the Z Offsets
 */
void TAlignment::initialiseDetectorAlignment(TSettings::alignmentMode mode) {
    if (align == NULL && mode == TSettings::transparentMode){
        cout<<"Loading transparent analysis mode..."<<flush;
        loadDetectorAlignment(mode);
    }
    else if (align ==NULL) {
        align = new TDetectorAlignment();
        align->setVerbosity(verbosity-4>0?verbosity-4:0);
        htmlAlign->setAlignment(align);
        cout << "TAlignment::Align::Detector alignment did not exist, so created new DetectorAlignment" << endl;
        align->SetZOffset(0, detectorD0Z);
        align->SetZOffset(1, detectorD1Z);
        align->SetZOffset(2, detectorD2Z);
        align->SetZOffset(3, detectorD3Z);
        align->SetZOffset(4, detectorDiaZ);
    }
}

void TAlignment::loadDetectorAlignment(TSettings::alignmentMode mode){
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
//				if(verbosity)align->Print();
                alignmentFile->Close();
//				align = (TDetectorAlignment*)align->Clone();
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
    if(mode == TSettings::transparentMode){
        cout<<" Loaded Alignment for transparent Alignment. "<<endl;
        align->Print();
    }

}


void TAlignment::createTransparentEventVectors(UInt_t nEvents, UInt_t startEvent) {
    initialiseDetectorAlignment(TSettings::transparentMode);
    cout<<"nEvents: "<<nEvents<<", startEvent: "<<startEvent<<endl;
    if (nEvents == 0) nEvents = eventReader->GetEntries() - startEvent;
    if (nEvents + startEvent > eventReader->GetEntries()) nEvents = eventReader->GetEntries() - startEvent;
    int noHitDet = 0;
//	int falseClusterSizeDet=0;
//	int noHitDia=0;
//	int falseClusterSizeDia = 0;
    int nCandidates = 0;
//	int nScreened = 0;
    int nNotInFidCut=0;
    int nNotGoodEta = 0;
    int nAnalysedEvents = 0;
    Float_t minEtaDif = settings->getMinimalAbsoluteEtaValue();
    cout << "CREATING VECTOR OF VALID EVENTS TRANSPARENT..."<<minEtaDif << endl;
    for (nEvent = startEvent; nEvent < nEvents + startEvent; nEvent++) {
        while(telescopeAlignmentEvent.size()<events.size()){ //do not use events which fullfill all critera from the alignment of the DUT
            telescopeAlignmentEvent.push_back(1);
        }
        TRawEventSaver::showStatusBar(nEvent - startEvent, nEvents, 1000);
        if(!settings->useForAlignment(nEvent,nEvents))
            break;
        eventReader->LoadEvent(nEvent);
        nAnalysedEvents++;
        if (!eventReader->isValidTrack()) {
            noHitDet++;
            continue;
        }
        Float_t fidCutX = eventReader->getFiducialValueX();
        Float_t fidCutY = eventReader->getFiducialValueY();
        if(!settings->isInAlignmentFiducialRegion(fidCutX,fidCutY)){
            nNotInFidCut++;
            continue;
        }

//			float fiducialValueX=eventReader->getFiducialValueX();
//			float fiducialValueY=eventReader->getFiducialValueY();
        Int_t subjectPlane = TPlaneProperties::getDiamondPlane();
        UInt_t subjectDetector = TPlaneProperties::getDetDiamond();
        vector< UInt_t > refPlanes =  TPlaneProperties::getSiliconPlaneVector();
        TPositionPrediction* predictedPosition = eventReader->predictPosition(subjectPlane,refPlanes,false);
//			cout<<"\n"<<endl;
//			predictedPosition->Print(2);
        Float_t predX = predictedPosition->getPositionX();
        Float_t predY = predictedPosition->getPositionY();
        Float_t metricPosInDetSystem = eventReader->getPositionInDetSystem(subjectDetector,predX,predY);
        Float_t centerPos = settings->convertMetricToChannelSpace(subjectDetector,metricPosInDetSystem);
        if(centerPos <  0)
            cout<<nEvent<<" "<<predX<<"/"<<predY<<": "<<metricPosInDetSystem<<" -->"<<centerPos<<endl;
//			cout<<nEvent<<" "<<predX<<"/"<<predY<<"\t"<<metricPosInDetSystem<< " - "<< centerPos<<endl;
        TCluster cluster = TTransparentAnalysis::makeTransparentCluster(eventReader,settings,subjectDetector,centerPos,10);
//			cluster.Print(1);
        TEvent *event = eventReader->getEvent();
        TEvent* clonedEvent = (TEvent*)event->Clone();
        TPlane* plane = clonedEvent->getPlanePointer(subjectPlane);
//			cout<<"\n\n***\n";
//			event->Print(1);
        vector<TCluster> vecClus(1,cluster);
        plane->SetXClusters(vecClus);
//			cout<<"\n";
//			clonedEvent->Print(1);
        for(UInt_t det=0;det<TPlaneProperties::getNSiliconDetectors();det++){
            Float_t clusPos = eventReader->getCluster(det,0).getPosition(settings->doCommonModeNoiseCorrection());
            if(clusPos<0||clusPos>=TPlaneProperties::getNChannels(det)||eventReader->getNClusters(det)!=1){
                cout<<"Do not take event clusPos is not valid...."<<endl;
                continue;
            }
        }
        Float_t eta = cluster.getEta(true);
        if (eta<minEtaDif || eta > 1- minEtaDif){
            nNotGoodEta++;
            if(clonedEvent) delete clonedEvent;
            continue;
        }
        nCandidates++;

        this->events.push_back(*clonedEvent);
        this->fiducialValueX.push_back(fidCutX);
        this->fiducialValueY.push_back(fidCutY);
        telescopeAlignmentEvent.push_back(0);

//			if (eventReader->useForAnalysis())
//				break;
    }

    cout << "nAnalysedEvents " << setw(7) << nAnalysedEvents << endl;
    cout << " -     noHitDet " << setw(7) << noHitDet << endl;
    cout << " - nNotInFidCut " << setw(7) << nNotInFidCut << endl;
    cout << " -  nNotGoodEta " << setw(7) << nNotGoodEta << endl;
    cout << " =  nCandidates " << setw(7) << nCandidates << endl;
}


/**
 * Fills a vector of TEvents which are candidates to be used in the Alignment
 * @param nEvents
 * @param startEvent
 */
void TAlignment::createEventVectors(UInt_t nEvents, UInt_t startEvent,enumDetectorsToAlign detAlign) {
    initialiseDetectorAlignment();
    if (nEvents == 0) nEvents = eventReader->GetEntries() - startEvent;
    if (nEvents + startEvent > eventReader->GetEntries()) nEvents = eventReader->GetEntries() - startEvent;
    int noHitDet = 0;
    UInt_t invalidClusterSize = 0;
//	int falseClusterSizeDet=0;
//	int noHitDia=0;
    int falseClusterSizeDia = 0;
    int nCandidates = 0;
    int nScreened = 0;
    int nNotInFidCut=0;
    cout << "CREATING VECTOR OF VALID EVENTS..." << endl;
    UInt_t nMaxSiliconEvents = 10000;
    UInt_t nTelescopeAlignmentEvents = 0;
    cout<<"DO Alignment: "<<detAlign<<endl;
    for (nEvent = startEvent; nEvent < nEvents + startEvent; nEvent++) {
        while(telescopeAlignmentEvent.size()<events.size()){ //do not use events which fullfill all critera from the alignment of the DUT
            telescopeAlignmentEvent.push_back(1);
            nTelescopeAlignmentEvents++;
        }
        if(events.size()%(int)100==0||(nEvent - startEvent)%(int)1000==0)
            cout<<TString::Format("\r%6d %6d",(int)events.size(),(int)(events.size()-nTelescopeAlignmentEvents))<<flush;
        TRawEventSaver::showStatusBar(nEvent - startEvent, nEvents, 1000);
        if(!settings->useForAlignment(nEvent,nEvents))
            return;
        eventReader->LoadEvent(nEvent);
        if (!eventReader->isValidTrack()) {
            noHitDet++;
            continue;
        }
        if (eventReader->isDetMasked()) {
            nScreened++;
            continue;
        }
        if (!eventReader->hasSmallSiliconClusterSizes(3)){
            invalidClusterSize++;
            continue;
        }


        float fidCutX=eventReader->getFiducialValueX();
        float fidCutY=eventReader->getFiducialValueY();

        if (events.size()<nMaxSiliconEvents && (detAlign == silAlignment || detAlign == bothAlignment)){
            this->fiducialValueX.push_back(fidCutX);
            this->fiducialValueY.push_back(fidCutY);
            this->events.push_back(*eventReader->getEvent());
        }
        if (eventReader->getNDiamondClusters() != 1) {
            falseClusterSizeDia++;
            continue;
        }


        if(!settings->isInAlignmentFiducialRegion(fidCutX,fidCutY)){
            if(verbosity>10)cout<<nEvent<<"\tevent not in correct fiducial region "<<fidCutX<<"/"<<fidCutY<<"-->"<<settings->getSelectionFidCuts()->getFiducialCutIndex(fidCutX,fidCutY)<<endl;
            nNotInFidCut++;
            continue;
        }
        else
            if(verbosity>8)cout<<nEvent<<"\tevent in alignment fiducial region "<<fidCutX<<"/"<<fidCutY<<"-->"<<settings->getSelectionFidCuts()->getFiducialCutIndex(fidCutX,fidCutY)<<endl;

        if(nEvent==startEvent&&verbosity>4)
            cout<<"\nEvent\tvalid\tnClus\tmasked\tFidCut\tAlign"<<endl;
        if(verbosity>20)
            cout<<nEvent<<"\t"<<eventReader->isValidTrack()<<"\t"<<eventReader->getNDiamondClusters()
            <<"\t"<<eventReader->isDetMasked()<<"\t"<<eventReader->IsInFiducialCut()<<"\t"<<eventReader->useForAlignment()<<endl;
        if (eventReader->useForAlignment()) {
            bool bBreak = false;
            for(UInt_t det=0;det<TPlaneProperties::getNDetectors()&&!bBreak;det++){
                Float_t clusPos = eventReader->getCluster(det,0).getPosition(settings->doCommonModeNoiseCorrection());
                if(clusPos<0||clusPos>=TPlaneProperties::getNChannels(det)||eventReader->getNClusters(det)!=1){
                    bBreak = true;
                    cout<<"Do not take event clusPos is not valid...."<<endl;
                }
            }
            if(bBreak){
                continue;
            }
            nCandidates++;
            if (events.size()>=nMaxSiliconEvents|| ! (detAlign == silAlignment || detAlign == bothAlignment)){
                this->fiducialValueX.push_back(fidCutX);
                this->fiducialValueY.push_back(fidCutY);
                this->events.push_back(*eventReader->getEvent());
            }
            telescopeAlignmentEvent.push_back(0);
        }
//		if (eventReader->useForAnalysis()) {
//			cout << "\nFound first Event for Analysis ->BREAK" << endl;
//			break;
//		}
    }
    if(verbosity||true){
        cout<<"\nCreated Vector of Events with one and only one Hit in Silicon  + one diamond Cluster + in Fiducial Cut Area"<<endl;
        cout<<"first Event: "<<startEvent<<"\t last Event: "<<nEvent<<"\t total Events: "<<endl;//todo
        cout<<"Cut Flow:\n";
        cout<<"\tTotal Events looked at:            "<<setw(7)<<nEvent-startEvent<<"\n";
        cout<<"\tEvents with no Silicon Hit:        "<<setw(7)<<noHitDet<<"\n";
        cout<<"\tInvalid Silicon Clustersize:       "<<setw(7)<<invalidClusterSize<<"\n";
        cout<<"\tSil Track not in Fid Cut:          "<<setw(7)<<nNotInFidCut<<"\n";
        cout<<"\tNo of Diamond Clust. != 1:         "<<setw(7)<<falseClusterSizeDia<<"\n";
        cout<<"\t                                   "<<"-------\n";
        cout<<"\tTotal Events for Dia Alignment:    "<<setw(7)<<nCandidates<<" = "<<(Float_t)nCandidates/(Float_t)(nEvent-startEvent)*100.<<"%\n"<<endl;
        cout<<"\tAdd. Events for Silicon Alignment: "<<setw(7)<<nTelescopeAlignmentEvents<<"\n";
    }
    align->addEventIntervall(startEvent, nEvent);
    align->setNUsedEvents((UInt_t) events.size());
    align->setAlignmentTrainingTrackFraction(settings->getAlignment_training_track_fraction());
    cout<<"AlignmentTrainingTrackFraction: "<< align->getAlignmentTrainingTrackFraction() <<endl;
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
    cout << "Alignment of ";
    switch (detToAlign){
        case diaAlignment:
            cout <<" Diamond."<<endl;break;
        case silAlignment:
            cout<<" Silicon."<<endl;break;
        case bothAlignment:
            cout<<" Silicon and Diamond."<<endl;break;
    }

    if (verbosity>2) {
        cout << "\n\n\nTAlignment::Align:Starting \"" << histSaver->GetPlotsPath() << "\"" << endl;
        cout << "\t\t" << events.size() << "\t";
        cout << "\t\t " << eventReader << " ." << endl;
    }

    if(this->mode == TSettings::transparentMode){
        cout<< "\n***Transparent Alignment mode***" << endl;
        cout<< " Det to align:: "<<detToAlign<<endl;
        cout<< " ResetAlignmen: "<<settings->resetAlignment()<<endl;
        cout<< " AlignmentMode of diamond: "<<diaCalcMode<<endl;
    }
    if(detToAlign!=diaAlignment&&settings->resetAlignment())
        initialiseDetectorAlignment();
    else
        loadDetectorAlignment();
    vector<Float_t> test1,test2;

    //create an TTrack object and set the eta distributions.
    if (myTrack == NULL) {
        if(verbosity>2)cout << "TAlignment::Align::create new TTrack" << endl;
        myTrack = new TTrack(align,settings);
        if(verbosity>2)cout << "TAlignment::Align::created new TTrack" << endl;
        for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++){
            TH1F* etaInt = eventReader->getEtaIntegral(det);
            if(etaInt==0){char t;cout<<"eta Int ==0"<<det<<"\tPress a key and enter to confirm"<<endl;cin >>t;}
            myTrack->setEtaIntegral(det,etaInt );
        }

        myTrack->setVerbosity(verbosity-4>0?verbosity-4:0);
    }
    if(!myTrack){
        cerr<<"could not create my Track ----> EXIT"<<endl;
        exit(-1);
    }
    if (false)
        UpdateResolutions(test1,test2);

    if (events.size() == 0) createEventVectors(nEvents, startEvent,detToAlign);

    myTrack->setDetectorAlignment(align);

    /*
     * distinguish here between diamond and silicon alignment
     */
    if((detToAlign==silAlignment||detToAlign==bothAlignment)){
        AlignSiliconPlanes();
    }
    this->saveAlignment(this->mode);
    if(detToAlign==diaAlignment||detToAlign==bothAlignment) {
        if(settings->resetAlignment())
            align->ResetAlignment(TPlaneProperties::getDiamondPlane());
        AlignDiamondPlane();
    }
//	align->PrintResults((UInt_t) 1);
    cout<<"Done with alignment of ";
    if (detToAlign==diaAlignment)cout<<"diamond"<<endl;
    else if(detToAlign==silAlignment)cout<<"silicon"<<endl;
    else if(detToAlign==bothAlignment)cout<<"silicon and diamond"<<endl;
    return (1);
}

/**
 *
 */
void TAlignment::AlignSiliconPlanes() {
    silicon_chi2_cut = false;
    cout<<"Alignment of Silicon Planes. max. Alignment Steps: "<<nAlignSteps<<endl;
    nAlignmentStep = -1;
    if(!align->isPreAligned() || bPlotAll){
        std::cout << "Detector is NOT pre-aligned" << endl;
        if (verbosity)cout << "\nCheck Detector Alignment:" << endl;
        resPlane1 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 1, 0, 3, true);
        resPlane2 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 2, 0, 3, true);
        resPlane3 = CheckDetectorAlignment(TPlaneProperties::XY_COR, 3, 1, 2, true);
    }
    else{
        std::cout << "Detector is pre-aligned" << endl;
        if(verbosity)
            cout<<"Detectors are already pre alignend do not create new Prealigment plots."<<endl;
    }


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
        //nAlignmentStep-1>=nAlignSteps||bPlotAll,
        siliconAlignmentStep(true,false);
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
        yPhiOff3 = TMath::Abs(align->GetLastPhiYOffset(3));
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

    settings->diamondPattern.hasInvalidIntervals();
    silicon_chi2_cut = true;
    //create ReferencePlane Vector: using all Silicon Planes for Alignment 0,1,2,3
    UInt_t diaPlane = 4;
    vector<UInt_t> vecRefPlanes;
    for (UInt_t i = 0; i < 4; i++)
        if (i != diaPlane) vecRefPlanes.push_back(i);
    nDiaAlignmentStep = -1;
    //checking Residual
    settings->diamondPattern.hasInvalidIntervals();
    cout <<  "check Strip detector alignment "<<endl;
    TResidual resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false, true);
    bool diaAlignmentDone = false;
    settings->diamondPattern.hasInvalidIntervals();
    for (nDiaAlignmentStep = 0; (nDiaAlignmentStep < nDiaAlignSteps) && (!diaAlignmentDone||nDiaAlignmentStep<2)	 ; nDiaAlignmentStep++) {
        cout << "\n\n ****" << nDiaAlignmentStep << " of " << nDiaAlignSteps << " Steps..." << endl;
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
    cout<<"\n\n\nCheckDIAMONDStripDetectorAlignment"<<endl;
    nDiaAlignmentStep = nDiaAlignSteps;
    resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, true, true, resDia);
    cout<<"\n\n\nCheckDIAMONDStripDetectorAlignment: CHI2"<<endl;
    CheckStripDetectorAlignmentChi2(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false, true, settings->getAlignment_chi2());
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
    if( res.getUsedTracks()/(Float_t)(vecXLabPredMetric.size()) < 0.1){
        cout<<"Something is wrong used less than 10% of the tracks for alignment: "<<res.getUsedTracks()<<" Tracks out of "<<vecXLabPredMetric.size()<<endl;
        if(verbosity>2&&verbosity%2==1){
            char t;
            cin >> t;
        }
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
    if(verbosity>7&&verbosity%2==1){
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

TString TAlignment::GetReferencePlaneString(vector<UInt_t> *vecRefPlanes){
    if (verbosity){
        cout<<"Ref Planes: [";
        for(UInt_t refPlane1=0;refPlane1<vecRefPlanes->size();refPlane1++){
            cout<<" "<<vecRefPlanes->at(refPlane1);
        }
        cout<<"]"<<endl;
    }
    for(UInt_t refPlane1=0;refPlane1<vecRefPlanes->size()-1;refPlane1++){
        for(UInt_t refPlane2=refPlane1+1;refPlane2<vecRefPlanes->size();refPlane2++){
            if(vecRefPlanes->at(refPlane1)==vecRefPlanes->at(refPlane2)){
                if (verbosity>1) cout<<" Removing "<<refPlane1<<" "<<refPlane2<<":"<< vecRefPlanes->at(refPlane1)<<" from"<<vecRefPlanes->size()<<"-";
                vecRefPlanes->erase(vecRefPlanes->begin()+refPlane1);
                if (verbosity>1) cout<<vecRefPlanes->size()<<endl;
            }
        }
    }
    TString refPlaneString;
    for (UInt_t i = 0; i < vecRefPlanes->size(); i++)
        if (i == 0)
            refPlaneString = TString::Format("%d",vecRefPlanes->at(i));
        else if (i + 1 < vecRefPlanes->size())
            refPlaneString.Append(TString::Format("_%d",vecRefPlanes->at(i)));
        else
            refPlaneString.Append(TString::Format("_and_%d",vecRefPlanes->at(i)));
    return refPlaneString;
}


TResidual TAlignment::Residual(alignmentMode aligning, TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        vector<UInt_t> vecRefPlanes,bool bAlign, bool bPlot, TResidual resOld,
        TCluster::calculationMode_t mode,resCalcMode calcMode,Float_t maxChi2){

    cout<<"\n\n[TAlignment::Residual]"<<endl;
    TString refPlaneString = GetReferencePlaneString(&vecRefPlanes);

    bool isStripAlignment = (aligning == singleStrip);
    bool verb = verbosity||isStripAlignment;
    if (verb){
        cout << "TAlignment::get";
        if(isStripAlignment)cout<<"SingleStrip";
        cout<<"Residual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor);
        cout << " with " << refPlaneString << ", plotting: " << bPlot<< ", using "<<events.size()<<" Events\t";
        if (resOld.isTestResidual()) cout<<" having a TestResidual";
        cout<<"\tClusterPos: ";
        switch (mode) {
            case TCluster::maxValue: cout<<"MaxValue";break;
            case TCluster::chargeWeighted: cout<<"Charge Weigthed";break;
            case TCluster::highest2Centroid: cout<<"Highest2Centroid";break;
            case TCluster::eta: cout<<"Eta";break;
            case TCluster::corEta: cout<<"corEta";break;
            case TCluster::highest2CentroidNoSmallHits: cout<<"highest2CentroidNoSmallHits";
            default:
                cout<<"Unknown";
        }
            cout<<"\tResidualCalcMode: ";
            switch (calcMode){
                case chi2CalcMode: cout<<"Chi2 CalcMode";break;
                case normalCalcMode: cout<<"normal CalcMode";break;
                case resolutionCalcMode: cout<<"resolution CalcMode";break;
                default:
                    cout<<"Unknown CalcMode";break;
            }
        cout<<endl;

        if(verb&&resOld.isTestResidual())resOld.Print(1);
    }
    clearMeasuredVectors();

    UInt_t nEvents = events.size();
    UInt_t nUsedEvents=0;
    UInt_t nNotUsedEvents=0;

    TPositionPrediction* predictedPostionMetric = 0;
    Float_t xPredSigma,yPredSigma;
    Float_t xLabPredictedMetric, yLabPredictedMetric;
    Float_t xDetPredictedMetric, yDetPredictedMetric;

    Float_t xDelta,yDelta, xLabMeasuredMetric,yLabMeasuredMetric, xDetMeasuredMetric, yDetMeasuredMetric, resxtest, resytest, xPhi,yPhi,chi2x,chi2y,eta;
    Int_t clusterSize;
    bool isTelescopeAlignment = TPlaneProperties::isSiliconPlane(subjectPlane) && TPlaneProperties::AreAllSiliconPlanes(vecRefPlanes);
    if (isTelescopeAlignment){
        cout<<" is Telescope Alignment with "<<nEvents<<endl;
    }
    else
        cout<<" is Diamond Alignment with "<<nEvents<<endl;
    bool cmnCorrected = settings->doCommonModeNoiseCorrection();
    Float_t pw_sil = settings->getPitchWidth(0,0);
    Float_t maxXLabMetric = TPlaneProperties::getNChannelsSilicon() * pw_sil;
    Float_t maxYLabMetric = TPlaneProperties::getNChannelsSilicon() * pw_sil;
    if (silicon_chi2_cut){
        Float_t chi2 = settings->getAlignment_chi2();
        cout<<"Applying Chi2 cut on silicon for Alignment of diamond: "<<chi2<<endl;
        if (maxChi2 > chi2)
            maxChi2 = chi2;
    }
    Int_t nNonTelescopeAlignmentEvents = 0;
    Int_t nInvalidxLabPredictedMetric = 0;
    Int_t nInvalidyLabPredictedMetric = 0;
    Int_t nInvalidXYLabMeasurement = 0;
    Int_t nCalcModeCut = 0;
    Int_t nTelescoepAlignmentEvent = 0;
    Int_t nInvalidPredictedRelPos = 0;
    Int_t nInvalidMeasuredRelPos = 0;
    Int_t nIgnoreForStripAlignment = 0;
    Int_t nInvalidSiliconChi2Cut = 0;
    vector<Float_t> calcModeCutX;
    vector<Float_t> calcModeCutY;
    for (UInt_t nEvent = 0; nEvent < nEvents; nEvent++) {
        TRawEventSaver::showStatusBar(nEvent, nEvents);
        if (!isTelescopeAlignment&&telescopeAlignmentEvent[nEvent]){
            nNonTelescopeAlignmentEvents++;
            continue;
        }
//        cout<<"get Event"<<endl;
        myTrack->setEvent(&events.at(nEvent));
        if (verbosity > 8) cout << "Event no.: " << nEvent << endl;

        if (verbosity > 5) cout<<"Predict Position: "<<endl;
//        cout<<"\tpredict pos"<<endl;
        bool bVerbose = verbosity >7||(isStripAlignment && verbosity);
        predictedPostionMetric = myTrack->predictPosition(subjectPlane, vecRefPlanes, cmnCorrected, silCalcMode, bVerbose);
        xLabPredictedMetric=predictedPostionMetric->getPositionX();
        yLabPredictedMetric=predictedPostionMetric->getPositionY();
        xDetPredictedMetric =myTrack->getPositionInDetSystem(2*subjectPlane,xLabPredictedMetric,yLabPredictedMetric);;
        yDetPredictedMetric =myTrack->getPositionInDetSystem(2*subjectPlane+1,xLabPredictedMetric,yLabPredictedMetric);;

        if(verbosity>6) predictedPostionMetric->Print(1);

        xPredSigma = predictedPostionMetric->getSigmaX();
        yPredSigma = predictedPostionMetric->getSigmaY();
        xPhi = predictedPostionMetric->getPhiX()*360./TMath::TwoPi();
        yPhi = predictedPostionMetric->getPhiY()*360./TMath::TwoPi();
        chi2x = predictedPostionMetric->getChi2X();
        chi2y = predictedPostionMetric->getChi2Y();
//        cout<<"\t measure posution"<<endl;
        if(isStripAlignment){
            xLabMeasuredMetric = myTrack->getStripXPosition(subjectPlane, yLabPredictedMetric,cmnCorrected,diaCalcMode);
            yLabMeasuredMetric = INVALID_POSITION;
        }
        else{
            if(mode!=TCluster::corEta)
                cerr<<"[Residual]: Silcion plane but wrong mode: "<<mode<<endl;
            xLabMeasuredMetric = myTrack->getPositionInLabFrame(TPlaneProperties::X_COR, subjectPlane, cmnCorrected,mode,myTrack->getEtaIntegral(subjectPlane*2));
            yLabMeasuredMetric = myTrack->getPositionInLabFrame(TPlaneProperties::Y_COR, subjectPlane, cmnCorrected,isStripAlignment?diaCalcMode:mode,myTrack->getEtaIntegral(subjectPlane*2+1));
        }

        TCluster::calculationMode_t currentMode = isStripAlignment?diaCalcMode:mode;
        xDetMeasuredMetric  = myTrack->getXMeasuredClusterPositionMetricSpace(subjectPlane, cmnCorrected, currentMode,myTrack->getEtaIntegral(subjectPlane*2));
        currentMode = mode;
        yDetMeasuredMetric  = myTrack->getYMeasuredClusterPositionMetricSpace(subjectPlane, cmnCorrected, currentMode,myTrack->getEtaIntegral(subjectPlane*2+1));

        if(verbosity>5&&!isStripAlignment) cout<< "\tLabMeasMetric: "<<xLabMeasuredMetric<<"/"<<yLabMeasuredMetric<<endl;
        if(!isStripAlignment&&(xLabMeasuredMetric<-400e3||yLabMeasuredMetric <-400e3)){
            nInvalidXYLabMeasurement ++;
            continue;
        }


        xDelta = xLabMeasuredMetric - xLabPredictedMetric;    //X_OBS-X_Pred
        yDelta = yLabMeasuredMetric - yLabPredictedMetric;    //Y_OBS-Y_Pred

        resxtest = TMath::Abs(TMath::Abs(xDelta - resOld.getXMean()) / resOld.getXSigma());
        resytest = TMath::Abs(TMath::Abs(yDelta - resOld.getYMean()) / resOld.getYSigma());

        clusterSize = myTrack->getClusterSize(subjectPlane*2+cor==TPlaneProperties::X_COR?0:1,0);

        if(isStripAlignment)
            eta = events.at(nEvent).getCluster(subjectPlane,cor==TPlaneProperties::XY_COR?TPlaneProperties::X_COR:cor,0).getEta(cmnCorrected);
        else
            eta = -1;

        Int_t subjectDet = 2*subjectPlane;
        Float_t predHitPosDetCh = myTrack->inChannelDetectorSpace(subjectDet,xDetPredictedMetric);
        TCluster clus = events[nEvent].getCluster(subjectPlane,TPlaneProperties::X_COR,0);
        Float_t channelPosXMeas = clus.getPosition(cmnCorrected,TCluster::maxValue);

        Float_t xDetPredictedMetric = myTrack->getPositionInDetSystem(subjectDet,xLabPredictedMetric,yLabPredictedMetric);
        Float_t relHitPosMeasuredMetric = myTrack->getRelativeHitPosition(subjectDet,xDetMeasuredMetric);
        Float_t relHitPosPredictedMetric = myTrack->getRelativeHitPosition(subjectDet,xDetPredictedMetric);


        bool useEvent=false;
        switch (calcMode){
            case normalCalcMode:
                if(cor==TPlaneProperties::XY_COR||isStripAlignment)  useEvent = resxtest < res_keep_factor && resytest < res_keep_factor;
                else if(cor==TPlaneProperties::X_COR)        useEvent = resxtest < res_keep_factor;
                else if(cor==TPlaneProperties::Y_COR)   useEvent = resytest < res_keep_factor;
                if (!useEvent){
                    calcModeCutX.push_back(resxtest);
                    calcModeCutY.push_back(resytest);
                }
                break;
            case chi2CalcMode:
                if(cor==TPlaneProperties::XY_COR || isStripAlignment)  useEvent = chi2x < maxChi2 && chi2y < maxChi2;
                else if(cor==TPlaneProperties::X_COR)        useEvent = chi2x < maxChi2;
                else if(cor==TPlaneProperties::Y_COR)   useEvent = chi2y < maxChi2;
                if (!useEvent){
                    calcModeCutX.push_back(chi2x);
                    calcModeCutY.push_back(chi2y);
                }
                break;
            case resolutionCalcMode:
                useEvent = true;
                break;
        }
        if (!useEvent) nCalcModeCut++;

        if(useEvent&& isStripAlignment){
            useEvent = useEvent && !telescopeAlignmentEvent[nEvent];
            if (!useEvent) nTelescoepAlignmentEvent++;
        }

        if(isStripAlignment&& false){
            cout<<TString::Format("%6d %6.1f - %6.1f = % 6.2f,",nEvent,xLabMeasuredMetric,xLabPredictedMetric,xDelta);
            cout<<TString::Format("  %6.2f",xDetMeasuredMetric);
            cout<<TString::Format("   %5.2f, %5.2f  %6.2f",predHitPosDetCh,channelPosXMeas,channelPosXMeas-predHitPosDetCh);
            cout<<TString::Format("  %6.1f/%5.1f",relHitPosPredictedMetric,relHitPosMeasuredMetric);
            cout<<TString::Format("  %+08.2f",xDetMeasuredMetric-xDetPredictedMetric);
            cout<<" ->"<<useEvent<<endl;

        }
        if (useEvent && abs(xLabPredictedMetric) > maxXLabMetric){
            cout<<" Invalid xLabPredictedMetric: "<<xLabPredictedMetric<<" / "<<maxXLabMetric<<endl;
            nInvalidxLabPredictedMetric++;
            useEvent = false;
        }
        if (useEvent && abs(yLabPredictedMetric) > maxYLabMetric){
//            cout<<" Invalid yLabPredictedMetric: "<<yLabPredictedMetric<<" / "<<maxYLabMetric<<endl;
            useEvent = false;
            nInvalidyLabPredictedMetric++;
        }
        if (useEvent && nDiaAlignmentStep > 1 && relHitPosPredictedMetric == N_INVALID){
            cout<<" Invalid relHitPosPredictedMetric"<<nDiaAlignmentStep<<" "<<relHitPosPredictedMetric<<"/"<<subjectDet<<"/"<<xDetPredictedMetric<<endl;
            useEvent = false;
            nInvalidPredictedRelPos++;
        }
        if (useEvent && relHitPosMeasuredMetric == N_INVALID){
            cout<<" Invalid relHitPosMeasuredMetric"<<xDetMeasuredMetric<<"/"<<yDetMeasuredMetric<<" -- "<<relHitPosPredictedMetric<<"/"<<subjectDet<<"/"<<relHitPosMeasuredMetric<<endl;
            useEvent = false;
            nInvalidMeasuredRelPos++;
        }
        if (useEvent && settings->IgnoreStripForAlignment(subjectDet,predHitPosDetCh)){
//            cout<<" Ignoring Strip "<<subjectDet<<"/"<<predHitPosDetCh<<" for alignment"<<endl;
            useEvent = false;
            nIgnoreForStripAlignment++;
        }
        if (useEvent && silicon_chi2_cut){
            useEvent = chi2x < maxChi2 && chi2y < maxChi2;
            if (!useEvent){
                predictedPostionMetric->Print();
                cout<<"Invalid Silicon Chi2 Cut: "<<chi2x<<"/"<<chi2y<<" with cut on "<<maxChi2<<endl;
                nInvalidSiliconChi2Cut++;
            }
        }

        vecXDetRelHitPosPredMetricAll.push_back(relHitPosPredictedMetric);
        vecXDetRelHitPosMeasMetricAll.push_back(relHitPosMeasuredMetric);
        vecDeltaXMetricAll.push_back(xDelta);
        vecDeltaYMetricAll.push_back(yDelta);
        vecUsedEventAll.push_back(useEvent);

        if (useEvent) {
            if(verbosity>4)cout<<nEvent<<" add to vector"<<endl;

            vecXLabPredMetric.push_back(xLabPredictedMetric);
            vecYLabPredMetric.push_back(yLabPredictedMetric);

            vecXLabMeasMetric.push_back(xLabMeasuredMetric);
            vecYLabMeasMetric.push_back(yLabMeasuredMetric);

            vecXDetPredMetric.push_back(xDetPredictedMetric);
            vecYDetPredMetric.push_back(yDetPredictedMetric);

            vecXDetPredChannel.push_back(predHitPosDetCh);

            vecXDetMeasMetric.push_back(xDetMeasuredMetric);
            vecYDetMeasMetric.push_back(yDetMeasuredMetric);

            vecXLabDeltaMetric.push_back(xDelta);
            vecYLabDeltaMetric.push_back(yDelta);

            vecXFidValue.push_back(fiducialValueX.at(nEvent));
            vecYFidValue.push_back(fiducialValueY.at(nEvent));
            vecXPhi.push_back(xPhi);
            vecYPhi.push_back(yPhi);
            vecXResPrediction.push_back(xPredSigma);
            vecYResPrediction.push_back(yPredSigma);
            vecXChi2.push_back(chi2x);
            vecYChi2.push_back(chi2y);
//            int det = subjectPlane*2+cor==TPlaneProperties::X_COR?0:1;
            vecEta.push_back(eta);
            vecClusterSize.push_back(clusterSize);
            nUsedEvents++;
            if (verbosity > 4) cout << "Measured: " <<xDetMeasuredMetric << " / " << yDetMeasuredMetric << endl;
            if (verbosity > 4) cout << "Observed: " << xLabMeasuredMetric << " / " << yLabMeasuredMetric << endl;
            if (verbosity > 4) cout << "Predicted: " << predictedPostionMetric->getPositionX() << " / " << predictedPostionMetric->getPositionY() << endl;
            if (verbosity > 4) cout << "Predicted: " << xLabPredictedMetric << " / " << yLabPredictedMetric << endl;
            if (verbosity > 4) cout << "Delta:    " << xDelta << " / " << yDelta << endl;
            if (verbosity > 4) cout << "ResTest:  " << resxtest << " / " << resytest << "\n\n" << endl;
        }
        else{
            nNotUsedEvents++;
        }
        if(predictedPostionMetric){predictedPostionMetric->Delete();predictedPostionMetric=0;}
        if (verbosity > 3) cout << xDelta << " " << yDelta << endl;
    }

    cout<<"using "<<vecXLabDeltaMetric.size() <<"/"<<nUsedEvents<<" Events, of "<<nUsedEvents+nNotUsedEvents<<" "<<(Float_t)nUsedEvents/(Float_t)(nUsedEvents+nNotUsedEvents)*100.<<endl;
    if (verbosity > 2)
        cout << vecXLabDeltaMetric.size() << " " << vecYLabDeltaMetric.size() << " " << vecXLabPredMetric.size() << " " << vecYLabPredMetric.size() << endl;

    if(nUsedEvents==0){
        cout<<"ResidualCalcMode: ";
        switch (calcMode){
            case chi2CalcMode: cout<<"Chi2 CalcMode";break;
            case normalCalcMode: cout<<"normal CalcMode";break;
            case resolutionCalcMode: cout<<"resolution CalcMode";break;
            default:
                cout<<"Unknown CalcMode";break;
        }
        cout << endl;
        cout << "nNonTelescopeAlignmentEvents: "<< nNonTelescopeAlignmentEvents << endl;
        cout << "nInvalidxLabPredictedMetric:  "<< nInvalidxLabPredictedMetric << endl;
        cout << "nInvalidyLabPredictedMetric:  "<< nInvalidyLabPredictedMetric << endl;
        cout << "nInvalidXYLabMeasurement:     "<< nInvalidXYLabMeasurement << endl;
        cout << "nCalcModeCut:                 "<< nCalcModeCut << endl;
        if (nCalcModeCut > 0){
            TH2F* h = histSaver->CreateScatterHisto("hCalcModeCutProblem",calcModeCutY,calcModeCutX,100,100);
            histSaver->SaveHistogram(h);
            delete h;
            TGraph g = histSaver->CreateDipendencyGraph("gCalcModeCutProblem",calcModeCutY,calcModeCutX);
            std::sort(calcModeCutY.begin(),calcModeCutY.end());
            std::sort(calcModeCutX.begin(),calcModeCutX.end());
            histSaver->SaveGraph(&g,"gCalcModeCutProblem");
            cout<<"MinX: "<<calcModeCutX.front()<<" | MaxX: "<<calcModeCutX.back()<<endl;
            cout<<"MinY: "<<calcModeCutY.front()<<" | MaxY: "<<calcModeCutY.back()<<endl;
        }


        cout << "nTelescoepAlignmentEvent:     "<< nTelescoepAlignmentEvent << endl;
        cout << "nInvalidPredictedRelPos:      "<< nInvalidPredictedRelPos << endl;
        cout << "nInvalidMeasuredRelPos:       "<< nInvalidMeasuredRelPos << endl;
        cout << "nIgnoreForStripAlignment:     "<< nIgnoreForStripAlignment << endl;
        cout << "nInvalidSiliconChi2Cut:       "<< nInvalidSiliconChi2Cut << endl;
        cout << "TAlignment::get";
        if(isStripAlignment)cout<<"SingleStrip";
        cout<<"Residual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor);
        cout << " with " << refPlaneString << ", plotting: " << bPlot<< ", using "<<events.size()<<" Events\t";

        cerr<< "cannot calculate Residual/Alignment for 0 Events. BREAK!"<<endl;
        align->PrintResults(1);
        cerr<< "cannot calculate Residual/Alignment for 0 Events. BREAK!"<<endl;

        exit(-1);
    }
    if (verbosity > 2)
        cout << vecXLabDeltaMetric.size() << " " << vecYLabDeltaMetric.size() << " " << vecXLabPredMetric.size() << " "
        << vecYLabPredMetric.size()  <<" "<<vecXDetMeasMetric.size()<<" "<<vecYDetMeasMetric.size()<< endl;
    //first estimate residuals widths
    TResidual res;
    res.setVerbosity(verbosity>4?verbosity-4:0);
    res.setResKeepFactor(res_keep_factor);
    res.calculateResidual(cor, &vecXLabPredMetric, &vecXLabDeltaMetric, &vecYLabPredMetric, &vecYLabDeltaMetric);
    if(isStripAlignment)
        res.Print();
    if(isStripAlignment&& calcMode == chi2CalcMode && bAlign){
        align->setNDiamondAlignmentEvents((UInt_t) vecXLabDeltaMetric.size());
        align->setDiaChi2(settings->getAlignment_chi2());
    }
    bool bChi2 = calcMode==chi2CalcMode;
    bool bUpdateResolutions = calcMode==chi2CalcMode&&bPlot;
    if(isStripAlignment){
        bUpdateResolutions = bAlign;
        this->CreatePlots(cor, subjectPlane,(string) refPlaneString, bPlot, bAlign,false);
    }
    this->CreatePlots(cor, subjectPlane, (string)refPlaneString, bPlot,bUpdateResolutions,bChi2);
    clearMeasuredVectors();
    return res;
}



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
    if (TPlaneProperties::isDiamondPlane(subjectPlane))
    cout<<"[TAlignment::CheckStripDetectorAlignmentChi2] "<<subjectPlane<<" Align: "<<bAlign<<"/"<<bPlot<<" Chi2: "<<maxChi2<<endl;
    TResidual res;
    res = getStripResidual(cor, subjectPlane, vecRefPlanes, bAlign, true,res,diaCalcMode,chi2CalcMode, maxChi2);
////    if (verbosity)
//        res.Print();
    return res;
}
/**
 *
 */
void TAlignment::saveAlignment(TSettings::alignmentMode mode) {
    settings->goToAlignmentRootDir();
//	if(mode == TSettings::normalMode)
    string fileName = settings->getAlignmentFilePath(mode);
    TFile *alignmentFile = new TFile(fileName.c_str(), "RECREATE");
    cout << "TAlignment:saveAlignment(): path: \"" << fileName << "\"" << endl;
    if(mode==TSettings::transparentMode&&verbosity%2==1&&verbosity>6){
        cout<<"Please confirm...."<<flush;
        char t; cin >>t;
    }
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
    UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
    bool isTelescopeAlignment = TPlaneProperties::isSiliconPlane(subjectPlane)&&TPlaneProperties::AreAllSiliconPlanes(vecRefPlanes);

    for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
        if (!isTelescopeAlignment && telescopeAlignmentEvent[nEvent])
            continue;
        TRawEventSaver::showStatusBar(nEvent, events.size());
        myTrack->setEvent(&events.at(nEvent));
        Float_t sumDeltaX = 0;
        Float_t sumDeltaY = 0;
        predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, settings->doCommonModeNoiseCorrection(), TCluster::corEta, false);
        chi2X = predictedPosition->getChi2X();
        chi2Y = predictedPosition->getChi2Y();
        xPhi = predictedPosition->getPhiX()*360./TMath::TwoPi();
        yPhi = predictedPosition->getPhiY()*360./TMath::TwoPi();
        xPredSigma = predictedPosition->getSigmaX();
        yPredSigma = predictedPosition->getSigmaY();
        for (subjectPlane = 0; subjectPlane < 4; subjectPlane++) {
            if (subjectPlane != 0) {
                predictedPosition->Delete();
                predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, settings->doCommonModeNoiseCorrection(), TCluster::corEta, false);
            }
            Float_t deltaX = myTrack->getXPositionMetric(subjectPlane,settings->doCommonModeNoiseCorrection(),TCluster::corEta);
            Float_t deltaY = myTrack->getYPositionMetric(subjectPlane,settings->doCommonModeNoiseCorrection(),TCluster::corEta);

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
    TH1F *histoChi2X = histSaver->CreateDistributionHisto(histName.str(), vecXChi2, 1000, HistogrammSaver::manual,0,20);
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
    TH1F *histoChi2Y = histSaver->CreateDistributionHisto(histName.str(), vecYChi2, 1000, HistogrammSaver::manual,0,20);
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
//	TResidual res = CheckDetectorAlignment(TPlaneProperties::XY_COR,3,0,0,true);
//	CheckDetectorAlignment(TPlaneProperties::XY_COR,3,0,0,true,res);
}


void TAlignment::LoadResolutionFromSettingsFile(){
    cout<<" Loading Resolution From Settings File..."<<flush;
    for(UInt_t pl = 0; pl < TPlaneProperties::getNSiliconPlanes(); pl++){
        align->setXResolution(settings->GetDefaultResolutionX(pl),pl);
        align->setYResolution(settings->GetDefaultResolutionY(pl),pl);
    }
    cout<<"DONE"<<endl;
}


void TAlignment::SetResolutionsWithUserInput() {
    cout<<"The measured Residuals are: "<<endl;
    cout<<"\tPlane"<<"\t"<<"X-COR"<<"\t"<<"Y-COR"<<endl;
    for(int i=0;i<4;i++)
        cout<<"\t"<<setw(5)<<i<<"\t"<<setw(5)<<gausFitValuesX.at(i).second<<"\t"<<setw(5)<<gausFitValuesY.at(i).second<<endl;
    cout<<endl;
    cout<<"UserInput: "<<settings->isUseUserResolutionInput()<<endl;
    if (!settings->isUseUserResolutionInput()){
        LoadResolutionFromSettingsFile();
        return;
    }

    cout<<"Please enter the resolutions of each plane:"<<endl;
    for(int i=0;i<4;i++){
        inputResolution(i,TPlaneProperties::X_COR);
    }
    cout<<"\n And now the resolutions of the Y Detectors"<<endl;
    for(int i=0;i<4;i++){
        inputResolution(i,TPlaneProperties::Y_COR);
    }
}

void TAlignment::inputResolution(UInt_t plane, TPlaneProperties::enumCoordinate cor){
    float resolution=-1;
    bool validValue = false;
//	cin.ignore(100,'\n');
    while (!validValue){
        cout<<"Enter the resolution of detector D"<<plane<<TPlaneProperties::getCoordinateString(cor)<<":\t"<<flush;
        cin>>resolution;
        string checkRes="";
        while (checkRes.size()!=1|| checkRes.find_first_of("jJyYnN")!=0){
            cout<<"Is the resolution for detector D"<<plane<<TPlaneProperties::getCoordinateString(cor)<<": "<<resolution<<" mum correct? [y,Y,n,N]"<<flush;
            cin>>checkRes;
            cin.clear();
            cout<<(checkRes.size()!=1)<<" "<<(checkRes.find_first_of("jJyYnN")!=0)<<endl;
            cout<<((checkRes.size()!=1&& checkRes.find_first_of("jJyYnN")!=0))<<endl;
        }
        if( checkRes.find_first_of("nN")!=string::npos){//not ok
            continue;
        }
        else{
            cout<<"Setting resolution of D"<<plane<<TPlaneProperties::getCoordinateString(cor)<<": "<<resolution<<endl;
            align->setResolution(resolution,plane,cor);
            validValue=true;
        }
    }
}


/**
 * Takes the current resolutions and the current Residuals and tries to converge to new resolutions
 *
 */
void TAlignment::UpdateResolutions(vector<Float_t> residuals, vector<Float_t> resolutions) {
    cout<<"TAlignment::UpdateResolutions"<<endl;
    vector<Float_t> newResolutions;
    newResolutions.resize(resolutions.size());
    vector<Float_t> zPositions;
    zPositions.push_back(-0.2);
    zPositions.push_back(-0);
    zPositions.push_back(+19);
    zPositions.push_back(+19.9);

    residuals.clear();
    residuals.push_back(2);
    residuals.push_back(2);
    residuals.push_back(2.1);
    residuals.push_back(2.1);

    resolutions.clear();
    resolutions.push_back(1.5);
    resolutions.push_back(1.5);
    resolutions.push_back(2);
    resolutions.push_back(2);

    if (residuals.size()!=resolutions.size()||resolutions.size()!=zPositions.size()){
        cerr << "[TAlignment::UpdateResolutions]: Vectors are from differnet size: "<<residuals.size()<<" "<<resolutions.size()<<" "<<zPositions.size()<<endl;
        exit(-1);
    }
    newResolutions = resolutions;
    for(int it = 0; it < 50; it++){
        cout<<"\n***** iteration "<< it<<"/50"<<endl;
        for(UInt_t subject = 0; subject <  residuals.size();subject++){
            vector<Float_t> res,zPos;
            for(UInt_t k=0;k<residuals.size();k++){
                if(k==subject)
                    continue;
                res.push_back(resolutions[k]);
                zPos.push_back(zPositions[k]);
            }
            Float_t position = zPositions[subject];
//			cout<<subject<<" "<<position<<" "<<res.size()<<" "<<zPos.size()<<endl;
            Float_t trackResolution =  myTrack->calculateTrackResolution(position,zPos,res);
            Float_t residual = residuals[subject];
            if(trackResolution>residual){
                cout<<subject<<"  Track Resolution > residual...."<<endl;
                continue;
            }
            Float_t newResolution = TMath::Sqrt(residual*residual - trackResolution*trackResolution);
            Float_t resolution = resolutions[subject];
            if(newResolution/resolution<.5)
                newResolution = resolution *.8;
            else if (newResolution/resolution >1.5)
                newResolution = resolution *1.2;


            cout<<"old Resolution in "<<subject<<": "<<resolutions.at(subject)<<"("<<newResolution/resolution*100.<<") --> "<< newResolution<<endl;
            newResolutions[subject]=newResolution;
        }
        resolutions = newResolutions;
    }
//	cout<<"Resolution is: "<<resolution<<endl;
//    char t; cin>>t;
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

    //do first convergence step
    for (UInt_t plane = 0; plane < 4; plane++) {
        vector<UInt_t> vecRefPlanes;
        for (UInt_t i = 0; i < 4; i++)
            if (i != plane) vecRefPlanes.push_back(i);
        TResidual res = getResidual(TPlaneProperties::XY_COR, plane, vecRefPlanes,true,TResidual(),getClusterCalcMode(plane),chi2CalcMode,maxChi2);
    }
    //get user input
    if (settings->isUseUserResolutionInput())
        SetResolutionsWithUserInput();
    else
        LoadResolutionFromSettingsFile();
    //do at maximum 5 convergence steps...
    int nSilicionResolutionMaxSteps=5;
    nSilicionResolutionMaxSteps=0;
    for(int nSiliconResolutionStep=0;nSiliconResolutionStep<nSilicionResolutionMaxSteps;nSiliconResolutionStep++){
        for (UInt_t plane = 0; plane < 4; plane++) {
            vector<UInt_t> vecRefPlanes;
            for (UInt_t i = 0; i < 4; i++)
                if (i != plane) vecRefPlanes.push_back(i);
            TResidual res = getResidual(TPlaneProperties::XY_COR, plane, vecRefPlanes,true,TResidual(),getClusterCalcMode(plane),chi2CalcMode,maxChi2);
        }
    }

}




TString TAlignment::GetPlotPreName(UInt_t subjectPlane){
    TString preName="";
    if (subjectPlane == 4) {
        preName  = "hDiamond_";
        if (nDiaAlignmentStep == -1)
            preName.Append("PreAlignment");
        else if (nDiaAlignmentStep == nDiaAlignSteps)
            preName.Append("PostAlignment");
        else
            preName.Append(TString::Format("%d_Step", nDiaAlignmentStep));
    } else {
        preName = "hSilicon_";
        if (nAlignmentStep == -1)
            preName.Append("PreAlignment");
        else if (nAlignmentStep == nAlignSteps){
            preName.Append("PostAlignment");
        }
        else
            preName.Append(TString::Format("%d_Step",nAlignmentStep));
    }
    return preName;
}

TString TAlignment::GetPlotPostName(bool bChi2){
    TString postName="";
    if(bChi2) postName=TString::Format("with_Chi2_cut");
    return postName;
}
TString TAlignment::GetPlotPostTitle(bool bChi2){
    TString postName="";
    if(bChi2) postName=TString::Format("_on_%.0f",settings->getAlignment_chi2());
    return postName;
}
TString TAlignment::GetPlotPostTitle(TString postName){
    return GetPlotPostTitle(postName.BeginsWith("with_Chi2_cut"));
}

Float_t TAlignment::CreateSigmaOfPredictionXPlots(TPlaneProperties::enumCoordinate cor,UInt_t subjectPlane, TString preName,TString postName,TString refPlaneString, bool bPlot){
    Float_t xPredictionSigma=0;
    if (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) {                                                                     //SigmaOfPredictionX
        TString histName = preName;
        histName.Append(TString::Format("_SigmaOfPredictionX_Plane_%d_with_",subjectPlane));
        histName.Append(refPlaneString + postName);
        TH1F* histo = histSaver->CreateDistributionHisto((string)histName, vecXResPrediction, 512, HistogrammSaver::threeSigma);
        if (!histo)
            cerr<<"Could not CreateDistributionHisto: "<<histName<<endl;
        else{
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->Draw("goff");
            histo->GetXaxis()->SetTitle("Delta X/ #mum");
            histo->GetYaxis()->SetTitle("Number of entries");
            xPredictionSigma =  histSaver->GetMean(vecXResPrediction);
            cout<<"X-Prediction of sigma (*100): "<<xPredictionSigma*100<<"\thisto:"<<histo->GetMean()*100<<endl;
            if(bPlot) histSaver->SaveHistogram(histo);
            delete histo;
        }
    }
    return xPredictionSigma;
}



void TAlignment::CreateDistributionPlotDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, Float_t xPredictionSigma) {
    if(!(cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))
        return;

    TString histName = preName;
    histName.Append(TString::Format("_Distribution_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName);
    if (verbosity >4) cout<<histName<<endl;
    TH1F* histo = 0;
    bool verb = false;;

    if(TPlaneProperties::isDiamondPlane(subjectPlane)&&(nDiaAlignmentStep == nDiaAlignSteps)){
        verb = true;
        if (verb) cout<<"\n[TAlignment::CreateDistributionPlotDeltaX]"<<histName<<endl;
        Float_t pitchWidth = settings->getDiamondPitchWidth();
        Float_t xmin = -1.5*pitchWidth;
        Float_t xmax =  1.5*pitchWidth;
        histo = histSaver->CreateDistributionHisto((string)histName, vecXLabDeltaMetric, 512, HistogrammSaver::manual,xmin,xmax);
        if (verb)cout<<"Pitch width "<< pitchWidth<<" "<<xmin<<"-"<<xmax<<endl;
    }
    else
        histo = histSaver->CreateDistributionHisto((string)histName, vecXLabDeltaMetric, 512, HistogrammSaver::threeSigma);
    if (!histo){
        cerr<<"Could not CreateDistributionHisto: "<<histName<<endl;
        return;
    }
    histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
    histo->Draw("");
    Float_t sigma = histo->GetRMS();
    Float_t fitWidth = sigma * 1.5;
    Float_t mean = histo->GetMean();
    cout << "Alignment for plane " << subjectPlane << endl;
    if(verb) cout<<"bins: "<<histo->GetNbinsX()<<"\tRange: "<<histo->GetXaxis()->GetXmin()<<"-"<<histo->GetXaxis()->GetXmax()<<endl;
    if(verb) cout<<"mean: "<<mean<<"+/-"<<sigma<<" / "<<fitWidth<<endl;
    TF1* fitX=0;
    if(TPlaneProperties::isDiamondPlane(subjectPlane) && this->mode == TSettings::transparentMode){
        fitWidth = 3 * sigma;
        fitX = new TF1("doubleGausFit","gaus(0)+gaus(3)",mean-fitWidth,mean+fitWidth);
//              fitX->SetParLimits(0,0,4*sigma);
        fitX->SetParLimits(1,-4*sigma,4*sigma);
        fitX->SetParLimits(2,0,4*sigma);
//              fitX->SetParLimits(3,0,4*sigma);
        fitX->SetParLimits(4,-4*sigma,4*sigma);
        fitX->SetParLimits(5,0,4*sigma);
        fitX->SetParNames("C_{0}","#mu_{0}","#sigma_{0}","#C_{1}","#mu_{1}","#sigma_{1}");
        fitX->SetParameter(1,mean);
        fitX->SetParameter(2,sigma);
        fitX->SetParameter(4,mean);
        fitX->SetParameter(5,sigma);
        if (verb) cout<<"Double Gaus Fit"<<endl;
    }
    else if(TPlaneProperties::isDiamondPlane(subjectPlane)){
        Float_t xmin = mean-1* fitWidth;
        Float_t xmax = mean+1* fitWidth;
        fitWidth = 3*sigma;
        fitX = new TF1("fit","[0]*TMath::Sqrt(TMath::Pi()/2)*[1]*(TMath::Erf(([2]+[3]-x)/TMath::Sqrt(2)/[1])+TMath::Erf(([3]-[2]+x)/TMath::Sqrt(2)/[1]))",mean-fitWidth,mean+fitWidth);
        fitX->FixParameter(3,settings->getDiamondPitchWidth()/2);//TODO
        fitX->SetParLimits(1,0,sigma);
        fitX->SetParNames("Integral","sigma of Gaus","position");
        fitX->SetParameter(2,0);
        fitX->SetParameter(1,0.1);
        if (verb) cout<<"Box Fit"<<endl;
    }
    else{
        fitX = new TF1("fitGausX", "gaus", mean - fitWidth , mean + fitWidth);
        fitX->SetParameter(1,mean);
        fitX->SetParameter(2,sigma);
        if (verb) cout<<"Gaus Fit"<<endl;
    }
    if (this->mode==TSettings::transparentMode){
        cout<<"Fitfunction: "<<fitX->GetName()<<" "<<fitX->GetTitle()<<endl;
    }
    if (verb) cout<<"FIT: "<<mean-fitWidth<<"-"<<mean+fitWidth<<endl;
    if(verb) cout<<"bins: "<<histo->GetNbinsX()<<"\tRange: "<<histo->GetXaxis()->GetXmin()<<"-"<<histo->GetXaxis()->GetXmax()<<endl;
    if(verb) cout<<"mean: "<<mean<<"+/-"<<sigma<<" / "<<fitWidth<<endl;
    if (fitX)
        fitX->SetNpx(1000);
    histo->Fit(fitX, "Q", "",mean-fitWidth, mean+fitWidth);
    Float_t xRes=0;


    if(TPlaneProperties::isDiamondPlane(subjectPlane) && mode == TSettings::transparentMode){
        //check TODO
    }
    else if(TPlaneProperties::isDiamondPlane(subjectPlane)){
        xRes = TMath::Max(fitX->GetParameter(1),fitX->GetParError(2));
        mean = fitX->GetParameter(2);
    }
    else{
        mean = fitX->GetParameter(1);
        xRes = fitX->GetParameter(2);
        gausFitValuesX.at(subjectPlane)= make_pair(mean,xRes);
    }
    if (verb) cout<<"Fit: "<<mean<<"+/-"<<xRes<<endl;

    if(xRes>0&&bUpdateResolution&&histo->GetEntries()>0){
        if(TPlaneProperties::isDiamondPlane(subjectPlane)){
            cout << "set Resolution via Gaus fit for diamond: " << xRes*100 << " with " << vecXLabDeltaMetric.size() << " Events" << endl;
            if (xRes >0)
                align->setXResolution(xRes,subjectPlane);
            align->setXMean(mean,subjectPlane);
        }
        else{
            cout << "\n\nset Resolution via Gaus-Fit: " << xRes*100 << " with " << vecXLabDeltaMetric.size() << " Events" << endl;
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
            cout<< "xDet: "<<xres2*100 <<" = Sqrt("<<xRes*100<<"^2 - "<<xPredictionSigma*100<<"^2)"<<endl;
            if (xres2 > 0)
                align->setXResolution(xres2, subjectPlane);
            align->setXMean(mean,subjectPlane);
        }
        if(xRes<fitWidth/3&&TPlaneProperties::isDiamondPlane(subjectPlane))
            xRes = fitWidth/3;
        histo->GetXaxis()->SetRangeUser(mean-4 * xRes, mean + 4 * xRes);
    }
    histo->GetXaxis()->SetTitle("Delta X / #mum");
    histo->GetYaxis()->SetTitle("Number of entries");
    if (verb)cout<<"save"<<histo->GetName()<<endl;
    if(verb) cout<<"bins: "<<histo->GetNbinsX()<<"\tRange: "<<histo->GetXaxis()->GetXmin()<<"-"<<histo->GetXaxis()->GetXmax()<<endl;
//    if (verb) { cout<<"Press a key"<<endl; char t; cin>>t;}
    if (bPlot) histSaver->SaveHistogram(histo);
    delete fitX;
    delete histo;
    histo =0;
}

void TAlignment::CreateScatterPlotPredYvsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName+TString::Format("_ScatterPlot_YPred_vs_DeltaX_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
	TString histName2 = preName+TString::Format("_ScatterPlot_YPred_vs_DeltaX_v2_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    TH2F *histo=0;
	TH2F *histo2 = 0;
    if(isSiliconPostAlignment){
        Float_t xmin = -1e9;
        Float_t xmax = +1e9;
        Float_t ymin = -50;
        Float_t ymax = 50;
        histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric,vecYLabPredMetric,256,512,xmin,xmax,ymin,ymax);
        Float_t mean = histo->GetMean(2);
        Float_t sigma = histo->GetRMS(2);
        delete histo;
        histo = 0;
        ymin = mean - 3 * sigma;
        ymax = mean + 3 * sigma;
        histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric,vecYLabPredMetric,256,512,xmin,xmax,ymin,ymax);
        if(!histo || (histo->GetEntries() == 0))
            cout<<"isSiliconPostAlignment: "<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<mean<<"+/-"<<sigma<<endl;
		Float_t mean_tmp = round(histo->GetMean(2)*10.)/10.;
        histo2 = histSaver->CreateScatterHisto((string)histName2, vecXLabDeltaMetric, vecYLabPredMetric, 100, 700, 0., 1e4, mean_tmp - 35., mean_tmp + 35., 0., false);
    }
    else
        histo = histSaver->CreateScatterHisto((string)histName, vecXLabDeltaMetric,vecYLabPredMetric, 256);
		Float_t mean_tmp = round(histo->GetMean(2)*10.)/10.;
        histo2 = histSaver->CreateScatterHisto((string)histName2, vecXLabDeltaMetric, vecYLabPredMetric, 100, 700, 0., 1e4, mean_tmp - 35., mean_tmp + 35., 0., false);

    if (histo->GetEntries() == 0)
        cout<<"entries == 0";
    histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
    if(histo && histo->GetEntries()){
        histo->GetXaxis()->SetTitle("Y Predicted / #mum");
        histo->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    else{
        std::cout<<"Something is wrong CreateScatterPlotPredYvsDeltaX"<<endl;
        char t;
        cin>>t;
    }
    histo2->SetTitle(histo2->GetTitle()+GetPlotPostTitle(postName));
    if(histo2 && histo2->GetEntries()){
        histo2->GetXaxis()->SetTitle("Y Predicted / #mum");
        histo2->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram(histo2);
        delete histo2;
    }
    else{
        std::cout << "Something is wrong CreateScatterPlotPredYvsDeltaX " << histName2 << endl;
		cout << histo2->GetEntries() << endl;
        char t;
        cin>>t;
    }
    histName.Append( "_graph");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecYLabPredMetric);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("predicted Y position / #mum");
    graph.GetYaxis()->SetTitle("delta X / #mum");
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}

void TAlignment::CreateScatterPlotMeasXvsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_XMeasured_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TH2F* histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric, vecXDetMeasMetric,256);
    histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
    if(!histo)
        cerr<<"Could not CreateScatterHisto: "<<histName<<endl;
    else{
        histo->GetXaxis()->SetTitle("X Measured / #mum");
        histo->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXDetMeasMetric);
//      if(&graph==0)
//          cerr<<"Could not CreateDipendencyGraph: "<<histName.str()<<endl;
//      else{
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("measured X  / #mum");
    graph.GetYaxis()->SetTitle("delta X / #mum");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}


void TAlignment::CreateScatterPlotPredXChvsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_XPredictedCht_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TH2F* histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric, vecXDetPredChannel,256);
    if(!histo)
        cerr<<"Could not CreateScatterHisto: "<<histName<<endl;
    else{
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle("X Predicted_{Det} / CH");
        histo->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXDetPredChannel);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("X Predicted_{Det} / CH");
    graph.GetYaxis()->SetTitle("Delta X / #mum");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}

void TAlignment::CreateScatterPlotPredXDetvsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_XPredictedDet_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TH2F* histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric, vecXDetPredMetric,256);
    if(!histo)
        cerr<<"Could not CreateScatterHisto: "<<histName<<endl;
    else{
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle("X Predicted_{Det} / #mum");
        histo->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXDetPredMetric);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("X Predicted_{Det} / #mum");
    graph.GetYaxis()->SetTitle("Delta X / #mum");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}

void TAlignment::CreateScatterChi2vsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {
    bool verb  = true;
    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_Chi2X_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TString xTitle = "Chi2X";
    TString yTitle = "Delta X / #mum";
    TH2F *histo=0;
    histo = histSaver->CreateScatterHisto((string)histName, vecXLabDeltaMetric,vecXChi2, 512);
    if(verb)cout<<" HISTO: "<<histo<<endl;
    if(histo){
        if(verb)cout<<" save: "<<histo->GetName()<<"\n entries: "<<histo->GetEntries()<<endl;
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXChi2);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;

    histName = preName + TString::Format("_ScatterPlot_Chi2Y_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    xTitle = "Chi2Y";
    yTitle = "Delta X / #mum";
    histo=0;
    histo = histSaver->CreateScatterHisto((string)histName, vecXLabDeltaMetric,vecYChi2, 512);
    if(verb)cout<<" HISTO: "<<histo<<endl;
    if(histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        if(verb)cout<<" save: "<<histo->GetName()<<"\n entries: "<<histo->GetEntries()<<endl;
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecYChi2);
    graph.Draw("APL");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}
void TAlignment::CreateScatterPlotPredXvsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot, bool bUpdateResolution, bool isSiliconPostAlignment) {

    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_XPred_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    bool verb = histName.BeginsWith("hDiamond_PostAlignment_ScatterPlot_XPred_vs_DeltaX_Plane_4");
    TString xTitle = "X Predicted / #mum";
    TString yTitle = "Delta X / #mum";
    TH2F *histo=0;
    if (verb){
        std::cout<< "[CreateScatterPlotPredXvsDeltaX]: "<<histName<< " bPlot:"<<bPlot<<" bUpdateResolution:"<<bUpdateResolution<<" isSiliconPostAlignment:"<<isSiliconPostAlignment<<endl;
    }
    if(isSiliconPostAlignment){
        if(verb)cout<<"isSiliconPostAlignment"<<endl;
        Float_t xmin = -1e9;
        Float_t xmax = +1e9;
        Float_t ymin = -50;
        Float_t ymax = 50;
        histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric,vecXLabPredMetric,256,512,xmin,xmax,ymin,ymax);
        if (verb) cout<<"x: "<<xmin <<"-"<<xmax<<"\ty: "<<ymin<<"-"<<ymax<<endl;
        Float_t mean = histo->GetMean(2);
        Float_t sigma = histo->GetRMS(2);
        delete histo;
        ymin = mean - 3 * sigma;
        ymax = mean + 3 * sigma;
        histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric,vecXLabPredMetric,256,512,xmin,xmax,ymin,ymax);
        if (verb) cout<<"x: "<<xmin <<"-"<<xmax<<"\ty: "<<ymin<<"-"<<ymax<<endl;
    }
    else
        histo = histSaver->CreateScatterHisto((string)histName, vecXLabDeltaMetric,vecXLabPredMetric, 512);
    if(verb)cout<<" HISTO: "<<histo<<endl;
    if(histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        if(verb)cout<<" save: "<<histo->GetName()<<"\n entries: "<<histo->GetEntries()<<endl;
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXLabPredMetric);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
//    if (verb){cout<<"Press a key to continue."<<flush;char t;cin>>t;}
}


void TAlignment::CreateScatterPlotEtaVsDeltaX(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot){

    if (!(bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) && mode == TSettings::transparentMode))
        return;
    TString histName = preName;
    histName.Append(TString::Format("_ScatterPlot_Eta_vs_DeltaX_Plane_%d_with_",subjectPlane)+refPlaneString+postName);
    TH2F *histo = histSaver->CreateScatterHisto((string)histName, vecEta, vecXLabDeltaMetric, 256);
//    histo.Draw("goff");
    if(histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetYaxis()->SetTitle("#eta");
        histo->GetXaxis()->SetTitle("Delta X / #mum");
        Float_t minX = histo->GetMean(1) - histo->GetRMS(1)*5;
        Float_t maxX = histo->GetMean(1) + histo->GetRMS(1)*5;
        histo->GetXaxis()->SetRangeUser(minX,maxX);
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecEta, vecXLabDeltaMetric);
    graph.Draw("APL");
    graph.GetYaxis()->SetTitle("#eta");
    graph.GetXaxis()->SetTitle("delta X / #mum");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    TGraph *gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr,(string) histName);
    if(gr) delete gr;
}

void TAlignment::CreateRelHitPosXPredDetMetricVsUseEventPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot){
    if(cor==TPlaneProperties::Y_COR)
        return;
    if(verbosity>2)
        cout<<"[CreateRelHitPosXPredDetMetricVsUseEventPlot] subjectPlane:"<<subjectPlane<<" "<<vecXDetRelHitPosPredMetricAll.size()<<"/"<<vecUsedEventAll.size()<<endl;
    TString histName = preName + TString::Format("_ScatterPlot_RelHitPosXPredDet_vs_UseEvent_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    if(vecXDetRelHitPosPredMetricAll.size()==0 ||(vecXDetRelHitPosPredMetricAll.size()!=vecUsedEventAll.size())){
        cout<<histName<< "\t number of entries: "<<vecXDetRelHitPosPredMetricAll.size()<<"/"<<vecUsedEventAll.size()<<endl;
        return;
    }
    TH2F *histo = histSaver->CreateScatterHisto((string)histName, vecXDetRelHitPosPredMetricAll, vecUsedEventAll, 2,256,-1,2,-100,100);
//    histo.Draw("goff");
    if(histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->Draw("goff");
        histo->GetYaxis()->SetTitle("rel pred hit Pos /#mum");
        histo->GetXaxis()->SetTitle("use event");
        histSaver->SaveHistogram(histo);
    }
    vector<Float_t> relHitPosUsed,relHitPosNotUsed;
    for(UInt_t i=0;i<vecXDetRelHitPosPredMetricAll.size()&& i < vecUsedEventAll.size();i++){
        if(vecUsedEventAll[i])
            relHitPosUsed.push_back(vecXDetRelHitPosPredMetricAll[i]);
        else
            relHitPosNotUsed.push_back(vecXDetRelHitPosPredMetricAll[i]);
    }
    //TODO
    Int_t bin = histo->GetXaxis()->FindBin(1.0);
    TH1D* hProjUsed = histo->ProjectionY((TString)(histName+(TString)"_UsedEvents"),bin,bin);
    if(hProjUsed) hProjUsed->SetLineColor(kGreen);
    bin = histo->GetXaxis()->FindBin(0.0);
    TH1D* hProjNotUsed = histo->ProjectionY((TString)(histName+(TString)"_NotUsedEvents"),bin,bin);
    if(hProjNotUsed) hProjNotUsed->SetLineColor(kRed);

    TString name = preName+(TString)("_StackRelHitPosXPred_")+refPlaneString+postName;
    THStack* stack = new THStack(name,TString::Format(";Rel. pred. hit pos X,plane %d /#mum;number of entries",subjectPlane));
    if(hProjUsed)
        stack->Add(hProjUsed);
    if(hProjNotUsed)
        stack->Add(hProjNotUsed);
    stack->Draw("nostack");
    stack->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
    if(hProjNotUsed||hProjUsed)
        if (stack->GetXaxis())
            stack->GetXaxis()->SetTitle("rel pred hit pos /#mum");

    stack->SetName(name+(TString)("_noStack"));
    histSaver->SaveStack(stack,"nostack");

    if(hProjUsed) hProjUsed->SetFillColor(kGreen);
    if(hProjNotUsed) hProjNotUsed->SetFillColor(kRed);
    stack->SetName("hStackRelPredHitPos");
    histSaver->SaveStack(stack,"");

    stack->SetName(name+(TString)("_hist"));
    histSaver->SaveStack(stack,"hist");

    delete stack;
    if (histo)
        delete histo;
}

void TAlignment::CreateRelHitPosXMeasDetMetricVsUseEventPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot){
    if (cor == TPlaneProperties::Y_COR)
        return;
    if(verbosity>3)cout<<"[CreateRelHitPosMeasDetMetricVsUseEventPlot]"<<vecXDetRelHitPosMeasMetricAll.size()<<"/"<<vecUsedEventAll.size()<<endl;
    if(vecXDetRelHitPosPredMetricAll.size()==0)
        return;
    TString corName =TPlaneProperties::X_COR;
    TString histName = preName;
    histName.Append(TString::Format("_ScatterPlot_RelHitPosMeasDet_vs_UseEvent_Plane_%d",subjectPlane)+corName+(TString)"_with_" + refPlaneString + postName);
    TH2F *histo = histSaver->CreateScatterHisto((string)histName, vecXDetRelHitPosMeasMetricAll, vecUsedEventAll, 2,256);
//    histo.Draw("goff");
    if(histo){
        histo->Draw("goff");
        histo->GetYaxis()->SetTitle("rel. meas. hit pos/#mum");
        histo->GetXaxis()->SetTitle("use event");
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    vector<Float_t> relHitPosUsed,relHitPosNotUsed;
    for(UInt_t i=0;i<vecXDetRelHitPosMeasMetricAll.size()&&vecUsedEventAll.size();i++){
        if(vecUsedEventAll[i])
            relHitPosUsed.push_back(vecXDetRelHitPosMeasMetricAll[i]);
        else
            relHitPosNotUsed.push_back(vecXDetRelHitPosMeasMetricAll[i]);
    }
    TH1F* hProjUsed = histSaver->CreateDistributionHisto((string)(histName+(TString)"_UsedEvents"),relHitPosUsed,256);
    if(hProjUsed) hProjUsed->SetLineColor(kGreen);
    TH1F* hProjNotUsed = histSaver->CreateDistributionHisto((string)(histName+(TString)"_NotUsedEvents"),relHitPosNotUsed,256);
    if(hProjNotUsed) hProjNotUsed->SetLineColor(kRed);

    TString name = preName+ TString::Format("_StackRelMeasHitPos_Plane_%d",subjectPlane) + corName + (TString)("_with_") +refPlaneString+postName;
    THStack* stack = new THStack(name,TString::Format("Rel. meas. hit pos, plane %d",subjectPlane)+corName+(TString)" /#mum");
    if(hProjUsed)
        stack->Add(hProjUsed);
    if(hProjNotUsed)
        stack->Add(hProjNotUsed);
    stack->Draw("goff");
    if(hProjNotUsed||hProjUsed)
        stack->GetXaxis()->SetTitle("rel meas hit pos");

    stack->SetName(name+TString("_nostack"));
    histSaver->SaveStack(stack,"nostack");

    if(hProjUsed) hProjUsed->SetFillColor(kGreen);
    if(hProjNotUsed) hProjNotUsed->SetFillColor(kRed);
//        stack->SetName("hStackRelMeasHitPos");
//        histSaver->SaveStack(stack,"");

    stack->SetName(name+TString("_hist"));
    histSaver->SaveStack(stack,"hist");

    delete stack;
}


void TAlignment::CreateRelHitPosMeasXPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot){
    if(!bPlot)  return;
    if (!(cor == TPlaneProperties::X_COR || cor == TPlaneProperties::XY_COR)) return;
    vector<Float_t >vecRelPos;
    UInt_t subjectDet = TPlaneProperties::getDetNo(cor,subjectPlane);
    for(UInt_t i = 0; i < vecXLabMeasMetric.size(); i++){
        Float_t xDetMeasMetric = vecXDetMeasMetric.at(i);
        Float_t relHitPos = myTrack->getRelativeHitPosition(subjectDet,xDetMeasMetric);
//        if(TPlaneProperties::isDiamondPlane(subjectPlane))
//            channelPos = myTrack->getXPositionInStripDetSystem(subjectPlane,xLabMeasMetric,vecYLabPredMetric[i]);
//        else
//            channelPos = myTrack->getPositionInDetSystem(subjectDet,xLabMeasMetric,yLabMeasMetric);
//        Float_t relPos = channelPos-(int)(channelPos+.5);
//        if(subjectPlane == 4)cout<<i<<" "<<xLabMeasMetric<<"-->"<<channelPos<<" "<<relPos<<"\n";
        vecRelPos.push_back(relHitPos);
    }
    TString histName = preName + TString::Format("_RelHitPosMeasX_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    TH1F* histo = histSaver->CreateDistributionHisto((string)histName, vecRelPos, 256);
    if(histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle("relative Hit Position_{observed} / ch");
        histo->GetYaxis()->SetTitle("number of entries");
        histSaver->SaveHistogram(histo);
        delete histo;
    }
}

void TAlignment::CreateFidValueXVsDeltaX(TPlaneProperties::enumCoordinate cor,
        UInt_t subjectPlane, TString preName, TString postName,
        TString refPlaneString, bool bPlot, bool bUpdateResolution,
        bool isSiliconPostAlignment) {

    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_FidCutX_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TString xTitle = "Fiducial Value X / ch";
    TString yTitle = "delta X /#mum";
    TH2F* histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric, vecXFidValue,256);
    if(!histo)
        cerr<<"Could not CreateScatterHisto: "<<histName<<endl;
    else{
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecXFidValue);
//      if(&graph==0)
//          cerr<<"Could not CreateDipendencyGraph: "<<histName.str()<<endl;
//      else{
    graph.Draw("APL");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;
}

void TAlignment::CreateFidValueYVsDeltaX(TPlaneProperties::enumCoordinate cor,
        UInt_t subjectPlane, TString preName, TString postName,
        TString refPlaneString, bool bPlot, bool bUpdateResolution,
        bool isSiliconPostAlignment) {

    if (!bPlot) return;
    if (!((cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR))) return;
    TString histName = preName + TString::Format("_ScatterPlot_FidCutY_vs_DeltaX_Plane_%d_with_",subjectPlane) +refPlaneString+postName;
    TString xTitle = "Fiducial Value Y / ch";
    TString yTitle = "delta X /#mum";
    TH2F* histo = histSaver->CreateScatterHisto((string)histName,vecXLabDeltaMetric, vecYFidValue,256,256,0,256);
    if(!histo)
        cerr<<"Could not CreateScatterHisto: "<<histName<<endl;
    else{
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram(histo);
        delete histo;
    }
    histName.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)histName, vecXLabDeltaMetric, vecYFidValue);
    graph.Draw("APL");
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    TGraph* gr = (TGraph*) graph.Clone();
    histSaver->SaveGraph(gr, (string)histName);
    if(gr) delete gr;

}

Float_t TAlignment::CreateSigmaOfPredictionYPlots(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    Float_t yPredictionSigma = 0;
    if (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR) {                                                                      //DistributionPlot DeltaX
        TString name =  preName + TString::Format("_SigmaOfPredictionY_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
        TString xTitle = "sigma of prediction Y/ #mum";
        TString yTitle = "number of entries #";
        TH1F* histo = histSaver->CreateDistributionHisto((string)name, vecYResPrediction, 512, HistogrammSaver::threeSigma);
        if (!histo)
            cerr<<"Could not CreateDistributionHisto: "<<name<<endl;
        else{
            histo->Draw("goff");
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->GetXaxis()->SetTitle(xTitle);
            histo->GetYaxis()->SetTitle(yTitle);
            yPredictionSigma =  histSaver->GetMean(vecYResPrediction);
            cout<<"Y-Prediction of sigma (*100): "<<yPredictionSigma*100<<"\thisto:"<<histo->GetMean()*100<<endl;
            if (bPlot) histSaver->SaveHistogram(histo);
            delete histo;
        }
    }
    return yPredictionSigma;
}

void TAlignment::CreateScatterPlotPredXvsDeltaY(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    TString name = preName + TString::Format("_ScatterPlot_XPred_vs_DeltaY_Plane_%d_with_",subjectPlane )+refPlaneString+postName;
	TString name2 = preName + TString::Format("_ScatterPlot_XPred_vs_DeltaY_v2_Plane_%d_with_",subjectPlane )+refPlaneString+postName;
    if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterPlot DeltaY vs Xpred

        TString xTitle = "X predicted /#mum";
        TString yTitle = "Delta Y /#mum";
        if(verbosity>3) cout<<"Save: "<<name<<" "<<flush;
        TH2F *histo = histSaver->CreateScatterHisto((string)name, vecYLabDeltaMetric,  vecXLabPredMetric, 256);
		int     xbins_v2 = 100;
		Float_t  xmin_v2 = 0.;
		Float_t  xmax_v2 = 1e4;
		int     ybins_v2 = 700;
		Float_t  ymin_v2 = -35.;
		Float_t  ymax_v2 = 35.;
		if (!isSiliconPostAlignment) {
			Float_t mean_v2 = round(histo->GetMean(2)*10.)/10.;
			ymin_v2 += mean_v2;
			ymax_v2 += mean_v2;
		}
		TH2F *histo2 = histSaver->CreateScatterHisto((string)name2, vecYLabDeltaMetric,  vecXLabPredMetric, xbins_v2, ybins_v2, xmin_v2, xmax_v2, ymin_v2, ymax_v2, 0., false);
        if(!histo)
            cerr<<"Could not create "<<name<<endl;
        else if (histo->GetEntries()==0)
            cerr<<"histo has 0 entries"<<vecYLabDeltaMetric.size()<<"/"<<vecXLabPredMetric.size()<<endl;
        else{
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->GetXaxis()->SetTitle(xTitle);
            histo->GetYaxis()->SetTitle(yTitle);
            histSaver->SaveHistogram((TH2F*) histo->Clone());
            delete histo;
        }
		if (histo2) {
			histo2->SetTitle(histo2->GetTitle()+GetPlotPostTitle(postName));
			histo2->GetXaxis()->SetTitle(xTitle);
			histo2->GetYaxis()->SetTitle(yTitle);
			histSaver->SaveHistogram((TH2F*) histo2->Clone());
			delete histo2;
		}

        name.Replace(0,1,"g");
        TGraph graph = histSaver->CreateDipendencyGraph((string)name, vecYLabDeltaMetric, vecXLabPredMetric);
        graph.Draw("APL");
        graph.GetXaxis()->SetTitle(xTitle);
        graph.GetYaxis()->SetTitle(yTitle);
        graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
        histSaver->SaveGraph((TGraph*) graph.Clone(), (string)name);
        if(verbosity>3)cout<<" DONE"<<endl;
    }
}

void TAlignment::CreateScatterPlotPredYvsDeltaY(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    TString name = preName + TString::Format("_ScatterPlot_YPred_vs_DeltaY_Plane_%d_with_",subjectPlane )+refPlaneString+postName;
    bool verb = name.BeginsWith("hSilicon_PostAlignment_ScatterPlot_YPred_vs_DeltaY_Plane_1");
    if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterPlot DeltaY vs Xpred
        TString xTitle = "Y predicted /#mum";
        TString yTitle = "Delta Y /#mum";
        if(verbosity>3) cout<<"Save: "<<name<<" "<<flush;
        TH2F *histo = histSaver->CreateScatterHisto((string)name, vecYLabDeltaMetric,  vecYLabPredMetric, 256);
        if(!histo)
            cerr<<"Could not create "<<name<<endl;
        else if (histo->GetEntries()==0)
            cerr<<"histo has 0 entries"<<vecYLabDeltaMetric.size()<<"/"<<vecXLabPredMetric.size()<<endl;
        else{
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->GetXaxis()->SetTitle(xTitle);
            histo->GetYaxis()->SetTitle(yTitle);
            histSaver->SaveHistogram((TH2F*) histo->Clone());
            delete histo;
        }

        name.Replace(0,1,"g");
        TGraph graph = histSaver->CreateDipendencyGraph((string)name, vecYLabDeltaMetric, vecYLabPredMetric);
        graph.Draw("APL");
        graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
        graph.GetXaxis()->SetTitle(xTitle);
        graph.GetYaxis()->SetTitle(yTitle);
        histSaver->SaveGraph((TGraph*) graph.Clone(), (string)name);
        if(verbosity>3)cout<<" DONE"<<endl;
    }
    if(verb){
        cout<<"[CreateScatterPlotPredXvsDeltaY] Pres a key"<<endl;
//        char t;
//        cin>>t;
    }

}

void TAlignment::CreateRelHitPosVsChi2Plots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString){
    if (cor == TPlaneProperties::XY_COR){
        CreateRelHitPosVsChi2Plots(TPlaneProperties::X_COR,subjectPlane,preName,postName,refPlaneString);
        CreateRelHitPosVsChi2Plots(TPlaneProperties::Y_COR,subjectPlane,preName,postName,refPlaneString);
        return;
    }
    vector<Float_t> vecRelPos;
    TString histName;
    vecRelPos.clear();
    vector<Float_t> *v;
    if(cor == TPlaneProperties::Y_COR)
        v = &vecYLabPredMetric;
    else
        v = &vecXLabPredMetric;

    for(UInt_t i = 0; i < v->size(); i++){
        Float_t yPred = v->at(i);
        int subjectDet = subjectPlane*2+1;
        Float_t channelPos = myTrack->inChannelDetectorSpace(subjectDet,yPred);
        Float_t relPos = channelPos-(int)(channelPos+.5);
        vecRelPos.push_back(relPos);
    }

    histName = preName +TString("_RelHitPos")+(TString)TPlaneProperties::getCoordinateString(cor)+  TString::Format("_Plane_%d_with_",subjectPlane) + refPlaneString+postName;
    TH1F* histo = histSaver->CreateDistributionHisto((string)histName, vecRelPos, 512);
    if(histo){
        histo->GetXaxis()->SetTitle("relative Hit Position / ch");
        histo->GetYaxis()->SetTitle("number of entries");
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    histName = preName + (TString)"_relHitPos"+(TString)TPlaneProperties::getCoordinateString(cor);
    histName.Append(TString::Format("vsChi2X_Plane_%d_with",subjectPlane)+refPlaneString+postName);
    TString yTitle = TString("relative Hit Position") +(TString)TPlaneProperties::getCoordinateString(cor) + TString(" / ch");

    TH2F* hist = histSaver->CreateScatterHisto((string)histName, vecRelPos,vecXChi2);
    if(hist){
        hist->GetYaxis()->SetTitle(yTitle);
        hist->GetXaxis()->SetTitle("#chi^{2}_{X}");
        hist->GetZaxis()->SetTitle("number of entries #");
        histSaver->SaveHistogram(hist);
        if(hist)delete hist;
    }
    histName = preName + (TString)"_relHitPos"+(TString)TPlaneProperties::getCoordinateString(cor);
    histName.Append(TString::Format("vsChi2Y_Plane_%d_with",subjectPlane)+refPlaneString+postName);
    hist = histSaver->CreateScatterHisto((string)histName, vecRelPos,vecYChi2);
    if(hist){
        hist->GetYaxis()->SetTitle(yTitle);
        hist->GetXaxis()->SetTitle("#chi^{2}_{Y}");
        hist->GetZaxis()->SetTitle("number of entries #");
        histSaver->SaveHistogram(hist);
        if(hist)delete hist;
    }
}

void TAlignment::CreateResidualVsAngle(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot){
    if(!bPlot)  return;

    TString histName = preName + TString::Format("_DeltaX_Vs_PhiX_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    TH2F* hist;
    if (cor!=TPlaneProperties::Y_COR){
        hist = histSaver->CreateScatterHisto((string)histName,this->vecXLabDeltaMetric,this->vecXPhi);
        if(hist){
            hist->GetXaxis()->SetTitle("Angle #Phi_X");
            hist->GetYaxis()->SetTitle("Delta X");
            hist->GetZaxis()->SetTitle("number of entries");
            histSaver->SaveHistogram(hist);
            delete hist;
        }
        histName = preName + TString::Format("_DeltaX_Vs_PhiY_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
        hist = histSaver->CreateScatterHisto((string)histName,this->vecXLabDeltaMetric,this->vecYPhi);
        if(hist){
            hist->GetXaxis()->SetTitle("Angle #Phi_Y");
            hist->GetYaxis()->SetTitle("Delta X");
            hist->GetZaxis()->SetTitle("number of entries");
            histSaver->SaveHistogram(hist);
            delete hist;
        }
    }

    if (cor!=TPlaneProperties::Y_COR)
            return;
    histName = preName + TString::Format("_DeltaY_Vs_PhiX_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    hist = histSaver->CreateScatterHisto((string)histName,this->vecYLabDeltaMetric,this->vecXPhi);
    if(hist){
        hist->GetXaxis()->SetTitle("Angle #PhiX");
        hist->GetYaxis()->SetTitle("Delta Y");
        hist->GetZaxis()->SetTitle("number of entries");
        histSaver->SaveHistogram(hist);
        delete hist;
    }
    histName = preName + TString::Format("_DeltaY_Vs_PhiY_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    hist = histSaver->CreateScatterHisto((string)histName,this->vecYLabDeltaMetric,this->vecYPhi);
    if(hist){
        hist->GetYaxis()->SetTitle("Angle #Phi_Y");
        hist->GetXaxis()->SetTitle("Delta Y");
        hist->GetZaxis()->SetTitle("number of entries");
        histSaver->SaveHistogram(hist);
        delete hist;
    }
}

void TAlignment::CreateRelHitPosPredXPlot(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,TString preName, TString postName, TString refPlaneString,bool bPlot){
    if(!bPlot)  return;
    if (!(cor == TPlaneProperties::X_COR || cor == TPlaneProperties::XY_COR)) return;
    vector<Float_t >vecRelPos;
    UInt_t subjectDet = TPlaneProperties::getDetNo(cor,subjectPlane);
    for(UInt_t i = 0; i < vecXLabPredMetric.size(); i++){
        Float_t xDetMeasMetric = vecXDetPredMetric.at(i);
        Float_t relHitPos = myTrack->getRelativeHitPosition(subjectDet,xDetMeasMetric);
//        if(TPlaneProperties::isDiamondPlane(subjectPlane))
//            channelPos = myTrack->getXPositionInStripDetSystem(subjectPlane,xLabMeasMetric,vecYLabPredMetric[i]);
//        else
//            channelPos = myTrack->getPositionInDetSystem(subjectDet,xLabMeasMetric,yLabMeasMetric);
//        Float_t relPos = channelPos-(int)(channelPos+.5);
//        if(subjectPlane == 4)cout<<i<<" "<<xLabMeasMetric<<"-->"<<channelPos<<" "<<relPos<<"\n";
        vecRelPos.push_back(relHitPos);
    }
    TString histName = preName + TString::Format("_RelHitPosMeasX_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    TH1F* histo = histSaver->CreateDistributionHisto((string)histName, vecRelPos, 256);
    if(histo){
        histo->GetXaxis()->SetTitle("relative Hit Position_{observed} / ch");
        histo->GetYaxis()->SetTitle("number of entries");

        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histSaver->SaveHistogram(histo);
        delete histo;
    }

    histName = preName + TString::Format("_RelHitPosMeasX_Vs_xDetMeasMetric_Plane_%d_with_",subjectPlane)+refPlaneString+postName;
    TH2F* histo2 = histSaver->CreateScatterHisto((string)histName, vecRelPos,vecXDetPredMetric);
    if(histo2){
        histo2->GetYaxis()->SetTitle("relative Hit Position_{observed} / #mum");
        histo2->GetXaxis()->SetTitle("predicted Hit Position_{observed}  / #mum");
        histo2->GetZaxis()->SetTitle("number of entries");
        histo2->SetTitle(histo2->GetTitle()+GetPlotPostTitle(postName));
        histSaver->SaveHistogram(histo2);
        delete histo2;
    }
}

void TAlignment::CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, string refPlaneString, bool bPlot, bool bUpdateResolution, bool bChi2) {
    if (!bPlot && !bUpdateResolution) return;
    if (bPlot)
        if(verbosity>0)cout<<"Save Histograms: "<<  vecXLabDeltaMetric.size() << " " << vecYLabDeltaMetric.size() << " " << vecXLabPredMetric.size() << " " << vecYLabPredMetric.size() << " " << vecXLabMeasMetric.size() << " " << vecYLabMeasMetric.size() << endl;
    // define preName
    TString preName  = GetPlotPreName(subjectPlane);;
    TString postName = GetPlotPostName(bChi2);
    bool isSiliconPostAlignment = (subjectPlane!=4)&&(nAlignmentStep == nAlignSteps);
    isSiliconPostAlignment  = isSiliconPostAlignment  || ((subjectPlane==4) && nDiaAlignSteps == nDiaAlignmentStep);
    if (verbosity>2){
        std::cout << "Create Plots: pre:\""<<preName<<"\"  post\""<<postName<<"\"" << endl;
        std::cout << "   Update Resolution: "<<bUpdateResolution << endl;
        std::cout << "   SubjectPlane "<<subjectPlane << endl;
        std::cout << "   nAlignmentStep: "<<nAlignmentStep << endl;
        std::cout << "   nAlignSteps: " << nAlignSteps << endl;
        std::cout << "   nDiaAlignSteps: " << nDiaAlignSteps << endl;
        std::cout << "   nDiaAlignmentStep " << nDiaAlignmentStep << endl;
        std::cout << "   IsSiliconPostAlignment: "<<isSiliconPostAlignment << endl;;
    }

    stringstream histName;
    if(verbosity<=2){
        cout << "\nCreatePlots with " << preName << " " << (subjectPlane!=4?nAlignmentStep:nDiaAlignmentStep) <<" Step" << flush;
        if (bUpdateResolution )
            cout << "\twith Alignment Resolution Update\n" << endl;
        else
            cout << endl;
    }

    CreateRelHitPosXPredDetMetricVsUseEventPlot( cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateRelHitPosXMeasDetMetricVsUseEventPlot( cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    Float_t xPredictionSigma =-1;
    if(bUpdateResolution)
        xPredictionSigma = CreateSigmaOfPredictionXPlots( cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateDistributionPlotDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,xPredictionSigma);

    CreateScatterPlotPredYvsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterChi2vsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotPredXvsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotMeasXvsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotPredXDetvsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotPredXChvsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotEtaVsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateRelHitPosMeasXPlot(cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateRelHitPosPredXPlot(cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateResidualVsAngle(cor,subjectPlane,preName,postName,refPlaneString,bPlot);
    CreateFidValueXVsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);
    CreateFidValueYVsDeltaX(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);
    CreateScatterPlotDeltaXvsChi2X(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);
    CreateScatterPlotDeltaYvsChi2Y(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);
    Float_t yPredictionSigma = CreateSigmaOfPredictionYPlots(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);

    CreateDistributionPlotDeltaY(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,yPredictionSigma);
    CreateScatterPlotPredXvsDeltaY(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);
    CreateScatterPlotPredYvsDeltaY(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);

    CreateScatterPlotObsXvsObsY(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);

    CreateAngularDistributionPlot(cor,subjectPlane,preName,postName,refPlaneString,bPlot,bUpdateResolution,isSiliconPostAlignment);//,xPredictionSigma);

    if (bPlot && subjectPlane == 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    //DeltaX vs ClusterSize
        histName.str("");
        histName << preName << "_ScatterPlot_ClusterSize_vs_DeltaX_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName;
        if(verbosity>3) cout<<"Save: "<<histName.str()<<" "<<flush;
        TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXLabDeltaMetric, vecClusterSize, 512);
        histo->Draw("goff");
        histo->GetXaxis()->SetTitle("Cluster Size");
        histo->GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveHistogram((TH2F*) histo->Clone());
        histName << "_graph";
        TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXLabDeltaMetric, vecClusterSize);
        graph.Draw("APL");
        graph.GetXaxis()->SetTitle("Cluster Size");
        graph.GetYaxis()->SetTitle("Delta X / #mum");
        histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
        delete histo;
        if(verbosity>3)cout<<"DONE"<<endl;
    }

    if(bPlot){
        vector<Float_t > vecRelPos;

        vecRelPos.clear();

        if(cor == TPlaneProperties::X_COR || cor == TPlaneProperties::XY_COR){
            for(UInt_t i = 0; i < vecXLabPredMetric.size(); i++){
//                Float_t xLabPredictedMetric = vecXLabPredMetric.at(i);
//                Float_t yLabPredictedMetric = vecYLabPredMetric.at(i);
                Float_t xPositionMeasuredMetric = vecXLabMeasMetric[i];
                int subjectDet = subjectPlane*2;
                Float_t channelPos = myTrack->inChannelDetectorSpace(subjectDet,xPositionMeasuredMetric);
                Float_t relPos = channelPos-(int)(channelPos+.5);
                if(subjectPlane == 4 && verbosity > 6)cout<<i<<" "<<xPositionMeasuredMetric<<"-->"<<channelPos<<" "<<relPos<<"\n";
                vecRelPos.push_back(relPos);
            }

            histName.str("");
            histName.clear();
            histName << preName << "_RelHitPosX_Plane_" << subjectPlane << "_with_" << refPlaneString<<postName;
            if(subjectPlane == 4)
                cout<<histName.str()<<endl;
            TH1F* histo = histSaver->CreateDistributionHisto(histName.str(), vecRelPos, 512);
            if(histo){
                histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(bChi2));
                histo->GetXaxis()->SetTitle("relative Hit Position / ch");
                histo->GetYaxis()->SetTitle("number of entries");
                histSaver->SaveHistogram(histo);
                delete histo;
            }
        }

        CreateRelHitPosVsChi2Plots(cor,subjectPlane,preName,postName,refPlaneString);
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
            TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, settings->doCommonModeNoiseCorrection(), TCluster::corEta, false);
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

    for (nEvent = 0; nEvent < events.size(); nEvent++) {
        TRawEventSaver::showStatusBar(nEvent, events.size());
        myTrack->setEvent(&events.at(nEvent));

        for (UInt_t subjectPlane = 0; subjectPlane < TPlaneProperties::getNSiliconPlanes(); subjectPlane++) {
            vector<UInt_t> refPlanes;
            for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
                if (subjectPlane != refPlane) refPlanes.push_back(refPlane);
            bool isTelescopeAlignment = TPlaneProperties::isSiliconPlane(subjectPlane) && TPlaneProperties::AreAllSiliconPlanes(refPlanes);
            if(!isTelescopeAlignment&&telescopeAlignmentEvent[nEvent])
                continue;
            TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, settings->doCommonModeNoiseCorrection(), TCluster::corEta, false);

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
    vecXLabPredMetric.clear();
    vecYLabPredMetric.clear();
    vecXLabMeasMetric.clear();
    vecYLabMeasMetric.clear();
    vecXLabDeltaMetric.clear();
    vecYLabDeltaMetric.clear();
    vecXDetPredMetric.clear();
    vecYDetPredMetric.clear();
    vecXDetPredChannel.clear();
    vecXFidValue.clear();
    vecYFidValue.clear();
    vecXChi2.clear();
    vecYChi2.clear();
    vecXPhi.clear();
    vecYPhi.clear();
    vecXDetMeasMetric.clear();
    vecClusterSize.clear();
    vecYDetMeasMetric.clear();
    vecXResPrediction.clear();
    vecYResPrediction.clear();
    vecEta.clear();
    vecUsedEventAll.clear();
    vecXDetRelHitPosPredMetricAll.clear();
    vecXDetRelHitPosMeasMetricAll.clear();
    vecDeltaXMetricAll.clear();
    vecDeltaYMetricAll.clear();

}

void TAlignment::CreateDistributionPlotDeltaY
(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, Float_t yPredictionSigma) {


    if (subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {        //DistributionPlot DeltaY
        TString name = preName + TString::Format("_Distribution_DeltaY_Plane_%d_with_",subjectPlane) + refPlaneString + postName;
        TString xTitle = "Delta Y /#mum";
        TString yTitle = "number fo entries #";
        if(verbosity>3)cout<<"Save: "<<name<<flush;

        TH1F* histo = (TH1F*) histSaver->CreateDistributionHisto((string)name, vecYLabDeltaMetric, 512, HistogrammSaver::threeSigma);
        bool verb = name.BeginsWith("hSilicon_PostAlignment_Distribution_DeltaY_Plane_0");
        if(!histo)
            cerr<<"Could not CreateDistributionHisto: "<<name<<endl;
        else{
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->Draw("goff");
            histo->GetXaxis()->SetTitle(xTitle);
            histo->GetYaxis()->SetTitle(yTitle);
            Float_t sigma = histo->GetRMS();
            Float_t fitWidth = sigma *1.5;
            Float_t mean = histo->GetMean();
            if (verb){
                cout<<"Name: "<<name<<endl;
                cout<<"mean "<<mean<<endl;
                cout<<"sigma "<<sigma<<endl;
                cout<<"fitWidth "<<fitWidth<<endl;
            }
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

            gausFitValuesY.at(subjectPlane)= make_pair(mean,yRes);
            if (bUpdateResolution&&histo->GetEntries()>0 && yRes > 0) {
                cout << "\n\nset Y-Resolution via Gaus-Fit: " << yRes*100 << " with " << vecYLabDeltaMetric.size() << " Events" << endl;
//                cout << "yRes: "<<yRes*100<<endl;
//                cout << "yPre: "<<yPredictionSigma*100<<endl;
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
                Float_t xmin = mean - 4 * yRes;
                Float_t xmax = mean + 4 * yRes;
                Int_t binLow = histo->GetXaxis()->FindBin(xmin);
                Int_t binUp  = histo->GetXaxis()->FindBin(xmax);
                if (binUp-binLow>5)
                    histo->GetXaxis()->SetRangeUser(mean-4 * yRes,mean +4 * yRes);
                if (verb){
                    cout<<"mean "<<mean<<endl;
                    cout<<"yRes "<<yRes<<endl;
                    cout<<"yres2 "<<yres2<<endl;
                    cout<<"yPredictionSigma "<<yPredictionSigma<<endl;
                    cout<<"xmin "<<xmin<<endl;
                    cout<<"xmax "<<xmax<<endl;
                    cout<<"binLow "<<binLow<<endl;
                    cout<<"binUp "<<binUp<<endl;
                }

            }
            if (bPlot) histSaver->SaveHistogram(histo);
            if (verb){
//                std::cout<<"Press a key"<<endl;
//                char t; cin >>t;
            }
            if(verbosity>3)cout<<" DONE"<<endl;
            delete fitGausY;
            delete histo;
        }
    }

}

void TAlignment::CreateScatterPlotObsXvsObsY(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR)) {   //ScatterHisto XObs vs YObs
        TString name = preName + TString::Format("_ScatterPlot_XObs_vs_YObs_Plane_%d_with_",subjectPlane) + refPlaneString +postName;
        if(verbosity>3) cout<<"Save: "<<name<<" "<<flush;
        TH2F *histo = histSaver->CreateScatterHisto((string)name, vecYLabMeasMetric, vecXLabMeasMetric,512);
        if(!histo)
            cerr<<"Could not create CreateScatterHisto: "<<name<<endl;
        else{
            histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
            histo->GetXaxis()->SetTitle("XObs / #mum");
            histo->GetYaxis()->SetTitle("YObs / #mum");
            histSaver->SaveHistogram(histo);    //,histName.str());
            delete histo;
        }
        if(verbosity>3)cout<<"DONE"<<endl;
    }
}

void TAlignment::CreateScatterPlotDeltaXvsChi2X(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot)
        return;
    if (cor == TPlaneProperties::Y_COR)
        return;

    TString name = preName+ TString::Format("_ScatterPlot_DeltaX_vs_Chi2X_Plane_%d_with_",subjectPlane)+refPlaneString + postName;
    TString xTitle = "#chi^{2}_{X}";
    TString yTitle ="Sum of Delta X / #mum";
    if(verbosity>3) cout<<"Save: "<<(string)name<<" "<<flush;
    TH2F *histo = histSaver->CreateScatterHisto((string)name, vecXLabDeltaMetric, vecXChi2, 256);
//    histo->Draw("goff");
    if (histo){
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram((TH2F*) histo->Clone());
    }
    name.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)name, vecXLabDeltaMetric, vecXChi2);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    histSaver->SaveGraph((TGraph*) graph.Clone(), (string)name);
    delete histo;
    if(verbosity>3)cout<<"DONE"<<endl;
}

void TAlignment::CreateScatterPlotDeltaYvsChi2Y(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (!bPlot)
        return;
    if (cor == TPlaneProperties::X_COR)
        return;
    TString name = preName+ TString::Format("_ScatterPlot_DeltaY_vs_Chi2Y_Plane_%d_with_",subjectPlane)+refPlaneString + postName;
    TString xTitle = "#chi^{2}_{Y}";
    TString yTitle ="Sum of Delta Y / #mum";
    if(verbosity>3) cout<<"Save: "<<(string)name<<" "<<flush;
    TH2F *histo = histSaver->CreateScatterHisto((string)name, vecYLabDeltaMetric, vecYChi2, 256);
    histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetYaxis()->SetTitle(yTitle);

    histSaver->SaveHistogram((TH2F*) histo->Clone());
    name.Replace(0,1,"g");
    TGraph graph = histSaver->CreateDipendencyGraph((string)name, vecYLabDeltaMetric, vecYChi2);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle(xTitle);
    graph.GetYaxis()->SetTitle(yTitle);
    graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
    histSaver->SaveGraph((TGraph*) graph.Clone(), (string)name);
    delete histo;
    if(verbosity>3)cout<<"DONE"<<endl;
}

void TAlignment::CreateAngularDistributionPlot(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString, bool bPlot,
        bool bUpdateResolution, bool isSiliconPostAlignment) {
    if (bPlot ){      //hAngularDistribution
        TString name  = preName + TString::Format("_AngularDistribution_for_Plan_%d_with_",subjectPlane) + refPlaneString + postName;
        TString xTitle = "#Phi_{X} / degree";
        TString yTitle = "#Phi_{Y} / degree";
        if(verbosity>3) cout<<"Save: "<<name<<" "<<flush;
        TH2F *histo = histSaver->CreateScatterHisto((string)name,vecXPhi,vecYPhi,512);
        histo->SetTitle(histo->GetTitle()+GetPlotPostTitle(postName));
        histo->GetXaxis()->SetTitle(xTitle);
        histo->GetYaxis()->SetTitle(yTitle);
        histSaver->SaveHistogram((TH2F*) histo->Clone());
        name.Replace(0,1,"g");
        TGraph graph = histSaver->CreateDipendencyGraph((string)name, vecYPhi, vecXPhi);
        graph.Draw("APL");
        graph.GetXaxis()->SetTitle(xTitle);
        graph.GetYaxis()->SetTitle(yTitle);
        graph.SetTitle(graph.GetTitle()+GetPlotPostTitle(postName));
        histSaver->SaveGraph((TGraph*)graph.Clone(),(string)name);
        delete histo;
        if(verbosity>3)cout<<"DONE"<<endl;
    }

}

void TAlignment::CreateChi2DistributionPlots(
        TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane,
        TString preName, TString postName, TString refPlaneString) {
}
