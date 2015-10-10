/*
 * TAnalysisOf3DShortAnalysis.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DShortAnalysis.hh>

TAnalysisOf3DShortAnalysis::TAnalysisOf3DShortAnalysis(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna) {
    // TODO Auto-generated constructor stub
    this->settings = settings;
    this->histSaver = histSaver;
    useCMN = true;
    verbosity = 10;//settings->getVerbosity();

    PulseHeightBins = 256;
    PulseHeightMin = 1;
    PulseHeightMax = 2800;
    PulseHeightMinMeanCharge = 1;
    PulseHeightMaxMeanCharge = 1500;
//    histos = new THistogramManager(histSaver);

    vecEdgePredX.resize(settings->get3dEdgeFidCuts()->getNFidCuts());
    vecEdgePredY.resize(settings->get3dEdgeFidCuts()->getNFidCuts());
    vecEdgePulseHeight.resize(settings->get3dEdgeFidCuts()->getNFidCuts());

    this->bTransAna = bTransAna;
    if (this->bTransAna)
        appendix ="_trans";
    else
        appendix ="";
}

TAnalysisOf3DShortAnalysis::~TAnalysisOf3DShortAnalysis() {
    // TODO Auto-generated destructor stub
//    delete histos;
}

void TAnalysisOf3DShortAnalysis::setEventReader(TTracking* eventReader) {
    this->eventReader = eventReader;
}

void TAnalysisOf3DShortAnalysis::initHistos() {
    cout<<"[TAnalysisOf3dDiamonds::initialiseShortAnalysisHistos()]"<<endl;
    //Universal histograms
    TString name = "hRelativeChargeTwoClustersX";
    name.Append(appendix);
    Float_t xmax = settings->get3dMetallisationFidCuts()->getXHigh();
    Float_t ymax = settings->get3dMetallisationFidCuts()->getYHigh();
    //histos->addHistogram("TH1F",name,name+TString(";X;rel. ph: ph_{clus1}/ph_{clus2}"),"",1024,0,xmax,100,0,20);
    hRelativeChargeTwoClustersX = new TH2F(name,name,1024,0,xmax,100,0,20);
    hRelativeChargeTwoClustersX->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersX->GetYaxis()->SetTitle("rel. ph: ph_{clus1}/ph_{clus2}");

    name ="hFidCutsVsMeanCharge";
    name.Append(appendix);
    hFidCutsVsMeanCharge = new TProfile2D(name,name, 256,-.3*xmax,xmax*1.3,256,-.3*ymax,ymax*1.3);
    hFidCutsVsMeanCharge->GetXaxis()->SetTitle("X / #mum");
    hFidCutsVsMeanCharge->GetYaxis()->SetTitle("Y / #mum");
    hFidCutsVsMeanCharge->GetZaxis()->SetTitle("avrg total Charge of clusters - nClusters <= 2");

    name = "hRelativeChargeTwoClustersY";
    name.Append(appendix);
    hRelativeChargeTwoClustersY = new TH2F(name,name,1024,0,ymax,100,0,20);
    hRelativeChargeTwoClustersY->GetXaxis()->SetTitle("Y");
    hRelativeChargeTwoClustersY->GetYaxis()->SetTitle("rel. ph: ph_{clus1}/ph_{clus2}");

    name="hRelativeChargeTwoClustersXY";
    name.Append(appendix);
    hRelativeChargeTwoClustersXY = new TProfile2D(name,name,128,0,xmax,128,0,ymax);
    hRelativeChargeTwoClustersXY->GetXaxis()->SetTitle("X");
    hRelativeChargeTwoClustersXY->GetYaxis()->SetTitle("Y");
    hRelativeChargeTwoClustersXY->GetZaxis()->SetTitle("avrg. rel. ph: ph{clus1}/ph_{clus2}");

    name ="hShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    hShortAnalysis2TotalChargeXY = new TProfile2D(name,name,128,0,xmax,128,0,ymax);
    hShortAnalysis2TotalChargeXY->GetXaxis()->SetTitle("X");
    hShortAnalysis2TotalChargeXY->GetYaxis()->SetTitle("Y");
    hShortAnalysis2TotalChargeXY->GetZaxis()->SetTitle("total charge: ph_{clus_1{}} + ph_{clus_{2}} /ADC");

    name = "hShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    hTotalAvrgChargeXY = new TProfile2D(name,name, xmax/10.,0,xmax,ymax/10.,0,ymax);
    hTotalAvrgChargeXY->Draw();
    hTotalAvrgChargeXY->GetXaxis()->SetTitle("#it{x} position / #mum");
    hTotalAvrgChargeXY->GetYaxis()->SetTitle("#it{y} position / #mum");
    hTotalAvrgChargeXY->GetZaxis()->SetTitle("avrg signal / ADC");
    hTotalAvrgChargeXY->GetZaxis()->SetTitleOffset(1.2);

    name = "hRelatviveNumberOfMultipleClusterEvents";
    name.Append(appendix);
    hRelatviveNumberOfMultipleClusterEvents = new TProfile(name,name,3,.5,3.5);
    hRelatviveNumberOfMultipleClusterEvents->GetXaxis()->SetTitle("predicted pattern");
    hRelatviveNumberOfMultipleClusterEvents->GetYaxis()->SetTitle("#_{multiple cluster events}/#_{total Events}");

    name = "hRelatviveNumberOfMultipleClusterEventsSamePattern";
    name.Append(appendix);
    hRelatviveNumberOfMultipleClusterEventsSamePattern = new TProfile(name,name,3,.5,3.5);
    hRelatviveNumberOfMultipleClusterEventsSamePattern->GetXaxis()->SetTitle("predicted pattern");
    hRelatviveNumberOfMultipleClusterEventsSamePattern->GetYaxis()->SetTitle("#_{multiple cluster events}/#_{total Events}");

    //hNumberofClusters
    name = "hNumberofClusters";
    hNumberofClusters = new TH1F(name,name,4,0,4);
    name.Append(appendix);
    hNumberofClusters->GetXaxis()->SetTitle("Number of Clusters");
    hNumberofClusters->GetYaxis()->SetTitle("Number of Entries #");

    //hEventsvsChannelCombined
    name = "hEventsvsChannelCombined";
    name.Append(appendix);
    hEventsvsChannelCombined = new TH1F(name,name,100,0,100);
    hEventsvsChannelCombined->GetXaxis()->SetTitle("Channel");
    hEventsvsChannelCombined->GetYaxis()->SetTitle("Number of Entries #");

    //hDoubleClusterPos
    name="hDoubleClusterPos";
    name.Append(appendix);
    hDoubleClusterPos = new TH1F(name,name,80,20,100);
    hDoubleClusterPos->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos->GetYaxis()->SetTitle("Number of Entries #");

    //hDoubleClusterPos0
    name ="hDoubleClusterPos0";
    name.Append(appendix);
    hDoubleClusterPos0 = new TH1F(name,name,80,20,100);
    hDoubleClusterPos0->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos0->GetYaxis()->SetTitle("Number of Entries #");
    hDoubleClusterPos0->SetFillColor(2);

    //hDoubleClusterPos1
    name = "hDoubleClusterPos1";
    name.Append(appendix);
    hDoubleClusterPos1 = new TH1F(name,name,80,20,100);
    hDoubleClusterPos1->GetXaxis()->SetTitle("HighestPH Channel Hit");
    hDoubleClusterPos1->GetYaxis()->SetTitle("Number of Entries #");
    hDoubleClusterPos1->SetFillColor(3);

    //hLandauCluster1
    name = "hLandauCluster1";
    name.Append(appendix);
    hLandauCluster1 = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauCluster1->GetXaxis()->SetTitle("Number of Clusters");
    hLandauCluster1->GetYaxis()->SetTitle("Number of Entries #");
    hLandauCluster1->SetFillColor(2);

    //hLandauCluster2
    name = "hLandauCluster2";
    name.Append(appendix);
    hLandauCluster2 = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauCluster2->GetXaxis()->SetTitle("Number of Clusters");
    hLandauCluster2->GetYaxis()->SetTitle("Number of Entries #");
    hLandauCluster2->SetFillColor(3);

    //hLandauDoubleCombined
    name = "hLandauDoubleCombined";
    name.Append(appendix);
    hLandauDoubleCombined = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauDoubleCombined->GetXaxis()->SetTitle("Number of Clusters");
    hLandauDoubleCombined->GetYaxis()->SetTitle("Number of Entries #");

    if(verbosity) cout<<"loop over patterns"<<endl;
    for(UInt_t i=0; i<settings->diamondPattern.getNIntervals(); i++){
        if(verbosity) cout<<"Loop: "<<i<<" "<<settings->diamondPattern.getNIntervals()<<endl;
        pair<int,int> channels =settings->diamondPattern.getPatternChannels(i+1);
        //hLandau
        name = TString::Format("hLandau_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        if(verbosity>1) cout<<"Create "<<name<<endl;
        hLandau.push_back(new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax));
        if(hLandau.back()){
            hLandau.back()->GetXaxis()->SetTitle("PH of diamond cluster");
            hLandau.back()->GetYaxis()->SetTitle("number of entries #");
        }
        else
            cerr<<"hLandau:'"<<name<<"' wasn't created correctly"<<endl;

        //hPHvsChannel
        name = TString::Format("hEventsvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hEventsvsChannel.push_back(new TH1F(name,name,100,0,100));
        if(hEventsvsChannel.back()){
            hEventsvsChannel.back()->GetXaxis()->SetTitle("HighestPH [ch]");
            hEventsvsChannel.back()->GetYaxis()->SetTitle("No. Events");
        }

        //hPHvsChannel
        name = TString::Format("hPHvsChannel_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hPHvsChannel.push_back(new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,100,0,100));
        if(hPHvsChannel.back()){
            hPHvsChannel.back()->Draw();
            hPHvsChannel.back()->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->GetXaxis()->SetLimits(PulseHeightMin,PulseHeightMax);
            hPHvsChannel.back()->SetAxisRange(PulseHeightMin, PulseHeightMax, "X");
            hPHvsChannel.back()->GetYaxis()->SetRangeUser(channels.first-1,channels.second+1);
            hPHvsChannel.back()->GetXaxis()->SetTitle("cluster charge /ADC");
            hPHvsChannel.back()->GetYaxis()->SetTitle("XPos /ch");
            //hPHvsChannel.back()->GetXaxis()->SetRange(PulseHeightMin,PulseHeightMax);
        }

        //hPHvsChannel.at(i)->SetMaximum(3000);
        //hPHvsChannel.at(i)->SetMinimum(0);

        //hHitandSeedCount
        name = TString::Format("hHitandSeedCount_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hHitandSeedCount.push_back(new TH2F(name,name,10,-.5,9.5,10,-.5,9.5));
        hHitandSeedCount.back()->GetXaxis()->SetTitle("Hit Count");
        hHitandSeedCount.back()->GetYaxis()->SetTitle("Seed Count");

        //hChi2XChi2Y
        name = TString::Format("hChi2X_vs_chi2Y_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hChi2XChi2Y.push_back(new TH2F(name,name,60,0,30,60,0,30));
        hChi2XChi2Y.back()->GetXaxis()->SetTitle("#chi^{2}_{X}");
        hChi2XChi2Y.back()->GetYaxis()->SetTitle("#chi^{2}_{Y}");

        //hFidCutXvsFidCutY
        name = TString::Format("hFidCutXvsFidCutY_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hFidCutXvsFidCutY.push_back(new TH2F(name,name,160,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),120,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutY.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutY.at(i)->GetYaxis()->SetTitle("FidCutY");

        //hFidCutXvsFidCutYvsCharge     For TH2D

        name = TString::Format("hFidCutXvsFidCutYvsCharge_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        name.Append(appendix);
        hFidCutXvsFidCutYvsCharge.push_back(new TH2D(name,name,
                213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
                160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutYvsCharge.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutYvsCharge.at(i)->GetYaxis()->SetTitle("FidCutY");
        hFidCutXvsFidCutYvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        name = "hFidCutXvsFidCutYvsEvents";
        name.Append(appendix);
        hFidCutXvsFidCutYvsEvents.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone(name));

        //hFidCutXvsFidCutYvsMeanCharge
        name = "hFidCutXvsFidCutYvsMeanCharge";
        name.Append(appendix);
        hFidCutXvsFidCutYvsMeanCharge.push_back((TH2D*)hFidCutXvsFidCutYvsCharge.at(i)->Clone(name));

        //hXdetvsYdetvsCharge       For TH2D     -------> Error is around here!!!!!
        Int_t DiamondPattern = i+1;
        Float_t xLow = settings->get3dMetallisationFidCuts()->getXLow(DiamondPattern);
        Float_t xHigh = settings->get3dMetallisationFidCuts()->getXHigh(DiamondPattern);
        Float_t yLow = settings->get3dMetallisationFidCuts()->getYLow(DiamondPattern);
        Float_t yHigh = settings->get3dMetallisationFidCuts()->getYHigh(DiamondPattern);
        if(verbosity) cout<<"("<<xLow<<"-"<<xHigh<<", "<<yLow<<"-"<<yHigh<<")"<<endl;
        Float_t xDiv = (xHigh - xLow)/5;
        Float_t yDiv = (yHigh - yLow)/5;
        name = TString::Format("hXdetvsYdetvsCharge_pattern_%d_ch_%02d_to_%02d",i,channels.first,channels.second);
        //name = TString::Format("hXdetvsYdetvsCharge%02d_%02d%",channels.first,channels.second);
        name.Append(appendix);
        hXdetvsYdetvsCharge.push_back(new TH2D(name,name,xDiv,xLow,xHigh,yDiv,yLow,yHigh));
        hXdetvsYdetvsCharge.at(i)->GetXaxis()->SetTitle("X (um)");
        hXdetvsYdetvsCharge.at(i)->GetYaxis()->SetTitle("Y (um)");
        hXdetvsYdetvsCharge.at(i)->GetZaxis()->SetTitle("Charge ADC");

        //hFidCutXvsFidCutYvsEvents
        name ="hXdetvsYdetvsEvents";
        name.Append(appendix);
        hXdetvsYdetvsEvents.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone(name));

        //hFidCutXvsFidCutYvsMeanCharge
        name ="hXdetvsYdetvsMeanCharge";
        name.Append(appendix);
        hXdetvsYdetvsMeanCharge.push_back((TH2D*)hXdetvsYdetvsCharge.at(i)->Clone(name));
    }
    if(verbosity>2) cout<<"DONE"<<endl;
    //hFidCutXvsFidCutYvsMeanChargeAllDetectors
    hFidCutXvsFidCutYvsMeanChargeAllDetectors = (TH2D*)hFidCutXvsFidCutYvsCharge.at(0)->Clone("hFidCutXvsFidCutYvsMeanChargeAllDetectors");

    for(int i=0;i<7;i++){
        cout<<" "<<i<<"/7"<<endl;
        //hFidCutXvsFidCutYClusters For TH2D    {0,1,1_1Seed,2_FirstCluster,2_SecondCluster,3}
        name = TString::Format("hFidCutXvsFidCutYClusters%i",i);
        name.Append(appendix);
        hFidCutXvsFidCutYClusters.push_back(new TH2D(name,name,213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh()));
        hFidCutXvsFidCutYClusters.at(i)->GetXaxis()->SetTitle("FidCutX");
        hFidCutXvsFidCutYClusters.at(i)->GetYaxis()->SetTitle("FidCutY");
        hFidCutXvsFidCutYClusters.at(i)->GetZaxis()->SetTitle("Events");
    }

    name ="hShortAnalysis2ClusterHitPattern_1stCluster";
    name.Append(appendix);
    int nFidCuts =settings->get3dMetallisationFidCuts()->getNFidCuts();
    hShortAnalysis2ClusterHitPattern_1stCluster =  new TH2F(name,name,
            nFidCuts,.5,nFidCuts+.5,
            5,-nFidCuts+.5,nFidCuts-.5);
    hShortAnalysis2ClusterHitPattern_1stCluster->GetXaxis()->SetTitle("predicteded hit pattern");
    hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()->SetTitle("cluster_{1}-hit-pattern - predicted-hit-pattern");


    cout<<"hShortAnalysis2ClusterHitPattern_2ndCluster"<<endl;
    name = "hShortAnalysis2ClusterHitPattern_2ndCluster";
    name.Append(appendix);
    hShortAnalysis2ClusterHitPattern_2ndCluster = (TH2F*)hShortAnalysis2ClusterHitPattern_1stCluster->Clone(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->SetTitle(name);
    hShortAnalysis2ClusterHitPattern_2ndCluster->GetYaxis()->SetTitle("cluster_{2}-hit-pattern - predicted-hit-pattern");

    cout<<"Landau"<<endl;
    name = "hLandau3D";
    name.Append(appendix);
    hLandau3DWithColumns = new TH1F(name,"3D",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DWithColumns->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DWithColumns->GetYaxis()->SetTitle("number of entries #");
    hLandau3DWithColumns->SetLineColor(kBlack);

    cout<<"hLandau3DWO"<<endl;
    name = "hLandau3DWO";
    name.Append(appendix);
    cout<<PulseHeightBins<< " "<<PulseHeightMin<<" "<<PulseHeightMax<<" "<<name<<endl;
    hLandau3DWithoutColumns = new TH1F(name,"3D Phantom",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DWithoutColumns->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DWithoutColumns->GetYaxis()->SetTitle("number of entries #");
    hLandau3DWithoutColumns->SetLineColor(kRed);
    hLandau3DWithoutColumns->SetLineStyle(7);

    cout<<"end hLandau3DWO"<<endl;

    cout<<"begin hLandau3DWO_subset"<<endl;
    TString name2="hLandau3DWO_subset";
    name2.Append(appendix);
    name = "hLandau3DWO_subset";
    name.Append(appendix);
    hLandau3DWithoutColumns_subset = new TH1F(name,"3D Phantom, central Region",PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandau3DWithoutColumns_subset->GetXaxis()->SetTitle("charge / ADC");
    hLandau3DWithoutColumns_subset->GetYaxis()->SetTitle("number of entries #");
    hLandau3DWithoutColumns_subset->SetLineColor(kRed);
    hLandau3DWithoutColumns_subset->SetLineWidth(2);
    cout<<"end init"<<endl;


    name = "hLandau3DWithoutColumnsFidCutXvsFidCutY";
    hLandau3DWithoutColumnsFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandau3DWithoutColumnsFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandau3DWithoutColumnsFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");


    name = "hLandau3DWithColumnsFidCutXvsFidCutY";
    hLandau3DWithColumnsFidCutXvsFidCutY = new TH2F(name,name,
            213,settings->getSi_avg_fidcut_xlow(),settings->getSi_avg_fidcut_xhigh(),
            160,settings->getSi_avg_fidcut_ylow(),settings->getSi_avg_fidcut_yhigh());
    hLandau3DWithColumnsFidCutXvsFidCutY->GetXaxis()->SetTitle("FidCutX");
    hLandau3DWithColumnsFidCutXvsFidCutY->GetYaxis()->SetTitle("FidCutY");
}

void TAnalysisOf3DShortAnalysis::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred,
        Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y) {
    xPred = x_pred;
    yPred = y_pred;
    xFid = x_fid;
    yFid = y_fid;
    xChi2 = chi2x;
    yChi2 = chi2y;
    if(!bTransAna){
        if(chi2x>settings->getChi2Cut3D_X()||chi2y>settings->getChi2Cut3D_Y())
            return;

        hNumberofClusters->Fill(eventReader->getNDiamondClusters());
        ClusterPlots(eventReader->getNDiamondClusters());

        Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(x_pred,y_pred);
        FillEdgeAlignmentHistos();
        if (predictedDetector <1 || predictedDetector > settings->getNDiaDetectorAreas())
            return;
        switch (eventReader->getNDiamondClusters()) {
            case 1:
                Analyse1Cluster();
                hRelatviveNumberOfMultipleClusterEvents->Fill(predictedDetector,0);
                hRelatviveNumberOfMultipleClusterEventsSamePattern->Fill(predictedDetector,0);
                break;
            case 2:
                Analyse2Cluster();
            default:
                hRelatviveNumberOfMultipleClusterEvents->Fill(predictedDetector,1);
        }
    }
    else
        Analyse1Cluster();
}

void TAnalysisOf3DShortAnalysis::FillEdgeAlignmentHistos(){
    Int_t nClusters = eventReader->getNDiamondClusters();
    if (nClusters ==0 || nClusters >2)
        return;
    Float_t charge = 0;
    for (UInt_t i = 0; i<nClusters; i++)
        charge += eventReader->getCluster(TPlaneProperties::getDetDiamond(),i).getPositiveCharge(useCMN);
    hFidCutsVsMeanCharge->Fill(xPred,yPred,charge);
    FillEdgeDistributions(charge);
}

void TAnalysisOf3DShortAnalysis::ClusterPlots(int nClusters) {
    if(nClusters==0){
        if(hFidCutXvsFidCutYClusters[0])
            hFidCutXvsFidCutYClusters.at(0)->Fill(xFid,yFid,1);
    }
    if(nClusters==1){
        if(hFidCutXvsFidCutYClusters[1])
            hFidCutXvsFidCutYClusters.at(1)->Fill(xFid,yFid,1);
    }
    if(nClusters==1){
        if(HitCount==0&&SeedCount==1){
            if(hFidCutXvsFidCutYClusters[2])
                hFidCutXvsFidCutYClusters.at(2)->Fill(xFid,yFid,1);
        }
    }
    if(nClusters==2){
        if(hFidCutXvsFidCutYClusters[3])
            hFidCutXvsFidCutYClusters.at(3)->Fill(xFid,yFid,1);
        TCluster diamondCluster0 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
        TCluster diamondCluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);

        hDoubleClusterPos->Fill(diamondCluster0.getHighestSignalChannel());
        hDoubleClusterPos->Fill(diamondCluster1.getHighestSignalChannel());
        if(diamondCluster0.getHighestSignalChannel()==85||diamondCluster1.getHighestSignalChannel()==85){
            hDoubleClusterPos0->Fill(diamondCluster0.getHighestSignalChannel());
            hDoubleClusterPos0->Fill(diamondCluster1.getHighestSignalChannel());
            hLandauCluster1->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
            hFidCutXvsFidCutYClusters.at(4)->Fill(xFid,yFid,1);
        }
        if(diamondCluster0.getHighestSignalChannel()==55||diamondCluster1.getHighestSignalChannel()==55){
            hDoubleClusterPos1->Fill(diamondCluster0.getHighestSignalChannel());
            hDoubleClusterPos1->Fill(diamondCluster1.getHighestSignalChannel());
            hLandauCluster2->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
            hFidCutXvsFidCutYClusters.at(5)->Fill(xFid,yFid,1);
        }
        if((!diamondCluster0.getHighestSignalChannel()==55&&!diamondCluster1.getHighestSignalChannel()==55)||(!diamondCluster0.getHighestSignalChannel()==85&&!diamondCluster1.getHighestSignalChannel()==85)){
            hLandauDoubleCombined->Fill((diamondCluster0.getPositiveCharge(useCMN)+diamondCluster1.getPositiveCharge(useCMN)));
        }

    }
    if(nClusters==3){
        hFidCutXvsFidCutYClusters.at(6)->Fill(xFid,yFid,1);
    }
}

void TAnalysisOf3DShortAnalysis::FillMeanChargeVector(
        Float_t clusterCharge) {
    vecPredDetX.push_back(xPred);
    vecPredDetY.push_back(yPred);
    vecPulseHeight.push_back(clusterCharge);

    Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPred,yPred);
    if (predictedDetector == 2){
        hLandau3DWithoutColumns->Fill(clusterCharge);
        hLandau3DWithoutColumnsFidCutXvsFidCutY->Fill(xFid,yFid);
        if(settings->centralRegion3DnH->IsInFiducialCut(xPred,yPred))
            hLandau3DWithoutColumns_subset->Fill(clusterCharge);
    }
    else if (predictedDetector == 3){
        hLandau3DWithColumns->Fill(clusterCharge);
        hLandau3DWithColumnsFidCutXvsFidCutY->Fill(xFid,yFid);
    }
}

void TAnalysisOf3DShortAnalysis::HitandSeedCount(TCluster* nCluster) {
    int Hit=0;int Seed=0;
    for (UInt_t i=0;i<nCluster->getClusterSize();i++){
        if(nCluster->isHit(i)) Hit++;
        if(nCluster->isSeed(i)) Seed++;
    }
    HitCount=(Hit-Seed);
    SeedCount=Seed;
}

void TAnalysisOf3DShortAnalysis::Analyse1Cluster(UInt_t clusterNo){
    if(!bTransAna){
        clusteredCluster = eventReader->getCluster(TPlaneProperties::getDetDiamond(),clusterNo);
        diamondCluster = &clusteredCluster;
        if(diamondCluster->isSaturatedCluster())
            return;
        if( !settings->diamondPattern.isValidCluster(diamondCluster)){
            //      cerr <<" Cluster is invalid: ";
            //      diamondCluster->Print(1);
            return;
        }
        HitandSeedCount(diamondCluster);
        Int_t clusterSize = diamondCluster->size()-2;
        vecClusterSize.push_back(clusterSize);
    }
    else{
        if(!isTransparentCluster)
            return;
        diamondCluster = transparentCluster;
    }
    //Edge Finding
    //    ShortAnalysis_FillEdgeDistributions(diamondCluster->getPositiveCharge(false));
    Float_t clusterCharge = diamondCluster->getPositiveCharge(useCMN,true);
    FillMeanChargeVector(clusterCharge);
    hTotalAvrgChargeXY->Fill(xPred,yPred,clusterCharge);
    //Universal PHvsChannel Plot
    for(UInt_t i=0; i < settings->diamondPattern.getNIntervals();i++){
        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        //cout<<"Diamond pattern: "<<i<<" Channels: "<<channels.first<<"-"<<channels.second<<endl;
        Int_t HighestSignalChannel = diamondCluster->getHighestSignalChannel();

        if(HighestSignalChannel<=channels.second && HighestSignalChannel>=channels.first){

            if(!bTransAna){
            }
            if(!settings->getSelectionFidCuts()->getFidCut(i+1)->IsInFiducialCut(xFid,yFid))
                return;

            Float_t charge = diamondCluster->getPositiveCharge(false);
            hFidCutXvsFidCutYvsCharge.at(i)->Fill(xFid,yFid,charge);
            hFidCutXvsFidCutYvsEvents.at(i)->Fill(xFid,yFid,1);

            hEventsvsChannel[i]->Fill(diamondCluster->getHighestSignalChannel());
            hPHvsChannel[i]->Fill(charge,diamondCluster->getHighestSignalChannel());
            hLandau[i]->Fill(charge);
            //            vecPHDiamondHit[i]->push_back(charge);
            hFidCutXvsFidCutY[i]->Fill(xFid,yFid);

            if(!bTransAna){
                hHitandSeedCount[i]->Fill(HitCount,SeedCount);
                hChi2XChi2Y[i]->Fill(xChi2, yChi2);
            }
        }
    }       //End of for diamond patterns
}

void TAnalysisOf3DShortAnalysis::Analyse2Cluster(){
    TCluster cluster1 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),0);
    TCluster cluster2 = eventReader->getCluster(TPlaneProperties::getDetDiamond(),1);
    Float_t pos1 = cluster1.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t pos2 = cluster2.getPosition(TCluster::highest2CentroidNoSmallHits);
    Float_t ph1 = cluster1.getPositiveCharge(useCMN);
    Float_t ph2 = cluster2.getPositiveCharge(useCMN);
    Int_t pattern1 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos1));
    Int_t pattern2 = settings->diamondPattern.getPatternOfHit(settings->diamondPattern.convertChannelToMetric(pos2));
    Int_t predictedArea = settings->getSelectionFidCuts()->getFidCutRegion(xFid,yFid);
    Int_t predictedDetector = settings->get3dMetallisationFidCuts()->getFidCutRegion(xPred,yPred);
    if (pattern1 == pattern2)
        hRelatviveNumberOfMultipleClusterEventsSamePattern->Fill(predictedDetector,1);
    if(pattern1<0)pattern1=-2;
    if(pattern2<0)pattern2=-2;
    pattern1++;
    pattern2++;
    if (!settings->diamondPattern.isValidChannelPosition(pos1)){
        Analyse1Cluster(1);//todo how can we use them for long analysis
        return;
    }
    if (!settings->diamondPattern.isValidChannelPosition(pos2)){
        Analyse1Cluster(0);//todo how can we use them for long analysis
        return;
    }
    if (pos1>128||pos1<0||pos2>128||pos2<0){
        return;
    }
    FillMeanChargeVector(ph1+ph2);

    Int_t delta1 = pattern1 - predictedDetector;
    Int_t delta2 = pattern2 - predictedDetector;
    hShortAnalysis2ClusterHitPattern_1stCluster->Fill(predictedDetector,delta1);
    hShortAnalysis2ClusterHitPattern_2ndCluster->Fill(predictedDetector,delta2);

    vecPH_Cluster1.push_back(ph1);
    vecPH_Cluster2.push_back(ph2);
    vecCh_Cluster1.push_back(pos1);
    vecCh_Cluster2.push_back(pos2);
    if (predictedDetector == pattern1||predictedDetector==pattern2){
        if(predictedDetector == 2|| true){
            Double_t relCharge = ph2/ph1;
            hRelativeChargeTwoClustersX->Fill(xPred,relCharge);
            hRelativeChargeTwoClustersY->Fill(yPred,relCharge);
            hRelativeChargeTwoClustersXY->Fill(xPred,yPred,relCharge);
            hShortAnalysis2TotalChargeXY->Fill(xPred,yPred,ph1+ph2);
        }
    }
    hTotalAvrgChargeXY->Fill(xPred,yPred,ph1+ph2);
}

/**
 * Cell Labeling
 *          +-----------+-----------+
 *          +           +           +
 *          +     1     +     3     +       ^
 *          +           +           +     Y |
 *          +-----------+-----------+       |
 *          +           +           +       |
 *          +     0     +     2     +       |
 *          +           +           +       |
 *          +-----------+-----------+       +-------->
 *                                                 X
 * @param xDet
 * @param yDet
 * @return
 */

void TAnalysisOf3DShortAnalysis::FillEdgeDistributions(Float_t clusterCharge){
    for(UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts();i++ ){
        TFiducialCut* fidCut = settings->get3dEdgeFidCuts()->getFidCut(i+1);
        if(!fidCut)
            continue;
        if(fidCut->IsInFiducialCut(xFid,yFid)){
            vecEdgePredX[i].push_back(xPred);
            vecEdgePredY[i].push_back(yPred);
            vecEdgePulseHeight[i].push_back(clusterCharge);
        }
    }
}

void TAnalysisOf3DShortAnalysis::saveHistos(TH1F* hLandauStrip) {
    SaveMeanChargeVector();
    Save2ClusterPlots();
    vector<Float_t> xPred;
    vector<Float_t> yPred;
    vector<Float_t> charge;
    histSaver->SaveHistogram(hRelativeChargeTwoClustersX);
    histSaver->SaveHistogram(hRelativeChargeTwoClustersY);
    histSaver->SaveHistogram(hFidCutsVsMeanCharge);
    TString name = "cRelativeChargeTwoClustersXY";
    name.Append(appendix);
    TCanvas *c1 = new TCanvas(name,name);
    c1->cd();
    hRelativeChargeTwoClustersXY->Draw("colz");
    hRelativeChargeTwoClustersXY->GetZaxis()->SetRangeUser(0,20);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    histSaver->SaveCanvas(c1);
    delete c1;

    name = "cShortAnalysis2TotalChargeXY";
    name.Append(appendix);
    c1 = new TCanvas(name,name);
    c1->cd();
    hShortAnalysis2TotalChargeXY->Draw("colz");
    hShortAnalysis2TotalChargeXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    histSaver->SaveCanvas(c1);
    delete c1;

    histSaver->SaveHistogram(hRelatviveNumberOfMultipleClusterEventsSamePattern);
    histSaver->SaveHistogram(hRelatviveNumberOfMultipleClusterEvents);
    if(hRelatviveNumberOfMultipleClusterEventsSamePattern) delete hRelatviveNumberOfMultipleClusterEventsSamePattern;
    if(hRelatviveNumberOfMultipleClusterEvents) delete hRelatviveNumberOfMultipleClusterEvents;

    name = "cTotalAvrgChargeXY";
    name.Append(appendix);
    if (settings->IsPaperMode())
        gStyle->SetCanvasDefW(gStyle->GetCanvasDefH()*2);
    c1 = new TCanvas(name,name);
    c1->cd();
    c1->SetRightMargin(.2);
    c1->SetObjectStat(false);
    hTotalAvrgChargeXY->Draw("colz");
    //    hTotalAvrgChargeXY->GetZaxis()->SetTitleOffset(1.2);
    Float_t xmin = hTotalAvrgChargeXY->GetXaxis()->GetXmin();
    Float_t xmax = hTotalAvrgChargeXY->GetXaxis()->GetXmax();
    Float_t deltax = xmax - xmin;
    Float_t ymin = hTotalAvrgChargeXY->GetYaxis()->GetXmin();
    Float_t ymax = hTotalAvrgChargeXY->GetYaxis()->GetXmax();
    Float_t deltay = ymax-ymin;
    xmin = xmin -.05*deltax;
    xmax = xmax +.05*deltax;
    ymin = ymin -.05*deltay;
    ymax = ymax +.05*deltay;
    TH1F* histo = c1->DrawFrame(xmin,ymin,xmax,ymax, hTotalAvrgChargeXY->GetTitle());
    histo->GetXaxis()->SetTitle(hTotalAvrgChargeXY->GetXaxis()->GetTitle());
    histo->GetYaxis()->SetTitle(hTotalAvrgChargeXY->GetYaxis()->GetTitle());
    hTotalAvrgChargeXY->SetObjectStat(false);
    hTotalAvrgChargeXY->Draw("colz same");
    hTotalAvrgChargeXY->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    gPad->Update();
    c1->Update();
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1,false);
    TCutG* centralRegion = settings->centralRegion3DnH->GetFiducialAreaCut();
    centralRegion->SetLineColor(kRed);
    centralRegion->Draw("same");
    histSaver->SaveCanvas(c1);

    delete c1;

    if (settings->IsPaperMode())
        gStyle->SetCanvasDefW(gStyle->GetCanvasDefH());

    for(UInt_t i = 0; i < vecEdgePredX.size(); i++){
        xPred.insert(  xPred.end(), vecEdgePredX[i].begin(), vecEdgePredX[i].end());
        yPred.insert(  yPred.end(), vecEdgePredY[i].begin(), vecEdgePredY[i].end());
        charge.insert(charge.end(), vecEdgePulseHeight[i].begin(), vecEdgePulseHeight[i].end());
    }
    //    SaveEdgeFittingDistributions();
    //a.end(), b.begin(), b.end());
    histSaver->SaveHistogram(histSaver->CreateScatterHisto("hEdgeFittingPredictedPosition",yPred,xPred),false);
    TH3F* hEdgeFittingCharge = histSaver->Create3DHisto("hEdgeFittingCharge",xPred,yPred,charge);
    TH2F* hEdgeFittingAvrgCharge = (TH2F*) hEdgeFittingCharge->Project3DProfile("yx");
    hEdgeFittingAvrgCharge->SetName("hEdgeFittingAvrgCharge");
    hEdgeFittingAvrgCharge->SetTitle("hEdgeFittingAvrgCharge");
    histSaver->SaveHistogram(hEdgeFittingAvrgCharge);
    cout<<"vecEdgePredX: "<<vecEdgePredX.size()<<endl;
    if(verbosity%2==1){char t; cin>>t;}
    for(int i = 0; i < vecEdgePredX.size(); i++){
        cout<<"Edge no"<<i<<" "<<vecEdgePredX[i].size()<<" "<<vecEdgePredY[i].size()<<" "<<vecEdgePulseHeight[i].size()<<" "<<endl;
        name = "hEdgeFittingAvrgCharge_";
        name.Append(settings->getEdgePositionName(i));
        TH2F* hEdgeFittingAvrgCharge;
        if(settings->getEdgePositionType(i) == TPlaneProperties::X_COR)
            hEdgeFittingAvrgCharge = histSaver->CreateScatterHisto((string)name,vecEdgePulseHeight[i],vecEdgePredX[i],200);
        else
            hEdgeFittingAvrgCharge = histSaver->CreateScatterHisto((string)name,vecEdgePulseHeight[i],vecEdgePredY[i],200);

        hEdgeFittingAvrgCharge->GetYaxis()->SetTitle("Pulse Height /ADC");
        TString title = "predicted Position ";
        title.Append(TPlaneProperties::getCoordinateString(settings->getEdgePositionType(i)).c_str());
        title.Append(" / #mum");
        hEdgeFittingAvrgCharge->GetXaxis()->SetTitle(title);//"predicted Position X / #mum");
        histSaver->SaveHistogram(hEdgeFittingAvrgCharge);
        TH1F* hEdgeFittingAvrgCharge_pfx = (TH1F*)hEdgeFittingAvrgCharge->ProfileX();
        if(hEdgeFittingAvrgCharge_pfx){
            hEdgeFittingAvrgCharge_pfx->GetYaxis()->SetTitle("avrg. Charge / ADC");
            TCutG *cut = this->settings->getEdgePosition(i);
            name = "c";
            name.Append(hEdgeFittingAvrgCharge_pfx->GetName());
            TCanvas *c1 = new TCanvas(name,name);
            hEdgeFittingAvrgCharge_pfx->Draw();
            if(cut)cut->Draw();
            histSaver->SaveCanvas(c1);
            cout<<"Saved "<<c1->GetName()<<endl;
            delete c1;
        }
        else
            cout<<" Cannot create ProfileX" << endl;
        TH1D* histo1st_py;
        TH1D* histo2nd_py;
        histSaver->SaveHistogram(hShortAnalysis2ClusterHitPattern_1stCluster);
        histSaver->SaveHistogram(hShortAnalysis2ClusterHitPattern_2ndCluster);
        for(int i = 0; i <=hShortAnalysis2ClusterHitPattern_1stCluster->GetNbinsX();i++){
            TString extension = TString::Format("_pattern%d",i);
            if (i==0)
                extension = "_all";
            name = hShortAnalysis2ClusterHitPattern_1stCluster->GetName();
            name.Append(extension);
            if(verbosity>3) cout<<name<<endl;
            if(i==0)
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name);
            else{
                //              int bin = hShortAnalysis2ClusterHitPattern_1stCluster->GetYaxis()
                histo1st_py= hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionY(name,i,i);
            }
            name = hShortAnalysis2ClusterHitPattern_2ndCluster->GetName();
            name.Append(extension);
            if(verbosity>3)cout<<name<<endl;

            if(i==0)
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name);
            else
                histo2nd_py = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionY(name,i,i);
            name = "h2ClusterAnalysis_ClusterPatterns";
            name.Append(extension);
            histSaver->SaveTwoHistos((string)name,histo1st_py,histo2nd_py);
            histSaver->SaveHistogram(histo1st_py);
            histSaver->SaveHistogram(histo2nd_py);
            delete histo1st_py;
            delete histo2nd_py;
        }
        name = "c2ClusterAnalysis_ClusterPatterns";
        name.Append(appendix);
        TH1D* histo1st_px = hShortAnalysis2ClusterHitPattern_1stCluster->ProjectionX();//name);
        TH1D* histo2nd_px = hShortAnalysis2ClusterHitPattern_2ndCluster->ProjectionX();//name);
        histSaver->SaveTwoHistos((string)name,histo1st_px,histo2nd_px);
        histSaver->SaveHistogram(histo1st_px);
        histSaver->SaveHistogram(histo2nd_px);
        //      if(histo1st_px) delete histo1st_px;
        //      if(histo2nd_px) delete histo2nd_px;
    }

    //  char t; cin>>t;
    //hNumberofClusters
    histSaver->SaveHistogram(hNumberofClusters);
    for(UInt_t i=0;i<settings->diamondPattern.getNIntervals();i++){
        hEventsvsChannelCombined->Add(hEventsvsChannel.at(i));
    }
    histSaver->SaveHistogram(hEventsvsChannelCombined);
    //  vector<TH1*> hLandauSorted;

    for(UInt_t i=0;i<settings->diamondPattern.getNIntervals();i++){

        pair<int,int> channels = settings->diamondPattern.getPatternChannels(i+1);
        name = "c_";
        name.Append(hLandau[i]->GetName());

        Float_t factor = hLandau[i]->GetBinContent(hLandau[i]->GetMaximumBin());
        factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
        name.Append("_normalized");
        histSaver->SaveTwoHistosNormalized((string)name,hLandau[i],hLandauStrip);

        Float_t max = hHitandSeedCount[i]->GetBinContent(hHitandSeedCount[i]->GetMaximumBin());
        hHitandSeedCount[i]->Scale(1./max);
        histSaver->SaveHistogram(hHitandSeedCount[i]);

        name = "c_"+(TString)hPHvsChannel[i]->GetName();
        TCanvas *c1 = new TCanvas(name,name);
        name = "h"+(TString)hPHvsChannel[i]->GetName();

        Int_t min = channels.first-1;
        max = channels.second+1;
        Int_t bins = max-min;
        TH2F* histo = new TH2F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax,bins,min,max);
        c1->cd();
        histo->Draw();
        hPHvsChannel[i]->Draw("goff");
        hPHvsChannel[i]->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->GetXaxis()->SetLimits(PulseHeightMin,PulseHeightMax);
        //hPHvsChannel[i]->GetXaxis()->SetRange(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->GetXaxis()->SetRangeUser(PulseHeightMin,PulseHeightMax);
        hPHvsChannel[i]->Draw("colzsame");
        histSaver->SaveCanvas(c1);
        histSaver->SaveHistogram(hPHvsChannel[i]);
        //histSaver->SaveHistogram(hPHvsPredictedXPos.at(i));
        Float_t maxChi2 = 12;
        hChi2XChi2Y[i]->Draw("colz");
        hChi2XChi2Y[i]->GetXaxis()->SetRangeUser(0,maxChi2);
        hChi2XChi2Y[i]->GetYaxis()->SetRangeUser(0,maxChi2);
        histSaver->SaveHistogram(hChi2XChi2Y[i]);
        histSaver->SaveHistogram(hFidCutXvsFidCutY.at(i));
        //histSaver->SaveHistogram(hPHvsPredictedChannel.at(i));
        //histSaver->SaveHistogram(hFidCutXvsFidCutYvsCharge.at(i));

        //hFidCutXvsFidCutYvsMeanCharge
        //           ptrCanvasMean.at(i)->cd();
        TCanvas *cc = new TCanvas();
        *hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
        hFidCutXvsFidCutYvsMeanCharge.at(i)->SetEntries(hFidCutXvsFidCutYvsEvents.at(i)->Integral());
        hFidCutXvsFidCutYvsMeanCharge.at(i)->Draw("COLZ");
        hFidCutXvsFidCutYvsMeanCharge.at(i)->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
        TString hName  = TString::Format("cFidCutXvsFidCutYvsMeanCharge_%d_%d",channels.first,channels.second);
        hName.Append(appendix);
        cc->SetName(hName);
        histSaver->SaveCanvas(cc);
        delete cc;

        /*//hXdetvsYdetvsEvents
           ptrCanvasXdetvsYdetMeanCharge.push_back(new TCanvas());
           ptrCanvasXdetvsYdetMeanCharge.at(i)->cd();
           hFidCutXvsFidCutYvsMeanCharge.at(i) = (*hFidCutXvsFidCutYvsCharge.at(i)/(*hFidCutXvsFidCutYvsEvents.at(i)));
           hXdetvsYdetvsEvents.at(i)->SetEntries(hXdetvsYdetvsEvents.at(i)->Integral());
           hXdetvsYdetvsEvents.at(i)->Draw("COLZ");
           hName  = TString::Format("cXdetvsYdetMeanCharge_%d_%d",channels.first,channels.second);
           ptrCanvasXdetvsYdetMeanCharge.at(i)->SetName(hName);
           histSaver->SaveCanvas(ptrCanvasXdetvsYdetMeanCharge[i]);
         */

    } //End of for loop
    TCanvas* cCombinedMeanCharge = new TCanvas();
    cCombinedMeanCharge->cd();

    // no Fiducial Cuts Drawn
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(0));
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(1));
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Add(hFidCutXvsFidCutYvsMeanCharge.at(2));
    name = "hFidCutXvsFidCutYvsMeanChargeAllDetectorsNoFidDrawn";
    name.Append(appendix);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    // Selection Fiducial Cuts Drawn
    name = "hFidCutXvsFidCutYvsMeanChargeAllDetectors";
    name.Append(appendix);
    settings->getSelectionFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    // Edge F
    cCombinedMeanCharge->Clear();
    name = "hFidCutXvsFidCutYvsMeanChargeAllEdges";
    name.Append(appendix);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->SetTitle(name);
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
    hFidCutXvsFidCutYvsMeanChargeAllDetectors->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);

    settings->get3dEdgeFidCuts()->DrawFiducialCutsToCanvas(cCombinedMeanCharge,true);
    cCombinedMeanCharge->SetName(name);
    histSaver->SaveCanvas(cCombinedMeanCharge);

    for ( UInt_t i = 0; i < settings->get3dEdgeFidCuts()->getNFidCuts(); i++){
        cCombinedMeanCharge->Clear();
        hFidCutXvsFidCutYvsMeanChargeAllDetectors->Draw("COLZ");
        settings->get3dEdgeFidCuts()->getFidCut(i+1)->DrawFiducialCutToCanvas(cCombinedMeanCharge,true);
        TString name = "hFidCutXvsFidCutYvsMeanCharge_";
        name.Append(settings->getEdgePositionName(i));
        name.Append(appendix);
        cCombinedMeanCharge->SetName(name);
        histSaver->SaveCanvas(cCombinedMeanCharge);
    }

    for( UInt_t i=0; i < hFidCutXvsFidCutYClusters.size(); i++){
        if (hFidCutXvsFidCutYClusters[i])
            hFidCutXvsFidCutYClusters[i]->SetName(TString::Format("hFidCutXvsFidCutYClusters_%d",i));
        histSaver->SaveHistogram((TH2F*)hFidCutXvsFidCutYClusters[i]);
    }


    //TODO FIX ME
    //       if(settings->do3dShortAnalysis()){
    //           TString name = TString::Format("hTransparentAnalysisTransparentChargeWithoutEdgeBadCellsComparison_DiamondPattern2_to_ClusterSize1");
    //           hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]->SetLineColor(kRed);
    //           histSaver->SaveTwoHistos((string)name,hLandau[1],hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]);
    //
    //           name = TString::Format("hTransparentAnalysisTransparentChargeWithoutEdgeBadCellsComparison_DiamondPattern2_to_ClusterSize1_normalized");
    //           histSaver->SaveTwoHistosNormalized((string)name,hLandau[1],hTransparentAnalysisTransparentChargeBadCellsWithoutEdge[0]);
    //       }

}

void TAnalysisOf3DShortAnalysis::SaveMeanChargeVector() {
    cout<<"SaveMeanChargeVector"<<endl;
    cout<<"vecPredDetX_ShortAna    "<<vecPredDetX.size()<<endl;
    cout<<"vecPredDetY_ShortAna    "<<vecPredDetY.size()<<endl;
    cout<<"vecPulseHeight_ShortAna "<<vecPulseHeight.size()<<endl;

    if (bTransAna)
        cout<<"Create Project3dProfile for hChargeDistribution3D ..."<<flush;

    TProfile2D* hMeanCharge = histSaver->CreateProfile2D("hChargeDistribution3D",vecPredDetX,vecPredDetY,vecPulseHeight,1024,1024);
    if (!hMeanCharge){
        cerr<<" hChargDistribution3D: was not created:"<<endl;
        return;
    }
    else if (hMeanCharge->GetEntries() == 0){
        cerr<<" hChargDistribution3D: number of entries is 0"<<endl;
        for (int i = 0; i<vecPredDetX.size()&&i<100;i++)
            cout<<TString::Format("%3d %5.1f/%5.1f --> %6.1f",i,vecPredDetX.at(i),vecPredDetY.at(i),vecPulseHeight.at(i))<<endl;
        return;

    }
    hMeanCharge->GetXaxis()->SetTitle("#it{X} / #mum");
    hMeanCharge->GetYaxis()->SetTitle("#it{Y}/ #mum");
    hMeanCharge->GetYaxis()->SetTitleOffset(1.4);
    hMeanCharge->GetZaxis()->SetTitleOffset(1.3);
    cout<<"\t[done]"<<endl;
    //  if(!hMeanCharge)
    //      return;
    //  else
    //      cerr<<" hChargDistribution3D_pfyx: was not created:"<<endl;
    hMeanCharge->GetZaxis()->SetTitle("Avrg. pulse height /ADC");
    hMeanCharge->SetTitle("Avrg. pulse height in detector system");
    hMeanCharge->GetZaxis()->SetRangeUser(PulseHeightMinMeanCharge,PulseHeightMaxMeanCharge);
    TString name = "hAvrgPulseHeigthDetSystem";
    name.Append(appendix);
    hMeanCharge->SetName(name);
    histSaver->SaveHistogram(hMeanCharge,false);

    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer";
    name.Append(appendix);
    TCanvas *c1 = new TCanvas(name, name);
    c1->cd();
    c1->SetRightMargin(.15);
    hMeanCharge->Draw("colz");
    c1->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hMeanCharge->GetListOfFunctions()->FindObject("palette");
    //    palette->SetY2NDC(0.7);
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
    settings->DrawMetallisationGrid(c1,3);
    histSaver->SaveCanvas(c1);
    c1->Clear();

    TH2D* hMeanChargeEntries = hMeanCharge->ProjectionXY(name,"B");
    hMeanChargeEntries->Draw("colz");
    settings->get3dMetallisationFidCuts()->DrawFiducialCutsToCanvas(c1);
    settings->DrawMetallisationGrid(c1,3);
    name = "cAvrgPulseHeigthDetSystemEntries_MetalizationLayer";
    name.Append(appendix);
    c1->SetName(name);
    if (hMeanChargeEntries)
        delete hMeanChargeEntries;

    TFiducialCut *fidCut3dWithColumns = settings->get3dMetallisationFidCuts()->getFidCut(3);
    cout<<"3d FidCut: "<<endl;
    fidCut3dWithColumns->Print(1);
    cout<<endl;
    Float_t xmin = fidCut3dWithColumns->GetXLow();
    Float_t xmax = fidCut3dWithColumns->GetXHigh();
    Float_t deltaX = TMath::Abs(.05*(xmax-xmin));
    Float_t ymin = fidCut3dWithColumns->GetYLow();
    Float_t ymax = fidCut3dWithColumns->GetYHigh();
    Float_t deltaY = TMath::Abs(.05*(ymax-ymin));

    hMeanCharge->GetXaxis()->SetRangeUser(xmin-deltaX,xmax+deltaX);
    hMeanCharge->GetYaxis()->SetRangeUser(ymin-deltaY,ymax+deltaY);
    hMeanCharge->Draw("colz");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom";
    name.Append(appendix);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom_rebinned";
    TProfile2D* hMeanCharge3D = histSaver->CreateProfile2D("hChargeDistribution3D_3D",
            vecPredDetX,vecPredDetY,vecPulseHeight,
            settings->getNColumns3d()*4,settings->getNRows3d()*4,
            xmin,xmax,ymin,ymax
    );
    hMeanCharge3D->Draw("colz");
    c1->Update();
    name = "cAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom_rebinned";
    name.Append(appendix);
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    name = "hEntriesAvrgPulseHeigthDetSystem_MetalizationLayer_Zoom_rebinned";
    name.Append(appendix);
    hMeanChargeEntries = hMeanCharge3D->ProjectionXY(name,"B");
    hMeanChargeEntries->Draw("colz");
    name.Replace(0,1,"c");
    c1->SetName(name);
    histSaver->SaveCanvas(c1);

    if(hMeanChargeEntries)
        delete hMeanChargeEntries;
    if (hMeanCharge)
        delete hMeanCharge;
    if(hMeanCharge3D)
        delete hMeanCharge3D;
}

void TAnalysisOf3DShortAnalysis::SaveEdgeDistributions() {
}

void TAnalysisOf3DShortAnalysis::Save2ClusterPlots() {
    TH2F * hPH = histSaver->CreateScatterHisto("hPulseHeightComparision2Clusters",vecPH_Cluster2,vecPH_Cluster1,
            PulseHeightBins,PulseHeightBins,PulseHeightMin,PulseHeightMax,PulseHeightMin-1000,PulseHeightMax-1000);
    hPH->GetXaxis()->SetTitle("Pulse height cluster no. 1");
    hPH->GetYaxis()->SetTitle("Pulse height cluster no. 2");
    histSaver->SaveHistogram(hPH);
    delete hPH;

    TH2F* hCh = histSaver->CreateScatterHisto("hChannelComparision2Clusters",vecCh_Cluster2,vecCh_Cluster1,128,128,0,128,0,128);
    hCh->GetXaxis()->SetTitle("Channel no for cluster no. 1");
    hCh->GetYaxis()->SetTitle("Channel no for cluster no. 2");
    hCh->Draw("colz");
    //  Float_t xmax = hCh->GetZaxis()->GetXmax();
    //  Float_t xmin = hCh->GetZaxis()->GetXmin();
    //  hCh->GetZaxis()->SetRangeUser(xmin+.1*(xmax-xmin),xmax);
    histSaver->SaveHistogram(hCh,false);
    delete hCh;
}
