/**
 * @file TAnalysisOfAsymmetricEta.cpp
 *
 * @date May 13, 2013
 * @author bachmair
 * @description
 */

#include "../include/TAnalysisOfAsymmetricEta.hh"

TAnalysisOfAsymmetricEta::TAnalysisOfAsymmetricEta(TSettings *settings) {
	// TODO Auto-generated constructor stub
	if (settings!=0)
		this->settings=settings;
	else
		exit(-1);

	alpha=0;
	histSaver=new HistogrammSaver(settings);
	histSaver->SetVerbosity(0);
	maxTriesPeakfinding = 30;
	maxTriesAlpha = 120;
	settings->goToClusterAnalysisDir();
	histSaver->SetPlotsPath(settings->getClusterAnalysisDir().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	verbosity = settings->getVerbosity();

	hAsymmetricEta2D=0;
	bestAlpha = 0;
}

TAnalysisOfAsymmetricEta::~TAnalysisOfAsymmetricEta() {
	// TODO Auto-generated destructor stub
}


void TAnalysisOfAsymmetricEta::FillEtaDistribution(TH2F* histo, Float_t correctionAlpha, TString hName){
	if (verbosity>4){
		cout<<" TAnalysisOfAsymmetricEta::FillEtaDistribution"<<settings->getNDiamonds()<<" "
			<<settings->getNDiaDetectorAreas()<<endl;
		char t;
		if (verbosity%2==1 && det >=8	){
			cout<<"Press a key"<<endl;
			cin>>t;
		}
	}
	//Fill histo with Asymmetric Cluster data..
	hAsymmetricEta2D->Reset();
	hAsymmetricEta2D->Clear();
	TString name = hName;
	name.Append(TString::Format("_%05d",(int)(TMath::Abs(alpha)*10000)));
	hAsymmetricEta2D->SetName(name);
	name = hName;
	name.Append(TString::Format("_%03.2f",alpha));
	hAsymmetricEta2D->SetTitle(name);
	if(verbosity>4)
		cout<<" Fill " << histo->GetName() <<" with " << clusters.size() << "Entries."<< endl;
	Int_t toSmallEta = 0;
	Int_t toBigEta = 0;
	Int_t isMasked = 0;
	Int_t toBigCluster = 0;
	Int_t correct = 0;
	Int_t total = 0;
	for (UInt_t i=0;i<clusters.size();i++){
		total++;
		TCluster  clus = clusters.at(i);
		if(clus.getClusterSize()>=settings->getMaxAllowedClustersize(det)){
			toBigCluster++;
			continue;
		}
		if (settings->isMaskedCluster(det,clus)){
			isMasked ++;
		    continue;
		}
		TCluster newClus = clus.getCrossTalkCorrectedCluster(correctionAlpha);
		Int_t leftChannel,oldChannel;
		Float_t newEta = newClus.getEta(leftChannel);
		Int_t nDia = settings->getDiaDetectorAreaOfChannel(leftChannel)+1;
		Float_t oldEta = clus.getEta(oldChannel);
		Int_t oldDia = settings->getDiaDetectorAreaOfChannel(oldChannel)+1;
//		if (det >=8 && leftChannel!=oldChannel){
		if (det >=8 && false)
		    if (newClus.getHighestChannelNumber() == clus.getHighestChannelNumber()){}
		    else if(nDia==oldDia){}
		    else{
		    cout<<"\n\n"<<i<<" nDia: "<<nDia<<"/"<<oldDia<<" Channel: "<<leftChannel<<"/"<<oldChannel<<" eta: "<<newEta<<"/"<<oldEta<<endl;
		    newClus.Print();
		    clus.Print();
		}
//        if (verbosity>4)
		if(newEta>0&&newEta<1){
			histo->Fill(newEta,nDia);
			correct++;
			if (det >= 8 && nDia==0){
			    cout<<"\n\n"<<i<<" nDia: "<<nDia<<"/"<<oldDia<<" Channel: "<<leftChannel<<"/"<<oldChannel<<" eta: "<<newEta<<"/"<<oldEta<<endl;
			               newClus.Print();
			}
			//				cout<<i<<" nDia: "<<nDia<<" eta: "<<newEta<<endl;
		}
		else{
			if (newEta<=0)
				toSmallEta++;
			else
				toBigEta++;
		}
	}
	if (verbosity >4){
		cout<<"histo: "<<histo->GetName()<<" "<<histo->GetEntries()<<endl;
		cout<<"CutFlow"<<endl;
		cout<<"\tCluster Size wrong: "<<toBigEta<<endl;
		cout<<"\tIs Masked Cluster:  "<<isMasked<<endl;
		cout<<"\tto Small Eta:       "<<toSmallEta<<endl;
		cout<<"\tto Big   Eta:       "<<toBigEta<<endl;
		cout<<"\tCorrect Size:       "<<correct<<endl;

	}
}

void TAnalysisOfAsymmetricEta::setClusters(vector<TCluster> clusters) {
	for(UInt_t i = 0; i< clusters.size();i++){
		this->clusters.push_back(clusters.at(i));
	}
}

TH1F* TAnalysisOfAsymmetricEta::getProjection(){
	//	cout<<"[TAnalysisOfAsymmetricEta::getProjection] "<<hAsymmetricEta2D<<endl;
	if(hAsymmetricEta2D==0){
		cout<<"Cannot get Projection of invalid histo."<<endl;
		exit(-1);
	}
//	cout<<"hAsymmetricEta2D Entries: "<<hAsymmetricEta2D->GetEntries()<<endl;
	//	cout<<"GetProjection"<<endl;
	if(settings->getVerbosity()>4)cout<<"GetProjection"<<endl;
	TH1F* hProjection = 0;
	TString hName = "hAsymmetricEta_";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	hName.Append("_");
	TString hTitle = "AsymmetricEta, ";
	hTitle.Append(TPlaneProperties::getStringForDetector(det).c_str());
	hTitle.Append(", #alpha = ");
	if( TPlaneProperties::isDiamondDetector(det)){
		if(verbosity>4)cout<<"Diamond Detector"<<endl;
		int nDias = settings->getNDiamonds();
		int nDia = settings->getAnalysedDiamond();
		if(settings->getVerbosity()>4)cout<<"nDias: "<<nDias<<endl;
		if(settings->getVerbosity()>4)cout<<"nDia: "<<nDia<<endl;
		int bin;
		if(settings->getNDiamonds()!=1){
			if (settings->getAnalysedDiamond()==TSettings::leftDia){//todo
				if (verbosity>4)cout<<"LEFT DIAMOND ANALYSED..."<<1.0<<endl;
				bin = hAsymmetricEta2D->GetYaxis()->FindBin(1.0);
				if(settings->getVerbosity()>4)cout<<"bin:"<<bin<<endl;
				TString name = hName;//+"_left";//TString::Format("%s_left",hName);
				name.Append(TString::Format("%06d_left",(int)(TMath::Abs(alpha)*100000)));
				if(settings->getVerbosity()>4)cout<<"Title: "<<flush;
				hTitle.Append(TString::Format("%02.2f %%, left",alpha*100));
				if(settings->getVerbosity()>4)cout<<hTitle<<endl;
				if(settings->getVerbosity()>4)cout<<name<<endl;
				hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(name,bin,bin);
				if(settings->getVerbosity()>4)cout<<"lDia:"<<hProjection<<endl;
			}
			else if(settings->getAnalysedDiamond()==TSettings::rightDia){//todo
				if (verbosity>4)cout<<"RIGHT DIAMOND ANALYSED..."<<2.0<<endl;
				bin = hAsymmetricEta2D->GetYaxis()->FindBin(2.0);
				if(settings->getVerbosity()>4)cout<<"bin:"<<bin<<endl;
				TString name = hName;//+"_left";//TString::Format("%s_left",hName);
				name.Append(TString::Format("%06d_right",(int)(TMath::Abs(alpha)*100000)));
				if(settings->getVerbosity()>4)cout<<"Title: "<<flush;
				hTitle.Append(TString::Format("%02.2f %%, right",alpha*100));
				if(settings->getVerbosity()>4)cout<<hTitle<<endl;
				if(settings->getVerbosity()>4)cout<<name<<endl;
				hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(name,bin,bin);
				if(alpha==0)
					histSaver->SaveHistogram(hProjection);
				if(settings->getVerbosity()>4)cout<<"rDia:"<<hProjection<<endl;
			}
			else {//todo
				if (verbosity>4) cout<<"NONE OF TWO DIAMONDs ANALYSED: "<<1.0<<endl;
				bin = hAsymmetricEta2D->GetYaxis()->FindBin(1.0);
				if(settings->getVerbosity()>4)cout<<"bin:"<<bin<<endl;
				TString name = hName;//+"_left";//TString::Format("%s_left",hName);
				name.Append(TString::Format("%06d_all",(int)(TMath::Abs(alpha)*100000)));
				if(settings->getVerbosity()>4)cout<<"Title: "<<flush;
				hTitle.Append(TString::Format("%02.2f %%, all",alpha*100));
				if(settings->getVerbosity()>4)cout<<hTitle<<endl;
				if(settings->getVerbosity()>4)cout<<name<<endl;
				hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(name);//,bin,bin);
				if(settings->getVerbosity()>4)cout<<"nDia:"<<hProjection<<endl;

			}
		}
		else{// has only one diamond
			if (verbosity>4) cout<<"ONLY ONE DIAMOND - TAKE ALL"<<endl;
			hName.Append("all_px");
			hTitle.Append(TString::Format("%02.2f %%, all",alpha*100));
			hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(hName);
			if(settings->getVerbosity()>4)cout<<"1Dia"<<hProjection<<endl;
		}
		if (verbosity%2==1&&verbosity>4){
			cout<<"Press a key to continue"<<endl;
			char t;
			cin>>t;
		}

	}
	else{// is Silicon Detector
		if(settings->getVerbosity()>4)cout<<"Silicon Detector "<<endl;
		hName.Append("_px");
		hTitle.Append(TString::Format("%02.2f %%",alpha*100));
		hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(hName);
		if(settings->getVerbosity()>4)cout<<"Sil: "<<hProjection<<endl;
	}
	if(hProjection)
		hProjection->SetTitle(hTitle);
	else
		if(settings->getVerbosity()>4)cout<<"hProjection is not Valid."<<endl;
	if(settings->getVerbosity()>4)cout<<"returning "<<hProjection<<": "<<hProjection->GetTitle()<<" "<<hProjection->GetEntries()<<endl;;
	return hProjection;
}


UInt_t TAnalysisOfAsymmetricEta::analyse() {
	cout<<"ANALYSE ASYMMETRIC ETA SAMPLE - CROSS TALK CORRECTION - Detector "<<this->det<<endl;
	//	cout<<"Press a key..."<<flush;
	//	char t; cin>>t;

	TString hName ="hAsymmetricEtaCorrected";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	int nDiamonds = settings->getNDiamonds();
	cout<<"nDiamonds = "<<nDiamonds<<endl;
	if (det>=8)
		settings->diamondPattern.showPatterns();
	//	hAsymmetricEta2D = new TH2F(hName,hName,512,0,1,nDiamonds+3,-1.5,nDiamonds+1.5);
	//	TH1F* hAsymmetricEta = new TH1F(hName,hName,512,0,1);
	alpha= 0;
	UInt_t nTries =0;

	TH1F* hAsymmetricEta_px=0;
	bool valid=false;
	Reset();
	TH2F* hAsymmetricEta2DNoCorrection=0;
	if(hAsymmetricEta2D){
		delete hAsymmetricEta2D;
		hAsymmetricEta2D=0;
	}
	if(verbosity>4)cout<<"hName"<<endl;
	hAsymmetricEta2D = new TH2F(hName,hName,512,0,1,nDiamonds+2,-0.5,nDiamonds+1.5);
	if(hAsymmetricEta2D ==0){
		cout<<"histo with name: '"<<hName<<"' was not created...."<<endl;//hAsymmetricEta2D<<
	}
	while (!valid && nTries < maxTriesAlpha){
		if(verbosity>4)cout<<"\n\nnew try no. "<<nTries << " / " << maxTriesAlpha <<" "<<hAsymmetricEta2D << " "<<hAsymmetricEta2D->GetName();

		FillEtaDistribution(hAsymmetricEta2D,alpha,hName);

		if(alpha==0) {
			hAsymmetricEta2DNoCorrection = (TH2F*)hAsymmetricEta2D->Clone();
			saveAsymmetricEtaPerArea(hAsymmetricEta2D,hAsymmetricEta2DNoCorrection,"hAsymmetricEta_",alpha);
		}

		if(hAsymmetricEta_px) {
			delete hAsymmetricEta_px;
			hAsymmetricEta_px=0;
		}
		if(verbosity>4)cout<<hAsymmetricEta2D<<endl;
		hAsymmetricEta_px = getProjection();
		Float_t skewness  = checkConvergence(hAsymmetricEta_px,nTries);
		Float_t mean = hAsymmetricEta_px->GetMean();
		valid = (skewness==0);
		if(!valid ) updateAlpha(skewness,mean);

		nTries++;
	}
	checkBadConvergence(nTries);
	createConvergencePlot();
//	cout<<" BEST ALPHA: "<<bestAlpha<<endl;
	if (nTries==maxTriesAlpha){
		cout<<TString::Format("Not converged,best Try %d, best alpha: %02.2f %%",bestTry,bestAlpha *100)<<endl;
		hName= "hAsymmetricEtaBestAlpha_";
		hName.Append(TPlaneProperties::getDetectorNameString(det).c_str());
		FillEtaDistribution(hAsymmetricEta2D,bestAlpha,hName);
		if(hAsymmetricEta_px){
			delete hAsymmetricEta_px;
			hAsymmetricEta_px=0;
		}
		hAsymmetricEta_px = getProjection();
	}
	else
		cout<<TString::Format("Convergered after %d/%d, best alpha: %02.2f %%",nTries,maxTriesAlpha,bestAlpha *100)<<endl;
	hName ="hAsymmetricEta_";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	hName.Append(TString::Format("%05d",(int)(TMath::Abs(bestAlpha)*10000)));
	hAsymmetricEta_px->SetName(hName);
	hName ="hAsymmetricEta ";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	hName.Append(TString::Format(" %2.2f %%",bestAlpha*100));
	hAsymmetricEta_px->SetTitle(hName);

	//	hAsymmetricEta->Smooth(3);
	TCanvas c1;
	c1.cd();
	hAsymmetricEta_px->Draw("");
	c1.Update();
	TPaveStats* ps = (TPaveStats *)hAsymmetricEta_px->GetListOfFunctions()->FindObject("stats");
	c1.Update();
	if(ps){
		ps->SetX1NDC(0.35);
		ps->SetX2NDC(0.55);
		ps->SetY1NDC(.5);
		ps->SetY2NDC(.8);
	}
	else
		cout<<"cannot find ps"<<endl;
	histSaver->SaveHistogram(hAsymmetricEta_px,false,false,true);;
	histSaver->SaveHistogram(hAsymmetricEta_px,false,false,true);;
	if(verbosity>4)cout<<"nDiamonds: "<<settings->getNDiamonds()<<endl;
	hName = "hAsymmetricEtaFinal_";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	saveAsymmetricEtaPerArea(hAsymmetricEta2D,hAsymmetricEta2DNoCorrection,hName,bestAlpha);
	TString name = "hDiamondHitNo";
	histSaver->SaveHistogram((TH1F*)hAsymmetricEta2D->ProjectionY(name));
	//	histSaver->SaveHistogram(hAsymmetricEta);
	cout<<flush;
	//	char t;
	//	cin>>t;
	return nTries;
}

Float_t TAnalysisOfAsymmetricEta::checkConvergence(TH1F* histo, UInt_t nTries){
	if (verbosity>4){
		cout<<"TAnalysisOfAsymmetricEta::checkConvergence: nTries "<<nTries<<"\t"<<histo->GetName()<<endl;
		cout<<"entries: "<<histo->GetEntries()<<endl;
	}
	int leftHalf = 0;
	int bin=0;
	for (bin=0; histo->GetBinLowEdge(bin)<.5;bin++){
		leftHalf+=histo->GetBinContent(bin);
	}
	int rightHalf = 0;
	for (bin=bin; bin<histo->GetNbinsX();bin++){
		rightHalf+=histo->GetBinContent(bin);
	}
	Float_t mean = histo->GetMean();
	if(verbosity>4)cout<<"Analysising: "<<histo->GetName()<<endl;
	if(verbosity>4)cout<<"Asymmetric eta with a charge share of "<<alpha*100<<"%"<<endl;
	if(verbosity>4)cout<<"nentries: "<<histo->GetEntries()<<endl;
	Float_t skewness = histo->GetSkewness();
	if(verbosity>4)cout<<"skewness*100: "<<skewness*100<<endl;
	if(verbosity>4)cout<<"left: "<<leftHalf<<"\tright: "<<rightHalf<<" "<<(Float_t)leftHalf/(Float_t)rightHalf*100<<endl;

	bool valid = TMath::Abs(mean-.5)<.005;
	valid = valid || TMath::Abs(skewness)<0.005;
	if(verbosity>4)cout<<"Mean "<<mean*100<<" "<<valid<<endl;
	TPolyMarker* pm = FindPeaks(histo,2);
	Float_t p1=0;
	Float_t p2=0;
	if(pm){
		for (int i=0;i< pm->GetN();i++){
			Float_t peakPos = pm->GetX()[i];
			if (peakPos>.5)
				peakPos = TMath::Abs(peakPos-1);
			if(i==0)
				p1 = peakPos;
			else if(i==1)
				p2 = peakPos;
			if(verbosity>4)cout<<"\tPeak - "<<i<<"\t"<<peakPos*100.<<": "<<pm->GetY()[i]<<"\n";
		}
		if (pm->GetN()==2){
			if(verbosity>4)cout<< (TMath::Abs(p1/p2-1)*100);
			valid = valid && TMath::Abs(p1/p2-1)*100<8;
			if(verbosity>4)cout<<" "<<valid;
		}
		else
			valid=false;
		if(verbosity>4)cout<<" "<<valid<<"\n";
	}
	else
		valid = false;
	valid = TMath::Abs(mean-.5)<.001;
	valid = valid || TMath::Abs(skewness)<0.0001;
	alphaValues.push_back(alpha*10);
	Float_t value = (Float_t)leftHalf/(Float_t)rightHalf-1;
	if(value!=value) value = -1;
	rightLefts.push_back(value);
	if(skewness!=skewness)skewness= -1;
	skewnesses.push_back(skewness);
	if(mean!=mean) mean = -1.5;
	means.push_back(mean-.5);
	vecTries.push_back(nTries);
	if (det>=8){
		cout<<"\tAlpha: "<<alphaValues.back()<<endl;
		cout<<"\tMean:  "<<means.back()<<endl;
		cout<<"\tSkew:  "<<skewnesses.back()<<endl;
		cout<<"\tL/RA:  "<<rightLefts.back()<<endl;
	}
	if(valid)
		skewness=0;
	return skewness;
}

bool TAnalysisOfAsymmetricEta::updateAlpha(Float_t skewness, Float_t mean){
	if(verbosity>4||true){
		cout<<"TAnalysisOfAsymmetricEta::updateAlpha for "<<det<<"\tMean:"<<mean<<" Skew:"<<skewness<<endl;
		cout<<"Alpha: "<<alpha<<endl;
	}
	if (alpha == 0){
			if(mean- .5>0)
				alpha = -.01;
			else
				alpha = +.01;
	}
	else {
		if(det == 6 || det ==2) 
			skewness*=-1;
		//				if( TMath::Abs(mean-.5)>.05){
			//					if(alpha<0)
		//						alpha *= .5/mean*2;
		//					else
		//						alpha *= mean/.5/2;
		//				}
		//				else
		if (skewness < 0 ){
			if(alpha<0)
				alpha*=1.01;
			else
				alpha *=.95;
		}
		else{//skewness<0
			if (alpha<0)
				alpha*=.95;
			else
				alpha*=1.01;
		}
	}
	if (verbosity>4||true)
		cout<<"Alpha: "<<alpha<<endl;
	return true;
}
void TAnalysisOfAsymmetricEta::createConvergencePlot(){
		TGraph *gr1 = new TGraph(vecTries.size(),&vecTries[0],&alphaValues[0]);
		gr1->SetTitle("#alpha Values (10%)");
		gr1->SetLineColor(kBlack);
		gr1->SetMarkerStyle(20);
		gr1->SetFillColor(kWhite);
		TGraph *gr2 = new TGraph(vecTries.size(),&vecTries[0],&means.at(0));
		gr2->SetTitle("mean value");
		gr2->SetLineColor(kRed);
		gr2->SetMarkerStyle(21);
		gr2->SetFillColor(kWhite);
		TGraph *gr3 = new TGraph(vecTries.size(),&vecTries[0],&rightLefts[0]);
		gr3->SetTitle("rightLeft asymmetry");
		gr3->SetLineColor(kBlue);
		gr3->SetMarkerStyle(22);
		gr3->SetFillColor(kWhite);
		TGraph *gr4 = new TGraph(vecTries.size(),&vecTries[0],&skewnesses[0]);
		gr4->SetTitle("skewness");
		gr4->SetLineColor(kGreen);
		gr4->SetMarkerStyle(23);
		gr4->SetFillColor(kWhite);
		TString gTitle = TString::Format("convergence of crosstalk #alpha (#alpha_{final} = %02.2f %%) for ",alpha*100);
		gTitle.Append(TPlaneProperties::getStringForDetector(det).c_str());

		TMultiGraph *mg = new TMultiGraph("mg",gTitle);
		mg->Add(gr1,"lp");
		mg->Add(gr2,"lp");
		mg->Add(gr3,"lp");
		mg->Add(gr4,"lp");
		TString cName = "gCrossTalkConvergence_";
		cName.Append(TPlaneProperties::getStringForDetector(det).c_str());
		TCanvas* c2 = new TCanvas(cName);
		c2->cd();
		mg->Draw("APL");
		TCutG * cut = new TCutG("cutBestTry",2);
		cut->SetPoint(0,bestTry,-9999);
		cut->SetPoint(1,bestTry,+9999);
		cut->SetLineColor(kRed);
		cut->SetLineStyle(2);

		TLegend *leg = c2->BuildLegend();
		leg->SetFillColor(kWhite);

		cut->Draw("same");
		leg->Draw();
		c2->SetGridy();
		c2->Draw();
		if (verbosity>4)
		cout<<"save: "<<c2->GetName()<<endl;
		histSaver->SaveCanvas(c2);
		histSaver->SaveCanvas(c2);
		delete c2;
}

bool TAnalysisOfAsymmetricEta::checkBadConvergence(UInt_t nTries){
	if(nTries == maxTriesAlpha){
		cout<<"Asymmetric Eta Sample of detector "<<det<<" did not converge, searching for best correction factor."<<endl;
		// take 'best' alpha where skweness is the smallest
		Float_t minSkewness = 1e9;
		Float_t newBestAlpha = 0;
		bestTry = 0;
		if(skewnesses.size()>0){
			minSkewness = TMath::Abs(skewnesses.at(0));
			UInt_t i=0;
			for( i = 0; i < skewnesses.size(); i++){
				Float_t skewness = TMath::Abs(skewnesses.at(i));
				if(skewness<minSkewness){
					minSkewness = skewness;
					newBestAlpha = alphaValues.at(i)/10;
//					cout<<"Set new alpha in step "<< i<<":"<<newBestAlpha<<endl;
					bestTry = i;
				}
			}
			cout<<"Found best Alpha in step no "<< bestTry <<"/"<<nTries << ": "<< newBestAlpha *10 <<"% with skewness of "<<minSkewness*100<<"%"<<endl;
			alpha = newBestAlpha;

		}
		this->bestAlpha = newBestAlpha;
		if(verbosity%2==1&&verbosity>6){cout<<"Press a key"<<endl;char t; cin>>t;}
		return true;
	}
	else{
		bestAlpha = alpha;
		bestTry = nTries;
	}
	return false;
}

void TAnalysisOfAsymmetricEta::saveAsymmetricEtaPerArea(TH2F* histo,  TH2F* histoNoCor, TString startname, Float_t alpha){
	histSaver->SaveHistogram(hAsymmetricEta2D,false);
	//	TString name = "hAsymmetricEtaFinal_";
	int maxDia =(int)settings->getNDiamonds()+2;
	if(TPlaneProperties::isSiliconDetector(this->det))
		maxDia = 0;
	if(!startname.Contains(TPlaneProperties::getStringForDetector(det)))
		startname.Append(TPlaneProperties::getStringForDetector(det));

	if(startname.Contains("hAsymmetricEtaFinal")){
		TH2F* histo2 = (TH2F*)histo->Clone(startname);
		histSaver->SaveHistogram(histo2);
		delete histo2;
	}

	for ( Int_t dia=-1; histo->GetYaxis()->GetBinCenter(dia)<=maxDia; dia++ ){
		TString histName = startname;
		if(!startname.Contains("hAsymmetricEtaFinal"))
			histName.Append(TString::Format("%05d_",(int)(alpha*100000)));
		if(dia==-1)
			histName.Append("All");
		else
			histName.Append(TString::Format("Area%d",(int)histo->GetYaxis()->GetBinCenter(dia)));

		TH1F* hist = 0;
		TH1F* histNoCor = 0;
		if(dia==-1){
			hist = (TH1F*)histo->ProjectionX(histName);
			if(histoNoCor) histNoCor =  (TH1F*)histoNoCor->ProjectionX();
		}
		else{
			hist = (TH1F*)histo->ProjectionX(histName,dia,dia);
			if (histoNoCor) histNoCor =  (TH1F*)histoNoCor->ProjectionX("_px",dia,dia);
		}
		if(histoNoCor && alpha!=0) {
			histNoCor->SetLineColor(kGray);
//			histNoCor->SetLineStyle(2);
//			histNoCor->Smooth(4);
		}
		//			cout<<"BinCenter of "<<dia<<": "<<histo->GetYaxis()->GetBinCenter(dia)<<endl;
		TString title = startname;
		if(dia==-1)
			title.Append(TString::Format("%02.3f_All",alpha*100));
		else
			title.Append(TString::Format("%02.3f Area%d",alpha*100,(int)histo->GetYaxis()->GetBinCenter(dia)));
		int nbins = 512;
		TString histName2 = histName;
		TString title2 = title;
		histName2.Append("inverted");
		title.Append(" inverted");
//		hist->Smooth(4);
		TH1F* histInverted = new TH1F(histName2,title2,nbins,0,1);
		for(int bin = 1; bin <= nbins;bin++){
			histInverted->SetBinContent(bin,hist->GetBinContent(nbins+1-bin));
		}
		histInverted->SetLineColor(kBlue-7);
//		if(hist)cout<<hist->GetName()<<": "<<hist<<endl;
//		if(histInverted)cout<<histInverted->GetName()<<": "<<histInverted<<endl;
//		if(histNoCor) cout<<histNoCor->GetName()<<": " << histNoCor<<endl;
//		//			cout<<title<<" "<<hist<<endl;
		if(hist){
			hist->SetTitle(title);
			string name = "c_";
			name.append(histName);
			if(hist->GetEntries()>0){
				TCanvas *c1 = new TCanvas(name.c_str());
				c1->cd();
				hist->Draw();
				if(histInverted && alpha!=0) histInverted->Draw("same");
				hist->Draw("same");
				histSaver->SaveCanvas(c1);
				if(histNoCor&& histInverted && alpha!=0) {
					name.insert(1,"1");
					name.append("_diff");
					c1->Clear();
					c1->SetName(name.c_str());
					c1->SetTitle(name.c_str());
					name = "stack_";
					name.append(histName);
					THStack *stack = new THStack(name.c_str(),hist->GetTitle());
					stack->Add(histNoCor,"");//->Draw("AXIS");
					stack->Add(histInverted,"");//->Draw("same");
					stack->Add(hist,"");//->Draw("same");
					stack->Draw("nostack");
//					histNoCor->Draw("same");
					histSaver->SaveCanvas(c1);
				}
			}
			else{
				TCanvas *c1 = new TCanvas(name.c_str());
				histSaver->SaveCanvas(c1);
			}
//				histSaver->SaveTwoHistos(name,hist,histInverted);
			if(hist) delete hist;

		}
		if(histInverted) delete histInverted;
		if(histNoCor) delete histNoCor;
	}
}

void TAnalysisOfAsymmetricEta::Reset() {
	alphaValues.clear();
	means.clear();
	skewnesses.clear();
	rightLefts.clear();
	vecTries.clear();
}

TPolyMarker* TAnalysisOfAsymmetricEta::FindPeaks(TH1F* histo, int nPeaks, Float_t sigma, TString option, Float_t threshold){
	histo->Smooth(2);
	return FindPeaks(maxTriesPeakfinding,histo,nPeaks,sigma,option,threshold);
}

TPolyMarker* TAnalysisOfAsymmetricEta::FindPeaks(UInt_t nTries,TH1F* histo, int nPeaks, Float_t sigma, TString option, Float_t threshold){
	if (!histo)
		return NULL;
	int n = histo->ShowPeaks(sigma,option,threshold);
	TList *functions = histo->GetListOfFunctions();
	TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");

	if(nPeaks != n && nTries != 0){
		if(verbosity>8)cout<<"Peakfinding: " << nTries <<"/"<<maxTriesPeakfinding<<" "<< n <<" / "<<nPeaks<< " Peaks, thrs. " <<threshold<<", sigma: " << sigma <<endl;
		if (n>nPeaks){
			if(verbosity>10)cout<<"increase threshold"<<endl;
			threshold*=1.1;
		}
		else if(n<nPeaks){
			threshold/=.9;
			if(verbosity>10)cout<<"decrease threshold"<<endl;
		}
		return FindPeaks(nTries-1,histo,nPeaks,sigma,option,threshold);
	}
	else
		if (nTries==0)
			if(verbosity>6)cout<<"Couldn't find "<< nPeaks << " peaks, having: "<<n<<" Peaks"<<endl;

	return pm;



}
