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
	histSaver=new HistogrammSaver();

	maxTriesPeakfinding = 20;
	maxTriesAlpha = 60;
	settings->goToClusterAnalysisDir();
	histSaver->SetPlotsPath(settings->getClusterAnalysisDir().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	verbosity = settings->getVerbosity();

	hAsymmetricEta2D=0;
}

TAnalysisOfAsymmetricEta::~TAnalysisOfAsymmetricEta() {
	// TODO Auto-generated destructor stub
}


void TAnalysisOfAsymmetricEta::FillEtaDistribution(TH2F* histo){
	//Fill histo with Asymmetric Cluster data..
	cout<<" Fill " << histo->GetName() <<" with " << clusters.size() << "Entries."<< endl;
	for (UInt_t i=0;i<clusters.size();i++){
		TCluster  clus = clusters.at(i);
		if(clus.getClusterSize()>=settings->getMaxAllowedClustersize(det))
			continue;
		TCluster newClus = clus.getCrossTalkCorrectedCluster(alpha);
		Int_t leftChannel;
		Float_t newEta = newClus.getEta(leftChannel);
		Int_t nDia = settings->getDiaDetectorAreaOfChannel(leftChannel)+1;
		if(newEta>0&&newEta<1){
			histo->Fill(newEta,nDia);
			//				cout<<i<<" nDia: "<<nDia<<" eta: "<<newEta<<endl;
		}
	}
	cout<<"histo: "<<histo->GetName()<<" "<<histo->GetEntries()<<endl;
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
		if(settings->getVerbosity()>4)cout<<"Diamond Detector"<<endl;
		int nDias = settings->getNDiamonds();
		int nDia = settings->getAnalysedDiamond();
		if(settings->getVerbosity()>4)cout<<"nDias: "<<nDias<<endl;
		if(settings->getVerbosity()>4)cout<<"nDia: "<<nDia<<endl;
		int bin;
		if(settings->getNDiamonds()!=1){
			if (true){//settings->getAnalysedDiamond()==TSettings::leftDia||true){//todo
				bin = hAsymmetricEta2D->GetYaxis()->FindBin(0.0);
				if(settings->getVerbosity()>4)cout<<"bin:"<<bin<<endl;
				TString name = hName;//+"_left";//TString::Format("%s_left",hName);
				name.Append(TString::Format("%06d_left",(int)(TMath::Abs(alpha)*100000)));
				if(settings->getVerbosity()>4)cout<<"Title: "<<flush;
				hTitle.Append(TString::Format("%02.2f %%, left",alpha*100));
				if(settings->getVerbosity()>4)cout<<hTitle<<endl;
				if(settings->getVerbosity()>4)cout<<name<<endl;
				hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(name,bin+1,bin+1);
				if(settings->getVerbosity()>4)cout<<"lDia:"<<hProjection<<endl;
			}
			else if(settings->getAnalysedDiamond()==TSettings::rightDia){//todo
			}
			else {//todo

			}
		}
		else{// has only one diamond
			hName.Append("all_px");
			hTitle.Append(TString::Format("%02.2f %%, all",alpha*100));
			hProjection = (TH1F*) hAsymmetricEta2D->ProjectionX(hName);
			if(settings->getVerbosity()>4)cout<<"1Dia"<<hProjection<<endl;
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
	//	hAsymmetricEta2D = new TH2F(hName,hName,512,0,1,nDiamonds+3,-1.5,nDiamonds+1.5);
	//	TH1F* hAsymmetricEta = new TH1F(hName,hName,512,0,1);
	alpha= 0;
	Float_t mean = 0;
	UInt_t nTries =0;
	TPolyMarker *pm = 0;

	TH1F* hAsymmetricEta_px=0;
	bool valid=false;
	Reset();
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
		if(verbosity>4)cout<<" "<<hAsymmetricEta2D->GetEntries()<<endl;
		hAsymmetricEta2D->Reset();
		hAsymmetricEta2D->Clear();
		TString name = hName;
		name.Append(TString::Format("_%05d",(int)(TMath::Abs(alpha)*10000)));
		hAsymmetricEta2D->SetName(name);
		name = hName;
		name.Append(TString::Format("_%03.2f",alpha));
		hAsymmetricEta2D->SetTitle(name);
		FillEtaDistribution(hAsymmetricEta2D);
		//		cout<<"save asym. etaDists"<<endl;
		if(alpha==0) saveAsymmetricEtaPerArea(hAsymmetricEta2D,"hAsymmetricEta_",alpha);
		//		cout<<"delete projection"<<endl;
		if(hAsymmetricEta_px) {
			delete hAsymmetricEta_px;
			hAsymmetricEta_px=0;
		}
		//		cout<<"get Projection"<<endl;
		if(verbosity>4)cout<<hAsymmetricEta2D<<endl;
		hAsymmetricEta_px = getProjection();
		//		cout<<"find new alpha"<<endl;
		int leftHalf = 0;
		int bin=0;
		for (bin=0; hAsymmetricEta_px->GetBinLowEdge(bin)<.5;bin++){
			leftHalf+=hAsymmetricEta_px->GetBinContent(bin);
		}
		int rightHalf = 0;
		for (bin=bin; bin<hAsymmetricEta_px->GetNbinsX();bin++){
			rightHalf+=hAsymmetricEta_px->GetBinContent(bin);
		}
		mean = hAsymmetricEta_px->GetMean();
		if(verbosity>4)cout<<"Analysising: "<<hAsymmetricEta_px->GetName()<<endl;
		if(verbosity>4)cout<<"Asymmetric eta with a charge share of "<<alpha*100<<"%"<<endl;
		if(verbosity>4)cout<<"nentries: "<<hAsymmetricEta_px->GetEntries()<<endl;
		Float_t skewness = hAsymmetricEta_px->GetSkewness();
		if(verbosity>4)cout<<"skewness*100: "<<skewness*100<<endl;
		if(verbosity>4)cout<<"left: "<<leftHalf<<"\tright: "<<rightHalf<<" "<<(Float_t)leftHalf/(Float_t)rightHalf*100<<endl;

		valid = TMath::Abs(mean-.5)<.005;
		valid = valid || TMath::Abs(skewness)<0.005;
		if(verbosity>4)cout<<"Mean "<<mean*100<<" "<<valid<<endl;
		pm = FindPeaks(hAsymmetricEta_px,2);
		Float_t p1,p2;
		if(pm){
			for (int i=0;i< pm->GetN();i++){
				Float_t peakPos = pm->GetX()[i];
				if (peakPos>.5)
					peakPos = TMath::Abs(peakPos-1);
				if(i==0)
					p1 = peakPos;
				else if(i==1)
					p2 = peakPos;
				if(verbosity>4)cout<<"\t"<<i<<"\t"<<peakPos*100.<<": "<<pm->GetY()[i]<<"\n";
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
		rightLefts.push_back((Float_t)leftHalf/(Float_t)rightHalf-1);
		skewnesses.push_back(skewness);
		means.push_back(mean-.5);
		vecTries.push_back(nTries);
		if(!valid){
			if (alpha ==0){
				if(det == 6 || det ==2){
					if(mean- .5>0)
						alpha = -.01;
					else
						alpha = +.01;
				}
				else{
					if(mean- .5>0)
						alpha = .01;
					else
						alpha = -.01;
				}
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
				if (skewness > 0 ){
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
		}
		nTries++;
	}

	if (true)
	{//convergence plot
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
		TString gTitle = TString::Format("convergence of crosstalk #alpha (#alpha_{final} = %02.2f %% for ",alpha*100);
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
		TLegend *leg = c2->BuildLegend();
		leg->SetFillColor(kWhite);
		leg->Draw();
		c2->SetGridy();
		c2->Draw();
		cout<<"save: "<<c2->GetName()<<endl;
		histSaver->SaveCanvas(c2);
		histSaver->SaveCanvas(c2);
		delete c2;
	}
	cout<<"Convergered after "<<nTries<<"/"<<maxTriesAlpha<<" best alpha: "<<alpha<<endl;
	hName ="hAsymmetricEta_";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	hName.Append(TString::Format("%05d",(int)(TMath::Abs(alpha)*10000)));
	hAsymmetricEta_px->SetName(hName);
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
	cout<<"nDiamonds: "<<settings->getNDiamonds()<<endl;
	hName = "hAsymmetricEtaFinal_";
	hName.Append(TPlaneProperties::getStringForDetector(det).c_str());
	saveAsymmetricEtaPerArea(hAsymmetricEta2D,hName,alpha);
	TString name = "hDiamondHitNo";
	histSaver->SaveHistogram((TH1F*)hAsymmetricEta2D->ProjectionY(name));
	//	histSaver->SaveHistogram(hAsymmetricEta);
	cout<<flush;
	//	char t;
	//	cin>>t;
	return nTries;
}




void TAnalysisOfAsymmetricEta::saveAsymmetricEtaPerArea(TH2F* histo,  TString startname, Float_t alpha){
	histSaver->SaveHistogram(hAsymmetricEta2D);
	//	TString name = "hAsymmetricEtaFinal_";
	int maxDia =(int)settings->getNDiamonds()+2;
	if(TPlaneProperties::isSiliconDetector(this->det))
		maxDia = 0;
	if(!startname.Contains(TPlaneProperties::getStringForDetector(det)))
		startname.Append(TPlaneProperties::getStringForDetector(det));
	for ( Int_t dia=-1; histo->GetYaxis()->GetBinCenter(dia)<maxDia; dia++ ){
		TString histName = startname;
		if(dia==-1)
			histName.Append(TString::Format("%05d_All",(int)(alpha*100000)));
		else
			histName.Append(TString::Format("%05d_Area%d",(int)(alpha*100000),(int)histo->GetYaxis()->GetBinCenter(dia)));

		TH1F* hist;
		if(dia==-1)
			hist = (TH1F*)histo->ProjectionX(histName);
		else
			hist = (TH1F*)histo->ProjectionX(histName,dia,dia);
		//			cout<<"BinCenter of "<<dia<<": "<<histo->GetYaxis()->GetBinCenter(dia)<<endl;
		TString title = startname;
		if(dia==-1)
			title.Append(TString::Format("%02.3f_All",alpha*100));
		else
			title.Append(TString::Format("%02.3f_Area%d",alpha*100,(int)histo->GetYaxis()->GetBinCenter(dia)));
		//			cout<<title<<" "<<hist<<endl;
		if(hist){
			hist->SetTitle(title);
			histSaver->SaveHistogram(hist,false,false,true);
			delete hist;
		}
		else
			cout<<"hist==0 "<<histName<<endl;
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
