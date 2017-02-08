/*
 * TAvrgChargePerBinMonteCarlo.cpp
 *
 *  Created on: Sep 17, 2015
 *      Author: bachmair
 */

#include <TAvrgChargePerBinMonteCarlo.hh>

TAvrgChargePerBinMonteCarlo::TAvrgChargePerBinMonteCarlo(
        HistogrammSaver* histSaver) {
    this->histSaver=histSaver;
}

TAvrgChargePerBinMonteCarlo::~TAvrgChargePerBinMonteCarlo() {
    // TODO Auto-generated destructor stub
}

void TAvrgChargePerBinMonteCarlo::doMonteCarlo(TProfile2D* profOverlay,
        TH1F* hLandauOfOverlay) {
    if(!profOverlay || !hLandauOfOverlay){
        cerr<<"[TAnalysisOf3dDiamonds::DoMonteCarloOfAvrgChargePerBinInOverlay] one of the histograms is not defined! ";
        cout<<profOverlay<<" "<<hLandauOfOverlay<<endl;
        return;
    }
    this->hLandauOfOverlay = hLandauOfOverlay;
    init_histos(profOverlay);
    do_loop();
    save_histos();
}

void TAvrgChargePerBinMonteCarlo::init_histos(TProfile2D* profOverlay) {
    hBinContents = histSaver->GetBinContentHisto(profOverlay);
    TString name = profOverlay->GetName()+(TString)"_Entries";
    Int_t max = hBinContents->GetBinContent(hBinContents->GetMaximumBin())+1;
    Int_t min = hBinContents->GetBinContent(hBinContents->GetMinimumBin())+1;
    max -=5;
    min -=5;
    Int_t bins = max-min;
    hEntries = new TH1F(name,name,bins,min,max);
    for(Int_t binx = 1;binx < hBinContents->GetNbinsX();binx++)
        for(Int_t biny = 1;biny < hBinContents->GetNbinsY();biny++){
            Int_t bin = hBinContents->GetBin(binx,biny);
            hEntries->Fill(hBinContents->GetBinContent(bin));
        }

    name = "hMonteCarloAvrgChargePerBin";
    hMonteCarloAvrgChargePerBin = new TH1D(name,name,PulseHeightBins,PulseHeightMinMeanCharge-100,PulseHeightMaxMeanCharge+100);
    hMonteCarloAvrgChargePerBin->GetXaxis()->SetTitle("mean charge per bin_{MonteCarlo}");
    hMonteCarloAvrgChargePerBin->GetYaxis()->SetTitle("number of entries #");

    name =TString::Format("hMonteCarloNBelowCut_%d",(int)cut);;
    hNumberOfEntriesBelowCut = new TH1D(name,name,100,0,100);
    hNumberOfEntriesBelowCut->GetXaxis()->SetTitle(TString::Format("number of  entries below %.f",cut));
    hNumberOfEntriesBelowCut->GetYaxis()->SetTitle("number of  entries #");

    name =TString::Format("hMonteCarloRelativeBelowCut_%d",(int)cut);;
    hRelativeNumberBelowCut = new TH1D(name,name,100,0,1);
    hRelativeNumberBelowCut->GetXaxis()->SetTitle(TString::Format("rel. number of  entries below %.f",cut));
    hRelativeNumberBelowCut->GetYaxis()->SetTitle("number of  entries #");
}


void TAvrgChargePerBinMonteCarlo::save_histos() {
    histSaver->SaveHistogram(hMonteCarloAvrgChargePerBin,false,true,true);
    histSaver->SaveHistogram(hNumberOfEntriesBelowCut,false,true,true);
    histSaver->SaveHistogram(hRelativeNumberBelowCut,false,true,true);
    histSaver->SaveHistogram(hEntries);
    histSaver->SaveHistogram(hBinContents,false,false,(TString)"colzText");
    delete hMonteCarloAvrgChargePerBin;
    delete hBinContents;
    delete hEntries;
    delete hNumberOfEntriesBelowCut;
    delete hRelativeNumberBelowCut;
}

void TAvrgChargePerBinMonteCarlo::do_loop() {
    Float_t mean =0;
    Int_t entries = 0;
    int nBelowCut = 0;
    for(UInt_t nMC = 0; nMC < 1e6; nMC++){
        entries = hEntries->GetRandom();
        mean = 0;
        nBelowCut = 0;
        for(Int_t entry =0; entry < entries; entry++){
            Float_t charge =  hLandauOfOverlay->GetRandom();
            mean += charge;
            if(charge<cut)
                nBelowCut ++;
        }
        mean /= (Float_t)entries;
        hNumberOfEntriesBelowCut->Fill(nBelowCut);
        hRelativeNumberBelowCut->Fill((Float_t)nBelowCut/(Float_t)entries);
        //        cout<<TString::Format("%6d --> %6.1f (%d)",nMC,mean,entries)<<endl;;
        hMonteCarloAvrgChargePerBin->Fill(mean);
    }
}
