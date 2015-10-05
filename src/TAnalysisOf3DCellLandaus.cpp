/*
 * TAnalysisOf3DCellLandaus.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DCellLandaus.hh>

TAnalysisOf3D_CellLandaus::TAnalysisOf3D_CellLandaus(TSettings *settings,HistogrammSaver *histSaver) {
    // TODO Auto-generated constructor stub
    this->settings = settings;
    this->histSaver = histSaver;
    PulseHeightBins = 1024;
    PulseHeightMin = 0;
    PulseHeightMax = 2500;
    verbosity = 0;
    hLandauStrip = 0;
}

TAnalysisOf3D_CellLandaus::~TAnalysisOf3D_CellLandaus() {
    // TODO Auto-generated destructor stub
}

void TAnalysisOf3D_CellLandaus::initQuarterCellLandaus() {
    TString appendix  = "";
    hQuarterCellsLandau.resize(settings->GetNCells3d());
    hQuarterCellsClusterSize.resize(settings->GetNCells3d());
    TString name;
    for (UInt_t cell =0; cell < settings->GetNCells3d(); cell ++){
        //        hQuarterCellsLandau.push_back(vector<TH1F*>());
        hQuarterCellsLandau[cell].resize(settings->getNQuarters3d());
        hQuarterCellsClusterSize[cell].resize(settings->getNQuarters3d());
        for(UInt_t quarter=0;quarter<settings->getNQuarters3d();quarter++){
            name = TString::Format("hQuaterCellsLandau_%d_%d",cell,quarter);
            name.Append(appendix);
            if(verbosity>1)cout<<"Create "<<name<<endl;
            hQuarterCellsLandau[cell][quarter] =
                    new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
            name = TString::Format("hQuarterCellsClusterSize_%d_%d",cell,quarter);
            name.Append(appendix);
            hQuarterCellsClusterSize[cell].push_back(new TH1F(name,name,6,0,6));
        }
    }
}

void TAnalysisOf3D_CellLandaus::addEvent(int cellNo, int quarterNo, float charge) {
    if(0<=cellNo && cellNo < hQuarterCellsLandau.size()){
        if (0<=quarterNo && quarterNo < hQuarterCellsLandau[cellNo].size())
            hQuarterCellsLandau[cellNo][quarterNo]->Fill(charge);
    }
}

TH1F* TAnalysisOf3D_CellLandaus::getCellLandau(unsigned int cell) {
    TString appendix = "";
    TString name = TString::Format("hCellsLandau_%d",cell);
    name.Append(appendix);
    if(verbosity>1)cout<<"Create "<<name<<endl;
    TH1F* hCellLandau = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hCellLandau->Reset();
    for (UInt_t quarter =0; quarter< hQuarterCellsLandau[cell].size();quarter++)
        hCellLandau->Add(hQuarterCellsLandau[cell][quarter]);
    return hCellLandau;
}

void TAnalysisOf3D_CellLandaus::SaveGoodAndBadCellLandaus() {
    cout<<"[TAnalysisOf3D_CellLandaus::SaveGoodAndBadCellLandaus]"<<endl;
    TString appendix ="";
    if (settings->do3dTransparentAnalysis())
        appendix ="_trans";
    TList *listGoodCells = new TList;
    TList *listBadCells = new TList;
    string plots_path = histSaver->GetPlotsPath();
    histSaver->SetPlotsPath(plots_path+(string)"/CellLandaus/");

    for(UInt_t column=0;column<settings->getNColumns3d();column++){
        for(UInt_t row=0;row<settings->getNRows3d();row++){
            Int_t cell = settings->get3DCellNo((int)column,row);
            //hCellNumbering->SetBinContent(column+1,row+1,cell); //This should be a clone of the 2D Cell Mean Charge Plot, Wait till Felix has finished.
            TH1F* hCellLandau = getCellLandau(cell);
            histSaver->SaveHistogram(hCellLandau);
            for(UInt_t i=0; i<settings->getBadCells3D().size(); i++)
                if(cell==settings->getBadCells3D().at(i)){
                    //                    Int_t Entries = hLandauBadCells->GetEntries();
                    //                    hLandauBadCells->Add(hCellsLandau.at(cell),1);    //Not working for some reason, ask Felix
                    //                    hLandauBadCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                    listBadCells->Add(hCellLandau);
                    //                    cout<<"\tbad"<<endl;
                }
            for(UInt_t i=0; i<settings->getGoodCellRegions3d().size(); i++){
                for(UInt_t j=0; j<settings->getGoodCellRegions3d().at(i).size(); j++)
                    if(cell==settings->getGoodCellRegions3d().at(i).at(j)){
                        //                        Int_t Entries = hLandauGoodCells->GetEntries();
                        //                        hLandauGoodCells->Add(hCellsLandau.at(cell)); //Not working for some reason, ask Felix
                        //                        hLandauGoodCells->SetEntries(Entries+hCellsLandau.at(cell)->GetEntries());
                        listGoodCells->Add(hCellLandau);
                        //                        cout<<"\tgood"<<endl;
                    }
            }
        }
    }
    histSaver->SetPlotsPath(plots_path);
    cout<<"List Good Cells: "<<listGoodCells->GetEntries()<<endl;
    listGoodCells->Print();
    cout<<"\nList Bad Cells: "<<listBadCells->GetEntries()<<endl;
    listBadCells->Print();

    TString name = "hLandauBadCells";
    name.Append(appendix);
    TH1F* hLandauBadCells = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauBadCells->SetTitle("Landau of bad cells");
    hLandauBadCells->GetXaxis()->SetTitle("pulse height /adc");
    hLandauBadCells->GetXaxis()->SetTitle("number of entries #");
    hLandauBadCells->Reset();
    hLandauBadCells->Merge(listBadCells);
    histSaver->SaveHistogram(hLandauBadCells);

    name = "hLandauGoodCells";
    name.Append(appendix);
    TH1F* hLandauGoodCells = new TH1F(name,name,PulseHeightBins,PulseHeightMin,PulseHeightMax);
    hLandauGoodCells->SetTitle("Landau of good cells");
    hLandauGoodCells->GetXaxis()->SetTitle("pulse height /adc");
    hLandauGoodCells->GetXaxis()->SetTitle("number of entries #");
    hLandauGoodCells->Reset();
    hLandauGoodCells->Merge(listGoodCells);
    //    histSaver->SaveHistogram(hLandauGoodCells);


    cout<<"hLandauGoodCells: "<<hLandauGoodCells->GetEntries()<<endl;
    cout<<"hLandauStrip:     "<<hLandauStrip->GetEntries()<<endl;

    Float_t factor = hLandauGoodCells->GetBinContent(hLandauGoodCells->GetMaximumBin());
    factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
    name = "cLandauGoodCells";
    name.Append(appendix);
    histSaver->SaveTwoHistos(name,hLandauGoodCells,hLandauStrip,factor,"right");

    name = "cLandauGoodCellsNormalized";
    name.Append(appendix);
    histSaver->SaveTwoHistosNormalized(name,hLandauGoodCells,hLandauStrip,1,"right");

    factor = hLandauBadCells->GetBinContent(hLandauBadCells->GetMaximumBin());
    factor/= (Float_t) hLandauStrip->GetBinContent(hLandauStrip->GetMaximumBin());
    name = "cLandauBadCells";
    name.Append(appendix);
    histSaver->SaveTwoHistos(name,hLandauBadCells,hLandauStrip,factor,"right");
    name = "cLandauBadCellsNormalized";
    name.Append(appendix);
    histSaver->SaveTwoHistosNormalized(name,hLandauBadCells,hLandauStrip,1,"right");
}
