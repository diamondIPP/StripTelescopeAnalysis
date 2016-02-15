/*
 * TAnalysisOfGoodCells.hh
 *
 *  Created on: Dec 13, 2013
 *      Author: bachmair
 */

#ifndef TANALYSISOFGOODCELLS_HH_
#define TANALYSISOFGOODCELLS_HH_
#include "TSettings.class.hh"
#include "TCluster.hh"
#include "HistogrammSaver.class.hh"
#include "TH1F.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include <map>
using namespace std;
class TAnalysisOfAnalysisDifferences {
public:
    TAnalysisOfAnalysisDifferences(TSettings* settings,HistogrammSaver* histSaver,TString extension="");
    void Analysis();
//    void SetTransparentMap(std::map<Int_t,TCluster>* transparentMap);
//    void SetClusteredMap(std::map<Int_t,TCluster>* clusteredMap);
    virtual ~TAnalysisOfAnalysisDifferences();
    std::map<Int_t, TCluster>* getClusteredMap() const;
    void setClusteredMap (std::map<Int_t, TCluster>* clusteredMap);
    std::map<Int_t, TCluster>* getTransparentMap() const;
    void setTransparentMap( std::map<Int_t, TCluster>* transparentMap);
    std::map<Int_t,std::pair<Float_t,Float_t> >* getPredictedPositions() const;
    void setPredictedPositions( std::map<Int_t, pair<Float_t, Float_t> >* predictedPositionMap);
    void setStripHistogram(TH1F* histo);
    void setNegativeChargeCut(Float_t value){negChargeCut = value;}
    Float_t getNegativeChargeCut(){return negChargeCut;}
    void set3DPhantomLandau(TH1F* hPhantomLandau);
private:
    TH1F* hPhantomLandau;
    void InitHistograms();
    void InitTransparentHistos();
    void InitClusteredHistos();
    void InitSameHistos();
    void SaveHistograms();
    void SaveTransparentClusteredComparison();
    void SaveComparisonPlots(TString name,TH1* histo,bool includePhantom=false);
    void LoopOverBothMaps();
    void UpdatePredictedPosition();
    void AnalyseSameEvent();
    void AnalyseTransparentEvent();
    void AnalyseOnlyTransparentEvent();
    void AnalyseOnlyClusteredEvent();
//    bool hasNegativeCharge(std::map<Int_t,TCluster>::iterator it,Int_t &pos,Float_t& charge);

private:
    TH1F* stripHisto;
    TSettings *settings;
    string oldPlotPath;
    Int_t verbosity;
    Float_t negChargeCut;
    HistogrammSaver *histSaver;
    std::map<Int_t,TCluster>* transparentMap;
    std::map<Int_t,TCluster>* clusteredMap;
    map<Int_t,pair<Float_t, Float_t> > *predictedPositions;

private:
    TString extension;
    Int_t nSameEvents;
    Int_t nOnlyClustered;
    Int_t nOnlyTransparent;
    map<Int_t,TCluster>::iterator itClustered;
    map<Int_t,TCluster>::iterator itTransparent;
    map<Int_t, std::pair<Float_t,Float_t> >::iterator itPredicted;
private:
    TH1F* hChargeDifference;
    map<TString,TH1*> mapHistos;
};

#endif /* TANALYSISOFGOODCELLS_HH_ */
