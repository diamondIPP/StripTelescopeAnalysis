/*
 * TAnalysisOf3DGoodCellsLandau.hh
 *
 *  Created on: Oct 6, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DGOODCELLSLANDAU_HH_
#define TANALYSISOF3DGOODCELLSLANDAU_HH_
#include <fstream>
#include <iostream>
#include <iostream>
#include <algorithm>    // std::sort
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TAnalysisOfAnalysisDifferences.hh"

#include "TTracking.hh"
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"

class TAnalysisOf3DGoodCellsLandau {
    public:
        TAnalysisOf3DGoodCellsLandau(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna=false);
        virtual ~TAnalysisOf3DGoodCellsLandau();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y);
        void initHistos();
        void saveHistos(TH1F* hLandauStrip, TH1F* hLandauPhantom);
        void setEventReader(TTracking* eventReader){this->eventReader=eventReader;}
        void setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster){
            this->isTransparentCluster=isTransparentCluster;this->transparentCluster=transparentCluster;}
        void setClusteredCluster(TCluster* cluster){this->clusteredCluster = cluster;}
        void setDiamondCluster(TCluster* diamondCluster){this->diamondCluster=diamondCluster;}
        void setValidClusterAnalysis(bool val){validClusteredAnalysis=val;}
        void setValidTransparentAnalysis(bool val){validTransparentAnalysis=val;}
        void setEventNo(UInt_t nEvent){this->nEvent=nEvent;}
    private:
        TCluster* clusteredCluster;
        bool validClusteredAnalysis, validTransparentAnalysis;
        void FillGoodCellsLandaus(Float_t charge);
        Int_t PulseHeightBins, PulseHeightMin, PulseHeightMax,PulseHeightMaxMeanCharge,PulseHeightMinMeanCharge;
        UInt_t nEvent;
        bool useCMN;
        bool isTransparentCluster;
        Int_t verbosity;
        Float_t chi2x, chi2y, predx, predy, fidx,fidy;
        TCluster *transparentCluster;
        TCluster* diamondCluster;
        TTracking* eventReader;
        HistogrammSaver* histSaver;
        TSettings* settings;
        UInt_t subjectDetector;
        UInt_t subjectPlane;
        TString appendix;
        bool bTransAna;
        map<Int_t, TCluster> mapClusteredAnalysisGoodCells;
        map<Int_t, TCluster> mapTransparentAnalysisGoodCells;
        map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsGoodCells;
//        map<Int_t, TCluster> mapClusteredAnalysisAllCells;
//        map<Int_t, TCluster> mapTransparentAnalysisAllCells;
//        map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsAllCells;
        TH1F* hLandauGoodCellsWithoutEdges;
        TH1F* hLandauGoodCellsWithoutColumns;
        TH1F* hLandauGoodCells;

        TProfile2D* hPulseHeightVsDetectorHitPostionXYGoodCells;
};

#endif /* TANALYSISOF3DGOODCELLSLANDAU_HH_ */
