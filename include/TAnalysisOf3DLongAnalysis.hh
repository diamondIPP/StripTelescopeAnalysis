/*
 * TAnalysisOf3DLongAnalysis.hh
 *
 *  Created on: Sep 18, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DLONGANALYSIS_HH_
#define TANALYSISOF3DLONGANALYSIS_HH_
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

#include "TTracking.hh"
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"
#include "TAnalysisOf3DResolutionStudies.hh"
#include "TAnalysisOf3DGoodCellsLandau.hh"

class TAnalysisOf3DLongAnalysis {
    public:
        TAnalysisOf3DLongAnalysis(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna);
        virtual ~TAnalysisOf3DLongAnalysis();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid,Float_t chi2x, Float_t chi2y);
        void initHistos();
        void saveHistos(TH1F* hLandauStrip,TH1F* hLandauPhantom);
        void setEventReader(TTracking* eventReader){this->eventReader = eventReader;}
        void setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster);
        void setDiamondCluster(TCluster* diamondCluster){this->diamondCluster=diamondCluster;}
        void setEventNo(UInt_t nEvent){this->nEvent=nEvent;}
    private:
        void checkClusteredAnalysis();
        void checkTransparentAnalysis();
        void FillResolutionPlots();
    private:
        TAnalysisOf3DResolutionStudies* resolutionStudy;
        TAnalysisOf3DGoodCellsLandau* goodCellsLandau;
        vector<TH1F*> vecHResolutionPerCell_maxValue;
        vector<TH1F*> vecHResolutionPerCell_chargeWeighted;
        vector<TH1F*> vecHResolutionPerCell_highest2Centroid;
        map<Int_t, TCluster> mapClusteredAnalysisGoodCells;
        map<Int_t, TCluster> mapTransparentAnalysisGoodCells;
        map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsGoodCells;
        map<Int_t, TCluster> mapClusteredAnalysisAllCells;
        map<Int_t, TCluster> mapTransparentAnalysisAllCells;
        map<Int_t, pair<Float_t,Float_t> > mapPredictedPositionsAllCells;
        TH2D* hNegativeChargePosition;
        TH1F* hNegativeCharges;
        TH2F* hInvalidCellNo;
        TH2F* hInvalidCluster;
        TH2F* hValidEventsDetSpace;
        TProfile2D* hPulseHeightVsDetectorHitPostionXY;
        TH1F* hLandauGoodCellsWithoutEdges;
        TH1F* hLandauGoodCellsWithoutColumns;
        TH1F* hLandauGoodCells;

    private:
        UInt_t nEvent;
        Int_t PulseHeightBins, PulseHeightMin, PulseHeightMax,PulseHeightMaxMeanCharge,PulseHeightMinMeanCharge;
        UInt_t subjectDetector;
        UInt_t subjectPlane;
        bool useCMN;
        bool isTransparentCluster;
        Int_t verbosity;
        Float_t chi2x, chi2y, predx, predy, fidx,fidy;
        HistogrammSaver* histSaver;
        TSettings* settings;
        bool bTransAna;
        TString appendix;
        TTracking* eventReader;
        TCluster *transparentCluster;
        TCluster clusteredCluster;
        TCluster *diamondCluster;
        bool validClusteredAnalysis;
        bool validTransparentAnalysis;

};

#endif /* TANALYSISOF3DLONGANALYSIS_HH_ */
