/*
 * TAnalysisOf3DStripAnalysis.hh
 *
 *  Created on: Oct 7, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DSTRIPANALYSIS_HH_
#define TANALYSISOF3DSTRIPANALYSIS_HH_
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

class TAnalysisOf3DStrip {
    public:
        TAnalysisOf3DStrip(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna=false);
        virtual ~TAnalysisOf3DStrip();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y);
        void initHistos();
        void saveHistos();
        void setEventReader(TTracking* eventReader){this->eventReader=eventReader;}
        void setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster){
            this->isTransparentCluster=isTransparentCluster;this->transparentCluster=transparentCluster;}
        void setClusteredCluster(TCluster* cluster){this->clusteredCluster = cluster;}
        void setDiamondCluster(TCluster* diamondCluster){this->diamondCluster=diamondCluster;}
        void setValidClusterAnalysis(bool val){validClusteredAnalysis=val;}
        void setValidTransparentAnalysis(bool val){validTransparentAnalysis=val;}
        void setEventNo(UInt_t nEvent){this->nEvent=nEvent;}
        TH1F* getLandau(){return hLandauStrip;}
    private:
        void fillPlots();
        TCluster* clusteredCluster;
        bool validClusteredAnalysis, validTransparentAnalysis;
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
        TH1F* hLandauStrip;
        TH2F* hLandauStripFidCutXvsFidCutY;

        TH2F* hPHvsChannelStrip;
        TH2F* hPHvsPredictedXPosStrip;
};

#endif /* TANALYSISOF3DSTRIPANALYSIS_HH_ */
