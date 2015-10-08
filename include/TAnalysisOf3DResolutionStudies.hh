/*
 * TAnalysisOf3DResolutionStudies.hh
 *
 *  Created on: Oct 6, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DRESOLUTIONSTUDIES_HH_
#define TANALYSISOF3DRESOLUTIONSTUDIES_HH_
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
#include "THStack.h"
#include "TPaletteAxis.h"
#include "TTracking.hh"

//#include "THistogramManager.hh"
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"

class TAnalysisOf3DResolutionStudies {
    public:
        TAnalysisOf3DResolutionStudies(TSettings *settings,HistogrammSaver *histSaver,bool bTransAna=false);
        virtual ~TAnalysisOf3DResolutionStudies();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y);
        void initHistos();
        void saveHistos();
        void setEventReader(TTracking* eventReader);
        void setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster){
            this->isTransparentCluster=isTransparentCluster;this->transparentCluster=transparentCluster;}
        void setDiamondCluster(TCluster* diamondCluster){this->diamondCluster=diamondCluster;}
    private:
        void fillResolutionPlots();
        void CreateResolutionPlots(vector<TH1F*>*vec,TString kind);

        vector<TH1F*> vecHResolutionPerCell_maxValue;
        vector<TH1F*> vecHResolutionPerCell_chargeWeighted;
        vector<TH1F*> vecHResolutionPerCell_highest2Centroid;
    private:
        TString appendix;
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
        bool bTransAna;
};

#endif /* TANALYSISOF3DRESOLUTIONSTUDIES_HH_ */
