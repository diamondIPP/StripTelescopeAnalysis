/*
 * TAnalysisOf3DResolution.hh
 *
 *  Created on: Mar 7, 2016
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DRESOLUTION_HH_
#define TANALYSISOF3DRESOLUTION_HH_

#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <deque>
#include <algorithm>    // std::min_element, std::max_element

#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile2D.h"

#include "HistogrammSaver.class.hh"
#include "TADCEventReader.hh"
#include "TTracking.hh"
#include "TSettings.class.hh"
#include "TCluster.hh"

class TAnalysisOf3DResolution {
    public:
        TAnalysisOf3DResolution(TSettings *settings,  HistogrammSaver *histSaver, TString appendix);
        virtual ~TAnalysisOf3DResolution();
        void Fill(TCluster* clus,Float_t xPredDet,Float_t yPredDet);
    private:
        typedef std::map< TString, std::vector<TH1*> > TCellHistoMap;
        void initialiseHistos();
        void saveHistos();
        void deleteHistos();
        void SaveResolutionPlots(vector<TH1F*>*vec,TString kind);
        TString appendix;
        TSettings *settings;
        HistogrammSaver *histSaver;
        Int_t subjectDetector;
        Float_t maxsnr;
        TH2F* hAdjacentChannels_SNR;
        TH2F* hAdjacentSNR_vs_cellNo;
        TH2F* hAdjacentChannels_Signal;
        vector<TH1F*> vecHResolutionPerCell_maxValue;
        vector<TH1F*> vecHResolutionPerCell_chargeWeighted;
        vector<TH1F*> vecHResolutionPerCell_highest2Centroid;
        vector<TH1F*> vecHResolutionPerCell_h2C_WithCut;
        vector<TH2*> vecHResolutionPerCell_maxValue_vs_SNR;
        vector<TH2*> vecHResolutionPerCell_chargeWeighted_vs_SNR;
        vector<TH2*> vecHResolutionPerCell_highest2Centroid_vs_SNR;
        vector<TH2*> vecHResolutionPerCell_h2C_WithCut_vs_SNR;
        vector<TH2*> vecHResolutionPerCell_maxValue_vs_PredHit;
        vector<TH2*> vecHResolutionPerCell_maxValue_vs_PredHitY;
        vector<TH2*> vecHResolutionPerCell_chargeWeighted_vs_PredHit;
        vector<TH2*> vecHResolutionPerCell_chargeWeighted_vs_PredHitY;
        vector<TH2*> vecHResolutionPerCell_highest2Centroid_vs_PredHit;
        vector<TH2*> vecHResolutionPerCell_highest2Centroid_vs_PredHitY;
        vector<TH2*> vecHResolutionPerCell_h2C_WithCut_vs_PredHit;
        vector<TH2*> vecHResolutionPerCell_h2C_WithCut_vs_PredHitY;
        TCellHistoMap cellHistos;
        bool useCMN;
};

#endif /* TANALYSISOF3DRESOLUTION_HH_ */
