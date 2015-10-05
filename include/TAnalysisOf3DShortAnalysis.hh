/*
 * TAnalysisOf3DShortAnalysis.hh
 *
 *  Created on: Sep 18, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DSHORTANALYSIS_HH_
#define TANALYSISOF3DSHORTANALYSIS_HH_
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

class TAnalysisOf3DShortAnalysis {
    public:
        TAnalysisOf3DShortAnalysis(TSettings *settings,HistogrammSaver *histSaver);
        virtual ~TAnalysisOf3DShortAnalysis();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid, Float_t chi2x, Float_t chi2y);
        void initHistos();
        void saveHistos(TH1F* hLandauStrip);
        void setEventReader(TTracking* eventReader);
        void setTransparentCluster(bool isTransparentCluster, TCluster* transparentCluster){
            this->isTransparentCluster=isTransparentCluster;this->transparentCluster=transparentCluster;}
        void setDiamondCluster(TCluster* diamondCluster){this->diamondCluster=diamondCluster;}
    private:
//        THistogramManager histos;
        Int_t PulseHeightBins, PulseHeightMin, PulseHeightMax,PulseHeightMaxMeanCharge,PulseHeightMinMeanCharge;
        vector <Float_t> vecPredDetX,vecPredDetY,vecPulseHeight,vecClusterSize;
        vector <Float_t> vecPH_Cluster1,vecPH_Cluster2,vecCh_Cluster2,vecCh_Cluster1;
        Float_t xPred,yPred;
        Float_t xFid,yFid;
        Float_t xChi2, yChi2;

        TH1F* hLandau3DWithoutColumns;
        TH2F* hLandau3DWithoutColumnsFidCutXvsFidCutY;
        TH1F* hLandau3DWithoutColumns_subset;
        TH2F* hLandau3DWithColumnsFidCutXvsFidCutY;
        TH1F* hLandau3DWithColumns;
        TCluster *transparentCluster;
        TCluster clusteredCluster;
        TCluster *diamondCluster;
        vector< vector <Float_t> > vecEdgePredX,vecEdgePredY,vecEdgePulseHeight;
        int verbosity;
    private:
        void FillEdgeAlignmentHistos();
        void HitandSeedCount(TCluster* nCluster);//
        void Analyse1Cluster(UInt_t clusterNo=0);//
        void Analyse2Cluster(); //

        void ClusterPlots(int nClusters);//, float nfiducialValueX, float nfiducialValueY);//exists
        void FillEdgeDistributions(Float_t clusterCharge);//exists
        void FillMeanChargeVector(Float_t clusterCharge);//exists
        void SaveMeanChargeVector(); // exists
        void SaveEdgeDistributions(); //missing
        void Save2ClusterPlots(); //exists
        bool isTransparentCluster;
        bool useCMN;
        int HitCount, SeedCount;
        TTracking* eventReader;
        HistogrammSaver* histSaver;
        TSettings* settings;
        TH2F* hRelativeChargeTwoClustersX;
        TH2F* hRelativeChargeTwoClustersY;
        TProfile2D* hFidCutsVsMeanCharge;
        TProfile2D* hRelativeChargeTwoClustersXY;
        TProfile2D* hShortAnalysis2TotalChargeXY;
        TProfile2D* hTotalAvrgChargeXY;
        TProfile *hRelatviveNumberOfMultipleClusterEvents;
        TProfile *hRelatviveNumberOfMultipleClusterEventsSamePattern;
        TH1F* hNumberofClusters;
        TH1F* hEventsvsChannelCombined;
        TH1F* hDoubleClusterPos;
        TH1F* hDoubleClusterPos0;
        TH1F* hDoubleClusterPos1;
        TH1F* hLandauCluster1;
        TH1F* hLandauCluster2;
        TH1F* hLandauDoubleCombined;
        std::vector<TH1F*> hLandau;
        std::vector<TH1F*> hEventsvsChannel;
        std::vector<TH2F*> hPHvsChannel;
        std::vector<TH2F*> hHitandSeedCount;
        std::vector<TH2F*> hChi2XChi2Y;
        std::vector<TH2F*> hFidCutXvsFidCutY;
        std::vector<TH2D*> hFidCutXvsFidCutYvsCharge;
        std::vector<TH2D*> hFidCutXvsFidCutYvsEvents;
        std::vector<TH2D*> hFidCutXvsFidCutYvsMeanCharge;
        std::vector<TH2D*> hXdetvsYdetvsCharge;
        std::vector<TH2D*> hXdetvsYdetvsEvents;
        std::vector<TH2D*> hXdetvsYdetvsMeanCharge;
        TH2D* hFidCutXvsFidCutYvsMeanChargeAllDetectors;
        std::vector<TH2D*> hFidCutXvsFidCutYClusters;
        TH2F* hShortAnalysis2ClusterHitPattern_1stCluster;
        TH2F* hShortAnalysis2ClusterHitPattern_2ndCluster;
};

#endif /* TANALYSISOF3DSHORTANALYSIS_HH_ */
