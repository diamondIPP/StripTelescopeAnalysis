/*
 * TAnalysisOf3DCellLandaus.hh
 *
 *  Created on: Sep 22, 2015
 *      Author: bachmair
 */

#ifndef TANALYSISOF3DCELLLANDAUS_HH_
#define TANALYSISOF3DCELLLANDAUS_HH_
#include "TH1F.h"
#include <vector>
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"
class TAnalysisOf3D_CellLandaus {
    public:
        TAnalysisOf3D_CellLandaus(TSettings *settings,HistogrammSaver *histSaver);
        virtual ~TAnalysisOf3D_CellLandaus();
        TH1F* getCellLandau(unsigned int cell);
        void setLandauStrip(TH1F* histo){hLandauStrip=histo;}
        void addEvent(int cellNo,int quarterNo,float charge);
    private:
        TSettings *settings;
        HistogrammSaver *histSaver;
        std::vector<TH1F*> hCellsLandaus;
        std::vector< std::vector<TH1F*> > hQuarterCellsLandau;
        std::vector< std::vector<TH1F*> > hQuarterCellsClusterSize;
        void initQuarterCellLandaus();
        void LongAnalysis_InitGoodCellsLandaus();
        void LongAnalysis_FillGoodCellsLandaus(Float_t charge);
        void LongAnalysis_SaveGoodCellsLandaus();
        void SaveGoodAndBadCellLandaus();
        int PulseHeightBins;
        float PulseHeightMin,PulseHeightMax;
        TH1F* hLandauStrip;
        unsigned int verbosity;
};

#endif /* TANALYSISOF3DCELLLANDAUS_HH_ */
