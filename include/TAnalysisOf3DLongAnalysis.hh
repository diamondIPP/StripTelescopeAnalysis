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

#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"

class TAnalysisOf3D_LongAnalysis {
    public:
        TAnalysisOf3D_LongAnalysis(TSettings *settings,HistogrammSaver *histSaver);
        virtual ~TAnalysisOf3D_LongAnalysis();
        void addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid);
        void initHistos();
        void saveHistos();
    private:
        HistogrammSaver* histSaver;
        TSettings* settings;

};

#endif /* TANALYSISOF3DLONGANALYSIS_HH_ */
