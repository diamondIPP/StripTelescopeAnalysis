/*
 * TAvrgChargePerBinMonteCarlo.hh
 *
 *  Created on: Sep 17, 2015
 *      Author: bachmair
 */

#ifndef TAVRGCHARGEPERBINMONTECARLO_HH_
#define TAVRGCHARGEPERBINMONTECARLO_HH_
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
#include "TCluster.hh"

class TAvrgChargePerBinMonteCarlo {
    public:
        TAvrgChargePerBinMonteCarlo(HistogrammSaver *histSaver);
        virtual ~TAvrgChargePerBinMonteCarlo();
        void doMonteCarlo(
                TProfile2D* profOverlay, TH1F* hLandauOfOverlay);
    private:
        void init_histos(TProfile2D* profOverlay);
        void save_histos();
        void do_loop();
        TH1D* hMonteCarloAvrgChargePerBin;
        TH1D* hNumberOfEntriesBelowCut;
        TH1D* hRelativeNumberBelowCut;
        TH1F* hEntries;
        TH2D* hBinContents;
        TH1F* hLandauOfOverlay;
    public:

        Float_t getCut() const {
            return cut;
        }

        void setCut(Float_t cut) {
            this->cut = cut;
        }

        Int_t getPulseHeightBins() const {
            return PulseHeightBins;
        }

        void setPulseHeightBins(Int_t pulseHeightBins) {
            PulseHeightBins = pulseHeightBins;
        }

        Int_t getPulseHeightMax() const {
            return PulseHeightMax;
        }

        void setPulseHeightMax(Int_t pulseHeightMax) {
            PulseHeightMax = pulseHeightMax;
        }

        Int_t getPulseHeightMaxMeanCharge() const {
            return PulseHeightMaxMeanCharge;
        }

        void setPulseHeightMaxMeanCharge(Int_t pulseHeightMaxMeanCharge) {
            PulseHeightMaxMeanCharge = pulseHeightMaxMeanCharge;
        }

        Int_t getPulseHeightMin() const {
            return PulseHeightMin;
        }

        void setPulseHeightMin(Int_t pulseHeightMin) {
            PulseHeightMin = pulseHeightMin;
        }

        Int_t getPulseHeightMinMeanCharge() const {
            return PulseHeightMinMeanCharge;
        }

        void setPulseHeightMinMeanCharge(Int_t pulseHeightMinMeanCharge) {
            PulseHeightMinMeanCharge = pulseHeightMinMeanCharge;
        }

    private:
        Float_t cut;
        HistogrammSaver *histSaver;
        Int_t PulseHeightBins, PulseHeightMin, PulseHeightMax,PulseHeightMaxMeanCharge,PulseHeightMinMeanCharge;
};

#endif /* TAVRGCHARGEPERBINMONTECARLO_HH_ */
