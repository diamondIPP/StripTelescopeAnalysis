/*
 * TAnalysisOf3DLongAnalysis.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: bachmair
 */

#include <TAnalysisOf3DLongAnalysis.hh>

TAnalysisOf3D_LongAnalysis::TAnalysisOf3D_LongAnalysis(TSettings* settings,HistogramSaver* histSaver) {
    // TODO Auto-generated constructor stub
    this->histSaver = histSaver;
    this->settings = settings;

}

TAnalysisOf3D_LongAnalysis::~TAnalysisOf3D_LongAnalysis() {
    // TODO Auto-generated destructor stub
}

void TAnalysisOf3D_LongAnalysis::addEvent(TCluster* cluster, Float_t x_pred, Float_t y_pred, Float_t x_fid, Float_t y_fid){

}

void TAnalysisOf3D_LongAnalysis::initHistos() {
}

void TAnalysisOf3D_LongAnalysis::saveHistos() {
}
