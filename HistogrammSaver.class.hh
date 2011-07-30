/*
 * HistogrammSaver.class.hh
 *
 *  Created on: 29.07.2011
 *      Author: Felix Bachmair
 */

#ifndef HISTOGRAMMSAVER_CLASS_HH_
#define HISTOGRAMMSAVER_CLASS_HH_
#include "TH1F.h"
#include "TCanvas.h"
class HistogrammSaver {
public:
	HistogrammSaver();
	virtual ~HistogrammSaver();
    void SaveHistogram(TH1F* histo);
    void SaveHistogram(TH2F* histo);
    void SaveHistogramPNG(TH1F* histo);
    void SaveHistogramPNG(TH2F* histo);
    void SaveHistogramROOT(TH1F* histo);
    void SaveHistogramROOT(TH2F* histo);
    void SaveHistogramPDF(TH1F* histo);
    void SaveHistogramPDF(TH2F* histo);
};

#endif /* HISTOGRAMMSAVER_CLASS_HH_ */
