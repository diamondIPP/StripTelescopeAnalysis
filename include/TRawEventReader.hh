/*
 * TRawEventReader.hh
 *
 *  Created on: 01.11.2011
 *      Author: bachmair
 */

#ifndef TRAWEVENTREADER_HH_
#define TRAWEVENTREADER_HH_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include <bitset>

//ROOT Class Headers
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TMath.h"
using namespace TMath;
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TImage.h"
//#include "TROOT.h"
#include "TDatime.h"
#include "TObjArray.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TText.h"


#include "RZEvent.struct.hh" //the header file that is connected to the Diamond/telescope data
#include "TDetector_Data.hh"
#include "TSettings.class.hh"
#include "TSystem.h"
using namespace std;

class TRawEventReader {
public:
	TRawEventReader(TSettings* settings);
	virtual ~TRawEventReader();
	int ReadRawEvent(int EventNumber, bool verbose = 0);
	TDetector_Data getPlane(int det,UInt_t diaInput);
    TDetector_Data getD0X() const;
    void setD0X(TDetector_Data d0X);
    TDetector_Data getD0Y() const;
    TDetector_Data getD1X() const;
    TDetector_Data getD1Y() const;
    TDetector_Data getD2X() const;
    TDetector_Data getD2Y() const;
    TDetector_Data getD3X() const;
    TDetector_Data getD3Y() const;
    TDetector_Data getDia0() const;
    TDetector_Data getDia1() const;
    TDetector_Data getDia(UInt_t diaInput) const;
    void setD0Y(TDetector_Data d0Y);
    void setD1X(TDetector_Data d1X);
    void setD1Y(TDetector_Data d1Y);
    void setD2X(TDetector_Data d2X);
    void setD2Y(TDetector_Data d2Y);
    void setD3X(TDetector_Data d3X);
    void setD3Y(TDetector_Data d3Y);
    void setDia0(TDetector_Data dia0);
    void setDia1(TDetector_Data dia1);
private:
    TSettings *settings;
    Int_t EventsPerFile;
    std::string current_rz_filename;
    ifstream current_rz_file;
    int run_number;
    RZEvent rzEvent;
    int verbosity;
public://TODO: Get them in the private area...
    TDetector_Data D0X;
    TDetector_Data D1X;
    TDetector_Data D2X;
    TDetector_Data D3X;
    TDetector_Data D0Y;
    TDetector_Data D1Y;
    TDetector_Data D2Y;
    TDetector_Data D3Y;
    TDetector_Data Dia0;
    TDetector_Data Dia1;
private:
    //This function swaps the endianess of the read in data for a signed 32-bit integer.
    void endian_swap(int & x)
    {
        x = (x >> 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x << 24);
    }

    //This function is overloaded to swap the endianness of the read in data for an unsigned 32-bit integer.
    void uendian_swap(unsigned int & x)
    {
        x = (x >> 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x << 24);
    }

    //This function swaps the endianess of the read in data for a signed 16-bit integer.
    void short_endian_swap(short int & x)
    {
        x = (x >> 8) | (x << 8);
    }

    //This function swaps the endianness of the read in data for an unsigned 16-bit integer.
    void ushort_endian_swap(unsigned short int & x)
    {
        x = (x >> 8) | (x << 8);
    }

};

#endif /* TRAWEVENTREADER_HH_ */
