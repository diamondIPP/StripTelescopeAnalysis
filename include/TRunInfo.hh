/*
 * TRunInfo.hh
 *
 *  Created on: Jul 17, 2012
 *      Author: bachmair
 */

#ifndef TRUNINFO_H_
#define TRUNINFO_H_
#include <fstream>
#include <iostream>
#include "TROOT.h"
using namespace std;
class TRunInfo {
public:
  TRunInfo();
  virtual ~TRunInfo();
private:
  UInt_t nRunNumber;
  UInt_t nVerbosity;
  UInt_t nEvents;
  UInt_t nStartEvent;
  std::string RunDescription;
  bool bPedestalAnalysis;
  bool bClusterAnalysis;
  bool bSelectionAnalysis;
  bool bAlignment;
  bool bAlignmentAnalysis;
  bool bTransparentAnalysis;
public:
  void setParameters(UInt_t nRunNo,std::string sRunDes,UInt_t nVeb,UInt_t NEvents,UInt_t nStartEvent,bool bPedAna,bool bClusAna,bool bSelAna,bool bAlign,bool bAlignAna,bool bTransAna){
    setRunNumber(nRunNo);
    setRunDescription(sRunDes);
    setVerbosity(nVeb);
    setEvents(NEvents);
    setStartEvent(nStartEvent);
    setPedestalAnalysis(bPedAna);
    setClusterAnalysis(bClusAna);
    setSelectionAnalysis(bSelAna);
    setAlignment(bAlign);
    setAlignmentAnalysis(bAlignAna);
    setTransparentAnalysis(bTransAna);
  }
  void setTransparentAnalysis(bool bTransAna){bTransparentAnalysis=bTransAna;}
  void setRunDescription(std::string rundescribtion){RunDescription=rundescribtion;}
  void setAlignment(bool alignment){ bAlignment = alignment;}
  void setAlignmentAnalysis(bool alignmentAnalysis){ bAlignmentAnalysis = alignmentAnalysis;}
  void setClusterAnalysis(bool clusterAnalysis){ bClusterAnalysis = clusterAnalysis;}
  void setEvents(UInt_t events){ nEvents = events;}
  void setPedestalAnalysis(bool pedestalAnalysis){ bPedestalAnalysis = pedestalAnalysis;}
  void setRunNumber(UInt_t runNumber){ nRunNumber = runNumber;}
  void setSelectionAnalysis(bool selectionAnalysis){ bSelectionAnalysis = selectionAnalysis;}
  void setStartEvent(UInt_t startEvent){ nStartEvent = startEvent;}
  void setVerbosity(UInt_t verbosity){ nVerbosity = verbosity;}

  UInt_t getEvents() const{return nEvents;}
  UInt_t getRunNumber() const{return nRunNumber;}
  UInt_t getStartEvent() const{return nStartEvent;}
  UInt_t getVerbosity() const{return nVerbosity;}
  std::string getRunDescription() const{return RunDescription;}
  bool doAlignment() const{return bAlignment;}
  bool doAlignmentAnalysis() const{return bAlignmentAnalysis;}
  bool doClusterAnalysis() const{return bClusterAnalysis;}
  bool doPedestalAnalysis() const{return bPedestalAnalysis;}
  bool doSelectionAnalysis() const{return bSelectionAnalysis;}
  bool doTransparentAnalysis() const{return bTransparentAnalysis;}
  void Print(){
    cout << endl << endl << endl << endl;
    cout << "====================================" << endl;
    cout << "==> Settings of RUN "<<getRunNumber()<<".." << endl;
    cout << "====================================" << endl << endl;
    cout << "NEVENTS: " << getEvents() << endl;
    cout << "RUNDESCRIPTION: " << getRunDescription() << endl;
    cout << "VERBOSITY: " << getVerbosity() << endl;
    cout << "INITIAL_EVENT: " << getStartEvent() << endl;
    cout << "DO_PEDESTALANALYSIS: "<<doPedestalAnalysis()<<endl;
    cout << "DO_CLUSTERANALYSIS: "<<doClusterAnalysis()<<endl;
    cout << "DO_SELECTIONANALYSIS: "<<doSelectionAnalysis()<<endl;
    cout << "DO_ALIGNMENT: " << doAlignment() << endl;
    cout << "DO_ALIGNMENTANALYSIS: "<<doAlignmentAnalysis()<<endl;
    cout << "DO_TRANSPARENTANALYSIS: "<<doTransparentAnalysis()<<endl;
    cout << endl << endl;
  };
};

#endif /* TRUNINFO_H_ */
