#include <fstream>
#include <iostream>

bool RunListOK;
std::string inputDir;
std::string outputDir;
std::string runListPath="RunList.ini";

int ReadRunList();

class RunInfo
{
public:
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
public:
    void setParameters(UInt_t nRunNo,std::string sRunDes,UInt_t nVeb,UInt_t NEvents,UInt_t nStartEvent,bool bPedAna,bool bClusAna,bool bSelAna,bool bAlign,bool bAlignAna){
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
    }
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
};
std::vector<RunInfo> RunParameters;
