#ifdef __CINT__
#include <deque>
#include <utility>
#include <vector>
#include <map>

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class deque<bool>+;

#pragma link C++ class vector<bool>+;
#pragma link C++ class vector<Double_t>+;
#pragma link C++ class vector<Float_t>+;

#pragma link C++ class map<UInt_t,UInt_t>+;
#pragma link C++ class map<TString,Float_t>+;
#pragma link C++ class map<TString,Int_t>+;
#pragma link C++ class map<TString,TString>+;
#pragma link C++ class map<TString, map<TString,TString> >+;

#pragma link C++ class pair< unsigned short, float>+;
#pragma link C++ class pair< UShort_t, Float_t>+;
#pragma link C++ class pair< int, Float_t>+;
#pragma link C++ class pair<UInt_t,UInt_t>+;
#pragma link C++ class pair<Float_t ,Float_t>+;

#pragma link C++ class deque< pair<int,Float_t> >+;
#pragma link C++ class deque< pair< UShort_t, Float_t> >+;
#pragma link C++ class vector< pair< Float_t, Float_t> > +;

#pragma link C++ class TResults+;
#pragma link C++ class TCluster+;
#pragma link C++ class vector<TCluster>+;
#pragma link C++ class vector< vector<TCluster> >+;
#pragma link C++ class vector<TCluster*>+;
#pragma link C++ class vector< vector<TCluster*> >+;
#pragma link C++ class TPlane+;
#pragma link C++ class vector<TPlane>+;
#pragma link C++ class TEvent+;
#pragma link C++ class TPlaneProperties+;
#pragma link C++ class TDiamondPattern+;
#pragma link C++ class TChannelMapping+;
#pragma link C++ class ChannelScreen+;
#pragma link C++ class TSettings+;
#pragma link C++ class TFiducialCut+;

#pragma link C++ class vector<TFiducialCut*>+;
#pragma link C++ class TFidCutRegions+;

#pragma link C++ class TDetectorAlignment+;
#endif
