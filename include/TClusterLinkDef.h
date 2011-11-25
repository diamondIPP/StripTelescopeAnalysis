#ifdef __CINT__
#include <deque>
#include <pair>
#include <vector>
#pragma link C++ class pair< UShort_t, Float_t>+;
#pragma link C++ class pair< int, Float_t>+;
#pragma link C++ class deque< pair<int,Float_t> >+;
#pragma link C++ class deque< pair< UShort_t, Float_t> >+;
#pragma link C++ class TCluster;
#pragma link C++ class vector<TCluster>;
#pragma link C++ class vector< vector<TCluster> >;
#pragma link C++ class vector<TCluster*>;
#pragma link C++ class vector< vector<TCluster*> >;
#endif
