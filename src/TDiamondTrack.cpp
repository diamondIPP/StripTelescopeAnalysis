/*
 * TDiamondTrack.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TDiamondTrack.hh"


TDiamondTrack::TDiamondTrack(int event, TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3)
{
   D0 = Det0;
   D1 = Det1;
   D2 = Det2;
   D3 = Det3;
	event_number = event;
	FakeTrack = false;
}

TDiamondTrack::TDiamondTrack(int event, TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3, TDetectorPlane Dia)
{
   D0 = Det0;
   D1 = Det1;
   D2 = Det2;
   D3 = Det3;
   D4 = Dia;
	event_number = event;
	FakeTrack = false;
}


TDiamondTrack::~TDiamondTrack() {}

TDetectorPlane const TDiamondTrack::GetD(Int_t detector)
{
   if(detector == 0) return D0;
   if(detector == 1) return D1;
   if(detector == 2) return D2;
   if(detector == 3) return D3;
   if(detector == 4) return D4;
   else {
      cout<<"TDiamondTrack::GetD("<<detector<<"): detector "<<detector<<" is not a valid detector"<<endl;
      return TDetectorPlane();
   }
}

Float_t TDiamondTrack::GetDetectorHitPosition(Int_t det /*det = 0 to 8*/) {
   switch(det) {
      case 0: return D0.GetX(); break;
      case 1: return D0.GetY(); break;
      case 2: return D1.GetX(); break;
      case 3: return D1.GetY(); break;
      case 4: return D2.GetX(); break;
      case 5: return D2.GetY(); break;
      case 6: return D3.GetX(); break;
      case 7: return D3.GetY(); break;
      case 8: return D4.GetX(); break;
      case 9: return D4.GetY(); break;
   }
   return -1;
}

void TDiamondTrack::SetDetectorHitPosition(Int_t det /*det = 0 to 8*/, Float_t pos) {
   switch(det) {
      case 0: D0.SetX(pos); break;
      case 1: D0.SetY(pos); break;
      case 2: D1.SetX(pos); break;
      case 3: D1.SetY(pos); break;
      case 4: D2.SetX(pos); break;
      case 5: D2.SetY(pos); break;
      case 6: D3.SetX(pos); break;
      case 7: D3.SetY(pos); break;
      case 8: D4.SetX(pos); break;
      case 9: D4.SetY(pos); break;
      }
}
