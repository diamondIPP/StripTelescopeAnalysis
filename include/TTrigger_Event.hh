/*
 * TTrigger_Event.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TTRIGGER_EVENT_HH_
#define TTRIGGER_EVENT_HH_
#include "TMath.h"
#include <iostream>
#include "TDetector_Data.hh"

using namespace std;

class TTrigger_Event {

   public:
      TTrigger_Event();
      virtual ~TTrigger_Event();
      //TDetector_Data GetDetector(TDetector_Data Detector);
      void SetD0X(TDetector_Data Detector);
      void SetD0Y(TDetector_Data Detector);
      void SetD1X(TDetector_Data Detector);
      void SetD1Y(TDetector_Data Detector);
      void SetD2X(TDetector_Data Detector);
      void SetD2Y(TDetector_Data Detector);
      void SetD3X(TDetector_Data Detector);
      void SetD3Y(TDetector_Data Detector);
      void SetDia0(TDetector_Data Detector);
      void SetDia1(TDetector_Data Detector);
      void SetAny(TDetector_Data Detector);
      TDetector_Data GetD0X() const {return D0X;}
      TDetector_Data GetD0Y() const {return D0Y;}
      TDetector_Data GetD1X() const {return D1X;}
      TDetector_Data GetD1Y() const {return D1Y;}
      TDetector_Data GetD2X() const {return D2X;}
      TDetector_Data GetD2Y() const {return D2Y;}
      TDetector_Data GetD3X() const {return D3X;}
      TDetector_Data GetD3Y() const {return D3Y;}
      TDetector_Data GetDia0() const {return Dia0;}
      TDetector_Data GetDia1() const {return Dia1;}
      TDetector_Data GetAny() const {return Any;}

   protected:
      TDetector_Data D0X;
      TDetector_Data D0Y;
      TDetector_Data D1X;
      TDetector_Data D1Y;
      TDetector_Data D2X;
      TDetector_Data D2Y;
      TDetector_Data D3X;
      TDetector_Data D3Y;
      TDetector_Data Dia0;
      TDetector_Data Dia1;
      TDetector_Data Any;

};

#endif /* TTRIGGER_EVENT_HH_ */
