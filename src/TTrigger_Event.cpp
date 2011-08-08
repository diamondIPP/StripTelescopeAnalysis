/*
 * TTrigger_Event.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TTrigger_Event.hh"
using namespace std;


TTrigger_Event::TTrigger_Event() {
}

TTrigger_Event::~TTrigger_Event() {
}

void TTrigger_Event::SetD0X(TDetector_Data Detector)
{
   D0X = Detector;
}
void TTrigger_Event::SetD0Y(TDetector_Data Detector)
{
   D0Y = Detector;
}
void TTrigger_Event::SetD1X(TDetector_Data Detector)
{
   D1X = Detector;
}
void TTrigger_Event::SetD1Y(TDetector_Data Detector)
{
   D1Y = Detector;
}
void TTrigger_Event::SetD2X(TDetector_Data Detector)
{
   D2X = Detector;
}
void TTrigger_Event::SetD2Y(TDetector_Data Detector)
{
   D2Y = Detector;
}
void TTrigger_Event::SetD3X(TDetector_Data Detector)
{
   D3X = Detector;
}
void TTrigger_Event::SetD3Y(TDetector_Data Detector)
{
   D3Y = Detector;
}
void TTrigger_Event::SetDia0(TDetector_Data Detector)
{
   Dia0 = Detector;
}
void TTrigger_Event::SetDia1(TDetector_Data Detector)
{
   Dia1 = Detector;
}
void TTrigger_Event::SetAny(TDetector_Data Detector)
{
   Any = Detector;
}
