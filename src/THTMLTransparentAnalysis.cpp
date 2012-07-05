/*
 * THTMLTransparentAnalysis.cpp
 *
 *  Created on: Jul 5, 2012
 *      Author: bachmair
 */

#include "THTMLTransparentAnalysis.hh"

THTMLTransparentAnalysis::THTMLTransparentAnalysis(TSettings* settings):THTMLGenerator(settings) {
  // TODO Auto-generated constructor stub
  this->setFileName("transparentAnalysis.html");
  this->setMainPath("../");
  this->setSubdirPath("transparentAnalysis/");
  this->setTitle("Transparent Analysis");


}

THTMLTransparentAnalysis::~THTMLTransparentAnalysis() {
  // TODO Auto-generated destructor stub
}

