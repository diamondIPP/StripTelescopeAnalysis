/*
 * THTMLAllignment.cpp
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#include "THTMLAlignment.hh"

THTMLAlignment::THTMLAlignment(TSettings *settings):THTMLGenerator(settings) {

  this->setFileName("alignment.html");
  this->setSubdirPath("/alignment/");
  this->setTitle("Alignment");


}

THTMLAlignment::~THTMLAlignment() {
	// TODO Auto-generated destructor stub
}

