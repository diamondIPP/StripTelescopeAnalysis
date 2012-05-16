/*
 * THTMLSelection.cpp
 *
 *  Created on: May 16, 2012
 *      Author: bachmair
 */

#include "../include/THTMLSelection.hh"

THTMLSelection::THTMLSelection(TSettings *settings):THTMLGenerator(settings) {
	// TODO Auto-generated constructor stub
	this->setFileName("selection.html");
	this->setSubdirPath("/selections/");
	this->setTitle("Selection - Cut Flow");
}

THTMLSelection::~THTMLSelection() {
	// TODO Auto-generated destructor stub
}

void THTMLSelection::createCutFlowTable(std::vector<int> vecCutFlow	)
{
}



