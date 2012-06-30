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

void THTMLAlignment::createContent()
{
  createOverviewTable();
  createPostDiamondOverview();
  createPostSiliconOverview();
  createPreDiamondOverview();
  createPreSiliconOverview();
}



void THTMLAlignment::createPostSiliconOverview()
{

  stringstream sectionContent;
  sectionContent<<"<h1>Post Alignment: Silicon</h1>\n";
  this->addSection("Post Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createPostDiamondOverview()
{
  stringstream sectionContent;
  sectionContent<<"<h1>Post Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



void THTMLAlignment::createPreSiliconOverview()
{
  stringstream sectionContent;
  sectionContent<<"<h1>Pre Alignment: Silicon</h1>\n";
  this->addSection("Pre Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createOverviewTable()
{
  stringstream sectionContent;
    sectionContent<<"<h1>Alignment Overview</h1>\n";
    this->addSection("Alignment Overview",sectionContent.str());
}

void THTMLAlignment::createPreDiamondOverview()
{
  stringstream sectionContent;
  sectionContent<<"<h1>Pre Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



