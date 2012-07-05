/*
 * THTMLAllignment.cpp
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#include "THTMLAlignment.hh"

THTMLAlignment::THTMLAlignment(TSettings *settings):THTMLGenerator(settings) {

  this->setFileName("alignment.html");
  this->setSubdirPath("alignment/");
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
  vector< vector< string> > table;
  table.resize(6);
  table.at(0).push_back("");
  table.at(0).push_back("Res X");
  table.at(0).push_back("Res Y");
  table.at(0).push_back("Mean X");
  table.at(0).push_back("Mean Y");

  for(UInt_t plane=0;plane<TPlaneProperties::getNSiliconPlanes();plane++){
    table.at(plane+1).push_back(combineToString((string)"plane ",plane));
    table.at(plane+1).push_back(combineToString((string)"",alignment->getXResolution(plane)));
    table.at(plane+1).push_back(combineToString((string)"",alignment->getYResolution(plane)));
    table.at(plane+1).push_back(combineToString((string)"",alignment->getXMean(plane)));
    table.at(plane+1).push_back(combineToString((string)"",alignment->getYMean(plane)));
  }
  table.at(5).push_back("Diamond");
  table.at(5).push_back(combineToString((string)"",alignment->getXResolution(4)));
  table.at(5).push_back("");
  table.at(5).push_back(combineToString((string)"",alignment->getXMean(4)));
  table.at(5).push_back("");
  sectionContent<<"<h1>Alignment Overview</h1>\n";
  sectionContent<<this->createTable(table)<<endl;
  this->addSection("Alignment Overview",sectionContent.str());
}

void THTMLAlignment::createPreDiamondOverview()
{
  stringstream sectionContent;
  sectionContent<<"<h1>Pre Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



