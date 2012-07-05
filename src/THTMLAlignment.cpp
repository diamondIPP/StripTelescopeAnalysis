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
//  sectionContent<<"<h1>Post Alignment: Silicon</h1>\n";
  this->addSection("Post Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createPostDiamondOverview()
{
  stringstream sectionContent;
//  sectionContent<<"<h1>Post Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



void THTMLAlignment::createPreSiliconOverview()
{
  stringstream sectionContent;
//  sectionContent<<"<h1>Pre Alignment: Silicon</h1>\n";
  this->addSection("Pre Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createOverviewTable()
{
  stringstream sectionContent;
  vector< vector< string> > table;
  table.resize(6);
  table.at(0).push_back("");
  table.at(0).push_back("Res X [&#956m]");
  table.at(0).push_back("Res Y [&#956m]");
  table.at(0).push_back("Mean X [&#956m]");
  table.at(0).push_back("Mean Y [&#956m]");
  if(alignment!=0){
    for(UInt_t plane=0;plane<TPlaneProperties::getNSiliconPlanes();plane++){
      table.at(plane+1).push_back(center(combineToString((string)"Plane ",plane)));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getXResolution(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getYResolution(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getXMean(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getYMean(plane)*TPlaneProperties::getStripDistance())));
    }
    table.at(5).push_back("Diamond");
    table.at(5).push_back(center(combineToString((string)"",alignment->getXResolution(4)*TPlaneProperties::getStripDistance())));
    table.at(5).push_back(center("--"));
    table.at(5).push_back(center(combineToString((string)"",alignment->getXMean(4)*TPlaneProperties::getStripDistance())));
    table.at(5).push_back(center("--"));
  }
//  sectionContent<<"<h1>Alignment Overview</h1>\n";
  sectionContent<<this->createTable(table)<<endl;
  this->addSection("Alignment Overview",sectionContent.str());
}

void THTMLAlignment::createPreDiamondOverview()
{
  stringstream sectionContent;
//  sectionContent<<"<h1>Pre Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



