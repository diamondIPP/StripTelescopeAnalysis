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
  sectionContent<<"<h3>X Range</h3><br>\n";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaX_-_Plane_0_with_1_2_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaX_-_Plane_1_with_0_2_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaX_-_Plane_2_with_0_1_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaX_-_Plane_3_with_0_1_and_2with_Chi2_cut_on_1","png",24)<<" <br<br>\n";
  sectionContent<<"<h3>X Range</h3><br>\n";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaY_-_Plane_0_with_1_2_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaY_-_Plane_1_with_0_2_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaY_-_Plane_2_with_0_1_and_3with_Chi2_cut_on_1","png",24)<<" ";
  sectionContent<<putImage(".","hSilicon_PostAlignment_DistributionPlot_DeltaY_-_Plane_3_with_0_1_and_2with_Chi2_cut_on_1","png",24)<<" <br<br>\n";

  this->addSection("Post Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createPostDiamondOverview()
{
  stringstream sectionContent;
  UInt_t nDiamondAlignmentEvents = alignment->getDiamondAlignmentEvents();
  UInt_t nUsedEvents = alignment->getNUsedEvents();
  Float_t percentage = (Float_t)nDiamondAlignmentEvents/(Float_t)nUsedEvents*100;
  sectionContent<<"For the diamond Alignment "<<nDiamondAlignmentEvents<<" of "<<nUsedEvents <<" ("<<setprecision(2)<<percentage<<"%) fullfill a  Chi2 cut at "<<alignment->getDiaChi2()<<".<br><br>\n";
  sectionContent<<"The diamond is aligned with a digital resoltuion convoluted with a gaus of "<<setprecision(2)<<alignment->getXResolution(4)*TPlaneProperties::getStripDistance()<<" &#956m";
  sectionContent<<" (pure digital resolution: "<<setprecision(2)<<1./TMath::Sqrt(12)*TPlaneProperties::getStripDistance()<<"&#956m)<br><br>\n\n";
  sectionContent<<center(putImage(".","hDiamond_PostAlignment_DistributionPlot_DeltaX_-_Plane_4_with_0_1_2_and_3","png",40))<<"<br>\n";
  sectionContent<<putImage(".","hDiamond_PostAlignment_ScatterPlot_XMeasured_vs_DeltaX_-_Plane_4_with_0_1_2_and_3","png",30)<<" ";
  sectionContent<<putImage(".","hDiamond_PostAlignment_ScatterPlot_XPred_vs_DeltaX_-_Plane_4_with_0_1_2_and_3","png",30)<<" ";
  sectionContent<<putImage(".","hDiamond_PostAlignment_ScatterPlot_YPred_vs_DeltaX_-_Plane_4_with_0_1_2_and_3","png",30)<<"<br>\n";
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
  table.at(0).push_back(" Res X [&#956m] ");
  table.at(0).push_back(" Res Y [&#956m] ");
  table.at(0).push_back(" ");
  table.at(0).push_back(" Mean X [&#956m] ");
  table.at(0).push_back(" Mean Y [&#956m] ");
  if(alignment!=0){
    for(UInt_t plane=0;plane<TPlaneProperties::getNSiliconPlanes();plane++){
      table.at(plane+1).push_back(center(combineToString((string)"Plane ",plane)));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getXResolution(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getYResolution(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back("");
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getXMean(plane)*TPlaneProperties::getStripDistance())));
      table.at(plane+1).push_back(center(combineToString((string)" ",alignment->getYMean(plane)*TPlaneProperties::getStripDistance())));
    }
    table.at(5).push_back("Diamond");
    table.at(5).push_back(center(combineToString((string)"",alignment->getXResolution(4)*TPlaneProperties::getStripDistance())));
    table.at(5).push_back(center("--"));
    table.at(5).push_back("");
    table.at(5).push_back(center(combineToString((string)"",alignment->getXMean(4)*TPlaneProperties::getStripDistance())));
    table.at(5).push_back(center("--"));
  }
//  sectionContent<<"<h1>Alignment Overview</h1>\n";

  sectionContent<<"<p> Alignent of RUN "<<settings->getRunNumber()<<"<br>\n ";
  sectionContent<<" Made on "<<alignment->getLastUpdateTimeAsString()<<"<br>\n";
  sectionContent<<" Used "<<alignment->getNUsedEvents()<<"Events for the Alignemnt procedure.<br>\n";
  sectionContent<<this->createTable(table)<<endl;
  this->addSection("Alignment Overview",sectionContent.str());
}

void THTMLAlignment::createPreDiamondOverview()
{
  stringstream sectionContent;
//  sectionContent<<"<h1>Pre Alignment: Diamond</h1>\n";
  this->addSection("Post Alignment Diamond",sectionContent.str());
}



