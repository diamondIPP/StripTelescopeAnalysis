/*
 * THTMLAllignment.cpp
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#include "THTMLAlignment.hh"

THTMLAlignment::THTMLAlignment(TSettings *settings):THTMLGenerator(settings) {

	this->setFileName("alignment.html");
	this->setSubdirPath("alignment");
	this->setTitle("Alignment");
	alignment =0;


}

THTMLAlignment::~THTMLAlignment() {
	// TODO Auto-generated destructor stub
}

void THTMLAlignment::createContent()
{
	createOverviewTable();
	createPostDiamondOverview();
	createPostSiliconOverview();
	createChi2Overview();
	createPreDiamondOverview();
	createPreSiliconOverview();
}



void THTMLAlignment::createPostSiliconOverview()
{

	stringstream sectionContent;

	sectionContent<<"<h3>X Resolution</h3>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaX_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n ";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaX_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n ";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaX_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n ";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaX_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<"<br><br>\n";

	sectionContent<<"<h3>Y Resolution</h3>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaY_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaY_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaY_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_Distribution_DeltaY_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n<br<br>\n";
	sectionContent<<"<br><br>\n";

	sectionContent<<"<h3>X_Obs vs. Y_Obs</h3>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XObs_vs_YObs_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_cut_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XObs_vs_YObs_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_cut_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XObs_vs_YObs_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",24)<<"\n ";//_cut_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XObs_vs_YObs_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",24)<<"\n ";//_cut_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<"<br><br>\n";

	sectionContent<<"<h3>Y_Pred vs DeltaX</h3>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<"<br><br>\n";

	sectionContent<<"<h3>X_Pred vs DeltaY</h3>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",24)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",24)<<"\n";
	sectionContent<<"<br><br>\n";

	this->addSection("Post Alignment Silicon",sectionContent.str());
}



void THTMLAlignment::createPostDiamondOverview()
{
	stringstream sectionContent;
	UInt_t nDiamondAlignmentEvents = alignment!=0?alignment->getDiamondAlignmentEvents():0;
	UInt_t nUsedEvents = alignment!=0?alignment->getNUsedEvents():0;
	Float_t percentage = (Float_t)nDiamondAlignmentEvents/(Float_t)nUsedEvents*100;
	sectionContent<<"For the diamond Alignment "<<nDiamondAlignmentEvents<<" of "<<nUsedEvents <<" ("<<setprecision(4)<<percentage<<"%) fullfill a  Chi2 cut at ";
	sectionContent<<(float)((alignment!=0)?alignment->getDiaChi2():-1.)<<".<br><br>\n";
	sectionContent<<"The diamond is aligned with a digital resoltuion convoluted with a gaus of "<<setprecision(4)<<(float)(alignment!=0?alignment->getXResolution(4):-1)<<" &#956m";
	sectionContent<<" (pure digital resolution: "<<setprecision(4)<<1./TMath::Sqrt(12)*TPlaneProperties::getStripDistance()<<"&#956m)<br><br>\n\n";
	sectionContent<<center(putImage(".",(string)"hDiamond_PostAlignment_Distribution_DeltaX_Plane_4_with_0_1_2_and_3","png",40))<<"<br>\n";
	sectionContent<<putImage(".",(string)"hDiamond_PostAlignment_ScatterPlot_XMeasured_vs_DeltaX_Plane_4_with_0_1_2_and_3","png",33)<<" ";
	sectionContent<<putImage(".",(string)"hDiamond_PostAlignment_ScatterPlot_XPred_vs_DeltaX_Plane_4_with_0_1_2_and_3","png",33)<<" ";
	sectionContent<<putImage(".",(string)"hDiamond_PostAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_4_with_0_1_2_and_3","png",33)<<"<br>\n";
	this->addSection("Post Alignment Diamond",sectionContent.str());
}



void THTMLAlignment::createPreSiliconOverview()
{//
	stringstream sectionContent;
	sectionContent<<"<h3>X Resolution</h3><br>\n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaX_Plane_1_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaX_Plane_2_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaX_Plane_3_with_1_and_2","png",33)<<" \n";
	sectionContent<<" <br<br>\n";
	sectionContent<<"<h3>Y Resolution</h3><br>\n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaY_Plane_1_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaY_Plane_2_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_Distribution_DeltaY_Plane_3_with_1_and_2","png",33)<<" \n";
	sectionContent<<" <br<br>\n";
	sectionContent<<"<h3>XObs vs YObs Range</h3><br>\n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XObs_vs_YObs_Plane_1_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XObs_vs_YObs_Plane_2_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XObs_vs_YObs_Plane_3_with_1_and_2","png",33)<<" \n";
	sectionContent<<" <br<br>\n";
	sectionContent<<"<h3>XPred vs DeltaY</h3><br>\n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_1_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_2_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_XPred_vs_DeltaY_Plane_3_with_1_and_2","png",33)<<" \n";
	sectionContent<<" <br<br>\n";
	sectionContent<<"<h3>YPred vs DeltaX</h3><br>\n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_1_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_2_with_0_and_3","png",33)<<" \n";
	sectionContent<<putImage(".",(string)"hSilicon_PreAlignment_ScatterPlot_YPred_vs_DeltaX_Plane_3_with_1_and_2","png",33)<<" \n";

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
			table.at(plane+1).push_back(center(combineToString((string)" ",(float)(alignment!=0?alignment->getXResolution(plane):0))));
			table.at(plane+1).push_back(center(combineToString((string)" ",(float)(alignment!=0?alignment->getYResolution(plane):0))));
			table.at(plane+1).push_back("");
			table.at(plane+1).push_back(center(combineToString((string)" ",(float)(alignment!=0?alignment->getXMean(plane):0))));
			table.at(plane+1).push_back(center(combineToString((string)" ",(float)(alignment!=0?alignment->getYMean(plane):0))));
		}
		table.at(5).push_back("Diamond");
		table.at(5).push_back(center(combineToString((string)"",alignment!=0?alignment->getXResolution(4):0)));
		table.at(5).push_back(center("--"));
		table.at(5).push_back("");
		table.at(5).push_back(center(combineToString((string)"",alignment!=0?alignment->getXMean(4):0)));
		table.at(5).push_back(center("--"));
	}
	//  sectionContent<<"<h1>Alignment Overview</h1>\n";

	sectionContent<<"<p> Alignent of RUN "<<settings->getRunNumber()<<"<br>\n ";
	sectionContent<<" Made on "<<(string)(alignment!=0?alignment->getLastUpdateTimeAsString():"0 ")<<"<br>\n";
	sectionContent<<" Used "<<(int)(alignment!=0?alignment->getNUsedEvents():-1)<<"Events for the Alignemnt procedure.<br>\n";

	sectionContent<<"chi2 cut: "<<settings->getAlignment_chi2()<<"<br>\n";
	sectionContent<<"Used TrainingMethod '";
	if(settings->getTrainingMethod()==TSettings::enumFraction)sectionContent<<"fraction' with "<<settings->getAlignment_training_track_fraction()*100<<"%. <br>\n";
	else sectionContent<<"number' with "<<settings->getAlignmentTrainingTrackNumber()<<" used Events.<br>\n";
	sectionContent<<this->createTable(table)<<endl;
	this->addSection("Alignment Overview",sectionContent.str());
}

void THTMLAlignment::createPreDiamondOverview()
{
	stringstream sectionContent;
	//  sectionContent<<"<h1>Pre Alignment: Diamond</h1>\n";
	this->addSection("Pre Alignment Diamond",sectionContent.str());
}



void THTMLAlignment::createChi2Overview(){
	stringstream sectionContent;
	sectionContent<<putImage(".",(string)"hPostAlignment_Chi2X_Distribution","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"hPostAlignment_Chi2Y_Distribution","png",45)<<"<br> \n";
	sectionContent<<putImage(".",(string)"hPostAlignment_Chi2X_vs_SumDeltaX","png",45)<<"\n";
	sectionContent<<putImage(".",(string)"hPostAlignment_Chi2Y_vs_SumDeltaY","png",45)<<"\n";
	sectionContent<<"<h3>XPred vs DeltaY</h3><br>\n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_DeltaX_vs_Chi2X_Plane_0_with_1_2_and_3with_Chi2_cut"),"png",33)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",33)<<" \n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_DeltaX_vs_Chi2X_Plane_1_with_0_2_and_3with_Chi2_cut"),"png",33)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",33)<<" \n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_DeltaX_vs_Chi2X_Plane_2_with_0_1_and_3with_Chi2_cut"),"png",33)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",33)<<" \n";
	sectionContent<<putImage(".",TString::Format("hSilicon_PostAlignment_ScatterPlot_DeltaX_vs_Chi2X_Plane_3_with_0_1_and_2with_Chi2_cut"),"png",33)<<"\n ";//_on_%1.f",settings->getAlignment_chi2()),"png",33)<<" \n";
	sectionContent<<" <br<br>\n";

	this->addSection("Chi2 Distributions",sectionContent.str());
}
