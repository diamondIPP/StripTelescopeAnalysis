/*
 * THTMLAllignment.cpp
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#include "THTML3DAnalysis.hh"

THTML3DAnalysis::THTML3DAnalysis(TSettings *settings):THTMLGenerator(settings) {

	this->setFileName("3dAnalysis.html");
	this->setSubdirPath("3dDiamondAnalysis/");
	this->setTitle("3DAnalysis");


}

THTML3DAnalysis::~THTML3DAnalysis() {
	// TODO Auto-generated destructor stub
}

void THTML3DAnalysis::createContent()
{
    createOverviewPlots();
    createDistributionComparePlots();
    createResolutionPlots();
}

void THTML3DAnalysis::createOverviewPlots(){
	stringstream sectionContent;
	sectionContent<<putImage(".",(string)"cAvrgPulseHeightDetSystem_MetalizationLayer","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"cAvrgPulseHeightDetSystem_MetalizationLayer_Zoom_rebinned","png",45)<<"<br> \n";
	sectionContent<<putImage(".",(string)"cProfRebinned","png",45)<<"\n";
	sectionContent<<putImage(".",(string)"cPulseHeightVsDetectorHitPostionXY_rebinned","png",45)<<"\n";
    sectionContent<<putImage(".",(string)"cPulseHeightVsDetectorHitPostionXY_rebinned_trans","png",45)<<"\n";

	sectionContent<<putImage(".",(string)"cPulseHeightVsDetectorHitPostionXYGoodCells","png",45)<<"\n";
	sectionContent<<putImage(".",(string)"cTotalAvrgChargeXY","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"cTotalAvrgChargeXY_trans","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"hPulseHeightVsDetectorHitPostionXY_rebinned","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"hPulseHeightVsDetectorHitPostionXY_rebinned_trans","png",45)<<" \n";

	sectionContent<<putImage(".",(string)"hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"hTransparentAnalysisTransparentSmallChargeInFirstChannelGoodCells_grid","png",45)<<" \n";
	sectionContent<<putImage(".",(string)"chNegativeChargePositionGrid_trans","png",45)<<" \n";
	sectionContent<<" <br<br>\n";

	this->addSection("Overview Plots",sectionContent.str());
}

void THTML3DAnalysis::createDistributionComparePlots(){

    stringstream sectionContent;
    sectionContent<<putImage(".",(string)"cNoNegativeCharge_PulseHeight_Comparision_good","png",45)<<" \n";
    sectionContent<<putImage(".",(string)"c_stackPulseHeights_GoodCells_good","png",45)<<" \n";
    sectionContent<<putImage(".",(string)"c_stack_NoNegativeCharge_PulseHeight_Comparision_good","png",45)<<" \n";
//
    sectionContent<<putImage(".",(string)"c_sAllPulseHeightDistributions_scaled","png",45)<<" \n";
    sectionContent<<putImage(".",(string)"c_sAllPulseHeightDistributions_scaled_trans","png",45)<<" \n";
    sectionContent<<" <br<br>\n";
    this->addSection("Distribution Compare",sectionContent.str());
}

void THTML3DAnalysis::createResolutionPlots(){
    stringstream sectionContent;
    sectionContent<<putImage(".",(string)"hResolutionAllCells_maxValue_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionAllCells_highest2Centroid_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionAllCells_chargeWeighted_trans","png",32)<<" <br>\n";
    sectionContent<<putImage(".",(string)"hResolutionAllButBadCells_maxValue_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionAllButBadCells_highest2Centroid_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionAllButBadCells_chargeWeighted_trans","png",32)<<" <br>\n";
    sectionContent<<putImage(".",(string)"hResolutionGoodCells_maxValue_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionGoodCells_highest2Centroid_trans","png",32)<<" \n";
    sectionContent<<putImage(".",(string)"hResolutionGoodCells_chargeWeighted_trans","png",32)<<" \n";
    sectionContent<<" <br<br>\n";
    this->addSection("Resolution Plots",sectionContent.str());
}
