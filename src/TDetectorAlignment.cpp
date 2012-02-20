/*
 * TDetectorAlignment.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TDetectorAlignment.hh"
using namespace std;
ClassImp(TDetectorAlignment);
/*
 *   XY		 YX					   XY		YX
 *   ||      ||				       ||       ||
 *   01   	 23					   45	    67
 *   0		  1						2 		3
 *   D0Z	D1Z
 *   0
 *
 *   |||||
 *   -----
 *   PLANE1: X OFF,YOFF,PHIX,PHIY
 */

TDetectorAlignment::TDetectorAlignment(){
	this->SetName("currentAlignment");
	this->SetTitle("currentAlignment");
	verbosity=0;
   nDetectors = 5;

   for(Int_t i=0; i<5; i++) {
      det_x_offset[i] = 0;
      det_y_offset[i] = 0;
      det_z_offset[i] = 0;
      det_phix_offset[i] = 0;
      det_phiy_offset[i] = 0;
      xResolution[i]=0.3;
      yResolution[i]=0.3;
   }
   nUsedEvents=0;
   alignmentTrainingTrackFraction=0;
   intervallBeginEventNo.clear();
   intervallEndEventNo.clear();
   nEvents=0;
   diaTime=TDatime(1995,0,0,0,0,0);
   silTime=TDatime(1995,0,0,0,0,0);
}

void TDetectorAlignment::addEventIntervall(UInt_t first, UInt_t last){
	if(first>last)
		return;
	intervallBeginEventNo.push_back(first);
	intervallEndEventNo.push_back(last);
	nEvents+=(last-first);
}


int TDetectorAlignment::getVerbosity() const
{
    return verbosity;
}

void TDetectorAlignment::AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::addPhiXOffset of Plane "<<plane<<": "<<addPhiXOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetPhiXOffset[plane].push_back(addPhiXOffset);
		det_phix_offset[plane]+=addPhiXOffset;
	}
	if(verbosity)cout<<det_phix_offset[plane]<<endl;
}
void TDetectorAlignment::AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::addPhiYOffset of Plane "<<plane<<": "<<addPhiYOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetPhiYOffset[plane].push_back(addPhiYOffset);
		det_phiy_offset[plane]+=addPhiYOffset;
	}
	if(verbosity)cout<<det_phiy_offset[plane]<<endl;
}

void TDetectorAlignment::AddToXOffset(UInt_t plane, Float_t addXOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::AddToXOffset of Plane "<<plane<<": "<<addXOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetXOffset[plane].push_back(addXOffset);
		det_x_offset[plane]+=addXOffset;
	}
	if(verbosity)cout<<det_x_offset[plane]<<endl;
}

void TDetectorAlignment::AddToYOffset(UInt_t plane, Float_t addYOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::AddToYOffset of Plane "<<plane<<": "<<addYOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetYOffset[plane].push_back(addYOffset);
		det_y_offset[plane]+=addYOffset;
	}
	if(verbosity)cout<<det_y_offset[plane]<<endl;
}

void TDetectorAlignment::AddToZOffset(UInt_t plane, Float_t addZOffset)
{
	if(plane<6){
		vecDetZOffset[plane].push_back(addZOffset);
		det_z_offset[plane]+=addZOffset;
	}
}

std::string TDetectorAlignment::PrintXOffset(UInt_t plane,UInt_t level)
{
	stringstream output;
	output<<TCluster::Intent(level)<<"X-Offset:" <<det_x_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetXOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		output<<" "<<vecHis.at(i);
	output<<endl;
	return output.str();
}

std::string TDetectorAlignment::PrintYOffset(UInt_t plane,UInt_t level)
{
	stringstream output;
	output<<TCluster::Intent(level)<<"Y-Offset:" <<det_y_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetYOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		output<<" "<<vecHis.at(i);
	output<<endl;
	return output.str();
}

std::string TDetectorAlignment::PrintZOffset(UInt_t plane,UInt_t level)
{
	stringstream output;
	output<<TCluster::Intent(level)<<"Z-Offset:" <<det_z_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetZOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		output<<" "<<vecHis.at(i);
	output<<endl;
	return output.str();
}

std::string TDetectorAlignment::PrintPhiXOffset(UInt_t plane, UInt_t level)
{
	stringstream output;
	output<<TCluster::Intent(level)<<"PhiX-Offset:" <<det_phix_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetPhiXOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		output<<" "<<vecHis.at(i);
	output<<endl;
	return output.str();
}

std::string TDetectorAlignment::PrintPhiYOffset(UInt_t plane, UInt_t level)
{
	stringstream output;
	output<<TCluster::Intent(level)<<"PhiY-Offset:" <<det_phiy_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetPhiYOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		output<<" "<<vecHis.at(i);
	output<<endl;
	return output.str();
}

std::string  TDetectorAlignment::PrintResolution(UInt_t plane, UInt_t level){
	stringstream output;
	output<<TCluster::Intent(level)<<"Resolution of "<<plane<<" Plane in X: "<<this->getXResolution(plane)<<" and Y:"<<this->getYResolution(plane)<<endl;
	return output.str();
}

std::string TDetectorAlignment::PrintResults(UInt_t level)
{
	stringstream output;
	output<<"\n"<<TCluster::Intent(level)<<"Alignment of Run "<<runNumber<<"\n";
	output<<TCluster::Intent(level+1)<<"Silicon alignment made on "<<silTime.AsString()<<"\n";
	output<<TCluster::Intent(level+1)<<"Diamond alignment made on "<<diaTime.AsString()<<"\n";
	output<<TCluster::Intent(level+1)<<"used "<<nUsedEvents<<" from "<<nEvents<<" with a training Fraction of "<<setw(5)<<this->alignmentTrainingTrackFraction*100<<"% (";
	output<<(Float_t)nUsedEvents/alignmentTrainingTrackFraction/nEvents*100<<"%)"<<"\n";
	output<<TCluster::Intent(level+1)<<"used "<<this->nDiamondAlignmentEvents<<" from "<<nUsedEvents;
	output<<"for diamond Alignment ("<<(Float_t)nDiamondAlignmentEvents/(Float_t)nUsedEvents*100<<"%) at a max Chi2 cut at "<<setw(4)<<this->diaChi2<<"\n";
	for(UInt_t i=0;i<this->intervallBeginEventNo.size()&&i<intervallEndEventNo.size();i++){
		output<<TCluster::Intent(level+2)<<"Intervall "<<i<<" of "<<intervallBeginEventNo.size()<<": [ ";
		output<<intervallBeginEventNo.at(i)<<" , "<<intervallEndEventNo.at(i)<<" ]"<<"\n";
	}

	output<<"\nAlignment Results"<<"\n";

	for(UInt_t plane=0;plane<this->nDetectors;plane++){
		output<<TCluster::Intent(level+1)<<" Plane "<< plane<<endl;
		output<<PrintXOffset(plane,level+2);
		output<<PrintYOffset(plane,level+2);
		output<<PrintZOffset(plane,level+2);
//		cout<<endl;
		output<<PrintPhiXOffset(plane,level+2);
		output<<PrintPhiYOffset(plane,level+2);
		output<<PrintResolution(plane,level+2);
		output<<endl;
	}
	cout<<output.str()<<flush;
	return output.str();
}

Double_t TDetectorAlignment::getXResolution(UInt_t plane)
{
	if(plane<6)
    return xResolution[plane];
	return -9999;
}

void TDetectorAlignment::setXResolution(Double_t xres,UInt_t plane)
{
	printf("Set X-Resolution of Plane %d to %2.6f\n",plane,xres);
	if(plane<6)
		xResolution[plane] = xres;
}

Double_t TDetectorAlignment::getYResolution(UInt_t plane)
{
	if(plane<6)
    return yResolution[plane];
	return -9999;
}

void TDetectorAlignment::setYResolution(Double_t resolution,UInt_t plane)
{
	printf("Set Y-Resolution of Plane %d to %2.6f\n",plane,resolution);
	if(plane<6)
		yResolution[plane] = resolution;
}

void TDetectorAlignment::setVerbosity(int verbosity)
{
	if(verbosity!=this->verbosity){
		cout<<"TDetectorAlignment::setVerbosity from "<<getVerbosity()<<" to "<<verbosity<<endl;
		this->verbosity = verbosity;
	}
}


