/**
 * @file TCellAnalysisClass.hh
 *
 * @date Jun 20, 2013
 * @author bachmair
 * @description
 */

#ifndef TCELLANALYSISCLASS_HH_
#define TCELLANALYSISCLASS_HH_

#include <vector>

#include "TTree.h"
#include "TH1.h"
#include "TSettings.class.hh"
/*
 *
 */
class TCellAnalysisClass {
public:
	TCellAnalysisClass(TSettings *settings);
	virtual ~TCellAnalysisClass();
	void addEvent(Float_t xPred,Float_t yPred, Int_t cell, UInt_t quarter, Float_t relCellPosX, Float_t relCellPosY, TCluster clus);
	TH1* getHistogram(string histoName,string varexp,string selection, string drawOption);
	int getEntries(string selection=""){return cellAnalysisTree->GetEntries(selection.c_str());}
	TTree* cellAnalysisTree;
private:
	void initialiseBranchAddresses();
	/**
	 * quarterFloat/cellFloat is a definition of a variable type which contians a Float for each entry in a three/two dimensional matrix
	 * first row
	 * second column
	 * (third quarter)
	 */
	typedef std::vector< std::vector <std::vector<Float_t> > > quarterFloat_t;
	typedef std::vector <std::vector<Float_t> > cellFloat;


	quarterFloat_t vecAvrgPH;
	cellFloat vecAvrgCellPH;
	TSettings *settings;

	Float_t ph;
	Float_t xPred;
	Float_t yPred;
	Int_t nCell;
	Int_t nColumn;
	Int_t nRow;
	Int_t nQuarter;
	Float_t relCellPosX;
	Float_t relCellPosY;
	Int_t clusterSize;
	TCluster* pCluster;

};

#endif /* TCELLANALYSISCLASS_HH_ */
