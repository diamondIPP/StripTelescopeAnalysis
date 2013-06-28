//settings

#ifndef TSETTINGS_CLASS_HH
#define TSETTINGS_CLASS_HH


//C++ standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>

//ROOT libraries
#include "Rtypes.h"


#include "ChannelScreen.hh"
#include "TChannelMapping.hh"
#include "TObject.h"
#include "TCluster.hh"
#include "TPlaneProperties.hh"
#include "TDiamondPattern.hh"
//#include "TSettings.class.hh"
#include "TSystem.h"
#include "TFile.h"
#include "TFiducialCut.hh"
#include "TFidCutRegions.hh"
#include "TRunInfo.hh"
#include <string>
#include <sys/stat.h>

class TSettings:public TObject {
private:
	std::string runDescription;
	std::string outputDir;
	std::string inputDir;
public:
	static bool existsDirectory(std::string dir);
	TSettings(TRunInfo* runInfo);
	TSettings(UInt_t runNumber=0);
	TSettings(std::string fileName,UInt_t runNumber=0);
	std::string getAbsoluteOuputPath(bool withRunDescribtion=0);
	std::string getAbsoluteInputPath(){return inputDir;};//todo
	std::string getRawTreeFilePath();
	std::string getPedestalTreeFilePath();
	std::string getClusterTreeFilePath();
	std::string getAlignmentFilePath();
	std::string getSelectionTreeFilePath();
	std::string getSelectionAnalysisPath(){return this->getAbsoluteOuputPath(true).append("/selectionAnalysis/");};

	std::string getSelectionPath(){return this->getAbsoluteOuputPath(true).append("/selectionss/");}
	std::string getEtaDistributionPath(Int_t step=-1);
	std::string get3dDiamondTreeFilePath();
	std::string get3dDiamondAnalysisPath(){return this->getAbsoluteOuputPath(true).append("/3dDiamondAnalysis/");};
	bool doCommonModeNoiseCorrection() const {return DO_CMC;}
	void goToRawTreeDir();
	void goToClusterTreeDir(){goToDir(this->getAbsoluteOuputPath(false));}
	void goToSelectionTreeDir();
	void goToOutputDir();
	void goToPedestalTreeDir(){goToDir(this->getAbsoluteOuputPath(false));}
	void goToAlignmentRootDir(){goToDir(this->getAbsoluteOuputPath(false));}


	void goTo3dDiamondTreeDir();
	void goTo3dDiamondAnalysisDir(){goToDir(this->get3dDiamondAnalysisPath());}
	void goToPedestalAnalysisDir();
	void goToClusterAnalysisDir();
	void goToSelectionDir(){goToDir(this->getAbsoluteOuputPath(true).append("/selections/"));}
	void goToSelectionAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/selectionAnalysis/"));}
	void goToAlignmentDir(){goToDir(this->getAbsoluteOuputPath(true).append("/alignment/"));}
	void goToAlignmentAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/anaAlignmnet/"));}
	void goToTransparentAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/transparentAnalysis/"));}
	std::string getTransparentAnalysisDir(){return this->getAbsoluteOuputPath(true).append("/transparentAnalysis/");}
	std::string getToPedestalAnalysisDir(){return this->getAbsoluteOuputPath(true).append("/pedestalAnalysis/");}
	std::string getAlignmentDir(){return this->getAbsoluteOuputPath(true).append("/alignment/");};
	std::string getAlignmentAnalysisFilePath(){return this->getAbsoluteOuputPath(true).append("/anaAlignmnet/");};
	bool isSpecialAnalysis(){return getRunDescription().at(0)!='0';};

private:
	void goToDir(std::string dir);
	void setVerbosity(int verb){this->verbosity=verb;cout<<"Set Verbosity to: "<<verbosity<<endl;}
	void checkSettings();
public:
	virtual ~TSettings();
	void setFidCut(TFiducialCut* fidcut);
	void saveSettings();
	void loadSettingsFromRootFile();
	void compareSettings();
	void createSettingsRootFile();
	void setRunDescription(std::string runDescription);
	void setOutputDir(std::string ouputDir){this->outputDir=outputDir;}
	void setInputDir (std::string inputDir){this->inputDir=inputDir;};
	bool is3dDiamond(){return b3dDiamond;};
	bool b3dDiamond;
	std::string getInputDir()const {return inputDir;};
	std::string getOutputDir()const {return outputDir;};
	enum enumAlignmentTrainingMethod{enumFraction=0, enumEvents=1};
	std::string getRunDescription() const {return runDescription;};
	Float_t getPHinSigmaPlotFactor() const{return 0.8;}
	Float_t getClusterSeedFactor(UInt_t det,UInt_t ch);
	Float_t getClusterHitFactor(UInt_t det,UInt_t ch);
	Float_t getAlignment_chi2() const;
	void setAlignment_chi2(Float_t alignment_chi2);
	Float_t getTransparentChi2() const{return transparentChi2;}
	void setTransparentChi2(Float_t chi2){transparentChi2=chi2;}
	float getFix_dia_noise() const;
	Int_t getIter_Size() const;
	Int_t getPedestalSildingLength(){return getIter_Size();};
	Int_t getTaylor_speed_throttle() const;
	Int_t getDia_input() const;
	Float_t getDi_Pedestal_Hit_Factor() const;
	Float_t getSi_Pedestal_Hit_Factor() const;
	Int_t getDO_CMC() const;
	Int_t getCMN_cut() const;
	Float_t getDi_Cluster_Hit_Factor() const;
	Float_t getDi_Cluster_Seed_Factor() const;
	Float_t getSi_Cluster_Hit_Factor() const;
	Float_t getSi_Cluster_Seed_Factor() const;
	Float_t getSi_avg_fidcut_xhigh() const;
	Float_t getSi_avg_fidcut_xlow() const;
	Float_t getSi_avg_fidcut_yhigh() const;
	Float_t getSi_avg_fidcut_ylow() const;
	Int_t getClosePlotsOnSave() const;
	Int_t getIndexProduceSwitch() const;
	Float_t getPulse_height_di_max() const;
	Int_t getPulse_height_num_bins() const;
	Float_t getPulse_height_si_max() const;
	Float_t getPulse_height_max(UInt_t det) const;
	Int_t getSaveAllFilesSwitch() const;
	Float_t getSnr_distribution_di_max() const;
	Float_t getSnr_distribution_si_max() const;
	Float_t getEta_hiq_slice_hi() const;
	Float_t getEta_hiq_slice_low() const;
	Float_t getEta_lowq_slice_hi() const;
	Float_t getEta_lowq_slice_low() const;
	std::vector<Float_t> getAlignment_phi_offsets() const;
	std::vector<Float_t> getAlignment_x_offsets() const;
	std::vector<Float_t> getAlignment_y_offsets() const;
	std::vector<Float_t> getAlignment_z_offsets() const;
	Int_t getEtavsq_n_landau_slices() const;
	Int_t getSnr_plots_enable() const;
	Float_t getAlignment_training_track_fraction() const;
	ChannelScreen getDet_channel_screen(int i);
	bool isDet_channel_screened(UInt_t det,UInt_t ch);
	std::vector<int> getDet_channel_screen_channels(int i) const;
	std::vector<int> getDet_channel_screen_regions(int i) const;
	bool getAlternativeClustering() const;
	bool getUseAutoFidCut() const;
	Float_t getAutoFidCutPercentage() const{return 0.4;};//todo in settingsfile adden
	UInt_t getAutoFidCutEvents()const {return 40000;};
	UInt_t getNDiamonds()const{return nDiamonds;};
	std::string getFileName() const;
	bool getSingle_channel_analysis_enable();
	Int_t getSingle_channel_analysis_eventwindow();
	Int_t getCMN_corr_low();
	Int_t getCMN_corr_high();
	std::vector<int> getSingle_channel_analysis_channels();
	float getStore_threshold();
	UInt_t getRunNumber();

	void setAlignment_training_track_fraction(Float_t alignment_training_track_fraction);
	void setFix_dia_noise(float fix_dia_noise);
	void setIter_Size(Int_t Iter_Size);
	void setTaylor_speed_throttle(Int_t Taylor_speed_throttle);
	void setDi_Pedestal_Hit_Factor(Float_t Di_Pedestal_Hit_Factor);
	void setDia_input(Int_t dia_input);
	void setSi_Pedestal_Hit_Factor(Float_t Si_Pedestal_Hit_Factor);
	void setDO_CMC(Int_t DO_CMC);
	void setCMN_cut(Int_t CMN_cut);
	void setDi_Cluster_Hit_Factor(Float_t Di_Cluster_Hit_Factor);
	void setDi_Cluster_Seed_Factor(Float_t Di_Cluster_Seed_Factor);
	void setSi_Cluster_Hit_Factor(Float_t Si_Cluster_Hit_Factor);
	void setSi_Cluster_Seed_Factor(Float_t Si_Cluster_Seed_Factor);
	void setSi_avg_fidcut_xhigh(Float_t si_avg_fidcut_xhigh);
	void setSi_avg_fidcut_xlow(Float_t si_avg_fidcut_xlow);
	void setSi_avg_fidcut_yhigh(Float_t si_avg_fidcut_yhigh);
	void setSi_avg_fidcut_ylow(Float_t si_avg_fidcut_ylow);
	void setClosePlotsOnSave(Int_t ClosePlotsOnSave);
	void setIndexProduceSwitch(Int_t IndexProduceSwitch);
	void setPulse_height_di_max(Float_t pulse_height_di_max);
	void setPulse_height_num_bins(Int_t pulse_height_num_bins);
	void setPulse_height_si_max(Float_t pulse_height_si_max);
	void setSaveAllFilesSwitch(Int_t SaveAllFilesSwitch);
	void setSnr_distribution_di_max(Float_t snr_distribution_di_max);
	void setSnr_distribution_si_max(Float_t snr_distribution_si_max);
	void setEta_hiq_slice_hi(Float_t eta_hiq_slice_hi);
	void setEta_hiq_slice_low(Float_t eta_hiq_slice_low);
	void setEta_lowq_slice_hi(Float_t eta_lowq_slice_hi);
	void setEta_lowq_slice_low(Float_t eta_lowq_slice_low);
	void setAlignment_phi_offsets(std::vector<Float_t> alignment_phi_offsets);
	void setAlignment_x_offsets(std::vector<Float_t> alignment_x_offsets);
	void setAlignment_y_offsets(std::vector<Float_t> alignment_y_offsets);
	void setAlignment_z_offsets(std::vector<Float_t> alignment_z_offsets);
	void setEtavsq_n_landau_slices(Int_t etavsq_n_landau_slices);
	void setSnr_plots_enable(Int_t snr_plots_enable);
	void setDet_channel_screen(int i,ChannelScreen Det_channel_screen);
	void setDet_channel_screen_channels(int i,std::vector<int> Det_channel_screen_channels);
	void setDet_channel_screen_regions(int i,std::vector<int> Det_channel_screen_regions);
	void setAlternativeClustering(bool AlternativeClustering);
	void setUseAutoFidCut(bool UseAutoFidCut);
	void setSingle_channel_analysis_enable(bool singleChannelAnalysisEnable);
	void setSingle_channel_analysis_eventwindow(Int_t singleChannelAnalysisEventWindow);
	void setCMN_corr_low(Int_t CMN_corr_low);
	void setCMN_corr_high(Int_t CMN_corr_high);
	void setStore_threshold(float storeThreshold);
	void setNDiamonds(UInt_t nDia);
	Int_t getPlotChannelOn() const;
	void setPlotChannelOn(Int_t plotChannelOn);
	Int_t getMakeBufferPlots() const;
	Int_t getPlotDiamond() const;
	void setMakeBufferPlots(Int_t makeBufferPlots);
	void setPlotDiamond(Int_t plotDiamond);
	Int_t getEventPrintHex() const;
	Int_t getMakeDiamondPlots() const;
	Int_t getMakeHits2D() const;
	Int_t getMakeNoise2D() const;
	Int_t getMakePedRmsTree() const;
	Int_t getMakePullDist() const;
	Int_t getSingleChannel2000plots() const;
	void setEventPrintHex(Int_t eventPrintHex);
	void setMakeDiamondPlots(Int_t makeDiamondPlots);
	void setMakeHits2D(Int_t makeHits2D);
	void setMakeNoise2D(Int_t makeNoise2D);
	void setMakePedRmsTree(Int_t makePedRmsTree);
	void setMakePullDist(Int_t makePullDist);
	void setSingleChannel2000plots(Int_t singleChannel2000plots);
	UInt_t getPlottedChannel() const;
	void setPlottedChannel(UInt_t plottedChannel);
	Int_t getMaxBufferPlots() const;
	void setMaxBufferPlots(Int_t maxBufferPlots);
	Float_t getRmsSigmaDifferenceCut() const;
	void setRmsSigmaDifferenceCut(Float_t rmsSigmaDifferenceCut);
	Int_t getHighRmsCut() const;
	Float_t getRmsCut() const;
	void setHighRmsCut(Int_t highRmsCut);
	void setRmsCut(Float_t rmsCut);
	Float_t getMaxNoise2D() const;
	Int_t getSingleTrack2D() const;
	Int_t getSingleTrack2DmaxClusterSize() const;
	Int_t getZoomDiamondPlots() const;
	void setMaxNoise2D(Float_t maxNoise2D);
	void setSingleTrack2D(Int_t singleTrack2D);
	void setSingleTrack2DmaxClusterSize(Int_t singleTrack2DmaxClusterSize);
	void setZoomDiamondPlots(Int_t zoomDiamondPlots);
	Float_t getRes_keep_factor();
	UInt_t GetSiliconAlignmentSteps(){return siliconAlignmentSteps;}
	UInt_t GetDiamondAlignmentSteps(){return diamondAlignmentSteps;}
	bool doAllAlignmentPlots(){return bDoAllAlignmentPlots;}
	enumAlignmentTrainingMethod getTrainingMethod() const;
	void setTrainingMethod(enumAlignmentTrainingMethod trainingMethod);
	void Print();
	UInt_t getDetChannelNo(UInt_t vaCh);
	UInt_t getVaChannelNo(UInt_t detChNo);
	Int_t getVerbosity();
	bool useForAlignment(UInt_t eventNumber, UInt_t nEvents=0);
	bool isInAlignmentFiducialRegion(Float_t, Float_t);
	UInt_t getAlignmentTrainingTrackNumber() const {return alignment_training_track_number;}
	Float_t getAlignmentPrecisionOffset()const{return alignmentPrecision_Offset;}
	Float_t getAlignmentPrecisionAngle()const{return alignmentPrecision_Angle;}
	bool resetAlignment() const{return bResetAlignment;};
	//	void setAlignmentTrainingTrackNumber(UInt_t alignmentTrainingTrackNumber);
	Int_t getNDiaDetectorAreas(){return vecDiaDetectorAreasInChannel.size();}
	TFidCutRegions* getSelectionFidCuts(){return fidCutsSelection;}
	TFidCutRegions* get3dFidCuts(){return fidCuts3D;};
	Float_t getMinDiamondChannel();
	Float_t getMaxDiamondChannel();
	std::pair< Int_t , Int_t > getDiaDetectorArea(Int_t n);
	bool isInDiaDetectorArea(Int_t ch,Int_t area);
	bool isClusterInDiaDetectorArea(TCluster cluster, Int_t area);
	int getDiaDetectorAreaOfChannel(Int_t ch);
	bool isDiaDetectorAreaBorderChannel(UInt_t ch);
	bool isMaskedCluster(UInt_t det, TCluster cluster,bool checkAdjacentChannels=true);
	bool hasBorderSeed(UInt_t det, TCluster cluster);
	bool hasBorderHit(UInt_t det, TCluster cluster);
	Float_t getSiliconPitchWidth(){return this->pitchWidthSil;}
	Float_t getDiamondPitchWidth(){return this->pitchWidthDia;}
	Float_t convertChannelToMetric(UInt_t det, Float_t channel);
	Float_t convertMetricToChannelSpace(UInt_t det, Float_t metricValue);
	void PrintPatterns(int k=0);
	Float_t getChi2Cut3D(){return chi2Cut3D;}
	Float_t get3DYOffset(){return yOffset3D;};//todo
	vector<int> getxEdgeFicucialRegion(){return xEdgeFicucialRegion;};
	vector<int> getyEdgeFicucialRegion(){return yEdgeFicucialRegion;};
	vector<int> getstripAnalysisFidCut(){return stripAnalysisFidCut;};
	vector<int> get3DnHAnalysisFidCut(){return TDnHAnalysisFidCut;};
	vector<int> get3DwHAnalysisFidCut(){return TDwHAnalysisFidCut;};
private:
	TFidCutRegions* fidCutsSelection;
	TFidCutRegions* fidCuts3D;
protected:
	float store_threshold;
private:
	bool isStandardSelectionFidCut,isStandard3dFidCut;
	void checkAlignmentFidcuts();
	void SetFileName(std::string fileName);
	void LoadSettings();
	void DefaultLoadDefaultSettings();
	void ParseStringArray(std::string key, std::string value, std::vector<std::string> &vec);
	void ParseFloatArray(std::string key, std::string value, std::vector<float> & vec);
	void ParseIntArray(std::string key, std::string value, std::vector<int> & vec);
	void ParseRegionArray(std::string key, std::string value, std::vector< std::pair<Int_t, Int_t> > &vec);
	void ParsePattern(std::string key, std::string value);
	void ParseFidCut(std::string key, std::string value, TFidCutRegions* fidCutRegions,bool &isStandardFidCut);
	void ParseScreenedChannelArray(std::string key, std::string value, std::vector<int> & vec);
	std::pair< std::string,std::string > ParseRegionString(std::string key, string value);
	bool ParseFloat(std::string key, std::string value,float  &output);
	Float_t ParseFloat(std::string key, std::string value){float output;ParseFloat(key,value,output);return output;}
	Int_t ParseInt(std::string key, std::string value){cout<<"Parse: "<<value<<endl;Int_t output;ParseInt(key,value,output);return output;}
	Int_t ParseInt(std::string value){return  (int)strtod(value.c_str(),0);}
	bool ParseInt(std::string key, std::string value, int &output);
	bool ParseInt(std::string key, std::string value, UInt_t &output);
	bool ParseBool(std::string key, std::string value, bool &output);
	void ParseCellArray(std::string key, std::string value, vector<Int_t> &cells);
	void Parse(std::string key, std::string value, std::vector<float> & vec){ ParseFloatArray(key,value,vec);}
	void Parse(std::string key, std::string value, std::vector<int> & vec){ ParseIntArray(key,value,vec);}
	bool Parse(std::string key, std::string value, bool &output){return ParseBool(key,value,output);}
	bool Parse(std::string key, std::string value, int &output){return ParseInt(key,value,output);}
	bool Parse(std::string key, std::string value, UInt_t &output){return ParseInt(key,value,output);}
	bool Parse(std::string key, std::string value, float &output){return ParseFloat(key,value,output);}
	pair<char,int> ParseCellPosition(std::string value);

private:
	std::string path;
	std::string fileName;
	TSystem *sys;
	TFile *settingsFile;
private:
	Float_t chi2Cut3D;
	Float_t transparentChi2;
	std::vector< std::pair<Int_t,Int_t> > vecDiaDetectorAreasInChannel;
	bool bResetAlignment;
	Float_t alignmentPrecision_Offset;
	Float_t alignmentPrecision_Angle;
	bool bDoAllAlignmentPlots;
	UInt_t siliconAlignmentSteps;
	UInt_t diamondAlignmentSteps;
	Int_t nDiamonds;
	Int_t SaveAllFilesSwitch;
	Int_t ClosePlotsOnSave;
	Int_t IndexProduceSwitch;
	float fix_dia_noise;
	Int_t Iter_Size;
	Int_t Taylor_speed_throttle;
	Int_t dia_input;
	Float_t Si_Pedestal_Hit_Factor;
	Float_t Di_Pedestal_Hit_Factor;
	bool single_channel_analysis_enable;
	Int_t single_channel_analysis_eventwindow;
	std::vector<int> single_channel_analysis_channels;
	Int_t DO_CMC;
	Int_t CMN_cut;
	Int_t CMN_corr_low;
	Int_t CMN_corr_high;
	Float_t Si_Cluster_Seed_Factor;
	Float_t Si_Cluster_Hit_Factor;
	Float_t Di_Cluster_Seed_Factor;
	Float_t Di_Cluster_Hit_Factor;
	Float_t si_avg_fidcut_xlow;
	Float_t si_avg_fidcut_xhigh;
	Float_t si_avg_fidcut_ylow;
	Float_t si_avg_fidcut_yhigh;
	Int_t pulse_height_num_bins;
	Float_t pulse_height_si_max;
	Float_t pulse_height_di_max;
	Float_t snr_distribution_si_max;
	Float_t snr_distribution_di_max;
	Float_t eta_lowq_slice_low;
	Float_t eta_lowq_slice_hi;
	Float_t eta_hiq_slice_low;
	Float_t eta_hiq_slice_hi;
	Int_t etavsq_n_landau_slices;
	Int_t snr_plots_enable;

	std::vector<Float_t> alignment_x_offsets;
	std::vector<Float_t> alignment_y_offsets;
	std::vector<Float_t> alignment_phi_offsets;
	std::vector<Float_t> alignment_z_offsets;
	std::vector<Int_t> alignmentFidCuts;

	Float_t alignment_training_track_fraction;
	UInt_t alignment_training_track_number;
	enumAlignmentTrainingMethod trainingMethod;
	std::vector<int> Det_channel_screen_channels[9];
	std::vector<int> Det_channel_screen_regions[9];
	ChannelScreen Det_channel_screen[9];
	bool dia_x_aligned;
	bool eta_correction;
	Float_t alignment_chi2;
	bool UseAutoFidCut;
	bool AlternativeClustering;
	Int_t plotChannel_on;
	Float_t res_keep_factor;
	Int_t plotDiamond; //make Buffer Noise plots for the diamond instead
	Int_t makeBufferPlots; //make Buffer Plot whenever sigma and rms differ by rms_sigma_difference_cut
	Int_t SingleChannel2000plots; //make SC_Pedestal plots for all silicon detectors and channels
	Int_t makeDiamondPlots; //make DC_Pedestal plots for all diamond channels
	Int_t makeHits2D; //make 2D histogram of hits and seeds
	Int_t makeNoise2D; //make 2D histogram of noise per channel
	Int_t makePullDist; //make pull distribution
	Int_t makePedRMSTree; //make .root file of pedestal and rms values
	Int_t eventPrintHex; //print hex (should match .rz data)
	UInt_t plottedChannel;
	Int_t maxBufferPlots;
	Float_t rms_sigma_difference_cut;
	Int_t high_rms_cut; //cut on absolute rms value instead of comparing to Gaussian
	Float_t rms_cut; //value to use if high_rms_cut

	Int_t zoomDiamondPlots; //zoom in on DC_Pedestal (100 event / window)

	Int_t singleTrack2D; //plot single tracks only in 2D hits histogram
	Int_t singleTrack2DmaxClusterSize; //max size of clusters in silicon track (cluster = Di_Hit_Factor hits; no check for seeds/shoulders)

	Float_t maxNoise2D; //highest noise value plotted in 2D noise histogram
	UInt_t runNumber;
	std::vector<Float_t>clusterHitFactors;
	std::vector<Float_t>clusterSeedFactors;

	vector<Float_t> vecClusterSeedFactorsDia;
	vector<Float_t> vecClusterHitFactorsDia;
	TChannelMapping *diamondMapping;
	Float_t pitchWidthSil;
	Float_t pitchWidthDia;
private:
	//Filter tracks not in good fiducial region w/o bad strips
	Int_t align_sil_fid_xlow;
	Int_t align_sil_fid_xhi;
	Int_t align_sil_fid_ylow;
	Int_t align_sil_fid_yhi;
public:
	TDiamondPattern diamondPattern;
private:
	int verbosity;
	Float_t silPitchWidth;
	Float_t diaPitchWidth;
	Float_t diaOffsetMetricSpace;
	Float_t diaStartingChannel;
private:
	Float_t yOffset3D;
	vector<int> xEdgeFicucialRegion;
	vector<int> yEdgeFicucialRegion;
	vector<int> stripAnalysisFidCut;
	vector<int> TDnHAnalysisFidCut;
	vector<int> TDwHAnalysisFidCut;

	int nRows3d;
	int nColumns3d;
	vector<Int_t> badCells3d;
	vector<Int_t> goodCells3d;
	vector<Int_t> deadCell3d;
	int XmetalisationStart3d;
	int XmetalisationEnd3d;
	int YmetalisationEnd3d;
public:
	void setNRows3d(int nRows){nRows3d=nRows;};
	int getNRows3d(){return nRows3d;};
	int getNColumns3d(){return nColumns3d;};
	int getNQuarters3d(){return 4;}
	void setNColumns3d(int nColumns){nColumns3d=nColumns;};
	void setXMetalisationStart3d(int nXMetalisationStart3d){XmetalisationStart3d=nXMetalisationStart3d;};
	int getXMetalisationStart3d(){return XmetalisationStart3d;};
	void setXMetalisationEnd3d(int nXMetalisationEnd3d){XmetalisationEnd3d=nXMetalisationEnd3d;};
	int getXMetalisationEnd3d(){return XmetalisationEnd3d;};
	void setYMetalisationEnd3d(int nYMetalisationEnd3d){YmetalisationEnd3d=nYMetalisationEnd3d;};
	int getYMetalisationEnd3d(){return YmetalisationEnd3d;};
	//void setNColumns3d(int nColumns){nColumns3d=nColumns;};
	vector<Int_t> getGoodCells3D(){return goodCells3d;};
	vector<Int_t> getBadCells3D(){return badCells3d;};
	vector<Int_t> getDeadCell3D(){return deadCell3d;};
	int get3DCellNo(char row, int column);
	int get3DCellNo(pair<char,int> pos){return get3DCellNo(pos.first,pos.second);};
	UInt_t get3dWithHolesDiamondPattern(){return 3;};

	ClassDef(TSettings,4);
};
#endif
