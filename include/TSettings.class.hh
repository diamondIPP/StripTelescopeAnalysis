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
#include "TCutG.h"
#include "TObject.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMacro.h"
#include "TProfile2D.h"
#include "ChannelScreen.hh"
#include "TChannelMapping.hh"
#include "TCluster.hh"
#include "TPlaneProperties.hh"
#include "TDiamondPattern.hh"
//#include "TSettings.class.hh"

#include "TFiducialCut.hh"
#include "TFidCutRegions.hh"
#include "TRunInfo.hh"
#include <string>
#include <sys/stat.h>
class TSettings:public TObject {
public:
	enum alignmentMode {normalMode=0,transparentMode=1};
	enum enumRunDescription {leftDia=1,rightDia=2,allDia=0,unknown=-1};
private:
	std::string runDescription;
	std::string outputDir;
	std::string inputDir;
public:
	TMacro* GetMacro(){return macro;}
	static UInt_t GetColor(int no);
	static bool existsDirectory(std::string dir);
	TSettings(TRunInfo* runInfo);
	TSettings(UInt_t runNumber=0);
	TSettings(std::string fileName,UInt_t runNumber=0);
	virtual ~TSettings();
    TString GetCurrentBegin(){return currentBegin;}
    TString GetCurrentEnd(){return currentEnd;}
	std::string getAbsoluteOuputPath(bool withRunDescribtion=0);
	std::string getAbsoluteInputPath(){return inputDir;};//todo
	std::string getRawTreeFilePath();
	std::string getPedestalTreeFilePath();
	std::string getClusterTreeFilePath();
	std::string getAlignmentFilePath(TSettings::alignmentMode mode = TSettings::normalMode);
	std::string getSelectionTreeFilePath();
	std::string getSelectionAnalysisPath(){return this->getAbsoluteOuputPath(true).append("/selectionAnalysis/");};
	std::string getCrossTalkFactorsFileName();
	std::string getSelectionPath(){return this->getAbsoluteOuputPath(true).append("/selectionss/");}
	std::string getEtaDistributionPath(Int_t step=-1);
	std::string get3dDiamondTreeFilePath();
	std::string get3dDiamondAnalysisPath();
	std::string getResultsRootFilePath();
	bool doCommonModeNoiseCorrection() const {return DO_CMC;}
	void goToRawTreeDir();
	void goToClusterTreeDir(){goToDir(this->getAbsoluteOuputPath(false));}
	void goToSelectionTreeDir();
	void goToOutputDir();
	void goToPedestalTreeDir(){goToDir(this->getAbsoluteOuputPath(false));}
	void goToAlignmentRootDir(alignmentMode mode = normalMode){goToDir(this->getAbsoluteOuputPath(false));}
	Float_t getNegativeChargeCut() const{return negativeChargeCut;}
	Float_t getNegativeChargeCutStrip() const{return negativeChargeCutStrip;}
	void setNegativeChargeCut(Float_t value){negativeChargeCut  = value;}


	void goTo3dDiamondTreeDir();
	void goTo3dDiamondAnalysisDir(){goToDir(this->get3dDiamondAnalysisPath());}
	void goToPedestalAnalysisDir();
	void goToClusterAnalysisDir();
	void goToSelectionDir(){goToDir(this->getAbsoluteOuputPath(true).append("/selections/"));}
	void goToSelectionAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/selectionAnalysis/"));}
	void goToAlignmentDir(alignmentMode mode =normalMode);
	void goToAlignmentAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/anaAlignmnet/"));}
	void goToTransparentAnalysisDir(){goToDir(this->getAbsoluteOuputPath(true).append("/transparentAnalysis/"));}
	std::string getTransparentAnalysisDir(TSettings::alignmentMode mode = TSettings::normalMode);

	std::string getClusterAnalysisDir(){return this->getAbsoluteOuputPath(true).append("/clustering/");}
	std::string getToPedestalAnalysisDir(){return this->getAbsoluteOuputPath(true).append("/pedestalAnalysis/");}
	std::string getAlignmentDir(TSettings::alignmentMode mode = normalMode);
	std::string getAlignmentAnalysisFilePath(){return this->getAbsoluteOuputPath(true).append("/anaAlignmnet/");};
	bool isSpecialAnalysis(){return getRunDescription().at(0)!='0';};

	enumRunDescription getAnalysedDiamond();
	TCutG* GetCell(Int_t nCell,TString name);
	void DrawMetallisationGrid(TCanvas* nCanvas,  int DiamondPattern);
	vector< TH1*> sortHistosByPulseHeight(vector<TH1*> vec);
	static bool SorterForPulseHeightOfHisto(TH1* a, TH1* b);
private:
	void goToDir(std::string dir);
	void setVerbosity(int verb){this->verbosity=verb;cout<<"Set Verbosity to: "<<verbosity<<endl;}
	void checkSettings();
public:
	bool IsValidChannel(UInt_t det, UInt_t channel){return TPlaneProperties::IsValidChannel(det,channel);};
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
	bool doTransparentAlignment()const {return bTransparentAlignment;}
	int getAlignmentMode()const {return detectorsToAlign;}
	int getDetectorsToAlign()const{return detectorsToAlign;}
	std::string getInputDir()const {return inputDir;};
	std::string getOutputDir()const {return outputDir;};
	enum enumAlignmentTrainingMethod{enumFraction=0, enumEvents=1};
	std::string getRunDescription() const {return runDescription;};
	Int_t getRepeaterCard(){return repeaterCardNo;}
	Float_t getPHinSigmaPlotFactor() const{return 0.8;}
	Float_t getClusterSeedFactor(UInt_t det,UInt_t ch);
	Float_t getClusterHitFactor(UInt_t det,UInt_t ch);
	Float_t getAlignment_chi2() const;
	bool IgnoreStripForAlignment(UInt_t det,Float_t predHitPosDetCh);
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
	Int_t   getNoise_si_num_bins() const;
	Float_t getNoise_si_max     () const;
	Int_t   getNoise_di_num_bins() const;
	Float_t getNoise_di_max     () const;
	Int_t getSaveAllFilesSwitch() const;
	Int_t getSnr_distribution_num_bins() const;
	Float_t getSnr_distribution_di_max() const;
	Float_t getSnr_distribution_si_max() const;
	Float_t getEta_hiq_slice_hi() const;
	Float_t getEta_hiq_slice_low() const;
	Float_t getEta_lowq_slice_hi() const;
	Float_t getEta_lowq_slice_low() const;
	UInt_t getMaxAllowedClustersize(UInt_t det);

	void setAsymmetricSample(bool asymmetricSample) {bAsymmetricSample = asymmetricSample;}

	std::vector<Float_t> getAlignment_phi_offsets() const;
	std::vector<Float_t> getAlignment_x_offsets() const;
	std::vector<Float_t> getAlignment_y_offsets() const;
	std::vector<Float_t> getAlignment_z_offsets() const;
	Int_t getEtavsq_n_landau_slices() const;
	Int_t getSnr_plots_enable() const;
	Float_t getAlignment_training_track_fraction() const;
	ChannelScreen getDet_channel_screen(int i);
	bool isDet_channel_screened(UInt_t det,UInt_t ch);
	bool IsMasked(UInt_t det, UInt_t ch){return isDet_channel_screened(det,ch);}
	bool IsMasked(UInt_t det, Float_t ch){return (IsMasked(det,(UInt_t)ch)||IsMasked(det,(UInt_t)ch+1));}
	UInt_t getNoisePlotChannel(){return 48;}//todo: variable in setttings file
	std::vector<int> getDet_channel_screen_channels(int i) const;
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

	bool isAsymmetricSample(){return bAsymmetricSample;};

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
	void setNoise_si_num_bins(Int_t   noise_si_num_bins);
	void setNoise_si_max     (Float_t noise_si_max     );
	void setNoise_di_num_bins(Int_t   noise_di_num_bins);
	void setNoise_di_max     (Float_t noise_di_max     );
	void setSaveAllFilesSwitch(Int_t SaveAllFilesSwitch);
	void setSnr_distribution_num_bins(Float_t snr_distribution_num_bins);
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
//	bool isInTransparentFiducialRegion(Float_t, Float_t);
	bool isInRoughFiducialCut(Float_t,Float_t);
	UInt_t getAlignmentTrainingTrackNumber() const {return alignment_training_track_number;}
	Float_t getAlignmentPrecisionOffset()const{return alignmentPrecision_Offset;}
	Float_t getAlignmentPrecisionAngle()const{return alignmentPrecision_Angle;}
	UInt_t getAlignmentEvents(UInt_t nEvents);
	bool resetAlignment() const{return bResetAlignment;};
	//	void setAlignmentTrainingTrackNumber(UInt_t alignmentTrainingTrackNumber);
	UInt_t getNDiaDetectorAreas(){return this->diamondPattern.getNPatterns();}
	TFidCutRegions* getSelectionFidCuts(){return fidCutsSelection;}
	TFidCutRegions* get3dEdgeFidCuts(){return fidCuts3DEdge;};
	TFidCutRegions* get3dMetallisationFidCuts(){return fidCuts3DMetallisation;};
	Float_t getMinDiamondChannel();
	Float_t getMaxDiamondChannel();
	std::pair< Int_t , Int_t > getDiaDetectorArea(Int_t n);
	bool isInDiaDetectorArea(Int_t ch,Int_t area);
	bool isClusterInDiaDetectorArea(TCluster cluster, Int_t area);
	bool isClusterInDiaDetectorArea(TCluster* cluster, Int_t area);
	bool isValidCellNo(UInt_t cellNo){return (cellNo<nRows3d*nColumns3d);}
	int getDiaDetectorAreaOfChannel(Int_t ch, UInt_t verbosity = 0);
	bool isDiaDetectorAreaBorderChannel(UInt_t ch);
	bool isMaskedCluster(UInt_t det, TCluster cluster,bool checkAdjacentChannels=true);
	bool checkAdjacentChannelsMasked(){return false;}//todo
	bool hasBorderSeed(UInt_t det, TCluster cluster);
	bool hasBorderHit(UInt_t det, TCluster cluster);
	Float_t getSiliconPitchWidth(){return this->pitchWidthSil;}
	Float_t getDiamondPitchWidth(){return this->pitchWidthDia;}
	Float_t getPitchWidth(UInt_t det, UInt_t area =0);
	Float_t convertChannelToMetric(UInt_t det, Float_t channel);
	Float_t convertMetricToChannelSpace(UInt_t det, Float_t metricValue);
	void PrintPatterns(int k=0);
	Float_t getChi2Cut3D(){return chi2Cut3D;}
	Float_t getChi2Cut3D_X(){return chi2Cut3D;}
	Float_t getChi2Cut3D_Y(){return chi2Cut3D;}
	Float_t getMinimalAbsoluteEtaValue(){return minAbsEtaVal;}
	void setMinimalAbsoluteEtaValue(Float_t minAbsEtaVal) {this->minAbsEtaVal = minAbsEtaVal;}
	bool isUseUserResolutionInput() const {return bUseUserResolutionInput;}
	void setUseUserResolutionInput(bool useUserResolutionInput) {bUseUserResolutionInput = useUserResolutionInput;}
	Float_t GetDefaultResolutionX(UInt_t plane){if(plane<9)return alignment_resolutions.at(plane*2);return 20;}
	Float_t GetDefaultResolutionY(UInt_t plane){if(plane<8)return alignment_resolutions.at(plane*2+1);return 20;};
	//	Float_t GetDefaultResolution(TPlaneProperties::enumCoordinate cor, UInt_t plane);//todo
	//	int getAreaOfInterest(){return inde
	Float_t get3DYOffset(){return yOffset3D;};//todo
	vector<int> getstripAnalysisFidCut(){return stripAnalysisFidCut;};//todo: ????
	vector<int> get3DnHAnalysisFidCut(){return TDnHAnalysisFidCut;};//todo: ????
	vector<int> get3DwHAnalysisFidCut(){return TDwHAnalysisFidCut;};//todo: ????
	TString  getAlignmentFidCuts(){TString output; for (UInt_t i=0;i<alignmentFidCuts.size();i++) output.Append(TString::Format("%d, ",i));return output;}
	TString getDiamond(){return diamondName;}
	Int_t getVoltage(){return voltage;}
    const std::vector<int>& getDiaChannelNoisy() const {return Dia_channel_noisy;}
    const std::vector<int>& getDiaChannelNotConnected() const {return Dia_channel_not_connected;}
    bool IsNotConnectedChannel(Int_t ch);
    bool IsNoisyChannel(Int_t ch);
    UInt_t getNEvents(){return nEvents;};
    const std::pair<float,float> get3DOverlayRange(){return OverlayRange3d;}
    bool RerunSelection(){return bRerunSelection;}
private:

	TString diamondName;
	Int_t voltage;
	TFidCutRegions* fidCutsSelection;
	TFidCutRegions* fidCuts3DEdge;
	TFidCutRegions* fidCuts3DMetallisation;
	std::pair<float,float> OverlayRange3d;
	TMacro *macro;
protected:
	float store_threshold;
private:
	Float_t minAbsEtaVal;
	bool isStandardSelectionFidCut,isStandardArea,isStandard3dEdgeFidCut,isStandard3dMetallisationFidCut;;
	void CheckAlignmentFidcuts();
	void SetFileName(std::string fileName);
	void LoadSettings();
	void DefaultLoadDefaultSettings();
	void ParseStringArray(std::string key, std::string value, std::vector<std::string> &vec);
	void ParseFloatArray(std::string key, std::string value, std::vector<float> & vec);
	void ParseFloatPair(std::string key, std::string value, std::pair<float,float> & p);
	void ParseIntArray(std::string key, std::string value, std::vector<int> & vec);
	void ParseRegionArray(std::string key, std::string value, std::vector< std::pair<Int_t, Int_t> > &vec);
	void ParsePattern(std::string key, std::string value);
    void ParseFidCut(std::string key, std::string value, TFiducialCut* fidcut);
	void ParseFidCutRegion(std::string key, std::string value, TFidCutRegions* fidCutRegions,bool &isStandardFidCut);
	void ParseScreenedChannelArray(std::string key, std::string value, std::vector<int> & vec);
	std::pair< std::string,std::string > ParseRegionString(std::string key, string value);
	bool ParseFloat(std::string key, std::string value,float  &output);
	Float_t ParseFloat(std::string key, std::string value){float output;ParseFloat(key,value,output);return output;}
	Int_t ParseInt(std::string key, std::string value){Int_t output;ParseInt(key,value,output);return output;}
	Int_t ParseInt(std::string value){return  (int)strtod(value.c_str(),0);}
	bool ParseInt(std::string key, std::string value, int &output);
	bool ParseInt(std::string key, std::string value, UInt_t &output);
	bool ParseBool(std::string key, std::string value, bool &output);
	bool ParseTString(std::string key, std::string value, TString &output);
    bool ParseString(std::string key, std::string value, string &output);
	void ParseCellArray(std::string key, std::string value, vector<Int_t> &cells);
	void ParseCellArray(string key, string value, vector<int> &vecCells, vector< vector<int> > &PtrvecCells);
	void Parse(std::string key, std::string value, std::vector<float> & vec){ ParseFloatArray(key,value,vec);}
	void Parse(std::string key, std::string value, std::vector<int> & vec){ ParseIntArray(key,value,vec);}
	bool Parse(std::string key, std::string value, bool &output){return ParseBool(key,value,output);}
	bool Parse(std::string key, std::string value, int &output){return ParseInt(key,value,output);}
	bool Parse(std::string key, std::string value, UInt_t &output){return ParseInt(key,value,output);}
	bool Parse(std::string key, std::string value, float &output){return ParseFloat(key,value,output);}
	bool Parse(std::string key, std::string value, TString &output){return ParseTString(key,value,output);}
	bool Parse(std::string key, std::string value, string &output){return ParseString(key,value,output);}
	pair<char,int> ParseCellPosition(std::string value);
	void LoadDefaultResolutions();


private:
	std::string path;
	std::string fileName;
	TSystem *sys;
	TFile *settingsFile;
private:
	bool bRerunSelection;
	bool bUseUserResolutionInput;
	bool bTransparentAlignment;
	int detectorsToAlign;
	bool bAsymmetricSample;
	Float_t chi2Cut3D;
	Float_t chi2Cut3D_X;
	Float_t chi2Cut3D_Y;
	Float_t transparentChi2;
	std::vector< std::pair<Int_t,Int_t> > vecDiaDetectorAreasInChannel;
	std::vector< std::vector<Int_t> > vecAlignmentIgnoreChannels;
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
	Int_t snr_distribution_num_bins;
	Float_t snr_distribution_si_max;
	Float_t snr_distribution_di_max;
	Int_t   noise_si_num_bins;
	Float_t noise_si_max;
	Int_t   noise_di_num_bins;
	Float_t noise_di_max;
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
	std::vector<Float_t> alignment_resolutions;
	std::vector<Int_t> alignmentFidCuts;

	Float_t alignment_training_track_fraction;
	UInt_t alignment_training_track_number;
	enumAlignmentTrainingMethod trainingMethod;
	std::vector<int> Det_channel_screen_channels[9];
	std::vector<int> Dia_channel_not_connected;
	std::vector<int> Dia_channel_noisy;
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
    TString currentBegin;
    TString currentEnd;
    Int_t repeaterCardNo;
    bool bPaperMode;
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
	Int_t nEvents;
	Float_t yOffset3D;
//	vector<int> xEdgeFicucialRegion;
//	vector<int> yEdgeFicucialRegion;
	vector<int> stripAnalysisFidCut;
	vector<int> TDnHAnalysisFidCut;
	vector<int> TDwHAnalysisFidCut;

	int b3dShortAnalysis;
	int b3dLongAnalysis;
	int b3dTransparentAnalysis;
	UInt_t nRows3d;
	UInt_t nColumns3d;
	Float_t columnRadius;//in mum
	Int_t MaxOverlayClusterSize;
	Float_t CentralColumnOverlayXLow;
	Float_t CentralColumnOverlayXHigh;
	Float_t CentralColumnOverlayYLow;
	Float_t CentralColumnOverlayYHigh;
	Float_t CentralColumnOverlayXBins;
	Float_t CentralColumnOverlayYBins;
	Float_t BiasColumnOverlayXLow;
	Float_t BiasColumnOverlayXHigh;
	Float_t BiasColumnOverlayYLow;
	Float_t BiasColumnOverlayYHigh;
	Float_t OverlayOffsetX;
	Float_t OverlayOffsetY;
	Float_t OverlayColumnPulseHeightCut;
	vector<Int_t> badCells3d;
	vector<Int_t> badCells3dnH;
	vector<Int_t> goodCells3d;
	vector< vector<Int_t> > goodCellRegions3d;
	vector<Int_t> deadCell3d;
	vector<TPlaneProperties::enumCoordinate> vecEdgePositionType;
	vector<TString> vecEdgePositionName;
	vector<Int_t> vecEdgePositionDetector;
	Float_t cellHeight;
	void CheckEdgeFidcuialCuts();
	Float_t minimumEdgeDistance;
	Float_t adcToElectronConversion;
	TString PathExtension3d;
	Float_t resolutionSNR;
	Float_t lowResponseThreshold;
	Float_t negativeChargeCut;
	Float_t negativeChargeCutStrip;
	pair<Float_t,Float_t> responseWindow;
public:
	Float_t getLowResponseThreshold(){return lowResponseThreshold;}
	TH2F* GetOverlayHisto(TString name,Int_t pattern=2,UInt_t nbinsx=15,UInt_t nbinsy=15);
	TProfile2D* GetOverlayProfile(TString name,Int_t pattern=2,UInt_t nbinsx=15,UInt_t nbinsy=15);
	std::pair<Float_t,Float_t> getResponseWindow(){return responseWindow;}
    TFiducialCut* centralRegion3DnH;
	Float_t GetCellHeight(){return cellHeight;}
	Float_t GetResolutionSNR(){return resolutionSNR;}
	Float_t GetCellWidth(int det,int pattern){return getPitchWidth(det,pattern);}
	TString getEdgePositionName(UInt_t i){if(i<vecEdgePositionName.size()) return vecEdgePositionName[i];return "";}//todo make it safe
	TPlaneProperties::enumCoordinate getEdgePositionType(UInt_t i){if(i<vecEdgePositionType.size()) return vecEdgePositionType[i];return TPlaneProperties::UNKOWN_COR;}
	TCutG* getEdgePosition(UInt_t i);
	int do3dShortAnalysis() {return b3dShortAnalysis;}
	int do3dLongAnalysis() {return b3dLongAnalysis;}
	int do3dTransparentAnalysis() {return b3dTransparentAnalysis;}
	void set3dTransparentAnalysis(int trans){b3dTransparentAnalysis = trans;}
    void PrintCellPosition(UInt_t cell, int diamondPattern);
    std::pair<Float_t, Float_t> getCellPositionX(UInt_t cell, int diamondPattern);
    std::pair<Float_t, Float_t> getCellPositionY(UInt_t cell, int diamondPattern);
    std::pair<Float_t, Float_t> getCellPositionX(UInt_t column,UInt_t row, int diamondPattern);
    std::pair<Float_t, Float_t> getCellPositionY(UInt_t column,UInt_t row,int diamondPattern);
	void setNRows3d(UInt_t nRows){nRows3d=nRows;};
	UInt_t getNRows3d(){return nRows3d;}; // Horizontally
	void setNColumns3d(UInt_t nColumns){nColumns3d=nColumns;};
	UInt_t getNColumns3d(){return nColumns3d;};//Vertical
	UInt_t GetNCells3d(){return nRows3d * nColumns3d;};
	UInt_t getNQuarters3d(){return 4;}
	//void setNColumns3d(int nColumns){nColumns3d=nColumns;};
	//vector<Int_t> getGoodCells3D(){return goodCells3d;};		//Changed by Iain 5/9/13
	vector< vector<Int_t> > getGoodCellRegions3d(){return goodCellRegions3d;}
	vector<Int_t> getBadCells3D(){return badCells3d;};
	vector<Int_t> getBadCells3DnH(){return badCells3dnH;};
	vector<Int_t> getDeadCell3D(){return deadCell3d;}
	pair<Float_t, Float_t> getAllGoodCellsXpos();
	pair<Float_t, Float_t> getAllGoodCellsYpos();
	bool isDeadCell(UInt_t nDiamondPattern, Int_t cellNo);
	bool isBadCell(UInt_t nDiamondPattern, Int_t cellNo);
	bool isBadCell(UInt_t nDiamondPattern, Float_t xDet, Float_t yDet);
	bool IsGoodCell(UInt_t nDiamondPattern, Int_t cellNo);
	int get3DQuarterNo( int column, int row,int quarter);
	int get3DCellNo(char column, int row);
	int get3DCellNo(int column, int row);
	int get3DCellNo(pair<char,int> pos){return get3DCellNo((char)pos.first,pos.second);};
	Int_t getCellNo(Float_t xDet,Float_t yDet);
	pair<int,int> getCellAndQuarterNo(Float_t xDet, Float_t yDet);

	Int_t getMaxOverlayClusterSize(){return MaxOverlayClusterSize;}
	Float_t getCentralColumnOverlayXLow(){return CentralColumnOverlayXLow;}
	Float_t getCentralColumnOverlayXHigh(){return CentralColumnOverlayXHigh;}
	Float_t getCentralColumnOverlayYLow(){return CentralColumnOverlayYLow;}
	Float_t getCentralColumnOverlayYHigh(){return CentralColumnOverlayYHigh;}
	Float_t getCentralColumnOverlayXBins(){return CentralColumnOverlayXBins;}
	Float_t getCentralColumnOverlayYBins(){return CentralColumnOverlayYBins;}

	Float_t getBiasColumnOverlayXLow(){return BiasColumnOverlayXLow;}
		Float_t getBiasColumnOverlayXHigh(){return BiasColumnOverlayXHigh;}
		Float_t getBiasColumnOverlayYLow(){return BiasColumnOverlayYLow;}
		Float_t getBiasColumnOverlayYHigh(){return BiasColumnOverlayYHigh;}

	Float_t getOverlayOffsetX(){return OverlayOffsetX;}
	Float_t getOverlayOffsetY(){return OverlayOffsetY;}
	Float_t getOverlayColumnPulseHeightCut(){return OverlayColumnPulseHeightCut;}

	inline Int_t getRowOfCell(Int_t cellNo){return cellNo % getNRows3d();}
	inline Int_t getColumnOfCell(Int_t cellNo){return cellNo / getNRows3d();}
	char getColumnChar(Int_t column);
	char getColumnCharOfCell(Int_t cellNo);
	pair<Float_t,Float_t> getRelativePositionInCell(Float_t xPredDet,Float_t yPredDet);
	UInt_t get3dWithHolesDiamondPattern(){return 3;};
	bool IsWithInTheColumnRadius(Float_t relCellPosX,Float_t relCellPosY);
	bool IsOnTheEdgeOfCell(Float_t relCellPosX,Float_t relCellPosY,Float_t minDistanceToEdge);
	bool IsOnTheEdgeOfCell(Float_t relCellPosX,Float_t relCellPosY){return IsOnTheEdgeOfCell(relCellPosX,relCellPosY,minimumEdgeDistance);}
	bool IsOnTheEdgeOfCell(pair<Float_t,Float_t> relCellPos, Float_t minDistanceToEdge){return  IsOnTheEdgeOfCell(relCellPos.first,relCellPos.second,minDistanceToEdge);}
	bool IsOnTheEdgeOfCell(pair<Float_t,Float_t> relCellPos){return IsOnTheEdgeOfCell(relCellPos,minimumEdgeDistance);}
	Float_t GetMinimumEdgeDistance(){return minimumEdgeDistance;};
	bool IsPaperMode(){return bPaperMode;}
	Float_t getAdcToElectronConversion(){return 1;}

	ClassDef(TSettings,6)

};
#endif
