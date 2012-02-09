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

//ROOT libraries
#include "Rtypes.h"


#include "ChannelScreen.hh"


class TSettings {
public:
	TSettings(std::string fileName,UInt_t runNumber=0);
	virtual ~TSettings();

	enum enumAlignmentTrainingMethod{enumFraction, enumEvents};

	Float_t getClusterSeedFactor(UInt_t det);
	Float_t getClusterHitFactor(UInt_t det);
	Float_t getAlignment_chi2() const;
	void setAlignment_chi2(Float_t alignment_chi2);
	float getFix_dia_noise() const;
	Int_t getIter_Size() const;
	Int_t getTaylor_speed_throttle() const;
	Float_t getDi_Pedestal_Hit_Factor() const;
	Int_t getDia_input() const;
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
    enumAlignmentTrainingMethod getTrainingMethod() const;
    void setTrainingMethod(enumAlignmentTrainingMethod trainingMethod);
    void Print();
protected:
    float store_threshold;
private:
    void SetFileName(std::string fileName);
    void LoadSettings();
    void DefaultLoadDefaultSettings();
    void ParseFloatArray(std::string value, std::vector<float> & vec);
    void ParseIntArray(std::string value, std::vector<int> & vec);
private:
    std::string fileName;
private:
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
private:
    //Filter tracks not in good fiducial region w/o bad strips
    Int_t align_sil_fid_xlow;
    Int_t align_sil_fid_xhi;
    Int_t align_sil_fid_ylow;
    Int_t align_sil_fid_yhi;
private:
    int verbosity;
};
#endif
