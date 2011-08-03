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
using namespace std;


class TSettings {
public:
	TSettings(string fileName);
	virtual ~TSettings();

	Float_t getAlignment_chi2() const;
	void setAlignment_chi2(Float_t alignment_chi2);
	float getFix_dia_noise() const;
	Int_t getIter_Size() const;
	Int_t getTaylor_speed_throttle() const;
	Float_t getDi_Pedestal_Hit_Factor() const;
	Int_t getDia_input() const;
	Float_t getSi_Pedestal_Hit_Factor() const;
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
	Int_t getSaveAllFilesSwitch() const;
	Float_t getSnr_distribution_di_max() const;
	Float_t getSnr_distribution_si_max() const;
	Float_t getEta_hiq_slice_hi() const;
	Float_t getEta_hiq_slice_low() const;
	Float_t getEta_lowq_slice_hi() const;
	Float_t getEta_lowq_slice_low() const;
	vector<Float_t> getAlignment_phi_offsets() const;
	vector<Float_t> getAlignment_x_offsets() const;
	vector<Float_t> getAlignment_y_offsets() const;
	vector<Float_t> getAlignment_z_offsets() const;
	Int_t getEtavsq_n_landau_slices() const;
	Int_t getSnr_plots_enable() const;
	Float_t getAlignment_training_track_fraction() const;
	ChannelScreen getDet_channel_screen(int i);
	vector<int> getDet_channel_screen_channels(int i) const;
	vector<int> getDet_channel_screen_regions(int i) const;
	bool getAlternativeClustering() const;
	bool getUseAutoFidCut() const;

	void setAlignment_training_track_fraction(Float_t alignment_training_track_fraction);
	void setFix_dia_noise(float fix_dia_noise);
	void setIter_Size(Int_t Iter_Size);
	void setTaylor_speed_throttle(Int_t Taylor_speed_throttle);
	void setDi_Pedestal_Hit_Factor(Float_t Di_Pedestal_Hit_Factor);
	void setDia_input(Int_t dia_input);
	void setSi_Pedestal_Hit_Factor(Float_t Si_Pedestal_Hit_Factor);
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
	void setAlignment_phi_offsets(vector<Float_t> alignment_phi_offsets);
	void setAlignment_x_offsets(vector<Float_t> alignment_x_offsets);
	void setAlignment_y_offsets(vector<Float_t> alignment_y_offsets);
	void setAlignment_z_offsets(vector<Float_t> alignment_z_offsets);
	void setEtavsq_n_landau_slices(Int_t etavsq_n_landau_slices);
	void setSnr_plots_enable(Int_t snr_plots_enable);
	void setDet_channel_screen(int i,ChannelScreen Det_channel_screen);
	void setDet_channel_screen_channels(int i,vector<int> Det_channel_screen_channels);
	void setDet_channel_screen_regions(int i,vector<int> Det_channel_screen_regions);
	void setAlternativeClustering(bool AlternativeClustering);
	void setUseAutoFidCut(bool UseAutoFidCut);

protected:
private:
	void SetFileName(string fileName);
	void LoadSettings();
	void DefaultLoadDefaultSettings();
	void ParseFloatArray(string value, vector<float> & vec);
	void ParseIntArray(string value, vector<int> & vec);
private:
	string fileName;
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

	Int_t CMN_cut;
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
	vector<Float_t> alignment_x_offsets;
	vector<Float_t> alignment_y_offsets;
	vector<Float_t> alignment_phi_offsets;
	vector<Float_t> alignment_z_offsets;
	Float_t alignment_training_track_fraction;
	vector<int> Det_channel_screen_channels[9];
	vector<int> Det_channel_screen_regions[9];
	ChannelScreen Det_channel_screen[9];
    bool dia_x_aligned;
    bool eta_correction;
	Float_t alignment_chi2;
	bool UseAutoFidCut;
	bool AlternativeClustering;

private:
    //Filter tracks not in good fiducial region w/o bad strips
     Int_t align_sil_fid_xlow;
     Int_t align_sil_fid_xhi;
     Int_t align_sil_fid_ylow;
     Int_t align_sil_fid_yhi;

private: int verbosity;

};
#endif
