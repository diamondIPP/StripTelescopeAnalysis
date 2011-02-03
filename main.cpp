/*
 * main.cpp
 * Diamond
 *
 * Created by Lukas Baeni on 19.01.11.
 */

//#include "TROOT.h"
//#include "Clustering.class.cpp"
#include "SlidingPedestal.class.cpp"
#include "Clustering.class.cpp"
#include <fstream>
#include <iostream>
#include "main.h"

using namespace std;

int main () {
	cout << "hello world!" << endl;
	initVariables();
	RunListOK = ReadRunList();
	
	cout << "RunListOK? " << RunListOK << endl;
	cout << "RUNNUMBER: " << RUNNUMBER << " NEVENTS: " << NEVENTS << " RUNDESCRIPTION: " << RUNDESCRIPTION << endl;
	
	cout << "Runnumbers ";
	for (int i = 0; i < RunParameters.size(); i++) {
		cout << RunParameters[i].RunNumber;
		if (i+1 < RunParameters.size()) cout << ", ";
	}
	cout << " will be analysed.." << endl;
	
	if (!RunListOK) return 0;
	
	for (int i = 0; i < RunParameters.size(); i++) {
		RunParameters[i].GetParameters();
		cout << endl << endl << endl << endl;
		cout << "====================================" << endl;
		cout << "==> Starting analysis.." << endl;
		cout << "====================================" << endl << endl;
		cout << "RUNNUMBER: " << RUNNUMBER << endl << "NEVENTS: " << NEVENTS << endl << "RUNDESCRIPTION: " << RUNDESCRIPTION << endl;
		cout << endl << endl << endl;
		
		if (DO_SLIDINGPEDESTAL) {
			SlidingPedestal sl(RUNNUMBER,RUNDESCRIPTION);
			sl.Slide(NEVENTS,INITIAL_EVENT,HIT_OCCUPANCY);
		}
//		return 0;
		Clustering cl(RUNNUMBER,RUNDESCRIPTION);
		if (DO_ALIGNMENT) {
			cl.Align(PLOTS, ALTERNATIVECLUSTERING);
		}
		else {
			cl.ClusterRun(PLOTS,ALTERNATIVECLUSTERING);
		}

		
	}
	return 0;
}

void initVariables() {
	RUNDESCRIPTION = "";
	NEVENTS = 10000;
	INITIAL_EVENT = 1000;
	HIT_OCCUPANCY = 0;
	PLOTS = 1;
	ALTERNATIVECLUSTERING = 0;
}

int ReadRunList() {
	RunInfo run;
	char RunDescription[200];
	int NEvents, Initial_Event;
	RunParameters.clear();
	cout << endl << "reading runlist.." << endl;
	ifstream file("RunList.ini");
	if (!file) {
		cout << "An error has encountered while trying to open RunList.ini" << endl;
		return 0;
	}
	else cout << "RunList.ini" << " successfully opened." << endl << endl;
	
	while (!file.eof()) {
		
		//get next line
		string line;
		getline(file,line);
		
		//check if comment or empty line
		if ((line.substr(0, 1) == ";") || (line.substr(0, 1) == "#") || (line.substr(0, 1) == "/") || line.empty()) {
			continue;
		}
		
		
		sscanf(line.c_str(), "%d %s %d %d %d %d", &RUNNUMBER, RunDescription, &NEvents, &Initial_Event, &DO_SLIDINGPEDESTAL, &DO_ALIGNMENT);
		if (NEvents != 0) NEVENTS = NEvents;
		if (Initial_Event != 0) INITIAL_EVENT = Initial_Event;
		cout << "RunDescription Char: " << RunDescription[0] << endl;
		if (RunDescription[0] != '0') RUNDESCRIPTION = RunDescription;
		
		run.SetParameters();
		RunParameters.push_back(run);
	}
}

/*void Clustering::LoadSettings() {
	
	cout<<endl<<"Overriding default settings with settings in Settings.ini"<<endl<<endl;
	
	ifstream file(settings_file.c_str());
	if(!file) {
		cout << "An error has encountered while trying to open file " << settings_file << endl;
		cout << "Keeping default settings; no channels will be screened." << endl;
		return;
	}
	else cout << settings_file << " successfully opened." << endl << endl;
	
	
	while(!file.eof()) {
		
		//get next line
		string line;
		getline(file,line);
		
		//check if comment or empty line
		if ((line.substr(0, 1) == ";") || (line.substr(0, 1) == "#") || (line.substr(0, 1) == "/") || line.empty()) {
			continue;
		}
		
		//find the index of first '=' character on the line
		string::size_type offsetl = line.find_first_of('=');
		string::size_type offsetr = line.find_first_of('=');
		
		//extract the key (LHS of the ini line)
		string key = line.substr(0, offsetl);
		
		//trim spaces from key
		while(line.at(offsetl-1)==' ') {
			offsetl--;
		}
		key = line.substr(0, offsetl);
		
		//extract the value (RHS of the ini line)
		string value = line.substr(offsetr+1, line.length()-(offsetr+1));
		
		//trim spaces from value
		while(line.at(offsetr+1)==' ') {
			offsetr++;
		}
		value = line.substr(offsetr+1, line.length()-(offsetr+1));
		
		//trim end ';' from end of key if found
		if(value.find_first_of(';')!=string::npos) {
			value = line.substr(offsetr+1, value.find_first_of(';'));//line.length()-(offsetr+1)-1);
		}
		
		//cant switch on strings so use if statements
		if(key=="SaveAllFilesSwitch") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			SaveAllFilesSwitch = (int)strtod(value.c_str(),0);
		}
		if(key=="ClosePlotsOnSave") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ClosePlotsOnSave = (int)strtod(value.c_str(),0);
		}
		if(key=="IndexProduceSwitch") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			IndexProduceSwitch = (int)strtod(value.c_str(),0);
		}
		if(key=="snr_plots_enable") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			snr_plots_enable = (int)strtod(value.c_str(),0);
		}
		if(key=="fix_dia_noise") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			fix_dia_noise = (int)strtod(value.c_str(),0);
		}
		if(key=="store_threshold") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			store_threshold = (float)strtod(value.c_str(),0);
		}
		if(key=="CMN_cut") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			CMN_cut = (int)strtod(value.c_str(),0);
		}
		if(key=="Iter_Size") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Iter_Size = (int)strtod(value.c_str(),0);
		}
		if(key=="Taylor_speed_throttle") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Taylor_speed_throttle = (int)strtod(value.c_str(),0);
		}
		if(key=="dia_input") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			dia_input = (int)strtod(value.c_str(),0);
		}
		if(key=="Si_Pedestal_Hit_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Si_Pedestal_Hit_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="Di_Pedestal_Hit_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Di_Pedestal_Hit_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="Si_Cluster_Seed_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Si_Cluster_Seed_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="Di_Cluster_Seed_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Di_Cluster_Seed_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="Si_Cluster_Hit_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Si_Cluster_Hit_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="Di_Cluster_Hit_Factor") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			Di_Cluster_Hit_Factor = (float)strtod(value.c_str(),0);
		}
		if(key=="eta_lowq_slice_low") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			eta_lowq_slice_low = (float)strtod(value.c_str(),0);
		}
		if(key=="eta_lowq_slice_hi") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			eta_lowq_slice_hi = (float)strtod(value.c_str(),0);
		}
		if(key=="eta_hiq_slice_low") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			eta_hiq_slice_low = (float)strtod(value.c_str(),0);
		}
		if(key=="eta_hiq_slice_hi") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			eta_hiq_slice_hi = (float)strtod(value.c_str(),0);
		}
		if(key=="etavsq_n_landau_slices") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			etavsq_n_landau_slices = (float)strtod(value.c_str(),0);
		}
		if(key=="alignment_x_offsets") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseFloatArray(value,alignment_x_offsets);
		}
		if(key=="alignment_y_offsets") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseFloatArray(value,alignment_y_offsets);
		}
		if(key=="alignment_phi_offsets") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseFloatArray(value,alignment_phi_offsets);
		}
		if(key=="alignment_z_offsets") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseFloatArray(value,alignment_z_offsets);
		}
		if(key=="alignment_training_track_fraction") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			alignment_training_track_fraction = (float)strtod(value.c_str(),0);
		}
		if(key=="D0X_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[0]);
		}
		if(key=="D0Y_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[1]);
		}
		if(key=="D1X_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[2]);
		}
		if(key=="D1Y_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[3]);
		}
		if(key=="D2X_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[4]);
		}
		if(key=="D2Y_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[5]);
		}
		if(key=="D3X_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[6]);
		}
		if(key=="D3Y_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[7]);
		}
		if(key=="Dia_channel_screen_channels") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_channels[8]);
		}
		if(key=="D0X_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[0]);
		}
		if(key=="D0Y_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[1]);
		}
		if(key=="D1X_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[2]);
		}
		if(key=="D1Y_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[3]);
		}
		if(key=="D2X_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[4]);
		}
		if(key=="D2Y_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[5]);
		}
		if(key=="D3X_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[6]);
		}
		if(key=="D3Y_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[7]);
		}
		if(key=="Dia_channel_screen_regions") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			ParseIntArray(value,Det_channel_screen_regions[8]);
		}
		if(key=="si_avg_fidcut_xlow") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			si_avg_fidcut_xlow = (int)strtod(value.c_str(),0);
		}
		if(key=="si_avg_fidcut_xhigh") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			si_avg_fidcut_xhigh = (int)strtod(value.c_str(),0);
		}
		if(key=="si_avg_fidcut_ylow") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			si_avg_fidcut_ylow = (int)strtod(value.c_str(),0);
		}
		if(key=="si_avg_fidcut_yhigh") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			si_avg_fidcut_yhigh = (int)strtod(value.c_str(),0);
		}
		if(key=="pulse_height_num_bins") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			pulse_height_num_bins = (int)strtod(value.c_str(),0);
		}
		if(key=="pulse_height_si_max") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			pulse_height_si_max = (int)strtod(value.c_str(),0);
		}
		if(key=="pulse_height_di_max") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			pulse_height_di_max = (int)strtod(value.c_str(),0);
		}
		if(key=="snr_distribution_si_max") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			snr_distribution_si_max = (int)strtod(value.c_str(),0);
		}
		if(key=="snr_distribution_di_max") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			snr_distribution_di_max = (int)strtod(value.c_str(),0);
		}
		if (key == "alignment_chi2") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			alignment_chi2 = (Float_t)strtod(value.c_str(),0);
		}
	}
	
	file.close();
	cout<<endl<<"Finished importing settings from "<<settings_file<<endl<<endl;
}*/