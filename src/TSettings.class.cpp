/*
 * TSettings.class.cpp
 *
 *  Created on: 02.08.2011
 *      Author: Felix Bachmair
 *
 *
 *      todo: Write copy constructor
 */

#include "TSettings.class.hh"
ClassImp(TSettings);
using namespace std;

TSettings::TSettings(UInt_t runNumber){
	verbosity=1;
	if(verbosity)
		cout<<"TSettings:Create TSettings-member with file:\""<<fileName<<"\""<<endl;
	diamondMapping=0;
	DefaultLoadDefaultSettings();
//	diamondMapping.PrintMapping();
	SetFileName("SETTINGS.new.ini");
	this->runNumber=runNumber;
  sys = gSystem;
	path = sys->pwd();
	createSettingsRootFile();
}

TSettings::TSettings(string fileName,UInt_t runNumber){
	verbosity=1;
	if(verbosity)
		cout<<"TSettings:Create TSettings-member with file:\""<<fileName<<"\""<<endl;
	diamondMapping=0;
	DefaultLoadDefaultSettings();
	this->runNumber=runNumber;
	sys = gSystem;
	path = sys->pwd();
  SetFileName(fileName);
	createSettingsRootFile();
}

TSettings::~TSettings(){

  saveSettings();
  settingsFile->Close();
}

void TSettings::SetFileName(string newFileName){
	if(verbosity)
		cout<<"TSettings::SetFileName:\""<<newFileName<<"\""<<endl;
	fileName=newFileName;
	LoadSettings();
}

void TSettings::LoadSettings(){
	if (verbosity)
		cout << endl << "TSettings::Overriding default settings with settings from \"" <<this->fileName<<"\""<< endl << endl;

	ifstream file(fileName.c_str());
	if(!file) {
		cout << "TSettings::An error has encountered while trying to open file " << fileName << endl;
		cout << "TSettings::Keeping default settings; no channels will be screened." << endl;
		return;
	}
	else cout <<"TSettings::"<< fileName << " successfully opened." << endl << endl;


	while(!file.eof()) {

		//get next line
		string line;
		getline(file,line);

		//check if comment or empty line: comments: ";","#","/"
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
	    if(key=="single_channel_analysis_channels") {
			 cout << key.c_str() << " = " << value.c_str() << endl;
			 ParseIntArray(value,single_channel_analysis_channels);
		  }
		if(key=="single_channel_analysis_enable") {
		         cout << key.c_str() << " = " << value.c_str() << endl;
		         single_channel_analysis_enable = (int)strtod(value.c_str(),0);
		      }
		  if(key=="single_channel_analysis_eventwindow") {
			 cout << key.c_str() << " = " << value.c_str() << endl;
			 single_channel_analysis_eventwindow = (int)strtod(value.c_str(),0);
		  }
		/*if(key=="store_threshold") {//TODO It's needed in settings reader
	         cout << key.c_str() << " = " << value.c_str() << endl;
	        store_threshold = (float)strtod(value.c_str(),0);
	      }*/
		  if(key=="CMN_corr_low") {
		  		         cout << key.c_str() << " = " << value.c_str() << endl;
		  		         CMN_corr_low = (int)strtod(value.c_str(),0);
		  		      }
		  if(key=="CMN_corr_high") {
		         cout << key.c_str() << " = " << value.c_str() << endl;
		         CMN_corr_high = (int)strtod(value.c_str(),0);
		      }
		if(key=="CMN_cut") {
		    cout << key.c_str() << " = " << value.c_str() << endl;
		    CMN_cut = (int)strtod(value.c_str(),0);
		}
		if(key=="DO_CMC") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			DO_CMC = (int)strtod(value.c_str(),0);
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
			cout << key.c_str() << " = " << value.c_str() << " size:";
			ParseIntArray(value,Det_channel_screen_channels[8]);
			cout<<Det_channel_screen_channels[8].size()<<endl;
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
		if (key == "UseAutoFidCut") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			UseAutoFidCut = (bool)strtod(value.c_str(),0);
		}
		if (key == "AlternativeClustering") {
			cout << key.c_str() << " = " << value.c_str() << endl;
			AlternativeClustering = (bool)strtod(value.c_str(),0);
		}
		 if(key=="store_threshold") {
			 cout << key.c_str() << " = " << value.c_str() << endl;
			 store_threshold = (float)strtod(value.c_str(),0);
		  }
	      if(key=="plotChannel_on") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         plotChannel_on = (int)strtod(value.c_str(),0);
	      }
	      if(key=="SingleChannel2000plots") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  SingleChannel2000plots = (int)strtod(value.c_str(),0);
	      }
	      if(key=="makeDiamondPlots") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  makeDiamondPlots = (int)strtod(value.c_str(),0);
	      }
	      if(key=="makeHits2D") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  makeHits2D = (int)strtod(value.c_str(),0);
	      }
	      if(key=="makeNoise2D") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  makeNoise2D = (int)strtod(value.c_str(),0);
	      }
	      if(key=="makePullDist") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  makePullDist = (int)strtod(value.c_str(),0);
	      }
	      if(key=="makePedRMSTree") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  makePedRMSTree = (int)strtod(value.c_str(),0);
	      }
	      if(key=="eventPrintHex") {
	    	  cout << key.c_str() << " = " << value.c_str() << endl;
	    	  eventPrintHex = (int)strtod(value.c_str(),0);
	      }

	      if(key=="plottedChannel") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         //plottedChannel = (int)strtod(value.c_str(),0);
	      }

	      if(key=="high_rms_cut") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         high_rms_cut = (int)strtod(value.c_str(),0);
	      }
	      if(key=="rms_cut") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         rms_cut = (float)strtod(value.c_str(),0);
	      }
	      if(key=="zoomDiamondPlots") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         zoomDiamondPlots = (int)strtod(value.c_str(),0);
	      }
	      if(key=="singleTrack2D") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         singleTrack2D = (int)strtod(value.c_str(),0);
	      }
	      if(key=="singleTrack2DmaxClusterSize") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         singleTrack2DmaxClusterSize = (int)strtod(value.c_str(),0);
	      }
	      if(key=="maxNoise2D") {
	         cout << key.c_str() << " = " << value.c_str() << endl;
	         maxNoise2D = (float)strtod(value.c_str(),0);
	      }
	      if(key=="clusterHitFactors") {
	    	  cout<<key<< " = "<< value.c_str() <<endl;
	    	 ParseFloatArray(value,clusterHitFactors);
	      }
	      if(key=="clusterSeedFactors") {
	    	  cout<<key<< " = "<< value.c_str() <<endl;
	    	  ParseFloatArray(value,clusterSeedFactors);
	      }
	      if(key=="diamondMapping") {
	    	  cout<<key<<" = "<<value.c_str()<<endl;
	    	  std::vector<int>vecDiaMapping;
	    	  ParseIntArray(value,vecDiaMapping);
	    	  if(diamondMapping==0)
	    		  delete diamondMapping;
	    	  diamondMapping=new TChannelMapping(vecDiaMapping);
	    	  diamondMapping->PrintMapping();
	    	  cout<<diamondMapping<<endl;
	    		getDetChannelNo(0);
	      }
	}

	file.close();

	for(int det=0; det<9; det++) {
		this->Det_channel_screen[det].ScreenChannels(this->getDet_channel_screen_channels(det));
		//this->getDet_channel_screen(det).ScreenRegions(this->getDet_channel_screen_regions(det));
		cout<<"Detector "<<det<<" screened channels: ";
		this->getDet_channel_screen(det).PrintScreenedChannels();
		cout<<endl;
	}

	for(int det=0;det<9;det++){
		cout<<"analyse detector "<<det<< " with "<<getClusterSeedFactor(det)<<"/"<<getClusterHitFactor(det)<<endl;
	}

	cout<<endl<<"TSettings::Finished importing settings from "<<fileName<<endl<<endl;
}

void TSettings::DefaultLoadDefaultSettings(){
	if(verbosity)
		cout<<"TSettings::LoadDefaultSettings"<<endl;
	//default general settings
	SaveAllFilesSwitch = 1; //1 for save files, 0 for don't
	ClosePlotsOnSave = 1;
	IndexProduceSwitch = 1;
	store_threshold=2;

	//default pedestal settings
	fix_dia_noise = -1;//7.7; // fix_dia_noise<0 disables diamond noise-fixing
	dia_input = 0; // 1 for 2006 and 0 for the rest
	DO_CMC = 1;
	CMN_cut = 4;  //Should be less than or equal to CMN_coor_high
	Iter_Size = 500; //buffer size
	Taylor_speed_throttle = 1000; //# of events to recalculate RMS the old way; set to 1 to disable
	CMN_corr_high=7;
	CMN_corr_low=3;

	alignment_chi2=1.0;
	alignment_training_track_fraction=0.25;
	alignment_training_track_number=10000;
	trainingMethod=enumFraction;

	//default clustering settings
	snr_plots_enable = 0;

	Di_Cluster_Seed_Factor = 10;
	Di_Cluster_Hit_Factor = 7;

	Si_Cluster_Seed_Factor = 5;
	Si_Cluster_Hit_Factor = 3;

	Si_Pedestal_Hit_Factor = 5;
	Di_Pedestal_Hit_Factor = 5;

	si_avg_fidcut_xlow = 90;
	si_avg_fidcut_xhigh = 165;
	si_avg_fidcut_ylow = 80;
	si_avg_fidcut_yhigh = 160;

	pulse_height_num_bins = 300;
	pulse_height_si_max = 300;
	pulse_height_di_max = 3000;

	snr_distribution_si_max = 2500;
	snr_distribution_di_max = 2500;
	UseAutoFidCut = 0;
	AlternativeClustering = 0;

	//Hi/low eta slices
	eta_lowq_slice_low = 600;
	eta_lowq_slice_hi = 700;
	eta_hiq_slice_low = 1200;
	eta_hiq_slice_hi = 1500;

	//Number of slices (<1 to disable)
	etavsq_n_landau_slices = 0;

	plotChannel_on = 0; //make RMS Difference plot for all detectors, and Buffer Noise plots for D0X
	plotDiamond = 1; //make Buffer Noise plots for the diamond instead
	makeBufferPlots = 0; //make Buffer Plot whenever sigma and rms differ by rms_sigma_difference_cut
	//NOTE: only works if plotChannel_on = 1 and plottedChannel < 256
	SingleChannel2000plots = 0; //make SC_Pedestal plots for all silicon detectors and channels
	makeDiamondPlots = 0; //make DC_Pedestal plots for all diamond channels
	makeHits2D = 0; //make 2D histogram of hits and seeds
	makeNoise2D = 0; //make 2D histogram of noise per channel
	makePullDist = 0; //make pull distribution
	makePedRMSTree = 0; //make .root file of pedestal and rms values
	eventPrintHex = 10000; //print hex (should match .rz data)

	res_keep_factor=2;

	maxBufferPlots = 100;
	rms_sigma_difference_cut = 0.3;

	high_rms_cut = 1; //cut on absolute rms value instead of comparing to Gaussian
	rms_cut = 20.; //value to use if high_rms_cut


	zoomDiamondPlots = 0; //zoom in on DC_Pedestal (100 event / window)

	singleTrack2D = 1; //plot single tracks only in 2D hits histogram
	singleTrack2DmaxClusterSize = 2; //max size of clusters in silicon track (cluster = Di_Hit_Factor hits; no check for seeds/shoulders)

	maxNoise2D = 20.; //highest noise value plotted in 2D noise histogram
	single_channel_analysis_enable=false;
	   //default settings
	single_channel_analysis_eventwindow=5000; // Number of events to put in each histogram
	plottedChannel=256; //256 = enter channel on run. also, set to 256 and type 256 to turn off buffer noise plots
	UInt_t nDiaChannels=128;
	diamondMapping=new TChannelMapping(nDiaChannels);
	getDetChannelNo(0);
	cout<<"Print DefaultMapping:"<<endl;
//	diamondMapping.PrintMapping();
	cout<<"DONE"<<endl;
}

void TSettings::ParseFloatArray(string value, vector<float> &vec) {
	int index=0;
	string::size_type offset1 = value.find_first_of('{')+1;
	string::size_type offset2 = value.find_first_of(',');

	vec.push_back((float)strtod(value.substr(offset1,offset2-offset1).c_str(),0));

	//cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(offset1,offset2-offset1)="<<value.substr(offset1,offset2-offset1)<<"\tstrtod(value.substr(offset1,offset2-offset1),0)="<<strtod(value.substr(offset1,offset2-offset1).c_str(),0)<<endl;//check

	value = value.substr(offset2+1,value.length()-(offset2+1));

	while(value.length()>2) {
		offset2 = TMath::Min(value.find_first_of(','),value.find_first_of('}'));
		vec.push_back((float)strtod(value.substr(0,offset2).c_str(),0));
		index++;
		//cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(0,offset2)="<<value.substr(0,offset2)<<"\tstrtod(value.substr(0,offset2),0)="<<strtod(value.substr(0,offset2).c_str(),0)<<endl;//check
		if(value.find_first_of(';')<2) break;
		value = value.substr(offset2+1,value.length()-(offset2+1));
	}
}

/**
 * TODO Ueberarbeiten so dass beliebige werte fuer jeden detector eingetragen werden koennen
 * @param det
 * @return
 */
Float_t TSettings::getClusterSeedFactor(UInt_t det){
//	cout<<"get Cluster Seed Factor: "<<det<<" "<<clusterSeedFactors.size()<<endl;
	if(det<clusterSeedFactors.size())
		return clusterSeedFactors.at(det);
	if(det==8)
		return getDi_Cluster_Seed_Factor();
	return getSi_Cluster_Seed_Factor();
}

Float_t TSettings::getClusterHitFactor(UInt_t det){

	if(clusterHitFactors.size()>det)
		return clusterHitFactors.at(det);
	if(det==8)
		return getDi_Cluster_Hit_Factor();
	return getSi_Cluster_Hit_Factor();
}

Int_t TSettings::getMakeBufferPlots() const
{
    return makeBufferPlots;
}

Int_t TSettings::getPlotDiamond() const
{
    return plotDiamond;
}

void TSettings::setMakeBufferPlots(Int_t makeBufferPlots)
{
    this->makeBufferPlots = makeBufferPlots;
}

void TSettings::setPlotDiamond(Int_t plotDiamond)
{
    this->plotDiamond = plotDiamond;
}

Int_t TSettings::getEventPrintHex() const
{
    return eventPrintHex;
}

Int_t TSettings::getMakeDiamondPlots() const
{
    return makeDiamondPlots;
}

Int_t TSettings::getMakeHits2D() const
{
    return makeHits2D;
}

Int_t TSettings::getMakeNoise2D() const
{
    return makeNoise2D;
}

Int_t TSettings::getMakePedRmsTree() const
{
    return makePedRMSTree;
}

Int_t TSettings::getMakePullDist() const
{
    return makePullDist;
}

Int_t TSettings::getSingleChannel2000plots() const
{
    return SingleChannel2000plots;
}

void TSettings::setEventPrintHex(Int_t eventPrintHex)
{
    this->eventPrintHex = eventPrintHex;
}

void TSettings::setMakeDiamondPlots(Int_t makeDiamondPlots)
{
    this->makeDiamondPlots = makeDiamondPlots;
}

void TSettings::setMakeHits2D(Int_t makeHits2D)
{
    this->makeHits2D = makeHits2D;
}

void TSettings::setMakeNoise2D(Int_t makeNoise2D)
{
    this->makeNoise2D = makeNoise2D;
}

void TSettings::setMakePedRmsTree(Int_t makePedRmsTree)
{
    makePedRMSTree = makePedRmsTree;
}

void TSettings::setMakePullDist(Int_t makePullDist)
{
    this->makePullDist = makePullDist;
}

void TSettings::setSingleChannel2000plots(Int_t singleChannel2000plots)
{
    SingleChannel2000plots = singleChannel2000plots;
}

UInt_t TSettings::getPlottedChannel() const
{
    return plottedChannel;
}

void TSettings::setPlottedChannel(UInt_t plottedChannel)
{
    this->plottedChannel = plottedChannel;
}

Int_t TSettings::getMaxBufferPlots() const
{
    return maxBufferPlots;
}

void TSettings::setMaxBufferPlots(Int_t maxBufferPlots)
{
    this->maxBufferPlots = maxBufferPlots;
}

Float_t TSettings::getRmsSigmaDifferenceCut() const
{
    return rms_sigma_difference_cut;
}

void TSettings::setRmsSigmaDifferenceCut(Float_t rmsSigmaDifferenceCut)
{
    rms_sigma_difference_cut = rmsSigmaDifferenceCut;
}

Int_t TSettings::getHighRmsCut() const
{
    return high_rms_cut;
}

Float_t TSettings::getRmsCut() const
{
    return rms_cut;
}

void TSettings::setHighRmsCut(Int_t highRmsCut)
{
    high_rms_cut = highRmsCut;
}

void TSettings::setRmsCut(Float_t rmsCut)
{
    rms_cut = rmsCut;
}

Float_t TSettings::getMaxNoise2D() const
{
    return maxNoise2D;
}

Int_t TSettings::getSingleTrack2D() const
{
    return singleTrack2D;
}

Int_t TSettings::getSingleTrack2DmaxClusterSize() const
{
    return singleTrack2DmaxClusterSize;
}

Int_t TSettings::getZoomDiamondPlots() const
{
    return zoomDiamondPlots;
}

void TSettings::setMaxNoise2D(Float_t maxNoise2D)
{
    this->maxNoise2D = maxNoise2D;
}

void TSettings::setSingleTrack2D(Int_t singleTrack2D)
{
    this->singleTrack2D = singleTrack2D;
}

void TSettings::setSingleTrack2DmaxClusterSize(Int_t singleTrack2DmaxClusterSize)
{
    this->singleTrack2DmaxClusterSize = singleTrack2DmaxClusterSize;
}

void TSettings::setZoomDiamondPlots(Int_t zoomDiamondPlots)
{
    this->zoomDiamondPlots = zoomDiamondPlots;
}

bool TSettings::isDet_channel_screened(UInt_t det, UInt_t ch)
{
	if(det<9&&ch<256)
		return this->Det_channel_screen[det].CheckChannel(ch);
	else
		return true;
}

TSettings::enumAlignmentTrainingMethod TSettings::getTrainingMethod() const
{
    return trainingMethod;
}

void TSettings::setTrainingMethod(enumAlignmentTrainingMethod trainingMethod)
{
    this->trainingMethod = trainingMethod;
}

UInt_t TSettings::getDetChannelNo(UInt_t vaCh)
{
	UInt_t detCh;
	detCh = this->diamondMapping->getDetChannelNo(vaCh);
	return detCh;
}

UInt_t TSettings::getVaChannelNo(UInt_t detChNo)
{

	UInt_t vaCh;
	vaCh = this->diamondMapping->getVAChannelNo(detChNo);
	return vaCh;
}

//todo
void TSettings::saveSettings()
{
  cout<<"SAVE SETTINGS TO ROOT FILE"<<endl;
  settingsFile->cd();
//  this->Write();
}


//todo
void TSettings::compareSettings()
{
  cout<<"compareSettings"<<endl;
}

//todo
void TSettings::createSettingsRootFile()
{
  stringstream name;
  name << path<< "/"<<this->runNumber<<"/Settings."<<this->runNumber<<".root";
//  settingsFile = TFile::Open(name.str().c_str());
  cout<<"Open settings from root file: "<<endl;
  settingsFile = new TFile(name.str().c_str(),"Update");
  if(settingsFile->IsZombie()){
    cout<<"file does not exist create new one!"<<endl;
    delete settingsFile;
    settingsFile = new TFile(name.str().c_str(),"RECREATE");
  }
  else{
    compareSettings();
  }


}

void TSettings::Print()
{
}

void TSettings::loadSettingsFromRootFile()
{
}

void TSettings::ParseIntArray(string value, vector<int> &vec) {

	int index=0;
	string::size_type offset1 = value.find_first_of('{')+1;
	string::size_type offset2 = value.find_first_of(',');

	vec.push_back((int)strtod(value.substr(offset1,offset2-offset1).c_str(),0));

	//cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(offset1,offset2-offset1)="<<value.substr(offset1,offset2-offset1)<<"\tstrtod(value.substr(offset1,offset2-offset1),0)="<<strtod(value.substr(offset1,offset2-offset1).c_str(),0)<<endl;//check

	value = value.substr(offset2+1,value.length()-(offset2+1));

	while(value.length()>2) {
		offset2 = TMath::Min(value.find_first_of(','),value.find_first_of('}'));
		vec.push_back((int)strtod(value.substr(0,offset2).c_str(),0));
		index++;
		//cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(0,offset2)="<<value.substr(0,offset2)<<"\tstrtod(value.substr(0,offset2),0)="<<strtod(value.substr(0,offset2).c_str(),0)<<endl;//check
		if(value.find_first_of(';')<2) break;
		value = value.substr(offset2+1,value.length()-(offset2+1));
	}
}


Float_t TSettings::getAlignment_chi2() const
{
	return alignment_chi2;
}

void TSettings::setAlignment_chi2(Float_t alignment_chi2)
{
	this->alignment_chi2 = alignment_chi2;
}

float TSettings::getFix_dia_noise() const
{
	return fix_dia_noise;
}

void TSettings::setFix_dia_noise(float fix_dia_noise)
{
	this->fix_dia_noise = fix_dia_noise;
}
Int_t TSettings::getIter_Size() const
{
	return Iter_Size;
}

void TSettings::setIter_Size(Int_t Iter_Size)
{
	this->Iter_Size = Iter_Size;
}

Int_t TSettings::getTaylor_speed_throttle() const
{
	return Taylor_speed_throttle;
}

void TSettings::setTaylor_speed_throttle(Int_t Taylor_speed_throttle)
{
	this->Taylor_speed_throttle = Taylor_speed_throttle;
}

Float_t TSettings::getDi_Pedestal_Hit_Factor() const
{
	return Di_Pedestal_Hit_Factor;
}

Int_t TSettings::getDia_input() const
{
	return dia_input;
}

Float_t TSettings::getSi_Pedestal_Hit_Factor() const
{
	return Si_Pedestal_Hit_Factor;
}

void TSettings::setDi_Pedestal_Hit_Factor(Float_t Di_Pedestal_Hit_Factor)
{
	this->Di_Pedestal_Hit_Factor = Di_Pedestal_Hit_Factor;
}

void TSettings::setDia_input(Int_t dia_input)
{
	this->dia_input = dia_input;
}

void TSettings::setSi_Pedestal_Hit_Factor(Float_t Si_Pedestal_Hit_Factor)
{
	this->Si_Pedestal_Hit_Factor = Si_Pedestal_Hit_Factor;
}


Int_t TSettings::getCMN_cut() const
{
	return CMN_cut;
}

Int_t TSettings::getDO_CMC() const
{
	return DO_CMC;
}

Float_t TSettings::getDi_Cluster_Hit_Factor() const
{
	return Di_Cluster_Hit_Factor;
}

Float_t TSettings::getDi_Cluster_Seed_Factor() const
{
	return Di_Cluster_Seed_Factor;
}

Float_t TSettings::getSi_Cluster_Hit_Factor() const
{
	return Si_Cluster_Hit_Factor;
}

Float_t TSettings::getSi_Cluster_Seed_Factor() const
{
	return Si_Cluster_Seed_Factor;
}

Float_t TSettings::getSi_avg_fidcut_xhigh() const
{
	return si_avg_fidcut_xhigh;
}

Float_t TSettings::getSi_avg_fidcut_xlow() const
{
	return si_avg_fidcut_xlow;
}

Float_t TSettings::getSi_avg_fidcut_yhigh() const
{
	return si_avg_fidcut_yhigh;
}

Float_t TSettings::getSi_avg_fidcut_ylow() const
{
	return si_avg_fidcut_ylow;
}

void TSettings::setCMN_cut(Int_t CMN_cut)
{
	this->CMN_cut = CMN_cut;
}

void TSettings::setDO_CMC(Int_t DO_CMC)
{
	this->DO_CMC = DO_CMC;
}

void TSettings::setDi_Cluster_Hit_Factor(Float_t Di_Cluster_Hit_Factor)
{
	this->Di_Cluster_Hit_Factor = Di_Cluster_Hit_Factor;
}

void TSettings::setDi_Cluster_Seed_Factor(Float_t Di_Cluster_Seed_Factor)
{
	this->Di_Cluster_Seed_Factor = Di_Cluster_Seed_Factor;
}

void TSettings::setSi_Cluster_Hit_Factor(Float_t Si_Cluster_Hit_Factor)
{
	this->Si_Cluster_Hit_Factor = Si_Cluster_Hit_Factor;
}

void TSettings::setSi_Cluster_Seed_Factor(Float_t Si_Cluster_Seed_Factor)
{
	this->Si_Cluster_Seed_Factor = Si_Cluster_Seed_Factor;
}

void TSettings::setSi_avg_fidcut_xhigh(Float_t si_avg_fidcut_xhigh)
{
	this->si_avg_fidcut_xhigh = si_avg_fidcut_xhigh;
}

void TSettings::setSi_avg_fidcut_xlow(Float_t si_avg_fidcut_xlow)
{
	this->si_avg_fidcut_xlow = si_avg_fidcut_xlow;
}

void TSettings::setSi_avg_fidcut_yhigh(Float_t si_avg_fidcut_yhigh)
{
	this->si_avg_fidcut_yhigh = si_avg_fidcut_yhigh;
}

void TSettings::setSi_avg_fidcut_ylow(Float_t si_avg_fidcut_ylow)
{
	this->si_avg_fidcut_ylow = si_avg_fidcut_ylow;
}

bool TSettings::getSingle_channel_analysis_enable()
{
	return this->single_channel_analysis_enable;
}

Int_t TSettings::getSingle_channel_analysis_eventwindow()
{
	return this->single_channel_analysis_eventwindow;
}

void TSettings::setSingle_channel_analysis_enable(bool singleChannelAnalysisEnable)
{
	this->single_channel_analysis_enable=singleChannelAnalysisEnable;
}

void TSettings::setSingle_channel_analysis_eventwindow(Int_t singleChannelAnalysisEventWindow)
{
	this->single_channel_analysis_eventwindow=singleChannelAnalysisEventWindow;
}

Int_t TSettings::getClosePlotsOnSave() const
{
	return ClosePlotsOnSave;
}

Int_t TSettings::getIndexProduceSwitch() const
{
	return IndexProduceSwitch;
}

Float_t TSettings::getPulse_height_di_max() const
{
	return pulse_height_di_max;
}

Int_t TSettings::getPulse_height_num_bins() const
{
	return pulse_height_num_bins;
}

Float_t TSettings::getPulse_height_si_max() const
{
	return pulse_height_si_max;
}

Float_t TSettings::getPulse_height_max(UInt_t det) const
{
	if (det < 8) return this->getPulse_height_si_max();
	if (det == 8) return this->getPulse_height_di_max();
	return -1;
}

Int_t TSettings::getSaveAllFilesSwitch() const
{
	return SaveAllFilesSwitch;
}

Float_t TSettings::getSnr_distribution_di_max() const
{
	return snr_distribution_di_max;
}

Float_t TSettings::getSnr_distribution_si_max() const
{
	return snr_distribution_si_max;
}

void TSettings::setClosePlotsOnSave(Int_t ClosePlotsOnSave)
{
	this->ClosePlotsOnSave = ClosePlotsOnSave;
}

void TSettings::setIndexProduceSwitch(Int_t IndexProduceSwitch)
{
	this->IndexProduceSwitch = IndexProduceSwitch;
}

void TSettings::setPulse_height_di_max(Float_t pulse_height_di_max)
{
	this->pulse_height_di_max = pulse_height_di_max;
}

void TSettings::setPulse_height_num_bins(Int_t pulse_height_num_bins)
{
	this->pulse_height_num_bins = pulse_height_num_bins;
}

void TSettings::setPulse_height_si_max(Float_t pulse_height_si_max)
{
	this->pulse_height_si_max = pulse_height_si_max;
}

void TSettings::setSaveAllFilesSwitch(Int_t SaveAllFilesSwitch)
{
	this->SaveAllFilesSwitch = SaveAllFilesSwitch;
}

void TSettings::setSnr_distribution_di_max(Float_t snr_distribution_di_max)
{
	this->snr_distribution_di_max = snr_distribution_di_max;
}

void TSettings::setSnr_distribution_si_max(Float_t snr_distribution_si_max)
{
	this->snr_distribution_si_max = snr_distribution_si_max;
}

Float_t TSettings::getEta_hiq_slice_hi() const
{
	return eta_hiq_slice_hi;
}

Float_t TSettings::getEta_hiq_slice_low() const
{
	return eta_hiq_slice_low;
}

Float_t TSettings::getEta_lowq_slice_hi() const
{
	return eta_lowq_slice_hi;
}

Float_t TSettings::getEta_lowq_slice_low() const
{
	return eta_lowq_slice_low;
}

void TSettings::setEta_hiq_slice_hi(Float_t eta_hiq_slice_hi)
{
	this->eta_hiq_slice_hi = eta_hiq_slice_hi;
}

void TSettings::setEta_hiq_slice_low(Float_t eta_hiq_slice_low)
{
	this->eta_hiq_slice_low = eta_hiq_slice_low;
}

void TSettings::setEta_lowq_slice_hi(Float_t eta_lowq_slice_hi)
{
	this->eta_lowq_slice_hi = eta_lowq_slice_hi;
}

void TSettings::setEta_lowq_slice_low(Float_t eta_lowq_slice_low)
{
	this->eta_lowq_slice_low = eta_lowq_slice_low;
}

vector<Float_t> TSettings::getAlignment_phi_offsets() const
{
	return alignment_phi_offsets;
}

vector<Float_t> TSettings::getAlignment_x_offsets() const
{
	return alignment_x_offsets;
}

vector<Float_t> TSettings::getAlignment_y_offsets() const
{
	return alignment_y_offsets;
}

vector<Float_t> TSettings::getAlignment_z_offsets() const
{
	return alignment_z_offsets;
}

Int_t TSettings::getEtavsq_n_landau_slices() const
{
	return etavsq_n_landau_slices;
}

Int_t TSettings::getSnr_plots_enable() const
{
	return snr_plots_enable;
}

void TSettings::setAlignment_phi_offsets(vector<Float_t> alignment_phi_offsets)
{
	this->alignment_phi_offsets = alignment_phi_offsets;
}

void TSettings::setAlignment_x_offsets(vector<Float_t> alignment_x_offsets)
{
	this->alignment_x_offsets = alignment_x_offsets;
}

void TSettings::setAlignment_y_offsets(vector<Float_t> alignment_y_offsets)
{
	this->alignment_y_offsets = alignment_y_offsets;
}

void TSettings::setAlignment_z_offsets(vector<Float_t> alignment_z_offsets)
{
	this->alignment_z_offsets = alignment_z_offsets;
}

void TSettings::setEtavsq_n_landau_slices(Int_t etavsq_n_landau_slices)
{
	this->etavsq_n_landau_slices = etavsq_n_landau_slices;
}

void TSettings::setSnr_plots_enable(Int_t snr_plots_enable)
{
	this->snr_plots_enable = snr_plots_enable;
}


Float_t TSettings::getAlignment_training_track_fraction() const
{
	return alignment_training_track_fraction;
}

ChannelScreen  TSettings::getDet_channel_screen(int i) {
	return Det_channel_screen[i];
}

vector<int> TSettings::getDet_channel_screen_channels(int i) const
{
	return Det_channel_screen_channels[i];
}

vector<int> TSettings::getDet_channel_screen_regions(int i) const
{
	return Det_channel_screen_regions[i];
}

void TSettings::setAlignment_training_track_fraction(Float_t alignment_training_track_fraction)
{
	this->alignment_training_track_fraction = alignment_training_track_fraction;
}

void TSettings::setDet_channel_screen(int i, ChannelScreen Det_channel_screen)
{
	this->Det_channel_screen[i] = Det_channel_screen;
}

void TSettings::setDet_channel_screen_channels(int i, vector<int> Det_channel_screen_channels)
{
	this->Det_channel_screen_channels[i] = Det_channel_screen_channels;
}

void TSettings::setDet_channel_screen_regions(int i, vector<int> Det_channel_screen_regions)
{
	this->Det_channel_screen_regions[i] = Det_channel_screen_regions;
}
bool TSettings::getAlternativeClustering() const
{
    return AlternativeClustering;
}

bool TSettings::getUseAutoFidCut() const
{
    return UseAutoFidCut;
}

void TSettings::setAlternativeClustering(bool AlternativeClustering)
{
    this->AlternativeClustering = AlternativeClustering;
}

void TSettings::setUseAutoFidCut(bool UseAutoFidCut)
{
    this->UseAutoFidCut = UseAutoFidCut;
}

string TSettings::getFileName() const
{
	return this->fileName;
}



Int_t TSettings::getCMN_corr_low()
{
	return this->CMN_corr_low;
}

void TSettings::setCMN_corr_low(Int_t CMN_corr_low)
{
	this->CMN_corr_low=CMN_corr_low;
}

Int_t TSettings::getCMN_corr_high()
{
	return this->CMN_corr_high;
}

void TSettings::setCMN_corr_high(Int_t CMN_corr_high)
{
	this->CMN_corr_high=CMN_corr_high;
}


std::vector<int> TSettings::getSingle_channel_analysis_channels()
{
	return this->single_channel_analysis_channels;
}

float TSettings::getStore_threshold()
{
	return this->store_threshold;
}


void TSettings::setStore_threshold(float storeThreshold)
{
	this->store_threshold=storeThreshold;
}
Int_t TSettings::getPlotChannelOn() const
{
    return plotChannel_on;
}

void TSettings::setPlotChannelOn(Int_t plotChannelOn)
{
    this->plotChannel_on = plotChannelOn;
}

UInt_t TSettings::getRunNumber(){
	return this->runNumber;
}

Float_t TSettings::getRes_keep_factor(){
	return this->res_keep_factor;
}


void Print(){

}
