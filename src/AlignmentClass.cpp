/*
 * AlignmentClass.cpp
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#include "AlignmentClass.hh"

AlignmentClass::AlignmentClass() {
	// TODO Auto-generated constructor stub
	nAlignSteps =2;
	histSaver= new HistogrammSaver(0);

	/*histSaver->SetRunNumber(run_number);
	histSaver->SetNumberOfEvents((unsigned int)PedTree->GetEntries());
	histSaver->SetPlotsPath(plots_path);*/
}

AlignmentClass::~AlignmentClass() {
	// TODO Auto-generated destructor stub
}

int AlignmentClass::Align(bool plots, bool CutFakeTracksOn){
	if(tracks.size()==0) {
		cout<<"Clustering::Align: No tracks found; need to CallClustering::ClusterRun first..."<<endl;
		return -1;
		//ClusterRun(plots); // doesn't use alternative clustering
	}

	if (tracks.size() == 0) {
		cout << "No tracks available. Alignment not possible. (tracks.size() = " << tracks.size() << ")" << endl;
		return 0;
	}

	/*	int count_true = 0;
	int count_false = 0;
	for (int i = 0; i < tracks_mask.size(); i++) {
		cout << "track " << i << " has mask " << tracks_mask[i] << endl;
		if (tracks_mask[i]) count_true++;
		if (!tracks_mask[i]) count_false++;
	}
	cout << count_true << " tracks are masked as true and " << count_false << " tracks as false." << endl;
	return;*/

	//*string plots_path_save = plots_path;
	vector<TDiamondTrack> alignment_tracks = tracks;
	//*vector<bool> alignment_tracks_mask = tracks_mask;
	//*vector<TDiamondTrack> alignment_tracks_fidcut = tracks_fidcut;
	//*vector<bool> alignment_tracks_fidcut_mask = tracks_fidcut_mask;
	ostringstream plots_path_alignment, plots_path_alignment_CutFakeTracks;
	plots_path_alignment << plots_path << "/alignment/";
	plots_path_alignment_CutFakeTracks << plots_path << "/alignment_CutFakeTracks/";
	plots_path = plots_path_alignment.str();

	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
	for (int alignStep = 0; alignStep < nAlignSteps; alignStep++) {
		sys->mkdir(plots_path.c_str());
		//		TDetectorAlignment* align = new TDetectorAlignment(plots_path, tracks, tracks_mask);
		TDetectorAlignment* align = new TDetectorAlignment(plots_path, alignment_tracks, alignment_tracks_mask);

		Int_t nPasses = 10;
		Double_t plot_width_factor = 3; // scales the widths of the plots; range is a 3*width of distribution centered on mean

		align->PlotAngularDistribution(); //look at angular distribution of tracks
		align->PlotCartesianDistribution(); //look at slope distribution of tracks

		string prename = "alignment_PrealignmentResiduals";
		string postname = "alignment_PostalignmentResiduals";

		// generate residuals before alignment
		align->CheckDetectorAlignmentXYPlots(0, 1, 3, prename);
		align->CheckDetectorAlignmentXYPlots(1, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(2, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(3, 0, 2, prename);

		// itterative alignment loop
		for(int i=0; i<nPasses; i++) {
			cout << "\n\nPass " << i+1 << endl<< endl;
			//xy alignment
			align->CheckDetectorAlignmentXY(0, 1, 3);
			align->AlignDetectorXY(1, 0, 3);
			align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);
			//align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);
			//align->AlignDetectorXY(2, 0, 3);
			align->AlignDetectorXY(3, 0, 2);
			//align->AlignDetectorXY(2, 0, 3);
			//align->AlignDetectorXY(1, 0, 3);

			//phi alignment: this chunk of code causes seg faulting in code at top of loop!
			//align->AlignDetectorPhi(1, 0, 3, false, false);
			//align->AlignDetectorPhi(2, 0, 3, false, false);
			//align->AlignDetectorPhi(3, 0, 2, false, false);

			//phi alignment: this chunk of code causes seg faulting in code at top of loop!
			//align->AlignDetectorZ(1, 0, 3, false, false);
			//align->AlignDetectorZ(2, 0, 3, false, false);
			//align->AlignDetectorZ(3, 0, 2, false, false);
		}

		cout<<endl;
		cout<<endl;
		cout<<"Checking final residuals"<<endl;
		cout<<endl;

		// generate residuals after alignment
		align->CheckDetectorAlignmentXYPlots(0, 1, 3, postname);
		align->CheckDetectorAlignmentXYPlots(1, 0, 3, postname);
		align->CheckDetectorAlignmentXYPlots(2, 0, 3, postname);
		align->CheckDetectorAlignmentXYPlots(3, 0, 2, postname);

		cout<<endl;


		//Now align the diamond

		//load fidcut tracks w/ 1 diamond cluster
		//		align->LoadTracks(tracks_fidcut, tracks_fidcut_mask);
		align->LoadTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask);

		//check that the silicon is still aligned for these tracks_fidcut
		cout<<"Check that the telescope alignment still holds for fidcut tracks w/ single diamond cluster"<<endl;
		align->CheckDetectorAlignmentXY(0, 1, 3);
		align->CheckDetectorAlignmentXY(1, 0, 3);
		align->CheckDetectorAlignmentXY(2, 0, 3);
		align->CheckDetectorAlignmentXY(3, 0, 2);

		// generate residuals before alignment
		align->CheckDetectorAlignmentXYPlots(4, 1, 2, prename);

		// itterative alignment loop
		for(int i=0; i<5; i++) {
			cout << "\n\nPass " << i+1 << endl<< endl;
			//xy alignment
			align->AlignDetectorXY(4, 1, 2);
			//align->AlignDetectorXY(4, 0, 3);
			//align->AlignDetectorXY(4, 1, 3);
			//align->AlignDetectorXY(4, 0, 2);
			//align->AlignDetectorXY(4, 1, 2);
		}

		cout<<endl;
		cout<<endl;
		cout<<"Checking final diamond residuals"<<endl;
		cout<<endl;

		// generate residuals after alignment
		align->CheckDetectorAlignmentXYPlots(4, 1, 2, postname);


		//report results in a file

		ostringstream alignment_summary_path;
		alignment_summary_path << plots_path << "Alignment_Summary.txt";
		ofstream alignment_summary(alignment_summary_path.str().c_str());

		alignment_summary << "Alignment summary " << endl;
		alignment_summary << "----------------- " << endl << endl;

		alignment_summary << "Offsets (multiples of 50um):" << endl << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_x_offset[" << det << "] = " << align->det_x_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_y_offset[" << det << "] = " << align->det_y_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_z_offset[" << det << "] = " << align->det_z_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_phix_offset[" << det << "] = " << align->det_phix_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_phiy_offset[" << det << "] = " << align->det_phiy_offset[det] <<endl;
		}
		alignment_summary << endl;

		alignment_summary << "Resolutions (multiples of 50um):" << endl << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_x_resolution[" << det << "] = " << align->det_x_resolution[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << "det_y_resolution[" << det << "] = " << align->det_y_resolution[det] <<endl;
		}


		alignment_summary << endl << endl;
		alignment_summary << "Alignment summary (for pasting into a spread sheet) " << endl;
		alignment_summary << "--------------------------------------------------- " << endl << endl;

		alignment_summary << "Offsets (multiples of 50um):" << endl << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_x_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_y_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_z_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_phix_offset[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_phiy_offset[det] <<endl;
		}
		alignment_summary << endl;

		alignment_summary << "Resolutions (multiples of 50um):" << endl << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_x_resolution[det] <<endl;
		}
		alignment_summary << endl;
		for(int det=0; det<5; det++) {
			alignment_summary << align->det_y_resolution[det] <<endl;
		}
		alignment_summary.close();

		cout << "Intrinsic silicon resolution " << align->GetSiResolution() << " strips or " << align->GetSiResolution() * 50 << "um" << endl;

		if (CutFakeTracksOn && alignStep < 1) {
			align->CutFakeTracks(alignment_tracks, alignment_tracks_mask, alignment_chi2, CutFakeTracksOn, true);
			align->CutFakeTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, alignment_chi2, CutFakeTracksOn, true);
			plots_path = plots_path_alignment_CutFakeTracks.str();
		}
		if (!CutFakeTracksOn || alignStep == 1) {
			dia_offset.clear();
			//			dia_offset.push_back(align->GetXOffset(5));
			//			dia_offset.push_back(align->GetYOffset(5));
			//			dia_offset.push_back(align->GetZOffset(5));
			//			cout << "align->GetXOffset(4) = " << align->GetXOffset(4) << endl;
			//			cout << "align->GetXOffset(5) = " << align->GetXOffset(5) << endl;
			//			TransparentClustering(alignment_tracks, alignment_tracks_mask, align);
			this->TransparentClustering(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, align);
			break;
		}
	} // end alignment loop
	/*
   //Plot out the offsets
   for(Int_t plane=1; plane<4; plane++) {
      align->PlotXOffsetHistory(plane);
      align->PlotYOffsetHistory(plane);
      align->PlotZOffsetHistory(plane);
      align->PlotPhiOffsetHistory(plane);
   }

   //Print out the offsets
   for(Int_t plane=1; plane<4; plane++) {
      cout<<"Detector D"<<plane<<":\t";
      cout<<align->GetXOffset(plane)<<"\t";
      cout<<align->GetYOffset(plane)<<"\t";
      cout<<align->GetZOffset(plane)<<"\t";
      cout<<align->GetPhiOffset(plane)<<endl;
   }
	 */
}

void AlignmentClass::TransparentAnalysis(vector<TDiamondTrack> & tracks, vector<bool> & tracks_mask, TDetectorAlignment *align, bool verbose)
{ cout << "Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
int event;
vector<Float_t> x_positions, y_positions, z_positions, par, par_y;
Float_t diamond_hit_position = 0;
Float_t diamond_hit_y_position = 0;
int diamond_hit_channel = 0;
int diamond_secondhit_channel = 0;
int current_channel, current_sign;
Float_t cluster_adc = 0;
Float_t transp_eta = 0;
Float_t firstchannel_adc, secondchannel_adc;
vector<int> event_numbers;
Double_t diamond_x_offset = align->GetXOffset(3);
Double_t diamond_y_offset = align->GetYOffset(3);
Double_t diamond_phi_offset = align->GetPhiXOffset(3);
Double_t diamond_phi_y_offset = align->GetPhiYOffset(3);
Double_t diamond_z_position = 19.625; // TODO: is the z position of the diamond always 10.2??
Float_t eff_diamond_hit_position = 0;
int eff_diamond_hit_channel = 0;
Float_t hit_diff = 0;

cout << "diamond x offset: " << diamond_x_offset << endl;
cout << "diamond y offset: " << diamond_y_offset << endl;
cout << "diamond phi offset: " << diamond_phi_offset << endl;
cout << "diamond phi y offset: " << diamond_phi_y_offset << endl;
cout << "diamond hit factor: " << Di_Cluster_Hit_Factor << endl;

cout << "init histograms for transparent clustering.." << endl;
for (int i = 0; i < 10; i++) {
	ostringstream histoname_landau, histoname_eta;
	histoname_landau << "PulseHeight_D3X_" << (i+1) << "HitTransparClusters_8HitsFidcut";
	cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
	histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),pulse_height_num_bins,-0.5,pulse_height_di_max+0.5);
	//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
	//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
}
histo_transparentclustering_landau_mean = new TH1F("PulseHeightMeanVsChannels_Det?_TranspAna_","PulseHeightMeanVsChannels_Det?_TranspAna",10,0.5,10.5);
ostringstream histoname_eta;
histoname_eta << "Eta_D3X_2CentroidHits_TransparClusters";
histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_D3X_TransparClusters","DiffEstEffHit_D3X_TransparClusters", 200, -5.,5.);
histo_transparentclustering_hitdiff_scatter = new TH2F("DiffEstEffHit_Scatter_Dia_TransparClusters","DiffEstEffHit_Scatter_Dia_TransparClusters", 200, -5.,5.,128,0,127);
histo_transparentclustering_2Channel_PulseHeight = new TH1F("PulseHeight_D3X_2Channel_TranspCluster_8HitsFidcut","PulseHeight_D3X_2Channel_TranspCluster_8HitsFidcut",pulse_height_num_bins,-0.5,pulse_height_di_max+0.5);
cout << " done." << endl;

verbose = true;

PedFile = new TFile(pedfile_path.c_str());
PedTree = (TTree*)PedFile->Get("PedTree");
if (!PedTree)
{
	cerr << "PedTree not found!" << endl;
}

cout << PedTree->GetEntries() << " events in PedTree " << endl;

//Event Header Branches
PedTree->SetBranchAddress("RunNumber",&run_number);
PedTree->SetBranchAddress("EventNumber",&event_number);
PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
PedTree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
PedTree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);

//Telescope Data Branches
PedTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
PedTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
PedTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
PedTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
PedTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
PedTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
PedTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
PedTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
PedTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
PedTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
PedTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
PedTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
PedTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
PedTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
PedTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
PedTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
PedTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
PedTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
PedTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
PedTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
PedTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
PedTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
PedTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
PedTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
PedTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
PedTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
PedTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
PedTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
PedTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
PedTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
PedTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
PedTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
PedTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
PedTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
PedTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
PedTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
PedTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
PedTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
PedTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
PedTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
PedTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
PedTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
PedTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
PedTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
PedTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);

PedTree->GetEvent(0);
cout<< "Loaded first event in PedTree: "<<event_number<<endl;
cout<< "RunNumber is: "<<run_number<<endl;
cout<< "StoreThreshold is: "<<store_threshold<<endl;

event_numbers.clear();
for (int j = 0; j < PedTree->GetEntries(); j++) {
	PedTree->GetEvent(j);
	//		cout << endl << endl << endl;
	//		cout << "event " << j << " in PedTree has event_number: " << event_number << endl;
	event_numbers.push_back(event_number);
	//		for (int blablabla = 0; blablabla < 254; blablabla++) {
	//			cout << "Dia_ADC = " << Dia_ADC[blablabla] << ",\tDet_PedMean = " << Det_PedMean[8][blablabla] << endl;
	//			cout << "Collected charge in channel " << blablabla << " of diamond: " << Dia_ADC[blablabla]-Det_PedMean[8][blablabla] << endl;
	//		}
}
//	return;

// loop over tracks
for (int i = 0; i < tracks.size(); i++) {
	if (verbose) cout << " -- starting transparent clustering for track " << i << endl;

	// check if track is masked
	if (tracks[i].FakeTrack) {
		if (verbose) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
		continue;
	}

	// get event number for track
	if (verbose) cout << "Getting event number.. ";
	event = tracks[i].GetEventNumber();
	if (verbose) cout << " -> track " << i << " corresponds to event " << event << endl;

	// check if event number is valid
	if (event < 0) {
		cout << "Track " << i << " has no event number. Skipping this track.." << endl;
		continue;
	}

	// load data (apply offset)
	align->LoadData(tracks[i]);

	// read out x, y, z positions
	x_positions.clear();
	y_positions.clear();
	z_positions.clear();
	for (int det = 0; det < 3; det++) {
		x_positions.push_back(align->track_holder.GetD(det).GetX());
		y_positions.push_back(align->track_holder.GetD(det).GetY());
		z_positions.push_back(align->track_holder.GetD(det).GetZ());
		cout << "Det"<< det << " has z position " << align->track_holder.GetD(det).GetZ() << endl;
	}

	// read out effictive diamond hit position
	eff_diamond_hit_position = align->track_holder.GetD(3).GetX();

	// fit track
	par.clear();
	align->LinTrackFit(z_positions, x_positions, par);
	if (verbose) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;

	// fit y position
	par_y.clear();
	align->LinTrackFit(z_positions, y_positions, par_y);
	if (verbose) cout << "linear fit of track:\tpar_y[0] = " << par_y[0] << ",\tpar_y[1] = " << par_y[1] << endl;

	// estimate hit position in diamond
	//		diamond_z_position = align->track_holder.GetD(4).GetZ();
	diamond_hit_position = par[0] + par[1] * diamond_z_position;
	diamond_hit_y_position = par_y[0] + par_y[1] * diamond_z_position;
	diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
	diamond_hit_y_position = diamond_hit_y_position + diamond_y_offset;
	diamond_hit_position = (diamond_hit_position - 64) * TMath::Cos(diamond_phi_offset) - (diamond_hit_y_position - 64) * TMath::Sin(diamond_phi_offset) + 64; // add the tilt correction
	diamond_hit_position += 0.5; // added 0.5 to take the middle of the channel instead of the edge
	diamond_hit_channel = (int)diamond_hit_position;
	if (verbose) cout << "z position of diamond is " << diamond_z_position << endl;

	// difference between estimated and effective hit in diamond
	hit_diff = TMath::Abs(eff_diamond_hit_position - diamond_hit_position);
	//		cout << "effective hit position in diamond:\t" << eff_diamond_hit_position << "\testimated position in diamond:\t" << diamond_hit_position << endl;

	//get event
	if (verbose) cout << "getting event " << event << ".." << endl;
	PedTree->GetEvent(event);
	if (verbose) cout << "event_number = " << event_number << endl;

	// find biggest hit in diamond
	eff_diamond_hit_channel = 0;
	for (int j = 0; j < 128; j++) {
		if (Det_ADC[6][j]-Det_PedMean[6][j] > (Det_ADC[6][eff_diamond_hit_channel]-Det_PedMean[6][eff_diamond_hit_channel])) {
			eff_diamond_hit_channel = j;
		}
	}
	histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
	histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position, diamond_hit_y_position);
	cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;

	// cluster diamond channels around estimated hit position
	//		cout << " --" << endl;
	if (verbose) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
	if (verbose) cout << "Det_ADC[6] = " << Det_ADC[6][diamond_hit_channel] << ",\tDet_PedMean = " << Det_PedMean[6][diamond_hit_channel] << endl;
	if (verbose) cout << "Collected charge in channel " << (int)diamond_hit_position << " of D3X: " << Det_ADC[6][(int)diamond_hit_position]-Det_PedMean[6][(int)diamond_hit_position] << endl;

	int sign;

	if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
	else sign = -1;

	// calculate eta for the two closest two channels to the estimated hit position
	diamond_secondhit_channel = diamond_hit_channel - sign;
	firstchannel_adc = Det_ADC[6][diamond_hit_channel]-Det_PedMean[6][diamond_hit_channel];
	secondchannel_adc = Det_ADC[6][diamond_secondhit_channel]-Det_PedMean[6][diamond_secondhit_channel];
	if (sign == 1) transp_eta = firstchannel_adc / (firstchannel_adc + secondchannel_adc);
	else transp_eta = secondchannel_adc / (firstchannel_adc + secondchannel_adc);
	histo_transparentclustering_eta->Fill(transp_eta);

	// fill pulse height histogram
	histo_transparentclustering_2Channel_PulseHeight->Fill(firstchannel_adc+secondchannel_adc);

	if (verbose) cout << "clusters for track " << i << ":" << endl;
	// loop over different cluster sizes
	for (int j = 0; j < 10; j++) {
		cluster_adc = 0;
		current_channel = diamond_hit_channel;
		if (verbose) cout << "selected channels for " << j+1 << " hit transparent cluster: ";
		current_sign = sign;
		// sum adc for n channel cluster
		for (int channel = 0; channel <= j; channel++) {
			current_channel = current_channel + current_sign * channel;
			current_sign = (-1) * current_sign;
			if (verbose) cout << current_channel;
			if (verbose) if (channel < j) cout << ", ";
			if (current_channel > 0 && current_channel < 128 /* && Dia_ADC[current_channel]-Det_PedMean[8][current_channel] > Di_Cluster_Hit_Factor*Det_PedWidth[8][current_channel]*/)
				cluster_adc = cluster_adc + Det_ADC[6][current_channel]-Det_PedMean[6][current_channel];
		}
		if (verbose) cout << endl;
		if (verbose) cout << "total charge of cluster: " << cluster_adc << endl;
		if (verbose) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
		if (current_channel <= 0 || current_channel >= 128) break;
		histo_transparentclustering_landau[j]->Fill(cluster_adc);
	} // end loop over cluster sizes
} // end loop over tracks

// save histograms
for (int i = 0; i < 10; i++) {
	histSaver->SaveHistogram(histo_transparentclustering_landau[i]);
}
histSaver->SaveHistogram(histo_transparentclustering_eta);
histSaver->SaveHistogram(histo_transparentclustering_hitdiff);
histSaver->SaveHistogram(histo_transparentclustering_hitdiff_scatter);
histSaver->SaveHistogram(histo_transparentclustering_2Channel_PulseHeight);

PedFile->Close();
}



void AlignmentClass::TransparentClustering(vector<TDiamondTrack> & tracks, vector<bool> & tracks_mask, TDetectorAlignment *align, bool verbose)
{
	cout << "Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
	int event;
	vector<Float_t> x_positions, y_positions, z_positions, par, par_y;
	Float_t diamond_hit_position = 0;
	Float_t diamond_hit_y_position = 0;
	int diamond_hit_channel = 0;
	int diamond_secondhit_channel = 0;
	int current_channel, current_sign;
	Float_t cluster_adc = 0;
	Float_t transp_eta = 0;
	Float_t firstchannel_adc, secondchannel_adc;
	vector<int> event_numbers;
	Double_t diamond_x_offset = align->GetXOffset(4);
	Double_t diamond_y_offset = align->GetYOffset(4);
	Double_t diamond_phi_offset = align->GetPhiXOffset(4);
	Double_t diamond_phi_y_offset = align->GetPhiYOffset(4);
	Double_t diamond_z_position = 10.2; // TODO: is the z position of the diamond always 10.2??
	Float_t eff_diamond_hit_position = 0;
	int eff_diamond_hit_channel = 0;
	Float_t hit_diff = 0;

	cout << "diamond x offset: " << diamond_x_offset << endl;
	cout << "diamond y offset: " << diamond_y_offset << endl;
	cout << "diamond phi offset: " << diamond_phi_offset << endl;
	cout << "diamond phi y offset: " << diamond_phi_y_offset << endl;
	cout << "diamond hit factor: " << Di_Cluster_Hit_Factor << endl;

	cout << "init histograms for transparent clustering.." << endl;
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_Dia_" << (i+1) << "ChannelsTransparAna_8HitsFidcut";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;


		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),pulse_height_num_bins,-0.5,pulse_height_di_max+0.5);
		//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
		//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	}

	histo_transparentclustering_landau_mean = new TH1F("PulseHeightMeanVsChannels_Dia_TranspAna_","PulseHeightMeanVsChannels_Dia_TranspAna",10,0.5,10.5);
	ostringstream histoname_eta;
	histoname_eta << "Eta_Dia_2CentroidHits_TransparClusters";
	histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
	histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_Dia_TransparClusters","DiffEstEffHit_Dia_TransparClusters", 200, -5.,5.);
	histo_transparentclustering_hitdiff_scatter = new TH2F("DiffEstEffHit_Scatter_Dia_TransparClusters","DiffEstEffHit_Scatter_Dia_TransparClusters", 200, -5.,5.,128,0,127);
	histo_transparentclustering_2Channel_PulseHeight = new TH1F("PulseHeight_Dia_2Channel_TranspCluster_8HitsFidcut","PulseHeight_Dia_2Channel_TranspCluster_8HitsFidcut",pulse_height_num_bins,-0.5,pulse_height_di_max+0.5);
	cout << " done." << endl;

	//	verbose = true;

	PedFile = new TFile(pedfile_path.c_str());
	PedTree = (TTree*)PedFile->Get("PedTree");
	if (!PedTree)
	{
		cerr << "PedTree not found!" << endl;
	}

	cout << PedTree->GetEntries() << " events in PedTree " << endl;

	//Event Header Branches
	PedTree->SetBranchAddress("RunNumber",&run_number);
	PedTree->SetBranchAddress("EventNumber",&event_number);
	PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
	PedTree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
	PedTree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);

	//Telescope Data Branches
	PedTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
	PedTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
	PedTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
	PedTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
	PedTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
	PedTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
	PedTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
	PedTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
	PedTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
	PedTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
	PedTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
	PedTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
	PedTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
	PedTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
	PedTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
	PedTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
	PedTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
	PedTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
	PedTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
	PedTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
	PedTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
	PedTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
	PedTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
	PedTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
	PedTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
	PedTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
	PedTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
	PedTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
	PedTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
	PedTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
	PedTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
	PedTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
	PedTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
	PedTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
	PedTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
	PedTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
	PedTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
	PedTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
	PedTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
	PedTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
	PedTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
	PedTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
	PedTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
	PedTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
	PedTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);

	PedTree->GetEvent(0);
	cout<< "Loaded first event in PedTree: "<<event_number<<endl;
	cout<< "RunNumber is: "<<run_number<<endl;
	cout<< "StoreThreshold is: "<<store_threshold<<endl;

	event_numbers.clear();
	for (int j = 0; j < PedTree->GetEntries(); j++) {
		PedTree->GetEvent(j);
		//		cout << endl << endl << endl;
		//		cout << "event " << j << " in PedTree has event_number: " << event_number << endl;
		event_numbers.push_back(event_number);
		//		for (int blablabla = 0; blablabla < 254; blablabla++) {
		//			cout << "Dia_ADC = " << Dia_ADC[blablabla] << ",\tDet_PedMean = " << Det_PedMean[8][blablabla] << endl;
		//			cout << "Collected charge in channel " << blablabla << " of diamond: " << Dia_ADC[blablabla]-Det_PedMean[8][blablabla] << endl;
		//		}
	}
	//	return;

	// loop over tracks
	for (int i = 0; i < tracks.size(); i++) {
		if (verbose) cout << " -- starting transparent clustering for track " << i << endl;

		// check if track is masked
		if (tracks[i].FakeTrack) {
			if (verbose) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
			continue;
		}

		// get event number for track
		if (verbose) cout << "Getting event number.. ";
		event = tracks[i].GetEventNumber();
		if (verbose) cout << " -> track " << i << " corresponds to event " << event << endl;

		// check if event number is valid
		if (event < 0) {
			cout << "Track " << i << " has no event number. Skipping this track.." << endl;
			continue;
		}

		// load data (apply offset)
		align->LoadData(tracks[i]);

		// read out x, y, z positions
		x_positions.clear();
		y_positions.clear();
		z_positions.clear();
		for (int det = 0; det < 4; det++) {
			x_positions.push_back(align->track_holder.GetD(det).GetX());
			y_positions.push_back(align->track_holder.GetD(det).GetY());
			z_positions.push_back(align->track_holder.GetD(det).GetZ());
			cout << "Det"<< det << " has z position " << align->track_holder.GetD(det).GetZ() << endl;
		}

		// read out effictive diamond hit position
		eff_diamond_hit_position = align->track_holder.GetD(4).GetX();

		// fit track
		par.clear();
		align->LinTrackFit(z_positions, x_positions, par);
		if (verbose) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;

		// fit y position
		par_y.clear();
		align->LinTrackFit(z_positions, y_positions, par_y);
		if (verbose) cout << "linear fit of track:\tpar_y[0] = " << par_y[0] << ",\tpar_y[1] = " << par_y[1] << endl;

		// estimate hit position in diamond
		//		diamond_z_position = align->track_holder.GetD(4).GetZ();
		diamond_hit_position = par[0] + par[1] * diamond_z_position;
		diamond_hit_y_position = par_y[0] + par_y[1] * diamond_z_position;
		diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
		diamond_hit_y_position = diamond_hit_y_position + diamond_y_offset;
		diamond_hit_position = (diamond_hit_position - 64) * TMath::Cos(diamond_phi_offset) - (diamond_hit_y_position - 64) * TMath::Sin(diamond_phi_offset) + 64; // add the tilt correction
		diamond_hit_position += 0.5; // added 0.5 to take the middle of the channel instead of the edge
		diamond_hit_channel = (int)diamond_hit_position;
		if (verbose) cout << "z position of diamond is " << diamond_z_position << endl;

		// difference between estimated and effective hit in diamond
		hit_diff = TMath::Abs(eff_diamond_hit_position - diamond_hit_position);
		//		cout << "effective hit position in diamond:\t" << eff_diamond_hit_position << "\testimated position in diamond:\t" << diamond_hit_position << endl;

		//get event
		if (verbose) cout << "getting event " << event << ".." << endl;
		PedTree->GetEvent(event);
		if (verbose) cout << "event_number = " << event_number << endl;

		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int j = 0; j < 128; j++) {
			if (Dia_ADC[j]-Det_PedMean[8][j] > (Dia_ADC[eff_diamond_hit_channel]-Det_PedMean[8][eff_diamond_hit_channel])) {
				eff_diamond_hit_channel = j;
			}
		}
		histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
		histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position, diamond_hit_y_position);
		cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;

		// cluster diamond channels around estimated hit position
		//		cout << " --" << endl;
		if (verbose) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbose) cout << "Dia_ADC = " << Dia_ADC[diamond_hit_channel] << ",\tDet_PedMean = " << Det_PedMean[8][diamond_hit_channel] << endl;
		if (verbose) cout << "Collected charge in channel " << (int)diamond_hit_position << " of diamond: " << Dia_ADC[(int)diamond_hit_position]-Det_PedMean[8][(int)diamond_hit_position] << endl;

		int sign;

		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;

		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = Dia_ADC[diamond_hit_channel]-Det_PedMean[8][diamond_hit_channel];
		secondchannel_adc = Dia_ADC[diamond_secondhit_channel]-Det_PedMean[8][diamond_secondhit_channel];
		if (sign == 1) transp_eta = firstchannel_adc / (firstchannel_adc + secondchannel_adc);
		else transp_eta = secondchannel_adc / (firstchannel_adc + secondchannel_adc);
		histo_transparentclustering_eta->Fill(transp_eta);

		// fill pulse height histogram
		histo_transparentclustering_2Channel_PulseHeight->Fill(firstchannel_adc+secondchannel_adc);

		if (verbose) cout << "clusters for track " << i << ":" << endl;
		// loop over different cluster sizes
		for (int j = 0; j < 10; j++) {
			cluster_adc = 0;
			current_channel = diamond_hit_channel;
			if (verbose) cout << "selected channels for " << j+1 << " hit transparent cluster: ";
			current_sign = sign;
			// sum adc for n channel cluster
			for (int channel = 0; channel <= j; channel++) {
				current_channel = current_channel + current_sign * channel;
				current_sign = (-1) * current_sign;
				if (verbose) cout << current_channel;
				if (verbose) if (channel < j) cout << ", ";
				if (current_channel > 0 && current_channel < 128 /* && Dia_ADC[current_channel]-Det_PedMean[8][current_channel] > Di_Cluster_Hit_Factor*Det_PedWidth[8][current_channel]*/)
					cluster_adc = cluster_adc + Dia_ADC[current_channel]-Det_PedMean[8][current_channel];
			}
			if (verbose) cout << endl;
			if (verbose) cout << "total charge of cluster: " << cluster_adc << endl;
			if (verbose) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
			if (current_channel <= 0 || current_channel >= 128) break;
			histo_transparentclustering_landau[j]->Fill(cluster_adc);
		} // end loop over cluster sizes
	} // end loop over tracks

	// save histograms
	for (int i = 0; i < 10; i++) {
		histo_transparentclustering_landau_mean->SetBinContent(i+1,histo_transparentclustering_landau[i]->GetMean()); // plot pulse hight means into a histogram
		histSaver->SaveHistogram(histo_transparentclustering_landau[i]);
	}
	histSaver->SaveHistogram(histo_transparentclustering_landau_mean);
	histSaver->SaveHistogram(histo_transparentclustering_eta);
	histSaver->SaveHistogram(histo_transparentclustering_hitdiff);
	histSaver->SaveHistogram(histo_transparentclustering_hitdiff_scatter);
	histSaver->SaveHistogram(histo_transparentclustering_2Channel_PulseHeight);

	PedFile->Close();
}




