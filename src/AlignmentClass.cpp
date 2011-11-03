/*
 * AlignmentClass.cpp
 *
 *  Created on: 31.07.2011
 *      Author: Felix Bachmair
 */

#include "AlignmentClass.hh"
using namespace std;

AlignmentClass::AlignmentClass(string fileName,UInt_t nEventNumber) {
	verbosity=1;
	PedFileName=fileName;
	align =NULL;
	eventReader=NULL;
	histSaver=NULL;
	histSaverAlignment=NULL;
	histSaverAlignmentFakedTracks=NULL;
	settings=NULL;
	dianoise_sigma[0]=1;
	dianoise_sigma[1]=1;
	nAlignSteps =2;
	histSaver= new HistogrammSaver(0);
	if (verbosity)cout<<"AlignmentClass::AlignmentClass: eventReader:loading:"<<endl;
	if (eventReader==NULL) eventReader = new TADCEventReader(PedFileName);
	if (eventReader==NULL || !eventReader->isOK()){
		cout<<"AlignmentClass::eventReader could not be initialized... EXIT PROGRAM"<<endl;
		exit (-1);
	}
	if (verbosity)cout<<"AlignmentClass::AlignmentClass: initialising Histogramm Savers:"<<endl;
	histSaverAlignment = new HistogrammSaver(0);
	histSaverAlignmentFakedTracks = new HistogrammSaver(0);
	if (verbosity)cout<<"AlignmentClass::AlignmentClass: eventReader:get Event:"<<nEventNumber<<endl;
	eventReader->GetEvent(nEventNumber);
	verbosity=3;
	if (verbosity)cout<<"AlignmentClass::AlignmentClass: eventReader:"<<eventReader->isOK()<<" "<<eventReader<<" entries:"<<eventReader->GetEntries()<<endl;
}

AlignmentClass::~AlignmentClass() {
	if(eventReader!=NULL) delete eventReader;
	if(histSaver!=NULL) delete histSaver;
	if (settings!=NULL) delete settings;
}

void AlignmentClass::SetSettings(TSettings *newSettings){
	if (verbosity)cout<<"AlignmentClass::SetSettings: getting Settings."<<flush;
	if(settings!=NULL) delete settings;
	settings = new TSettings(newSettings->getFileName());
	if (verbosity)cout<<"...DONE."<<endl;
}

void AlignmentClass::SetPlotsPath(string plotsPath){
	histSaver->SetPlotsPath(plotsPath);
	ostringstream plots_path_alignment, plots_path_alignment_CutFakeTracks;
	plots_path_alignment << plotsPath << "/alignment/";
	plots_path_alignment_CutFakeTracks << plotsPath << "/alignment_CutFakeTracks/";
	histSaverAlignment->SetPlotsPath(plots_path_alignment.str());
	histSaverAlignmentFakedTracks->SetPlotsPath(plots_path_alignment_CutFakeTracks.str());
}

bool AlignmentClass::GetEvent(UInt_t nEvent){
	return eventReader->GetEvent(nEvent);
}

void AlignmentClass::SetTracks(vector<TDiamondTrack> tracks)
{
	this->alignment_tracks=tracks;
}

void AlignmentClass::SetAlignment_tracks_mask(vector<bool> tracks_mask) {
	this->alignment_tracks_mask=tracks_mask;
}

void AlignmentClass::SetAlignment_tracks_fidcut(vector<TDiamondTrack> alignment_tracks_fidcut){
	this-> alignment_tracks_fidcut= alignment_tracks_fidcut;
}

void AlignmentClass::SetAlignment_tracks_fidcut_mask(vector<bool> alignment_tracks_fidcut_mask){
	this->alignment_tracks_fidcut_mask=alignment_tracks_fidcut_mask;
}

int AlignmentClass::Align(bool plots, bool CutFakeTracksOn){
	if (verbosity){
		cout<<"AlignentClass::Align:Starting \""<<histSaver->GetPlotsPath()<<"\""<<endl;
		cout<<"\t\t"<<alignment_tracks.size()<<" "<<alignment_tracks_mask.size()<< " ";
		cout<<		  alignment_tracks_fidcut.size()<<" "<<alignment_tracks_fidcut_mask.size()<<endl;
		cout << "\t\t "<<eventReader<<" ."<<endl;
	}
	if(alignment_tracks.size()==0) {
		cout<<"AlignmentClass::Align: No tracks found; need to CallClustering::ClusterRun first..."<<endl;
		return -1;
		//ClusterRun(plots); // doesn't use alternative clustering
	}

	if (alignment_tracks.size() == 0) {
		cout << "AlignmentClass::Align:No tracks available. Alignment not possible. (tracks.size() = " << alignment_tracks.size() << ")" << endl;
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

	//*vector<bool> alignment_tracks_mask = tracks_mask;
	//*vector<TDiamondTrack> alignment_tracks_fidcut = tracks_fidcut;
	//*vector<bool> alignment_tracks_fidcut_mask = tracks_fidcut_mask;

	if(align==NULL){
		align = new TDetectorAlignment(histSaver->GetPlotsPath(), alignment_tracks, alignment_tracks_mask);
		cout<<"AlignmentClass::Align::Detectoralignment did not exist, so created new DetectorAlignment"<<endl;
	}

	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
	for (int alignStep = 0; alignStep < nAlignSteps; alignStep++) {
		//TDetectorAlignment* align = new TDetectorAlignment(plots_path, tracks, tracks_mask);
		if (verbosity) cout<<"AlignmentClass::Align:start with alignmentStep no. "<<alignStep+1 << " of " <<nAlignSteps<<endl;
		Int_t nPasses = 10;
		Double_t plot_width_factor = 3; // scales the widths of the plots; range is a 3*width of distribution centered on mean

		align->PlotAngularDistribution(); //look at angular distribution of tracks
		align->PlotCartesianDistribution(); //look at slope distribution of tracks

		string prename = "alignment_PrealignmentResiduals";
		string postname = "alignment_PostalignmentResiduals";

		// generate residuals before alignment
		if(verbosity)cout<<"AlignmentClass::Align::CheckDetectorAlignment"<<endl;
		align->LoadTracks(this->alignment_tracks, this->alignment_tracks_mask);
		align->CheckDetectorAlignmentXYPlots(0, 1, 3, prename);
		align->CheckDetectorAlignmentXYPlots(1, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(2, 0, 3, prename);
		align->CheckDetectorAlignmentXYPlots(3, 0, 2, prename);

		// itterative alignment loop
		if(verbosity)cout<<"AlignmentClass::Align::alignmentloop"<<endl;

		align->setVerbosity(0);
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
		cout<<"AlignmentClass::Align:Checking final residuals"<<endl;
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
		cout<<"AlignmentClass::Align:Check that the telescope alignment still holds for fidcut tracks w/ single diamond cluster"<<endl;
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
		cout<<"AlignmentClass::Align:Checking final diamond residuals"<<endl;
		cout<<endl;

		// generate residuals after alignment
		align->CheckDetectorAlignmentXYPlots(4, 1, 2, postname);


		this->createAlignmentSummary();

		cout << "AlignmentClass::Align:Intrinsic silicon resolution " << align->GetSiResolution() << " strips or " << align->GetSiResolution() * 50 << "um" << endl;

		if (CutFakeTracksOn && alignStep < 1) {
			align->CutFakeTracks(alignment_tracks, alignment_tracks_mask, settings->getAlignment_chi2(), CutFakeTracksOn, true);
			align->CutFakeTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, settings->getAlignment_chi2(), CutFakeTracksOn, true);
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
//			this->TransparentClustering(tracks_transparent, alignment_tracks_fidcut_mask, align);
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
{
	cout << "AlignmentClass::TransparentAnalysis:Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
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
	cout << "diamond hit factor: " << settings->getDi_Cluster_Hit_Factor() << endl;

	cout << "AlignmentClass::TransparentAnalysis:init histograms for transparent clustering.." << endl;
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_D3X_" << (i+1) << "HitTransparClusters_8HitsFidcut";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
		//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
		//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	}
	histo_transparentclustering_landau_mean = new TH1F("PulseHeightMeanVsChannels_Det?_TranspAna_","PulseHeightMeanVsChannels_Det?_TranspAna",10,0.5,10.5);
	ostringstream histoname_eta;
	histoname_eta << "Eta_D3X_2CentroidHits_TransparClusters";
	histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
	histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_D3X_TransparClusters","DiffEstEffHit_D3X_TransparClusters", 200, -5.,5.);
	histo_transparentclustering_hitdiff_scatter = new TH2F("DiffEstEffHit_Scatter_Dia_TransparClusters","DiffEstEffHit_Scatter_Dia_TransparClusters", 200, -5.,5.,128,0,127);
	histo_transparentclustering_2Channel_PulseHeight = new TH1F("PulseHeight_D3X_2Channel_TranspCluster_8HitsFidcut","PulseHeight_D3X_2Channel_TranspCluster_8HitsFidcut",settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
	cout << " done." << endl;

	verbose = true;



	event_numbers.clear();
	for (int j = 0; j < eventReader->GetEntries(); j++) {
		eventReader->GetEvent(j);
		//		cout << endl << endl << endl;
		//		cout << "event " << j << " in PedTree has eventReader->getEvent_number(): " << eventReader->getEvent_number() << endl;
		event_numbers.push_back(eventReader->getEvent_number());
		//		for (int blablabla = 0; blablabla < 254; blablabla++) {
		//			cout << "eventReader->Dia_ADC = " << eventReader->Dia_ADC[blablabla] << ",\teventReader->Det_PedMean = " << eventReader->Det_PedMean[8][blablabla] << endl;
		//			cout << "Collected charge in channel " << blablabla << " of diamond: " << eventReader->Dia_ADC[blablabla]-eventReader->Det_PedMean[8][blablabla] << endl;
		//		}
	}
	//	return;

	// loop over tracks
	for (int i = 0; i < tracks.size(); i++) {
		if (verbosity>1) cout << " -- starting transparent clustering for track " << i << endl;

		// check if track is masked
		if (tracks[i].FakeTrack) {
			if (verbosity>2) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
			continue;
		}

		// get event number for track
		if (verbosity>2) cout << "Getting event number.. ";
		event = tracks[i].GetEventNumber();
		if (verbosity>2) cout << " -> track " << i << " corresponds to event " << event << endl;

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
		if (verbosity>2) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;

		// fit y position
		par_y.clear();
		align->LinTrackFit(z_positions, y_positions, par_y);
		if (verbosity>2) cout << "linear fit of track:\tpar_y[0] = " << par_y[0] << ",\tpar_y[1] = " << par_y[1] << endl;

		// estimate hit position in diamond
		//		diamond_z_position = align->track_holder.GetD(4).GetZ();
		diamond_hit_position = par[0] + par[1] * diamond_z_position;
		diamond_hit_y_position = par_y[0] + par_y[1] * diamond_z_position;
		diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
		diamond_hit_y_position = diamond_hit_y_position + diamond_y_offset;
		diamond_hit_position = (diamond_hit_position - 64) * TMath::Cos(diamond_phi_offset) - (diamond_hit_y_position - 64) * TMath::Sin(diamond_phi_offset) + 64; // add the tilt correction
		diamond_hit_position += 0.5; // added 0.5 to take the middle of the channel instead of the edge
		diamond_hit_channel = (int)diamond_hit_position;
		if (verbosity>2) cout << "z position of diamond is " << diamond_z_position << endl;

		// difference between estimated and effective hit in diamond
		hit_diff = TMath::Abs(eff_diamond_hit_position - diamond_hit_position);
		//		cout << "effective hit position in diamond:\t" << eff_diamond_hit_position << "\testimated position in diamond:\t" << diamond_hit_position << endl;

		//get event
		if (verbosity>2) cout << "getting event " << event << ".." << endl;
		eventReader->GetEvent(event);
		if (verbosity>2) cout << "eventReader->getEvent_number() = " << eventReader->getEvent_number() << endl;

		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int j = 0; j < 128; j++) {
			if (eventReader->getDet_ADC(6,j)-eventReader->getDet_PedMean(6,j) > (eventReader->getDet_ADC(6,eff_diamond_hit_channel)-eventReader->getDet_PedMean(6,eff_diamond_hit_channel))) {
				eff_diamond_hit_channel = j;
			}
		}
		histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
		histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position, diamond_hit_y_position);
		cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;

		// cluster diamond channels around estimated hit position
		//		cout << " --" << endl;
		if (verbosity>2) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbosity>2) cout << "evenDet_PedMeantReader->Det_ADC[6] = " << eventReader->getDet_ADC(6,diamond_hit_channel) << ",\teventReader->Det_PedMean = " << eventReader->getDet_PedMean(6,diamond_hit_channel) << endl;
		if (verbosity>2) cout << "Collected charge in channel " << (int)diamond_hit_position << " of D3X: " << eventReader->getDet_ADC(6,(int)diamond_hit_position)-eventReader->getDet_PedMean(6,(int)diamond_hit_position) << endl;

		int sign;

		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;

		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = eventReader->getDet_ADC(6,diamond_hit_channel)-eventReader->getDet_PedMean(6,diamond_hit_channel);
		secondchannel_adc = eventReader->getDet_ADC(6,diamond_secondhit_channel)-eventReader->getDet_PedMean(6,diamond_secondhit_channel);
		if (sign == 1) transp_eta = firstchannel_adc / (firstchannel_adc + secondchannel_adc);
		else transp_eta = secondchannel_adc / (firstchannel_adc + secondchannel_adc);
		histo_transparentclustering_eta->Fill(transp_eta);

		// fill pulse height histogram
		histo_transparentclustering_2Channel_PulseHeight->Fill(firstchannel_adc+secondchannel_adc);

		if (verbosity>2) cout << "AlignmentClass::TransparentAnalysis:clusters for track " << i << ":" << endl;
		// loop over different cluster sizes
		for (int j = 0; j < 10; j++) {
			cluster_adc = 0;
			current_channel = diamond_hit_channel;
			if (verbosity>2) cout << "selected channels for " << j+1 << " hit transparent cluster: ";
			current_sign = sign;
			// sum adc for n channel cluster
			for (int channel = 0; channel <= j; channel++) {
				current_channel = current_channel + current_sign * channel;
				current_sign = (-1) * current_sign;
				if (verbosity>2) cout << current_channel;
				if (verbosity>2) if (channel < j) cout << ", ";
				if (current_channel > 0 && current_channel < 128 /* && eventReader->Dia_ADC[current_channel]-eventReader->Det_PedMean[8][current_channel] > Di_Cluster_Hit_Factor*eventReader->Det_PedWidth[8][current_channel]*/)
					cluster_adc = cluster_adc + eventReader->getDet_ADC(6,current_channel)-eventReader->getDet_PedMean(6,current_channel);
			}
			if (verbosity>2) cout << endl;
			if (verbosity>2) cout << "total charge of cluster: " << cluster_adc << endl;
			if (verbosity>2) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
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

}


void AlignmentClass::TransparentClustering(vector<TDiamondTrack> & tracks, vector<bool> & tracks_mask, TDetectorAlignment *align, bool verbose)
{
	if(verbosity) cout << "AlignmentClass::TransparentClustering::Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
	cout << "AlignmentClass::TransparentClustering::get Event Reader: "<<eventReader<<" ."<<endl;
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
	bool HasMaskedChannel = 0;
	Float_t chi2X = 0;
	Float_t chi2Y = 0;
	Float_t SiRes = align->GetSiResolution();
	
	cout << "D0: x offset: " << align->GetXOffset(0) << "\ty offset: " << align->GetYOffset(0) << "\tphix offset: " << align->GetPhiXOffset(0) << "\tphiy offset: " << align->GetPhiYOffset(0) << endl;
	cout << "D1: x offset: " << align->GetXOffset(1) << "\ty offset: " << align->GetYOffset(1) << "\tphix offset: " << align->GetPhiXOffset(1) << "\tphiy offset: " << align->GetPhiYOffset(1) << endl;
	cout << "D2: x offset: " << align->GetXOffset(2) << "\ty offset: " << align->GetYOffset(2) << "\tphix offset: " << align->GetPhiXOffset(2) << "\tphiy offset: " << align->GetPhiYOffset(2) << endl;
	cout << "D3: x offset: " << align->GetXOffset(3) << "\ty offset: " << align->GetYOffset(3) << "\tphix offset: " << align->GetPhiXOffset(3) << "\tphiy offset: " << align->GetPhiYOffset(3) << endl;
	cout << "D4: x offset: " << align->GetXOffset(4) << "\ty offset: " << align->GetYOffset(4) << "\tphix offset: " << align->GetPhiXOffset(4) << "\tphiy offset: " << align->GetPhiYOffset(4) << endl;
	cout << "silicon resolution: " << SiRes << endl;
	cout << "diamond x offset: " << diamond_x_offset << endl;
	cout << "diamond y offset: " << diamond_y_offset << endl;
	cout << "diamond phi offset: " << diamond_phi_offset << endl;
	cout << "diamond phi y offset: " << diamond_phi_y_offset << endl;
	cout << "diamond hit factor: " << settings->getDi_Cluster_Hit_Factor() << endl;

	if (verbosity) cout << "AlignmentClass::TransparentClustering::init histograms for transparent clustering.." << endl;
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_Dia_" << (i+1) << "ChannelsTransparAna_8HitsFidcut";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;


		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
		//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
		//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	}

	histo_transparentclustering_landau_mean = new TH1F("PulseHeightMeanVsChannels_Dia_TranspAna","PulseHeightMeanVsChannels_Dia_TranspAna",11,-0.5,10.5);
	histo_transparentclustering_SNR_vs_channel = new TH1F("SignificanceVsChannels_Dia_TranspAna","SignificanceVsChannels_Dia_TranspAna",11,-0.5,10.5);
	ostringstream histoname_eta;
	histoname_eta << "Eta_Dia_2CentroidHits_TransparClusters";
	histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
	histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_Dia_TransparClusters","DiffEstEffHit_Dia_TransparClusters", 200, -5.,5.);
	histo_transparentclustering_hitdiff_scatter = new TH2F("DiffEstEffHit_Scatter_Dia_TransparClusters","DiffEstEffHit_Scatter_Dia_TransparClusters", 200, -5.,5.,256,0,255);
    histo_transparentclustering_2Channel_PulseHeight = new TH1F("PulseHeight_Dia_2LargestHitsIn10Strips_TranspCluster_8HitsFidcut","PulseHeight_Dia_2LargestHitsIn10Strips_TranspCluster_8HitsFidcut",settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_residuals, histoname_residuals_scatter, histoname_residuals_largest_hit, histoname_residuals_largest_hit_scatter;
		histoname_residuals << "TranspAnaResidualsDia" << (i+1) << "StripsChargeWeighted";
		histoname_residuals_scatter << "TranspAnaResidualsDia" << (i+1) << "StripsChargeWeightedVsY";
		histoname_residuals_largest_hit << "TranspAnaResidualsDiaLargestHitIn" << (i+1) << "Strips";
		histoname_residuals_largest_hit_scatter << "TranspAnaResidualsDiaLargestHitIn" << (i+1) << "StripsVsY";
		histo_transparentclustering_residuals[i] = new TH1F(histoname_residuals.str().c_str(),histoname_residuals.str().c_str(), 200, -2.5, 2.5);
		histo_transparentclustering_residuals_scatter[i] = new TH2F(histoname_residuals_scatter.str().c_str(),histoname_residuals_scatter.str().c_str(), 200, -2.5, 2.5,256,0,255);
		histo_transparentclustering_residuals_largest_hit[i] = new TH1F(histoname_residuals_largest_hit.str().c_str(),histoname_residuals_largest_hit.str().c_str(), 200, -2.5, 2.5);
		histo_transparentclustering_residuals_largest_hit_scatter[i] = new TH2F(histoname_residuals_largest_hit_scatter.str().c_str(),histoname_residuals_largest_hit_scatter.str().c_str(), 200, -2.5, 2.5,256,0,255);
	}
	histo_transparentclustering_residuals_2largest_hits = new TH1F("TranspAnaResidualsDia2LargestHitsIn10StripsChargedWeighted","TranspAnaResidualsDia2LargestHitsIn10StripsChargedWeighted",200, -2.5, 2.5);
	histo_transparentclustering_residuals_2largest_hits_scatter = new TH2F("TranspAnaResidualsDia2LargestHitsIn10StripsChargedWeightedVsY","TranspAnaResidualsDia2LargestHitsIn10StripsChargedWeightedVsY",200, -2.5, 2.5, 256, 0, 255);
	histo_transparentclustering_chi2X = new TH1F("TranspAnaChi2X","TranspAnaChi2X",200,0,20);
	histo_transparentclustering_chi2Y = new TH1F("TranspAnaChi2Y","TranspAnaChi2Y",200,0,20);
	cout << " done." << endl;


	//Telescope Data Branches
	event_numbers.clear();
	eventReader=new TADCEventReader(PedFileName);
	if(verbosity) {
		cout << "AlignmentClass::TransparentClustering::get Event numbers: "<<eventReader<<" ."<<flush;
		cout<<eventReader->isOK()<<flush;
		cout <<". "<< eventReader->GetEntries()<<endl;
	}
	/*creating list of eventnumbers: event_numbers*/
	if(verbosity>=3)cout<<"create list of eventnumbers...";
	for (int j = 0; j < eventReader->GetEntries(); j++) {

		if(verbosity>=4) cout << "AlignmentClass::TransparentClustering::get Event"<<j<<flush;
		if(!eventReader->GetEvent(j)) continue;
		//		cout << endl << endl << endl;
		//		cout << "event " << j << " in PedTree has eventReader->getEvent_number(): " << eventReader->getEvent_number() << endl;

		if(verbosity>=4) cout << " push back "<<flush;
		event_numbers.push_back(eventReader->getEvent_number());
		if(verbosity>=4) cout << " done. "<<endl;
		//		for (int blablabla = 0; blablabla < 254; blablabla++) {
		//			cout << "eventReader->Dia_ADC = " << eventReader->Dia_ADC[blablabla] << ",\teventReader->Det_PedMean = " << eventReader->Det_PedMean[8][blablabla] << endl;
		//			cout << "Collected charge in channel " << blablabla << " of diamond: " << eventReader->Dia_ADC[blablabla]-eventReader->Det_PedMean[8][blablabla] << endl;
		//		}
	}
	if(verbosity>=3)cout<<"DONE"<<endl;
	//	return;

	// loop over tracks
	if(verbosity) cout << "AlignmentClass::TransparentClustering::loop over tracks"<<endl;
	for (int i = 0; i < tracks.size(); i++) {
		if (verbosity>=3) cout << " -- starting transparent clustering for track " << i << endl;

		// check if track is masked
		if (tracks[i].FakeTrack) {
			if (verbosity>=3) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
			continue;
		}

		// get event number for track
		if (verbosity>=3) cout << "Getting event number.. ";
		event = tracks[i].GetEventNumber();
		if (verbosity>=3) cout << " -> track " << i << " corresponds to event " << event << endl;

		// check if event number is valid
		if (event < 0) {
			if(verbosity>=3) cout << "Track " << i << " has no event number. Skipping this track.." << endl;
			continue;
		}

		// load data (apply offset)
		align->LoadData(tracks[i]);

		// read out x, y, z positions
		x_positions.clear();
		y_positions.clear();
		z_positions.clear();
		if (verbosity>=3) cout << "Det\tx\ty\tz"<< endl;
		for (int det = 0; det < 4; det++) {
			x_positions.push_back(align->track_holder.GetD(det).GetX());
			y_positions.push_back(align->track_holder.GetD(det).GetY());
			z_positions.push_back(align->track_holder.GetD(det).GetZ());
			if (verbosity>=3) cout << det<<"\t" << align->track_holder.GetD(det).GetX()<<"\t" << align->track_holder.GetD(det).GetY()<<"\t" << align->track_holder.GetD(det).GetZ() << endl;
		}

		// read out effictive diamond hit position
		eff_diamond_hit_position = align->track_holder.GetD(4).GetZ();
		//if (verbosity>=3) cout << det<<"\t" << align->track_holder.GetD(det).GetX()<<"\t" << align->track_holder.GetD(det).GetY()<<"\t" << align->track_holder.GetD(det).GetZ() << endl;.GetX();

		// fit track
		par.clear();
		chi2X = 0;
		chi2X = align->LinTrackFit(z_positions, x_positions, par, SiRes);
		if (verbosity>=3) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;
		histo_transparentclustering_chi2X->Fill(chi2X);

		// fit y position
		par_y.clear();
		chi2Y = 0;
        chi2Y = align->LinTrackFit(z_positions, y_positions, par_y, SiRes);
        if (verbosity>=3) cout << "linear fit of track:\tpar_y[0] = " << par_y[0] << ",\tpar_y[1] = " << par_y[1] << endl;
		histo_transparentclustering_chi2Y->Fill(chi2Y);
		
//		if (chi2X > 10 || chi2Y > 10) {
//			continue;
//		}

		// estimate hit position in diamond
		//		diamond_z_position = align->track_holder.GetD(4).GetZ();
		diamond_hit_position = par[0] + par[1] * diamond_z_position;
		diamond_hit_y_position = par_y[0] + par_y[1] * diamond_z_position;
		diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
		diamond_hit_y_position = diamond_hit_y_position + diamond_y_offset;
		diamond_hit_position = (diamond_hit_position - 64) * TMath::Cos(diamond_phi_offset) - (diamond_hit_y_position - 64) * TMath::Sin(diamond_phi_offset) + 64; // add the tilt correction
		diamond_hit_position += 0.5; // added 0.5 to take the middle of the channel instead of the edge
		diamond_hit_channel = (int)diamond_hit_position;
		if (verbosity>=3) cout << "z position of diamond is " << diamond_z_position << endl;
		
		if (diamond_hit_channel < 7 || diamond_hit_channel > 120) HasMaskedChannel = true;
		else {
			for (int j = diamond_hit_channel-5; j < diamond_hit_channel+6; j++) {
				if (!settings->getDet_channel_screen(8).CheckChannel((int)eventReader->getDet_Channels(8,j))) HasMaskedChannel = true;
			}
		}
		if (false) {//(HasMaskedChannel) {
			cout << "estimated hit position in diamond ist close to masked channel --> skipping this track.." << endl;
			continue;
		}

		// difference between estimated and effective hit in diamond
		hit_diff = TMath::Abs(eff_diamond_hit_position - diamond_hit_position);
		//		cout << "effective hit position in diamond:\t" << eff_diamond_hit_position << "\testimated position in diamond:\t" << diamond_hit_position << endl;

		//get event
		if (verbosity>=3) cout << "getting event " << event << ".." << endl;
		eventReader->GetEvent(event);
		if (verbosity>=3) cout << "eventReader->getEvent_number() = " << eventReader->getEvent_number() << endl;

		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int j = 0; j < 128; j++) {
			if (eventReader->getDia_ADC(j)-eventReader->getDet_PedMean(8,j) > (eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel))) {
				eff_diamond_hit_channel = j;
			}
		}
		
		// charge weighted effictive hit position (2 strips)
		if (false) {
			if (eventReader->getDia_ADC(eff_diamond_hit_channel-1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel-1) < eventReader->getDia_ADC(eff_diamond_hit_channel+1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel+1)) {
				Float_t q1, q2;
				q1 = eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel);
				q2 = eventReader->getDia_ADC(eff_diamond_hit_channel+1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel+1);
				eff_diamond_hit_position = eff_diamond_hit_channel + 0.5;//q2 / (q1+q2);
			}
			else {
				Float_t q1, q2;
				q1 = eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel);
				q2 = eventReader->getDia_ADC(eff_diamond_hit_channel-1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel-1);
				eff_diamond_hit_position = eff_diamond_hit_channel - 0.5;//q2 / (q1+q2);
			}
		}
		
		// charge weighted effictive hit position (3 strips)
		if (true) {
			Float_t q1,q2,q3;
			q1 = eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel);
			q2 = eventReader->getDia_ADC(eff_diamond_hit_channel-1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel-1);
			q3 = eventReader->getDia_ADC(eff_diamond_hit_channel+1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel+1);
			eff_diamond_hit_position = eff_diamond_hit_channel - q2/(q1+q2+q3) + q3/(q1+q2+q3);
		}
		
		
		//		eff_diamond_hit_position =- diamond_x_offset;
		
		histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_position + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
        histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_position + 0.5 - diamond_hit_position, diamond_hit_y_position);
		if (verbosity>=3) cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;

		// cluster diamond channels around estimated hit position
		//		cout << " --" << endl;
		if (verbosity>=3) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbosity>=3) cout << "eventReader->Dia_ADC = " << eventReader->getDia_ADC(diamond_hit_channel) << ",\teventReader->Det_PedMean = " << eventReader->getDet_PedMean(8,diamond_hit_channel) << endl;
		if (verbosity>=3) cout << "Collected charge in channel " << (int)diamond_hit_position << " of diamond: " << eventReader->getDia_ADC((int)diamond_hit_position)-eventReader->getDet_PedMean(8,(int)diamond_hit_position) << endl;

		int sign;

		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;

		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = eventReader->getDia_ADC(diamond_hit_channel)-eventReader->getDet_PedMean(8,diamond_hit_channel);
		secondchannel_adc = eventReader->getDia_ADC(diamond_secondhit_channel)-eventReader->getDet_PedMean(8,diamond_secondhit_channel);
		if (sign == 1) transp_eta = firstchannel_adc / (firstchannel_adc + secondchannel_adc);
		else transp_eta = secondchannel_adc / (firstchannel_adc + secondchannel_adc);
		histo_transparentclustering_eta->Fill(transp_eta);

		// fill pulse height histogram
//		histo_transparentclustering_2Channel_PulseHeight->Fill(firstchannel_adc+secondchannel_adc);

		if (verbosity>=3) cout << "clusters for track " << i << ":" << endl;
		Float_t charge_mean, charge_weighted_position, charge_weighted_2largest_position;
		int dia_largest_hit, dia_second_largest_hit;
		// loop over different cluster sizes
		for (int j = 0; j < 10; j++) {
			cluster_adc = 0;
			current_channel = diamond_hit_channel;
			if (verbosity>=3) cout << "selected channels for " << j+1 << " hit transparent cluster: ";
			current_sign = sign;
//			Float_t qtot = 0;
			charge_mean = 0;
			charge_weighted_position = 0;
			dia_largest_hit = 0;
			// sum adc for n channel cluster
			for (int channel = 0; channel <= j; channel++) {
				current_channel = current_channel + current_sign * channel;
				current_sign = (-1) * current_sign;
				if (verbosity>=3) cout << current_channel;
				if (verbosity>=3) if (channel < j) cout << ", ";
				if (current_channel > 0 && current_channel < 128 /* && eventReader->Dia_ADC[current_channel]-eventReader->Det_PedMean[8][current_channel] > Di_Cluster_Hit_Factor*eventReader->Det_PedWidth[8][current_channel]*/) {
					cluster_adc = cluster_adc + eventReader->getDia_ADC(current_channel)-eventReader->getDet_PedMean(8,current_channel);
					//					qtot += Dia_ADC[current_channel]-Det_PedMean[8][current_channel];
					charge_mean += (eventReader->getDia_ADC(current_channel)-eventReader->getDet_PedMean(8,current_channel)) * (current_channel+0.5);
					//					cout << "Dia_ADC["<<current_channel<<"]-Det_PedMean[8]["<<current_channel<<"] = " << Dia_ADC[current_channel]-Det_PedMean[8][current_channel];
					//					cout << "\tcharge_mean = " << charge_mean << endl;
					if (dia_largest_hit == 0 || (eventReader->getDia_ADC(current_channel)-eventReader->getDet_PedMean(8,current_channel)) > (eventReader->getDia_ADC(dia_largest_hit)-eventReader->getDet_PedMean(8,dia_largest_hit))) {
						dia_second_largest_hit = dia_largest_hit;
						dia_largest_hit = current_channel;
					}
					else {
						if (channel > 0) {
							if (dia_second_largest_hit == 0 || (eventReader->getDia_ADC(current_channel)-eventReader->getDet_PedMean(8,current_channel)) > (eventReader->getDia_ADC(dia_second_largest_hit)-eventReader->getDet_PedMean(8,dia_second_largest_hit))) {
								dia_second_largest_hit = current_channel;
							}
						}
					}
				}
			}
			charge_weighted_position = charge_mean / cluster_adc;
			if (j == 9) {
				Float_t q1,q2;
				q1 = eventReader->getDia_ADC(dia_largest_hit)-eventReader->getDet_PedMean(8,dia_largest_hit);
				q2 = eventReader->getDia_ADC(dia_second_largest_hit)-eventReader->getDet_PedMean(8,dia_second_largest_hit);
				charge_weighted_2largest_position = 1/(q1+q2) * (q1*(Float_t)dia_largest_hit + q2*(Float_t)dia_second_largest_hit);
				histo_transparentclustering_residuals_2largest_hits->Fill(charge_weighted_2largest_position + 0.5 - diamond_hit_position);
				histo_transparentclustering_residuals_2largest_hits_scatter->Fill(charge_weighted_2largest_position + 0.5 - diamond_hit_position, diamond_hit_y_position);
				histo_transparentclustering_2Channel_PulseHeight->Fill(q1+q2);
			}
			histo_transparentclustering_residuals[j]->Fill(charge_weighted_position - diamond_hit_position);
			histo_transparentclustering_residuals_scatter[j]->Fill(charge_weighted_position - diamond_hit_position, diamond_hit_y_position);
			histo_transparentclustering_residuals_largest_hit[j]->Fill(dia_largest_hit + 0.5 - diamond_hit_position);
			histo_transparentclustering_residuals_largest_hit_scatter[j]->Fill(dia_largest_hit + 0.5 - diamond_hit_position, diamond_hit_y_position);
			if (verbose) cout << endl;
			if (verbose) cout << "total charge of " << j+1 << " strips: " << cluster_adc << "\tcharge_weighted_position: " << charge_weighted_position << endl;
			if (verbose) cout << "total charge of cluster: " << cluster_adc << endl;
			if (verbose) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
			if (current_channel <= 0 || current_channel >= 128) break;
			histo_transparentclustering_landau[j]->Fill(cluster_adc);
		} // end loop over cluster sizes
	} // end loop over tracks
	
	Float_t noise_multiple_channels = dianoise_sigma[1];

	// save histograms
	for (int i = 0; i < 10; i++) {
        histo_transparentclustering_landau_mean->SetBinContent(i+2,histo_transparentclustering_landau[i]->GetMean()); // plot pulse hight means into a histogram
		histo_transparentclustering_SNR_vs_channel->SetBinContent(i+2,histo_transparentclustering_landau[i]->GetMean()/(TMath::Sqrt(i+1)*dianoise_sigma[1]));
		//		noise_multiple_channels = noise_multiple_channels * TMath::Sqrt(2);
		histSaver->SaveHistogram(histo_transparentclustering_landau[i]);
		if (i>0) histSaver->SaveHistogram(histo_transparentclustering_residuals[i],1);
		else histSaver->SaveHistogram(histo_transparentclustering_residuals[i],0);
		histSaver->SaveHistogram(histo_transparentclustering_residuals_scatter[i]);
		histSaver->SaveHistogram(histo_transparentclustering_residuals_largest_hit[i],0);
		histSaver->SaveHistogram(histo_transparentclustering_residuals_largest_hit_scatter[i]);
	}
	histSaver->SaveHistogram(histo_transparentclustering_residuals_2largest_hits,1);
	histSaver->SaveHistogram(histo_transparentclustering_residuals_2largest_hits_scatter);
	histo_transparentclustering_landau_mean->Scale(1/histo_transparentclustering_landau[9]->GetMean());
	
	histSaver->SaveHistogram(histo_transparentclustering_chi2X);
	histSaver->SaveHistogram(histo_transparentclustering_chi2Y);
	
	histSaver->SaveHistogram(histo_transparentclustering_hitdiff,1);
	
	
    histSaver->SaveHistogram(histo_transparentclustering_landau_mean);
	histSaver->SaveHistogram(histo_transparentclustering_SNR_vs_channel);
	histSaver->SaveHistogram(histo_transparentclustering_eta);
	//	histSaver->SaveHistogram(histo_transparentclustering_hitdiff);
    histSaver->SaveHistogram(histo_transparentclustering_hitdiff_scatter);
    histSaver->SaveHistogram(histo_transparentclustering_2Channel_PulseHeight);
}

void AlignmentClass::createAlignmentSummary(){
	//report results in a file

			ostringstream alignment_summary_path;
			alignment_summary_path << histSaver->GetPlotsPath() << "Alignment_Summary.txt";
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
}




