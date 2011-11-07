/*
 * TTransparentClustering.cpp
 *
 *  Created on: 04.11.2011
 *      Author: bachmair
 */

#include "../include/TTransparentClustering.hh"

TTransparentClustering::TTransparentClustering(string PedFileName) {
	// TODO Auto-generated constructor stub
	diamond_z_position = 10.2;
	verbosity=0;
	settings=NULL;
	histSaver=NULL;
	eventReader=NULL;
	currentAlignment=NULL;
	dianoise_sigma[0]=1;
	dianoise_sigma[1]=1;//TODO: Figure out what are correct values... how do they have to beset???

	this->eventNumberDifference=1500;


	histSaver= new HistogrammSaver(0);
	eventReader = new TADCEventReader(PedFileName);
	for (int i = 0; i < 10; i++) {
		histo_transparentclustering_landau[i] = NULL;
	}
}

TTransparentClustering::~TTransparentClustering() {
	// TODO Auto-generated destructor stub
}

void TTransparentClustering::MakeTransparentClustering()
{

	cout<<"\n\n\n\n\n*********************************************************************************"<<endl;
	cout<<"*********************************************************************************" <<endl;
	cout<<"****** TTransparentClustering::MakeTransparentClustering ************************"<<endl;
	cout<<"*********************************************************************************"<<endl;
	cout<<"*********************************************************************************\n\n"<<endl;

	int event;
	diamond_hit_position = 0;
	diamond_hit_y_position = 0;
	diamond_hit_channel = 0;
	diamond_secondhit_channel = 0;
	int current_channel, current_sign;
	Float_t cluster_adc = 0;
	Float_t transp_eta = 0;
	Float_t firstchannel_adc, secondchannel_adc;
	eff_diamond_hit_position = 0;
	eff_diamond_hit_channel = 0;
	bool HasMaskedChannel = 0;
	chi2X = 0;
	chi2Y = 0;


	cout << "D0: x offset: " << currentAlignment->GetXOffset(0) << "\ty offset: " << currentAlignment->GetYOffset(0) << "\tphix offset: " << currentAlignment->GetPhiXOffset(0) << "\tphiy offset: " << currentAlignment->GetPhiYOffset(0) << endl;
	cout << "D1: x offset: " << currentAlignment->GetXOffset(1) << "\ty offset: " << currentAlignment->GetYOffset(1) << "\tphix offset: " << currentAlignment->GetPhiXOffset(1) << "\tphiy offset: " << currentAlignment->GetPhiYOffset(1) << endl;
	cout << "D2: x offset: " << currentAlignment->GetXOffset(2) << "\ty offset: " << currentAlignment->GetYOffset(2) << "\tphix offset: " << currentAlignment->GetPhiXOffset(2) << "\tphiy offset: " << currentAlignment->GetPhiYOffset(2) << endl;
	cout << "D3: x offset: " << currentAlignment->GetXOffset(3) << "\ty offset: " << currentAlignment->GetYOffset(3) << "\tphix offset: " << currentAlignment->GetPhiXOffset(3) << "\tphiy offset: " << currentAlignment->GetPhiYOffset(3) << endl;
	cout << "D4: x offset: " << currentAlignment->GetXOffset(4) << "\ty offset: " << currentAlignment->GetYOffset(4) << "\tphix offset: " << currentAlignment->GetPhiXOffset(4) << "\tphiy offset: " << currentAlignment->GetPhiYOffset(4) << endl;
	cout << "silicon resolution: " << SiRes << endl;
	cout << "diamond x offset: " << diamond_x_offset << endl;
	cout << "diamond y offset: " << diamond_y_offset << endl;
	cout << "diamond phi offset: " << diamond_phi_offset << endl;
	cout << "diamond phi y offset: " << diamond_phi_y_offset << endl;
	//cout << "diamond hit factor: " << settings->getDi_Cluster_Hit_Factor() << endl;

	initHistogramms();

	createEventNumberList();


	//Telescope Data Branches



	// loop over tracks
	if(verbosity) cout << "AlignmentClass::TransparentClustering::loop over tracks"<<endl;
	for (int i = 0; i < tracks.size(); i++) {
		if (verbosity>=4) cout << " -- starting transparent clustering for track " << i << endl;

		// check if track is masked
		if (tracks[i].FakeTrack) {
			if (verbosity>=2) cout << "Clustering::TransparentClustering: Track " << i << " is masked as fake track and skipped." << endl;
			continue;
		}

		// get event number for track
		if (verbosity>=3) cout << "Getting event number.. ";
		event = tracks[i].GetEventNumber();
		if (verbosity>=3) cout << " -> track " << i << " corresponds to event " << event << endl;

		// check if event number is valid
		if (event < 0) {
			if(verbosity>=2) cout << "Track " << i << " has no event number. Skipping this track.." << endl;
			continue;
		}

		// load data (apply offset)
		this->currentAlignment->LoadData(tracks[i]);

		LoadXYZPositions();

		GetEffectiveDiamondPosition();

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

		//get event
		if (verbosity>=3) cout << "getting event " << event << ".." << endl;
		if (event>eventReader->GetEntries()){
			cout<<"event is bigger than number of events in tree: take next event..."<<endl;
			continue;
		}
		eventReader->GetEvent(event);
		if (verbosity>=3) cout << "eventReader->getEvent_number() = " << eventReader->getEvent_number() << endl;



		CalculateEffektiveHitPosition();
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

		if (verbosity>=3) cout << "clusters for track " << i << ":" <<endl;
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
						}//if channel >0
					}//else (dia_largest_hit == 0 ||
				}//IF 128>current_channel>0&&c
				else {
					cout<<"current_channel not between 0 and 128"<<endl;
				}
			}//end for channel
			if(cluster_adc!=0){
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
				if (verbosity) cout << endl;
				if (verbosity) cout << "total charge of " << j+1 << " strips: " << cluster_adc << "\tcharge_weighted_position: " << charge_weighted_position << endl;
				if (verbosity) cout << "total charge of cluster: " << cluster_adc << endl;
				if (verbosity) cout << "histo_transparentclustering_landau[" << j << "] address: " << histo_transparentclustering_landau[j] << endl;
				if (current_channel <= 0 || current_channel >= 128) break;
				histo_transparentclustering_landau[j]->Fill(cluster_adc);
			}
			else{
				cout<<"cluster_adc == 0 ==> BREAK"<<endl;
				continue;
			}

		} // end loop over cluster sizes
	} // end loop over tracks

	Float_t noise_multiple_channels = dianoise_sigma[1];

}




Double_t TTransparentClustering::getDiamondPhiOffset() const
{
    return diamond_phi_offset;
}

Double_t TTransparentClustering::getDiamondPhiYOffset() const
{
    return diamond_phi_y_offset;
}

Double_t TTransparentClustering::getDiamondXOffset() const
{
    return diamond_x_offset;
}

Double_t TTransparentClustering::getDiamondYOffset() const
{
    return diamond_y_offset;
}

Double_t TTransparentClustering::getDiamondZPosition() const
{
    return diamond_z_position;
}

void TTransparentClustering::setDiamondPhiOffset(Double_t diamondPhiOffset)
{
	if(verbosity) cout<<" TTransparentClustering::setDiamondPhiOffset::"<<diamondPhiOffset<<endl;
    diamond_phi_offset = diamondPhiOffset;
}

void TTransparentClustering::setDiamondPhiYOffset(Double_t diamondPhiYOffset)
{
	if(verbosity) cout<<" TTransparentClustering::setDiamondPhiYOffset::"<<diamondPhiYOffset<<endl;
    diamond_phi_y_offset = diamondPhiYOffset;
}

void TTransparentClustering::setDiamondXOffset(Double_t diamondXOffset)
{
	if(verbosity) cout<<" TTransparentClustering::setDiamondXOffset::"<<diamondXOffset<<endl;
    diamond_x_offset = diamondXOffset;
}

void TTransparentClustering::setDiamondYOffset(Double_t diamondYOffset)
{
	if(verbosity) cout<<" TTransparentClustering::setDiamondYOffset::"<<diamondYOffset<<endl;
    diamond_y_offset = diamondYOffset;
}

void TTransparentClustering::setDiamondZPosition(Double_t diamondZPosition)
{
	if(verbosity) cout<<" TTransparentClustering::diamondZPosition::"<<diamondZPosition<<endl;
    diamond_z_position = diamondZPosition;
}

Float_t TTransparentClustering::getSiRes() const
{
    return SiRes;
}

void TTransparentClustering::setSiRes(Float_t siRes)
{
	if(verbosity)
		cout<<"TTransparentClustering::setSiRes::"<<SiRes<<endl;
    SiRes = siRes;
}

TDetectorAlignment *TTransparentClustering::getCurrentAlignment() const
{
    return currentAlignment;
}

void TTransparentClustering::setCurrentAlignment(TDetectorAlignment *currentAlignment)
{
    this->currentAlignment = currentAlignment;
}

void TTransparentClustering::initHistogramms()
{
	if (verbosity) cout << "AlignmentClass::TransparentClustering::init histograms for transparent clustering.." << endl;

	if(settings==NULL) {
		cout<<"Settings not yet set... Please set first"<<endl;
		exit(-1);
	}

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
	if (verbosity) cout << " done." << endl;
}







void TTransparentClustering::SetSettings(TSettings *newSettings){
	if (verbosity)cout<<"TTransparentClustering::SetSettings: getting Settings."<<flush;
	if(settings!=NULL) delete settings;
	settings = new TSettings(newSettings->getFileName());
	if (verbosity)cout<<"...DONE."<<endl;
}

void TTransparentClustering::SetPlotsPath(string plotsPath){
	histSaver->SetPlotsPath(plotsPath);
}

vector<TDiamondTrack> TTransparentClustering::getTracks() const
{
    return tracks;
}

void TTransparentClustering::setTracks(vector<TDiamondTrack> tracks){
	if(verbosity)cout<< "TTransparentClustering::setTracks"<<endl;
//	this->tracks.clear();
//	for(int i=0;i<tracks.size();i++)
//		this->tracks.push_back(tracks.at(i));
	this->tracks=tracks;
}

void TTransparentClustering::createEventNumberList()
{
	/*creating list of eventnumbers: event_numbers*/
	if(verbosity) {
		cout << "AlignmentClass::TransparentClustering::get Event numbers: "<<eventReader<<" ."<<flush;
		cout<<eventReader->isOK()<<flush;
		cout <<". "<< eventReader->GetEntries()<<endl;
	}
	event_numbers.clear();
	for (int j = 0; j < eventReader->GetEntries(); j++) {

		if(verbosity>=4) cout << "AlignmentClass::TransparentClustering::get Event"<<j<<flush;
		if(!eventReader->GetEvent(j)) continue;
		if(verbosity>=4) cout << " push back "<<flush;
		event_numbers.push_back(eventReader->getEvent_number());
		if(verbosity>=4) cout << " done. "<<endl;
	}
	if(verbosity>=3)cout<<"DONE"<<endl;

}

void TTransparentClustering::LoadXYZPositions()
{

	// read out x, y, z positions
	x_positions.clear();
	y_positions.clear();
	z_positions.clear();
	if (verbosity>=3) cout << "Det\tx\ty\tz"<< endl;
	for (int det = 0; det < 4; det++) {
		x_positions.push_back(currentAlignment->track_holder.GetD(det).GetX());
		y_positions.push_back(currentAlignment->track_holder.GetD(det).GetY());
		z_positions.push_back(currentAlignment->track_holder.GetD(det).GetZ());
		if (verbosity>=3) cout << det<<"\t" << currentAlignment->track_holder.GetD(det).GetX()<<"\t" << currentAlignment->track_holder.GetD(det).GetY()<<"\t" << currentAlignment->track_holder.GetD(det).GetZ() << endl;
	}

}

Int_t TTransparentClustering::getVerbosity() const
{
    return verbosity;
}

void TTransparentClustering::setVerbosity(Int_t verbosity)
{
    this->verbosity = verbosity;
}

void TTransparentClustering::GetEffectiveDiamondPosition()
{

	// read out effictive diamond hit position
	eff_diamond_hit_position = currentAlignment->track_holder.GetD(4).GetZ();
	//if (verbosity>=3) cout << det<<"\t" << currentAlignment->track_holder.GetD(det).GetX()<<"\t" << currentAlignment->track_holder.GetD(det).GetY()<<"\t" << currentAlignment->track_holder.GetD(det).GetZ() << endl;.GetX();

	// fit track
	par.clear();
	chi2X = 0;
	chi2X = currentAlignment->LinTrackFit(z_positions, x_positions, par, SiRes);
	if (verbosity>=3) cout << "linear fit of track:\tpar[0] = " << par[0] << ",\tpar[1] = " << par[1] << endl;
	histo_transparentclustering_chi2X->Fill(chi2X);

	// fit y position
	par_y.clear();
	chi2Y = 0;
    chi2Y = currentAlignment->LinTrackFit(z_positions, y_positions, par_y, SiRes);
    if (verbosity>=3) cout << "linear fit of track:\tpar_y[0] = " << par_y[0] << ",\tpar_y[1] = " << par_y[1] << endl;
	histo_transparentclustering_chi2Y->Fill(chi2Y);

//		if (chi2X > 10 || chi2Y > 10) {
//			continue;
//		}

	// estimate hit position in diamond
	//		diamond_z_position = currentAlignment->track_holder.GetD(4).GetZ();
	diamond_hit_position = par[0] + par[1] * diamond_z_position;
	diamond_hit_y_position = par_y[0] + par_y[1] * diamond_z_position;
	diamond_hit_position = diamond_hit_position + diamond_x_offset; // add offset
	diamond_hit_y_position = diamond_hit_y_position + diamond_y_offset;
	diamond_hit_position = (diamond_hit_position - 64) * TMath::Cos(diamond_phi_offset) - (diamond_hit_y_position - 64) * TMath::Sin(diamond_phi_offset) + 64; // add the tilt correction
	diamond_hit_position += 0.5; // added 0.5 to take the middle of the channel instead of the edge
	diamond_hit_channel = (int)diamond_hit_position;
	if (verbosity>=3) cout << "z position of diamond is " << diamond_z_position << endl;


	histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_position + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
    histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_position + 0.5 - diamond_hit_position, diamond_hit_y_position);
	if (verbosity>=3) cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;

}

void TTransparentClustering::SaveHistogramms()
{

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

void TTransparentClustering::CalculateEffektiveHitPosition()
{
	// find biggest hit in diamond
	eff_diamond_hit_channel = 0;
	for (int j = 0; j < 128; j++) {
		if (eventReader->getDia_ADC(j)-eventReader->getDet_PedMean(8,j) > (eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel))) {
			eff_diamond_hit_channel = j;
		}
	}
	// charge weighted effictive hit position (3 strips)ENABLED
	if (true) {
		Float_t q1,q2,q3;
		q1 = eventReader->getDia_ADC(eff_diamond_hit_channel)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel);
		q2 = eventReader->getDia_ADC(eff_diamond_hit_channel-1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel-1);
		q3 = eventReader->getDia_ADC(eff_diamond_hit_channel+1)-eventReader->getDet_PedMean(8,eff_diamond_hit_channel+1);
		eff_diamond_hit_position = eff_diamond_hit_channel - q2/(q1+q2+q3) + q3/(q1+q2+q3);
	}

	// charge weighted effictive hit position (2 strips)DISABLED
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




}






