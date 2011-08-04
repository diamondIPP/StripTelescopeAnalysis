//machinery for clustering
//2010-07-11 Goal: would like to multithread by splitting a run into several windows, then on each window calc the ped, cluster, track, and at each stage combine windows to fill histograms. Can have histograms talk to each other by making trees of data to histogram or writing histograms.root and appending there. 
//2010-07-23 Import settings from Setting.ini
//2010-07-27 Clustering analysis done for the most part
//           Decided to save clusters as arrays of basic datatypes (can't save classes to trees)
//           Putting off data reduction to PedestalAnalyze and ClusterAnalyze 
//2010-08-01 Started major histogramming efforts
//2010-08-02 deleting histograms CRASHES root on second call to the constructor so no more deleting
//           draws histograms and prints cutflow stats
//2010-08-04 Fixed settings array parsing bug (first element in array always was zero)
//           Better channel screening
//2010-08-05 Added run stats to histograms
//2010-08-14 Skip analysis for CMN and ZeroDivisor events (flagged in Sliding Pedestal)
//2010-08-15 Created HTML summary page
//           Fixed a slew of problems with the HTML summary page
//           Fixed some cuts (now skip bad cluster events in silicon and clusters w/ bad seeds in diamond)
//           Various other fixes
//           TODO: figure out why D1 has no coincident clusters for runs 14109 and 14110 while it works for 14107
//2010-08-19 Created LumpyCluster flag (if cluster is not monotonically decreasing from the peak)
//           TODO: create plots for cut clusters
//           TODO: create plots for high/low etas
//           TODO: write cut cluster event numbers to csv files
//2010-09-02 Eta vs Q plots and saturated cluster removal and sat clus hit occ plots
//2010-09-18 TODO: Plan for tracking is as follows
//                    -> align adhoc without tracking
//                    -> track with residuals' widths as uncertainties
//                    -> align with tracking (this makes sure tracks are really tracks; thereby improving residuals)
//2010-11-12 Blinded the alignment resolution determination
//2010-12-01 Created AutoFidCut() to Clustering.class.cpp
//2011-03-07 XCode 4.0 Test max AND second Xcode 4 Test!
//2011-03-28 Xcode 4 check by Lukas!
//2011-04-26 Check by Lukas

//Class and stuct definitions
//#include "ChannelScreen.hh" //Channel Screen Class
#include "Clustering.class.hh"

Clustering::Clustering(unsigned int RunNumber, string RunDescription) {
	verbosity=1;
	if(verbosity)cout<<"Clustering::Clustering"<<endl;
	eventReader=NULL;
	settings =NULL;
	alignment=NULL;

	//Alignment

	//Telescope geometry (looks like this is the same geom for all runs)
	//Note: if alignment gives large residuals, try switching wide to compact geometry or vs

	//wide geometry; edges: 0, 2.40, 9, 18, 20.40 (si modules 2.40cm wide, x/y planes spaced 2mm, D0/D1 interspacing 9mm, dia module 1.9cm wide)
	detectorD0Z = 0.725; // by definition
	detectorD1Z = 1.625; // by definition
	detectorD2Z = 18.725; // by definition
	detectorD3Z = 19.625; // by definition
	detectorDiaZ = 10.2; // by definition
	//compact geometry; edges: 0, 2.40, 6, 12, 14.40 (si modules 2.40cm wide, x/y planes spaced 2mm, D0/D1 interspacing 9mm, dia module 1.9cm wide)
	//Double_t detectorD0Z = 0.725; // by definition
	//Double_t detectorD1Z = 1.625; // by definition
	//Double_t detectorD2Z = 12.725; // by definition
	//Double_t detectorD3Z = 13.625; // by definition
	//Double_t detectorDiaZ = 7.2; // by definition


	//default paths
	sys = gSystem;
	if(verbosity)cout<<"Clustering::Clustering: Set plotspath \""<<flush;
	ostringstream plotspath;
	plotspath << sys->pwd() << "/plots-" << RunNumber;
	if(RunDescription=="") plotspath << "/";
	else plotspath << "-" << RunDescription << "/";
	plots_path = plotspath.str();
	if(verbosity)cout<<plots_path<<"\""<<endl;
	png_file_char = plotspath.str();
	C_file_char = plotspath.str();
	root_file_char = plotspath.str();
	//make plots dir
	sys->mkdir(plots_path.c_str());

	//Load Settings
	if(verbosity)cout<<"Clustering::Clustering:setSettingpath \""<<flush;
	ostringstream settingspath;
	settingspath << sys->pwd() << "/Settings." << RunNumber;
	if(RunDescription=="") settingspath << ".ini";
	else settingspath << "-" << RunDescription << ".ini";
	settings_file = settingspath.str();
	if(verbosity)cout<<settings_file<<"\""<<endl;
	settings = new TSettings(settingspath.str());

	// create file histograms.root
	if(verbosity)cout<<"Clustering::Clustering:SetRoothistoPath"<<flush;
	ostringstream root_histo_file_path;
	root_histo_file_path << plotspath.str().c_str() << "histograms.root";
	if(verbosity) cout << "creating \"" << root_histo_file_path.str().c_str() << "\" .." << endl;
	TFile fplots(root_histo_file_path.str().c_str(),"RECREATE");;

	//SetPedfilePath
	if (verbosity) cout<<"Clustering::Clustering:Get PedTree"<<endl;
	ostringstream pedfilepath;
	pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
	if(RunDescription=="")
		pedfilepath << ".root";
	else
		pedfilepath << "-" << RunDescription << ".root";
	pedfile_path = pedfilepath.str();
	PedFile = new TFile(pedfilepath.str().c_str());
	PedTree = (TTree*)PedFile->Get("PedTree");
	if (!PedTree)
	{
		cerr << "Clustering::Clustering:PedTree not found!" << endl;
	}

	if (verbosity) cout<<"Clustering::Clustering:get eventReader: "<<flush;
	eventReader = new TADCEventReader(PedTree);
	if (verbosity) cout << eventReader->GetEntries() << " events in EventReader" << endl;
	current_event = 0;

	//initialize histograms
	InitializeHistograms();

	//Initialize FidCutRegions
	for(int i=0;i<4;i++) {
		FCR[i] = new FidCutRegion(i);
		FCR[i]->SetAllValuesZero();
	}

	//initialize counters
	total_events = 0, totalsurviving_events = 0, singlesitrack_events = 0, singlesitrack_1diamondclus_events = 0, singlesitrack_fidcut_events = 0, singlesitrack_fidcut_1diamondclus_events = 0, CMNEvents = 0, ZeroDivisorEvents = 0;
	for(int i=0; i<4; i++) detectorxycluster_events[i] = 0;
	for(int d=0; d<10; d++) {
		goldengatecluster_events[d] = 0;
		badchannelcluster_events[d] = 0;
		lumpycluster_events[d] = 0;
		saturatedcluster_events[d] = 0;
	}
	//alignment counters
	counter_alignment_tracks = 0;
	counter_alignment_fidcut_tracks = 0;
	counter_alignment_only_tracks = 0;
	counter_alignment_only_fidcut_tracks = 0;
	counter_alignment_tracks_zero_suppressed = 0;

	//create alignment class
	alignment = new AlignmentClass(PedTree,0);
	alignment->SetSettings(settings);
	alignment->SetPlotsPath(plots_path);

	histSaver= new HistogrammSaver(verbosity);
	histSaver->SetRunNumber(eventReader->run_number);
	histSaver->SetNumberOfEvents((unsigned int)eventReader->GetEntries());
	histSaver->SetPlotsPath(plots_path);
}

Clustering::~Clustering() {
   PedFile->Close();
   //delete PedTree;
   //delete PedFile; // delete this line if root segfaults

   //DeleteHistograms();
}


//take current event and cluster
void Clustering::ClusterEvent(bool verbose) {
   eventReader->GetEvent(current_event);
   if(verbose) cout<<endl<<endl;
   if(current_event%1000==0) cout<<"Clustering::ClusterEvent(): current_event = "<<current_event<<endl;

   clustered_event.Clear();
   clustered_event.SetEventNumber(current_event);

   vector<int> hits, cluster, badchannelclusterflags, goldengateclusterflags, lumpyclusterflags, saturatedclusterflags;
   vector< vector<int> > clusters;
   int previouschan, currentchan, hasaseed, hasmasked, previousseed, isgoldengate, islumpy, peakchan, hassaturated;
   float peakchan_psadc, currentchan_psadc, previouschan_psadc;

   for(int det=0; det<9; det++) {
      hits.clear();
      cluster.clear();
      clusters.clear();
      
      //look for hits
      if(verbose) cout<<endl<<endl<<"Detector "<<det<<" hits: ";
      for(int i=0; i<(int)eventReader->Det_NChannels[det]; i++) {
         if(det<8)
        	 if(eventReader->Det_ADC[det][i]-eventReader->Det_PedMean[det][i] > settings->getSi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][i]) {
        		 hits.push_back(i);
        		 if(verbose) {
        			 cout<<(int)eventReader->Det_Channels[det][i];
        			 if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][i]))
        				 cout<<"(masked)";
        			 cout<<", ";
        		 }
        	 }
         if(det==8)
        	 if(eventReader->Dia_ADC[i]-eventReader->Det_PedMean[det][i] > settings->getDi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][i]) {
            hits.push_back(i);
            if(verbose) {
               cout<<(int)eventReader->Det_Channels[det][i];
               if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][i]))
            	   cout<<"(masked)";
               cout<<", ";
            }
         }
      }
      if(verbose) {
         cout<<endl<<"Channels screened: ";
         for(int i=0; i<256; i++) if(!settings->getDet_channel_screen(det).CheckChannel(i)) cout<<i<<", ";
         cout<<endl<<hits.size()<<" hits found in "<<(int)eventReader->Det_NChannels[det]<<" saved channels"<<endl;
      }
      if(hits.size()==0) {
         if(verbose) cout<<"No hits found so skipping to next detector."<<endl;
         continue;
      }
      hits.push_back(-1);
      if(verbose) cout<<"pushed back -1 as hit for determining end of hits"<<endl;

      //now look for contiguous regions in hits and save clustered hits with seeds
      if(verbose) cout<<"before clear(): cluster.size()="<<cluster.size()<<" and cluster.size()="<<cluster.size()<<endl;
      cluster.clear();
      clusters.clear();
      goldengateclusterflags.clear();
      badchannelclusterflags.clear();
      lumpyclusterflags.clear();
      saturatedclusterflags.clear();
      if(verbose) cout<<"after clear(): cluster.size()="<<cluster.size()<<" and cluster.size()="<<cluster.size()<<endl;
      previouschan=-1;
      for(uint i=0; i<hits.size(); i++) {
         currentchan = eventReader->Det_Channels[det][hits[i]];
         if(verbose) {
            if(hits[i]==-1) cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or end of hits"<<endl;
            else {
               cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or channel "<<currentchan<<" (";
               if(det==8)
            	   currentchan_psadc = eventReader->Dia_ADC[hits[i]]-eventReader->Det_PedMean[det][hits[i]];
               if(det<8)
            	   currentchan_psadc = eventReader->Det_ADC[det][hits[i]]-eventReader->Det_PedMean[det][hits[i]];
               cout<<currentchan_psadc<<" psadc, "<<eventReader->Det_PedWidth[det][hits[i]]<<" pedrms, "<<currentchan_psadc/eventReader->Det_PedWidth[det][hits[i]]<<" snr)"<<endl;
            }
         }
         //build a cluster of hits
         if((previouschan==-1 || currentchan==previouschan+1) && hits[i]!=-1) {
            if(verbose) cout<<"adding channel to cluster"<<endl;
            cluster.push_back(hits[i]);
            previouschan = currentchan;
         }
         //found end of cluster so search current cluster for a seed
         else {
            if(hits[i]!=-1) i--;
            if(verbose) cout<<"found end of cluster; looking for seed:"<<endl;
            hasaseed=0;
            hasmasked=0;
            previousseed=-1;
            isgoldengate=0;
            islumpy=0;
            hassaturated=0;
            peakchan=-1;
            peakchan_psadc=-5000;
            
            //require no masked channels adjacent to the cluster
            if((int)eventReader->Det_Channels[det][cluster[0]]==0 || (int)eventReader->Det_Channels[det][cluster[cluster.size()-1]]==255) {
               hasmasked = 1;
               if(verbose) cout<< "Cluster is up against edge of detector; flagging cluster as bad channel cluster." << endl;
            }
            else if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][cluster[0]]-1)
               || !settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][cluster[cluster.size()-1]]+1)) {
               hasmasked = 1;
               if(verbose) cout<< "Channel(s) adjacent to the cluster is masked; flagging cluster as bad channel cluster." << endl;
            }
            
            for(uint j=0; j<cluster.size(); j++) {
               currentchan = cluster[j];
               if(verbose) cout<<"cluster["<<j<<"]="<<cluster[j]<<" or channel "<<(int)eventReader->Det_Channels[det][currentchan];
               if(det<8) {
                  currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];
                  //check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
                  if(cluster.size()>2) if(currentchan_psadc > peakchan_psadc) {
                     peakchan_psadc = currentchan_psadc;
                     peakchan = j;
                  }
                  //check for seed
                  if(currentchan_psadc > settings->getSi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][currentchan]) {
                     hasaseed=1;
                     if(verbose) cout<<" is a seed";
                     if(previousseed!=-1 && eventReader->Det_Channels[det][currentchan]!=previousseed+1) {
                        isgoldengate=1;
                        if(verbose) cout<<" (goldengate cluster)";
                     }
                     else previousseed=eventReader->Det_Channels[det][currentchan];
                     //check to see if the hit saturated the adc
                     if(eventReader->Det_ADC[det][currentchan]>254) {
                        hassaturated=1;
                        if(verbose) cout<<" (saturated adc)";
                     }
                  }
                  //require no masked hits in silicon
                  if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][currentchan])) {
                	  hasmasked=1; if(verbose) cout<< " (masked hit)";
                  }
               }
               if(det==8) {
                  currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];
                  //check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
                  if(cluster.size()>2) if(currentchan_psadc > peakchan_psadc) {
                     peakchan_psadc = currentchan_psadc;
                     peakchan = j;
                  }
                  //check for seed
                  if(currentchan_psadc > settings->getDi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][currentchan]) {
                     hasaseed=1;
                     if(verbose) cout<<" is a seed";
                     if(previousseed!=-1 && eventReader->Det_Channels[det][currentchan]!=previousseed+1) {
                        isgoldengate=1;
                        if(verbose) cout<<" (goldengate cluster)";
                     }
                     else previousseed=eventReader->Det_Channels[det][currentchan];
                     //check to see if the hit saturated the adc
                     if(eventReader->Dia_ADC[currentchan]>4094) {
                        hassaturated=1;
                        if(verbose) cout<<" (saturated adc)";
                     }
                     //require no masked seeds in diamond
                     if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][currentchan])) {
                    	 hasmasked=1; if(verbose) cout<< " (masked seed)";
                     }
                  }
               }
               if(verbose) cout<<endl;
            }
            //Lumpy cluster check (can't have lumpy cluster with less than 3 hits)
            if(cluster.size()>2) {
               if(verbose) cout<<"Found peak at cluster hit "<<peakchan<<" or channel index "<<cluster[peakchan]<<" with PSADC of "<<peakchan_psadc<<endl;
               //now scan right of peak to see if monotonically decreasing
               previouschan_psadc = peakchan_psadc;
               for(uint j=peakchan+1; j<cluster.size(); j++) {
                  currentchan = cluster[j];
                  if(det<8) currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];
                  if(det==8) currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];
                  if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
                  if(currentchan_psadc>previouschan_psadc) {
                	  islumpy = 1;
                	  if(verbose)
                		  cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)eventReader->Det_Channels[det][currentchan]<<")"<<endl;
                  }
                  previouschan_psadc = currentchan_psadc;
               }
               //now scan left of peak to see if monotonically decreasing
               previouschan_psadc = peakchan_psadc;
               for(int j=peakchan-1; j>=0; j--) {
                  currentchan = cluster[j];
                  if(det<8) currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];
                  if(det==8) currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];
                  if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
                  if(currentchan_psadc>previouschan_psadc) {
                	  islumpy = 1;
                	  if(verbose)
                		  cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)eventReader->Det_Channels[det][currentchan]<<")"<<endl;
                  }
                  previouschan_psadc = currentchan_psadc;
               }
            }
            //if there's a seed in the cluster, save it
            if(hasaseed) {
               clusters.push_back(cluster);
               badchannelclusterflags.push_back(hasmasked);
               goldengateclusterflags.push_back(isgoldengate);
               lumpyclusterflags.push_back(islumpy);
               saturatedclusterflags.push_back(hassaturated);
               if(verbose) {
                  cout<<"storing cluster"<<endl;
                  if(hasmasked) cout<<"Flagged as BadChannelCluster!"<<endl;
                  if(isgoldengate) cout<<"Flagged as GoldenGateCluster!"<<endl;
                  if(islumpy) cout<<"Flagged as LumpyCluster!"<<endl;
                  if(hassaturated) cout<<"Flagged as SaturatedCluster!"<<endl;
               }
            }
            cluster.clear(); //start storing a new cluster
            previouschan=-1; //save next channel
            //move on to the next cluster
         }//end else (found end of candidate cluster)
      }//end loop over hits in detector
      
      if(clusters.size()==0) {
         if(verbose) cout<<"No clusters found so skipping to next detector."<<endl;
         continue;
      }
      
      //now that we have lists of channels belonging to clusters, let's create Cluster objects to store the data
      Cluster* current_cluster = 0;
      int highest_index, nexthighest_index;
      float highest_psadc, nexthighest_psadc, current_psadc;
      if(clusters.size()>0) {
         if(verbose) cout<<"det="<<det<<"\tclusters.size()="<<clusters.size()<<endl;
         for(uint i=0; i<clusters.size(); i++) {
            //add cluster to a different list in current event depending on flags
            current_cluster = clustered_event.AddCluster(det,badchannelclusterflags[i]); //||goldengateclusterflags[i]||lumpyclusterflags[i]||saturatedclusterflags[i]);
            //flag cluster
            current_cluster->FlagBadChannelCluster(badchannelclusterflags[i]);
            current_cluster->FlagGoldenGateCluster(goldengateclusterflags[i]);
            current_cluster->FlagLumpyCluster(lumpyclusterflags[i]);
            current_cluster->FlagSaturatedCluster(saturatedclusterflags[i]);
            //calculate transparent eta while saving each channel to cluster (note that due to "zero supression" in the saved pedestal info, not all clusters will have an eta and will be reported as -1)
            highest_index = -1; nexthighest_index = -1;
            highest_psadc = 0; nexthighest_psadc = 0; current_psadc = 0;
            for(uint j=0; j<clusters[i].size(); j++) { //save each cluster one channel at a time
               currentchan = clusters[i][j];
               if(det<8) {
                  current_cluster->AddHit(eventReader->Det_Channels[det][currentchan], eventReader->Det_ADC[det][currentchan], eventReader->Det_PedMean[det][currentchan], eventReader->Det_PedWidth[det][currentchan]);
                  current_psadc = eventReader->Det_ADC[det][currentchan] - eventReader->Det_PedMean[det][currentchan];
               }
               if(det==8) {
                  current_cluster->AddHit(eventReader->Det_Channels[det][currentchan], eventReader->Dia_ADC[currentchan], eventReader->Det_PedMean[det][currentchan], eventReader->Det_PedWidth[det][currentchan]);
                  current_psadc = eventReader->Dia_ADC[currentchan] - eventReader->Det_PedMean[det][currentchan];
               }
               if(verbose)
            	   cout<<"Detector "<<det<<": cluster_saved/all_clusters="<<i+1<<"/"<<clusters.size()<<"\tchan_to_save="<<j+1<<"/"<<clusters[i].size()<<"\tcurrentchan="<<(int)eventReader->Det_Channels[det][currentchan]<<endl;
               //locate highest seed
               if(current_psadc>highest_psadc) {
                  highest_index = currentchan;
                  highest_psadc = current_psadc;
               }
            }//end loop over channels in a cluster to save
            //check which channel has next highest psadc
            if(highest_index<(int)eventReader->Det_NChannels[det]-1)
            	if(eventReader->Det_Channels[det][highest_index+1]==eventReader->Det_Channels[det][highest_index]+1 &&
            			settings->getDet_channel_screen(det).CheckChannel(eventReader->Det_Channels[det][highest_index+1])) {
               //first if the next channel is available and an ok channel, just assume it has the nexthighest_psadc
               nexthighest_index = highest_index+1;
               if(det==8) nexthighest_psadc = eventReader->Dia_ADC[nexthighest_index] - eventReader->Det_PedMean[det][nexthighest_index];
               else nexthighest_psadc = eventReader->Det_ADC[det][nexthighest_index] - eventReader->Det_PedMean[det][nexthighest_index];
            }
            if(highest_index>0)
            	if(eventReader->Det_Channels[det][highest_index-1]==eventReader->Det_Channels[det][highest_index]-1 &&
            			settings->getDet_channel_screen(det).CheckChannel(eventReader->Det_Channels[det][highest_index-1])) {
               //now if the previous channel is available and an ok channel, check whether it has a higher psadc than the next one
               if(det==8) current_psadc = eventReader->Dia_ADC[highest_index-1] - eventReader->Det_PedMean[det][highest_index-1];
               else current_psadc = eventReader->Det_ADC[det][highest_index-1] - eventReader->Det_PedMean[det][highest_index-1];
               if(current_psadc>nexthighest_psadc) {
                  nexthighest_index = highest_index-1;
                  nexthighest_psadc = current_psadc;
               }
            }
            //calculate eta and centroid for highest 2 channels (used later for alignment and tracking)
            if(nexthighest_index>-1) {
               //eta
               if(highest_index>nexthighest_index) current_cluster->SetEta(highest_psadc/(highest_psadc+nexthighest_psadc));
               else current_cluster->SetEta(nexthighest_psadc/(highest_psadc+nexthighest_psadc));

               //centroid
               current_cluster->highest2_centroid=(eventReader->Det_Channels[det][highest_index]*highest_psadc+eventReader->Det_Channels[det][nexthighest_index]*nexthighest_psadc)/(highest_psadc+nexthighest_psadc);

               /*
               if(totalsurviving_events%1000==0) {
                  cout<<"Event "<<&current_event-1<<"; detector "<<det<<":"<<endl;
                  cout<<"highest_index = "<<highest_index<<"\thighest_psadc="<<highest_psadc<<endl;
                  cout<<"nexthighest_index = "<<nexthighest_index<<"\tnexthighest_psadc="<<nexthighest_psadc<<endl;
                  cout<<"highest_psadc/(highest_psadc+nexthighest_psadc) = "<<highest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
                  cout<<"nexthighest_psadc/(highest_psadc+nexthighest_psadc) = "<<nexthighest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
               }
               */
            }
         }//end loop over clusters
      }
   }//end loop over detectors

   //get ready for next event
   current_event++;

}

// alternative function to ClusterEvent
// channels are searched for a seed. clusters are build around located seeds.
void Clustering::ClusterEventSeeds(bool verbose) {
	eventReader->GetEvent(current_event);
	if(verbose) cout<<endl<<endl;
	if(current_event%1000==0) cout<<"Clustering::ClusterEvent(): current_event = "<<current_event<<endl;

	clustered_event.Clear();
	clustered_event.SetEventNumber(current_event);

	vector<int> hits, cluster, badchannelclusterflags, goldengateclusterflags, lumpyclusterflags, saturatedclusterflags, seeds, cluster_channels;
	vector< vector<int> > clusters;
	int previouschan, currentchan, hasaseed, hasmasked, previousseed, isgoldengate, islumpy, peakchan, hassaturated;
	float peakchan_psadc, currentchan_psadc, previouschan_psadc;

	// -- loop over detectors
	for(int det=0; det<9; det++) {
		hits.clear();
		cluster.clear();
		clusters.clear();
		seeds.clear();

		// -- loop over channels and look for seeds
		if (verbose) cout << endl << endl << "Detector " << det << " seeds: ";
		for (int i = 0; i < (int)eventReader->Det_NChannels[det]; i++) {
			if (det < 8 && eventReader->Det_ADC[det][i]-eventReader->Det_PedMean[det][i] > settings->getSi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][i]) {
				seeds.push_back(i);
				if (verbose) {
					cout<<(int)eventReader->Det_Channels[det][i];
					if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][i]))
						cout<<"(masked)";
					cout<<", ";
				}
			}
			if (det == 8 && eventReader->Dia_ADC[i]-eventReader->Det_PedMean[det][i] > settings->getDi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][i]) {
				seeds.push_back(i);
				if (verbose) {
					cout<<(int)eventReader->Det_Channels[det][i];
					if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][i]))
						cout<<"(masked)";
					cout<<", ";
				}
			}
		}	// end loop over channels

		goldengateclusterflags.clear();
		badchannelclusterflags.clear();
		lumpyclusterflags.clear();
		saturatedclusterflags.clear();
		cluster_channels.clear();

		/*		for (int i=0; i < seeds.size()-1; i++) {
		 if (TMath::Abs(seeds[i]-seeds[i+1]) > 1) {
		 goldengateclusterflags.push_back(seeds[i]);
		 }
		 }*/

		if (seeds.size() == 0) {
//			cout << "No seeds found so skipping to next detector." << endl;
			continue;
		}

		// detect cluster for every seed
		for (int i = 0; i < seeds.size(); i++) {
			bool channel_already_used = false;
			int first_hit = 0;
			int last_hit = 0;
			hasaseed=1;
			hasmasked=0;
			previousseed=-1;
			isgoldengate=0;
			islumpy=0;
			hassaturated=0;
			peakchan=-1;
			peakchan_psadc=-5000;
			cluster.clear();

			// check if channel is already used in another cluster
			for (int j = 0; j < cluster_channels.size(); j++) {
				if (seeds[i] == cluster_channels[j]) channel_already_used = true;
			}
			if (channel_already_used) continue;

			if (det < 8) {
				// look for first hit of cluster
				first_hit = seeds[i];
				for (int j = seeds[i]-1; j >= 0 && eventReader->Det_ADC[det][j]-eventReader->Det_PedMean[det][j] > settings->getSi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][j] && (int)eventReader->Det_Channels[det][j]-(int)eventReader->Det_Channels[det][j-1] == 1; j--) {
					first_hit = j;
				}

				// look for last hit of cluster and write hits to cluster
				for (int j = first_hit; j < (int)eventReader->Det_NChannels[det] && eventReader->Det_ADC[det][j]-eventReader->Det_PedMean[det][j] > settings->getSi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][j]; j++) {
					if (cluster.size() > 0 && (int)eventReader->Det_Channels[det][j]-(int)eventReader->Det_Channels[det][j-1] != 1) break;
					cluster.push_back(j);
					cluster_channels.push_back(j);
				}


/*				for (int j = seeds[i]-1; j >= 0 && EventReader->Det_ADC[det][j]-EventReader->Det_PedMean[det][j] > Si_Cluster_Hit_Factor*EventReader->Det_PedWidth[det][j]; j--) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (EventReader->Det_ADC[det][j]-EventReader->Det_PedMean[det][j] > EventReader->Det_ADC[det][j+1]-EventReader->Det_PedMean[det][j+1])
							islumpy = 1;
					}
				}
				for (int j = seeds[i]+1; j >= 0 && EventReader->Det_ADC[det][j]-EventReader->Det_PedMean[det][j] > Si_Cluster_Hit_Factor*EventReader->Det_PedWidth[det][j]; j++) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (EventReader->Det_ADC[det][j]-EventReader->Det_PedMean[det][j] > EventReader->Det_ADC[det][j-1]-EventReader->Det_PedMean[det][j-1])
							islumpy = 1;
					}
				}*/
			}
			if (det == 8) {
				// look for first hit of the cluster
				first_hit = seeds[i];
				for (int j = seeds[i]-1; j >= 0 && eventReader->Dia_ADC[j]-eventReader->Det_PedMean[det][j] > settings->getDi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][j] && (int)eventReader->Det_Channels[det][j]-(int)eventReader->Det_Channels[det][j-1] == 1; j--) {
					first_hit = j;
				}

				// look for last hit of cluster and write hits to cluster
				for (int j = first_hit; j < (int)eventReader->Det_NChannels[det] && eventReader->Dia_ADC[j]-eventReader->Det_PedMean[det][j] > settings->getDi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][j]; j++) {
					if (cluster.size() > 0 && (int)eventReader->Det_Channels[det][j]-(int)eventReader->Det_Channels[det][j-1] != 1) break;
					cluster.push_back(j);
					cluster_channels.push_back(j);
				}


/*				for (int j = seeds[i]-1; j >= 0 && EventReader->Dia_ADC[j]-EventReader->Det_PedMean[det][j] > settings->getDi_Cluster_Hit_Factor()*EventReader->Det_PedWidth[det][j]; j--) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (EventReader->Dia_ADC[j]-EventReader->Det_PedMean[det][j] > EventReader->Dia_ADC[j+1]-EventReader->Det_PedMean[det][j+1])
							islumpy = 1;
					}
				}
				for (int j = seeds[i]+1; j >= 0 && EventReader->Dia_ADC[j]-EventReader->Det_PedMean[det][j] > settings->getDi_Cluster_Hit_Factor()*EventReader->Det_PedWidth[det][j]; j++) {
					if (find(cluster_channels.begin(), cluster_channels.end(), j) == cluster_channels.end()) {
						cluster.push_back(j);
						if (EventReader->Dia_ADC[j]-EventReader->Det_PedMean[det][j] > EventReader->Dia_ADC[j-1]-EventReader->Det_PedMean[det][j-1])
							islumpy = 1;
					}
				}*/
			}
			if (cluster.size() > 0) {

				//require no masked channels adjacent to the cluster
				if((int)eventReader->Det_Channels[det][cluster[0]]==0 || (int)eventReader->Det_Channels[det][cluster[cluster.size()-1]]==255) {
					hasmasked = 1;
					if(verbose) cout<< "Cluster is up against edge of detector; flagging cluster as bad channel cluster." << endl;
				}
				else if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][cluster[0]]-1)
						|| !settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][cluster[cluster.size()-1]]+1)) {
					hasmasked = 1;
					if(verbose) cout<< "Channel(s) adjacent to the cluster is masked; flagging cluster as bad channel cluster." << endl;
				}

				// -- look for peak, ckeck golden gate, check saturation and require no masked hits in silicon
				for(uint j=0; j<cluster.size(); j++) {
					currentchan = cluster[j];
					if(verbose) cout<<"cluster["<<j<<"]="<<cluster[j]<<" or channel "<<(int)eventReader->Det_Channels[det][currentchan];
					if(det<8) {
						currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];

						//check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
//						if(cluster.size()>2) {
							if(currentchan_psadc > peakchan_psadc) {
								peakchan_psadc = currentchan_psadc;
								peakchan = j;
							}
//						}

						// check for golden gate
						if(currentchan_psadc > settings->getSi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][currentchan]) {
							if(verbose) cout<<" is a seed";
							if(previousseed!=-1 && eventReader->Det_Channels[det][currentchan]!=previousseed+1) {
								isgoldengate=1;
								if(verbose) cout<<" (goldengate cluster)";
							}
							previousseed=eventReader->Det_Channels[det][currentchan];
						}

						//check to see if the hit saturated the adc
						if(eventReader->Det_ADC[det][currentchan]>254) {
							hassaturated=1;
							if(verbose) cout<<" (saturated adc)";
						}

						//require no masked hits in silicon
						if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked hit)";}
					}
					if(det==8) {
						currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];

						//check for peak (for lumpy cluster check; can only have lumpy cluster for >2 hits)
//						if(cluster.size()>2)
						if(currentchan_psadc > peakchan_psadc) {
							peakchan_psadc = currentchan_psadc;
							peakchan = j;
						}

						//check for golden gate
						if(currentchan_psadc > settings->getDi_Cluster_Seed_Factor()*eventReader->Det_PedWidth[det][currentchan]) {
							if(verbose) cout<<" is a seed";
							if(previousseed!=-1 && eventReader->Det_Channels[det][currentchan]!=previousseed+1) {
								isgoldengate=1;
								if(verbose) cout<<" (goldengate cluster)";
							}
							previousseed=eventReader->Det_Channels[det][currentchan];
							//require no masked seeds in diamond
							if(!settings->getDet_channel_screen(det).CheckChannel((int)eventReader->Det_Channels[det][currentchan])) {
								hasmasked=1;
								if(verbose) cout<< " (masked seed)";
							}
						}

						//check to see if the hit saturated the adc
						if(eventReader->Dia_ADC[currentchan]>4094) {
							hassaturated=1;
							if(verbose) cout<<" (saturated adc)";
						}
					}
					if(verbose) cout<<endl;
				}

				//Lumpy cluster check (can't have lumpy cluster with less than 3 hits)
				if(cluster.size()>2) {
					if(verbose) cout<<"Found peak at cluster hit "<<peakchan<<" or channel index "<<cluster[peakchan]<<" with PSADC of "<<peakchan_psadc<<endl;
					//now scan right of peak to see if monotonically decreasing
					previouschan_psadc = peakchan_psadc;
					for(uint j=peakchan+1; j<cluster.size(); j++) {
						currentchan = cluster[j];
						if(det<8) currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];
						if(det==8) currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];
						if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
						if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)eventReader->Det_Channels[det][currentchan]<<")"<<endl;}
						previouschan_psadc = currentchan_psadc;
					}
					//now scan left of peak to see if monotonically decreasing
					previouschan_psadc = peakchan_psadc;
					for(int j=peakchan-1; j>=0; j--) {
						currentchan = cluster[j];
						if(det<8) currentchan_psadc = eventReader->Det_ADC[det][currentchan]-eventReader->Det_PedMean[det][currentchan];
						if(det==8) currentchan_psadc = eventReader->Dia_ADC[currentchan]-eventReader->Det_PedMean[det][currentchan];
						if(verbose) cout<<"LumpyClusterCheck: clusterhit="<<j<<"\tchanindex="<<currentchan<<"\tcurrentchan_psadc="<<currentchan_psadc<<"\tpreviouschan_psadc="<<previouschan_psadc<<endl;
						if(currentchan_psadc>previouschan_psadc) {islumpy = 1; if(verbose) cout<<"LumpyClusterCheck: Cluster is lumpy (clusterhit "<<j<<", index "<<currentchan<<", or chan "<<(int)eventReader->Det_Channels[det][currentchan]<<")"<<endl;}
						previouschan_psadc = currentchan_psadc;
					}
				}

				//if there's a seed in the cluster, save it
				if(hasaseed) {
					clusters.push_back(cluster);
					badchannelclusterflags.push_back(hasmasked);
					goldengateclusterflags.push_back(isgoldengate);
					lumpyclusterflags.push_back(islumpy);
					saturatedclusterflags.push_back(hassaturated);
					if(verbose) {
						cout<<"storing cluster"<<endl;
						if(hasmasked) cout<<"Flagged as BadChannelCluster!"<<endl;
						if(isgoldengate) cout<<"Flagged as GoldenGateCluster!"<<endl;
						if(islumpy) cout<<"Flagged as LumpyCluster!"<<endl;
						if(hassaturated) cout<<"Flagged as SaturatedCluster!"<<endl;
					}
				}
			}
		} // end loop over seeds
		
		//now that we have lists of channels belonging to clusters, let's create Cluster objects to store the data
		Cluster* current_cluster = 0;
		int highest_index, nexthighest_index;
		float highest_psadc, nexthighest_psadc, current_psadc;
		if(clusters.size()>0) {
			if(verbose) cout<<"det="<<det<<"\tclusters.size()="<<clusters.size()<<endl;
			for(uint i=0; i<clusters.size(); i++) {
				//add cluster to a different list in current event depending on flags
				current_cluster = clustered_event.AddCluster(det,badchannelclusterflags[i]); //||goldengateclusterflags[i]||lumpyclusterflags[i]||saturatedclusterflags[i]);
				//flag cluster
				current_cluster->FlagBadChannelCluster(badchannelclusterflags[i]);
				current_cluster->FlagGoldenGateCluster(goldengateclusterflags[i]);
				current_cluster->FlagLumpyCluster(lumpyclusterflags[i]);
				current_cluster->FlagSaturatedCluster(saturatedclusterflags[i]);
				//calculate transparent eta while saving each channel to cluster (note that due to "zero supression" in the saved pedestal info, not all clusters will have an eta and will be reported as -1)
				highest_index = -1; nexthighest_index = -1;
				highest_psadc = 0; nexthighest_psadc = 0; current_psadc = 0;
				for(uint j=0; j<clusters[i].size(); j++) { //save each cluster one channel at a time
					currentchan = clusters[i][j];
					if(det<8) {
						current_cluster->AddHit(eventReader->Det_Channels[det][currentchan], eventReader->Det_ADC[det][currentchan], eventReader->Det_PedMean[det][currentchan], eventReader->Det_PedWidth[det][currentchan]);
						current_psadc = eventReader->Det_ADC[det][currentchan] - eventReader->Det_PedMean[det][currentchan];
					}
					if(det==8) {
						current_cluster->AddHit(eventReader->Det_Channels[det][currentchan], eventReader->Dia_ADC[currentchan], eventReader->Det_PedMean[det][currentchan], eventReader->Det_PedWidth[det][currentchan]);
						current_psadc = eventReader->Dia_ADC[currentchan] - eventReader->Det_PedMean[det][currentchan];
					}
					if(verbose) cout<<"Detector "<<det<<": cluster_saved/all_clusters="<<i+1<<"/"<<clusters.size()<<"\tchan_to_save="<<j+1<<"/"<<clusters[i].size()<<"\tcurrentchan="<<(int)eventReader->Det_Channels[det][currentchan]<<endl;
					//locate highest seed
					if(current_psadc>highest_psadc) {
						highest_index = currentchan;
						highest_psadc = current_psadc;
					}
				}//end loop over channels in a cluster to save
				//check which channel has next highest psadc
				if(highest_index<(int)eventReader->Det_NChannels[det]-1)
					if(eventReader->Det_Channels[det][highest_index+1]==eventReader->Det_Channels[det][highest_index]+1 &&
							settings->getDet_channel_screen(det).CheckChannel(eventReader->Det_Channels[det][highest_index+1])) {
					//first if the next channel is available and an ok channel, just assume it has the nexthighest_psadc
					nexthighest_index = highest_index+1;
					if(det==8) nexthighest_psadc = eventReader->Dia_ADC[nexthighest_index] - eventReader->Det_PedMean[det][nexthighest_index];
					else nexthighest_psadc = eventReader->Det_ADC[det][nexthighest_index] - eventReader->Det_PedMean[det][nexthighest_index];
				}
				if(highest_index>0)
					if(eventReader->Det_Channels[det][highest_index-1]==eventReader->Det_Channels[det][highest_index]-1 &&
							settings->getDet_channel_screen(det).CheckChannel(eventReader->Det_Channels[det][highest_index-1])) {
					//now if the previous channel is available and an ok channel, check whether it has a higher psadc than the next one
					if(det==8)
						current_psadc = eventReader->Dia_ADC[highest_index-1] - eventReader->Det_PedMean[det][highest_index-1];
					else
						current_psadc = eventReader->Det_ADC[det][highest_index-1] - eventReader->Det_PedMean[det][highest_index-1];
					if(current_psadc>nexthighest_psadc) {
						nexthighest_index = highest_index-1;
						nexthighest_psadc = current_psadc;
					}
				}
				//calculate eta and centroid for highest 2 channels (used later for alignment and tracking)
				if(nexthighest_index>-1) {
					//eta
					if(highest_index>nexthighest_index) current_cluster->SetEta(highest_psadc/(highest_psadc+nexthighest_psadc));
					else current_cluster->SetEta(nexthighest_psadc/(highest_psadc+nexthighest_psadc));

					//centroid
					current_cluster->highest2_centroid=(eventReader->Det_Channels[det][highest_index]*highest_psadc+eventReader->Det_Channels[det][nexthighest_index]*nexthighest_psadc)/(highest_psadc+nexthighest_psadc);

					/*
					 if(totalsurviving_events%1000==0) {
					 cout<<"Event "<<&current_event-1<<"; detector "<<det<<":"<<endl;
					 cout<<"highest_index = "<<highest_index<<"\thighest_psadc="<<highest_psadc<<endl;
					 cout<<"nexthighest_index = "<<nexthighest_index<<"\tnexthighest_psadc="<<nexthighest_psadc<<endl;
					 cout<<"highest_psadc/(highest_psadc+nexthighest_psadc) = "<<highest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
					 cout<<"nexthighest_psadc/(highest_psadc+nexthighest_psadc) = "<<nexthighest_psadc/(highest_psadc+nexthighest_psadc)<<endl;
					 }
					 */
				}
			}//end loop over clusters
		}

	}	// end loop over detectors
	current_event++;
	// -- end of lukas' code

}


void Clustering::InitializeHistograms() {

	if (verbosity) cout<<"Clustering::InitializeHistogramms "<<endl;
	{//Determine which histos to make
		   for(int cut=0; cut<2; cut++) histoswitch_dianoise[cut]=1; //if there's not a hit, fill the bin; second index: no fidcut, fidcut
		   for(int det=0; det<9; det++) for(int hits=0; hits<11; hits++) for(int cut=0; cut<3; cut++) {
			   if(det==8 && (cut==0||cut==1||(cut==2&&hits==10))) histoswitch_clusterocc[det][hits][cut] = 1; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
			   else if(det<8 && (cut==0||cut==1||cut==2) && (hits==0||hits==1||hits==2||hits==3||hits==4||hits==10)) histoswitch_clusterocc[det][hits][cut] = 1;
			   else histoswitch_clusterocc[det][hits][cut] = 0;
		   }
		   for(int det=0; det<9; det++) for(int hits=0; hits<11; hits++) for(int cut=0; cut<3; cut++) for(int snr=0; snr<2; snr++) {
			   if(det==8 && (cut==0||cut==1) && snr==0) histoswitch_landau[det][hits][cut][snr] = 1; //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
			   else if(det<8 && (cut==1||cut==2) && snr==0 && (hits==0||hits==1||hits==2||hits==3||hits==4||hits==10)) histoswitch_landau[det][hits][cut][snr] = 1;
			   else if(hits==5 && cut==0 && snr==0) histoswitch_landau[det][hits][cut][snr] = 1;
			   else histoswitch_landau[det][hits][cut][snr] = 0;
		   }
		   for(int det=0; det<9; det++) for(int cut=0; cut<2; cut++) histoswitch_clusterfreq[det][cut] = 1;
		   for(int det=0; det<9; det++) for(int metric=0; metric<3; metric++) for(int cut=0; cut<3; cut++) {
			   if(metric!=2 && det==8 && (cut==0||cut==1)) histoswitch_clustersizefreq[det][metric][cut] = 1; //second index: nhits, nseeds, standard deviation; third index: no fidcut, fidcut, no track req
			   else if(metric!=2 && det<8 && cut==0) histoswitch_clustersizefreq[det][metric][cut] = 1;
			   else histoswitch_clustersizefreq[det][metric][cut] = 0;
		   }
		   for(int det=0; det<5; det++) for(int cut=0; cut<3; cut++) histoswitch_scatter[det][cut] = 1; //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut
		   for(int det=0; det<9; det++) for(int type=0; type<4; type++) for(int cut=0; cut<3; cut++) for(int chip=0; chip<2; chip++) {
			   if(det==8 && (type==0||type==1) && (cut==0||cut==1)) histoswitch_eta[det][type][cut][chip] = 1; //second index: transparent,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req; fourth index: 1st readout chip, 2nd readout chip
			   else if(det<8 && (type==0||type==1) && (cut==0||cut==1)) histoswitch_eta[det][type][cut][chip] = 1;
			   else histoswitch_eta[det][type][cut][chip] = 0;
		   }
		   for(int det=0; det<9; det++) for(int type=0; type<3; type++) for(int chip=0; chip<2; chip++) {
			   if(1) histoswitch_etaintegral[det][type][chip] = 1; //second index: transparent,2hit,low,high,low(q-dist),hi(q-dist); third index: track but no fidcut, track and fidcut, no track req
			   else histoswitch_etaintegral[det][type][chip] = 0;
		   }
		   histoswitch_trackfreq_vs_time = 0; //show how many tracks at a given time
		   for(int det=0; det<9; det++) for(int etatype=0; etatype<2; etatype++) for(int chip=0; chip<2; chip++) histoswitch_eta_vs_Q[det][etatype][chip] = 1;
		   for(int det=0; det<9; det++) for(int frac=0; frac<2; frac++) {
			   if(det<8 && frac==0) histoswitch_hitocc_saturated[det][frac] = 1;
			   else if(det==8) histoswitch_hitocc_saturated[det][frac] = 1;
			   else histoswitch_hitocc_saturated[det][frac] = 0;
		   }
		   for(int det=0; det<9; det++) for(int trackcut=0; trackcut<3; trackcut++) { //second index: no fidcut, fidcut, no track req
			   if(trackcut==0) histoswitch_clusters_average[det][trackcut] = 1;
			   else histoswitch_clusters_average[det][trackcut] = 0;
		   }
	   }

   for(int i=0; i<2; i++) histo_dianoise[i] = 0;
    for (int i = 0; i < 8; i++) histo_detnoise[i] = 0;
   for(int i=0; i<5; i++) for(int j=0; j<3; j++) histo_scatter[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) histo_clusterocc[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) histo_landau[i][j][k][l] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) histo_clusterfreq[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) histo_clustersizefreq[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<6; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) histo_eta[i][j][k][l] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<2; k++) histo_etaintegral[i][j][k] = 0;
   histo_trackfreq_vs_time = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) histo_eta_vs_Q[i][j][k] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) histo_hitocc_saturated[i][j] = 0;
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) histo_clusters_average[i][j] = 0;
	histo_scatter_autofidcut = 0;

	string histo_name, det_name, trackcut_name, nhitsopt_name, etaopt_name, clussizeopt_name, chip_name;

	   histo_dianoise[0] = new TH1F("noise_diamond","noise_diamond",80,-40,40);
	   histo_dianoise[1] = new TH1F("noise_diamond_fidcutevents","noise_diamond_fidcutevents",80,-40,40);
	    for (int det = 0; det < 8; det++) {
	        switch(det) {
	            case 0: det_name = "_D0X"; break;
	            case 1: det_name = "_D0Y"; break;
	            case 2: det_name = "_D1X"; break;
	            case 3: det_name = "_D1Y"; break;
	            case 4: det_name = "_D2X"; break;
	            case 5: det_name = "_D2Y"; break;
	            case 6: det_name = "_D3X"; break;
	            case 7: det_name = "_D3Y"; break;
	        }
	        histo_name = "Noise" + det_name;
	        histo_detnoise[det] = new TH1F(histo_name.c_str(),histo_name.c_str(),160,-20,20);
	    }

	   for(int det=0; det<5; det++) {

	      switch(det) {
	      case 0: det_name = "_D0"; break;
	      case 1: det_name = "_D1"; break;
	      case 2: det_name = "_D2"; break;
	      case 3: det_name = "_D3"; break;
	      }

	      for(int trackcut=0; trackcut<3; trackcut++) {
	         if(det==4) histo_name = "Silicon8HitsScatter_Average";
	         else histo_name = "Silicon8HitsScatter" + det_name;
	         if(trackcut==1) histo_name += "_1DiamondCluster";
	         if(trackcut==2) histo_name += "_Fidcut";
	         histo_scatter[det][trackcut] = new TH2F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5,256,-0.5,255.5); //first index: d0,d1,d2,d3,<all>; second index: si-track, si-track + diamond, fidcut
	      }
	   }
		
		histo_scatter_autofidcut = new TH2F("AutoFidCut_Silicon8HitsScatter_Average_1DiamondCluster","AutoFidCut_Silicon8HitsScatter_Average_1DiamondCluster",256,-0.5,255.5,256,-0.5,255.5);

	   for(int det=0; det<9; det++) {

	      switch(det) {
	         case 0: det_name = "_D0X"; break;
	         case 1: det_name = "_D0Y"; break;
	         case 2: det_name = "_D1X"; break;
	         case 3: det_name = "_D1Y"; break;
	         case 4: det_name = "_D2X"; break;
	         case 5: det_name = "_D2Y"; break;
	         case 6: det_name = "_D3X"; break;
	         case 7: det_name = "_D3Y"; break;
	         case 8: det_name = "_Dia"; break;
	      }


	      for(int chip=0; chip<2; chip++) {
	         if(chip) chip_name = "_RightChip";
	         else chip_name = "_LeftChip";

	         histo_name = "EtaVsQ_TransparentEta_8HitsFidcut" + det_name + chip_name;
	         if(det<8) histo_eta_vs_Q[det][0][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_si_max(),101,-0.005,1.005);
	         if(det==8) histo_eta_vs_Q[det][0][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_di_max(),101,-0.005,1.005);

	         histo_name = "EtaVsQ_2HitEta_8HitsFidcut" + det_name + chip_name;
	         if(det<8) histo_eta_vs_Q[det][1][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_si_max(),101,-0.005,1.005);
	         if(det==8) histo_eta_vs_Q[det][1][chip] = new TH2F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),0,settings->getPulse_height_di_max(),101,-0.005,1.005);

	      }
	      for(int frac=0; frac<2; frac++) {
	         histo_name = "HitOccupancy_SaturatedChannels";
	         if(frac==1) histo_name += "_RelativeFrequency";
	         histo_name += det_name;
	         if(det<8) histo_hitocc_saturated[det][frac] = new TH1F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5);
	         if(det==8) histo_hitocc_saturated[det][frac] = new TH1F(histo_name.c_str(),histo_name.c_str(),128,-0.5,127.5);
	      }

	      for(int cut=0; cut<2; cut++) {
	         histo_name = "ClusterFrequency" + det_name;
	         if(cut) histo_name += "_CutClusters";
	         histo_clusterfreq[det][cut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-0.5,10.5);
	      }

	      for(int trackcut=0; trackcut<3; trackcut++) {

	         switch(trackcut) {
	            case 0: trackcut_name = "_8Hits"; break;
	            case 1: trackcut_name = "_8HitsFidcut"; break;
	            case 2: trackcut_name = "_No8Hits"; break;
	         }

	         histo_name = "ClustersAverage" + det_name + trackcut_name;
	         histo_clusters_average[det][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-5.5,5.5);

	         for(int clussizeopt=0; clussizeopt<3; clussizeopt++) {
	            switch(clussizeopt) {
	               case 0: clussizeopt_name = "_NHits"; break;
	               case 1: clussizeopt_name = "_NSeeds"; break;
	               case 2: clussizeopt_name = "_StdDev"; break;
	            }
	            histo_name = "ClusterSizeFrequency" + det_name + clussizeopt_name + trackcut_name;
	            histo_clustersizefreq[det][clussizeopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),11,-0.5,10.5); //second index: nhits, nseeds, 2nd moment; third index: no fidcut, fidcut, no track req
	         }

	         for(int etaopt=0; etaopt<6; etaopt++) {
	            switch(etaopt) {
	               case 0: etaopt_name = "_2Largest"; break;
	               case 1: etaopt_name = "_2HitCluster"; break;
	               case 2: etaopt_name = "_LowCharge"; break;
	               case 3: etaopt_name = "_HighCharge"; break;
	               case 4: etaopt_name = "_LowCharge_Landau"; break;
	               case 5: etaopt_name = "_HighCharge_Landau"; break;
	            }

	            for(int chip=0; chip<2; chip++) {
	               if(chip) chip_name = "_RightChip";
	               else chip_name = "_LeftChip";
	               histo_name = "Eta" + det_name + chip_name + etaopt_name + trackcut_name;
	               histo_eta[det][etaopt][trackcut][chip] = new TH1F(histo_name.c_str(),histo_name.c_str(),1001,-0.0005,1.0005); //second index: 2highest,2hit,low,high; third index: track but no fidcut, track and fidcut, no track req
	            }
	         }

	         if(trackcut==0) for(int etaopt=0; etaopt<3; etaopt++) {
	            switch(etaopt) {
	               case 0: etaopt_name = "_2Largest_LeftToRight"; break;
	               case 1: etaopt_name = "_2Largest_RightToLeft"; break;
	               case 2: etaopt_name = "_2Largest_Average"; break;
	            }

	            for(int chip=0; chip<2; chip++) {
	               if(chip) chip_name = "_RightChip";
	               else chip_name = "_LeftChip";
	               histo_name = "EtaIntegral" + det_name + chip_name + etaopt_name + trackcut_name;
	               histo_etaintegral[det][etaopt][chip] = new TH1F(histo_name.c_str(),histo_name.c_str(),1001,-0.0005,1.0005); //second index: left to right, right to left, average
	            }
	         }

	         for(int nhitsopt=0; nhitsopt<11; nhitsopt++) {

	            switch(nhitsopt) {
	               case 0: nhitsopt_name = "_1HitClusters"; break;
	               case 1: nhitsopt_name = "_2HitClusters"; break;
	               case 2: nhitsopt_name = "_3HitClusters"; break;
	               case 3: nhitsopt_name = "_4HitClusters"; break;
	               case 4: nhitsopt_name = "_5orMoreHitClusters"; break;
	               case 5: nhitsopt_name = "_1-2HitClusters"; break;
	               case 6: nhitsopt_name = "_1-3HitClusters"; break;
	               case 7: nhitsopt_name = "_1-4HitClusters"; break;
	               case 8: nhitsopt_name = "_2and3HitClusters"; break;
	               case 9: nhitsopt_name = "_3and4HitClusters"; break;
	               case 10: nhitsopt_name = "_AllHitClusters"; break;
	            }

	            histo_name = "ClusterCentroidOccupancy" + det_name + nhitsopt_name + trackcut_name;
	            if(det==8) histo_clusterocc[det][nhitsopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),128,-0.5,127.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req
	            else histo_clusterocc[det][nhitsopt][trackcut] = new TH1F(histo_name.c_str(),histo_name.c_str(),256,-0.5,255.5);//-0.5/5,255.5/5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req

	            for(int snr=0; snr<2; snr++) {
	               if(snr) {
	                  histo_name = "SNRDistribution" + det_name + nhitsopt_name + trackcut_name;
	                  if(det==8) histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getSnr_distribution_di_max()+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
	                  else histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getSnr_distribution_si_max()+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
	               }
	               else {
	                  histo_name = "PulseHeight" + det_name + nhitsopt_name + trackcut_name;
	                  if(det==8) histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
	                  else histo_landau[det][nhitsopt][trackcut][snr] = new TH1F(histo_name.c_str(),histo_name.c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_si_max()+0.5); //second index: 1,2,3,4,>=5,1&2,1-3,1-4,2&3,3&4,all; third index: no fidcut, fidcut, no track req; fourth index: psadc, snr
	               }
	            }
	         }
	      }
	   }

	   /*	cout << "init new histos..";
	      	for (int i = 0; i < 5; i++) {
	      		ostringstream histoname_landau, histoname_eta;
	      		histoname_landau << "PulseHeight_Dia_" << (i+1) << "HitTransparClusters";
	      		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
	      		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),2000,0.,2000.);
	      		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
	      		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	      		histo_transparentclustering_eta[i] = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),2000,0.,2000.);
	      	}
	      	cout << " done." << endl;*/


	   	histo_afc_x = new TH1F("histo_afc_x","histo_afc_x",256,0,255); // for afc (max, 19.11.2010)
	   	histo_afc_y = new TH1F("histo_afc_y","histo_afc_y",256,0,255); // for afc (max, 19.11.2010)
	   	histo_afc_clone = new TH2F("histo_afc_clone","histo_afc_clone",256,0,255,256,0,255);

	   	histo_afc_x_cut = new TH1F("histo_afc_x_cut","histo_afc_x_cut",256,0,255);
	   	histo_afc_y_cut = new TH1F("histo_afc_y_cut","histo_afc_y_cut",256,0,255);

	       histo_afc_unit_histo_1f = new TH1F("histo_afc_unit_histo_1f","histo_afc_unit_histo_1f",256,0,255);
	   	histo_afc_unit_histo_2f = new TH2F("histo_afc_unit_histo_2f","histo_afc_unit_histo_2f",256,0,255,256,0,255);
	   	histo_afc_region_1_mask = new TH2F("histo_afc_region_1_mask","histo_afc_region_1_mask",256,0,255,256,0,255);

}

void Clustering::DeleteHistograms() {
   for(int i=0; i<2; i++) delete histo_dianoise[i];
    for (int i = 0; i < 8; i++) delete histo_detnoise[i];
   for(int i=0; i<2; i++) delete histofit_dianoise[i];
   for(int i=0; i<5; i++) for(int j=0; j<3; j++) delete histo_scatter[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) delete histo_clusterocc[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<11; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) delete histo_landau[i][j][k][l];
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) delete histo_clusterfreq[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) delete histo_clustersizefreq[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<6; j++) for(int k=0; k<3; k++) for(int l=0; l<2; l++) delete histo_eta[i][j][k][l];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) for(int k=0; k<2; k++) delete histo_etaintegral[i][j][k];
   delete histo_trackfreq_vs_time;
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) delete histo_eta_vs_Q[i][j][k];
   for(int i=0; i<9; i++) for(int j=0; j<2; j++) delete histo_hitocc_saturated[i][j];
   for(int i=0; i<9; i++) for(int j=0; j<3; j++) delete histo_clusters_average[i][j];
	delete histo_scatter_autofidcut;
}

void Clustering::DrawHistograms() {
	if (verbosity>1)cout<<"Clustering::DrawHistograms()"<<endl;
   string tempname = "tempname";
   for(int cut=0; cut<2; cut++) if(histoswitch_dianoise[cut]) {
      tempname += "0";
      histofit_dianoise[cut] = new TF1(tempname.c_str(),"gaus",histo_dianoise[cut]->GetMean()-histo_dianoise[cut]->GetRMS(),histo_dianoise[cut]->GetMean()+histo_dianoise[cut]->GetRMS());
      histofit_dianoise[cut]->SetLineColor(kBlue);
      histo_dianoise[cut]->Fit(tempname.c_str(),"r"); // fit option "r" restricts the range of the fit
      histSaver->SaveHistogramPNG(histo_dianoise[cut]);
   }
    
    for (int det = 0; det < 8; det++) {
        histSaver->SaveHistogram(histo_detnoise[det]);
    }
   
   for(int det=0; det<5; det++)
      for(int trackcut=0; trackcut<3; trackcut++)
         if(histoswitch_scatter[det][trackcut]) histSaver->SaveHistogramPNG(histo_scatter[det][trackcut]);
   
   gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
   for(int det=0; det<9; det++) {
      
      for(int frac=0; frac<2; frac++)
         if(histoswitch_hitocc_saturated[det][frac]) histSaver->SaveHistogramPNG(histo_hitocc_saturated[det][frac]);
      
      for(int goodcut=0; goodcut<2; goodcut++)
         if(histoswitch_clusterfreq[det][goodcut]) histSaver->SaveHistogramPNG(histo_clusterfreq[det][goodcut]);
         
      histo_clusters_average[det][0]->SetNormFactor(histo_clusters_average[det][0]->Integral(-100,100)/singlesitrack_events);

  	if (verbosity>2)cout<<"Clustering::DrawHistograms()::forchip"<<endl;
      for(int chip=0; chip<2; chip++) {
         for(int etatype=0; etatype<2; etatype++) {
            if(histoswitch_eta_vs_Q[det][etatype][chip]) histSaver->SaveHistogramPNG(histo_eta_vs_Q[det][etatype][chip]);
         }
         
         //agregate eta plots

     	if (verbosity>2)cout<<"Clustering::DrawHistograms()::for chip::generate eta plots"<<endl;
         TCanvas etacanvas("etacanvas");
         if(histoswitch_etaintegral[det][2][chip]) histo_etaintegral[det][2][chip]->Draw();
         if(histoswitch_etaintegral[det][0][chip]) histo_etaintegral[det][0][chip]->Draw("same");
         if(histoswitch_etaintegral[det][1][chip]) histo_etaintegral[det][1][chip]->Draw("same");
         if (verbosity>2)cout<<"Clustering::DrawHistograms()::for chip::generate eta plots::draw pt"<<endl;
         //TODO ERSATZ FINDEN
         //pt->Draw();
         ostringstream plot_filename;
         if(chip) plot_filename << plots_path << "EtaIntegrals_D" << det << "_RightChip.png";
         else plot_filename << plots_path << "EtaIntegrals_D" << det << "_LeftChip.png";
         if (verbosity>2)cout<<"Clustering::DrawHistograms()::for chip::generate eta plots::print canvas"<<endl;
         etacanvas.Print(plot_filename.str().c_str());
         
         //individual eta plots_canvas
         if (verbosity)cout<<"Clustering::DrawHistograms()::for chip::generate eta plots::individual eta plots_canvas"<<endl;
         if(histo_etaintegral[det][0][chip]) histSaver->SaveHistogramPNG(histo_etaintegral[det][0][chip]);
         if(histo_etaintegral[det][2][chip]) histSaver->SaveHistogramPNG(histo_etaintegral[det][2][chip]);
            
      }
  	if (verbosity>2)cout<<"Clustering::DrawHistograms()::for trackcut"<<endl;
      
      for(int trackcut=0; trackcut<3; trackcut++) {

         if(histoswitch_clusters_average[det][trackcut]) histSaver->SaveHistogramPNG(histo_clusters_average[det][trackcut]);
         
         for(int clussizeopt=0; clussizeopt<3; clussizeopt++)
            if(histoswitch_clustersizefreq[det][clussizeopt][trackcut]) histSaver->SaveHistogramPNG(histo_clustersizefreq[det][clussizeopt][trackcut]);
         
         for(int etaopt=0; etaopt<6; etaopt++) {
            if(det==8) if(histoswitch_eta[det][etaopt][trackcut][0]) histSaver->SaveHistogramPNG(histo_eta[det][etaopt][trackcut][0]);
            if(det<8) for(int chip=0; chip<2; chip++)
               if(histoswitch_eta[det][etaopt][trackcut][chip]) histSaver->SaveHistogramPNG(histo_eta[det][etaopt][trackcut][chip]);
         }
         
         for(int nhitsopt=0; nhitsopt<11; nhitsopt++) {
            //hitoccs
            if(histoswitch_clusterocc[det][nhitsopt][trackcut]) 
            	histSaver->SaveHistogramPNG(histo_clusterocc[det][nhitsopt][trackcut]);
            
            //pulse heights
            for(int snr=0; snr<2; snr++) 
               if(histoswitch_landau[det][nhitsopt][trackcut][snr]) histSaver->SaveHistogram(histo_landau[det][nhitsopt][trackcut][snr]);
         }
      }
   }
}

void Clustering::BookHistograms() {
   
   TDetectorPlane D0, D1, D2, D3, Dia;
   D0.SetZ(detectorD0Z);
   D1.SetZ(detectorD1Z);
   D2.SetZ(detectorD2Z);
   D3.SetZ(detectorD3Z);
   Dia.SetZ(detectorDiaZ);
   bool mask;
   
   //cut flow counters
   //Int_t total_events, goldengatecluster_events, badchannelcluster_events, singlesitrack_events, singlesitrack_1diamondclus_events, singlesitrack_fidcut_events, singlesitrack_fidcut_1diamondclus_events;
   
   //count total events
   total_events++;
   
   //count bad events
   //if(clustered_event.HasSaturatedCluster()) saturatedcluster_events[9]++;
   
   
   //cut events containing the following clusters
   //silicon: saturated, lumpy, goldengate
   //diamond: saturated,
   if(eventReader->CMNEvent_flag || eventReader->ZeroDivisorEvent_flag || clustered_event.HasSaturatedCluster() || clustered_event.HasSaturatedCluster(8) || clustered_event.HasLumpyCluster() || clustered_event.HasGoldenGateCluster() /*|| clustered_event.HasBadChannelCluster()*/) {//|| !clustered_event.HasOneSiliconTrack()) {//clustered_event.HasBadChannelCluster()) {
      if(eventReader->CMNEvent_flag) CMNEvents++;
      if(eventReader->ZeroDivisorEvent_flag) ZeroDivisorEvents++;
      
      //count cut types of events here
      
      return;
   }
   
   //count various types of events surviving the above cuts
   if(clustered_event.HasGoldenGateCluster()) goldengatecluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasGoldenGateCluster(det)) goldengatecluster_events[det]++;
   if(clustered_event.HasBadChannelCluster()) badchannelcluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasBadChannelCluster(det)) badchannelcluster_events[det]++;
   if(clustered_event.HasLumpyCluster()) lumpycluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasLumpyCluster(det)) lumpycluster_events[det]++;
   if(clustered_event.HasSaturatedCluster()) saturatedcluster_events[9]++;
   for(int det=0; det<9; det++) if(clustered_event.HasSaturatedCluster(det)) {
      saturatedcluster_events[det]++;
      for(uint clus=0; clus<clustered_event.GetNCutClusters(det); clus++)
         for(int frac=0; frac<2; frac++)
            if(clustered_event.GetCutCluster(det,clus)->IsSaturatedCluster()) for(int hit=0; hit<clustered_event.GetCutCluster(det,clus)->GetNHits(); hit++) {
               if(det<8)
            	   if(clustered_event.GetCutCluster(det,clus)->GetADC(hit)>4094 &&
            			   settings->getDet_channel_screen(det).CheckChannel(clustered_event.GetCutCluster(det,clus)->GetChannel(hit))) histo_hitocc_saturated[det][frac]->Fill(clustered_event.GetCutCluster(det,clus)->GetChannel(hit));
               if(det==8)
            	   if(clustered_event.GetCutCluster(det,clus)->GetADC(hit)>4094 &&
            			   settings->getDet_channel_screen(det).CheckChannel(clustered_event.GetCutCluster(det,clus)->GetChannel(hit)) && clustered_event.GetCutCluster(det,clus)->GetChannel(hit)>63) histo_hitocc_saturated[det][frac]->Fill(clustered_event.GetCutCluster(det,clus)->GetChannel(hit)-63);
            }
   }
   totalsurviving_events++;
   
   Float_t centroid, charge, snr, stddev, eta, eta2ch = -1;
   Int_t nhits, nseeds, clus_peak;
   Cluster* current_cluster = 0;
   
   //simple tracking: require one and only one cluster in each silicon plane
   bool one_and_only_one = clustered_event.HasOneSiliconTrack();
   if(one_and_only_one) singlesitrack_events++;
   //else if((current_event-1)%1000==0) clustered_event.Print();
   
   //look to see if there are coincident hits in the silicon planes
   for(int d=0; d<4; d++)
      if(clustered_event.HasCoincidentClustersXY(d)) detectorxycluster_events[d]++;
    
    // silicon noise plots
    for (int det = 0; det < 8; det++) {
        for (int i = 0; i < eventReader->Det_NChannels[det]; i++) {
            float psadc = eventReader->Det_ADC[det][i] - eventReader->Det_PedMean[det][i];
            if (settings->getDet_channel_screen(det).CheckChannel(eventReader->Det_Channels[det][i]) && psadc<settings->getSi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[det][i]) {
                histo_detnoise[det]->Fill(psadc);
            }
        }
    }
   
   //if a track, check whether it's in the silicon fiducial region
   bool fiducial_track = 0;
   Float_t si_avg_x=0, si_avg_y=0, det_centroid[8];
   if(one_and_only_one) {
      for(int det=0; det<4; det++) {
         si_avg_x += clustered_event.GetCluster(2*det,0)->Get1stMoment();
         si_avg_y += clustered_event.GetCluster(2*det+1,0)->Get1stMoment();
      }
      si_avg_x = si_avg_x/4;
      si_avg_y = si_avg_y/4;
      
      if(	si_avg_x>settings->getSi_avg_fidcut_xlow() &&
    		si_avg_x<settings->getSi_avg_fidcut_xhigh() &&
    		si_avg_y>settings->getSi_avg_fidcut_ylow() &&
    		si_avg_y<settings->getSi_avg_fidcut_yhigh()) {
         fiducial_track=1;
         singlesitrack_fidcut_events++;
      }
      
      //silicon track clusters' centroids
      for(int det=0; det<8; det++)
         det_centroid[det] = clustered_event.GetCluster(det,0)->Get1stMoment();
      
      //silicon track scatter for each plane
      for(int det=0; det<4; det++)
         histo_scatter[det][0]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
      //silicon track average scatter
      histo_scatter[4][0]->Fill(si_avg_x,si_avg_y);
      
      //scatter of tracks in the diamond
      if(clustered_event.GetNClusters(8)==1) {
         //scatter of tracks in diamond in each plane
         for(int det=0; det<4; det++)
            histo_scatter[det][1]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
         //scatter of tracks in diamond averaged over silicon planes
         histo_scatter[4][1]->Fill(si_avg_x,si_avg_y);
		  //check usual events
//		  if (si_avg_x < 80) {
//			  EventMonitor(current_event-1);
//		  }
         //count diamond track events
         singlesitrack_1diamondclus_events++;
         //count fidcut diamond track events
         if(fiducial_track) singlesitrack_fidcut_1diamondclus_events++;
      }
      
      //scatter of diamond tracks in the silicon fiducial region
      if(fiducial_track && clustered_event.GetNClusters(8)==1) {
         //scatter of diamond tracks in the silicon fiducial region
         for(int det=0; det<4; det++) 
            histo_scatter[det][2]->Fill(det_centroid[2*det],det_centroid[2*det+1]);
         //scatter of diamond tracks in the silicon fiducial region averaged over silicon planes
         histo_scatter[4][2]->Fill(si_avg_x,si_avg_y);
      }
      
      //diamond noise plots
      float psadc;
      for(uint i=0; i<eventReader->Det_NChannels[8]; i++) {
         psadc=eventReader->Dia_ADC[i]-eventReader->Det_PedMean[8][i];
         if(settings->getDet_channel_screen(8).CheckChannel(eventReader->Det_Channels[8][i]) && psadc<settings->getDi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[8][i]) {
            histo_dianoise[0]->Fill(psadc);
            if(fiducial_track) histo_dianoise[1]->Fill(psadc);
         }
      }
      
      //save tracks for alignment
      
      bool zero_suppressed_track = 0;
      for(int det=0; det<8; det++)
         if(clustered_event.GetCluster(det,0)->highest2_centroid==-1)
            zero_suppressed_track=1;
         
      if(zero_suppressed_track) //skip events where pedestal calc zero-suppression prevents eta correction
         counter_alignment_tracks_zero_suppressed++;
      
      else { //if track is ok, then add to list of alignment tracks
         
         //now tracking with centroid of the seed and the highest neighbor (for eta correction)
         D0.SetX(clustered_event.GetCluster(0,0)->highest2_centroid);
         D1.SetX(clustered_event.GetCluster(2,0)->highest2_centroid);
         D2.SetX(clustered_event.GetCluster(4,0)->highest2_centroid);
         D3.SetX(clustered_event.GetCluster(6,0)->highest2_centroid);
         
         D0.SetY(clustered_event.GetCluster(1,0)->highest2_centroid);
         D1.SetY(clustered_event.GetCluster(3,0)->highest2_centroid);
         D2.SetY(clustered_event.GetCluster(5,0)->highest2_centroid);
         D3.SetY(clustered_event.GetCluster(7,0)->highest2_centroid);

         //now tracks in diamond
         if(clustered_event.GetNClusters(8)==1 && fiducial_track) {
            if(clustered_event.GetCluster(8,0)->highest2_centroid!=-1) {
               Dia.SetX(clustered_event.GetCluster(8,0)->highest2_centroid);
               TDiamondTrack track_fidcut = TDiamondTrack(clustered_event.event_number,D0,D1,D2,D3,Dia);
               tracks_fidcut.push_back(track_fidcut);
               counter_alignment_fidcut_tracks++;
               mask = rand.Uniform()<settings->getAlignment_training_track_fraction();
               tracks_fidcut_mask.push_back(mask);
               counter_alignment_only_fidcut_tracks += mask;
            }
         }
         else
            Dia.SetX(-1);
         
         TDiamondTrack track = TDiamondTrack(clustered_event.event_number,D0,D1,D2,D3);
         //track.SetEventNumber(Silicon_tracks[t].Event_number);
         tracks.push_back(track);
         counter_alignment_tracks++;
         mask = rand.Uniform()<settings->getAlignment_training_track_fraction();
         tracks_mask.push_back(mask);
         counter_alignment_only_tracks += mask;
         
      }//end else
      
   }//end if one_and_only_one
   
   uint nclusters;
   //loop over detectors
   for(int det=0; det<9; det++) {
      
      //cluster frequency
      nclusters=clustered_event.GetNClusters(det);
      histo_clusterfreq[det][0]->Fill(nclusters);
      histo_clusterfreq[det][1]->Fill(clustered_event.GetNCutClusters(det));
      
      //loop over all *good* clusters
      for(uint clus=0; clus<nclusters; clus++) {
         
         current_cluster = clustered_event.GetCluster(det,clus);
         
         //if(current_cluster->IsBadChannelCluster()) continue; // skip bad clusters
         
         nhits = current_cluster->GetNHits();
         nseeds = current_cluster->GetNSeeds();
         stddev = current_cluster->GetStandardDeviation();
         centroid = current_cluster->Get1stMoment();
         charge = current_cluster->GetCharge();
         snr = current_cluster->GetSNR();
         eta = current_cluster->GetEta();
         //if(totalsurviving_events%1000==0) cout<<"Event "<<current_event-1<<": Transparent eta for detector "<<det<<" is "<<eta<<endl;
         if(nhits==2) eta2ch = (current_cluster->GetADC(1)-current_cluster->GetPedMean(1))/current_cluster->GetCharge();
         clus_peak = current_cluster->GetPeakHit();
         
         //all good clusters
         //------------
         
         //cluster size frequency
         histo_clustersizefreq[det][0][2]->Fill(nhits);
         histo_clustersizefreq[det][1][2]->Fill(nseeds);
         histo_clustersizefreq[det][2][2]->Fill(stddev);
         
         //transparent eta
         if(eta!=-1) {
            if(centroid<128) histo_eta[det][0][2][0]->Fill(eta); //first chip
            else histo_eta[det][0][2][1]->Fill(eta); //second chip
         }
         
         //cluster property frequency plots
         if(nhits==1) {
            histo_clusterocc[det][0][2]->Fill(centroid); // fill 1 cluster histo
            histo_landau[det][0][2][0]->Fill(charge);
            histo_landau[det][0][2][1]->Fill(snr);
            histo_clusterocc[det][5][2]->Fill(centroid); // fill 1&2 cluster histo
            histo_landau[det][5][2][0]->Fill(charge);
            histo_landau[det][5][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==2) {
            if(centroid<128) histo_eta[det][1][2][0]->Fill(eta2ch); // fill 2hit cluster eta for first chip
            else histo_eta[det][1][2][1]->Fill(eta2ch); // fill 2hit cluster eta for second chip
            histo_clusterocc[det][1][2]->Fill(centroid); // fill 2 cluster histo
            histo_landau[det][1][2][0]->Fill(charge);
            histo_landau[det][1][2][1]->Fill(snr);
            histo_clusterocc[det][5][2]->Fill(centroid); // fill 1&2 cluster histo
            histo_landau[det][5][2][0]->Fill(charge);
            histo_landau[det][5][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][8][2]->Fill(centroid); // fill 2&3 cluster histo
            histo_landau[det][8][2][0]->Fill(charge);
            histo_landau[det][8][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==3) {
            histo_clusterocc[det][2][2]->Fill(centroid); // fill 3 cluster histo
            histo_landau[det][2][2][0]->Fill(charge);
            histo_landau[det][2][2][1]->Fill(snr);
            histo_clusterocc[det][6][2]->Fill(centroid); // fill 1-3 cluster histo
            histo_landau[det][6][2][0]->Fill(charge);
            histo_landau[det][6][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][8][2]->Fill(centroid); // fill 2&3 cluster histo
            histo_landau[det][8][2][0]->Fill(charge);
            histo_landau[det][8][2][1]->Fill(snr);
            histo_clusterocc[det][9][2]->Fill(centroid); // fill 3&4 cluster histo
            histo_landau[det][9][2][0]->Fill(charge);
            histo_landau[det][9][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits==4) {
            histo_clusterocc[det][3][2]->Fill(centroid); // fill 4 cluster histo
            histo_landau[det][3][2][0]->Fill(charge);
            histo_landau[det][3][2][1]->Fill(snr);
            histo_clusterocc[det][7][2]->Fill(centroid); // fill 1-4 cluster histo
            histo_landau[det][7][2][0]->Fill(charge);
            histo_landau[det][7][2][1]->Fill(snr);
            histo_clusterocc[det][9][2]->Fill(centroid); // fill 3&4 cluster histo
            histo_landau[det][9][2][0]->Fill(charge);
            histo_landau[det][9][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         if(nhits>4) {
            histo_clusterocc[det][4][2]->Fill(centroid); // fill >=5 cluster histo
            histo_landau[det][4][2][0]->Fill(charge);
            histo_landau[det][4][2][1]->Fill(snr);
            histo_clusterocc[det][10][2]->Fill(centroid); // fill all cluster histo
            histo_landau[det][10][2][0]->Fill(charge);
            histo_landau[det][10][2][1]->Fill(snr);
         }
         
         //cut on clean tracks
         //-------------------
         
         if(one_and_only_one) {
            
            //cluster size frequency
            histo_clustersizefreq[det][0][0]->Fill(nhits);
            histo_clustersizefreq[det][1][0]->Fill(nseeds);
            histo_clustersizefreq[det][2][0]->Fill(stddev);
            
            //transparent eta
            if(eta!=-1) {
               if(centroid<128) histo_eta[det][0][0][0]->Fill(eta); //first chip
               else histo_eta[det][0][0][1]->Fill(eta); //second chip
            }
            
            //stack clusters
            for(int h=0; h<nhits; h++)
               histo_clusters_average[det][0]->Fill(h-clus_peak,current_cluster->GetPSADC(h));
            
            //cluster property frequency plots
            if(nhits==1) {
               histo_clusterocc[det][0][0]->Fill(centroid); // fill 1 cluster histo
               histo_landau[det][0][0][0]->Fill(charge);
               histo_landau[det][0][0][1]->Fill(snr);
               histo_clusterocc[det][5][0]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][0][0]->Fill(charge);
               histo_landau[det][5][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==2) {
               if(centroid<128) histo_eta[det][1][0][0]->Fill(eta2ch); // fill 2hit cluster eta
               else histo_eta[det][1][0][1]->Fill(eta2ch); // fill 2hit cluster eta
               histo_clusterocc[det][1][0]->Fill(centroid); // fill 2 cluster histo
               histo_landau[det][1][0][0]->Fill(charge);
               histo_landau[det][1][0][1]->Fill(snr);
               histo_clusterocc[det][5][0]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][0][0]->Fill(charge);
               histo_landau[det][5][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][8][0]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][0][0]->Fill(charge);
               histo_landau[det][8][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==3) {
               histo_clusterocc[det][2][0]->Fill(centroid); // fill 3 cluster histo
               histo_landau[det][2][0][0]->Fill(charge);
               histo_landau[det][2][0][1]->Fill(snr);
               histo_clusterocc[det][6][0]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][0][0]->Fill(charge);
               histo_landau[det][6][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][8][0]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][0][0]->Fill(charge);
               histo_landau[det][8][0][1]->Fill(snr);
               histo_clusterocc[det][9][0]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][0][0]->Fill(charge);
               histo_landau[det][9][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits==4) {
               histo_clusterocc[det][3][0]->Fill(centroid); // fill 4 cluster histo
               histo_landau[det][3][0][0]->Fill(charge);
               histo_landau[det][3][0][1]->Fill(snr);
               histo_clusterocc[det][7][0]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][0][0]->Fill(charge);
               histo_landau[det][7][0][1]->Fill(snr);
               histo_clusterocc[det][9][0]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][0][0]->Fill(charge);
               histo_landau[det][9][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
            if(nhits>4) {
               histo_clusterocc[det][4][0]->Fill(centroid); // fill >=5 cluster histo
               histo_landau[det][4][0][0]->Fill(charge);
               histo_landau[det][4][0][1]->Fill(snr);
               histo_clusterocc[det][10][0]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][0][0]->Fill(charge);
               histo_landau[det][10][0][1]->Fill(snr);
            }
         }
         
         //cut fiducial tracks
         //-------------------
         
         if(fiducial_track && (det<8 || (det==8 && nclusters==1))) {
            
            //cluster size frequency
            histo_clustersizefreq[det][0][1]->Fill(nhits);
            histo_clustersizefreq[det][1][1]->Fill(nseeds);
            histo_clustersizefreq[det][2][1]->Fill(stddev);
			 
//			 cout << "Clustersize from event " << current_event << " filled in hist_clustersizefreq[8][*][1].." << endl;
            
            //transparent eta
            if(eta!=-1) {
               if(centroid<128) {
                  histo_eta[det][0][1][0]->Fill(eta);
                  histo_eta_vs_Q[det][0][0]->Fill(charge,eta); //fill eta vs charge plot
               }
               else {
                  histo_eta[det][0][1][1]->Fill(eta);
                  histo_eta_vs_Q[det][0][1]->Fill(charge,eta); //fill eta vs charge plot
               }
            }
            
            
            //cluster property frequency plots
            if(nhits==1) {
               histo_clusterocc[det][0][1]->Fill(centroid); // fill 1 cluster histo
               histo_landau[det][0][1][0]->Fill(charge);
               //histo_landau[det][0][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][0][1][1]->Fill(snr);
               histo_clusterocc[det][5][1]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][1][0]->Fill(charge);
               //histo_landau[det][5][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][5][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==2) {
               if(centroid<128) {
                  histo_eta[det][1][1][0]->Fill(eta2ch); // fill 2hit cluster eta
                  histo_eta_vs_Q[det][1][0]->Fill(charge,eta2ch); //fill eta vs charge plot
               }
               else {
                  histo_eta[det][1][1][1]->Fill(eta2ch); // fill 2hit cluster eta
                  histo_eta_vs_Q[det][1][1]->Fill(charge,eta2ch); //fill eta vs charge plot
               }
               histo_clusterocc[det][1][1]->Fill(centroid); // fill 2 cluster histo
               histo_landau[det][1][1][0]->Fill(charge);
               histo_landau[det][1][1][1]->Fill(snr);
               histo_clusterocc[det][5][1]->Fill(centroid); // fill 1&2 cluster histo
               histo_landau[det][5][1][0]->Fill(charge);
               //histo_landau[det][5][1][0]->Fill(current_cluster->GetTotalADC());
               histo_landau[det][5][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][8][1]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][1][0]->Fill(charge);
               histo_landau[det][8][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==3) {
               histo_clusterocc[det][2][1]->Fill(centroid); // fill 3 cluster histo
               histo_landau[det][2][1][0]->Fill(charge);
               histo_landau[det][2][1][1]->Fill(snr);
               histo_clusterocc[det][6][1]->Fill(centroid); // fill 1-3 cluster histo
               histo_landau[det][6][1][0]->Fill(charge);
               histo_landau[det][6][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][8][1]->Fill(centroid); // fill 2&3 cluster histo
               histo_landau[det][8][1][0]->Fill(charge);
               histo_landau[det][8][1][1]->Fill(snr);
               histo_clusterocc[det][9][1]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][1][0]->Fill(charge);
               histo_landau[det][9][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits==4) {
               histo_clusterocc[det][3][1]->Fill(centroid); // fill 4 cluster histo
               histo_landau[det][3][1][0]->Fill(charge);
               histo_landau[det][3][1][1]->Fill(snr);
               histo_clusterocc[det][7][1]->Fill(centroid); // fill 1-4 cluster histo
               histo_landau[det][7][1][0]->Fill(charge);
               histo_landau[det][7][1][1]->Fill(snr);
               histo_clusterocc[det][9][1]->Fill(centroid); // fill 3&4 cluster histo
               histo_landau[det][9][1][0]->Fill(charge);
               histo_landau[det][9][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
            if(nhits>4) {
               histo_clusterocc[det][4][1]->Fill(centroid); // fill >=5 cluster histo
               histo_landau[det][4][1][0]->Fill(charge);
               histo_landau[det][4][1][1]->Fill(snr);
               histo_clusterocc[det][10][1]->Fill(centroid); // fill all cluster histo
               histo_landau[det][10][1][0]->Fill(charge);
               histo_landau[det][10][1][1]->Fill(snr);
            }
         }
      }//end loop over clusters
   }//end loop over detectors
}

void Clustering::GenerateHTML() {
	if (verbosity>1)cout<<"Clustering::GenerateHTML()"<<endl;
   if(eventReader==NULL) return;
   int section, subsection;
   
   //summary page
   if (verbosity>2)cout<<"Clustering::GenerateHTML():sethtml_plotspath"<<endl;
   string html_summary_path = plots_path + "index.html";
   ofstream html_summary(html_summary_path.c_str());
   if (verbosity>2)cout << "Clustering::GenerateHTML():Run summary HTML created at " << html_summary_path << endl;
   
   html_summary << "<html>" << endl; 
   if (verbosity>2)cout << "Clustering::GenerateHTML():start html" <<eventReader<< endl;

   html_summary << "<title>Run "<<eventReader->run_number<<" analysis results - Summary</title>" << endl;
   if (verbosity>2)cout << "Clustering::GenerateHTML():done with eventReader" << endl;
   html_summary << "<body bgcolor=\"White\">" << endl;
   html_summary << "<font face=\"Arial,Sans\" size=2>"<<endl;
   html_summary << "<a name=\"top\"></a>" << endl;

   if (verbosity>2)cout<<"Clustering::GenerateHTML():start summary"<<endl;
   html_summary << "<b>Summary</b>||"
                << "<a href=\"d8.html\">Diamond</a>||"
                << "<a href=\"d0.html\">D0X</a>||"
                << "<a href=\"d1.html\">D0Y</a>||"
                << "<a href=\"d2.html\">D1X</a>||"
                << "<a href=\"d3.html\">D1Y</a>||"
                << "<a href=\"d4.html\">D2X</a>||"
                << "<a href=\"d5.html\">D2Y</a>||"
                << "<a href=\"d6.html\">D3X</a>||"
                << "<a href=\"d7.html\">D3Y</a><p>"<<endl;
   html_summary << "<center><font color=\"#000000\"><h1>Run "<<eventReader->run_number<<" analysis results - Summary</h1></font></center>"<<endl;
   html_summary << "<hr size=\"10\" Color=\"#ffff33\">" << endl; 
   html_summary << "Results from " << dateandtime.GetMonth() << "/ " << dateandtime.GetDay() << "/" << dateandtime.GetYear() << " at " << dateandtime.GetHour() << ":" << dateandtime.GetMinute() << ":" << dateandtime.GetSecond() << ")<p>" << endl;
   
   if (verbosity>2)cout<<"Clustering::GenerateHTML():Settings for pedestal calc"<<endl;
   html_summary<<"<h3>Settings for pedestal calculation</h3>"<<endl;
   if(settings->getDia_input()==0) html_summary<<"Dia0 or sirocco input 4 selected for analysis<br>"<<endl;
   if(settings->getDia_input()==1) html_summary<<"Dia1 or sirocco input 5 selected for analysis<br>"<<endl;
   html_summary<<"Buffer size of "<<settings->getIter_Size()<<" events used for pedestal analysis.<br>"<<endl;
   if(settings->getFix_dia_noise()>=0) html_summary<<"Noise for all diamond channels FIXED at "<<settings->getFix_dia_noise()<<"<br>"<<endl;
   if(eventReader!=NULL)html_summary<<"Silicon channel adc and pedestal info saved for psadc fluctuations >"<<eventReader->store_threshold<<"sigma<br>"<<endl;
   html_summary<<"For pedestal analysis, a silicon hit is >"<< settings->getSi_Pedestal_Hit_Factor() <<"sigma and a diamond hit is >"<< settings->getDi_Pedestal_Hit_Factor() <<"sigma<br>"<<endl;
   html_summary<<"Events with diamond common mode noise fluctations >"<<settings->getCMN_cut()<<"sigma flagged for removal from subsequent clustering analysis<br>"<<endl;
   if(settings->getTaylor_speed_throttle()==1) html_summary<<"Taylor's RMS speed tweak disabled<br>"<<endl;
   else html_summary<<"Taylor's RMS speed tweak enabled; RMS recalculated normally every "<<settings->getTaylor_speed_throttle()<<" events<br>"<<endl;
   
   html_summary<<"<h3>Settings for clustering</h3>"<<endl;
   html_summary<<"For silicon cluster analysis, a hit is >"<< settings->getSi_Cluster_Hit_Factor() <<"sigma and a seed is >"<< settings->getSi_Cluster_Seed_Factor() <<"sigma<br>"<<endl;
   html_summary<<"For diamond cluster analysis, a hit is >"<< settings->getDi_Cluster_Hit_Factor() <<"sigma and a seed is >"<< settings->getDi_Cluster_Seed_Factor() <<"sigma<br>"<<endl;
   html_summary<<"The following diamond channels were screened in the clustering analysis: ";
   for(int i=0; i<128; i++) {
      if(!settings->getDet_channel_screen(8).CheckChannel(i)) html_summary<<i<<",";
   }
   html_summary<<"<br>"<<endl;
   
   html_summary<<"<h3>Cut flow</h3>"<<endl;
   //html_summary<<"<ul>"<<endl;//1
   html_summary<<"<b>Events slated for analysis = "<<total_events<<"</b>"<<endl;
   html_summary<<"<ul>"<<endl;//2
   html_summary<<"<li><b>Events cut from analysis= "<<total_events-totalsurviving_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//3a
   html_summary<<"<li>CMNEvents = "<<CMNEvents<<"</li>"<<endl;
   html_summary<<"<li>ZeroDivisorEvents = "<<ZeroDivisorEvents<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//3a
   html_summary<<"<li><b>Events surviving above cuts analyzed = "<<totalsurviving_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//3b
   html_summary<<"<li>Golden gate cluster events = "<<goldengatecluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4a
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" golden gate cluster events = "<<goldengatecluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia golden gate cluster events = "<<goldengatecluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4a
   html_summary<<"<li>Bad channel cluster events = "<<badchannelcluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4b
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" bad channel cluster events = "<<badchannelcluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia bad channel cluster events = "<<badchannelcluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4b
   html_summary<<"<li>Lumpy cluster events = "<<lumpycluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4c
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" lumpy cluster events = "<<lumpycluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia lumpy cluster events = "<<lumpycluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4c
   html_summary<<"<li>Saturated cluster events = "<<saturatedcluster_events[9]<<"</li>"<<endl;
   html_summary<<"<ul>"<<endl;//4d
   for(int d=0; d<8; d++) {
      html_summary<<"<li>D"<<d/2;
      if(d%2) html_summary<<"Y";
      else html_summary<<"X";
      html_summary<<" saturated cluster events = "<<saturatedcluster_events[d]<<"</li>"<<endl;
   }
   html_summary<<"<li>Dia saturated cluster events = "<<saturatedcluster_events[8]<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//4d
   for(int d=0; d<4; d++)
      html_summary<<"<li>detectorxycluster_events["<<d<<"] = "<<detectorxycluster_events[d]<<"</li>"<<endl;
   html_summary<<"<li><b>Events with a simple <u>track</u> or <i>one and only one cluster in each silicon</i> = "<<singlesitrack_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//4c
   html_summary<<"<li>Events with single track in silicon and single cluster in diamond = "<<singlesitrack_1diamondclus_events<<"</li>"<<endl;
   html_summary<<"<li><b>Events with single track in silicon fiducial region = "<<singlesitrack_fidcut_events<<"</b></li>"<<endl;
   html_summary<<"<ul>"<<endl;//5
   html_summary<<"<li><b>Events with single track in silicon fiducial region and single cluster in diamond = "<<singlesitrack_fidcut_1diamondclus_events<<"</b></li>"<<endl;
   if(singlesitrack_fidcut_1diamondclus_events!=histo_clusterocc[8][10][1]->GetEntries()) cout<<"Clustering::GenerateHTML: We've got a problem here. "<<singlesitrack_fidcut_1diamondclus_events<<" is not equal to "<<histo_clusterocc[8][10][1]->GetEntries()<<endl;
   html_summary<<"<ul>"<<endl;//6
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 1-hit cluster in diamond = "<<(float)histo_clusterocc[8][0][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 2-hit cluster in diamond = "<<(float)histo_clusterocc[8][1][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 3-hit cluster in diamond = "<<(float)histo_clusterocc[8][2][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single 4-hit cluster in diamond = "<<(float)histo_clusterocc[8][3][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"<li>Fraction of events with single track in silicon fiducial region and single >4-hit cluster in diamond = "<<(float)histo_clusterocc[8][4][1]->GetEntries()/histo_clusterocc[8][10][1]->GetEntries()<<"</li>"<<endl;
   html_summary<<"</ul>"<<endl;//6
   html_summary<<"</ul>"<<endl;//5
   html_summary<<"</ul>"<<endl;//4c
   html_summary<<"</ul>"<<endl;//3b
   html_summary<<"</ul>"<<endl;//2
   //html_summary<<"</ul>"<<endl;//1
   html_summary<<"<hr size=\"5\">" << endl;
   
   
   // Plot index
   if (verbosity>2)cout<<"Clustering::GenerateHTML():plot index"<<endl;
   section=0;
   subsection=1;
   html_summary << "<a name=\"plot_index\"></a>" << endl;
   //html_summary << "<font color=\"f00000\"><a href=\"silicon.html\">For silicon plots click here.</a></font><p>" << endl;
   html_summary << "<p><h2>Plot Index</h2></p>" << endl;
   html_summary << "<p><a href=\"#Pedestal_Data_jump\">("<<subsection++<<") Diamond pedestal plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Hitocc_Data_jump\">("<<subsection++<<") Diamond track occupancy plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Scatter_Data_jump\">("<<subsection++<<") Silicon scatter plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Pulse_Height_jump\">("<<subsection++<<") Diamond Pulse Height Landau Plots </a></p>" << endl;
   html_summary << "<p><a href=\"#Pulse_Height_fid_jump\">("<<subsection++<<") Diamond Pulse Height Landau Plots with Silicon Fiducial Cut </a></p>" << endl;
   html_summary << "<p><a href=\"#Charge_Dist_jump\">("<<subsection++<<") Adjacent Channel Cluster Asymmetry Plots </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   
   // Diamond pedestal plots 
   if (verbosity>2)cout<<"Clustering::GenerateHTML():diamond plots"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Pedestal_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Diamond pedestal plots</h1>" << endl;
   html_summary << "<p><h3>("<<section<<"."<<subsection++<<") Diamond Pedestal Mean</h3> Diamond pedestal means by Channel</p>" << endl;
   html_summary << "<p><img SRC=\"Pedestal_Values.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<") Diamond Pedestal RMS (Initial)</h3> A diamond pedestal RMS by channel for initial values" << endl;
   html_summary << "<p><img SRC=\"Channel_Noise_Initial.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<") Diamond Pedestal RMS (Final)</h3> A diamond pedestal RMS by channel for final calculation (from running calculation)" << endl;
   html_summary << "<p><img SRC=\"Channel_Noise_Final.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the raw ADC values verses event number (a zoomed in section)." << endl;
   html_summary << "<p><img SRC=\"Raw_ADC_vs_Event_zoom.png\" BORDER=\"1\"></p>";
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the raw ADC values versus event number with CMN events removed." << endl;
   html_summary << "<p><img SRC=\"Raw_ADC_vs_Event_CMN_cut_zoom.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the pedestal subtracted ADC values vesus event number (zoom)" << endl;
   html_summary << "<p><img SRC=\"PS_ADC_vs_Event_zoom.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   
   // Diamond track occupancy plots
   if (verbosity>2)cout<<"Clustering::GenerateHTML():diamond occupancy plots"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Hitocc_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Diamond occupancy plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A display of the number of raw hits (>"<<settings->getDi_Pedestal_Hit_Factor()<<"sigma) in each channel that the pedestal analysis found." << endl;
   html_summary << "<p><img SRC=\"hit_occup_can_dia.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all events in the diamond." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all single cluster events in the diamond corresponding to <b>tracks</b> in the silicon." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Locations of clusters' centroids for all single cluster events in the diamond corresponding to tracks in the silicon <b>fiducial</b> region." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
   
   // Silicon scatter plots
   if (verbosity>2)cout<<"Clustering::GenerateHTML()::Silicon scatter plots"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Scatter_Data_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Silicon scatter plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A Scatter Plot of the Average X values for the 1 and 2 silicon planes versus an average of the y values for the 1 and 2 silicon planes for events that go through all 8 silicon planes" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the above but restricted to events where the event had a single diamond cluster." << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the above but restricted to events where the event had a single diamond cluster and is in the fiducial cut on the diamond.  A cut on the average value of the X and Y channels of the silicon is made for tracks through all 8 planes and the diamond.</p>" << endl;
   html_summary << "<p>The range in silicon space is: Channels " << settings->getSi_avg_fidcut_xlow() << " to " << settings->getSi_avg_fidcut_xhigh() << " on the x-axis and " << settings->getSi_avg_fidcut_ylow() << " to " << settings->getSi_avg_fidcut_yhigh() << " on the y-axis</p>" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
    
   // Full diamond landaus
   if (verbosity>2)cout<<"Clustering::GenerateHTML()::diamond landaus"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Pulse_Height_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Full Diamond Pulse Height Landau Plots</h1>" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of diamond cluster size or number of hits (>"<<settings->getDi_Cluster_Hit_Factor()<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of number of seeds (>"<<settings->getDi_Cluster_Seed_Factor()<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>all channel size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][10][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][0][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][1][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][2][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>4 channel clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][3][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>>4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][4][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>both 1 and 2 hit clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][5][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, and 3 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][6][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, 3 and 4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][7][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 and 3 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][8][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 and 4 hit size clusters for the entire diamond</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][9][0][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
     
   // Fidcut diamond landaus
   if (verbosity>2)cout<<"Clustering::GenerateHTML()::fidcut diamond landaus"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Pulse_Height_fid_jump\"></a>" << endl;
   html_summary << "<h1>("<<section<<") Fiducial Cut Pulse Height and Noise</h1>" << endl;
   html_summary << "<p>The following plots are taken from a cut in the silicon planes. A cut on the average value of the X and Y channels of the silicon is made for tracks through all 8 planes and the diamond.</p>" << endl;
   html_summary << "<p>The range in silicon space is: Channels " << settings->getSi_avg_fidcut_xlow() ;
   html_summary << " to " << settings->getSi_avg_fidcut_xhigh() << " on the x-axis and " << settings->getSi_avg_fidcut_ylow() << " to " << settings->getSi_avg_fidcut_yhigh() << " on the y-axis</p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<p><h3>("<<section<<"."<<subsection++<<")</h3> A scatter plot of the X and Y channels of the silicon surrounding the diamond with the fiducial cut.</p>" << endl;
   html_summary << "<p><img SRC=\""<<histo_scatter[4][2]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Projection of the silicon fiducial region in the diamond." << endl;
   html_summary << "<p><img SRC=\""<<histo_clusterocc[8][10][1]->GetTitle()<<".png\" BORDER=\"1\"></p>";
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of diamond cluster size or number of hits (>"<<settings->getDi_Cluster_Hit_Factor()<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][0][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of number of seeds (>"<<settings->getDi_Cluster_Seed_Factor()<<"sigma) in each diamond cluster" << endl;
   html_summary << "<p><img SRC=\""<<histo_clustersizefreq[8][1][1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>all channel size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][10][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][0][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][1][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][2][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>4 channel clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][3][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>>4 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][4][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>both 1 and 2 hit clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][5][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, and 3 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][6][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>1, 2, 3 and 4 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][7][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>2 and 3 hit size clusters for the silicon fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][8][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A pulse height distribution of <b>3 and 4 hit size clusters for the fiducial</b> channel range (minus noisy channels)" << endl;
   html_summary << "<p><img SRC=\""<<histo_landau[8][9][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> A distribution of non-hit ADC values (<b>noise</b>) for all channels when there was a track in the silicon fiducial region.  The width of the distribution describes the noise in the diamond signal." << endl;
   html_summary << "<p><img SRC=\""<<histo_dianoise[1]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<p><a href=\"#plot_index\">Back to plot index </a></p>" << endl;
   html_summary << "<hr size=\"5\">" << endl;
      
   // Charge sharing plots
   if (verbosity>2)cout<<"Clustering::GenerateHTML()::charge sharing plots"<<endl;
   section++;
   subsection=1;
   html_summary << "<a name=\"Charge_Dist_jump\"></a>" << endl;
   html_summary << "<p><h1>("<<section<<") Adjacent Channel Cluster Asymmetry Analysis</h1></p>" << endl;
   html_summary << "<p>This section shows the left to right distribution of charge for one, two, and three channel clusters and for the transparency hits. The main analysis is looking at the distribution of charge for two channel clusters, that is, which channel tends to have more charge deposited than the other. When ordered 0-128, the channels to the left are the lower number indexed channels and the higher are the right. The plot is made by taking the amount of charge in the left channel and dividing by the total. The histogram fills in bins less than 0.5 for left-biased clusters and greater than 0.5 for right-biased clusters. For one channel clusters, the two channels used are the actual cluster, and then the next highest adjacent channel. For three channel clusters, the two largest adjacent charge values are used. For the transparency, any hit above 5 sigma is used as a seed and the next highest adjacent channel is also take.</p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Two Channel Cluster Distribution" << endl;
   html_summary << "<p><img SRC=\""<<histo_eta[8][1][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Transparency Distribution" << endl;
   html_summary << "<p><img SRC=\""<<histo_eta[8][0][1][0]->GetTitle()<<".png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   /*
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> One Channel Cluster and Next Highest Channel Distribution" << endl;
   html_summary << "<p><img SRC=\"One_Channel_Cluster_and_next_highest_channel_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster Two Highest Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_two_highest_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster High to Low Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_hitolow_Distribution.png\" BORDER=\"1\"></p>" << endl;
   html_summary << "<hr size=\"1\">" << endl;
   html_summary << "<h3>("<<section<<"."<<subsection++<<")</h3> Three Channel Cluster Two Lowest Channels Distribution" << endl;
   html_summary << "<p><img SRC=\"Three_Channel_Cluster_twolow_Distribution.png\" BORDER=\"1\"></p>" << endl;
   */
   html_summary << "<hr size=\"5\">" << endl;

   html_summary << "</font>"<<endl;
   html_summary << "</HTML>" << endl;
   if (verbosity>2)cout<<"Clustering::GenerateHTML()::close html_summary"<<endl;
   html_summary.close();
   
}



/******************88
Note that the following ClusterRun is a flexible script;
for example, we could draw histograms and generate html
for the first 100k events, then change plot directory
and plot the second 100k events, 
********************/

//sequentially cluster all events
void Clustering::ClusterRun(bool plots) {
	if(verbosity)	cout<<"Clustering::ClusterRun() \""<<plots_path<<"\""<<endl;
   //Start Timer
   TStopwatch Watch;
   Watch.Start();
   
   //record large clusters
   string large_clusters_filename = plots_path + "large_clusters.txt";
   ofstream large_clusters_file(large_clusters_filename.c_str());
   large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
   uint nclusters, nhits;
   Cluster* current_cluster;
   
   //histograms of hit positions to find out if the eta correction helped
   TH1F* histo_hitpos[9][2];
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         ostringstream histotitle;
         histotitle << "Hit_Position_Interstrip_D" << det << "_chip" << chip;
         histo_hitpos[det][chip] = new TH1F(histotitle.str().c_str(), histotitle.str().c_str(), 101, -0.005, 1.005);
      }
   }
	

	current_event = 0;

   //loop over events
	if(verbosity)	cout<<"Clustering::ClusterRun()::StartLoop"<<endl;
   for(uint e=0; e<eventReader->GetEntries(); e++) {
   //for(uint e=0; e<10000; e++) {
//	   bool AlternativeClustering = false;
	   if (!settings->getAlternativeClustering()) ClusterEvent();
	   else ClusterEventSeeds();
      if(e%10000==0) clustered_event.Print();
      BookHistograms();
	   
	   // -- output for events flagged by ClusterEventSeeds() but not by ClusterEvent() (TODO: REMOVE THIS SECTION)
//	   if (event_number == 1706 || event_number == 3609 || event_number == 7509 || event_number == 9734) {
//		   cout << "event " << event_number << " processed.." << endl;
//		   cout << "GoldenGate: " << clustered_event.HasGoldenGateCluster(-1) << endl;
//		   cout << "BadChannel: " << clustered_event.HasBadChannelCluster(-1) << endl;
//		   cout << "Lumpy: " << clustered_event.HasLumpyCluster(-1) << endl;
//		   cout << "Saturated: " << clustered_event.HasSaturatedCluster(-1) << endl;
//		   cout << "OneSiliconTrack: " << clustered_event.HasOneSiliconTrack() << endl;
//	   }// -- end of flag output
      
      //record large clusters
      //if a track, check whether it's in the silicon fiducial region
      if(clustered_event.HasOneSiliconTrack()) {
         
         //loop over detectors to record large clusters to a txt file
         for(int det=0; det<9; det++) {
            nclusters = clustered_event.GetNClusters(det);
            //loop over all *good* clusters
            for(uint clus=0; clus<nclusters; clus++) {
               current_cluster = clustered_event.GetCluster(det,clus);
               nhits = current_cluster->GetNHits();
               if(nhits>8) {
                  large_clusters_file << eventReader->event_number << "\t" << det << "\t" << clus << "/" << nclusters << "\t" << nhits << "\t" << current_cluster->GetNSeeds();
                  if(current_cluster->IsBadChannelCluster()) large_clusters_file << " (bad chan cluster)"; // skip bad clusters
                  large_clusters_file << endl;
               }
            }
         }//end large cluster record
         
      }
   }
   
   //make the hitoccupancy plots for the saturated hits
   for(int det=0; det<9; det++) {
      for(int bin=0; bin<=histo_hitocc_saturated[det][1]->GetNbinsX(); bin++) {
         histo_hitocc_saturated[det][1]->SetBinContent(bin,(float)histo_hitocc_saturated[det][1]->GetBinContent(bin)/total_events);
      }
      for(int frac=0; frac<2; frac++)
         if(det==8) histo_hitocc_saturated[det][frac]->GetXaxis()->SetRangeUser(1,63);
   }
   
   //generate the integrated etas
   float leftintegral, rightintegral;
   int nbins;
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         nbins = histo_eta[det][0][0][chip]->GetNbinsX();
         for(int bin=1; bin<=nbins; bin++) {
            leftintegral = histo_eta[det][0][0][chip]->Integral(1,bin)/histo_eta[det][0][0][chip]->Integral(1,nbins);
            rightintegral = histo_eta[det][0][0][chip]->Integral(nbins-bin+1,nbins)/histo_eta[det][0][0][chip]->Integral(1,nbins);
            histo_etaintegral[det][0][chip]->SetBinContent(bin,leftintegral);
            histo_etaintegral[det][1][chip]->SetBinContent(bin,rightintegral);
            histo_etaintegral[det][2][chip]->SetBinContent(bin,(leftintegral+rightintegral)/2.);
         }
      }
   }
   
   //eta correct hits
   double hit_position;
   for(uint t=0; t<tracks.size(); t++) {
      for(int det=0; det<8; det++) {
         hit_position = tracks[t].GetDetectorHitPosition(det);
         if(0&&t%1000==0) {
            cout<<histo_etaintegral[det][2][0]->GetTitle()<<endl;
            cout<<hit_position<<endl;
            cout<<int(hit_position)<<endl;
            cout<<int(100*(hit_position-int(hit_position)))<<endl;
            cout<<int(100*(hit_position-int(hit_position))+0.5)<<endl;
            cout<<histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
            cout<<histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
         }
         //left to right integrated eta
         if(hit_position<128) {
            nbins = histo_etaintegral[det][0][0]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            histo_hitpos[det][0]->Fill(hit_position-int(hit_position));
         }
         else {
            nbins = histo_etaintegral[det][0][1]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            histo_hitpos[det][1]->Fill(hit_position-int(hit_position));
         }
         //symmetrzed
         //if(hit_position<128) hit_position = int(hit_position) + histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //else hit_position = int(hit_position) + histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         tracks[t].SetDetectorHitPosition(det, hit_position);
      }
   }
	cout<<"Clustering::ClusterRun():for tracks_fidcut"<<endl;
   for(uint t=0; t<tracks_fidcut.size(); t++) {
      for(int det=0; det<9; det++) {
         hit_position = tracks_fidcut[t].GetDetectorHitPosition(det);
         if(0&&t%1000==0) {
            cout<<histo_etaintegral[det][2][0]->GetTitle()<<endl;
            cout<<hit_position<<endl;
            cout<<int(hit_position)<<endl;
            cout<<int(100*(hit_position-int(hit_position)))<<endl;
            cout<<int(100*(hit_position-int(hit_position))+0.5)<<endl;
            cout<<histo_etaintegral[det][2][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
            cout<<histo_etaintegral[det][2][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5))<<endl;
         }
         //left to right integrated eta
         //if(hit_position<128) hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //else hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
         //symmetrzed
         if(hit_position<128) {
            nbins = histo_etaintegral[det][0][0]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][0]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            if(det==8) histo_hitpos[det][0]->Fill(hit_position-int(hit_position));
         }
         else {
            nbins = histo_etaintegral[det][0][1]->GetNbinsX();
            //hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(100*(hit_position-int(hit_position))+0.5));
            hit_position = int(hit_position) + histo_etaintegral[det][0][1]->GetBinContent(int(nbins*(hit_position-int(hit_position))+0.5));
            if(det==8) histo_hitpos[det][1]->Fill(hit_position-int(hit_position));
         }
         tracks_fidcut.at(t).SetDetectorHitPosition(det, hit_position);
      }
   }
   
//	AutoFidCut(); //added 2010-12-01 (max)
	
   //Draw plots and generate HTML
   if(plots) {
	   if (verbosity)cout<<"Clustering::ClusterRun():draw plots"<<endl;
      DrawHistograms();
	   if (verbosity)cout<<"Clustering::ClusterRun():generateHtml"<<endl;
      
      
      GenerateHTML();
      ////DeleteHistograms(); //this causes segfaults for some unknown reason
   }
   TCanvas tempcanv("tempcanv");
   for(int det=0; det<9; det++) {
      for(int chip=0; chip<2; chip++) {
         ostringstream tempsstream;
         tempsstream << plots_path << histo_hitpos[det][chip]->GetTitle() << ".png";
         histo_hitpos[det][chip]->Draw();
         tempcanv.Print(tempsstream.str().c_str());
      }
   }
   
   large_clusters_file.close();
   PedFile->Close(); // somehow histo_etaintegral[0][2][0] is deleted by this line?!
   
   
   cout<<"Cut flow is as follows:"<<endl;
   cout<<"total_events = "<<total_events<<endl;
   cout<<"[ ]CMNEvents = "<<CMNEvents<<endl;
   cout<<"[ ]ZeroDivisorEvents = "<<ZeroDivisorEvents<<endl;
   cout<<"[a]totalsurviving_events = "<<totalsurviving_events<<endl;
   cout<<"   [ ]goldengatecluster_events[9] = "<<goldengatecluster_events[9]<<endl;
   cout<<"   [ ]badchannelcluster_events[9] = "<<badchannelcluster_events[9]<<endl;
   cout<<"   [ ]lumpycluster_events[9] = "<<lumpycluster_events[9]<<endl;
   cout<<"   [ ]saturatedcluster_events[9] = "<<saturatedcluster_events[9]<<endl;
   for(int d=0; d<4; d++)
      cout<<"   [ ]detectorxycluster_events["<<d<<"] = "<<detectorxycluster_events[d]<<endl;
   cout<<"   [b]singlesitrack_events = "<<singlesitrack_events<<endl;
   cout<<"      [ ]singlesitrack_1diamondclus_events = "<<singlesitrack_1diamondclus_events<<endl;
   cout<<"      [c]singlesitrack_fidcut_events = "<<singlesitrack_fidcut_events<<endl;
   cout<<"         [ ]singlesitrack_fidcut_1diamondclus_events = "<<singlesitrack_fidcut_1diamondclus_events<<endl;
         
   Watch.Stop();
   Watch.Print("u");
   cout<<"tracks.size() = "<<tracks.size()<<endl;
   cout<<"tracks_fidcut.size() = "<<tracks_fidcut.size()<<endl;
   cout<<"counter_alignment_tracks_zero_suppressed = "<<counter_alignment_tracks_zero_suppressed<<endl;
   cout<<"counter_alignment_tracks = "<<counter_alignment_tracks<<endl;
   cout<<"counter_alignment_only_tracks = "<<counter_alignment_only_tracks<<endl;
   cout<<"counter_alignment_only_tracks/float(counter_alignment_tracks) = "<<counter_alignment_only_tracks/float(counter_alignment_tracks)<<endl;
   cout<<"counter_alignment_fidcut_tracks = "<<counter_alignment_fidcut_tracks<<endl;
   cout<<"counter_alignment_only_fidcut_tracks = "<<counter_alignment_only_fidcut_tracks<<endl;
   cout<<"counter_alignment_only_fidcut_tracks/float(counter_alignment_fidcut_tracks) = "<<counter_alignment_only_fidcut_tracks/float(counter_alignment_fidcut_tracks)<<endl;
   cout<<"singlesitrack_fidcut_events = "<<singlesitrack_fidcut_events<<endl;
   
}


void Clustering::HistCleaner(int regid, TH2F* histo) {
    
    
	cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
	cout << endl << "HistCleaner called on histogram " << histo->GetName() << "." << endl;
	
	int afc_current_histo_x_colsums[256];
	int afc_current_histo_y_colsums[256];
	
	TH1F* afc_current_histo_x = new TH1F("afc_current_histo_x","afc_current_histo_x",256,0,255);
	TH1F* afc_current_histo_y = new TH1F("afc_current_histo_y","afc_current_histo_y",256,0,255);
	
	for(int i=0;i<256;i++) {
		afc_current_histo_x_colsums[i] = 0;
		for(int j=0;j<256;j++) {
			afc_current_histo_x_colsums[i] = afc_current_histo_x_colsums[i] + histo->GetBinContent(i,j);
			afc_current_histo_x->SetBinContent(i,afc_current_histo_x_colsums[i]);
		}
		
	}
	
	//calculate x bin values for each y row and count # of non-zero bins and set the corresponding histogram bins
	for(int j=0;j<256;j++) {
		afc_current_histo_y_colsums[j] = 0;
		for(int i=0;i<256;i++) {
			afc_current_histo_y_colsums[j] = afc_current_histo_y_colsums[j] + histo->GetBinContent(i,j);
			afc_current_histo_y->SetBinContent(j,afc_current_histo_y_colsums[j]);
		}
		
	}
	
	int xpos = 2;
	int ypos = 2;
    
	cout << endl;
	cout << "Results for histo " << histo->GetName() << " with region id " << regid << " are:" << endl << endl;
    
	
	for(int i=xpos;i<256;i++) {
		if( (afc_current_histo_x->GetBinContent(i-1) == 0) && (afc_current_histo_x->GetBinContent(i) > 0) ) {
			//FCR[z]->SetValueXLow(i+5);
			cout << "Setting " << histo->GetName() << ".ValueXLow to: " << i+5 << endl;
			for(int j=i+5;j<256;j++) {
				if ( (afc_current_histo_x->GetBinContent(j-1) > 0) && (afc_current_histo_x->GetBinContent(j) == 0) ) {
                    //					FCR[z]->SetValueXHigh(j-5);
					cout << "Setting " << histo->GetName() << ".ValueXHigh to: " << j-5 << endl;					
					xpos=j+5;
					break;
				}} break;
			//break;
		}}
    
	for(int j=ypos;j<256;j++) {
		if( (afc_current_histo_y->GetBinContent(j-1) == 0) && (afc_current_histo_y->GetBinContent(j) > 0) ) {
			//FCR[z]->SetValueXLow(i+5);
			cout << "Setting " << histo->GetName() << ".ValueYLow to: " << j+5 << endl;
			for(int k=j+5;k<256;k++) {
				if ( (afc_current_histo_y->GetBinContent(k-1) > 0) && (afc_current_histo_y->GetBinContent(k) == 0) ) {
					//					FCR[z]->SetValueXHigh(j-5);
					cout << "Setting " << histo->GetName() << ".ValueYHigh to: " << k-5 << endl;					
					ypos=k+5;
					break;
				}} break;
			//break;
		}}
	
	
	afc_current_histo_x->Delete();
	afc_current_histo_y->Delete();
    
	cout << endl << "Finished HistCleaner on " << histo->GetName() << "." << endl;  
	cout << "HistCleaner exiting..." << endl << endl;
}




void Clustering::AutoFidCut() {
		
	histo_afc_unit_histo_1f->Reset(); // reset 1f unit histogram (all bins zero)
	histo_afc_unit_histo_2f->Reset(); // reset 2f unit histogram (all bins zero)
    
	//SaveHistogramPDF(histo_afc_unit_histo_1f);
	//SaveHistogramPDF(histo_afc_unit_histo_2f);
	
	//define sum and counter vectors for x columns and y rows
	int histo_afc_col_sums_x[256];
	int histo_afc_col_sums_y[256];
	
	
    if(verbosity>0) {
    cout << "\n";
	cout << "AutoFidCut: I'm the AutoFidCut function.\n\n";
	cout << " START AutoFidCut \n\n";	
    }    

    if(verbosity>0) {    
	// produce necessary plots to detect fidcut region
	cout << endl << endl << "-- produce plot for AutoFidCut().." << endl;
	cout << "running over " << eventReader->GetEntries() << " events.." << endl;
    }
    
	//record large clusters
//	string large_clusters_filename = plots_path + "large_clusters.txt";
//	ofstream large_clusters_file(large_clusters_filename.c_str());
//	large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
//	uint nclusters, nhits;
//	Cluster* current_cluster;
	
	current_event = 0;
	for (uint e=0; e<eventReader->GetEntries(); e++) { // maybe it's not necessary to run over all events for the AutoFidCut?!
		if (!settings->getAlternativeClustering()) ClusterEvent();
		else ClusterEventSeeds();
		if (e%10000==0) clustered_event.Print();
		
		// current_cluster = 0;
		
		// -- produce scatter plot for AutoFidCut
		bool one_and_only_one = clustered_event.HasOneSiliconTrack();
		Float_t si_avg_x=0, si_avg_y=0;
		if (eventReader->CMNEvent_flag || eventReader->ZeroDivisorEvent_flag || clustered_event.HasSaturatedCluster() || clustered_event.HasSaturatedCluster(8) || clustered_event.HasLumpyCluster() || clustered_event.HasGoldenGateCluster() || clustered_event.HasBadChannelCluster())
			continue;
		if (one_and_only_one) {
			for (int det=0; det<4; det++) {
				si_avg_x += clustered_event.GetCluster(2*det,0)->Get1stMoment();
				si_avg_y += clustered_event.GetCluster(2*det+1,0)->Get1stMoment();
			}
			si_avg_x = si_avg_x/4;
			si_avg_y = si_avg_y/4;
			
			if (clustered_event.GetNClusters(8)==1)
				histo_scatter_autofidcut->Fill(si_avg_x,si_avg_y);
		}
	}
    
    if(verbosity>0) {
	cout << "plot production for AutoFidCut: done." << endl << endl;
    }
	
    histSaver->SaveHistogram(histo_scatter_autofidcut);
	
	// -- end of scatter plot production
    
    //calculate y bin values for each x column and count # of non-zero bins and set the corresponding histogram bins
    for(int i=0;i<256;i++) {
        histo_afc_col_sums_x[i] = 0;
        for(int j=0;j<256;j++) {
            histo_afc_col_sums_x[i] = histo_afc_col_sums_x[i] + histo_scatter_autofidcut->GetBinContent(i,j);
            histo_afc_x->SetBinContent(i,histo_afc_col_sums_x[i]);
        }
        
    }
    
    //calculate x bin values for each y row and count # of non-zero bins and set the corresponding histogram bins
    for(int j=0;j<256;j++) {
        histo_afc_col_sums_y[j] = 0;
        for(int i=0;i<256;i++) {
            histo_afc_col_sums_y[j] = histo_afc_col_sums_y[j] + histo_scatter_autofidcut->GetBinContent(i,j);
            histo_afc_y->SetBinContent(j,histo_afc_col_sums_y[j]);
        }
        
    }
    
    //algorithm for finding numeric differences between maxima of the x and y and a certain cut factor (afc_cut_factor)
    
    //get maximums for x and y histogram
    int histo_afc_x_max = (int)histo_afc_x->GetMaximum();
    int histo_afc_y_max = (int)histo_afc_y->GetMaximum();
    
    //set cut factor for pre-fidcut
    float afc_cut_factor = 0.4;
    
    //write maxima of x and y histograms to stdout
    if(verbosity>0) {
    cout << "histo_afc_x->GetMaximum() is:\t" << (int)histo_afc_x->GetMaximum() << "\n";
    cout << "histo_afc_y->GetMaximum() is:\t" << (int)histo_afc_y->GetMaximum() << "\n";
    }
        
    //set x and y fidcut histogram bin contents according to afc_cut_factor and maxima (set to zero if values are less than afc_cut_factor*maximum, keep if values are greater)
    for(int i=0;i<256;i++) {
        if( (histo_afc_x->GetBinContent(i)) < (afc_cut_factor*histo_afc_x_max)) {
            histo_afc_x_cut->SetBinContent(i,0);
        } else {
            histo_afc_x_cut->SetBinContent(i,histo_afc_x->GetBinContent(i));
        }
    }
    
    for(int j=0;j<256;j++) {
        if( (histo_afc_y->GetBinContent(j)) < (afc_cut_factor*histo_afc_y_max)) {
            histo_afc_y_cut->SetBinContent(j,0);
        } else {
            histo_afc_y_cut->SetBinContent(j,histo_afc_y->GetBinContent(j));
        }
        
    }	
    
    
    //TH2F histo_afc_scatter_firstfidcut = new TH2F("histo_afc_scatter_firstfidcut",256,0,255,256,0,255);
    TH2F *histo_afc_scatter_firstfidcut = (TH2F*)histo_scatter_autofidcut->Clone();
    histo_afc_scatter_firstfidcut->SetName("histo_afc_scatter_firstfidcut");
    
    int afc_width_min = 20; //variable for defining minimum block width to be recognized as partial fid cut region
    
    int afc_width_counter = 0; //temporary counter for the width of the current block
    int afc_last_start = 0; //temporary variable to store last start position in block finder algorithm
    int afc_last_end = 0; //temporary variable to store last end position in block finder algorithn
    
    int afc_block_count_x = 0;
    int afc_block_count_y = 0;
    
    //loop for y direction to count projected consecutive blocks with a block width over afc_width_min
    for(int j=1;j<256;j++) {
        if ( (histo_afc_x_cut->GetBinContent(j) > 0) && (histo_afc_x_cut->GetBinContent(j-1) == 0) ) {
            afc_last_start = j;
            afc_width_counter++;
            afc_block_count_x++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 1 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        } else if ( (histo_afc_x_cut->GetBinContent(j) > 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) ) {
            afc_width_counter++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 2 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) && (afc_width_counter < afc_width_min) ) {
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 3 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            for(int k=afc_last_start;k<=j;k++) {
                histo_afc_x_cut->SetBinContent(k,0);
            }
            afc_block_count_x--;
            afc_width_counter = 0;
            continue;
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) && (histo_afc_x_cut->GetBinContent(j-1) > 0) && (afc_width_counter > afc_width_min) ) {
            afc_width_counter = 0;
            afc_block_count_x++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 5 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            
        } else if ( (histo_afc_x_cut->GetBinContent(j) == 0) ) {
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 0 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        }
        else {
            afc_width_counter++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_x_cut->GetBinContent(j) << " is case 4 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        }
    }
    
    //loop for y direction to count projected consecutive blocks with a block width over afc_width_min    
    for(int j=1;j<256;j++) {
        if ( (histo_afc_y_cut->GetBinContent(j) > 0) && (histo_afc_y_cut->GetBinContent(j-1) == 0) ) {
            afc_last_start = j;
            afc_width_counter++;
            afc_block_count_y++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 1 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        } else if ( (histo_afc_y_cut->GetBinContent(j) > 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) ) {
            afc_width_counter++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 2 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) && (afc_width_counter < afc_width_min) ) {
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 3 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            for(int k=afc_last_start;k<=j;k++) {
                histo_afc_y_cut->SetBinContent(k,0);
            }
            afc_block_count_y--;
            afc_width_counter = 0;
            continue;
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) && (histo_afc_y_cut->GetBinContent(j-1) > 0) && (afc_width_counter > afc_width_min) ) {
            afc_width_counter = 0;
            afc_block_count_y++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 5 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            
        } else if ( (histo_afc_y_cut->GetBinContent(j) == 0) ) {
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 0 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
            continue;
        }
        else {
            afc_width_counter++;
            if(verbosity>=2) cout << "column:\t" << j << " has value: " << histo_afc_y_cut->GetBinContent(j) << " is case 4 width afc_last_start= " << afc_last_start << " and widthcounter= " << afc_width_counter << "\n";
        }
    }
    
	
	
    //divide the detector plane in a maximum of four quadrants
    
    if(verbosity>=1) {
    cout << "#-#-#�#�#�#�#�#�#�#�#�#�#�#�#�#" << endl;
    cout << "#" << endl;
    cout << "# afc region recognition started..." << endl;
    cout << "#" << endl;
	}
	
    int afc_x_divider = 0;
    int afc_y_divider = 0;
    bool one_and_only_afc_region_x = 1;
    bool one_and_only_afc_region_y = 1;
    bool one_and_only_afc_region = 0;
	
    if(verbosity>=2) cout << "# x divider function called." << endl;
    for(int i=2;i<256;i++) {			
        if( (histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0) ) {
            for(int j=i;j<256;j++) {
                if ( (histo_afc_x_cut->GetBinContent(j-1) > 0) && (histo_afc_x_cut->GetBinContent(j) == 0) ) {
                    afc_x_divider = j+5;
                    break;
                }
            }
            break;
        }
    }
    
    for(int k=afc_x_divider-4;k<256;k++) {
        if(histo_afc_x_cut->GetBinContent(k)!=0) {
            one_and_only_afc_region_x = 0;
        }
    }
	
    if (one_and_only_afc_region_x) {
        afc_x_divider = 0;
        if(verbosity>=1) cout << "# this run has one and only one afc region in x direction. afc_x_divider set to 0." << endl;
    } else {
        if(verbosity>=1) cout << "# this run has more than one afc region in x direction." << endl;
    }
	
    
    if(verbosity>=1) cout << "#" << endl << "# y divider function called." << endl;
    for(int i=2;i<256;i++) {
        if( (histo_afc_y_cut->GetBinContent(i-1) == 0) && (histo_afc_y_cut->GetBinContent(i) > 0) ) {
            for(int j=i;j<256;j++) {
                if ( (histo_afc_y_cut->GetBinContent(j-1) > 0) && (histo_afc_y_cut->GetBinContent(j) == 0) ) {
                    afc_y_divider = j+5;
                    break;
                }
            }
            break;
        }
    }
    
    for(int k=afc_y_divider-4;k<256;k++) {
        if(histo_afc_y_cut->GetBinContent(k)!=0) {
            one_and_only_afc_region_y = 0;
        }
    }
	
    if (one_and_only_afc_region_y) {
        afc_y_divider = 0;
        if(verbosity>=1) cout << "# this run has one and only one afc region in y direction. afc_y_divider set to 0." << endl;
    } else {
        if(verbosity>=1) cout << "# this run has more than one afc region in y direction." << endl;
    }
	
    if( one_and_only_afc_region_x && one_and_only_afc_region_y ) {
        one_and_only_afc_region = 1;
        if(verbosity>=1) cout << "#" << endl << "#" << endl << "# this run has one and only one afc region." << endl;
    } else {
        if(verbosity>=1) cout << "#" << endl << "#" << endl << "# this run has more than one afc region." << endl;
    }
	
    if(verbosity>=1) cout << "#" << endl << "# afc region recognition ended..." << endl << "#" << endl << "#-#-#�#�#�#�#�#�#�#�#�#�#�#�#�#" << endl;
	
    //end divider code
	
	if(verbosity>=1) { 
    cout << endl;
    cout << "# # # # # # # # # # # # # # # # # # # # # # #" << endl;
    cout << "#" << endl;
    cout << "# Famous afc_x_divider has value: " << afc_x_divider << endl;
    cout << "#" << endl;
    cout << "# Famous afc_y_divider has value: " << afc_y_divider << endl;
    cout << "#" << endl;
    cout << "# # # # # # # # # # # # # # # # # # # # # # #" << endl;	
    cout << endl;
    }
	
    for(int z=1;z<2;z++) {
        for(int u=0;u<afc_x_divider;u++) {
            for(int l=1;l<afc_y_divider;l++) {
                histo_afc_unit_histo_2f->SetBinContent(u,l,1);
                
            }
        }
    }
    
    //		for(int z=1;z<2;z++) {
    //			for(int u=afc_x_divider;u<256;u++) {
    //				for(int l=1;l<256;l++) {
    //					histo_afc_unit_histo_2f->SetBinContent(u,l,0);
    //				
    //				}
    //			}
    //		}
	
	
	
    histo_afc_unit_histo_2f->SetName("histo_afc_unit_histo_2f_r1");
    
    histSaver->SaveHistogramPDF(histo_afc_unit_histo_2f);
	
    (*histo_afc_region_1_mask) = (*histo_afc_unit_histo_2f)*(*histo_afc_scatter_firstfidcut);
	
    histo_afc_region_1_mask->SetName("histo_afc_region_1_mask");
	
    histSaver->SaveHistogramPDF(histo_afc_region_1_mask);
    if(verbosity>=2) cout << endl << endl << "afc_region_1 saved" << endl << endl;
	
    // run HistCleaner on 
    if(verbosity>=2) cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
    HistCleaner(9999, histo_afc_region_1_mask);
    if(verbosity>=2) cout << endl << "Finished HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
	
    //		SaveHistogramPDF(histo_afc_region_1_mask);
    //		cout << endl << endl << "afc_region_1 saved" << endl << endl;
	
    //end divider stuff
	
    // set histogram bin value for histo_afc_scatter_firstfidcut
    for(int i=0;i<256;i++) {
        if(histo_afc_x_cut->GetBinContent(i) == 0) {
            for(int k=0;k<256;k++) {
                histo_afc_scatter_firstfidcut->SetBinContent(i,k,0);
            }
        }
    }
    
    for(int j=0;j<256;j++) {
        if(histo_afc_y_cut->GetBinContent(j) == 0) {
            for(int k=0;k<256;k++) {
                histo_afc_scatter_firstfidcut->SetBinContent(k,j,0);
            }
        }
    }
    
	//define counter variables for region counters
    int afc_region_counter_x = 0;
    int afc_region_counter_y = 0;
	
    //count number of consecutive regions in x direction
    for(int i=2;i<256;i++) {
        if( ((histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0)) || ((histo_afc_x_cut->GetBinContent(i-1) > 0) && (histo_afc_x_cut->GetBinContent(i) == 0))) {	
            afc_region_counter_x++;
        }
    }
	
    //count number of consecutive regions in y direction
    for(int i=2;i<256;i++) {
        if( ((histo_afc_y_cut->GetBinContent(i-1) == 0) && (histo_afc_y_cut->GetBinContent(i) > 0)) || ((histo_afc_y_cut->GetBinContent(i-1) > 0) && (histo_afc_y_cut->GetBinContent(i) == 0))) {	
			afc_region_counter_y++;
        }
    }
	
    if(verbosity>=1) {
    cout << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    cout << "It is most likely that this detector setup has " << (afc_region_counter_x/2)*(afc_region_counter_y/2) << " active diamond regions, i guess." << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    }
	
    int xpos = 2;
    int ypos = 2;
	
    for(int z=1;z<((afc_region_counter_x/2)*(afc_region_counter_y/2))+1;z++) {
        
		for(int i=xpos;i<256;i++) {
			if( (histo_afc_x_cut->GetBinContent(i-1) == 0) && (histo_afc_x_cut->GetBinContent(i) > 0) ) {
				FCR[z]->SetValueXLow(i+5);
				if(verbosity>=2) cout << "Setting FCR[" << z << "].ValueXLow to: " << i+5 << endl;
				for(int j=i+5;j<256;j++) {
					if ( (histo_afc_x_cut->GetBinContent(j-1) > 0) && (histo_afc_x_cut->GetBinContent(j) == 0) ) {
                        FCR[z]->SetValueXHigh(j-5);
                        if(verbosity>=2) cout << "Setting FCR[" << z << "].ValueXHigh to: " << j-5 << endl;					
                        xpos=j+5;
                        break;
					}} break;
                //break;
			}}
        
		
    }
	
    if(verbosity>=2) {
    cout << endl;
    cout << "#####################################################################" << endl;
    cout << endl;
    cout << "FCR[1]->ValueXLow is: " << FCR[1]->GetValueXLow() << endl;
    cout << "FCR[1]->ValueXHigh is: " << FCR[1]->GetValueXHigh() << endl;
    cout << endl;
    cout << "#####################################################################" << endl;
    }
    
    int afc_regions_count = (afc_block_count_x/2)*(afc_block_count_y/2);
    int afc_regions_safety_delta = 5;
    
    if(verbosity>=2) {
    cout << "\nafc_block_count_x is:\t" << afc_block_count_x << "\n";
    cout << "afc_block_count_y is:\t" << afc_block_count_y << "\n";
    cout << "That means, we have " << afc_block_count_x/2 << " block(s) in x- and " << afc_block_count_y/2 << " block(s) in y-direction. so: " << afc_regions_count << " in total.\n\n";
    }
    
    
    //take afc_block_counts_x and afc_block_counts_y and populate FCRs
    
    
    for(int i=0;i<afc_regions_count;i++) {
        FCR[i]->SetAllValuesZero();
        
        
        //		FCR[i]->SetValueXLow(some_var+afc_regions_safety_delta);
        //		FCR[i]->SetValueXHigh(some_var-afc_regions_safety_delta);
        //		FCR[i]->SetValueYLow(some_var+afc_regions_safety_delta);
        //		FCR[i]->SetValueYHigh(some_var-afc_regions_safety_delta);
    }
    
    histo_afc_scatter_firstfidcut->SetTitle("First AutoFidcut based on Scatter Plot");
    histo_afc_scatter_firstfidcut->SetTitleFont(42);
    
    histo_afc_region_1_mask->Reset();
	
    (*histo_afc_region_1_mask) = (*histo_afc_unit_histo_2f)*(*histo_afc_scatter_firstfidcut);
	
    histo_afc_region_1_mask->SetName("histo_afc_region_1_mask_LATEST");
    
    int afc_region_start = 0;
	
    if(afc_y_divider == 0) {
        afc_region_start = 1;
    }
	
	for(int i=1;i<5;i++) {
		
		ostringstream phname;
		phname << "histo_current_2f_" << i;
		
		TH2F* afc_current_histo_2f = new TH2F(phname.str().c_str(),phname.str().c_str(),256,0,255,256,0,255);
		
		histo_afc_unit_histo_2f->Reset();
		if(i==1) {
			for(int a=2;a<afc_x_divider;a++) {
				for(int b=1;b<afc_y_divider;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
        
		if(i==2) {
			for(int a=afc_x_divider;a<256;a++) {
				for(int b=1;b<afc_y_divider;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		if(i==3) {
			for(int a=2;a<afc_x_divider;a++) {
				for(int b=afc_y_divider;b<256;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		if(i==4) {
			for(int a=afc_x_divider;a<256;a++) {
				for(int b=afc_y_divider;b<256;b++) {
					histo_afc_unit_histo_2f->SetBinContent(a,b,1);
				}
			}
		}
		
		
		(*afc_current_histo_2f) = (*histo_afc_scatter_firstfidcut)*(*histo_afc_unit_histo_2f); 
		
		afc_current_histo_2f->SetName(phname.str().c_str());
		
        if(verbosity>=1) {
		cout << "TENGTARATENG" << endl;
		cout << "FOR LOOP FOR FIDCUT REGIONS. NAME OF CURRENT HISTO IS: " << afc_current_histo_2f->GetName() << "." << endl;
        }
		HistCleaner(i, afc_current_histo_2f);
		
		
		histSaver->SaveHistogramPDF(afc_current_histo_2f);
		
		afc_current_histo_2f->Delete();
		
	}
	
	
    
    //		cout << endl << "Calling HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;
    HistCleaner(1, histo_afc_region_1_mask);
    //		 cout << endl << "Finished HistCleaner on " << histo_afc_region_1_mask->GetName() << "." << endl;  
	
	
    histSaver->SaveHistogramPDF(histo_afc_region_1_mask);
    
    HistCleaner(1, histo_afc_scatter_firstfidcut);
	
    FCR[0]->GetAllValues();
    
    FCR[3]->SetAllValuesZero();
    FCR[3]->SetValueXLow(1);
    FCR[3]->SetValueYHigh(2);
    FCR[3]->GetAllValues();
    
    histSaver->SaveHistogramPDF(histo_afc_x);
    histSaver->SaveHistogramPDF(histo_afc_y);
    histSaver->SaveHistogramPDF(histo_afc_x_cut);
    histSaver->SaveHistogramPDF(histo_afc_y_cut);
    histSaver->SaveHistogramPDF(histo_afc_scatter_firstfidcut);
    
    
    histSaver->SaveHistogramPNG(histo_afc_x);
    histSaver->SaveHistogramPNG(histo_afc_y);
    histSaver->SaveHistogramPNG(histo_afc_x_cut);
    histSaver->SaveHistogramPNG(histo_afc_y_cut);
    histSaver->SaveHistogramPNG(histo_afc_scatter_firstfidcut);
    
    if(verbosity>=0) {
    cout << "\n";
    cout << "FINISHED AutoFidCut \n\n\n";
    }
}	


void Clustering::Alignment(bool plots, bool CutFakeTracksOn){
	//Align(plots, CutFakeTracksOn);

	alignment->GetEvent(current_event);
	alignment->SetTracks(tracks);
	alignment->SetAlignment_tracks_fidcut(tracks_fidcut);
	alignment->SetAlignment_tracks_fidcut_mask(tracks_fidcut_mask);
	if( alignment->Align(0,CutFakeTracksOn)==-1){
		if(verbosity) cout<<"Clustering::Alignment::No Clusters build yet, do ClusterRun"<<endl;
		ClusterRun(0);
		if(verbosity) cout<<"Clustering::Alignment::track size: "<<tracks.size()<<endl;
		alignment->SetTracks(tracks);
		alignment->SetAlignment_tracks_mask(tracks_mask);
		alignment->SetAlignment_tracks_fidcut(tracks_fidcut);
		alignment->SetAlignment_tracks_fidcut_mask(tracks_fidcut_mask);
		if(alignment->Align(0,CutFakeTracksOn)==-1) cout<<"clustering::Alignment:Clustering did not worked"<<endl;
		exit (-1);
	}
}

void Clustering::Align(bool plots, bool CutFakeTracksOn) {
	
	if(tracks.size()==0) {
		cout<<"Clustering::Align: No tracks found; calling Clustering::ClusterRun first..."<<endl;
		ClusterRun(plots); // doesn't use alternative clustering
	}
	
	if (tracks.size() == 0) {
		cout << "No tracks available. Alignment not possible. (tracks.size() = " << tracks.size() << ")" << endl;
		return;
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
	
	string plots_path_save = plots_path;
	vector<TDiamondTrack> alignment_tracks = tracks;
	vector<bool> alignment_tracks_mask = tracks_mask;
	vector<TDiamondTrack> alignment_tracks_fidcut = tracks_fidcut;
	vector<bool> alignment_tracks_fidcut_mask = tracks_fidcut_mask;
	ostringstream plots_path_alignment, plots_path_alignment_CutFakeTracks;
	plots_path_alignment << plots_path << "/alignment/";
	plots_path_alignment_CutFakeTracks << plots_path << "/alignment_CutFakeTracks/";
	plots_path = plots_path_alignment.str();
	
	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
	for (int alignStep = 0; alignStep < 2; alignStep++) {
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
			align->CutFakeTracks(alignment_tracks, alignment_tracks_mask, settings->getAlignment_chi2(), CutFakeTracksOn, true);
			align->CutFakeTracks(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, settings->getAlignment_chi2(), CutFakeTracksOn, true);
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
			TransparentClustering(alignment_tracks_fidcut, alignment_tracks_fidcut_mask, align);
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
	plots_path = plots_path_save;
}

// -- transparent analysis for any plane
void Clustering::TransparentAnalysis(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose) {
    cout << "Clustering::TransparentAnalysis:Starting transparent clustering with " << tracks.size() << " tracks.." << endl;
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
	
	cout << "init histograms for transparent clustering.." << endl;
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
		//PedTree->GetEvent(j);
		eventReader->GetEvent(j);
        //		cout << endl << endl << endl;
        //		cout << "event " << j << " in PedTree has event_number: " << event_number << endl;
		event_numbers.push_back(eventReader->event_number);
        //		for (int blablabla = 0; blablabla < 254; blablabla++) {
        //			cout << "EventReader->Dia_ADC = " << EventReader->Dia_ADC[blablabla] << ",\tEventReader->Det_PedMean = " << EventReader->Det_PedMean[8][blablabla] << endl;
        //			cout << "Collected charge in channel " << blablabla << " of diamond: " << EventReader->Dia_ADC[blablabla]-EventReader->Det_PedMean[8][blablabla] << endl;
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
		eventReader->GetEvent(event);
		if (verbose) cout << "event_number = " << eventReader->event_number << endl;
		
		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int j = 0; j < 128; j++) {
			if (eventReader->Det_ADC[6][j]-eventReader->Det_PedMean[6][j] > (eventReader->Det_ADC[6][eff_diamond_hit_channel]-eventReader->Det_PedMean[6][eff_diamond_hit_channel])) {
				eff_diamond_hit_channel = j;
			}
		}
		histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
        histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position, diamond_hit_y_position);
		cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;
		
		// cluster diamond channels around estimated hit position
        //		cout << " --" << endl;
		if (verbose) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbose) cout << "EventReader->Det_ADC[6] = " << eventReader->Det_ADC[6][diamond_hit_channel] << ",\tEventReader->Det_PedMean = " << eventReader->Det_PedMean[6][diamond_hit_channel] << endl;
		if (verbose) cout << "Collected charge in channel " << (int)diamond_hit_position << " of D3X: " << eventReader->Det_ADC[6][(int)diamond_hit_position]-eventReader->Det_PedMean[6][(int)diamond_hit_position] << endl;
		
		int sign;
		
		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;
		
		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = eventReader->Det_ADC[6][diamond_hit_channel]-eventReader->Det_PedMean[6][diamond_hit_channel];
		secondchannel_adc = eventReader->Det_ADC[6][diamond_secondhit_channel]-eventReader->Det_PedMean[6][diamond_secondhit_channel];
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
				if (current_channel > 0 && current_channel < 128 /* && EventReader->Dia_ADC[current_channel]-EventReader->Det_PedMean[8][current_channel] > settings->getDi_Cluster_Hit_Factor()*EventReader->Det_PedWidth[8][current_channel]*/)
					cluster_adc = cluster_adc + eventReader->Det_ADC[6][current_channel]-eventReader->Det_PedMean[6][current_channel];
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

// -- take track after alignment, point into diamond and produce a cluster with the surrounding channels
void Clustering::TransparentClustering(vector<TDiamondTrack> &tracks, vector<bool> &tracks_mask, TDetectorAlignment *align, bool verbose) {
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
	cout << "diamond hit factor: " << settings->getDi_Cluster_Hit_Factor() << endl;
	
	cout << "init histograms for transparent clustering.." << endl;
	for (int i = 0; i < 10; i++) {
		ostringstream histoname_landau, histoname_eta;
		histoname_landau << "PulseHeight_Dia_" << (i+1) << "ChannelsTransparAna_8HitsFidcut";
		cout << "histoname_landau: " << histoname_landau.str().c_str() << endl;
		histo_transparentclustering_landau[i] = new TH1F(histoname_landau.str().c_str(),histoname_landau.str().c_str(),settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
//		histoname_eta << "Eta_Dia_" << (i+1) << "HitTransparClusters";
//		cout << "histoname_eta: " << histoname_eta.str().c_str() << endl;
	}
    histo_transparentclustering_landau_mean = new TH1F("PulseHeightMeanVsChannels_Dia_TranspAna_","PulseHeightMeanVsChannels_Dia_TranspAna",10,0.5,10.5);
	ostringstream histoname_eta;
	histoname_eta << "Eta_Dia_2CentroidHits_TransparClusters";
	histo_transparentclustering_eta = new TH1F(histoname_eta.str().c_str(),histoname_eta.str().c_str(),100,0.,1.);
	histo_transparentclustering_hitdiff = new TH1F("DiffEstEffHit_Dia_TransparClusters","DiffEstEffHit_Dia_TransparClusters", 200, -5.,5.);
	histo_transparentclustering_hitdiff_scatter = new TH2F("DiffEstEffHit_Scatter_Dia_TransparClusters","DiffEstEffHit_Scatter_Dia_TransparClusters", 200, -5.,5.,128,0,127);
    histo_transparentclustering_2Channel_PulseHeight = new TH1F("PulseHeight_Dia_2Channel_TranspCluster_8HitsFidcut","PulseHeight_Dia_2Channel_TranspCluster_8HitsFidcut",settings->getPulse_height_num_bins(),-0.5,settings->getPulse_height_di_max()+0.5);
	cout << " done." << endl;
	
//	verbose = true;
	
	
	eventReader->GetEvent(0);
	cout<< "Loaded first event in EventReader: "<<eventReader->event_number<<endl;
	cout<< "RunNumber is: "<<eventReader->run_number<<endl;
	cout<< "StoreThreshold is: "<<eventReader->store_threshold<<endl;
	
	event_numbers.clear();
	for (int j = 0; j < eventReader->GetEntries(); j++) {
		eventReader->GetEvent(j);
//		cout << endl << endl << endl;
//		cout << "event " << j << " in PedTree has event_number: " << event_number << endl;
		event_numbers.push_back(eventReader->event_number);
//		for (int blablabla = 0; blablabla < 254; blablabla++) {
//			cout << "EventReader->Dia_ADC = " << EventReader->Dia_ADC[blablabla] << ",\tEventReader->Det_PedMean = " << EventReader->Det_PedMean[8][blablabla] << endl;
//			cout << "Collected charge in channel " << blablabla << " of diamond: " << EventReader->Dia_ADC[blablabla]-EventReader->Det_PedMean[8][blablabla] << endl;
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
		eventReader->GetEvent(event);
		if (verbose) cout << "event_number = " << eventReader->event_number << endl;
		
		// find biggest hit in diamond
		eff_diamond_hit_channel = 0;
		for (int j = 0; j < 128; j++) {
			if (eventReader->Dia_ADC[j]-eventReader->Det_PedMean[8][j] > (eventReader->Dia_ADC[eff_diamond_hit_channel]-eventReader->Det_PedMean[8][eff_diamond_hit_channel])) {
				eff_diamond_hit_channel = j;
			}
		}
		histo_transparentclustering_hitdiff->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position); // added 0.5 to eff_diamond_hit_channel to take the middle of the channel instead of the edge
        histo_transparentclustering_hitdiff_scatter->Fill(eff_diamond_hit_channel + 0.5 - diamond_hit_position, diamond_hit_y_position);
		cout << "effective diamond hit channel: " << eff_diamond_hit_channel << endl;
		
		// cluster diamond channels around estimated hit position
//		cout << " --" << endl;
		if (verbose) cout << "track " << i << " has an estimated hit position at " << diamond_hit_position << " (channel " << diamond_hit_channel << ")" << endl;
		if (verbose) cout << "EventReader->Dia_ADC = " << eventReader->Dia_ADC[diamond_hit_channel] << ",\tEventReader->Det_PedMean = " << eventReader->Det_PedMean[8][diamond_hit_channel] << endl;
		if (verbose) cout << "Collected charge in channel " << (int)diamond_hit_position << " of diamond: " << eventReader->Dia_ADC[(int)diamond_hit_position]-eventReader->Det_PedMean[8][(int)diamond_hit_position] << endl;
		
		int sign;
		
		if (diamond_hit_position - diamond_hit_channel < 0.5) sign = 1;
		else sign = -1;
		
		// calculate eta for the two closest two channels to the estimated hit position
		diamond_secondhit_channel = diamond_hit_channel - sign;
		firstchannel_adc = eventReader->Dia_ADC[diamond_hit_channel]-eventReader->Det_PedMean[8][diamond_hit_channel];
		secondchannel_adc = eventReader->Dia_ADC[diamond_secondhit_channel]-eventReader->Det_PedMean[8][diamond_secondhit_channel];
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
				if (current_channel > 0 && current_channel < 128 /* && EventReader->Dia_ADC[current_channel]-EventReader->Det_PedMean[8][current_channel] > settings->getDi_Cluster_Hit_Factor()*EventReader->Det_PedWidth[8][current_channel]*/)
					cluster_adc = cluster_adc + eventReader->Dia_ADC[current_channel]-eventReader->Det_PedMean[8][current_channel];
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

// --
void Clustering::EventMonitor(int CurrentEvent) {
	current_event = current_event - 1;
	
	//record large clusters
	string large_clusters_filename = plots_path + "large_clusters.txt";
	ofstream large_clusters_file(large_clusters_filename.c_str());
	large_clusters_file << "event\tdet\tclus\tnhits\tnseeds" << endl;
	uint nclusters, nhits;
	Cluster* current_cluster;
	
	//histograms of hit positions to find out if the eta correction helped
	TH1F* histo_hitpos[9][2];
	for(int det=0; det<9; det++) {
		for(int chip=0; chip<2; chip++) {
			ostringstream histotitle;
			histotitle << "Hit_Position_Interstrip_D" << det << "_chip" << chip;
			histo_hitpos[det][chip] = new TH1F(histotitle.str().c_str(), histotitle.str().c_str(), 101, -0.005, 1.005);
		}
	}
	
	if (!settings->getAlternativeClustering()) ClusterEvent();
	else ClusterEventSeeds();
//	if(e%10000==0) clustered_event.Print();
//	BookHistograms();
	
	//record large clusters
	//if a track, check whether it's in the silicon fiducial region
	if(clustered_event.HasOneSiliconTrack()) {
		
		//loop over detectors to record large clusters to a txt file
		for(int det=0; det<9; det++) {
            nclusters = clustered_event.GetNClusters(det);
            //loop over all *good* clusters
            for(uint clus=0; clus<nclusters; clus++) {
				current_cluster = clustered_event.GetCluster(det,clus);
				nhits = current_cluster->GetNHits();
				if(nhits>8) {
					large_clusters_file << eventReader->event_number << "\t" << det << "\t" << clus << "/" << nclusters << "\t" << nhits << "\t" << current_cluster->GetNSeeds();
					if(current_cluster->IsBadChannelCluster()) large_clusters_file << " (bad chan cluster)"; // skip bad clusters
					large_clusters_file << endl;
				}
            }
		}//end large cluster record
		
	}
	
	ostringstream histo_diamond_name;
	histo_diamond_name << "Event" << current_event << "DiaClusters";
	TH2F* histo_clusters[5];
	TH1F* histo_diamond;
	histo_diamond = new TH1F(histo_diamond_name.str().c_str(),histo_diamond_name.str().c_str(),256,-0.5,255.5);
	Float_t si_avg_x=0, si_avg_y=0;
	
	for (int det = 0; det < 5; det++) {
		ostringstream histo_clusters_name;
		histo_clusters_name << "Event" << current_event << "Det" << det << "Clusters";
		histo_clusters[det] = new TH2F(histo_clusters_name.str().c_str(),histo_clusters_name.str().c_str(),256,-0.5,255.5,256,-0.5,255.5);
		Float_t si_x=0, si_y=0;
		Float_t cluster_size = 0;
		//		int bin_x = 0, bin_y = 0, bin = 0;
		si_x = clustered_event.GetCluster(2*det, 0)->Get1stMoment();
		if (det < 4) {
			si_y = clustered_event.GetCluster(2*det+1, 0)->Get1stMoment();
			cluster_size = (clustered_event.GetCluster(2*det, 0)->GetNHits() + clustered_event.GetCluster(2*det+1, 0)->GetNHits()) / 2.;
			si_avg_x += si_x;
			si_avg_y += si_y;
		}
//		bin_x = histo_clusters[det]->
//		bin = histo_clusters[det]->GetBin(
		if (det < 4) {
			histo_clusters[det]->SetBinContent(si_x+1,si_y+1,cluster_size);
		}
		else {
			si_avg_x = si_avg_x / 4;
			si_avg_y = si_avg_y / 4;
			histo_clusters[det]->SetBinContent(si_x+1,si_avg_y+1,clustered_event.GetCluster(2*det, 0)->GetNHits());
			
			for (int j = 0; j < eventReader->Det_NChannels[8]; j++) {
				if (eventReader->Dia_ADC[j]-eventReader->Det_PedMean[8][j] > settings->getDi_Cluster_Hit_Factor()*eventReader->Det_PedWidth[8][j] || true) {
//					histo_diamond->SetBinContent(j+1,EventReader->Dia_ADC[j]-EventReader->Det_PedMean[8][j]);
					histo_diamond->SetBinContent((int)eventReader->Det_Channels[8][j],eventReader->Dia_ADC[j]-eventReader->Det_PedMean[8][j]);
				}
			}
			histSaver->SaveHistogram(histo_diamond);
//			histo_diamond->SetBinContent(si_x+1,clustered_event.GetCluster(2*det, 0)->GetNHits());
		}

		histSaver->SaveHistogram(histo_clusters[det]);
	}
	
	current_event++;
}

void Clustering::SetRunParameters(int reg, FidCutRegion current_region, bool MultiRegions) {
	settings->setSi_avg_fidcut_xlow(current_region.GetValueXLow());
	settings->setSi_avg_fidcut_xhigh(current_region.GetValueXHigh());
	settings->setSi_avg_fidcut_ylow(current_region.GetValueYLow());
	settings->setSi_avg_fidcut_yhigh(current_region.GetValueYHigh());
	
	// TODO: alternative plots path only for more than 1 region!!
	
	if (MultiRegions) {
		
		ostringstream plotspath;
		plotspath << plots_path.c_str() << "-FidCutRegion" << reg;
		
		//	plotspath << sys->pwd() << "/plots-" << RunNumber;
		//	if(RunDescription=="") plotspath << "/";
		//	else plotspath << "-" << RunDescription << "/";
		
		plotspath << "-FidCutRegion" << reg;
		
		plots_path = plotspath.str();
		png_file_char = plotspath.str();
		C_file_char = plotspath.str();
		root_file_char = plotspath.str();
		//make plots dir
		sys->mkdir(plots_path.c_str());
	}
}
