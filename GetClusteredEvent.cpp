//GetClusteredEvent
//2010-09-16 Fixed event lookup

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include "Cluster.class.cpp"
#include "ClusteredEvent.class.cpp"
#include "ChannelScreen.h"


void ParseIntArray(string value, vector<int> &vec) {
   int index=0;
   string::size_type offset1 = value.find_first_of('{');
   string::size_type offset2 = value.find_first_of(',');
   vec.push_back((int)strtod(value.substr(offset1,offset2).c_str(),0));
   value = value.substr(offset2+1,value.length()-(offset2+1));
   while(value.length()>2) {
      offset2 = TMath::Min(value.find_first_of(','),value.find_first_of('}'));
      vec.push_back((int)strtod(value.substr(0,offset2).c_str(),0));
      index++;
      //cout<<"vec["<<index<<"]="<<vec[index]<<"\tvalue.length()="<<value.length()<<"\tvalue.substr(0,offset2)="<<value.substr(0,offset2)<<"\tstrtod(value.substr(0,offset2),0)="<<strtod(value.substr(0,offset2).c_str(),0)<<endl;//check
      if(value.find_first_of(',')>value.find_first_of('}')) break;
      value = value.substr(offset2+1,value.length()-(offset2+1));
   }
}

//get current event and cluster
ClusteredEvent GetClusteredEvent(unsigned int RunNumber, unsigned int EventNumber, string RunDescription = "", bool verbose = 0) {
   
   //settings
   Int_t SaveAllFilesSwitch; //1 for save files, 0 for don't
   Int_t ClosePlotsOnSave;
   Int_t IndexProduceSwitch;

   Float_t Si_Cluster_Seed_Factor;
   Float_t Si_Cluster_Hit_Factor;
   Float_t Di_Cluster_Seed_Factor;
   Float_t Di_Cluster_Hit_Factor;
      
   vector<int> single_channel_analysis_channels;
      
   //Channels to Screen
   vector<int> Det_channel_screen_input[9];
   ChannelScreen Det_channel_screen[9];
      
   //paths
   string png_file_char;
   string C_file_char;
   string root_file_char;
   TSystem* sys;
   string settings_file;
      
   //processed event storage; read from pedtree
   UInt_t run_number;
   UInt_t event_number;
   Float_t store_threshold;
   UInt_t Det_NChannels[9];
   UChar_t Det_Channels[9][256];
   UChar_t Det_ADC[8][256];
   UShort_t Dia_ADC[256];
   Float_t Det_PedMean[9][256];
   Float_t Det_PedWidth[9][256];
      
   //output clusters
   ClusteredEvent* clustered_event;
      
   //io
   UInt_t current_event;
   TFile *PedFile;
   TTree *PedTree;
   
   Si_Cluster_Seed_Factor = 5;
   Si_Cluster_Hit_Factor = 3;
   Di_Cluster_Seed_Factor = 5;
   Di_Cluster_Hit_Factor = 3;
   
   //default paths
   sys = gSystem;
   
   ostringstream settingspath;
   settingspath << sys->pwd() << "/Settings." << RunNumber;
   if(RunDescription=="") settingspath << ".ini";
   else settingspath << "-" << RunDescription << ".ini";
   settings_file = settingspath.str();
   
   ostringstream pedfilepath;
   pedfilepath << sys->pwd() << "/Pedestal." << RunNumber;
   if(RunDescription=="") pedfilepath << ".root";
   else pedfilepath << "-" << RunDescription << ".root";
   
   //---------------------------------
   //Load settings
   
   cout<<endl<<"Overriding default settings with settings in Settings.ini"<<endl<<endl;
   
   ifstream file(settings_file.c_str());
   if(!file) cout << "An error has encountered while trying to open file " << settings_file << endl;
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
      if(key=="D0X_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[0]);
      }
      if(key=="D0Y_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[1]);
      }
      if(key=="D1X_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[2]);
      }
      if(key=="D1Y_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[3]);
      }
      if(key=="D2X_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[4]);
      }
      if(key=="D2Y_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[5]);
      }
      if(key=="D3X_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[6]);
      }
      if(key=="D3Y_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[7]);
      }
      if(key=="Dia_channel_screen_input") {
         cout << key.c_str() << " = " << value.c_str() << endl;
         ParseIntArray(value,Det_channel_screen_input[8]);
      }
   }
   
   file.close();
   cout<<endl<<"Finished importing settings from Settings.ini"<<endl<<endl;
   
   //Load settings
   //---------------------------------
   
   //screen channels
   for(int det=0; det<9; det++) 
      Det_channel_screen[det].ScreenChannels(Det_channel_screen_input[det]);
   
   
   //Get pedestal subtracted event
   
   PedFile = new TFile(pedfilepath.str().c_str());
   PedTree = (TTree*)PedFile->Get("PedTree");
   if (!PedTree)
   {
      cerr << "PedTree not found!" << endl;
   }
      
      
   //Event Header Branches
   PedTree->SetBranchAddress("RunNumber",&run_number);
   PedTree->SetBranchAddress("EventNumber",&event_number);
   PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
   
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
      
   
   //print stats
   cout << PedTree->GetEntries() << " events in PedTree " << endl;
   PedTree->GetEvent(0);
   cout<<"First event in PedTree is "<<event_number<<endl;
   PedTree->GetEvent(PedTree->GetEntries()-1);
   cout<<"Last event in PedTree is "<<event_number<<endl;
     
   //try to locate event
   PedTree->GetEvent(0);
   for(int guessed_event = EventNumber-event_number; guessed_event>0; guessed_event--) {
      PedTree->GetEvent(guessed_event);
      //cout<<"event_number="<<event_number<<endl;
      if(event_number==EventNumber) break;
   }
   cout<<"Event "<<event_number<<" loaded"<<endl;
      
   PedFile->Close();
   
   
   
   //create new clustered event
   clustered_event = new ClusteredEvent(event_number, Si_Cluster_Hit_Factor, Si_Cluster_Seed_Factor, Di_Cluster_Hit_Factor, Di_Cluster_Seed_Factor);
   
   vector<int> hits, cluster, badchannelclusterflags, goldengateclusterflags;
   vector< vector<int> > clusters;
   int previouschan, currentchan, hasaseed, hasmasked, previousseed, isgoldengate;
   
   for(int det=0; det<9; det++) {
      hits.clear();
      cluster.clear();
      clusters.clear();
      
      //look for hits
      if(verbose) cout<<endl<<endl<<"Detector "<<det<<" hits: ";
      for(int i=0; i<(int)Det_NChannels[det]; i++) {
         if(det<8) if(Det_ADC[det][i]-Det_PedMean[det][i] > Si_Cluster_Hit_Factor*Det_PedWidth[det][i]) {
            hits.push_back(i);
            if(verbose) {
               cout<<(int)Det_Channels[det][i];
               if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
               cout<<", ";
            }
         }
         if(det==8) if(Dia_ADC[i]-Det_PedMean[det][i] > Di_Cluster_Hit_Factor*Det_PedWidth[det][i]) {
            hits.push_back(i);
            if(verbose) {
               cout<<(int)Det_Channels[det][i];
               if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][i])) cout<<"(masked)";
               cout<<", ";
            }
         }
      }
      if(verbose) {
         cout<<endl<<"Channels screened: ";
         for(int i=0; i<256; i++) if(!Det_channel_screen[det].CheckChannel(i)) cout<<i<<", ";
         cout<<endl<<hits.size()<<" hits found in "<<(int)Det_NChannels[det]<<" saved channels"<<endl;
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
      if(verbose) cout<<"after clear(): cluster.size()="<<cluster.size()<<" and cluster.size()="<<cluster.size()<<endl;
      previouschan=-1;
      for(uint i=0; i<hits.size(); i++) {
         currentchan = Det_Channels[det][hits[i]];
         if(verbose) {
            if(hits[i]==-1) cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or end of hits"<<endl;
            else {
               cout<<"examining hit "<<i<<" at channel index "<<hits[i]<<" or channel "<<currentchan<<" (";
               if(det==8) cout<<Dia_ADC[hits[i]]-Det_PedMean[det][hits[i]]<<" psadc)"<<endl;
               if(det<8) cout<<Det_ADC[det][hits[i]]-Det_PedMean[det][hits[i]]<<" psadc)"<<endl;
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
            isgoldengate=0;
            previousseed=-1;
            for(uint j=0; j<cluster.size(); j++) {
               currentchan = cluster[j];
               if(verbose) cout<<"cluster["<<j<<"]="<<cluster[j]<<" or channel "<<(int)Det_Channels[det][currentchan];
               if(det<8) if(Det_ADC[det][currentchan]-Det_PedMean[det][currentchan] > Si_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
                  hasaseed=1; 
                  if(verbose) cout<<" is a seed";
                  if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
                     isgoldengate=1;
                     cout<<" (goldengate cluster)";
                  }
                  else previousseed=Det_Channels[det][currentchan];
               }
               if(det==8) if(Dia_ADC[currentchan]-Det_PedMean[det][currentchan] > Si_Cluster_Seed_Factor*Det_PedWidth[det][currentchan]) {
                  hasaseed=1; 
                  if(verbose) cout<<" is a seed";
                  if(previousseed!=-1 && Det_Channels[det][currentchan]!=previousseed+1) {
                     isgoldengate=1;
                     cout<<" (goldengate cluster)";
                  }
                  else previousseed=Det_Channels[det][currentchan];
               }
               if(!Det_channel_screen[det].CheckChannel((int)Det_Channels[det][currentchan])) {hasmasked=1; if(verbose) cout<< " (masked channel)";}
               if(verbose) cout<<endl;
            }
            if(hasaseed) {
               clusters.push_back(cluster); 
               badchannelclusterflags.push_back(hasmasked); 
               goldengateclusterflags.push_back(isgoldengate); 
               if(verbose) {
                  cout<<"storing cluster"<<endl;
                  if(hasmasked) cout<<"flagging BadChannelCluster"<<endl;
                  if(isgoldengate) cout<<"flagging GoldenGateCluster"<<endl;
               }
            } //if there's a seed in the cluster, save it
            cluster.clear(); //start storing a new cluster
            previouschan=-1; //save next channel
            //move on to the next cluster
         }
      }
      
      if(verbose) if(clusters.size()==0)  {
         cout<<"No clusters found so skipping to next detector."<<endl;
         continue;
      }
      
      //now that we have lists of channels belonging to clusters, let's create Cluster objects to store the data
      Cluster* current_cluster = 0;
      if(clusters.size()>0) {
         for(uint i=0; i<clusters.size(); i++) {
            current_cluster = clustered_event->AddCluster(det);
            current_cluster->FlagBadChannelCluster(badchannelclusterflags[i]);
            current_cluster->FlagGoldenGateCluster(goldengateclusterflags[i]);
            for(uint j=0; j<clusters[i].size(); j++) { //save each cluster one channel at a time
               currentchan = clusters[i][j];
               if(det<8) {
                  current_cluster->AddHit(Det_Channels[det][currentchan], Det_ADC[det][currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
               }
               if(det==8) {
                  current_cluster->AddHit(Det_Channels[det][currentchan], Dia_ADC[currentchan], Det_PedMean[det][currentchan], Det_PedWidth[det][currentchan]);
               }
               if(verbose) cout<<"Detector "<<det<<": cluster_saved/all_clusters="<<i+1<<"/"<<clusters.size()<<"\tchan_to_save="<<j+1<<"/"<<clusters[i].size()<<"\tcurrentchan="<<(int)Det_Channels[det][currentchan]<<endl;
            }//end loop over channels in a cluster to save
         }//end loop over clusters
      }
   }//end loop over detectors
   
   return *clustered_event;
}


