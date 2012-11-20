//Separate data reduction from data processing
//2010-07-27 Start

#include "ClusterAnalyze.hh"

ClusterAnalyze::ClusterAnalyze(unsigned int RunNumber, string RunDescription = "") {
   run_number = RunNumber;
   
   //gROOT->ProcessLine("#include <vector>"); // so root can store vectors to trees
   //gROOT->ProcessLine("#include <Cluster.class.cpp>"); // so root can store vectors to trees
   //gROOT->ProcessLine(".L Cluster.class.cpp+"); // so root can store vectors to trees
   
   //default settings
   SaveAllFilesSwitch = 1; //1 for save files, 0 for don't
   ClosePlotsOnSave = 1;
   IndexProduceSwitch = 1;
   
   Si_Cluster_Seed_Factor = 5;
   Si_Cluster_Hit_Factor = 3;
   Di_Cluster_Seed_Factor = 5;
   Di_Cluster_Hit_Factor = 3;
   

   //default paths
   sys = gSystem;
   
   ostringstream plotspath;
   plotspath << sys->pwd() << "/plots-" << RunNumber;
   if(RunDescription=="") plotspath << "/";
   else plotspath << "-" << RunDescription << "/";
   png_file_char = plotspath.str();
   C_file_char = plotspath.str();
   root_file_char = plotspath.str();
   
   ostringstream settingspath;
   settingspath << sys->pwd() << "/Settings." << RunNumber;
   if(RunDescription=="") settingspath << ".ini";
   else settingspath << "-" << RunDescription << ".ini";
   settings_file = settingspath.str();
   
   ostringstream clusterfilepath;
   clusterfilepath << sys->pwd() << "/Cluster." << RunNumber;
   if(RunDescription=="") clusterfilepath << ".root";
   else clusterfilepath << "-" << RunDescription << ".root";
   
   
   LoadSettings();
   
   
   //store output pedestal subtracted data
   ClusterFile = new TFile(clusterfilepath.str().c_str(),"recreate");
   ClusterTree = new TTree("ClusterTree","A Tree for Clustered Data");
   //PedTree->SetMaxTreeSize(8000000000); 
   //PedTree->SetMaxTreeSize(Long64_t(TMath::Power(2,63))); // http://lists.healthgrid.org/archives/gate-users/2004-December/000396.html
   PedTree->SetMaxTreeSize(Long64_t(TMath::Power(2,45)));  // Gives ~1TB (actually 4TB) max file size
   //PedTree->Branch("EventBranch","PSEvent",&store_event,24000,0); //address must be address of a pointer to a stored event
   
   //Event Header Branches
   ClusterTree->SetBranchAddress("RunNumber",&run_number,"RunNumber/i");
   ClusterTree->SetBranchAddress("EventNumber",&event_number,"EventNumber/i");
   ClusterTree->SetBranchAddress("StoreThreshold",&store_threshold,"StoreThreshold/F");

   //Telescope Data Branches
   ClusterTree->SetBranchAddress("D0X_NClusters",&Det_NClusters[0]);
   ClusterTree->SetBranchAddress("D0Y_NClusters",&Det_NClusters[1]);
   ClusterTree->SetBranchAddress("D1X_NClusters",&Det_NClusters[2]);
   ClusterTree->SetBranchAddress("D1Y_NClusters",&Det_NClusters[3]);
   ClusterTree->SetBranchAddress("D2X_NClusters",&Det_NClusters[4]);
   ClusterTree->SetBranchAddress("D2Y_NClusters",&Det_NClusters[5]);
   ClusterTree->SetBranchAddress("D3X_NClusters",&Det_NClusters[6]);
   ClusterTree->SetBranchAddress("D3Y_NClusters",&Det_NClusters[7]);
   ClusterTree->SetBranchAddress("Dia_NClusters",&Det_NClusters[8]);
   ClusterTree->SetBranchAddress("D0X_Cluster_NChannels",&Det_Cluster_NChannels[0]);
   ClusterTree->SetBranchAddress("D0Y_Cluster_NChannels",&Det_Cluster_NChannels[1]);
   ClusterTree->SetBranchAddress("D1X_Cluster_NChannels",&Det_Cluster_NChannels[2]);
   ClusterTree->SetBranchAddress("D1Y_Cluster_NChannels",&Det_Cluster_NChannels[3]);
   ClusterTree->SetBranchAddress("D2X_Cluster_NChannels",&Det_Cluster_NChannels[4]);
   ClusterTree->SetBranchAddress("D2Y_Cluster_NChannels",&Det_Cluster_NChannels[5]);
   ClusterTree->SetBranchAddress("D3X_Cluster_NChannels",&Det_Cluster_NChannels[6]);
   ClusterTree->SetBranchAddress("D3Y_Cluster_NChannels",&Det_Cluster_NChannels[7]);
   ClusterTree->SetBranchAddress("Dia_Cluster_NChannels",&Det_Cluster_NChannels[8]);
   ClusterTree->SetBranchAddress("D0X_Cluster_Channels",&Det_Cluster_Channels[0]);
   ClusterTree->SetBranchAddress("D0Y_Cluster_Channels",&Det_Cluster_Channels[1]);
   ClusterTree->SetBranchAddress("D1X_Cluster_Channels",&Det_Cluster_Channels[2]);
   ClusterTree->SetBranchAddress("D1Y_Cluster_Channels",&Det_Cluster_Channels[3]);
   ClusterTree->SetBranchAddress("D2X_Cluster_Channels",&Det_Cluster_Channels[4]);
   ClusterTree->SetBranchAddress("D2Y_Cluster_Channels",&Det_Cluster_Channels[5]);
   ClusterTree->SetBranchAddress("D3X_Cluster_Channels",&Det_Cluster_Channels[6]);
   ClusterTree->SetBranchAddress("D3Y_Cluster_Channels",&Det_Cluster_Channels[7]);
   ClusterTree->SetBranchAddress("Dia_Cluster_Channels",&Det_Cluster_Channels[8]);
   ClusterTree->SetBranchAddress("D0X_Cluster_ADC",&Det_Cluster_ADC[0]);
   ClusterTree->SetBranchAddress("D0Y_Cluster_ADC",&Det_Cluster_ADC[1]);
   ClusterTree->SetBranchAddress("D1X_Cluster_ADC",&Det_Cluster_ADC[2]);
   ClusterTree->SetBranchAddress("D1Y_Cluster_ADC",&Det_Cluster_ADC[3]);
   ClusterTree->SetBranchAddress("D2X_Cluster_ADC",&Det_Cluster_ADC[4]);
   ClusterTree->SetBranchAddress("D2Y_Cluster_ADC",&Det_Cluster_ADC[5]);
   ClusterTree->SetBranchAddress("D3X_Cluster_ADC",&Det_Cluster_ADC[6]);
   ClusterTree->SetBranchAddress("D3Y_Cluster_ADC",&Det_Cluster_ADC[7]);
   ClusterTree->SetBranchAddress("Dia_Cluster_ADC",&Dia_Cluster_ADC);
   ClusterTree->SetBranchAddress("D0X_Cluster_PedMean",&Det_Cluster_PedMean[0]);
   ClusterTree->SetBranchAddress("D0Y_Cluster_PedMean",&Det_Cluster_PedMean[1]);
   ClusterTree->SetBranchAddress("D1X_Cluster_PedMean",&Det_Cluster_PedMean[2]);
   ClusterTree->SetBranchAddress("D1Y_Cluster_PedMean",&Det_Cluster_PedMean[3]);
   ClusterTree->SetBranchAddress("D2X_Cluster_PedMean",&Det_Cluster_PedMean[4]);
   ClusterTree->SetBranchAddress("D2Y_Cluster_PedMean",&Det_Cluster_PedMean[5]);
   ClusterTree->SetBranchAddress("D3X_Cluster_PedMean",&Det_Cluster_PedMean[6]);
   ClusterTree->SetBranchAddress("D3Y_Cluster_PedMean",&Det_Cluster_PedMean[7]);
   ClusterTree->SetBranchAddress("Dia_Cluster_PedMean",&Det_Cluster_PedMean[8]);
   ClusterTree->SetBranchAddress("D0X_Cluster_PedWidth",&Det_Cluster_PedWidth[0]);
   ClusterTree->SetBranchAddress("D0Y_Cluster_PedWidth",&Det_Cluster_PedWidth[1]);
   ClusterTree->SetBranchAddress("D1X_Cluster_PedWidth",&Det_Cluster_PedWidth[2]);
   ClusterTree->SetBranchAddress("D1Y_Cluster_PedWidth",&Det_Cluster_PedWidth[3]);
   ClusterTree->SetBranchAddress("D2X_Cluster_PedWidth",&Det_Cluster_PedWidth[4]);
   ClusterTree->SetBranchAddress("D2Y_Cluster_PedWidth",&Det_Cluster_PedWidth[5]);
   ClusterTree->SetBranchAddress("D3X_Cluster_PedWidth",&Det_Cluster_PedWidth[6]);
   ClusterTree->SetBranchAddress("D3Y_Cluster_PedWidth",&Det_Cluster_PedWidth[7]);
   ClusterTree->SetBranchAddress("Dia_Cluster_PedWidth",&Det_Cluster_PedWidth[8]);
      
   current_event = 0;
   ClusterTree->GetEvent(current_event);
   cout<< "Loaded first event in ClusterTree: "<<event_number<<endl;
   cout<< "RunNumber is: "<<run_number<<endl;
   cout<< "StoreThreshold is: "<<store_threshold<<endl;
   
}

ClusterAnalyze::~ClusterAnalyze() {
   delete PedTree;
   delete PedFile;
   delete ClusterTree;
   delete ClusterFile;
}


void ClusterAnalyze::LoadSettings() {
   
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
}

void ClusterAnalyze::ParseIntArray(string value, vector<int> &vec) {
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

void ClusterAnalyze::Cluster_Clear(int det) {
   //reset channel counter
   Det_Cluster_NChannels[det]=0;
}

void ClusterAnalyze::ClearClusters() {
   //reset channel counter
   for(int i=0; i<9; i++) Det_Cluster_NChannels[i]=0;
}

float ClusterAnalyze::Cluster_GetCharge(int det) {
   double totalcharge = 0;
   for(int i=0; i<Det_Cluster_NChannels[det]; i++) {
      if(det<8) if(Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i]>Si_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) totalcharge += Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i];
      if(det==8) if(Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i]>Di_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) totalcharge += Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i];
   }
   return totalcharge;
}

int ClusterAnalyze::Cluster_GetNHits(int det) {
   int nhits = 0;
   for(int i=0; i<Det_Cluster_NChannels[det]; i++) {
      if(det<8) if(Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i]>Si_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) nhits++;
      if(det==8) if(Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i]>Di_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) nhits++;
   }
   return nhits;
}

int ClusterAnalyze::Cluster_GetNSeeds(int det) {
   int nseeds = 0;
   for(int i=0; i<Det_Cluster_NChannels[det]; i++) {
      if(det<8) if(Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i]>Si_Cluster_Seed_Factor*Det_Cluster_PedWidth[det][i]) nseeds++;
      if(det==8) if(Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i]>Di_Cluster_Seed_Factor*Det_Cluster_PedWidth[det][i]) nseeds++;
   }
   return nseeds;
}

float ClusterAnalyze::Cluster_Get1stMoment(int det) {
   double firstmoment = 0;
   double psadc = 0;
   for(int i=0; i<Det_Cluster_NChannels[det]; i++) {
      if(det<8) {
         psadc = Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i];
         if(psadc>Si_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) firstmoment += psadc*Det_Cluster_Channels[det][i];
      }
      if(det==8) {
         psadc = Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i];
         if(psadc>Di_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) firstmoment += psadc*Det_Cluster_Channels[det][i];
      }
   }
   return firstmoment/Cluster_GetCharge();
}

float ClusterAnalyze::Cluster_Get2ndMoment(int det) {
   double secondmoment = 0;
   double psadc = 0;
   for(int i=0; i<Det_Cluster_NChannels[det]; i++) {
      if(det<8) {
         psadc = Det_Cluster_ADC[det][i]-Det_Cluster_PedMean[det][i];
         if(psadc>Si_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) secondmoment += psadc*Det_Cluster_Channels[det][i]*Det_Cluster_Channels[det][i];
      }
      if(det==8) {
         psadc = Dia_Cluster_ADC[i]-Det_Cluster_PedMean[det][i];
         if(psadc>Di_Cluster_Hit_Factor*Det_Cluster_PedWidth[det][i]) secondmoment += psadc*Det_Cluster_Channels[det][i]*Det_Cluster_Channels[det][i];
      }
   }
   return secondmoment/Cluster_GetCharge();
}

float ClusterAnalyze::Cluster_Get1stMomentEtaCorrected(int det) {
   return 0;
}

float ClusterAnalyze::Cluster_GetPositionLocal(int det) {
   return 0;
}

float ClusterAnalyze::Cluster_GetPositionGlobal(int det) {
   return 0;
}

//sequentially cluster all events
void ClusterAnalyze::AnalyzeRun() {
   
   for(uint e=0; e<ClusterTree->GetEntries(); e++) {
      BookHistograms();
   }
   
   PedFile->Close();
   ClusterTree->Write();
   ClusterFile->Write();
   ClusterTree->Print();
   cout << "Number of Events in ClusterTree: " << ClusterTree->GetEntries() << endl;
   ClusterFile->Close();
}
