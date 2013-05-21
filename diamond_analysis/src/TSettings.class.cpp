/*
 * TSettings.class.cpp
 *
 *  Created on: 02.08.2011
 *      Author: Felix Bachmair
 *
 *
 *      todo: Write copy constructor
 */

#include "../include/TSettings.class.hh"
ClassImp(TSettings);
using namespace std;


TSettings::TSettings(TRunInfo *runInfo)
{
	cout<<"TSettings TRunInfo"<<endl;
	//  verbosity=runInfo->getVerbosity();
	setVerbosity(runInfo->getVerbosity());
	diamondMapping=0;
	fidCutsSelection = new TFidCutRegions();
	fidCuts3D = new TFidCutRegions();
	DefaultLoadDefaultSettings();
	this->runNumber=runInfo->getRunNumber();
	sys = gSystem;
	setRunDescription(runInfo->getRunDescription());
	stringstream fileNameStr;
	fileNameStr<<path<<"/"<<runInfo->getRunSettingsDir()<<"/settings."<<runInfo->getRunNumber();
	if (runInfo->getRunDescription().at(0)!='0')
		fileNameStr<<"-"<<runInfo->getRunDescription();
	fileNameStr<<".ini";
	cout<<"get Settingsfile: fileName =\""<<fileNameStr.str()<<"\""<<endl;
	int fileExist = existsDirectory(fileNameStr.str());
	cout<< "File Exists: "<<fileExist<<endl;
	cout<<endl;
	if(!fileExist){
		cout<<"Settingsfile: "<<fileNameStr.str()<<" does not exists"<<endl;
		cout<<"finish on press char:"<<flush;
		exit(-1);
	}

	if(getVerbosity())
		cout<<"TSettings:Create TSettings-member with file:\""<<fileNameStr.str()<<"\""<<endl;
	outputDir=runInfo->outputDir;//TODO:getOutputDir();
	cout<<runInfo->getInputDir()<<endl;
	setInputDir(runInfo->getInputDir());
	SetFileName(fileNameStr.str());
	checkSettings();
}

TSettings::TSettings(UInt_t runNumber){
	setVerbosity(0);//standard is 0
	if(verbosity)
		cout<<"TSettings:Create TSettings-member with file:\""<<fileName<<"\""<<endl;
	diamondMapping=0;
	fidCutsSelection = new TFidCutRegions();
	fidCuts3D = new TFidCutRegions();
	DefaultLoadDefaultSettings();
	SetFileName("SETTINGS.new.ini");
	this->runNumber=runNumber;
	sys = gSystem;
	path = sys->pwd();
	runDescription="";
	checkSettings();
}

TSettings::TSettings(string fileName,UInt_t runNumber){
	setVerbosity(0);//standard is 0
	if(verbosity)
		cout<<"TSettings:Create TSettings-member with file:\""<<fileName<<"\""<<endl;
	diamondMapping=0;
	fidCutsSelection = new TFidCutRegions();
	fidCuts3D = new TFidCutRegions();
	DefaultLoadDefaultSettings();
	this->runNumber=runNumber;
	sys = gSystem;
	path = sys->pwd();
	runDescription="";
	SetFileName(fileName);
	checkSettings();
	//	createSettingsRootFile();
}

TSettings::~TSettings(){

	//  saveSettings();
	//  settingsFile->Close();
	cout<<"delete Settings"<<endl;
}

void TSettings::checkSettings(){
	cout<<"Check Settings..."<<endl;
	if (!fidCutsSelection){
		fidCutsSelection = new TFidCutRegions();
	}
	fidCutsSelection->SetName("standard-FidCuts");
	if (isStandardSelectionFidCut==true){
		fidCutsSelection->Reset();
		cout<<"Standard FidCut: "<<endl;
		fidCutsSelection->SetName("standard-FidCuts");
		fidCutsSelection->addFiducialCut(getSi_avg_fidcut_xlow(),getSi_avg_fidcut_xhigh(),getSi_avg_fidcut_ylow(),getSi_avg_fidcut_yhigh());
	}
	if(!fidCutsSelection)
		fidCutsSelection =  new TFidCutRegions();
	fidCutsSelection->SetName("3D-FidCuts");

	if (isStandard3dFidCut==true){
			fidCuts3D->Reset();
			fidCutsSelection->SetName("3D-FidCuts");
			fidCuts3D->addFiducialCut(-1e9,1e9,-1e9,1e9);
		}
//	diamondPattern.Print();
	cout<<"NDiamond Patterns: "<<diamondPattern.getNPatterns()<<endl;

	this->checkAlignmentFidcuts();
	cout<<"Settings seems to be ok."<<endl;
	if(verbosity>5&&verbosity%2==1){char t; cin>>t;}
}

void TSettings::checkAlignmentFidcuts(){
	for(UInt_t i=0;i<this->alignmentFidCuts.size();i++){
		Int_t fidCutIndex = alignmentFidCuts[i];
		TFiducialCut* fidCut = this->getSelectionFidCuts()->getFidCut(fidCutIndex);
		if(fidCut==0){
			cout<<endl;
			this->getSelectionFidCuts()->Print(1);
			cout<<"fidCutIndex "<<fidCutIndex<<" seems not exist... ===> EXIT, Press a key and enter to confirm"<<endl;
			char t;
			cin>>t;
			exit(-1);
		}
		else{
			cout<<"fidCut No:"<<fidCutIndex<<endl;
			fidCut->Print(1);
		}

	}
	//exit(-1);
}

bool TSettings::existsDirectory(std::string dir){
	struct stat sta;
	int retVal = stat(dir.c_str(),&sta);
	return (retVal>=0);
}

std::string TSettings::get3dDiamondTreeFilePath(){
        stringstream path;
        path<<getAbsoluteOuputPath(false);
        path<<"selectionData."<<getRunNumber();
        if(this->isSpecialAnalysis())
            path<<"-"<<getRunDescription();
        path<<".root";
        return path.str();
}
void TSettings::goTo3dDiamondTreeDir(){
 goToDir(this->getAbsoluteOuputPath(false));
}

void TSettings::goToOutputDir(){
 goToDir(this->getOutputDir());
}



std::string TSettings::getRawTreeFilePath()
{
	stringstream path;
	path<<getAbsoluteOuputPath(false);
	path<<"rawData."<<getRunNumber()<<".root";
	return path.str();
}

std::string TSettings::getPedestalTreeFilePath()
{
	stringstream path;
	path<<getAbsoluteOuputPath(false);
	path<<"pedestalData."<<getRunNumber()<<".root";
	return path.str();
}

std::string TSettings::getClusterTreeFilePath()
{
	stringstream path;
	path<<getAbsoluteOuputPath(false);
	path<<"clusterData."<<getRunNumber()<<".root";
	return path.str();
}

std::string TSettings::getAlignmentFilePath(TSettings::alignmentMode mode)
{

	stringstream path;
	path<<getAbsoluteOuputPath(false);
	if(mode == TSettings::transparentMode)
		path << "alignment-trans.";
	else
		path<<"alignment.";
	path<<getRunNumber();
	if(this->isSpecialAnalysis())
		path<<"-"<<getRunDescription();
	path<<".root";
	return path.str();
}

std::string TSettings::getEtaDistributionPath(Int_t step){
	stringstream output;
	output<<this->getAbsoluteOuputPath(false)<<"/etaCorrection";
	if(step>=0)
		output<<"_Step"<<step;
	output<<"."<<runNumber<<".root";
	return output.str();
}

std::string TSettings::getSelectionTreeFilePath()
{
	stringstream path;
	path<<getAbsoluteOuputPath(false);
	path<<"selectionData."<<getRunNumber();
	if(this->isSpecialAnalysis())
		path<<"-"<<getRunDescription();
	path<<".root";
	return path.str();
}

std::string TSettings::getCrossTalkFactorsFileName(){
	string hName = this->getAbsoluteOuputPath(false);
	hName.append(TString::Format("crossTalkCorrectionFactors.%d",runNumber));
	if(isSpecialAnalysis())
		hName.append("-"+getRunDescription());
	hName.append(".txt");
	return hName;
}

void TSettings::goToDir(std::string dir){
	if(this->getVerbosity()>3)cout<<"\ncurrent Dir: "<<sys->pwd()<<endl;
	if(this->getVerbosity()>3)cout<<"goTo Dir: "<<dir<<endl;
	if(!existsDirectory(dir))
		sys->mkdir( dir.c_str(),true);
	sys->cd(dir.c_str());
	if(this->getVerbosity()>3)cout<<"new Dir: "<<sys->pwd()<<endl;
}
void TSettings::goToRawTreeDir(){
	goToDir(this->getAbsoluteOuputPath(false));
}

void TSettings::goToSelectionTreeDir(){
	goToDir(this->getAbsoluteOuputPath(false));
}

/*void TSettings::goToOutputDir(){
	goToDir(this->getOutputDir());
}*/

void TSettings::goToPedestalAnalysisDir(){
	goToDir(this->getAbsoluteOuputPath(isSpecialAnalysis()).append("/pedestalAnalysis/"));
}

void TSettings::goToClusterAnalysisDir(){
	goToDir(this->getAbsoluteOuputPath(isSpecialAnalysis()).append("/clustering/"));
}

void TSettings::goToAlignmentDir(alignmentMode mode) {
	if(mode == transparentMode)
		goToDir(this->getAbsoluteOuputPath(true).append("/alignment2/"));
	else
		goToDir(this->getAbsoluteOuputPath(true).append("/alignment/"));
}


std::string TSettings::getAlignmentDir(TSettings::alignmentMode mode) {
	if(mode == transparentMode)
		return this->getAbsoluteOuputPath(true).append("/alignment2/");
	else
		return this->getAbsoluteOuputPath(true).append("/alignment/");
}

std::string TSettings::getTransparentAnalysisDir(alignmentMode mode){
	if (mode == transparentMode)
		return this->getAbsoluteOuputPath(true).append("/transparentAnalysis2/");
	return this->getAbsoluteOuputPath(true).append("/transparentAnalysis/");
}

TSettings::enumRunDescription TSettings::getAnalysedDiamond(){
	cout<<"rundes: "<<getRunDescription()<<endl;
	if(!isSpecialAnalysis())
		return allDia;
	if(getRunDescription()=="left"||getRunDescription()=="1")
		return leftDia;
	if(getRunDescription()=="right"||getRunDescription()=="2")
		return rightDia;
	return unknown;
}


void TSettings::setRunDescription(std::string runDescription)
{
	this->runDescription = runDescription;
}

std::string TSettings::getAbsoluteOuputPath(bool withRunDescribtion){
	if(this->getVerbosity()>3)cout<<"Absolute PATH:"<<withRunDescribtion<<" "<<(int)(runDescription.at(0)!='0')<<endl;
	stringstream output;
	output<<this->getOutputDir();
	output<<"/"<<runNumber<<"/";
	if(withRunDescribtion&&runDescription.at(0)!='0'){
		output<<runDescription<<"/";
		sys->MakeDirectory(output.str().c_str());
	}
	if(this->getVerbosity()>3)cout<<"OUTPUT: "<<output.str()<<endl;
	return output.str();
}
void TSettings::SetFileName(string newFileName){
	if(getVerbosity())
		cout<<"TSettings::SetFileName:\""<<newFileName<<"\""<<endl;
	fileName=newFileName;
	LoadSettings();
	if (fidCutsSelection){
		fidCutsSelection->Print();
	}
}

void TSettings::LoadSettings(){
	if (getVerbosity())
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
		if(key == "asymmetricSample") Parse(key,value,bAsymmetricSample);
		if(key == "SaveAllFilesSwitch") Parse(key,value,SaveAllFilesSwitch);
		if(key == "siliconAlignmentSteps")Parse(key,value,siliconAlignmentSteps);
		if(key == "ClosePlotsOnSave") Parse(key,value,ClosePlotsOnSave);
		if(key == "IndexProduceSwitch")Parse(key,value,IndexProduceSwitch);
		if(key == "snr_plots_enable") Parse(key,value,snr_plots_enable);
		if(key == "fix_dia_noise") Parse(key,value,fix_dia_noise);
		if(key == "single_channel_analysis_channels") Parse(key, value,single_channel_analysis_channels);
		if(key == "single_channel_analysis_enable") Parse(key,value,single_channel_analysis_enable);
		if(key == "single_channel_analysis_eventwindow") Parse(key,value,single_channel_analysis_eventwindow);
		if(key == "CMN_corr_low") Parse(key,value,CMN_corr_low);
		if(key == "CMN_corr_high")Parse(key,value,CMN_corr_high);
		if(key == "resetAlignment")Parse(key,value,bResetAlignment);
		if(key == "CMN_cut")Parse(key,value,CMN_cut);
		if(key == "DO_CMC")Parse(key,value,DO_CMC);
		if(key == "CMN_cut") Parse(key,value,CMN_cut);
		if(key == "Iter_Size")Parse(key,value,Iter_Size);
		if(key == "Taylor_speed_throttle") Parse(key,value,Taylor_speed_throttle);
		if(key == "dia_input") Parse(key,value,dia_input);
		if(key == "alignment_training_track_fraction") Parse(key,value,alignment_training_track_fraction);
		if(key == "alignment_training_track_number") Parse(key,value,alignment_training_track_number);
		if(key == "alignment_training_method"){
			cout << key.c_str() << " = " << value.c_str() << endl;
			int method = (int)strtod(value.c_str(),0);
			if(method >=0&&method<=2 )
				setTrainingMethod((enumAlignmentTrainingMethod)method);
			else
				cerr<<"Not a valid Input for alignment Training Method : "<<method<<endl;
		}
		if(key == "alignment_training_fidcuts") {
			ParseIntArray(key,value,alignmentFidCuts);
			if(verbosity){
				for(UInt_t i = 0; i<alignmentFidCuts.size();i++){
					cout<<TCluster::Intent(1)<<"Region "<<alignmentFidCuts.at(i)<<endl;
				}
			}
		}
		if (key == "res_keep_factor") {ParseFloat(key,value,res_keep_factor);}
		if (key == "MinimalAbsoluteEtaValue") ParseFloat(key,value,minAbsEtaVal);
		if(key == "Si_Pedestal_Hit_Factor") ParseFloat(key,value,Si_Pedestal_Hit_Factor);
		if(key == "Di_Pedestal_Hit_Factor") ParseFloat(key,value,Di_Pedestal_Hit_Factor);
		if(key == "Si_Cluster_Seed_Factor") ParseFloat(key,value,Si_Cluster_Seed_Factor);
		if(key == "Di_Cluster_Seed_Factor") ParseFloat(key,value,Di_Cluster_Seed_Factor);
		if(key == "Si_Cluster_Hit_Factor") ParseFloat(key,value,Si_Cluster_Hit_Factor);
		if(key == "Di_Cluster_Hit_Factor") ParseFloat(key,value,Di_Cluster_Hit_Factor);
		if(key == "eta_lowq_slice_low") ParseFloat(key,value,eta_lowq_slice_low);
		if(key == "eta_lowq_slice_hi") ParseFloat(key,value,eta_lowq_slice_hi);
		if(key == "eta_hiq_slice_low") ParseFloat(key,value,eta_hiq_slice_low);
		if(key == "eta_hiq_slice_hi") ParseFloat(key,value,eta_hiq_slice_hi);
		if(key == "etavsq_n_landau_slices") ParseInt(key,value,etavsq_n_landau_slices);
		if(key == "alignment_x_offsets") ParseFloatArray(key,value,alignment_x_offsets);
		if(key == "alignment_y_offsets") ParseFloatArray(key, value,alignment_y_offsets);
		if(key == "alignment_phi_offsets") ParseFloatArray(key,value,alignment_phi_offsets);
		if(key == "alignment_z_offsets") ParseFloatArray(key,value,alignment_z_offsets);
		if(key == "D0X_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[0]);
		if(key == "D0Y_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[1]);
		if(key == "D1X_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[2]);
		if(key == "D1Y_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[3]);
		if(key == "D2X_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[4]);
		if(key == "D2Y_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[5]);
		if(key == "D3X_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[6]);
		if(key == "D3Y_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[7]);
		if(key == "Dia_channel_screen_channels") ParseScreenedChannelArray(key,value,Det_channel_screen_channels[8]);
		if(key == "D0X_channel_screen_regions")  ParseScreenedChannelArray(key,value,Det_channel_screen_regions[0]);
		if(key == "D0Y_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[1]);
		if(key == "D1X_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[2]);
		if(key == "D1Y_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[3]);
		if(key == "D2X_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[4]);
		if(key == "D2Y_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[5]);
		if(key == "D3X_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[6]);
		if(key == "D3Y_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[7]);
		if(key == "Dia_channel_screen_regions") ParseScreenedChannelArray(key,value,Det_channel_screen_regions[8]);
		if(key == "chi2Cut3D") ParseFloat(key,value,chi2Cut3D);
		if(key == "si_avg_fidcut_xlow") {ParseFloat(key,value,si_avg_fidcut_xlow);};
		if(key == "si_avg_fidcut_xhigh") ParseFloat(key,value,si_avg_fidcut_xhigh);
		if(key == "si_avg_fidcut_ylow") ParseFloat(key,value,si_avg_fidcut_ylow);
		if(key == "si_avg_fidcut_yhigh") ParseFloat(key,value,si_avg_fidcut_yhigh);
		if(key == "selectionFidCut") {if (!fidCutsSelection) fidCutsSelection=new TFidCutRegions();ParseFidCut(key,value,fidCutsSelection,isStandardSelectionFidCut);}
		if(key == "3dFitCut"){if (!fidCuts3D) fidCuts3D=new TFidCutRegions();ParseFidCut(key,value,fidCuts3D,isStandard3dFidCut);}
		if(key == "pulse_height_num_bins") ParseInt(key,value,pulse_height_num_bins);
		if(key == "pulse_height_si_max") ParseFloat(key,value,pulse_height_si_max);
		if(key == "pulse_height_di_max")  ParseFloat(key,value,pulse_height_di_max);
		if(key == "snr_distribution_si_max")  Parse(key,value,snr_distribution_si_max);
		if(key == "snr_distribution_di_max")  Parse(key,value,snr_distribution_di_max);
		if (key == "alignment_chi2") Parse(key,value,alignment_chi2);
		if (key == "UseAutoFidCut") Parse(key,value,UseAutoFidCut);
		if (key == "nDiamonds")this->setNDiamonds(ParseInt(key,value));
		if (key == "AlternativeClustering") Parse(key,value,AlternativeClustering);
		if(key == "store_threshold") Parse(key,value,store_threshold);
		if(key == "plotChannel_on") Parse(key,value,plotChannel_on);
		if(key == "SingleChannel2000plots") Parse(key,value,SingleChannel2000plots);
		if(key == "makeDiamondPlots")  Parse(key,value,makeDiamondPlots);
		if(key == "alignmentPrecision_Offset") Parse(key,value,alignmentPrecision_Offset);
		if(key == "alignmentPrecision_Angle") Parse(key,value,alignmentPrecision_Angle);
		if(key == "makeHits2D")  Parse(key,value,makeHits2D);
		if(key == "makeNoise2D")  Parse(key,value,makeNoise2D);
		if(key == "makePullDist")  Parse(key,value,makePullDist);
		if(key == "makePedRMSTree")  Parse(key,value,makePedRMSTree);
		if(key == "eventPrintHex")  Parse(key,value,eventPrintHex);
		if(key == "plottedChannel")  Parse(key,value,plottedChannel);
		if(key == "high_rms_cut")  Parse(key,value,high_rms_cut);
		if(key == "rms_cut")  Parse(key,value,rms_cut);
		if(key == "zoomDiamondPlots")  Parse(key,value,zoomDiamondPlots);
		if(key == "singleTrack2D")  Parse(key,value,singleTrack2D);
		if(key == "singleTrack2DmaxClusterSize")  Parse(key,value,singleTrack2DmaxClusterSize);
		if(key == "maxNoise2D")  Parse(key,value,maxNoise2D);
		if(key == "clusterHitFactors") ParseFloatArray(key, value,clusterHitFactors);
		if(key == "clusterSeedFactors") ParseFloatArray(key, value,clusterSeedFactors);
		if(key == "doAllAlignmentPlots") Parse(key,value,bDoAllAlignmentPlots);
		if(key == "pitchWidthDia") Parse(key,value,pitchWidthDia);
		if(key == "pitchWidthSil") Parse(key,value,pitchWidthSil);
		if(key == "diamondPattern") ParsePattern(key,value);
		if(key == "TransparentAlignment") ParseBool(key,value,bTransparentAlignment);
		if(key == "diamondMapping") {
			cout<<key<<" = "<<value.c_str()<<endl;
			std::vector<int>vecDiaMapping;
			ParseIntArray(key, value,vecDiaMapping);
			if(diamondMapping==0)
				delete diamondMapping;
			diamondMapping=new TChannelMapping(vecDiaMapping);
			diamondMapping->PrintMapping();
			cout<<diamondMapping<<endl;
			getDetChannelNo(0);
		}
		if(key == "Dia_DetectorChannels") {
			cout<<key<<" = "<<value.c_str()<<endl;
//			vector<string> vecDetectorChannelString;
//			ParseStringArray(key, value,vecDetectorChannelString);
			diamondPattern.resetPattern();
			ParseRegionArray(key, value,vecDiaDetectorAreasInChannel);
			Int_t detChannel = -1;
			for(UInt_t i=0;i<vecDiaDetectorAreasInChannel.size();i++){
				cout<<i<<" "<<vecDiaDetectorAreasInChannel[i].first<<" "<<vecDiaDetectorAreasInChannel[i].second<<endl;
				if (vecDiaDetectorAreasInChannel[i].second< detChannel){
					cout<<"this Definition of DetectorChannels doesn't work, please update"<<endl;
					exit(-1);
				}
			}
		}
		if(key == "Dia_ClusterSeedFactors"){
			ParseFloatArray(key, value,vecClusterSeedFactorsDia);
			if((Int_t)vecClusterSeedFactorsDia.size()!=getNDiaDetectorAreas()){
				cerr<<"The number of defined ClusterSeedFactors for the diamond Areas does not fit with the number of defined areas:\t"<<flush;
				cerr<<vecClusterSeedFactorsDia.size()<<" "<<getNDiaDetectorAreas()<<endl;
				exit(-1);
			}
			cout<<key<<endl;
			for(UInt_t i=0;i<vecClusterSeedFactorsDia.size();i++)
				cout<<i<<"\t"<<getDiaDetectorArea(i).first<<"-"<<getDiaDetectorArea(i).second<<": "<<vecClusterSeedFactorsDia.at(i)<<endl;
		}
		if(key == "Dia_ClusterHitFactors"){
			ParseFloatArray(key, value,vecClusterHitFactorsDia);
			if((Int_t)vecClusterHitFactorsDia.size()!=getNDiaDetectorAreas()){
				cerr<<"The number of defined ClusterHitFactors for the diamond Areas does not fit with the number of defined areas:\t"<<flush;
				cerr<<vecClusterHitFactorsDia.size()<<" "<<getNDiaDetectorAreas()<<endl;
				exit(-1);
			}
			cout<<key<<endl;
			for(UInt_t i=0;i<vecClusterHitFactorsDia.size();i++)
				cout<<i<<"\t"<<getDiaDetectorArea(i).first<<"-"<<getDiaDetectorArea(i).second<<": "<<vecClusterHitFactorsDia.at(i)<<endl;
		}
        if(key == "is3dDiamond"){
           cout<<key<<" =" <<value.c_str()<<endl;
           b3dDiamond = (bool)strtod(value.c_str(),0);
        }
		/*if(key == "store_threshold") {//TODO It's needed in settings reader
	         cout << key.c_str() << " = " << value.c_str() << endl;
	        store_threshold = (float)strtod(value.c_str(),0);
	      }*/
	}

//	for(UInt_t ch=0;ch<TPlaneProperties::getNChannelsDiamond();ch++)
//		cout<<setw(3)<<ch<<":"<<setw(3)<<getDiaDetectorAreaOfChannel(ch)<<"\t"<<getClusterSeedFactor(TPlaneProperties::getDetDiamond(),ch)
//			<<"--"<<getClusterHitFactor(TPlaneProperties::getDetDiamond(),ch)<<endl;
//	char t;
//	cin>>t;
	file.close();

	for(int det=0; det<9; det++) {
		this->Det_channel_screen[det].ScreenChannels(this->getDet_channel_screen_channels(det));
		//this->getDet_channel_screen(det).ScreenRegions(this->getDet_channel_screen_regions(det));
		cout<<"Detector "<<det<<" screened channels: ";
		this->getDet_channel_screen(det).PrintScreenedChannels();
		cout<<endl;
	}

	for(int det=0;det<9;det++){
		cout<<"analyse detector "<<det<< " with "<<getClusterSeedFactor(det,0)<<"/"<<getClusterHitFactor(det,0)<<endl;
	}
	checkSettings();
}

void TSettings::LoadDefaultResolutions(){
	alignment_resolutions.resize(9);
	//X
	alignment_resolutions[0] = 2.3;
	alignment_resolutions[2] = 2.2;
	alignment_resolutions[4] = 1.6;
	alignment_resolutions[6] = 1.7;

	//Y
	alignment_resolutions[1] = 1.6;
	alignment_resolutions[3] = 1.5;
	alignment_resolutions[5] = 1.9;
	alignment_resolutions[7] = 2.0;

	alignment_resolutions[8] = 10.0;//todo
}
void TSettings::DefaultLoadDefaultSettings(){

	if(getVerbosity())
		cout<<"TSettings::LoadDefaultSettings"<<endl;
	//default general settings
	isStandardArea=true;
	isStandardSelectionFidCut=true;
	runDescription="";
	SaveAllFilesSwitch = 1; //1 for save files, 0 for don't
	ClosePlotsOnSave = 1;
	IndexProduceSwitch = 1;
	store_threshold=2;
	//default pedestal settings
	fix_dia_noise = -1;//7.7; // fix_dia_noise<0 disables diamond noise-fixing
	dia_input = 0; // 1 for 2006 and 0 for the rest
	DO_CMC = 0;
	CMN_cut = 4;  //Should be less than or equal to CMN_coor_high
	Iter_Size = 500; //buffer size
	Taylor_speed_throttle = 1000; //# of events to recalculate RMS the old way; set to 1 to disable
	CMN_corr_high=7;
	CMN_corr_low=3;
	bDoAllAlignmentPlots = false;


	res_keep_factor=2;
	alignmentPrecision_Offset = 0.01;
	alignmentPrecision_Angle = 0.001;
	alignment_chi2 = 4.0;
	alignment_training_track_fraction = 0.25;
	alignment_training_track_number = 100000;
	trainingMethod=enumEvents;
	bResetAlignment=false;

	//default clustering settings
	snr_plots_enable = 0;
	bTransparentAlignment=false;

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
	transparentChi2 = 5;
	UseAutoFidCut = 0;
	nDiamonds=1;
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


	maxBufferPlots = 100;
	rms_sigma_difference_cut = 0.3;

	high_rms_cut = 1; //cut on absolute rms value instead of comparing to Gaussian
	rms_cut = 20.; //value to use if high_rms_cut

	siliconAlignmentSteps=5;
	diamondAlignmentSteps=5;
	zoomDiamondPlots = 0; //zoom in on DC_Pedestal (100 event / window)
	pitchWidthDia = 50; // in mum
	pitchWidthSil = 50; // in mum

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
    b3dDiamond =false;
    silPitchWidth=50;//in um
    diaPitchWidth=50;//in um
    diaOffsetMetricSpace=0;
    diaStartingChannel=0;
	cout<<"Print DefaultMapping:"<<endl;
    //	cout<<"Print DefaultMapping:"<<endl;
	//	diamondMapping.PrintMapping();
	diamondPattern.loadStandardPitchWidthSettings();
	alignmentFidCuts.clear();
	alignmentFidCuts.push_back(1);
	cout<<"DONE"<<endl;
	isStandardSelectionFidCut=true;
	isStandardArea = true;
	chi2Cut3D=4.0;
	bAsymmetricSample=false;
	minAbsEtaVal = .05;
	bUseUserResolutionInput = false;
	LoadDefaultResolutions();

	checkSettings();
}



/**
 * function which parses an string of the format  '{XX,XX,XX,XX,XX}' to an
 * vector of strings
 */
void TSettings::ParseStringArray(string key, string value, vector<string> &vec){
	cout << key.c_str() << " = " << value.c_str() << endl;

	if(value.find('{')==string::npos||value.find('}')==string::npos){
		cerr<<"the string \'"<<value<<"\' cannot be parsed as a float array since bracket is missing"<<endl;
		exit(-1);
	}
	string::size_type beginning = value.find_first_of('{')+1;
	string::size_type ending = value.find_last_of('}');
	string analyseString = value.substr(beginning,ending-beginning);
	//  cout<<"analyze: \'"<<analyseString<<"\'"<<endl;
	size_t i;
	while((i=analyseString.find(','))!=string::npos){
		string data = analyseString.substr(0,i);
		vec.push_back(data);
		analyseString = analyseString.substr(i+1);
		//    cout<<"analyseString: \'"<<analyseString<<"\'"<<endl;
	}
	string data = analyseString.substr(0,i);
	vec.push_back(data);
}

bool TSettings::ParseFloat(string key, string value, float &output){
	cout << key.c_str() << " = " << value.c_str() << endl;
	output = (float)strtod(value.c_str(),0);
	return true;
}

bool TSettings::ParseInt(string key, string value, int &output){
	cout << key.c_str() << " = " << value.c_str() << endl;
	output = (int)strtod(value.c_str(),0);
	return true;
}

bool TSettings::ParseInt(string key, string value, UInt_t &output){
	cout << key.c_str() << " = " << value.c_str() << endl;
	output = (UInt_t)strtod(value.c_str(),0);
	return true;
}

bool TSettings::ParseBool(string key, string value, bool &output){
	cout << key.c_str() << " = " << value.c_str() << endl;
	output = (bool)strtod(value.c_str(),0);
	return true;
}

void TSettings::ParseFloatArray(string key, string value, vector<float> &vec) {
	cout << key.c_str() << " = " << value.c_str() << endl;
	std::vector <std::string> stringArray;
	ParseStringArray(key, value,stringArray);
	vec.clear();
	//  cout<<value<<" --> Array length: "<<stringArray.size()<<endl;
	for(UInt_t i=0;i<stringArray.size();i++)
		vec.push_back((float)strtod(stringArray.at(i).c_str(),0));
}

void TSettings::ParseIntArray(string key, string value, vector<int> &vec) {
	cout << key.c_str() << " = " << value.c_str() << endl;
	std::vector <std::string> stringArray;
	ParseStringArray(key, value,stringArray);
	vec.clear();
	//    cout<<value<<" --> Array length: "<<stringArray.size()<<endl;
	for(UInt_t i=0;i<stringArray.size();i++)
		vec.push_back((int)strtod(stringArray.at(i).c_str(),0));
}

void TSettings::ParseScreenedChannelArray(std::string key, std::string value, std::vector< int> &vec){
//	cout << key.c_str() << " = " << value.c_str() << endl;
	std::vector <std::string> stringArray;
	ParseStringArray(key, value,stringArray);
	vec.clear();
	for(UInt_t i=0;i<stringArray.size();i++){
		std::string val= stringArray.at(i);
		cout<<i<<"/"<<stringArray.size()<<"\tval: '"<<val<<"'"<<endl;
		if(val.find('-')!=string::npos){
			std::pair< std::string,std::string > region = ParseRegionString(key, stringArray[i]);
			Int_t begin = ParseInt(region.first.c_str());
			Int_t end = ParseInt(region.second.c_str());
			if(verbosity>5){
				cout<<"screen channel "<<begin<<"-"<<end;
				if(verbosity%2==1){
					cout<<"\tPress a key and enter to continue."<<flush;
					char t; cin>>t;}
				else
					cout<<endl;
			}
			for(Int_t i=begin;i<=end;i++)
				vec.push_back(i);
		}
		else
			vec.push_back(ParseInt(val));
	}

//	this->ParseIntArray(key,value,vec);
}

void TSettings::ParseRegionArray(string key, string value, std::vector< std::pair<Int_t,Int_t> > &vec){
	cout << "Parsing: "<<key.c_str() << " = " << value.c_str() << endl;
	std::vector <std::string> stringArray;
	ParseStringArray(key, value,stringArray);
	vec.clear();
	for(UInt_t i=0;i<stringArray.size();i++){
		std::pair< std::string,std::string > region = ParseRegionString(key, stringArray[i]);
		Int_t begin = (int)strtod(region.first.c_str(),0);
		Int_t end = (int)strtod(region.second.c_str(),0);
		cout<<i<<" "<<begin<<"-"<<end<<endl;
		if(begin<end){
			Float_t pos = begin*getDiamondPitchWidth();
			cout<<pos<<","<<begin<<"-"<<end<<endl;
			diamondPattern.addPattern(getDiamondPitchWidth(),pos,begin,end);
		}
	}

}


void TSettings::ParsePattern(std::string key, std::string value){
	cout<< "\nParse Pattern"<<endl;
	vector<Float_t> vecEntries;
	ParseFloatArray(key,value,vecEntries);
	if(vecEntries.size()==4){
		Float_t pos = vecEntries[0];
		Float_t pitchWidth = vecEntries[1];
		UInt_t firstCh = vecEntries[2];
		UInt_t lastCh = vecEntries[3];
		cout<<"Position: "<<pos<<"\tpw:"<<pitchWidth<<" first: "<<firstCh<<" last: lastCh"<<endl;
		if(diamondPattern.isStandardPitchWidth())
			diamondPattern.clear();
		diamondPattern.addPattern(pitchWidth,pos,firstCh,lastCh);
		cout<<"new Area of Interest: no. "<<diamondPattern.getNPatterns()-1<<" "<<firstCh<<"-"<<lastCh<<endl;
		PrintPatterns();
		isStandardArea=false;
	}
	else
		cout<<"vecEntries.size(): "<<vecEntries.size()<<endl;
	if (verbosity>5&&verbosity%2==1){
		cout<<"Press a key and enter: "<<flush;
		char t;
		cin>>t;
	}

}


void TSettings::ParseFidCut(std::string key, std::string value, TFidCutRegions* fidCutRegions,bool &isStandardFidCut){
	cout<< "\nParse FidCut: "<<value<<endl;

	if (fidCutRegions==0){
		cerr<<"TSettings::ParseFidCut: Couldn't Parse Since fidCutRegions ==0 "<<fidCutRegions<<endl;
		return;
	}
	std::vector <std::string> stringArray;
	ParseStringArray(key, value,stringArray);
	if(stringArray.size()==2){
		std::pair< std::string,std::string > region = ParseRegionString(key, stringArray[0]);
		Float_t beginX = (int)strtod(region.first.c_str(),0);
		Float_t endX = (int)strtod(region.second.c_str(),0);
		region = ParseRegionString(key, stringArray[1]);
		Float_t beginY = (int)strtod(region.first.c_str(),0);
		Float_t endY = (int)strtod(region.second.c_str(),0);
		if(beginX<endX&&beginY<endY){
			if(isStandardFidCut){
				fidCutRegions->Reset();
				isStandardFidCut=false;
			}
			fidCutRegions->addFiducialCut(beginX,endX,beginY,endY);
		}
		else
			cerr<<"TSettings::ParseFidCut: Cannot Parse FidCut - entries are wrong, "<<beginX<<"-"<<endX<<", "<<beginY<<"-"<<endY<<endl;
	}
	else
		cerr<<"TSettings::ParseFidCut: Cannot Parse FidCut - Size of vecEntries does not fit: "<<stringArray.size()<<endl;
	if (verbosity>5 &&verbosity%2==1){
		cout<<"Press a key and enter: "<<flush;
		char t;
		cin>>t;
	}
}

std::pair< std::string,std::string > TSettings::ParseRegionString(string key, string value){
	cout << key.c_str() << " = " << value.c_str() << endl;
	size_t index = value.find_first_of(":-");
	string beginString,endString;
	if(index!=string::npos){
		beginString = value.substr(0,index);
		endString = value.substr(index+1);
	}
	else{
		beginString = "";
		endString = "";
	}
	return make_pair(beginString,endString);
}
/**
 * TODO Ueberarbeiten so dass beliebige werte fuer jeden detector eingetragen werden koennen
 * @param det
 * @return
 */
Float_t TSettings::getClusterSeedFactor(UInt_t det,UInt_t ch){
//	cout<<"get Cluster Seed Factor: "<<det<<" "<<clusterSeedFactors.size()<<endl;
	if(TPlaneProperties::isDiamondDetector(det)){
		Int_t area = getDiaDetectorAreaOfChannel(ch);
//		cout<<"Diamond: "<<det<<":"<<ch<<"--->"<<area<<endl;
		if ((Int_t)vecClusterSeedFactorsDia.size()>area&&area>-1)
					return vecClusterSeedFactorsDia.at(area);
		else{
			if(det<clusterSeedFactors.size())
					return clusterSeedFactors.at(det);
			return getDi_Cluster_Seed_Factor();
		}
	}
	if(det<clusterSeedFactors.size())
			return clusterSeedFactors.at(det);
	return getSi_Cluster_Seed_Factor();
}

Float_t TSettings::getClusterHitFactor(UInt_t det,UInt_t ch){

	if(TPlaneProperties::isDiamondDetector(det)){
		Int_t area = getDiaDetectorAreaOfChannel(ch);
//		cout<<det<<":"<<ch<<"--->"<<area<<endl;
		if ((Int_t)vecClusterHitFactorsDia.size()>area&&area>-1)
			return vecClusterHitFactorsDia.at(area);
		else {
			if(clusterHitFactors.size()>det)
				return clusterHitFactors.at(det);
			return getDi_Cluster_Hit_Factor();
		}
	}
	if(clusterHitFactors.size()>det)
		return clusterHitFactors.at(det);
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
	if(det<TPlaneProperties::getNDetectors()&&ch<TPlaneProperties::getNChannels(det))
		return this->Det_channel_screen[det].isScreened(ch);
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
	//  cout<<"SAVE SETTINGS TO ROOT FILE"<<endl;
	//  settingsFile->cd();
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
	//  stringstream name;
	//  name << path<< "/"<<this->runNumber<<"/Settings."<<this->runNumber<<".root";
	////  settingsFile = TFile::Open(name.str().c_str());
	//  cout<<"Open settings from root file: "<<endl;
	//  settingsFile = new TFile(name.str().c_str(),"Update");
	//  if(settingsFile->IsZombie()){
	//    cout<<"file does not exist create new one!"<<endl;
	//    delete settingsFile;
	//    settingsFile = new TFile(name.str().c_str(),"RECREATE");
	//  }
	//  else{
	//    compareSettings();
	//  }


}

//todo
UInt_t TSettings::getMaxAllowedClustersize(UInt_t det){
	if(TPlaneProperties::isDiamondDetector(det))
		return 8+2;//two because of outer additional channels
	else
		return 4+2;//two because of outer additional channels

}

/**
 * @todo
 * @param k
 */
void TSettings::PrintPatterns(int k){

}

/**
 * @todo
 */
void TSettings::Print()
{
}

void TSettings::loadSettingsFromRootFile()
{
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


void TSettings::setFidCut(TFiducialCut *fidCut){
	if(fidCut==0) return;
	setSi_avg_fidcut_xhigh(fidCut->GetXHigh());
	setSi_avg_fidcut_yhigh(fidCut->GetYHigh());
	setSi_avg_fidcut_xlow(fidCut->GetXLow());
	setSi_avg_fidcut_ylow(fidCut->GetYLow());
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
	trainingMethod = enumFraction;
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

void TSettings::setNDiamonds(UInt_t nDia){
	this->nDiamonds=nDia;
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
//
//void TSettings::setAlignmentTrainingTrackNumber(){
//	alignment_training_track_number = alignmentTrainingTrackNumber;
//	trainingMethod =  enumEvents;
//}


bool TSettings::useForAlignment(UInt_t eventNumber, UInt_t nEvents) {
	if(getTrainingMethod()==enumEvents)
		return eventNumber<=getAlignmentTrainingTrackNumber();
	else{
		Float_t fraction = (Float_t)eventNumber/(Float_t)nEvents;
		return fraction<=getAlignment_training_track_fraction();
	}
	return false;
}

bool TSettings::isInAlignmentFiducialRegion(Float_t xVal,Float_t yVal){

	Int_t fidCutRegion = this->getSelectionFidCuts()->getFiducialCutIndex(xVal,yVal);
	if(verbosity>6)cout<<" isInAlignmentFiducialRegion\t"<<fidCutRegion<<flush;
	for(UInt_t i=0; i < alignmentFidCuts.size();i++)
		if(alignmentFidCuts.at(i)==fidCutRegion){
			if(verbosity>6)cout<<"\tTrue"<<endl;
			return true;
		}
	if(verbosity>6)cout<<"\tFalse"<<endl;
	return false;
}

Int_t TSettings::getVerbosity(){
	return this->verbosity;
}

bool TSettings::isInDiaDetectorArea(Int_t ch,Int_t area){
	if(area<getNDiaDetectorAreas()){
		pair<Int_t, Int_t > interval = getDiaDetectorArea(area);
		return interval.first <= ch && ch <= interval.second;
	}
	else return false;
}
std::pair< Int_t , Int_t > TSettings::getDiaDetectorArea(int n){
		return diamondPattern.getInterval(n);
}
Float_t TSettings::getMinDiamondChannel(){
	if(getNDiaDetectorAreas()>0){
		Float_t min = getDiaDetectorArea(0).first;
		for(Int_t area = 1; area < getNDiaDetectorAreas(); area++)
			min = TMath::Min(min,(Float_t)getDiaDetectorArea(area).first);
		return min;
	}
	return 0;
}

Float_t TSettings::getMaxDiamondChannel(){
	if(getNDiaDetectorAreas()>0){
		Float_t max = getDiaDetectorArea(0).second;
		for(Int_t area = 1; area < getNDiaDetectorAreas(); area++)
			max = TMath::Max(max,(Float_t)getDiaDetectorArea(area).second);
		return max;
		}
	return TPlaneProperties::getNChannelsDiamond();

}
int TSettings::getDiaDetectorAreaOfChannel(Int_t ch){
	for(Int_t area = 0; area < getNDiaDetectorAreas(); area++)
		if(isInDiaDetectorArea(ch,area))
			return area;
//	cout<<"cannot find ch "<<ch<<" in "<<getNDiaDetectorAreas()<<""<<endl;
	if(verbosity>5)
		diamondPattern.showPatterns();
	return -1;
}

bool TSettings::isDiaDetectorAreaBorderChannel(UInt_t ch){
//		cout<<"check isDiaDetectorAreaBorderChannel "<<ch<<"\t"<<flush;
	for(Int_t area =0;area< getNDiaDetectorAreas();area++){
//		cout<<area<<":";
		UInt_t leftBorder = getDiaDetectorArea(area).first;
		UInt_t rightBorder = getDiaDetectorArea(area).second;
//		cout<<(ch==leftBorder)<<(ch==rightBorder)<<"\t"<<flush;
		if ( ch == leftBorder || ch == rightBorder ){
//				cout<<"\tchannel "<<ch<<" is a boarder channel of Area "<<area<<endl;
			return true;
		}
	}
	if (verbosity>5)
		cout<<"channel "<<ch<<" is NOT a boarder channel "<<endl;
	return false;
}

bool TSettings::isMaskedCluster(UInt_t det, TCluster cluster,bool checkAdjacentChannels){
	bool isMasked = false;
	for(UInt_t i=0;i<cluster.getClusterSize()&&!isMasked;i++){
		int channelNo = cluster.getChannel(i);
		bool isScreened = this->isDet_channel_screened(det,channelNo);
		bool isAdjacentToCluster = !cluster.isHit(i);
		if(checkAdjacentChannels)
			isMasked = isMasked || isScreened;
		else if (isAdjacentToCluster){
				continue;
			}
			else{
				isMasked = isMasked || isScreened;
			}
	}
//	cout<<"==>"<<isMasked<<endl;
	return isMasked;
}

bool TSettings::hasBorderSeed(UInt_t det, TCluster cluster){
	UInt_t clSize = cluster.getClusterSize();
	if(verbosity>10)cout<<"\n\ncheck has Border seed "<<det<<"  -- "<<clSize<<endl;
	if(verbosity>10)cluster.Print(1);
	for (UInt_t clPos=0;clPos < clSize;clPos++){
		UInt_t ch = cluster.getChannel(clPos);
//		cout<<"\nCheck ClPos "<<clPos<<"/"<<clSize<<", channel "<<ch<<flush;
		if (cluster.isSeed(clPos)){
//			cout<<", isSeed, "<<flush;
			if(TPlaneProperties::isSiliconDetector(det)){
				if( ch == 0 || ch == TPlaneProperties::getNChannels(det)-1)
					return true;
			}
			else if (TPlaneProperties::isDiamondDetector(det)){
				if(isDiaDetectorAreaBorderChannel(ch))
					return true;
			}
		}
		else
			cout<<flush;
	}
	return false;
}
bool TSettings::hasBorderHit(UInt_t det, TCluster cluster){
	UInt_t clSize = cluster.getClusterSize();
	if(verbosity>10) cout<<"\n\ncheck has Border seed "<<det<<"  -- "<<clSize<<endl;
	if(verbosity>10) cluster.Print(1);
	for (UInt_t clPos=0;clPos < clSize;clPos++){
		UInt_t ch = cluster.getChannel(clPos);
//		cout<<"\nCheck ClPos "<<clPos<<"/"<<clSize<<", channel "<<ch<<flush;
		if (cluster.isHit(clPos)){
//			cout<<", isSeed, "<<flush;
			if(TPlaneProperties::isSiliconDetector(det)){
				if( ch == 0 || ch == TPlaneProperties::getNChannels(det)-1)
					return true;
			}
			else if (TPlaneProperties::isDiamondDetector(det)){
				if(isDiaDetectorAreaBorderChannel(ch))
					return true;
			}
		}
		//for symmetry reasons do not accept Seeds which are next to a boarder channel
		else if(cluster.isSeed(clPos)){
			if(TPlaneProperties::isSiliconDetector(det)){
				if( ch == 1|| ch == TPlaneProperties::getNChannels(det)-2)
					return true;
			}
			else if (TPlaneProperties::isDiamondDetector(det)){
				if(isDiaDetectorAreaBorderChannel(ch-1) || isDiaDetectorAreaBorderChannel(ch+1))
					return true;
			}
		}
		else
			cout<<flush;
	}
	return false;
}

Float_t TSettings::getPitchWidth(UInt_t det, UInt_t area){
	if (TPlaneProperties::isDiamondDetector(det) ){
		if (this->getNDiaDetectorAreas()>area)
			return diamondPattern.getPitchWidth(area);
		else
			return getDiamondPitchWidth();
	}
	else
		return this->getSiliconPitchWidth();
}
Float_t TSettings::convertChannelToMetric(UInt_t det, Float_t channel){
	if (TPlaneProperties::isDiamondDetector(det))
		return diamondPattern.convertChannelToMetric(channel);
	return channel*getSiliconPitchWidth();
}

Float_t TSettings::convertMetricToChannelSpace(UInt_t det, Float_t metricValue){
	Float_t channelPosition  = N_INVALID;
	channelPosition =  metricValue/this->getSiliconPitchWidth();
	if(TPlaneProperties::isDiamondDetector(det)){
			channelPosition = this->diamondPattern.convertMetricToChannel(metricValue);
	}
	return channelPosition;
}
