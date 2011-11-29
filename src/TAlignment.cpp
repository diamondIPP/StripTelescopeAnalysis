/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(int runNumber) {
	cout<<"**********************************************************"<<endl;
	cout<<"*************TAlignment::TAlignment***********************"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<runNumber<<endl;
	// TODO Auto-generated constructor stub
	sys = gSystem;
	stringstream settingsFilePath;
	settingsFilePath<<sys->pwd()<<"/Settings."<<runNumber<<".ini";
	cout<<settingsFilePath.str()<<endl;
	settings = new TSettings("settingsFilePath.str()");

	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());
	gErrorIgnoreLevel=1001;
	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<"clusterData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str());
	eventReader->checkADC();
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("alignment");
	sys->cd("alignment");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	initialiseHistos();
	cout<<"end initialise"<<endl;
	alignmentPercentage=1;
	Float_t stripSize= 50./10000.;//mu m
	detectorD0Z = 0.725/stripSize; // by definition in cm
	detectorD1Z = 1.625/stripSize; // by definition in cm
	detectorD2Z = 18.725/stripSize; // by definition in cm
	detectorD3Z = 19.625/stripSize; // by definition in cm
	detectorDiaZ = 10.2/stripSize; // by definition in cm
	this->runNumber=runNumber;
	verbosity=2;
	nAlignSteps=1;

	align=NULL;

}

TAlignment::~TAlignment() {
	// TODO Auto-generated destructor stub
	saveHistos();
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TAlignment::createVectors(){
	createVectors(eventReader->GetEntries());
}
void TAlignment::createVectors(UInt_t nEvents){
	int noHitDet=0;
	int falseClusterSizeDet=0;
	int noHitDia=0;
	int falseClusterSizeDia=0;
	int nCandidates=0;
	int nScreened=0;
	bool isScreened=false;
	cout<<"ANALYSE VECTORS...."<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		//check if every plane has exactly one cluster
		bool candidate=true;
		for(UInt_t det=0;det<8&&candidate;det++){
			if(eventReader->getCluster()->at(det).size()!=1){
				candidate=false;
				noHitDet++;
				break;
			}
			if(candidate){
				TCluster cluster = eventReader->getCluster()->at(det).at(0);
				if(cluster.size()>2){
					candidate=false;
					falseClusterSizeDet++;
					break;
				}
			}
		}//for det
		if(candidate&&eventReader->getCluster()->at(8).size()!=1){
			candidate=false;
			//cout<<"dia size:"<<eventReader->getCluster()->at(8).size()<<endl;
			noHitDia++;
		}

		//events with candidate=true areevents which have exactly one cluster in each plane
		// and the cluster size is 2
		if (candidate){
			isScreened=eventReader->getCluster()->at(8).at(0).isScreened();
			if(isScreened){
				nScreened++;
				candidate=false;
			}
		}
		if (candidate){
			nCandidates++;
			addEventToTracks();
		}
	}
	cout<<"\n\nDetAnalysed "<<nEvents<<": "<<nEvents-noHitDet-falseClusterSizeDet<<" Candidates while "<< noHitDet<<" Events have not exactly one Cluster and "<<falseClusterSizeDet<<" Events have wrong cluster size"<<endl;
	cout<<"\n\nDiaAnalysed "<<nEvents-noHitDet-falseClusterSizeDet<<": "<<nCandidates<<" Candidates while "<< noHitDia<<" Events have not exactly one Cluster and "<<falseClusterSizeDia<<" Events have wrong cluster size"<<endl;
	cout<<"EVENTS SCREENED:"<<nScreened<<endl;
	cout<<tracks.size()<<" "<<tracks_masked.size()<<" "<<tracks_fidcut.size()<<" "<<tracks_masked_fidcut.size()<<endl;

	for(UInt_t trackNo=0;trackNo<tracks.size();trackNo++){
		Float_t xPos=tracks.at(trackNo).GetDetectorHitPosition(0);
		Float_t yPos=tracks.at(trackNo).GetDetectorHitPosition(1);
		cout<<trackNo<<" "<<tracks.at(trackNo).GetEventNumber()<<flush;
		for(int no=1;no<4;no++){
			Float_t deltaX = tracks.at(trackNo).GetDetectorHitPosition(no*2)-xPos;
			cout<<" "<<deltaX;
			hXPositionDifference[no-1]->Fill(deltaX);
			hXXPositionDifference[no-1]->Fill(deltaX,tracks.at(trackNo).GetDetectorHitPosition(no*2));
			hXYPositionDifference[no-1]->Fill(deltaX,tracks.at(trackNo).GetDetectorHitPosition(no*2+1));
		}
		cout<<"\ty:  ";
		for(int no=1;no<4;no++){
			Float_t deltaY = tracks.at(trackNo).GetDetectorHitPosition(no*2+1)-yPos;
			cout<<" "<<deltaY;
			hYPositionDifference[no-1]->Fill(deltaY);
			hYXPositionDifference[no-1]->Fill(deltaY,tracks.at(trackNo).GetDetectorHitPosition(no*2));
			hYYPositionDifference[no-1]->Fill(deltaY,tracks.at(trackNo).GetDetectorHitPosition(no*2+1));
		}
		cout<<endl;
	}
}

void TAlignment::initialiseHistos(){
	for(int no = 0;no <4;no++){
		stringstream histName;
		histName<<"ScatterPlot_VectorCreation_8Hits_Plane_"<<no;
		this->hScatterPosition[no] = new TH2F(histName.str().c_str(),histName.str().c_str(),256,0,255,256,0,255);
	}
	for(int no=0;no<3;no++){
		//*******************XPOS******************************************
		stringstream histName;
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no+1)*2);
		this->hXPositionDifference[no]=new TH1F(histName.str().c_str(),histName.str().c_str(),2048,-256,256);

		histName.str("");
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no)*2+2)<<"_vsX";
		this->hXXPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		histName.str("");
		histName<<"xPos_Difference_D0X-"<<TADCEventReader::getStringForPlane((no)*2+2)<<"_vsY";
		this->hXYPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		//*********************YPOS***********************************************

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3);
		this->hYPositionDifference[no]=new TH1F(histName.str().c_str(),histName.str().c_str(),2048,-256,256);

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3)<<"_vsX";
		this->hYXPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

		histName.str("");
		histName<<"yPos_Difference_D0Y-"<<TADCEventReader::getStringForPlane((no)*2+3)<<"_vsY";
		this->hYYPositionDifference[no]=new TH2F(histName.str().c_str(),histName.str().c_str(),2048,-256,256,2048,0,256);

	}
}

void TAlignment::saveHistos(){

	for(int no=0;no<4;no++){
		this->histSaver->SaveHistogramPNG(this->hScatterPosition[no]);
		delete this->hScatterPosition[no];
	}
	for(int no=0;no<3;no++){
		//** XPOS***
		float mean = hXPositionDifference[no]->GetMean();
		float sigma= hXPositionDifference[no]->GetRMS();
		hXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		this->histSaver->SaveHistogramPNG(hXPositionDifference[no]);
		delete hXPositionDifference[no];
		hXXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXXPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		hXYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		this->histSaver->SaveHistogramPNG(hXXPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hXYPositionDifference[no]);
		delete this->hXXPositionDifference[no];
		delete this->hXYPositionDifference[no];

		//** YPOS***

		mean = hYPositionDifference[no]->GetMean();
		sigma= hYPositionDifference[no]->GetRMS();
		hYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hYXPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hYYPositionDifference[no]->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);
		hXYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		hYYPositionDifference[no]->GetYaxis()->SetRangeUser(50,170);
		this->histSaver->SaveHistogramPNG(hYPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hYXPositionDifference[no]);
		this->histSaver->SaveHistogramPNG(hYYPositionDifference[no]);
		delete hYPositionDifference[no];
		delete this->hYXPositionDifference[no];
		delete this->hYYPositionDifference[no];
	}
}

void TAlignment::addEventToTracks()
{
	TDetectorPlane D0;
	TDetectorPlane D1;
	TDetectorPlane D2;
	TDetectorPlane D3;
	TDetectorPlane Dia;
	D0.SetZ(detectorD0Z);
	D1.SetZ(detectorD1Z);
	D2.SetZ(detectorD2Z);
	D3.SetZ(detectorD3Z);
	Dia.SetZ(detectorDiaZ);

	D0.SetX(eventReader->getCluster()->at(0).at(0).getPosition());
	D1.SetX(eventReader->getCluster()->at(2).at(0).getPosition());
	D2.SetX(eventReader->getCluster()->at(4).at(0).getPosition());
	D3.SetX(eventReader->getCluster()->at(4).at(0).getPosition());

	D0.SetY(eventReader->getCluster()->at(1).at(0).getPosition());
	D1.SetY(eventReader->getCluster()->at(3).at(0).getPosition());
	D2.SetY(eventReader->getCluster()->at(5).at(0).getPosition());
	D3.SetY(eventReader->getCluster()->at(7).at(0).getPosition());

	Dia.SetX(eventReader->getCluster()->at(8).at(0).getPosition());

	TDiamondTrack newDiamondTrack(nEvent,D0,D1,D2,D3,Dia);

	Float_t xPos=eventReader->getCluster()->at(0).at(0).getPosition();
	Float_t yPos=eventReader->getCluster()->at(1).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(0,xPos);
	newDiamondTrack.SetDetectorHitPosition(1,yPos);
	hScatterPosition[0]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(2).at(0).getPosition();
	yPos=eventReader->getCluster()->at(3).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(2,xPos);
	newDiamondTrack.SetDetectorHitPosition(3,yPos);
	hScatterPosition[1]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(4).at(0).getPosition();
	yPos=eventReader->getCluster()->at(5).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(4,xPos);
	newDiamondTrack.SetDetectorHitPosition(5,yPos);
	hScatterPosition[2]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(6).at(0).getPosition();
	yPos=eventReader->getCluster()->at(7).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(6,xPos);
	newDiamondTrack.SetDetectorHitPosition(7,yPos);
	hScatterPosition[3]->Fill(xPos,yPos);

	xPos=eventReader->getCluster()->at(8).at(0).getPosition();
	newDiamondTrack.SetDetectorHitPosition(8,xPos);

	tracks.push_back(newDiamondTrack);
	tracks_masked.push_back((rand.Rndm()>this->alignmentPercentage));
	tracks_fidcut.push_back(newDiamondTrack);
	tracks_masked_fidcut.push_back((rand.Rndm()>this->alignmentPercentage));

}


int TAlignment::Align()
{
	if(tracks.size()==0)
		createVectors();

	if (verbosity){
		cout<<"\n\n\nTAlignment::Align:Starting \""<<histSaver->GetPlotsPath()<<"\""<<endl;
		cout<<"\t\t"<<tracks.size()<<" "<<tracks_masked.size()<< " ";
		cout<<		  tracks_fidcut.size()<<" "<<tracks_masked_fidcut.size()<<endl;
		cout << "\t\t "<<eventReader<<" ."<<endl;
	}
	if(tracks.size()==0) {
		cout<<"TAlignment::Align: No tracks found; need to CallClustering::ClusterRun first..."<<endl;
		return -1;
		//ClusterRun(plots); // doesn't use alternative clustering
	}

	if (tracks.size() == 0) {
		cout << "TAlignment::Align:No tracks available. Alignment not possible. (tracks.size() = " << tracks.size() << ")" << endl;
		return 0;
	}

	if(align==NULL){
		align = new TDetectorAlignment(histSaver->GetPlotsPath(), tracks, tracks_masked);
		cout<<"TAlignment::Align::Detectoralignment did not exist, so created new DetectorAlignment"<<endl;
	}

	// now start the telescope alignment!
	// alignment loop: align, cut fake tracks, align again (if CutFakeTracksOn is set true)
	for (int alignStep = 0; alignStep < nAlignSteps; alignStep++) {
		//TDetectorAlignment* align = new TDetectorAlignment(plots_path, tracks, tracks_mask);
		align->LoadTracks(tracks,tracks_masked);
		if (verbosity) cout<<"TAlignment::Align:start with alignmentStep no. "<<alignStep+1 << " of " <<nAlignSteps<<endl;
		doDetAlignmentStep();
		cout << "AlignmentClass::Align:Intrinsic silicon resolution " << align->GetSiResolution() << " strips or " << align->GetSiResolution() * 50 << "um" << endl;
		doDiaAlignmentStep();
	}
	return 1;
}


void TAlignment::doDetAlignmentStep(){
	Int_t nPasses = 10;
	//Double_t plot_width_factor = 3; // scales the widths of the plots; range is a 3*width of distribution centered on mean
	if(verbosity)cout<<"TAlignment::doAlignmentStep::PlotAngularDistribution"<<endl;

	align->PlotAngularDistribution(); //look at angular distribution of tracks
	align->PlotCartesianDistribution(); //look at slope distribution of tracks

	string prename = "alignment_PrealignmentResiduals";
	string postname = "alignment_PostalignmentResiduals";

	// generate residuals before alignment
	if(verbosity)cout<<"TAlignment::Align::CheckDetectorAlignment"<<endl;
	//align->LoadTracks(this->tracks, this->tracks_masked);
	align->CheckDetectorAlignmentXYPlots(0, 1, 3, prename);
	align->CheckDetectorAlignmentXYPlots(1, 0, 3, prename);
	align->CheckDetectorAlignmentXYPlots(2, 0, 3, prename);
	align->CheckDetectorAlignmentXYPlots(3, 0, 2, prename);
	// itterative alignment loop
	if(verbosity)cout<<"TAlignment::doDetAlignmentStep::alignmentloop"<<endl;

	align->setVerbosity(verbosity-1);
	for(int i=0; i<nPasses; i++) {
		cout << "\n\nPass " << i+1 << endl<< endl;
		//xy alignment
		align->CheckDetectorAlignmentXY(0, 1, 3);
		align->AlignDetectorXY(1, 0, 3,verbosity>1);
		align->AlignDetectorXY(2, 0, 3,verbosity>1);
		//align->AlignDetectorXY(1, 0, 3);
		//align->AlignDetectorXY(2, 0, 3);
		//align->AlignDetectorXY(1, 0, 3);
		//align->AlignDetectorXY(2, 0, 3);
		align->AlignDetectorXY(3, 0, 2,verbosity>1);
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
	cout<<"TAlignment::doDetAlignmentStep:Checking final Silicon residuals"<<endl;
	cout<<endl;

	// generate residuals after alignment
	align->CheckDetectorAlignmentXYPlots(0, 1, 3, postname);
	align->CheckDetectorAlignmentXYPlots(1, 0, 3, postname);
	align->CheckDetectorAlignmentXYPlots(2, 0, 3, postname);
	align->CheckDetectorAlignmentXYPlots(3, 0, 2, postname);

	cout<<endl;
}

void TAlignment::doDiaAlignmentStep()
{
	align->LoadTracks(tracks_fidcut, tracks_masked_fidcut);

	//check that the silicon is still aligned for these tracks_fidcut
	cout<<"TAlignment::doDiaAlignmentStep:Check that the telescope alignment still holds for fidcut tracks w/ single diamond cluster"<<endl;
	align->CheckDetectorAlignmentXY(0, 1, 3);
	align->CheckDetectorAlignmentXY(1, 0, 3);
	align->CheckDetectorAlignmentXY(2, 0, 3);
	align->CheckDetectorAlignmentXY(3, 0, 2);

	string prename = "alignment_PrealignmentResiduals";
	string postname = "alignment_PostalignmentResiduals";

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
	cout<<"TAlignment::doDiaAlignmentStep:Checking final diamond residuals"<<endl;
	cout<<endl;

	// generate residuals after alignment
	align->CheckDetectorAlignmentXYPlots(4, 1, 2, postname);

}



