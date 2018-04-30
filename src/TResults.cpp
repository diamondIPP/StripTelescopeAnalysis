/*
 * TResults.cpp
 *
 *  Created on: May 29, 2012
 *      Author: bachmair
 */

#include "../include/TResults.hh"

ClassImp(TResults);

using namespace std;

TResults::TResults(UInt_t runnumber){
    path = gSystem->pwd();
    this->runnumber = runnumber;
    TString name =  TString::Format("results_%d",this->runnumber);
    this->SetName(name);
//    initialiseResults();
    rootFileName = TString::Format("results.%d.root",this->runnumber);
}

TResults::TResults(TSettings *settings) {
    path = gSystem->pwd();
    initialiseResults(settings);
    setResultsFromSettings(settings);
    cout<<"New TResults with settings "<<settings;
    if (settings) cout << "\t run: "<<settings->getRunNumber()<<endl;
    else cout<<"\t run: XXXXXXXX "<<endl;
    openResults(settings);
    textFileName = settings->getAbsoluteOuputPath(true);
    textFileName.Append(TString::Format("/results_%d.txt",this->runnumber));

    resultsFileName = settings->getAbsoluteOuputPath(true);
    resultsFileName.Append(TString::Format("/results_%d.res",this->runnumber));
    if (settings)
        rootFileName = settings->getResultsRootFilePath();
    else
        rootFileName = TString::Format("results.%d.root",this->runnumber);
    //  this->settings=settings;
}

TResults::~TResults() {
    // TODO Auto-generated destructor stub
    cout<<"delete Results of run "<<runnumber<<endl;

    writeFiles();
}

TResults::TResults(const   TResults& rhs){//copy constructor
//    initialiseResults();
    inheritOldResults(rhs);
}
//
//TResults::TResults &operator=(const   TResults &src){ //class assignment function
//
//}


void TResults::inheritOldResults(const TResults & rhs)
{
    cout<<"InheritOldResults"<<endl;
    this->seedSigma.clear();
    for(UInt_t det=0;det<rhs.seedSigma.size();det++)this->seedSigma.push_back(rhs.seedSigma[det]);
    this->hitSigma.clear();
    for(UInt_t det=0;det<rhs.hitSigma.size();det++)this->hitSigma.push_back(rhs.hitSigma[det]);
    this->noise.clear();
    for(UInt_t det=0;det<rhs.noise.size();det++)this->noise.push_back(rhs.noise[det]);
    signalFeedOverCorrection.clear();
    for(UInt_t det=0;det<rhs.signalFeedOverCorrection.size();det++)
        signalFeedOverCorrection.push_back(signalFeedOverCorrection[det]<-1?-1:signalFeedOverCorrection[det]);

    diaCMCNoise = rhs.diaCMCNoise;
    CMN = rhs.CMN;
    mean2outOf10_normal = rhs.mean2outOf10_normal;
    mean2outOf10_trans = rhs.mean2outOf10_trans;
    meanNoutOfN_normal.clear();
    for(UInt_t cl=0;cl<rhs.meanNoutOfN_normal.size();cl++)this->meanNoutOfN_normal.push_back(rhs.meanNoutOfN_normal[cl]);
    meanNoutOfN_trans.clear();
    for(UInt_t cl=0;cl<rhs.meanNoutOfN_trans.size();cl++)this->meanNoutOfN_trans.push_back(rhs.meanNoutOfN_trans[cl]);

    mp2outOf10_normal = rhs.mp2outOf10_normal;
    mp2outOf10_trans = rhs.mp2outOf10_trans;
    width2outOf10_normal = rhs.width2outOf10_normal;
    width2outOf10_trans = rhs.width2outOf10_trans;
    gSigma2outOf10_normal = rhs.gSigma2outOf10_normal;
    gSigma2outOf10_trans = rhs.gSigma2outOf10_trans;
    this->lastUpdate = rhs.lastUpdate;

    mean_clustered_normal = rhs.mean_clustered_normal;
    mp_clustered_normal = rhs.mp_clustered_normal;
    width_clustered_normal = rhs.width_clustered_normal;
    gSigma_clustered_normal = rhs.gSigma_clustered_normal;

    mean_clustered_trans = rhs.mean_clustered_trans;
    mp_clustered_trans = rhs.mp_clustered_trans;
    width_clustered_trans = rhs.width_clustered_trans;
    gSigma_clustered_trans = rhs.gSigma_clustered_trans;
    repeaterCard = rhs.repeaterCard;
    maskedChannels.clear();
    for (set<Int_t>::iterator it = rhs.maskedChannels.begin(); it !=  rhs.maskedChannels.end();it++)
        maskedChannels.insert(*it);

    noisyChannels.clear();
    for (set<Int_t>::iterator it = rhs.noisyChannels.begin(); it !=  rhs.noisyChannels.end();it++)
        noisyChannels.insert(*it);;
    notConnectedChannels.clear();
    for (set<Int_t>::iterator it = rhs.notConnectedChannels.begin(); it !=  rhs.notConnectedChannels.end();it++)
        notConnectedChannels.insert(*it);;
//    std::sort(maskedChannels.begin(),maskedChannels.end());

    IntegerMap.clear();
    cout<<"IntegerMap.size() "<<rhs.IntegerMap.size()<<endl;
    for(map<TString,Int_t>::const_iterator it = rhs.IntegerMap.begin();it!=rhs.IntegerMap.end();it++)
        IntegerMap[it->first] = it->second;

    FloatMap.clear();
    cout<<"FloatMap.size(): "<<rhs.FloatMap.size()<<endl;
    for(map<TString,Float_t>::const_iterator it = rhs.FloatMap.begin();it!=rhs.FloatMap.end();it++)
        FloatMap[it->first] = it->second;

    StringMap.clear();
    cout<<"FloatMap.size()"<<rhs.FloatMap.size()<<endl;
    for(map<TString,TString>::const_iterator it = rhs.StringMap.begin();it!=rhs.StringMap.end();it++)
        StringMap[it->first] = it->second;

    cout<<"initOldResults: Copy keyList"<<endl;
    keyList.clear();

    cout<<"rhs.keyList.size(): "<<rhs.keyList.size()<<endl;
    int i =0;
    map<TString, map<TString,TString> >::const_iterator it1;
    int k = 0;
    for(it1 = rhs.keyList.begin();it1!=rhs.keyList.end()&&i<rhs.keyList.size();it1++){
        cout<<i<<endl;
        TString section = it1->first;
        cout<<"Create Section:\t\""<<section<<"\""<<endl;
        cout<<section.IsAlnum()<<" "<<section.IsAscii()<<endl;
        i++;
        if (!section.IsAscii()){cout<<"continue"<<endl; continue;}
            keyList[section] = map<TString,TString>();
            map<TString,TString>::const_iterator it2;
            int size = it1->second.size();
            int j =0;
            for (it2 = (it1->second).begin(); it2!=(it1->second).end() && j < size;it2++){
                TString key = it2->first;
                cout<<"Add "<<section<<"\t\t\""<<key<<"\""<<endl;
                if (!key.IsAscii()){cout<<"continue"<<endl; continue;}
                keyList[section][key] = it2->second;
                j++;
                k++;
            }
        }
    cout<<"added "<<k<< " keys in "<< i << " section..."<<endl;
}

void TResults::initialiseResults(TSettings *settings){
    repeaterCard = -1;
    runnumber = -1;
    seedSigma.resize(TPlaneProperties::getNDetectors(),-1);
    hitSigma.resize(TPlaneProperties::getNDetectors(),-1);
    noise.resize(TPlaneProperties::getNDetectors(),-1);
//    meanNoutOfN_normal.resize(TPlaneProperties::getMaxTransparentClusterSize(TPlaneProperties::getDetDiamond()),-1); // DA
    meanNoutOfN_normal.resize(settings->getMaxTransparentClusterSize()); // DA
//    meanNoutOfN_trans.resize(TPlaneProperties::getMaxTransparentClusterSize(TPlaneProperties::getDetDiamond()),-1); // DA
    meanNoutOfN_trans.resize(settings->getMaxTransparentClusterSize()); // DA
    doubleGaus1_normal = -1;
    doubleGaus2_normal = -1;
    doubleGaus1_trans = -1;
    doubleGaus2_trans = -1;
    singleGausFWTM_normal = -1;
    singleGausFWTM_trans = -1;

    singleGausShort_normal = -1;
    singleGausShort_trans = -1;
    singleGaus_normal = -1;
    singleGaus_trans = -1;
    signalFeedOverCorrection.resize(TPlaneProperties::getNDetectors(),-1.0);

    mean_clustered_normal = -1;
    mp_clustered_normal = -1;
    width_clustered_normal = -1;
    gSigma_clustered_normal = -1;

    mean_clustered_trans = -1;
    mp_clustered_trans = -1;
    width_clustered_trans = -1;
    gSigma_clustered_trans = -1;

    diaCMCNoise = -1;
    CMN = -1;
    nAllEvents = -1;
    nNoSiliconHit = -1;
    nOneAndOnlyOneSiliconNotFiducialCut = -1;
    nValidSiliconTrack =-1;
    nNoDiamondHit =-1;
    nExactlyOneDiamondHit =-1;
    nUseForAlignment =-1;
    nUseForAnalysis = -1;
    maskedChannels.clear();
    noisyChannels.clear();
    notConnectedChannels.clear();
}


void TResults::openResults(TSettings *settings){
    cout<< settings<<endl;
    //	if (!settings)return;
    cout<<"open Results with settings "<<settings<<"\t run: "<<settings->getRunNumber()<<flush;
    rootFileName = settings->getResultsRootFilePath();
    //  this->Settings = *settings;
    runnumber = settings->getRunNumber();
    cout<<" "<<runnumber<<endl;
    //	std::stringstream resultsFile;
    //	resultsFile<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
    cout<<((string)(rootFileName))<<endl;

    TFile *file =  new TFile(rootFileName,"READ");
    TResults *oldResults;

    if(file->IsZombie()) {
        cout << "FIle does not exists, create new File!"<<endl;
        delete file;
        oldResults = new TResults(runnumber);
    }
    else{
        file->GetListOfKeys()->Print();
        stringstream name;
        name << "results_"<<runnumber;
        cout<<"Name of key: \""<<name.str()<<"\""<<endl;
        TIter next(file->GetListOfKeys());
        TKey *key;
        oldResults = 0;
        bool foundResult = false;
        while ((key=(TKey*)next())) {
            TString className = key->GetClassName();
            if (className.Contains("TResults")){
                if (foundResult == true && oldResults){
                    cout <<" found a second Results key in the file..."<<oldResults->GetName()<<" "<<key->GetName()<<endl;
                }
                else{
                    cout<<"FOUND a valid class: "<<key->GetName()<<endl;
                    oldResults =  (TResults*)key->ReadObj();
                    foundResult = true;
                }
            }
        }
        //		oldResults = (TResults*)file->FindGet(name.str().c_str());
        cout<<"old Results: "<<oldResults<<flush;
        if (oldResults) cout<<" "<<oldResults->GetName()<<flush;

        if(oldResults==0){
            cerr<< "Something is wrong, results does not exists..."<<endl;
            return;
        }
        cout<<oldResults->IsZombie()<<endl;
        //		cout<<"LAST UPDATE ON "<<oldResults->getLastUpdateDate().AsString()<<endl;
        this->inheritOldResults(*oldResults);
        setResultsFromSettings(settings);

    }
    Print();
}
TString TResults::getChannelsStringList(std::set<Int_t> channelSet){
    TString output = "[";
    for (set<Int_t>::iterator it = channelSet.begin(); it !=  channelSet.end();it++)
            output.Append(TString::Format("%d, ",*it));
    //    }
        output = output.Strip(TString::kBoth,' ');
        output = output.Strip(TString::kBoth,',');
        output.Append("]");
        return output;
}
TString TResults::getMaskedChannels() {
    return getChannelsStringList(maskedChannels);
}

//TString TResults::getMaskedChannels() {
//    TString output = "[";
////    for (UInt_t i  = 0; i< maskedChannels.size();i++){
//
//    for (set<Int_t>::iterator it = maskedChannels.begin(); it !=  maskedChannels.end();it++)
//        output.Append(TString::Format("%d, ",*it));
////    }
//    output = output.Strip(TString::kBoth,' ');
//    output = output.Strip(TString::kBoth,',');
//    output.Append("]");
//    return output;
//}

void TResults::setResultsFromSettings(TSettings* settings){
    updated();;
    //	cout<<"setResultsFromSettings"<<endl;
    hitSigma.resize(TPlaneProperties::getNDetectors(),-1);
    seedSigma.resize(TPlaneProperties::getNDetectors(),-1);
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
        seedSigma.at(det)=settings->getClusterSeedFactor(det,0);
        hitSigma.at(det)=settings->getClusterHitFactor(det,0);
    }
    runDescription = settings->getRunDescription();
    UInt_t diaDet = TPlaneProperties::getDetDiamond();
    for (UInt_t ch=0; ch< TPlaneProperties::getNChannels(diaDet);ch++){
        if (settings->IsMasked(diaDet,ch))
            maskedChannels.insert(ch);
        if (settings->IsNoisyChannel(ch))
            noisyChannels.insert(ch);
        if (settings->IsNotConnectedChannel(ch))
            notConnectedChannels.insert(ch);
    }
    std::string rundes = settings->getRunDescription();
    std::cout<<"get diamondChannels: "<<rundes<<std::endl;
    diamondChannels = settings->diamondPattern.getIntervalOfDiamond(rundes);
    repeaterCard = settings->getRepeaterCard();
    diamondName = settings->getDiamond();
    voltage = settings->getVoltage();
    StringMap["currentBegin"] = settings->GetCurrentBegin();
    StringMap["currentEnd"] = settings->GetCurrentEnd();
//    std::sort(maskedChannels.begin(),maskedChannels.end());
    //	cout<<runDescription<<endl;
    //	char t;
    //	cin>>t;
}

void TResults::saveResults(TString name){
    updated();;
    //	std::stringstream fileName;
    //	fileName<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
    //	cout<<"save File:"<<fileName.str()<<endl;
    //	std::stringstream name;
    //	name <<"results_"<<runnumber;
    //	this->SetName(name.str().c_str());
    //	TFile *file =  new TFile(fileName.str().c_str(),"RECREATE");
    //	file->cd();
    //	this->Write();
    //	this->Write(fileName.str().c_str());
    cout<< "Save "<<name<<endl;
    //	this->Write(name);
    //	TFile *file =  new TFile(name.Append(".bak"),"RECREATE");
    //	this->Write();
    //	file->Close();
    writeFiles();
}

void TResults::Print(){
    //	getLastUpdateDate().Print();
    //	cout<<"\t"<<�getLastUpdateDate().AsString()<<endl;
    cout<<getLastUpdateDate().AsString()<<endl;
    cout<<"\tdet\tseed  \thit  \tnoise"<<endl;
    for(UInt_t det=0; det<TPlaneProperties::getNDetectors();det++){
        cout<<"\t"<<det<<"\t"<<std::setw(4)<<this->seedSigma[det]<<"\t"<<std::setw(4)<<this->hitSigma[det]<<" \t"<<std::setw(4)<<this->noise[det]<<endl;
    }
    //	alignment.Print();

    //  settings->Print();z
    cout << "\tmean2outOf10_normal   = " << mean2outOf10_normal << "\n"
            << "\tmp2outOf10_normal     = " << mp2outOf10_normal << "\n"
            << "\tmp2outOf10_normal     = " << mp2outOf10_normal << "\n"
            << "\tgSigma2outOf10_normal = " << gSigma2outOf10_normal<<endl;
    cout<<"DONE"<<endl;
    writeFiles();
}

void TResults::setSignalFeedOverCorrection(UInt_t det, Float_t correction){
    updated();;
    if(det>=TPlaneProperties::getNDetectors())
        return;
    if(signalFeedOverCorrection.size() <= det)
        signalFeedOverCorrection.resize(det+1 , -1);
    signalFeedOverCorrection[det] = correction;
    if (false) cout<<"Set Results: signal feed over correction of det "<<det<<": "<<correction*100<<" %"<<endl;
}

void TResults::setNoise(UInt_t det,Float_t detNoise){
    updated();;
    if(det>=TPlaneProperties::getNDetectors())
        return;
    if(noise.size()<=det)
        noise.resize(det+1,0);
    noise[det]=detNoise;
    if (false) cout<<"Set Results: Noise of det "<<det<<": "<<detNoise<<endl;
}


void TResults::setAlignment(TDetectorAlignment* newAlignment)
{
    updated();;
    this->alignment= *newAlignment;
    if (false) cout<<"TResults:SetAlignment"<<endl;
}


int TResults::getRunNumber() const {
    return runnumber;
}

void TResults::setRunNumber(UInt_t runNumber) {
    updated();;
    this->runnumber = runNumber;
}


void TResults::setDoubleGaussianResolution(Float_t gaus1,Float_t gaus2, TSettings::alignmentMode mode){
    updated();;
    if (gaus2 > gaus1){
        Float_t val = gaus1;
        gaus1 = gaus2;
        gaus2 = val;
    }
    if (mode == TSettings::normalMode){
        doubleGaus1_normal = gaus1;
        doubleGaus2_normal = gaus2;
    }
    else if ( mode == TSettings::transparentMode){
        doubleGaus1_trans = gaus1;
        doubleGaus2_trans = gaus2;
    }
}

void TResults::setSingleGaussianFixedResolution(Float_t gaus,TSettings::alignmentMode mode){
    updated();;
    if (mode == TSettings::normalMode){
        singleGausFixed_normal = gaus;
    }
    else if ( mode == TSettings::transparentMode){
        singleGausFixed_trans = gaus;
    }
}

void TResults::setSingleGaussianResolution(Float_t gaus,TSettings::alignmentMode mode){
    updated();;
    if (mode == TSettings::normalMode){
        singleGaus_normal = gaus;
    }
    else if ( mode == TSettings::transparentMode){
        singleGaus_trans = gaus;
    }
}
void TResults::setSingleGaussianShortResolution(Float_t gaus,TSettings::alignmentMode mode){
    updated();;
    if (mode == TSettings::normalMode){
        singleGausShort_normal = gaus;
    }
    else if ( mode == TSettings::transparentMode){
        singleGausShort_trans = gaus;
    }
}

void TResults::setSingleGaussianFWTMResolution(Float_t gaus,TSettings::alignmentMode mode){
    updated();;
    if (mode == TSettings::normalMode){
        singleGausFWTM_normal = gaus;
    }
    else if ( mode == TSettings::transparentMode){
        singleGausFWTM_trans = gaus;
    }
}


void TResults::setPH_clustered(Float_t mean, Float_t mp, Float_t width,
        Float_t gSigma, TSettings::alignmentMode mode) {
    updated();;
    if (false) cout<<"SET PH_clustered "<<mean<<" "<<mp<<" "<<width<<" "<<gSigma<<" "<<mode<<endl;
    if (mode == TSettings::transparentMode){
        mean_clustered_trans = mean;
        mp_clustered_trans = mp;
        width_clustered_trans = width;
        gSigma_clustered_trans = gSigma;
        if (false) cout<<"SET PH_clustered trans. "<<mean_clustered_trans<<" "<<mp_clustered_trans<<" "<<width_clustered_trans<<" "<<gSigma_clustered_trans<<endl;
    }
    else{
        mean_clustered_normal = mean;
        mp_clustered_normal = mp;
        width_clustered_normal = width;
        gSigma_clustered_normal = gSigma;
        if (false) cout<<"SET PH_clustered norm. "<<mean_clustered_normal<<" "<<mp_clustered_normal<<" "<<width_clustered_normal<<" "<<gSigma_clustered_normal<<endl;
    }
    writeFiles();
}

void TResults::setPH_2outOf10(Float_t mean, Float_t mp, Float_t width, Float_t gSigma, TSettings::alignmentMode mode){
    updated();;
    if (mode == TSettings::normalMode){
        mean2outOf10_normal = mean;
        mp2outOf10_normal = mp;
        width2outOf10_normal = width;
        gSigma2outOf10_normal = gSigma;
    }
    else if ( mode == TSettings::transparentMode){
        mean2outOf10_trans = mean;
        mp2outOf10_trans = mp;
        width2outOf10_trans = width;
        gSigma2outOf10_trans = gSigma;
    }
}

void TResults::setPH_NoutOfN(vector<Float_t> vecPHNoutOfN, TSettings::alignmentMode mode){
    updated();;
    if(meanNoutOfN_trans.size()<vecPHNoutOfN.size())
        meanNoutOfN_trans.resize(vecPHNoutOfN.size(),-1);
    if(meanNoutOfN_normal.size()<vecPHNoutOfN.size())
        meanNoutOfN_normal.resize(vecPHNoutOfN.size(),-1);

    if(mode == TSettings::normalMode){
        for (UInt_t i = 0; i< vecPHNoutOfN.size();i++)
            meanNoutOfN_normal.at(i) = vecPHNoutOfN.at(i);
    }
    else if (mode == TSettings::transparentMode){
        for (UInt_t i = 0; i< vecPHNoutOfN.size();i++)
            meanNoutOfN_trans.at(i) = vecPHNoutOfN.at(i);
    }

}

void TResults::setCMN(Float_t cmn){
    updated();;
    this->CMN = cmn;
}

Float_t TResults::getAvergDiamondCorrection(){
    UInt_t det = TPlaneProperties::getDetDiamond();
    if ( signalFeedOverCorrection.size()>det){
        return signalFeedOverCorrection[det];
    }
    return -100;
}
std::pair<Float_t, Float_t> TResults::getAvergSiliconCorrection(){
    UInt_t nSilDetectors = 0;
    Float_t mean = 0;
    Float_t sigma2 = 0;
    for( UInt_t det = 0; det < TPlaneProperties::getNDetectors() && det < signalFeedOverCorrection.size(); det++){
        if(TPlaneProperties::isSiliconDetector(det) && signalFeedOverCorrection[det] != -100){
            nSilDetectors ++;
            Float_t correction = signalFeedOverCorrection[det]*100.;
            mean += correction;
            sigma2 += correction * correction;
            //			cout<<"det: "<<correction*100.
            //					<<" "<<mean/nSilDetectors<<"+/-"<<sigma2<<" "<<correction*100.*correction*100.<<endl;
        }
    }
    mean = mean/(Float_t) nSilDetectors;
    sigma2 = sigma2/(Float_t) nSilDetectors;
    //	cout<< mean*1e2 <<" "<<sigma2*1e4<<endl;
    //	cout<< mean *mean *1e4<<" "<<sigma2*1e4<<endl;
    sigma2 = sigma2 - mean * mean;
    mean/=100.;
    sigma2/=100;
    if(mean>100 || mean < -100){
        mean = -1;
        sigma2 = -1;
    }

    if (false) cout << "avrg Sil Feed over correction: " << mean*100 << " +/- " << sigma2*100 <<" % in "<< nSilDetectors<< " Dectectors." << endl;
    return make_pair(mean,sigma2);
}


string TResults::emptyString(UInt_t nChars,char character){
    string str = "";
    str.resize(nChars,character);
    return str;
}

string TResults::createSection(TString sectionName,map<TString, TString> results) {
    stringstream output;
    output <<"["<<sectionName<<"]\n";
    map<TString,TString>::iterator it;
    Int_t length1 = 0;
    Int_t length2 = 0;
    for (it = results.begin();it!=results.end();it++){
        length1 = TMath::Max(it->first.Length(),length1);
        length2 = TMath::Max(it->second.Length(),length2);
    }
    for (it = results.begin();it!=results.end();it++){
        Int_t nChars1 = length1 - (it->first.Length());
        Int_t nChars2 = length2 - (it->second.Length());
        output<<it->first.Strip(TString::kBoth)<<":"<<emptyString(nChars1)<<"\t"<<it->second.Strip(TString::kBoth)<<"\n";
    }
    output<<"\n";
    return output.str();
}




void TResults::createOutputResultFile(){
    UInt_t det  = 8;
    ofstream myfile;
    myfile.open (resultsFileName, ios::out |ios::trunc);

    std::map<TString,TString> results;
    results["RunNo"] = TString::Format("%7d\t",runnumber);
    results["descr."] = runDescription;
    results["dia"] = diamondName;
    results["RepeaterCardNo"] = TString::Format("%d",repeaterCard);
    results["corrected"] = TString::Format("%d",(runnumber>1e5));
    if (voltage == 0)
        results["Voltage"] = "UNKOWN";
    else
        results["Voltage"] = TString::Format("%+d",voltage);
    results["SVN_REV"] = SVN_REV;
    results["lastUpdate"] = lastUpdate.AsString();
    results["maskedChannels"] =  getChannelsStringList(maskedChannels);
    results["noisyChannels"] = getChannelsStringList(noisyChannels);
    results["notConnectedChannels"] = getChannelsStringList(notConnectedChannels);
    results["channels_diamond"] = TString::Format("%d - %d ",diamondChannels.first,diamondChannels.second);
    if (StringMap.count("currentBegin"))
        results["currentBegin"] = StringMap["currentBegin"];
    else
        results["currentBegin"] = "??";
    if (StringMap.count("currentEnd"))
        results["currentEnd"] = StringMap["currentEnd"];
    else
        results["currentEnd"] = "??";
    myfile << createSection("RunInfo",results);

    results.clear();
    results["avrgNoise_sil"] =  TString::Format("%6.2f\t",getAvrgSilNoise().first);
    results["sigNoise_sil"] =   TString::Format("%6.2f\t",getAvrgSilNoise().second);
    results["Noise_dia"] =      TString::Format("%6.2f\t",noise[det]);
    results["CMC_Noise_dia"] =  TString::Format("%6.2f\t",diaCMCNoise);
    results["CM_Noise_dia"] =   TString::Format("%6.2f\t",CMN);
    myfile << createSection("Noise",results);

    results.clear();
    results["CorSil"] = TString::Format("%+6.2f\t",getAvergSiliconCorrection().first*100);
    results["sigSil"] = TString::Format("%+6.2f\t",getAvergSiliconCorrection().second*100);
    results["CorDia"] = TString::Format("%+6.2f\t",getAvergDiamondCorrection()*100);
    myfile << createSection("Feed_Through_Correction",results);

    results.clear();
    results["mean2outOf10_clustered"] =     TString::Format("%+7.2f\t",mean_clustered_normal);;
    results["mp2outOf10_clustered"] =       TString::Format("%+7.2f\t",mp_clustered_normal);
    results ["width2outOf10_clustered"] =   TString::Format("%+7.2f\t",width_clustered_normal);
    results ["gSigma2outOf10_clustered"] =  TString::Format("%+7.2f\t",gSigma_clustered_normal);
    myfile << createSection("Landau_clustered",results);

    results.clear();
    results["mean2outOf10_normal"] =    TString::Format("%+7.2f\t",mean2outOf10_normal);;
    results["mp2outOf10_normal"] =      TString::Format("%+7.2f\t",mp2outOf10_normal);
    results ["width2outOf10_normal"] =  TString::Format("%+7.2f\t",width2outOf10_normal);
    results ["gSigma2outOf10_normal"] = TString::Format("%+7.2f\t",gSigma2outOf10_normal);
    results["m2/2_normal"] =            TString::Format("%+7.2f\t",meanNoutOfN_normal[1]);
    results ["m4/4_normal"] =           TString::Format("%+7.2f\t",meanNoutOfN_normal[3]);
    myfile << createSection("Landau_normal",results);

    results.clear();
    results["mean2outOf10_trans"] =     TString::Format("%+7.2f\t",mean2outOf10_trans);;
    results["mp2outOf10_trans"] =       TString::Format("%+7.2f\t",mp2outOf10_trans);
    results ["width2outOf10_trans"] =   TString::Format("%+7.2f\t",width2outOf10_trans);
    results ["gSigma2outOf10_trans"] =  TString::Format("%+7.2f\t",gSigma2outOf10_trans);
    results["m2/2_trans"] =             TString::Format("%+7.2f\t",meanNoutOfN_trans[1]);
    results ["m4/4_trans"] =            TString::Format("%+7.2f\t",meanNoutOfN_trans[3]);
    myfile << createSection("Landau_trans",results);

    results.clear();
    results["singleGausFixed_normal"] = TString::Format("%8.2f\t",singleGausFixed_normal);
    results["singleGausShort_normal"] = TString::Format("%8.2f\t",singleGausShort_normal);
    results["singleGaus_normal"] =      TString::Format("%8.2f\t",singleGaus_normal);
    results["singleGausFWTM_normal"] =  TString::Format("%8.2f\t",singleGausFWTM_normal);
    results["doubleGaus1_normal"] =     TString::Format("%8.2f\t",doubleGaus1_normal);
    results["doubleGaus2_normal"] =     TString::Format("%8.2f\t",doubleGaus2_normal);
    myfile << createSection("Resolution_normal",results);

    results.clear();
    results["singleGausFixed_trans"] =  TString::Format("%8.2f\t",singleGausFixed_trans);
    results["singleGausShort_trans"] =  TString::Format("%8.2f\t",singleGausShort_trans);
    results["singleGaus_trans"] =       TString::Format("%8.2f\t",singleGaus_trans);
    results["singleGausFWTM_trans"] =   TString::Format("%8.2f\t",singleGausFWTM_trans);
    results["doubleGaus1_trans"] =      TString::Format("%8.2f\t",doubleGaus1_trans);
    results["doubleGaus2_trans"] =      TString::Format("%8.2f\t",doubleGaus2_trans);
    myfile << createSection("Resolution_trans",results);


    for (std::map<TString,map<TString,TString> >::iterator it = keyList.begin();it!=keyList.end();it++){
        results.clear();
        TString section = it->first;
        for (std::map<TString,TString>::iterator it2 = it->second.begin();it2!=it->second.end();it2++){
            TString key = it2->first;
            TString mapKey = section+(TString)"_"+key;
            TString type = it2->second;
            TString answer;
            if (type=="Float"){
                answer = TString::Format("%f",FloatMap[mapKey]);
            }
            else if(type=="Int"){
                answer = TString::Format("%d",IntegerMap[mapKey]);
            }
            else if(type == "String"){
                answer = StringMap[mapKey];
            }
            results[it2->first] = answer;
        }
        myfile << createSection(section,results);
    }

    results.clear();
    for (std::map<TString, Int_t>::iterator it = IntegerMap.begin();it!=IntegerMap.end();it++)
        if(it->first.First('_')<=0)
            results[it->first] = TString::Format("%d",it->second);
    for (std::map<TString, Float_t>::iterator it = FloatMap.begin();it!=FloatMap.end();it++)
        if(it->first.First('_')<=0)
            results[it->first] = TString::Format("%f",it->second);
    for (std::map<TString, TString>::iterator it = StringMap.begin();it!=StringMap.end();it++)
        if(it->first.First('_')<=0)
            results[it->first] = it->second;
    myfile << createSection("additional",results);
    myfile.close();
}

void TResults::createOutputTextFile(){
    //	cout<<"CREATE OUPTUT TEXT FILE"<<endl;
    ofstream myfile;
    myfile.open (textFileName, ios::out	|ios::trunc);
    UInt_t det  = 8;
    std::stringstream resultStream;
    std::stringstream header1Stream;
    std::stringstream header2Stream;
    header1Stream <<"#RunNo"<<emptyString(1)<<"\t";
    header2Stream <<"#RunNo"<<emptyString(1)<<"\t";
    resultStream << TString::Format("%7d\t",runnumber);

    header1Stream <<emptyString(TMath::Max((Int_t)runDescription.Length(),6))<<"\t";
    header2Stream <<"descr."<<emptyString(runDescription.Length()-6>0?runDescription.Length()-6:0)<<"\t";
    resultStream << runDescription <<emptyString(6-runDescription.Length()>0?6-runDescription.Length():0)<<"\t";

    header1Stream <<"Noise\t";
    header2Stream <<"Noise"<<"\t";
    resultStream << TString::Format("%5.2f\t",noise[det]);

    header1Stream <<emptyString(6)<<"\t";
    header2Stream <<"CMCnoi\t";
    resultStream << TString::Format("%6.2f\t",diaCMCNoise);

    header1Stream <<emptyString(5)<<"\t";
    header2Stream <<"CMN"<<emptyString(2)<<"\t";
    resultStream << TString::Format("%5.2f\t",CMN);

    header1Stream <<"Corre.\t";
    header2Stream <<"CorSil\t";
    resultStream << TString::Format("%+6.2f\t",getAvergSiliconCorrection().first*100);

    header1Stream <<emptyString(6)<<"\t";
    header2Stream <<"sigSil\t";
    resultStream << TString::Format("%+6.2f\t",getAvergSiliconCorrection().second*100);

    header1Stream <<emptyString(6)<<"\t";
    header2Stream <<"CorDia\t";
    resultStream << TString::Format("%+6.2f\t",getAvergDiamondCorrection()*100);

    header1Stream <<"clust. \t";
    header2Stream <<"m_clus \t";
    resultStream << TString::Format("%+7.2f\t",mean_clustered_normal);

    header1Stream <<"       \t";
    header2Stream <<"mp_clus\t";
    resultStream << TString::Format("%+7.2f\t",mp_clustered_normal);

    header1Stream <<"       \t";
    header2Stream <<"w_clus \t";
    resultStream << TString::Format("%+7.2f\t",width_clustered_normal);

    header1Stream <<"       \t";
    header2Stream <<"gs_clus\t";
    resultStream << TString::Format("%+7.2f\t",gSigma_clustered_normal);

    header1Stream <<"normal \t";
    header2Stream <<"m2/10  \t";
    resultStream << TString::Format("%+7.2f\t",mean2outOf10_normal);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"mp2/10 \t";
    resultStream << TString::Format("%+7.2f\t",mp2outOf10_normal);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"w2/10  \t";
    resultStream << TString::Format("%+7.2f\t",width2outOf10_normal);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"sig2/10\t";
    resultStream << TString::Format("%+7.2f\t",gSigma2outOf10_normal);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"m2/2   \t";
    resultStream << TString::Format("%+7.2f\t",meanNoutOfN_normal[1]);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"m4/4   \t";
    resultStream << TString::Format("%+7.2f\t",meanNoutOfN_normal[3]);

    header1Stream <<"trans  \t";
    header2Stream <<"m2/10  \t";
    resultStream << TString::Format("%+7.2f\t",mean2outOf10_trans);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"mp2/10 \t";
    resultStream << TString::Format("%+7.2f\t",mp2outOf10_trans);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"w2/10  \t";
    resultStream << TString::Format("%+7.2f\t",width2outOf10_trans);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"sig2/10\t";
    resultStream << TString::Format("%+7.2f\t",gSigma2outOf10_trans);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"m2/2   \t";
    resultStream << TString::Format("%+7.2f\t",meanNoutOfN_trans[1]);

    header1Stream <<emptyString(7)<<"\t";
    header2Stream <<"m4/4   \t";
    resultStream << TString::Format("%+7.2f\t",meanNoutOfN_trans[3]);

    header1Stream <<"RES_NORM\t";
    header2Stream <<"Res_DG1n\t";
    resultStream << TString::Format("%8.2f\t",doubleGaus1_normal);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_DG2n\t";
    resultStream << TString::Format("%8.2f\t",doubleGaus2_normal);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGSn\t";
    resultStream << TString::Format("%8.2f\t",singleGausShort_normal);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGNn\t";
    resultStream << TString::Format("%8.2f\t",singleGaus_normal);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGFn\t";
    resultStream << TString::Format("%8.2f\t",singleGausFixed_normal);

    header1Stream <<"RES_TRAN\t";
    header2Stream <<"Res_DG1t\t";
    resultStream << TString::Format("%8.2f\t",doubleGaus1_trans);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_DG2t\t";
    resultStream << TString::Format("%8.2f\t",doubleGaus2_trans);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGSt\t";
    resultStream << TString::Format("%8.2f\t",singleGausShort_trans);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGNt\t";
    resultStream << TString::Format("%8.2f\t",singleGausShort_trans);

    header1Stream <<emptyString(8)<<"\t";
    header2Stream <<"Res_SGFt\t";
    resultStream << TString::Format("%8.2f\t",singleGausFixed_trans);

    header1Stream<<emptyString(5)<<"\t";
    header2Stream<<"REV"<<emptyString(2)<<"\t";
    resultStream<<SVN_REV;

    myfile<<header1Stream.str()<<endl;
    myfile<<header2Stream.str()<<endl;
    myfile<<resultStream.str()<<endl;
    myfile.close();
    //	cout<<"CREATE OUPTUT TEXT FILE DONE "<<endl;
}

std::pair<Float_t,Float_t> TResults::getAvrgSilNoise() {
    Float_t mean=0;
    Float_t mean2=0;
    Int_t i = 0;
    for (Int_t det = 0; det < 8 && det< noise.size()-1; det++){
        mean += noise[det];
        mean2+= noise[det] * noise[det];
        i++;
    }
    mean /= (Float_t)i;
    mean2 /= (Float_t)i;
    Float_t sigma = TMath::Sqrt(mean2-mean*mean);
    return make_pair(mean,sigma);
}


void TResults::setFloatValue(TString section, TString key, Float_t value) {
    FloatMap[section+(TString)"_"+key] = value;
    addKey(section,key,"Float");
}

void TResults::setIntValue(TString section, TString key, Int_t value) {
    IntegerMap[section+(TString)"_"+key] = value;
    addKey(section,key,"Int");
}

void TResults::setStringValue(TString section, TString key, TString value) {
    StringMap[section+(TString)"_"+key] = value;
    addKey(section,key,"String");
}

void TResults::addKey(TString section, TString key, TString type) {
    if (false) cout <<"addKey: "<<section<<" - "<<key<<": "<<type<<endl;
    if(keyList.count(section) ==  0)
        keyList[section] = map<TString,TString>();
    keyList[section][key] = type;
}
