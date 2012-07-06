/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(TSettings* settings) {
  cout << "\n\n\n**********************************************************" << endl;
  cout << "*************TAlignment::TAlignment***********************" << endl;
  cout << "**********************************************************" << endl;

  sys = gSystem;
  setSettings(settings);
  runNumber = settings->getRunNumber();
  cout << runNumber << endl;
  stringstream runString;
  runString.str("");
  runString << settings->getRunNumber();

  sys->MakeDirectory(runString.str().c_str());
  sys->cd(runString.str().c_str());
  htmlAlign = new THTMLAlignment(settings);
  stringstream filepath;
  filepath.str("");
  filepath << sys->pwd();
  filepath << "/selectionData." << runNumber << ".root";
  cout << "OPEN eventReader with file \"" << filepath.str() << endl;
  eventReader = new TADCEventReader(filepath.str(), settings->getRunNumber());
  //eventReader->checkADC();
  histSaver = new HistogrammSaver();
  sys->MakeDirectory("alignment");
  sys->cd("alignment");
  stringstream plotsPath;
  plotsPath << sys->pwd() << "/";
  histSaver->SetPlotsPath(plotsPath.str().c_str());
  histSaver->SetRunNumber(runNumber);
  histSaver->SetNumberOfEvents(eventReader->GetEntries());
  htmlAlign->setFileGeneratingPath(sys->pwd());
  sys->cd("..");
  cout << "end initialise" << endl;
  alignmentPercentage = settings->getAlignment_training_track_fraction();
  Float_t stripSize = 1.;    // 50./10000.;//mu m
  detectorD0Z = 0.725 / stripSize;    // by definition in cm
  detectorD1Z = 1.625 / stripSize;    // by definition in cm
  detectorD2Z = 18.725 / stripSize;    // by definition in cm
  detectorD3Z = 19.625 / stripSize;    // by definition in cm
  detectorDiaZ = 10.2 / stripSize;    // by definition in cm
  verbosity = 1;
  res_keep_factor = settings->getRes_keep_factor();
  cout << "Res Keep factor is set to " << res_keep_factor << endl;
  align = NULL;
  myTrack = NULL;
  nAlignmentStep = -1;
  nAlignSteps = 5;
  nDiaAlignmentStep = -1;
  nDiaAlignSteps = 5;

  diaCalcMode = TCluster::maxValue;    //todo
  silCalcMode = TCluster::corEta;    //todo

  plotAll = true;
  results=0;

}

TAlignment::~TAlignment() {
  cout << "TAlignment deconstructor" << endl;
  htmlAlign->createContent();
  htmlAlign->generateHTMLFile();
  if(htmlAlign!=0)delete htmlAlign;
  if (results!=0)results->setAlignment(this->align);
  if (myTrack) delete myTrack;
  if (histSaver) delete histSaver;
  //	if(eventReader)delete eventReader;
  sys->cd("..");
}

void TAlignment::setSettings(TSettings* settings) {
  this->settings = settings;
}

/**
 * initialise the variable align and set the Z Offsets
 */
void TAlignment::initialiseDetectorAlignment() {
  if (align == NULL) {
    align = new TDetectorAlignment();
    htmlAlign->setAlignment(align);
    cout << "TAlignment::Align::Detectoralignment did not exist, so created new DetectorAlignment" << endl;
    align->SetZOffset(0, detectorD0Z);
    align->SetZOffset(1, detectorD1Z);
    align->SetZOffset(2, detectorD2Z);
    align->SetZOffset(3, detectorD3Z);
    align->SetZOffset(4, detectorDiaZ);
  }
}
/**
 *
 * @param nEvents
 * @param startEvent
 */
void TAlignment::createEventVectors(UInt_t nEvents, UInt_t startEvent) {
  initialiseDetectorAlignment();
  if (nEvents == 0) nEvents = eventReader->GetEntries() - startEvent;
  if (nEvents + startEvent > eventReader->GetEntries()) nEvents = eventReader->GetEntries() - startEvent;
  int noHitDet = 0;
  //	int falseClusterSizeDet=0;
  //int noHitDia=0;
  int falseClusterSizeDia = 0;
  int nCandidates = 0;
  int nScreened = 0;
  cout << "CREATING VECTOR OF VALID EVENTS..." << endl;

  for (nEvent = startEvent; nEvent < nEvents + startEvent; nEvent++) {
    TRawEventSaver::showStatusBar(nEvent - startEvent, nEvents, 100);
    eventReader->LoadEvent(nEvent);
    if (!eventReader->isValidTrack()) {
      noHitDet++;
      continue;
    }
    if (eventReader->isDetMasked()) {
      nScreened++;
      continue;
    }
    if (eventReader->getNDiamondClusters() != 1) {
      falseClusterSizeDia++;
      continue;
    }
    if (eventReader->useForAlignment()) {
      nCandidates++;
      this->events.push_back(*eventReader->getEvent());
    }
    if (eventReader->useForAnalysis()) {
      cout << "\nFound first Event for Analyis ->BREAK" << endl;
      break;
    }
  }
  align->addEventIntervall(startEvent, nEvent);
  align->setNUsedEvents((UInt_t) events.size());
  align->setAlignmentTrainingTrackFraction(settings->getAlignment_training_track_fraction());
  align->setRunNumber(settings->getRunNumber());
}

/**
 * @brief main function of TAlignment, creates eventVector and makes alignemnt
 *
 * @param nEvents	number of Events form which the alignment events are taken
 * @param startEvent first event number
 * @return 1 if everything worked else return 0
 *
 * @todo work on return values
 */
int TAlignment::Align(UInt_t nEvents, UInt_t startEvent) {
  if (verbosity) {
    cout << "\n\n\nTAlignment::Align:Starting \"" << histSaver->GetPlotsPath() << "\"" << endl;
    cout << "\t\t" << events.size() << "\t";
    cout << "\t\t " << eventReader << " ." << endl;
  }

  initialiseDetectorAlignment();
  if (events.size() == 0) createEventVectors(nEvents, startEvent);
  if (myTrack == NULL) {
    cout << "TAlignment::Align::create new TTrack" << endl;
    myTrack = new TTrack(align);
    cout << "TAlignment::Align::created new TTrack" << endl;
    for (UInt_t det = 0; det < TPlaneProperties::getNDetectors(); det++)
      myTrack->setEtaIntergral(det, eventReader->getEtaIntegral(det));
  }
  alignSiliconPlanes();
  AlignDiamondPlane();
  align->PrintResults((UInt_t) 1);

  this->saveAlignment();

  return (1);
}

/**
 *
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 */
void TAlignment::alignSiliconPlanes() {
  nAlignmentStep = -1;
  //getChi2Distribution(10000);
  cout << "Check Detector Alignment" << endl;

//  CheckDetectorAlignment(TPlaneProperties::XY_COR, 1, 0, 0, true);
//  CheckDetectorAlignment(TPlaneProperties::XY_COR, 2, 0, 0, true);
//  CheckDetectorAlignment(TPlaneProperties::XY_COR, 3, 0, 0, true);

  CheckDetectorAlignment(TPlaneProperties::XY_COR, 1, 0, 3, true);
  CheckDetectorAlignment(TPlaneProperties::XY_COR, 2, 0, 3, true);
  CheckDetectorAlignment(TPlaneProperties::XY_COR, 3, 1, 2, true);
//  cout << "ALign Detector with Plane 0" << endl;
//  alignDetector(TPlaneProperties::XY_COR, 1, 0, 0, false);
//  alignDetector(TPlaneProperties::XY_COR, 2, 0, 0, false);
  alignDetector(TPlaneProperties::XY_COR, 3, 0, 0, false);
  alignDetector(TPlaneProperties::XY_COR, 3, 0, 0, false);

  if (verbosity) cout << endl;
  cout << "STart with Alignment Steps" << endl;
  bool bAlignmentGood = true;
  for (nAlignmentStep = 0; nAlignmentStep < nAlignSteps; nAlignmentStep++) {    //||(nAlignmentStep<10&&!bAlignmentGood);nAlignmentStep++){
    bAlignmentGood = true;
    bAlignmentGood = siliconAlignmentStep(false||plotAll);
  }
  getFinalSiliconAlignmentResuluts();
}

bool TAlignment::siliconAlignmentStep(bool bPlot) {
 TPlaneProperties::enumCoordinate coordinates = TPlaneProperties::XY_COR;//((nAlignmentStep < nAlignSteps) ? TPlaneProperties::XY_COR : TPlaneProperties::Y_COR);

  cout << "\n\n\nALIGNMENT STEP:\t" << nAlignmentStep + 1 << " of " << nAlignSteps << "\n\n" << endl;
  bool isAlignmentDone = true;
//  if (nAlignmentStep < 2) {
//    for (UInt_t plane = 0; plane < TPlaneProperties::getNSiliconPlanes(); plane++) {
//      vector<UInt_t> vecRefPlanes;
//      for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
//        if (refPlane != plane) vecRefPlanes.push_back(refPlane);
//      TResidual res = CheckDetectorAlignment(TPlaneProperties::XY_COR, plane, vecRefPlanes, false);
//      res = CheckDetectorAlignment(coordinates, plane, vecRefPlanes, false, res);
//      align->setXResolution(res.getXSigma(), plane);
//      align->setYResolution(res.getYSigma(), plane);
//    }
//    cout << endl;
//  }

  UInt_t subjectPlane1 = 1;
  UInt_t refPlane1_1 = 0, refPlane1_2 = 3;
  cout << "\nAlign Plane " << subjectPlane1 << " with Plane " << refPlane1_1 << " and " << refPlane1_2 << endl;
  TResidual resPlane1 = CheckDetectorAlignment(TPlaneProperties::XY_COR, subjectPlane1, refPlane1_1, refPlane1_2, false);
  if (verbosity) cout << "\nnCheckAlignmentPlane1:" << resPlane1.isTestResidual() << "  " << resPlane1.getXMean() << " +/- " << resPlane1.getXSigma() << endl;
  alignDetector(coordinates, subjectPlane1, refPlane1_1, refPlane1_2, bPlot || plotAll, resPlane1);
  Float_t xOff = align->GetLastXOffset(1);
  Float_t yOff = align->GetLastYOffset(1);
  if (xOff > 0.01 || yOff > 0.01) isAlignmentDone = false;

  UInt_t subjectPlane2 = 2;
  UInt_t refPlane2_1 = 0, refPlane2_2 = 3;
  cout << "\nAlign Plane " << subjectPlane2 << " with Plane " << refPlane2_1 << " and " << refPlane2_2 << endl;
  TResidual resPlane2 = CheckDetectorAlignment(TPlaneProperties::XY_COR, subjectPlane2, refPlane2_1, refPlane2_2, false);
  if (verbosity) cout << "CheckAlignmentPlane2:" << resPlane2.isTestResidual() << ", " << resPlane2.getXMean() << " +/- " << resPlane2.getXSigma() << endl;
  alignDetector(coordinates, subjectPlane2, refPlane2_1, refPlane2_2, bPlot || plotAll, resPlane2);
  xOff = align->GetLastXOffset(2);
  yOff = align->GetLastYOffset(2);
  if (xOff > 0.01 || yOff > 0.01) isAlignmentDone = false;
//
  UInt_t subjectPlane3 = 3;
  UInt_t refPlane3_1 = 0, refPlane3_2 = 0;
  cout << "\nAlign Plane " << subjectPlane3 << " with Plane " << refPlane3_1 << " and " << refPlane3_2 << endl;
  vector<UInt_t> vecRefPlanes;
  vecRefPlanes.push_back(refPlane3_1);
  vecRefPlanes.push_back(refPlane3_2);
  //		vecRefPlanes.push_back(2);
  TResidual resPlane3 = CheckDetectorAlignment(TPlaneProperties::XY_COR, subjectPlane3, vecRefPlanes, false);
  if (verbosity) cout << "CheckAlignment Plane3:" << resPlane3.isTestResidual() << " ," << resPlane3.getXMean() << " +/- " << resPlane3.getXSigma() << endl;
  resPlane3 = alignDetector(coordinates, subjectPlane3, vecRefPlanes, bPlot || plotAll, resPlane3);
  xOff = align->GetLastXOffset(3);
  yOff = align->GetLastYOffset(3);
  if (xOff > 0.01 || yOff > 0.01) isAlignmentDone = false;


  //	align->PrintResults(0);
  if (isAlignmentDone)
    cout << "AlignmentStep was successfully done" << endl;
  else
    cout << "AlignmentStep was NOT the finally alignment step" << endl;
  cout << endl;
  return (isAlignmentDone);
}

void TAlignment::AlignDiamondPlane() {
  verbosity = 1;
  cout << "\n\n\n*******************************************************" << endl;
  cout << "*******************************************************" << endl;
  cout << "***************** Align Diamond ***********************" << endl;
  cout << "*******************************************************" << endl;
  cout << "*******************************************************\n" << endl;

  //create ReferencePlane Vector:
  UInt_t diaPlane = 4;
  vector<UInt_t> vecRefPlanes;
  for (UInt_t i = 0; i < 4; i++)
    if (i != diaPlane) vecRefPlanes.push_back(i);
  nDiaAlignmentStep = -1;
  TResidual resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false, true);
  for (nDiaAlignmentStep = 0; nDiaAlignmentStep < nDiaAlignSteps; nDiaAlignmentStep++) {
    cout << "\n\n " << nDiaAlignmentStep << " of " << nDiaAlignSteps << " Steps..." << endl;
    alignStripDetector(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, false, resDia);
    resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes);
  }
  verbosity = 0;
  resDia = CheckStripDetectorAlignment(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, true, true, resDia);
  CheckStripDetectorAlignmentChi2(TPlaneProperties::X_COR, diaPlane, vecRefPlanes, true, true, settings->getAlignment_chi2());
  align->setDiamondDate();
}


/**
 *
 * @param cor
 * @param subjectPlane
 * @param refPlane1
 * @param refPlane2
 * @param bPlot
 * @param resOld
 * @return
 */
TResidual TAlignment::alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld) {
  vector<UInt_t> vecRefPlanes;
  vecRefPlanes.push_back(refPlane1);
  if (refPlane1 != refPlane2) vecRefPlanes.push_back(refPlane2);
  return (alignDetector(cor, subjectPlane, vecRefPlanes, bPlot, resOld));
}

TResidual TAlignment::alignDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {
  if (true) {    //some outputs
    if (verbosity) cout << "TAlignment::AlignDetector::\taligning Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << vecRefPlanes.size() << " Planes: ";
    for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
      if (verbosity) cout << vecRefPlanes.at(i) << " ";
    if (verbosity) cout << "\t";
    if (bPlot) if (verbosity) cout << " plots are created\t";
    if (resOld.isTestResidual()) if (verbosity)
      cout << "resOld is a testResidual" << endl;
    else {
      if (verbosity) cout << endl;
      if (verbosity) resOld.Print(2);
    }
    for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
      if (vecRefPlanes.at(i) == subjectPlane) cerr << "Plane " << subjectPlane << " is used as a reference Plane and as a subject Plane" << endl;
  }

  //get Residual
  TResidual res = this->getResidual(cor, subjectPlane, vecRefPlanes, bPlot, resOld);
  Float_t x_offset = res.getXOffset();
  Float_t phix_offset = res.getPhiXOffset();
  Float_t y_offset = res.getYOffset();
  Float_t phiy_offset = res.getPhiYOffset();
  if (subjectPlane == 3) {
    cout << "Xoff: " << x_offset << "\tPhiXoff: " << phix_offset << "\tYoff: " << y_offset << "\tphiY:" << phiy_offset << endl;
  }
  printf("Correction Values: X: %2.6f,  PhiX: %2.6f,   Y: %2.6f,  PhiY: %2.6f\n", x_offset, phix_offset, y_offset, phiy_offset);

  //save corrections to alignment
  if (vecRefPlanes.size() == 1) {
    if (TPlaneProperties::X_COR == cor || TPlaneProperties::XY_COR == cor) align->AddToXOffset(subjectPlane, x_offset);
    if (TPlaneProperties::Y_COR == cor || TPlaneProperties::XY_COR == cor) align->AddToYOffset(subjectPlane, y_offset);
  } else {
    if (subjectPlane == 4 && nDiaAlignmentStep == 0) {
      align->AddToXOffset(subjectPlane, x_offset);
      return (res);
    }
    if (TPlaneProperties::X_COR == cor || TPlaneProperties::XY_COR == cor) {
      if (x_offset != N_INVALID) align->AddToXOffset(subjectPlane, x_offset);
      if (phiy_offset != N_INVALID) align->AddToPhiYOffset(subjectPlane, phiy_offset);
    }
    if (TPlaneProperties::Y_COR == cor || TPlaneProperties::XY_COR == cor) {
      if (y_offset != N_INVALID) align->AddToYOffset(subjectPlane, y_offset);
      if (phix_offset != N_INVALID) align->AddToPhiXOffset(subjectPlane, phix_offset);
    }
  }
  this->myTrack->setDetectorAlignment(align);
  return (res);
}

TResidual TAlignment::alignStripDetector(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {
  if (verbosity) {
    cout << "TAlignment::AlignStripDetector::\taligning Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << vecRefPlanes.size() << " Planes: ";

    for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
      cout << vecRefPlanes.at(i) << " ";
    cout << "\t";
    if (bPlot) cout << " plots are created\t";
    if (resOld.isTestResidual())
      cout << "resOld is a testResidual" << endl;
    else {
      cout << endl;
      resOld.Print(2);
    }
    for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
      if (vecRefPlanes.at(i) == subjectPlane) cerr << "Plane " << subjectPlane << " is used as a reference Plane and as a subject Plane" << endl;
  }
  //get residual
  TResidual res = this->getStripResidual(cor, subjectPlane, vecRefPlanes, false, bPlot, resOld);
  Float_t x_offset = res.getXOffset();
  Float_t phix_offset = res.getPhiYOffset();
  Float_t y_offset = res.getYOffset();
  Float_t phiy_offset = res.getPhiXOffset();

  if (verbosity||true) cout << "Correction Values: X:" << x_offset << ",  PhiX: " << phix_offset << ",   Y: " << y_offset << ",  PhiY: " << phiy_offset << "\n" << endl;

  if (subjectPlane == 4 && nDiaAlignmentStep == 0) {
    align->AddToXOffset(subjectPlane, x_offset);
    return (res);
  }
  if (TPlaneProperties::X_COR == cor || TPlaneProperties::XY_COR == cor) {
    align->AddToXOffset(subjectPlane, x_offset);
    align->AddToPhiXOffset(subjectPlane, phix_offset);
  }
  if (TPlaneProperties::Y_COR == cor || TPlaneProperties::XY_COR == cor) {
    align->AddToYOffset(subjectPlane, y_offset);
    align->AddToPhiYOffset(subjectPlane, phiy_offset);
  }

  return res;
}

TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld, TCluster::calculationMode_t mode) {
  vector<UInt_t> vecRefPlanes;
  vecRefPlanes.push_back(refPlane1);
  if (refPlane1 != refPlane2) vecRefPlanes.push_back(refPlane2);
  return getResidual(cor, subjectPlane, vecRefPlanes, bPlot, resOld, mode);
}

TResidual TAlignment::getResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld, TCluster::calculationMode_t mode) {
  for(UInt_t refPlane1=0;refPlane1<vecRefPlanes.size()-1;refPlane1++){
    for(UInt_t refPlane2=refPlane1+1;refPlane2<vecRefPlanes.size();refPlane2++){
        if(vecRefPlanes.at(refPlane1)==vecRefPlanes.at(refPlane2))
          vecRefPlanes.erase(vecRefPlanes.begin()+refPlane1);
    }
  }
  stringstream refPlaneString;
  for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
    if (i == 0)
      refPlaneString << vecRefPlanes.at(i);
    else if (i + 1 < vecRefPlanes.size())
      refPlaneString << "_" << vecRefPlanes.at(i);
    else
      refPlaneString << "_and_" << vecRefPlanes.at(i);
  if (verbosity) cout << "TAlignment::getResidual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << refPlaneString.str() << ", plotting: " << bPlot << "  with " << alignmentPercentage << "\t" << resOld.isTestResidual() << endl;
  clearMeasuredVectors();
  Float_t xDelta, yDelta;
  Float_t xPositionObserved, yPositionObserved;
  Float_t xMeasured, yMeasured;
  Float_t xPredicted, yPredicted;
  Float_t resxtest, resytest;
  Float_t xPhi,yPhi;
  TPositionPrediction* predictedPostion = 0;
  for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));
    xPositionObserved = myTrack->getPosition(TPlaneProperties::X_COR, subjectPlane, mode);
    yPositionObserved = myTrack->getPosition(TPlaneProperties::Y_COR, subjectPlane, mode);
    predictedPostion = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
    xPredicted = predictedPostion->getPositionX();
    yPredicted = predictedPostion->getPositionY();
    xPhi = predictedPostion->getPhiX()*360./TMath::TwoPi();
    yPhi = predictedPostion->getPhiY()*360./TMath::TwoPi();
    xDelta = xPositionObserved - xPredicted;    //X_OBS-X_Pred
    yDelta = yPositionObserved - yPredicted;    //Y_OBS-Y_Pred
    resxtest = TMath::Abs(xDelta - resOld.getXMean()) / resOld.getXSigma();
    resytest = TMath::Abs(yDelta - resOld.getYMean()) / resOld.getYSigma();
    if (verbosity > 3) cout << nEvent << endl;
    //if(verbosity>3)	predictedPostion->Print();
    if (verbosity > 3) cout << "Measured: " << myTrack->getXMeasured(subjectPlane, mode) << "/" << myTrack->getYMeasured(subjectPlane, mode) << endl;
    if (verbosity > 3) cout << "Observed: " << xPositionObserved << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "Predicted: " << predictedPostion->getPositionX() << " / " << predictedPostion->getPositionY() << endl;
    if (verbosity > 3) cout << "Delta:    " << xDelta << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "ResTest:  " << resxtest << " / " << resytest << "\n\n" << endl;

    if (resxtest < res_keep_factor && resytest < res_keep_factor) {
      vecXObs.push_back(xPositionObserved);
      vecYObs.push_back(yPositionObserved);
      vecXDelta.push_back(xDelta);
      vecYDelta.push_back(yDelta);
      vecXPred.push_back(xPredicted);
      vecYPred.push_back(yPredicted);
      vecXMeasured.push_back(xMeasured);
      vecYMeasured.push_back(yMeasured);
      vecXPhi.push_back(xPhi);
      vecYPhi.push_back(yPhi);
      vecXChi2.push_back(predictedPostion->getChi2X());
      vecYChi2.push_back(predictedPostion->getChi2Y());
    }
    if (verbosity > 3) cout << xDelta << " " << yDelta << endl;
    predictedPostion->Delete();
  }

  if (verbosity > 2) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << endl;
  //first estimate residuals widths
  TResidual res;
  res.setResKeepFactor(res_keep_factor);
  res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta, resOld);
  this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot);

  clearMeasuredVectors();
  return res;

}
/**
 * @brief creates element TResidual to adjust the alignment
 *
 * creates a vector of pedicted X positions, predicted Y positions, delta X and delta Y
 * and use the function calculateResidual to get the residual with this vectors
 * @param	cor coordindate for which the residual is calculated
 * @param	subjectPlane plane for which the residual should be calculated
 * @param	refPlane1 first reference plane
 * @param	refPlane2 second reference plane
 * @param	bPlot	variable to create plots or not
 *
 * @return
 */
TResidual TAlignment::getStripResidual(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld, TCluster::calculationMode_t mode) {

  stringstream refPlaneString;
  for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
    if (i == 0)
      refPlaneString << vecRefPlanes.at(i);
    else if (i + 1 < vecRefPlanes.size())
      refPlaneString << "_" << vecRefPlanes.at(i);
    else
      refPlaneString << "_and_" << vecRefPlanes.at(i);
  if (verbosity) cout << "TAlignment::getStripResidual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << refPlaneString.str() << ", plotting: " << bPlot << "  with " << alignmentPercentage << "\t" << resOld.isTestResidual() << endl;
  clearMeasuredVectors();

  Float_t deltaX, deltaY;
  Float_t xPositionObserved;
  Float_t yPositionObserved;
  TPositionPrediction* predictedPostion = 0;
  Float_t resxtest, resytest;
  Float_t xPhi,yPhi;
  for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));
    predictedPostion = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, verbosity > 1);
    xPositionObserved = myTrack->getStripXPosition(subjectPlane, predictedPostion->getPositionY(), TCluster::maxValue);
    yPositionObserved = myTrack->getPosition(TPlaneProperties::Y_COR, subjectPlane, TCluster::maxValue);
    deltaX = xPositionObserved - predictedPostion->getPositionX();    //X_OBS-X_Pred
    deltaY = yPositionObserved - predictedPostion->getPositionY();    //Y_OBS-Y_Pred
    resxtest = TMath::Abs(deltaX - resOld.getXMean()) / resOld.getXSigma();
    resytest = TMath::Abs(deltaY - resOld.getYMean()) / resOld.getYSigma();
    xPhi = predictedPostion->getPhiX();
    yPhi = predictedPostion->getPhiY();
    if (verbosity > 3) cout << "Event no.: " << nEvent << endl;
    //if(verbosity>3)	predictedPostion->Print();
    if (verbosity > 3) events.at(nEvent).getPlane(subjectPlane).Print();
    if (verbosity > 3) cout << "Measured: " << myTrack->getXMeasured(subjectPlane) << "/" << myTrack->getYMeasured(subjectPlane) << endl;
    if (verbosity > 3) cout << "Observed: " << xPositionObserved << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "Predicted: " << predictedPostion->getPositionX() << " / " << predictedPostion->getPositionY() << endl;
    if (verbosity > 3) cout << "Delta:    " << deltaX << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "ResTest:  " << resxtest << " / " << resytest << "\n\n" << endl;

    if (resxtest < res_keep_factor && resytest < res_keep_factor) {
      vecXObs.push_back(xPositionObserved);
      vecYObs.push_back(yPositionObserved);
      vecXDelta.push_back(deltaX);
      vecYDelta.push_back(deltaY);
      vecXPred.push_back(predictedPostion->getPositionX());
      vecYPred.push_back(predictedPostion->getPositionY());
      vecXPhi.push_back(xPhi);
      vecYPhi.push_back(yPhi);

      vecXChi2.push_back(predictedPostion->getChi2X());
      vecYChi2.push_back(predictedPostion->getChi2Y());
    }
    if (verbosity > 3) cout << deltaX << " " << deltaY << endl;
    predictedPostion->Delete();
  }

  if (verbosity > 2) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << endl;

  //first estimate residuals widths
  TResidual res;
  res.setResKeepFactor(res_keep_factor);
  res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta, resOld);
  this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot, bAlign);
  clearMeasuredVectors();
  return res;

}

TResidual TAlignment::getStripResidualChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2, TCluster::calculationMode_t mode) {

  TResidual resOld;
  stringstream refPlaneString;
  for (UInt_t i = 0; i < vecRefPlanes.size(); i++)
    if (i == 0)
      refPlaneString << vecRefPlanes.at(i);
    else if (i + 1 < vecRefPlanes.size())
      refPlaneString << "_" << vecRefPlanes.at(i);
    else
      refPlaneString << "_and_" << vecRefPlanes.at(i);
  if (verbosity) cout << "TAlignment::getStripResidual of Plane " << subjectPlane << TPlaneProperties::getCoordinateString(cor) << " with " << refPlaneString.str() << ", plotting: " << bPlot << "  with " << alignmentPercentage << "\t" << resOld.isTestResidual() << endl;
  clearMeasuredVectors();

  Float_t deltaX, deltaY;
  Float_t xPositionObserved;
  Float_t yPositionObserved;
  Float_t xMeasured, yMeasured;
  TPositionPrediction* predictedPosition = 0;
  Float_t resxtest, resytest;
  Float_t chi2x, chi2y;
  Float_t xPhi,yPhi;
  for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));
    predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, nEvent < 1);
    chi2x = predictedPosition->getChi2X();
    chi2y = predictedPosition->getChi2Y();
    xPositionObserved = myTrack->getStripXPosition(subjectPlane, predictedPosition->getPositionY(), TCluster::maxValue);
    yPositionObserved = myTrack->getPosition(TPlaneProperties::Y_COR, subjectPlane, TCluster::maxValue);
    deltaX = xPositionObserved - predictedPosition->getPositionX();    //X_OBS-X_Pred
    deltaY = yPositionObserved - predictedPosition->getPositionY();    //Y_OBS-Y_Pred
    xMeasured = myTrack->getMeasured(TPlaneProperties::X_COR, subjectPlane, TCluster::maxValue);
    yMeasured = myTrack->getMeasured(TPlaneProperties::Y_COR, subjectPlane, TCluster::maxValue);
    resxtest = TMath::Abs(deltaX - resOld.getXMean()) / resOld.getXSigma();
    resytest = TMath::Abs(deltaY - resOld.getYMean()) / resOld.getYSigma();
    xPhi=predictedPosition->getPhiX();
    yPhi=predictedPosition->getPhiY();
    if (verbosity > 3) cout << "Event no.: " << nEvent << endl;
    //if(verbosity>3)	predictedPostion->Print();
    if (verbosity > 3) events.at(nEvent).getPlane(subjectPlane).Print();
    if (verbosity > 3) cout << "Measured: " << myTrack->getXMeasured(subjectPlane) << "/" << myTrack->getYMeasured(subjectPlane) << endl;
    if (verbosity > 3) cout << "Observed: " << xPositionObserved << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "Predicted: " << predictedPosition->getPositionX() << " / " << predictedPosition->getPositionY() << endl;
    if (verbosity > 3) cout << "Delta:    " << deltaX << " / " << yPositionObserved << endl;
    if (verbosity > 3) cout << "Chi2:     " << chi2x << " / " << chi2y << endl;
    if (verbosity > 3) cout << "ResTest:  " << resxtest << " / " << resytest << endl;

    if (chi2x < settings->getAlignment_chi2() && chi2y < settings->getAlignment_chi2()) {
      if (verbosity > 3) cout << "Take event for Residual" << endl;
      vecXObs.push_back(xPositionObserved);
      vecYObs.push_back(yPositionObserved);
      vecXDelta.push_back(deltaX);
      vecYDelta.push_back(deltaY);
      vecXPred.push_back(predictedPosition->getPositionX());
      vecYPred.push_back(predictedPosition->getPositionY());
      vecXMeasured.push_back(xMeasured);
      vecYMeasured.push_back(yMeasured);
      vecXPhi.push_back(xPhi);
      vecYPhi.push_back(yPhi);
      vecXChi2.push_back(predictedPosition->getChi2X());
      vecYChi2.push_back(predictedPosition->getChi2Y());
    } else if (verbosity > 3) cout << "through event away..." << endl;
    if (verbosity > 3) cout << "\n\n" << endl;
    //			if(verbosity>3)	cout<< deltaX<<" "<<deltaY<<endl;
    predictedPosition->Delete();
  }
  if (bAlign) {
    align->setNDiamondAlignmentEvents((UInt_t) vecXDelta.size());
    align->setDiaChi2(settings->getAlignment_chi2());
  }

  if (verbosity) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << " " << vecXObs.size() << " " << vecYObs.size() << endl;
  if (verbosity) cout << "used " << vecXDelta.size() << " Events of " << events.size() << " with a Chi2 < " << settings->getAlignment_chi2() << endl;
  //first estimate residuals widths
  TResidual res;
  res.setResKeepFactor(res_keep_factor);
  res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta);
  this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot, bAlign);

  clearMeasuredVectors();
  return res;
}
/**
 * calculateResidual if there is no residuals calculatet to cut on the res_keep_factor;
 * this funcition opens the other calculateResidual function but uses a TResidual res
 * which is using all items...
 */
TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, UInt_t refPlane1, UInt_t refPlane2, bool bPlot, TResidual resOld) {
  vector<UInt_t> vecRefPlanes;
  vecRefPlanes.push_back(refPlane1);
  vecRefPlanes.push_back(refPlane2);
  return CheckDetectorAlignment(cor, subjectPlane, vecRefPlanes, bPlot, resOld);
}

TResidual TAlignment::CheckDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPlot, TResidual resOld) {
  if (verbosity) cout << "\n\nTAlignment::checkDetectorAlignment\n\t check " << TPlaneProperties::getCoordinateString(cor) << " coordinate of Plane " << subjectPlane << " with " << vecRefPlanes.size() << " Planes" << endl;
  //	if(refPlane1==subjectPlane || refPlane2==subjectPlane){
  //		return TResidual(true);
  //	}
  int verb = verbosity;
  verbosity = 0;
  TResidual res = getResidual(cor, subjectPlane, vecRefPlanes, false, resOld);
  if (verbosity) cout << endl;
  res.SetTestResidual(false);
  res = getResidual(cor, subjectPlane, vecRefPlanes, bPlot, res);
  verbosity = verb;
  if (verbosity) res.Print();
  return res;
}

TResidual TAlignment::CheckStripDetectorAlignment(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, TResidual resOld) {
  int verb = verbosity;
  verbosity = 0;
  TResidual res = getStripResidual(cor, subjectPlane, vecRefPlanes, false, false, resOld);
  if (verbosity) cout << endl;
  res.SetTestResidual(false);
  res = getStripResidual(cor, subjectPlane, vecRefPlanes, bAlign, bPlot, res);
  verbosity = verb;
  if (verbosity) res.Print();
  return res;
}

TResidual TAlignment::CheckStripDetectorAlignmentChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bAlign, bool bPlot, Float_t maxChi2) {
  int verb = verbosity;
  verbosity = 1;

  //	getStripResidual(cor,subjectPlane,vecRefPlanes,false,false,maxChi2);
  if (verbosity) cout << endl;
  TResidual res = getStripResidualChi2(cor, subjectPlane, vecRefPlanes, bAlign, bPlot, maxChi2);
  if (verbosity) res.Print();
  verbosity = verb;
  return res;
}
/**
 *
 */
void TAlignment::saveAlignment() {
  stringstream fileName;
  fileName << "alignment." << settings->getRunNumber() << ".root";
  TFile *alignmentFile = new TFile(fileName.str().c_str(), "RECREATE");
  cout << "TAlignment:saveAlignment(): path: \"" << sys->pwd() << "\", file Name:\"" << fileName.str() << "\"" << endl;
  alignmentFile->cd();
  align->SetName("alignment");
  stringstream title;
  title << "alignment of Run " << settings->getRunNumber();
  align->SetTitle(title.str().c_str());
  align->Write();
  alignmentFile->Write();
  alignmentFile->Close();
}

/**
 *
 */
void TAlignment::getChi2Distribution(Float_t maxChi2) {
  vecXChi2.clear();
  vecYChi2.clear();
  //	if(verbosity)
  cout << "TAlignment::getChi2Distribution" << endl;
  //	Float_t xPositionObserved,yPositionObserved,deltaX,deltaY,resxtest,resytest;
  TPositionPrediction* predictedPosition = 0;
  vector<UInt_t> vecRefPlanes;

  for (UInt_t i = 0; i < 4; i++)
    vecRefPlanes.push_back(i);
  //	UInt_t oldVerbosity=myTrack->getVerbosity();
  //	myTrack->setVerbosity(4);
  vector<Float_t> vecSumDeltaX;
  vector<Float_t> vecSumDeltaY;
  for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));
    Float_t sumDeltaX = 0;
    Float_t sumDeltaY = 0;
    Float_t chi2X = 0;
    Float_t chi2Y = 0;
    UInt_t subjectPlane = 0;
    predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
    chi2X = predictedPosition->getChi2X();
    chi2Y = predictedPosition->getChi2Y();
    if (predictedPosition->getChi2X() < maxChi2 && predictedPosition->getChi2Y() < maxChi2) {
      for (subjectPlane = 0; subjectPlane < 4; subjectPlane++) {
        if (subjectPlane != 0) {
          predictedPosition->Delete();
          predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
        }
        Float_t deltaX = myTrack->getXPosition(subjectPlane);
        ;
        Float_t deltaY = myTrack->getYPosition(subjectPlane);
        ;

        deltaX -= predictedPosition->getPositionX();
        deltaY -= predictedPosition->getPositionY();
        sumDeltaX += TMath::Abs(deltaX);
        sumDeltaY += TMath::Abs(deltaY);
      }    //for loop over subjectPlane
      vecXChi2.push_back(chi2X);
      vecYChi2.push_back(chi2Y);
      vecSumDeltaX.push_back(sumDeltaX);
      vecSumDeltaY.push_back(sumDeltaY);
    }    //end if chi2x<maxCi && chi2y<maxChi
    predictedPosition->Delete();
  }    //end for loop over nEvent

  //	myTrack->setVerbosity(oldVerbosity);
  stringstream histName;
  //Chi2X Distribution
  histName.str("");
  if (nAlignmentStep == -1)
    histName << "hPreAlignment";
  else if (nAlignmentStep >= nAlignSteps - 1)
    histName << "hPostAlignment";
  else
    histName << "h_" << nAlignmentStep << "_Step";
  histName << "_Chi2X_Distribution";
  TH1F histoChi2X = histSaver->CreateDistributionHisto(histName.str(), vecXChi2, 4096, HistogrammSaver::positiveSigma);
  histoChi2X.GetXaxis()->SetTitle("Chi^2/NDF of X plane");
  histoChi2X.GetYaxis()->SetTitle("number of entries");

  histSaver->SaveHistogram(&histoChi2X);

  //Chi2Y Distribution
  histName.str("");
  if (nAlignmentStep == -1)
    histName << "hPreAlignment";
  else if (nAlignmentStep >= nAlignSteps - 1)
    histName << "hPostAlignment";
  else
    histName << "h_" << nAlignmentStep << "_Step";
  histName << "_Chi2Y_Distribution";
  TH1F histoChi2Y = histSaver->CreateDistributionHisto(histName.str(), vecYChi2, 4096, HistogrammSaver::positiveSigma);
  histoChi2Y.GetXaxis()->SetTitle("Chi^2/NDF of Y plane");
  histoChi2Y.GetYaxis()->SetTitle("number of entries");
  histSaver->SaveHistogram(&histoChi2Y);

  histName.str("");
  if (nAlignmentStep == -1)
    histName << "hPreAlignment";
  else if (nAlignmentStep >= nAlignSteps - 1)
    histName << "hPostAlignment";
  else
    histName << "h_" << nAlignmentStep << "_Step";
  histName << "_Chi2X_vs_SumDeltaX";
  cout << "CREATE: " << histName.str() << endl;
  TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecSumDeltaX, vecXChi2);
  histo->GetYaxis()->SetTitle("Sum of Delta X");
  histo->GetYaxis()->SetTitle("Chi2 in X");
  histSaver->SaveHistogram(histo);

  histName.str("");
  if (nAlignmentStep == -1)
    histName << "hPreAlignment";
  else if (nAlignmentStep >= nAlignSteps - 1)
    histName << "hPostAlignment";
  else
    histName << "h_" << nAlignmentStep << "_Step";
  histName << "_Chi2Y_vs_SumDeltaY";
  cout << "CREATE PLOT: \"" << histName.str() << "\"" << endl;
  TH2F *histo1 = histSaver->CreateScatterHisto(histName.str(), vecSumDeltaY, vecYChi2);
  histName << "_graph";
  TGraph graph1 = histSaver->CreateDipendencyGraph(histName.str(), vecSumDeltaY, vecYChi2);
  graph1.Draw("AP");
  graph1.GetYaxis()->SetTitle("Sum of Delta Y");
  graph1.GetXaxis()->SetTitle("Chi2 in Y");

  histo1->GetYaxis()->SetTitle("Sum of Delta Y");
  histo1->GetXaxis()->SetTitle("Chi2 in Y");
  histSaver->SaveHistogram(histo1);
  histSaver->SaveGraph(&graph1, histName.str(), "AP");
}

TResidual TAlignment::calculateResidualWithChi2(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, Float_t maxChi2, bool bAlign, bool bPlot) {
  if (verbosity) cout << "\n\nTAlignment::calculateResidualWithChi2" << endl;
  clearMeasuredVectors();
  stringstream refPlaneString;
  for (UInt_t i = 0; i < vecRefPlanes.size(); i++) {
    if (i == 0)
      refPlaneString << vecRefPlanes.at(i);
    else if (i + 1 < vecRefPlanes.size())
      refPlaneString << "_" << vecRefPlanes.at(i);
    else
      refPlaneString << "_and_" << vecRefPlanes.at(i);
  }
  refPlaneString << "with_Chi2_cut_on_" << maxChi2;

  //	Float_t xPositionObserved,yPositionObserved,deltaX,deltaY,resxtest,resytest;
  TPositionPrediction* predictedPosition = 0;

  //	UInt_t oldVerbosity=myTrack->getVerbosity();
  //	myTrack->setVerbosity(4);
  Float_t chi2x, chi2y, xPositionObserved, yPositionObserved, deltaX, deltaY;
  ;    //_
  for (UInt_t nEvent = 0; nEvent < events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));
    predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta, false);
    chi2x = predictedPosition->getChi2X();
    chi2y = predictedPosition->getChi2Y();
    xPositionObserved = myTrack->getPosition(TPlaneProperties::X_COR, subjectPlane);
    yPositionObserved = myTrack->getPosition(TPlaneProperties::Y_COR, subjectPlane);
    predictedPosition = myTrack->predictPosition(subjectPlane, vecRefPlanes, TCluster::corEta);
    if (verbosity > 3) predictedPosition->Print();
    if (verbosity > 3) cout << xPositionObserved << " / " << yPositionObserved << endl;
    deltaX = xPositionObserved - predictedPosition->getPositionX();    //X_OBS-X_Pred
    deltaY = yPositionObserved - predictedPosition->getPositionY();    //Y_OBS-Y_Pred

    if (chi2x < maxChi2 && chi2y < maxChi2) {
      vecXObs.push_back(xPositionObserved);
      vecYObs.push_back(yPositionObserved);
      vecXDelta.push_back(deltaX);
      vecYDelta.push_back(deltaY);
      vecXPred.push_back(predictedPosition->getPositionX());
      vecYPred.push_back(predictedPosition->getPositionY());
      vecXChi2.push_back(predictedPosition->getChi2X());
      vecYChi2.push_back(predictedPosition->getChi2Y());
    }

    predictedPosition->Delete();
  }
  cout << "use " << vecXDelta.size() << " of " << events.size() << " Events ( " << (Float_t) vecXDelta.size() / (Float_t) events.size() * 100. << " %)" << endl;

  TResidual res;
  res.setResKeepFactor(res_keep_factor);
  res.calculateResidual(cor, &vecXPred, &vecXDelta, &vecYPred, &vecYDelta);
  res.SetTestResidual(false);
  this->CreatePlots(cor, subjectPlane, refPlaneString.str(), bPlot, bAlign);

  clearMeasuredVectors();

  if (verbosity) res.Print(1);
  return res;
}

void TAlignment::getFinalSiliconAlignmentResuluts() {
  Float_t maxChi2 = settings->getAlignment_chi2();
  cout << "get Final Silicon Resolution with a  maximum chi2 of " << maxChi2 << endl;
  //set Resolutions
  //	setDetectorResolution(maxChi2);
  setSiliconDetectorResolution(maxChi2);

  //	vector<UInt_t>vecRefPlanes;
  //	for(UInt_t plane=0;plane<4;plane++){
  //		vecRefPlanes.clear();
  //		for(UInt_t i=0;i<4;i++)
  //			if(i!=plane)vecRefPlanes.push_back(i);
  //			TResidual res0=calculateResidualWithChi2(TPlaneProperties::XY_COR,plane,vecRefPlanes,maxChi2,false,true);
  //	}

  getChi2Distribution(15);
}

void TAlignment::setSiliconDetectorResolution(Float_t maxChi2) {
  //get  something like a aprroximate Sigma with calculating the residual
  for (UInt_t plane = 0; plane < 4; plane++) {
    vector<UInt_t> vecRefPlanes;
    for (UInt_t i = 0; i < 4; i++)
      if (i != plane) vecRefPlanes.push_back(i);
    TResidual res = CheckDetectorAlignment(TPlaneProperties::XY_COR, plane, vecRefPlanes, false);
    res = CheckDetectorAlignment(TPlaneProperties::XY_COR, plane, vecRefPlanes, false, res);
    align->setXResolution(res.getXSigma(), plane);
    align->setYResolution(res.getYSigma(), plane);
  }

  for (UInt_t plane = 0; plane < 4; plane++) {
    vector<UInt_t> vecRefPlanes;
    for (UInt_t i = 0; i < 4; i++)
      if (i != plane) vecRefPlanes.push_back(i);
    TResidual res = calculateResidualWithChi2(TPlaneProperties::XY_COR, plane, vecRefPlanes, maxChi2, true, true);
  }
}

void TAlignment::CreatePlots(TPlaneProperties::enumCoordinate cor, UInt_t subjectPlane, string refPlaneString, bool bPlot, bool bUpdateAlignment, bool bChi2) {
  if (bPlot) if (verbosity) cout << vecXDelta.size() << " " << vecYDelta.size() << " " << vecXPred.size() << " " << vecYPred.size() << " " << vecXObs.size() << " " << vecYObs.size() << endl;
  if (!bPlot && !bUpdateAlignment) return;
  stringstream preName;
  if (subjectPlane == 4) {
    preName << "hDiamond_";
    if (nDiaAlignmentStep == -1)
      preName << "PreAlignment";
    else if (nDiaAlignmentStep == nDiaAlignSteps)
      preName << "PostAlignment";
    else
      preName << nDiaAlignmentStep << "Step";
  } else {
    preName << "hSilicon_";
    if (nAlignmentStep == -1)
      preName << "PreAlignment";
    else if (nAlignmentStep == nAlignSteps)
      preName << "PostAlignment";
    else
      preName << nAlignmentStep << "_Step";
  }
  stringstream histName;
  cout << "\nCreatePlots with " << preName.str() << " " << nAlignmentStep << flush;
  if (bUpdateAlignment)
    cout << "\twith Alignment Resolution Update\n" << endl;
  else
    cout << endl;

  if (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) {    //DistributionPlot DeltaX
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_DistributionPlot_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
     //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH1F* histo = (TH1F*) histSaver->CreateDistributionHisto(histName.str(), vecXDelta, 512, HistogrammSaver::threeSigma).Clone();
    int ntries=0;
    while ((histo->GetBinContent(histo->GetMaximumBin())/histo->GetEntries())<0.05&&ntries<3){//todo change hardcoding
      histo->Rebin();ntries++;
    }
    TF1* fitGausX = new TF1("fitGaus", "gaus", -1, 1);
    if (bUpdateAlignment) {
      cout << "Alignment for plane" << subjectPlane << endl;
      histo->Draw("goff");
      histo->Fit(fitGausX, "Q", "", 0.5, 0.5);
      Float_t xRes = fitGausX->GetParameter(2);
      cout << "set Resolution via Gaus-Fit: " << xRes << " with " << vecXDelta.size() << " Events" << endl;
      align->setXResolution(xRes, subjectPlane);
      histo->GetXaxis()->SetRangeUser(-5 * xRes, +5 * xRes);
    }
    histo->GetXaxis()->SetTitle("Delta X in Channels");
    histo->GetYaxis()->SetTitle("Number of entries");
    if (bPlot) histSaver->SaveHistogram(histo);
    delete fitGausX;
    delete histo;
  }

  if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    //DistributionPlot DeltaX vs Ypred
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_YPred_vs_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateDipendencyHisto(histName.str(), vecYPred, vecXDelta, 256);
//    histo.Draw("goff");
    histo->GetXaxis()->SetTitle("Y Predicted");
    histo->GetYaxis()->SetTitle("Delta X");

    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram(histo);
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecYPred);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("predicted Y position in ChannelNo.");
    graph.GetYaxis()->SetTitle("delta X in Channel No.");
    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

  if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) && subjectPlane == TPlaneProperties::getDiamondPlane()) {    //DistributionPlot DeltaX vs Xpred
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_XPred_vs_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateDipendencyHisto(histName.str(), vecXPred, vecXDelta, 256);
    histo->GetXaxis()->SetTitle("X Predicted");
    histo->GetYaxis()->SetTitle("Delta X");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram(histo);
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXPred);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("predicted X position in ChannelNo.");
    graph.GetYaxis()->SetTitle("delta X in Channel No.");
    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

  if (bPlot && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR) && subjectPlane == TPlaneProperties::getDiamondPlane()) {    //DistributionPlot DeltaX vs Xpred
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_XMeasured_vs_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F* histo = histSaver->CreateDipendencyHisto(histName.str(), vecXMeasured, vecXDelta, 256);
//    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("X Measured");
    histo->GetYaxis()->SetTitle("Delta X");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram(histo);
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXMeasured);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("measured X position in ChannelNo.");
    graph.GetYaxis()->SetTitle("delta X in Channel No.");
    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

  if (subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //DistributionPlot DeltaY
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_DistributionPlot_DeltaY";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH1F* histo = (TH1F*) histSaver->CreateDistributionHisto(histName.str(), vecYDelta, 512, HistogrammSaver::threeSigma).Clone();

    histo->Draw("goff");
    Float_t sigma = histo->GetRMS();
    Float_t fitWidth = sigma *1.5;
    Float_t mean = histo->GetMean();
    TF1* fitGausY = new TF1("fitGaus", "gaus", mean-fitWidth*2, mean+fitWidth*2);
    histo->Fit(fitGausY, "Q", "", mean-fitWidth, mean+fitWidth);
    Float_t yRes = fitGausY->GetParameter(2);
    Float_t skewness = histo->GetSkewness();
    cout<<"\nSKEWNESS of "<<histName.str()<<": "<<skewness<<endl;
    if (bUpdateAlignment) {
      align->setYResolution(yRes, subjectPlane);
      histo->GetXaxis()->SetRangeUser(-5 * yRes, +5 * yRes);
    }
    histo->GetXaxis()->SetTitle("Delta Y in Channels");
    histo->GetYaxis()->SetTitle("Number of entries");
    if (bPlot) histSaver->SaveHistogram(histo);
    delete fitGausY;
    delete histo;
  }

  if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterPlot DeltaY vs Xpred
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_XPred_vs_DeltaY";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateDipendencyHisto(histName.str(), vecXPred, vecYDelta, 256);
//    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("X Predicted");
    histo->GetYaxis()->SetTitle("Delta Y");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecYDelta, vecXPred);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Predicted X position in channel no");
    graph.GetYaxis()->SetTitle("Delta Y in channel no");

    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

  if (bPlot && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterHisto XObs vs YObs
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_YPred_vs_YObs";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXObs, vecYObs);
    histo->GetXaxis()->SetTitle("XObs");
    histo->GetYaxis()->SetTitle("YObs");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());    //,histName.str());
    delete histo;
  }

  if (bPlot && nAlignmentStep > -1 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterHisto DeltaX vs Chi2X
    histName.str("");
    histName.clear();
    histName << preName.str();
    histName << "_ScatterPlot_DeltaX_vs_Chi2X";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecXDelta, vecXChi2, 256);
//    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("Delta X");
    histo->GetYaxis()->SetTitle("Chi2 X");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXChi2);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Chi^2 per NDF");
    graph.GetYaxis()->SetTitle("Delta X in channel no");
    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }
  if (bPlot && nAlignmentStep > -1 && subjectPlane < 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::Y_COR)) {    //ScatterHisto DeltaY vs Chi2Y
    histName.str("");
    histName << preName.str();
    histName << "_ScatterPlot_DeltaX_vs_Chi2X";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateScatterHisto(histName.str(), vecYDelta, vecYChi2, 256);
//    histo->Draw("goff");
    histo->GetYaxis()->SetTitle("Delta Y");
    histo->GetXaxis()->SetTitle("Chi2 Y");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecYDelta, vecYChi2);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Chi^2 per NDF");
    graph.GetYaxis()->SetTitle("Delta Y in channel no");
    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }
  if (bPlot && subjectPlane == 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    //predX vs deltaX
    histName.str("");
    histName << preName.str();
    histName << "_ScatterPlot_XPred_vs_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateDipendencyHisto(histName.str(), vecXPred, vecXDelta, 512);
//    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("X Predicted");
    histo->GetYaxis()->SetTitle("Delta X");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecXPred);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Predicted X position in channel no");
    graph.GetYaxis()->SetTitle("Delta X in channel no");

    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

  if (bPlot && subjectPlane ==4){
    histName.str("");histName.clear();
    histName<<"hAngularDistribution";
    TH2F *histo = histSaver->CreateScatterHisto(histName.str(),vecXPhi,vecYPhi,512);
//    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("Phi X in Degree");
    histo->GetYaxis()->SetTitle("Phi Y in Degree");
    TH1D* hProj=histo->ProjectionX();
    int binxMin=0;
    for(binxMin=0;binxMin<hProj->GetNbinsX();binxMin++)if(hProj->GetBinContent(binxMin))break;
    int binxMax;
    for(binxMax=hProj->GetNbinsX();binxMax>0;binxMax--)if(hProj->GetBinContent(binxMax))break;
    histo->GetXaxis()->SetRangeUser(hProj->GetBinLowEdge(binxMin-1),hProj->GetBinLowEdge(binxMax+1));
    delete hProj;
    hProj=histo->ProjectionY();
    int binyMin=0;
    for(binyMin=0;binyMin<hProj->GetNbinsX();binyMin++)if(hProj->GetBinContent(binyMin))break;
    int binyMax;
    for(binyMax=hProj->GetNbinsX();binyMax>0;binyMax--)if(hProj->GetBinContent(binyMax))break;

    histo->GetYaxis()->SetRangeUser(hProj->GetBinLowEdge(binyMin-1),hProj->GetBinLowEdge(binyMax+1));
    delete hProj;
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecClusterSize);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Phi Y in Degree");
    graph.GetYaxis()->SetTitle("Phi Y in Degree");
    histSaver->SaveGraph((TGraph*)graph.Clone(),histName.str());
    delete histo;
  }
  if (bPlot && subjectPlane == 4 && (cor == TPlaneProperties::XY_COR || cor == TPlaneProperties::X_COR)) {    //DeltaX vs ClusterSize
    histName.str("");
    histName << preName.str();
    histName << "_ScatterPlot_ClusterSize_vs_DeltaX";
    histName << "_-_Plane_" << subjectPlane << "_with_" << refPlaneString;
    ;    //<<"_with"<<refPlane1<<"_and_"<<refPlane2;
    TH2F *histo = histSaver->CreateDipendencyHisto(histName.str(), vecClusterSize, vecXDelta, 512);
    histo->Draw("goff");
    histo->GetXaxis()->SetTitle("Cluster Size");
    histo->GetYaxis()->SetTitle("Delta X");
    histSaver->SaveHistogram((TH2F*) histo->Clone());
    histName << "_graph";
    TGraph graph = histSaver->CreateDipendencyGraph(histName.str(), vecXDelta, vecClusterSize);
    graph.Draw("APL");
    graph.GetXaxis()->SetTitle("Cluster Size");
    graph.GetYaxis()->SetTitle("Delta X in channel no");

    histSaver->SaveGraph((TGraph*) graph.Clone(), histName.str());
    delete histo;
  }

}

void TAlignment::DoEtaCorrectionSilicon(UInt_t correctionStep) {
  cout << "****************************************" << endl;
  cout << "Do Eta correction for all silicon planes " << correctionStep << endl;
  cout << "****************************************" << endl;
  vector<TH1F*> histoStripDistribution, histoStripDistributionFlattned;
  vector<TH1F *> vecHEta;
  stringstream fileName;
  fileName << "etaIntegral_Step" << correctionStep << "." << settings->getRunNumber() << ".root";
  TFile* correctedEtaFile = new TFile(fileName.str().c_str(), "RECREATE");
  cout << "create Histos..." << endl;
  for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
    stringstream histoTitle;
    histoTitle << "hPredictedStripPosition" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
    histoStripDistribution.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, -0.501, 0.501));
    histoTitle << "_flattend";
    histoStripDistributionFlattned.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, -0.501, 0.501));
    histoTitle.str("");
    histoTitle.clear();
    histoTitle << "hCorrectedEtaDistribution" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
    vecHEta.push_back(new TH1F(histoTitle.str().c_str(), histoTitle.str().c_str(), 128, 0, 1));
  }

  cout << "fill first strip hit histo" << events.size() << endl;

  for (nEvent = 0; nEvent < this->events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));

    for (UInt_t subjectPlane = 0; subjectPlane < TPlaneProperties::getNSiliconPlanes(); subjectPlane++) {
      //			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
      //				continue;
      vector<UInt_t> refPlanes;
      for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
        if (subjectPlane != refPlane) refPlanes.push_back(refPlane);
      TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, TCluster::corEta, false);
      Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane * 2, pred->getPositionX(), pred->getPositionY());
      Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane * 2 + 1, pred->getPositionX(), pred->getPositionY());
      UInt_t stripMiddleX = (UInt_t) (predictedStripPositionX + 0.5);
      Float_t deltaX = predictedStripPositionX - stripMiddleX;
      UInt_t stripMiddleY = (UInt_t) (predictedStripPositionY + 0.5);
      Float_t deltaY = predictedStripPositionY - stripMiddleY;
      //			cout<<nEvent<<": "<<subjectPlane<<"Fill "<<deltaX<<" "<<deltaY<<endl;
      histoStripDistribution.at(subjectPlane * 2)->Fill(deltaX);
      histoStripDistribution.at(subjectPlane * 2 + 1)->Fill(deltaY);
    }
  }
  vector<UInt_t> vecMinEntries;
  cout << "Minimal Entries in a bin of historgram:" << endl;
  for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
    TH1F* histo = histoStripDistribution.at(det);
    Int_t minBin = histo->GetMinimumBin();
    Int_t nMinEntries = histo->GetBinContent(minBin);
    cout << endl << det << ": " << minBin << " " << nMinEntries << "\t";
    vecMinEntries.push_back(nMinEntries);
  }
  cout << "\n\n" << endl;

  vector<UInt_t> vecEventNo[9];
  cout << "create flattened strip hit histo " << events.size() << endl;
  for (nEvent = 0; nEvent < this->events.size(); nEvent++) {
    TRawEventSaver::showStatusBar(nEvent, events.size());
    myTrack->setEvent(&events.at(nEvent));

    for (UInt_t subjectPlane = 0; subjectPlane < TPlaneProperties::getNSiliconPlanes(); subjectPlane++) {
      //			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
      //				continue;
      vector<UInt_t> refPlanes;
      for (UInt_t refPlane = 0; refPlane < TPlaneProperties::getNSiliconPlanes(); refPlane++)
        if (subjectPlane != refPlane) refPlanes.push_back(refPlane);
      TPositionPrediction* pred = myTrack->predictPosition(subjectPlane, refPlanes, TCluster::corEta, false);

      Float_t predictedStripPositionX = myTrack->getPositionInDetSystem(subjectPlane * 2, pred->getPositionX(), pred->getPositionY());
      Float_t predictedStripPositionY = myTrack->getPositionInDetSystem(subjectPlane * 2 + 1, pred->getPositionX(), pred->getPositionY());
      UInt_t stripMiddleX = (UInt_t) (predictedStripPositionX + 0.5);
      Float_t deltaX = predictedStripPositionX - stripMiddleX;
      //			histoStripDistribution.at(subjectPlane*2)->Fill(deltaX);
      UInt_t stripMiddleY = (UInt_t) (predictedStripPositionY + 0.5);
      Float_t deltaY = predictedStripPositionY - stripMiddleY;
      //			histoStripDistribution.at(subjectPlane*2)->Fill(deltaY);

      Int_t binX = histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaX);
      Int_t binY = histoStripDistributionFlattned.at(subjectPlane)->FindBin(deltaY);
      if (histoStripDistributionFlattned.at(subjectPlane * 2)->GetBinContent(binX) < vecMinEntries.at(subjectPlane * 2)) {
        vecEventNo[subjectPlane * 2].push_back(nEvent);
        histoStripDistributionFlattned.at(subjectPlane * 2)->Fill(deltaX);
      }
      if (histoStripDistributionFlattned.at(subjectPlane * 2 + 1)->GetBinContent(binY) < vecMinEntries.at(subjectPlane * 2 + 1)) {
        vecEventNo[subjectPlane * 2 + 1].push_back(nEvent);
        histoStripDistributionFlattned.at(subjectPlane * 2 + 1)->Fill(deltaY);
      }
    }
  }
  for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
    cout << "save histogram: " << det << "  " << histoStripDistributionFlattned.at(det)->GetTitle() << "  " << histoStripDistribution.at(det)->GetTitle() << endl;
    histSaver->SaveHistogram(histoStripDistributionFlattned.at(det));
    correctedEtaFile->Add(histoStripDistributionFlattned.at(det));
    histSaver->SaveHistogram(histoStripDistribution.at(det));
    correctedEtaFile->Add(histoStripDistribution.at(det));
  }

  cout << "\n\ncreate eta correction histo" << endl;
  for (UInt_t det = 0; det < TPlaneProperties::getNSiliconDetectors(); det++) {
    for (UInt_t i = 0; i < vecEventNo[det].size(); i++) {
      nEvent = vecEventNo[det].at(i);
      eventReader->LoadEvent(nEvent);

      if (!eventReader->useForAlignment() && !eventReader->useForAnalysis()) continue;
      Float_t eta = eventReader->getCluster(det, 0).getEta();
      vecHEta.at(det)->Fill(eta);

    }

    stringstream histName;
    histName << "hEtaIntegral" << "_step" << correctionStep << "_" << TPlaneProperties::getStringForDetector(det);
    ;
    UInt_t nBins = vecHEta.at(det)->GetNbinsX();
    TH1F *histo = new TH1F(histName.str().c_str(), histName.str().c_str(), nBins, 0, 1);
    Int_t entries = vecHEta.at(det)->GetEntries();
    entries -= vecHEta.at(det)->GetBinContent(0);
    entries -= vecHEta.at(det)->GetBinContent(nBins + 1);
    Int_t sum = 0;
    for (UInt_t bin = 1; bin < nBins + 1; bin++) {
      Int_t binContent = vecHEta.at(det)->GetBinContent(bin);
      sum += binContent;
      Float_t pos = vecHEta.at(det)->GetBinCenter(bin);
      histo->Fill(pos, (Float_t) sum / (Float_t) entries);
    }
    cout << "save " << vecHEta.at(det)->GetTitle() << " " << vecHEta.at(det)->GetEntries() << endl;
    histSaver->SaveHistogram(vecHEta.at(det));
    correctedEtaFile->Add(vecHEta.at(det));
    cout << "save " << histo->GetTitle() << " " << histo->GetEntries() << endl;
    histSaver->SaveHistogram(histo);
    correctedEtaFile->Add(histo->Clone());
  }
  correctedEtaFile->Write();
  correctedEtaFile->Close();
}

/**
 * @brief: Print Informations for all Events smaller than maxEvent, starting with startEvent
 * @param maxEvent  maximum Event No to print event
 * @param:  startEvent  first Event where to start with printing Event Information
 *
 */
void TAlignment::PrintEvents(UInt_t maxEvent, UInt_t startEvent) {
  if (maxEvent == 0) maxEvent = eventReader->GetEntries();
  if (maxEvent > eventReader->GetEntries()) maxEvent = eventReader->GetEntries();
  cout << "\n\n\n\n\n\nPrintEVENTS" << maxEvent << "\n\n\n" << flush;
  for (UInt_t event = startEvent; event < maxEvent; event++) {
    eventReader->LoadEvent(event);
    cout << event << ":\t" << flush;
    eventReader->getEvent()->Print(1);
  }
  cout << "\n\n\n\n\n\n" << flush;

}

void TAlignment::clearMeasuredVectors() {
  vecXPred.clear();
  vecYPred.clear();
  vecXObs.clear();
  vecYObs.clear();
  vecXDelta.clear();
  vecYDelta.clear();
  vecXChi2.clear();
  vecYChi2.clear();
  vecXPhi.clear();
  vecYPhi.clear();
  vecXMeasured.clear();
  vecYMeasured.clear();

}

