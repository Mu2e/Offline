/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: IlcDCHReconstructor.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for DCH reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "IlcDCHReconstructor.h"
//#include "IlcRunLoader.h"
//#include "IlcRun.h"
//#include "IlcRawReader.h"
#include "IlcDCHtracker.h"
#include "IlcDCHParam.h"
//#include "IlcDCH.h"

//ClassImp(IlcDCHReconstructor)


IlcDCHRecoParam *    IlcDCHReconstructor::fgkRecoParam =0;  // reconstruction parameters
Int_t    IlcDCHReconstructor::fgStreamLevel     = 0;        // stream (debug) level


IlcDCHReconstructor::IlcDCHReconstructor()/*: IlcReconstructor()*/ {
  //
  // default constructor
  //
  if (!fgkRecoParam) {
    //IlcError("The Reconstruction parameters nonitiilczed - Used default one");
    fgkRecoParam = IlcDCHRecoParam::GetHighFluxParam();
  }
}

////_____________________________________________________________________________
//void IlcDCHReconstructor::Reconstruct(IlcRunLoader* runLoader) const
//{
//// reconstruct clusters
//
//  IlcLoader* loader = runLoader->GetLoader("DCHLoader");
//  if (!loader) {
//     Error("Reconstruct", "DCH loader not found");
//    return;
//  }
//  loader->LoadRecPoints("recreate");
//  loader->LoadDigits("read");
//
//  IlcDCHParam* param = GetDCHParam(runLoader);
//  if (!param) return;
////   IlcDCHclusterer clusterer(param);
//  Int_t nEvents = runLoader->GetNumberOfEvents();
//
//  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
//    runLoader->GetEvent(iEvent);
//
//    TTree* treeClusters = loader->TreeR();
//    if (!treeClusters) {
//      loader->MakeTree("R");
//      treeClusters = loader->TreeR();
//    }
//    TTree* treeDigits = loader->TreeD();
//    if (!treeDigits) {
//      Error("Reconstruct", "Can't get digits tree !");
//      return;
//    }
//
////     clusterer.SetInput(treeDigits);
////     clusterer.SetOutput(treeClusters);
////     clusterer.Digits2Clusters();
//
//    loader->WriteRecPoints("OVERWRITE");
//  }
//
//  loader->UnloadRecPoints();
//  loader->UnloadDigits();
//}
//
////_____________________________________________________________________________
//void IlcDCHReconstructor::Reconstruct(IlcRunLoader* runLoader,
//				      IlcRawReader* rawReader) const
//{
//// reconstruct clusters from raw data
//
//  IlcLoader* loader = runLoader->GetLoader("DCHLoader");
//  if (!loader) {
//    // Error("Reconstruct", "DCH loader not found");
//    return;
//  }
//  loader->LoadRecPoints("recreate");
//
//  IlcDCHParam* param = GetDCHParam(runLoader);
//  if (!param) {
//    IlcWarning("Loading default DCH parameters !");
//    param = new IlcDCHParam;
//  }
////   IlcDCHclusterer clusterer(param);
//
//  TString option = GetOption();
//  //  if (option.Contains("PedestalSubtraction"))
//  //  clusterer.SetPedSubtraction(kTRUE);
//  //if (option.Contains("OldRCUFormat"))
//  //  clusterer.SetOldRCUFormat(kTRUE);
//
//  Int_t iEvent = 0;
//  while (rawReader->NextEvent()) {
//    runLoader->GetEvent(iEvent++);
//
//    TTree* treeClusters = loader->TreeR();
//    if (!treeClusters) {
//      loader->MakeTree("R");
//      treeClusters = loader->TreeR();
//    }
//
////     clusterer.SetOutput(treeClusters);
////     clusterer.Digits2Clusters(rawReader);
//
//    loader->WriteRecPoints("OVERWRITE");
//  }
//
//  loader->UnloadRecPoints();
//
//}
//
////_____________________________________________________________________________
//IlcTracker* IlcDCHReconstructor::CreateTracker(IlcRunLoader* runLoader) const
//{
//// create a DCH tracker
//
//  IlcDCHParam* param = GetDCHParam(runLoader);
//  if (!param) {
//    IlcWarning("Loading default DCH parameters !");
//    param = new IlcDCHParam;
//  }
//  //  param->ReadGeoMatrices();
//  return new IlcDCHtracker(param);
//}
//
////_____________________________________________________________________________
//void IlcDCHReconstructor::FillESD(IlcRunLoader* /*runLoader* /,
//				  IlcESD* /*esd* /) const
//{
//// make PID
//
////   Double_t parDCH[] = {47., 0.10, 10.};
////   IlcDCHpidESD tpcPID(parDCH);
////   tpcPID.MakePID(esd);
//}
//
//
////_____________________________________________________________________________
//IlcDCHParam* IlcDCHReconstructor::GetDCHParam(IlcRunLoader* runLoader) const
//{
//// get the DCH parameters
//
//  TDirectory* saveDir = gDirectory;
//
//  IlcRun* ilcRun = runLoader->GetIlcRun();
//  if (!ilcRun) {
//    IlcError("Couldn't get IlcRun object");
//    return 0;
//  }
//  IlcDCH* dch = (IlcDCH*) ilcRun->GetDetector("DCH");
//  if (!dch) {
//    IlcError("Couldn't get DCH detector");
//    return 0;
//  }
//  IlcDCHParam* param = dch->GetParam();
//
//  //runLoader->CdGAFile();
//  //  IlcDCHParam* param = (IlcDCHParam*) gDirectory->Get("75x40_100x60_150x60");
//  // if (!param) Error("GetDCHParam", "no DCH parameters found");
//
//  saveDir->cd();
//  return param;
//}
