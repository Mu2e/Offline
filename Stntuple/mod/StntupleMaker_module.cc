///////////////////////////////////////////////////////////////////////////////
// Class StntupleMaker : fills Stntuple (P.Murat)
// ------------------------------------------
// order of the data blocks is essential - they are filled according to the
// order in which they are declared...
//
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <cstdio>

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TH1.h"
#include "TString.h"
#include "TProfile.h"
#include "TFolder.h"
#include "TSystem.h"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnErrorLogger.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"

//  #include "Stntuple/obj/TStnTriggerBlock.hh"

#include "Stntuple/mod/StntupleGlobals.hh"
#include "Stntuple/mod/StntupleMaker_module.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/mod/StntupleGlobals.hh"

// ClassImp(StntupleMaker)

static const char rcsid[] = "$Name:  $";

void stntuple_get_version(char*& ver, char*& test);
namespace mu2e {
//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
StntupleMaker::StntupleMaker(fhicl::ParameterSet const& PSet): 
  StntupleModule        (PSet,"StntupleMaker")
  , fProcessName        (PSet.get<std::string> ("processName"    ,"PROD"          ))
  , fMakeCalData        (PSet.get<int>         ("makeCalData"    ,        0       ))
  , fMakeClusters       (PSet.get<int>         ("makeClusters"   ,        1       ))
  , fMakeStrawData      (PSet.get<int>         ("makeStrawData"  ,        0       ))
  , fMakeTracks         (PSet.get<int>         ("makeTracks"     ,        1       ))
  , fMakeTracksUem      (PSet.get<int>         ("makeTracksUem"  ,        0       ))
  , fMakeTracksDmm      (PSet.get<int>         ("makeTracksDmm"  ,        0       ))
  , fMakeTracksUmm      (PSet.get<int>         ("makeTracksUmm"  ,        0       ))
  , fMakeTrigger        (PSet.get<int>         ("makeTrigger"    ,        0       ))
  , fMakeGenp           (PSet.get<int>         ("makeGenp"       ,        1       ))
  , fMakeSimp           (PSet.get<int>         ("makeSimp"       ,        1       ))
  , fMakeVirtualHits    (PSet.get<int>         ("makeVirtualHits",        0       ))

  , fG4ModuleLabel      (PSet.get<std::string> ("g4ModuleLabel"      , "g4run"        ))
  , fStrawHitMaker      (PSet.get<std::string> ("strawHitMaker"      , "makeSH"       ))
  , fTrkPatRecDem       (PSet.get<std::string> ("trkPatRecDem"       , "trkPatRecDem" ))
  , fTrkPatRecUem       (PSet.get<std::string> ("trkPatRecUem"       , "trkPatRecUem" ))
  , fTrkPatRecDmm       (PSet.get<std::string> ("trkPatRecDmm"       , "trkPatRecDmm" ))
  , fTrkPatRecUmm       (PSet.get<std::string> ("trkPatRecUmm"       , "trkPatRecUmm" ))
  , fCaloCrystalHitMaker(PSet.get<std::string> ("caloCrystalHitMaker", "CaloCrystalHitsMaker"))
  , fCaloClusterMaker   (PSet.get<std::string> ("caloClusterMaker"   , "makeCaloCluster"))
  , fTrkExtrapol        (PSet.get<std::string> ("trkExtrapol"        , "trkExtrapol"  ))
  , fTrkCalMatch        (PSet.get<std::string> ("trkCalMatch"        , "trackCaloMatching"))
  , fPidDem             (PSet.get<std::string> ("pidDem"             , "undefined"    ))

  , fTrkPatRecDem2       (PSet.get<std::string> ("trkPatRecDem2"       , "trkPatRecDemHack" ))
  , fTrkPatRecUem2       (PSet.get<std::string> ("trkPatRecUem2"       , "trkPatRecUemHack" ))
  , fCaloCrystalHitMaker2(PSet.get<std::string> ("caloCrystalHitMaker2", "CaloCrystalHitsMakerHack"))
  , fCaloClusterMaker2   (PSet.get<std::string> ("caloClusterMaker2"   , "makeCaloClusterHack"))
  , fTrkExtrapol2        (PSet.get<std::string> ("trkExtrapol2"        , "trkExtrapolHack"    ))
  , fTrkCalMatch2        (PSet.get<std::string> ("trkCalMatch2"        , "caloMatchingHack"   ))
  , fPidDem2             (PSet.get<std::string> ("pidDem2"             , "undefined"          ))

  , fMinTActive         (PSet.get<double>      ("minTActive"     ,   710.         ))
  , fMinECrystal        (PSet.get<double>      ("minECrystal"    ,    0.1         ))
{

  char  *ver, *text;
  stntuple_get_version(ver,text);

  fVersion      = new TNamed(ver,text);
  TModule::fFolder->Add(fVersion);

  fgTimeOffsets = new SimParticleTimeOffset(PSet.get<fhicl::ParameterSet>("TimeOffsets"));

}


//------------------------------------------------------------------------------
StntupleMaker::~StntupleMaker() {
  delete fgTimeOffsets;
  delete fVersion;
}


//------------------------------------------------------------------------------
bool StntupleMaker::beginRun(art::Run& aRun) {

  static int first_begin_run = 1;

  THistModule::beforeBeginRun(aRun);

  if (first_begin_run) {
//-----------------------------------------------------------------------------
// if we runnning stnmaker_prod.exe, save revision of the TCL file in STNTUPLE
//-----------------------------------------------------------------------------
    first_begin_run = 0;
    const char* c = gSystem->Getenv("STNMAKER_PROD_TCL");
    if (c) TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL",c));
    else   TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL","unknown"));
  }

  THistModule::afterBeginRun(aRun);

  return 1;
}

//------------------------------------------------------------------------------
bool StntupleMaker::endRun(art::Run& aRun ) {
  THistModule::beforeEndRun(aRun);
  THistModule::afterEndRun (aRun);
  return 1;
}


//------------------------------------------------------------------------------
void StntupleMaker::endJob() {

  THistModule::beforeEndJob();
  THistModule::afterEndJob ();

}

//------------------------------------------------------------------------------
void StntupleMaker::beginJob() {

  int split_mode, compression_level, buffer_size;

  THistModule::beforeBeginJob();

  // create data blocks and branches

  fgStntupleFolder->Add(new TNamed("ProcessName"     ,fProcessName));

					// for the moment do it by hands...
					// create default branches to go into 
					// STNTUPLE

  split_mode        = THistModule::SplitLevel();
  compression_level = THistModule::CompressionLevel();
  buffer_size       = THistModule::BufferSize();
//-----------------------------------------------------------------------------
// calorimeter hit data
// this is not RAW hit data yet...
//-----------------------------------------------------------------------------
  if (fMakeCalData) {
    TStnDataBlock* cal_data;

    cal_data = AddDataBlock("CalDataBlock","TCalDataBlock",
			    StntupleInitMu2eCalDataBlock,
			    buffer_size,
			    split_mode,
			    compression_level);
    if (cal_data) {
      cal_data->AddCollName("mu2e::CaloCrystalHitCollection",fCaloCrystalHitMaker.data(),"");
    }
  }
//-----------------------------------------------------------------------------
// straw hit data
//-----------------------------------------------------------------------------
  if (fMakeStrawData) {
    TStnDataBlock* straw_data;

    straw_data = AddDataBlock("StrawDataBlock","TStrawDataBlock",
			      StntupleInitMu2eStrawDataBlock,
			      buffer_size,
			      split_mode,
			      compression_level);
    if (straw_data) {
      straw_data->AddCollName("mu2e::StrawHitCollection",fStrawHitMaker.data(),"");
    }
  }
//-----------------------------------------------------------------------------
// 2 trigger data branches: one for the trigger data, 
// another one for the emulated trigger data
// always NON-SPLIT
//-----------------------------------------------------------------------------
  // if (fMakeTrigger.value()) {
  //   TStnTriggerBlock *data;
  //   data = (TStnTriggerBlock*) AddDataBlock("TriggerBlock","TStnTriggerBlock",
  // 					    StntupleInitTriggerBlock,
  // 					    buffer_size.value(),
  // 					    -1,
  // 					    compression_level.value());
  //   data->SetL3Source(fL3Source.value());
  // }
//-----------------------------------------------------------------------------
// track branches: for ROOT v3 to use streamers one has to specify split=-1
//-----------------------------------------------------------------------------
  if (fMakeTracks) {
    TStnDataBlock* track_data;
					// always store defTracks for the 
					// default process in the "TrackBlock"

    track_data = AddDataBlock("TrackBlock",
			      "TStnTrackBlock",
			      StntupleInitMu2eTrackBlock,
			      buffer_size,
			      split_mode,
			      compression_level);

    SetResolveLinksMethod("TrackBlock",StntupleInitMu2eTrackBlockLinks);

    if (track_data) {
      track_data->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecDem.data()    ,"DownstreameMinus");
      track_data->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
      track_data->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol.data()     ,"");
      //      track_data->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch.data()     ,"");
      track_data->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCalMatch.data()     ,"");
      track_data->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
      track_data->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
      track_data->AddCollName("mu2e::PIDProductCollection"          ,fPidDem.data()          ,"");
      track_data->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");
    }
  }

//-----------------------------------------------------------------------------
// UEM: Upstream Electron Minus
//-----------------------------------------------------------------------------
  if (fMakeTracksUem) {
    TStnDataBlock* track_data_uem;

    track_data_uem = AddDataBlock("TrackBlockUem",
				  "TStnTrackBlock",
				  StntupleInitMu2eTrackBlock,
				  buffer_size,
				  split_mode,
				  compression_level);

    SetResolveLinksMethod("TrackBlockUem",StntupleInitMu2eTrackBlockLinks);

    if (track_data_uem) {
      track_data_uem->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");
      track_data_uem->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecUem.data()    ,"UpstreameMinus");
      track_data_uem->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
      track_data_uem->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
      track_data_uem->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
//-----------------------------------------------------------------------------
// for the moment, don't run extrapolation and PID for these
//-----------------------------------------------------------------------------
//       track_data_uem->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::PIDProductCollection"          ,fPidDem2.data()         ,"");
    }
  }
//-----------------------------------------------------------------------------
// DMM: Downstream Muon Minus
//-----------------------------------------------------------------------------
  if (fMakeTracksDmm) {
    TStnDataBlock* track_data_dmm;

    track_data_dmm = AddDataBlock("TrackBlockDmm",
				  "TStnTrackBlock",
				  StntupleInitMu2eTrackBlock,
				  buffer_size,
				  split_mode,
				  compression_level);

    SetResolveLinksMethod("TrackBlockDmm",StntupleInitMu2eTrackBlockLinks);

    if (track_data_dmm) {
      track_data_dmm->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecDmm.data()    ,"DownstreammuMinus");
      track_data_dmm->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");
      track_data_dmm->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
      track_data_dmm->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
      track_data_dmm->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
//-----------------------------------------------------------------------------
// for the moment, don't run extrapolation and PID for these
//-----------------------------------------------------------------------------
//       track_data_uem->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::PIDProductCollection"          ,fPidDem2.data()         ,"");
    }
  }
//-----------------------------------------------------------------------------
// UMM: Upstream Muon Minus
//-----------------------------------------------------------------------------
  if (fMakeTracksUmm) {
    TStnDataBlock* track_data_umm;

    track_data_umm = AddDataBlock("TrackBlockUmm",
				  "TStnTrackBlock",
				  StntupleInitMu2eTrackBlock,
				  buffer_size,
				  split_mode,
				  compression_level);

    SetResolveLinksMethod("TrackBlockUmm",StntupleInitMu2eTrackBlockLinks);

    if (track_data_umm) {
      track_data_umm->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecUmm.data()    ,"UpstreammuMinus");
      track_data_umm->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");
      track_data_umm->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
      track_data_umm->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
      track_data_umm->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
//-----------------------------------------------------------------------------
// for the moment, don't run extrapolation and PID for these
//-----------------------------------------------------------------------------
//       track_data_uem->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch2.data()    ,"");
//       track_data_uem->AddCollName("mu2e::PIDProductCollection"          ,fPidDem2.data()         ,"");
    }
  }
//-----------------------------------------------------------------------------
// clusters 
//-----------------------------------------------------------------------------
  if (fMakeClusters) {
    TStnDataBlock* cluster_data;

    cluster_data = AddDataBlock("ClusterBlock",
				"TStnClusterBlock",
				StntupleInitMu2eClusterBlock,
				buffer_size,
				split_mode,
				compression_level);

    SetResolveLinksMethod("ClusterBlock",StntupleInitMu2eClusterBlockLinks);

    if (cluster_data) {
      cluster_data->AddCollName("mu2e::CaloClusterCollection",fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    }
  }
//-----------------------------------------------------------------------------
// generator particles 
//-----------------------------------------------------------------------------
  if (fMakeGenp) {
    TStnDataBlock* genp_data;

    genp_data = AddDataBlock("GenpBlock",
			     "TGenpBlock",
			     StntupleInitMu2eGenpBlock,
			     buffer_size,
			     split_mode,
			     compression_level);

    //    SetResolveLinksMethod("GenpBlock",StntupleInitMu2eClusterBlockLinks);

    if (genp_data) {
      genp_data->AddCollName("mu2e::GenParticleCollection","","");
    }
  }
//-----------------------------------------------------------------------------
// simulated particles 
//-----------------------------------------------------------------------------
  if (fMakeSimp) {
    TStnDataBlock* simp_data;

    simp_data = AddDataBlock("SimpBlock",
			     "TSimpBlock",
			     StntupleInitMu2eSimpBlock,
			     buffer_size,
			     split_mode,
			     compression_level);

    //    SetResolveLinksMethod("GenpBlock",StntupleInitMu2eClusterBlockLinks);

    if (simp_data) {
      simp_data->AddCollName("mu2e::SimParticleCollection",""                   ,"");
      simp_data->AddCollName("mu2e::StrawHitCollection"   ,fStrawHitMaker.data(),"");
      simp_data->AddCollName("mu2e::StepPointMCCollection",fG4ModuleLabel.data(),"");
    }
  }
//-----------------------------------------------------------------------------
// virtual detectior hits
//-----------------------------------------------------------------------------
  if (fMakeVirtualHits) {
    TStnDataBlock* virtual_data;

    virtual_data = AddDataBlock("VdetBlock",
				"TVdetDataBlock",
				StntupleInitMu2eVirtualDataBlock,
				buffer_size,
				split_mode,
				compression_level);

    if (virtual_data) {
      virtual_data->AddCollName("mu2e::StepPointMCCollection", fG4ModuleLabel.data(),"virtualdetector");
    }
  }  


  THistModule::afterEndJob();
}

//_____________________________________________________________________________
bool StntupleMaker::filter(AbsEvent& AnEvent) {

  // when execution comes here al the registered data blocks are already
  // initialized with the event data. Left: variables in the data blocks
  // which depend on the variable defined in other blocks, like track number
  // for a muon or an electron - the idea is that these are defined during the
  // 2nd loop in FillStntupleModule, where ResolveLinks methods are called
  // for each data block

//-----------------------------------------------------------------------------
// connect to the error reporting facility
//-----------------------------------------------------------------------------
  TStnErrorLogger* logger = Event()->GetErrorLogger();
//   logger->Connect("Report(Int_t, const char*)",
// 		  "StntupleModule",
// 		  this,
// 		  "LogError(const char*)");
//-----------------------------------------------------------------------------
// disconnect from the error reporting signal and return back to AC++
//-----------------------------------------------------------------------------
//   logger->Disconnect("Report(Int_t,const char*)",
// 		     this,"LogError(Int_t,const char*)");

  return 1;
}


//------------------------------------------------------------------------------
void StntupleMaker::GetDefTrackCollName(char* Name) {
  // put in a working kludge first

  strcpy(Name,"default");
}


// //_____________________________________________________________________________
// int StntupleMaker::InitCalDataBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eCalDataBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitHeaderBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eHeaderBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitTrackBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTrackBlock(Block,event,mode);
// }

//_____________________________________________________________________________
// int StntupleMaker::InitTriggerBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTriggerBlock(Block,event,mode);
// }

} // end namespace mu2e

using mu2e::StntupleMaker;

DEFINE_ART_MODULE(StntupleMaker);
