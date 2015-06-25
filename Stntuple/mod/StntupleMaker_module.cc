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

#include "fhiclcpp/ParameterSet.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TNamed.h"
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

#include "Stntuple/mod/StntupleModule.hh"
#include "Stntuple/mod/StntupleGlobals.hh"
// #include "Stntuple/mod/StntupleMaker_module.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

#include "KalmanTests/inc/TrkFitDirection.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/mod/StntupleGlobals.hh"

using namespace std; 

// ClassImp(StntupleMaker)

static const char rcsid[] = "$Name:  $";

void stntuple_get_version(char*& ver, char*& test);
namespace mu2e {
class StntupleMaker : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
					// process name, default - PROD
  std::string      fProcessName;
					// switches for individual branches
  int              fMakeCalData;
  int              fMakeClusters;
  int              fMakeStrawData;
  int              fMakeTracks;
  int              fMakeTrigger;
  int              fMakeGenp;
  int              fMakeSimp;
  int              fMakeVirtualHits;
//-----------------------------------------------------------------------------
// module parameters
// generator
//-----------------------------------------------------------------------------
  std::string              fG4ModuleLabel;
  std::string              fGeneratorModuleLabel;   // defines collection to save, default: "" (all)
  std::string              fMakeStrawHitModuleLabel;

  std::vector<std::string> fTrackBlockName;
  std::vector<std::string> fTrkRecoModuleLabel;
  std::vector<std::string> fTrkExtrapolModuleLabel;
  std::vector<std::string> fTrkCaloMatchModuleLabel;
  std::vector<std::string> fPidModuleLabel;

  std::vector<int>         fFitParticle;
  std::vector<int>         fFitDirection;

  std::string              fCaloCrystalHitMaker;
  std::string              fCaloClusterMaker;
  
  double                   fMinTActive  ;  // start of the active window
  double                   fMinECrystal ;  // 

  TNamed*                  fVersion;
//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  StntupleMaker(fhicl::ParameterSet const& pset);

  ~StntupleMaker();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  void GetDefTrackCollName(char* Name);

					// ****** setters

//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run& ARun);
  virtual bool endRun  (art::Run& ARun);
  virtual void beginJob();
  virtual void endJob  ();
  virtual bool filter  (AbsEvent& event);

  //  ClassDef(StntupleMaker,0)
};



//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
StntupleMaker::StntupleMaker(fhicl::ParameterSet const& PSet): 
  StntupleModule        (PSet,"StntupleMaker")
  , fProcessName        (PSet.get<string> ("processName"    ,"PROD"          ))
  , fMakeCalData        (PSet.get<int>         ("makeCalData"    ,        0       ))
  , fMakeClusters       (PSet.get<int>         ("makeClusters"   ,        1       ))
  , fMakeStrawData      (PSet.get<int>         ("makeStrawData"  ))
  , fMakeTracks         (PSet.get<int>         ("makeTracks"     ))
  , fMakeTrigger        (PSet.get<int>         ("makeTrigger"    ,        0       ))
  , fMakeGenp           (PSet.get<int>         ("makeGenp"       ,        1       ))
  , fMakeSimp           (PSet.get<int>         ("makeSimp"       ,        1       ))
  , fMakeVirtualHits    (PSet.get<int>         ("makeVirtualHits",        0       ))
  
  , fG4ModuleLabel          (PSet.get<string>        ("g4ModuleLabel"          ))
  , fGeneratorModuleLabel   (PSet.get<string>        ("generatorModuleLabel"   ))
  , fMakeStrawHitModuleLabel(PSet.get<string>        ("makeStrawHitModuleLabel"))
  , fTrackBlockName         (PSet.get<vector<string>>("trackBlockName"         ))
  , fTrkRecoModuleLabel     (PSet.get<vector<string>>("trkRecoModuleLabel"     ))
  , fTrkExtrapolModuleLabel (PSet.get<vector<string>>("trkExtrapolModuleLabel" ))
  , fTrkCaloMatchModuleLabel(PSet.get<vector<string>>("trkCaloMatchModuleLabel"))
  , fPidModuleLabel         (PSet.get<vector<string>>("pidModuleLabel"         ))
  
  , fFitParticle            (PSet.get<vector<int>>        ("fitParticle"       ))
  , fFitDirection           (PSet.get<vector<int>>        ("fitDirection"      ))

  , fCaloCrystalHitMaker(PSet.get<string> ("caloCrystalHitsMaker"))
  , fCaloClusterMaker   (PSet.get<string> ("caloClusterMaker"    ))
  
  , fMinTActive         (PSet.get<double>      ("minTActive"     ))
  , fMinECrystal        (PSet.get<double>      ("minECrystal"    ))
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

  string      _iname1;	// data instance name

  THistModule::beforeBeginJob();

  // create data blocks and branches

  fgStntupleFolder->Add(new TNamed("ProcessName"     ,fProcessName));

					// for the moment do it by hands...
					// create default branches to go into 
					// STNTUPLE

  split_mode        = THistModule::SplitLevel();
  compression_level = THistModule::CompressionLevel();
  buffer_size       = THistModule::BufferSize();

  //  _iname1           = fFitDirection1.name() + fFitParticle1.name();
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
      straw_data->AddCollName("mu2e::StrawHitCollection",fMakeStrawHitModuleLabel.data(),"");
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
  TStnDataBlock *track_data;
  const char    *block_name;
  int            nblocks;
  string         iname;

  if (fMakeTracks) {
    nblocks = fTrackBlockName.size();

    for (int i=0; i<nblocks; i++) {
					// always store defTracks for the 
					// default process in the "TrackBlock"

      block_name = fTrackBlockName[i].data();
      track_data = AddDataBlock(block_name,
				"TStnTrackBlock",
				StntupleInitMu2eTrackBlock,
				buffer_size,
				split_mode,
				compression_level);

      SetResolveLinksMethod(block_name,StntupleInitMu2eTrackBlockLinks);

      if (track_data) {
	TrkFitDirection fit_dir = TrkFitDirection((TrkFitDirection::FitDirection)fFitDirection[i]);
	TrkParticle     part    = TrkParticle((TrkParticle::type)fFitParticle[i]);
	iname                   = fit_dir.name() + part.name();

	track_data->AddCollName("mu2e::KalRepCollection"              ,fTrkRecoModuleLabel[i].data()     ,iname.data());
	track_data->AddCollName("mu2e::StrawHitCollection"            ,fMakeStrawHitModuleLabel.data()   ,"");
	track_data->AddCollName("mu2e::StrawDigiMCCollection"         ,fMakeStrawHitModuleLabel.data()   ,"StrawHitMC");
	track_data->AddCollName("mu2e::PtrStepPointMCVectorCollection",fMakeStrawHitModuleLabel.data()   ,"StrawHitMCPtr");
	track_data->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapolModuleLabel [i].data(),"");
	track_data->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data()          ,"");
	track_data->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCaloMatchModuleLabel[i].data(),"");
	track_data->AddCollName("mu2e::PIDProductCollection"          ,fPidModuleLabel[i].data()         ,"");
	track_data->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()             ,"");
      }
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
      cluster_data->AddCollName("mu2e::CaloClusterCollection",fCaloClusterMaker.data(),"");
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

    genp_data->AddCollName("mu2e::GenParticleCollection",fGeneratorModuleLabel.data(),"");
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
      simp_data->AddCollName("mu2e::StrawHitCollection"   ,fMakeStrawHitModuleLabel.data(),"");
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
//  TStnErrorLogger* logger = Event()->GetErrorLogger();
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
