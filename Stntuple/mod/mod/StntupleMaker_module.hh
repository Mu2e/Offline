//--------------------------------------------------------------------------
// Class: StntupleMaker 
//
// Environment: Software developed for the CDF at FNAL.
//
// Copyright Information: 
//	Copyright (C) 1997		Fermilab
//
//  revision history:
//  -----------------
//------------------------------------------------------------------------
#ifndef Stntuple_mod_StntupleMaker
#define Stntuple_mod_StntupleMaker

#include "TNamed.h"

#include "Stntuple/mod/StntupleModule.hh"
namespace mu2e {
class StntupleMaker : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
					// process name, default - PROD
  std::string   fProcessName;

  int           fMakeCalData;
  int           fMakeClusters;
  int           fMakeStrawData;
  int           fMakeTracks;
  int           fMakeTracks2;
  int           fMakeTrigger;
  int           fMakeGenp;
  int           fMakeSimp;
//-----------------------------------------------------------------------------
// module parameters
//-----------------------------------------------------------------------------
  std::string   fG4ModuleLabel;
  std::string   fStrawHitMaker;
  std::string   fTrkPatRecDem;
  std::string   fTrkPatRecUem;
  std::string   fCaloCrystalHitMaker;
  std::string   fCaloClusterMaker;
  std::string   fTrkExtrapol;
  std::string   fTrkCalMatch;
  std::string   fPidDem;
  
  std::string   fStrawHitMaker2;
  std::string   fTrkPatRecDem2;
  std::string   fTrkPatRecUem2;
  std::string   fCaloCrystalHitMaker2;
  std::string   fCaloClusterMaker2;
  std::string   fTrkExtrapol2;
  std::string   fTrkCalMatch2;
  std::string   fPidDem2;
  
  double        fMinTActive  ;  // start of the active window
  double        fMinECrystal ;  // 

  TNamed*       fVersion;
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

  // static int InitCalDataBlock (TStnDataBlock* Block);
  // static int InitElectronBlock(TStnDataBlock* Block);
  // static int InitTrackBlock   (TStnDataBlock* Block);
  //  static int InitTriggerBlock (TStnDataBlock* Block);

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
} // end namespace mu2e

#endif
