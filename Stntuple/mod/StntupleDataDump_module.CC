//--------------------------------------------------------------------------
// Description:
// -----------
// Class StntupleDataDump : books tree and does other initializations 
//                            for STNTUPLE
//
// Nov 23 2000 P.Murat
//------------------------------------------------------------------------
#include <string>
#include <cstdio>

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TH1.h"
#include "TProfile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTime.h"

// Framework includes.

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include <Stntuple/obj/TStnDBManager.hh>

#include <Stntuple/obj/TStnDataBlock.hh>
#include <Stntuple/obj/TStnNode.hh>
#include <Stntuple/obj/TStnEvent.hh>
#include <Stntuple/obj/TStnErrorLogger.hh>
#include <Stntuple/alg/TStntuple.hh>

#ifndef Stntuple_mod_InitStntuple
#define Stntuple_mod_InitStntuple

class TObjArray;

#include "Stntuple/mod/StntupleModule.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

using namespace std;

namespace mu2e {

class InitStntuple : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:

  string                   fMakeStrawHitModuleLabel;

  Int_t                    fLastRun;
  Float_t                  fSumInstLum; //! avg inst lum for evaluating
  Int_t                    fnLum;       //! exe speed
  Float_t                  fCpuSpeed;   //! MHz of CPU

//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  InitStntuple(fhicl::ParameterSet const& Pset);

  ~InitStntuple();
				        // ****** accessors
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run&   r);
  virtual bool filter  (AbsEvent&   e);
  virtual bool endRun (art::Run&   r);
  virtual void beginJob();
  virtual void endJob  ();
					// ****** functions of the module

  Int_t     ProcessNewRun      (art::Run* ARun);
};

//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
StntupleDataDump::StntupleDataDump(fhicl::ParameterSet const& Pset): 
  StntupleModule   (Pset,"StntupleDataDump"),
  fMakeStrawHitsModuleLabel(PSet.get<string>       ("makeStrawHitModuleLabel"))
{
//-----------------------------------------------------------------------------
// dont create subdirectories for the modules: they will have different folders
//-----------------------------------------------------------------------------
  fLastRun      = -1;
}


//------------------------------------------------------------------------------
StntupleDataDump::~StntupleDataDump() {
  // do not need to delete anything
}


//------------------------------------------------------------------------------
void StntupleDataDump::beginJob() {
}


//------------------------------------------------------------------------------
bool StntupleDataDump::beginRun(art::Run&  aRun) {
  // fetch calibration constants for a new run


  return 1;
}

//------------------------------------------------------------------------------
bool StntupleDataDump::filter(AbsEvent& AnEvent) {

  int rc = initData();
  return 1;
}

//------------------------------------------------------------------------------
bool StntupleDataDump::endRun(art::Run& ARun) {
  return true;
}

//-----------------------------------------------------------------------------
void StntupleDataDump::endJob() {
}

} //end namespace mu2e
using mu2e::StntupleDataDump;

DEFINE_ART_MODULE(StntupleDataDump);
