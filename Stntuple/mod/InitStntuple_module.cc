//--------------------------------------------------------------------------
// Description:
// -----------
// Class InitStntuple : books tree and does other initializations 
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

#include "Stntuple/mod/InitStntuple_module.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"

// ClassImp(InitStntuple)

namespace mu2e {
//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
InitStntuple::InitStntuple(fhicl::ParameterSet const& Pset): 
  StntupleModule   (Pset,"InitStntuple")
{
//-----------------------------------------------------------------------------
// dont create subdirectories for the modules: they will have different folders
//-----------------------------------------------------------------------------
  THistModule::fgMakeSubdirs = 0;

  fLastRun      = -1;
}


//------------------------------------------------------------------------------
InitStntuple::~InitStntuple() {
  // do not need to delete anything
}


//------------------------------------------------------------------------------
void InitStntuple::beginJob() {

  THistModule::beforeBeginJob();

  // book the tree, for ROOT 3 kludge split mode: set it always to
  // "non-split,old"
  // header block, however is always written in split mode

  fgTree   = new TTree("STNTUPLE", "STNTUPLE");

  AddDataBlock("HeaderBlock","TStnHeaderBlock",
	       StntupleInitMu2eHeaderBlock,
	       THistModule::BufferSize(),
	       0, // 99,                          // fSplitMode.value(), always split
	       THistModule::CompressionLevel());

  SetResolveLinksMethod("HeaderBlock",StntupleInitMu2eHeaderBlockLinks);

  fnLum = 0;
  fSumInstLum = 0.0;

  FILE* pipe;
  pipe = gSystem->OpenPipe(
     "cat /proc/cpuinfo | grep MHz | tail -1 | awk '{print $4}'","r");
  fscanf(pipe,"%f",&fCpuSpeed);
  gSystem->ClosePipe(pipe);

  THistModule::afterBeginJob();

//  return AppResult::OK;
}


//------------------------------------------------------------------------------
Int_t InitStntuple::ProcessNewRun(art::Run* ARun) {

  // InitRunSummary     ();
  // InitTriggerTable   ();
  // InitDeadList       ();

  TStntuple::Init(ARun->run());

  return 0;
}

//------------------------------------------------------------------------------
bool InitStntuple::beginRun(art::Run&  aRun) {
  // fetch calibration constants for a new run

  THistModule::beforeBeginRun(aRun);

  int runnum  = aRun.run();
  //  int mc_flag = 1; // env->monteFlag();

  if (runnum != fLastRun) {
    ProcessNewRun(&aRun);
    fLastRun = runnum;
  }

  THistModule::afterBeginRun(aRun);

  return 1;
}

//------------------------------------------------------------------------------
bool InitStntuple::filter(AbsEvent& AnEvent) {
  // event entry point: initialize all the registered data blocks with the event 
  // data
  // it is assumed that the tree itself is filled in FillStntupleModule
  // order in which branches are filled may be important - for example,
  // it is better to fill track branch before the electron branch,
  // missing Et branch logically is the last one
  // assume that InitStntuple is executed before any other STNTUPLE-related
  // module
  // it decides whether we are about to close the file and to open a new one

  THistModule::beforeEvent(AnEvent);
//-----------------------------------------------------------------------------
// connect to the error reporting facility
//-----------------------------------------------------------------------------
//  TStnErrorLogger* logger = Event()->GetErrorLogger();

  //  printf(">>> InitStntuple::filter: Couldn't connect the ErrorLogger\n");

//   logger->Connect("Report(Int_t, const char*)",
// 		  "StntupleModule",
// 		  this,
// 		  "LogError(const char*)");
//-----------------------------------------------------------------------------
// initialization
//-----------------------------------------------------------------------------
  unsigned long etime = (unsigned long)(gSystem->Now());
  Event()->Init(&AnEvent,0);
  etime = (unsigned long)(gSystem->Now()) - etime;

  //compute avg inst lum
  TStnHeaderBlock* fHeaderBlock = 
    (TStnHeaderBlock*) Event()->GetDataBlock("HeaderBlock");
  if(fHeaderBlock) {
    float ilum = fHeaderBlock->InstLum()*1.0e-30;
    if(ilum>0.1 && ilum < 10000.0) {
      fSumInstLum += ilum;
      fnLum++;
    }
    // store the cpu speed in units of int(GHz*5)
    int speed = int(fCpuSpeed/1000.0*5.0);
    if(speed>255) speed = 255;
    //store event filling time in s*10, TTime is in ms
    int ietime = int(float(etime)/100.0);
    if(ietime>=(1<<24)) ietime=((1<<24)-1);
    fHeaderBlock->fCpu = (ietime<<8 | speed);
  }
//-----------------------------------------------------------------------------
// disconnect from the error reporting signal and return back to AC++
//-----------------------------------------------------------------------------
//   logger->Disconnect("Report(Int_t,const char*)",
// 		     this,"LogError(const char*)");

  THistModule::afterEvent(AnEvent);

  return 1;
}

//------------------------------------------------------------------------------
bool InitStntuple::endRun(art::Run& ARun) {
  THistModule::beforeEndRun(ARun);
  THistModule::afterEndRun (ARun);
  return true;
}

//-----------------------------------------------------------------------------
void InitStntuple::endJob() {

  THistModule::beforeEndJob();

  if(fnLum>0) {
    std::cout <<"InitStntuple::endJob avg inst lum = "
	 <<fSumInstLum/fnLum << " e30 " << std::endl;
  } else {
    std::cout <<"InitStntuple::endJob avg inst lum = "
	 << " undefined " << std::endl;
  }

  THistModule::afterEndJob();
}

} //end namespace mu2e
using mu2e::InitStntuple;

DEFINE_ART_MODULE(InitStntuple);
