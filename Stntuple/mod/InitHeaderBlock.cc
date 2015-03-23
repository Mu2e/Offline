///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TFolder.h"
#include "TNamed.h"

#include <Stntuple/obj/TStnEvent.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/mod/InitStntupleDataBlocks.hh>
#include <Stntuple/mod/StntupleUtilities.hh>

void stntuple_get_version(char*& ver, char*& test);

//_____________________________________________________________________________
Int_t StntupleInitMu2eHeaderBlock(TStnDataBlock* block, AbsEvent* AnEvent, int mode) 
{
  // Run II version, section number is defined

  static TFolder*         fol              = NULL;
  static TNamed*          processName      = NULL;

  if (! fol) {
//-----------------------------------------------------------------------------
//  initialize local static variables
//-----------------------------------------------------------------------------
    fol     = (TFolder*) gROOT->GetRootFolder()->FindObject("Stntuple");
    processName      = (TNamed*) fol->FindObject("ProcessName");
  }

  TStnHeaderBlock* header = (TStnHeaderBlock*) block;

  header->fEventNumber   = AnEvent->event();
  header->fRunNumber     = AnEvent->run();
  header->fSectionNumber = AnEvent->subRun();
  header->fMcFlag        = 1.; // gblEnv->monteFlag();
  header->fInstLum       = -1;

  header->fNTracks       = -1;
//-----------------------------------------------------------------------------
//  instantaneous luminosity
//-----------------------------------------------------------------------------
  char  *ver, *text;
  stntuple_get_version(ver,text);
  header->fStnVersion = ver;

  return 0;
}

//_____________________________________________________________________________
Int_t StntupleInitMu2eHeaderBlockLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // Mu2e version, section number is defined

  Int_t  ev_number, rn_number;

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (! Block->Initialized(ev_number,rn_number)) return -1;

					// do not do initialize links 2nd time

  if (Block->LinksInitialized()) return 0;

  TStnHeaderBlock* header = (TStnHeaderBlock*) Block;
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  header->fLinksInitialized = 1;

  return 0;
}


