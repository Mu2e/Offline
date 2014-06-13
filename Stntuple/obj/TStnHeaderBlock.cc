///////////////////////////////////////////////////////////////////////////////
//  always written out in split mode, so Streamer should never be called
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <iostream>
#include <iomanip>

#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

ClassImp(TStnHeaderBlock)

//______________________________________________________________________________
void TStnHeaderBlock::ReadV50(TBuffer &R__b)
{
  // first version of the header block had label V50 on it
  //
  R__b >> fVersion;
  R__b >> fEventNumber;
  R__b >> fRunNumber;
  R__b >> fInstLum;
  R__b >> fGoodRun;
  R__b >> fBrCode;
  R__b >> fGoodTrig;
  R__b >> fTrigWord;
//-----------------------------------------------------------------------------
//  variables not written out in V50
//-----------------------------------------------------------------------------
  fMcFlag        =  0;
  fInstLum       = -1;
//-----------------------------------------------------------------------------
//  set variables added in V51 to the default values
//-----------------------------------------------------------------------------
  fNTracks       = -1;
}

//______________________________________________________________________________
void TStnHeaderBlock::Streamer(TBuffer &R__b)
{
   // Stream an object of class TStnHeaderBlock: should never be called

  int nwi     = ((Int_t*  )&fInstLum   )-&fVersion;
  int nwf     = ((Float_t*)&fLastNumber)-&fInstLum;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if (R__v == 50) { 
      ReadV50(R__b);
    }
    else {
//-----------------------------------------------------------------------------
//  current version : 51
//-----------------------------------------------------------------------------
      R__b.ReadFastArray(&fVersion,nwi);
      R__b.ReadFastArray(&fInstLum,nwf);
    }
  } 
  else {
    R__b.WriteVersion(TStnHeaderBlock::IsA());
    R__b.WriteFastArray(&fVersion,nwi);
    R__b.WriteFastArray(&fInstLum,nwf);
  }
}

//------------------------------------------------------------------------------
//  init MET object from a UC standard Ntuple
//------------------------------------------------------------------------------
TStnHeaderBlock::TStnHeaderBlock() {
  fLastNumber    = -1;
  fEventNumber   = -1;
  fRunNumber     = -1;
  fLastRunNumber = -1;
}


//_____________________________________________________________________________
TStnHeaderBlock::~TStnHeaderBlock() {
}


//_____________________________________________________________________________
void TStnHeaderBlock::Clear(Option_t* opt) {

  fLinksInitialized = 0;
}


//_____________________________________________________________________________
void TStnHeaderBlock::Print(Option_t* opt) const {
  // don't print header for the same event 2 times in a raw
  static char f_last_opt[1000] = {0};

  TStnHeaderBlock* block = (TStnHeaderBlock*) this;

  if ((fLastNumber            == fEventNumber) && 
      (fLastRunNumber         == fRunNumber  ) &&
      (strcmp(opt,f_last_opt) == 0           )    ) return;

  block->fLastNumber    = fEventNumber;
  block->fLastRunNumber = fRunNumber;
  strncpy(f_last_opt,opt,1000);

  printf(" *** Run,Event: %6i,%-9i Rs: %5i : %s\n", 
	 fRunNumber,fEventNumber,fSectionNumber,opt);
}



