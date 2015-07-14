#ifndef STNTUPLE_TStnPidBlock
#define STNTUPLE_TStnPidBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE track block- prototype for MU2E
//  cloned from CDF TStnPidBlock
//
//  Author:    Pavel Murat (CDF/FNAL)
//  Date:      Mar 07 2013
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale TStnPidBlock
//  will be written in split mode, however the TClonesArray itself may not
//  be split (if one needs tracks, he normally needs all the tracks, not only
//  track Pt)
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TStnPid.hh"

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"

class TStnEvent;
class TStnPid;

class TStnPidBlock: public TStnDataBlock {

  friend Int_t StntupleInitMu2ePidBlock     (TStnDataBlock*, AbsEvent*, int);
  friend Int_t StntupleInitMu2ePidBlockLinks(TStnDataBlock*, AbsEvent*, int);
public:
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
  Int_t          fNTracks;
  TClonesArray*  fListOfTrackPid;
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnPidBlock();
  virtual ~TStnPidBlock();
					// ****** accessors

  Int_t           NTracks       () { return fNTracks; }
  TClonesArray*   ListOfTrackPid() { return fListOfTrackPid; }

  TStnPid*      Pid(int i) { 
    return (TStnPid*) fListOfTrackPid->UncheckedAt(i); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TStnPid* NewPid() {
    TStnPid* t = new ((*fListOfTrackPid)[fNTracks]) TStnPid(fNTracks);
    fNTracks++;
    return t;
  }
//-----------------------------------------------------------------------------
// overloaded functions of  TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* option = "") const;

  ClassDef(TStnPidBlock,1)
};

#endif
