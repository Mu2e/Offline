#ifndef STNTUPLE_TStnTrackBlock
#define STNTUPLE_TStnTrackBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE track block- prototype for MU2E
//  cloned from CDF TStnTrackBlock
//
//  Author:    Pavel Murat (CDF/FNAL)
//  Date:      Mar 07 2013
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale TStnTrackBlock
//  will be written in split mode, however the TClonesArray itself may not
//  be split (if one needs tracks, he normally needs all the tracks, not only
//  track Pt)
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TStnTrack.hh"

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"

class TStnEvent;
class TStnTrack;

class TStnTrackBlock: public TStnDataBlock {

  friend Int_t StntupleInitMu2eTrackBlock     (TStnDataBlock*, AbsEvent*, int);
  friend Int_t StntupleInitMu2eTrackBlockLinks(TStnDataBlock*, AbsEvent*, int);
public:
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
  Int_t          fNTracks;
  TClonesArray*  fListOfTracks;
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnTrackBlock();
  virtual ~TStnTrackBlock();
					// ****** accessors

  Int_t           NTracks     () { return fNTracks; }
  TClonesArray*   ListOfTracks() { return fListOfTracks; }

  TStnTrack*      Track(int i) { 
    return (TStnTrack*) fListOfTracks->UncheckedAt(i); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TStnTrack* NewTrack() {
    TStnTrack* t = new ((*fListOfTracks)[fNTracks]) TStnTrack(fNTracks);
    fNTracks++;
    return t;
  }
//-----------------------------------------------------------------------------
// overloaded functions of  TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* option = "") const;

  ClassDef(TStnTrackBlock,1)
};

#endif
