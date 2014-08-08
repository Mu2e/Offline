#ifndef STNTUPLE_TVdetDataBlock
#define STNTUPLE_TVdetDataBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TVdetHitData.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

class TVdetDataBlock: public TStnDataBlock {
  friend Int_t StntupleInitMu2eVirtualDataBlock(TStnDataBlock* block, 
					      AbsEvent*      event, 
					      int            mode);
public:
  Int_t          fNHits;	        // number of hits in the virtual detectors
  TClonesArray*  fListOfHits;		// list of hits
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TVdetDataBlock();
  virtual ~TVdetDataBlock();
					// ****** accessors
  Int_t          NHits    () { return fNHits; }
  TVdetHitData* Hit (int i) { return (TVdetHitData*) fListOfHits->UncheckedAt(i); }
  
  TClonesArray* GetListOfHits () { return fListOfHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        //Create hit, increse number of hits

  TVdetHitData* NewHit() { return new ((*fListOfHits)[fNHits++]) TVdetHitData(); } 
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TVdetDataBlock,1)	// virtual data block
};


#endif


