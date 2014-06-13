#ifndef STNTUPLE_TStrawDataBlock
#define STNTUPLE_TStrawDataBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStrawHitData.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

class TStrawDataBlock: public TStnDataBlock {
  friend Int_t StntupleInitMu2eStrawDataBlock(TStnDataBlock* block, 
					      AbsEvent*      event, 
					      int            mode);
public:
  Int_t          fNHits;	        // number of hits in the straw tracker
  TClonesArray*  fListOfHits;		// list of hits
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStrawDataBlock();
  virtual ~TStrawDataBlock();
					// ****** accessors
  Int_t          NHits    () { return fNHits; }
  TStrawHitData* Hit (int i) { return (TStrawHitData*) fListOfHits->UncheckedAt(i); }
  
  TClonesArray* GetListOfHits () { return fListOfHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        //Create hit, increse number of hits

  TStrawHitData* NewHit() { return new ((*fListOfHits)[fNHits++]) TStrawHitData(); } 
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TStrawDataBlock,1)	// straw tracker data block
};


#endif


