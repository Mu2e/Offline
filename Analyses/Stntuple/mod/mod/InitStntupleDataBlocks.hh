#ifndef Stntuple_mod_InitStntupleDataBlocks_hh
#define Stntuple_mod_InitStntupleDataBlocks_hh

#include "TObject.h" 

#include "Stntuple/obj/AbsEvent.hh"

class TStnDataBlock;

Int_t StntupleInitMu2eStrawDataBlock(TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eCalDataBlock  (TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eClusterBlock  (TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eHeaderBlock   (TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eTrackBlock    (TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eGenpBlock     (TStnDataBlock* blk, AbsEvent* evt, int mode);
Int_t StntupleInitMu2eSimpBlock     (TStnDataBlock* blk, AbsEvent* evt, int mode);

					// block-to-block link resolution

Int_t StntupleInitMu2eHeaderBlockLinks (TStnDataBlock* Block, AbsEvent* AnEvent, int Mode);
Int_t StntupleInitMu2eTrackBlockLinks  (TStnDataBlock* blk  , AbsEvent* evt    , int mode);
Int_t StntupleInitMu2eClusterBlockLinks(TStnDataBlock* blk  , AbsEvent* evt    , int mode);


#endif
