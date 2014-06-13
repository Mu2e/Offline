#ifndef Stntuple_mod_StntupleInitDataBlocks_hh
#define Stntuple_mod_StntupleInitDataBlocks_hh


#include "TObject.h" 

#include "Stntuple/obj/AbsEvent.hh"


class TAnaTrackBlock;
class TStnDataBlock;

Int_t InitMu2eTrackBlock  (TStnDataBlock* Block, AbsEvent* Evt, Int_t Mode);
Int_t InitMu2eClusterBlock(TStnDataBlock* Block, AbsEvent* Evt, Int_t Mode);

#endif
