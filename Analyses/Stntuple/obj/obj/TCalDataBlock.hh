#ifndef TCalDataBlock_hh
#define TCalDataBlock_hh

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnCrystal.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "TCalHitData.hh"

class TCalDataBlock : public TStnDataBlock {
  friend Int_t StntupleInitMu2eCalDataBlock(TStnDataBlock*, AbsEvent*, int);
public:
  // this is version v2
  int            fNHits;		// number of hit crystals
  int            fNDisks;               // 
  int            fNCrystals  [4];	// 
  float          fRMin       [4];	// as a temporary measure, store 
  float          fRMax       [4];
  float          fZ0         [4];       // 
  float          fCrystalSize;
  float          fMinFraction;		// min fr of the included crystal area
  float          fWrapperThickness;     //
  float          fShellThickness;        //

  TClonesArray*  fListOfCalHitData;	// list of crystal hit data 
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TCalDataBlock();
  virtual ~TCalDataBlock();
					// ****** accessors

  Int_t         NHits           () { return fNHits; }

  TCalHitData*  CalHitData(int I) { 
    return (TCalHitData*) fListOfCalHitData->UncheckedAt(I);
  }

  TClonesArray* GetListOfCalHitData () { return fListOfCalHitData ; } 
  int           NDisks() { return fNDisks; }

  float         RMin(int I) { return fRMin[I]; }
  float         RMax(int I) { return fRMax[I]; }
  float         Z0  (int I) { return fZ0  [I]; }

  int           NCrystals(int I) { return fNCrystals[I]; }

  float         CrystalSize() { return fCrystalSize; }
  float         MinFraction() { return fMinFraction; }
  float         WrapperThickness() { return fWrapperThickness; }
  float         ShellThickness  () { return fShellThickness;   }

  int           DiskNumber   (int ID);
  int           CrystalRadius(int ID);
//-----------------------------------------------------------------------------
// these routines do loops - be careful using them if you need performance
//-----------------------------------------------------------------------------
//   TCalHitData*  GetTowerByKey(Int_t Key) {
//     return Tower(TCalHitData::IEta(Key), TCalHitData::IPhi(Key));
//   }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCalHitData*  NewCalHitData() { 
    return new ((*fListOfCalHitData)[fNHits++]) TCalHitData();
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void        ReadV1(TBuffer& R__b);
  void        ReadV2(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void        Clear(Option_t* opt="");
  void        Print(Option_t* opt="") const;

  ClassDef(TCalDataBlock,2)
};

#endif
