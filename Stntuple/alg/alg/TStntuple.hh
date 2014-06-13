#ifndef STNTUPLE_ALG_TStntuple_hh
#define STNTUPLE_ALG_TStntuple_hh

#include <math.h>
#include "TClonesArray.h"
#include "TVector2.h"
#include "TVector3.h"

class TStnTrack;
class TStnCluster;
class TStnElectron;

class TStntuple: public TObject {
protected:
  static TStntuple*  fgInstance;

  static Int_t             fgRunNumber;
  static Float_t           fgEventVertex;

  class  Cleaner {
  public:
    Cleaner ();
    ~Cleaner();
  };
  friend class Cleaner;
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
                                // don't use constructor! use Instance()
                                // instead (if ever necessary...)
  TStntuple();
  virtual ~TStntuple();

  static TStntuple*  Instance();

  static Int_t     Init(Int_t RunNumber);

  static double DioWeightAl(double P);
  static double DioWeightTi(double P);
//-----------------------------------------------------------------------------
// print routines - sometimes it is not possible to do it from a single block
//-----------------------------------------------------------------------------
  // static int  PrintElectron(TStnElectron*       Ele,
  // 			    TStnElectronBlock*  ElectronBlock,
  // 			    TCalDataBlock*      CalDataBlock);
  ClassDef(TStntuple,0)
};

#endif

