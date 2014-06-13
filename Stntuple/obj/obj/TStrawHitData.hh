//-----------------------------------------------------------------------------
//  2014-01-26 P.Murat - TStrawHitData
//-----------------------------------------------------------------------------
#ifndef TStrawHitData_hh
#define TStrawHitData_hh

#include <math.h>
#include "TMath.h"
#include "TObject.h"

class TStrawHitData : public TObject {
protected:
  int     fIndex;
  float   fTime;
  float   fDt;
  float   fEnergy;
public:
                                        // ****** constructors and destructors
  TStrawHitData();
  virtual ~TStrawHitData();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Index () { return fIndex; }
  float   Time  () { return fTime; }
  float   Dt    () { return fDt; }
  float   Energy() { return fEnergy; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    Set(int Index, float Time, float Dt, float EnergyDep) {
    fIndex = Index; fTime = Time; fDt = Dt; fEnergy = EnergyDep;
  }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;

  ClassDef (TStrawHitData,1)
};

#endif
