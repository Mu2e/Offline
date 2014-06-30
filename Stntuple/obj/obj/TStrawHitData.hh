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
  int     fPdgCode;
  int     fMotherPdgCode;
  int     fGeneratorCode;
  int     fSimID;
  float   fTime;
  float   fDt;
  float   fEnergy;
  float   fMcMomentum;			// MC particle momentum
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

  int     PdgCode      () { return fPdgCode;       }
  int     MotherPdgCode() { return fMotherPdgCode; }
  int     GeneratorCode() { return fGeneratorCode; }
  int     SimID        () { return fSimID;         }

  float   McMomentum   () { return fMcMomentum;    }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    Set(int Index, float Time, float Dt, float EnergyDep,
	      int PdgID, int MotherPdgID, int GenCode, int SimID, 
	      float McMomentum);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);

  ClassDef (TStrawHitData,2)
};

#endif
