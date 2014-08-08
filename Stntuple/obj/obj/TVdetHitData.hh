//-----------------------------------------------------------------------------
//  2014-07-15 TVdetHitData
//-----------------------------------------------------------------------------
#ifndef TVdetHitData_hh
#define TVdetHitData_hh

#include <math.h>
#include "TMath.h"
#include "TObject.h"

class TVdetHitData : public TObject {
protected:
  int     fIndex;
  int     fPdgCode;
  int     fGeneratorCode;
  float   fTime;
  float   fMass;
  float   fEnergyKin;
  float   fEnergy;
  float   fMcMomentum;			// MC particle momentum
  float   fMcMomentumX;			// MC particle momentum X - component
  float   fMcMomentumY;			// MC particle momentum Y - component
  float   fMcMomentumZ;			// MC particle momentum Z - component
  float   fMcPositionX;			// MC particle position X - component
  float   fMcPositionY;			// MC particle position Y - component
  float   fMcPositionZ;			// MC particle position Z - component
  
public:
                                        // ****** constructors and destructors
  TVdetHitData();
  virtual ~TVdetHitData();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Index () { return fIndex; }
  float   Time  () { return fTime; }
  float   Mass  () { return fMass; }

  float   EnergyKin () { return fEnergyKin; }
  float   Energy    () { return fEnergy; }

  int     PdgCode      () { return fPdgCode;       }
  int     GeneratorCode() { return fGeneratorCode; }

  float   McMomentum   () { return fMcMomentum;    }
  float   McMomentumX   () { return fMcMomentumX;    }
  float   McMomentumY   () { return fMcMomentumY;    }
  float   McMomentumZ   () { return fMcMomentumZ;    }

  float   McPositionX   () { return fMcPositionX;    }
  float   McPositionY   () { return fMcPositionY;    }
  float   McPositionZ   () { return fMcPositionZ;    }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    Set(int Index, float Time,  float Mass, float EnergyKin, 
	      float Energy,
	      int PdgID, int GenCode,
	      float McMomentum, 
	      float McMomentumX, float McMomentumY, float McMomentumZ,
	      float McPositionX, float McPositionY, float McPositionZ);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);

  ClassDef (TVdetHitData,2)
};

#endif
