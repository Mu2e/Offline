//-----------------------------------------------------------------------------
//  2014-01-26 P.Murat: Mu2e calorimeter hit data
//-----------------------------------------------------------------------------
#ifndef TCalHitData_hh
#define TCalHitData_hh

#include "TObject.h"

class TCalHitData : public TObject {
public:

protected: 
  int        fID;         // hit ID, cods disk,  x1, x2
  int        fNChannels;  // number of readout channels, 1 or 2 (kludge)
  float      fTime; 
  float      fEnergy;
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCalHitData();
  virtual ~TCalHitData();
					// ****** initialization
  //  static void InitStaticVariables();
//-----------------------------------------------------------------------------
// static methods
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     ID    () { return fID;     }
  int     NChannels() { return fNChannels; }
  float   Time  () { return fTime;   }
  float   Energy() { return fEnergy; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void Set(int ID, int NChannels, float Time,  float Energy) {
    fID = ID; fNChannels = NChannels; fTime = Time; fEnergy = Energy;
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;

  ClassDef (TCalHitData,1)
};

#endif
