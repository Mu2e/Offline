#ifndef STNTUPLE_TGenpBlock
#define STNTUPLE_TGenpBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE GENP (generator-level) particle block
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 25 2000
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale 
//  TGenpBlock will be written in split mode, however the TClonesArray 
//  itself may not be split
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "TArrayI.h"
#include "TClonesArray.h"

#include "Stntuple/base/TStnArrayI.hh"
#include "TStnDataBlock.hh"
#include "TGenParticle.hh"

class TStnEvent;

class TGenpBlock: public TStnDataBlock {
  friend Int_t StntupleInitMu2eGenpBlock  (TStnDataBlock*, TStnEvent*, int);
protected:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNParticles;		// total # of particles
  TClonesArray*  fListOfParticles;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TGenpBlock();
  virtual ~TGenpBlock();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NParticles        () { return fNParticles; }

					// `i'-th particle in the global list
  TGenParticle*   Particle(int i) { 
    return (TGenParticle*) fListOfParticles->UncheckedAt(i); 
  }
//-----------------------------------------------------------------------------
//  modifiers
//-----------------------------------------------------------------------------
					// currently: add LAST particle to the
					// LAST interaction

  TGenParticle*  NewParticle(Int_t PdgCode, Int_t GeneratorID, 
			     int   m1, Int_t m2, Int_t d1, Int_t d2, 
			     float px, Float_t py, Float_t pz, Float_t e,
			     float vx, Float_t vy, Float_t vz, Float_t t,
			     float ProperTime);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
//  I/O and schema evolution
//-----------------------------------------------------------------------------
  ClassDef(TGenpBlock,1)		// GENP block: output of MC event generators
};

#endif


