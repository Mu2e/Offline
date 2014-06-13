///////////////////////////////////////////////////////////////////////////////
// Mu2e Sim Particle - store only primary ones
///////////////////////////////////////////////////////////////////////////////
#ifndef STNTUPLE_TSimParticle
#define STNTUPLE_TSimParticle

#include "TMath.h"
#include "TObject.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

class TSimParticle : public TObject {
public:
  int             fParentID;
  int             fPdgCode;
  int             fCreationCode;
  int             fStartVolumeIndex;
  int             fTerminationCode;
  int             fEndVolumeIndex;
  int             fNStrawHits;

  float           fMomTargetEnd;
  float           fMomTrackerFront;		// entrance to ST

  TLorentzVector  fStartPos;
  TLorentzVector  fStartMom;

public:
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
				// ****** constructors and destructor
  TSimParticle();

  TSimParticle(Int_t ID, Int_t ParentID, Int_t PdgCode, 
	       int CreationCode, int TerminationCode,
	       int StartVolumeIndex, int EndVolumeIndex,
	       Float_t px, Float_t py, Float_t pz, Float_t e,
	       Float_t vx, Float_t vy, Float_t vz, Float_t t);

  virtual ~TSimParticle();
//-----------------------------------------------------------------------------
// init methods
//-----------------------------------------------------------------------------
  int Init(Int_t id, Int_t ParentID, Int_t PdgCode, 
	   int CreationCode, int TerminationCode,
	   int StartVolumeIndex, int EndVolumeIndex,
	   Float_t px, Float_t py, Float_t pz, Float_t e,
	   Float_t vx, Float_t vy, Float_t vz, Float_t t);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int    Number    () const { return GetUniqueID(); }
  int    NStrawHits() const { return fNStrawHits; }
//------------------------------------------------------------------------------
//  missing TParticle accessors and setters
//------------------------------------------------------------------------------
  void     SetNStrawHits(int N) { fNStrawHits = N; }

  void     SetMomTargetEnd   (double P) { fMomTargetEnd    = P; }
  void     SetMomTrackerFront(double P) { fMomTrackerFront = P; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void     Print(Option_t* opt = "") const;

  ClassDef(TSimParticle,1)
};

#endif


