//-----------------------------------------------------------------------------
// TSimParticle: STNTUPLE description of the generator-level MC particle
// It inherits from ROOT's TParticle and adds more functions to it
// particle number in the original list is saved in TObject::fUniqueID 
// Nov 25 2000 P.Murat (CDF/FNAL)
//----------------------------------------------------------------------------- 
#include <limits.h>	// UINT_MAX
#include "TDatabasePDG.h"
#include "Stntuple/obj/TSimParticle.hh"

ClassImp(TSimParticle)

//-----------------------------------------------------------------------------
void TSimParticle::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = ((int*  ) &fMomTargetEnd) - &fParentID;
  nwf = ((float*) &fStartPos    ) - &fMomTargetEnd ;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TObject::Streamer(R__b);

    R__b.ReadFastArray(&fParentID   ,nwi);
    R__b.ReadFastArray(&fMomTargetEnd,nwf);

    fStartPos.Streamer(R__b);
    fStartMom.Streamer(R__b);
  }
  else {
    R__b.WriteVersion(TSimParticle::IsA());
    TObject::Streamer(R__b);

    R__b.WriteFastArray(&fParentID   ,nwi);
    R__b.WriteFastArray(&fMomTargetEnd,nwf);

    fStartPos.Streamer(R__b);
    fStartMom.Streamer(R__b);
  }
}

//_____________________________________________________________________________
TSimParticle::TSimParticle() {
  SetUniqueID(UINT_MAX);
}

//_____________________________________________________________________________
TSimParticle::TSimParticle(Int_t ID, Int_t ParentID, Int_t PdgCode, 
			   int CreationCode, int TerminationCode,
			   int StartVolumeIndex, int EndVolumeIndex,
			   Float_t px, Float_t py, Float_t pz, Float_t e,
			   Float_t vx, Float_t vy, Float_t vz, Float_t t):
  TObject(),
  fStartPos(vx,vy,vz,t),
  fStartMom(px,py,pz,e)
  
{
  SetUniqueID(ID);
  fParentID         = ParentID;
  fPdgCode          = PdgCode;
  fCreationCode     = CreationCode;
  fTerminationCode  = TerminationCode;
  fStartVolumeIndex = StartVolumeIndex;
  fEndVolumeIndex   = EndVolumeIndex;
  fNStrawHits       = 0;
  fMomTargetEnd     = -1.;
  fMomTrackerFront  = -1.;
}

//-----------------------------------------------------------------------------
TSimParticle::~TSimParticle() {
}

//_____________________________________________________________________________
int  TSimParticle::Init(Int_t ID, Int_t ParentID, Int_t PdgCode, 
			int CreationCode, int TerminationCode,
			int StartVolumeIndex, int EndVolumeIndex,
			Float_t px, Float_t py, Float_t pz, Float_t e,
			Float_t vx, Float_t vy, Float_t vz, Float_t t) 
{
  SetUniqueID(ID);
  fParentID         = ParentID;
  fPdgCode          = PdgCode;
  fCreationCode     = CreationCode;
  fTerminationCode  = TerminationCode;
  fStartVolumeIndex = StartVolumeIndex;
  fEndVolumeIndex   = EndVolumeIndex;
  fNStrawHits       = 0;
  fMomTargetEnd     = -1.;
  fMomTrackerFront  = -1.;

  fStartPos.SetXYZT(vx,vy,vz,t);
  fStartMom.SetXYZT(px,py,pz,e);

  return 0;
}

//_____________________________________________________________________________
void TSimParticle::Print(Option_t* Opt) const {

  TString opt = Opt;
  if ((opt == "banner") || (opt == "")) {
				// print banner
    printf("   i name                   PDG  isthep  im1  im2  id1  id2      px");
    printf("      py       pz       e        vx       vy        vz       t\n");
  }

  TDatabasePDG* db = TDatabasePDG::Instance();

  TParticlePDG* pdg = db->GetParticle(fPdgCode);

  if ((opt == "data") || (opt == "")) {
    printf("%4i",Number());
    if (pdg) printf(" %-19s",pdg->GetName());
    else          printf(" %-19s","*** unknown ***");
    printf("%7i"  ,fPdgCode);
    printf("%9.3f",fStartMom.Px());
    printf("%9.3f",fStartMom.Py());
    printf("%9.3f",fStartMom.Pz());
    printf("%9.3f",fStartMom.Energy());
    printf("%9.3f",fStartPos.X());
    printf("%9.3f",fStartPos.Y());
    printf("%9.3f",fStartPos.X());
    printf("%9.3f",fStartPos.T());
    printf("\n");
  }
}


