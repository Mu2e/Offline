//-----------------------------------------------------------------------------
// TGenParticle: STNTUPLE description of the generator-level MC particle
// It inherits from ROOT's TParticle and adds more functions to it
// particle number in the original list is saved in TObject::fUniqueID 
// Nov 25 2000 P.Murat (CDF/FNAL)
//----------------------------------------------------------------------------- 
#include <limits.h>	// UINT_MAX
#include "TDatabasePDG.h"
#include "Stntuple/obj/TGenParticle.hh"

ClassImp(TGenParticle)

//_____________________________________________________________________________
TGenParticle::TGenParticle() {
  SetUniqueID(UINT_MAX);
}


// //_____________________________________________________________________________
// TGenParticle::TGenParticle(GENP_PARTICLE* mc) {
// }

//_____________________________________________________________________________
TGenParticle::TGenParticle(Int_t ip, Int_t PdgCode, Int_t GeneratorID, 
			   Int_t jm1, Int_t jm2, Int_t jd1, Int_t jd2, 
			   Float_t px, Float_t py, Float_t pz, Float_t e,
			   Float_t vx, Float_t vy, Float_t vz, Float_t t,
			   float ProperTime) :
  TParticle(PdgCode,GeneratorID,jm1,jm2,jd1,jd2,px,py,pz,e,vx,vy,vz,t)
{
  SetUniqueID(ip);
  fProperTime = ProperTime;
}

//_____________________________________________________________________________
Int_t TGenParticle::Init(Int_t ip, Int_t idhep, Int_t istdhep,
			 Int_t jm1, Int_t jm2, Int_t jd1, Int_t jd2, 
			 Float_t px, Float_t py, Float_t pz, Float_t e,
			 Float_t vx, Float_t vy, Float_t vz, Float_t t,
			 float   ProperTime) 
{
  SetUniqueID(ip);
  fPdgCode     = idhep;
  fStatusCode  = istdhep;
  fMother[0]   = jm1;
  fMother[1]   = jm2;
  fDaughter[0] = jd1;
  fDaughter[1] = jd2;
  fWeight      = 1.;

  fParticlePDG = TDatabasePDG::Instance()->GetParticle(idhep);

  fPx          = px;
  fPy          = py;
  fPz          = pz;
  fE           = e;
  fVx          = vx;
  fVy          = vy;
  fVz          = vz;
  fVt          = t;
  fPolarTheta  = 0;
  fPolarPhi    = 0;
  
  if (fParticlePDG) {
    fCalcMass    = fParticlePDG->Mass();
  } 
  else {
    Double_t a2 = fE*fE -fPx*fPx -fPy*fPy -fPz*fPz;
    if (a2 >= 0) fCalcMass =  TMath::Sqrt(a2);
    else         fCalcMass = -TMath::Sqrt(-a2);
  }

  fProperTime = ProperTime;
  return 0;
}

//_____________________________________________________________________________
TGenParticle::~TGenParticle() {
}

//_____________________________________________________________________________
void TGenParticle::Print(Option_t* Opt) const {

  TString opt = Opt;
  if ((opt == "banner") || (opt == "")) {
				// print banner
    printf("   i name                   PDG  isthep  im1  im2  id1  id2      px");
    printf("      py       pz       e        vx       vy        vz       t\n");
  }

  //  TDatabasePDG* db = TDatabasePDG::Instance();

  TParticlePDG* pdg = ((TParticle*)this)->GetPDG();

  if ((opt == "data") || (opt == "")) {
    printf("%4i",Number());
    if (pdg) printf(" %-19s",pdg->GetName());
    else          printf(" %-19s","*** unknown ***");
    printf("%7i"   ,GetPdgCode());
    printf("%6i"   ,GetStatusCode());
    printf("  ");
    printf("%5i"   ,GetMother(0));
    printf("%5i"   ,GetMother(1));
    printf("%5i"   ,GetDaughter(0));
    printf("%5i"   ,GetDaughter(1));
    printf("%11.3f",Px());
    printf("%11.3f",Py());
    printf("%11.3f",Pz());
    printf("%11.3f",Energy());
    printf("%11.3f",Vx());
    printf("%11.3f",Vy());
    printf("%11.3f",Vz());
    printf("%11.3f",T ());
    printf("\n");
  }
}


