
#include "Validation/inc/TValPar.hh"

//ClassImp(TValPar)

//_____________________________________________________________________________
void TValPar::SetIndependent(Int_t x) {
  fIndep = x; 
  // if the loose and tight were set by user, don't overwrite
  if(!fChanged) {
    if(fIndep==0) {
      fLoose = 0.99; 
      fTight = 0.999;
    } else {
      fLoose = 0.001; 
      fTight = 0.01;
    }
  }
}

//_____________________________________________________________________________
void TValPar::Clear(Option_t* Opt) {
  fMode = 0;
  fScale1 = 1.0;
  fScale2 = 1.0;
  fUnder = 0;
  fOver = 0;
  fLoose = 0.99;
  fTight = 0.999;
  fChanged = false;
}

