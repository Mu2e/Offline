
#include "Validation/inc/TValHist.hh"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"

//ClassImp(TValHist)

//_____________________________________________________________________________
void TValHist::ClearB(Option_t* Opt) {
  fKsProb = 0.0;
  fFrProb = 0.0;
  fDiff = true;
  fStatus = 10;
  fFontScale = 1.0;
}

