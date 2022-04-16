
#include "Offline/Validation/inc/TValHist.hh"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

// ClassImp(TValHist)

//_____________________________________________________________________________
void TValHist::ClearB(Option_t* Opt) {
  fKsProb = 0.0;
  fFrProb = 0.0;
  fDiff = true;
  fStatus = fCantCompare;
  fFontScale = 1.0;
  fEmpty = false;
}
