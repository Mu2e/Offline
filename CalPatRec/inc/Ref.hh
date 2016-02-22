///////////////////////////////////////////////////////////////////////////////
// proxy
///////////////////////////////////////////////////////////////////////////////

#ifndef CalPatRec_Ref
#define CalPatRec_Ref

#include "TNamed.h"

namespace mu2e {
  class HelixFitHack;
  class KalFitHack;
};

class Ref : public TNamed {

public:
  mu2e::HelixFitHack*  fHelixFit;
  mu2e::KalFitHack*    fSeedFit;
  mu2e::KalFitHack*    fKalFit;

  Ref();
  Ref(const char* Name, const char* Title, 
      mu2e::HelixFitHack* HelixFit, mu2e::KalFitHack* SeedFit, mu2e::KalFitHack* KalFit);

  ~Ref();

  void PlotXY  (int ISet = -1);
  void PlotZPhi(int ISet = -1);

  //  ClassDef(Ref,0)
};

#endif
