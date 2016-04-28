///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "CalPatRec/inc/Ref.hh"

#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/KalFitHack.hh"

//-----------------------------------------------------------------------------
Ref::Ref():TNamed() {
  fHelixFit = NULL;
  fSeedFit  = NULL;
  fKalFit   = NULL;
}

//-----------------------------------------------------------------------------
Ref::Ref(const char* Name, const char* Title, 
	 mu2e::HelixFitHack* HelixFit, mu2e::KalFitHack* SeedFit, mu2e::KalFitHack* KalFit):
  TNamed(Name,Title) 
{
  fHelixFit = HelixFit;
  fSeedFit  = SeedFit;
  fKalFit   = KalFit;
}
//-----------------------------------------------------------------------------
Ref::Ref(const char* Name, const char* Title, 
	 mu2e::HelixFitHack* HelixFit):
  TNamed(Name,Title) 
{
  fHelixFit = HelixFit;
  fSeedFit  = 0;
  fKalFit   = 0;
}
//-----------------------------------------------------------------------------
Ref::Ref(const char* Name, const char* Title, 
	 mu2e::KalFitHack*KalFit ):
  TNamed(Name,Title) 
{
  fHelixFit = 0;
  fSeedFit  = 0;
  fKalFit   = KalFit;
}

//-----------------------------------------------------------------------------
Ref::~Ref() {
}

//-----------------------------------------------------------------------------
void Ref::PlotXY(int ISet) {
  fHelixFit->plotXY(ISet);
}

//-----------------------------------------------------------------------------
void Ref::PlotZPhi(int ISet) {
  fHelixFit->plotZPhi(ISet);
}

