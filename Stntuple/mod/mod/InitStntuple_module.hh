//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//	Class InitStntuple: initializations for STNTUPLE maker
//
// Environment:
//	Software developed for the CDF at FNAL.
//
// Nov 23 2000 P.Murat
//
//  revision history:
//  -----------------
//------------------------------------------------------------------------

#ifndef Stntuple_mod_InitStntuple
#define Stntuple_mod_InitStntuple

class TObjArray;

#include "Stntuple/mod/StntupleModule.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

namespace mu2e {
class InitStntuple : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:

  Int_t                    fLastRun;
  Float_t                  fSumInstLum; //! avg inst lum for evaluating
  Int_t                    fnLum;       //! exe speed
  Float_t                  fCpuSpeed;   //! MHz of CPU

//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  InitStntuple(fhicl::ParameterSet const& Pset);

  ~InitStntuple();
				        // ****** accessors
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run&   r);
  virtual bool filter  (AbsEvent&   e);
  virtual bool endRun (art::Run&   r);
  virtual void beginJob();
  virtual void endJob  ();
					// ****** functions of the module

  Int_t     ProcessNewRun      (art::Run* ARun);
  // Int_t     InitTriggerTable   ();
  // Int_t     InitRunSummary     ();
					// ****** setters
  //  ClassDef(InitStntuple,0)
};
} // end namespace mu2e

#endif
