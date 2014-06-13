//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//	Class FillTreeModule: histogrammming for Geant3 module
//
// Environment:
//	Software developed for the CDF at FNAL.
//
// Aug 31 2000 P.Murat: module to fill Stntuple tree
//
// Copyright Information:
//	Copyright (C) 1997		Fermilab
//
//  revision history:
//  -----------------
//------------------------------------------------------------------------

#ifndef Stntuple_mod_FillStntuple
#define Stntuple_mod_FillStntuple

#include "Stntuple/mod/StntupleModule.hh"

namespace mu2e {

class FillStntuple : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
  Int_t             fLastRun;		// last run with events
//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  FillStntuple(fhicl::ParameterSet const& PSet);

  ~FillStntuple();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  Int_t     ProcessNewRun      (art::Run* aRun);
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run&   r);
  virtual bool endRun  (art::Run&   r);
  virtual bool filter  (AbsEvent&   e);

  //  ClassDef(FillStntuple,0)
};

}
#endif





