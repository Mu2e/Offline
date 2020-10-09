#ifndef Mu2eUtilities_rm48_hh
#define Mu2eUtilities_rm48_hh

//
//  Adapter to make CLHEP::RandFlat look like the cernlib rn48.
//
//
//  Original author Rob Kutschke.
//
//  The Daya Bay cosmic code, written in fortran, calls the
//  cernlib random number generator rm48, this is a flat
//  distribution on ]0,1[.
//  The code in this file and rm48.cc is an adapter to
//  allow calls to rm48 to call through to an instance of
//  CLHEP::RandFlat.  The existence of the instance
//  is managed in CosmicDYB.cc.
//
//  There is a potential problem.  If other codes also want
//  to use rm48, this implementation will break the rule
//  that each module gets its own random engine.

#include "CLHEP/Random/RandFlat.h"

namespace CLHEP { class RandFlat; }

// Called by the Daya Bay comsic code.
extern "C" {
  void rm48_ ( double *v, int* n);
}

namespace mu2e {

  // Called by CosmicDYB.
  void setRm48Distribution( CLHEP::RandFlat& dist);

}

#endif /* Mu2eUtilities_rm48_hh */
