#ifndef RM48_HH
#define RM48_HH

//
//  Adapter to make RandFlat look like the cernlib rn48.
//
//  $Id: rm48.hh,v 1.1 2010/03/15 21:27:24 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/03/15 21:27:24 $
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

// Called by the Daya Bay comsic code.
extern "C" {
  void rm48_ ( double *v, int* n);
}

namespace mu2e {

  // Called by CosmicDYB.
  void setRm48Distribution( CLHEP::RandFlat& dist);

}

#endif
