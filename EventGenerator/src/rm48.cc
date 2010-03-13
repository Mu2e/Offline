//
//  Adapter to make RandFlat look like the cernlib rn48.
//
//  $Id: rm48.cc,v 1.1 2010/03/13 00:14:44 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2010/03/13 00:14:44 $
//
//  Original author Rob Kutschke.
//
//  See rm48.hh for details.

#include "EventGenerator/src/rm48.hh"
#include "CLHEP/Random/RandFlat.h"

static CLHEP::RandFlat* distribution(0);

void rm48_ ( double *v, int *n){
  
  //if (!distribution){
  // throw here.
  //}

  for ( int i=0; i<*n; ++i){
    *(v++) = distribution->shoot();
  }

}

void setRm48Distribution( CLHEP::RandFlat& dist){
  distribution = &dist;
}
