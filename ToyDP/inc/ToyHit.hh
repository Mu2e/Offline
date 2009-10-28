#ifndef ToyDP_ToyHit_hh
#define ToyDP_ToyHit_hh

//
// Define a hit in the CTracker.
// The CTracker is a toy detector that exists to illustrate
// the framework.
//
// $Id: ToyHit.hh,v 1.2 2009/10/28 13:41:16 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/28 13:41:16 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  struct ToyHit {
    ToyHit();
    ToyHit( double x, double y, double z, double ph);
    ToyHit( CLHEP::Hep3Vector const& pos, double ph);

    virtual ~ToyHit();

    CLHEP::Hep3Vector _position;
    double _pulseheight;
  };

} // end namespace mu2e

#endif
