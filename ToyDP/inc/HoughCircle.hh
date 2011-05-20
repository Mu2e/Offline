#ifndef ToyDP_HoughCircle_hh
#define ToyDP_HoughCircle_hh

//
// A candidate hough circle - i.e., a geometric cicle and some information
// about hough peak quality (eventually)
//
// $Id: HoughCircle.hh,v 1.5 2011/05/20 20:18:24 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:18:24 $
//
// Original author Peter Shanahan
//

#include "CLHEP/Vector/TwoVector.h"

namespace mu2e {

  // first version.  We'll undoubtably want to add more later
  class HoughCircle {

    public:
      HoughCircle();
      HoughCircle( double x, double y, double radius, unsigned nstraw);
      HoughCircle( CLHEP::Hep2Vector const& cent, double radius, unsigned nstraw);
      virtual ~HoughCircle();

      const CLHEP::Hep2Vector& Center() const {return _center;}
      const double Radius() const {return _radius;}
      const unsigned NStraws() const {return _nStraws;}

    private:
      CLHEP::Hep2Vector _center;
      double _radius;
      unsigned _nStraws;
  };

} // end namespace mu2e

#endif /* ToyDP_HoughCircle_hh */
