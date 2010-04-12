#ifndef ToyDP_HoughCircle_hh
#define ToyDP_HoughCircle_hh

//
// A candidate hough circle - i.e., a geometric cicle and some information
// about hough peak quality (eventually)
//
// $Id: HoughCircle.hh,v 1.1 2010/04/12 18:12:28 shanahan Exp $
// $Author: shanahan $
// $Date: 2010/04/12 18:12:28 $
//
// Original author Peter Shanahan
//

#include "CLHEP/Vector/TwoVector.h"

namespace mu2e {

  // first version.  We'll undoubtably want to add more later
  class HoughCircle {

    public:
      HoughCircle();
      HoughCircle( double x, double y, double radius, uint32_t nstraw);
      HoughCircle( CLHEP::Hep2Vector const& cent, double radius, uint32_t nstraw);
      virtual ~HoughCircle();

      const CLHEP::Hep2Vector& Center() const {return _center;};
      const double Radius() const {return _radius;};
      const uint32_t NStraws() const {return _nStraws;};

    private: 
      CLHEP::Hep2Vector _center;
      double _radius;
      uint32_t _nStraws;
  };

} // end namespace mu2e

#endif // ToyDP_HoughCircle_hh
