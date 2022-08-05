//
// Shape description for TrkExt
//
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtShape_HH
#define TrkExtShape_HH

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TrkExtShape {

    public:
      TrkExtShape(double boundaryLimit = 0.001) ;
      ~TrkExtShape() { }

      virtual bool contains (CLHEP::Hep3Vector& p)  =0 ;
      virtual void initialize (void)  =0 ;
      CLHEP::Hep3Vector  intersection (const CLHEP::Hep3Vector & x1, const CLHEP::Hep3Vector & x2) ;

    protected:

      double _limit;

  };

} // end namespace mu2e


#endif
