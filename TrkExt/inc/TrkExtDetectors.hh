//
// Detector description for TrkExt
//
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtDetectors_HH
#define TrkExtDetectors_HH

#include <string>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/TrkExt/inc/TrkExtShape.hh"
#include "Offline/TrkExt/inc/TrkExtToyDS.hh"
#include "Offline/TrkExt/inc/TrkExtProtonAbsorber.hh"
#include "Offline/TrkExt/inc/TrkExtStoppingTarget.hh"

namespace mu2e {

  namespace TrkExtDetectorList {
    enum Enum {
      Undefined = -1,
      ToyDS = 0,
      ProtonAbsorber = 1,
      StoppingTarget = 2
    };
  }

  class TrkExtDetectors : public TrkExtShape
  {

    public:
      TrkExtDetectors() ;
      ~TrkExtDetectors() { }

      void initialize () ;
      bool contains (CLHEP::Hep3Vector& p) ;
      TrkExtDetectorList::Enum volumeId(CLHEP::Hep3Vector &xx) ;
      double limit() { return _limit; }
      double mostProbableEnergyLoss (const CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid = TrkExtDetectorList::Undefined) ;
      double meanEnergyLoss (const CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid = TrkExtDetectorList::Undefined) ;
      double scatteringAngle (const CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid) ;
      CLHEP::Hep3Vector intersection (const CLHEP::Hep3Vector & x1, const CLHEP::Hep3Vector & x2) ;


    private:

      TrkExtToyDS _ds;
      TrkExtProtonAbsorber _pa;
      TrkExtStoppingTarget _st;

  };



} // end namespace mu2e


#endif
