//
// Detector description for TrkExt
//
//  $Id: TrkExtDetectors.hh,v 1.1 2012/08/04 00:22:09 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:22:09 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtDetectors_HH
#define TrkExtDetectors_HH

#include <string>
#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtShape.hh"
#include "TrkExt/inc/TrkExtToyDS.hh"
#include "TrkExt/inc/TrkExtProtonAbsorber.hh"
#include "TrkExt/inc/TrkExtStoppingTarget.hh"

namespace mu2e {

  namespace TrkExtDetectorList {
    enum Enum {
      Undefined = -1,
      ToyDS = 0,
      ProtonAbsorber = 1,
      StoppingTarget = 2 
    };
  };

  class TrkExtDetectors : public TrkExtShape
  {

  public:
    TrkExtDetectors() ;
    ~TrkExtDetectors() { }

    void initialize () ;
    bool contains (CLHEP::Hep3Vector& p) ;
    TrkExtDetectorList::Enum volumeId(CLHEP::Hep3Vector &xx) ;
    double limit() { return _limit; }
    double mostProbableEnergyLoss (CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid = TrkExtDetectorList::Undefined) ;
    double meanEnergyLoss (CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid = TrkExtDetectorList::Undefined) ;
    double scatteringAngle (CLHEP::Hep3Vector& p, double ds, TrkExtDetectorList::Enum volid) ;
    CLHEP::Hep3Vector intersection (CLHEP::Hep3Vector & x1, CLHEP::Hep3Vector & x2) ;


  private:

    TrkExtToyDS _ds;
    TrkExtProtonAbsorber _pa;
    TrkExtStoppingTarget _st;

  };



} // end namespace mu2e


#endif
