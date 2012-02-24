// Geometry of the production target. This also defines the proton beam direction.
// 
// Andrei Gaponenko, 2011

#ifndef PRODUCTIONTARGET_HH
#define PRODUCTIONTARGET_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {
    
  class ProductionTargetMaker;

  class ProductionTarget : public Detector {
  public: 

    // cylinder parameters
    double rOut() const { return _rOut; }
    double halfLength() const { return _halfLength; }

    // in mu2e coordinates
    const CLHEP::Hep3Vector& position() const { return _prodTargetPosition; }

    // this is used to transorm particle momentum and position from 
    // the PrimaryProtonGun frame to the Mu2e frame
    const CLHEP::HepRotation& protonBeamRotation() const { return _protonBeamRotation; }

    // "passive" rotation, used for placing the production target.
    // This is the inverse of protonBeamRotation.
    const CLHEP::HepRotation& productionTargetRotation() const { return _protonBeamInverseRotation; }

    //----------------------------------------------------------------
  private: 
    friend class ProductionTargetMaker;

    // Private ctr: the class should be only obtained via ProductionTargetFNAL::ProductionTargetMaker.
    ProductionTarget(double rOut, double halfLength, double rotX, double rotY, const CLHEP::Hep3Vector& position);

    CLHEP::HepRotation _protonBeamRotation;

    // can't return by const ref if invert on the fly so need to store redundant data
    CLHEP::HepRotation _protonBeamInverseRotation;

    CLHEP::Hep3Vector _prodTargetPosition;
    double _rOut;
    double _halfLength;
  };
}

#endif/*PRODUCTIONTARGET_HH*/
