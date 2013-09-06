// Geometry of the production target. This also defines the proton beam direction.
//
// Andrei Gaponenko, 2011

#ifndef PRODUCTIONTARGET_HH
#define PRODUCTIONTARGET_HH

#include <vector>
#include <map>

#include "art/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeomPrimitives/inc/Polycone.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProductionTargetMaker;

  class ProductionTarget : virtual public Detector {
  public:

    // cylinder parameters
    double rOut() const { return _rOut; }
    double halfLength() const { return _halfLength; }
    double envHalfLength() const { return _envelHalfLength; }

    // in mu2e coordinates
    const CLHEP::Hep3Vector& position() const { return _prodTargetPosition; }

    // this is used to transorm particle momentum and position from
    // the PrimaryProtonGun frame to the Mu2e frame
    const CLHEP::HepRotation& protonBeamRotation() const { return _protonBeamRotation; }

    // "passive" rotation, used for placing the production target.
    // This is the inverse of protonBeamRotation.
    const CLHEP::HepRotation& productionTargetRotation() const { return _protonBeamInverseRotation; }

    Polycone const * const getHubsRgtPtr() const {return _pHubsRgtParams.get();}
    Polycone const * const getHubsLftPtr() const {return _pHubsLftParams.get();}
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsRgt() const { return _anchoringPntsRgt; }
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsLft() const { return _anchoringPntsLft; }

    //----------------------------------------------------------------

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    ProductionTarget( ProductionTarget const & );
    ProductionTarget const & operator= ( ProductionTarget const & );

  private:
    friend class ProductionTargetMaker;

    // Private ctr: the class should be only obtained via ProductionTargetFNAL::ProductionTargetMaker.
    ProductionTarget(double rOut, double halfLength, double rotX, double rotY, const CLHEP::Hep3Vector& position);

    CLHEP::HepRotation _protonBeamRotation;

    // can't return by const ref if invert on the fly so need to store redundant data
    CLHEP::HepRotation _protonBeamInverseRotation;// FIXME: should be transient

    std::unique_ptr<Polycone> _pHubsRgtParams;
    std::unique_ptr<Polycone> _pHubsLftParams;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsRgt;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsLft;

    CLHEP::Hep3Vector _prodTargetPosition;
    double _rOut;
    double _halfLength;
    double _envelHalfLength;

    // Needed for persistency
    template<class T> friend class art::Wrapper;
    ProductionTarget() {}
  };
}

#endif/*PRODUCTIONTARGET_HH*/
