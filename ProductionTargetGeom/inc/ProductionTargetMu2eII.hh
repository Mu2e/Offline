// Geometry of the Mu2e-II production target. 
// This also defines the proton beam direction.
//
// Based on Andrei Gaponenko, 2011 ProductionTarget.hh
// Created by Michael MacKenzie, 2020


#ifndef PRODUCTIONTARGETMU2EII_HH
#define PRODUCTIONTARGETMU2EII_HH

#include <map>

#include "canvas/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeomPrimitives/inc/Polycone.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProductionTargetMaker;

  class ProductionTargetMu2eII : virtual public Detector {
  public:

    // cylinder parameters
    int    version() const { return _version; }
    // double rOut() const { return _rOut; }
    // double halfLength() const { return _halfLength; }
    // double envHalfLength() const { return _envelHalfLength; }

    bool   isRotating() const {return _isRotating;}
    bool   isConveyor() const {return _isConveyor;}


    // in mu2e coordinates
    const CLHEP::Hep3Vector& position() const { return _prodTargetPosition; }

    // this is used to transorm particle momentum and position from
    // the PrimaryProtonGun frame to the Mu2e frame
    const CLHEP::HepRotation& protonBeamRotation() const { return _protonBeamRotation; }

    // "passive" rotation, used for placing the production target.
    // This is the inverse of protonBeamRotation.
    const CLHEP::HepRotation& productionTargetRotation() const { return _protonBeamInverseRotation; }

    //Conveyor parameters
    double      conveyorBallRadius   () const { return _conveyorBallRadius  ; }
    double      conveyorBallGap      () const { return _conveyorBallGap     ; }
    std::string conveyorBallMaterial () const { return _conveyorBallMaterial; }
    int         conveyorNBalls       () const { return _conveyorNBalls      ; }
    const std::vector<double>& conveyorBallXs() const { return _conveyorBallXs; }
    const std::vector<double>& conveyorBallYs() const { return _conveyorBallYs; }
    const std::vector<double>& conveyorBallZs() const { return _conveyorBallZs; }


    //Rotating parameters
    int rotatingNRods () const { return _rotatingNRods; }

    
    ~ProductionTargetMu2eII() {}


    std::string targetCoreMaterial()         const  {return _targetCoreMaterial;}
    std::string targetVacuumMaterial()       const  {return _targetVacuumMaterial;}



    //----------------------------------------------------------------

  private:
    friend class ProductionTargetMaker;

    // Private ctr: the class should be only obtained via ProductionTargetFNAL::ProductionTargetMaker.
    ProductionTargetMu2eII(std::string type, int version);

    CLHEP::HepRotation _protonBeamRotation;

    // can't return by const ref if invert on the fly so need to store redundant data
    CLHEP::HepRotation _protonBeamInverseRotation;// FIXME: should be transient


    CLHEP::Hep3Vector _prodTargetPosition;

    std::string _type;
    int    _version;
    bool   _isRotating;
    bool   _isConveyor;

    std::string _targetCoreMaterial;
    std::string _targetVacuumMaterial;

    //conveyor parameters
    int         _conveyorNBalls;
    double      _conveyorBallRadius;
    double      _conveyorBallGap;
    std::string _conveyorBallMaterial;
    std::vector<double> _conveyorBallXs; //ball centers
    std::vector<double> _conveyorBallYs;
    std::vector<double> _conveyorBallZs;

    //rotating paramerters
    double _rotatingNRods;
    
    //beam parameters
    double _beamRotX;
    double _beamRotY;
    double _beamRotZ;

    
    // Needed for persistency
    template<class T> friend class art::Wrapper;
    ProductionTargetMu2eII() {}
  };
}

#endif/*PRODUCTIONTARGETMU2UII_HH*/
