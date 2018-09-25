// Geometry of the production target. This also defines the proton beam direction.
//
// Andrei Gaponenko, 2011
// Update to v2 in 2018 by David Norvil Brown

#ifndef PRODUCTIONTARGET_HH
#define PRODUCTIONTARGET_HH

#include <map>

#include "canvas/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeomPrimitives/inc/Polycone.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProductionTargetMaker;

  class ProductionTarget : virtual public Detector {
  public:

    // cylinder parameters
    int    version() const { return _version; }
    double rOut() const { return _rOut; }
    double halfLength() const { return _halfLength; }
    double envHalfLength() const { return _envelHalfLength; }

    // Parameters for version 2 for stiffening parts
    double nFins       () const   { return _nFins;    }
    double finHeight   () const   { return _finHeight;    }
    double finThickness() const   { return _finThickness; }
    double hubDistUS() const      { return _hubDistUS;    }
    double hubDistDS() const      { return _hubDistDS;    }
    double hubAngleUS() const     { return _hubAngleUS;   }
    double hubAngleDS() const     { return _hubAngleDS;   }
    double hubOverhangUS() const  { return _hubOverhangUS;}
    double hubOverhangDS() const  { return _hubOverhangDS;}
    double hubLenUS() const       { return _hubLenUS;     }
    double hubLenDS() const       { return _hubLenDS;     }

    void   setHubLenUS(double& aLen ) { _hubLenUS = aLen; }
    void   setHubLenDS(double& aLen ) { _hubLenDS = aLen; }


    // in mu2e coordinates
    const CLHEP::Hep3Vector& position() const { return _prodTargetPosition; }

    // this is used to transorm particle momentum and position from
    // the PrimaryProtonGun frame to the Mu2e frame
    const CLHEP::HepRotation& protonBeamRotation() const { return _protonBeamRotation; }

    // "passive" rotation, used for placing the production target.
    // This is the inverse of protonBeamRotation.
    const CLHEP::HepRotation& productionTargetRotation() const { return _protonBeamInverseRotation; }

    const Polycone * getHubsRgtPtr() const {return _pHubsRgtParams/*.get()*/;}
    const Polycone * getHubsLftPtr() const {return _pHubsLftParams/*.get()*/;}
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsRgt() const { return _anchoringPntsRgt; }
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsLft() const { return _anchoringPntsLft; }

    ~ProductionTarget() {
           if (_pHubsRgtParams!=NULL) delete _pHubsRgtParams;
           if (_pHubsLftParams!=NULL) delete _pHubsLftParams;
    }

    //----------------------------------------------------------------

  private:
    friend class ProductionTargetMaker;

    // Private ctr: the class should be only obtained via ProductionTargetFNAL::ProductionTargetMaker.
    ProductionTarget(int version, double rOut, double halfLength, double rotX,
		     double rotY, const CLHEP::Hep3Vector& position, 
		     int    nFins,
		     double finHeight, double finThickness, 
		     double hubDistUS, double hubDistDS,
		     double hubAngleUS, double hubAngleDS,
		     double hubOverhangUS, double hubOverhangDS );

    CLHEP::HepRotation _protonBeamRotation;

    // can't return by const ref if invert on the fly so need to store redundant data
    CLHEP::HepRotation _protonBeamInverseRotation;// FIXME: should be transient

//    std::unique_ptr<Polycone> _pHubsRgtParams;
//    std::unique_ptr<Polycone> _pHubsLftParams;
    Polycone * _pHubsRgtParams;
    Polycone * _pHubsLftParams;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsRgt;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsLft;

    CLHEP::Hep3Vector _prodTargetPosition;

    int    _version;
    double _rOut;
    double _halfLength;
    double _envelHalfLength;
    
    // version 1+ parameters
    int    _nFins;
    double _finHeight;
    double _finThickness;
    double _hubDistUS;
    double _hubDistDS;
    double _hubAngleUS;
    double _hubAngleDS;
    double _hubOverhangUS;
    double _hubOverhangDS;
    double _hubLenUS;
    double _hubLenDS;

    // Needed for persistency
    template<class T> friend class art::Wrapper;
    ProductionTarget():_pHubsRgtParams(NULL), _pHubsLftParams(NULL) {}
  };
}

#endif/*PRODUCTIONTARGET_HH*/
