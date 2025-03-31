#ifndef BeamlineGeom_TorusSection_hh
#define BeamlineGeom_TorusSection_hh

//
// Class to represent the transport solenoid
//
#include <array>
#include <string>

#include "Offline/GeomPrimitives/inc/Torus.hh"
#include "Offline/BeamlineGeom/inc/TSSection.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TorusSection : public TSSection {

  friend class BeamlineMaker;

  public:

    constexpr static std::size_t _npar_data{5};

    // fixme: improve  _materialName initialization
    TorusSection() :
      _rTorus(0.),_rIn(0.),_rOut(0.),
      _phiBegin(0.),_deltaPhi(0.)
    {
      fillData();
    }


    TorusSection(double rTorus, double rIn, double rOut, double phi0, double dPhi,
                 CLHEP::Hep3Vector const & origin,
                 CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                 std::string const & materialName = ""):
      TSSection(origin, rotation, materialName),
      _rTorus(rTorus),_rIn(rIn),_rOut(rOut),
      _phiBegin(phi0),_deltaPhi(dPhi)
    {
      fillData();
    }

    TorusSection(const Torus& torus) :
      TSSection(torus.originInMu2e(), torus.rotation(), torus.materialName()),
      _rTorus(torus.torusRadius()),
      _rIn(torus.innerRadius()),
      _rOut(torus.outerRadius()),
      _phiBegin(torus.phi0()),
      _deltaPhi(torus.deltaPhi())
    {
      fillData();
    }

    void set(double rTorus, double rIn, double rOut, double phi0, double dPhi,
             CLHEP::Hep3Vector  const & origin,
             CLHEP::HepRotation  const & rotation = CLHEP::HepRotation(),
             std::string const & materialName = "") {
      _rTorus  =rTorus;
      _rIn     =rIn;
      _rOut    =rOut;
      _phiBegin=phi0;
      _deltaPhi=dPhi;
      _origin=origin;
      _rotation=rotation;
      _materialName=materialName;
      fillData();
    }

    double torusRadius() const {return _rTorus;}
    double rIn()         const {return _rIn;   }
    double rOut()        const {return _rOut;  }
    double phiStart()    const {return _phiBegin; }
    double deltaPhi()    const {return _deltaPhi; }
    const std::array<double,_npar_data>& getParameters() const { return _data; }

  private:

    // All dimensions in mm.
    double _rTorus;
    double _rIn;
    double _rOut;

    double _phiBegin;
    double _deltaPhi;

    std::array<double,_npar_data> _data;

    void fillData() {
      _data[0] = _rIn;
      _data[1] = _rOut;
      _data[2] = _rTorus;
      _data[3] = _phiBegin;
      _data[4] = _deltaPhi;
    }

};

}
#endif /* BeamlineGeom_TorusSection_hh */
