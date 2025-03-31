#ifndef BeamlineGeom_ConeSection_hh
#define BeamlineGeom_ConeSection_hh

//
// Class to represent the transport solenoid
//
#include <string>

#include "Offline/BeamlineGeom/inc/TSSection.hh"
#include "Offline/GeomPrimitives/inc/Cone.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class ConeSection : public TSSection {

    friend class BeamlineMaker;

  public:

    constexpr static std::size_t _npar_data = 7;

    // fixme: improve  _materialName initialization
    ConeSection() : TSSection(),
                    _rIn1(0.), _rOut1(0.),  _rIn2(0.), _rOut2(0.), _halfZ(0.),
                    _phi0(0.), _deltaPhi(CLHEP::twopi)
    {
      fillData();
    }

    ConeSection(double rIn1, double rOut1, double rIn2, double rOut2, double halfZ,
                double phi0, double deltaPhi,
                CLHEP::Hep3Vector const & origin,
                CLHEP::HepRotation  const & rotation = CLHEP::HepRotation(),
                std::string const & materialName = "") :
      TSSection(origin, rotation, materialName),
      _rIn1(rIn1), _rOut1(rOut1),  _rIn2(rIn2), _rOut2(rOut2), _halfZ(halfZ),
      _phi0(phi0), _deltaPhi(deltaPhi)
    {
      fillData();
    }

    ConeSection( const Cone& cone )  :
      TSSection(cone.originInMu2e(),  cone.rotation(), cone.materialName()),
      _rIn1  ( cone.innerRadius1() ),
      _rOut1 ( cone.outerRadius1() ),
      _rIn2  ( cone.innerRadius2() ),
      _rOut2 ( cone.outerRadius2() ),
      _halfZ ( cone.zHalfLength()  ),
      _phi0  ( cone.phi0()         ),
      _deltaPhi ( cone.deltaPhi()  )
    {
      fillData();
    }

    void set(double rIn1, double rOut1, double rIn2, double rOut2, double halfZ,
             double phi0, double deltaPhi,
             CLHEP::Hep3Vector const & origin, CLHEP::HepRotation const & rotation=CLHEP::HepRotation(),
             std::string const & materialName = "" ) {
      _rIn1=rIn1;
      _rOut1=rOut1;
      _rIn2=rIn2;
      _rOut2=rOut2;
      _halfZ=halfZ;
      _phi0=phi0;
      _deltaPhi=deltaPhi;
      _origin=origin;
      _rotation=rotation;
      _materialName=materialName;
      fillData();
    }

    double rIn1()  const { return _rIn1; }
    double rOut1() const { return _rOut1; }
    double rIn2()  const { return _rIn2; }
    double rOut2() const { return _rOut2; }
    double getHalfLength() const { return _halfZ; }
    double phiStart()    const {return _phi0; }
    double deltaPhi()    const {return _deltaPhi; }
    const std::array<double,_npar_data>& getParameters() const { return _data; }

  private:

    double _rIn1;
    double _rOut1;
    double _rIn2;
    double _rOut2;
    double _halfZ;
    double _phi0;
    double _deltaPhi;

    std::array<double,_npar_data> _data;

    void fillData() {
      _data[0] = _rIn1;
      _data[1] = _rOut1;
      _data[2] = _rIn2;
      _data[3] = _rOut2;
      _data[4] = _halfZ;
      _data[5] = _phi0;
      _data[6] = _deltaPhi;
    }

  };

}
#endif /* BeamlineGeom_ConeSection_hh */
