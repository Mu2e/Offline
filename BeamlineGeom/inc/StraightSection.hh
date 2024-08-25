#ifndef BeamlineGeom_StraightSection_hh
#define BeamlineGeom_StraightSection_hh

//
// Class to represent the transport solenoid
//
#include <string>

#include "Offline/BeamlineGeom/inc/TSSection.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class StraightSection : public TSSection {

  friend class BeamlineMaker;

  public:

    StraightSection() : TSSection(),
                        _rIn(0.), _rOut(0.), _halfZ(0.), _diffZ(0.)
    {}

    StraightSection(double rIn, double rOut, double halfZ,
                    CLHEP::Hep3Vector const & origin,
                    CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                    std::string const & materialName = "", double diffZ=0.) :
      TSSection(origin, rotation, materialName),
      _rIn(rIn), _rOut(rOut), _halfZ(halfZ), _diffZ(diffZ)
    {}

    explicit StraightSection( const Tube& tube )  :
      TSSection(tube.originInMu2e(),  tube.rotation(), tube.materialName()),
      _rIn  ( tube.innerRadius() ),
      _rOut ( tube.outerRadius() ),
      _halfZ( tube.zHalfLength() ),
      _diffZ(0.)
    {}

    void set(double rIn, double rOut, double halfZ,
             CLHEP::Hep3Vector const & origin, CLHEP::HepRotation const & rotation=CLHEP::HepRotation(),
             std::string const & materialName = "", double diffZ = 0.) {
      _rIn=rIn;
      _rOut=rOut;
      _halfZ=halfZ;
      _origin=origin;
      _rotation=rotation;
      _materialName=materialName;
      _diffZ=diffZ;
    }

    double rIn() const { return _rIn; }
    double rOut() const { return _rOut; }
    double getHalfLength() const { return _halfZ; }
    double getLengthDiff() const { return _diffZ; }
    void   setLengthDiff(double zDiff) { _diffZ=zDiff; }

  private:

    double _rIn;
    double _rOut;
    double _halfZ;
    double _diffZ; //amount to end before high Z end

  };

}
#endif /* BeamlineGeom_StraightSection_hh */
