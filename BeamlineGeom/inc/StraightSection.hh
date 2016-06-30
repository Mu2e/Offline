#ifndef BeamlineGeom_StraightSection_hh
#define BeamlineGeom_StraightSection_hh

//
// Class to represent the transport solenoid
//
#include <string>

#include "BeamlineGeom/inc/TSSection.hh"
#include "GeomPrimitives/inc/Tube.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class StraightSection : public TSSection {

  friend class BeamlineMaker;

  public:

    StraightSection() : TSSection(),
      _rIn(0.), _rOut(0.), _halfZ(0.) 
    {}

    StraightSection(double rIn, double rOut, double halfZ, 
                    CLHEP::Hep3Vector const & origin, 
                    CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                    std::string const & materialName = "") :
      TSSection(origin, rotation, materialName),
      _rIn(rIn), _rOut(rOut), _halfZ(halfZ)
    {}

    explicit StraightSection( const Tube& tube )  :
      TSSection(tube.originInMu2e(),  tube.rotation(), tube.materialName()),
      _rIn  ( tube.innerRadius() ),
      _rOut ( tube.outerRadius() ),
      _halfZ( tube.zHalfLength() ) 
    {}

    ~StraightSection(){}

    void set(double rIn, double rOut, double halfZ, 
             CLHEP::Hep3Vector const & origin, CLHEP::HepRotation const & rotation=CLHEP::HepRotation(),
             std::string const & materialName = "") {
      _rIn=rIn;
      _rOut=rOut;
      _halfZ=halfZ;
      _origin=origin;
      _rotation=rotation;
      _materialName=materialName;
    }

    double rIn() const { return _rIn; }
    double rOut() const { return _rOut; }
    double getHalfLength() const { return _halfZ; }

  private:

    double _rIn;
    double _rOut;
    double _halfZ;

  };

}
#endif /* BeamlineGeom_StraightSection_hh */
