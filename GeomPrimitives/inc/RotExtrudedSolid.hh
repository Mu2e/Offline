#ifndef ROTEXTRUDEDSOLID_HH
#define ROTEXTRUDEDSOLID_HH

#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class RotExtrudedSolid {
  public:
    RotExtrudedSolid() : yHalfThickness_() {}

    RotExtrudedSolid( const std::string& name,
                      const std::string& mat,
                      const CLHEP::Hep3Vector& offset,
                      double yht,
                      const std::vector<CLHEP::Hep2Vector>& v,
                      CLHEP::HepRotation const & rotation = CLHEP::HepRotation() )
      : name_(name)
      , material_(mat)
      , offsetFromMu2eOrigin_(offset)
      , yHalfThickness_(yht)
      , vertices_(v),
      _rotation    ( rotation     )
    {}

    const std::string& getName()     const { return name_; }
    const std::string& getMaterial() const { return material_; }

    const CLHEP::Hep3Vector& getOffsetFromMu2eOrigin() const { return offsetFromMu2eOrigin_; }

    double getYhalfThickness() const { return yHalfThickness_; }

    const std::vector<CLHEP::Hep2Vector>& getVertices() const { return vertices_; }

    CLHEP::HepRotation const & getRotation()     const { return _rotation; }

    CLHEP::Hep2Vector& modifyVertex( int i ) { return vertices_.at( (std::size_t)i); }

  private:
    std::string name_;
    std::string material_;
    CLHEP::Hep3Vector offsetFromMu2eOrigin_;
    double yHalfThickness_;
    std::vector<CLHEP::Hep2Vector> vertices_;
    CLHEP::HepRotation             _rotation; // wrt to parent volume
  };

}

#endif/*ROTEXTRUDEDSOLID_HH*/
