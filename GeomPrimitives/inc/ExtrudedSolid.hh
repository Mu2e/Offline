#ifndef EXTRUDEDSOLID_HH
#define EXTRUDEDSOLID_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

namespace mu2e {

  class ExtrudedSolid {
  public:
    ExtrudedSolid() : yHalfThickness_() {}

    ExtrudedSolid( const std::string& name,
                   const std::string& mat,
                   const CLHEP::Hep3Vector& offset,
                   double yht,
                   const std::vector<CLHEP::Hep2Vector>& v )
      : name_(name)
      , material_(mat)
      , offsetFromMu2eOrigin_(offset)
      , yHalfThickness_(yht)
      , vertices_(v)
    {}

    const std::string& getName()     const { return name_; }
    const std::string& getMaterial() const { return material_; }

    const CLHEP::Hep3Vector& getOffsetFromMu2eOrigin() const { return offsetFromMu2eOrigin_; }

    double getYhalfThickness() const { return yHalfThickness_; }

    const std::vector<CLHEP::Hep2Vector>& getVertices() const { return vertices_; }

    CLHEP::Hep2Vector& modifyVertex( int i ) { return vertices_.at( (std::size_t)i); }

  private:
    std::string name_;
    std::string material_;
    CLHEP::Hep3Vector offsetFromMu2eOrigin_;
    double yHalfThickness_;
    std::vector<CLHEP::Hep2Vector> vertices_;
  };

}

#endif/*EXTRUDEDSOLID_HH*/
