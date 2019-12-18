#ifndef GeomPrimitives_GenericTrap_hh
#define GeomPrimitives_GenericTrap_hh
//
// Copied from Torus description, with some method naming to match ExtrudedSolid
//
// "Typical" Generic Trapezoid object
//
// $Id: GenericTrap.hh,v 1.1 2013/08/07 20:20:07 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/08/07 20:20:07 $
//
// Original author: Kyle Knoepfel
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "cetlib_except/exception.h"

namespace mu2e {

  class GenericTrap {

  public:

    GenericTrap(std::string const & name,
		std::string const & materialName,
		CLHEP::Hep3Vector const & originInMu2e,
		double const dz, std::vector<CLHEP::Hep2Vector> const & vertices,
		CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
		) :
      _name(name),
      _materialName( materialName ),
      _originInMu2e( originInMu2e ),
      _dz(dz),
      _vertices(vertices),
      _rotation    ( rotation     )
    {
      if(_vertices.size() != 8)
	throw cet::exception("GEOM")<< "mu2e::GenericTrap: Need exactly 8 vertices, given "
				    << _vertices.size() << "\n";

    }
    
    std::string const & getName() const { return _name; }
    std::string const & getMaterial() const { return _materialName; }
    
    CLHEP::Hep3Vector const  & getOffsetFromMu2eOrigin() const { return _originInMu2e; }

    double getYhalfThickness() const { return _dz; }
    std::vector<CLHEP::Hep2Vector> getVertices() const { return _vertices; }
    
    CLHEP::HepRotation const & getRotation()     const { return _rotation; }

    CLHEP::Hep2Vector& modifyVertex( int i ) { return _vertices.at( (std::size_t)i); }
    
    // Genreflex can't do persistency of vector<GenericTrap> without a default constructor
    GenericTrap() {}

  private:

    std::string                    _name;
    std::string                    _materialName;
    CLHEP::Hep3Vector              _originInMu2e;
    double                         _dz; 
    std::vector<CLHEP::Hep2Vector> _vertices;
    CLHEP::HepRotation             _rotation; // wrt to parent volume

  };

}

#endif/*GeomPrimitives_GenericTrap_hh*/
