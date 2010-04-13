#ifndef CRYSTAL_HH
#define CRYSTAL_HH
// $Id: Crystal.hh,v 1.3 2010/04/13 17:14:55 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:14:55 $

// original authors Julie Managan and Robert Bernstein

namespace mu2e{
  namespace calorimeter{

//
// C++ includes
#include <vector>

//
// Mu2e includes
#include "CalorimeterGeom/inc/CrystalId.hh"
#include "CalorimeterGeom/inc/CrystalIndex.hh"
#include "CalorimeterGeom/inc/CrystalDetail.hh"

//
//Other includes
#include "CLHEP/Vector/ThreeVector.h"


    class Crystal{

      friend class Layer;
      friend class Device;
      friend class Calorimeter;
      friend class CalorimeterMaker;


    public:

      // A free function, returning void, that takes a const Crystal& as an argument.
      typedef void (*CrystalFunction)( const Crystal& c);

      Crystal();

      // Constructor using wire tangents.
      Crystal( const CrystalId& id,
	       CrystalIndex index,
	       const CLHEP::Hep3Vector& c,
	       const CrystalDetail* detail,
	       int    detailIndex,
	       double wtx = 0.,
	       double wty = 0.
	       );

      // Constructor using wire unit vector.
      Crystal( const CrystalId& id,
	       CrystalIndex index,
	       const CLHEP::Hep3Vector& c,
	       const CrystalDetail* detail,
	       int detailIndex,
	       CLHEP::Hep3Vector const& t
	       );
  
      ~Crystal ();

      const CrystalId& Id() const { return _id;}
      CrystalIndex Index() const { return _index;}

      const CrystalDetail& getDetail() const { return *_detail;}

      const std::vector<const Crystal*>& nearestNeighbours() const{
	return _nearest;
      }
  
      // Compiler generated copy and assignment constructors
      // should be OK.
  
      // Nominal mid-point of the crystal, ignoring sag.
      const CLHEP::Hep3Vector& midPoint() const {return _cNominal;}

      // Unit vector in the nominal direction of the wire.
      const CLHEP::Hep3Vector& direction() const { return _w;}

      // Half length of the crystal.
      double zhalfLength() const { return _detail->zhalfLength();}

      int hack;
  
    protected:

      // Identifier
      CrystalId _id;

      // Index into the array of all crystals.
      CrystalIndex _index;

      // Nominal mid-point of the crystal, ignoring sag.
      CLHEP::Hep3Vector _cNominal;

      // Detailed description of a crystal.
      const CrystalDetail* _detail;
      int _detailIndex;

      // Unit vector along the wire direction.
      CLHEP::Hep3Vector _w;

      // Nearest neighbours.
      std::vector<const Crystal *> _nearest;

      std::vector<CrystalId> _nearestById;
      std::vector<CrystalIndex> _nearestByIndex;

    };
  } //namespace calorimeter
} //namespace mu2e

#endif
