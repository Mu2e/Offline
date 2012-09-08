#ifndef CalorimeterGeom_Vane_hh
#define CalorimeterGeom_Vane_hh

//
// Hold information about one vane in the calorimter.
//

//
// $Id: Vane.hh,v 1.9 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author R, Bernstein and Rob Kutschke
//

#include <vector>
//#include <iostream>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

    class Vane{

      friend class VaneCalorimeter;
      friend class VaneCalorimeterMaker;

    public:

      Vane():_id(-1){}
      Vane( int & id ):_id(id){}
      ~Vane(){}

      // Compiler generated copy and assignment constructors
      // should be OK.

      int id() const { return _id;}

      // Get position in the global Mu2e frame
      CLHEP::Hep3Vector const& getOrigin() const      { return _origin; }
      CLHEP::Hep3Vector const& getOriginLocal() const { return _originLocal; }
      CLHEP::Hep3Vector const& getSize() const        { return _size; }
      CLHEP::HepRotation * getRotation() const {
        return const_cast<CLHEP::HepRotation *>(&_rotation);
      }

    protected:

      int _id;
      CLHEP::Hep3Vector _origin;
      CLHEP::Hep3Vector _originLocal;
      CLHEP::Hep3Vector _size;
      CLHEP::HepRotation _rotation;

    };

} //namespace mu2e

#endif /* CalorimeterGeom_Vane_hh */
