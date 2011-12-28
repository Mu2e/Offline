//
// Hold information about one Collimator in Extinction Monitor.
//
// $Id: ExtMonUCICol.hh,v 1.2 2011/12/28 00:25:05 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/28 00:25:05 $
//

#ifndef ExtMonUCICol_hh
#define ExtMonUCICol_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <algorithm>
#include <vector>

namespace mu2e {

  namespace ExtMonUCI {

    class ExtMonCol{

    friend class ExtMon;
    friend class ExtMonMaker;

    public:

      ExtMonCol():_id(-1) {}
      ExtMonCol(int& id):_id(id) {}
      ~ExtMonCol() {}

      // Formatted string embedding the id of the collimator.
      std::string name( std::string const& base ) const;

      // Accept the compiler generated destructor, copy constructor and assignment operators

      const int& getId() const { return _id;}

      const CLHEP::Hep3Vector&  origin()      const { return _origin; }
      const CLHEP::Hep3Vector&  originLocal() const { return _originLocal; }
      const CLHEP::HepRotation& rotation()    const { return _rotation; }

      const std::vector<double>& paramsOuter() const { return _paramsOuter; }
      const std::vector<double>& paramsInner() const { return _paramsInner; }

    protected:

      int _id;

      CLHEP::Hep3Vector _position1;
      CLHEP::Hep3Vector _position2;

      CLHEP::Hep3Vector   _origin;
      CLHEP::Hep3Vector   _originLocal;
      CLHEP::HepRotation  _rotation;

      std::vector<double> _paramsOuter;
      std::vector<double> _paramsInner;
    };

  } // namespace ExtMonUCI

} // namespace mu2e

#endif /* ExtMonCol_hh */
