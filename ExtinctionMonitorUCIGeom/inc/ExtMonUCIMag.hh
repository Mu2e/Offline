//
// Hold information about one Magnet in Extinction Monitor.
//
// $Id: ExtMonUCIMag.hh,v 1.3 2012/09/05 00:44:59 youzy Exp $
// $Author: youzy $
// $Date: 2012/09/05 00:44:59 $

#ifndef ExtMonUCIMag_hh
#define ExtMonUCIMag_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <algorithm>
#include <vector>

namespace mu2e {

  namespace ExtMonUCI {

    class ExtMonMag{

    friend class ExtMon;
    friend class ExtMonMaker;

    public:

      ExtMonMag():_id(-1) {}
      ExtMonMag(int& id):_id(id) {}
      ~ExtMonMag() {}

      // Formatted string embedding the id of the collimator.
      std::string name( std::string const& base ) const;

      // Accept the compiler generated destructor, copy constructor and assignment operators

      const int& getId() const { return _id;}

      const double& psi() const { return _psi; }
      const CLHEP::Hep3Vector&  origin()      const { return _origin; }
      const CLHEP::Hep3Vector&  originLocal() const { return _originLocal; }
      const CLHEP::HepRotation& rotation()    const { return _rotation; }

      const std::vector<double>& paramsOuter() const { return _paramsOuter; }
      const std::vector<double>& paramsInner() const { return _paramsInner; }

    protected:

      int _id;

      CLHEP::Hep3Vector _position1;
      CLHEP::Hep3Vector _position2;
      double _psi;

      CLHEP::Hep3Vector   _origin;
      CLHEP::Hep3Vector   _originLocal;
      CLHEP::HepRotation  _rotation;

      std::vector<double> _paramsOuter;
      std::vector<double> _paramsInner;
    };

  } // namespace ExtMonUCI

} // namespace mu2e

#endif /* ExtMonMag_hh */
