//
// Hold information about one shielding block in Extinction Monitor.
//
// $Id: ExtMonUCIShd.hh,v 1.1 2012/02/16 20:25:46 youzy Exp $
// $Author: youzy $
// $Date: 2012/02/16 20:25:46 $

#ifndef ExtMonUCIShd_hh
#define ExtMonUCIShd_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <algorithm>
#include <vector>

namespace mu2e {

  namespace ExtMonUCI {

    class ExtMonShd{

    friend class ExtMon;
    friend class ExtMonMaker;

    public:

      ExtMonShd():_id(-1) {}
      ExtMonShd(int& id):_id(id) {}
      ~ExtMonShd() {}

      // Formatted string embedding the id of the shd.
      std::string name( std::string const& base ) const;

      // Accept the compiler generated destructor, copy constructor and assignment operators

      const int& getId() const { return _id;}
      void setId(int id) { _id = id; }

      const CLHEP::Hep3Vector&  origin()      const { return _origin; }
      const CLHEP::Hep3Vector&  originLocal() const { return _originLocal; }
      const CLHEP::HepRotation& rotation()    const { return _rotation; }

      const std::vector<double>& params() const { return _params; }

    protected:

      int _id;

      CLHEP::Hep3Vector   _origin;
      CLHEP::Hep3Vector   _originLocal;
      CLHEP::HepRotation  _rotation;

      std::vector<double> _params;
    };

  } // namespace ExtMonUCI

} // namespace mu2e

#endif /* ExtMonShd_hh */
