//
// Hold information about one Tof (Station/Segment) in Extinction Monitor.
//
// $Id: ExtMonUCITof.hh,v 1.2 2011/12/28 00:25:05 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/28 00:25:05 $

#ifndef ExtMonUCITof_hh
#define ExtMonUCITof_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <algorithm>
#include <vector>

namespace mu2e {

  namespace ExtMonUCI {

    class ExtMonTof{

    friend class ExtMon;
    friend class ExtMonMaker;

    public:

      ExtMonTof():_id(-1) {}
      ExtMonTof(int& id):_id(id) {}
      ExtMonTof(int& station, int& segment):_station(station),_segment(segment) {}
      ~ExtMonTof() {}

      // Formatted string embedding the id of the tof.
      std::string name( std::string const& base ) const;

      // Accept the compiler generated destructor, copy constructor and assignment operators

      const int& getStation()  { return _station; }
      const int& getSegment()  { return _segment; }
      const int& getId() const { return _id;}
      void setId(int id) { _id = id; }

      const CLHEP::Hep3Vector&  origin()      const { return _origin; }
      const CLHEP::Hep3Vector&  originLocal() const { return _originLocal; }
      const CLHEP::HepRotation& rotation()    const { return _rotation; }

      const std::vector<double>& params() const { return _params; }

    protected:

      int _station;
      int _segment;
      int _id;

      CLHEP::Hep3Vector   _origin;
      CLHEP::Hep3Vector   _originLocal;
      CLHEP::HepRotation  _rotation;

      std::vector<double> _params;
    };

  } // namespace ExtMonUCI

} // namespace mu2e

#endif /* ExtMonTof_hh */
