#ifndef TrackerGeom_Straw_hh
#define TrackerGeom_Straw_hh
//
// Hold information about one straw in a tracker.
//
// Original author Rob Kutschke
//

#include <vector>
#include <array>
#include "cetlib_except/exception.h"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Tracker;

  class Straw{
    public:
      using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

      Straw():_hlen(-1.0) {}; // default constructor should be deleted, it results in a non-fuctional object FIXME!

      // Constructor using wire tangents.
      Straw( const StrawId& id,
          const xyzVec& c,
          double halflen,
          double wtx = 0.,
          double wty = 0.
          );

      // Constructor using wire unit vector.
      Straw( const StrawId& id,
          const xyzVec& c,
          xyzVec const& t,
          double halflen
          );
      // aligned tracker constructor takes the end positions
      Straw (const StrawId& id,
          const xyzVec& calwireend, const xyzVec& hvwireend,
          const xyzVec& calstrawend, const xyzVec& hvstrawend);

      Straw (const StrawId& id,
          const xyzVec wireends[2],
          const xyzVec strawends[2]);

      // Accept the compiler copy constructor and assignment operators


      const StrawId& id() const { return _id;}

      // Formatted string embedding the id of the straw.
      std::string name( std::string const& base ) const;

      // Half length: same for straws and wires
      float halfLength() const { return _hlen; }

      // end positions
      //
      xyzVec wireEnd(StrawEnd const& end ) const { return wireEnd(end.end()); }
      xyzVec wireEnd(StrawEnd::End end ) const { return end == StrawEnd::cal ? _wmid + _hlen*_wdir : _wmid - _hlen*_wdir; }
      xyzVec strawEnd(StrawEnd const& end ) const { return strawEnd(end.end()); }
      xyzVec strawEnd(StrawEnd::End end ) const { return end == StrawEnd::cal ? _smid + _hlen*_sdir : _smid - _hlen*_sdir; }

      bool operator==(const Straw& other) const {
        return _id == other.id();
      }
      bool operator>(const Straw& other) const {
        return _id > other.id();
      }
      bool operator<(const Straw& other) const {
        return _id < other.id();
      }

      // aligned interface: these provide more refined information.  They call down
      // to the unaligned information in case this is the unaligned tracker
      // These take the position along the straw (WRT the center) as argument, accounting for misalignment,
      // and (eventually) gravitational sag and electrostatic displacement.
      xyzVec wireDirection(float upos=0.0) const { return _wdir; }
      xyzVec strawDirection(float upos=0.0) const { return _sdir; }
      xyzVec wirePosition(float upos=0.0) const { return _wmid + upos*_wdir; }
      xyzVec strawPosition(float upos=0.0) const { return _smid + upos*_sdir; }
      // define the origin as the wire position in the middle
      xyzVec const& origin() const { return _wmid; }

      // deprecated interface
      xyzVec const& direction() const { return _wdir; }
      xyzVec const& getMidPoint() const {return _wmid; }
      xyzVec const& getDirection() const { return _wdir; }
    private:
      // Identifier
      StrawId _id;
      // wire and straw midpoints.
      xyzVec _wmid, _smid;
      // wire and straw directions
      xyzVec _wdir, _sdir;
      // half length
      float _hlen;
  };

}  //namespace mu2e
#endif /* TrackerGeom_Straw_hh */
