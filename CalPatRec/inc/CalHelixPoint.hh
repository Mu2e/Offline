#ifndef CalHelixPoint_HH
#define CalHelixPoint_HH

// tracker
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

//-----------------------------------------------------------------------------
// utility struct
//-----------------------------------------------------------------------------
  struct CalHelixPoint {
    static double       _efac;
    static StrawHitFlag _useflag;	// flag bits to define use

    size_t             _ind;		// straw hit index
    CLHEP::Hep3Vector  _pos;		// position
    double             _phi;	        // ambiguity-resolved phi angle
    StrawHitFlag       _flag;		// flag
    int                _used;           // = 1 if the strawhit is used by another track-candidate
    CLHEP::Hep3Vector  _wdir;		// wire direction
    CLHEP::Hep3Vector  _sdir;           // straw radial direction, perp to Z and wire direction
					// errors are asymmetric; along the wire is given by time division, 
					// perp to the wire by the straw size/sqrt(12)
    const Straw*       _straw;          // pointer to the straw
    const StrawHit*    _strawhit;       // pointer to the strawHit
    double             _perr;
    double             _rerr;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
					// initialize some variables on construction
    CalHelixPoint():_phi(0.0) {}

    CalHelixPoint(size_t ind, const StrawHit& sh, const StrawHitPosition& shp, const Straw& straw, const StrawHitFlag& flag);

    CalHelixPoint(size_t ind, const CLHEP::Hep3Vector& pos, const CLHEP::Hep3Vector& wdir, double werr, double serr);
 
    bool         use       () const;
    int          isUsed    () { return _used;}
    bool         stereo    () const;
    bool         isOutlier () const;
    bool         isCalosel () const;
    void         setOutlier();
    void         setUse    (bool use);

    double       x         () const { return _pos.x(); }
    double       y         () const { return _pos.y(); }
    double       z         () const { return _pos.z(); }
  };

};
#endif

