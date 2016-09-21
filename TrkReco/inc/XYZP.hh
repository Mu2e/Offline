//
//  Information about position for helix finding
//  Original author: Dave Brown (LBNL)  2014
//
#ifndef TrkReco_XYZP_HH
#define TrkReco_XYZP_HH
#include "CLHEP/Vector/ThreeVector.h"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
namespace mu2e {
  class StrawHit;
  class StrawHitPosition;
  class Straw;
// utility struct; value plus error
  struct VALERR {
    double _val;
    double _err;
  };
  struct FZ {
    VALERR _phi;
    double _z;
  };
  class XYZP;
  typedef std::vector<XYZP> XYZPVector;
// utility struct
  struct XYZP {
// straw hit index
    int _ind;
// position
    CLHEP::Hep3Vector _pos;
// ambiguity-resolved azimuth WRT the center of the helix
    double _phi;
// flag
    StrawHitFlag _flag;
// wire direction
    CLHEP::Hep3Vector _wdir;
// straw radial direction, perp to Z and wire direction
    CLHEP::Hep3Vector _sdir;
// errors are asymmetric; along the wire is given by time division, perp to the wire by the straw size/sqrt(12)
    double _perr,_rerr;
// initialize some variables on construction
    XYZP():_ind(-1),_phi(0.0),_perr(1000.0),_rerr(1000.0){}
    XYZP(size_t index,StrawHit const& sh, StrawHitPosition const& shp,Straw const& straw);
    XYZP(size_t ind,CLHEP::Hep3Vector const& pos, CLHEP::Hep3Vector const& wdir, double werr, double serr);
    XYZP(CLHEP::Hep3Vector const& pos, double size);
// radial position information
    virtual void rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const;
    virtual void finfo(CLHEP::Hep3Vector const& center, VALERR& phi) const;
    bool use() const;
    bool stereo() const;
    void setOutlier();
    void setUse(bool use);
// flag bits to define use
    static StrawHitFlag _useflag, _dontuseflag;
    static int _debug;
     // fill function; make it static to keep the namespace clean
    static void fillXYZP(StrawHitCollection const& shcol,
	StrawHitPositionCollection const& shpos, std::vector<size_t>, XYZPVector& xyzp);
  };
}

#endif
