//
// The parameters of a Polycone
//
// $Id: Polycone.cc,v 1.5 2013/04/30 14:56:57 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/04/30 14:56:57 $
//
// Original author KLG
//

#include "GeomPrimitives/inc/Polycone.hh"

#include <iterator>
#include <algorithm>

#include "cetlib_except/exception.h"

namespace mu2e {

  Polycone::Polycone(const std::vector<double>& zPlanes,
                     const std::vector<double>& rInner,
                     const std::vector<double>& rOuter,
                     const CLHEP::Hep3Vector& originInMu2e,
                     const std::string& materialName,
                     double phiStart,
                     double phiTotal):
    _params ( zPlanes,
              rInner,
              rOuter,
              phiStart,
              phiTotal ),
    _originInMu2e(originInMu2e),
    _materialName(materialName)
  {
    // In our use the polycone shape parameters can come directly from
    // user inputs.  Therefore an assert() check is not appropriate,
    // we'll throw an exception to report problems.
    if(rInner.size() != zPlanes.size()) {
      throw cet::exception("GEOM")<<"mu2e::Polycone: (rInner.size() = "<<rInner.size()
                                  <<" does not match zPlanes.size() = "<<zPlanes.size()
                                  <<"\n";
    }
    if(rOuter.size() != zPlanes.size()) {
      throw cet::exception("GEOM")<<"mu2e::Polycone: (rOuter.size() = "<<rOuter.size()
                                  <<" does not match zPlanes.size() = "<<zPlanes.size()
                                  <<"\n";
    }
  };

  std::ostream& operator<<(std::ostream& os, const Polycone& p) {
    os<<"Polycone(numZPlanes = "<<p.numZPlanes()<<", z={";
    std::copy(p.zPlanes().begin(), p.zPlanes().end(), std::ostream_iterator<double>(os, ", "));
    os<<"}, rIn={";
    std::copy(p.rInner().begin(), p.rInner().end(), std::ostream_iterator<double>(os, ", "));
    os<<"}, rOut={";
    std::copy(p.rOuter().begin(), p.rOuter().end(), std::ostream_iterator<double>(os, ", "));
    os<<"}"
      <<", originInMu2e="<<p.originInMu2e()
      <<", phi0="<<p.phi0()
      <<", phiTotal="<<p.phiTotal()
      <<", materialName="<<p.materialName()
      <<")";
    return os;
  }

}
