//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.cc,v 1.5 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG; somewhat based on Rob Kutschke StrawDetail
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"

using namespace std;

namespace mu2e 
{
  CRSScintillatorBarDetail::CRSScintillatorBarDetail(std::string const& materialName,
                                                     std::vector<double> const& halfLengths,
                                                     std::vector<int> const& localToWorld) :
    _materialName(materialName),
    _halfLengths(halfLengths),
    _localToWorld(localToWorld)
  {
  }

  CLHEP::Hep3Vector CRSScintillatorBarDetail::toWorld(const CLHEP::Hep3Vector &localPosition, 
                                                      const CLHEP::Hep3Vector &barPosition) const
  {
    CLHEP::Hep3Vector worldPosition=barPosition;
    for(int i=0; i<3; i++) worldPosition[_localToWorld[i]]+=localPosition[i];
    return worldPosition;
  }

  CLHEP::Hep3Vector CRSScintillatorBarDetail::toLocal(const CLHEP::Hep3Vector &worldPosition, 
                                                      const CLHEP::Hep3Vector &barPosition) const
  {
    CLHEP::Hep3Vector localPosition;
    CLHEP::Hep3Vector tmp=worldPosition-barPosition;
    for(int i=0; i<3; i++) localPosition[i]=tmp[_localToWorld[i]];
    return localPosition;
  }

  CLHEP::Hep3Vector CRSScintillatorBarDetail::toLocalNormalized(const CLHEP::Hep3Vector &worldPosition, 
                                                                const CLHEP::Hep3Vector &barPosition) const
  {
    CLHEP::Hep3Vector localPosition;
    CLHEP::Hep3Vector tmp=worldPosition-barPosition;
    for(int i=0; i<3; i++) localPosition[i]=tmp[_localToWorld[i]]/_halfLengths[_localToWorld[i]];
    return localPosition;
  }

  bool CRSScintillatorBarDetail::isInside(const CLHEP::Hep3Vector &worldPosition, 
                                          const CLHEP::Hep3Vector &barPosition) const
  {
    CLHEP::Hep3Vector tmp=worldPosition-barPosition;
    for(int i=0; i<3; i++) 
    {
      if(abs(tmp[i]/_halfLengths[i])>1.0) return false;
    }
    return true;
  }

} // namespace mu2e
