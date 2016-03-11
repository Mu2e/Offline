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
                                                     std::vector<int> const& localToWorld,
                                                     std::string const& CMBmaterialName,
                                                     double CMBoffset, double CMBhalfThickness,
                                                     bool CMBside0, bool CMBside1) :
    _materialName(materialName),
    _halfLengths(halfLengths),
    _localToWorld(localToWorld),
    _CMBmaterialName(CMBmaterialName),
    _CMBoffset(CMBoffset),
    _CMBhalfThickness(CMBhalfThickness),
    _CMBside0(CMBside0),
    _CMBside1(CMBside1)
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


  /********************/
  // counter motherboard section

  CLHEP::Hep3Vector CRSScintillatorBarDetail::getCMBPosition(int side, const CLHEP::Hep3Vector &barPosition) const
  {
    int CMBcoordinate=_localToWorld[2];
    CLHEP::Hep3Vector CMBdifference;
    CMBdifference[CMBcoordinate]=_halfLengths[CMBcoordinate]+_CMBoffset;

    if(side==0) CMBdifference*=-1; //side==0 is for the negative side of the counter, and side==1 is for the positive side of the counter
    CLHEP::Hep3Vector CMBposition=barPosition+CMBdifference;
    return CMBposition;
  }

  std::vector<double> CRSScintillatorBarDetail::getCMBHalfLengths() const
  {
    int CMBcoordinate=_localToWorld[2];
    std::vector<double> CMBhalfLengths=_halfLengths;
    CMBhalfLengths[CMBcoordinate]=_CMBhalfThickness;
    return CMBhalfLengths;
  }

  bool CRSScintillatorBarDetail::hasCMB(int side) const
  {
    switch(side)
    {
      case 0: return _CMBside0;
      case 1: return _CMBside1;
    }
    return false;
  }
} // namespace mu2e
