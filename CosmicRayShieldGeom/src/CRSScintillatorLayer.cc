//
// Representation of one Scintillator Layer in CosmicRayShield
//
//
// $Id: CRSScintillatorLayer.cc,v 1.6 2013/10/25 05:06:33 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/10/25 05:06:33 $
//
// Original author KLG based on Rob Kutschke's Layer
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#ifndef __CINT__

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e 
{

  CRSScintillatorLayer::CRSScintillatorLayer():
    _id(CRSScintillatorLayerId())
  {}

  CRSScintillatorLayer::CRSScintillatorLayer(CRSScintillatorLayerId const& id):
    _id(id)
  {}

  void CRSScintillatorLayer::getDimensions(std::vector<double> &halflengths, CLHEP::Hep3Vector &center) const
  {

    double min[3], max[3];
    for(int i=0; i<3; i++)
    {
      min[i]=NAN; max[i]=NAN;
    }

    std::vector<const CRSScintillatorBar*>::const_iterator ibar;
    for(ibar=_bars.begin(); ibar!=_bars.end(); ++ibar) 
    {
      const CRSScintillatorBar &bar = **ibar; 
      const CLHEP::Hep3Vector &barPosition = bar.getPosition();
      const std::vector<double> &barHalfLengths = bar.getHalfLengths();
      for(int i=0; i<3; i++)
      {
        double minPositionI = barPosition[i]-barHalfLengths[i];
        double maxPositionI = barPosition[i]+barHalfLengths[i];
        if(minPositionI<min[i] || std::isnan(min[i])) min[i]=minPositionI;
        if(maxPositionI>max[i] || std::isnan(max[i])) max[i]=maxPositionI;
      }
    }

    for(int i=0; i<3; i++)
    {
      halflengths[i]=(max[i]-min[i])/2.0;
      center[i]=(max[i]+min[i])/2.0;
    }

  }

  string CRSScintillatorLayer::name( string const& base ) const
  {
    ostringstream os;
    os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber() << "_"
       << _id.getLayerNumber();
    return os.str();
  }

} // namespace mu2e
#endif

