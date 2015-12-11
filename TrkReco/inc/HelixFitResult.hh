//
// Helix fit output 
//
// $Id: HelixFit.hh,v 1.8 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
#ifndef HelixFitResult_HH
#define HelixFitResult_HH
#include "TrkReco/inc/HelixDef.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
namespace mu2e 
{
// output struct
  struct HelixFitResult {
    HelixDef _hdef; // must copy by value as references can't be re-assigned
// fit status
    TrkErrCode _fit; // error code from last fit
// circle parameters; the z center is ignored.
    CLHEP::Hep3Vector _center;
    double _radius;
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)), fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
    double _dfdz, _fz0;
    HelixFitResult(TrkDef const& tdef) : _hdef(tdef),  _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
    HelixFitResult(HelixDef const& hdef) : _hdef(hdef),  _fit(TrkErrCode::fail),_radius(-1.0),_dfdz(0.0),_fz0(0.0) {}
    HelixFitResult& operator =(HelixFitResult const& other);
 };
}
#endif

