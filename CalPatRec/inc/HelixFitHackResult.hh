//
#ifndef HelixFitHackResult_HH
#define HelixFitHackResult_HH

#include "CalPatRec/inc/LsqSums4.hh"
#include "CalPatRec/inc/HelixDefHack.hh"
#include "TrkBase/TrkErrCode.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  struct HelixFitHackResult {
    HelixDefHack       _hdef;         // must copy by value as references can't be re-assigned
    TrkErrCode         _fit;	      // fit status code from last fit
//-----------------------------------------------------------------------------
// circle parameters; the z center is ignored.
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxy;
    ::LsqSums4         _srphi;

    CLHEP::Hep3Vector  _center;
    double             _radius;

    double             _chi2;
//-----------------------------------------------------------------------------
// 2015-02-06 P.Murat: fit with non-equal weights - XY-only
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxyw;
    CLHEP::Hep3Vector  _cw;
    double             _rw;
    double             _chi2w;
//-----------------------------------------------------------------------------
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)), 
// fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
//-----------------------------------------------------------------------------
    double             _dfdz;
    double             _fz0;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    HelixFitHackResult();

    HelixFitHackResult(TrkDef const& tdef) : 
      _hdef(tdef),  
      _fit(TrkErrCode::fail),
      _radius(-1.0),
      _dfdz(0.0),
      _fz0(0.0) {
    }

    HelixFitHackResult(HelixDefHack const& hdef) : 
      _hdef(hdef),  
      _fit(TrkErrCode::fail),
      _radius(-1.0),
      _dfdz(0.0),
      _fz0(0.0) {
    }

    HelixFitHackResult(const HelixFitHackResult& hfitRes);

    HelixFitHackResult& operator =(HelixFitHackResult const& other);

    bool fitIsValid() { return _sxy.qn() > 0; }

    bool weightedFitIsValid() { return _sxyw.qn() > 0; }

    void print(const char* Title);

    void clear();
  };

};
#endif

