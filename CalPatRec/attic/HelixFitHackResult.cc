//

#include "CalPatRec/inc/HelixFitHackResult.hh"


namespace mu2e {

//-----------------------------------------------------------------------------
// HelixFitHackResult
//-----------------------------------------------------------------------------
  HelixFitHackResult::HelixFitHackResult() {
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  HelixFitHackResult::HelixFitHackResult(const HelixFitHackResult & hfitRes) :
    _hdef(hfitRes._hdef),
    _fit(hfitRes._fit),
    
    _sxy(hfitRes._sxy),
    _srphi(hfitRes._srphi),
    _center(hfitRes._center),
    _radius(hfitRes._radius),
    _chi2(hfitRes._chi2),
    
    _sxyw(hfitRes._sxyw),
    _cw(hfitRes._cw),
    _rw(hfitRes._rw),
    _chi2w(hfitRes._chi2w),
      
    _dfdz(hfitRes._dfdz),
    _fz0(hfitRes._fz0) 
  { }


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void HelixFitHackResult::clear() {
    _hdef.init();
    _fit.setFailure(1,"failure");
    
    _sxy.clear();
    _srphi.clear();
    _chi2   = -1.;
    _radius = -1.;
    
    _sxyw.clear();
    _rw   = -1.;
    _chi2w = -1.;
      
    _dfdz = -1.e6;
    _fz0  = -1.e6;
  }



//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  HelixFitHackResult& HelixFitHackResult::operator = (HelixFitHackResult const& other) {
    if (this != &other){
      _hdef   = other._hdef;
      _fit    = other._fit;
					// equal weights - XY and RPhi fits
      _sxy    = other._sxy;
      _srphi  = other._srphi;
      _chi2   = other._chi2;
      _center = other._center;
      _radius = other._radius;
					// non-equal weights - XY fit only
      _sxyw   = other._sxyw;
      _cw     = other._cw;
      _rw     = other._rw;
      _chi2w  = other._chi2w;

      _dfdz   = other._dfdz;
      _fz0    = other._fz0;
    }
    return *this;
  }

//-----------------------------------------------------------------------------
  void HelixFitHackResult::print(const char* Title) {

    printf(" HelixFitHackResult::print: %s\n",Title);

    printf(" _sxy (N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
	   _sxy.qn(),_sxy.x0(),_sxy.y0(),_sxy.radius(),_sxy.chi2DofCircle());

    printf(" center, radius, chi2: %8.3f %8.3f %8.3f %10.2f\n", _center.x(),_center.y(),_radius,_chi2);

    printf(" _sxyw(N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
	   _sxyw.qn(),_sxyw.x0(),_sxyw.y0(),_sxyw.radius(),_sxyw.chi2DofCircle());

    printf(" cw, rw, chi2w: %8.3f %8.3f %8.3f %10.2f\n", _cw.x(),_cw.y(),_rw,_chi2w);

    printf(" _srphi(phi0, df/dz, chi2: %3.0f %8.3f %10.5f %10.2f\n",
	   _sxyw.qn(),_srphi.phi0(),_srphi.dfdz(),_srphi.chi2DofLine());

    printf(" _dfdz, _fz0, : %10.4f %10.4f\n", _dfdz,_fz0);

  }

};
