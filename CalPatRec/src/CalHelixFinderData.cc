//

#include "RecoDataProducts/inc/TimeCluster.hh"
#include "CalPatRec/inc/CalHelixFinderData.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

//-----------------------------------------------------------------------------
// CalHelixFinderData
//-----------------------------------------------------------------------------
  CalHelixFinderData::CalHelixFinderData() {
    _helix = NULL;
  }

//-----------------------------------------------------------------------------
  CalHelixFinderData::CalHelixFinderData(const CalHelixFinderData& Data) {
    _timeCluster    = Data._timeCluster;
    _timeClusterPtr = Data._timeClusterPtr;
					// the only pointer owned 
    if (Data._helix) _helix = Data._helix->clone();
    else             _helix = NULL;

    _goodhits       = Data._goodhits;
    _tpart          = Data._tpart;
    _fdir           = Data._fdir;
    _shcol          = Data._shcol;
    _shpos          = Data._shpos;
    _shfcol         = Data._shfcol;
    _fit            = Data._fit;
    _sxy            = Data._sxy;
    _srphi          = Data._srphi;
    _center         = Data._center;
    _radius         = Data._radius;
    _chi2           = Data._chi2;
    _sxyw           = Data._sxyw;
    _cw             = Data._cw;
    _rw             = Data._rw;
    _chi2w          = Data._chi2w;
    _dfdz           = Data._dfdz;
    _fz0            = Data._fz0;
    _diag           = Data._diag;
  }

//-----------------------------------------------------------------------------
  CalHelixFinderData::~CalHelixFinderData() {
    if (_helix) delete _helix;
  }

//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
  void CalHelixFinderData::clearTempVariables() {

    _timeCluster    = NULL;
    _timeClusterPtr = art::Ptr<TimeCluster>();

    _goodhits.clear();
    
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
  void CalHelixFinderData::print(const char* Title) {

    printf(" CalHelixFinderData::print: %s\n",Title);

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
