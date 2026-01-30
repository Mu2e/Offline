//

#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/CalPatRec/inc/CalHelixFinderData.hh"
#include "Offline/BTrkLegacy/inc/HelixParams.hh"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

//-----------------------------------------------------------------------------
// CalHelixFinderData
//-----------------------------------------------------------------------------
  CalHelixFinderData::CalHelixFinderData() {
    _helix = NULL;
    _goodhits.reserve(kNMaxChHits);
    _chHitsToProcess. reserve(kNMaxChHits);
  }

//-----------------------------------------------------------------------------
  // CalHelixFinderData::CalHelixFinderData(const CalHelixFinderData& Data) {
  //   _timeCluster    = Data._timeCluster;
  //   _timeClusterPtr = Data._timeClusterPtr;
  //                                         // the only pointer owned
  //   if (Data._helix) _helix = Data._helix->clone();
  //   else             _helix = NULL;
  //   _helicity       = Data._helicity;

  //   _chHitsToProcess= Data._chHitsToProcess;
  //   _zFace          = Data._zFace;
  //   _phiPanel       = Data._phiPanel;

  //   _goodhits       = Data._goodhits;
  //   _tpart          = Data._tpart;
  //   _fdir           = Data._fdir;
  //   _chcol          = Data._chcol;

  //   _shfcol         = Data._shfcol;
  //   _fit            = Data._fit;
  //   _sxy            = Data._sxy;
  //   _szphi          = Data._szphi;
  //   _center         = Data._center;
  //   _radius         = Data._radius;

  //   _dfdz           = Data._dfdz;
  //   _fz0            = Data._fz0;
  //   _diag           = Data._diag;


  //   _nXYSh          = Data._nXYSh;
  //   _nZPhiSh        = Data._nZPhiSh;
  //   _nStrawHits     = Data._nStrawHits;
  //   _nComboHits     = Data._nComboHits;
  //   _nFiltPoints    = Data._nFiltPoints;
  //   _nFiltStrawHits = Data._nFiltStrawHits;

  //   _helixChi2      = Data._helixChi2;

  //   _seedIndex      = Data._seedIndex;
  //   _candIndex      = Data._candIndex;

  //   //copy the info relative to the  panels
  //   _oTracker       = Data._oTracker;

  //   _hitsUsed       = Data._hitsUsed;
  // }

//-----------------------------------------------------------------------------
  CalHelixFinderData::~CalHelixFinderData() {
    if (_helix) delete _helix;
  }

//-----------------------------------------------------------------------------
  void CalHelixFinderData::orderID(ChannelID* X, ChannelID* O) {
    if (X->Panel % 2 == 0) X->Face = 0;
    else                   X->Face = 1; // define original face

    O->Station = X->Station; // stations already ordered
    O->Plane   = X->Plane;   // planes already ordered, but not necessary for ordered construct

    if (X->Station % 2 == 0) {
      if (X->Plane == 0) O->Face = 1 - X->Face;
      else               O->Face = X->Face + 2;
    }
    else {
      if (X->Plane == 0) O->Face = X->Face;
      else               O->Face = 3 - X->Face; // order face
    }

    O->Panel = int(X->Panel/2);                // order panel

    // int n = X->Station + X->Plane + X->Face;   // pattern has no intrinsic meaning, just works
    // if (n % 2 == 0) O->Layer = 1 - X->Layer;
    // else            O->Layer = X->Layer;       // order layer
  }

//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
  void CalHelixFinderData::clearTempVariables() {

    _timeCluster    = NULL;
    _timeClusterPtr = art::Ptr<TimeCluster>();

    _chHitsToProcess.clear();

    _goodhits.clear();

    _fit.setFailure(1,"failure");

    _sxy.clear();
    _szphi.clear();

    _radius = -1.;

    _dfdz = -1.e6;
    _fz0  = -1.e6;

    _nFiltPoints    = 0;
    _nFiltStrawHits = 0;

    _nXYSh       = 0;
    _nZPhiSh     = 0;

    _nStrawHits  = 0;
    _nComboHits  = 0;

    _helixChi2   = 1e10;

    _seedIndex   = HitInfo_t();
    _candIndex   = HitInfo_t();

    //clear the panel-based structure
    for (int f=0; f<StrawId::_ntotalfaces; ++f) {
      FaceZ_t*  facez  = &_oTracker[f];
      facez->bestFaceHit = -1;
      facez->idChBegin   = -1;
      facez->idChEnd     = -1;

      for (int p=0; p<FaceZ_t::kNPanels; ++p) {
        PanelZ_t* panelz = &facez->panelZs[p];
        panelz->idChBegin = -1;
        panelz->idChEnd   = -1;
      }
    }

    _hitsUsed =  {0};
  }

 void CalHelixFinderData::clearTimeClusterInfo() {

    _timeCluster    = NULL;
    _timeClusterPtr = art::Ptr<TimeCluster>();

    _chHitsToProcess.clear();

    _nFiltPoints    = 0;
    _nFiltStrawHits = 0;

    //clear the panel-based structure
    for (int f=0; f<StrawId::_ntotalfaces; ++f) {
      FaceZ_t*  facez  = &_oTracker[f];
      facez->bestFaceHit = -1;
      facez->idChBegin   = -1;
      facez->idChEnd     = -1;

      for (int p=0; p<FaceZ_t::kNPanels; ++p) {
        PanelZ_t* panelz = &facez->panelZs[p];
        panelz->idChBegin = -1;
        panelz->idChEnd   = -1;
      }
    }

    _hitsUsed =  {0};
  }

void CalHelixFinderData::clearHelixInfo() {

    _goodhits.clear();

    _fit.setFailure(1,"failure");

    _sxy.clear();
    _szphi.clear();

    _radius = -1.;

    _dfdz = -1.e6;
    _fz0  = -1.e6;

    _nXYSh       = 0;
    _nZPhiSh     = 0;

    _nStrawHits  = 0;
    _nComboHits  = 0;

    _helixChi2   = 1e10;

    _seedIndex   = HitInfo_t();
    _candIndex   = HitInfo_t();

  }
//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
   void CalHelixFinderData::clearResults() {

    _goodhits.clear();

    _chHitsToProcess.clear();

    _fit.setFailure(1,"failure");

    _sxy.clear();
    _szphi.clear();
    //    _chi2   = -1.;
    _radius = -1.;

    _dfdz = -1.e6;
    _fz0  = -1.e6;


    _nXYSh       = 0;
    _nZPhiSh     = 0;

    _nStrawHits  = 0;
    _nComboHits  = 0;

    _helixChi2   = 1e10;

    // _seedIndex   = SeedInfo_t(-1,-1);
    // _candIndex   = SeedInfo_t(-1,-1);

    _seedIndex   = HitInfo_t();
    _candIndex   = HitInfo_t();

    _hitsUsed    = {0};
  }


//-----------------------------------------------------------------------------
  void CalHelixFinderData::print(const char* Title) {

    printf(" CalHelixFinderData::print: %s\n",Title);

    printf(" _sxy (N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
           _sxy.qn(),_sxy.x0(),_sxy.y0(),_sxy.radius(),_sxy.chi2DofCircle());

    printf(" center, radius, chi2: %8.3f %8.3f %8.3f %10.2f\n", _center.x(),_center.y(),_radius,_sxy.chi2DofCircle());

    // printf(" _sxyw(N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
    //            _sxyw.qn(),_sxyw.x0(),_sxyw.y0(),_sxyw.radius(),_sxyw.chi2DofCircle());

    // printf(" cw, rw, chi2w: %8.3f %8.3f %8.3f %10.2f\n", _cw.x(),_cw.y(),_rw,_chi2w);

    printf(" _szphi(phi0, df/dz, chi2: %3.0f %8.3f %10.5f %10.2f\n",
           _sxy.qn(),_szphi.phi0(),_szphi.dfdz(),_szphi.chi2DofLine());

    printf(" _dfdz, _fz0, : %10.4f %10.4f\n", _dfdz,_fz0);

  }

}
