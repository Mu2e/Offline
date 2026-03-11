//

#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrkReco/inc/RobustHelixFinderData.hh"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {

  //-----------------------------------------------------------------------------
  // RobustHelixFinderData
  //-----------------------------------------------------------------------------
  RobustHelixFinderData::RobustHelixFinderData() {
    _chHitsToProcess.reserve(kNMaxChHits);
    _chHitsWPos     .reserve(kNMaxChHits);
  }


  //-----------------------------------------------------------------------------
  void RobustHelixFinderData::orderID(RobustHelixFinderData::ChannelID* X, RobustHelixFinderData::ChannelID* O) {
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
  void RobustHelixFinderData::clearTempVariables() {

    _timeCluster    = NULL;
    _timeClusterPtr = art::Ptr<TimeCluster>();

    _chHitsToProcess.clear();
    _chHitsWPos.clear();

    _nStrawHits = 0;
    _nComboHits = 0;

    _sxy.clear();
    _szphi.clear();
    _nFiltComboHits = 0;
    _nFiltStrawHits = 0;

    _nXYSh       = 0;
    _nZPhiSh     = 0;

    _nStrawHits  = 0;
    _nComboHits  = 0;

    _nFiltComboHits = 0;
    _nFiltStrawHits = 0;

    //clear the panel-based structure
    for (int f=0; f<StrawId::_ntotalfaces; ++f) {
      FaceZ_t* facez = &_oTracker[f];
      facez->bestFaceHit = -1;
      facez->idChBegin   = -1;
      facez->idChEnd     = -1;
      for (int p=0; p<FaceZ_t::kNPanels; ++p){
        PanelZ_t* panelz  = &facez->panelZs[p];
        panelz->idChBegin = -1;
        panelz->idChEnd   = -1;
      }
    }
  }

  //-----------------------------------------------------------------------------
  // don't clear the diagnostics part.
  //-----------------------------------------------------------------------------
  void RobustHelixFinderData::clearResults() {

    // _goodhits.clear();

    _sxy.clear();
    _szphi.clear();

    _nXYSh       = 0;
    _nZPhiSh     = 0;

    _nStrawHits  = 0;
    _nComboHits  = 0;
  }


  //-----------------------------------------------------------------------------
  void RobustHelixFinderData::print(const char* Title) {

    printf(" RobustHelixFinderData::print: %s\n",Title);

    // printf(" _sxy (N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
    //     _sxy.qn(),_sxy.x0(),_sxy.y0(),_sxy.radius(),_sxy.chi2DofCircle());

    printf(" center, radius, chi2: %8.3f %8.3f %8.3f %10.2f\n", _hseed._helix.centerx(),
        _hseed._helix.centery(),_hseed._helix._radius,_sxy.chi2DofCircle());

    // printf(" _sxyw(N, X0, Y0, R, chi2: %3.0f %8.3f %8.3f %8.3f %10.2f)\n",
    //     _sxyw.qn(),_sxyw.x0(),_sxyw.y0(),_sxyw.radius(),_sxyw.chi2DofCircle());

    // printf(" cw, rw, chi2w: %8.3f %8.3f %8.3f %10.2f\n", _cw.x(),_cw.y(),_rw,_chi2w);

    // printf(" _szphi(phi0, df/dz, chi2: %3.0f %8.3f %10.5f %10.2f\n",
    //     _sxy.qn(),_szphi.phi0(),_szphi.dfdz(),_szphi.chi2DofLine());

    printf(" _dfdz, _fz0, : %10.4f %10.4f\n", 1./_hseed._helix._lambda,_hseed._helix._fz0);

  }

}
