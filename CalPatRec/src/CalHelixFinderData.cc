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
    _chcol          = Data._chcol;
    // _shpos          = Data._shpos;
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
    
    _nPoints        = Data._nPoints;
    _nFiltPoints    = Data._nFiltPoints;
    _helixChi2      = Data._helixChi2;

    _seedIndex      = Data._seedIndex;
    _candIndex      = Data._candIndex;

    //copy the info relative to the  panels
    for (int p=0; p<kNTotalPanels; ++p){
      _oTracker[p] = Data._oTracker[p];//PanelZ_t(Data._oTracker[p]);
    }    
    
    _hitsUsed       = Data._hitsUsed;
  }

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
// set the Helix points as tested according to the last check made
//-----------------------------------------------------------------------------
  void CalHelixFinderData::setTestHelixPoints(){
    // PanelZ_t*   panel(0);
    // for (int p=0; p<kNTotalPanels; ++p){
    //   panel  = &_oTracker[p];
    //   int nhits = panel->fHitData.size();
    //   for (int i=0; i<nhits; ++i){
    // 	panel->fHitData.at(i)._tested = panel->fHitData.at(i)._used;
    //   }
    // }
  }
//-----------------------------------------------------------------------------
// set the Helix points as used according to the last test made
//-----------------------------------------------------------------------------
  void CalHelixFinderData::markHelixPoints(){
    // PanelZ_t*   panel(0);
    // int         counter(0);
    // for (int p=0; p<kNTotalPanels; ++p){
    //   panel  = &_oTracker[p];
    //   int nhits = panel->fHitData.size();
    //   for (int i=0; i<nhits; ++i){
    // 	panel->fHitData.at(i)._used       = panel->fHitData.at(i)._tested;
    // 	panel->fHitData.at(i)._dzFromSeed = panel->fHitData.at(i)._testDzFromSeed;
    // 	panel->fHitData.at(i)._drFromPred = panel->fHitData.at(i)._testDrFromPred;
    // 	if (panel->fHitData.at(i)._used == 1) ++counter;
    //   }
    // }
    // //update the value of the number of good points associated with the helix candidate
    // _nPoints = counter;
  }
//-----------------------------------------------------------------------------
// reset the "tested" flag of the Helix points
//-----------------------------------------------------------------------------
  void CalHelixFinderData::resetTestHelixPoints(){
    // PanelZ_t*   panel(0);
    // for (int p=0; p<kNTotalPanels; ++p){
    //   panel  = &_oTracker[p];
    //   int nhits = panel->fHitData.size();
    //   for (int i=0; i<nhits; ++i){
    // 	panel->fHitData.at(i)._tested     = 0;
    // 	panel->fHitData.at(i)._dzFromSeed = 0;
    // 	panel->fHitData.at(i)._drFromPred = 0;
    //   }
    // }
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

    _nFiltPoints = 0;
    _nPoints     = 0;
    _helixChi2   = 1e10;

    _seedIndex   = SeedInfo_t(-1,-1);
    _candIndex   = SeedInfo_t(-1,-1);

    //clear the panel-based structure
    for (int p=0; p<kNTotalPanels; ++p) {
      PanelZ_t* panelz = &_oTracker[p];
      panelz->fNHits = 0;
      panelz->fHitData.clear() ;
    }

    _hitsUsed =  {0};
  }

//-----------------------------------------------------------------------------
// don't clear the diagnostics part.
//-----------------------------------------------------------------------------
   void CalHelixFinderData::clearResults() {

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

    _nPoints     = 0;
    _helixChi2   = 1e10;

    _seedIndex   = SeedInfo_t(-1,-1);
    _candIndex   = SeedInfo_t(-1,-1);
    
    _hitsUsed    = {0};
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
