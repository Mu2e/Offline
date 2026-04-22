//
#ifndef CalHelixFinderData_HH
#define CalHelixFinderData_HH

#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"

#include "Offline/BTrkLegacy/inc/TrkErrCode.hh"
#include "Offline/BTrkLegacy/inc/TrkParticle.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/TrkReco/inc/TrkFaceData.hh"
#include "Offline/BTrkLegacy/inc/HelixParams.hh"

#include <array>


namespace mu2e {

  class TimeCluster;
  //  class Panel;

  struct ChannelID {
    int Station;
    int Plane;
    int Face;
    int Panel;
    //      int Layer;
  };

  //---------------------------------------------------------------------------
  // output struct
  //-----------------------------------------------------------------------------
  class CalHelixFinderData {
  public:


    // enum {
    //   // kNStations      = 20,
    //   // kNTotalFaces    = 80,
    //   // kNTotalPanels   = 240,
    //   kNMaxHitsPerPanel= 20
    //   // kNFaces         =  4,
    //   // kNPanelsPerFace =  3
    // };

    enum { kMaxResidIndex = 5000 };

    constexpr static uint16_t        kNMaxChHits = 5000;

    struct Diag_t {
      int       loopId_4;    // fData[4]
      float    radius_5;    // fData[5]
      float    phi0_6;
      float    chi2_circle;
      float    z0_6;
      float    rdfdz_7;
      float    dfdz_8;
      int       n_rescued_points_9;    // fData[9]
      float    dz_10;
      int       n_active_11;           // fData[11]
      float    chi2_dof_circle_12;
      float    chi2_dof_line_13;
      float    radius_14;
      float    chi2_dof_circle_15;
      int       n_rescued_points_16;    // fData[16]
      float    dfdzres_17;
      float    dfdzres_18;
      float    dfdzres_19;
      float    dr_20;
      float    dr_21;
      float    dphi0res_22;
      float    dphi0res_23;
      float    dphi0res_24;

      int       nStationPairs;

      float    dfdz;
      float    dfdz_scaled;
      float    chi2_line;
      int       n_active;
      float    resid[kMaxResidIndex];
      float    dist [kMaxResidIndex];
      float    dz   [kMaxResidIndex];

      float    dr;
      float    straw_mean_radius;
      float    chi2d_helix;

      int      nLoops;
      float    meanHitRadialDist;

      float    nHitsRatio;
      float    eDepAvg;
    };

    const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
    art::Ptr<TimeCluster>             _timeClusterPtr;

    HelixTraj*                        _helix;
    Helicity                          _helicity;

    std::vector<int>                  _goodhits;

    HitInfo_t                         _seedIndex;
    HitInfo_t                         _candIndex;

    int                               _nStrawHits;
    int                               _nComboHits;

    int                               _nXYSh;
    int                               _nZPhiSh;

    int                               _nFiltPoints;     //ComboHits from the TimeCluster + DeltaFinder filtering
    int                               _nFiltStrawHits;  //StrawHits from the TimeCluster + DeltaFinder filtering

    float                            _helixChi2;

    TrkParticle                       _tpart;
    TrkFitDirection                   _fdir;

    const ComboHitCollection*         _chcol;
    // const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _shfcol;

    TrkErrCode                        _fit;            // fit status code from last fit
//-----------------------------------------------------------------------------
// circle parameters; the z center is ignored.
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxy;
    ::LsqSums4         _szphi;

    XYZVectorF             _center;
    float             _radius;

    //    float             _chi2;
//-----------------------------------------------------------------------------
// 2015-02-06 P.Murat: fit with non-equal weights - XY-only
//-----------------------------------------------------------------------------
    // ::LsqSums4         _sxyw;
    // XYZVectorF             _cw;
    // float             _rw;
    // float             _chi2w;
//-----------------------------------------------------------------------------
// Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)),
// fz0 is the phi value of the particle where it goes through z=0
// note that dfdz has a physical ambiguity in q*zdir.
//-----------------------------------------------------------------------------
    float             _dfdz;
    float             _fz0;
//-----------------------------------------------------------------------------
// diagnostics, histogramming
//-----------------------------------------------------------------------------
    Diag_t             _diag;
//-----------------------------------------------------------------------------
// structure used to organize thei strawHits for the pattern recognition
//-----------------------------------------------------------------------------
//    PanelZ_t                                           _oTracker[kNTotalPanels];
//    std::array<int,kNTotalPanels*kNMaxHitsPerPanel>     _hitsUsed;
    std::array<FaceZ_t,StrawId::_ntotalfaces>             _oTracker;
    std::array<float,StrawId::_ntotalfaces>               _zFace;
    std::array<float,StrawId::_nupanels>                  _phiPanel;

    std::vector<ComboHit>                                 _chHitsToProcess;
    std::array<int,kNMaxChHits>                           _hitsUsed;

    //    std::array<int,StrawId::_nupanels*PanelZ_t::kNMaxPanelHits>  _hitsUsed;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    CalHelixFinderData();
    ~CalHelixFinderData();
    CalHelixFinderData(const CalHelixFinderData&) = default;
    CalHelixFinderData& operator=(const CalHelixFinderData &) = default;

    const ComboHitCollection*         chcol () { return _chcol ; }
    // const StrawHitPositionCollection* shpos () { return _shpos ; }
    const StrawHitFlagCollection*     shfcol() { return _shfcol; }

    bool          fitIsValid        () { return _sxy.qn() > 0; }
    bool          weightedFitIsValid() { return _sxy.qn() > 0; }
    int           maxIndex          () { return kMaxResidIndex; }
    HelixTraj*    helix             () { return _helix;        }

    int           nGoodHits         () { return _goodhits.size(); }

    void          orderID           (ChannelID* X, ChannelID* O);

    void          print(const char* Title);
    void          clearTimeClusterInfo();
    void          clearHelixInfo();
    void          clearTempVariables();
    void          clearResults();

  };

}
#endif

