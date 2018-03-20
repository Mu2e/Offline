//
#ifndef CalHelixFinderData_HH
#define CalHelixFinderData_HH

#include "TObject.h"
#include "CalPatRec/inc/LsqSums4.hh"
#include "CalPatRec/inc/CalHelixPoint.hh"

#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

class HelixTraj;

namespace mu2e {

  class TimeCluster;
  class Panel;

  struct ChannelID {
    int Station;
    int Plane; 
    int Face; 
    int Panel; 
    //      int Layer;
  };

  struct HitData_t {
    const StrawHit*         fHit;
    const StrawHitPosition* fPos;
    const Straw*            fStraw;
    int                     fSeedNumber;
    int                     fNSecondHits;
    int                     fDeltaIndex;
    int                     fIsOutlier;
    float                   fChi2Min;
    float                   fSigW;     // cached resolution along the wire
    float                   fRMid;
    float                   fDr;	// work variable

    HitData_t(const StrawHit* Hit, const StrawHitPosition* Pos, const Straw* aStraw, float SigW) {
      fHit         = Hit; 
      fPos         = Pos; 
      fStraw       = aStraw; 
      fChi2Min     = 1.1e10; 
      fSigW        = SigW; 
      fSeedNumber  = -1; 
      fNSecondHits = -1;
      fDeltaIndex  = -1;
      fIsOutlier   = -1;
      fRMid        = fStraw->getMidPoint().perp();
      fDr          = 1.1e10;
    }

    int Used     () { return (fChi2Min < 1.e10) ; }
    int isOutlier() { return (fIsOutlier > 0)   ; }
  };
    
  
  struct PanelZ_t {
    int                              fNHits;      // guess, total number of hits per panel
    std::vector<CalHelixPoint>       fHitData;
    const Panel*                     fPanel;      // backward pointer to the tracker panel
    double                           wx;          // direction cosines of the wires, assumed to be all the same
    double                           wy;      
    double                           phi;         // phi angle of the wire
    double                           z;           // 

    PanelZ_t    (){
      fNHits  = 0;
      fPanel  = NULL;
    }

    PanelZ_t    (const PanelZ_t&Copy){
      fNHits = Copy.fNHits;  
      fPanel = Copy.fPanel;  
      wx     = Copy.wx;      
      wy     = Copy.wy;      
      phi    = Copy.phi;     
      z      = Copy.z;       
      int nhits = Copy.fHitData.size();
      for (int i=0; i<nhits; ++i){
	fHitData.push_back(CalHelixPoint(Copy.fHitData.at(i)));
      }
    }
  }; 

  struct SeedInfo_t {
    int  Panel;	 
    int  PanelHitIndex;
    SeedInfo_t(){
      Panel         = -1;	 
      PanelHitIndex = -1;
    }
    SeedInfo_t(int P, int H){
      Panel         = P;	 
      PanelHitIndex = H;
    }
  };
  //---------------------------------------------------------------------------
  // output struct
  //-----------------------------------------------------------------------------
  class CalHelixFinderData : public TObject {
  public:
   

    enum {
      kNStations      = 20,
      kNTotalFaces    = 80,
      kNTotalPanels   = 240,
      kNFaces         =  4,
      kNPanelsPerFace =  3
    };
    
    enum { kMaxResidIndex = 500 };
    
    struct Diag_t {
      int       loopId_4;    // fData[4]
      double    radius_5;    // fData[5]
      double    phi0_6;
      double    chi2_circle;
      double    z0_6;
      double    rdfdz_7;
      double    dfdz_8;
      int       n_rescued_points_9;    // fData[9]
      double    dz_10;
      int       n_active_11;           // fData[11]
      double    chi2_dof_circle_12;
      double    chi2_dof_line_13;
      double    radius_14;
      double    chi2_dof_circle_15;
      int       n_rescued_points_16;    // fData[16]
      double    dfdzres_17;
      double    dfdzres_18;
      double    dfdzres_19;
      double    dr_20;
      double    dr_21;
      double    dphi0res_22;
      double    dphi0res_23;
      double    dphi0res_24;

      int       nStationPairs;
      
      double    dfdz;
      double    dfdz_scaled;
      double    chi2_line;
      int       n_active;
      double    resid[kMaxResidIndex];
      double    dist [kMaxResidIndex];
      double    dz   [kMaxResidIndex];

      double    dr;
      double    straw_mean_radius;
      double    chi2d_helix;

    };
    
    const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
    art::Ptr<TimeCluster>             _timeClusterPtr;

    HelixTraj*                        _helix;

    std::vector<StrawHitIndex>        _goodhits;

    SeedInfo_t                        _seedIndex;
    SeedInfo_t                        _candIndex;

    int                               _nPoints;      // n(hits)
    int                               _nFiltPoints;

    double                            _helixChi2;

    TrkParticle                       _tpart;
    TrkFitDirection                   _fdir;

    const StrawHitCollection*         _shcol;
    const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _shfcol;
    
    TrkErrCode                        _fit;	    // fit status code from last fit
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
// diagnostics, histogramming
//-----------------------------------------------------------------------------
    Diag_t             _diag;
//-----------------------------------------------------------------------------
// structure used to organize thei strawHits for the pattern recognition
//-----------------------------------------------------------------------------
//    PanelZ_t             _oTracker[kNStations][kNFaces][kNPanelsPerFace];
    PanelZ_t             _oTracker[kNTotalPanels];
    std::vector<int>     _hitsUsed[kNTotalPanels];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    CalHelixFinderData();
    ~CalHelixFinderData();

    CalHelixFinderData(const CalHelixFinderData& Data);

    // CalHelixFinderData& operator =(CalHelixFinderData const& other);

    const StrawHitCollection*         shcol () { return _shcol ; }
    const StrawHitPositionCollection* shpos () { return _shpos ; }
    const StrawHitFlagCollection*     shfcol() { return _shfcol; }

    bool          fitIsValid        () { return _sxy.qn() > 0; }
    bool          weightedFitIsValid() { return _sxyw.qn() > 0; }
    int           maxIndex          () { return kMaxResidIndex; }
    HelixTraj*    helix             () { return _helix;        }

    int           nGoodHits         () { return _goodhits.size(); }

    void          orderID           (ChannelID* X, ChannelID* O);

    void          print(const char* Title);
    void          clearTempVariables();
    void          clearResults();

    void          markHelixPoints        ();

    void          setTestHelixPoints     ();
    void          resetTestHelixPoints   ();

  };

};
#endif

