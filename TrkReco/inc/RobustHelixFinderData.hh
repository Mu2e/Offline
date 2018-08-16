//
#ifndef RobustHelixFinderData_HH
#define RobustHelixFinderData_HH

#include "Mu2eUtilities/inc/LsqSums4.hh"
// #include "CalPatRec/inc/CalHelixPoint.hh"

#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"

#include "TrkReco/inc/TrkFaceData.hh"

#include <array>

class HelixTraj;

namespace mu2e {

  class TimeCluster;
  class FaceZ_t;
  //  class Panel;

  //---------------------------------------------------------------------------
  // output struct
  //-----------------------------------------------------------------------------
  class RobustHelixFinderData {
  public:
    
    enum { kMaxResidIndex = 500 };
    
    struct ChannelID {
      int Station;
      int Plane; 
      int Face; 
      int Panel; 
      //      int Layer;
    };

    struct Diag_t {
      
      int       nChPPanel;

      double    resid[kMaxResidIndex];
      double    dist [kMaxResidIndex];
      double    dz   [kMaxResidIndex];
      
      int       circleFitCounter;
      int       nrescuedhits;

      double    dr;
      double    chi2d_helix;
      
      double    chi2dXY;

      int       ntriple_0;    //number of triplets used in the first call of RobustHelix::fitCircle
      double    radius_0;     //radius resulting from the first call of RobustHelix::fitCircle
  
      int       nshsxy_0;
      double    rsxy_0;
      double    chi2dsxy_0;

      int       nshsxy_1;
      double    rsxy_1;
      double    chi2dsxy_1;

      int       nfz0counter;

      int       nshszphi;
      double    lambdaszphi;
      double    chi2dszphi;

      int       nshszphi_0;
      double    lambdaszphi_0;
      double    chi2dszphi_0;

      int       nshszphi_1;
      double    lambdaszphi_1;
      double    chi2dszphi_1;


      int       ntriple_1;    //number of triplets used in the first call of RobustHelix::fitCircle
      double    radius_1;     //radius resulting from the first call of RobustHelix::fitCircle
      
      int       ntriple_2;
      double    radius_2;

      double    lambda_0;
      double    lambda_1;

      int       xyniter;
      int       fzniter;
      int       niter;

    };
    
    const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
    art::Ptr<TimeCluster>             _timeClusterPtr;

    //    HelixTraj*                        _helix;

    HelixSeed                         _hseed;

    std::vector<StrawHitIndex>        _goodhits;

    // SeedInfo_t                        _seedIndex;
    // SeedInfo_t                        _candIndex;

    int                               _nStrawHits;      
    int                               _nComboHits;    

    int                               _nXYSh;
    int                               _nZPhiSh;
  
    int                               _nFiltComboHits;  //ComboHits from the TimeCluster + DeltaFinder filtering 
    int                               _nFiltStrawHits;  //StrawHits from the TimeCluster + DeltaFinder filtering 

    double                            _helixChi2;

    TrkParticle                       _tpart;
    TrkFitDirection                   _fdir;

    const ComboHitCollection*         _chcol;
    // const StrawHitPositionCollection* _shpos;
    const StrawHitFlagCollection*     _chfcol;
    
    TrkErrCode                        _fit;	    // fit status code from last fit
//-----------------------------------------------------------------------------
// circle parameters; the z center is ignored.
//-----------------------------------------------------------------------------
    ::LsqSums4         _sxy;
    ::LsqSums4         _szphi;

    XYZVec             _center;
    double             _radius;

    double             _chi2;
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
    std::array<FaceZ_t,StrawId::_ntotalfaces>            _oTracker;
    // std::array<int,kNTotalPanels*kNMaxHitsPerPanel>     _hitsUsed;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    RobustHelixFinderData();
    ~RobustHelixFinderData();

    // RobustHelixFinderData(const RobustHelixFinderData& Data);

    // RobustHelixFinderData& operator =(RobustHelixFinderData const& other);

    const ComboHitCollection*         chcol () { return _chcol ; }
    // const StrawHitPositionCollection* shpos () { return _shpos ; }
    const StrawHitFlagCollection*     chfcol() { return _chfcol; }

    // bool          fitIsValid        () { return _sxy.qn() > 0; }
    // bool          weightedFitIsValid() { return _sxy.qn() > 0; }
    int           maxIndex          () { return kMaxResidIndex; }
    // HelixTraj*    helix             () { return _helix;        }

    int           nGoodHits         () { return _goodhits.size(); }

    void          orderID           (ChannelID* X, ChannelID* O);

    void          print(const char* Title);
    void          clearTempVariables();
    void          clearResults();

  };

};
#endif

