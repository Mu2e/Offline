//
#ifndef RobustHelixFinderData_HH
#define RobustHelixFinderData_HH

#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
// #include "CalPatRec/inc/CalHelixPoint.hh"

#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"

#include "Offline/TrkReco/inc/TrkFaceData.hh"

#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//c++
#include <array>

class HelixTraj;

using namespace ROOT::Math::VectorUtil;

namespace mu2e {
  typedef ROOT::Math::XYVectorF  XYVec;
  // struct for weighted positions
  class XYWVec : public XYVec {
    public :
      XYWVec(XYZVectorF pos, int face, float weight=1.0) : XYVec(pos.x(),pos.y()), _face(face), _weight(weight){}
      float weight() const { return _weight; }
      int   face() const { return _face; }

    private :
      int   _face;
      float _weight; // weight for this position
  };

  class TimeCluster;
  class FaceZ_t;
  //  class Panel;

  //---------------------------------------------------------------------------
  // output struct
  //-----------------------------------------------------------------------------
  class RobustHelixFinderData {
    public:

      enum { kMaxResidIndex = 500 };

      constexpr static uint16_t        kNMaxChHits = 150;

      struct ChannelID {
        int Station;
        int Plane;
        int Face;
        int Panel;
        //      int Layer;
      };

      struct Diag_t {

        int       nShFitCircle;
        int       nChFitCircle;

        int       nShFitXY;
        int       nChFitXY;

        int       nChPPanel;
        int       nChHits;

        float    resid[kMaxResidIndex];
        float    dist [kMaxResidIndex];
        float    rwdot[kMaxResidIndex];
        float    dz   [kMaxResidIndex];

        int      nXYCh;
        int      nZPhiCh;

        int       circleFitCounter;
        int       nrescuedhits;

        float    dr;
        float    chi2d_helix;

        float    chi2dXY;
        float    chi2dZPhi;

        int       ntriple_0;    //number of triplets used in the first call of RobustHelix::fitCircle
        float    radius_0;     //radius resulting from the first call of RobustHelix::fitCircle

        int       nshsxy_0;
        float    rsxy_0;
        float    chi2dsxy_0;

        int       nshsxy_1;
        float    rsxy_1;
        float    chi2dsxy_1;

        int       nfz0counter;

        int       nshszphi;
        float    lambdaszphi;
        float    chi2dszphi;

        int       nshszphi_0;
        float    lambdaszphi_0;
        float    chi2dszphi_0;

        int       nshszphi_1;
        float    lambdaszphi_1;
        float    chi2dszphi_1;


        int       ntriple_1;    //number of triplets used in the first call of RobustHelix::fitCircle
        float    radius_1;     //radius resulting from the first call of RobustHelix::fitCircle

        int       ntriple_2;
        float    radius_2;

        float    lambda_0;
        float    lambda_1;

        int       xyniter;
        int       fzniter;
        int       niter;

        int       nLoops;
        float     meanHitRadialDist;
        int       nHitsLoopFailed;


      };

      const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
      art::Ptr<TimeCluster>             _timeClusterPtr;

      //    HelixTraj*                        _helix;

      HelixSeed                         _hseed;

      //    std::vector<StrawHitIndex>        _goodhits;

      // SeedInfo_t                        _seedIndex;
      // SeedInfo_t                        _candIndex;

      int                               _nStrawHits;
      int                               _nComboHits;

      int                               _nXYSh;
      int                               _nZPhiSh;

      int                               _nXYCh;
      int                               _nZPhiCh;

      int                               _nFiltComboHits;  //ComboHits from the TimeCluster + DeltaFinder filtering
      int                               _nFiltStrawHits;  //StrawHits from the TimeCluster + DeltaFinder filtering

      // double                            _helixChi2;

      // TrkParticle                       _tpart;
      // TrkFitDirection                   _fdir;

      const ComboHitCollection*         _chcol;

      //    TrkErrCode                        _fit;     // fit status code from last fit
      //-----------------------------------------------------------------------------
      // circle parameters; the z center is ignored.
      //-----------------------------------------------------------------------------
      ::LsqSums4         _sxy;
      ::LsqSums4         _szphi;

      // XYZVectorF             _center;
      // double             _radius;

      // double             _chi2;
      //-----------------------------------------------------------------------------
      // Z parameters; dfdz is the slope of phi vs z (=-sign(1.0,qBzdir)/(R*tandip)),
      // fz0 is the phi value of the particle where it goes through z=0
      // note that dfdz has a physical ambiguity in q*zdir.
      //-----------------------------------------------------------------------------
      // double             _dfdz;
      // double             _fz0;
      //-----------------------------------------------------------------------------
      // diagnostics, histogramming
      //-----------------------------------------------------------------------------
      Diag_t             _diag;
      //-----------------------------------------------------------------------------
      // structure used to organize thei strawHits for the pattern recognition
      //-----------------------------------------------------------------------------
      std::array<FaceZ_t,StrawId::_ntotalfaces>            _oTracker;

      std::vector<ComboHit>                                _chHitsToProcess;
      std::vector<XYWVec>                                  _chHitsWPos;
      // std::array<int,kNTotalPanels*kNMaxHitsPerPanel>     _hitsUsed;
      //-----------------------------------------------------------------------------
      // functions
      //-----------------------------------------------------------------------------
      RobustHelixFinderData();

      const ComboHitCollection*         chcol () { return _chcol ; }

      // bool          fitIsValid        () { return _sxy.qn() > 0; }
      // bool          weightedFitIsValid() { return _sxy.qn() > 0; }
      int           maxIndex          () { return kMaxResidIndex; }
      // HelixTraj*    helix             () { return _helix;        }

      // int           nGoodHits         () { return _goodhits.size(); }

      void          orderID           (ChannelID* X, ChannelID* O);

      void          print(const char* Title);
      void          clearTempVariables();
      void          clearResults();

  };

}
#endif

