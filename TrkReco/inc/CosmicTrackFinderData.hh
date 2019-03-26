//Author: S Middleton 2018
#ifndef CosmicTrackFinderData_HH
#define CosmicTrackFinderData_HH


#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "RecoDataProducts/inc/CosmicTrackSeed.hh" //CHANGE SEED
#include "Mu2eUtilities/inc/BuildMatrixSums.hh"
#include "TrkReco/inc/TrkFaceData.hh"

#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//c++
#include <array>

//class HelixTraj;

using namespace ROOT::Math::VectorUtil;

namespace mu2e {
  typedef ROOT::Math::XYVectorF  XYVec;
  // struct for weighted positions
  class XYWVec : public XYVec {
  public :
    XYWVec(XYZVec pos, int face, float weight=1.0) : XYVec(pos.x(),pos.y()), _face(face), _weight(weight){}
    float weight() const { return _weight; }
    int   face() const { return _face; }

  private :
    int   _face;
    float _weight; // weight for this position
  };

  class TimeCluster;
 
  //  class Panel;

  //---------------------------------------------------------------------------
  // output struct
  //-----------------------------------------------------------------------------
  class CosmicTrackFinderData {
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
   
      int       nShFit; //after fit
      int       nChFit; //after fit

      int       nChPPanel;
      int       nChHits;

      float    hit_residualX[kMaxResidIndex];
      float    hit_residualY[kMaxResidIndex];
      float    hit_pullX[kMaxResidIndex];
      float    hit_pullY[kMaxResidIndex];
      
      //float    dz   [kMaxResidIndex];
      
      float    recon_eff;
      int      CosmicTrackFitCounter;
      float    chi2d_track;
      
      unsigned      niters;    
    };
   
    const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
    art::Ptr<TimeCluster>             _timeClusterPtr;
    const art::Event*                 event;
    CosmicTrackSeed                   _tseed;
    int                               _nStrawHits;      
    int                               _nComboHits;    
    int                               _nXYSh; //SH at start
    int                               _nXYCh; //CH at start
    int                               _nFiltComboHits;  //ComboHits from the TimeCluster filtering 
    int                               _nFiltStrawHits;  //StrawHits from the TimeCluster filtering 
    //const CosmicTrackSeedCollection*         _stcol;
    const ComboHitCollection*         _chcol;
    const StrawHitCollection*         _shcol;
    const TimeClusterCollection*      _tccol;
    ::BuildMatrixSums         _S;

//-----------------------------------------------------------------------------
// diagnostics, histogramming
//-----------------------------------------------------------------------------
    Diag_t             _diag;
//-----------------------------------------------------------------------------
// structure used to organize the strawHits for the pattern recognition
//-----------------------------------------------------------------------------
    std::array<FaceZ_t,StrawId::_ntotalfaces>            _oTracker;//array of faces, length of number of faces
    std::vector<ComboHit>                                _chHitsToProcess;
    std::vector<XYWVec>                                  _chHitsWPos;
    
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    CosmicTrackFinderData();
    ~CosmicTrackFinderData();

    //const CosmicTrackSeedCollection*   stcol(){return      _stcol;}
    const ComboHitCollection*         chcol () { return _chcol ; }
    const StrawHitCollection*         shcol () {return _shcol;}
    //const StrawHitFlagCollection*     shfcol() { return _shfcol; }
    const TimeClusterCollection*         tccol () { return _tccol ; }
    int           maxIndex          () { return kMaxResidIndex; }
    void          orderID           (ChannelID* X, ChannelID* O);
    void          deleteTrack ();
    void          print(const char* Title);
    void          clearTempVariables();
    void          clearResults();

  };

};
#endif

