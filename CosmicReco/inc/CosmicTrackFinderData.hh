//Author: S Middleton 
//Purpose: Store details of cosmic track
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
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh" //CHANGE SEED
#include "Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"
#include "TrkReco/inc/TrkFaceData.hh"

#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//c++
#include <array>

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
  class CosmicTrackFinderData {
  public:
    
    enum { kMaxResidIndex = 500 };

    constexpr static uint16_t        kNMaxChHits = 150;
    
    struct ChannelID {
      int Station;
      int Plane; 
      int Face; 
      int Panel; 
    };

   
    const art::Event*                 event;
    const art::Run*		      run;
    
    const TimeCluster*                _timeCluster;     // hides vector of its time cluster straw hit indices
    art::Ptr<TimeCluster>             _timeClusterPtr;
    
    CosmicTrackSeed                   _tseed;
    int                               _nStrawHits;      
    int                               _nComboHits;    
    
    int                               _nXYCh; //CH at start
    int                               _nFiltComboHits;  //ComboHits from the TimeCluster filtering 
    int                               _nFiltStrawHits;  //StrawHits from the TimeCluster filtering 
    
    const ComboHitCollection*         _chcol;
    const StrawHitCollection*         _shcol;
    const TimeClusterCollection*      _tccol;
    const StrawDigiMCCollection*      _mccol;
    ::BuildLinearFitMatrixSums         _S;//USED?

//-----------------------------------------------------------------------------
// diagnostics, histogramming
//-----------------------------------------------------------------------------
    std::array<FaceZ_t,StrawId::_ntotalfaces>         _oTracker;//array of faces, length of number of faces
    ComboHitCollection                                _chHitsToProcess;
    std::vector<XYWVec>                               _chHitsWPos;
    std::vector<StrawDigiMC>		              _mcDigisToProcess;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    CosmicTrackFinderData();
    ~CosmicTrackFinderData();

    const ComboHitCollection*         chcol () { return _chcol ; }
    const StrawHitCollection*         shcol () {return _shcol;}
    const TimeClusterCollection*         tccol () { return _tccol ; }
    int           maxIndex          () { return kMaxResidIndex; }
    void          orderID           (ChannelID* X, ChannelID* O);
    void          deleteTrack ();
    void          print(const char* Title);
    void          clearMCVariables();
    void          clearTempVariables();
    void          clearResults();
    
  };

};
#endif

