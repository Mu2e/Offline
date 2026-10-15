#ifndef STMGeom_STM_hh
#define STMGeom_STM_hh

// Stopping Target Monitor (STM) Object
//
// Author: Anthony Palladino
// Update: Haichuan Cao August 2023
// Update: Yongyi Wu Sept-Oct 2026, hand-stacked downstream shielding


#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/STMGeom/inc/STMDownstreamEnvelope.hh"
#include "Offline/STMGeom/inc/PermanentMagnet.hh"
#include "Offline/STMGeom/inc/TransportPipe.hh"
#include "Offline/STMGeom/inc/SupportTable.hh"
#include "Offline/STMGeom/inc/STMCollimator.hh"
#include "Offline/STMGeom/inc/GeDetector.hh"
#include "Offline/STMGeom/inc/ShieldPipe.hh"
#include "Offline/STMGeom/inc/STM_SSC.hh"
#include "Offline/STMGeom/inc/SSCSupport.hh"
#include "Offline/STMGeom/inc/LeadBrick.hh"
#include "Offline/STMGeom/inc/BrickWall.hh"
#include "Offline/STMGeom/inc/SSCFrontShield.hh"
#include "Offline/STMGeom/inc/FrontShieldingRight.hh"
#include "Offline/STMGeom/inc/FrontShieldingLeft.hh"
#include "Offline/STMGeom/inc/RightWall.hh"
#include "Offline/STMGeom/inc/LeftWall.hh"
#include "Offline/STMGeom/inc/TopWall.hh"
#include "Offline/STMGeom/inc/BottomWall.hh"
#include "Offline/STMGeom/inc/HPGeDetector.hh"
#include "Offline/STMGeom/inc/LaBrDetector.hh"
#include "Offline/STMGeom/inc/FrontShielding.hh"
#include "Offline/STMGeom/inc/LeftShielding.hh"
#include "Offline/STMGeom/inc/RightShielding.hh"
#include "Offline/STMGeom/inc/TopShielding.hh"
#include "Offline/STMGeom/inc/BottomShielding.hh"
#include "Offline/STMGeom/inc/InnerShielding.hh"
#include "Offline/STMGeom/inc/BackShielding.hh"
#include "Offline/STMGeom/inc/ElectronicShielding.hh"
#include "Offline/STMGeom/inc/STM_Absorber.hh"


namespace mu2e {

  class STMMaker;

  class STM : virtual public Detector {

  public:

    ~STM() override = default;

    // delete automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    STM( STM const& ) = delete;
    STM( STM const&& ) = delete;
    STM& operator= ( STM const& ) = delete;
    STM& operator= ( STM&&      ) = delete;


    STMDownstreamEnvelope  const * getSTMDnStrEnvPtr()       const { return _pSTMDnStrEnvParams.get(); }
    PermanentMagnet  const * getSTMMagnetPtr()               const { return _pSTMMagnetParams.get(); }
    TransportPipe    const * getSTMTransportPipePtr()        const { return _pSTMTransportPipeParams.get(); }
    SupportTable     const * getSTMMagnetSupportTablePtr()   const { return _pSTMMagnetSupportTableParams.get(); }
    STMCollimator    const * getSTMFOVCollimatorPtr()        const { return _pSTMFOVCollimatorParams.get(); }
    STMCollimator    const * getSTMSSCollimatorPtr()         const { return _pSTMSSCollimatorParams.get(); }
    SupportTable     const * getSTMDetectorSupportTablePtr() const { return _pSTMDetectorSupportTableParams.get(); }
    GeDetector       const * getSTMDetector1Ptr()            const { return _pSTMDetector1Params.get(); }
    GeDetector       const * getSTMDetector2Ptr()            const { return _pSTMDetector2Params.get(); }
    ShieldPipe       const * getSTMShieldPipePtr()           const { return _pSTMShieldPipeParams.get(); }
    STM_SSC          const * getSTM_SSCPtr()                 const { return _pSTM_SSCParams.get(); }
    SSCSupport       const * getSSCSupportPtr()              const { return _pSSCSupportParams.get(); }
    // The standard lead bricks, shared by every structure that stacks them.
    LeadBrick        const * getLeadBrickPtr()               const { return _pLeadBrickParams.get(); }
    SSCFrontShield   const * getSSCFrontShieldPtr()          const { return _pSSCFrontShieldParams.get(); }
    FrontShieldingRight const * getFrontShieldingRightPtr()  const { return _pFrontShieldingRightParams.get(); }
    FrontShieldingLeft  const * getFrontShieldingLeftPtr()   const { return _pFrontShieldingLeftParams.get(); }
    RightWall           const * getRightWallPtr()            const { return _pRightWallParams.get(); }
    LeftWall            const * getLeftWallPtr()             const { return _pLeftWallParams.get(); }
    TopWall             const * getTopWallPtr()              const { return _pTopWallParams.get(); }
    BottomWall          const * getBottomWallPtr()           const { return _pBottomWallParams.get(); }

    HPGeDetector     const * getHPGeDetectorPtr()            const { return _pSTMHPGeDetectorParams.get(); }
    LaBrDetector     const * getLaBrDetectorPtr()            const { return _pSTMLaBrDetectorParams.get(); }

    FrontShielding   const * getFrontShieldingPtr()        const { return _pSTMFrontShieldingParams.get(); }
    LeftShielding    const * getLeftShieldingPtr()         const { return _pSTMLeftShieldingParams.get(); }
    RightShielding   const * getRightShieldingPtr()        const { return _pSTMRightShieldingParams.get(); }
    TopShielding     const * getTopShieldingPtr()          const { return _pSTMTopShieldingParams.get(); }
    BottomShielding  const * getBottomShieldingPtr()       const { return _pSTMBottomShieldingParams.get(); }
    InnerShielding   const * getInnerShieldingPtr()        const { return _pSTMInnerShieldingParams.get(); }
    BackShielding    const * getBackShieldingPtr()        const { return _pSTMBackShieldingParams.get(); }

    ElectronicShielding   const * getElectronicShieldingPtr()        const { return _pSTMElectronicShieldingParams.get(); }
    STM_Absorber          const * getSTM_AbsorberPtr()               const { return _pSTMSTM_AbsorberParams.get(); }


    CLHEP::Hep3Vector const & originInMu2e() const { return _originInMu2e; };

    // True when the downstream shielding is the updated, hand-stacked version
    // This is one switch for the whole downstream geometry
    bool handstacked() const { return _handstacked; }


  private:

    friend class STMMaker;

    // The class should only be constructed via STM::STMMaker.
    STM(){};

    std::unique_ptr<STMDownstreamEnvelope>  _pSTMDnStrEnvParams;
    std::unique_ptr<PermanentMagnet>  _pSTMMagnetParams;
    std::unique_ptr<TransportPipe>    _pSTMTransportPipeParams;
    std::unique_ptr<SupportTable>     _pSTMMagnetSupportTableParams;
    std::unique_ptr<STMCollimator>    _pSTMFOVCollimatorParams;
    std::unique_ptr<STMCollimator>    _pSTMSSCollimatorParams;
    std::unique_ptr<SupportTable>     _pSTMDetectorSupportTableParams;
    std::unique_ptr<GeDetector>       _pSTMDetector1Params;
    std::unique_ptr<GeDetector>       _pSTMDetector2Params;
    std::unique_ptr<ShieldPipe>       _pSTMShieldPipeParams;
    std::unique_ptr<STM_SSC>          _pSTM_SSCParams;
    std::unique_ptr<SSCSupport>       _pSSCSupportParams;
    std::unique_ptr<LeadBrick>        _pLeadBrickParams;
    std::unique_ptr<SSCFrontShield>   _pSSCFrontShieldParams;
    std::unique_ptr<FrontShieldingRight> _pFrontShieldingRightParams;
    std::unique_ptr<FrontShieldingLeft>  _pFrontShieldingLeftParams;
    std::unique_ptr<RightWall>           _pRightWallParams;
    std::unique_ptr<LeftWall>            _pLeftWallParams;
    std::unique_ptr<TopWall>             _pTopWallParams;
    std::unique_ptr<BottomWall>          _pBottomWallParams;

    std::unique_ptr<HPGeDetector>     _pSTMHPGeDetectorParams;
    std::unique_ptr<LaBrDetector>     _pSTMLaBrDetectorParams;

    std::unique_ptr<FrontShielding>   _pSTMFrontShieldingParams;
    std::unique_ptr<LeftShielding>    _pSTMLeftShieldingParams;
    std::unique_ptr<RightShielding>   _pSTMRightShieldingParams;
    std::unique_ptr<TopShielding>     _pSTMTopShieldingParams;
    std::unique_ptr<BottomShielding>  _pSTMBottomShieldingParams;
    std::unique_ptr<InnerShielding>   _pSTMInnerShieldingParams;
    std::unique_ptr<BackShielding>    _pSTMBackShieldingParams;

    std::unique_ptr<ElectronicShielding>   _pSTMElectronicShieldingParams;
    std::unique_ptr<STM_Absorber>          _pSTMSTM_AbsorberParams;

    CLHEP::Hep3Vector   _originInMu2e;

    bool                _handstacked = false;
  };

}

#endif/*STMGeom_STM_hh*/
