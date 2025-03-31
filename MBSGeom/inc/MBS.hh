#ifndef MBSGeom_MBS_hh
#define MBSGeom_MBS_hh

// Muon Beam Stop Object
//
//
// Original author KLG
//
// Modified 2015/10/30 by dnbrow01 to reflect current design of MBS
// as found in docdb-1351 v7.  Call this Version 2.
// 2016/08/05 dnbrow01:  Add holes for services access.  Call this version 3
// C++ includes
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class MBSMaker;

  class MBS : virtual public Detector {

  public:

    ~MBS() override = default;

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    MBS( MBS const& ) = delete;
    MBS( MBS&&      ) = delete;
    MBS& operator=( MBS const& ) = delete;
    MBS& operator=( MBS&&      ) = delete;

    // Volume names a per figure 8.29 in the CDR <== somewhat outdated now
    // MBSM is mother volume
    // BSTS is stainless steel pipe in three sections
    // SPBS is the outer HDPE shield that doesn't cover the bottom
    // BSTC is the upstream inner HDPE
    // BSBS is the downstream inner HDPE
    // CLV2 is the HDPE end plug on the downstream side.
    // First is the mother volume
    Polycone const * getMBSMPtr()     const { return _pMBSMParams.get();     }
    Polycone const * getBSTSPtr()     const { return _pBSTSParams.get();     }
    Tube     const * getSPBSSup1Ptr() const { return _pSPBSSup1Params.get(); }
    Tube     const * getSPBSSup2Ptr() const { return _pSPBSSup2Params.get(); }
    Tube     const * getSPBSLPtr()    const { return _pSPBSLParams.get();    }
    Tube     const * getSPBSCPtr()    const { return _pSPBSCParams.get();    }
    Tube     const * getSPBSRPtr()    const { return _pSPBSRParams.get();    }
    Polycone const * getBSTCPtr()     const { return _pBSTCParams.get();     }
    Polycone const * getBSBSPtr()     const { return _pBSBSParams.get();     }
    Polycone const * getCLV2Ptr()     const { return _pCLV2Params.get();     }
    Tube     const * getCLV2ABSPtr()  const { return _pCLV2ABSParams.get();  }
    Tube     const * getCalRingShieldPtr()  const { return _pCalShieldRingParams.get();  }
    int      const   getVersion()     const { return _version;               }

    CLHEP::Hep3Vector const & originInMu2e() const { return _originInMu2e; };

    // DS3 vacuum volume is minimally extended to contain the following MBS envelope:
    double getEnvelopeZmax() const { return /*_pCLV2Params->originInMu2e().z() +*/ _zMax; }
    double getEnvelopeRmax() const { return _rMax; }
    double getEnvelopeRmin() const { return _rMin; }
    double getEnvelopeTLng() const { return _totLength; }

    std::vector<CLHEP::Hep3Vector> getHoleCentersInSteel() const
    { return _holeCentersInSteel; }
    std::vector<CLHEP::Hep3Vector> getHoleCentersInUpstreamPoly() const
    { return _holeCentersInUpstreamPoly; }
    std::vector<CLHEP::Hep3Vector> getHoleCentersInDownstreamPoly() const
    { return _holeCentersInDownstreamPoly; }
    double getHoleXDimInSteel()    const { return _holeXDimInSteel; }
    double getHoleYDimInSteel()    const { return _holeYDimInSteel; }
    double getHoleZDimInSteel()    const { return _holeZDimInSteel; }
    double getHoleXDimInUpPoly()   const { return _holeXDimInUpPoly; }
    double getHoleYDimInUpPoly()   const { return _holeYDimInUpPoly; }
    double getHoleZDimInUpPoly()   const { return _holeZDimInUpPoly; }
    double getHoleXDimInDownPoly() const { return _holeXDimInDownPoly; }
    double getHoleYDimInDownPoly() const { return _holeYDimInDownPoly; }
    double getHoleZDimInDownPoly() const { return _holeZDimInDownPoly; }


  private:

    friend class MBSMaker;

    // The class should only be constructed via MBS::MBSMaker.
    MBS(){};

    // several concentric components

    std::unique_ptr<Polycone> _pMBSMParams; // mother envelope
    std::unique_ptr<Polycone> _pBSTSParams; // Stainless pipe (1-3 sections)
    std::unique_ptr<Tube>     _pSPBSSup1Params;
    std::unique_ptr<Tube>     _pSPBSSup2Params;
    std::unique_ptr<Tube>     _pSPBSLParams;
    std::unique_ptr<Tube>     _pSPBSCParams;// Outer HDPE
    std::unique_ptr<Tube>     _pSPBSRParams;
    std::unique_ptr<Polycone> _pBSTCParams; // Inner HDPE upstream
    std::unique_ptr<Polycone> _pBSBSParams; // Inner HDPE downstream
    std::unique_ptr<Polycone> _pCLV2Params; // HDPE end plug
    std::unique_ptr<Tube>     _pCLV2ABSParams;
    std::unique_ptr<Tube>     _pCalShieldRingParams; //Shield to protect the calorimeter

    double _rMax = 0;
    double _rMin = 0;
    double _zMax = 0;
    double _totLength = 0;

    CLHEP::Hep3Vector   _originInMu2e;

    int    _version = 0;
    std::vector<CLHEP::Hep3Vector> _holeCentersInSteel;
    std::vector<CLHEP::Hep3Vector> _holeCentersInUpstreamPoly;
    std::vector<CLHEP::Hep3Vector> _holeCentersInDownstreamPoly;
    double _holeXDimInSteel = 0.;
    double _holeYDimInSteel = 0.;
    double _holeZDimInSteel = 0.;
    double _holeXDimInUpPoly = 0.;
    double _holeYDimInUpPoly = 0.;
    double _holeZDimInUpPoly = 0.;
    double _holeXDimInDownPoly = 0.;
    double _holeYDimInDownPoly = 0.;
    double _holeZDimInDownPoly = 0.;

  };

}

#endif/*MBSGeom_MBS_hh*/
