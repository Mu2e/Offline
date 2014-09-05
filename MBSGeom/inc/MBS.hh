#ifndef MBSGeom_MBS_hh
#define MBSGeom_MBS_hh

// Muon Beam Stop Object
//
// $Id: MBS.hh,v 1.5 2014/09/05 14:42:20 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/05 14:42:20 $
//
// Original author KLG
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class MBSMaker;

  class MBS : virtual public Detector {

  public:

    // Volume names a per figure 8.29 in the CDR
    Polycone const * getMBSMPtr()     const { return _pMBSMParams.get();     }
    Tube     const * getBSTSPtr()     const { return _pBSTSParams.get();     }
    Tube     const * getSPBSSup1Ptr() const { return _pSPBSSup1Params.get(); }
    Tube     const * getSPBSSup2Ptr() const { return _pSPBSSup2Params.get(); }
    Tube     const * getSPBSLPtr()    const { return _pSPBSLParams.get();    }
    Tube     const * getSPBSCPtr()    const { return _pSPBSCParams.get();    }
    Tube     const * getSPBSRPtr()    const { return _pSPBSRParams.get();    }
    Polycone const * getBSTCPtr()     const { return _pBSTCParams.get();     }
    Polycone const * getBSBSPtr()     const { return _pBSBSParams.get();     }
    Polycone const * getCLV2Ptr()     const { return _pCLV2Params.get();     }

    CLHEP::Hep3Vector const & originInMu2e() const { return _originInMu2e; };

    // DS3 vacuum volume is minimally extended to contain the following MBS envelope:
    double getEnvelopeZmax() const { return /*_pCLV2Params->originInMu2e().z() +*/ _zMax; }
    double getEnvelopeRmax() const { return _rMax; }
    double getEnvelopeRmin() const { return _rMin; }
    double getEnvelopeTLng() const { return _totLength; }

  private:

    friend class MBSMaker;

    // The class should only be constructed via MBS::MBSMaker.
    MBS(){};

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    MBS( MBS const & );
    MBS const & operator= ( MBS const & );

    // several concentric components

    std::unique_ptr<Polycone> _pMBSMParams;
    std::unique_ptr<Tube>     _pBSTSParams;
    std::unique_ptr<Tube>     _pSPBSSup1Params;
    std::unique_ptr<Tube>     _pSPBSSup2Params;
    std::unique_ptr<Tube>     _pSPBSLParams;
    std::unique_ptr<Tube>     _pSPBSCParams;
    std::unique_ptr<Tube>     _pSPBSRParams;
    std::unique_ptr<Polycone> _pBSTCParams;
    std::unique_ptr<Polycone> _pBSBSParams;
    std::unique_ptr<Polycone> _pCLV2Params;

    double _rMax;
    double _rMin;
    double _zMax;
    double _totLength;

    CLHEP::Hep3Vector   _originInMu2e;

  };

}

#endif/*MBSGeom_MBS_hh*/
