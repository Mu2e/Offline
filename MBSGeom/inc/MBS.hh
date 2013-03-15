#ifndef MBSGeom_MBS_hh
#define MBSGeom_MBS_hh

// Muon Beam Stop Object
//
// $Id: MBS.hh,v 1.2 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
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

    Tube const * const getBSTSPtr() const {return _pBSTSParams.get();}
    Tube const * const getSPBSPtr() const {return _pSPBSParams.get();}
    Tube const * const getBSTCPtr() const {return _pBSTCParams.get();}
    Tube const * const getBSBSPtr() const {return _pBSBSParams.get();}
    Tube const * const getCLV2Ptr() const {return _pCLV2Params.get();}

    CLHEP::Hep3Vector const & originInMu2e() const { return _originInMu2e; };

  private:

    friend class MBSMaker;

    // The class should only be constructed via MBS::MBSMaker.
    MBS(){};

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    MBS( MBS const & );
    MBS const & operator= ( MBS const & );

    // several concentric components

    std::unique_ptr<Tube> _pBSTSParams;
    std::unique_ptr<Tube> _pSPBSParams;
    std::unique_ptr<Tube> _pBSTCParams;
    std::unique_ptr<Tube> _pBSBSParams;
    std::unique_ptr<Tube> _pCLV2Params;

    CLHEP::Hep3Vector   _originInMu2e;

  };

}

#endif/*MBSGeom_MBS_hh*/
