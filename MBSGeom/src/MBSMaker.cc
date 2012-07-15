//
// Construct and return MBS
//
// $Id: MBSMaker.cc,v 1.2 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//
// Original author KLG
//
// Notes
// see mu2e-doc-1351 for naming convetions etc...

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "MBSGeom/inc/MBSMaker.hh"
#include "MBSGeom/inc/MBS.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  MBSMaker::MBSMaker(SimpleConfig const & _config,
                     double solenoidOffset)
  {

    // if( ! _config.getBool("hasMBS",false) ) return;

    // create an empty MBS
    _mbs = auto_ptr<MBS>(new MBS());

    // access its object through a reference

    MBS & mbs = *_mbs.get();

    parseConfig(_config);

    // now create the specific components

    CLHEP::Hep3Vector _BSTSOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_BSTSZ);
    mbs._originInMu2e = _BSTSOffsetInMu2e;

    mbs._pBSTSParams = std::auto_ptr<Tube>
      (new Tube(_BSTSMaterialName,
                _BSTSOffsetInMu2e,
                _BSTSInnerRadius,
                _BSTSOuterRadius,
                _BSTSHLength));

    CLHEP::Hep3Vector _SPBSOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_SPBSZ);

    mbs._pSPBSParams = std::auto_ptr<Tube>
      (new Tube(_SPBSMaterialName,
                _SPBSOffsetInMu2e,
                _SPBSInnerRadius,
                _SPBSOuterRadius,
                _SPBSHLength));

    CLHEP::Hep3Vector _BSTCOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_BSTCZ);

    mbs._pBSTCParams = std::auto_ptr<Tube>
      (new Tube(_BSTCMaterialName,
                _BSTCOffsetInMu2e,
                _BSTCInnerRadius,
                _BSTCOuterRadius,
                _BSTCHLength));

    CLHEP::Hep3Vector _BSBSOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_BSBSZ);

    mbs._pBSBSParams = std::auto_ptr<Tube>
      (new Tube(_BSBSMaterialName,
                _BSBSOffsetInMu2e,
                _BSBSInnerRadius,
                _BSBSOuterRadius,
                _BSBSHLength));

    CLHEP::Hep3Vector _CLV2OffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_CLV2Z);

    mbs._pCLV2Params = std::auto_ptr<Tube>
      (new Tube(_CLV2MaterialName,
                _CLV2OffsetInMu2e,
                _CLV2InnerRadius,
                _CLV2OuterRadius,
                _CLV2HLength));

  }

  void MBSMaker::parseConfig( SimpleConfig const & _config ){

    _verbosityLevel                   = _config.getInt("mbs.verbosityLevel");

    _BSTSInnerRadius   = _config.getDouble("mbs.BSTSInnerRadius");
    _BSTSOuterRadius   = _config.getDouble("mbs.BSTSOuterRadius");
    _BSTSHLength       = _config.getDouble("mbs.BSTSHLength");
    _BSTSMaterialName  = _config.getString("mbs.BSTSMaterialName");
    _BSTSZ             = _config.getDouble("mbs.BSTSZ");
    _SPBSInnerRadius   = _config.getDouble("mbs.SPBSInnerRadius");
    _SPBSOuterRadius   = _config.getDouble("mbs.SPBSOuterRadius");
    _SPBSHLength       = _config.getDouble("mbs.SPBSHLength");
    _SPBSMaterialName  = _config.getString("mbs.SPBSMaterialName");
    _SPBSZ             = _BSTSZ + _BSTSHLength - _SPBSHLength;
    _BSTCInnerRadius   = _config.getDouble("mbs.BSTCInnerRadius");
    _BSTCOuterRadius   = _config.getDouble("mbs.BSTCOuterRadius");
    _BSTCHLength       = _config.getDouble("mbs.BSTCHLength");
    _BSTCMaterialName  = _config.getString("mbs.BSTCMaterialName");
    _BSTCZ             = _BSTSZ - _BSTSHLength + _BSTCHLength;
    _BSBSInnerRadius   = _config.getDouble("mbs.BSBSInnerRadius");
    _BSBSOuterRadius   = _config.getDouble("mbs.BSBSOuterRadius");
    _BSBSHLength       =_BSTSHLength -_BSTCHLength;
    _BSBSMaterialName  = _config.getString("mbs.BSBSMaterialName");
    _BSBSZ             = _BSTSZ + _BSTSHLength - _BSBSHLength;
    _CLV2InnerRadius   = _config.getDouble("mbs.CLV2InnerRadius");
    _CLV2OuterRadius   = _config.getDouble("mbs.CLV2OuterRadius");
    _CLV2HLength       = _config.getDouble("mbs.CLV2HLength");
    _CLV2MaterialName  = _config.getString("mbs.CLV2MaterialName");
    _CLV2Z             = _BSTSZ + _BSTSHLength - _CLV2HLength;

  }

} // namespace mu2e
