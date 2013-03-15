#ifndef MBSGeom_MBSMaker_hh
#define MBSGeom_MBSMaker_hh
//
// Class to construct and return MBS
//
// $Id: MBSMaker.hh,v 1.2 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author KLG
//

#include <memory>
#include <string>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class MBS;
  class SimpleConfig;
  class Tube;

  class MBSMaker {

  public:

    MBSMaker( SimpleConfig const & config,
              double solenoidOffset);

    void parseConfig( SimpleConfig const & _config );

    // This is deprecated and will go away soon.
    // Still needed for root graphics version.
    const MBS& getMBS() const { return *_mbs;}

    // This is the accessor that will remain.
    std::unique_ptr<MBS> getMBSPtr() { return std::move(_mbs); }

  private:

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    MBSMaker( MBSMaker const & );
    MBSMaker const & operator= ( MBSMaker const & );

    std::unique_ptr<MBS> _mbs;

    int    _verbosityLevel;

    double _BSTSInnerRadius;
    double _BSTSOuterRadius;
    double _BSTSHLength;
    std::string _BSTSMaterialName;
    double _BSTSZ;
    double _SPBSInnerRadius;
    double _SPBSOuterRadius;
    double _SPBSHLength;
    std::string _SPBSMaterialName;
    double _SPBSZ;
    double _BSTCInnerRadius;
    double _BSTCOuterRadius;
    double _BSTCHLength;
    std::string _BSTCMaterialName;
    double _BSTCZ;
    double _BSBSInnerRadius;
    double _BSBSOuterRadius;
    double _BSBSHLength;
    std::string _BSBSMaterialName;
    double _BSBSZ;
    double _CLV2InnerRadius;
    double _CLV2OuterRadius;
    double _CLV2HLength;
    std::string _CLV2MaterialName;
    double _CLV2Z;

  };

}  //namespace mu2e

#endif /* MBSGeom_MBSMaker_hh */
