#ifndef MBSGeom_MBSMaker_hh
#define MBSGeom_MBSMaker_hh
//
// Class to construct and return MBS
//
// $Id: MBSMaker.hh,v 1.3 2013/07/02 15:57:07 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/07/02 15:57:07 $
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

    double _MBSCZ;

    double _BSTSInnerRadius;
    double _BSTSOuterRadius;
    double _BSTSHLength;
    std::string _BSTSMaterialName;
    double _BSTSZ;

    double _SPBSSup1InnerRadius;
    double _SPBSSup1OuterRadius;
    double _SPBSSup1HLength;
    std::string _SPBSSup1MaterialName;
    double _SPBSSup1Z;
    double _SPBSSup2InnerRadius;
    double _SPBSSup2OuterRadius;
    double _SPBSSup2HLength;
    std::string _SPBSSup2MaterialName;
    double _SPBSSup2Z;

    double _SPBSLInnerRadius;
    double _SPBSLOuterRadius;
    double _SPBSLHLength;
    std::string _SPBSLMaterialName;
    double _SPBSLZ;
    double _SPBSCInnerRadius;
    double _SPBSCOuterRadius;
    double _SPBSCHLength;
    double _SPBSCminAngle;
    double _SPBSCmaxAngle;
   std::string _SPBSCMaterialName;
    double _SPBSCZ;
    double _SPBSRInnerRadius;
    double _SPBSROuterRadius;
    double _SPBSRHLength;
    std::string _SPBSRMaterialName;
    double _SPBSRZ;

    std::vector<double> _BSTCInnerRadii;
    std::vector<double> _BSTCOuterRadii;
    std::vector<double> _BSTCLengths;
    std::string _BSTCMaterialName;
    double _BSTCZ;
    std::vector<double> _BSBSInnerRadii;
    std::vector<double> _BSBSOuterRadii;
    std::vector<double> _BSBSLengths;
    std::string _BSBSMaterialName;
    double _BSBSZ;
    std::vector<double> _CLV2InnerRadii;
    std::vector<double> _CLV2OuterRadii;
    std::vector<double> _CLV2Lengths;
    std::string _CLV2MaterialName;
    double _CLV2Z;

  };

}  //namespace mu2e

#endif /* MBSGeom_MBSMaker_hh */
