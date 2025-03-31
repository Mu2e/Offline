#ifndef MBSGeom_MBSMaker_hh
#define MBSGeom_MBSMaker_hh
//
// Class to construct and return MBS
//
//
// Original author KLG
//
// Update by D. No. Brown to Version 2, representing version in
// WBS 5.8 v7 (Mu2e-doc-1351 v7)

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class MBS;
  class SimpleConfig;
  class Tube;

  class MBSMaker {

  public:

    MBSMaker( SimpleConfig const & config,
              double solenoidOffset);
    ~MBSMaker() = default;

    // delete automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    MBSMaker( MBSMaker const & ) = delete;
    MBSMaker( MBSMaker&&       ) = delete;
    MBSMaker& operator=( MBSMaker const & ) = delete;
    MBSMaker& operator=( MBSMaker&&       )  = delete;

    void parseConfig( SimpleConfig const & _config );

    // This is deprecated and will go away soon.
    // Still needed for root graphics version.
    const MBS& getMBS() const { return *_mbs;}

    // This is the accessor that will remain.
    std::unique_ptr<MBS> getMBSPtr() { return std::move(_mbs); }

  private:

    std::unique_ptr<MBS> _mbs;

    int    _verbosityLevel;
    int    _MBSVersion;

    double _MBSCZ;

    double _BSTSInnerRadius;  // Version 1 only
    double _BSTSOuterRadius;  // Versions 1 and 2
    std::vector<double> _BSTSInnRadii; // Version 2 only
    std::vector<double> _BSTSOutRadii; // Version 2 only
    std::vector<double> _BSTSZLengths; // Version 2 only
    double _BSTSHLength;  // Versions 1 and 2
    std::string _BSTSMaterialName;  // Versions 1 and 2
    double _BSTSZ;  // Versions 1 and 2

    // Support Ring 1 - Version 1 only
    double _SPBSSup1InnerRadius;
    double _SPBSSup1OuterRadius;
    double _SPBSSup1HLength;
    std::string _SPBSSup1MaterialName;
    double _SPBSSup1Z;

    // Support Ring 2 - Version 1 only
    double _SPBSSup2InnerRadius;
    double _SPBSSup2OuterRadius;
    double _SPBSSup2HLength;
    std::string _SPBSSup2MaterialName;
    double _SPBSSup2Z;

    // Upstream outer HDPE - Version 1 only
    double _SPBSLInnerRadius;
    double _SPBSLOuterRadius;
    double _SPBSLHLength;
    std::string _SPBSLMaterialName;
    double _SPBSLZ;

    // Midstream outer HDPE - Versions 1 and 2
    double _SPBSCInnerRadius;
    double _SPBSCOuterRadius;
    double _SPBSCHLength;
    double _SPBSCminAngle;
    double _SPBSCmaxAngle;
   std::string _SPBSCMaterialName;
    double _SPBSCZ;

    // Downstream outer HDPE - Version 1 only
    double _SPBSRInnerRadius;
    double _SPBSROuterRadius;
    double _SPBSRHLength;
    std::string _SPBSRMaterialName;
    double _SPBSRZ;

    // Inner HDPE upstream - Versions 1 and 2
    std::vector<double> _BSTCInnerRadii;
    std::vector<double> _BSTCOuterRadii;
    std::vector<double> _BSTCLengths;
    std::string _BSTCMaterialName;
    double _BSTCZ;

    // Inner HDPE downstream - Versions 1 and 2
    std::vector<double> _BSBSInnerRadii;
    std::vector<double> _BSBSOuterRadii;
    std::vector<double> _BSBSLengths;
    std::string _BSBSMaterialName;
    double _BSBSZ;

    // HDPE end plug - Versions 1 and 2
    std::vector<double> _CLV2InnerRadii;
    std::vector<double> _CLV2OuterRadii;
    std::vector<double> _CLV2Lengths;
    std::string _CLV2MaterialName;
    double _CLV2Z;

    // Absorber
    bool _CLV2AbsBuild;
    std::string _CLV2AbsMaterialName;
    double _CLV2AbsHLength;

    // Holes in MBS
    int _nHolesSt;
    int _nHolesUP;
    int _nHolesDP;
    double _BSTSHoleXDim;
    double _BSTSHoleYDim;
    double _BSTSHoleZDim;
    double _upPolyHoleXDim;
    double _upPolyHoleYDim;
    double _upPolyHoleZDim;
    double _downPolyHoleXDim;
    double _downPolyHoleYDim;
    double _downPolyHoleZDim;
    std::vector<CLHEP::Hep3Vector> _BSTSHoleCenters;
    std::vector<CLHEP::Hep3Vector> _upPolyHoleCenters;
    std::vector<CLHEP::Hep3Vector> _downPolyHoleCenters;


    //Calorimeter shield definitions
    double _CalShieldRingInnerRadius;
    double _CalShieldRingOuterRadius;
    double _CalShieldRingHLength;
    std::string _CalShieldRingMaterialName;
    double _CalShieldRingZ;

  };

}  //namespace mu2e

#endif /* MBSGeom_MBSMaker_hh */
