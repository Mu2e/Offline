//
// Construct and return ProtonAbsorber
//
//
// Original author MyeongJae Lee
//
// All in mu2e coordinate system.
//
// To modify the size of MECO-style (conical) proton absorber, use :
// bool protonabsorber.isShorterCone = true;
// and modify protonabsorber.halfLength and protonabsorber.distFromTargetEnd
// Modifing the radii protobabsorber.OutRaduis[0,1] is not recommanded.
//
// Updated by David Norvil Brown to allow for a "barrel of slats" geometry
// for the outer proton absorber.
#include <iostream>
#include <iomanip>
#include <cmath>


// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/BeamlineGeom/inc/ProtonAbsorber.hh"
#include "Offline/GeometryService/inc/ProtonAbsorberMaker.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  ProtonAbsorberMaker::ProtonAbsorberMaker(SimpleConfig const& _config, const DetectorSolenoid& ds, const StoppingTarget& target)
    : _ds(&ds), _target(&target)
  {
    BuildIt (_config);

    int verbosity(_config.getInt("protonabsorber.verbosityLevel",0));
    if ( verbosity > 0 ) PrintConfig();

  }

  void ProtonAbsorberMaker::BuildIt ( SimpleConfig const& _config)
  {
    _IPAVersion = _config.getInt("protonabsorber.version", 1);
    Build(_config);
  }

  void ProtonAbsorberMaker::Build ( SimpleConfig const& _config)
  {

    //////////////////////////////////////
    // All variables read from geometry txt files
    //////////////////////////////////////

    // geometry and material
    // not recommended to modify following variables
    double r1out0        = _config.getDouble("protonabsorber.OutRadius0", 335.2);
    double r2out1        = _config.getDouble("protonabsorber.OutRadius1", 380.0);
    // new variable for the length of PA region
    double pabsZHalfLen  = _config.getDouble("protonabsorber.halfLength", 1250.0);
    double thick             = _config.getDouble("protonabsorber.thickness", 0.5);
    double distFromTargetEnd = _config.getDouble("protonabsorber.distFromTargetEnd", 0);
    // region for PA
    const double pabsRegionLength = 2500.;
    std::string materialName = _config.getString("protonabsorber.materialName", "Polyethylene092");
    if (distFromTargetEnd > pabsRegionLength || distFromTargetEnd <0) {
      distFromTargetEnd = 0;
    }

    int oPAversion              = _config.getInt("protonabsorber.outerPAversion",1);
    int oPAnSides = 0;
    double oPAslotWidth  = 0.0;
    double oPAslotLength = 0.0;
    double oPAslotOffset = 0.0;
    double oPASlatHeight = 0.0;
    double oPASlatWidth  = 0.0;
    double oPASlatLength = 0.0;
    if ( oPAversion > 1 ) {
      oPAnSides = _config.getInt("protonabsorber.outerPAnSides",0);
      oPAslotWidth = _config.getDouble("protonabsorber.outerPAslotWidth",0.0);
      oPAslotLength = _config.getDouble("protonabsorber.outerPAslotLength",0.0);
      oPAslotOffset = _config.getDouble("protonabsorber.outerPAslotOffset",0.0);
      oPASlatHeight = _config.getDouble("protonabsorber.oPASupportSlatHeight",0.0);
      oPASlatWidth = _config.getDouble("protonabsorber.oPASupportSlatWidth",0.0);
      oPASlatLength = _config.getDouble("protonabsorber.oPASupportSlatLength",0.0);
    }
    double oPAin0               = _config.getDouble("protonabsorber.outerPAInnerRadius0", 436.0);
    double oPAin1               = _config.getDouble("protonabsorber.outerPAInnerRadius1", 720.0);
    double oPAhl                = _config.getDouble("protonabsorber.outerPAHalfLength", 2250.0);
    double oPAthick             = _config.getDouble("protonabsorber.outerPAThickness", 10.0);
    double oPAzcenter           = _config.getDouble("protonabsorber.outerPAZCenter", 6250.0);
    std::string oPAmaterialName = _config.getString("protonabsorber.outerPAMaterialName", "Polyethylene092");
    double ds23split = _ds->vac_zLocDs23Split();


    // TS and DS geometry for locating the center of Proton Absorber
    double solenoidOffset = _config.getDouble("mu2e.solenoidOffset");

    // adding virtual detector before and after target
    double vdHL = 0.;
    if (_config.getBool("hasVirtualDetector", true) ) {
      vdHL = _config.getDouble("vd.halfLength");
    }

    // subtract virtual detector thickness from the larger outer
    // radius of the proton absorber
    r2out1 -= 2.*vdHL;



    /////////////////////////////
    // Geometry related variables
    /////////////////////////////

    // z of target end in mu2e coordinate
    // we add space for the virtual detector here
    double targetEnd = _target->centerInMu2e().z() + 0.5*_target->cylinderLength() + 2.*vdHL;;

    // distance from target end to ds2-ds3 boundary
    double targetEndToDS2End = _ds->vac_zLocDs23Split() - targetEnd;

    //////////////////////////////////////
    // Decide which pabs will be turned on
    //////////////////////////////////////

    bool pabs1 = true;
    bool pabs2 = true;
    bool opa1 = _config.getBool("protonabsorber.outerPA", false);
    bool opa2 = true;

    // if pabs starts from DS3 region
    if (distFromTargetEnd > targetEndToDS2End) {
      pabs1 = false;
    }

    // if pabs is short enouhg to locate at DS2 region only
    if (distFromTargetEnd + pabsZHalfLen*2.< targetEndToDS2End) {
      pabs2 = false;
    }
    if (!pabs1 && !pabs2) {
      return;
    }

    // if opa2 is short enough to locate at DS2 region only
    if (oPAzcenter + oPAhl< ds23split) {
      opa2 = false;
    }

    //////////////
    // Half length
    //////////////

    double pabs1halflen = 0, pabs2halflen = 0, pabs3halflen = 0, pabs4halflen = 0;
    if (pabs1) {
      pabs1halflen = (targetEndToDS2End - distFromTargetEnd)*0.5;
      if (pabs1halflen > pabsZHalfLen ) pabs1halflen = pabsZHalfLen;
    }
    if (pabs2) {
      pabs2halflen = pabsZHalfLen - pabs1halflen;
    }

    pabs3halflen = oPAhl;
    if (opa2) {
      pabs3halflen = oPAhl - (oPAzcenter + oPAhl - ds23split)*0.5;
      pabs4halflen = (oPAzcenter + oPAhl - ds23split)*0.5;
    }


    /////////
    // Offset
    /////////

    double pabs1ZOffset = targetEnd + 0.01;
    double pabs2ZOffset = targetEnd + targetEndToDS2End + 0.01;
    if (pabs1) {
      pabs1ZOffset = targetEnd + distFromTargetEnd + pabs1halflen;
    }
    if (pabs2) {
      if (pabs1) {
        pabs2ZOffset = targetEnd + targetEndToDS2End + pabs2halflen;
      }
      else {
        pabs2ZOffset = targetEnd + distFromTargetEnd + pabs2halflen;
      }
    }


    double pabs3ZOffset = oPAzcenter;
    double pabs4ZOffset = 0;
    if (opa2) {
      pabs3ZOffset = pabs3ZOffset  - pabs4halflen;
      pabs4ZOffset = ds23split + pabs4halflen;
    }

    /////////
    // Radius
    /////////

    double pabs1rIn0 = 0, pabs1rIn1 = 0, pabs1rOut0 = 0, pabs1rOut1 = 0;
    double pabs2rIn0 = 0, pabs2rIn1 = 0, pabs2rOut0 = 0, pabs2rOut1 = 0;

    if (pabs1) {
      pabs1rOut0 = (r2out1 - r1out0)/pabsRegionLength*(pabs1ZOffset-pabs1halflen-targetEnd) + r1out0;
      pabs1rOut1 = (r2out1 - r1out0)/pabsRegionLength*(pabs1ZOffset+pabs1halflen-targetEnd) + r1out0;
      pabs1rIn0  = pabs1rOut0 - thick;
      pabs1rIn1  = pabs1rOut1 - thick;
    }
    if (pabs2) {
      pabs2rOut0 = (r2out1 - r1out0)/pabsRegionLength*(pabs2ZOffset-pabs2halflen-targetEnd) + r1out0;
      pabs2rOut1 = (r2out1 - r1out0)/pabsRegionLength*(pabs2ZOffset+pabs2halflen-targetEnd) + r1out0;
      pabs2rIn0  = pabs2rOut0 - thick;
      pabs2rIn1  = pabs2rOut1 - thick;
    }

    double oPAout0 = oPAin0 + oPAthick;
    double oPAout1 = oPAin1 + oPAthick;
    double pabs3rIn0 = oPAin0, pabs3rIn1 = oPAin1, pabs3rOut0 = oPAout0, pabs3rOut1 = oPAout1;
    double pabs4rIn0 = 0, pabs4rIn1 = 0, pabs4rOut0 = 0, pabs4rOut1 = 0;

    if (opa2) {

      // pabs3rOut1 = pabs3halflen/(oPAhl - pabs4halflen)*oPAout1;
      // pabs3rIn1  = pabs3rOut1 - thick;

      pabs3rOut1 = oPAout0 + (oPAout1 - oPAout0)*pabs3halflen/oPAhl;
      pabs3rIn1  = pabs3rOut1 - oPAthick;

      pabs4rOut0 = pabs3rOut1;
      pabs4rIn0  = pabs4rOut0 - oPAthick;
      pabs4rOut1 = oPAout1;
      pabs4rIn1  = pabs4rOut1 - oPAthick;
    }

    /////////
    // Build
    /////////

    _pabs = unique_ptr<ProtonAbsorber>(new ProtonAbsorber());
    _pabs->_oPAVersion = oPAversion;
    _pabs->_nOPASupportSlats = 3; //default for versions before version 3
    _pabs->_nOPASupportSlatTypes = 1; //default for versions before version 3
    if(oPAversion > 2) {
      _pabs->_nOPASupportSlats     = _config.getInt("protonabsorber.nOPASupportSlats");
      _config.getVectorDouble("protonabsorber.oPASupportSlatAngles",_pabs->_oPASlatAngles,_pabs->_nOPASupportSlats);
      _config.getVectorInt   ("protonabsorber.oPASupportSlatTypes",_pabs->_oPASlatTypes,_pabs->_nOPASupportSlats);
      _pabs->_nOPASupportSlatTypes = _config.getInt("protonabsorber.nOPASupportSlatTypes");
      _config.getVectorString("protonabsorber.oPASupportSlatMaterials",_pabs->_oPASlatMaterials,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatHeights",_pabs->_oPASlatHeights,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatWidths",_pabs->_oPASlatWidths,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatLengths",_pabs->_oPASlatLengths,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatSideThicknesses",_pabs->_oPASlatSideThicknesses,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatTopThicknesses",_pabs->_oPASlatTopThicknesses,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatFillParameter1",_pabs->_oPASlatFillParameter1,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatFillParameter2",_pabs->_oPASlatFillParameter2,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatFillParameter3",_pabs->_oPASlatFillParameter3,_pabs->_nOPASupportSlatTypes);
      _config.getVectorDouble("protonabsorber.oPASupportSlatFillParameter4",_pabs->_oPASlatFillParameter4,_pabs->_nOPASupportSlatTypes);
      _config.getVectorString("protonabsorber.oPASupportSlatFillMaterials",_pabs->_oPASlatFillMaterials,_pabs->_nOPASupportSlatTypes);
      //cross supports
      _pabs->_nCrossSupports       = _config.getInt   ("protonabsorber.nCrossSupports");
      _pabs->_crossSupportMaterial = _config.getString("protonabsorber.crossSupportMaterial");
      _config.getVectorDouble("protonabsorber.crossSupportThicknesses",_pabs->_crossSupportThicknesses,_pabs->_nCrossSupports);
      _config.getVectorDouble("protonabsorber.crossSupportWidth"      ,_pabs->_crossSupportWidth      ,_pabs->_nCrossSupports);
      _config.getVectorInt   ("protonabsorber.crossSupportOneIndex"   ,_pabs->_crossSupportOneIndex   ,_pabs->_nCrossSupports);
      _config.getVectorInt   ("protonabsorber.crossSupportTwoIndex"   ,_pabs->_crossSupportTwoIndex   ,_pabs->_nCrossSupports);
      _config.getVectorDouble("protonabsorber.crossSupportPhis"       ,_pabs->_crossSupportPhis       ,_pabs->_nCrossSupports);
      _config.getVectorDouble("protonabsorber.crossSupportHeights"    ,_pabs->_crossSupportHeights    ,_pabs->_nCrossSupports);
      _config.getVectorDouble("protonabsorber.crossSupportRadii"      ,_pabs->_crossSupportRadii      ,_pabs->_nCrossSupports);

    } else {
      std::string opaSupportMaterial = _config.getString("protonabsorber.oPASupportMaterialName");
      for(int iSlat = 0; iSlat < _pabs->_nOPASupportSlatTypes; ++iSlat) {
        _pabs->_oPASlatMaterials.push_back(opaSupportMaterial);
        _pabs->_oPASlatHeights  .push_back(oPASlatHeight);
        _pabs->_oPASlatWidths   .push_back(oPASlatWidth);
        _pabs->_oPASlatLengths  .push_back(oPASlatLength);
      }
      _pabs->_oPASlatAngles.push_back(0.);
      _pabs->_oPASlatAngles.push_back(120.);
      _pabs->_oPASlatAngles.push_back(240.);
      _pabs->_oPASlatTypes.push_back(0);
      _pabs->_oPASlatTypes.push_back(0);
      _pabs->_oPASlatTypes.push_back(0);
    }
    _pabs->_oPAslotWidth  = oPAslotWidth;
    _pabs->_oPAslotLength = oPAslotLength;
    _pabs->_oPAslotOffset = oPAslotOffset;
    CLHEP::Hep3Vector pabs1Offset(-1.*solenoidOffset, 0.0, pabs1ZOffset);
    CLHEP::Hep3Vector pabs2Offset(-1.*solenoidOffset, 0.0, pabs2ZOffset);
    CLHEP::Hep3Vector pabs3Offset(-1.*solenoidOffset, 0.0, pabs3ZOffset);
    CLHEP::Hep3Vector pabs4Offset(-1.*solenoidOffset, 0.0, pabs4ZOffset);

    // std::cout<< "Outer proton Z-end in Mu2e: "   << pabs4ZOffset+pabs4halflen << std::endl;
    // std::cout<< "Outer proton Z-start in Mu2e: " << pabs3ZOffset-pabs3halflen << std::endl;

    _pabs->_parts.push_back( ProtonAbsorberPart( 0, pabs1Offset, pabs1rOut0, pabs1rIn0, pabs1rOut1, pabs1rIn1, pabs1halflen, materialName));
    _pabs->_parts.push_back( ProtonAbsorberPart( 1, pabs2Offset, pabs2rOut0, pabs2rIn0, pabs2rOut1, pabs2rIn1, pabs2halflen, materialName));
    if ( oPAversion > 1 ) {
      _pabs->_parts.push_back( ProtonAbsorberPart( 2, pabs3Offset, pabs3rOut0, pabs3rIn0, pabs3rOut1, pabs3rIn1, pabs3halflen, oPAnSides, oPAmaterialName));
      _pabs->_parts.push_back( ProtonAbsorberPart( 3, pabs4Offset, pabs4rOut0, pabs4rIn0, pabs4rOut1, pabs4rIn1, pabs4halflen, oPAnSides, oPAmaterialName));
    } else {
      _pabs->_parts.push_back( ProtonAbsorberPart( 2, pabs3Offset, pabs3rOut0, pabs3rIn0, pabs3rOut1, pabs3rIn1, pabs3halflen, oPAmaterialName));
      _pabs->_parts.push_back( ProtonAbsorberPart( 3, pabs4Offset, pabs4rOut0, pabs4rIn0, pabs4rOut1, pabs4rIn1, pabs4halflen, oPAmaterialName));
    }

    // global variables
    (_pabs->_pabs1flag) = pabs1;
    (_pabs->_pabs2flag) = pabs2;
    (_pabs->_materialName) = materialName;
    (_pabs->_vdHL) = vdHL;
    (_pabs->_distfromtargetend) = distFromTargetEnd;
    (_pabs->_halflength) = pabsZHalfLen;
    (_pabs->_thickness) = thick;
    // glaboal variable for OPA
    (_pabs->_oPA1flag) = opa1;
    (_pabs->_oPA2flag) = opa2;
    (_pabs->_oPAmaterialName) = oPAmaterialName;
    (_pabs->_oPAzcenter) = oPAzcenter;
    (_pabs->_oPAhalflength) = oPAhl;
    (_pabs->_oPAthickness) = oPAthick;

    ////////////////////////
    // Support Structure for OPA
    ////////////////////////
    (_pabs->_oPAnSupports) = _config.getInt("protonabsorber.oPASupportNSupportRing",0);
    if ( _pabs->_oPAnSupports > 0 ) {
      if(oPAversion > 2) //each support can have individual materials
        _config.getVectorString("protonabsorber.oPASupportMaterials",_pabs->_oPASupportMaterials,_pabs->_oPAnSupports);
      else {//assume one material for all supports
        std::string opaSupportMaterial = _config.getString("protonabsorber.oPASupportMaterialName");
        for(int iSup = 0; iSup < _pabs->_oPAnSupports; ++iSup) _pabs->_oPASupportMaterials.push_back(opaSupportMaterial);
      }
      _config.getVectorDouble("protonabsorber.oPASupportInnerRadii",_pabs->_oPASupportIR,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportOuterRadii",_pabs->_oPASupportOR,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportHalflength",_pabs->_oPASupportHL,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportZMidpoints",_pabs->_oPASupportZM,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportHasExtra",_pabs->_oPASupportHE,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportExtraRad",_pabs->_oPASupportXR,(_pabs->_oPAnSupports));
      _config.getVectorDouble("protonabsorber.oPASupportExtraDPhi",_pabs->_oPASupportPH,(_pabs->_oPAnSupports));
    }

    ////////////////////////
    // Support Structure for IPA
    ////////////////////////

    const std::size_t supportSets     = _config.getInt   ("protonabsorber.ipa.nSets"          ,0 );
    const std::size_t wiresPerSet     = _config.getInt   ("protonabsorber.ipa.nWiresPerSet"   ,0 );
    const double      wireRadius      = _config.getDouble("protonabsorber.ipa.wireRadius"     ,0.);
    const string      wireMaterial    = _config.getString("protonabsorber.ipa.wireMaterial"      );
    const double      wireAngleOffset = _config.getDouble("protonabsorber.ipa.wireAngleOffset",0.);

    _pabs->_buildSupports = _config.getBool("protonabsorber.ipa.buildSupports"   );
    _pabs->_ipaSupport.reset( new InnerProtonAbsSupport( supportSets, wiresPerSet ) );
    _pabs->_ipaSupport->setWireAngleOffset(wireAngleOffset);
    // Calculate zPosition spacing
    double zSpacing = 0;
    if (_IPAVersion == 1) {
      zSpacing = 2*pabs1halflen/(supportSets+1);
    }
    else if (_IPAVersion >= 2) {
      zSpacing = 2*pabs1halflen; // want one set at each end of the IPA so the spacing between the sets for v2
    }

    const double ipazstart = targetEnd + distFromTargetEnd;
    const double opazstart = _pabs->_oPAzcenter - _pabs->_oPAhalflength;

    // Will be used for v2 only
    const double wireRotation = _config.getDouble("protonabsorber.ipa.wireRotationToVertical", 45);
    const double wire_rotation_from_vertical = wireRotation*CLHEP::rad;

    //******************
    // Pion Degrader
    //******************

    _pabs->_degraderBuild = _config.getBool("degrader.build"  ,false);
    if ( _pabs->_degraderBuild ) {
      _pabs->_degraderVersion = _config.getInt   ("degrader.version"  ,2);
      _pabs->_degraderRot     = _config.getDouble("degrader.rotation");
      _pabs->_degraderZ0      = _config.getDouble("degrader.upstreamEdge.z");
      _pabs->_degraderFiltMaterial = _config.getString("degrader.filter.materialName");

      if (_pabs->_degraderVersion >= 3) {
        _pabs->_degraderFilter2Material   = _config.getString("degrader.filter2.materialName");
        _pabs->_degraderConverterMaterial = _config.getString("degrader.converter.materialName");
      }

      _pabs->_degraderFramMaterial = _config.getString("degrader.frame.materialName");
      _pabs->_degraderCowtMaterial = _config.getString("degrader.counterweight.materialName");
      _pabs->_degraderRodMaterial  = _config.getString("degrader.rod.materialName");
      _pabs->_degraderSuptMaterial = _config.getString("degrader.support.materialName","G4_Al");

      double rin = _config.getDouble("degrader.frame.rIn");
      _pabs->_degraderFrameDims.push_back(rin);
      double rout = _config.getDouble("degrader.frame.rOut");
      _pabs->_degraderFrameDims.push_back(rout);
      double hl  = _config.getDouble("degrader.frame.halfLength");
      _pabs->_degraderFrameDims.push_back(hl);
      double coffs = _config.getDouble("degrader.frame.centerOffFromPivot");
      _pabs->_degraderFrameDims.push_back(coffs);
      // Based on doc-32335, frame may not form complete circle, so allow
      // for arc angle.
      double arcAngle = _config.getDouble("degrader.frame.arcAngle",360.0);
      _pabs->_degraderFrameDims.push_back(arcAngle);

      // Pivot point on mount, in x and y
      double xPivot = _config.getDouble("degrader.pivot.x",-3576.5);
      _pabs->_degraderPivotPos.push_back(xPivot);
      double yPivot = _config.getDouble("degrader.pivot.y",-567.25);
      _pabs->_degraderPivotPos.push_back(yPivot);

      // parameters for filter itself
      rin = _config.getDouble("degrader.filter.rIn");
      _pabs->_degraderFilterDims.push_back(rin);
      rout = _config.getDouble("degrader.filter.rOut");
      _pabs->_degraderFilterDims.push_back(rout);
      hl  = _config.getDouble("degrader.filter.halfLength");
      _pabs->_degraderFilterDims.push_back(hl);
      _pabs->_degraderFilterDims.push_back(coffs);

      rin = _config.getDouble("degrader.counterweight.rIn");
      _pabs->_degraderCounterDims.push_back(rin);
      rout = _config.getDouble("degrader.counterweight.rOut");
      _pabs->_degraderCounterDims.push_back(rout);
      hl  = _config.getDouble("degrader.counterweight.halfLength");
      _pabs->_degraderCounterDims.push_back(hl);
      coffs = _config.getDouble("degrader.counterweight.centerOffFromPivot");
      _pabs->_degraderCounterDims.push_back(coffs);
      double width = _config.getDouble("degrader.rod.width");
      _pabs->_degraderRodDims.push_back(width);
      double depth = _config.getDouble("degrader.rod.depth");
      _pabs->_degraderRodDims.push_back(depth);

      // Allow for a support structure for stepper motor and axle for
      // pion degrader.  There will be four "arms" sticking out from the
      // OPA support upstream side and a plate mounted on them.  First
      // configure the arms.
      // First get the half-length dimensions of the arms
      double hlx = _config.getDouble("degrader.supportArm.dx",0.0);
      _pabs->_degraderSupportArmDims.push_back(hlx);
      double hly = _config.getDouble("degrader.supportArm.dy",0.0);
      _pabs->_degraderSupportArmDims.push_back(hly);
      double hlz = _config.getDouble("degrader.supportArm.dz",0.0);
      _pabs->_degraderSupportArmDims.push_back(hlz);
      // Now get the offsets from the pivot point (in x&y) and
      // from the z0 of the filter above in z.
      double osx = _config.getDouble("degrader.supportArm.offsetx",0.0);
      _pabs->_degraderSupportArmDims.push_back(osx);
      double osy = _config.getDouble("degrader.supportArm.offsety",0.0);
      _pabs->_degraderSupportArmDims.push_back(osy);
      double osz = _config.getDouble("degrader.supportArm.offsetz",0.0);
      _pabs->_degraderSupportArmDims.push_back(osz);
      // The plate will be at the end of the arms.  Give its dimension
      // in x, y, and z
      double hlpx = _config.getDouble("degrader.supportPlate.dx",0.0);
      _pabs->_degraderSupportPlateDims.push_back(hlpx);
      double hlpy = _config.getDouble("degrader.supportPlate.dy",0.0);
      _pabs->_degraderSupportPlateDims.push_back(hlpy);
      double hlpz = _config.getDouble("degrader.supportPlate.dz",0.0);
      _pabs->_degraderSupportPlateDims.push_back(hlpz);

      if (_pabs->_degraderVersion >= 3) {
//-----------------------------------------------------------------------------
// P.Murat: parameters of the second filter disk and the converter
//-----------------------------------------------------------------------------
        double rin2 = _config.getDouble("degrader.filter2.rIn");
        _pabs->_degraderFilter2Dims.push_back(rin2);
        double rout2 = _config.getDouble("degrader.filter2.rOut");
        _pabs->_degraderFilter2Dims.push_back(rout2);
        double hl2  = _config.getDouble("degrader.filter2.halfLength");
        _pabs->_degraderFilter2Dims.push_back(hl2);
//-----------------------------------------------------------------------------
// parameters of the second filter disk and the converter
// to minimize the energy losses may want to move the converter ring a bit forward
//-----------------------------------------------------------------------------
        double rin_conv  = _config.getDouble("degrader.converter.rIn");
        _pabs->_degraderConverterDims.push_back(rin_conv);
        double rout_conv = _config.getDouble("degrader.converter.rOut");
        _pabs->_degraderConverterDims.push_back(rout_conv);
        double hl_conv   = _config.getDouble("degrader.converter.halfLength");
        _pabs->_degraderConverterDims.push_back(hl_conv);
        double dz_conv   = _config.getDouble("degrader.converter.dz");
        _pabs->_degraderConverterDz = dz_conv;
      }

    } // end of the if for Pion Degrader building


    ////////////////////
    // End Rings for IPA

    // defaults chosen so that there are no end-rings for v1 of the IPA
    _pabs->_ipaSupport->_nEndRings = _config.getInt("protonabsorber.ipa.nEndRings", 0);
    const double endRingHalfLength = _config.getDouble("protonabsorber.ipa.endRingHalfLength", 0);
    const double endRingRadialLength = _config.getDouble("protonabsorber.ipa.endRingRadialLength", 0);
    const string endRingMaterial = _config.getString("protonabsorber.ipa.endRingMaterial", "DSVacuum");

    for (std::size_t iRing(0); iRing < _pabs->_ipaSupport->nEndRings(); iRing++) {
      const double zPosition = ipazstart + zSpacing*(iRing); // same z spacing as for the wire sets

      const double endRingOuterRadius = pabs1rIn0 + ( zPosition-ipazstart )*(pabs1rIn1-pabs1rIn0)/(2*pabs1halflen) - 0.01; // have the end ring on the inside of the IPA
      const double endRingInnerRadius = endRingOuterRadius - endRingRadialLength;

      _pabs->_ipaSupport->_endRingMap.push_back( Tube( endRingInnerRadius,
                                                       endRingOuterRadius,
                                                       endRingHalfLength,
                                                       CLHEP::Hep3Vector( -3904, 0., zPosition+endRingHalfLength ), // want the front end of the ring to be at the front of the IPA
                                                       CLHEP::HepRotation(), // put in identiy matrix and determine rotation later
                                                       0,
                                                       CLHEP::twopi,
                                                       endRingMaterial
                                                       )
                                                 );
    }

    /////////////////
    // Wires for IPA

    for ( std::size_t iS(0); iS < _pabs->_ipaSupport->nSets() ; iS++ ) {
      std::vector<Tube> wireVector;

      double zPosition = 0;
      if (_IPAVersion == 1) {
        zPosition = ipazstart + zSpacing*(iS+1);
      }
      else if (_IPAVersion >= 2) {
        zPosition = ipazstart + zSpacing*(iS);
      }

      const double wireOuterRadius = oPAin0 + ( zPosition-opazstart )*(oPAin1-oPAin0)/(2*_pabs->_oPAhalflength);
      const double wireInnerRadius = pabs1rOut0 + ( zPosition-ipazstart )*(pabs1rOut1-pabs1rOut0)/(2*pabs1halflen);
      double wireLength      = wireOuterRadius - wireInnerRadius;

      if (_IPAVersion >= 2) { // for v2, we want the wire to be angled from the vertical in order to provide longitudinal tension
        //        wireLength = wireLength / cos(wire_rotation_from_vertical);
        // However, because the OPA is conical the wires at one end need to be longer than the wires at the other
        double theta_opa = atan2(oPAin1 - oPAin0, 2*_pabs->_oPAhalflength);
        //if OPA entirely covers the region of the wires (as it should) -->
        //tan(angle of OPA descent) = delta r / delta z
        //r at OPA where wire connects = r at OPA at IPA edge / cos(wire angle from horizontal) / (tan(wire angle from horizontal) + tan(OPA angle))

        double wire_angle = wire_rotation_from_vertical;
        if(_IPAVersion > 2) {
          wire_angle *= CLHEP::degree; //fix bug in newer versions
          wireLength -= (0.1 + 2.*abs(wireRadius*tan(wire_angle))); //add buffer + wire radius effect
        }
        if (zPosition < ipazstart+pabs1halflen) {
          // if we're closer to the target
            wireLength = (wireLength * sin(CLHEP::pi - 2*wire_angle - theta_opa)) / sin(wire_angle + theta_opa);
        }
        else {
          // we're further from the target so the wire needs to be longer
            wireLength = (wireLength * sin(CLHEP::pi -2*wire_angle - theta_opa)) / sin(wire_angle - theta_opa);
        }
      }


      for ( std::size_t iW(0); iW < _pabs->_ipaSupport->nWiresPerSet() ; iW++ ) {
        wireVector.push_back( Tube( 0. ,
                                    wireRadius ,
                                    wireLength*0.5,
                                    CLHEP::Hep3Vector( -3904, 0., zPosition ),
                                    CLHEP::HepRotation(), // put in identiy matrix and determine rotation later
                                    0,
                                    CLHEP::twopi,
                                    wireMaterial
                                    )
                              );
      }
      _pabs->_ipaSupport->_supportWireMap.push_back( wireVector );
    }
  }


  void ProtonAbsorberMaker::PrintConfig ( ) {

    double pabs1z = (_pabs->part(0)).center().z();
    double pabs2z = (_pabs->part(1)).center().z();
    double pabs1hl = (_pabs->part(0)).halfLength();
    double pabs2hl = (_pabs->part(1)).halfLength();
    std::cout<<"ProtonAbsorberMaker Configuration -----------------"<<std::endl;
    std::cout<<" Dist from target end : " << _pabs->distanceFromTargetEnd() << std::endl;
    std::cout<<" vdHL : " << _pabs->virtualDetectorHalfLength() << std::endl;
    std::cout<<" pabs1len (full length) : " << pabs1hl*2. << std::endl;
    std::cout<<" pabs2len (full length) : " << pabs2hl*2. << std::endl;
    std::cout<<" pabs1 extent in Mu2e : " << pabs1z - pabs1hl << ", " << pabs1z << ", " << pabs1z+pabs1hl << std::endl;
    std::cout<<" pabs2 extent in Mu2e : " << pabs2z - pabs2hl << ", " << pabs2z << ", " << pabs2z+pabs2hl << std::endl;
    std::cout<<" [rIn, rOut] : " << std::endl;
    std::cout<<" =  [" << _pabs->part(0).innerRadiusAtStart() << ", " << _pabs->part(0).outerRadiusAtStart() <<"], " << std::endl
             <<"    [" << _pabs->part(0).innerRadiusAtEnd()   << ", " << _pabs->part(0).outerRadiusAtEnd()   <<"], " << std::endl
             <<"    [" << _pabs->part(1).innerRadiusAtStart() << ", " << _pabs->part(1).outerRadiusAtStart() <<"], " << std::endl
             <<"    [" << _pabs->part(1).innerRadiusAtEnd()   << ", " << _pabs->part(1).outerRadiusAtEnd()   <<"]"   << std::endl;
    std::cout<<" Material : " << _pabs->fillMaterial() <<  std::endl;
    std::cout<<" pabs1 is " << ( (_pabs->isAvailable(0)) ? "valid" : "invalid" ) << ", " ;
    std::cout<<" pabs2 is " << ( (_pabs->isAvailable(1)) ? "valid" : "invalid" ) << std::endl ;

    if (_pabs->_oPA1flag) {
      std::cout<<" oPAlen (full length) : " << _pabs->_oPAhalflength*2. << std::endl;
      std::cout<<" oPA extent in Mu2e : " << _pabs->_oPAzcenter - _pabs->_oPAhalflength << ", "
                                          << _pabs->_oPAzcenter << ", "
                                          << _pabs->_oPAzcenter + _pabs->_oPAhalflength << std::endl;
      std::cout<<" [rIn, rOut] : " << std::endl;
      std::cout<<" =  [" << _pabs->part(2).innerRadiusAtStart() << ", " << _pabs->part(2).outerRadiusAtStart() <<"], " << std::endl
               <<"    [" << _pabs->part(2).innerRadiusAtEnd()   << ", " << _pabs->part(2).outerRadiusAtEnd()   <<"], " << std::endl;
      std::cout<<" Material : " << _pabs->_oPAmaterialName <<  std::endl;
    }
    else {
      std::cout<<" oPA is invalid" << std::endl ;
    }

  }



} // namespace mu2e
