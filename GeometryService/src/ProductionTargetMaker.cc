#include "GeometryService/inc/ProductionTargetMaker.hh"
#include "cetlib_except/exception.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include <iostream>
#include <algorithm>

namespace mu2e {

  std::unique_ptr<ProductionTarget> ProductionTargetMaker::make(const SimpleConfig& c, double solenoidOffset) {


    if (c.getString("targetPS_model") == "MDC2018"){ 
      //     std::cout << "making Tier1 in maker" << std::endl;
      return std::move(makeTier1(c, solenoidOffset));
	} else 
      if (c.getString("targetPS_model") == "Hayman_v_2_0"){ 
	//	std::cout << " making Hayman in Maker" << std::endl;
	return std::move(makeHayman_v_2_0(c, solenoidOffset));
	  } else 
	{throw cet::exception("GEOM") << " illegal production target version specified = " << c.getInt("targetPS_version")  << std::endl;}
    return 0;

  }

  std::unique_ptr<ProductionTarget> ProductionTargetMaker::makeTier1(const SimpleConfig& c, double solenoidOffset){
    std::unique_ptr<ProductionTarget> tgtPS
      (new ProductionTarget(
			    c.getString("targetPS_model","NULL"),
			    c.getInt("targetPS_version",1),
			    c.getDouble("targetPS_rOut"),
                            c.getDouble("targetPS_halfLength"),
                            c.getDouble("targetPS_rotX") * CLHEP::degree,
                            c.getDouble("targetPS_rotY") * CLHEP::degree,
                            CLHEP::Hep3Vector(solenoidOffset,
                                              0,
                                              c.getDouble("productionTarget.zNominal")
                                              )
                            + c.getHep3Vector("productionTarget.offset", CLHEP::Hep3Vector(0,0,0)),
                            
			    c.getInt   ("targetPS_nFins"    ,    0  ),
			    c.getDouble("targetPS_finHeight",    0.0),
			    c.getDouble("targetPS_finThickness", 0.0),
			    c.getDouble("targetPS_hubDistanceUS",0.0),
			    c.getDouble("targetPS_hubDistanceDS",0.0),
			    c.getDouble("targetPS_hubAngleUS"   ,0.0),
			    c.getDouble("targetPS_hubAngleDS"   ,0.0),
			    c.getDouble("targetPS_hubOverhangUS",0.0),
			    c.getDouble("targetPS_hubOverhangDS",0.0) ) );
   double trgtMaxAngle = c.getDouble("targetPS_rotY");
    if (c.getDouble("targetPS_rotX")>trgtMaxAngle) { trgtMaxAngle=c.getDouble("targetPS_rotX"); }
    trgtMaxAngle *= CLHEP::deg;

    int    version = tgtPS->version();
    int    nSpokeperside = c.getInt("targetPS_Spoke_nsperside");
    double spokeSideDangle = c.getDouble("targetPS_Spoke_sideDangle");
    double spokeAnchordist = c.getDouble("targetPS_Spoke_anchordist");
    double spokeAngleStep = 360.0/((double) nSpokeperside);

    double Hub_thickness       = c.getDouble("targetPS_Hub_thickness",0.0);
    double Hub_hang_Length     = c.getDouble("targetPS_Hub_hang_Length",0.0);
    double Hub_overhang_Length = c.getDouble("targetPS_Hub_overhang_Length",0.0);
    double Hub_overhang_angle  = c.getDouble("targetPS_Hub_overhang_angle",0.0);

    if (version > 1) {
      Hub_overhang_Length = tgtPS->hubOverhangDS();
      Hub_overhang_angle  = tgtPS->hubAngleDS();
    }
    std::vector<double> HubRgtCornersZ, HubRgtCornersInnRadii, HubRgtCornersOutRadii;
    std::vector<double> HubLftCornersZ, HubLftCornersInnRadii, HubLftCornersOutRadii;

    //double Hub_hang_TotLength = Hub_hang_Length;
    if (Hub_overhang_angle>=0.0 && Hub_overhang_angle<90.0) {
            Hub_overhang_angle *= CLHEP::deg;
            double deltaRad = Hub_overhang_Length*tan(Hub_overhang_angle);
            double cosHub_overhang_angle = cos(Hub_overhang_angle);
            double appThick = Hub_thickness/cosHub_overhang_angle;
	    double actual_hang_Length = Hub_hang_Length;
	    if ( version > 1 ) actual_hang_Length = Hub_thickness/tan(Hub_overhang_angle);
	    tgtPS->setHubLenDS(actual_hang_Length);
            //Hub_hang_TotLength+=Hub_overhang_Length;
            HubRgtCornersZ.push_back(tgtPS->halfLength()+Hub_overhang_Length
				     - tgtPS->hubDistDS());
            HubRgtCornersInnRadii.push_back(tgtPS->rOut()+deltaRad);
            HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+appThick);

            HubRgtCornersZ.push_back(tgtPS->halfLength()-tgtPS->hubDistDS());
            HubRgtCornersInnRadii.push_back(tgtPS->rOut());
	    HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);


            HubRgtCornersZ.push_back(tgtPS->halfLength() - actual_hang_Length 
				     - tgtPS->hubDistDS());
            HubRgtCornersInnRadii.push_back(tgtPS->rOut());

	    if ( version > 1 ) {
	      HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+0.001);
	    } else {
	      HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);
	    }
	    // Finished downstream/right, now do upstream/left
	    if (version > 1) {
	      Hub_overhang_Length = tgtPS->hubOverhangUS();
	      Hub_overhang_angle  = tgtPS->hubAngleUS();
	      Hub_overhang_angle *= CLHEP::deg;
	      deltaRad = Hub_overhang_Length*tan(Hub_overhang_angle);
	      cosHub_overhang_angle = cos(Hub_overhang_angle);
	      appThick = Hub_thickness/cosHub_overhang_angle;
	      actual_hang_Length = Hub_hang_Length;
	      if ( version > 1 ) actual_hang_Length = Hub_thickness/tan(Hub_overhang_angle);
	      tgtPS->setHubLenUS(actual_hang_Length);
	    }

            HubLftCornersZ.push_back(-tgtPS->halfLength()-Hub_overhang_Length
				     + tgtPS->hubDistUS());
            HubLftCornersInnRadii.push_back(tgtPS->rOut()+deltaRad);
            HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+appThick);

            HubLftCornersZ.push_back(-tgtPS->halfLength() + tgtPS->hubDistUS());
            HubLftCornersInnRadii.push_back(tgtPS->rOut());
            HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);

            HubLftCornersZ.push_back(-tgtPS->halfLength()+actual_hang_Length 
				     + tgtPS->hubDistUS());
            HubLftCornersInnRadii.push_back(tgtPS->rOut());
	    if (version > 1 ) {
	      HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+0.1);
	    } else {
	      HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);
	    }

	    // Created both hubs


            double SpkAnchrPosZ = HubRgtCornersZ.at(0)-spokeAnchordist*cosHub_overhang_angle;
            double SpkAnchrPosRad = HubRgtCornersOutRadii.at(0)-spokeAnchordist*sin(Hub_overhang_angle);
            double iSpkAngleRgt(0.0), iSpkAngleLft(-spokeSideDangle), tmpAngle(0.0);
            for (int iSpk=0; iSpk<nSpokeperside; ++iSpk) {
	      // Switch back to Downstream version of #s
	      // if version 2 or higher...
	      if (version > 1) {
		Hub_overhang_Length = tgtPS->hubOverhangDS();
		Hub_overhang_angle  = tgtPS->hubAngleDS();
		Hub_overhang_angle *= CLHEP::deg;
		deltaRad = Hub_overhang_Length*tan(Hub_overhang_angle);
		cosHub_overhang_angle = cos(Hub_overhang_angle);
		appThick = Hub_thickness/cosHub_overhang_angle;
		SpkAnchrPosZ = HubRgtCornersZ.at(0)-spokeAnchordist*cosHub_overhang_angle;
		SpkAnchrPosRad = HubRgtCornersOutRadii.at(0)-spokeAnchordist*sin(Hub_overhang_angle);
	      }
                    tmpAngle = iSpkAngleRgt * CLHEP::deg;
                    std::map<double,CLHEP::Hep3Vector>::iterator pntPos;
                    pntPos = tgtPS->_anchoringPntsRgt.insert(
                                    std::pair<double,CLHEP::Hep3Vector>(
                                                    iSpkAngleRgt,
                                                    CLHEP::Hep3Vector(SpkAnchrPosRad*cos(tmpAngle),
                                                                    SpkAnchrPosRad*sin(tmpAngle),
                                                                    SpkAnchrPosZ)
                                    )
                             ).first;
                    pntPos->second.transform(tgtPS->protonBeamRotation());
		    // Switch back to Upstream version of #s
		    // if version 2 or higher...
		    if (version > 1) {
		      Hub_overhang_Length = tgtPS->hubOverhangUS();
		      Hub_overhang_angle  = tgtPS->hubAngleUS();
		      Hub_overhang_angle *= CLHEP::deg;
		      deltaRad = Hub_overhang_Length*tan(Hub_overhang_angle);
		      cosHub_overhang_angle = cos(Hub_overhang_angle);
		      appThick = Hub_thickness/cosHub_overhang_angle;
		      SpkAnchrPosZ = HubLftCornersZ.at(0)+spokeAnchordist*cosHub_overhang_angle;
		      SpkAnchrPosZ *= -1.0;
		      SpkAnchrPosRad = HubLftCornersOutRadii.at(0)-spokeAnchordist*sin(Hub_overhang_angle);
	      }
                    tmpAngle = iSpkAngleLft * CLHEP::deg;
                    pntPos = tgtPS->_anchoringPntsLft.insert(
                                    std::pair<double,CLHEP::Hep3Vector>(
                                                    iSpkAngleLft,
                                                    CLHEP::Hep3Vector(SpkAnchrPosRad*cos(tmpAngle),
                                                                    SpkAnchrPosRad*sin(tmpAngle),
                                                                    -SpkAnchrPosZ)
                                    )
                             ).first;
                    pntPos->second.transform(tgtPS->protonBeamRotation());
                    iSpkAngleRgt += spokeAngleStep;
                    iSpkAngleLft += spokeAngleStep;
            }

	    if ( version > 1 ) {
	      double maxHang = tgtPS->hubOverhangUS();
	      double tmpHang = tgtPS->hubOverhangDS();
	      if ( tmpHang > maxHang ) maxHang = tmpHang;
	      tgtPS->_envelHalfLength += maxHang + 2.0;
	    } else {
	      tgtPS->_envelHalfLength += Hub_overhang_Length;
	    }

    } else if (Hub_overhang_angle==90.0) {

           HubRgtCornersZ.push_back(tgtPS->halfLength());
           HubRgtCornersInnRadii.push_back(tgtPS->rOut());
           HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness+Hub_overhang_Length);

           HubRgtCornersZ.push_back(tgtPS->halfLength()-Hub_thickness);
           HubRgtCornersInnRadii.push_back(tgtPS->rOut());
           HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness+Hub_overhang_Length);

           HubRgtCornersZ.push_back(tgtPS->halfLength()-(Hub_thickness+0.000000001));
           HubRgtCornersInnRadii.push_back(tgtPS->rOut());
           HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);

           HubRgtCornersZ.push_back(tgtPS->halfLength()-Hub_hang_Length);
           HubRgtCornersInnRadii.push_back(tgtPS->rOut());
           HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);

           HubLftCornersZ.push_back(-tgtPS->halfLength());
           HubLftCornersInnRadii.push_back(tgtPS->rOut());
           HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness+Hub_overhang_Length);

           HubLftCornersZ.push_back(-tgtPS->halfLength()+Hub_thickness);
           HubLftCornersInnRadii.push_back(tgtPS->rOut());
           HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness+Hub_overhang_Length);

           HubLftCornersZ.push_back(-tgtPS->halfLength()+(Hub_thickness+0.000000001));
           HubLftCornersInnRadii.push_back(tgtPS->rOut());
           HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);

           HubLftCornersZ.push_back(-tgtPS->halfLength()+Hub_hang_Length);
           HubLftCornersInnRadii.push_back(tgtPS->rOut());
           HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);

           double SpkAnchrPosZ = HubRgtCornersZ.at(1);
           double SpkAnchrPosRad = HubRgtCornersOutRadii.at(1)-spokeAnchordist;
           double iSpkAngleRgt(0.0), iSpkAngleLft(-spokeSideDangle), tmpAngle(0.0);
           for (int iSpk=0; iSpk<nSpokeperside; ++iSpk) {
                   tmpAngle = iSpkAngleRgt * CLHEP::deg;
                   std::map<double,CLHEP::Hep3Vector>::iterator pntPos;
                   pntPos = tgtPS->_anchoringPntsRgt.insert(
                                   std::pair<double,CLHEP::Hep3Vector>(
                                                   iSpkAngleRgt,
                                                   CLHEP::Hep3Vector(SpkAnchrPosRad*cos(tmpAngle),
                                                                   SpkAnchrPosRad*sin(tmpAngle),
                                                                   SpkAnchrPosZ)
                                   )
                            ).first;
                   pntPos->second.transform(tgtPS->protonBeamRotation());
                   tmpAngle = iSpkAngleLft * CLHEP::deg;
                   pntPos = tgtPS->_anchoringPntsLft.insert(
                                   std::pair<double,CLHEP::Hep3Vector>(
                                                   iSpkAngleLft,
                                                   CLHEP::Hep3Vector(SpkAnchrPosRad*cos(tmpAngle),
                                                                   SpkAnchrPosRad*sin(tmpAngle),
                                                                   -SpkAnchrPosZ)
                                   )
                            ).first;
                   pntPos->second.transform(tgtPS->protonBeamRotation());
                   iSpkAngleRgt += spokeAngleStep;
                   iSpkAngleLft += spokeAngleStep;
           }

    } else {
            throw cet::exception("GEOM")
            << "Production Target, angle of the Hub overhang is not implemented \n";
    }

    tgtPS->_envelHalfLength *= cos(trgtMaxAngle);
    tgtPS->_envelHalfLength += HubRgtCornersOutRadii.at(0)*sin(trgtMaxAngle);

    bool   const _has_virtualDet              = c.getBool("targetPS.hasVD.backward",false) || c.getBool("targetPS.hasVD.forward",false);
    double const _sideVD_thickness            = 2.0*c.getDouble("vd.halfLength",0.0);
    if (_has_virtualDet) { tgtPS->_envelHalfLength+=_sideVD_thickness; }

    //const CLHEP::Hep3Vector& tgtPS_pos = tgtPS->position();
    const CLHEP::Hep3Vector center(0,0,0);
//     tgtPS->_pHubsRgtParams = std::unique_ptr<Polycone>
//       (new Polycone(HubRgtCornersZ,
//                       HubRgtCornersInnRadii,
//                       HubRgtCornersOutRadii,
//                       center,
//                       c.getString("targetPS_Hub_materialName")));
//
//     tgtPS->_pHubsLftParams = std::unique_ptr<Polycone>
//       (new Polycone(HubLftCornersZ,
//                       HubLftCornersInnRadii,
//                       HubLftCornersOutRadii,
//                       center,
//                       c.getString("targetPS_Hub_materialName")));
    tgtPS->_pHubsRgtParams = new Polycone(HubRgtCornersZ,
                      HubRgtCornersInnRadii,
                      HubRgtCornersOutRadii,
                      center,
                      c.getString("targetPS_Hub_materialName"));

    tgtPS->_pHubsLftParams = new Polycone(HubLftCornersZ,
                      HubLftCornersInnRadii,
                      HubLftCornersOutRadii,
                      center,
                      c.getString("targetPS_Hub_materialName"));


    return std::move(tgtPS);
  }

  std::unique_ptr<ProductionTarget> ProductionTargetMaker::makeHayman_v_2_0(const SimpleConfig& c, double solenoidOffset){

    //simple config does not have a method for returning a vector: all voids.  So I get it, stuff it in a 
    //temporary variable, and put that into the constructor.

    std::vector<double> startingSectionThickness;
    std::vector<int>    numberOfSegmentsPerSection;
    std::vector<double> thicknessOfSegmentPerSection;
    std::vector<double> heightOfRectangularGapPerSection;
    std::vector<double> thicknessOfGapPerSection;
    std::vector<double> finAngles;
    c.getVectorDouble("targetPS_startingSectionThickness",startingSectionThickness);
    c.getVectorInt("targetPS_numberOfSegmentsPerSection",numberOfSegmentsPerSection);
    c.getVectorDouble("targetPS_thicknessOfSegmentPerSection",thicknessOfSegmentPerSection);
    c.getVectorDouble("targetPS_heightOfRectangularGapPerSection",heightOfRectangularGapPerSection);
    c.getVectorDouble("targetPS_thicknessOfGapPerSection",thicknessOfGapPerSection);
    c.getVectorDouble("targetPS_finAngles",finAngles);
    for_each (finAngles.begin(),finAngles.end(),[]( double& elem){elem *= CLHEP::degree;});
    std::unique_ptr<ProductionTarget> tgtPS
      (new ProductionTarget(
			    c.getString("targetPS_model","NULL"),
			    c.getInt("targetPS_version"),
			    c.getDouble("targetPS_productionTargetMotherOuterRadius"),
			    c.getDouble("targetPS_productionTargetMotherHalfLength"),
			    c.getDouble("targetPS_rOut"),
                            c.getDouble("targetPS_halfHaymanLength"),
                            c.getDouble("targetPS_rotX") * CLHEP::degree,
                            c.getDouble("targetPS_rotY") * CLHEP::degree,
                            c.getDouble("targetPS_rotZ") * CLHEP::degree,
                            CLHEP::Hep3Vector(solenoidOffset,
                                              0,
                                              c.getDouble("productionTarget.zNominal")
                                              )
                            + c.getHep3Vector("productionTarget.offset"),
                            c.getString("targetPS_targetCoreMaterial"),
			    c.getString("targetPS_targetFinMaterial"),
			    c.getString("targetPS_targetVacuumMaterial"),
			    c.getString("targetPS_supportRingMaterial"),
			    c.getString("targetPS_spokeMaterial"),
			    c.getInt("targetPS_numberOfTargetSections"),
			    startingSectionThickness,
                            numberOfSegmentsPerSection,
			    thicknessOfSegmentPerSection,
			    heightOfRectangularGapPerSection,
			    thicknessOfGapPerSection,
			    c.getInt("targetPS_nHaymanFins"),
			    finAngles,
			    c.getDouble("targetPS_finThickness"),
			    c.getDouble("targetPS_finOuterRadius"),
			    c.getDouble("targetPS_supportRingLength"),
                            c.getDouble("targetPS_supportRingInnerRadius"),
                            c.getDouble("targetPS_supportRingOuterRadius"),
                            c.getDouble("targetPS_supportRingCutoutThickness"),
                            c.getDouble("targetPS_supportRingCutoutLength")
			    ));
    //check if we should configure supports
    //switch to using '.' instead of '_' to differentiate levels
    tgtPS->_supportsBuild = c.getBool("targetPS.supports.build", false);
    if(tgtPS->_supportsBuild) {
      //support wheel parameters
      tgtPS->_supportWheelRIn      = c.getDouble("targetPS.supports.wheel.rIn");
      tgtPS->_supportWheelROut     = c.getDouble("targetPS.supports.wheel.rOut");
      tgtPS->_supportWheelHL       = c.getDouble("targetPS.supports.wheel.halfLength");
      tgtPS->_supportWheelMaterial = c.getString("targetPS.supports.wheel.material");
      //number of support rods and wires
      tgtPS->_nSpokesPerSide = c.getInt("targetPS.supports.nSpokes");
      //features on the wheel near each support rod
      c.getVectorDouble("targetPS.supports.features.angles", tgtPS->_supportWheelFeatureAngles, tgtPS->_nSpokesPerSide); 
      c.getVectorDouble("targetPS.supports.features.arcs"  , tgtPS->_supportWheelFeatureArcs  , tgtPS->_nSpokesPerSide); 
      c.getVectorDouble("targetPS.supports.features.rIns"  , tgtPS->_supportWheelFeatureRIns  , tgtPS->_nSpokesPerSide); 
      //support wheel rods parameters
      c.getVectorDouble("targetPS.supports.rods.halfLength", tgtPS->_supportWheelRodHL, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.offset", tgtPS->_supportWheelRodOffset, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.radius", tgtPS->_supportWheelRodRadius, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.radialOffset", tgtPS->_supportWheelRodRadialOffset, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.wireOffset.downstream", tgtPS->_supportWheelRodWireOffsetD, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.wireOffset.upstream", tgtPS->_supportWheelRodWireOffsetU, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.angles", tgtPS->_supportWheelRodAngles, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.rods.tilts", tgtPS->_supportWheelRodTilts, tgtPS->_nSpokesPerSide);
      //support wire (spokes) parameters
      c.getVectorDouble("targetPS.supports.spokes.targetAngles.downstream", tgtPS->_spokeTargetAnglesD, tgtPS->_nSpokesPerSide);
      c.getVectorDouble("targetPS.supports.spokes.targetAngles.upstream", tgtPS->_spokeTargetAnglesU, tgtPS->_nSpokesPerSide);
      tgtPS->_spokeRadius = 0.5*c.getDouble("targetPS.supports.spokes.diameter");
      //override old format of the material if new syntax is found
      tgtPS->_spokeMaterial = c.getString("targetPS.supports.spokes.material", tgtPS->_spokeMaterial);
      //loop through support rods, update mother outer radius if they extend beyond it
      double rMother = tgtPS->_productionTargetMotherOuterRadius;
      for(int ispoke = 0; ispoke < tgtPS->_nSpokesPerSide; ++ispoke) {
	double rodR0   = tgtPS->_supportWheelRodRadialOffset[ispoke];
	double rodR    = tgtPS->_supportWheelRodRadius[ispoke];
	double rodTilt = tgtPS->_supportWheelRodTilts[ispoke]*CLHEP::degree;
	double rodHL   = tgtPS->_supportWheelRodHL[ispoke];
	double rodDZ   = tgtPS->_supportWheelRodOffset[ispoke];
	//maximum radius is on the side angled away from the center
	int side = (rodTilt > 0.) ? 1 : -1; //downstream side is higher for positive angles
	double rMax = rodR0 + (rodHL + side*rodDZ)*abs(sin(rodTilt)) + rodR*cos(rodTilt);
	if(rMax >= rMother) rMother = rMax + 0.001; //add a small buffer
      }
      tgtPS->_productionTargetMotherOuterRadius = rMother;
    }
    return std::move(tgtPS);
  }
 
 
} // namespace mu2e
