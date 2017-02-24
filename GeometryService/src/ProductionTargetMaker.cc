#include "GeometryService/inc/ProductionTargetMaker.hh"

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "ProductionTargetGeom/inc/ProductionTarget.hh"

#include <iostream>

namespace mu2e {

  std::unique_ptr<ProductionTarget> ProductionTargetMaker::make(const SimpleConfig& c, double solenoidOffset) {

    std::unique_ptr<ProductionTarget> tgtPS
      (new ProductionTarget(c.getDouble("targetPS_rOut"),
                            c.getDouble("targetPS_halfLength"),
                            c.getDouble("targetPS_rotX") * CLHEP::degree,
                            c.getDouble("targetPS_rotY") * CLHEP::degree,
                            CLHEP::Hep3Vector(solenoidOffset,
                                              0,
                                              c.getDouble("productionTarget.zNominal")
                                              )
                            + c.getHep3Vector("productionTarget.offset", CLHEP::Hep3Vector(0,0,0))
                            )
       );

    double trgtMaxAngle = c.getDouble("targetPS_rotY");
    if (c.getDouble("targetPS_rotX")>trgtMaxAngle) { trgtMaxAngle=c.getDouble("targetPS_rotX"); }
    trgtMaxAngle *= CLHEP::deg;

    int    nSpokeperside = c.getInt("targetPS_Spoke_nsperside");
    double spokeSideDangle = c.getDouble("targetPS_Spoke_sideDangle");
    double spokeAnchordist = c.getDouble("targetPS_Spoke_anchordist");
    double spokeAngleStep = 360.0/((double) nSpokeperside);

    double Hub_thickness = c.getDouble("targetPS_Hub_thickness",0.0);
    double Hub_hang_Length = c.getDouble("targetPS_Hub_hang_Length",0.0);
    double Hub_overhang_Length = c.getDouble("targetPS_Hub_overhang_Length",0.0);
    double Hub_overhang_angle = c.getDouble("targetPS_Hub_overhang_angle",0.0);
    std::vector<double> HubRgtCornersZ, HubRgtCornersInnRadii, HubRgtCornersOutRadii;
    std::vector<double> HubLftCornersZ, HubLftCornersInnRadii, HubLftCornersOutRadii;

    //double Hub_hang_TotLength = Hub_hang_Length;
    if (Hub_overhang_angle>=0.0 && Hub_overhang_angle<90.0) {
            Hub_overhang_angle *= CLHEP::deg;
            double deltaRad = Hub_overhang_Length*tan(Hub_overhang_angle);
            double cosHub_overhang_angle = cos(Hub_overhang_angle);
            double appThick = Hub_thickness/cosHub_overhang_angle;
            //Hub_hang_TotLength+=Hub_overhang_Length;
            HubRgtCornersZ.push_back(tgtPS->halfLength()+Hub_overhang_Length);
            HubRgtCornersInnRadii.push_back(tgtPS->rOut()+deltaRad);
            HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+appThick);

            HubRgtCornersZ.push_back(tgtPS->halfLength());
            HubRgtCornersInnRadii.push_back(tgtPS->rOut());
            HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);

            HubRgtCornersZ.push_back(tgtPS->halfLength()-Hub_hang_Length);
            HubRgtCornersInnRadii.push_back(tgtPS->rOut());
            HubRgtCornersOutRadii.push_back(HubRgtCornersInnRadii.back()+Hub_thickness);

            HubLftCornersZ.push_back(-tgtPS->halfLength()-Hub_overhang_Length);
            HubLftCornersInnRadii.push_back(tgtPS->rOut()+deltaRad);
            HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+appThick);

            HubLftCornersZ.push_back(-tgtPS->halfLength());
            HubLftCornersInnRadii.push_back(tgtPS->rOut());
            HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);

            HubLftCornersZ.push_back(-tgtPS->halfLength()+Hub_hang_Length);
            HubLftCornersInnRadii.push_back(tgtPS->rOut());
            HubLftCornersOutRadii.push_back(HubLftCornersInnRadii.back()+Hub_thickness);

            double SpkAnchrPosZ = HubRgtCornersZ.at(0)-spokeAnchordist*cosHub_overhang_angle;
            double SpkAnchrPosRad = HubRgtCornersOutRadii.at(0)-spokeAnchordist*sin(Hub_overhang_angle);
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

            tgtPS->_envelHalfLength += Hub_overhang_Length;

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
  } // make()

} // namespace mu2e
