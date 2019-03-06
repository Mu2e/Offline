//
// Construct VirtualDetectors
//
// $Id: VirtualDetectorMaker.cc,v 1.35 2014/09/16 21:57:43 jrquirk Exp $
// $Author: jrquirk $
//

#include <iostream>
#include <iomanip>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "GeometryService/inc/VirtualDetectorMaker.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/Collimator_TS1.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

//#include "G4Helper/inc/G4Helper.hh"
//#include "G4Helper/inc/VolumeInfo.hh"

using namespace std;
using namespace CLHEP;

namespace mu2e {

  std::unique_ptr<VirtualDetector> VirtualDetectorMaker::make(const SimpleConfig& c) {

    unique_ptr<VirtualDetector> vd(new VirtualDetector());

    art::ServiceHandle<GeometryService> geom;

    vd = unique_ptr<VirtualDetector>(new VirtualDetector());

    if(c.getBool("hasVirtualDetector",false)) {

      const double vdHL = c.getDouble("vd.halfLength",0.01*mm);
      vd->_halfLength = vdHL;

      // Add configurable control for verbosity.
      const int verbosityLevel = c.getInt("vd.verbosityLevel",0);

      // Need some data from other subsystems
      GeomHandle<Beamline> bg;
      double solenoidOffset = bg->solenoidOffset();

      // VD Coll1_In and Coll1_Out are at the front and back of
      // collimator 1, which is placed inside TS1.

      double coll1HL = bg->getTS().getColl1().halfLength();
      double coll1ColLen = bg->getTS().getColl1().collarHalfLength() * 2.0;
      const HepRotation *ts1rot = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS1,
                                                        TransportSolenoid::TSRadialPart::IN)->getRotation();
      Hep3Vector   ts1pos = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS1,
                                                  TransportSolenoid::TSRadialPart::IN)->getGlobal();

      Hep3Vector coll1pos = bg->getTS().getColl1().getLocal();
      double zDiff = bg->getTS().getColl1().collarMarginZ();


      Hep3Vector deltaZ1(0,0,coll1HL-vdHL);  
      Hep3Vector deltaZ1Collar(0,0,coll1HL - vdHL - zDiff);
      Hep3Vector deltaZ0Collar(0,0,coll1HL - coll1ColLen - zDiff + vdHL);

      vd->addVirtualDetector( VirtualDetectorId::Coll1_In,
                               ts1pos, ts1rot, coll1pos-deltaZ1);
      vd->addVirtualDetector( VirtualDetectorId::Coll1_Out,
                               ts1pos, ts1rot, coll1pos+deltaZ1);
      vd->addVirtualDetector( VirtualDetectorId::Coll1_pBarCollar_Out,
			      ts1pos, ts1rot, coll1pos+deltaZ1Collar);
      vd->addVirtualDetector( VirtualDetectorId::Coll1_pBarCollar_In,
			      ts1pos, ts1rot, coll1pos+deltaZ0Collar);

      //************************************************************
      // VD TS2_Bend, TS4_Bend are placed at the nominal beamline in
      // the centers of the bends of TS2 and TS4.  DNB (Lou) 17 Jan 2016

      // This gives us the center of curvature of the bend
      Hep3Vector   ts2pos = bg->getTS().getTSCryo(
            TransportSolenoid::TSRegion::TS2,
	    TransportSolenoid::TSRadialPart::IN)->getGlobal();

      // Now displace from center of curvature to nominal beamline
      double rTor = bg->getTS().torusRadius();

      // Because of rotation of the torus 90 degrees, have to switch put
      // desire z-displacement in negative y-direction...
      Hep3Vector displace(rTor*sin(45.0*CLHEP::deg),
			  -rTor*cos(45.0*CLHEP::deg),0.0);

      const HepRotation *ts2rot =
      	bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS2,
      			      TransportSolenoid::TSRadialPart::IN)->getRotation();

      vd->addVirtualDetector( VirtualDetectorId::TS2_Bend,
			      ts2pos, ts2rot, displace);

      Hep3Vector   ts4pos = bg->getTS().getTSCryo(
            TransportSolenoid::TSRegion::TS4,
	    TransportSolenoid::TSRadialPart::IN)->getGlobal();

      const HepRotation *ts4rot =
      	bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS4,
      			      TransportSolenoid::TSRadialPart::IN)->getRotation();

      vd->addVirtualDetector( VirtualDetectorId::TS4_Bend,
			      ts4pos, ts4rot, -displace );

      //***************************************************
      // VD Coll31_In, Coll31_Out, Coll32_In, Coll32_Out are placed
      // around two sections of collimator 3

      double coll31HL = bg->getTS().getColl31().halfLength();
      double coll32HL = bg->getTS().getColl32().halfLength();

      const HepRotation *ts3rot = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS3,
                                                        TransportSolenoid::TSRadialPart::IN)->getRotation();
      Hep3Vector   ts3pos = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS3,
                                                  TransportSolenoid::TSRadialPart::IN)->getGlobal();

      Hep3Vector coll31pos = bg->getTS().getColl31().getLocal();
      Hep3Vector coll32pos = bg->getTS().getColl32().getLocal();
      Hep3Vector deltaZ31(0,0,coll31HL-vdHL);
      Hep3Vector deltaZ32(0,0,coll32HL-vdHL);

      vd->addVirtualDetector( VirtualDetectorId::Coll31_In,
                               ts3pos, ts3rot, coll31pos-deltaZ31);
      vd->addVirtualDetector( VirtualDetectorId::Coll31_Out,
                               ts3pos, ts3rot, coll31pos+deltaZ31);
      vd->addVirtualDetector( VirtualDetectorId::Coll32_In,
                               ts3pos, ts3rot, coll32pos-deltaZ32);
      vd->addVirtualDetector( VirtualDetectorId::Coll32_Out,
                               ts3pos, ts3rot, coll32pos+deltaZ32);

      // VD Coll5_In, Coll5_Out are at the front and back of collimator
      // 5, which is placed inside TS5.

      double coll5HL = bg->getTS().getColl51().halfLength();

      const HepRotation *ts5rot = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS5,
                                                        TransportSolenoid::TSRadialPart::IN)->getRotation();
      Hep3Vector   ts5pos = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS5,
                                                        TransportSolenoid::TSRadialPart::IN)->getGlobal();
      Hep3Vector coll5pos = bg->getTS().getColl51().getLocal();


      Hep3Vector deltaZ5(0,0,coll5HL-vdHL);

      vd->addVirtualDetector( VirtualDetectorId::Coll5_In,
                               ts5pos, ts5rot, coll5pos-deltaZ5);
      vd->addVirtualDetector( VirtualDetectorId::Coll5_Out,
                               ts5pos, ts5rot, coll5pos+deltaZ5);

      // VD ST_In, ST_Out are placed inside DS2, just before and after
      // stopping target

      GeomHandle<StoppingTarget> target;
      GeomHandle<DetectorSolenoid> ds;

      const CLHEP::Hep3Vector ds2centerInMu2e(ds->position().x(), ds->position().y(), ds->vac_zLocDs23Split() - ds->vac_halfLengthDs2());
      const Hep3Vector targetOffset(target->centerInMu2e() - ds2centerInMu2e);
      Hep3Vector shift(0., 0., vdHL + target->cylinderLength()/2);

      vd->addVirtualDetector( VirtualDetectorId::ST_In,
                               ds2centerInMu2e, 0, targetOffset-shift);
      vd->addVirtualDetector( VirtualDetectorId::ST_Out,
                               ds2centerInMu2e, 0, targetOffset+shift);


      /*******new virtual detector for STM Upstream halfway between coll5Out and STIn   ****/

      if ( verbosityLevel > 0 ) {
	std::cout << "coll 5 " << coll5pos.z() << " " << deltaZ5.z() << " " << targetOffset.z() << " " << shift.z() << std::endl;
      }
      const Hep3Vector STMOffset(targetOffset.x()-shift.x(),targetOffset.y()-shift.y(), targetOffset.z()-shift.z() - 0.5*( (coll5pos.z()+deltaZ5.z()) - (targetOffset.z()-shift.z()) ));
      vd->addVirtualDetector( VirtualDetectorId::STMUpstream,
			      ds2centerInMu2e,0,STMOffset);



      if (c.getBool("hasTracker",false)){

        ostringstream vdName(VirtualDetectorId::name(VirtualDetectorId::TT_Mid));

        if(StrawId::_nplanes%2!=0){
          throw cet::exception("GEOM")
            << "This virtual detector " << vdName.str()
            << " can only be placed if the Tracker has an even number of planes \n";
        }

        Tracker const & tracker = *(GeomHandle<Tracker>());
        Hep3Vector ttOffset(-solenoidOffset,0.,tracker.z0());

        // VD TT_Mid is placed inside the tracker mother volume in the
        // middle of the tracker shifted by the half length of vd
        // VD TT_MidInner is placed inside the tracker at the same z position as
        // VD TT_Mid but from radius 0 to the inner radius of the tracker
        // mother volume. However, its mother volume is DS3Vacuum
        // which has a different offset. We will use the global offset
        // here (!) as DS is not in the geometry service yet


        Hep3Vector vdTTMidOffset(0.,0.,0.);
	// Version 4 adds brass rings in Tracker, have to move vd to the side
	if ( c.getBool("TrackerHasBrassRings",false) ) vdTTMidOffset.setZ(10.1);

        vd->addVirtualDetector( VirtualDetectorId::TT_Mid,
                                 ttOffset, 0, vdTTMidOffset);

        vd->addVirtualDetector( VirtualDetectorId::TT_MidInner,
                                 ttOffset, 0, vdTTMidOffset);

        //       if (verbosityLevel >0) {
        //         for ( int vdId=11; vdId<=12; ++vdId) {
        //           cout << __func__ << " VD " << vdId << " offsets L, G " <<
        //             vd->getLocal(vdId) << ", " <<
        //             vd->getGlobal(vdId) << endl;
        //         }
        //       }

        // Global position is in Mu2e coordinates; local position in the detector system.
        double zFrontGlobal = tracker.mother().position().z()-tracker.mother().tubsParams().zHalfLength()-vdHL;
        double zBackGlobal  = tracker.mother().position().z()+tracker.mother().tubsParams().zHalfLength()+vdHL;
        double zFrontLocal  = zFrontGlobal - tracker.z0();
        double zBackLocal   = zBackGlobal  - tracker.z0();

        Hep3Vector vdTTFrontOffset(0.,
                                   0.,
                                   zFrontLocal);

        // VD TT_FrontHollow is placed outside the tracker mother
        // volume in front of the tracker "outside" of the proton
        // absorber


        // formally VD TT_FrontHollow, TT_FrontPA are placed in DS3Vacuum, but it is a
        // complicated subtraction volume, so we pretend to place them in
        // the Tracker and rely on the global offsets in the mu2e
        // detector frame (note that their local offsets are wrt Tracker)

        vd->addVirtualDetector( VirtualDetectorId::TT_FrontHollow,
                                 ttOffset,
                                 0,
                                 vdTTFrontOffset);

        // we add another VD detector at the same Z "inside" the proton absorber

        if (c.getBool("hasProtonAbsorber",false)){

          vd->addVirtualDetector(  VirtualDetectorId::TT_FrontPA,
                                    ttOffset,
                                    0,
                                    vdTTFrontOffset);
        }

        Hep3Vector vdTTBackOffset(0.,
                                  0.,
                                  zBackLocal);

        vd->addVirtualDetector( VirtualDetectorId::TT_Back,
                                 ttOffset,
                                 0,
                                 vdTTBackOffset);

        // these next two detectors are also thin, but they are not disks but cylinders
        // placed on the inner and outer surface of the tracker envelope

        Hep3Vector vdTTOutSurfOffset(0.,0.,tracker.mother().position().z()-tracker.z0());

        vd->addVirtualDetector( VirtualDetectorId::TT_OutSurf,
                                 ttOffset, 0, vdTTOutSurfOffset);

        vd->addVirtualDetector( VirtualDetectorId::TT_InSurf,
                                 ttOffset, 0, vdTTOutSurfOffset);

      }

      if(geom->hasElement<ExtMonFNALBuilding>() && c.getBool("extMonFNAL.filter.vd.enabled", false)) {
        CLHEP::Hep3Vector vzero;

        // This detector will be placed on the face of beam dump
        // shielding.  Computing offsets here is inconvenient since we
        // don't have VolumeInfo for the parent. Just ignore them.
        vd->addVirtualDetector(VirtualDetectorId::EMFC1Entrance, vzero, 0, vzero);

        // Detector inside the ExtMonFNAL magnet room, on the face of the upstream wall
        vd->addVirtualDetector(VirtualDetectorId::EMFC1Exit, vzero, 0, vzero);

        // Detector inside the ExtMonFNAL magnet room, on the face of the downstream wall
        vd->addVirtualDetector(VirtualDetectorId::EMFC2Entrance, vzero, 0, vzero);

        // Detector inside the ExtMonFNAL detector room, on the face of the upstream wall
        vd->addVirtualDetector(VirtualDetectorId::EMFC2Exit, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
      }

      // This VD is related to PS
      if(true) {
        const CLHEP::Hep3Vector vzero(0,0,0);

        // This detector will be placed on the front (-Z direction) exit of PS.
        // Computing offsets in G4 since PS is not in GeomHandle.
        vd->addVirtualDetector(VirtualDetectorId::PS_FrontExit,
                                vzero, 0, vzero
                                );

        if (c.getBool("targetPS.hasVD.backward",false)) {

                GeomHandle<ProductionTarget> tgt;
                const CLHEP::Hep3Vector vBCKWOff(0,0,-tgt->halfLength()+c.getDouble("vd.halfLength",0.0));

                vd->addVirtualDetector(VirtualDetectorId::PT_Back,
                                        vzero, 0, vBCKWOff
                                        );
        }
        if (c.getBool("targetPS.hasVD.forward",false)) {

                GeomHandle<ProductionTarget> tgt;
                const CLHEP::Hep3Vector vFRNWOff(0,0,tgt->halfLength()-c.getDouble("vd.halfLength",0.0));

                vd->addVirtualDetector(VirtualDetectorId::PT_Front,
                                        vzero, 0, vFRNWOff
                                        );
        }

      }

      if(c.hasName("vd.ExtMonCommonPlane.z")) {
        // Position and half length of this detector are best computed
        // in one place.  Since the VirtualDetector data structure
        // does not store half size, we'll do the computations later.
        vd->addVirtualDetector(VirtualDetectorId::ExtMonCommonPlane, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
      }


      if(c.getBool("vd.ProtonBeamDumpCoreFace.enabled", false)) {
        // Position and half length of this detector are best computed
        // in one place.  Since the VirtualDetector data structure
        // does not store half size, we'll do the computations later.
        vd->addVirtualDetector(VirtualDetectorId::ProtonBeamDumpCoreFace, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
      }

      if(c.hasName("extMonFNAL.vd.enabled")) {
        throw cet::exception("BADCONFIG")<<"Error: geometry parameter extMonFNAL.vd.enabled is obsolete and should not be used.\n";
      }
      if(geom->hasElement<ExtMonFNAL::ExtMon>() && c.getBool("extMonFNAL.detector.vd.enabled", false)) {
        vd->addVirtualDetector(VirtualDetectorId::EMFDetectorUpEntrance, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFDetectorUpExit, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFDetectorDnEntrance, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFDetectorDnExit, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
      }

      if(geom->hasElement<ExtMonFNAL::ExtMon>() && c.getBool("extMonFNAL.box.vd.enabled", false)) {
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxFront, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxSW, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxBottom, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxBack, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxNE, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
        vd->addVirtualDetector(VirtualDetectorId::EMFBoxTop, CLHEP::Hep3Vector(), 0, CLHEP::Hep3Vector());
      }

      // VD Coll5_OutSurf is at the outer surfaceof the collimator
      // 5, which is placed inside TS5.

      vd->addVirtualDetector( VirtualDetectorId::Coll5_OutSurf,
                               ts5pos, ts5rot, coll5pos);

      // VD DSNeutronShieldExit is at the downstream part of the
      // aperture in the neutron shielding outside of the IFB/VPSP
      if ( c.getBool("vd.DSNeutronShieldExit.build", false ) ) {

//         const std::vector<CLHEP::Hep3Vector>& ENSCBcentersOfBoxes = enscendb->centersOfBoxes();

//         size_t nBox = ENSCBcentersOfBoxes.size();
//         size_t ib;
//         for(ib = 0; ib < nBox; ++ib) {
//           if ( enscendb->hasHole(ib) ) break;
//         }
//         int hID = enscendb->holeIndex(ib);
//         // locations are wrt HallAir
//         // for some reason the location has to be taken from the box and not the hole tbd
//         //        CLHEP::Hep3Vector holeLocation = enscendb->holeLocation(hID);
//        CLHEP::Hep3Vector holeLocation = ENSCBcentersOfBoxes[ib];
	CLHEP::Hep3Vector holeLocation(
				       c.getDouble("ExtShieldDownstream.detecHoleX")*CLHEP::mm,
				       c.getDouble("ExtShieldDownstream.detecHoleY")*CLHEP::mm,
				       c.getDouble("ExtShieldDownstream.detecHoleZ")*CLHEP::mm);
	double holeHalfLength = c.getDouble("ExtShieldDownstream.detecHoleHalflength")*CLHEP::mm;

	// End of bit added by Dave (Louisville) Brown

        GeomHandle<DetectorSolenoid> ds;
        CLHEP::Hep3Vector const & dsP ( ds->position() );
        CLHEP::Hep3Vector vdPositionInMu2e(dsP.x(), dsP.y(),
                                           holeLocation.z() + holeHalfLength + vd->getHalfLength());

        GeomHandle<Mu2eEnvelope> env;
        const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                       (env->xmax() + env->xmin())/2.,
                                                       (env->ymax() + env->ymin())/2.,
                                                       (env->zmax() + env->zmin())/2.
                                                       );

        vd->addVirtualDetector( VirtualDetectorId::DSNeutronShieldExit,
                                hallFormalCenterInMu2e,
                                0x0,
                                vdPositionInMu2e - hallFormalCenterInMu2e);


        if ( verbosityLevel > 0) {
           cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::DSNeutronShieldExit) << endl;
           cout << "               at local=" << vd->getLocal(VirtualDetectorId::DSNeutronShieldExit) << " global="<< vd->getGlobal(VirtualDetectorId::DSNeutronShieldExit) <<endl;
        }

      }



      if (c.getBool("hasDiskCalorimeter",true))
      {
	GeomHandle<DiskCalorimeter> cg;

	int vdIdDiskSurf = VirtualDetectorId::EMC_Disk_0_SurfIn;
	int vdIdDiskEdge = VirtualDetectorId::EMC_Disk_0_EdgeIn;
        int vdIdFEBEdge  = VirtualDetectorId::EMC_FEB_0_EdgeIn;
        int vdIdFEBSurf  = VirtualDetectorId::EMC_FEB_0_SurfIn;

        double crateHalfLength = cg->caloInfo().getDouble("crateZLength")/2.0;           
	double delta           = 2*vdHL+0.02;

        CLHEP::Hep3Vector parentInMu2e = cg->geomUtil().origin();

	for(size_t i=0; i<cg->nDisk(); ++i)
        {
           const CLHEP::Hep3Vector& sizeDisk = cg->disk(i).geomInfo().size();
           CLHEP::Hep3Vector posDiskLocal  = cg->disk(i).geomInfo().origin() -
           cg->geomUtil().origin();
           CLHEP::Hep3Vector posCrateLocal = posDiskLocal + CLHEP::Hep3Vector(0.0,0.0,cg->disk(i).geomInfo().crateDeltaZ());

           CLHEP::Hep3Vector  posFrontDisk = posDiskLocal - CLHEP::Hep3Vector (0,0,sizeDisk.z()/2.0+delta);
           CLHEP::Hep3Vector  posBackDisk  = posDiskLocal + CLHEP::Hep3Vector (0,0,sizeDisk.z()/2.0+delta);
           CLHEP::Hep3Vector  posInnerDisk = posDiskLocal;

           CLHEP::Hep3Vector  posFrontFEB  = posCrateLocal - CLHEP::Hep3Vector (0,0,crateHalfLength+delta);
           CLHEP::Hep3Vector  posBackFEB   = posCrateLocal + CLHEP::Hep3Vector (0,0,crateHalfLength+delta);
           CLHEP::Hep3Vector  posInnerFEB  = posCrateLocal;

          vd->addVirtualDetector( vdIdDiskSurf,
				  parentInMu2e,
				  0,
				  posFrontDisk);
	  ++vdIdDiskSurf;

	  vd->addVirtualDetector( vdIdDiskSurf,
				  parentInMu2e,
				  0,
				  posBackDisk);
	  ++vdIdDiskSurf;


	  vd->addVirtualDetector( vdIdDiskEdge,
				  parentInMu2e,
				  0,
				  posInnerDisk);
	  ++vdIdDiskEdge;

	  vd->addVirtualDetector( vdIdDiskEdge,
				  parentInMu2e,
				  0,
				  posInnerDisk);
	  ++vdIdDiskEdge;



	  vd->addVirtualDetector( vdIdFEBSurf,
				  parentInMu2e,
				  0,
				  posFrontFEB);
	  ++vdIdFEBSurf;

	  vd->addVirtualDetector( vdIdFEBSurf,
				  parentInMu2e,
				  0,
				  posBackFEB);
	  ++vdIdFEBSurf;


 	  vd->addVirtualDetector( vdIdFEBEdge,
				  parentInMu2e,
				  0,
				  posInnerFEB);
	  ++vdIdFEBEdge;

	  vd->addVirtualDetector( vdIdFEBEdge,
				  parentInMu2e,
				  0,
				  posInnerFEB);
	  ++vdIdFEBEdge;


        }

      }

      if ( c.getBool("hasSTM",false) ) {
        // temporary arrangements until MSTM is in GeometryService
        //const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
        //CLHEP::Hep3Vector parentPositionInMu2e = parent.centerInMu2e();
        GeomHandle<Mu2eEnvelope> env;
        const CLHEP::Hep3Vector parentPositionInMu2e(
                                                     (env->xmax() + env->xmin())/2.,
                                                     (env->ymax() + env->ymin())/2.,
                                                     (env->zmax() + env->zmin())/2.
                                                    );

        GeomHandle<DetectorSolenoid> ds;
        CLHEP::Hep3Vector const & dsP ( ds->position() );

        const CLHEP::Hep3Vector zeroVector(0.,0.,0.);

        GeomHandle<CosmicRayShield> CRS;
        const double z_crv_max = CRS->getSectorPosition("D").z() + (CRS->getSectorHalfLengths("D"))[2];

        //Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
        // Our reference z is the downstream edge of the CRV
        const CLHEP::Hep3Vector mstmReferencePositionInMu2e(dsP.x(),
                                              0.0,
                                              z_crv_max );
        const CLHEP::Hep3Vector mstmReferencePositionInParent = mstmReferencePositionInMu2e - parentPositionInMu2e;

        //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
        //CLHEP::Hep3Vector const& parentInMu2e = parent.centerInMu2e();
        //WARNING: This must be the same as in constructMSTM.cc for now.
        GeomHandle<Mu2eHall> hall;
        const double z_hall_inside_max = hall->getWallExtentz("dsArea",1)/CLHEP::mm;//the integer allows you to specify which side of which wall you want the z for: 1 = west side of east wall (i.e. the z of the inside surface of the east wall)

        const double mstmMotherHalfLength = (z_hall_inside_max - z_crv_max)/2.0;

        if ( c.getBool("vd.STMUpStr.build", false) ) {
          //place this VD 1 cm downstream of the CRS (Cosmic Ray Shield)
          // const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
          // const double yExtentLow      = c.getDouble("yOfFloorSurface.below.mu2eOrigin");
          // const double y_vd_center     = (yExtentLow + y_crv_max)/2.0;
          const double y_vd_center = 0.0;

          CLHEP::Hep3Vector vdPositionWRTparent     = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,y_vd_center, 1.0*mm-vd->_halfLength);

          vd->addVirtualDetector(VirtualDetectorId::STM_UpStr, //ID
                                 parentPositionInMu2e,         //reference position
                                 0x0,                          //rotation
                                 vdPositionWRTparent);         //placement w.r.t. reference


           if ( verbosityLevel > 0) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_UpStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_UpStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_UpStr) <<endl;
           }
        }

//         if ( c.getBool("vd.STMCRVShieldDnStr.build", false) ) {
//           //place this VD just downstream of the shield wall dnStr of the CRV (Cosmic Ray Veto)
//           // const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
//           // const double yExtentLow      = c.getDouble("yOfFloorSurface.below.mu2eOrigin");
//           // const double y_vd_center     = (yExtentLow + y_crv_max)/2.0;
//           const double y_vd_center = 0.0;
//
//           //CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,y_vd_center, -mstmMotherHalfLength+1.0*mm-vd->_halfLength);
//           const double z_offset =   c.getDouble("stm.crvshield.upStrSpace")
//                                   + 2.0*c.getDouble("stm.crvshield.halflength")
//                                   + 1.0*mm   // another 1mm gap
//                                   - vd->_halfLength;
//           CLHEP::Hep3Vector vdPositionWRTparent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,y_vd_center, z_offset);
//
//           vd->addVirtualDetector(VirtualDetectorId::STM_CRVShieldDnStr, //ID
//                                  parentPositionInMu2e,//mstmMotherPositionInMu2e,//reference position
//                                  0x0,                               //rotation
//                                  vdPositionWRTparent);    //vdPositionWRTmstmMother);//placement w.r.t. reference
//
//            if ( verbosityLevel > -1) {
//               cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_CRVShieldDnStr) << endl;
//               cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_CRVShieldDnStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_CRVShieldDnStr) <<endl;
//            }
//         }

        if ( c.getBool("vd.STMFieldOfViewCollDnStr.build", false) ) {
          //place this VD just downstream of the STM Field-Of-View Collimator
          // const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
          // const double yExtentLow      = c.getDouble("yOfFloorSurface.below.mu2eOrigin");
          // const double y_vd_center     = (yExtentLow + y_crv_max)/2.0;
          const double y_vd_center = 0.0;

          double z_offset =   c.getDouble("stm.magnet.UpStrSpace")
                                  + 2.0*c.getDouble("stm.magnet.halfLength")
                                  + c.getDouble("stm.FOVcollimator.UpStrSpace")
                                  + 2.0*c.getDouble("stm.FOVcollimator.halfLength")
                                  + 100.0*mm   // a 10cm gap
                                  - vd->_halfLength;
          if (c.getBool("stm.pipe.build", false)){
             z_offset +=  2.0*c.getDouble("stm.pipe.DnStrHalfLength");
          }
          CLHEP::Hep3Vector vdPositionWRTparent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,y_vd_center, z_offset);

          vd->addVirtualDetector(VirtualDetectorId::STM_FieldOfViewCollDnStr, //ID
                                 parentPositionInMu2e,    //reference position
                                 0x0,                     //rotation
                                 vdPositionWRTparent);    //placement w.r.t. reference


           if ( verbosityLevel > 0) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_FieldOfViewCollDnStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_FieldOfViewCollDnStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_FieldOfViewCollDnStr) <<endl;
           }
        }

        if ( c.getBool("vd.STMMagDnStr.build", false) ) {
          //place this VD downstream of the magnet
          // const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
          // const double yExtentLow      = c.getDouble("yOfFloorSurface.below.mu2eOrigin");
          // const double y_vd_center     = (yExtentLow + y_crv_max)/2.0;
          const double y_vd_center = 0.0;

          double z_offset =   c.getDouble("stm.magnet.UpStrSpace")
                            + 2.0*c.getDouble("stm.magnet.halfLength")
                            + vd->_halfLength;
          if (c.getBool("stm.pipe.build", false)){
             z_offset +=  2.0*c.getDouble("stm.pipe.DnStrHalfLength");
          }
          if (c.getBool("stm.FOVcollimator.build", false)){
             z_offset +=  0.5*c.getDouble("stm.FOVcollimator.UpStrSpace");
          }

          CLHEP::Hep3Vector vdPositionWRTparent     = mstmReferencePositionInParent +  CLHEP::Hep3Vector(0.0,y_vd_center, z_offset );

          vd->addVirtualDetector(VirtualDetectorId::STM_MagDnStr, //ID
                                 parentPositionInMu2e,
                                 0x0,                 //rotation
                                 vdPositionWRTparent);


          if ( verbosityLevel > 0) {
            cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_MagDnStr) << endl;
            cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_MagDnStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_MagDnStr) <<endl;
          }
        }

        if ( c.getBool("vd.STMSSCollUpStr.build", false) ) {
          const double mstmZAllowed =  c.getDouble("stm.z.allowed");
          const double mstmCollHalfLength =  c.getDouble("stm.SScollimator.halfLength");
          CLHEP::Hep3Vector mstmCollPositionInParent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,0.0,2.0*mstmMotherHalfLength) - CLHEP::Hep3Vector(0.0,0.0,mstmZAllowed) + CLHEP::Hep3Vector(0.0,0.0,mstmCollHalfLength);
          CLHEP::Hep3Vector vdPositionWRTparent     = mstmCollPositionInParent + CLHEP::Hep3Vector(0.0,0.0,-mstmCollHalfLength-10.0);

          vd->addVirtualDetector(VirtualDetectorId::STM_SpotSizeCollUpStr, //ID
                                 parentPositionInMu2e, //reference position
                                 0x0,                  //rotation
                                 vdPositionWRTparent); //placement w.r.t. reference


           if ( verbosityLevel > 0) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_SpotSizeCollUpStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_SpotSizeCollUpStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_SpotSizeCollUpStr) <<endl;
           }
        }

        if ( c.getBool("vd.STMCollDnStr.build", false) ) {
          const double mstmZAllowed =  c.getDouble("stm.z.allowed");
          const double mstmCollHalfLength =  c.getDouble("stm.SScollimator.halfLength");
          CLHEP::Hep3Vector mstmCollPositionInParent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,0.0,2.0*mstmMotherHalfLength) - CLHEP::Hep3Vector(0.0,0.0,mstmZAllowed) + CLHEP::Hep3Vector(0.0,0.0,mstmCollHalfLength);
          const double mstmCanUpStrSpace             = c.getDouble("stm.det1.can.UpStrSpace");

          CLHEP::Hep3Vector vdPositionWRTparent     = mstmCollPositionInParent + CLHEP::Hep3Vector(0.0,0.0,mstmCollHalfLength+0.5*mstmCanUpStrSpace+vd->_halfLength);

          vd->addVirtualDetector(VirtualDetectorId::STM_CollDnStr, //ID
                                 parentPositionInMu2e, //reference position
                                 0x0,                  //rotation
                                 vdPositionWRTparent); //placement w.r.t. reference


           if ( verbosityLevel > 0) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_CollDnStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_CollDnStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_CollDnStr) <<endl;
           }
        }

        if ( c.getBool("vd.STMDet1UpStr.build", false) ) {
          const double mstmZAllowed =  c.getDouble("stm.z.allowed");
          const double mstmCollHalfLength =  c.getDouble("stm.SScollimator.halfLength");
          CLHEP::Hep3Vector mstmCollPositionInParent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,0.0,2.0*mstmMotherHalfLength) - CLHEP::Hep3Vector(0.0,0.0,mstmZAllowed) + CLHEP::Hep3Vector(0.0,0.0,mstmCollHalfLength);
          const double mstmCanUpStrSpace            =  c.getDouble("stm.det1.can.UpStrSpace");
          const double mstmCanUpStrWindowHalfLength =  c.getDouble("stm.det1.can.UpStrWindowHalfLength");

          CLHEP::Hep3Vector vdPositionWRTparent     = mstmCollPositionInParent + CLHEP::Hep3Vector(c.getDouble("stm.det1.xoffset"),0.0,mstmCollHalfLength+mstmCanUpStrSpace+2.0*mstmCanUpStrWindowHalfLength+vd->_halfLength);

          vd->addVirtualDetector(VirtualDetectorId::STM_Det1UpStr,   //ID
                                 parentPositionInMu2e, //reference position
                                 0x0,                  //rotation
                                 vdPositionWRTparent); //placement w.r.t. reference


           if ( verbosityLevel > 0) {
             cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_Det1UpStr) << endl;
             cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_Det1UpStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_Det1UpStr) <<endl;
           }
        }

        if ( c.getBool("vd.STMDet2UpStr.build", false) ) {
          const double mstmZAllowed =  c.getDouble("stm.z.allowed");
          const double mstmCollHalfLength =  c.getDouble("stm.SScollimator.halfLength");
          CLHEP::Hep3Vector mstmCollPositionInParent = mstmReferencePositionInParent + CLHEP::Hep3Vector(0.0,0.0,2.0*mstmMotherHalfLength) - CLHEP::Hep3Vector(0.0,0.0,mstmZAllowed) + CLHEP::Hep3Vector(0.0,0.0,mstmCollHalfLength);
          const double mstmCanUpStrSpace            =  c.getDouble("stm.det2.can.UpStrSpace");
          const double mstmCanUpStrWindowHalfLength =  c.getDouble("stm.det2.can.UpStrWindowHalfLength");

          CLHEP::Hep3Vector vdPositionWRTparent     = mstmCollPositionInParent + CLHEP::Hep3Vector(c.getDouble("stm.det2.xoffset"),0.0,mstmCollHalfLength+mstmCanUpStrSpace+2.0*mstmCanUpStrWindowHalfLength+vd->_halfLength);

          vd->addVirtualDetector(VirtualDetectorId::STM_Det2UpStr,   //ID
                                 parentPositionInMu2e, //reference position
                                 0x0,                  //rotation
                                 vdPositionWRTparent); //placement w.r.t. reference


          if ( verbosityLevel > 0) {
            cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::STM_Det2UpStr) << endl;
            cout << "               at local=" << vd->getLocal(VirtualDetectorId::STM_Det2UpStr) << " global="<< vd->getGlobal(VirtualDetectorId::STM_Det2UpStr) <<endl;
          }
        }
      }

      if ( c.getBool("pbar.coll1In.build", false) ) {

         GeomHandle<Beamline> bl;
         TransportSolenoid const& ts = bl->getTS();
         CollimatorTS1 const& coll1  = ts.getColl1() ;

         double pbarTS1InHalfLength = c.getDouble("pbar.coll1In.halfLength");
         double pbarTS1InOffset = c.getDouble("pbar.coll1In.offset", 1.0);

         CLHEP::Hep3Vector pbarTS1InPos = coll1.getLocal();
	 if (verbosityLevel > 0){
	 std::cout << "starting coll1 position " << pbarTS1InPos << std::endl;
	 }
         CLHEP::Hep3Vector parentCenterInMu2e;
	 // 
	 // make the VD 1 mm upstream of the window; the window is much thinner.  Just add a throw to make sure...
	 double windowToVDOffset = 1.0 * CLHEP::mm;  //envisaging a day when this will be configurable
	 if (windowToVDOffset < pbarTS1InHalfLength) throw cet::exception("GEOM") << __func__ << "window thicker than pbarTS1InHalfLength" << std::endl;
	 CLHEP::Hep3Vector windowLocIn(0.,0.,0.);
 
	 if (pbarTS1InOffset >= 0.0) {
	   // use local when put in the TS1Vacuum
	   pbarTS1InPos = coll1.getLocal();
	   //
	   // put these together before you change pbarTS1InPos.z()
	   windowLocIn = pbarTS1InPos;
	   windowLocIn.setZ( pbarTS1InPos.z()  + 2.*vdHL - windowToVDOffset + pbarTS1InOffset);

	   pbarTS1InPos.setZ( pbarTS1InPos.z() - coll1.halfLength() + 2.*vdHL + pbarTS1InHalfLength + pbarTS1InOffset);
	   parentCenterInMu2e = ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getGlobal();
         }
         else { // pbarTS1InOffset < 0.0
	   // use global when put in the HallAir
	   Tube const & psVacuumParams  = GeomHandle<PSVacuum>()->vacuum();
	   pbarTS1InPos = ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getGlobal();
	   //
	   // put these together before you change pbarTS1InPos.z()
	   windowLocIn = pbarTS1InPos;
	   windowLocIn.setZ   ( pbarTS1InPos.z() - ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() - windowToVDOffset + pbarTS1InOffset );

	   pbarTS1InPos.setZ( pbarTS1InPos.z() - ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() - pbarTS1InHalfLength + pbarTS1InOffset);
	   if (verbosityLevel > 0){
	     std::cout << pbarTS1InPos.z() << " " << ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() << " " << pbarTS1InHalfLength << " " <<  pbarTS1InOffset << std::endl;
	   }
	   CLHEP::Hep3Vector psVacuumOriginInMu2e = psVacuumParams.originInMu2e();
	   pbarTS1InPos = pbarTS1InPos - psVacuumOriginInMu2e;
	   parentCenterInMu2e = psVacuumOriginInMu2e;
	   windowLocIn = windowLocIn - psVacuumOriginInMu2e;
         }
         CLHEP::Hep3Vector posPSPbarIn = pbarTS1InPos;
	 posPSPbarIn.setZ( pbarTS1InPos.z() - pbarTS1InHalfLength - vdHL );
	 if (verbosityLevel > 0){
	   cout << "posPSPbarIn, windowLocIn, and parent Center is psVacuumOrigin " << posPSPbarIn << " " << windowLocIn << " " << parentCenterInMu2e << endl;
	 }
         vd->addVirtualDetector(VirtualDetectorId::PSPbarIn, parentCenterInMu2e, 0, windowLocIn);


	 //
	 //floating VD
	 //      CLHEP::Hep3Vector posPSPbarOut = pbarTS1InPos;
	 //	 posPSPbarOut.setZ( pbarTS1InPos.z() + pbarTS1InHalfLength + vdHL );
         //      posPSPbarOut.setZ( pbarTS1InPos.z() + windowToVDOffset + vdHL );
	 //         vd->addVirtualDetector(VirtualDetectorId::PSPbarOut, parentCenterInMu2e, 0, posPSPbarOut);

	 CLHEP::Hep3Vector windowLocOut = windowLocIn;
	 windowLocOut.setZ(windowLocOut.z() + 2.*windowToVDOffset);
	 if (verbosityLevel > 0){
	   std::cout << "windowLocOut = " << windowLocOut << std::endl;
	 }
         vd->addVirtualDetector(VirtualDetectorId::PSPbarOut, parentCenterInMu2e, 0, windowLocOut);


	 if ( verbosityLevel > 0 ) {
	   cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::PSPbarIn) << endl;
	   cout << "               at local=" << vd->getLocal(VirtualDetectorId::PSPbarIn) << " global="<< vd->getGlobal(VirtualDetectorId::PSPbarIn) <<endl;
	   cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::PSPbarOut) << endl;
	   cout << "               at local=" << vd->getLocal(VirtualDetectorId::PSPbarOut) << " global="<< vd->getGlobal(VirtualDetectorId::PSPbarOut) <<endl;
	 }
      }

      if(c.getBool("vd.crv.build", false))
      {
        GeomHandle<CosmicRayShield> CRS;
        GeomHandle<Mu2eEnvelope> env;
        const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                       (env->xmax() + env->xmin())/2.,
                                                       (env->ymax() + env->ymin())/2.,
                                                       (env->zmax() + env->zmin())/2.
                                                       );

        for(int vdId=VirtualDetectorId::CRV_R; vdId<=VirtualDetectorId::CRV_U; vdId++)
        {
          std::string vdName;
          vdName = VirtualDetector::volumeName(vdId).back();
          const std::vector<double> crvHalfLengths = CRS->getSectorHalfLengths(vdName);
          const CLHEP::Hep3Vector vdDirection = c.getHep3Vector("crs.vdDirection"+vdName);

          CLHEP::Hep3Vector vdPosInHall = CRS->getSectorPosition(vdName) - hallFormalCenterInMu2e;
          for(int i=0; i<3; i++) vdPosInHall[i] += vdDirection[i]*crvHalfLengths[i] + vdDirection[i]*vd->_halfLength;

          vd->addVirtualDetector(vdId, //ID
                                 hallFormalCenterInMu2e,  //reference position
                                 0x0,                     //rotation
                                 vdPosInHall);            //placement w.r.t. reference



          if(verbosityLevel > 0)
          {
            cout << " Constructing " << VirtualDetector::volumeName(vdId) << endl;
            cout << "               at local=" << vd->getLocal(vdId) << " global="<< vd->getGlobal(vdId) <<endl;
          }
        }
      }

    } // if(hasVirtualDetector)

    return vd;

  } // make()

} // namespace mu2e
