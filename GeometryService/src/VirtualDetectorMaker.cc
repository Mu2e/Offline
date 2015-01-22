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
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxes.hh"

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

      // Need some data from other subsystems
      GeomHandle<Beamline> bg;
      double solenoidOffset = bg->solenoidOffset();

      // VD Coll1_In and Coll1_Out are at the front and back of
      // collimator 1, which is placed inside TS1.

      double coll1HL = bg->getTS().getColl1().halfLength();

      const HepRotation *ts1rot = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS1,
                                                        TransportSolenoid::TSRadialPart::IN)->getRotation();
      Hep3Vector   ts1pos = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS1,
                                                  TransportSolenoid::TSRadialPart::IN)->getGlobal();

      Hep3Vector coll1pos = bg->getTS().getColl1().getLocal();


      Hep3Vector deltaZ1(0,0,coll1HL-vdHL);

      vd->addVirtualDetector( VirtualDetectorId::Coll1_In,
                               ts1pos, ts1rot, coll1pos-deltaZ1);
      vd->addVirtualDetector( VirtualDetectorId::Coll1_Out,
                               ts1pos, ts1rot, coll1pos+deltaZ1);

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

      double coll5HL = bg->getTS().getColl5().halfLength();

      const HepRotation *ts5rot = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS5,
                                                        TransportSolenoid::TSRadialPart::IN)->getRotation();
      Hep3Vector   ts5pos = bg->getTS().getTSCryo(TransportSolenoid::TSRegion::TS5,
                                                        TransportSolenoid::TSRadialPart::IN)->getGlobal();
      Hep3Vector coll5pos = bg->getTS().getColl5().getLocal();


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

      std::cout << "coll 5 " << coll5pos.z() << " " << deltaZ5.z() << " " << targetOffset.z() << " " << shift.z() << std::endl;
      const Hep3Vector STMOffset(targetOffset.x()-shift.x(),targetOffset.y()-shift.y(), targetOffset.z()-shift.z() - 0.5*( (coll5pos.z()+deltaZ5.z()) - (targetOffset.z()-shift.z()) ));
      vd->addVirtualDetector( VirtualDetectorId::STMUpstream,
			      ds2centerInMu2e,0,STMOffset);



      if (c.getBool("hasTTracker",false)){

        ostringstream vdName(VirtualDetectorId::name(VirtualDetectorId::TT_Mid));

        if(c.getInt("ttracker.numDevices")%2!=0){
          throw cet::exception("GEOM")
            << "This virtual detector " << vdName
            << " can only be placed if the TTracker has an even number of devices \n";
        }

        TTracker const & ttracker = *(GeomHandle<TTracker>());
        Hep3Vector ttOffset(-solenoidOffset,0.,ttracker.z0());

        // VD TT_Mid is placed inside the ttracker mother volume in the
        // middle of the ttracker shifted by the half length of vd
        // VD TT_MidInner is placed inside the ttracker at the same z position as
        // VD TT_Mid but from radius 0 to the inner radius of the ttracker
        // mother volume. However, its mother volume is DS3Vacuum
        // which has a different offset. We will use the global offset
        // here (!) as DS is not in the geometry service yet

        Hep3Vector vdTTMidOffset(0.,0.,0.);

        vd->addVirtualDetector( VirtualDetectorId::TT_Mid,
                                 ttOffset, 0, vdTTMidOffset);

        vd->addVirtualDetector( VirtualDetectorId::TT_MidInner,
                                 ttOffset, 0, vdTTMidOffset);

        //       int static const verbosityLevel = 1;
        //       if (verbosityLevel >0) {
        //         for ( int vdId=11; vdId<=12; ++vdId) {
        //           cout << __func__ << " VD " << vdId << " offsets L, G " <<
        //             vd->getLocal(vdId) << ", " <<
        //             vd->getGlobal(vdId) << endl;
        //         }
        //       }

        // Global position is in Mu2e coordinates; local position in the detector system.
        double zFrontGlobal = ttracker.mother().position().z()-ttracker.mother().tubsParams().zHalfLength()-vdHL;
        double zBackGlobal  = ttracker.mother().position().z()+ttracker.mother().tubsParams().zHalfLength()+vdHL;
        double zFrontLocal  = zFrontGlobal - ttracker.z0();
        double zBackLocal   = zBackGlobal  - ttracker.z0();

        Hep3Vector vdTTFrontOffset(0.,
                                   0.,
                                   zFrontLocal);

        // VD TT_FrontHollow is placed outside the ttracker mother
        // volume in front of the ttracker "outside" of the proton
        // absorber


        // formally VD TT_FrontHollow, TT_FrontPA are placed in DS3Vacuum, but it is a
        // complicated subtraction volume, so we pretend to place them in
        // the TTracker and rely on the global offsets in the mu2e
        // detector frame (note that their local offsets are wrt TTracker)

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
        // placed on the inner and outer surface of the ttracker envelope

        Hep3Vector vdTTOutSurfOffset(0.,0.,ttracker.mother().position().z()-ttracker.z0());

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
        GeomHandle<ExtNeutShieldCendBoxes> enscendb;

        const std::vector<CLHEP::Hep3Vector>& ENSCBcentersOfBoxes = enscendb->centersOfBoxes();

        size_t nBox = ENSCBcentersOfBoxes.size();
        size_t ib;
        for(ib = 0; ib < nBox; ++ib) {
          if ( enscendb->hasHole(ib) ) break;
        }
        int hID = enscendb->holeIndex(ib);
        // locations are wrt HallAir
        // for some reason the location has to be taken from the box and not the hole tbd
        //        CLHEP::Hep3Vector holeLocation = enscendb->holeLocation(hID);
        CLHEP::Hep3Vector holeLocation = ENSCBcentersOfBoxes[ib];

        GeomHandle<DetectorSolenoid> ds;
        CLHEP::Hep3Vector const & dsP ( ds->position() );
        CLHEP::Hep3Vector vdPositionInMu2e(dsP.x(), dsP.y(),
                                           holeLocation.z() + enscendb->holeHalfLength(hID) + vd->getHalfLength());

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
        
        int static const verbosityLevel = 1;
        if ( verbosityLevel > 0) {
           cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::DSNeutronShieldExit) << endl;
           cout << "               at local=" << vd->getLocal(VirtualDetectorId::DSNeutronShieldExit) << " global="<< vd->getGlobal(VirtualDetectorId::DSNeutronShieldExit) <<endl;
        }

      }
                              

      //placing virtual detector around the calorimeter vanes
      if (c.getBool("hasVaneCalorimeter",true)){
	GeomHandle<VaneCalorimeter> cg;
	
	int vdIdFront = VirtualDetectorId::EMC_0_FrontIn;
	int vdIdEdge = VirtualDetectorId::EMC_0_EdgeIn;
	int vdIdSurf = VirtualDetectorId::EMC_0_SurfIn;	      

	const Hep3Vector FrontOffsetOut(0.0, 0.0, (cg->vane(0).size().z() + vdHL) );
	const Hep3Vector FrontOffsetIn(0.0, 0.0,  -(cg->vane(0).size().z() + vdHL) );

	const Hep3Vector EdgeOffsetOut(0.0,  (cg->vane(0).size().y() + vdHL), 0.0);
	const Hep3Vector EdgeOffsetIn(0.0,  -(cg->vane(0).size().y() + vdHL), 0.0);
	
	const Hep3Vector ROplaneOffsetOut( -(cg->vane(0).size().x() + vdHL), 0.0, 0.0);
	const Hep3Vector ROplaneOffsetIn( (cg->vane(0).size().x() + vdHL), 0.0, 0.0);

	for(int i=0; i<cg->nVane(); ++i){
	    Hep3Vector FrontOffsetInRot(0.0, 0.0, 0.0);
	    FrontOffsetInRot = (cg->vane(i).rotation().inverse())*FrontOffsetIn;
	   //  cout<<"vane origin ("<<i<<") = "<<cg->vane(i).origin() <<endl;
// 	    cout<<"vdIdFront = "<<vdIdFront<<endl;	  
// 	    cout<<" FrontOffsetInRot = "<<FrontOffsetInRot<<endl;
	  
	    vd->addVirtualDetector( vdIdFront,
				    cg->vane(i).origin(),
				    0,
				    FrontOffsetInRot);
	    ++vdIdFront;

	    Hep3Vector FrontOffsetOutRot(0.0, 0.0, 0.0);
	    FrontOffsetOutRot = (cg->vane(i).rotation().inverse())*FrontOffsetOut;
	  
	    // cout<<"vdIdFront = "<<vdIdFront<<endl;
	    // cout<<" FrontOffsetOutRot = "<<FrontOffsetOutRot<<endl;
	    

	    vd->addVirtualDetector( vdIdFront,
				    cg->vane(i).origin(),
				    0,
				    FrontOffsetOutRot);
	    ++vdIdFront;
	
	    Hep3Vector EdgeOffsetInRot(0.0, 0.0, 0.0);
	    EdgeOffsetInRot = (cg->vane(i).rotation().inverse())*EdgeOffsetIn;
// 	    cout<<"vane origin ("<<i<<") = "<<cg->vane(i).origin() <<endl;
// 	    cout<<"vdIdEdge = "<<vdIdEdge<<endl;	  
// 	    cout<<" EdgeOffsetInRot = "<<EdgeOffsetInRot<<endl;
	  
	    vd->addVirtualDetector( vdIdEdge,
				    cg->vane(i).origin(),
				    0,
				    EdgeOffsetInRot);
	    ++vdIdEdge;

	    Hep3Vector EdgeOffsetOutRot(0.0, 0.0, 0.0);
	    EdgeOffsetOutRot = (cg->vane(i).rotation().inverse())*EdgeOffsetOut;
	  
	 //    cout<<"vdIdEdge = "<<vdIdEdge<<endl;
// 	    cout<<" EdgeOffsetOutRot = "<<EdgeOffsetOutRot<<endl;
	    

	    vd->addVirtualDetector( vdIdEdge,
				    cg->vane(i).origin(),
				    0,
				    EdgeOffsetOutRot);
	    ++vdIdEdge;
	
	   

	    Hep3Vector ROplaneOffsetInRot(0.0, 0.0, 0.0);
	    ROplaneOffsetInRot = (cg->vane(i).rotation().inverse())*ROplaneOffsetIn;
	  //   cout<<"//---------------------------------------//"<<endl
// 		<<"//---------------------------------------//"<<endl
// 		<<"//---------------------------------------//"<<endl;
// 	    cout<<"vdIdSurf = "<<vdIdSurf<<endl;
// 	    cout<<" ROplaneOffsetInRot = "<<ROplaneOffsetInRot<<endl;
	    
	    vd->addVirtualDetector( vdIdSurf,
				    cg->vane(i).origin(),
				    0,
				    ROplaneOffsetInRot);
	    ++vdIdSurf;
	    
	    Hep3Vector ROplaneOffsetOutRot(0.0, 0.0, 0.0);
	    ROplaneOffsetOutRot = (cg->vane(i).rotation().inverse())*ROplaneOffsetOut;
	  //   cout<<"vdIdSurf = "<<vdIdSurf<<endl;	   
// 	    cout<<" ROplaneOffsetOutRot = "<<ROplaneOffsetOutRot<<endl;

	    vd->addVirtualDetector( vdIdSurf,
				    cg->vane(i).origin(),
				    0,
				    ROplaneOffsetOutRot);
	    ++vdIdSurf;
	  }
	    
      }
      
      if (c.getBool("hasDiskCalorimeter",true)){
	GeomHandle<DiskCalorimeter> cg;
	
	int vdIdSurf = VirtualDetectorId::EMC_Disk_0_SurfIn;
	int vdIdEdge = VirtualDetectorId::EMC_Disk_0_EdgeIn;

	Hep3Vector EdgeOffset(0.0, 0.0, 0.0);
	const Hep3Vector OffsetOut(0.0, 0.0, (cg->disk(0).size().z() + vdHL) );
	const Hep3Vector OffsetIn(0.0, 0.0,  -(cg->disk(0).size().z() + vdHL) );
	
	for(size_t i=0; i<cg->nDisk(); ++i){
	 
	 
	
// 	  cout<<"disk origin ("<<i<<") = "<<cg->disk(i).origin() <<endl;
// 	  cout<<"vdIdSurf = "<<vdIdSurf<<endl;	  
	  
	  vd->addVirtualDetector( vdIdSurf,
				  cg->disk(i).origin(),
				  0,
				  OffsetIn);
	  ++vdIdSurf;

	  // cout<<"vdIdSurf = "<<vdIdSurf<<endl;	  
	 
	  vd->addVirtualDetector( vdIdSurf,
				  cg->disk(i).origin(),
				  0,
				  OffsetOut);
	  ++vdIdSurf;
	
	//   cout<<"disk origin ("<<i<<") = "<<cg->disk(i).origin() <<endl;
// 	  cout<<"vdIdEdge = "<<vdIdEdge<<endl;	  
	  
	  vd->addVirtualDetector( vdIdEdge,
				  cg->disk(i).origin(),
				  0,
				  EdgeOffset);
	  ++vdIdEdge;

	  //	  cout<<"vdIdEdge = "<<vdIdEdge<<endl;

	  vd->addVirtualDetector( vdIdEdge,
				  cg->disk(i).origin(),
				  0,
				  EdgeOffset);
	  ++vdIdEdge;
	}
	    
      }

      if ( c.getBool("mstm.build", false) ) {
        // VirtualDetector in front of the MSTM detector crystal
        // temporary arrangements till MSTM is in GeometryService
        
	const double mstmUpStreamWallUpStrSpace = c.getDouble("mstm.wallUpStr.UpStrSpace");
        const double mstmUpStreamWallHalfLength = c.getDouble("mstm.wallUpStr.halfLength");
	const double mstmMagnetUpStrSpace       = c.getDouble("mstm.magnet.UpStrSpace");
        const double mstmMagnetHalfLength       = c.getDouble("mstm.magnet.halfLength");
        const double mstmColl1UpStrSpace        = c.getDouble("mstm.collimator1.UpStrSpace");
        const double mstmColl1HalfLength        = c.getDouble("mstm.collimator1.halfLength");
        //const double mstmPipe0HalfLength       = c.getDouble("mstm.pipe0.halfLength");
        //const double mstmColl1UpStrSpace       = c.getDouble("mstm.collimator1.UpStrSpace");
        //const double mstmColl1HalfLength       = c.getDouble("mstm.collimator1.halfLength");
        const double mstmShutterUpStrSpace     = c.getDouble("mstm.shutter.UpStrSpace");
        const double mstmShutterNumberSegments = c.getInt("mstm.shutter.numberSegments");
	double mstmShutterHalfLength           = 0.;
	for (int segment = 1; segment <= mstmShutterNumberSegments; ++segment) {
	  std::stringstream mstmShutterSegmentConfig;
	  mstmShutterSegmentConfig << "mstm.shutter.segment" << segment << ".halfLength";
	  mstmShutterHalfLength += c.getDouble(mstmShutterSegmentConfig.str());
	}
	const double mstmColl2UpStrSpace       = c.getDouble("mstm.collimator2.UpStrSpace");
        const double mstmColl2HalfLength       = c.getDouble("mstm.collimator2.halfLength");
        const double mstmPipe1UpStrSpace       = c.getDouble("mstm.pipe1.UpStrSpace");
	const double mstmPipe1HalfLength       = c.getDouble("mstm.pipe1.halfLength");
        const double mstmColl3UpStrSpace       = c.getDouble("mstm.collimator3.UpStrSpace");
        const double mstmColl3HalfLength       = c.getDouble("mstm.collimator3.halfLength");
        const double mstmCanUpStrSpace         = c.getDouble("mstm.can.UpStrSpace");        
        //const double mstmCanHalfLength         = c.getDouble("mstm.can.halfLength");
        //const double mstmCrystalHalfLength     = c.getDouble("mstm.crystal.halfLength");

        GeomHandle<ExtNeutShieldCendBoxes> enscendb;

        const std::vector<CLHEP::Hep3Vector>& ENSCBcentersOfBoxes = enscendb->centersOfBoxes();

        size_t nBox = ENSCBcentersOfBoxes.size();
        size_t ib;
        for(ib = 0; ib < nBox; ++ib) {
          // cout << __func__ << " " << ib << " ENSCBcentersOfBoxes: " << ENSCBcentersOfBoxes[ib] << endl;
          if ( enscendb->hasHole(ib) ) break;
        }
        int hID = enscendb->holeIndex(ib);
        // locations are wrt HallAir
        // for some reason the location has to be taken from the box and not the hole tbd
        //        CLHEP::Hep3Vector holeLocation = enscendb->holeLocation(hID);
        CLHEP::Hep3Vector holeLocation = ENSCBcentersOfBoxes[ib];

        GeomHandle<DetectorSolenoid> ds;
        CLHEP::Hep3Vector const & dsP ( ds->position() );
        CLHEP::Hep3Vector mstmReferencePositionInMu2e(dsP.x(), dsP.y(),
                                                      holeLocation.z() + enscendb->holeHalfLength(hID) + 2.0*vd->getHalfLength() );

        GeomHandle<Mu2eEnvelope> env;
        const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                       (env->xmax() + env->xmin())/2.,
                                                       (env->ymax() + env->ymin())/2.,
                                                       (env->zmax() + env->zmin())/2.
                                                       );        
        //mstmReferencePositionInMu2e  = mstmReferencePositionInMu2e - hallFormalCenterInMu2e;
      
        //WARNING: This must be changed so that we can use the VolumeInfo to get the 
        //         mother position
        //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
        //CLHEP::Hep3Vector const& parentInMu2e = parent.centerInMu2e();
        //WARNING: This must be the same as in constructMSTM.cc for now.
        const double mstmMotherHalfLength =   1.0 * c.getDouble("mstm.pipe0.halfLength")
                                    + 0.5 * c.getDouble("mstm.collimator1.UpStrSpace")
                                    + 1.0 * c.getDouble("mstm.collimator1.halfLength")
                                    + 0.5 * c.getDouble("mstm.collimator2.UpStrSpace")
                                    + 1.0 * c.getDouble("mstm.collimator2.halfLength")
                                    + 0.5 * c.getDouble("mstm.shutter.UpStrSpace")
                                    + 1.0 * 500.0 //over-estimate 1 meter for the shutter length
                                    + 1.0 * c.getDouble("mstm.pipe1.halfLength")
                                    + 0.5 * c.getDouble("mstm.collimator3.UpStrSpace")
                                    + 1.0 * c.getDouble("mstm.collimator3.halfLength")
                                    + 0.5 * c.getDouble("mstm.can.UpStrSpace")
                                    + 1.0 * c.getDouble("mstm.can.halfLength")
                                    + 0.5 * 1500.0; //make the Mother another 0.5m longer
        
        const CLHEP::Hep3Vector mstmMotherPositionInMu2e = mstmReferencePositionInMu2e + CLHEP::Hep3Vector(0.0, 0.0, mstmMotherHalfLength);      

        
        if ( c.getBool("vd.MSTMWallUpStr.build", false) ) {
           
	   //place this VD 1 cm in front of the concrete shielding wall of the STM area
           double z_offset =  1.0*mstmUpStreamWallUpStrSpace
                            - 1.0*cm
                            - vd->_halfLength;
                 
           CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,0.0,z_offset - mstmMotherHalfLength);        

           vd->addVirtualDetector(VirtualDetectorId::MSTM_WallUpStr, //ID
                                  mstmMotherPositionInMu2e,          //reference position
                                  0x0,                               //rotation
                                  vdPositionWRTmstmMother);          //placement w.r.t. reference
           
           
           int static const verbosityLevel = 1;
           if ( verbosityLevel > -1) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::MSTM_WallUpStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::MSTM_WallUpStr) << " global="<< vd->getGlobal(VirtualDetectorId::MSTM_WallUpStr) <<endl;
           }
        }
        
        if ( c.getBool("vd.MSTMColl1DnStr.build", false) ) {
           
           double z_offset =  mstmUpStreamWallUpStrSpace
                            + 2.0*mstmUpStreamWallHalfLength
                            + mstmMagnetUpStrSpace
                            + 2.0*mstmMagnetHalfLength
                            + mstmColl1UpStrSpace
                            + 2.0*mstmColl1HalfLength
                            + 0.5*mstmShutterUpStrSpace
                            - vd->_halfLength;
                 
           CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,0.0,z_offset - mstmMotherHalfLength);        

           vd->addVirtualDetector(VirtualDetectorId::MSTM_Coll1DnStr, //ID
                                  mstmMotherPositionInMu2e,           //reference position
                                  0x0,                                //rotation
                                  vdPositionWRTmstmMother);           //placement w.r.t. reference
           
           
           int static const verbosityLevel = 1;
           if ( verbosityLevel > -1) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::MSTM_Coll1DnStr) << endl;
           }
        }
        
        if ( c.getBool("vd.MSTMShutterDnStr.build", false) ) {
           
           double z_offset =  mstmUpStreamWallUpStrSpace
                            + 2.0*mstmUpStreamWallHalfLength
                            + mstmMagnetUpStrSpace
                            + 2.0*mstmMagnetHalfLength
                            + mstmColl1UpStrSpace
                            + 2.0*mstmColl1HalfLength
                            + mstmShutterUpStrSpace
                            + 2.0*mstmShutterHalfLength
                            + 0.5*mstmColl2UpStrSpace
                            - vd->_halfLength;
                 
           CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,0.0,z_offset - mstmMotherHalfLength);        

           vd->addVirtualDetector(VirtualDetectorId::MSTM_ShutterDnStr, //ID
                                  mstmMotherPositionInMu2e,           //reference position
                                  0x0,                                //rotation
                                  vdPositionWRTmstmMother);           //placement w.r.t. reference
           
           
           int static const verbosityLevel = 1;
           if ( verbosityLevel > -1) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::MSTM_ShutterDnStr) << endl;
           }
        }
        
        if ( c.getBool("vd.MSTMColl2DnStr.build", false) ) {
           
           double z_offset =  mstmUpStreamWallUpStrSpace
                 + 2.0*mstmUpStreamWallHalfLength
                 + mstmMagnetUpStrSpace
                 + 2.0*mstmMagnetHalfLength
                 + mstmColl1UpStrSpace
                 + 2.0*mstmColl1HalfLength
                 + mstmShutterUpStrSpace
                 + 2.0*mstmShutterHalfLength
                 + mstmColl2UpStrSpace
                 + 2.0*mstmColl2HalfLength
                 + 0.5*mstmPipe1UpStrSpace
                 - vd->_halfLength;
                 
           CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,0.0,z_offset - mstmMotherHalfLength);        
           
           vd->addVirtualDetector(VirtualDetectorId::MSTM_Coll2DnStr, //ID
                                  mstmMotherPositionInMu2e,           //reference position
                                  0x0,                                //rotation
                                  vdPositionWRTmstmMother);           //placement w.r.t. reference
           
           
           int static const verbosityLevel = 1;
           if ( verbosityLevel > -1) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::MSTM_Coll2DnStr) << endl;
           }
        }
        
        if ( c.getBool("vd.MSTMColl3DnStr.build", false) ) {
           
           double z_offset =  mstmUpStreamWallUpStrSpace
                 + 2.0*mstmUpStreamWallHalfLength
                 + mstmMagnetUpStrSpace
                 + 2.0*mstmMagnetHalfLength
                 + mstmColl1UpStrSpace
                 + 2.0*mstmColl1HalfLength
                 + mstmShutterUpStrSpace
                 + 2.0*mstmShutterHalfLength
                 + mstmColl2UpStrSpace
                 + 2.0*mstmColl2HalfLength
                 + mstmPipe1UpStrSpace
                 + 2.0*mstmPipe1HalfLength
                 + mstmColl3UpStrSpace
                 + 2.0*mstmColl3HalfLength
                 + 0.5*mstmCanUpStrSpace
                 //+ mstmCanHalfLength
                 //- mstmCrystalHalfLength 
                 - vd->_halfLength;
                 
           CLHEP::Hep3Vector vdPositionWRTmstmMother = CLHEP::Hep3Vector(0.0,0.0,z_offset - mstmMotherHalfLength);        

           vd->addVirtualDetector(VirtualDetectorId::MSTM_Coll3DnStr, //ID
                                  mstmMotherPositionInMu2e,           //reference position
                                  0x0,                                //rotation
                                  vdPositionWRTmstmMother);           //placement w.r.t. reference
           
           
           int static const verbosityLevel = 1;
           if ( verbosityLevel > -1) {
              cout << " Constructing " << VirtualDetector::volumeName(VirtualDetectorId::MSTM_Coll3DnStr) << endl;
              cout << "               at local=" << vd->getLocal(VirtualDetectorId::MSTM_Coll3DnStr) << " global="<< vd->getGlobal(VirtualDetectorId::MSTM_Coll3DnStr) <<endl;
           }
        }
        
        
      }

    } // if(hasVirtualDetector)

    return vd;

  } // make()

} // namespace mu2e
