//
// Construct and return an VirtualDetector.
//
// $Id: VirtualDetectorMaker.cc,v 1.14 2012/01/25 02:02:00 youzy Exp $
// $Author: youzy $
//

#include <iostream>
#include <iomanip>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "VirtualDetectorGeom/inc/VirtualDetectorMaker.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"
#include "TargetGeom/inc/Target.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

using namespace std;
using namespace CLHEP;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  VirtualDetectorMaker::VirtualDetectorMaker(SimpleConfig const& c)
  {
    art::ServiceHandle<GeometryService> geom;

    _vd = auto_ptr<VirtualDetector>(new VirtualDetector());

    if( ! c.getBool("hasVirtualDetector",false) ) return;

    const double vdHL = c.getDouble("vd.halfLength",0.01*mm);
    _vd->_halfLength = vdHL;

    // Need some data from other subsystems
    GeomHandle<Beamline> bg;
    double solenoidOffset = bg->solenoidOffset();

    // VD Coll1_In and Coll1_Out are at the front and back of
    // collimator 1, which is placed inside TS1.

    //double ts1HL   = bg->getTS().getTS1().getHalfLength();
    double coll1HL = bg->getTS().getColl1().getHalfLength();

    HepRotation *ts1rot = bg->getTS().getTS1().getRotation();
    Hep3Vector   ts1pos = bg->getTS().getTS1().getGlobal();
    Hep3Vector coll1pos = bg->getTS().getColl1().getLocal();
    Hep3Vector deltaZ1(0,0,coll1HL-vdHL);

    _vd->addVirtualDetector( VirtualDetectorId::Coll1_In, 
			     ts1pos, ts1rot, coll1pos-deltaZ1);
    _vd->addVirtualDetector( VirtualDetectorId::Coll1_Out,
			     ts1pos, ts1rot, coll1pos+deltaZ1);

    // VD Coll31_In, Coll31_Out, Coll32_In, Coll32_Out are placed
    // around two sections of collimator 3

    //double ts3HL   = bg->getTS().getTS1().getHalfLength();
    double coll31HL = bg->getTS().getColl31().getHalfLength();
    double coll32HL = bg->getTS().getColl32().getHalfLength();

    HepRotation *ts3rot = bg->getTS().getTS3().getRotation();
    Hep3Vector   ts3pos = bg->getTS().getTS3().getGlobal();

    Hep3Vector coll31pos = bg->getTS().getColl31().getLocal();
    Hep3Vector coll32pos = bg->getTS().getColl32().getLocal();
    Hep3Vector deltaZ31(0,0,coll31HL-vdHL);
    Hep3Vector deltaZ32(0,0,coll32HL-vdHL);

    _vd->addVirtualDetector( VirtualDetectorId::Coll31_In,
			     ts3pos, ts3rot, coll31pos-deltaZ31);
    _vd->addVirtualDetector( VirtualDetectorId::Coll31_Out,
			     ts3pos, ts3rot, coll31pos+deltaZ31);
    _vd->addVirtualDetector( VirtualDetectorId::Coll32_In,
			     ts3pos, ts3rot, coll32pos-deltaZ32);
    _vd->addVirtualDetector( VirtualDetectorId::Coll32_Out,
			     ts3pos, ts3rot, coll32pos+deltaZ32);

    // VD Coll5_In, Coll5_Out are at the front and back of collimator
    // 5, which is placed inside TS5.

    //double ts5HL   = bg->getTS().getTS5().getHalfLength();
    double coll5HL = bg->getTS().getColl5().getHalfLength();

    HepRotation *ts5rot = bg->getTS().getTS5().getRotation();
    Hep3Vector   ts5pos = bg->getTS().getTS5().getGlobal();
    Hep3Vector coll5pos = bg->getTS().getColl5().getLocal();
    Hep3Vector deltaZ5(0,0,coll5HL-vdHL);

    _vd->addVirtualDetector( VirtualDetectorId::Coll5_In,
			     ts5pos, ts5rot, coll5pos-deltaZ5);
    _vd->addVirtualDetector( VirtualDetectorId::Coll5_Out,
			     ts5pos, ts5rot, coll5pos+deltaZ5);

    // VD ST_In, ST_Out are placed inside DS2, just before and after
    // stopping target

    GeomHandle<Target> target;

    double z0    = target->cylinderCenter();
    double zHalf = target->cylinderLength()/2.0;

    double rTorus        = bg->getTS().torusRadius();
    double ts5HalfLength = bg->getTS().getTS5().getHalfLength();
    double ds2HalfLength = c.getDouble("toyDS2.halfLength");
    double ds2Z0         = rTorus + 2.*ts5HalfLength + ds2HalfLength;

    Hep3Vector ds2Offset(-solenoidOffset,0.,ds2Z0);
    Hep3Vector targetOffset(0.,0.,(12000.+z0-ds2Z0));
    Hep3Vector shift(0.,0.,zHalf+vdHL);

    _vd->addVirtualDetector( VirtualDetectorId::ST_In,
			     ds2Offset, 0, targetOffset-shift);
    _vd->addVirtualDetector( VirtualDetectorId::ST_Out,
			     ds2Offset, 0, targetOffset+shift);

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
      // mother volume. However, its mother volume is ToyDS3Vacuum
      // which has a different offset. We will use the global offset
      // here (!) as DS is not in the geometry service yet

      Hep3Vector vdTTMidOffset(0.,0.,0.);

      _vd->addVirtualDetector( VirtualDetectorId::TT_Mid,
			       ttOffset, 0, vdTTMidOffset);

      _vd->addVirtualDetector( VirtualDetectorId::TT_MidInner,
			       ttOffset, 0, vdTTMidOffset);

//       int static const verbosityLevel = 1;
//       if (verbosityLevel >0) {
//         for ( int vdId=11; vdId<=12; ++vdId) {
//           cout << __func__ << " VD " << vdId << " offsets L, G " << 
//             _vd->getLocal(vdId) << ", " << 
//             _vd->getGlobal(vdId) << endl;
//         }
//       }

      Hep3Vector vdTTFrontOffset(0.,
                                 0.,
                                 -ttracker.getTrackerEnvelopeParams().zHalfLength()-vdHL);

      // VD TT_FrontHollow is placed outside the ttracker mother
      // volume in front of the ttracker "outside" of the proton
      // absorber
      

      // formally VD TT_FrontHollow, TT_FrontPA are placed in ToyDS3Vacuum, but it is a
      // complicated subtraction volume, so we pretend to place them in
      // the TTracker and rely on the global offsets in the mu2e
      // detector frame (note that their local offsets are wrt TTracker)

      _vd->addVirtualDetector( VirtualDetectorId::TT_FrontHollow,
			       ttOffset,
			       0,
			       vdTTFrontOffset);

      // we add another VD detector at the same Z "inside" the proton absorber

      if (c.getBool("hasProtonAbsorber",false)){

	_vd->addVirtualDetector(  VirtualDetectorId::TT_FrontPA,
				 ttOffset, 
				 0,
				 vdTTFrontOffset);
      }

      Hep3Vector vdTTBackOffset(0.,
				0.,
				ttracker.getTrackerEnvelopeParams().zHalfLength()+vdHL);

      _vd->addVirtualDetector( VirtualDetectorId::TT_Back,
			       ttOffset,
			       0,
			       vdTTBackOffset);
    


    }

    // These VDs are related to the beam dump, which is always present
    GeomHandle<ProtonBeamDump> dump;
    if(true) {
      const CLHEP::Hep3Vector vzero(0,0,0);
      
      // This detector will be placed on the face of beam dump
      // shielding.  Computing offsets here is invonvenient since we
      // don't have VolumeInfo for the parent. Just ignore them.
      _vd->addVirtualDetector(VirtualDetectorId::EMFC1Entrance,
			      vzero, 0, vzero
			      );

      // Detector inside the ExtMonFNAL magnet pit, on the face of the upstream wall
      _vd->addVirtualDetector(VirtualDetectorId::EMFC1Exit,
			      dump->magnetPitCenterInEnclosure() + dump->enclosureCenterInMu2e(),
			      &dump->enclosureRotationInMu2e(),
			      /*local position in parent*/
			      CLHEP::Hep3Vector(0, 0, +dump->magnetPitHalfSize()[2] - vdHL)
			      );
      
      // Detector inside the ExtMonFNAL magnet pit, on the face of the downstream wall
      _vd->addVirtualDetector(VirtualDetectorId::EMFC2Entrance,
			      dump->magnetPitCenterInEnclosure() + dump->enclosureCenterInMu2e(),
			      &dump->enclosureRotationInMu2e(),
			      /*local position in parent*/
			      CLHEP::Hep3Vector(0, 0, -dump->magnetPitHalfSize()[2] + vdHL)
			      );
      
    }

    if(geom->hasElement<ExtMonFNAL::ExtMon>()) {
      // Detector inside the ExtMonFNAL detector room, on the face of the upstream wall
      GeomHandle<ExtMonFNAL::ExtMon> extmon;
      _vd->addVirtualDetector(VirtualDetectorId::EMFC2Exit,
			      extmon->roomCenterInMu2e(),
			      &dump->enclosureRotationInMu2e(),
			      /*local position in parent*/
			      CLHEP::Hep3Vector(0, 0, extmon->roomHalfSize()[2] - vdHL)
			      );
    }

    if(geom->hasElement<ExtMonUCI::ExtMon>()) {
      // Detector in front of the ExtMonUCI face shielding
      GeomHandle<ExtMonUCI::ExtMon> extmon;
      vector<double> params = extmon->envelopeParams();

      _vd->addVirtualDetector(VirtualDetectorId::EMIEntrance,
                              extmon->origin(),
                              0,
                              CLHEP::Hep3Vector(0, 0, extmon->envelopeParams()[2] - 0.5*vdHL)
                             );
    }

  }

  VirtualDetectorMaker::~VirtualDetectorMaker (){}

} // namespace mu2e
