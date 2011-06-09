//
// Construct and return an VirtualDetector.
//
// $Id: VirtualDetectorMaker.cc,v 1.8 2011/06/09 19:57:53 genser Exp $
// $Author: genser $
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
#include "BeamlineGeom/inc/Beamline.hh"
#include "TargetGeom/inc/Target.hh"
#include "TTrackerGeom/inc/TTracker.hh"

using namespace std;
using namespace CLHEP;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  VirtualDetectorMaker::VirtualDetectorMaker(SimpleConfig const& c)
  {
    _vd = auto_ptr<VirtualDetector>(new VirtualDetector());

    if( ! c.getBool("hasVirtualDetector",false) ) return;

    double vdHL = c.getDouble("vd.halfLength",0.1*mm);
    _vd->_halfLength = vdHL;

    // Need some data from other subsystems
    GeomHandle<Beamline> bg;
    double solenoidOffset = bg->solenoidOffset();

    // VD 1 and 2 are at the front and back of collimator 1, which is placed inside TS1.

    //double ts1HL   = bg->getTS().getTS1().getHalfLength();
    double coll1HL = bg->getTS().getColl1().getHalfLength();

    HepRotation *ts1rot = bg->getTS().getTS1().getRotation();
    Hep3Vector   ts1pos = bg->getTS().getTS1().getGlobal();
    Hep3Vector coll1pos = bg->getTS().getColl1().getLocal();
    Hep3Vector deltaZ1(0,0,coll1HL-vdHL);

    _vd->addVirtualDetector( 1, "Coll1_In",  ts1pos, ts1rot, coll1pos-deltaZ1);
    _vd->addVirtualDetector( 2, "Coll1_Out", ts1pos, ts1rot, coll1pos+deltaZ1);

    // VD 3-6 are placed around two sections of collimator 3

    //double ts3HL   = bg->getTS().getTS1().getHalfLength();
    double coll31HL = bg->getTS().getColl31().getHalfLength();
    double coll32HL = bg->getTS().getColl32().getHalfLength();

    HepRotation *ts3rot = bg->getTS().getTS3().getRotation();
    Hep3Vector   ts3pos = bg->getTS().getTS3().getGlobal();

    Hep3Vector coll31pos = bg->getTS().getColl31().getLocal();
    Hep3Vector coll32pos = bg->getTS().getColl32().getLocal();
    Hep3Vector deltaZ31(0,0,coll31HL-vdHL);
    Hep3Vector deltaZ32(0,0,coll32HL-vdHL);

    _vd->addVirtualDetector( 3, "Coll31_In",  ts3pos, ts3rot, coll31pos-deltaZ31);
    _vd->addVirtualDetector( 4, "Coll31_Out", ts3pos, ts3rot, coll31pos+deltaZ31);
    _vd->addVirtualDetector( 5, "Coll32_In",  ts3pos, ts3rot, coll32pos-deltaZ32);
    _vd->addVirtualDetector( 6, "Coll32_Out", ts3pos, ts3rot, coll32pos+deltaZ32);

    // VD 7 and 8 are at the front and back of collimator 5, which is placed inside TS5.

    //double ts5HL   = bg->getTS().getTS5().getHalfLength();
    double coll5HL = bg->getTS().getColl5().getHalfLength();

    HepRotation *ts5rot = bg->getTS().getTS5().getRotation();
    Hep3Vector   ts5pos = bg->getTS().getTS5().getGlobal();
    Hep3Vector coll5pos = bg->getTS().getColl5().getLocal();
    Hep3Vector deltaZ5(0,0,coll5HL-vdHL);

    _vd->addVirtualDetector( 7, "Coll5_In",  ts5pos, ts5rot, coll5pos-deltaZ5);
    _vd->addVirtualDetector( 8, "Coll5_Out", ts5pos, ts5rot, coll5pos+deltaZ5);

    // VD 9 and 10 are placed inside DS2, just before and after stopping target

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

    _vd->addVirtualDetector(  9, "ST_In",  ds2Offset, 0, targetOffset-shift);
    _vd->addVirtualDetector( 10, "ST_Out", ds2Offset, 0, targetOffset+shift);

    if (c.getBool("hasTTracker",false)){

      ostringstream vdName("TT_Mid");

      if(c.getInt("ttracker.numDevices")%2!=0){
        throw cet::exception("GEOM")
          << "This virtual detector " << vdName 
          << " can only be placed if the TTracker has an even number of devices \n";
      }

      TTracker const & ttracker = *(GeomHandle<TTracker>());
      Hep3Vector ttrackerOffset(-solenoidOffset,0.,ttracker.z0());

      // VD 11 is placed inside the ttracker mother volume in the middle of the ttracker
      // VD 12 is placed inside the ttracker at the same z position as
      // VD 11 but from radius 0 to the inner radius of the ttracker
      // mother volume. However, its mother volume is ToyDS3Vacuum
      // which has a different offset. We will use the global offset
      // here (!) as DS is not in the geometry service yet

      _vd->addVirtualDetector( 11, vdName.str(),  ttrackerOffset, 0, Hep3Vector(0.,0.,0.));
      _vd->addVirtualDetector( 12, "TT_MidInner", ttrackerOffset, 0, Hep3Vector(0.,0.,0.));

      vdName.str("TT_FrontHollow");

      Hep3Vector ttrackerFrontVDOffset(-solenoidOffset,
				       0.,
				       ttracker.z0()-ttracker.getTrackerEnvelopeParams().zHalfLength()-vdHL);

      // VD 13 is placed outside the ttracker mother volume in front of the ttracker
      // "outside" of the proton absorber
      
      _vd->addVirtualDetector( 13,
			       vdName.str(),
			       ttrackerFrontVDOffset,
			       0,
			       Hep3Vector(0.,0.,0.));

      // we add another VD detector at the same Z "inside" the proton absorber

      if (c.getBool("hasProtonAbsorber",false)){

	vdName.str("TT_FrontPA");
	_vd->addVirtualDetector( 14,
				 vdName.str(),
				 ttrackerFrontVDOffset, 
				 0,
				 Hep3Vector(0.,0.,0.));
      }
    }
  }

  VirtualDetectorMaker::~VirtualDetectorMaker (){}

} // namespace mu2e
