// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/DetectorSolenoidShieldingMaker.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"

// C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  std::unique_ptr<DetectorSolenoidShielding> DetectorSolenoidShieldingMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {

    std::unique_ptr<DetectorSolenoidShielding> dss ( new DetectorSolenoidShielding() );

    const CLHEP::Hep3Vector& dsPos = ds.position();

    // Vacuum pump spool piece
    dss->_vpspCryoSeal.reset( new Tube( c.getString("vpsp.material"),
                                        CLHEP::Hep3Vector( dsPos.x(),
                                                           dsPos.y(),
                                                           ds.cryoZMax()+c.getDouble("vpsp.cryoseal.z") ),
                                        c.getDouble("vpsp.cryoseal.rIn"),
                                        c.getDouble("vpsp.cryoseal.rOut"),
                                        c.getDouble("vpsp.cryoseal.halfLength") ) );

    dss->_dssTubes.push_back( dss->getVPSPCryoSeal() );

    dss->_vpspMain.reset( new Tube( c.getString("vpsp.material"),
                                    CLHEP::Hep3Vector( dsPos.x(),
                                                       dsPos.y(),
                                                       ds.cryoZMax()+c.getDouble("vpsp.main.z") ),
                                    c.getDouble("vpsp.main.rIn" ,ds.rIn1() ),
                                    c.getDouble("vpsp.main.rOut",ds.rIn2()),
                                    c.getDouble("vpsp.main.halfLength") ) );

    dss->_dssTubes.push_back( dss->getVPSPmain() );

    dss->_vpspEndSeal.reset( new Tube( c.getString("vpsp.material"),
                                       CLHEP::Hep3Vector( dsPos.x(),
                                                          dsPos.y(),
                                                          ds.cryoZMax()+c.getDouble("vpsp.endseal.z") ),
                                       c.getDouble("vpsp.endseal.rIn"),
                                       c.getDouble("vpsp.endseal.rOut"),
                                       c.getDouble("vpsp.endseal.halfLength") ) );

    dss->_dssTubes.push_back( dss->getVPSPendSeal() );

    dss->_vpspEndFlange.reset( new Tube( c.getString("vpsp.material"),
                                         CLHEP::Hep3Vector( dsPos.x(),
                                                            dsPos.y(),
                                                            ds.cryoZMax()+c.getDouble("vpsp.endflange.z") ),
                                         c.getDouble("vpsp.endflange.rIn"),
                                         c.getDouble("vpsp.endflange.rOut"),
                                         c.getDouble("vpsp.endflange.halfLength") ) );

    dss->_dssTubes.push_back( dss->getVPSPendFlange() );

    // Instrumentation bulkhead feedthrough
    dss->_ifbMain.reset( new Tube( c.getString("ifb.material"),
                                   CLHEP::Hep3Vector( dsPos.x(),
                                                      dsPos.y(),
                                                      ds.cryoZMax()+c.getDouble("ifb.main.z") ),
                                   c.getDouble("ifb.main.rIn"),
                                   c.getDouble("ifb.main.rOut"),
                                   c.getDouble("ifb.main.halfLength") ) );

    dss->_dssTubes.push_back( dss->getIFBmain() );

    dss->_ifbEndSeal.reset( new Tube( c.getString("ifb.material"),
                                      CLHEP::Hep3Vector( dsPos.x(),
                                                         dsPos.y(),
                                                         ds.cryoZMax()+c.getDouble("ifb.endseal.z") ),
                                      c.getDouble("ifb.endseal.rIn"),
                                      c.getDouble("ifb.endseal.rOut"),
                                      c.getDouble("ifb.endseal.halfLength") ) );

    dss->_dssTubes.push_back( dss->getIFBendSeal() );

    dss->_ifbEndPlug.reset( new Tube( c.getString("ifb.material"),
                                      CLHEP::Hep3Vector( dsPos.x(),
                                                         dsPos.y(),
                                                         ds.cryoZMax()+c.getDouble("ifb.endplug.z") ),
                                      c.getDouble("ifb.endplug.rIn"),
                                      c.getDouble("ifb.endplug.rOut"),
                                      c.getDouble("ifb.endplug.halfLength") ) );

    dss->_dssTubes.push_back( dss->getIFBendPlug() );

    dss->_ifbEndWindow.reset( new Tube( c.getString("ifb.endwindow.material", c.getString("ifb.material") ),
                                        CLHEP::Hep3Vector( dsPos.x(),
                                                           dsPos.y(),
                                                           ds.cryoZMax()+c.getDouble("ifb.endwindow.z") ),
                                        c.getDouble("ifb.endwindow.rIn"),
                                        c.getDouble("ifb.endwindow.rOut"),
                                        c.getDouble("ifb.endwindow.halfLength") ) );

    dss->_dssTubes.push_back( dss->getIFBendWindow() );
    /*
    dss->_ifbEndWindowFrameInside.reset( new Tube( c.getString("ifb.endwindowFrameInside.material", c.getString("ifb.material") ),
                                        CLHEP::Hep3Vector( dsPos.x(),
                                                           dsPos.y(),
                                                           ds.cryoZMax()+c.getDouble("ifb.endwindow.z") ),
                                        c.getDouble("ifb.endwindowFrameInside.rIn"),
                                        c.getDouble("ifb.endwindowFrameInside.rOut"),
                                        c.getDouble("ifb.endwindowFrameInside.halfLength") ) );

    dss->_dssTubes.push_back( dss->getIFBendWindowFrameInside() );
    */
    // Set names of Tubes
    dss->_vpspCryoSeal ->setName("VPSP_CryoSeal" );
    dss->_vpspMain     ->setName("VPSP_Main"     );
    dss->_vpspEndSeal  ->setName("VPSP_EndSeal"  );
    dss->_vpspEndFlange->setName("VPSP_EndFlange");

    dss->_ifbMain     ->setName("IFB_Main"     );
    dss->_ifbEndSeal  ->setName("IFB_EndSeal"  );
    dss->_ifbEndPlug  ->setName("IFB_EndPlug"  );
    dss->_ifbEndWindow->setName("IFB_EndWindow");

    return dss;

  } // make()

} // namespace mu2e
