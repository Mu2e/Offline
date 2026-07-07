// Geometry of the detector solenoid downstream shielding
//
// Original author: Kyle Knoepfel, 2013

#ifndef DETECTORSOLENOIDSHIELDING_HH
#define DETECTORSOLENOIDSHIELDING_HH

// #include "art/Persistency/Common/Wrapper.h"

#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include <memory>
#include <vector>

namespace mu2e {

  class DetectorSolenoidShielding : virtual public Detector {
  public:

    // VPSP
    const Tube* getVPSPCryoSeal()  const { return _vpspCryoSeal.get();  }
    const Tube* getVPSPmain()      const { return _vpspMain.get();      }
    const Tube* getVPSPendSeal()   const { return _vpspEndSeal.get();   }
    const Tube* getVPSPendFlange() const { return _vpspEndFlange.get(); }

    // IFB
    const Tube* getIFBmain()       const { return _ifbMain.get();      }
    const Tube* getIFBendSeal()    const { return _ifbEndSeal.get();   }
    const Tube* getIFBendPlug()    const { return _ifbEndPlug.get();   }
    const Tube* getIFBendWindow()  const { return _ifbEndWindow.get(); }
    const Tube* getIFBendWindowFrameInside()  const { return _ifbEndWindowFrameInside.get(); }
    const Tube* getIFBendWindowFrameOutside()  const { return _ifbEndWindowFrameOutside.get(); }

    const std::vector<const Tube*> getTubes() const { return _dssTubes; }

    //----------------------------------------------------------------
  private:
    friend class DetectorSolenoidShieldingMaker;

    // Private ctr: the class should be only obtained via
    // DetectorSolenoidShielding::DetectorSolenoidShieldingMaker.
    DetectorSolenoidShielding(){}

    // VPSP
    std::unique_ptr<Tube> _vpspCryoSeal;
    std::unique_ptr<Tube> _vpspMain;
    std::unique_ptr<Tube> _vpspEndSeal;
    std::unique_ptr<Tube> _vpspEndFlange;

    // IFB
    std::unique_ptr<Tube> _ifbMain;
    std::unique_ptr<Tube> _ifbEndSeal;
    std::unique_ptr<Tube> _ifbEndPlug;
    std::unique_ptr<Tube> _ifbEndWindow;
    std::unique_ptr<Tube> _ifbEndWindowFrameInside;
    std::unique_ptr<Tube> _ifbEndWindowFrameOutside;

    // Map for both
    std::vector<const Tube*> _dssTubes;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    DetectorSolenoidShielding() {}
  };
}

#endif/*DETECTORSOLENOIDSHIELDING_HH*/
