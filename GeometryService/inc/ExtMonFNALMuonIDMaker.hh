//Maker File for Muon ID Detector Box

//Jackson Waters 2018

#ifndef EXTMONFNALMUONIDMAKER_HH
#define EXTMONFNALMUONIDMAKER_HH

#include <memory>
#include <string>

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMuonID.hh"

namespace CLHEP { class Hep3Vector; }
namespace CLHEP { class HepRotation; }

namespace mu2e { class SimpleConfig;        }
namespace mu2e { class ExtMonFNALMuonID;     }
namespace mu2e { namespace ExtMonFNAL { class ExtMon; } }

namespace mu2e {

  class ExtMonFNALMuonIDMaker {
     public:

    static ExtMonFNALMuonID read(const SimpleConfig& c,
                                 const std::string& prefix,
                                 const CLHEP::HepRotation& muonIDInRotationInMu2e,
                                 const CLHEP::Hep3Vector& refTrajmuonIDEntranceInMu2e,
    double nominalMomentum);

    static std::unique_ptr<ExtMonFNALMuonID> make(const SimpleConfig& config);
    };

  };

#endif/*EXTMONFNALMUONIDMAKER_HH*/


