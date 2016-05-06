//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_module.cc,v 1.82 2014/03/24 21:39:01 gandr Exp $
// $Author: gandr $
// $Date: 2014/03/24 21:39:01 $
//
// Original author Rob Kutschke

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <stdexcept>

namespace mu2e {
  class G4 : public art::EDProducer {
  public:
    G4(const fhicl::ParameterSet&) {}
    virtual void produce(art::Event&) {
      throw std::runtime_error("\n\nG4_module is obsolete and not functional any more.  Please upgrade your job configuration to use Mu2eG4 module instead.\n\n\n");
    }
  };
} // End of namespace mu2e

using mu2e::G4;
DEFINE_ART_MODULE(G4);
