//
// A Producer Module that runs Geant4 and adds its output to the event.
// ******Meant for Geant4 Studies not for Mu2e Simulations**********
//
// $Id: Mu2eG4Study_module.cc,v 1.9 2013/12/17 21:49:06 genser Exp $
// $Author: genser $
// $Date: 2013/12/17 21:49:06 $
//
// Original author K. Genser, based on Rob's G4_module
//
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.
//


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <stdexcept>

namespace mu2e {

  class Mu2eG4Study : public art::EDProducer {
  public:
    Mu2eG4Study(fhicl::ParameterSet const& pSet) {
      throw std::runtime_error("\n\nMu2eG4Study module is very obsolete. Please upgrade your job configuration to use the modern mainstream Mu2eG4 module instead.\nStarting with v5_7_7 and still as of v6_2_3 Mu2eG4/fcl/g4study2.fcl and g4study2Calo_01.fcl are examples on how to do it.\n\n\n");
    }
    virtual void produce(art::Event& e) override {}
  };

} // End of namespace mu2e

using mu2e::Mu2eG4Study;
DEFINE_ART_MODULE(Mu2eG4Study);
