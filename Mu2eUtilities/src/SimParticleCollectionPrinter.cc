// Andrei Gaponenko, 2013

#include <utility>                            // for pair

#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "CLHEP/Vector/LorentzVector.h"       // for operator<<
#include "CLHEP/Vector/ThreeVector.h"         // for operator<<
#include "MCDataProducts/inc/ProcessCode.hh"  // for operator<<
#include "cetlib/map_vector.h"                // for operator<<

namespace mu2e {
  //================================================================
  SimParticleCollectionPrinter::SimParticleCollectionPrinter(const Config& conf)
    : prefix_(conf.prefix())
    , enabled_(conf.enabled())
    , primariesOnly_(conf.primariesOnly())
  {}

  //================================================================
  SimParticleCollectionPrinter::SimParticleCollectionPrinter()
    : prefix_("")
    , enabled_(false)
    , primariesOnly_(false)
  {}

  //================================================================
  std::ostream& SimParticleCollectionPrinter::print(std::ostream& os, const SimParticle& p) {
    os<<" id = "<<p.id()
      <<", parent="<<p.parentId()
      <<", pdgId="<<p.pdgId()
      <<", start="<<p.startPosition()
      <<", pstart="<<p.startMomentum()
      <<", end="<<p.endPosition()
      <<", pend="<<p.endMomentum()
      <<", nSteps = "<<p.nSteps()
      <<", preLastStepKE = "<<p.preLastStepKineticEnergy()
      <<", creationCode = "<<p.creationCode()
      <<", stoppingCode = "<<p.stoppingCode()
      <<", startG4Status = "<<p.startG4Status()
      <<", endG4Status = "<<p.endG4Status();
    return os;
  }

  //================================================================
  std::ostream& SimParticleCollectionPrinter::print(std::ostream& os, const SimParticleCollection& particles) const {
    if(enabled_) {
      for(const auto& i: particles) {
        const SimParticle& p{i.second};
        if(!primariesOnly_ || p.isPrimary()) {
          os<<prefix_;
          print(os, p);
          os<<std::endl;
        }
      }
    }

    return os;
  }

  //================================================================
}
