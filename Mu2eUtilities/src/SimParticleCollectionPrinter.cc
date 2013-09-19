// Andrei Gaponenko, 2013

#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {
  //================================================================
  SimParticleCollectionPrinter::SimParticleCollectionPrinter(const fhicl::ParameterSet& pset)
    : prefix_(pset.get<std::string>("prefix", ""))
    , enabled_(pset.get<bool>("enabled", true))
    , primariesOnly_(pset.get<bool>("primariesOnly", false))
  {}

  //================================================================
  std::ostream& SimParticleCollectionPrinter::print(std::ostream& os, const SimParticle& p) {
    os<<" id = "<<p.id()
      <<", parent="<<p.parentId()
      <<", pdgId="<<p.pdgId()
      <<", p="<<p.startMomentum()
      <<", start="<<p.startPosition()
      <<", end="<<p.endPosition()
      <<", stoppingCode = "<<p.stoppingCode();
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
  fhicl::ParameterSet SimParticleCollectionPrinter::defaultPSet() {
    fhicl::ParameterSet res;
    res.put<bool>("enabled", false);
    return res;
  }

  //================================================================
}
