// Look through the generated particles in an event and define which one was the 'primary' (signal-like)
// There should be at most one of these.  In case of none, a null primary is added to the event
// Original author: David Brown (LBNL) Jan 2019
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e 
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include <vector>
#include <iostream>
#include <string>
using CLHEP::HepLorentzVector;
using CLHEP::Hep3Vector;
namespace mu2e {
  class FindMCPrimary : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
	fhicl::Atom<art::InputTag> genPC{  Name("GenParticles"),
	  Comment("GenParticle collection containing the primary")};
	fhicl::Atom<art::InputTag> simPC{  Name("SimParticles"),
	  Comment("SimParticle collection")};
	fhicl::Sequence<std::string> genIDs { Name("PrimaryGenIds"),
	  Comment("Generator IDs of potential Primary Particles")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit FindMCPrimary(const Parameters& conf);
      void produce(art::Event& evt) override;

    private:
      int _debug;
      art::InputTag _gpc, _spc;
      std::vector<int> _pgenids; 
  };

  FindMCPrimary::FindMCPrimary(const Parameters& config )  : 
    _debug(config().debug()),
    _gpc(config().genPC()),
    _spc(config().simPC())
  {
    consumes<GenParticleCollection>(_gpc);
    consumes<SimParticleCollection>(_spc);
    produces <PrimaryParticle>(); 
    for(auto const& genid : config().genIDs()) {
      _pgenids.push_back(GenId::findByName(genid).id());
      if(_debug > 0)
	std::cout << "GenId " << genid << " value " << _pgenids.back() << " defined as primary " << std::endl;
    }
  }


  void FindMCPrimary::produce(art::Event& event) {
    // create output: by default, this is null
    GenParticle pgp(PDGCode::null, GenId::unknown,CLHEP::Hep3Vector(),
	CLHEP::HepLorentzVector(), 0.0,0.0);
    // find input
    auto gpch = event.getValidHandle<GenParticleCollection>(_gpc);
    auto const& gpc = *gpch;
    auto spch = event.getValidHandle<SimParticleCollection>(_spc);
    auto const& spc = *spch;
    // loop over input gen particles and save those whose GenIds match the list
    std::vector<GenParticleCollection::const_iterator> pgps;
    for(auto igp= gpc.begin(); igp != gpc.end(); igp++) {
      if(std::find(_pgenids.begin(),_pgenids.end(),igp->generatorId().id()) != _pgenids.end()){
	pgps.push_back(igp);
      }
    }
    // now examine the list of potential primaries
    if (pgps.size() == 1 ) {
      // exactly 1 primary: we're done
      pgp = *pgps.front();
      if(_debug > 1) std::cout << "Found single particle primary " << pgp << std::endl;
    } else if (pgps.size() >= 2){ 
      // internal RMC and RPC generates 2 'primary' particles.  Create the true primary
      // virtual photon from these
      if(pgps.size() == 2 && pgps[0]->generatorId() == pgps[1]->generatorId() &&
	  (pgps[0]->generatorId() == GenId::InternalRPC ||
	   pgps[0]->generatorId() == GenId::InternalRMC ) ) {
	// double -check consistency
	if(fabs(pgps[0]->time() - pgps[1]->time()) > 1e-6)
	  throw cet::exception("Simulation")<<"FindMCPrimary: RMC/RPC origins don't match" << std::endl;
	// Add the 4vectors of the 2 consistuent electrons
	HepLorentzVector momsum = pgps[0]->momentum() + pgps[1]->momentum();
	// update the primary; particle type is (offshell) photon
	pgp = GenParticle(PDGCode::gamma,pgps[0]->generatorId(),
	    pgps[0]->position(),momsum, pgps[0]->time(), pgps[0]->properTime());
	if(_debug > 1) std::cout << "Found 2-particle primary " << pgp << std::endl;

      } else {
	// for now, assume this is an error
	throw cet::exception("Simulation")<<"FindMCPrimary: Multiple Primaries" << std::endl;
      }
    } else if(_debug > 1)
      std::cout << "Found no particle primary " << pgp << std::endl;
    // now find the SimParticles associated with these GenParticles.  There is no forwards map,
    // so I just have to exhaustively search
    PrimaryParticle::SPPV simps;
    for(auto igp=gpc.begin(); igp != gpc.end(); ++igp) {
      // find the index for this iterator
      size_t gpindex = std::distance(gpc.begin(),igp);
      // create a ptr for this GenParticle
      auto gpp = art::Ptr<GenParticle>(gpch,gpindex);
      // find the sim particle giving this ptr as its GenParticle
      for(auto isp=spc.begin(); isp != spc.end(); ++isp) {
	if(isp->second.genParticle() == gpp){
      // create a ptr for this SimParticle and remember it
	  size_t spindex = std::distance(spc.begin(),isp);
	  auto spp = art::Ptr<SimParticle>(spch,spindex);
	  simps.push_back(spp);
	}
      }
    }
    // check consistency
    if(pgps.size() != simps.size())
      throw cet::exception("Simulation")<<"FindMCPrimary: GenParticle <-> SimParticle inconsistency" << std::endl;
    // create output object
    PrimaryParticle pp(pgp, simps);
    // put in event
    event.put(std::make_unique<PrimaryParticle>(pp));
  }
}
DEFINE_ART_MODULE(mu2e::FindMCPrimary)
