// Look through the SimParticles in an event and define which one was the 'primary' (signal-like)
// There should be at most one of these.  In case of none, a null primary is added to the event
// Original author: David Brown (LBNL) Jan 2019
// June 2021 rewritten for art-based muon stops D. Brown
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
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
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> simPC{  Name("SimParticles"), Comment("SimParticle collection")};
        fhicl::OptionalAtom<std::string> primaryProcess{ Name("PrimaryProcess"), Comment("Process code that produced the primary physics particle") };
        fhicl::OptionalSequence<std::string> primaryGenIds{ Name("PrimaryGenIds"), Comment("Generator codes that produced the primary physics particle") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit FindMCPrimary(const Parameters& conf);
      void produce(art::Event& evt) override;
    private:
      int _debug;
      art::InputTag _spc;
      ProcessCode _pcode;
      std::vector<GenId> _gcodes;
  };

  FindMCPrimary::FindMCPrimary(const Parameters& config )  :
    art::EDProducer{config},
    _debug(config().debug()),
    _spc(config().simPC())
    {
      std::string pcode;
      std::vector<std::string> gcodes;
      if( config().primaryProcess(pcode)) _pcode = ProcessCode::findByName(pcode);
      if( config().primaryGenIds(gcodes)){
        for(auto const& gcodename : gcodes) {
          _gcodes.emplace_back(GenId::findByName(gcodename));
        }
      }
      if(_pcode != ProcessCode::uninitialized && _gcodes.size() > 0 )
        throw cet::exception("Configuration") << "Both process and gen codes specified " << std::endl;
      else if(_pcode == ProcessCode::uninitialized && _gcodes.size() == 0)
        throw cet::exception("Configuration") << "Neither process or gen codes specified " << std::endl;
      consumes<SimParticleCollection>(_spc);
      produces <PrimaryParticle>();
      if(_debug > 0){
        if(_pcode != ProcessCode::uninitialized)
          std::cout << "Looking for primary SimParticles with code " << _pcode << std::endl;
        else if(_gcodes.size() > 0){
          std::cout << "Looking for primary GenParticles with code: " << std::endl;
          for( auto const& gcode : _gcodes) {
            std::cout << gcode << std::endl;
          }
        }
      }
    }

  void FindMCPrimary::produce(art::Event& event) {
    typedef std::vector<art::Ptr<SimParticle> > SPPV;
    SPPV sims;
    // find input
    auto spch = event.getValidHandle<SimParticleCollection>(_spc);
    auto const& spc = *spch;
    for(auto isp = spc.begin(); isp != spc.end(); ++isp){
      // first, search for particles coming from a primary process
      if(isp->second.creationCode() == _pcode){
        sims.emplace_back(spch, isp->first.asInt());
        if(_debug > 1)std::cout << "found primary SimParticle, code " << _pcode << std::endl;
      }
      // then, by creation code
      if(isp->second.genParticle().isNonnull()){
        for(auto const& gcode : _gcodes) {
          if(isp->second.genParticle()->generatorId() == gcode){
            sims.emplace_back(spch, isp->first.asInt());
            if(_debug > 1)std::cout << "found primary GenParticle, code " << gcode << std::endl;
          }
        }
      }
    }
    // if no primary was found, throw
    if(sims.size() == 0)throw cet::exception("Simulation") << " No Primary particle found " << std::endl;
    // check these all have the same creation code
    if(sims.size() > 1){
      auto isim = sims.begin()++;
      while(isim != sims.end()){
        if((*isim)->creationCode() != sims.front()->creationCode())
          throw cet::exception("Simulation")<<"PrimaryParticle: creation codes don't match" << std::endl;
        ++isim;
      }
    }
    // create output object; this checks for consistency
    PrimaryParticle pp(sims);
    event.put(std::make_unique<PrimaryParticle>(pp));
  }
}
DEFINE_ART_MODULE(mu2e::FindMCPrimary)
