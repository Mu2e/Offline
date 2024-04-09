//-----------------------------------------------------------------------------
// P.M. : keep is simple - generic generator module which uses plugins
//        to generate various signal processes
//        probably could do for pions as well
//-----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

#include "canvas/Utilities/InputTag.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

namespace mu2e {
//================================================================
  class EventGeneratorMu : public art::EDProducer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<int>           stopPdgCode   {Name("stopPdgCode"   ),Comment("PDG code of the stop"       )};
      fhicl::Atom<double>        lifetime      {Name("lifetime"      ),Comment("mu+ lifetime"               )};
      fhicl::Atom<double>        tmin          {Name("tmin"          ), Comment("T min"                     )};
      fhicl::Atom<double>        tmax          {Name("tmax"          ), Comment("T max"                     )};
      fhicl::Atom<art::InputTag> simpCollTag   {Name("simpCollTag"   ),Comment("input SimParticleCollection")};
      fhicl::DelegatedParameter  generator     {Name("generator"     ),Comment("generator tool"             )};
      fhicl::Atom<int>           verbosityLevel{Name("verbosityLevel"),Comment("verbosity Level"            )};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit EventGeneratorMu(const Parameters& conf);

    virtual void produce(art::Event& event) override;

//----------------------------------------------------------------
  private:
    PDGCode::type                              _stopPdgCode;
    double                                     _lifetime   ;          // mu+ lifetime
    double                                     _tmin       ;
    double                                     _tmax       ;
    const art::InputTag                        _simpCollTag;
    art::RandomNumberGenerator::base_engine_t& _eng        ;
    std::unique_ptr<CLHEP::RandExponential>    _randExp    ;
    std::unique_ptr<ParticleGeneratorTool>     _generator  ;
    int                                        _verbosityLevel;
  };

//================================================================
  EventGeneratorMu::EventGeneratorMu(const Parameters& conf) : EDProducer{conf}
    , _lifetime    (conf().lifetime()   )
    , _tmin        (conf().tmin       ())
    , _tmax        (conf().tmax       ())
    , _simpCollTag (conf().simpCollTag())
    , _eng         (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , _randExp     (std::make_unique<CLHEP::RandExponential>(_eng))
    , _verbosityLevel(conf().verbosityLevel())
  {
    _stopPdgCode = static_cast<PDGCode::type>(conf().stopPdgCode());

    if (_simpCollTag.empty()) {
      produces <mu2e::GenParticleCollection>();
      produces <mu2e::SimParticleCollection>();
    }
    else {
      consumes<SimParticleCollection>(_simpCollTag);
      produces<mu2e::StageParticleCollection>();
    }

    const auto pset = conf().generator.get<fhicl::ParameterSet>();

    _generator      = art::make_tool<ParticleGeneratorTool>(pset);
    _generator->finishInitialization(_eng,"wow");

    if(_verbosityLevel > 0) {
      mf::LogInfo log("EventGeneratorMu");
      log << ", lifetime = " << _lifetime << std::endl;
    }
  }

//-----------------------------------------------------------------------------
  void EventGeneratorMu::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
//-----------------------------------------------------------------------------
// practice of passing vectors by value is sprawling
//-----------------------------------------------------------------------------
    std::vector<art::Ptr<SimParticle>> list;

    if (not _simpCollTag.empty()) {
      art::Handle<SimParticleCollection> simpch;
      event.getByLabel<SimParticleCollection>(_simpCollTag,simpch);
      if (not simpch.isValid()) {
        printf("EventGeneratorMu::%s ERROR: can\'t find SimParticleCollection tag=%s, BAIL OUT\n",
               __func__,_simpCollTag.encode().data());
        return;
      }
//-----------------------------------------------------------------------------
// fill list of stopped particles with given PDG code
//-----------------------------------------------------------------------------
      for(auto i = simpch->begin(); i != simpch->end(); ++i) {
        const auto& inpart = i->second;
        if((inpart.pdgId() == _stopPdgCode) && (inpart.endMomentum().vect().mag2() == 0)) {
          list.emplace_back(simpch, i->first.asInt());
        }
      }
//-----------------------------------------------------------------------------
// for each stop (vertex) generate secondary particle(s)
// decay_time = stop_time + lifetime
// for the antiptoton annihilation, will set the _lifetime to zero
//-----------------------------------------------------------------------------
      for(const auto& stop: list) {
        double decay_time = stop->endGlobalTime()+_randExp->fire(_lifetime);
//-----------------------------------------------------------------------------
// re-package daughters into a list of "stage particles"
//-----------------------------------------------------------------------------
        auto daughters = _generator->generate();
        for(const auto& d: daughters) {
          output->emplace_back(stop                     ,
                               _generator->processCode(),
                               d.pdgId                  ,
                               stop->endPosition()      ,
                               d.fourmom                ,
                               decay_time               );
        }
      }
    }
    else {
//-----------------------------------------------------------------------------
// particle gun mode: need to produce the stop
// --------------------
// 1. create a GenParticle, store it in the event, and get a handle
//-----------------------------------------------------------------------------
      std::unique_ptr<GenParticleCollection> genp_collp(new GenParticleCollection);

      CLHEP::Hep3Vector       pos ;
      double                  tstop(0);
      CLHEP::HepLorentzVector mom(0,0,0,0) ;

      _generator->getXYZ(&pos);

      GenParticle genp(_stopPdgCode,_generator->genId(),pos,mom,tstop);

      genp_collp->push_back(genp);
      auto genpch = event.put(std::move(genp_collp));

      GenParticleCollection* genpc =  (GenParticleCollection*) genpch.product();
//-----------------------------------------------------------------------------
// 2. create a fake SimParticle at rest corresponding to the GenParticle
//    and store it in the event
//-----------------------------------------------------------------------------
      std::unique_ptr<SimParticleCollection> simp_collp(new SimParticleCollection);
      art::Ptr<GenParticle>::key_type k9(1);
      art::Ptr<GenParticle> genpPtr = art::Ptr<GenParticle>(genpch,k9);

      SimParticle simp(SimParticle::key_type(),
                       1,
                       art::Ptr<SimParticle>(), // stopped particle has no parent
                       _stopPdgCode,
                       genpPtr,
                       CLHEP::HepLorentzVector(pos,tstop),
                       mom,
                       0,
                       0,
                       0,
                       0,
                       _generator->processCode());

      simp.addEndInfo(pos,mom,tstop,0,0,0,_generator->processCode(),0.,0,0.);

      cet::map_vector_key o_my_gosh{1};
      simp_collp->insert(std::make_pair(o_my_gosh,simp));
      auto simpc_ph = event.put(std::move(simp_collp));

      art::Ptr<SimParticle>::key_type s9(1);
      art::Ptr<SimParticle> stop = art::Ptr<SimParticle>(simpc_ph,s9);
      list.push_back(stop);

      for(const auto& stop: list) {
        double decay_time = _tmin;
        if (_tmax > _tmin) {
          do { decay_time = _tmin + _randExp->fire(_lifetime); } while (decay_time > _tmax);
        }
//-----------------------------------------------------------------------------
// re-package daughters into a list of "stage particles"
//-----------------------------------------------------------------------------
        auto daughters = _generator->generate();
        for(const auto& d: daughters) {
          genpc->push_back(GenParticle(d.pdgId,_generator->genId(),stop->endPosition(),d.fourmom,decay_time));
        }
      }
    }

    if(_verbosityLevel >= 9) {
      std::cout<<"EventGeneratorMu output: "<<*output<<std::endl;
    }

    event.put(std::move(output));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EventGeneratorMu)
