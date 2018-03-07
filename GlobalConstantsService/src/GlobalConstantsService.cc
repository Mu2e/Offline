// To add a new entity type to the service:
//
//   1) add an entry to the makers_ map using either the generic or
//      your custom defined maker.
//
//   2) explicitly instantiate the getElement() method
//
// Andrei Gaponenko, 2012

#include "GlobalConstantsService/inc/GlobalConstantsService.hh"

#include <iostream>
#include <typeinfo>
#include <cassert>

#include "cetlib_except/exception.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {

  //================================================================
  namespace {
    // Generic maker for types that can be constructed from SimpleConfig
    template<class T> ConditionsEntity* makeT(const SimpleConfig& config, const std::string&, const std::string&) {
      return new T(config);
    }
  }

  //================================================================
  GlobalConstantsService::GlobalConstantsService(fhicl::ParameterSet const& pset,
                                                 art::ActivityRegistry&iRegistry)
    : configStatsVerbosity_((pset.get<int>("configStatsVerbosity", 0))),
      config_(pset.get<std::string>("inputFile"),
              pset.get<bool>("allowReplacement",     true),
              pset.get<bool>("messageOnReplacement", false),
              pset.get<bool>("messageOnDefault",     false)
             )
  {

    config_.printOpen(std::cout,"GlobalConstants");
    iRegistry.sPostEndJob.watch (this, &GlobalConstantsService::postEndJob );

    if ( pset.get<bool>("printConfig",false) ) {
      config_.print(std::cout,"GlobalConstants: ");
    }

    makers_[typeid(mu2e::ParticleDataTable).name()] = &makeT<mu2e::ParticleDataTable>;
    makers_[typeid(mu2e::PhysicsParams).name()]     = &makeT<mu2e::PhysicsParams>;


    //preloadAllEntities();
  }

  void GlobalConstantsService::postEndJob(){
    config_.printAllSummaries( std::cout, configStatsVerbosity_, "GlobalConstants: " );
  }

  //================================================================
  // For now the key and version arguments are ignored.
  template <class ENTITY>
  const ENTITY* GlobalConstantsService::getElement(const std::string& key, const std::string& version) const
  {
    ConditionsEntity *raw(0);

    std::string name = typeid(ENTITY).name();
    ConditionsMap::const_iterator it(entities_.find(name));
    if(it==entities_.end()) {
      FactoryMap::const_iterator im(makers_.find(name));
      if(im==makers_.end()) {
        throw cet::exception("GEOM")
          << "GlobalConstantsService: unknown entity type "<<name<< "\n";
      }
      else {
        std::pair<ConditionsMap::iterator, bool>
          insert_result =
          entities_.insert(ConditionsMap::value_type(name,
                                                     ConditionsMap::mapped_type(im->second(config_, key,version))));
        raw = insert_result.first->second.get();
      }
    }
    else {
      raw = dynamic_cast<ENTITY*>(it->second.get());
    }

    assert(raw);

    const ENTITY *res(dynamic_cast<ENTITY*>(raw));

    assert(res);

    return res;
  }

  // Explicitly instantiate getElement() for all types we handle rather than exposing all the code in the header.
  template const ParticleDataTable* GlobalConstantsService::getElement<ParticleDataTable>(const std::string&, const std::string&) const;
  template const PhysicsParams* GlobalConstantsService::getElement<PhysicsParams>(const std::string&, const std::string&) const;

  //================================================================
  void GlobalConstantsService::preloadAllEntities() {
    for(FactoryMap::const_iterator im=makers_.begin(); im != makers_.end(); ++im) {
      if(entities_.find(im->first) == entities_.end()) {
        entities_.insert(ConditionsMap::value_type(im->first,
                                                   ConditionsMap::mapped_type(im->second(config_, "", ""))));
      }
    }
  }

  //================================================================

} // end namespace mu2e
