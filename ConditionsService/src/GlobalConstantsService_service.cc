// To add a new entity type to the service:
//
//   1) add an entry to the makers_ map using either the generic or
//      your custom defined maker.
//
//   2) explicitly instantiate the getElement() method
//
// Andrei Gaponenko, 2012

#include "ConditionsService/inc/GlobalConstantsService.hh"

#include <typeinfo>
#include <cassert>

#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"

namespace mu2e {

  //================================================================
  namespace {
    // Generic maker for types that can be constructed from SimpleConfig
    template<class T> ConditionsEntity* makeT(const SimpleConfig& config, const std::string&, const std::string&) {
      return new T(config);
    }
  }

  //================================================================
  GlobalConstantsService::GlobalConstantsService(fhicl::ParameterSet const& iPS,
                                                 art::ActivityRegistry&iRegistry)
    : config_(iPS.get<std::string>("inputFile"))
  {
    mf::LogPrint log("GlobalConstants");
    log<<"GlobalConstantsService input file is: "<<config_.inputFile()<<"\n";

    makers_[typeid(mu2e::ParticleDataTable).name()] = &makeT<mu2e::ParticleDataTable>;
    makers_[typeid(mu2e::PhysicsParams).name()]     = &makeT<mu2e::PhysicsParams>;

    //preloadAllEntities();
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
        raw = im->second(config_, key,version);
        entities_.insert(std::make_pair(name, raw));
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
        ConditionsEntity *raw = im->second(config_, "", "");
        entities_.insert(std::make_pair(im->first, raw));
      }
    }
  }

  //================================================================

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::GlobalConstantsService);
