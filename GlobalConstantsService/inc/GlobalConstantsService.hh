// Andrei Gaponenko, 2012

#ifndef GlobalConstantsService_GlobalConstantsService_hh
#define GlobalConstantsService_GlobalConstantsService_hh

#include <string>

#include "boost/shared_ptr.hpp"
#include "boost/utility.hpp"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e {

  template <typename ENTITY> class GlobalConstantsHandle;

  class GlobalConstantsService : boost::noncopyable {

  public:
    GlobalConstantsService(const fhicl::ParameterSet&, art::ActivityRegistry&);

    // The class is intended to be used via GlobalConstantsHandle
    // No public accessors are provided.

    // Normally the entities are created on demand.
    // The preload method is only needed for
    //
    //   1) to quickly check that all defined entities can be loaded
    //
    //   2) in cases when timing is important, pre-loading allows to
    //      avoid the first access delay inside the event loop.

    void preloadAllEntities();

    void postEndJob();

  private:

    template <typename ENTITY> friend class GlobalConstantsHandle;

    typedef boost::shared_ptr<ConditionsEntity> ConditionsEntityPtr;
    typedef std::map<std::string,ConditionsEntityPtr> ConditionsMap;

    typedef ConditionsEntity* (*MakerPtr)(const SimpleConfig& config, const std::string& key, const std::string& version);
    typedef std::map<std::string, MakerPtr> FactoryMap;

    int  configStatsVerbosity_;
    SimpleConfig config_;
    FactoryMap makers_;
    mutable ConditionsMap entities_;

    template <class ENTITY> const ENTITY* getElement(const std::string& key, const std::string& version) const;
  };

}

DECLARE_ART_SERVICE(mu2e::GlobalConstantsService, SHARED)
#endif /* GlobalConstantsService_GlobalConstantsService_hh */
