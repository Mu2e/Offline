#ifndef DbExample_ConditionsService2_hh
#define DbExample_ConditionsService2_hh

//
//
// C++ include files
#include <string>

// Framework include files
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "cetlib_except/exception.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "DbExample/inc/ConditionsEntity2.hh"
#include "DbExample/inc/ConditionsCache.hh"

// Other external include files.
//#include "boost/shared_ptr.hpp"

namespace mu2e {

  class ConditionsService2 {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> fileName{Name("fileName"),
	  Comment("SimpleConfig file name")};
      fhicl::Atom<int> verbose{Name("verbose"),
	  Comment("verbosity 0 or 1"),0};
    };

    // this line is required by art to allow the command line help print
    typedef art::ServiceTable<Config> Parameters;

    ConditionsService2(Parameters const& config, 
		       art::ActivityRegistry& iRegistry);

    ConditionsCache::ptr getCache(std::string name) { 
      return _caches[name];
    }
    //void postBeginJob();

  private:

    // This is not copyable or assignable - private and unimplemented.
    ConditionsService2 const& operator=(ConditionsService2 const& rhs);
    ConditionsService2(ConditionsService2 const& rhs);

    Config _config;
    SimpleConfig _simpleConfig;
    int _verbose;
    std::map<std::string,ConditionsCache::ptr> _caches;
  };

}

DECLARE_ART_SERVICE(mu2e::ConditionsService2, LEGACY)
#endif /* ConditionsService2_ConditionsService2_hh */
