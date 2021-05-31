#ifndef DbService_DbService_hh
#define DbService_DbService_hh

#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Services/Registry/ServiceTable.h" 
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "cetlib_except/exception.h"
#include "DbService/inc/DbEngine.hh"



namespace mu2e {

  class DbService {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> purpose{Name("purpose"),
	  Comment("defines what tables to load"),"PRODUCTION"};
      fhicl::Atom<std::string> version{Name("version"),
	  Comment("the version of the purpose"),"v*_*_*"};
      fhicl::Atom<std::string> dbName{Name("dbName"),
	  Comment("which database to use"),"none"};
      fhicl::OptionalSequence<std::string> textFile{Name("textFile"),
	  Comment("list of text files containing override table data")};
      fhicl::Atom<int> verbose{Name("verbose"), 
	  Comment("verbose flag, 0 to 10"),0};
      fhicl::Atom<bool> saveCsv{Name("saveCsv"), 
	  Comment("save csv content in tables, default false"),false};
      fhicl::OptionalAtom<bool> fastStart{Name("fastStart"), 
	  Comment("read the DB immedatiately, not on first use")};
      fhicl::OptionalAtom<int> cacheLifetime{Name("cacheLifetime"), 
	  Comment("if >0, read IoV from cache, but renew each lifetime s")};
    };

    // this line is required by art to allow the command line help print
    typedef art::ServiceTable<Config> Parameters;

    explicit DbService(Parameters const& config, 
		       art::ActivityRegistry& iRegistry);
    ~DbService();

    // Functions registered for callbacks.
    void postBeginJob();
    void postEndJob();

    // how the DbHandle interacts with this service
    DbEngine& engine() {return _engine;}


  private:

    // This is not copyable or assignable - private and unimplemented.
    DbService const& operator=(DbService const& rhs);
    DbService(DbService const& rhs);

    Config _config;
    int _verbose;

    DbVersion _version;
    DbEngine _engine;

  };

}


DECLARE_ART_SERVICE(mu2e::DbService, SHARED)
#endif
