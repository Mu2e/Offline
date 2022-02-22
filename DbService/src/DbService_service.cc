
#include "Offline/DbService/inc/DbService.hh"
#include "Offline/DbService/inc/DbIdList.hh"
#include "Offline/DbTables/inc/DbUtil.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace mu2e {

  DbService::DbService(Parameters const& config,
		       art::ActivityRegistry& iRegistry):
    _config(config()),_verbose(config().verbose()),
    _version(config().purpose(),config().version()){

    // register callbacks
    iRegistry.sPostBeginJob.watch(this, &DbService::postBeginJob);
    iRegistry.sPostEndJob.watch (this, &DbService::postEndJob );

    if(_verbose>0) {
      std::cout << "DbService  " 
		<< config().purpose() << " " 
		<< config().version() << "   "
		<< config().dbName() << std::endl;
      std::vector<std::string> files;
      config().textFile(files);
      std::cout << "DbService  textFiles: " ;
      for(auto const& s : files) std::cout << " " << s;
      std::cout << std::endl;
    }



    if(_verbose>2) std::cout << "DbService::constructor " <<std::endl;

    // the engine which will read db, hold calibrations, deliver them
    _engine.setVerbose(_verbose);
    _engine.setSaveCsv(_config.saveCsv());
    if(_config.nearestMatch()) {
      _engine.setNearestMatch(true);
      mf::LogWarning warn("DbService");
      warn<<"WARNING: setting DbService nearestMatch true will probably result "
          <<"in approximate or unreproducible results\n";
    }

    DbIdList idList; // read file of db connection details
    _engine.setDbId( idList.getDbId(_config.dbName()) );
    _engine.setVersion( _version );

    // if there were text files containing calibrations,
    // then read them and tell the engine to let them override IOV
    std::vector<std::string> files;
    _config.textFile(files);
   ConfigFileLookupPolicy configFile;
   std::string fn;
    for(auto ss : files ) {
      if(_verbose>2) std::cout << "DbService::beginJob reading file "<<
		       ss <<std::endl;
      fn = configFile(ss);
      auto coll = DbUtil::readFile(fn,_config.saveCsv());
      if(_verbose>2) {
	for(auto const& lt : coll) {
	  std::cout << "  read table " << lt.table().name() <<std::endl;
	}
      }
      _engine.addOverride(coll);
    } // end loop over files

    int cacheLifetime = 0;
    if(_config.cacheLifetime(cacheLifetime)) {
      _engine.reader().setCacheLifetime(cacheLifetime);
    }

    int retryTimeout = 0;
    if(_config.retryTimeout(retryTimeout)) {
      _engine.reader().setTimeout(retryTimeout);
    }


    // service will start calling the database at the first event,
    // so the service can exist without the DB being contacted.  
    // fastStart overrides this and starts reading the DB imediately.
    bool fastStart = false;
    if(_config.fastStart(fastStart)) {
      if (fastStart) {
        // this prepares IOV infrastructure for efficient queries
        if(_verbose>2) std::cout << "DbService::beginJob initializing engine "
                                 <<std::endl;
        _engine.beginJob();
      }
    }

  }

  DbService::~DbService() {}

  /********************************************************/
  void DbService::postBeginJob(){

  }


  /********************************************************/
  void DbService::postEndJob(){
    // just print summaries according to verbosity
    _engine.endJob();
  }

}

DEFINE_ART_SERVICE(mu2e::DbService)

