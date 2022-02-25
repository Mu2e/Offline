#ifndef DbService_DbValTool_hh
#define DbService_DbValTool_hh

// This class takes a DbValCache, the container of the purpose, version, etc
// database heriarchy, performs some basic tasks, including flattening
// it into a DbSet, a list of IoVs.  These tasks are also used in 
// DbTool, so the code should be here, in one place
// its only member is a ref to the cache, so it can be easily created/destroyed

#include <string>
#include <chrono>
#include <curl/curl.h>
#include "Offline/DbTables/inc/DbValCache.hh"
#include "Offline/DbTables/inc/DbVersion.hh"
#include "Offline/DbTables/inc/DbSet.hh"


namespace mu2e {
  class DbValTool {
  public:

    DbValTool(DbValCache const& cache):_valcache(cache) {}

    // convert input strings to their database indexes
    // allow empty strings and return -1
    int findPid(std::string purpose) const;
    void findPidVid(std::string purpose, std::string version, 
                    int& pid, int& vid) const;
    // convert between tid and name, bool tells if found
    bool tidByName(std::string const& name, int& tid) const;
    bool nameByTid(int tid, std::string& name) const;
    void printSet(DbSet const& dbset) const;

    // given a purpose/version, which defines and entire calibration set,
    // fill a DbSet with the lists of IoV for each table
    void fillSetVer(DbVersion const& dver, DbSet& dbset) const;
    // fill a DbSet based on a list of groups
    void fillSetGid(std::vector<int> const& gids, DbSet& dbset) const;

  private:

    DbValCache const& _valcache;

  };

}

#endif
