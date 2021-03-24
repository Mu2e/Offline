#ifndef DbService_DbIdList_hh
#define DbService_DbIdList_hh

//
// Read the connections.txt file, usually from cvmfs,
// and hold the connections information.  
// connections are read from cvmfs in order to steer 
// running jobs in approximately real time.
// As needed, create new DbId's with the appropriate 
// connections info, for a given database name.
//

#include <string>
#include <vector>
#include "DbTables/inc/DbId.hh"

namespace mu2e {

  class DbIdList {
  public:
    DbIdList();
    // return the appropriate DbId given a database name
    DbId getDbId(const std::string& name="mu2e_conditions_prd");
  private:
    std::vector<DbId> _ids;
  };
}
#endif
