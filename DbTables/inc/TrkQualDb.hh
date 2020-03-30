#ifndef DbTables_TrkQualDb_hh
#define DbTables_TrkQualDb_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "DbTables/inc/MVAToolDb.hh"

namespace mu2e {
  
  class TrkQualDb : public MVAToolDb {
  public:
    typedef std::shared_ptr<TrkQualDb> ptr_t;
    typedef std::shared_ptr<const TrkQualDb> cptr_t;

    TrkQualDb():MVAToolDb("TrkQualDb") {
    }
  };
};
#endif
