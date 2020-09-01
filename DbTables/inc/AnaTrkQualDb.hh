#ifndef DbTables_AnaTrkQualDb_hh
#define DbTables_AnaTrkQualDb_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "DbTables/inc/MVAToolDb.hh"

namespace mu2e {
  
  class AnaTrkQualDb : public MVAToolDb {
  public:
    typedef std::shared_ptr<AnaTrkQualDb> ptr_t;
    typedef std::shared_ptr<const AnaTrkQualDb> cptr_t;

    AnaTrkQualDb():MVAToolDb("AnaTrkQualDb", "ana.trkqualdb") {
    }
  };
};
#endif
