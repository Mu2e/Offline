#ifndef DbTables_DbTableFactory_hh
#define DbTables_DbTableFactory_hh

#include "DbTables/inc/DbTable.hh"

namespace mu2e  {

  class DbTableFactory {
  public:
    static mu2e::DbTable::ptr_t newTable(std::string const& name);
  };

}

#endif
