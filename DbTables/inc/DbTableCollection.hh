#ifndef DataProducts_DbTableCollection_hh
#define DataProducts_DbTableCollection_hh

#include <iostream>
#include <iomanip>
#include "DbTables/inc/DbLiveTable.hh"

namespace mu2e {
    typedef std::vector<mu2e::DbLiveTable> DbTableCollection;
}

#endif
