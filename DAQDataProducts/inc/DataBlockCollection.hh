#ifndef DAQDataProducts_DataBlockCollection_hh
#define DAQDataProducts_DataBlockCollection_hh

//
// Define a type for a collection of DataBlock objects.
//
//
// Original author Tomonari Miyashita
//

#include <vector>

#include "DAQDataProducts/inc/DataBlock.hh"

namespace mu2e {
   typedef std::vector<mu2e::DataBlock> DataBlockCollection;
}

#endif /* DAQDataProducts_DataBlockCollection_hh */
