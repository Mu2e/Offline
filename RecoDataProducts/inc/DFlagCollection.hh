#ifndef RecoDataProducts_DFlagCollection_hh
#define RecoDataProducts_DFlagCollection_hh

//
// Define a type for a collection of DFlag objects.
//
// $Id: DFlagCollection.hh,v 1.1 2013/01/26 18:18:44 brownd Exp $
// $Author: brownd $
// $Date: 2013/01/26 18:18:44 $
//
// Original author David Brown
//

#include <vector>

#include "RecoDataProducts/inc/DFlag.hh"

namespace mu2e {
   typedef std::vector<mu2e::DFlag> DFlagCollection;
}

#endif /* RecoDataProducts_DFlagCollection_hh */
