#ifndef RecoDataProducts_StrawDigiCollection_hh
#define RecoDataProducts_StrawDigiCollection_hh

//
// Define a type for a collection of StrawDigi objects.
//
// $Id: StrawDigiCollection.hh,v 1.1 2013/12/07 19:50:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:50:42 $
//
// Original author David Brown 
//

#include <vector>
//#include <array>

#include "RecoDataProducts/inc/StrawDigi.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawDigi> StrawDigiCollection;
}

#endif /* DataProducts_StrawDigiCollection_hh */
