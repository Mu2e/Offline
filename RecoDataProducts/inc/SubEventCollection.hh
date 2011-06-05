#ifndef EarlyPatRec_SubEventCollection_hh
#define EarlyPatRec_SubEventCollection_hh

//
// Define a type for a collection of SubEvent objects.
//
// $Id: SubEventCollection.hh,v 1.1 2011/06/05 23:11:35 mf Exp $
// $Author: mf $
// $Date: 2011/06/05 23:11:35 $
//
// Original author Mark Fischler
//

#include <vector>

#include "RecoDataProducts/inc/SubEvent.hh"

namespace mu2e {
   typedef std::vector<mu2e::SubEvent> SubEventCollection;
}

#endif /* EarlyPatRec_SubEventCollection_hh */
