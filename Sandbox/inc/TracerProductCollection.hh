#ifndef Sandbox_TracerProductCollection_hh
#define Sandbox_TracerProductCollection_hh
//
// Access to a bunch of TracerProdcuts created on the heap.
// 
// $Id: TracerProductCollection.hh,v 1.1 2011/06/11 01:49:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 01:49:10 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/BarePointerCollection.hh"
#include "Sandbox/inc/TracerProduct.hh"

namespace mu2e {

  typedef BarePointerCollection<TracerProduct> TracerProductCollection;

} // namespace mu2e

#endif /* Sandbox_TracerProductCollection_hh */
