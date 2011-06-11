#ifndef Sandbox_TracerProductCollection_hh
#define Sandbox_TracerProductCollection_hh
//
// Access to a bunch of TracerProdcuts created on the heap.
// 
// $Id: TracerProductCollection.hh,v 1.2 2011/06/11 02:27:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 02:27:32 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/OwningPointerCollection.hh"
#include "Sandbox/inc/TracerProduct.hh"

namespace mu2e {

  typedef OwningPointerCollection<TracerProduct> TracerProductCollection;

} // namespace mu2e

#endif /* Sandbox_TracerProductCollection_hh */
