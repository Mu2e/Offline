#ifndef MCDataProducts_CRVHitCollection_hh
#define MCDataProducts_CRVHitCollection_hh

//
// Define a type for a collection of CRVHit objects.
// They are stored in time-ordered sets.
// Each set of CRVHits is associated with a CRV bar index (identifying a particular CRV counter).
//
// $Id: CRVHitCollection.hh,v 1.1 2014/08/07 01:33:41 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CRVHit.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,std::set<mu2e::CRVHit> > CRVHitCollection;
}

#endif /* MCDataProducts_CRVHitCollection_hh */
