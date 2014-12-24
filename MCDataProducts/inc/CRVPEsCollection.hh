#ifndef MCDataProducts_CRVPEsCollection_hh
#define MCDataProducts_CRVPEsCollection_hh

//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CRVPEs.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CRVPEs> CRVPEsCollection;
}

#endif /* MCDataProducts_CRVPEsCollection_hh */
