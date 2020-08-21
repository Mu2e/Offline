#ifndef MCDataProducts_CrvPhotonsCollection_hh
#define MCDataProducts_CrvPhotonsCollection_hh

//
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvPhotons> CrvPhotonsCollection;
}

#endif /* MCDataProducts_CrvPhotonsCollection_hh */
