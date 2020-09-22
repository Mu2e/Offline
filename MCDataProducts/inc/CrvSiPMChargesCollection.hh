#ifndef MCDataProducts_CrvSiPMChargesCollection_hh
#define MCDataProducts_CrvSiPMChargesCollection_hh

//
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvSiPMCharges> CrvSiPMChargesCollection;
}

#endif /* MCDataProducts_CrvSiPMChargesCollection_hh */
