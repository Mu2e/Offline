#ifndef MCDataProducts_CrvSiPMResponsesCollection_hh
#define MCDataProducts_CrvSiPMResponsesCollection_hh

//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CrvSiPMResponses.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvSiPMResponses> CrvSiPMResponsesCollection;
}

#endif /* MCDataProducts_CrvSiPMResponsesCollection_hh */
