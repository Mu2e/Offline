#ifndef MCDataProducts_CrvDigiCollection_hh
#define MCDataProducts_CrvDigiCollection_hh

//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "RecoDataProducts/inc/CrvDigi.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvDigi> CrvDigiCollection;
}

#endif /* RecoDataProducts_CrvDigiCollection_hh */
