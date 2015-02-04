#ifndef MCDataProducts_CrvWaveformsCollection_hh
#define MCDataProducts_CrvWaveformsCollection_hh

//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/CrvWaveforms.hh"

#include <map>
#include <set>

namespace mu2e 
{
  typedef std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvWaveforms> CrvWaveformsCollection;
}

#endif /* MCDataProducts_CrvWaveformsCollection_hh */
