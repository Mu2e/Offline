#ifndef RecoDataProducts_CrvCoincidence_hh
#define RecoDataProducts_CrvCoincidence_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e 
{
  class CrvCoincidence
  {
    public:

    CrvCoincidence() {}

    CrvCoincidence(const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses, int crvSectorType) : 
                   _crvRecoPulses(crvRecoPulses), _crvSectorType(crvSectorType) {}

    const std::vector<art::Ptr<CrvRecoPulse> > &GetCrvRecoPulses() const {return _crvRecoPulses;}
    int                                         GetCrvSectorType() const {return _crvSectorType;}

    private:

    std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
    int                                  _crvSectorType;
  };
}

#endif /* RecoDataProducts_CrvCoincidence_hh */
