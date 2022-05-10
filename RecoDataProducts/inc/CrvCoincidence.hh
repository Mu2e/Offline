#ifndef RecoDataProducts_CrvCoincidence_hh
#define RecoDataProducts_CrvCoincidence_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e
{
  class CrvCoincidence
  {
    public:

    CrvCoincidence() {}

    CrvCoincidence(const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses, int crvSectorType,
                   const std::vector<float> &slopes, const std::vector<int> &layers) :
                   _crvRecoPulses(crvRecoPulses), _crvSectorType(crvSectorType), _slopes(slopes), _layers(layers) {}

    const std::vector<art::Ptr<CrvRecoPulse> > &GetCrvRecoPulses() const {return _crvRecoPulses;}
    int                                         GetCrvSectorType() const {return _crvSectorType;}
    const std::vector<float>                          &GetSlopes() const {return _slopes;}
    const std::vector<int>                            &GetLayers() const {return _layers;}

    private:

    std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
    int                                  _crvSectorType;
    std::vector<float>                   _slopes; //width direction of counter / thickness direction of counter //slope=0 means straight through the module
    std::vector<int>                     _layers;
  };
  typedef std::vector<mu2e::CrvCoincidence> CrvCoincidenceCollection;
}

#endif /* RecoDataProducts_CrvCoincidence_hh */
