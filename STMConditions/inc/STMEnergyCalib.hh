#ifndef STMConditions_STMEnergyCalib_hh
#define STMConditions_STMEnergyCalib_hh

//
// Prodition to hold parameters to convert STM digis to energy
//

#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/STMConditions/inc/STMEnergyCorr.hh"
#include <map>

namespace mu2e {

class STMEnergyCalib : virtual public ProditionsEntity {
 public:
  typedef std::shared_ptr<STMEnergyCalib> ptr_t;
  typedef std::shared_ptr<const STMEnergyCalib> cptr_t;
  constexpr static const char* cxname = {"STMEnergyCalib"};

  typedef std::map<STMChannel, STMEnergyCorr> CalibMap;
  typedef std::map<STMChannel, float> PedestalMap;

  STMEnergyCalib(CalibMap const& cmap, PedestalMap const& pmap) :
      ProditionsEntity(cxname), _cmap(cmap), _pmap(pmap) {}

  // the calibration struct for a given channel
  STMEnergyCorr const& calib(STMChannel const& channel) const {
    return _cmap.at(channel);
  }

  float pedestal(STMChannel const& channel) const { return _pmap.at(channel); }

 private:
  CalibMap _cmap;
  PedestalMap _pmap;
};

}  // namespace mu2e

#endif
