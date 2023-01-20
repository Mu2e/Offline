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
  typedef std::map<STMChannel, float> SamplingFrequencyMap;

  STMEnergyCalib(CalibMap const& cmap, PedestalMap const& pmap, SamplingFrequencyMap const& sfmap) :
    ProditionsEntity(cxname), _cmap(cmap), _pmap(pmap), _sfmap(sfmap) {}

  // the calibration struct for a given channel
  STMEnergyCorr const& calib(STMChannel const& channel) const {
    return _cmap.at(channel);
  }

  float pedestal(STMChannel const& channel) const { return _pmap.at(channel); }
  float samplingFrequency(STMChannel const& channel) const { return _sfmap.at(channel); }

  void setCalib(CalibMap cmap) { _cmap = cmap; }
  void setPedestal(PedestalMap pmap) { _pmap = pmap; }

 private:
  CalibMap _cmap;
  PedestalMap _pmap;
  SamplingFrequencyMap _sfmap;
};

}  // namespace mu2e

#endif
