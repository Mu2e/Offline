#ifndef CRVConditions_CRVCalib_hh
#define CRVConditions_CRVCalib_hh

//
// Prodition to hold calibration constants for CRV SiPMS
//

#include "Offline/CRVConditions/inc/CRVCalibPar.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include <cstdint>
#include <vector>

namespace mu2e {

class CRVCalib : virtual public ProditionsEntity {
 public:
  typedef std::shared_ptr<CRVCalib> ptr_t;
  typedef std::shared_ptr<const CRVCalib> cptr_t;
  constexpr static const char* cxname = {"CRVCalib"};

  typedef std::vector<CRVCalibPar> CalibVec;

  CRVCalib(const CalibVec& cvec) : ProditionsEntity(cxname), _cvec(cvec) {}

  const CRVCalibPar& calib(std::uint16_t channel) const {
    return _cvec.at(channel);
  }
  float pedestal(std::uint16_t channel) const {
    return _cvec.at(channel).pedestal();
  }
  float pulseHeight(std::uint16_t channel) const {
    return _cvec.at(channel).pulseHeight();
  }
  float pulseArea(std::uint16_t channel) const {
    return _cvec.at(channel).pulseArea();
  }
  float timeOffset(std::uint16_t channel) const {
    return _cvec.at(channel).timeOffset();
  }

 private:
  // a vector of CRVCalibPar indexed by channel
  CalibVec _cvec;
};

}  // namespace mu2e

#endif
