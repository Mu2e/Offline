// Implementation of methods
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include <algorithm>

namespace mu2e {
  // Static methods to encode CAPHRI hit information
  unsigned short IntensityInfoCalo::encodeCaphriIndex(const int id) {
    const auto itr = std::find(CaloConst::_caphriId.begin(), CaloConst::_caphriId.end(), id);
    if(itr == CaloConst::_caphriId.end()) return -1; // not a CAPHRI crystal
    const unsigned short idx = std::distance(CaloConst::_caphriId.begin(), itr);
    return idx;
  }

  unsigned short IntensityInfoCalo::encodeCaphriEnergy(const double energy) {
    const unsigned short e_short = energy / CaloConst::_caphriEnergyUnits; // store in different units from MeV input
    return e_short;
  }

  unsigned short IntensityInfoCalo::encodeCaphriHit(const unsigned short energy, const unsigned short index) {
    return energy | index << CaloConst::_caphriIndexBits;
  }

  // Static methods to decode CAPHRI hit information
  int IntensityInfoCalo::decodeCaphriIndex(const unsigned short idx) {
    if(idx > 3) return -1;
    return CaloConst::_caphriId[idx];
  }

  double IntensityInfoCalo::decodeCaphriEnergy(const unsigned short e_short) {
    const double energy = e_short * CaloConst::_caphriEnergyUnits; // convert from stored units
    return energy;
  }

  void IntensityInfoCalo::decodeCaphriHit(const unsigned short encoded, unsigned short& energy, unsigned short& index) {
    energy = encoded & ~CaloConst::_caphriIndexMask;
    index = (encoded &  CaloConst::_caphriIndexMask) >> CaloConst::_caphriIndexBits;
  }
}
