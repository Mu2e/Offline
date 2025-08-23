// Implementation of methods
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include <algorithm>

namespace mu2e {

  // Store a Caphri hit
  bool IntensityInfoCalo::addCaphriHit(const double energy, const int ID) {

    // Map the crystal ID to an index between 0-3
    const auto itr = std::find(CaloConst::_caphriId.begin(), CaloConst::_caphriId.end(), ID);
    if(itr == CaloConst::_caphriId.end()) return false; // not a CAPHRI crystal
    const unsigned short index = std::distance(CaloConst::_caphriId.begin(), itr);
    if(index > 3) return false; // must be between 0-3

    const unsigned short e_short = energy / CaloConst::_caphriEnergyUnits; // store in different units from MeV input

    // encode the energy + index into a single short
    const unsigned short hit_info = e_short | index << CaloConst::_caphriIndexBits;

    // add it to the hit vector
    caphriHits_.push_back(hit_info);

    return true;
  }

  // Decode CAPHRI hit information
  void IntensityInfoCalo::getCaphriHit(const size_t ihit, double& energy, int& ID) const {
    ID = -1; energy = 0.;
    if(ihit >= caphriHits_.size()) return;

    const auto hit_info = caphriHits_[ihit];

    const unsigned short e_short =  hit_info & ~CaloConst::_caphriIndexMask;
    const unsigned short index   = (hit_info &  CaloConst::_caphriIndexMask) >> CaloConst::_caphriIndexBits;

    energy = e_short * CaloConst::_caphriEnergyUnits; // convert from stored units
    if(index < CaloConst::_caphriId.size()) ID = CaloConst::_caphriId[index];
  }
}
