// Implementation of methods
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "cetlib_except/exception.h"
#include <algorithm>

namespace mu2e {

  // Store a Caphri hit
  void IntensityInfoCalo::addCaphriHit(const double energy, const int ID) {

    // Map the crystal ID to an index between 0-3
    const auto itr = std::find(CaloConst::_caphriId.begin(), CaloConst::_caphriId.end(), ID);
    if(itr == CaloConst::_caphriId.end()) throw cet::exception("RECO") << "Crystal ID " << ID << " is not listed as a CAPHRI crystal!"; // not a CAPHRI crystal
    const unsigned short index = std::distance(CaloConst::_caphriId.begin(), itr);
    if(index > 3) throw cet::exception("RECO") << "CAPHRI crystal index " << index << "(crystal ID " << ID << ") is greater than 3!"; // must be between 0-3

    const unsigned short e_short = energy / caphriEnergyUnits_; // store in different units from MeV input to adjust the precision

    // encode the energy + index into a single short
    const unsigned short hit_info = e_short | index << caphriIndexBits_;

    // add it to the hit vector
    caphriHits_.push_back(hit_info);
  }

  // Decode CAPHRI hit information
  void IntensityInfoCalo::getCaphriHit(const size_t ihit, double& energy, int& ID) const {
    ID = -1; energy = 0.;
    if(ihit >= caphriHits_.size()) throw cet::exception("RECO") << "Hit index " << ihit << " is >= N(hits) = " << caphriHits_.size();

    const auto hit_info = caphriHits_[ihit];

    const unsigned short e_short =  hit_info & ~caphriIndexMask_;
    const unsigned short index   = (hit_info &  caphriIndexMask_) >> caphriIndexBits_;
    if(index >= CaloConst::_caphriId.size()) throw cet::exception("RECO") << "CAPHRI crystal index " << index << " is >= N(CAPHRI) = " << CaloConst::_caphriId.size();

    energy = e_short * caphriEnergyUnits_; // convert from stored units
    ID = CaloConst::_caphriId[index];
  }
}
