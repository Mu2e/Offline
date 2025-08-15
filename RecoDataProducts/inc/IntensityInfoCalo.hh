//
// Class to collect the info needed for monitoring POT / stop muons
//
//


#ifndef RecoDataProducts_IntensityInfoCalo_hh
#define RecoDataProducts_IntensityInfoCalo_hh

#include "Offline/DataProducts/inc/CaloConst.hh"
#include <vector>

namespace mu2e {

  class IntensityInfoCalo
  {
  public:
    IntensityInfoCalo() {}
    IntensityInfoCalo(/* unsigned short nCaloHits,*/ unsigned short caloEnergy, std::vector<unsigned short> caphriHits):
      // nCaloHits_(nCaloHits),
      caloEnergy_(caloEnergy),caphriHits_(caphriHits)
    {}


    // void setNCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
    void setNCaloHitsD0    (unsigned short tmp) {nCaloHitsD0_  = tmp;}
    void setNCaloHitsD1    (unsigned short tmp) {nCaloHitsD1_  = tmp;}
    void setCaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}
    void setCaphriHits     (std::vector<unsigned short> tmp) {caphriHits_  = tmp;}

    unsigned short nCaloHits    () const { return nCaloHitsD0_ + nCaloHitsD1_; }
    unsigned short nCaloHitsD0  () const { return nCaloHitsD0_ ; }
    unsigned short nCaloHitsD1  () const { return nCaloHitsD1_ ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }
    std::vector<unsigned short> caphriHits  () const { return caphriHits_ ; }
    size_t         nCaphriHits  () const { return caphriHits_.size(); }

    // Static methods to encode CAPHRI hit information
    static unsigned short encodeCaphriIndex(const int id) {
      const auto itr = std::find(CaloConst::_caphriId.begin(), CaloConst::_caphriId.end(), id);
      if(itr == CaloConst::_caphriId.end()) return -1; // not a CAPHRI crystal
      const unsigned short idx = std::distance(CaloConst::_caphriId.begin(), itr);
      return idx;
    }

    static unsigned short encodeCaphriEnergy(const double energy) {
      const unsigned short e_short = energy * 100.; // store in units of 0.01 MeV
      return e_short;
    }

    static unsigned short encodeCaphriHit(const unsigned short energy, const unsigned short index) {
      return energy | index << 14;
    }

    // Static methods to decode CAPHRI hit information
    static int decodeCaphriIndex(const unsigned short idx) {
      if(idx > 3) return -1;
      return CaloConst::_caphriId[idx];
    }

    static double decodeCaphriEnergy(const unsigned short e_short) {
      const double energy = e_short / 100.; // stored in units of 0.01 MeV
      return energy;
    }

    static void decodeCaphriHit(const unsigned short encoded, unsigned short& energy, unsigned short& index) {
      constexpr unsigned short mask = 0x3 << 14;
      energy = encoded & ~mask;
      index = (encoded & mask) >> 14;
    }

  private:
    // unsigned short  nCaloHits_    = 0;
    unsigned short  nCaloHitsD0_  = 0;
    unsigned short  nCaloHitsD1_  = 0;
    unsigned short  caloEnergy_   = 0;
    std::vector<unsigned short>  caphriHits_ = {};
  };

  typedef std::vector<mu2e::IntensityInfoCalo> IntensityInfosCalo;
}

#endif
