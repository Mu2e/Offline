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

    unsigned short nCaloHits    () const { return nCaloHitsD0_ + nCaloHitsD1_; }
    unsigned short nCaloHitsD0  () const { return nCaloHitsD0_ ; }
    unsigned short nCaloHitsD1  () const { return nCaloHitsD1_ ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }
    size_t         nCaphriHits  () const { return caphriHits_.size(); }

    // Methods to store/retrieve caphri hit info
    bool addCaphriHit(const double energy, const int ID);
    void getCaphriHit(const size_t ihit, double& energy, int& ID) const;

  private:
    // unsigned short  nCaloHits_    = 0;
    unsigned short  nCaloHitsD0_  = 0;
    unsigned short  nCaloHitsD1_  = 0;
    unsigned short  caloEnergy_   = 0;
    std::vector<unsigned short>  caphriHits_ = {};

    // For storing compact CAPHRI hit info
    constexpr static double caphriEnergyUnits_       = 0.01; // Store CAPHRI hit energies in units of 0.01 MeV
    constexpr static int    caphriIndexBits_         =   14; // Store the CAPHRI hit index 14 bits into the 16-bit word
  };

  typedef std::vector<mu2e::IntensityInfoCalo> IntensityInfosCalo;
}

#endif
