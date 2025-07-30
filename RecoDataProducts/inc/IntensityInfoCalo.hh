//
// Class to collect the info needed for monitoring POT / stop muons
//
//


#ifndef RecoDataProducts_IntensityInfoCalo_hh
#define RecoDataProducts_IntensityInfoCalo_hh
#include <vector>

namespace mu2e {

  class IntensityInfoCalo
  {
  public:
    IntensityInfoCalo() {}
    IntensityInfoCalo( unsigned short nCaloHits, unsigned short caloEnergy, std::vector<unsigned short> caphriHits):
      nCaloHits_(nCaloHits),caloEnergy_(caloEnergy),caphriHits_(caphriHits)
    {}


    void setNCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
    void setCaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}
    void setCaphriHits     (std::vector<unsigned short> tmp) {caphriHits_  = tmp;}

    unsigned short nCaloHits    () const { return nCaloHits_   ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }
    std::vector<unsigned short> caphriHits  () const { return caphriHits_ ; }

  private:
    unsigned short  nCaloHits_    = 0;
    unsigned short  caloEnergy_   = 0;
    std::vector<unsigned short>  caphriHits_ = {};
  };

  typedef std::vector<mu2e::IntensityInfoCalo> IntensityInfosCalo;
}

#endif
