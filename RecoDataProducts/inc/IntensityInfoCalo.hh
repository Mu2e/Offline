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
    IntensityInfoCalo( unsigned short nCaloHits, unsigned short caloEnergy, unsigned short nCaphriHits):
      nCaloHits_(nCaloHits),caloEnergy_(caloEnergy),nCaphriHits_(nCaphriHits)
    {}


    void setNCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
    void setCaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}
    void setNCaphriHits    (unsigned short tmp) {nCaphriHits_  = tmp;}

    unsigned short nCaloHits    () const { return nCaloHits_   ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }
    unsigned short nCaphriHits  () const { return nCaphriHits_ ; }

  private:
    unsigned short  nCaloHits_    = 0;
    unsigned short  caloEnergy_   = 0;
    unsigned short  nCaphriHits_  = 0;
  };

  typedef std::vector<mu2e::IntensityInfoCalo> IntensityInfosCalo;
}

#endif
