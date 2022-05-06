//
// Class to collect the info needed for monitoring POT / stop muons
//
//


#ifndef RecoDataProducts_IntensityInfoCalo_hh
#define RecoDataProducts_IntensityInfoCalo_hh

namespace mu2e {

  class IntensityInfoCalo
  {
  public:
    IntensityInfoCalo() {}
    IntensityInfoCalo( unsigned short nCaloHits, unsigned short caloEnergy):
      nCaloHits_(nCaloHits),caloEnergy_(caloEnergy)
    {}


    void setNCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
    void setCaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}

    unsigned short nCaloHits    () const { return nCaloHits_   ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }

  private:
    unsigned short  nCaloHits_    = 0;
    unsigned short  caloEnergy_   = 0;
  };
}

#endif
