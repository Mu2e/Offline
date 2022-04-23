//
// Class to collect the info needed for monitoring POT / stop muons
// bvitali May 2021
//


#ifndef RecoDataProducts_IntensityInfo_hh
#define RecoDataProducts_IntensityInfo_hh

namespace mu2e {

  class IntensityInfo
  {
  public:
    IntensityInfo() {}
    IntensityInfo(unsigned short nTrackerHits, unsigned short nCaloHits,
		  unsigned short nProtonTCs, unsigned short caloEnergy,
		  unsigned short testVariable):
      nTrackerHits_(nTrackerHits), nCaloHits_(nCaloHits),
      nProtonTCs_(nProtonTCs), caloEnergy_(caloEnergy),
      testVariable_(testVariable)
    {}

	   
    void setNTrackerHits   (unsigned short tmp) {nTrackerHits_ = tmp;}
    void setNCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
    void setNProtonTCs     (unsigned short tmp) {nProtonTCs_   = tmp;}
    void setCaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}
    void setTestVariable   (unsigned short tmp) {testVariable_ = tmp;} 

    unsigned short nTrackerHits () const { return nTrackerHits_; }
    unsigned short nCaloHits    () const { return nCaloHits_   ; }
    unsigned short nProtonTCs   () const { return nProtonTCs_  ; }
    unsigned short caloEnergy   () const { return caloEnergy_  ; }
    unsigned short testVariable () const { return testVariable_; }

  private:
    unsigned short  nTrackerHits_ = 0;
    unsigned short  nCaloHits_    = 0;            
    unsigned short  nProtonTCs_   = 0;        
    unsigned short  caloEnergy_   = 0;        
    unsigned short  testVariable_ = 0;            
  };
}

#endif
