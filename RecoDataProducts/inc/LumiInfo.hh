//
// Class to collect the info needed for monitoring POT / stop muons
// bvitali May 2021
//


#ifndef RecoDataProducts_LumiInfo_hh
#define RecoDataProducts_LumiInfo_hh

namespace mu2e {

    class LumiInfo
    {
        public:
            LumiInfo(unsigned short nTrackerHits, unsigned short nCaloHits,
                    unsigned short nProtonTCs, unsigned short caloEnergy,
                    unsigned short testVariable):
                nTrackerHits_(nTrackerHits), nCaloHits_(nCaloHits),
                nProtonTCs_(nProtonTCs), caloEnergy_(caloEnergy),
                testVariable_(testVariable)
            {}

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
