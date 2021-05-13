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
            LumiInfo(): nTrackerHits_(0), nCaloHits_(0),
                        nProtonTCs_(0), caloEnergy_(0),
                        testVariable_(0)
            {}

            LumiInfo(unsigned short nTrackerHits, unsigned short nCaloHits,
                    unsigned short nProtonTCs, unsigned short caloEnergy,
                    unsigned short testVariable):
                nTrackerHits_(nTrackerHits), nCaloHits_(nCaloHits),
                nProtonTCs_(nProtonTCs), caloEnergy_(caloEnergy),
                testVariable_(testVariable)
            {}

            void setnTrackerHits   (unsigned short tmp) {nTrackerHits_ = tmp;}
            void setnCaloHits      (unsigned short tmp) {nCaloHits_    = tmp;}
            void setnProtonTCs     (unsigned short tmp) {nProtonTCs_   = tmp;}
            void setcaloEnergy     (unsigned short tmp) {caloEnergy_   = tmp;}
            void settestVariable   (unsigned short tmp) {testVariable_ = tmp;}

            unsigned short nTrackerHits () const { return nTrackerHits_; }
            unsigned short nCaloHits    () const { return nCaloHits_;    }
            unsigned short nProtonTCs   () const { return nProtonTCs_;   }
            unsigned short caloEnergy   () const { return caloEnergy_;   }
            unsigned short testVariable () const { return testVariable_; }

        private:
            unsigned short  nTrackerHits_;
            unsigned short  nCaloHits_;            
            unsigned short  nProtonTCs_;        
            unsigned short  caloEnergy_;        
            unsigned short  testVariable_;            
    };
}

#endif
