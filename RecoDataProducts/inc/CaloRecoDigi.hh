#ifndef RecoDataProducts_CaloRecoDigi_hh
#define RecoDataProducts_CaloRecoDigi_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include <vector>


namespace mu2e
{
   class CaloRecoDigi
   {
      public:
        CaloRecoDigi():
           caloDigi_(),eDep_(0),eDepErr_(0),time_(0),timeErr_(0),chi2_(-1),ndf_(0),pileUp_(false)
        {}

        CaloRecoDigi(art::Ptr<CaloDigi> caloDigi, float eDep, float eDepErr, float time,
                        float timeErr, float chi2, unsigned ndf, bool pileUp) :
           caloDigi_(caloDigi),eDep_(eDep),eDepErr_(eDepErr),time_(time),timeErr_(timeErr),
           chi2_(chi2),ndf_(ndf),pileUp_(pileUp)
        {}

        const art::Ptr<CaloDigi>&  caloDigiPtr()   const {return caloDigi_;}
        int                        SiPMID()        const {return caloDigi_->SiPMID();}
        float                      energyDep()     const {return eDep_;}
        float                      energyDepErr()  const {return eDepErr_;}
        float                      time()          const {return time_;}
        float                      timeErr()       const {return timeErr_;}
        float                      chi2()          const {return chi2_;}
        unsigned                   ndf()           const {return ndf_;}
        bool                       pileUp()        const {return pileUp_;}

      private:
        art::Ptr<CaloDigi>  caloDigi_;
        float               eDep_;
        float               eDepErr_;
        float               time_;
        float               timeErr_;
        float               chi2_;
        unsigned            ndf_;
        bool                pileUp_;
   };

   using CaloRecoDigiCollection = std::vector<mu2e::CaloRecoDigi>;
}
#endif
