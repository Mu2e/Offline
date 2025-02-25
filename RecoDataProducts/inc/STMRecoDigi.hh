#ifndef RecoDataProducts_STMRecoDigi_hh
#define RecoDataProducts_STMRecoDigi_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include <vector>


namespace mu2e
{
   class STMRecoDigi
   {
      public:
        STMRecoDigi():
           stmDigi_(),eDep_(0),eDepErr_(0),time_(0),timeErr_(0),chi2_(-1),ndf_(0),pileUp_(false)
        {}

        STMRecoDigi(art::Ptr<STMDigi> stmDigi, float eDep, float eDepErr, float time,
                        float timeErr, float chi2, unsigned ndf, bool pileUp) :
           stmDigi_(stmDigi),eDep_(eDep),eDepErr_(eDepErr),time_(time),timeErr_(timeErr),
           chi2_(chi2),ndf_(ndf),pileUp_(pileUp)
        {}

        const art::Ptr<STMDigi>&  stmDigiPtr()   const {return stmDigi_;}
        int                        DetID()        const {return stmDigi_->DetID();}
        float                      energyDep()     const {return eDep_;}
        float                      energyDepErr()  const {return eDepErr_;}
        float                      time()          const {return time_;}
        float                      timeErr()       const {return timeErr_;}
        float                      chi2()          const {return chi2_;}
        unsigned                   ndf()           const {return ndf_;}
        bool                       pileUp()        const {return pileUp_;}

      private:
        art::Ptr<STMDigi>  stmDigi_;
        float               eDep_;
        float               eDepErr_;
        float               time_;
        float               timeErr_;
        float               chi2_;
        unsigned            ndf_;
        bool                pileUp_;
   };

   using STMRecoDigiCollection = std::vector<mu2e::STMRecoDigi>;
}
#endif
