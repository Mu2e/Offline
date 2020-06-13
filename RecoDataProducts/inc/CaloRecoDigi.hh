#ifndef RecoDataProducts_CaloRecoDigi_hh
#define RecoDataProducts_CaloRecoDigi_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include <vector>


namespace mu2e
{
  class CaloRecoDigi 
  {
     public:
       CaloRecoDigi(): 
          caloDigi_(),eDep_(0),eDepErr_(0),time_(0),timeErr_(0),chi2_(-1),ndf_(0),pileUp_(false)
       {} 

       CaloRecoDigi(art::Ptr<CaloDigi> caloDigi, double eDep, double eDepErr, double time, 
                       double timeErr, double chi2, int ndf, bool pileUp) :
	  caloDigi_(caloDigi),eDep_(eDep),eDepErr_(eDepErr),time_(time),timeErr_(timeErr),
	  chi2_(chi2),ndf_(ndf),pileUp_(pileUp)
       {}

       const art::Ptr<CaloDigi>&  caloDigiPtr()   const {return caloDigi_;}
             int                  ROid()          const {return caloDigi_->roId();} 
             double               energyDep()     const {return eDep_;} 
             double               energyDepErr()  const {return eDepErr_;} 
             double               time()          const {return time_;} 
             double               timeErr()       const {return timeErr_;} 
             double               chi2()          const {return chi2_;}   
             int                  ndf()           const {return ndf_;}   
             bool                 pileUp()        const {return pileUp_;}  


     private:
       art::Ptr<CaloDigi>  caloDigi_;   
       double              eDep_; 
       double              eDepErr_; 
       double              time_; 
       double              timeErr_; 
       double              chi2_; 
       int                 ndf_; 
       bool                pileUp_;
  };


  typedef std::vector<mu2e::CaloRecoDigi> CaloRecoDigiCollection;
}
#endif 
