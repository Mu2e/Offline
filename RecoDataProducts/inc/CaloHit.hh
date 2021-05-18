#ifndef RecoDataProducts_CaloHit_hh
#define RecoDataProducts_CaloHit_hh

#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "canvas/Persistency/Common/Ptr.h"

#include <vector>
#include <map>

namespace mu2e {

   class CaloHit
   {
       public:
          CaloHit(): crystalId_(-1),nSiPMs_(0),time_(0.),eDep_(0.),eDepErr_(0.),timeErr_(0.),recoCaloDigis_()
          {}

          CaloHit(int crystalId, int nSiPMs, float time, float timeErr, float eDep, float eDepErr, 
                  std::vector<art::Ptr<CaloRecoDigi>> &CaloRecoDigi)  :
                    crystalId_(crystalId),nSiPMs_(nSiPMs),time_(time),eDep_(eDep),
                    eDepErr_(eDepErr),timeErr_(timeErr),recoCaloDigis_(CaloRecoDigi)
          {}

          CaloHit(int crystalId, int nSiPMs, float time, float eDep)  :
	    crystalId_(crystalId),nSiPMs_(nSiPMs),time_(time),eDep_(eDep),
	    eDepErr_(0),timeErr_(0)
          {}

          int                                          crystalID      () const { return crystalId_; }	  
          int                                          nSiPMs         () const { return nSiPMs_;}	  
          float                                        time           () const { return time_;}		  
          float                                        timeErr        () const { return timeErr_;}	  
          float                                        energyDep      () const { return eDep_;} 	  
          float                                        energyDepErr   () const { return eDepErr_;} 	  
          float                                        energyDepTot   () const { return eDep_*nSiPMs_;}	  
          float                                        energyDepTotErr() const { return eDepErr_*nSiPMs_;}
          const std::vector<art::Ptr<CaloRecoDigi>>&   recoCaloDigis  () const { return recoCaloDigis_;}

          void     setCrystalID      (int  ID) { crystalId_ = ID; }	  
          void     setNSiPMs         (int   N) { nSiPMs_    = N;  }	  
          void     setTime           (float T) { time_      = T;  }		  
          void     setEDep           (float E) { eDep_      = E;  }  	  

        private:
          int    crystalId_;
          int    nSiPMs_;
          float  time_;             
          float  eDep_;        
          float  eDepErr_;        
          float  timeErr_;             
          std::vector<art::Ptr<CaloRecoDigi>> recoCaloDigis_;

   };
   
   using CaloHitPtr        = art::Ptr<CaloHit>;
   using CaloHitPtrVector  = std::vector<CaloHitPtr>;
   using CaloHitCollection = std::vector<mu2e::CaloHit> ;
   using CaloHitRemapping  = std::map<art::Ptr<CaloHit>,art::Ptr<CaloHit>>;
}

#endif
