#ifndef RecoDataProducts_CaloProtoCluster_hh
#define RecoDataProducts_CaloProtoCluster_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloHit.hh"
#include <vector>

namespace mu2e {


   class CaloProtoCluster {
       
       public:
           typedef art::Ptr<CaloHit>       CaloHitPtr;
           typedef std::vector<CaloHitPtr> CaloHitPtrVector;

	   CaloProtoCluster() : 
	      _time(0.),_timeErr(0.0),_energyDep(0.),_energyDepErr(0.0),_caloHitsPtrVector(),_isSplit(false)
	   {}

	   CaloProtoCluster(double time, double timeErr, double energy, double energyErr, 
	                    CaloHitPtrVector caloHit, bool isSplit) : 
	      _time(time),_timeErr(timeErr),_energyDep(energy),_energyDepErr(energyErr),
	      _caloHitsPtrVector(caloHit),_isSplit(isSplit)
	   {}


	         double              time()              const {return _time;}            
	         double              timeErr()           const {return _timeErr;}            
	         double              energyDep()         const {return _energyDep;}       
	         double              energyDepErr()      const {return _energyDepErr;}       
	   const CaloHitPtrVector&   caloHitsPtrVector() const {return _caloHitsPtrVector;}
	         bool                isSplit()           const {return _isSplit;} 



        private:
	   double            _time;       
	   double            _timeErr;       
	   double            _energyDep;  
	   double            _energyDepErr;  
	   CaloHitPtrVector  _caloHitsPtrVector;
	   bool              _isSplit;    

   };
   
   
   typedef std::vector<mu2e::CaloProtoCluster> CaloProtoClusterCollection;
} 

#endif 
