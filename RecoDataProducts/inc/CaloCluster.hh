#ifndef RecoDataProducts_CaloCluster_hh
#define RecoDataProducts_CaloCluster_hh
//
// Calorimeter cluster information
//
// Note: The size is set independently of the CaloHitPtrVector to work with fast clustering algorithm
//
#include "canvas/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloHit.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>
#include <ostream>


namespace mu2e {


   class CaloCluster {
       
       public:
            typedef art::Ptr<CaloHit>        CaloHitPtr;
            typedef std::vector<CaloHitPtr>  CaloHitPtrVector;

	    CaloCluster() : 
	       diskId_(-1),time_(0.),timeErr_(0.0),energyDep_(0.),energyDepErr_(0.),e1_(0.),e9_(0.),e25_(0.),
	       secondMoment_(0.),angle_(0),cog3Vector_(CLHEP::Hep3Vector(0,0,0)), caloHitsPtrVector_(),size_(0),isSplit_(false)
	    {}

	    CaloCluster(int iSection, double time, double timeErr, double energy, double energyErr, 
	                CaloHitPtrVector caloHits, unsigned size, bool isSplit) : 
	       diskId_(iSection),time_(time),timeErr_(timeErr),energyDep_(energy),energyDepErr_(energyErr),e1_(0),e9_(0),
	       e25_(0),secondMoment_(-1),angle_(0),cog3Vector_(CLHEP::Hep3Vector(0,0,0)),
	       caloHitsPtrVector_(caloHits),size_(size),isSplit_(isSplit)
	    {}


	          int                 diskId()            const{return diskId_;}       
	          int                 size()              const{return size_;}
	          double              time()              const{return time_;}            
	          double              timeErr()           const{return timeErr_;}            
	          double              energyDep()         const{return energyDep_;}       
	          double              energyDepErr()      const{return energyDepErr_;}       
	          double              e1()                const{return e1_;}       
	          double              e9()                const{return e9_;}       
	          double              e25()               const{return e25_;}       
	          double              secondMoment()      const{return secondMoment_;}       
	          double              angle()             const{return angle_;}       
	          bool                isSplit()           const{return isSplit_;} 
	    const CLHEP::Hep3Vector&  cog3Vector()        const{return cog3Vector_;}      
	    const CaloHitPtrVector&   caloHitsPtrVector() const{return caloHitsPtrVector_;}


	    //Setting parameters
	    void cog3Vector  (CLHEP::Hep3Vector cog3Vector)     {cog3Vector_ = cog3Vector;}
	    void energyRing  (double e1, double e9, double e25) {e1_ = e1;e9_ = e9;e25_ = e25;}
	    void secondMoment(double second)                    {secondMoment_ = second;}
	    void angle       (double angle)                     {angle_ = angle;}
         


	 private:
	    int               diskId_;   
	    double            time_;       
	    double            timeErr_;       
	    double            energyDep_;
	    double            energyDepErr_;
	    double            e1_;
	    double            e9_;
	    double            e25_;	      
	    double            secondMoment_;	      
	    double            angle_;	      
	    CLHEP::Hep3Vector cog3Vector_; 
	    CaloHitPtrVector  caloHitsPtrVector_;
	    unsigned          size_;
            bool              isSplit_;    

   };

   
   typedef std::vector<mu2e::CaloCluster> CaloClusterCollection;

} 

#endif 
