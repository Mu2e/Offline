#ifndef CaloMC_ShowerStepUtil_hh
#define CaloMC_ShowerStepUtil_hh
//
// Utility to summarize information from StepPointMCs into CaloShower objects.
//
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

     class ShowerStepUtil
     {
	   public:               
	       enum weight_type {none, energy};
	       
	       ShowerStepUtil(unsigned imax, weight_type type = weight_type::none) : 
	          imax_(imax),type_(type),n_(imax,0),eDepG4_(imax,0),eDepVis_(imax,0),
		  pIn_(imax,0),time_(imax,0),t0_(imax,0),x_(imax,0),y_(imax,0),z_(imax,0),
		  w_(imax,0),pos_(0,0,0)
	       {};

               void add(unsigned i, double eDepG4, double eDepVis, double time, double momentum, CLHEP::Hep3Vector& pos);
               void reset(unsigned i);
	       void printBucket(unsigned i);
	       
	       unsigned             nBuckets()              const {return imax_;}
	       unsigned             entries(unsigned i)     const {return n_.at(i);}
	       double               energyG4(unsigned i)    const {return eDepG4_.at(i);}
	       double               energyVis(unsigned i)   const {return eDepVis_.at(i);}
	       double               t0(unsigned i)          const {return t0_.at(i);}
               double               pIn(unsigned i)         const {return pIn_.at(i);}
               double               time(unsigned i)        const {return time_.at(i) / w_.at(i);}
               CLHEP::Hep3Vector&   pos(unsigned i);     


	   private:	       
	       unsigned              imax_;
	       weight_type           type_;
	       std::vector<unsigned> n_;
	       std::vector<double>   eDepG4_;
	       std::vector<double>   eDepVis_;
	       std::vector<double>   pIn_;
	       std::vector<double>   time_;
	       std::vector<double>   t0_;
	       std::vector<double>   x_;
	       std::vector<double>   y_;
	       std::vector<double>   z_;
	       std::vector<double>   w_;
	       CLHEP::Hep3Vector     pos_;
     };

} 

#endif
