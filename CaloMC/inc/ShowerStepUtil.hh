//
// Utility to summarize information from StepPointMCs into CaloShower objects.
//

#ifndef CaloMC_ShowerStepUtil_hh
#define CaloMC_ShowerStepUtil_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <vector>


namespace mu2e {


     class ShowerStepUtil {

	   public:
               
	       enum weight_type {none, energy};
	       
	       ShowerStepUtil(int imax, weight_type type = weight_type::none) : 
	          imax_(imax),type_(type),n_(imax,0),
		  edep_(imax,0),pIn_(imax,0),time_(imax,0),t0_(imax,0),
		  x_(imax,0),y_(imax,0),z_(imax,0),
		  x2_(imax,0),y2_(imax,0),z2_(imax,0),
		  xy_(imax,0),xz_(imax,0),yz_(imax,0),
		  w_(imax,0),w2_(imax,0),posIn_(imax,CLHEP::Hep3Vector(0,0,0)),
		  pos_(0,0,0),cov_(3,0)
	       {};


               void add(int i, double edep, double time, double momentum, CLHEP::Hep3Vector& pos);
               void init(int i, double time, double momentum, const CLHEP::Hep3Vector& posIn);
               void reset(int i);
	       void printBucket(int i);
	       
	       int                  nBuckets()         const {return imax_;}
	       int                  entries(int i)     const {return n_.at(i);}
	       double               energyDep(int i)   const {return edep_.at(i);}
	       double               t0(int i)          const {return t0_.at(i);}
               double               pIn(int i)         const {return pIn_.at(i);}
               double               time(int i)        const {return time_.at(i) / w_.at(i);}
               const CLHEP::Hep3Vector&   posIn(int i) const {return posIn_.at(i);}     
               CLHEP::Hep3Vector&   pos(int i);     
	       CLHEP::HepSymMatrix& covPos(int i);

	       
	       

	   private:
	       
	       int                   imax_;
	       weight_type           type_;

	       std::vector<int>      n_;
	       std::vector<double>   edep_;
	       std::vector<double>   pIn_;
	       std::vector<double>   time_;
	       std::vector<double>   t0_;
	       std::vector<double>   x_;
	       std::vector<double>   y_;
	       std::vector<double>   z_;
	       std::vector<double>   x2_;
	       std::vector<double>   y2_;
	       std::vector<double>   z2_;
	       std::vector<double>   xy_;
	       std::vector<double>   xz_;
	       std::vector<double>   yz_;
	       std::vector<double>   w_;
	       std::vector<double>   w2_;
	       std::vector<CLHEP::Hep3Vector> posIn_;

	       CLHEP::Hep3Vector     pos_;
	       CLHEP::HepSymMatrix   cov_;
     };

} 

#endif
