//
// Utility to summarize information from StepPointMCs into CaloShower objects.
//
// Original author B. Echenard
//


#ifndef HitMakers_CaloShowerStepMCUtil_hh
#define HitMakers_CaloShowerStepMCUtil_hh

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"


namespace mu2e {


     class CaloShowerStepMCUtil {

	   public:
               
	       enum weight_type {none, energy};
	       
	       CaloShowerStepMCUtil(int imax, weight_type type = weight_type::none) : 
	          _imax(imax),_type(type),_n(imax,0),
		  _edep(imax,0),_time(imax,0),_t0(imax,0),
		  _x(imax,0),_y(imax,0),_z(imax,0),
		  _x2(imax,0),_y2(imax,0),_z2(imax,0),
		  _xy(imax,0),_xz(imax,0),_yz(imax,0),
		  _w(imax,0),_w2(imax,0),_pos(0,0,0),_cov(3,0)
	       {};



               void add(int i, double edep, double time, CLHEP::Hep3Vector& pos);
               void setTinit(int i, double time) {_t0.at(i) = time;}
               void reset(int i);
	       void printBucket(int i);
	       
	       int                  nBuckets()     const {return _imax;}
	       int                  entries(int i) const {return _n[i];}
	       double               eDep(int i)    const {return _edep[i];}
	       double               t0(int i)      const {return _t0[i];}
               double               time(int i)    const {return _time.at(i) / _w.at(i);}
               CLHEP::Hep3Vector&   pos(int i);     
	       CLHEP::HepSymMatrix& covPos(int i);

	       
	       

	   private:
	       
	       int                   _imax;
	       weight_type           _type;

	       std::vector<int>      _n;
	       std::vector<double>   _edep;
	       std::vector<double>   _time;
	       std::vector<double>   _t0;
	       std::vector<double>   _x;
	       std::vector<double>   _y;
	       std::vector<double>   _z;
	       std::vector<double>   _x2;
	       std::vector<double>   _y2;
	       std::vector<double>   _z2;
	       std::vector<double>   _xy;
	       std::vector<double>   _xz;
	       std::vector<double>   _yz;
	       std::vector<double>   _w;
	       std::vector<double>   _w2;
	       CLHEP::Hep3Vector     _pos;
	       CLHEP::HepSymMatrix   _cov;
     };

} 

#endif
