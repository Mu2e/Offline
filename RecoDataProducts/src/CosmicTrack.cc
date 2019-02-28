#include "RecoDataProducts/inc/CosmicTrack.hh"

#include <vector>

using namespace std;

namespace mu2e{

	CosmicTrack::CosmicTrack() {
          	_Nhits = 0;

    		_par1=0.;
   		_par2=0.;
    		_par3=0.;
    		_par4=0.;

    		XYZVec _track_parameters(0,0,0);
    		XYZVec _track_equation(0,0,0);//r(t) expression

    		XYZVec _track_direction(0,0,0);//the "gradient" term
    		XYZVec _track_point0(0,0,0);//the "starting point" in fit line

   		_chisq=0;
    		_chisq_dof=0; 

   
	  
	 }

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear() {
	  _par1 = 0.;
   	  _par2 = 0.;
          _par3 = 0.;
          _par4 = 0.;

	  _Nhits = 0;
	  
	  _fit_residuals.erase(_fit_residuals.begin(),_fit_residuals.end());
	  _fit_residual_errors.erase(_fit_residual_errors.begin(),_fit_residual_errors.end());
	  
	}
	
     

	
}

