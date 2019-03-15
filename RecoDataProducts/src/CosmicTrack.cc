//S. Middleton, Feb 2019
#include "RecoDataProducts/inc/CosmicTrack.hh"


#include <vector>

using namespace std;

namespace mu2e{


	CosmicTrack::CosmicTrack() {
          	_Nhits = 0;

    		_a0=0.;
   		_a1=0.;
    		_b0=0.;
    		_b1=0.;
		_Sagitta = 0.;
   		XYZVec _track_parameters(0,0,0);
    		XYZVec _track_equation(0,0,0);//r(t) expression
    		XYZVec _track_direction(0,0,0);//the "gradient" term
    		XYZVec _track_point0(0,0,0);//the "starting point" in fit line

   		_chisq=0;
    		_chisq_dof=0; 

   
	  
	 }
    
	double GetSagitta(){
		//TODO ADD MATHS HERE
		return 1.;
	}

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear() {
	  _a0 = 0.;
   	  _a1 = 0.;
          _b0 = 0.;
          _b1 = 0.;

	  _Nhits = 0;
          _Sagitta = 0 ;
          
	  _hit_errorsTotal.erase(_hit_errorsTotal.begin(),_hit_errorsTotal.end());
	  
	  _fit_residualsX.erase(_fit_residualsX.begin(),_fit_residualsX.end());
	  _fit_residual_errorsX.erase(_fit_residual_errorsX.begin(),_fit_residual_errorsX.end());
	  
	  _fit_residualsY.erase(_fit_residualsY.begin(),_fit_residualsY.end());
	  _fit_residual_errorsY.erase(_fit_residual_errorsY.begin(),_fit_residual_errorsY.end());
	  
	}
	
     

	
}

