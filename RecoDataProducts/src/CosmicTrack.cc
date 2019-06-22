//S. Middleton, Feb 2019
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include <vector>

using namespace std;
namespace mu2e{


	CosmicTrack::CosmicTrack() {
          	_Nhits = 0;
    		
    		_inita0=0.;
		_inita1=0.;
		_initb0=0.;
		_initb1=0.;
		_Sagitta = 0.;
		
    		XYZVec _track_equation(0,0,0);//r(t) expression
    		XYZVec _track_direction(0,0,0);//the "gradient" term
    		XYZVec _track_position(0,0,0);//the "starting point" in fit line
		XYZVec _initial_track_direction(0,0,0);
   		
   			
	 }
    
	double GetSagitta(){
		return 1.;
	}

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear_parameters() {
	  
          _track_parameters.erase(_track_parameters.begin(), _track_parameters.end());
          
	}

	void CosmicTrack::clear_all() {
	  
	  _finalhit_errorsTotal.erase(_finalhit_errorsTotal.begin(),_finalhit_errorsTotal.end());
	  
	  _finalfit_residualsX.erase(_finalfit_residualsX.begin(),_finalfit_residualsX.end());
	  _finalfit_residual_errorsX.erase(_finalfit_residual_errorsX.begin(),_finalfit_residual_errorsX.end());
	  
	  _finalfit_residualsY.erase(_finalfit_residualsY.begin(),_finalfit_residualsY.end());
	  _finalfit_residual_errorsY.erase(_finalfit_residual_errorsY.begin(),_finalfit_residual_errorsY.end());
	  
	  _finalfit_pullsX.erase(_finalfit_pullsX.begin(),_finalfit_pullsX.end());
	  _finalfit_pullsY.erase(_finalfit_pullsY.begin(),_finalfit_pullsY.end());
	  
	  _cov.erase(_cov.begin(),_cov.end());
	 
	}
	
     

	
}

