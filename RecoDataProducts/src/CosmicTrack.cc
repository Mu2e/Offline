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
    		_inita0=0.;
		_inita1=0.;
		_initb0=0.;
		_initb1=0.;
		_Sagitta = 0.;
		
    		XYZVec _track_equation(0,0,0);//r(t) expression
    		XYZVec _track_direction(0,0,0);//the "gradient" term
    		XYZVec _track_position(0,0,0);//the "starting point" in fit line
		XYZVec _initial_track_direction(0,0,0);
   		_initchisq=0;
    		_initchisq_dof=0; 
		_finalchisq=0;
    		_finalchisq_dof=0; 	
	 }
    
	double GetSagitta(){
		return 1.;
	}

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear_parameters() {
	 _a0 = 0.;
   	  _a1 = 0.;
          _b0 = 0.;
          _b1 = 0.;  
          _track_parameters.erase(_track_parameters.begin(), _track_parameters.end());
          
	}

	void CosmicTrack::clear_all() {
	  //_a0 = 0.;
   	  //_a1 = 0.;
          //_b0 = 0.;
          //_b1 = 0.;
          //_inita0=0.;
	  //_inita1=0.;
          //_initb0=0.;
          //_initb1=0.;
          //_track_parameters.erase(_track_parameters.begin(), _track_parameters.end());
          //_initial_track_parameters.erase(_initial_track_parameters.begin(), _initial_track_parameters.end());
          
	  _Nhits = 0;
          _Sagitta = 0 ;
          
          _finalchisq = 0;
    	  _finalchisq_dof = 0;
    
    	  //_initchisq = 0;
    	  //_initchisq_dof = 0;
    
	  //_inithit_errorsTotal.erase(_inithit_errorsTotal.begin(),_inithit_errorsTotal.end());
	  _finalhit_errorsTotal.erase(_finalhit_errorsTotal.begin(),_finalhit_errorsTotal.end());
	  
	  //_initfit_residualsX.erase(_initfit_residualsX.begin(),_initfit_residualsX.end());
	  //_initfit_residual_errorsX.erase(_initfit_residual_errorsX.begin(),_initfit_residual_errorsX.end());
	  
	  _finalfit_residualsX.erase(_finalfit_residualsX.begin(),_finalfit_residualsX.end());
	  _finalfit_residual_errorsX.erase(_finalfit_residual_errorsX.begin(),_finalfit_residual_errorsX.end());
	  
	  //_initfit_residualsY.erase(_initfit_residualsY.begin(),_initfit_residualsY.end());
	  //_initfit_residual_errorsY.erase(_initfit_residual_errorsY.begin(),_initfit_residual_errorsY.end());
	  
	  _finalfit_residualsY.erase(_finalfit_residualsY.begin(),_finalfit_residualsY.end());
	  _finalfit_residual_errorsY.erase(_finalfit_residual_errorsY.begin(),_finalfit_residual_errorsY.end());
	  
	  //_res_MC_X.erase(_res_MC_X.begin(), _res_MC_X.end());
	  //_res_MC_Y.erase(_res_MC_Y.begin(), _res_MC_Y.end());
	  
	  _niters.erase(_niters.begin(),_niters.end());
	  
	  
	}
	
     

	
}

