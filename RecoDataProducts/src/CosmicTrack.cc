//S. Middleton, Aug 2019
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include <vector>

using namespace std;

TrackParams::TrackParams(){
	A0 = 0.;
	A1 = 0.;
	B0 = 0.;
	B1 = 0.;
	T0 = 0.;
} 

TrackCov::TrackCov(){
	sigA0 = 0.;
  	sigA1 = 0.;
  	sigB0 = 0.;
  	sigB1 = 0.;
} 

TrackAxes::TrackAxes(){
	_XDoublePrime.SetXYZ(0,0,0);
	_YDoublePrime.SetXYZ(0,0,0);
	_ZPrime.SetXYZ(0,0,0);

}
TrackEquation::TrackEquation(){
	Pos.SetXYZ(0,0,0);
	Dir.SetXYZ(0,0,0);
} 

TrackSeedDiag::TrackSeedDiag(){
   	FinalChiX = 0;
   	FinalChiY = 0;
   	FinalChiTot = 0;
   	
   	InitialChiX = 0;
   	InitialChiY = 0;
   	InitialChiTot = 0;
   	
	}
TrackDriftDiag::TrackDriftDiag(){
   	FinalChiX = 0;
   	FinalChiY = 0;
   	FinalChiTot = 0;

	}
namespace mu2e{

	CosmicTrack::CosmicTrack() {
          	_Nhits = 0;
    		
    		InitParams.A0 = 0;
          	InitParams.A1 = 0;
          	InitParams.B0 = 0;
          	InitParams.B1 = 0;
          	InitParams.T0 = 0;
          	
		_Sagitta = 0.;
		
    		//Pos.SetXYZ(0,0,0);
	        Direction.SetXYZ(0,0,0);
    		
   		
   			
	 }
    
	double GetSagitta(){
		return 1.;
	}

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear_parameters() {
          FitParams.A0 = 0;
          FitParams.A1 = 0;
          FitParams.B0 = 0;
          FitParams.B1 = 0;
          FitParams.T0 = 0;
	}
	
	
	void CosmicTrack::clear_diag(){
          Diag.FinalErrX.erase(Diag.FinalErrX.begin(),Diag.FinalErrX.end());
	  Diag.FinalErrY.erase(Diag.FinalErrY.begin(),Diag.FinalErrY.end());
	  Diag.FinalResidualsX.erase(Diag.FinalResidualsX.begin(),Diag.FinalResidualsX.end());
	  Diag.FinalResidualsY.erase(Diag.FinalResidualsY.begin(),Diag.FinalResidualsY.end());

	}

	void CosmicTrack::clear_errors() {  
	 
	  Diag.FinalErrX.erase(Diag.FinalErrX.begin(),Diag.FinalErrX.end());
	  Diag.FinalErrY.erase(Diag.FinalErrY.begin(),Diag.FinalErrY.end());
	  		  
	}
	void CosmicTrack::clear_cov() {  
	 
	  FitParams.Covarience.sigA0 = 0.;
	  FitParams.Covarience.sigA1 = 0.;
	  FitParams.Covarience.sigB0 = 0.;
	  FitParams.Covarience.sigB1 = 0.;		  
	}
	
	
	
}

