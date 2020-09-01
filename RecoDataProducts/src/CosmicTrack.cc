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

namespace mu2e{

	CosmicTrack::CosmicTrack() {
    		
    		InitParams.A0 = 0;
          	InitParams.A1 = 0;
          	InitParams.B0 = 0;
          	InitParams.B1 = 0;
          	InitParams.T0 = 0;
          	
	 }
    


	// Destructor
	CosmicTrack::~CosmicTrack() {}

	
	
	
}
