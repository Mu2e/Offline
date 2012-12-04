#include "IlcDCHrecoParam.h"

//ClassImp(IlcDCHrecoParam)

IlcDCHrecoParam::IlcDCHrecoParam(){
  fRoadY =1.;
  fRoadZ=3.;
  fSafetyZ=10.;
  fNearestByC=0.004;          
  fNearestByTgl=0.6;  
  fNearestByDistance2=325;
  fMinNOfClusetrs=7;
  fMinDensity=0.4;
  fSeedTrackSigma2=0.02;
  fCutTrackSigma2=0.6;
  fSeedBeamSigma2=4.*4.;
  fSeedBeamSigmaZ2=100.*100.;

  fMinDCAr=5;
  fMinDCAz=5;
  fScaleForQuilcty=3;
  fDeltaZ=60;
  fDistFromEdge=10;

  fMinRofV0=5; 
  fMaxRofV0=150; 
  fCutOnDistanceInV0=4; 
}
