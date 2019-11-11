#ifndef RecoDataProducts_CosmicTrackMCInfo_hh
#define RecoDataProducts_CosmicTrackMCInfo_hh

#include "TMath.h"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "DataProducts/inc/XYZVec.hh"

struct CosmicTrackMCInfo{

     double TrueTheta ;
     double TruePhi ;
     std::vector<double> TrueTimeResiduals;
     std::vector<double> TrueDOCA;
     std::vector<double> Ambig;
     TrackEquation TrueFitEquation;
     TrackAxes TrueTrackCoordSystem;
     TrackParams RawTrueParams;
     CosmicTrackMCInfo();
     CosmicTrackMCInfo(double theta, double phi, TrackEquation eqn, TrackAxes axis, TrackParams para) : TrueTheta(theta), TruePhi(phi), TrueFitEquation(eqn), TrueTrackCoordSystem(axis), RawTrueParams(para) {};

};

#endif
