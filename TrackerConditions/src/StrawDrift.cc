
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TMath.h"

#include "TrackerConditions/inc/StrawDrift.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


using namespace std;

namespace mu2e {
  

  //look up and return the average speed from vectors
  double StrawDrift::GetAverageSpeed(double dist) const {
    if (dist < _distances[1]){
      return _distances[1]/_times[_phiBins];
    }
    int index = min(int(_distances.size())-2,max(0,int(floor(dist/_deltaD))));
    double time = _times[index*_phiBins] + (dist - _distances[index])/(_deltaD) * (_times[(index+1)*_phiBins] - _times[index*_phiBins]);
    return dist/time;
  }
  
  double StrawDrift::GetInstantSpeedFromD(double dist) const {
    int index = min(int(_distances.size())-2,max(0,int(floor(dist/_deltaD))));
    return _instantSpeed[index] + (dist - _distances[index])/(_deltaD) * (_instantSpeed[index+1] - _instantSpeed[index]);
  }
  
  double StrawDrift::GetInstantSpeedFromT(double time) const
  {
    int lowerIndex = 0;
    //step through and find time larger than what is specified
    for (size_t i=0; i < (_distances.size() - 2); i++) {
      int fullIndex = i*(_phiBins) + 0; //map from a 2D index to a 1D index (at phi=0)
      if(time < _times[fullIndex]){
        lowerIndex = i;
        break;
      }
    }

    return _instantSpeed[lowerIndex] + (time - _times[lowerIndex*_phiBins])/(_times[(lowerIndex+1)*_phiBins]-_times[lowerIndex*_phiBins]) * (_instantSpeed[lowerIndex+1]-_instantSpeed[lowerIndex]);
  }
  
  //D2T for sims
  double StrawDrift::D2T(double distance, double phi) const {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down

    int index = min(int(_distances.size())-2,max(0,int(floor(distance/_deltaD))));
    double lowerTime = _times[index*_phiBins+lowerPhiIndex] + (distance - _distances[index])/(_deltaD) * (_times[(index+1)*_phiBins+lowerPhiIndex] - _times[index*_phiBins+lowerPhiIndex]);
    double upperTime = _times[index*_phiBins+upperPhiIndex] + (distance - _distances[index])/(_deltaD) * (_times[(index+1)*_phiBins+upperPhiIndex] - _times[index*_phiBins+upperPhiIndex]);
    if (phi == 0)
      return lowerTime;

    double lowerPhi = lowerPhiIndex * phiSliceWidth;
    return lowerTime + (reducedPhi - lowerPhi)/phiSliceWidth * (upperTime - lowerTime);
  }
  
  //T2D for reco
  double StrawDrift::T2D(double time, double phi) const {
    if (time < 0)
      return 0;
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down

    int index = 0;
    //step through and find time larger than what is specified
    for (size_t i=0; i < (_distances.size() - 2); i++) {
      int fullIndex = i*(_phiBins) + lowerPhiIndex; //map from a 2D index to a 1D index (at phi=0)
      if(time < _times[fullIndex]){
        index = i;
        break;
      }
      index = i;
    }

    double lowerDist = _distances[index] + (time - _times[index*_phiBins+lowerPhiIndex])/(_times[(index+1)*_phiBins+lowerPhiIndex]-_times[index*_phiBins+lowerPhiIndex]) * _deltaD;
    double upperDist = _distances[index] + (time - _times[index*_phiBins+upperPhiIndex])/(_times[(index+1)*_phiBins+upperPhiIndex]-_times[index*_phiBins+upperPhiIndex]) * _deltaD;
    if (phi == 0)
      return lowerDist;

    double lowerPhi = lowerPhiIndex * phiSliceWidth;
    return lowerDist + (reducedPhi - lowerPhi)/phiSliceWidth * (upperDist - lowerDist);
  }
  
  double StrawDrift::ConstrainAngle(double phi) const {
    if (phi < 0) {
      phi = -1.0*phi;
    }
    phi = fmod(phi,TMath::Pi());
    if (phi > TMath::Pi()/2.0) {
      phi = phi - 2.0*fmod(phi,TMath::Pi()/2.0);
    }
    return phi;
  }
 
  void StrawDrift::print(std::ostream& os) const {
    size_t n = _times.size();
    size_t nd = _distances.size();
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    os << endl << "StrawDrift parameters: "  << std::endl
       << "Times (size=" << n << ") Distances (size=" << nd << "): " << endl;
    os << "  effectiveSpeed = " << _distances[1]/_times[1*_phiBins] << " " 
       << _distances[1]/_times[1*_phiBins+1]
       << " ... " << _distances[nd-1]/_times[n-2] << " " 
       << _distances[nd-1]/_times[n-1] << endl;
    os << "  time = " << _times[0] << " " 
       << _times[1]
       << " ... " << _times[n-2] << " " 
       << _times[n-1] << endl;
    os << "  phi width= " << phiSliceWidth << endl;

    os << "distances = " << _distances[0] << " " << _distances[1]
       << " ... " << _distances[nd-2] << " " << _distances[nd-1] << endl;
    os << "instantSpeeds = " << _instantSpeed[0] << " "<<  _instantSpeed[1]
       << " ... " << _instantSpeed[nd-2] << " " << _instantSpeed[nd-1] << endl;
  }
}

