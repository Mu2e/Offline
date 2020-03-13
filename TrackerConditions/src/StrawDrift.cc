
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
  


    //step through and find distance larger than what is specified
  size_t StrawDrift::lowerDistanceBin(double dist) const {
    for (size_t i=0; i < (_distances.size() - 1); i++) {
      if(dist >= _distances[i]) return i;
    }
    return 0;
  }

  //look up and return the average speed from vectors
  double StrawDrift::GetAverageSpeed(double dist) const {
    return _averageSpeeds[ lowerDistanceBin(dist) ];
  }
  
  double StrawDrift::GetInstantSpeedFromD(double dist) const {
    return _instantSpeeds[ lowerDistanceBin(dist) ];
  }
  
  double StrawDrift::GetInstantSpeedFromT(double time) const
  {
    int lowerIndex = 0;
    //step through and find time larger than what is specified
    for (size_t i=0; i < (_distances.size() - 1); i++) {
      int fullIndex = i*(_phiBins) + 0; //map from a 2D index to a 1D index (at phi=0)
      if(time >= _D2Tinfos[fullIndex].time){
        lowerIndex = i;
        break;
      }
    }

    return _instantSpeeds[lowerIndex];
  }
  
  double StrawDrift::GetGammaFromD(double distance, double phi) const 
  {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, 
    // the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth);
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth);
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperGamma = 0;
    float lowerGamma = 0;
    for (size_t k=0; k < (_distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;//mapping from a 2D to a 1D index
      if (distance >= _D2Tinfos[fullIndex].distance){
        upperGamma = _D2Tinfos[fullIndex].gamma;//set the gamma associated with the higher index
        fullIndex = k*(_phiBins)+lowerPhiIndex; //mapping from a 2D to a 1D index
        lowerGamma = _D2Tinfos[fullIndex].gamma;//set the gamma associated with the lower index
        break;
      }
    }
    double Gamma = lowerGamma*lowerPhiWeight + upperGamma*upperPhiWeight;//compute the final gamma
    return Gamma;
  }
  
  double StrawDrift::GetGammaFromT(double time, double phi) const 
  {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperGamma = 0;
    float lowerGamma = 0;
    for (size_t k=0; k < (_distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;//mapping from a 2D to a 1D index
      if (time >= _D2Tinfos[fullIndex].time){
        upperGamma = _D2Tinfos[fullIndex].gamma;//set the gamma associated with the higher index
        fullIndex = k*(_phiBins)+lowerPhiIndex; //mapping from a 2D to a 1D index
        lowerGamma = _D2Tinfos[fullIndex].gamma;//set the gamma associated with the lower index
        break;
      }
    }
    double Gamma = lowerGamma*lowerPhiWeight + upperGamma*upperPhiWeight;//compute the final gamma
    return Gamma;
  }
  
  //look up and return the lorentz corrected r componenent of the average velocity
  double StrawDrift::GetEffectiveSpeed(double dist, double phi) const {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = fmod(phi,TMath::Pi()/2.0);
    if (reducedPhi < 0){
      reducedPhi += TMath::Pi()/2.0;
    };
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperSpeed = 0;
    float lowerSpeed = 0;
    float effectiveSpeed = 0;
    
    for (size_t k=0; k < (_distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (dist >= _D2Tinfos[fullIndex].distance){
        upperSpeed = _D2Tinfos[fullIndex].effectiveSpeed; //set the higher speed
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerSpeed = _D2Tinfos[fullIndex].effectiveSpeed; // set the lower speed
        break;
      }
    }
    effectiveSpeed = lowerSpeed*lowerPhiWeight + upperSpeed*upperPhiWeight;
    return effectiveSpeed;
  }
  
  //D2T for sims
  double StrawDrift::D2T(double distance, double phi) const {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperTime = 0;
    float lowerTime = 0;
    float time = 0;
    float gammaTest = 0.;
    for (size_t k=0; k < (_distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (distance >= _D2Tinfos[fullIndex].distance){
        upperTime = _D2Tinfos[fullIndex].time;//set the higher time
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerTime = _D2Tinfos[fullIndex].time;// set the lower time
        gammaTest = _D2Tinfos[fullIndex].gamma;//just another test
        break;
      }
    }
    time = lowerTime*lowerPhiWeight + upperTime*upperPhiWeight;//compute the final time
    gammaTest = 1.0*gammaTest;
    return time;
  }
  
  //T2D for reco
  double StrawDrift::T2D(double time, double phi) const {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = fmod(phi,TMath::Pi()/2.0);
    if (reducedPhi < 0){
      reducedPhi += TMath::Pi()/2.0;
    };
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperDist = 0;
    float lowerDist = 0;
    float distance = 0;
    for (size_t k=0; k < (_distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (time >= _D2Tinfos[fullIndex].time){
        upperDist = _D2Tinfos[fullIndex].distance; //set the higher distance
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerDist = _D2Tinfos[fullIndex].distance;//set the lower distance
        break;
      }
    }
    distance = lowerDist*lowerPhiWeight + upperDist*upperPhiWeight;//compute the final distance
    return distance;
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
    size_t n = _D2Tinfos.size();
    os << endl << "StrawDrift parameters: "  << std::endl
       << "DTInfo (size=" << n << "):" << endl;
    os << "  gamma = " << _D2Tinfos[0].gamma << " " 
       << _D2Tinfos[0].gamma
       << " ... " << _D2Tinfos[n-2].gamma << " " 
       << _D2Tinfos[n-1].gamma << endl;
    os << "  effectiveSpeed = " << _D2Tinfos[0].effectiveSpeed << " " 
       << _D2Tinfos[0].effectiveSpeed
       << " ... " << _D2Tinfos[n-2].effectiveSpeed << " " 
       << _D2Tinfos[n-1].effectiveSpeed << endl;
    os << "  time = " << _D2Tinfos[0].time << " " 
       << _D2Tinfos[0].time
       << " ... " << _D2Tinfos[n-2].time << " " 
       << _D2Tinfos[n-1].time << endl;
    os << "  phi = " << _D2Tinfos[0].phi << " " 
       << _D2Tinfos[0].phi
       << " ... " << _D2Tinfos[n-2].phi << " " 
       << _D2Tinfos[n-1].phi << endl;
    os << "  distance = " << _D2Tinfos[0].distance << " " 
       << _D2Tinfos[0].distance
       << " ... " << _D2Tinfos[n-2].distance << " " 
       << _D2Tinfos[n-1].distance << endl;
    os << "  instantaneousSpeed = " << _D2Tinfos[0].instantaneousSpeed << " " 
       << _D2Tinfos[0].instantaneousSpeed
       << " ... " << _D2Tinfos[n-2].instantaneousSpeed << " " 
       << _D2Tinfos[n-1].instantaneousSpeed << endl;

    n = _distances.size();
    os << "distances = " << _distances[0] << " " << _distances[1]
       << " ... " << _distances[n-2] << " " << _distances[n-1] << endl;
    n = _instantSpeeds.size();
    os << "instantSpeeds = " << _instantSpeeds[0] << " "<<  _instantSpeeds[1]
       << " ... " << _instantSpeeds[n-2] << " " << _instantSpeeds[n-1] << endl;
    n = _averageSpeeds.size();
    os << "averageSpeeds = " << _averageSpeeds[0] << " " << _averageSpeeds[1]
       << " ... " << _averageSpeeds[n-2] << " " << _averageSpeeds[n-1] << endl;

  }
}

