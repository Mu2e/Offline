#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TMath.h"


using namespace std;

namespace mu2e {


  //look up and return the average speed from vectors
  double StrawDrift::GetAverageSpeed(double dist) const {
    if (dist < _distances_dbins[1]){
      return _distances_dbins[1]/_times_dbins[_phiBins];
    }
    int index = min(int(_distances_dbins.size())-2,max(0,int(floor(dist/_deltaD))));
    double time = _times_dbins[index*_phiBins] + (dist - _distances_dbins[index])/(_deltaD) * (_times_dbins[(index+1)*_phiBins] - _times_dbins[index*_phiBins]);
    return dist/time;
  }

  double StrawDrift::GetInstantSpeedFromD(double dist) const {
    int index = min(int(_distances_dbins.size())-2,max(0,int(floor(dist/_deltaD))));
    return _instantSpeed_dbins[index] + (dist - _distances_dbins[index])/(_deltaD) * (_instantSpeed_dbins[index+1] - _instantSpeed_dbins[index]);
  }

  double StrawDrift::GetInstantSpeedFromT(double time) const
  {
    int lowerIndex = 0;
    //step through and find time larger than what is specified
    for (size_t i=0; i < (_distances_dbins.size() - 2); i++) {
      int fullIndex = i*(_phiBins) + 0; //map from a 2D index to a 1D index (at phi=0)
      if(time < _times_dbins[fullIndex]){
        lowerIndex = i;
        break;
      }
    }

    return _instantSpeed_dbins[lowerIndex] + (time - _times_dbins[lowerIndex*_phiBins])/(_times_dbins[(lowerIndex+1)*_phiBins]-_times_dbins[lowerIndex*_phiBins]) * (_instantSpeed_dbins[lowerIndex+1]-_instantSpeed_dbins[lowerIndex]);
  }

  //D2T for sims
  double StrawDrift::D2T(double distance, double phi) const {
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float fphi = foldPhi(phi);
    size_t phirange[2];
    phiRange(fphi,phirange);
    auto dindex = distIndex(distance);
    double lowerTime = _times_dbins[dindex*_phiBins+phirange[0]] + (distance - _distances_dbins[dindex])/(_deltaD) * (_times_dbins[(dindex+1)*_phiBins+phirange[0]] - _times_dbins[dindex*_phiBins+phirange[0]]);
    double upperTime = _times_dbins[dindex*_phiBins+phirange[1]] + (distance - _distances_dbins[dindex])/(_deltaD) * (_times_dbins[(dindex+1)*_phiBins+phirange[1]] - _times_dbins[dindex*_phiBins+phirange[1]]);
    if (phi == 0) return lowerTime;
    double lowerPhi = phirange[0] * _deltaPhi;
    return lowerTime + (fphi - lowerPhi)/_deltaPhi * (upperTime - lowerTime);
  }

  //T2D for reco
  double StrawDrift::T2D(double time, double phi, bool nonnegative) const {
    if (time < 0 && nonnegative) return 0;
    double fphi = foldPhi(phi);
    size_t phirange[2];
    phiRange(fphi,phirange);
    size_t tindex = timeIndex(time);
    size_t lowrange[2];
    lowrange[0] = tindex*_phiBins+phirange[0];
    lowrange[1] = (tindex+1)*_phiBins+phirange[0];
    double lowerDist = _distances_tbins[lowrange[0]] + (time - _times_tbins[tindex])/(_deltaT) * (_distances_tbins[lowrange[1]] - _distances_tbins[lowrange[0]]);
    if (fphi == 0) return lowerDist;
    size_t uprange[2];
    uprange[0] = tindex*_phiBins+phirange[1];
    uprange[1] = (tindex+1)*_phiBins+phirange[1];
    double upperDist = _distances_tbins[uprange[0]] + (time - _times_tbins[tindex])/(_deltaT) * (_distances_tbins[uprange[1]] - _distances_tbins[uprange[0]]);
    double lowerPhi = phirange[0] * _deltaPhi;
    return lowerDist + (fphi - lowerPhi)/_deltaPhi * (upperDist - lowerDist);
  }

  double StrawDrift::foldPhi(double phi) const {
    // fold phi onto [0,pi/2]
    phi = fabs(fmod(phi,M_PI));
    if (phi > M_PI_2) phi = M_PI-phi;
    return phi;
  }

  void StrawDrift::phiRange(double phi, size_t prange[2]) const {
    //for interpolation, define a high and a low index
    prange[0] = static_cast<size_t>(std::max(0,static_cast<int>(floor(phi/_deltaPhi))));
    prange[1] = static_cast<size_t>(std::min(static_cast<int>(_phiBins)-1,static_cast<int>(ceil(phi/_deltaPhi))));
  }

  size_t StrawDrift::timeIndex(double time) const {
    return static_cast<size_t>( std::min(static_cast<int>(_times_tbins.size())-2,std::max(0,static_cast<int>(floor(time/_deltaT)))));
  }

  size_t StrawDrift::distIndex(double dist) const {
    return static_cast<size_t>(std::min(static_cast<int>(_distances_dbins.size())-2,std::max(0,static_cast<int>(floor(dist/_deltaD)))));
  }

  void StrawDrift::indexRange(double time, double phi, size_t irange[2]) const {
    size_t phirange[2];
    phiRange(phi,phirange);
    size_t tindex = timeIndex(time);
    irange[0] = tindex*_phiBins+phirange[0];
    irange[1] = (tindex+1)*_phiBins+phirange[0];
  }

  void StrawDrift::print(std::ostream& os) const {
    size_t n = _times_dbins.size();
    size_t nd = _distances_dbins.size();
    os << endl << "StrawDrift parameters: "  << std::endl
      << "Times (size=" << n << ") Distances (size=" << nd << "): " << endl;
    os << "  effectiveSpeed = " << _distances_dbins[1]/_times_dbins[1*_phiBins] << " "
      << _distances_dbins[1]/_times_dbins[1*_phiBins+1]
      << " ... " << _distances_dbins[nd-1]/_times_dbins[n-2] << " "
      << _distances_dbins[nd-1]/_times_dbins[n-1] << endl;
    os << "  time = " << _times_dbins[0] << " "
      << _times_dbins[1]
      << " ... " << _times_dbins[n-2] << " "
      << _times_dbins[n-1] << endl;
    os << "  phi width= " << _deltaPhi << endl;

    os << "distances= " << _distances_dbins[0] << " " << _distances_dbins[1]
      << " ... " << _distances_dbins[nd-2] << " " << _distances_dbins[nd-1] << endl;
    os << "instantSpeeds = " << _instantSpeed_dbins[0] << " "<<  _instantSpeed_dbins[1]
      << " ... " << _instantSpeed_dbins[nd-2] << " " << _instantSpeed_dbins[nd-1] << endl;
  }
}

