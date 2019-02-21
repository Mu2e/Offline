#ifndef KALMAN_TRACKFIT_HH
#define KALMAN_TRACKFIT_HH

#include <string>
#include <map>

#include "RecoDataProducts/inc/KalmanTrack.hh"

#include "TMatrix.h"

namespace mu2e{
namespace Kalman {
  /*
  //Function to get residuals:
  State CalculateResidual(const State& state_1, const State& state_2);

  //Function to update chi squared:
  double CalculateChi2Update(const State& state);

  //Function to get Pull:
  State CalculatePull(unsigned int i) const;

  //Function to calculate predicted residual/pull 
  State CalculatePredictedResidual(unsigned int i) const;

  //Function to calculate filtered residual 
  State CalculateFilteredResidual(unsigned int i) const;

  //Function to calculate smoothed residual 
  State CalculateSmoothedResidual(unsigned int i) const;


  //Set Seed:
  void SetSeed(State seed);
  //Get Seed:
  State GetSeed() const { return _seed; }
  //Set Kalman Track:
  void SetTrack(KalTrack track);
  //Get Kalman Track
  KalTrack& GetTrack() { return _track; }
  const KalTrack& GetTrack() const { return _track; }

  //Calculate propagator matrix:
  virtual TMatrixD CalculatePropagator(const Hit& first, const Hit& second) = 0;
  //Function to Propoagte:
  void Propagate(const Hit& from, Hit& to) const;
  virtual TMatrixD CalculateProcessNoise(const Hit& start, const Hit& end) = 0;

  //Function to do Prediction:
  void Predict(Hit& hit) const 
  //Function to do filtering:
  void Filter(Hit& hit) const;
  //Function to do smoothing:
  void Smooth(Hit& hit) const;
  
  unsigned int GetDimension() const { return _dimension; }
  */
}
}
#endif

