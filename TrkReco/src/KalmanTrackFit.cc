#include "TrkReco/inc/KalmanTrackFit.hh"
#include "RecoDataProducts/inc/KalmanTrack.hh"

#include <iostream>
#include <sstream>

namespace mu2e {
namespace Kalman {
/*
  State KalmanTrackFit::CalculateResidual(const State& state1, const State& state2) {
    if (state1._dimension != state2._dimension) {
      throw "States have different dimensions Kalman::CalculateResidual()";
    }
    TMatrixD vector = state1._vector - state2._vector;
    TMatrixD covariance = state1._covariance + state2._covariance;

    State residual(vector, covariance);

    return residual;
  }

  double KalmanTrackFitCalculateChiSquaredUpdate(const State& state) {
    TMatrixD rT(TMatrixD::kTransposed, state._vector);
    TMatrixD RI(TMatrixD::kInverted, state._covariance);
    TMatrixD update = rT * RI * state._vector;
    return update(0, 0);
  }


  void KalmanTrackFit::Propagate(const Hit& first, Hit& second) const {
    
    propagator_matrix = this->CalculatePropagator(first, second);
    TMatrixD propT(TMatrixD::kTransposed, propagator_matrix);
    noise_matrix = this->CalculateProcessNoise(first, second);

    TMatrixD new_vec(GetDimension(), 1);
    new_vec = propagator_matrix * first.GetFiltered().GetVector();
    TMatrixD new_cov(GetDimension(), GetDimension());
    new_cov = propagator_matrix * first.GetFiltered().GetCovariance() * propT + noise_matrix;

    second.SetPredicted(State(new_vec, new_cov));
  }
  

  void KalmanTrackFit::Predict(Hit& hit) const {
    if (hit.HasData()) {
    }
  }
  void KalmanTrackFit::Filter(Hit& hit) const {
    if (hit.HasData()) {
      const State& data = hit.GetData();
      const State& predicted = hit.GetPredicted();

      //State measured = _measurements.at(hit.GetId())->Measure(predicted);
      //State pull = CalculateResidual(data, measured);
      //TMatrixD V = data.GetCovariance();
      //TMatrixD H = _measurements.at(tp.GetId())->GetMeasurementMatrix();
      //TMatrixD HT(TMatrixD::kTransposed, H);

      //TMatrixD cov_inv(TMatrixD::kInverted, pull.GetCovariance());

      //TMatrixD K = predicted.GetCovariance() * HT * cov_inv;
      //TMatrixD KT(TMatrixD::kTransposed, K);

      //TMatrixD gain = _identity_matrix - K*H;
      //TMatrixD gainT(TMatrixD::kTransposed, gain);
      //TMatrixD gain_constant = K * V * KT;

      //TMatrixD temp_vec = predicted.GetVector() + K * pull.GetVector();
      //TMatrixD temp_cov = gain * predicted.GetCovariance() * gainT + gain_constant;

      //tp.SetFiltered(State(temp_vec, temp_cov));
    //} else {
    //  tp.SetFiltered(tp.GetPredicted());
    //}
  }
}

void KalmanTrackFit::Smooth(Hit& hit) const {
    if (hit.HasData()) {

}
}
*/
}
}
