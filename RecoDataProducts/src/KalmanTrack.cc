#include "RecoDataProducts/inc/KalmanTrack.hh"

#include <iostream>
#include <sstream>
#include <exception>
#include <numeric>
#include <sstream>

using namespace std;

namespace mu2e {
namespace Kalman {

/*-----------STATE MEMBER FUNCTIONS-------------//
//----------------------------------------------*/

  State::State(unsigned int dim) :
    _dimension(dim),
    _vector(dim, 1),
    _covariance(dim, dim),
    _has_value(false) {
    }

  State::State(TMatrixD vector, TMatrixD covarience) :
    _dimension(vector.GetNrows()),
    _vector(vector),
    _covariance(covarience),
   _has_value(true) {
   if((_vector.GetNcols() != 1) ||
        (_covariance.GetNrows() != static_cast<int>(_dimension)) ||
        (_covariance.GetNcols() != static_cast<int>(_dimension)));
   {
   	throw "Inconsistancy between vector and covarience in Kalman::State::State()";
        
    }
   }//end has value

  State::State(const State& state) :
    _dimension(state._dimension),
    _vector(state._vector),
    _covariance(state._covariance),
    _has_value(state._has_value) {
   }

  //Copy over:
  State& State::operator=(const State& state) {
    this->_dimension = state._dimension;
    this->_vector = state._vector;
    this->_covariance = state._covariance;
    this->_has_value = state._has_value;
     
    return *this;
  }

  void State::SetVector(TMatrixD vector) {
    if((vector.GetNrows() != static_cast<int>(_dimension)) || (vector.GetNcols() != 1)){
    throw "ERROR: inconsistancy vector rows != dim";
    } 
    _vector = vector;
    _has_value = true;
  }

  void State::SetCovariance(TMatrixD covarience) {
    if((covarience.GetNrows() != static_cast<int>(_dimension)) ||
        (covarience.GetNcols() != static_cast<int>(_dimension)))
    { 
	throw "Covariance matrix has wrong dimensions in Kalman::State::SetCovariance()";
    }
    _covariance = covarience;
    _has_value = true;
  }

/*--------------------HIT MEMBER FUNCTIONS----------------------//
//--------------------------------------------------------------*/

  Hit::Hit(unsigned int dim, unsigned int data_dim, double position, int id) :
    _predicted(dim),
    _filtered(dim),
    _smoothed(dim),
    _data(data_dim),
    _hit_id(id),
    _position(position) {
    }

  Hit::Hit(unsigned int dim, double position, State data, int id) :
    _predicted(dim),
    _filtered(dim),
    _smoothed(dim),
    _data(data),
    _hit_id(id),
    _position(position) {
    }

  Hit::Hit(State predicted, State filtered, State smoothed, State data, double position, int id) :
    _predicted(predicted),
    _filtered(filtered),
    _smoothed(smoothed),
    _data(data),
    _hit_id(id),
    _position(position) {
    }

  Hit::Hit(const Hit& hit) :
    _predicted(hit._predicted),
    _filtered(hit._filtered),
    _smoothed(hit._smoothed),
    _data(hit._data),
    _hit_id(hit._hit_id),
    _position(hit._position) {
  }

  Hit& Hit::operator=(const Hit& hit) {
    if ((this->GetDimension()) != hit.GetDimension()) {
      throw "TrackPoint has wrong dimensions Kalman::Hit::operator=()";
    }
    this->_predicted = hit._predicted;
    this->_filtered = hit._filtered;
    this->_smoothed = hit._smoothed;
    this->_data = hit._data;

    this->_hit_id = hit._hit_id;
    this->_position = hit._position;

    return *this;
  }

  Hit& Hit::copy(Hit hit) { 
    if (this->GetDimension() != hit.GetDimension()) {
      throw "Hit has wrong dimensions Kalman::Hit::operator=()";
    }
    this->_predicted = hit._predicted;
    this->_filtered = hit._filtered;
    this->_smoothed = hit._smoothed;
    this->_data = hit._data;

    return *this;
  }


  State& Hit::GetFiltered() {
    if (_filtered.HasValue()) {
      return _filtered;
    } else {
      return _predicted;
    }
  }

  const State& Hit::GetFiltered() const {
    if (_filtered.HasValue()) {
      return _filtered;
    } else {
      return _predicted;
    }
  }

  State& Hit::GetSmoothed() {
    if (_smoothed.HasValue()) {
      return _smoothed;
    } else if (_filtered.HasValue()) {
      return _filtered;
    } else {
      return _predicted;
    }
  }

  const State& Hit::GetSmoothed() const {
    if (_smoothed.HasValue()) {
      return _smoothed;
    } else if (_filtered.HasValue()) {
      return _filtered;
    } else {
      return _predicted;
    }
  }

/*------------------KAL TRACK MEMBER FUNCTIONS------------//
//----------------------------------------------------------*/

  KalTrack::KalTrack(unsigned int dim, unsigned int length) :
    _dimension(dim) {
    for (unsigned int i = 0; i < length; ++i) {
      _track_vector.push_back(Hit(_dimension, _dimension, 0.0));
    }
  }

  void KalTrack::SetHit(unsigned int index, Hit hit) {
    if (hit.GetDimension() != this->_dimension) {
      throw "TrackPoint has wrong dimensions Kalman::KalTrack::SetHit()";
    }
    if (index >= _track_vector.size()) {
      throw  "Index out of bounds of track Kalman::KalTrack::SetHitPoint()";
    }
    _track_vector[index] = hit;
  }

  void KalTrack::AddHit(Hit hit) {
    if (hit.GetDimension() != this->_dimension) {
      throw "TrackPoint has wrong dimensions Kalman::KalTrack::AddHit()";
    }
    _track_vector.push_back(hit);
  }

  void KalTrack::DeleteHit(unsigned int index) {
    if (index >= _track_vector.size()) {
      throw  "Index out of bounds of track Kalman::KalTrack::SetState()";
    }
    _track_vector.erase(_track_vector.begin()+index);
  }


  void KalTrack::ResetTrack(const KalTrack& similar_track) {
    _track_vector.clear();

    for (unsigned int i = 0; i < similar_track.GetLength(); ++i) {
      Hit new_hit(_dimension, similar_track[i].GetData().GetDimension(),
                                                                   similar_track[i].GetPosition());
      new_hit.SetId(similar_track[i].GetId());
      _track_vector.push_back(new_hit);
    }
  }

}//End Kalman
} // namespace Mu2e


