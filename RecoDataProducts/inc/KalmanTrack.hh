#ifndef KALMAN_TRACK_HH
#define KALMAN_TRACK_HH
//S. Middleton 
//Feb 2019
// Multipurpose Kalman Filtering routine for discrete kalman filtering of hits to update track variables
//C++:
#include <vector>
#include <algorithm>
//ROOT:
#include "TMath.h"
#include "TMatrixD.h"
//Mu2e
#include "RecoDataProducts/inc/KalmanTrack.hh"

namespace mu2e {
namespace Kalman {

  // Class and Function Declarations
  class State;
  class Hit;
  class KalTrack;


  
  /*---------------------------------State-------------------------------------//
  //                       Basic Kalman State object
  /---------------------------------------------------------------------------*/
  class State{

     public:
	//Constructors:
        explicit State(unsigned int dimension);
        State(TMatrixD vector, TMatrixD covariance);
        //Copy Constructor:
        State(const State& st);

        //Assignment:
        State& operator=(const State& st);

	//Accessors:
        TMatrixD& GetVector() { return _vector; }
        const TMatrixD& GetVector() const { return _vector; }
	TMatrixD& GetCovariance() { return _covariance; }
        const TMatrixD& GetCovariance() const { return _covariance; }
        unsigned int GetDimension() const { return _dimension; }
	//Modiffiers:
        void SetVector(TMatrixD vec);
	void SetCovariance(TMatrixD cov);

	//Destructor=:
	virtual ~State(){}

	//Has Value Boolean:
	bool HasValue() const { return _has_value; }

        void SetHasValue(bool has_value) { _has_value = has_value; }
     private:
        unsigned int _dimension;
	TMatrixD _vector;
        TMatrixD _covariance;
	
        bool _has_value;

   };//End State Class
  /*-----------------------------------Hit------------------------------------//
    Stores info of single hit. This includes its position, data, predicted, filtered and smoothed states
  //--------------------------------------------------------------------------*/
  class Hit{
     public:

      //Constructors:
      Hit(unsigned int dimension, unsigned int measurement_dim, double position,int hit_id = 0);
      Hit(unsigned int dimension, double position, State data, int id = 0);
      Hit(State predicted, State filtered, State smoothed, State data, double position, int hit_id = 0 );

      //Copy Constructor: 
      Hit(const Hit& hit);

      //Assignment Operator:
      Hit& operator=(const Hit& hit);

      //Copt Operator:
      Hit& copy(Hit hit);

      //Destructor
       virtual ~Hit() {}
    
      //Accessors:
      unsigned int GetDimension() const { return this->_predicted.GetDimension(); }
      double GetPosition() const { return _position; }
      int GetId() const { return _hit_id; }
      State& GetPredicted() { return _predicted; }
      const State& GetPredicted() const { return _predicted; }
      State& GetFiltered();
      const State& GetFiltered() const;
      State& GetSmoothed();
      const State& GetSmoothed() const;
      State& GetData() { return _data; }
      const State& GetData() const { return _data; }

      //Modiffiers:
      void SetPosition(double new_pos) { _position = new_pos; }
      void SetPredicted(State state) { _predicted = state; }
      void SetFiltered(State state) { _filtered = state; }
      void SetSmoothed(State state) { _smoothed = state; }
      void SetData(State state) { _data = state; }
      void SetId(int new_id) { _hit_id = new_id; }

      //Check for Data:
      bool HasData() const { return _data.HasValue(); }

    private:

      State _predicted;
      State _filtered;
      State _smoothed;
      State _data;
      int _hit_id;
      double _position;
    };// End Hit Class


   /*-----------------------------------Kalman Track--------------------------//
                        Stores a list of filtered/smoothed hits
   //------------------------------------------------------------------------*/
   class KalTrack{
    typedef std::vector< Hit > HitArray;

    public:
      //Constructors:
      KalTrack(unsigned int dimension, unsigned int length = 0);

      //Accessors:
      Hit GetHit(unsigned int index) const { return _track_vector[index]; }
      //reference to the state at the specified index:
      Hit& operator[] (unsigned int index) { return _track_vector[index]; }
      //Get a const reference to the state at the specified index:
      const Hit& operator[] (unsigned int index) const { return _track_vector[index]; }
      unsigned int GetLength() const { return _track_vector.size(); }
      unsigned int GetDimension() const { return _dimension; }

      //Modiffiers:
      void SetHit(unsigned int index, Hit hit);
      void AddHit(Hit tp);
      void DeleteHit(unsigned int index);
      void ResetTrack(const KalTrack& similar_track);


    private:
      unsigned int _dimension;
      HitArray _track_vector;
  };

}//End Mu2e
}//End Kalman

#endif
