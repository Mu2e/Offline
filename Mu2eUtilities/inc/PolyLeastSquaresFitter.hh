
#ifndef _MU2E_UTILITIES_POLYLEASTSQUARESFITTER_HH
#define _MU2E_UTILITIES_POLYLEASTSQUARESFITTER_HH
// This Class is to enable the 3x3 fitting using matrices.
// Author: S. Middleton
// Date: Jan 2019

//c++
#include <vector>
//ROOT
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "Math/VectorUtil.h"
#include "TMatrixD.h"

namespace mu2e{

class PolyLeastSquaresFitter {
	 public:
	  // forward declaration of embedded class
	  class PolyCoeff;
	  //Constructors:
          PolyLeastSquaresFitter();

          PolyLeastSquaresFitter(int numberOfInputVariables,
                                   TMatrixD& polynomialCoefficients);
       
          //Destructor
          ~PolyLeastSquaresFitter();

          PolyLeastSquaresFitter* Fit(const std::vector< std::vector<double> >& points, const std::vector< std::vector<double> >& values,int polynomialOrder,const std::vector<double>& weights,TMatrixD point_errors);
   
        

	  std::vector<int> VecIndex(int index, int nInputVariables);
	  std::vector<int> PowIndex(int index, int nInputVariables);
	
	  void SetCoeffs(int pointDim, TMatrixD& coeff);
      
          double*  MakePolynomialVector(const double* point,
                                             double* polyVector)const;
          class PolyCoeff{
    		public:
                //Constructor
      		PolyCoeff(std::vector<int> inVariablesByVector,
                                            int outVariable, double coefficient){}
                
           };//end embedded class
	private:
	  int                                point_dimension_;
	  std::vector< std::vector<int> >    index_key_by_power_; // std::vector<int>[i_1] = j_1
	  std::vector< std::vector<int> >    index_key_by_vector_; // std::vector<int>[j_1] = i_1
	  TMatrixD&                    coefficient_matrix_;
	  std::vector<PolyCoeff> polynomial_vector_;
	  
     //};
};//end polyLLSFitter class
}//end namespace mu2e

#endif

