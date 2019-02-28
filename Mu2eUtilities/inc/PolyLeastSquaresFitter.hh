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
#include "Mu2eUtilities/inc/PolyCoeff.hh"

using namespace mu2e;

namespace PolyLeastSquaresFitter {
	 
         void PolyLeastSquaresFitter(int numberOfInputVariables,TMatrixD& polynomialCoefficients);


          void Fit(std::vector< std::vector<double> >& points, std::vector< std::vector<double> >& values,int polynomialOrder, std::vector<double>& weights,TMatrixD point_errors);
   
	  std::vector<int> VecIndex(int index, int nInputVariables);
	  std::vector<int> PowIndex(int index, int nInputVariables);
	
	  void SetCoeffs(int pointDim, TMatrixD& coeff);
	  //TMatrixD GetCoefficientsMatrix() const ();
          void PrintCoefficientsMatrix();
          double*  MakePolynomialVector(const double* point, double* polyVector, int pointdim, int ncoeffs);
          int Factorial(int val); 
          int Combinations(int N, int R) ;
          int GetNCoeffs(int pointDim,int polyorder);
          
	
}//~ namespace


#endif

