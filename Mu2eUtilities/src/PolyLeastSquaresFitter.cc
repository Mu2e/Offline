//Fitting Routines for Poly Fitting
//Author: S Middleton, Jan 2019
//Mu2e:
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"
#include "Mu2eUtilities/inc/PolyCoeff.hh"
#include "Mu2eUtilities/inc/PolyLeastSquaresFitter.hh"


//c++
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <exception>
//ROOT
#include "TMatrixD.h"
using namespace std;
using namespace mu2e;

namespace PolyLeastSquaresFitter{
//----------Calls to Initializes Fit Coefficents -----------------//


void PolyLeastSquaresFitter(int numberOfInputVariables,
                                   TMatrixD& polynomialCoefficients)  {
  
  SetCoeffs(numberOfInputVariables, polynomialCoefficients);
}


std::vector<int> VecIndex(int index, int nInputVariables) {
    //if index = 0, retruns empty:
    if (index == 0) return std::vector<int>();
    ///intiate indeices:
    std::vector<int> indices(1, -1); //1 input of value -1 i.e starting at -1
    nInputVariables--;//predecreement
    for (int i = 0; i < index;   ++i) {
	//if first index == nPoints
        if (indices.front() == nInputVariables) {
          indices = std::vector<int>(indices.size()+1, 0);//add on one
          indices.back()--;//take one back from end
        }
	//if last ==nPoints
        if (indices.back()  == nInputVariables) {
            int j = indices.size()-1;
            while (indices[j] == nInputVariables) j--;
            for (int k = indices.size()-1; k >= j; k--) {
                indices[k] = indices[j]+1;
            }
        } else {
          indices.back()++;//++1 onto last element
        }
    }
    //returns list of indicees
    return indices;
}

// ----------- Set Coeff Matrix--------//
void SetCoeffs(int pointDim, TMatrixD& dummy_coefficent_matrix) {
 
  int num_poly_coeff   = dummy_coefficent_matrix.GetNcols();
  
   std::vector< std::vector<int> > index_key_by_vector_ ;
   std::vector<PolyCoeff> polynomial_vector_ ;

  for (int i = 0; i < num_poly_coeff; ++i){
    index_key_by_vector_.push_back(VecIndex(i, pointDim));
  }
  for (int i = 0; i < dummy_coefficent_matrix.GetNrows(); ++i){
    for (int j = 0; j < num_poly_coeff; ++j){
      polynomial_vector_.push_back(
          PolyCoeff(index_key_by_vector_[j], i,dummy_coefficent_matrix[i][ j]));
          std::cout<<"end of set co"<<i<<j<<dummy_coefficent_matrix[i][ j]<<std::endl;
     }
  }
}



//-----Print Test-------//
void PrintCoefficientsMatrix(PolyCoeff coefficient_matrix_) {
  
}

TMatrixD GetCoefficientsAsMatrix(TMatrixD coefficient_matrix){
  return coefficient_matrix;
}


double*  MakePolynomialVector(const double* point,double* polyVector, int pointDim,int nCoeffs){
   std::vector< std::vector<int> >index_key_by_vector_ ;
  
  for (int i = 0; i < nCoeffs; ++i){
    index_key_by_vector_.push_back(VecIndex(i, pointDim));
  }

  for (size_t i = 0; i < index_key_by_vector_.size(); ++i) {//size = number of coeffs
    polyVector[i] = 1.; //rows =1
    for (size_t j = 0; j < index_key_by_vector_[i].size(); ++j)
      //*
      polyVector[i] *= point[ index_key_by_vector_[i][j] ];
  }
    return polyVector;
}

//n,m = m = n!/(m!(n-m)!):

int Factorial(int val) 
 { 
	    std::cout<<"Fact"<<std::endl;
	    int Result = 1; 
	    for(int i = 1; i <= val; i++) 
	    { 
		 Result *= i; 
	      } 
	    return Result;
} 

int Combinations(int N, int R) //ok for int-->long? 
{ 	    
	    std::cout<<"Com"<<std::endl;
            return (Factorial(N)) / ((Factorial(N-R)) * Factorial(R)); 
} 

int GetNCoeffs(int pointDim,int polyorder) {
    std::cout<<"Getting Coeffs..."<<std::endl;
    if (polyorder <= 0) {
	return 1;
    }
    int n = Combinations(pointDim, polyorder+1);
    return n+1;
    
}

void Fit(std::vector< std::vector<double> >& points,std::vector< std::vector<double> >& values, int polynomialOrder, std::vector<double>& weights,TMatrixD point_errors) {
  // Algorithm: We have F2 = sum_i ( f_k f_l) where f are polynomial terms;
  // FY = sum_i (f_)
  std::cout<<" IN Poly Fit ...."<<std::endl;
  int pointDim =static_cast<int> (points[0].size());
  int valueDim =static_cast<int> (values[0].size());
  int nPoints  = static_cast<int> (points.size());
  int nCoeffs =GetNCoeffs(pointDim, polynomialOrder); 
  std::cout<<"nceofs"<< nCoeffs<<std::endl;

  if (point_errors.GetNcols() == 0 &&
      point_errors.GetNrows()== 0) {
      point_errors(nCoeffs, nCoeffs);
  }
  
  // STEP 1 : create a dummy TMatrixD - fill with "1s": A (Nv x NCoeff):
  TMatrixD dummy_coefficent_matrix(valueDim, nCoeffs);
  for (int i = 0; i < valueDim; ++i){
    for (int j = 0; j < nCoeffs; ++j){
      std::cout<<"in dummy"<<i<<" "<<j<<std::endl;
      dummy_coefficent_matrix[i][j] = 1;
    }
  }
  std::cout<<" Got Dummy..."<<std::endl;
  //Step 2: Initialze the fit allowing for setting of the init coeffs
  PolyLeastSquaresFitter(pointDim, dummy_coefficent_matrix);

  //Step 3: Create the design matrix A NP x NCoeff
  std::vector<double> point_poly_vector(nCoeffs, 0);
  TMatrixD A(nPoints, nCoeffs);  // design matrix
  std::cout<<" Got Design..."<<std::endl;
  //Step 4: For each point make a polynomical vector and fill A
  for (int row = 0; row < nPoints; ++row) {
    MakePolynomialVector(&points[row][0],&point_poly_vector[0], pointDim, nCoeffs); 
    std::cout<<" Made PolyVec..."<<std::endl;                             
    for (int column = 0; column < nCoeffs; ++column) {
      //Fill A will the values from poly vector /TODO
      A[row][column] = point_poly_vector[column];
    }
  }

  //Step 5: Create the value matrix
  TMatrixD Y(nPoints, valueDim);  // value matrix
  for (int row = 0; row < nPoints; ++row) {
    for (int column = 0; column < valueDim; ++column) {
      //For each point add in value Y has number of columns = value dimensions.
      Y[row][column] = values[row][column];
    }
  }
  std::cout<<" Filled Value..."<<std::endl;
  //Step 6: Create the weight matrix (diagonal are the per-point/value weights), W=V_m i.e covariences in the measurements
  //std::vector<double>& weight_vector(nPoints);
  TMatrixD W(nPoints, nPoints);
  
  for (int index = 0; index < nPoints; ++index) {
    W[index][index] = weights[index];
  }
  //Define transpose of design matrix
  // Copy A to At- retain A
  TMatrixD At(A);                      
  // Transpose At (leaving A unchanged) - At is ( _nCoeffs x NP) matrix:
  At.T();

  // Step 7: Get F2 = A^T A, where A is the design matrix with linearly independent columns (F2=V_p) (Nc x Nc)
  std::cout<<"get f2"<<std::endl;
  TMatrixD F2(At * W * A);
  for (int i = 1; i < nCoeffs; ++i)
    for (int j = 1; j < nCoeffs; ++j)
      F2[i][j] -= point_errors[i][j]*nPoints;
  std::cout<<"get fy"<<std::endl;
  // Step 8:  Get Fy = A^T Y, where A is the design matrix and Y is the value matrix - Nc x NV
  TMatrixD Fy(At * W * Y);
  std::cout<<"inverting"<<std::endl;
  // Step 9 :  get coefficents F2^(-1) = (A^T A)^(-1), where A is the design matrix
  try {
    double* det = NULL;
    F2.Invert(det);
  } catch (exception& exc) {
    
    std::stringstream message;
    message << "Could not find least squares fit for data. Nested exception:"
            << std::endl;
   
  }

  // C = (A^T A)^(-1) A^T Y, where A is the design matrix and Y is the value matrix
  TMatrixT<double> coefficient_matrix = (F2 * Fy);
  TMatrixT<double> coefficient_matrix_transpose = coefficient_matrix.T();
  //Step 10: Reset the Coefficients in fitter:
  SetCoeffs(pointDim, coefficient_matrix_transpose);
  
  
}




}//end namespace PLS


