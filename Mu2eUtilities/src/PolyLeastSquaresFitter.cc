//Fitting Routines for Poly Fitting
//Author: S Middleton, Jan 2019
//Mu2e:
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"
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


PolyLeastSquaresFitter::PolyLeastSquaresFitter(int numberOfInputVariables,
                                   TMatrixD& polynomialCoefficients) : point_dimension_(numberOfInputVariables), index_key_by_power_(), index_key_by_vector_(), coefficient_matrix_(polynomialCoefficients) {
  SetCoeffs(numberOfInputVariables, polynomialCoefficients);
}


std::vector<int> PolyLeastSquaresFitter::VecIndex(int index, int nInputVariables) {
    if (index == 0) return std::vector<int>();
    std::vector<int> indices(1, -1);
    nInputVariables--;
    for (int i = 0; i < index;   ++i) {
        if (indices.front() == nInputVariables) {
          indices = std::vector<int>(indices.size()+1, 0);
          indices.back()--;
        }
        if (indices.back()  == nInputVariables) {
            int j = indices.size()-1;
            while (indices[j] == nInputVariables) j--;
            for (int k = indices.size()-1; k >= j; k--) {
                indices[k] = indices[j]+1;
            }
        } else {
          indices.back()++;
        }
    }
    return indices;
}

std::vector<int> PolyLeastSquaresFitter::PowIndex(int index, int nInputVariables) {
    std::vector<int> powerIndex(nInputVariables, 0);
    std::vector<int> vectorIndex = VecIndex(index, nInputVariables);
    for (size_t i = 0; i < vectorIndex.size(); ++i) {
        powerIndex[vectorIndex[i]]++;
    }
    return powerIndex;
}

void PolyLeastSquaresFitter::SetCoeffs(int pointDim, TMatrixD& dummy_coefficent_matrix) {
 
  point_dimension_     = pointDim;
  coefficient_matrix_  = dummy_coefficent_matrix;
  int num_poly_coeff   = coefficient_matrix_.GetNcols();
  index_key_by_power_  = std::vector< std::vector<int> >();
  index_key_by_vector_ = std::vector< std::vector<int> >();
  polynomial_vector_ = std::vector<PolyCoeff>();
  //for (int i = 0; i < num_poly_coeff; ++i){
  //  index_key_by_power.push_back(PowIndex(i, pointDim));
  //}
  for (int i = 0; i < num_poly_coeff; ++i){
    index_key_by_vector_.push_back(VecIndex(i, pointDim));
  }
  for (int i = 0; i < dummy_coefficent_matrix.GetNrows(); ++i){
    for (int j = 0; j < num_poly_coeff; ++j){
      polynomial_vector_.push_back(
          PolyLeastSquaresFitter::PolyCoeff(index_key_by_vector_[j], i,dummy_coefficent_matrix[i+1][ j+1]));
     }
  }
}

double*  PolyLeastSquaresFitter::MakePolynomialVector(const double* point,
                                             double* polyVector)
    const {
  for (size_t i = 0; i < index_key_by_vector_.size(); ++i) {
    polyVector[i] = 1.;
    for (size_t j = 0; j < index_key_by_vector_[i].size(); ++j)
      polyVector[i] *= point[ index_key_by_vector_[i][j] ];
  }
    return polyVector;
}



PolyLeastSquaresFitter* PolyLeastSquaresFitter::Fit(const std::vector< std::vector<double> >& points,const std::vector< std::vector<double> >& values, int polynomialOrder,const std::vector<double>& weights,TMatrixD point_errors) {
  // Algorithm: We have F2 = sum_i ( f_k f_l) where f are polynomial terms;
  // FY = sum_i (f_)

  int pointDim =static_cast<int> (points[0].size());
  int valueDim =static_cast<int> (values[0].size());
  int nPoints  = static_cast<int> (points.size());
  int nCoeffs;

  if (polynomialOrder <= 0) nCoeffs = 1;

  if (polynomialOrder > 0){
	  double num = 1, dem_1 =1, dem_2 =1;
	  for (int i = 0; i < polynomialOrder; ++i) {
		num *= (pointDim+i);
		dem_1 *=(i+1);
		dem_2 *=(pointDim-1);
	  }
	  for (int i = 0; i < polynomialOrder; ++i) {
		nCoeffs += num/(dem_1*dem_2);
	    }
  }

  if (point_errors.GetNcols() == 0 &&
      point_errors.GetNrows()== 0) {
      point_errors(nCoeffs, nCoeffs);
  }
  
  // create a dummy TMatrixD: A
  TMatrixD dummy_coefficent_matrix(valueDim, nCoeffs);
  for (int i = 0; i < valueDim; ++i){
    for (int j = 0; j < nCoeffs; ++j){
      dummy_coefficent_matrix[i+1][j+1] = 1;
    }
  }
  //Initialze the fit allowing for setting of the init coeffs
  PolyLeastSquaresFitter* Poly_Fit = new PolyLeastSquaresFitter(pointDim, dummy_coefficent_matrix);

  // Create the design matrix A
  std::vector<double> point_poly_vector(nCoeffs, 0);
  TMatrixD A(nPoints, nCoeffs);  // design matrix

  for (int row = 0; row < nPoints; ++row) {
    Poly_Fit->MakePolynomialVector(&points[row][0],&point_poly_vector[0]); 
                                
    for (int column = 0; column < nCoeffs; ++column) {
      A[row+1][column+1] = point_poly_vector[column];
    }
  }

  // Create the value matrix
  TMatrixD Y(nPoints, valueDim);  // value matrix
  for (int row = 0; row < nPoints; ++row) {
    for (int column = 0; column < valueDim; ++column) {
      Y[row+1][column+1] = values[row][column];
    }
  }

  // Create the weight matrix (diagonal are the per-point/value weights)
  //std::vector<double>& weight_vector;
  TMatrixD W(nPoints, nPoints);
  //if (weights.size() > 0) {
    
  //  weight_vector =  weights;
  //}
  for (int index = 0; index < nPoints; ++index) {
    W[index+1][index+1] = weights[index];
  }
  //Define transpose of design matrix
  TMatrixD At = A.T();

  // F2 = A^T A, where A is the design matrix with linearly independent columns
  TMatrixD F2(nCoeffs, nCoeffs);
  F2 = At * W * A;
  for (int i = 1; i < nCoeffs+1; ++i)
    for (int j = 1; j < nCoeffs+1; ++j)
      F2[i][j] -= point_errors[i][j]*nPoints;

  // Fy = A^T Y, where A is the design matrix and Y is the value matrix
  TMatrixD Fy(nCoeffs, valueDim);
  Fy = At * W * Y;

  // F2^(-1) = (A^T A)^(-1), where A is the design matrix
  
  
  try {
    double* det = NULL;
    F2.Invert(det);
  } catch (exception& exc) {
    
    std::stringstream message;
    message << "Could not find least squares fit for data. Nested exception:"
            << std::endl;
    //throw(exception& (Exception::recoverable,
    //             message.str(),
    //             "Error in PolynomialLeastSquaresFit -  cannot invert!!!"));
  }

  // C = (A^T A)^(-1) A^T Y, where A is the design matrix and Y is the value matrix
  TMatrixT<double> coefficient_matrix = (F2 * Fy);
  TMatrixT<double> coefficient_matrix_transpose = coefficient_matrix.T();
  Poly_Fit->SetCoeffs(pointDim, coefficient_matrix_transpose);
  std::cout<<Poly_Fit<<std::endl;
  return Poly_Fit;
}

