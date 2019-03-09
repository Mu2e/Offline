//Least Squares Fitting Routine
//Author: S Middleton, Dec 2018
//Mu2e:
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
#include "Mu2eUtilities/inc/LeastSquaresFitter.hh"

//c++
#include <vector>
#include <iostream>
#include <sstream>
#include <exception>
#include <numeric>
using namespace std;
using namespace mu2e;
namespace LeastSquaresFitter {

/*---------------For 2D X-Y type fit use:--------------------*/
	
	void xy_fit( const std::vector<double> &_x, const std::vector<double> &_y,
		        const std::vector<double> &_y_err, StraightTrack* line, TMatrixD& covariance) { 
	 
          
	  // Set up the matrices:

	  int n_points = static_cast<int>(_x.size());  // Number of measurements
	  // Represents the functional form: (N x 2) matrix
	  TMatrixD A(n_points, 2);    
	  // Covariance matrix of measurements: (N x N) matrix                 
	  TMatrixD V_m(n_points, n_points);   
	  // Measurements:( N x 1) matrix          
	  TMatrixD Y(n_points, 1);                     
          
	  for (int i = 0; i < static_cast<int>(_x.size()); ++i) {
            //fill first column with 1's:
	    A[i][0] = 1;
	    //fill second column with x co-ordinate of hit:
	    A[i][1] = _x[i];
	    // fill covarience matrix with "hit errors^2":
	    V_m[i][i] = (_y_err[i] * _y_err[i]);
	    //fill "column vector" y:
	    Y[i][0] = _y[i];
            
	  }

	  // Perform the inversions and multiplications which make up the least squares fit

	  double* det = NULL;                   // To hold the determinant
	  // Invert Covarience Matrix in place: (N x N) matrix still.
	  V_m.Invert(det);  
	  // Copy A to At                    
	  TMatrixD At(A);
           // Transpose At (leaving A unchanged). At is a (2 x N) matrix now:                       
	  At.T();   
          // The covariance matrix of the parameters of model i.e. m and c : (2 x 2) matrix:                     
	  TMatrixD V_p(At * V_m * A);
	  // Invert in place           
	  V_p.Invert(det);                      
	  covariance = V_p;
         
	  //  Find the least sqaures estimate of the parameters, P is (2 x 1) matix of c and m:
	  TMatrixD P(V_p * At * V_m * Y);       
          
	  // Extract the fit parameters and set as line parameters to allow plotting:
	  line->set_c_0(P[0][0]);
	  line->set_m_0(P[1][0]);
	  line->set_c_0_err(sqrt(V_p[0][0]));
	  line->set_m_0_err(sqrt(V_p[1][1]));

	  //Get Chi2 and Chi2/N and set:
         
	  //Calculate C an (Nx1) matrix:
	  TMatrixD C(Y - (A * P));
	  //Copy and keep:
	  TMatrixD Ct(C);
	  //make a transpose- Ct is (1 x N) matrix
	  Ct.T();
          //Get the result matrix - this is a (1 x 1) matrix i.e. it is just the chisq:
	  TMatrixD result(Ct * V_m * C);
	  //Set line parameters:
	  line->set_chisq(result[0][0]);
	  line->set_chisq_dof(result[0][0] / n_points);
	  

	  for(int j=0; j< static_cast<int>(_x.size()); j++){
	
		double residual_error = 1/sqrt(V_m[j][j]);
		line->set_fit_residuals(C[j][0]);
		
		line->set_fit_residual_errors(residual_error);
		
	   }
          
	} // ~linear_fit

/*----------For simple Chis-Squared Calc for estimator of y use:--------------*/

void xyz_fit( int _nCoeffs, const std::vector<double> &_measured_i, const std::vector<double> &_measured_j,const std::vector<double> &_measured_k,const std::vector<double> &_err, StraightTrack* line, TMatrixD& covariance) { 
	  
	  
	  // Set up the matrices
          // Number of hits (N):
	  int n_points = static_cast<int>(_measured_i.size());
	  line->set_N(n_points);  
          int n_values = 1; //here maping to one dependent varibale (e.g. y)

	  // Represents the functional form A (design matrix) is (NP x NCoeff) matrix:
	  TMatrixD A(n_points, _nCoeffs);    
	  // Covariance matrix of measurements ( NP x NP) matrix :                
	  TMatrixD V_m(n_points, n_points);    
	  //Set up Y (value matrix) as a (NP x NV) matrix, Y is the predicted value i.e. the one you are estimating :        
	  TMatrixD Y(n_points, n_values); 
                    
          //Fill A and Z:
	  for (int i = 0; i < static_cast<int>(_measured_i.size()); ++i) {
            //Fill first column with 1's:
	    A[i][0] = 1;
	    //Fill second column with x_i:
	    A[i][1] = _measured_i[i];
	    //Fill third column with y_i:
            A[i][2] = _measured_k[i];
           
            //Get Coverience matrix on measurement from calculated hit-error:
	    V_m[i][i] = (_err[i] * _err[i]);
            
	    //Fill Value Matrix:
	    Y[i][0] = _measured_j[i];
	    
            
	  }

	  // Perform the inversions and multiplications which make up the least squares fit
	  double* det = NULL; 
	  // Invert in place - ( NP x NP) matrix                  
	  V_m.Invert(det);                      
	  
           // Copy A to At- retain A
	  TMatrixD At(A);                      
	  
	  // Transpose At (leaving A unchanged) - At is ( _nCoeffs x NP) matrix:
	  At.T();
	  // The covariance matrix of the parameters of model(inv)  - Vp is (_nCoeffs x _nCoeffs) matrix F2:                         
	  TMatrixD V_p(At * V_m * A);           
         
          //Invert V_p:
          try {
    	       V_p.Invert(det); 
    
  	   } catch (exception& exc) {
    
		    std::stringstream message;
		    message << "cannot fit due to singular matrix error on inversion!" << std::endl;
          }                     
	  //FY is a NCoeff x NV matrix:
	  TMatrixD FY(_nCoeffs,n_values);
	  FY = At * V_m * Y;

	  covariance = V_p;
          // The least sqaures estimate of the parameters - P is a (_nCoeffs x NV) matrix:
	  TMatrixD P(V_p * FY );       
          
          
	  //Calculate the residuals:
	  TMatrixD C(Y - (A * P));

          // Calculate the fit chisq
	  TMatrixD Ct(C);
	  Ct.T();
	  TMatrixD result(Ct * V_m * C);

          
	  // Extract the Fit Parameters
	  line->set_c_0(P[0][0]);        
	  line->set_m_0(P[1][0]);  
	  //line->set_m_1(P[2][0]);
	  //Extract Covarience Parameters:
          
	  line->set_c_0_err(sqrt(V_p[0][0]));
	  line->set_m_0_err(sqrt(V_p[1][1]));
          //line->set_m_1_err(sqrt(V_p[2][2]));
      
	  //Extract Chi 2 Info:
	  line->set_chisq(result[0][0]);
	  line->set_chisq_dof(result[0][0] / n_points);

	  for(int i=0; i< n_points; i++){
	  	
		line->set_fit_residuals(C[i][0]);
		
	   }

          
	  
         
          
	} // ~linear_fit

std::vector <double> Means(std::vector<std::vector<double> > values, std::vector<double> weights) {
  if (values.size() < 1) {
    throw "No input values" ;
  }
  if (values[0].size() < 1) {
    throw"Dimension < 1" ;
  }
  if (weights.size() != values.size()) {
    weights = std::vector<double>(values.size(), 1.);
  }
  size_t npoints     = values.size();
  size_t dim         = values[0].size();
  std::vector<double>    means(dim, 0);
  double             totalWeight = 0;
  for (size_t x = 0; x < npoints; x++) {
    totalWeight += weights[x];
  }
  std::vector<double> normalized_weights;
  for (size_t x = 0; x < npoints; x++) {
    normalized_weights.push_back(weights[x] / totalWeight);
  }

  double mean;
  for (size_t i = 0; i < dim; ++i) {
    mean = 0.;
    for (size_t x = 0; x < npoints; ++x) {
      mean += values[x][i] * normalized_weights[x];
    }
    means[i+1] = mean;
  }

  return means;
}

//For 1D list of values:
double Mean(std::vector<double>  values, std::vector<double> weights) {
  
  if (weights.size() != values.size()) {
    weights = std::vector<double>(values.size(), 1.);
  }
  size_t npoints     = values.size();
 
  double  totalWeight = 0;
  for (size_t x = 0; x < npoints; x++) {
    totalWeight += weights[x];
  }
  std::vector<double> normalized_weights;
  for (size_t x = 0; x < npoints; x++) {
    normalized_weights.push_back(weights[x] / totalWeight);
  }

  double mean;
  
   mean = 0.;
   for (size_t x = 0; x < npoints; ++x) {
   mean += values[x] * normalized_weights[x];
   }
   

  return mean;
}

}//ends namespace
