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
          //line->set_cov(V_p);
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
	  	//(_x[i]*P[1][0]+P[0][0])*
                //TMatrixD D(V_m - A*V_p*At);
		//double fit_error = sqrt(D[i][i]*D[i][i]);
		//double hit_error = V_m[i][i];

		//double err_m_term = (_x[i]*P[1][0])*(_x[i]*P[1][0])*(((V_p[1][1])/(P[1][0]*P[1][0])));
	        //double err_m_term_1=(_x[i]*P[1][0])*(_x[i]*P[1][0])*(V_m[i][i]/_x[i])*(V_m[i][i]/_x[i]);
		//double cross_term = 2*sqrt((((V_p[1][1]))/(P[1][0]*P[1][0]))+(((V_m[i][i]/_x[i])*(V_m[i][i]/_x[i]))))*(_x[i]*P[1][0])*sqrt(V_p[0][0]);
		//double c_term = V_p[0][0];
		 
	 	//double fit_error =sqrt((V_p[0][0]) + (_x[i]*P[1][0])*(_x[i]*P[1][0])*(((V_p[1][1])/(P[1][0]*P[1][0]))+((V_m[i][i]/_x[i])*(V_m[i][i]/_x[i]))));
		//double residual_error = sqrt(fit_error*fit_error + hit_error*hit_error);
	        //double residual_error = sqrt(err_m_term+err_m_term_1+cross_term+c_term + hit_error*hit_error);
		
		double residual_error = 1/sqrt(V_m[j][j]);
		line->set_fit_residuals(C[j][0]);
		std::cout<<" res err  = " <<C[j][0]<< " / "<<residual_error<< " = "<<C[j][0]/residual_error<<std::endl;
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
	  line->set_m_1(P[2][0]);
	  //Extract Covarience Parameters:
          
	  line->set_c_0_err(sqrt(V_p[0][0]));
	  line->set_m_0_err(sqrt(V_p[1][1]));
          line->set_m_1_err(sqrt(V_p[2][2]));
      
	  //Extract Chi 2 Info:
	  line->set_chisq(result[0][0]);
	  line->set_chisq_dof(result[0][0] / n_points);

	 double sum_of_squares = 0;
         
	 for(int i=0; i< n_points; i++){
		sum_of_squares += C[i][0]*C[i][0];
	
          }
	  line->set_fit_residual_errors(sqrt(sum_of_squares));

	  for(int i=0; i< n_points; i++){
	  	
      
		
                //TMatrixD D(V_m - A*V_p*At);
		//double fit_error = sqrt(D[i][i]*D[i][i]);
		
		
		line->set_fit_residuals(C[i][0]);//C[i][0]/sqrt(sum_of_squares)); //scale to number of CHs in track? /n works
		
		
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





/*------------full fit for constrained chi-squared fit ---------*/
void minimise_chisq( int _nCoeffs, const std::vector<std::vector<double>> measurement, std::vector<double> &_err, StraightTrack* line, TMatrixD& covariance) { 
        
          
	





} // ~linear_fit






/*
void circle_fit(const double R_res_cut,  std::vector<double> &_x,  std::vector<double> &_y, TF2 &circle, CircleFit &Circle, TMatrixD& covariance){
// Number of measurements:
  int n_points = static_cast<int>(_x.size()); 
//Functional form now has 3 columns:
  TMatrixD A(n_points, 3); 
// Covariance matrix of measurements:                     
  TMatrixD C(n_points, n_points);  
//Kappa of circle:           
  TMatrixD K(n_points, 1); 
  for (int i = 0; i < static_cast<int>(_x.size()); ++i) {  
	//Get x and y position of hits
  	double x_i = x[i];
  	double y_i = y[i];
	//Do the maths:
  	A[i][0] = (x_i * x_i) + (y_i * y_i);
  	A[i][1] = x_i;
  	A[i][2] = y_i;
  	float dif = 1.;//WHAT IS THIS?
  	C[i][i] = (dif*dif);
  	K[i][0] = 1.;
        }

// Build the least squares fit:
  double* det = NULL; 
// Invert the measurement covariance matrix:            
  C.Invert(det);  
//Transpose A (but A still remains as used later):               
  TMatrixD At(A);                  
  At.T(); 
// The covariance matrix of the parameters of model (inv)                         
  TMatrixD V_p(At * C * A);      
  V_p.Invert(det);  
// The least sqaures estimate of the parameters                           
  TMatrixD P(V_p * At * C * K);  

// Extract the fit parameters
  double alpha, beta, gamma;
  alpha = P[0][0];
  beta = P[1][0];
  gamma = P[2][0];

// Convert the linear parameters into the circle center and radius
  double x0, y0, R;
  x0 = (-1*beta) / (2 * alpha);
  y0 = (-1*gamma) / (2 * alpha);
  if (((4 * alpha) + (beta * beta) + (gamma * gamma)) < 0){
    R = 0;
  }
  else{
    R = sqrt((4 * alpha) + (beta * beta) + (gamma * gamma)) / (2 * alpha);
  }
// Transform the covariance matrix to the same basis
  TMatrixD jacobian(3, 3);
  jacobian(0, 0) = beta / (2.0*alpha*alpha);
  jacobian(0, 1) = -1.0 / (2.0*alpha);
  jacobian(1, 0) = gamma / (2.0*alpha*alpha);
  jacobian(1, 2) = -1.0 / (2.0*alpha);
  jacobian(2, 0) = (-1.0/(2.0*alpha)) * (((beta*beta + gamma*gamma) / (2.0*alpha)) + 1) /sqrt(((beta*beta + gamma*gamma) / 4.0) + alpha);
  jacobian(2, 1) = (beta/(4.0*alpha*alpha)) /sqrt(((beta*beta + gamma*gamma)/(4.0*alpha*alpha)) + (1.0/alpha));
  jacobian(2, 2) = (gamma/(4.0*alpha*alpha)) /sqrt(((beta*beta + gamma*gamma)/(4.0*alpha*alpha)) + (1.0/alpha));
  TMatrixD jacobianT(3, 3);
  jacobianT.Transpose(jacobian);

  covariance = jacobian * V_p * jacobianT;

  R = fabs(R);

  if (R > R_res_cut){ //R_res_cut = tracker radius (look up)
     return false; 
  }

//Set circle parameters:
  circle.set_x0(x0);
  circle.set_y0(y0);
  circle.set_R(R);
  circle.set_alpha(alpha);
  circle.set_beta(beta);
  circle.set_gamma(gamma);

  // Calculate the fit chisq
  TMatrixD Z(K - (A * P));
  TMatrixD Zt(Z);
  Zt.T();
  TMatrixD result(Zt * C * Z);
  double chi2 = result[0][0];
  //circle.set_chisq(chi2);       // Left unreduced (not dividing by NDF)
 
}//end circle fit
*/
}//ends namespace
