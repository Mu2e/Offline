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
	  
          
	} // ~linear_fit

void xyz_fit( const std::vector<double> &_x, const std::vector<double> &_y,const std::vector<double> &_z,
		        const std::vector<double> &_err, StraightTrack* line, TMatrixD& covariance) { 
	  
	  
	  // Set up the matrices
          // Number of hits (N):
	  int n_points = static_cast<int>(_x.size());  
       
	  // Represents the functional form A is (N x 3) matrix:
	  TMatrixD A(n_points, 3);    
	  // Covariance matrix of measurements ( N x N) matrix :                
	  TMatrixD V_m(n_points, n_points);    
	  //Set up Y as a (N x 1) matrix :        
	  TMatrixD Y(n_points, 1);                     
          //Fill A and Z:
	  for (int i = 0; i < static_cast<int>(_x.size()); ++i) {
            //Fill first column with 1's:
	    A[i][0] = 1;
	    //Fill second column with x_i:
	    A[i][1] = _x[i];
	    //Fill third column with y_i:
            A[i][2] = _z[i];
           
            //Get Coverience matrix as hit-error:
	    V_m[i][i] = (_err[i] * _err[i]);
            
	    //Fill Y column vector:
	    Y[i][0] = _y[i];
	    
            
	  }

	  // Perform the inversions and multiplications which make up the least squares fit
	  double* det = NULL; 
	  // Invert in place - ( N x N) matrix                  
	  V_m.Invert(det);                      
	  
           // Copy A to At- retain A
	  TMatrixD At(A);                      
	  
	  // Transpose At (leaving A unchanged) - At is ( 3 x N ) matrix:
	  At.T();
	  // The covariance matrix of the parameters of model(inv)  - Vp is (3 x 3) matrix                              
	  TMatrixD V_p(At * V_m * A);           
         
          //Invert V_p:
          try {
    	       V_p.Invert(det); 
    
  	   } catch (exception& exc) {
    
		    std::stringstream message;
		    message << "cannot fit due to singular matrix error on inversion!" << std::endl;
          }                     
	  
	  
	  covariance = V_p;
          // The least sqaures estimate of the parameters - P is a (3 x 1) matrix:
	  TMatrixD P(V_p * At * V_m * Y);       
          
	  // Extract the Fit Parameters
	  line->set_c_0(P[0][0]);        
	  line->set_m_0(P[1][0]);  
	  line->set_m_1(P[2][0]);
	  //Extract Covarience Parameters:
          
	  line->set_c_0_err(sqrt(V_p[0][0]));
	  line->set_m_0_err(sqrt(V_p[1][1]));
	  std::cout<<line->get_m_0_err()<<std::endl;
          line->set_m_1_err(sqrt(V_p[2][2]));
          
	  // Calculate the fit chisq
	  TMatrixD C(Y - (A * P));
	  for(int i=0; i< n_points; i++){
	  	line->set_fit_residuals(C[i][0]);

		double fit_error = Y[i][0]*sqrt(((V_p[0][0])*(V_p[0][0]))+((V_p[1][1])*(V_p[1][1]))+((V_p[2][2])*(V_p[2][2])));
		
		line->set_fit_residual_errors(sqrt((_err[i]*_err[i]]) + (fit_error*fit_error)));
	   }
	  TMatrixD Ct(C);
	  Ct.T();
	  TMatrixD result(Ct * V_m * C);
	  line->set_chisq(result[0][0]);
	  line->set_chisq_dof(result[0][0] / n_points);
	  
         
          
	} // ~linear_fit
/*
void Plot_Fit(const std::vector<double> &_x, const std::vector<double> &_y, StraightTrack* line){
	TCanvas* c1=new TCanvas(); 
        TH2D* Event_Display = new TH2D("Event_Display", "Event_Display", 50,-1000., 1000., 50, -1000., 1000.);
        TGraph* hits = new TGraph(_x.size());
        TF1* fit = new TF1("fit", "[0]*x+[1]", -1000.,1000.);
        c1->Draw();
        
	for (int i = 0; i < static_cast<int>(_x.size()); ++i) {
            hits->SetPoint(i, _x[i], _y[i]);
        }
        
        
	fit->SetParameter(0,  line->get_m_0());
	fit->SetParameter(1,line->get_c_0());
        Event_Display->Draw("same");
	hits->SetMarkerSize(20);
        hits->Draw("Psame");
	fit->Draw("same");
	c1->Update();
        c1->SaveAs("new_event.root");
}
*/


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
