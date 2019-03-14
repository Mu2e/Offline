#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
////S. Middleton, Feb 2019 - Cosmic track class, will store info of fit result in terms of m and c and errors.
#include "TMath.h"
#include "TMatrixD.h"
#include "DataProducts/inc/XYZVec.hh"
#include<vector>
#include<bitset>
/*

Cosmic track is parameterised by a paramteric equatrion r(t) were t is a vector of arbituary parameters.

e.g. r(t) = 	   x0	      a
		[  y0 ] + t*[ b ]    assuming 3D with 4 parameters then: 
		   z0         c	      
(first part is point on line, second is direction)
Assume all pass through vertical (assumption for cosmics ...not true)
r(t) = 	   x0	      a
	[  1 ] +  t*[ 0 ]    assuming 3D with 4 parameters then: 
	   z0 	      b     

our parameters are then x0, z0, a and b
*/
using namespace std;

namespace mu2e {
  
  class CosmicTrack{
  public:
   
    // Constructors
    CosmicTrack();
    
    CosmicTrack(double par_1, double par_2, double par_3, double par_4) : _par1(par_1), _par2(par_2), _par3(par_3), _par4(par_4){};

    CosmicTrack(std::vector<double> track_parameters) : _track_parameters(track_parameters){};

    
    // Destructor
    ~CosmicTrack();


    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);

    // Accessors:
    //Number or hits in track
    int get_N() const { return _Nhits;}
    //Chi^2 info
    double get_chisq() const { return _chisq; }
    double get_chisq_dof() const { return _chisq_dof; }
    //track parameter info
    std::vector<double> get_track_parameters() const { return _track_parameters; } 
    
    double get_parameter(unsigned para_ID){ 
	if(_track_parameters.size() >=para_ID){
		return _track_parameters.at(para_ID);
	}throw "Error: parameter list not long enough";
     }
     
    XYZVec get_track_equation(){return _track_equation;}
    XYZVec get_track_direction(){ return _track_direction;}
    //curvature and momentum info
    double GetSagitta(){return _Sagitta;}
    //residuals info
    std::vector<double> get_fit_residuals() const { return _fit_residuals; }
    std::vector<double> get_fit_residual_errors() const { return _fit_residual_errors; }
    //hit errors info
    std::vector<double> get_hit_errorsX() const{ return _hit_errorsX;}
    std::vector<double> get_hit_errorsY() const{ return _hit_errorsY;}
    std::vector<double> get_hit_errorsTotal() const{ return _hit_errorsTotal;}
    
    // Modiffiers:
    void clear();
    void set_N(int N){_Nhits = N;}
    void set_parameters(double par1,double par2, double par3, double par4){ 
	_par1=par1;
	_par2=par2;
	_par3=par3;
	_par4=par4;
        _track_parameters.push_back(par1);
        _track_parameters.push_back(par2);
	_track_parameters.push_back(par3);
	_track_parameters.push_back(par4);

	}
    void set_parameters(unsigned para_ID, double par){
    	if(_track_parameters.size() >= para_ID){
		_track_parameters.at(para_ID) = par;
	}throw "Error: Parameter list not long enough";
	 
	
    }
    void set_track_direction(XYZVec direction){_track_direction = direction;}
    void set_chisq(double chisq) { _chisq = chisq; }
    void set_mom(XYZVec mom){_track_mommentum=mom;} 

    void set_fit_residuals(double residual) { _fit_residuals.push_back(residual); }
    void set_fit_residual_errors(double residual_err) { _fit_residual_errors.push_back(residual_err); }
    
    void set_hit_errorsX(double hit_errorX){_hit_errorsX.push_back(hit_errorX);}
    void set_hit_errorsY(double hit_errorY){_hit_errorsY.push_back(hit_errorY);}
    void set_hit_errorsTotal(double hit_errorTotal){_hit_errorsTotal.push_back(hit_errorTotal);}
    
    void set_cov(int r, int c, double element){ _cov_mat[r][c] = element;}
    
  private:
    
    int _Nhits;

    double _par1; //a0
    double _par2; //a1
    double _par3; //b0
    double _par4; //b1

    std::vector<double> _track_parameters; //FIXME NEED TO BE 4D vector.....

    double _Sagitta;//TODO
    XYZVec _track_mommentum;//TODO
    
    XYZVec _track_equation;//r(t) expression
    XYZVec _track_direction;//the "gradient" term
    XYZVec _track_point0;//the "starting point" in fit line
    
  
    double _chisq;
    double _chisq_dof;

    std::vector<double> _fit_residuals;
    std::vector<double> _fit_residual_errors;
    std::vector<double> _hit_errorsX;
    std::vector<double> _hit_errorsY;
    std::vector<double> _hit_errorsTotal;
    
    std::vector<std::vector<double> > _cov_mat;

    bool set_chosenOne;
   
  };
}

#endif

