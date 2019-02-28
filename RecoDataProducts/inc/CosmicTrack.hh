#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
//S. Middleton - Cosmic track class, will store info of fit result in terms of m and c and errors.
#include "TMath.h"
#include "TMatrixD.h"
#include "RecoDataProducts/inc/XYZVec.hh"
#include<vector>
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

    CosmicTrack(XYZVec track_parameters) : _track_parameters(track_parameters){};

    //CosmicTrack(XYZVec track_equation) : _track_equation(track_equation) {};
    /*
    CosmicTrack(XYZVec track_start, XYZVec track_dir) : _track_direction(track_dir), _track_point0(track_start) {
	
	_track_equation = track_start + _track_parameters.Cross(track_dir)


	};
    */
    // Destructor
    ~CosmicTrack();


    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);
    // Accessors:
    //int get_dim() const {return _dim;}
    int get_N() const { return _Nhits;}

    XYZVec get_parameters(){ 
	return _track_parameters;
	}

    double get_chisq() const { return _chisq; }
    double get_chisq_dof() const { return _chisq_dof; }
    
    XYZVec get_track_parameters() const { return _track_parameters; }
    XYZVec get_track_equation() {return _track_equation;}

    XYZVec get_track_direction(){ return _track_direction;}
   
    
    std::vector<double> get_fit_residuals() const { return _fit_residuals; }
    std::vector<double> get_fit_residual_errors() const { return _fit_residual_errors; }

    // Modiffiers:
    void clear();

    void set_N(int N){_Nhits = N;}

    void set_parameters(double par1,double par2, double par3, double par4){ 
	_par1=par1;
	_par2=par2;
	_par3=par3;
	_par4=par4;

	}

    void set_chisq(double chisq) { _chisq = chisq; }
    void set_chisq_dof(double chisq_dof) { _chisq_dof = chisq_dof; }

    void set_fit_residuals(double residual) { _fit_residuals.push_back(residual); }
    void set_fit_residual_errors(double residual_err) { _fit_residual_errors.push_back(residual_err); }

    void set_cov(int r, int c, double element){ _cov_mat[r][c] = element;}
    
  private:
  
    int _Nhits;

    double _par1;
    double _par2;
    double _par3;
    double _par4;

    XYZVec _track_parameters;
    XYZVec _track_equation;//r(t) expression

    XYZVec _track_direction;//the "gradient" term
    XYZVec _track_point0;//the "starting point" in fit line

    double _chisq;
    double _chisq_dof; 

    std::vector<double> _fit_residuals;
    std::vector<double> _fit_residual_errors;
    std::vector<std::vector<double> > _cov_mat;
    
    
  };
}

#endif

