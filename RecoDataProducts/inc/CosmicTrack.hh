#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
////S. Middleton, Feb 2019 - Cosmic track class, will store info of fit result in terms of m and c and errors.
#include "TMath.h"
#include "TMatrixD.h"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
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
    
    CosmicTrack(double par_1, double par_2, double par_3, double par_4) : _a0(par_1), _a1(par_2), _b0(par_3), _b1(par_4){};

    CosmicTrack(std::vector<double> track_parameters) : _track_parameters(track_parameters){};

    
    // Destructor
    ~CosmicTrack();


    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);

 //---------------Accessors:--------------//
    int get_N() const { return _Nhits;}
    double get_finalchisq() const { return _finalchisq; }
    double get_finalchisq_dof() const { return _finalchisq_dof; }
    double get_initchisq() const { return _initchisq; }
    double get_initchisq_dof() const { return _initchisq_dof; }
    
    double get_finalchisqX() const { return _finalchisqX; }
    double get_finalchisq_dofX() const { return _finalchisq_dofX; }
    double get_initchisqX() const { return _initchisqX; }
    double get_initchisq_dofX() const { return _initchisq_dofX; }
    
    double get_finalchisqY() const { return _finalchisqY; }
    double get_finalchisq_dofY() const { return _finalchisq_dofY; }
    double get_initchisqY() const { return _initchisqY; }
    double get_initchisq_dofY() const { return _initchisq_dofY; }
    
    std::vector<double> get_track_parameters() const { return _track_parameters; } 
    std::vector<double> get_initial_track_parameters() const { return _initial_track_parameters; } 
    
    double get_parameter(unsigned para_ID) const { 
	if(_track_parameters.size() >=para_ID){
		return _track_parameters.at(para_ID);
	}throw "Error: parameter list not long enough";
     }
     
    double get_track_length() const {return _track_length;}
    XYZVec get_track_equation() const{return _track_equation;}
    XYZVec get_track_direction() const{return _track_direction;}
    XYZVec get_track_position() const{return _track_position;}
    XYZVec get_initial_track_direction() const{ return _initial_track_direction;}
    
    XYZVec getXPrime() const { return _XPrime;}
    XYZVec getYPrime() const { return _YPrime;}
    XYZVec getZPrime() const { return _ZPrime;}
    
    double GetSagitta() const {return _Sagitta;}
    
    std::vector<double> get_init_fit_residualsX() const { return _initfit_residualsX; }
    std::vector<double> get_final_fit_residualsX() const { return _finalfit_residualsX; }
    std::vector<double> get_init_fit_residual_errorsX() const { return _initfit_residual_errorsX; }
    std::vector<double> get_final_fit_residual_errorsX() const { return _finalfit_residual_errorsX; }
    
    std::vector<double> get_final_fit_residualsY() const { return _finalfit_residualsY; }
    std::vector<double> get_final_fit_residual_errorsY() const { return _finalfit_residual_errorsY; }  
    std::vector<double> get_init_fit_residualsY() const { return _initfit_residualsY; }
    std::vector<double> get_init_fit_residual_errorsY() const { return _initfit_residual_errorsY; }
    //hit errors info
    std::vector<double> get_hit_errorsX() const{ return _hit_errorsX;}
    std::vector<double> get_hit_errorsY() const{ return _hit_errorsY;}
    std::vector<double> get_final_hit_errorsTotal() const{ return _finalhit_errorsTotal;}
    
    std::vector<double> get_init_hit_errorsTotal() const{ return _inithit_errorsTotal;}
    
    std::vector<int> get_iter(){return _niters;}
    std::vector<ComboHit*> get_outliers()const {return _outliers;}
//---------------------------------------//  
    
     void clear_all(); //clears track info and diagnostics
     void clear_parameters(); //clears just track info
     
//-------------Modiffiers---------------//
   
    void set_N(int N){_Nhits = N;}
    void set_track_length(double length) { _track_length = length; }
    void set_parameters(double a0,double a1, double b0, double b1){ 
	_a0=a0;
	_a1=a1;
	_b0=b0;
	_b1=b1;
        _track_parameters.push_back(a0);
        _track_parameters.push_back(a1);
	_track_parameters.push_back(b0);
	_track_parameters.push_back(b1);

	}
	
    void set_initial_parameters(double a0,double a1, double b0, double b1){ 
	_inita0=a0;
	_inita1=a1;
	_initb0=b0;
	_initb1=b1;
        _initial_track_parameters.push_back(a0);
        _initial_track_parameters.push_back(a1);
	_initial_track_parameters.push_back(b0);
	_initial_track_parameters.push_back(b1);

	}
	
    void set_parameters(unsigned para_ID, double par){
    	if(_track_parameters.size() >= para_ID){
		_track_parameters.at(para_ID) = par;
	}throw "Error: Parameter list not long enough";
	 
	
    }
    void set_track_direction(XYZVec direction){_track_direction = direction;}
    void set_track_position(double x, double y , double z){_track_position.SetXYZ(x,y,z);}
    void set_initial_track_direction(XYZVec direction){_initial_track_direction = direction;}
    
    void setXPrime(XYZVec XPrime){ _XPrime = XPrime;}
    void setYPrime(XYZVec YPrime){ _YPrime = YPrime;}
    void setZPrime(XYZVec ZPrime){_ZPrime = ZPrime;}
    
    void set_initchisq_dof(double initchisq_dof) { _initchisq_dof = initchisq_dof; }
    void set_finalchisq_dof(double finalchisq_dof) { _finalchisq_dof = finalchisq_dof; }
    
    void set_initchisq_dofX(double initchisq_dofX) { _initchisq_dofX = initchisq_dofX; }
    void set_finalchisq_dofX(double finalchisq_dofX) { _finalchisq_dofX = finalchisq_dofX; }
    
    void set_initchisq_dofY(double initchisq_dofY) { _initchisq_dofY = initchisq_dofY; }
    void set_finalchisq_dofY(double finalchisq_dofY) { _finalchisq_dofY = finalchisq_dofY; }
    
    void set_mom(XYZVec mom){_track_mommentum=mom;} 
    
    void set_init_fit_residualsX(double residual) { _initfit_residualsX.push_back(residual); }
    void set_final_fit_residualsX(double residual) { _finalfit_residualsX.push_back(residual); }
    void set_init_fit_residual_errorsX(double residual_err) { _initfit_residual_errorsX.push_back(residual_err); }
    void set_final_fit_residual_errorsX(double residual_err) { _finalfit_residual_errorsX.push_back(residual_err); }
    
    void set_init_fit_residualsY(double residual) { _initfit_residualsY.push_back(residual); }
    void set_final_fit_residualsY(double residual) { _finalfit_residualsY.push_back(residual); }  
    void set_init_fit_residual_errorsY(double residual_err) { _initfit_residual_errorsY.push_back(residual_err); }
    void set_final_fit_residual_errorsY(double residual_err) { _finalfit_residual_errorsY.push_back(residual_err); }
     
    void set_hit_errorsX(double hit_errorX){_hit_errorsX.push_back(hit_errorX);}
    void set_hit_errorsY(double hit_errorY){_hit_errorsY.push_back(hit_errorY);}
    
    void set_init_hit_errorsTotal(double hit_errorTotal){_inithit_errorsTotal.push_back(hit_errorTotal);}
    void set_final_hit_errorsTotal(double hit_errorTotal){_finalhit_errorsTotal.push_back(hit_errorTotal);}
    
    void set_cov(int r, int c, double element){ _cov_mat[r][c] = element;}
    
    void add_outlier(ComboHit* hit){_outliers.push_back(hit);}
    
    void set_niter(int iter){_niters.push_back(iter);}
  private:
    
    int _Nhits;

    double _a0; //a0
    double _a1; //a1
    double _b0; //b0
    double _b1; //b1
    
    double _inita0; //a0
    double _inita1; //a1
    double _initb0; //b0
    double _initb1;

    std::vector<double> _track_parameters; //FIXME NEED TO BE 4D vector.....
    std::vector<double> _initial_track_parameters; 
    double _track_length;
    double _Sagitta;//TODO
    XYZVec _track_mommentum;//TODO
    
    XYZVec _track_equation;//r(t) expression
    XYZVec _track_direction;//the "gradient" term
    XYZVec _track_position;//the "starting point" in fit line
    XYZVec _initial_track_direction; // the first estimate of the line
    //Get Track Co-ordinate System:
    XYZVec _XPrime;
    XYZVec _YPrime;
    XYZVec _ZPrime;
    
    double _finalchisq;
    double _finalchisq_dof;
    
    double _initchisq;
    double _initchisq_dof;
    
     double _finalchisqX;
    double _finalchisq_dofX;
    
    double _initchisqX;
    double _initchisq_dofX;
    
     double _finalchisqY;
    double _finalchisq_dofY;
    
    double _initchisqY;
    double _initchisq_dofY;
    
    std::vector<double> _initfit_residualsX;
    std::vector<double> _initfit_residual_errorsX;
    std::vector<double> _initfit_residualsY;
    std::vector<double> _initfit_residual_errorsY;
    
    std::vector<double> _finalfit_residualsX;
    std::vector<double> _finalfit_residual_errorsX;
    std::vector<double> _finalfit_residualsY;
    std::vector<double> _finalfit_residual_errorsY;
    
    std::vector<double> _inithit_errorsTotal;
    std::vector<double> _finalhit_errorsTotal;
    
    
    std::vector<double> _hit_errorsX;
    std::vector<double> _hit_errorsY;
    
    std::vector<std::vector<double> > _cov_mat;
    std::vector<ComboHit*> _outliers;
    
    std::vector<int> _niters;
    
  };
}

#endif

