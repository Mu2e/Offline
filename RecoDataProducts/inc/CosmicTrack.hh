#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
////S. Middleton, Feb 2019 - Cosmic track class, main purpose id to store diagnostics.
#include "TMath.h"
#include "TMatrixD.h"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include<vector>
#include<bitset>

using namespace std;

namespace mu2e {
  
  class CosmicTrack{
  public:
    // Constructors
    CosmicTrack();
    CosmicTrack(double par_1, double par_2, double par_3, double par_4) : _a0(par_1), _a1(par_2), _b0(par_3), _b1(par_4){};
    CosmicTrack(std::vector<double> track_parameters) : _track_parameters(track_parameters){};

    ~CosmicTrack();
    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);

 //---------------Accessors:--------------//
  double get_parameter(unsigned para_ID) const { 
	if(_track_parameters.size() >=para_ID){
		return _track_parameters.at(para_ID);
	}throw "Error: parameter list not long enough";
     }
  
    XYZVec get_track_direction() const{return _track_direction;}
    XYZVec get_track_position() const{return _track_position;}
    XYZVec get_initial_track_direction() const{ return _initial_track_direction;}
    
    XYZVec getXPrime() const { return _XPrime;}
    XYZVec getYPrime() const { return _YPrime;}
    XYZVec getZPrime() const { return _ZPrime;}
    
    XYZVec getinitXPrime() const { return _initXPrime;}
    XYZVec getinitYPrime() const { return _initYPrime;}
    XYZVec getinitZPrime() const { return _initZPrime;}
    
    double GetSagitta() const {return _Sagitta;}
     
    std::vector<double> get_track_parameters() const { return _track_parameters; } 
    std::vector<double> get_initial_track_parameters() const { return _initial_track_parameters; } 
    std::vector<double> get_true_track_parameters() const { return _true_track_parameters; }
    
    
    double get_true_phi() const { return true_track_angle_phi;}
    double get_fit_phi() const { return fit_track_angle_phi;}
    double get_true_theta() const { return true_track_angle_theta;}
    double get_fit_theta() const { return fit_track_angle_theta;}
    
    XYZVec get_true_track_direction() { return _true_track_direction;}
    std::vector<double> get_cov() const { return _cov;}
 
 //Diagnostics (MOVE THESE)
    int get_N() const { return _Nhits;}
    double get_finalchisq() const { return _finalchisq; }
    double get_finalchisq_dof() const { return _finalchisq_dof; }
    double get_true_finalchisq_dof() const { return _finalchisq_dofTrue; }
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
    
    double get_chi2_quant() const{ return _chi2_quant;}
   
    
    std::vector<double> get_init_fit_residualsX() const { return _initfit_residualsX; }
    std::vector<double> get_final_fit_residualsX() const { return _finalfit_residualsX; }
    std::vector<double> get_init_fit_residual_errorsX() const { return _initfit_residual_errorsX; }
    std::vector<double> get_final_fit_residual_errorsX() const { return _finalfit_residual_errorsX; }
    std::vector<double> get_init_fit_pullsX() const { return  _initfit_pullsX; }
    std::vector<double> get_final_fit_pullsX() const { return _finalfit_pullsX; }
    std::vector<double> get_init_fit_pullsY() const { return  _initfit_pullsY; }
    std::vector<double> get_final_fit_pullsY() const { return _finalfit_pullsY; }
    std::vector<double> get_final_fit_residualsY() const { return _finalfit_residualsY; }
    std::vector<double> get_final_fit_residual_errorsY() const { return _finalfit_residual_errorsY; }  
    std::vector<double> get_init_fit_residualsY() const { return _initfit_residualsY; }
    std::vector<double> get_init_fit_residual_errorsY() const { return _initfit_residual_errorsY; }
    //hit errors info
    std::vector<double> get_hit_errorsX() const{ return _hit_errorsX;}
    std::vector<double> get_hit_errorsY() const{ return _hit_errorsY;}
    std::vector<double> get_final_hit_errorsTotal() const{ return _finalhit_errorsTotal;} 
    std::vector<double> get_init_hit_errorsTotal() const{ return _inithit_errorsTotal;}
    int get_iter(){return _niters;}
    
    
//----------------END DIAGNOSTICS-----------------------//  
    
     void clear_all(); //clears track info and diagnostics
     void clear_parameters(); //clears just track info
     
//-------------Modiffiers of Track Parameters/Features---------------//
    void set_chi2_quant(double sum){_chi2_quant = sum;}
    void set_N(int N){_Nhits = N;}
    void set_parameters(double a0,double a1, double b0, double b1){ 
	_a0=a0;
	_a1=a1;
	_b0=b0;
	_b1=b1;
	_track_parameters.erase(_track_parameters.begin(), _track_parameters.end());
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
	
    void set_true_parameters(double a0,double a1, double b0, double b1){ 
	_truea0=a0;
	_truea1=a1;
	_trueb0=b0;
	_trueb1=b1;
        _true_track_parameters.push_back(a0);
        _true_track_parameters.push_back(a1);
	_true_track_parameters.push_back(b0);
	_true_track_parameters.push_back(b1);
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
    void setinitXPrime(XYZVec XPrime){ _initXPrime = XPrime;}
    void setinitYPrime(XYZVec YPrime){ _initYPrime = YPrime;}
    void setinitZPrime(XYZVec ZPrime){_initZPrime = ZPrime;}
    
    void set_true_track_direction(XYZVec direction){ _true_track_direction = direction;}
    void set_true_phi(double track_angle){true_track_angle_phi=track_angle;}
    void set_fit_phi(double track_angle){fit_track_angle_phi=track_angle;}
    void set_true_theta(double track_angle){true_track_angle_theta=track_angle;}
    void set_fit_theta(double track_angle){fit_track_angle_theta=track_angle;}
    void set_cov(double da0, double da1, double db0, double db1){
    	_cov.push_back(da0);
    	_cov.push_back(da1);
    	_cov.push_back(db0);
    	_cov.push_back(db1);
    }
    
    void set_mom(XYZVec mom){_track_mommentum=mom;}
    
    //_-----------Diganostics ---------------//
    void set_initchisq_dof(double initchisq_dof) { _initchisq_dof = initchisq_dof; }
    void set_finalchisq_dof(double finalchisq_dof) { _finalchisq_dof = finalchisq_dof; }
    void set_true_finalchisq_dof(double finalchisq_dof) { _finalchisq_dofTrue = finalchisq_dof; }
    void set_initchisq_dofX(double initchisq_dofX) { _initchisq_dofX = initchisq_dofX; }
    void set_finalchisq_dofX(double finalchisq_dofX) { _finalchisq_dofX = finalchisq_dofX; }
    void set_initchisq_dofY(double initchisq_dofY) { _initchisq_dofY = initchisq_dofY; }
    void set_finalchisq_dofY(double finalchisq_dofY) { _finalchisq_dofY = finalchisq_dofY; }
    
    void set_init_fit_residualsX(double residual) { _initfit_residualsX.push_back(residual); }
    void set_final_fit_residualsX(double residual) { _finalfit_residualsX.push_back(residual); }
    void set_init_fit_residual_errorsX(double residual_err) { _initfit_residual_errorsX.push_back(residual_err); }
    void set_final_fit_residual_errorsX(double residual_err) { _finalfit_residual_errorsX.push_back(residual_err); }
    void set_init_fit_residualsY(double residual) { _initfit_residualsY.push_back(residual); }
    void set_final_fit_residualsY(double residual) { _finalfit_residualsY.push_back(residual); }  
    void set_init_fit_residual_errorsY(double residual_err) { _initfit_residual_errorsY.push_back(residual_err); }
    void set_final_fit_residual_errorsY(double residual_err) { _finalfit_residual_errorsY.push_back(residual_err); }
    
    void set_init_pullsX(double residual) { _initfit_pullsX.push_back(residual); }
    void set_final_pullsX(double residual) { _finalfit_pullsX.push_back(residual); }
    void set_init_pullsY(double residual) { _initfit_pullsY.push_back(residual); }
    void set_final_pullsY(double residual) { _finalfit_pullsY.push_back(residual); }  
    
    void set_hit_errorsX(double hit_errorX){_hit_errorsX.push_back(hit_errorX);}
    void set_hit_errorsY(double hit_errorY){_hit_errorsY.push_back(hit_errorY);}
    void set_init_hit_errorsTotal(double hit_errorTotal){_inithit_errorsTotal.push_back(hit_errorTotal);}
    void set_final_hit_errorsTotal(double hit_errorTotal){_finalhit_errorsTotal.push_back(hit_errorTotal);}
    void set_niter(int iter){_niters= (iter);}
    //----------End Diagnostics ----------//
  private:
    int _Nhits;
    //Parameters:
    double _a0 = 0; //a0
    double _a1 = 0 ; //a1
    double _b0 = 0; //b0
    double _b1 = 0; //b1
    //true parameters:
    double _truea0 = 0; //a0
    double _truea1 = 0 ; //a1
    double _trueb0 = 0; //b0
    double _trueb1 = 0; //b1
    //Initial Estimates:
    double _inita0; //a0
    double _inita1; //a1
    double _initb0; //b0
    double _initb1;
    //Track info:
    double true_track_angle_phi;
    double fit_track_angle_phi;
    double true_track_angle_theta;
    double fit_track_angle_theta;
    
    std::vector<double> _track_parameters; //FIXME NEED TO BE 4D vector.....
    std::vector<double> _true_track_parameters;
    std::vector<double> _initial_track_parameters; 
   
    double _Sagitta;//TODO
    XYZVec _track_mommentum;//TODO
   
    XYZVec _track_equation;//r(t) expression
    XYZVec _track_direction;//the "gradient" term
    XYZVec _track_position;//the "starting point" in fit line
    XYZVec _initial_track_direction; // the first estimate of the line
    XYZVec _true_track_direction; // the first estimate of the line
    //Track Co-ordinate System:
    XYZVec _XPrime;
    XYZVec _YPrime;
    XYZVec _ZPrime;
    //Track Co-ordinate System:
    XYZVec _initXPrime;
    XYZVec _initYPrime;
    XYZVec _initZPrime;
    //Cov Matrix:
    std::vector<double> _cov;
    //Diagnostics:
    double _chi2_quant;
    double _finalchisq ;
    double _finalchisq_dof ;
    double _finalchisq_dofTrue ;
    double _initchisq ;
    double _initchisq_dof ;   
    double _finalchisqX ;
    double _finalchisq_dofX ;
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
    std::vector<double> _initfit_pullsY;
    std::vector<double> _initfit_pullsX;
    std::vector<double> _finalfit_pullsY;
    std::vector<double> _finalfit_pullsX;
    
    std::vector<double> _hit_errorsX;
    std::vector<double> _hit_errorsY;

    int _niters;
    
  };
}

#endif

