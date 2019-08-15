#ifndef RecoDataProducts_StraightTrack_hh
#define RecoDataProducts_StraightTrack_hh
//S. Middleton - straight track class, will store info of fit result in terms of m and c and errors.
//This is now almost obsolete - upgrade to comsic track
#include "TMath.h"
#include "TMatrixD.h"
#include<vector>
#include<bitset>

using namespace std;

namespace mu2e {
  
  class StraightTrack{
  public:

    // Constructors
    StraightTrack();
    StraightTrack(double c_0, double m_0);
    StraightTrack( int N, double c_0, double c_0_err, double m_0, double m_0_err,double chisq, double chisq_dof, std::vector<double> fit_residual, std::vector<double> fit_residual_error);

    // Destructor
    ~StraightTrack();

    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);
    // Accessors:
    int get_N() const { return _Nhits;}
    double get_c_0() const { return _c_0; }
    double get_c_0_err() const { return _c_0_err; }
   
    //double get_cov(int row, int column){ return _cov_diag[row][column];}

    double get_m_0() const { return _m_0; }
    double get_m_0_err() const { return _m_0_err; }

    double get_chisq() const { return _chisq; }
    double get_chisq_dof() const { return _chisq_dof; }

    std::vector<double> get_fit_residuals() const { return _fit_residuals; }
    std::vector<double> get_fit_residual_errors() const { return _fit_residual_errors; }

    // Modiffiers:
    void clear();

    void set_N(int N){_Nhits = N;}

    void set_c_0(double c_0) { _c_0 = c_0; }
    void set_c_0_err(double c_0_err) { _c_0_err = c_0_err; }
   

    void set_m_0(double m_0) { _m_0 = m_0; }
    void set_m_0_err(double m_0_err) { _m_0_err = m_0_err; }

    void set_chisq(double chisq) { _chisq = chisq; }
    void set_chisq_dof(double chisq_dof) { _chisq_dof = chisq_dof; }
 
    void set_fit_residuals(double residual) { _fit_residuals.push_back(residual); }
    void set_fit_residual_errors(double residual_err) { _fit_residual_errors.push_back(residual_err); }

    void set_cov(int r, int c, 
double element){ _cov_mat[r][c] = element;}
    

    void set_parameters(int N, double c_0, double c_0_err, double m_0, double m_0_err,double chisq, double chisq_dof,std::vector<double> fit_residual, std::vector<double> fit_residual_error);

    
    
  private:

    int _Nhits;
    double _c_0;
    double _c_0_err;
    double _m_0;//rel x
    double _m_0_err;
    double _chisq;
    double _chisq_dof; 

    std::vector<double> _fit_residuals;
    std::vector<double> _fit_residual_errors;
    std::vector<std::vector<double> > _cov_mat;
    
    
  };
}

#endif

