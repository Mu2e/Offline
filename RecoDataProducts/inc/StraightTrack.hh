#ifndef RecoDataProducts_StraightTrack_hh
#define RecoDataProducts_StraightTrack_hh
//S. Middleton - straight track class, will store info of fit result in terms of m and c and errors.
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
    StraightTrack(double c_0, double m_0, double m_1, double m_2);
    StraightTrack(int dim, int N, double c_0, double c_0_err, double m_0, double m_0_err,double m_1, double m_1_err, double m_2, double m_2_err, double chisq, double chisq_dof, std::bitset<32ul> leftHit, std::bitset<32ul> rightHit, std::vector<double> fit_residual, std::vector<double> fit_residual_error);

    // Destructor
    ~StraightTrack();


    std::vector<std::vector<double>> Initialize_Cov_Mat(int i, int j, int sizei, int sizej);
    // Accessors:
    //int get_dim() const {return _dim;}
    int get_N() const { return _Nhits;}
    double get_c_0() const { return _c_0; }
    double get_c_0_err() const { return _c_0_err; }
   
    //double get_cov(int row, int column){ return _cov_diag[row][column];}

    double get_m_0() const { return _m_0; }
    double get_m_0_err() const { return _m_0_err; }

    double get_m_1() const { return _m_1; }
    double get_m_1_err() const { return _m_1_err; }

    double get_m_2() const { return _m_2; }
    double get_m_2_err() const { return _m_2_err; }

    double get_chisq() const { return _chisq; }
    double get_chisq_dof() const { return _chisq_dof; }

    double get_u_chisq() const { return _uchisq; }
    double get_u_dof() const { return _udof; }

    double get_v_chisq() const { return _vchisq; }
    double get_v_dof() const { return _vdof; }
    
    std::bitset<32ul> getLHit(){ return _leftHit;}

    std::bitset<32ul> getRHit(){ return _rightHit;}
    
    std::vector<double> get_fit_residuals() const { return _fit_residuals; }
    std::vector<double> get_fit_residual_errors() const { return _fit_residual_errors; }

    // Modiffiers:
    void clear();

    //void set_dim(int dim){_dim = dim;}
    void set_N(int N){_Nhits = N;}

    void set_c_0(double c_0) { _c_0 = c_0; }
    void set_c_0_err(double c_0_err) { _c_0_err = c_0_err; }
   

    void set_m_0(double m_0) { _m_0 = m_0; }
    void set_m_0_err(double m_0_err) { _m_0_err = m_0_err; }

    void set_m_1(double m_1) { _m_1 = m_1; }
    void set_m_1_err(double m_1_err) { _m_1_err = m_1_err; }

    void set_m_2(double m_2) { _m_2 = m_2; }
    void set_m_2_err(double m_2_err) { _m_2_err = m_2_err; }

    void set_chisq(double chisq) { _chisq = chisq; }
    void set_chisq_dof(double chisq_dof) { _chisq_dof = chisq_dof; }
 
    void set_v_chisq(double chisq) { _vchisq = chisq; }
    void set_u_chisq(double chisq) { _uchisq = chisq; }

    void set_u_NDF(int udof){ _udof = udof;}
    void set_v_NDF(int vdof){ _vdof = vdof;}

    
    void set_fit_residuals(double residual) { _fit_residuals.push_back(residual); }
    void set_fit_residual_errors(double residual_err) { _fit_residual_errors.push_back(residual_err); }

    void set_cov(int r, int c, 
double element){ _cov_mat[r][c] = element;}
    
    void set_L(std::bitset<32ul>& L){_leftHit = L;}
    void set_R(std::bitset<32ul>& R){_rightHit = R;}

    void set_parameters(int dim, int N, double c_0, double c_0_err, double m_0, double m_0_err,double m_1, double m_1_err, double m_2, double m_2_err, double chisq, double chisq_dof, std::bitset<32ul> leftHit, std::bitset<32ul> rightHit, std::vector<double> fit_residual, std::vector<double> fit_residual_error);

    
    
  private:
    //int _dim;
    int _Nhits;
    double _c_0;
    double _c_0_err;
    double _m_0;//rel x
    double _m_0_err;
    double _m_1;//rel y
    double _m_1_err;
    double _m_2;//rel z
    double _m_2_err;
    double _chisq;
    double _chisq_dof; 
   
    double _uchisq;
    double _udof; 

    double _vchisq;
    double _vdof;

    std::bitset<32ul> _leftHit, _rightHit;

    std::vector<double> _fit_residuals;
    std::vector<double> _fit_residual_errors;
    std::vector<std::vector<double> > _cov_mat;
    
    
  };
}

#endif

