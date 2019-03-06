#include "RecoDataProducts/inc/StraightTrack.hh"
#include <iostream>
#include <sstream>
#include <exception>
#include <numeric>
#include <sstream>

using namespace std;

namespace mu2e{
	StraightTrack::StraightTrack() {
          //_dim = 3;
	  _Nhits = 0;
	  _c_0 = 0.0;
	  _c_0_err = 0.0;
	  _m_0 = 0.0;
	  _m_0_err = 0.0;
          _m_1 = 0.0;
          _m_1_err =0.0;
	  _m_2 = 0.0;
          _m_2_err =0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
	  
	  
	}

	StraightTrack::StraightTrack(double c_0,  double m_0, double m_1, double m_2) {
	  //_dim = 3;
	  _Nhits = 0;
	  _c_0 = c_0;
	  _c_0_err = 0.0;
	  _m_0 = m_0;
	  _m_0_err = 0.0;
          _m_1 = m_1;
	  _m_1_err = 0.0;
	  _m_2 = 0.0;
          _m_2_err =0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
         
	}

	StraightTrack::StraightTrack(int dim, int N, double c_0, double c_0_err, double m_0, double m_0_err,double m_1, double m_1_err, double m_2, double m_2_err, double chisq, double chisq_dof, std::bitset<32ul> leftHit, std::bitset<32ul> rightHit, std::vector<double> _fit_residual, std::vector<double> _fit_residual_error) {
          //_dim = 3;
	  _Nhits = N;
	  _c_0 = c_0;
	  _c_0_err = c_0_err;
	  _m_0 = m_0;
	  _m_0_err = m_0_err;
          _m_1 = m_1;
	  _m_1_err = m_1_err;
	  _m_2 = m_2;
	  _m_2_err = m_2_err;
	  _chisq = chisq;
	  _chisq_dof = chisq_dof;
          
        }
	
	std::vector<std::vector<double>>StraightTrack::Initialize_Cov_Mat(int i, int j, int size_i, int size_j){
	  _cov_mat[size_i][size_j] = { 0 };
	  return _cov_mat;
        }

	// Destructor
	StraightTrack::~StraightTrack() {}

	void StraightTrack::clear() {
	  //_dim = 3;
	  _Nhits = 0;
	  _c_0 = 0.0;
	  _c_0_err = 0.0;
	  _m_0 = 0.0;
	  _m_0_err = 0.0;
	  _m_1 = 0.0;
	  _m_1_err = 0.0;
	  _m_2 = 0.0;
	  _m_2_err = 0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
	  _fit_residuals.erase(_fit_residuals.begin(),_fit_residuals.end());
	  _fit_residual_errors.erase(_fit_residual_errors.begin(),_fit_residual_errors.end());
	  
	}
        
	void StraightTrack::set_parameters(int dim, int N, double c_0, double c_0_err, double m_0, double m_0_err,double m_1, double m_1_err, double m_2, double m_2_err, double chisq, double chisq_dof, std::bitset<32ul> leftHit, std::bitset<32ul> rightHit, std::vector<double> fit_residual, std::vector<double> fit_residual_error) {
	  //_dim = dim;
	  _Nhits = N;
	  _c_0 = c_0;
	  _c_0_err = c_0_err;
	  _m_0 = m_0;
	  _m_0_err = m_0_err;
          _m_1= m_1;
	  _m_1_err = m_1_err;
          _m_2= m_2;
	  _m_2_err = m_2_err;
	  _chisq = chisq;
	  _chisq_dof = chisq_dof;
	  _fit_residuals = fit_residual;
	  _fit_residual_errors= fit_residual_error;

	}
     

	
}

