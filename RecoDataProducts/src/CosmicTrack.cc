#include "RecoDataProducts/inc/CosmicTrack.hh"
//#include "TMatrixD.h"
using namespace std;

namespace mu2e{
	CosmicTrack::CosmicTrack() {
          
	  _Nhits = 0;
	  _c = 0.0;
	  _c_err = 0.0;
	  _m = 0.0;
	  _m_err = 0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
          
	}

	CosmicTrack::CosmicTrack(double c, double m) {
	  
	  _Nhits = 0;
	  _c = c;
	  _c_err = 0.0;
	  _m = m;
	  _m_err = 0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
          
	}

	CosmicTrack::CosmicTrack(int N, double c, double c_err, double m, double m_err,double chisq, double chisq_dof) {
 
	  _Nhits = N;
	  _c = c;
	  _c_err = c_err;
	  _m = m;
	  _m_err = m_err;
	  _chisq = chisq;
	  _chisq_dof = chisq_dof;
          
        }

	// Destructor
	CosmicTrack::~CosmicTrack() {}

	void CosmicTrack::clear() {
	
	  _Nhits = 0;
	  _c = 0.0;
	  _c_err = 0.0;
	  _m = 0.0;
	  _m_err = 0.0;
	  _chisq = 0.0;
	  _chisq_dof = 0.0;
	  
	}

	void CosmicTrack::set_parameters(int N, double c, double c_err, double m, double m_err,double chisq, double chisq_dof) {
	 
	  _Nhits = N;
	  _c = c;
	  _c_err = c_err;
	  _m = m;
	  _m_err = m_err;
	  _chisq = chisq;
	  _chisq_dof = chisq_dof;
	  
	 

	}
}

