#ifndef _MU2E_UTILITIES_LEASTSQUARESFITTER_HH
#define _MU2E_UTILITIES_LEASTSQUARESFITTER_HH
//Least Squares Fitting Routines for linear and circle fits
// Author: S. Middleton
// Date: Nov 2018
//c++
#include <vector>
//ROOT
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "Math/VectorUtil.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
using namespace mu2e;

namespace LeastSquaresFitter{
  	void xy_fit(const std::vector<double> &_x,  const std::vector<double> &_y,const std::vector<double> &_y_err, StraightTrack* line, TMatrixD& covariance);

	void xyz_fit(int _nCoeffs, const std::vector<double> &_x,  const std::vector<double> &_y,const std::vector<double> &_z, const std::vector<double> &_z_err, StraightTrack* line, TMatrixD& covariance);

	std::vector <double> Means(std::vector<std::vector<double> > values, std::vector<double> weights);
   
	double Mean(std::vector<double>  values, std::vector<double> weights);
	
        

} // ~namespace 

#endif

