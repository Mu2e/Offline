#ifndef _MU2E_UTILITIES_LIKLIHOODFUNCTIONS_HH
#define _MU2E_UTILITIES_LIKLIHOODFUNCTIONS_HH
// Author: S. Middleton 
// Date: July 2019
//Purpose: Will pass PDF function to Minuit 
#include "TrackerConditions/inc/StrawDrift.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"

//ROOT
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
//Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;
struct EndResult{
        public:
		std::vector<std::string> names;
		std::vector<double> bestfit;
		std::vector<double> bestfiterrors;
	
};

namespace LiklihoodFunctions {
	void DoFit(std::vector<double> times, std::vector<double> docas, std::vector<double> time_residuals);

}


#endif
