#ifndef _MU2E_UTILITIES_PDFFit_HH
#define _MU2E_UTILITIES_PDFFit_HH
// Author: S. Middleton 
// Date: July 2019
//Purpose: Class of PDFs for liklihood functions which will be sent to Minuit for cosmic track analysis in tracker
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

class PDFFit : public ROOT::Minuit2::FCNBase {
  public:
    std::vector<double> docas;
    std::vector<double> times;
    std::vector<double> time_residuals;
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    int nparams =4;
    double voltage = 1425.;
    double pdf_mint = -10;
    double pdf_mintau = 0.1;
    double pdf_mins = 0.1;
    double pdf_maxt = 5000;//40;
    double pdf_maxtau = 2000;//20.;
    double pdf_maxs = 10000;//10.;

    PDFFit(std::vector<double> &_docas, std::vector<double> &_times, std::vector<double> &_time_residuals, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k, double voltage) : docas(_docas), times(_times), time_residuals(_time_residuals), constraint_means(_constraint_means), constraints(_constraints) {};
    //Error definition of the function. MINUIT defines Parameter errors as the change in Parameter Value required to change the function Value by up.
    double Up() const { return 0.5; };
    //Logliklihood
    double operator() (const std::vector<double> &x) const;
    //PDF:
    void calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min=-1, double doca_max=1e20, bool dist=false) const;
    double TimeResidual(double doca, double time, double time_offset) const;
    //void calculate_full_pdf();

};


#endif
