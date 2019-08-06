#ifndef _MU2E_UTILITIES_PDFFit_HH
#define _MU2E_UTILITIES_PDFFit_HH
// Author: S. Middleton 
// Date: July 2019
//Purpose: Class of PDFs for liklihood functions which will be sent to Minuit for cosmic track analysis in tracker
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
//Tracker Details:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

//ROOT
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
//Minuit
#include <Minuit2/FCNBase.h>


using namespace mu2e;

class TimePDFFit : public ROOT::Minuit2::FCNBase {
  public:
    std::vector<double> docas;
    std::vector<double> combohit_times;
    std::vector<XYZVec> combohit_positions; 
    std::vector<Straw> straws;
    StrawResponse srep;
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    int nparams =5;
    double doca_min = 0;
    double doca_max = 1000;

    TimePDFFit(std::vector<double> &_combohit_times, std::vector<XYZVec> &_combohit_positions, std::vector<Straw> &_straws, StrawResponse _srep, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k) :  combohit_times(_combohit_times), combohit_positions (_combohit_positions), straws(_straws), srep(_srep), constraint_means(_constraint_means), constraints(_constraints) {};
    
    double Up() const { return 0.5; };
    
    double operator() (const std::vector<double> &x) const;
    double TimeResidual(Straw straw, double doca, double time, StrawResponse srep) const ;
    double calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1) const;

};

class PDFFit : public ROOT::Minuit2::FCNBase {
  public:
    std::vector<double> docas;
    std::vector<double> combohit_times;
    std::vector<XYZVec> combohit_positions; //vector of 3D positions
    std::vector<double> errorsX;
    std::vector<double> errorsY;
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    int nparams =5;
    double doca_min = 0;
    double doca_max = 1000;

    PDFFit(std::vector<double> &_combohit_times, std::vector<XYZVec> &_combohit_positions, std::vector<double> _errorsX, std::vector<double> _errorsY, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k) :  combohit_times(_combohit_times), combohit_positions (_combohit_positions), errorsX(_errorsX), errorsY(_errorsY), constraint_means(_constraint_means), constraints(_constraints) {};
   
    double Up() const { return 0.5; };
    double operator() (const std::vector<double> &x) const;
   
};


#endif
