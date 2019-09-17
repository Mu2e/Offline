#ifndef _MU2E_UTILITIES_PDFFit_HH
#define _MU2E_UTILITIES_PDFFit_HH

// Class of PDFs for liklihood functions which will be sent to Minuit for cosmic track analysis in tracker. 
// NB: This Code is based on Richie's prototype analysis code.

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
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
    
    //std::vector<double> time_residuals;
  
    ComboHitCollection chits;
    std::vector<Straw> straws;
    StrawResponse srep;
    CosmicTrack track;
    std::vector<double> docas;
     
    std::vector<double> constraint_means;
    std::vector<double> constraints;

    int k;


    XYZVec Direction;
    int nparams =5;
    double doca_min = -10;
    double doca_max = 10;

     TimePDFFit(ComboHitCollection _chits, std::vector<Straw> &_straws, StrawResponse _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k) :  chits(_chits), straws(_straws), srep(_srep), track(_track), constraint_means(_constraint_means), constraints(_constraints) , k(_k) {};
   
    double Up() const { return 0.5; };
    double operator() (const std::vector<double> &x) const;
    double TimeResidual(Straw straw, double doca, StrawResponse srep, double t0, ComboHit hit) const ;
    double calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit hit) const;
    
};

/*
class FullFit : public TimePDFFit {
  public:
    FullFit(ComboHitCollection _chits, std::vector<Straw> &_straws, StrawResponse _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k);
    double voltage=1425.;
    void calculate_full_pdf();

    double interpolatePDF(double time_residual, double sigma, double tau) const;

    double pdf_mint, pdf_mintau, pdf_mins;
    double pdf_maxt, pdf_maxtau, pdf_maxs;
    double pdf_deltat, pdf_deltatau, pdf_deltas;
    double *pdf_sigmas, *pdf_taus, *pdf_times;
    double *pdf;
    int k;
};
*/
class PDFFit : public ROOT::Minuit2::FCNBase {
  public:
    std::vector<double> docas;
    ComboHitCollection chits;
    std::vector<Straw> straws;
    StrawResponse srep;
   
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    int nparams =5;
    double doca_min = 0;
    double doca_max = 1000;

    PDFFit(ComboHitCollection _chits,  std::vector<Straw> &_straws, StrawResponse _srep, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k) :  chits(_chits) , constraint_means(_constraint_means), constraints(_constraints) {};
   
    double Up() const { return 0.5; };
    double operator() (const std::vector<double> &x) const;
    double TimeResidual(Straw straw, double doca, StrawResponse srep, double t0, ComboHit hit) const ;
    double calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit hit) const;
   
};


#endif
