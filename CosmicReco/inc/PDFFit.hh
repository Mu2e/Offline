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

class GaussianPDFFit : public ROOT::Minuit2::FCNBase {
  public:
    ComboHitCollection chits;
    std::vector<Straw> straws;
    StrawResponse srep;
    CosmicTrack track;
    std::vector<double> docas;
     
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    double sigma_t;
    int k;
    std::vector<double> DOCAs;
    std::vector<double> TimeResiduals;
    
    int nparams =5; 
    
     GaussianPDFFit(ComboHitCollection _chits, std::vector<Straw> &_straws, StrawResponse _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double _sigma_t, int _k) :  chits(_chits), straws(_straws), srep(_srep), track(_track), constraint_means(_constraint_means), constraints(_constraints) , sigma_t(_sigma_t), k(_k) {};
   
    double Up() const { return 0.5; };
    double operator() (const std::vector<double> &x) const;
    double TimeResidual(Straw straw, double doca, StrawResponse srep, double t0, ComboHit hit) const ;
    double calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1) const;
    double calculate_ambig(Straw const& straw, double a0, double a1, double b0, double b1) const;
    
};

class FullDriftFit : public GaussianPDFFit {
  public:
    FullDriftFit(ComboHitCollection _chits, std::vector<Straw> &_straws, StrawResponse _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double sigma_t, int _k);
int Factorial(int k);
    void CalculateFullPDF();
    double InterpolatePDF(double time_residual, double sigma, double tau) const;
    void DeleteArrays() const;
    double Min_t;
    double  Min_tau, Min_s;
    double Max_t, Max_tau, Max_s;
    double delta_T, delta_Tau, delta_S;
    double *pdf_sigmas, *pdf_taus, *pdf_times;
    double *pdf;
    int k;
    double operator() (const std::vector<double> &x) const;
};

#endif
