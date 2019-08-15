#ifndef _MU2E_UTILITIES_PDFFit_HH
#define _MU2E_UTILITIES_PDFFit_HH
// Author: S. Middleton 
// Date: July 2019
//Purpose: Class of PDFs for liklihood functions which will be sent to Minuit for cosmic track analysis in tracker
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
    
    std::vector<double> time_residuals;

    ComboHitCollection chits;
    std::vector<Straw> straws;
    StrawResponse srep;
    CosmicTrack track;
    std::vector<double> docas;
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    
    int nparams =5;
    double doca_min = -100;
    double doca_max = 100;

     TimePDFFit(ComboHitCollection _chits, std::vector<Straw> &_straws, StrawResponse _srep, CosmicTrack _track, std::vector<double> _docas, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k) :  chits(_chits), straws(_straws), srep(_srep), track(_track), docas(_docas), constraint_means(_constraint_means), constraints(_constraints) {};
   
    double Up() const { return 0.5; };//10
    double operator() (const std::vector<double> &x) const;
    double TimeResidual(Straw straw, double doca, StrawResponse srep, double t0, ComboHit hit) const ;
    double calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit hit, std::vector<double> distance) const;
    void SetDOCAList(double d) { docas.push_back(d);}
    std::vector<double> GetDOCAList() const{ return docas;}
};

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
