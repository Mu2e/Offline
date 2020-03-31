#ifndef _COSMIC_RECO_PDFFit_HH
#define _COSMIC_RECO_PDFFit_HH

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
		StrawResponse const& srep ;
		CosmicTrack track;
		std::vector<double> docas;

		std::vector<double> constraint_means;
		std::vector<double> constraints;
		double sigma_t;
		int k;

		const Tracker* tracker;
		std::vector<double> DOCAs;
		std::vector<double> TimeResiduals;

		int nparams =5; 

		GaussianPDFFit(ComboHitCollection _chits, StrawResponse const& _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double _sigma_t, int _k,  const Tracker* _tracker) :  chits(_chits),  srep(_srep), track(_track), constraint_means(_constraint_means), constraints(_constraints) , sigma_t(_sigma_t), k(_k), tracker(_tracker) {};

		double Up() const { return 0.5; };
		double operator() (const std::vector<double> &x) const;
		double TimeResidual(double doca, StrawResponse const& srep , double t0, ComboHit const& hit,  const Tracker* tracker) const ;
		double calculate_DOCA(ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker) const;
		double calculate_ambig(ComboHit const& chit, double a0, double a1, double b0, double b1,  const Tracker* tracker) const;

    };

    class FullDriftFit : public GaussianPDFFit {
    	public:
		FullDriftFit(ComboHitCollection _chits,  StrawResponse const& _srep, CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double sigma_t, int _k, const Tracker* _tracker);
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
