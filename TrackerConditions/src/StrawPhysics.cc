//
// StrawPhysics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawPhysics.cc,v 1.8 2014/03/16 15:12:54 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/16 15:12:54 $
//
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawPhysics.hh"
// boost
#include "boost/math/special_functions/gamma.hpp"
#include <math.h>
#include <algorithm>

using namespace std;
namespace mu2e {
  StrawPhysics::StrawPhysics(fhicl::ParameterSet const& pset) :
    _EIonize(pset.get<double>("EnergyPerIonization",1.6e-5)), // 100% Ar is between 27 ev/ionization and 100 ev/ionization, this number is tuned to reproduce the correct mean-free-path from the G4 energy deposits when setting the number of electrons/cluster to the distribution below
    _meanpath(pset.get<double>("MeanFreePath",0.357)), // mm, average distance between ionizations for STP Ar (Blum etal, table 1.1)
    _QIonize(pset.get<double>("ChargePerIonization",1.6e-7)), // e, pC
    _gasgain(pset.get<double>("GasGain",7.0e4)),
    _polyaA(pset.get<double>("PolyaA",1.25)), // A = 1/f = theta + 1.  A=1 -> exponential, A=infinity->delta-function
    _gslope(pset.get<double>("GainRMSSlope",0.809)), // slope of relative gain sigma on 1/sqrt(n)
    _nggauss(pset.get<unsigned>("NGainGauss",30)), // number of electrons/cluster to switch to a Gaussian model of the gain fluctuations
    _attlen{pset.get<double>("ShortAttentuationLength",50.0),pset.get<double>("LongAttentuationLength",27000.0)}, // from ATLAS TRT measurement
    _longfrac(pset.get<double>("LongAttenuationFraction",0.92)),
    _vprop(pset.get<double>("PropagationVelocity",273.0)), //mm/nsec
    _vdisp(pset.get<double>("PropagationVelocityDispersion",0.01)), //1/nsec
    _cdpoly(pset.get<vector<double> >("ClusterDriftPolynomial",vector<double>{0.0,16.0})), // linear term has units nanoseconds/mm
    _cdsigmapoly(pset.get<vector<double> >("ClusterDriftSigmaPolynomial",vector<double>{0.2})) // constant drift time error for now (ns)
  {
    // integrate the number of ionizations
    double ptot(0.0);
    std::vector<double> nProb = pset.get<vector<double> >("ProbPerCharge",vector<double>{0.656,0.15,0.064,0.035,0.0225,0.0155,0.0105,
      0.0081,0.0061, 0.0049, 0.0039, 0.0030, 0.0025, 0.0020, 0.0016, 0.0012, 0.00095, 0.00075}); // Blum, table 1.4
// renormalize these probs
    for(unsigned iprob=0;iprob< nProb.size(); ++iprob ) {
      ptot += nProb[iprob];
    }
    double norm = 1.0/ptot;
    double psum(0.0);
    _EAverage = _QAverage = 0.0;
    for(unsigned iprob=0;iprob< nProb.size(); ++iprob ) {
      double nprob = nProb[iprob]*norm;
      _intNProb.push_back(psum + nprob);
      psum += nprob;
      _EAverage += nprob*ionizationEnergy(iprob+1);
      _QAverage += nprob*ionizationCharge(iprob+1);
    }
    // now compute the average charge and energy per ionization from this distrubion
    
  }

  StrawPhysics::~StrawPhysics() {}

  // model gain fluctuations for a cluster using the Polya function.  This depends on the # of electrons
  double StrawPhysics::clusterGain(CLHEP::RandGaussQ& rgauss, CLHEP::RandFlat& rflat, unsigned nele) const {
    double gain(0.0);
    if(nele < _nggauss){
      // The theta parameter scales when nele>1
      double aval;
      if(nele == 1)
	aval = _polyaA;
      else 
	aval = nele*(_polyaA+1)-1;
      // use the inverse incomplete gamma function, which is directly related to the integral of the Polya function
      gain = strawGain()*boost::math::gamma_p_inv(aval,rflat.fire())/aval;
    } else {
// use a Gaussian approximation.
      double sigma = _gslope/sqrt(nele);
      gain = strawGain()*rgauss.fire(1.0,sigma);
    }
    return gain;
  }

  double StrawPhysics::driftDistanceToTime(double ddist, double phi) const {
    double retval(0.0);
    double arg(1.0);
    for(size_t ipow=0;ipow<_cdpoly.size();++ipow){
      retval += arg*_cdpoly[ipow];
      arg*=ddist;
    }
    return max(0.0,retval);
  }

  double StrawPhysics::driftTimeSpread(double ddist, double phi) const {
    double retval(0.0);
    double arg(1.0);
    for(size_t ipow=0;ipow<_cdsigmapoly.size();++ipow){
      retval += arg*_cdsigmapoly[ipow];
      arg*=ddist;
    }
    return max(0.0,retval);
  }
  
  double StrawPhysics::propagationAttenuation(double wdist) const {
    return  (1.0-_longfrac)*exp(-wdist/_attlen[0]) +
      _longfrac*exp(-wdist/_attlen[1]);
  }

  double StrawPhysics::propagationTime(double wdist) const {
    return  wdist/_vprop;
  }

  unsigned StrawPhysics::nePerIon(double urand) const {
  // sample the distribution for the number of ionizations.
    unsigned nele(_intNProb.size()+1);
    for(unsigned ie=0;ie< _intNProb.size(); ++ie){
      if(urand < _intNProb[ie]){
	// if the nput random number is below the integral for this number, we're done
	nele = ie+1;
	break;
      }
    }
    return nele;
  }

  unsigned StrawPhysics::nePerEIon(double EIon) const {
    // approximate assuming linear increase in ionization energy per electron.
    return std::min(static_cast<unsigned>(std::max(1.0,floor(0.5*(sqrt(1.0+8.0*EIon/_EIonize)-1.0)))),(unsigned)_intNProb.size());
  }

}

