
#include "Offline/TrackerConditions/inc/StrawPhysics.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

// boost
#include "boost/math/special_functions/gamma.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Matrix/Vector.h"
//--------------------

using namespace std;

namespace mu2e {


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
    // Find the maximum number of electrons this energy could produce
    unsigned ie(0);
    while(ie < _EIonize.size() && EIon > _EIonize[ie]){
      ++ie;
    }
    // require at least 1 electron
    return std::max((unsigned)1,ie);
  }


  // model gain fluctuations for a cluster using the Polya function.
  // This depends on the # of electrons
  double StrawPhysics::clusterGain(CLHEP::RandGaussQ& rgauss,
      CLHEP::RandFlat& rflat, unsigned nele) const {
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


  // get drift time
  double StrawPhysics::driftDistanceToTime(double ddist, double phi) const{
    if(_nonlindrift){
      return _strawDrift->D2T(ddist,phi);
    }
    else{
      return ddist/0.0625; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }


  // a debuging function
  double StrawPhysics::testStrawDrift(double ddist, double phi){
    double lorentzTime = _strawDrift->D2T(ddist,phi);
    double lorentzDist = _strawDrift->T2D(lorentzTime,phi);
    // cout <<"ddist - lorentzDist = "<<ddist - lorentzDist<<"\n";
    return ddist - lorentzDist;
  }


  double StrawPhysics::driftTimeSpread(double ddist) const {
    return _dtvar[0] + ddist*_dtvar[1];
  }


  double StrawPhysics::propagationTime(double wdist) const {
    return  wdist/_vprop;
  }


  double StrawPhysics::ionizationEnergy(unsigned nele) const {
    if(nele >=1 && nele <= _EIonize.size())
      return _EIonize[nele-1];
    else if(nele<1)
      return 0.0;
    else
      return _EIonize.back();
  }


  void StrawPhysics::print(std::ostream& os) const {
    os << endl << "StrawPhysics parameters: "  << std::endl
      << "EIonize = ";
    for(auto x : _EIonize ) os << x << " " ;
    os << " MeV " << std::endl;
    os << "Mean free path = " << _meanpath << " mm " << std::endl
      << "Ionization electron kinetic energy = " << _eKin << " MeV" << std::endl;
    os << "intNProb = ";
    for(auto x: _intNProb) os << x << " ";
    os << std::endl;
    os << "Average gas gain = " << _gasgain << std::endl
      << "Polya 'A' parameter = " << _polyaA << std::endl
      << "gslope = " << _gslope << std::endl
      << "nggauss = " << _nggauss << std::endl
      << "Signal propagagation velocity = " << _vprop << " mm/ns" <<  std::endl
      << "Average # of ionization electrons = " << _NAverage << std::endl
      << "Average electron ionization energy = " << _EAverage << " MeV" << std::endl;
    os << "cdpoly = ";
    for(auto x: _cdpoly) os << x << " " ;
    os << std::endl;
    os << "dtvar = ";
    for(auto x: _dtvar) os << x << " " ;
    os << std::endl;
    os << "nonlindrift = " << _nonlindrift << std::endl
      << "bz = " << _bz << std::endl;

  }
}


