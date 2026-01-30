// particle description used for tracking
// the implementation depends on the mu2e conditions service, but
// not the interface
#include "Offline/BTrkLegacy/TrkParticle.hh"
#include "Offline/BTrkLegacy/inc/ExternalInfo.hh"
#include "Offline/BTrkLegacy/inc/ParticleInfoInterface.hh"

#include <cmath>

TrkParticle::TrkParticle(TrkParticle::type ptype) : _type(ptype)
{}

TrkParticle::TrkParticle(TrkParticle const& other) : _type(other._type)
{}

TrkParticle::~TrkParticle() {}

TrkParticle&
TrkParticle::operator =(TrkParticle const& other) {
  if(this != &other){
    _type = other._type;
  }
  return *this;
}


double
TrkParticle::mass() const {
// avoid calling the Particle Data table on each call by buffering results in statics
  double retval(-1.0);
  switch (_type) {
    case e_minus: case e_plus: {
      static double e_mass = ExternalInfo::particleInfoInstance()->mass( e_plus );
      return e_mass;
    }

    case mu_minus: case mu_plus: {
      static double mu_mass = ExternalInfo::particleInfoInstance()->mass( mu_plus );
      return mu_mass;
    }

    case pi_minus: case pi_plus: {
      static double pi_mass = ExternalInfo::particleInfoInstance()->mass( pi_plus );
      return pi_mass;
    }

    case K_minus: case K_plus: {
      static double K_mass = ExternalInfo::particleInfoInstance()->mass( K_plus );
      return K_mass;
    }

    case anti_p_minus: case p_plus: {
      static double p_mass = ExternalInfo::particleInfoInstance()->mass( p_plus );
      return p_mass;
    }

    default: {
    }
  }
  return retval;
}

double
TrkParticle::charge() const {
// avoid calling the Particle Data table on each call by buffering results in statics
  double retval(0.0);
  switch (_type) {
    case e_minus: case mu_minus: case pi_minus: case K_minus: case anti_p_minus: {
      static double minus_charge = ExternalInfo::particleInfoInstance()->charge( e_minus );
      return minus_charge;
    }

    case e_plus: case mu_plus: case pi_plus: case K_plus: case p_plus: {
      static double plus_charge = ExternalInfo::particleInfoInstance()->charge( e_plus );
      return plus_charge;
    }

    default: {
    }
  }
  return retval;
}


std::string const&
TrkParticle::name() const {
// I can't use HepPDT for the names as those include Latex characters and so conflict
// with many standard software lexicons (like root).
  switch (_type) {
    case e_minus: {
      static std::string e_minus_name = "eMinus";
      return e_minus_name;
    }
    case e_plus: {
      static std::string e_plus_name = "ePlus";
      return e_plus_name;
    }
    case mu_minus: {
      static std::string mu_minus_name = "muMinus";
      return mu_minus_name;
    }
    case mu_plus: {
      static std::string mu_plus_name = "muPlus";
      return mu_plus_name;
    }
    case pi_minus: {
      static std::string pi_minus_name = "piMinus";
      return pi_minus_name;
    }
    case pi_plus: {
      static std::string pi_plus_name = "piPlus";
      return pi_plus_name;
    }
    case K_minus: {
      static std::string K_minus_name = "KMinus";
      return K_minus_name;
    }
    case K_plus: {
      static std::string K_plus_name = "KPlus";
      return K_plus_name;
    }
    case anti_p_minus: {
      static std::string anti_p_minus_name = "pMinus";
      return anti_p_minus_name;
    }
    case p_plus: {
      static std::string p_plus_name = "pPlus";
      return p_plus_name;
    }
    default: {
      static std::string unknown_name("unknown");
      return unknown_name;
    }
  }
}

double
TrkParticle::energy(double momentum) const {
  double m = mass();
  return sqrt(momentum*momentum + m*m);
}

double
TrkParticle::beta(double momentum) const {
  return fabs(momentum)/energy(momentum);
}

double
TrkParticle::betagamma(double momentum) const {
  return fabs(momentum)/mass();
}

double
TrkParticle::gamma(double momentum) const {
  return energy(momentum)/mass();
}
