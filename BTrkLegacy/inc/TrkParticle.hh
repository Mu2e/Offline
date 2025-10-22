// Class to define particle types which can be reconstructed as tracks.
// the type values are copied from the PDG codes.
#ifndef TrkParticle_hh
#define TrkParticle_hh
#include <string>
// for now just the 5 stable 'fundamental' particles
//  plus deuterons, tritons, alphas.  It should include nuclear fragments, Cascade, someday

class TrkParticle {
  public:
  enum type {
        e_minus = 11 ,
        e_plus = -11 ,
        mu_minus = 13 ,
        mu_plus = -13 ,
        pi_plus = 211 ,
        pi_minus = -211 ,
        K_plus = 321 ,
        K_minus = -321 ,
        p_plus = 2212 ,
        anti_p_minus = -2212,
	deuterium = 1000010020,
	tritium = 1000010030,
	He3 = 1000020030,
	He4 = 100002004
  };
// construct from a type
  TrkParticle(type ptype=e_minus);
  TrkParticle(TrkParticle const& other);
  ~TrkParticle();
// basic accessor
  type particleType() const { return _type; }
  bool operator == (TrkParticle const& other) const { return _type == other._type; }
  bool operator != (TrkParticle const& other) const { return ! this->operator==(other); }
  TrkParticle & operator =(TrkParticle const& other);
// return basic information
  double mass() const;
  double charge() const;
  std::string const& name() const;
// basic kinematics; provide momentum in CLHEP units
  double energy(double momentum) const; // return value in CLHEP units
  double beta(double momentum) const;
  double betagamma(double momentum) const;
  double gamma(double momentum) const;
  private:
  type _type;
};

#endif