
#ifndef ValGenParticle_HH_
#define ValGenParticle_HH_

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Validation/inc/ValId.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValGenParticle {
 public:
  ValGenParticle(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const GenParticleCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  ValId _id;
  TH1D* _hp;
  TH1D* _hlogp;
  TH1D* _hx;
  TH1D* _hxt;
  TH1D* _hy;
  TH1D* _hyt;
  TH1D* _hz;
  TH1D* _hzt;
  TH1D* _t;
  TH1D* _t2;
};
}  // namespace mu2e

#endif
