#ifndef TrkDiag_SPCount_HH
#define TrkDiag_SPCount_HH
#include "Offline/MCDataProducts/inc/SimParticle.hh"

// MC matching
namespace mu2e {
  struct SPCount {
    SPCount() : _count(0) {}
    SPCount(art::Ptr<SimParticle> const& spp) : _spp(spp), _count(1) {}
    void append(art::Ptr<SimParticle> const& sp) { if(sp == _spp)++_count; }
    bool operator ==(art::Ptr<SimParticle> const& sp) const { return _spp == sp; }
    art::Ptr<SimParticle> _spp;
    unsigned _count;
  };
  struct SPCountComp : public binary_function<SPCount, SPCount , bool> {
    bool operator() (SPCount a, SPCount b) { return a._count > b._count; }
  };
}
#endif
