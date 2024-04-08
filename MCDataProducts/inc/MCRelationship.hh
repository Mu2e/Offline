//
// MC relationships of some objects
// Dave Brown, LBNL, 20 Jun 2016
//
#ifndef MCRelationship_HH
#define MCRelationship_HH
#include <cstdint>
namespace art { template <typename product> class Ptr; }
namespace mu2e {
  class SimParticle;
  class StrawDigiMC;
  class MCRelationship {
    public:
    typedef art::Ptr<SimParticle> SPPtr;
    // values of relationship of 2 MC objects
    enum relation : int8_t {none=-1,same=0,daughter,mother,sibling,udaughter,umother,usibling};
    relation relationship() const { return (relation)_rel; }
    int8_t removal() const { return _rem; } // relationship generational distance
    // convenience operators
    bool operator ==(MCRelationship const& other ) const { return _rel == other._rel; }
    bool operator !=(MCRelationship const& other ) const { return _rel != other._rel; }
    bool operator < (MCRelationship const& other ) const {
    // note the sign flip: smaller removal = more related, to follow conventional meaning
      if(_rel != none && other._rel != none){
        return _rem > other._rem;
      } else if(_rel != none)
        return false;
      else
        return true;
    }
    bool operator > (MCRelationship const& other ) const {
      return ((!operator ==(other)) && (!operator < (other))); }
    bool operator ==(relation rval) const { return _rel == rval; }
    bool operator !=(relation rval) const { return _rel != rval; }
    // trivial constructor
    MCRelationship(relation rval=none) : _rel(rval), _rem(-1) {}
    // construct from SimParticles
    MCRelationship(SPPtr const& sppi,SPPtr const& sppj);
    // construct from StrawDigiMCs
    MCRelationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);
    // from the mixture
    MCRelationship(StrawDigiMC const& mcd, SPPtr const& spp);
    private:
    int8_t _rel = none; // relationship
    int8_t _rem = -1; // distance between relationship
  };
}
#endif
