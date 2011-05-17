#ifndef TTrackerGeom_ManifoldId_hh
#define TTrackerGeom_ManifoldId_hh

#include <ostream>

#include "TrackerGeom/inc/SectorId.hh"

namespace mu2e {
    
  struct ManifoldId{
    
  public:
    
    ManifoldId():
      _sid(),
      _manifold(-1){
    }
    
    ManifoldId(SectorId sector,
               int manifold
               ):
      _sid(sector),
      _manifold(manifold)
    {}
    
    ~ManifoldId(){}
    
    const SectorId getSectorId() const {
      return _sid;
    }
    
    const int getManifold() const {
      return _manifold;
    }
    
    bool operator== (const ManifoldId s) const{
      return (_sid == s._sid && _manifold == s._manifold);
    }
    
  private:
    
    // Member variables
    SectorId _sid;
    int _manifold;
    
  };
  
  inline std::ostream& operator<<(std::ostream& ost,
                                  const ManifoldId& m ){
    ost << m.getSectorId() << " " << m.getManifold();
    return ost;
  }
  
} // end namespace mu2e

#endif /* TTrackerGeom_ManifoldId_hh */
